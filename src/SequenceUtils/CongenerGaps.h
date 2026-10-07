// CongenerGaps.h - how far each species' copy of a marker gene lies from its congeners' copies, by alignment: per copy
// the distance to its nearest congener's copy (min) and to a typical one (median), as congener_gaps.tsv in the
// database. --build aligns each copy against the copies of its genus (all of them up to Settings::all, else its
// Settings::nearest nearest by their k-mer sketches and Settings::sample others drawn by a hash of the pair, the sample
// alone giving the median), with WFA2 and protal's read scores, both ends partly free (gene calls of congeners start and
// end in different places): the distance is the share of the aligned columns that differ (mismatches and gap bases),
// as a read's divergence is counted.
//
// A run compares each read's divergence on a copy with that copy's gap (Profiler.h, the "gaps" features): a strain of
// the species stays within the gap to its nearest congener, while a congener the database lacks lies about as far from
// the reference as its congeners do. The run's relative_* features compare a taxon only with the congeners present in
// the sample; this table holds every congener of the database (docs/claude/2026-10-07-error-read-signatures, section 6
// and the follow-up). A read covers part of a gene and the gene's divergence varies along it, so one read against the
// whole gene's gap is a noisy test; the features count over a taxon's reads.
#pragma once

#include "CongenerGapsTable.h"
#include "GeneConservation.h"
#include "WFA2Wrapper2.h"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace protal::congener_gaps {
    inline constexpr double kMaxDistance = 0.4;          // an alignment that would reach this is given up: "at least"
    inline constexpr double kFreeEnds = 0.15;            // of each copy's length, free at either end
    inline constexpr double kMinOverlap = 0.5;           // aligned columns needed, of the shorter copy's length
    inline constexpr size_t kSketchSize = 128;           // bottom sketch hashes for choosing the nearest congeners

    // Which congeners a copy is aligned against.
    struct Settings {
        size_t all = 24;      // a genus with at most this many other copies of the gene: all of them
        size_t nearest = 4;   // else the nearest by sketch distance (for the min)
        size_t sample = 16;   // and this many drawn by hash (the median's sample)
    };

    // The distance of two gene copies from an alignment with both ends partly free: of the aligned columns between
    // the first and last aligned base pair, the share that are mismatches or gap bases (X, I and D of WFA2's
    // operations). NaN when the alignment covers less than kMinOverlap of the shorter copy; kMaxDistance when the
    // copies are farther apart than that (the aligner gives up). `cigar` is scratch.
    inline double AlignedDistance(std::string const& a, std::string const& b, WFA2Wrapper2& aligner, std::string& cigar) {
        if (a.empty() || b.empty()) return std::nan("");
        if (a == b) return 0;
        int const free_a = static_cast<int>(kFreeEnds * static_cast<double>(a.size()));
        int const free_b = static_cast<int>(kFreeEnds * static_cast<double>(b.size()));
        // A mismatch costs 4 (protal's scores, 4/6/2): kMaxDistance of the longer copy's bases, as mismatches, and a
        // margin for the gaps that go with them.
        size_t const longer = std::max(a.size(), b.size());
        int const max_score = static_cast<int>(4.5 * kMaxDistance * static_cast<double>(longer)) + 64;
        aligner.Alignment(b, a, free_b, free_b, free_a, free_a, max_score);
        if (!aligner.Success()) return kMaxDistance;
        aligner.CigarInto(cigar);
        size_t begin = 0, end = cigar.size();
        while (begin < end && (cigar[begin] == 'I' || cigar[begin] == 'D')) begin++;
        while (end > begin && (cigar[end - 1] == 'I' || cigar[end - 1] == 'D')) end--;
        size_t matches = 0, differences = 0, a_bases = 0;
        for (size_t i = begin; i < end; i++) {
            char const op = cigar[i];
            if (op == 'M' || op == '=') matches++;
            else differences++;
            a_bases += op != 'I';
        }
        size_t const columns = matches + differences;
        size_t const shorter = std::min(a.size(), b.size());
        if (columns == 0 || static_cast<double>(std::min(a_bases, columns)) < kMinOverlap * static_cast<double>(shorter)) {
            return std::nan("");
        }
        return std::min(kMaxDistance, static_cast<double>(differences) / static_cast<double>(columns));
    }

    // Counts of a build's scan, for its log line.
    struct Stats {
        size_t genes = 0;           // genes with a genus of two copies or more
        size_t copies = 0;          // copies given a gap
        size_t alignments = 0;      // pairs aligned
        size_t far = 0;             // of them given up at kMaxDistance
        size_t unaligned = 0;       // or without enough overlap
        size_t sampled_genera = 0;  // genus-gene groups above Settings::all
    };

    namespace detail {
        inline uint64_t Mix(uint64_t x) {
            x += 0x9E3779B97F4A7C15ull;
            x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ull;
            x = (x ^ (x >> 27)) * 0x94D049BB133111EBull;
            return x ^ (x >> 31);
        }
        // A pair's draw for the median's sample: the same for (a, b) and (b, a), different per gene.
        inline uint64_t PairHash(uint32_t a, uint32_t b, uint32_t gene) {
            uint64_t const lo = std::min(a, b), hi = std::max(a, b);
            return Mix((lo << 32 | hi) ^ Mix(gene));
        }
        inline double Median(std::vector<double>& v) {
            std::sort(v.begin(), v.end());
            size_t const n = v.size();
            return n % 2 ? v[n / 2] : (v[n / 2 - 1] + v[n / 2]) / 2;
        }
    }

    // The gaps of one gene's copies (taxids[i] holds seqs[i]; one copy per taxid): each copy compared with the copies
    // of its genus (genus_of(taxid), 0 for none) as Settings says, on `threads` threads. Adds (taxid, Gap) for every
    // copy with at least one congener's copy aligned.
    template<typename GenusOf>
    void ScanGene(uint32_t geneid, std::vector<uint32_t> const& taxids, std::vector<std::string> const& seqs,
                  GenusOf&& genus_of, Settings const& settings, int threads,
                  std::vector<std::pair<uint32_t, Gap>>& out, Stats& stats) {
        threads = std::max(threads, 1);
        std::map<uint32_t, std::vector<uint32_t>> by_genus;  // genus -> copy indices (by taxid: taxids is sorted)
        for (size_t i = 0; i < taxids.size(); i++) {
            if (uint32_t const genus = genus_of(taxids[i])) by_genus[genus].push_back(static_cast<uint32_t>(i));
        }
        std::vector<std::pair<uint32_t, uint32_t>> pairs;            // copy indices, lower first
        std::vector<std::vector<uint32_t>> median_of(taxids.size());  // a sampled copy's median partners
        std::vector<uint8_t> sampled(taxids.size(), 0);
        std::vector<uint16_t> congeners(taxids.size(), 0);
        bool any = false;
        for (auto const& [genus, members] : by_genus) {
            size_t const n = members.size();
            if (n < 2) continue;
            any = true;
            for (uint32_t const i : members) congeners[i] = static_cast<uint16_t>(std::min<size_t>(n - 1, UINT16_MAX));
            if (n - 1 <= settings.all) {
                for (size_t x = 0; x < n; x++) {
                    for (size_t y = x + 1; y < n; y++) pairs.emplace_back(members[x], members[y]);
                }
                continue;
            }
            stats.sampled_genera++;
            std::vector<std::vector<uint32_t>> sketches(n);
            #pragma omp parallel for schedule(dynamic, 16) num_threads(threads)
            for (int64_t x = 0; x < static_cast<int64_t>(n); x++) {
                sketches[static_cast<size_t>(x)] = gene_conservation::BottomSketch(seqs[members[static_cast<size_t>(x)]], kSketchSize);
            }
            std::vector<std::vector<std::pair<uint32_t, uint32_t>>> chosen(n);
            #pragma omp parallel for schedule(dynamic, 8) num_threads(threads)
            for (int64_t x = 0; x < static_cast<int64_t>(n); x++) {
                size_t const xi = static_cast<size_t>(x);
                std::vector<std::pair<double, uint32_t>> by_distance;
                std::vector<std::pair<uint64_t, uint32_t>> by_hash;
                by_distance.reserve(n - 1);
                by_hash.reserve(n - 1);
                for (size_t y = 0; y < n; y++) {
                    if (y == xi) continue;
                    by_distance.emplace_back(gene_conservation::SketchDistance(sketches[xi], sketches[y]), members[y]);
                    by_hash.emplace_back(detail::PairHash(taxids[members[xi]], taxids[members[y]], geneid), members[y]);
                }
                size_t const k = std::min(settings.nearest, by_distance.size());
                std::partial_sort(by_distance.begin(), by_distance.begin() + static_cast<std::ptrdiff_t>(k), by_distance.end());
                size_t const s = std::min(settings.sample, by_hash.size());
                std::partial_sort(by_hash.begin(), by_hash.begin() + static_cast<std::ptrdiff_t>(s), by_hash.end());
                uint32_t const me = members[xi];
                std::vector<uint32_t> partners;
                for (size_t q = 0; q < s; q++) partners.push_back(by_hash[q].second);
                std::sort(partners.begin(), partners.end());
                median_of[me] = partners;
                for (size_t q = 0; q < k; q++) chosen[xi].emplace_back(std::min(me, by_distance[q].second), std::max(me, by_distance[q].second));
                for (uint32_t const p : partners) chosen[xi].emplace_back(std::min(me, p), std::max(me, p));
            }
            for (size_t x = 0; x < n; x++) {
                sampled[members[x]] = 1;
                pairs.insert(pairs.end(), chosen[x].begin(), chosen[x].end());
            }
        }
        if (!any) return;
        stats.genes++;
        std::sort(pairs.begin(), pairs.end());
        pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
        std::vector<double> distance(pairs.size());
        #pragma omp parallel num_threads(threads)
        {
            WFA2Wrapper2 aligner(4, 6, 2, 0);
            std::string cigar;
            #pragma omp for schedule(dynamic, 64)
            for (int64_t p = 0; p < static_cast<int64_t>(pairs.size()); p++) {
                auto const [i, j] = pairs[static_cast<size_t>(p)];
                distance[static_cast<size_t>(p)] = AlignedDistance(seqs[i], seqs[j], aligner, cigar);
            }
        }
        stats.alignments += pairs.size();
        // Every copy's aligned partners and distances, by partner.
        std::vector<std::vector<std::pair<uint32_t, double>>> found(taxids.size());
        for (size_t p = 0; p < pairs.size(); p++) {
            double const d = distance[p];
            if (std::isnan(d)) {
                stats.unaligned++;
                continue;
            }
            stats.far += d >= kMaxDistance;
            auto const [i, j] = pairs[p];
            found[i].emplace_back(j, d);
            found[j].emplace_back(i, d);
        }
        std::vector<double> values;
        for (size_t i = 0; i < taxids.size(); i++) {
            auto& list = found[i];
            if (list.empty()) continue;
            std::sort(list.begin(), list.end());
            double least = std::numeric_limits<double>::infinity();
            for (auto const& [_, d] : list) least = std::min(least, d);
            values.clear();
            if (sampled[i]) {
                for (auto const& [partner, d] : list) {
                    if (std::binary_search(median_of[i].begin(), median_of[i].end(), partner)) values.push_back(d);
                }
            } else {
                for (auto const& [_, d] : list) values.push_back(d);
            }
            if (values.empty()) continue;  // only nearest congeners aligned: no median sample
            out.emplace_back(taxids[i], Gap{ Scaled(least), Scaled(detail::Median(values)), congeners[i] });
            stats.copies++;
        }
    }
}
