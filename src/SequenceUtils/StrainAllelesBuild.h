// StrainAllelesBuild.h - how --build makes strain_alleles.tsv (StrainAlleles.h) from the full reference: every copy of a
// database species' gene from a genome that may give alleles (AlleleGenome, --allele_genome_share) is aligned against
// the representative's copy (WFA2, protal's scores, both ends partly free: gene calls start and end apart) and turned
// into edits; of a copy's distinct alleles, the Settings::candidates with the least hashes are kept (a sample that does
// not depend on the order the threads read them in), and of those Settings::alleles chosen greedily for the coverage of
// the sample (Select): each sampled allele stands for the strains near it, and the chosen ones together explain as much
// of the sample's edits, allele by allele, as a read of each would be scored with (Coverage). Until 2026-10-09 the
// choice was farthest first from the representative and each other (k-center), which kept the deep lineages and
// dropped a strain's near relatives: in a world grown on known trees, a strain's nearest allele genome was stored in
// half of its genes with eight allele genomes against nine tenths with two, and the models called the strains the
// alleles reach at 0.98 and the others at 0.3-0.6 (docs/claude/2026-10-09-ancestry-true-positive-test). An allele must
// stay on its species' side of the nearest congener: nearer the representative than the copy's nearest congener's copy
// is (congener_gaps.tsv), and within kMaxDivergence; one that is not is a misassigned genome or a congener's gene, whose
// reads the species would otherwise take.
#pragma once

#include "CongenerGapsTable.h"
#include "StrainAlleles.h"
#include "WFA2Wrapper2.h"

#include <algorithm>
#include <cstdint>
#include <mutex>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace protal::strain_alleles {
    inline constexpr double kFreeEnds = 0.15;      // of each copy's length, free at either end
    inline constexpr double kMinCover = 0.5;       // an allele covers at least this share of the representative's copy
    inline constexpr double kMaxDivergence = 0.1;  // and at most this share of its columns differ
    inline constexpr size_t kEndMatches = 10;      // an allele begins and ends with this many matches (Diff)

    struct Settings {
        size_t alleles = 4;      // kept per copy (--strain_alleles; 0: none)
        size_t candidates = 16;  // distinct alleles sampled per copy, of the least hashes, to choose from
        double share = 1.0;      // of the genomes, those that may give alleles (--allele_genome_share)
    };

    // Counts of a build's pass, for its log line.
    struct Stats {
        size_t records = 0;           // full-reference copies of the database's species and genes
        size_t outside_share = 0;     // of them from genomes that may not give alleles
        size_t identical = 0;         // identical to the representative's copy
        size_t duplicates = 0;        // a copy's sequence seen before
        size_t unaligned = 0;         // covering too little of the representative, too divergent, or past the aligner
        size_t beyond_congener = 0;   // as far from the representative as the nearest congener's copy, or farther
        size_t offered = 0;           // alleles offered to their copy's sample
        void Add(Stats const& o) {
            records += o.records;
            outside_share += o.outside_share;
            identical += o.identical;
            duplicates += o.duplicates;
            unaligned += o.unaligned;
            beyond_congener += o.beyond_congener;
            offered += o.offered;
        }
    };

    // The allele of `copy` against the representative's `rep`: its edits where it covers rep (between the alignment's
    // first and last runs of kEndMatches matches), nullopt if it covers less than kMinCover of rep, if more than
    // kMaxDivergence of its columns differ, or if the aligner gives up. `divergence` is that share; `ops` is scratch. A
    // base of the copy that is not ACGT is no edit (unknown).
    inline std::optional<Allele> Diff(std::string const& rep, std::string const& copy, WFA2Wrapper2& aligner, std::string& ops,
                                      double& divergence) {
        divergence = 0;
        if (rep.empty() || copy.empty() || rep.size() > UINT16_MAX) return std::nullopt;
        int const free_copy = static_cast<int>(kFreeEnds * static_cast<double>(copy.size()));
        int const free_rep = static_cast<int>(kFreeEnds * static_cast<double>(rep.size()));
        size_t const longer = std::max(rep.size(), copy.size());
        int const max_score = static_cast<int>(4.5 * kMaxDivergence * static_cast<double>(longer)) + 64;
        aligner.Alignment(copy, rep, free_copy, free_copy, free_rep, free_rep, max_score);
        if (!aligner.Success()) return std::nullopt;
        aligner.CigarInto(ops);
        // Each end trimmed to the first run of kEndMatches matches: ends-free alignment skips a prefix of one copy only,
        // so a copy that starts elsewhere *and* has a flank of its own gets the flank aligned as mismatches and gaps.
        size_t first = 0, last = ops.size();
        uint32_t ref = 0, query = 0;
        auto matches_from = [&ops](size_t i, size_t end) {
            if (i + kEndMatches > end) return false;
            for (size_t k = i; k < i + kEndMatches; k++) {
                if (ops[k] != 'M' && ops[k] != '=') return false;
            }
            return true;
        };
        while (first < last && !matches_from(first, last)) {
            ref += ops[first] != 'I';
            query += ops[first] != 'D';
            first++;
        }
        while (last > first && !matches_from(last - kEndMatches, last)) last--;
        Allele allele;
        allele.begin = static_cast<uint16_t>(ref);
        size_t columns = 0, differing = 0;
        for (size_t i = first; i < last;) {
            char const op = ops[i];
            size_t j = i + 1;
            while (j < last && ops[j] == op) j++;
            uint32_t const run = static_cast<uint32_t>(j - i);
            columns += run;
            switch (op) {
                case 'X':
                    for (uint32_t k = 0; k < run; k++) {
                        uint8_t const base = query + k < copy.size() ? BaseCode(copy[query + k]) : 4;
                        if (base <= 3) allele.edits.push_back(Edit::Substitution(static_cast<uint16_t>(ref + k), base));
                    }
                    differing += run;
                    ref += run;
                    query += run;
                    break;
                case 'I':
                    allele.edits.push_back(Edit::Indel(kInsertion, static_cast<uint16_t>(ref), static_cast<uint16_t>(std::min<uint32_t>(run, Edit::kMaxLength))));
                    differing += run;
                    query += run;
                    break;
                case 'D':
                    allele.edits.push_back(Edit::Indel(kDeletion, static_cast<uint16_t>(ref), static_cast<uint16_t>(std::min<uint32_t>(run, Edit::kMaxLength))));
                    differing += run;
                    ref += run;
                    break;
                default:  // M, =
                    ref += run;
                    query += run;
                    break;
            }
            i = j;
        }
        allele.end = static_cast<uint16_t>(std::min<uint32_t>(ref, UINT16_MAX));
        // By position, then kind: an insertion before a base and that base's substitution share a position, and the
        // alignment gives the insertion first (Edit::operator< puts the substitution first, as the table keeps them).
        std::sort(allele.edits.begin(), allele.edits.end());
        if (columns == 0 || allele.end <= allele.begin) return std::nullopt;
        divergence = static_cast<double>(differing) / static_cast<double>(columns);
        if (static_cast<double>(allele.end - allele.begin) < kMinCover * static_cast<double>(rep.size()) || divergence > kMaxDivergence) {
            return std::nullopt;
        }
        return allele;
    }

    // Edits two alleles do not share (the size of their symmetric difference), and an allele's distance from the
    // representative (its edits).
    inline size_t EditDistance(Allele const& a, Allele const& b) {
        size_t i = 0, j = 0, shared = 0;
        while (i < a.edits.size() && j < b.edits.size()) {
            if (a.edits[i] == b.edits[j]) {
                shared++;
                i++;
                j++;
            } else if (a.edits[i] < b.edits[j]) {
                i++;
            } else {
                j++;
            }
        }
        return a.edits.size() + b.edits.size() - 2 * shared;
    }

    // Edits two alleles share (both lists sorted).
    inline size_t SharedEdits(Allele const& a, Allele const& b) {
        return (a.edits.size() + b.edits.size() - EditDistance(a, b)) / 2;
    }

    // How much of the allele `a` the allele `c`, stored, would save a read of a's genome that spans the gene, as the
    // alignment scores count it (BestAllele: the edits of a that c shares, less the edits of c that a lacks, when
    // positive), as a share of a's edits in units of kCoverageUnit: the unit for c = a, 0 when c shares no more than
    // half of its own edits with a (the read keeps the representative) or when a has no edit.
    inline constexpr uint64_t kCoverageUnit = 1u << 20;

    inline uint64_t Coverage(Allele const& a, Allele const& c) {
        size_t const shared = SharedEdits(a, c);
        if (a.edits.empty() || 2 * shared <= c.edits.size()) return 0;
        return kCoverageUnit * static_cast<uint64_t>(2 * shared - c.edits.size()) / a.edits.size();
    }

    // Of a copy's sampled alleles (hash, allele), up to k chosen greedily for the sample's coverage: the sum over the
    // sampled alleles of the best Coverage a chosen allele gives each (every sampled allele stands for the strains near
    // it, with the same weight whatever its distance from the representative). Each step takes the allele that adds
    // most (of equals the least hash), until k or none adds any (one identical to a chosen allele). Integer units, so
    // that equal gains are equal on every machine.
    inline std::vector<Allele> Select(std::vector<std::pair<uint64_t, Allele>> sample, size_t k) {
        std::sort(sample.begin(), sample.end(), [](auto const& a, auto const& b) { return a.first < b.first; });
        size_t const n = sample.size();
        std::vector<uint64_t> coverage(n * n);  // of allele i by allele j stored
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) coverage[i * n + j] = Coverage(sample[i].second, sample[j].second);
        }
        std::vector<uint64_t> covered(n, 0);
        std::vector<uint8_t> taken(n, 0);
        std::vector<Allele> chosen;
        while (chosen.size() < k) {
            size_t best = n;
            uint64_t best_gain = 0;
            for (size_t j = 0; j < n; j++) {
                if (taken[j]) continue;
                uint64_t gain = 0;
                for (size_t i = 0; i < n; i++) {
                    uint64_t const c = coverage[i * n + j];
                    if (c > covered[i]) gain += c - covered[i];
                }
                if (gain > best_gain) {
                    best = j;
                    best_gain = gain;
                }
            }
            if (best == n) break;
            taken[best] = 1;
            chosen.push_back(sample[best].second);
            for (size_t i = 0; i < n; i++) covered[i] = std::max(covered[i], coverage[i * n + best]);
        }
        return chosen;
    }

    // The alleles offered per copy, from any thread: per copy the seen sequences (by hash, up to kSeen, so that the many
    // identical copies of a big species are aligned once) and the sample of the least hashes of its distinct alleles.
    class Collector {
    public:
        static constexpr size_t kShards = 1024;
        static constexpr size_t kSeen = 64;

        explicit Collector(size_t candidates) : m_candidates(std::max<size_t>(candidates, 1)), m_shards(kShards) {}

        static uint64_t Key(uint32_t taxid, uint32_t gene) { return static_cast<uint64_t>(taxid) << 32 | gene; }

        // Whether a copy's sequence (its hash) was offered before; records it if not.
        bool Seen(uint64_t key, uint64_t sequence_hash) {
            auto& shard = ShardOf(key);
            std::lock_guard<std::mutex> lock(shard.mutex);
            auto& slot = shard.slots[key];
            if (std::find(slot.seen.begin(), slot.seen.end(), sequence_hash) != slot.seen.end()) return true;
            if (slot.seen.size() < kSeen) slot.seen.push_back(sequence_hash);
            return false;
        }

        void Offer(uint64_t key, Allele allele) {
            uint64_t const hash = Fnv1a(AlleleText(allele));
            auto& shard = ShardOf(key);
            std::lock_guard<std::mutex> lock(shard.mutex);
            auto& kept = shard.slots[key].kept;
            for (auto const& [h, _] : kept) {
                if (h == hash) return;
            }
            if (kept.size() < m_candidates) {
                kept.emplace_back(hash, std::move(allele));
                return;
            }
            auto const largest = std::max_element(kept.begin(), kept.end(), [](auto const& a, auto const& b) { return a.first < b.first; });
            if (hash < largest->first) *largest = { hash, std::move(allele) };
        }

        // Rows for Table::FromRows: per copy with a sampled allele, the Select of k.
        std::vector<std::tuple<uint32_t, uint32_t, std::vector<Allele>>> Rows(size_t k) {
            std::vector<std::tuple<uint32_t, uint32_t, std::vector<Allele>>> rows;
            for (auto& shard : m_shards) {
                for (auto& [key, slot] : shard.slots) {
                    if (slot.kept.empty()) continue;
                    auto chosen = Select(std::move(slot.kept), k);
                    if (!chosen.empty()) rows.emplace_back(static_cast<uint32_t>(key >> 32), static_cast<uint32_t>(key & 0xffffffffu), std::move(chosen));
                }
                shard.slots.clear();
            }
            return rows;
        }

    private:
        struct Slot {
            std::vector<uint64_t> seen;
            std::vector<std::pair<uint64_t, Allele>> kept;
        };
        struct Shard {
            std::mutex mutex;
            std::unordered_map<uint64_t, Slot> slots;
        };
        Shard& ShardOf(uint64_t key) {
            uint64_t x = key + 0x9E3779B97F4A7C15ull;
            x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ull;
            x = (x ^ (x >> 27)) * 0x94D049BB133111EBull;
            return m_shards[(x ^ (x >> 31)) % kShards];
        }
        size_t m_candidates;
        std::vector<Shard> m_shards;
    };

    // One full-reference record of a database species' gene: its allele offered to `collector`, counted in `stats`.
    // `rep` is the representative's copy, `accession` the record's genome ("" if the full reference names none), `gap`
    // the copy's nearest congener's distance (-1: none known; 0, an identical congener's copy, admits no allele).
    // Thread-safe through the collector; `aligner` and the strings are the thread's.
    inline void AddRecord(Collector& collector, Stats& stats, Settings const& settings, uint32_t taxid, uint32_t gene,
                          std::string_view accession, std::string const& rep, std::string const& copy, double gap,
                          WFA2Wrapper2& aligner, std::string& ops) {
        stats.records++;
        if (!AlleleGenome(accession, settings.share)) {
            stats.outside_share++;
            return;
        }
        if (copy == rep) {
            stats.identical++;
            return;
        }
        uint64_t const key = Collector::Key(taxid, gene);
        if (collector.Seen(key, Fnv1a(copy))) {
            stats.duplicates++;
            return;
        }
        double divergence = 0;
        auto allele = Diff(rep, copy, aligner, ops, divergence);
        if (!allele || allele->edits.empty()) {
            (allele ? stats.identical : stats.unaligned)++;
            return;
        }
        if (gap >= 0 && divergence >= gap) {
            stats.beyond_congener++;
            return;
        }
        stats.offered++;
        collector.Offer(key, std::move(*allele));
    }
}
