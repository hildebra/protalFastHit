// AncestrySites.h - where a species' copy of a marker gene differs from its congeners' copies, and which side a read
// takes there. The sites are the species' own derived states: the positions where the congeners compared agree on one
// base (kConsensus of those compared at the position) and the species carries another. A strain of the species carries
// the species' base there; a species the database lacks that branched off the lineage below a site carries the
// congeners'. The reads of a false positive are mostly a novel congener's, at the same identity as a missed strain's,
// and this is the one place in the read where the two differ (docs/claude/2026-10-07-error-read-signatures). The
// index's unique k-mers see part of it (lu_per_kb); counting the sites with the reads' bases does not lose a site to
// every sequencing error around it.
//
// 0.7.9 compared the nearest congener alone. Half of the differences between two species are the congener's own
// derived states, where any third species agrees with the species, so a novel congener agreed at 0.68 instead of near
// 0 at r226 v17, and the residual errors were where that signal inverts (docs/claude/2026-10-08-r226-v17). The
// consensus over the congeners (one vote per species, kConsensus of those compared at the position, kMinCongeners of
// them at least) leaves the species' own derived sites; a large genus, whose nearest congener shares most of the
// species' history, gets as many as a small one; with fewer congeners the sites are the nearest's differences, as
// before.
//
// Two copies are compared along their shared 12-mers, without an alignment: the 12-mers unique to each copy pair them
// on their main diagonal, and the bases are compared on the stretches between two paired 12-mers that lie on it (an
// indel moves the diagonal and ends the stretch, so a stretch past it is not compared; a position no congener was
// compared at is no site). Per copy the sites are a sorted list of positions with the congeners' base: a few hundred
// bytes. They are computed once per run for each (species, gene) a record touches (Cache: the gene store holds every
// reference's copy; the congeners are the gene's nearest by alignment in congener_gaps.tsv and the species' nearest in
// species_neighbours.tsv, up to kCongeners) and counted per record from its CIGAR (M: the species' base at every
// site; X: the read's base looked up). The gaps features (CongenerGaps.h) take their nearest congener from the same
// table.
#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "CongenerGapsTable.h"
#include "SpeciesNeighbours.h"

namespace protal::ancestry {
    inline constexpr size_t kKmer = 12;            // the copies are paired by their shared unique 12-mers
    inline constexpr size_t kCongeners = 17;       // congeners compared at most: the gene's nearest by alignment and the
                                                   // species' 16 nearest (species_neighbours::kMaxNeighbours)
    inline constexpr double kConsensus = 0.9;      // of the congeners compared at a position, the share that must carry
                                                   // one base other than the species' for a site (three of three, ...,
                                                   // nine of ten)
    inline constexpr size_t kMinCongeners = 3;     // congeners compared at a position for the consensus to apply; with
                                                   // fewer the site is the nearest congener's difference, as in 0.7.9
    inline constexpr double kMinCompared = 0.5;    // of the species' copy compared, else the congener's copy is no use
    inline constexpr size_t kMaxLength = 65535;    // positions are 16-bit
    inline constexpr size_t kMaxCached = 200000;   // (species, gene) pairs kept; the cache starts over beyond
    inline constexpr double kUnknown = -1;         // the features of a taxon without sites

    // A base's code (A 0, C 1, G 2, T 3, any case), 4 for another letter.
    inline uint8_t Code(char c) {
        switch (c) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            default: return 4;
        }
    }

    // The sites of one copy: against one congener's copy (Compare: where the congener differs, with its base), or the
    // consensus over its congeners (Consensus: where kConsensus of those compared carry one other base, with that base).
    struct Sites {
        std::vector<uint16_t> positions;  // on the species' copy, 0-based, ascending
        std::vector<uint8_t> bases;       // the congeners' base at each (Code)
        uint32_t congener = 0;            // the nearest congener compared (the highest identity); 0: none
        uint16_t congeners = 0;           // congeners compared
        size_t compared = 0;              // positions compared with any congener
        float identity = 0;               // of the species' copy and the nearest congener's over the bases compared

        bool Empty() const { return positions.empty(); }

        // The sites in [begin, end) of the copy.
        size_t Count(size_t begin, size_t end) const {
            auto const lo = std::lower_bound(positions.begin(), positions.end(), begin);
            auto const hi = std::lower_bound(positions.begin(), positions.end(), end);
            return static_cast<size_t>(hi - lo);
        }

        // The index of the site at position p, or SIZE_MAX.
        size_t Find(size_t p) const {
            auto const it = std::lower_bound(positions.begin(), positions.end(), p);
            return it != positions.end() && *it == p ? static_cast<size_t>(it - positions.begin()) : SIZE_MAX;
        }
    };

    namespace detail {
        using Kmer = std::pair<uint32_t, uint32_t>;  // (value, position)

        // The 12-mers of a sequence as (value, position), every k-mer of A/C/G/T only.
        inline std::vector<Kmer> Kmers(std::string_view s) {
            std::vector<Kmer> out;
            if (s.size() < kKmer) return out;
            out.reserve(s.size());
            uint32_t value = 0, valid = 0;
            constexpr uint32_t mask = (1u << (2 * kKmer)) - 1;
            for (size_t i = 0; i < s.size(); i++) {
                uint8_t const c = Code(s[i]);
                if (c > 3) {
                    valid = 0;
                    value = 0;
                    continue;
                }
                value = ((value << 2) | c) & mask;
                if (++valid >= kKmer) out.emplace_back(value, static_cast<uint32_t>(i + 1 - kKmer));
            }
            return out;
        }

        // Keeps the k-mers that occur once; sorted by value.
        inline void KeepUnique(std::vector<Kmer>& kmers) {
            std::sort(kmers.begin(), kmers.end());
            std::vector<Kmer> unique;
            unique.reserve(kmers.size());
            for (size_t i = 0; i < kmers.size();) {
                size_t j = i;
                while (j < kmers.size() && kmers[j].first == kmers[i].first) j++;
                if (j == i + 1) unique.push_back(kmers[i]);
                i = j;
            }
            kmers.swap(unique);
        }

        // The unique 12-mers of a sequence, sorted by value: what a copy is compared with others by.
        inline std::vector<Kmer> UniqueKmers(std::string_view s) {
            auto kmers = Kmers(s);
            KeepUnique(kmers);
            return kmers;
        }
    }

    // The comparison of a copy with one congener's: where the congener differs (sites), and which positions of the
    // copy were compared at all (covered: one char per position, 1 when compared).
    struct Comparison {
        Sites sites;
        std::vector<char> covered;
    };

    // Compares `own` (its unique 12-mers in own_kmers, detail::UniqueKmers) with `other`: nullopt if fewer than
    // kMinCompared of own's bases could be compared (the copies do not pair: another gene, or too different), or own
    // is longer than kMaxLength.
    inline std::optional<Comparison> CompareCopies(std::string_view own, std::vector<detail::Kmer> const& own_kmers,
                                                   std::string_view other) {
        if (own.size() > kMaxLength || own.size() < kKmer || other.size() < kKmer) return std::nullopt;
        auto const& a = own_kmers;
        auto const b = detail::UniqueKmers(other);
        // The shared unique k-mers and their diagonals (position on own minus position on other).
        std::vector<std::pair<int64_t, uint32_t>> pairs;  // (diagonal, position on own)
        for (size_t i = 0, j = 0; i < a.size() && j < b.size();) {
            if (a[i].first < b[j].first) i++;
            else if (b[j].first < a[i].first) j++;
            else {
                pairs.emplace_back(static_cast<int64_t>(a[i].second) - static_cast<int64_t>(b[j].second), a[i].second);
                i++;
                j++;
            }
        }
        if (pairs.size() < 2) return std::nullopt;
        std::sort(pairs.begin(), pairs.end());
        // The main diagonal: the one most pairs lie on.
        int64_t diagonal = 0;
        size_t best = 0;
        for (size_t i = 0; i < pairs.size();) {
            size_t j = i;
            while (j < pairs.size() && pairs[j].first == pairs[i].first) j++;
            if (j - i > best) {
                best = j - i;
                diagonal = pairs[i].first;
            }
            i = j;
        }
        // The stretches between consecutive pairs on it are compared (the pairs are sorted by position within a
        // diagonal).
        Comparison c;
        c.covered.assign(own.size(), 0);
        uint32_t previous = 0;
        bool first = true;
        for (auto const& [d, p] : pairs) {
            if (d != diagonal) continue;
            if (!first) std::fill(c.covered.begin() + previous, c.covered.begin() + std::min<size_t>(own.size(), p + kKmer), 1);
            previous = p;
            first = false;
        }
        size_t n = 0;
        for (size_t i = 0; i < own.size(); i++) {
            if (!c.covered[i]) continue;
            int64_t const j = static_cast<int64_t>(i) - diagonal;
            uint8_t const x = Code(own[i]);
            uint8_t const y = j < 0 || j >= static_cast<int64_t>(other.size()) ? 4 : Code(other[static_cast<size_t>(j)]);
            if (x > 3 || y > 3) {
                c.covered[i] = 0;  // another letter on either copy, or beyond the other copy: not compared
                continue;
            }
            n++;
            if (x != y) {
                c.sites.positions.push_back(static_cast<uint16_t>(i));
                c.sites.bases.push_back(y);
            }
        }
        if (static_cast<double>(n) < kMinCompared * static_cast<double>(own.size())) return std::nullopt;
        c.sites.compared = n;
        c.sites.congeners = 1;
        c.sites.identity = static_cast<float>(1.0 - static_cast<double>(c.sites.positions.size()) / static_cast<double>(n));
        return c;
    }

    inline std::optional<Comparison> CompareCopies(std::string_view own, std::string_view other) {
        if (own.size() > kMaxLength || own.size() < kKmer) return std::nullopt;
        return CompareCopies(own, detail::UniqueKmers(own), other);
    }

    // The sites of `own` against one congener's copy `other`: where other differs, with its base (nullopt as
    // CompareCopies).
    inline std::optional<Sites> Compare(std::string_view own, std::string_view other) {
        auto c = CompareCopies(own, other);
        if (!c) return std::nullopt;
        return std::move(c->sites);
    }

    // The consensus sites of a copy of `length` over its comparisons with its congeners (one vote each): at the positions
    // compared with kMinCongeners congeners or more, those where kConsensus of them carry one base other than the
    // species' (three of three, nine of ten); at the positions compared with fewer (a small genus, or a stretch past an
    // indel in the others' copies), those where the nearest congener (the highest identity) differs, with its base, as
    // 0.7.9 had everywhere. congener and identity are the nearest congener's, compared counts the positions any
    // congener was compared at, congeners how many were. Empty without a comparison.
    inline Sites Consensus(size_t length, std::vector<Comparison> const& comparisons) {
        Sites out;
        if (comparisons.empty() || length > kMaxLength) return out;
        auto const nearest = std::max_element(comparisons.begin(), comparisons.end(), [](Comparison const& a, Comparison const& b) {
            return a.sites.identity < b.sites.identity;
        });
        std::vector<uint8_t> compared(length, 0);
        std::vector<std::array<uint8_t, 4>> votes(length, std::array<uint8_t, 4>{ 0, 0, 0, 0 });
        std::vector<int8_t> fallback(length, -1);  // the nearest congener's base where it differs
        for (auto const& c : comparisons) {
            for (size_t i = 0; i < length && i < c.covered.size(); i++) compared[i] += c.covered[i] != 0;
            for (size_t s = 0; s < c.sites.positions.size(); s++) {
                if (c.sites.positions[s] >= length) continue;
                votes[c.sites.positions[s]][c.sites.bases[s] & 3]++;
                if (&c == &*nearest) fallback[c.sites.positions[s]] = static_cast<int8_t>(c.sites.bases[s] & 3);
            }
        }
        for (size_t i = 0; i < length; i++) {
            if (compared[i] == 0) continue;
            out.compared++;
            if (compared[i] < kMinCongeners) {
                if (fallback[i] >= 0) {
                    out.positions.push_back(static_cast<uint16_t>(i));
                    out.bases.push_back(static_cast<uint8_t>(fallback[i]));
                }
                continue;
            }
            uint8_t best = 0;
            for (uint8_t b = 1; b < 4; b++) {
                if (votes[i][b] > votes[i][best]) best = b;
            }
            if (votes[i][best] > 0 && static_cast<double>(votes[i][best]) >= kConsensus * static_cast<double>(compared[i])) {
                out.positions.push_back(static_cast<uint16_t>(i));
                out.bases.push_back(best);
            }
        }
        out.congener = nearest->sites.congener;
        out.identity = nearest->sites.identity;
        out.congeners = static_cast<uint16_t>(std::min<size_t>(comparisons.size(), UINT16_MAX));
        return out;
    }

    // What a record covers of the sites: those it covers (M, = and X ops; a deletion skips its sites), those where
    // the read has the species' base (M and =, and X where the read's base is neither), and those where it has the
    // congeners'. A record without a sequence ("*") counts nothing.
    struct Counts {
        uint64_t sites = 0, agree = 0, congener = 0;
    };

    inline Counts Count(Sites const& sites, std::string const& cigar, size_t pos, std::string const& seq) {
        Counts c;
        if (sites.Empty() || seq.empty() || seq == "*") return c;
        size_t ref = pos == 0 ? 0 : pos - 1, query = 0, run = 0;
        for (char const op : cigar) {
            if (op >= '0' && op <= '9') {
                run = run * 10 + static_cast<size_t>(op - '0');
                continue;
            }
            if (op == 'M' || op == '=') {
                size_t const n = sites.Count(ref, ref + run);
                c.sites += n;
                c.agree += n;
                ref += run;
                query += run;
            } else if (op == 'X') {
                auto it = std::lower_bound(sites.positions.begin(), sites.positions.end(), ref);
                for (; it != sites.positions.end() && *it < ref + run; ++it) {
                    size_t const q = query + (*it - ref);
                    if (q >= seq.size()) break;
                    c.sites++;
                    if (Code(seq[q]) == sites.bases[static_cast<size_t>(it - sites.positions.begin())]) c.congener++;
                }
                ref += run;
                query += run;
            } else if (op == 'D' || op == 'N') {
                ref += run;
            } else if (op == 'I' || op == 'S') {
                query += run;
            }
            run = 0;
        }
        return c;
    }

    // The sites of every (species, gene) a run touches, computed once: the consensus (Consensus) over the congeners
    // that have the gene and whose copy pairs with the species' (CompareCopies): the gene's nearest congener by
    // alignment (congener_gaps.tsv) and the species' nearest congeners (species_neighbours.tsv, nearest first), up to
    // kCongeners. Thread-safe; a shared empty Sites for a species without one. `genomes` is the GenomeLoader (HasGene,
    // GetGeneOMP), a template so that this header needs none of it.
    class Cache {
    public:
        // gaps: the database's congener_gaps.tsv, or nullptr (the species' neighbours alone); neighbours may be empty (the
        // gene's nearest by the gaps alone).
        template<typename Genomes>
        std::shared_ptr<Sites const> Get(uint32_t taxid, uint32_t geneid, Genomes& genomes,
                                         species_neighbours::Table const& neighbours,
                                         congener_gaps::Table const* gaps = nullptr) {
            static std::shared_ptr<Sites const> const none = std::make_shared<Sites const>();
            bool const by_gaps = gaps != nullptr && !gaps->Empty();
            if (neighbours.Empty() && !by_gaps) return none;
            uint64_t const key = (static_cast<uint64_t>(taxid) << 32) | geneid;
            {
                std::lock_guard<std::mutex> lock(m_mutex);
                if (auto const it = m_sites.find(key); it != m_sites.end()) return it->second;
            }
            std::shared_ptr<Sites const> result = none;
            if (genomes.HasGene(taxid, geneid)) {
                auto const own = genomes.GetGeneOMP(taxid, geneid).Sequence();
                std::string_view const view = own.View();
                if (view.size() >= kKmer && view.size() <= kMaxLength) {
                    std::vector<uint32_t> candidates;
                    if (by_gaps) {
                        auto const gap = gaps->Find(taxid, geneid);
                        if (gap && gap->nearest != 0 && gap->nearest != taxid) candidates.push_back(gap->nearest);
                    }
                    for (auto const& n : neighbours.Of(taxid)) {
                        if (n.taxid != taxid && std::find(candidates.begin(), candidates.end(), n.taxid) == candidates.end()) {
                            candidates.push_back(n.taxid);
                        }
                    }
                    auto const own_kmers = detail::UniqueKmers(view);
                    std::vector<Comparison> comparisons;
                    for (uint32_t const congener : candidates) {
                        if (comparisons.size() >= kCongeners) break;
                        if (!genomes.HasGene(congener, geneid)) continue;
                        auto const other = genomes.GetGeneOMP(congener, geneid).Sequence();
                        auto comparison = CompareCopies(view, own_kmers, other.View());
                        if (!comparison) continue;
                        comparison->sites.congener = congener;
                        comparisons.push_back(std::move(*comparison));
                    }
                    if (!comparisons.empty()) result = std::make_shared<Sites const>(Consensus(view.size(), comparisons));
                }
            }
            std::lock_guard<std::mutex> lock(m_mutex);
            if (m_sites.size() >= kMaxCached) m_sites.clear();
            return m_sites.try_emplace(key, result).first->second;
        }

        size_t Size() const {
            std::lock_guard<std::mutex> lock(m_mutex);
            return m_sites.size();
        }

        void Clear() {
            std::lock_guard<std::mutex> lock(m_mutex);
            m_sites.clear();
        }

    private:
        mutable std::mutex m_mutex;
        std::unordered_map<uint64_t, std::shared_ptr<Sites const>> m_sites;
    };
}
