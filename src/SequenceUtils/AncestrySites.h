// AncestrySites.h - where a species' copy of a marker gene differs from its nearest congener's copy, and which side a
// read takes there. At those sites a strain of the species carries the species' base (its derived states, and the
// congener's); a species that branched off the lineage below some of them carries the congener's base at those. The
// reads of a false positive are mostly a novel congener's, at the same identity as a missed strain's, and this is
// the one place in the read where the two differ (docs/claude/2026-10-07-error-read-signatures). The index's unique
// k-mers see part of it (lu_per_kb); counting the sites with the reads' bases does not lose a site to every
// sequencing error around it.
//
// The two copies are compared along their shared 12-mers, without an alignment: the 12-mers unique to each copy pair
// them on their main diagonal, and the bases are compared on the stretches between two paired 12-mers that lie on it
// (an indel moves the diagonal and ends the stretch, so a stretch past it is not compared). Per copy the sites are a
// sorted list of positions with the congener's base: a few hundred bytes. They are computed once per run for each
// (species, gene) a record touches (Cache: the gene store holds every reference's copy, species_neighbours.tsv names
// the congeners) and counted per record from its CIGAR (M: the species' base at every site; X: the read's base
// looked up).
#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "SpeciesNeighbours.h"

namespace protal::ancestry {
    inline constexpr size_t kKmer = 12;            // the copies are paired by their shared unique 12-mers
    inline constexpr size_t kCongenersTried = 3;   // nearest congeners tried for a gene: the first whose copy pairs
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

    // The sites of one copy against its nearest congener's.
    struct Sites {
        std::vector<uint16_t> positions;  // on the species' copy, 0-based, ascending
        std::vector<uint8_t> bases;       // the congener's base at each (Code)
        uint32_t congener = 0;            // the congener compared; 0: none
        size_t compared = 0;              // bases compared
        float identity = 0;               // of the two copies over the bases compared

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
        // The 12-mers of a sequence as (value, position), every k-mer of A/C/G/T only.
        inline std::vector<std::pair<uint32_t, uint32_t>> Kmers(std::string_view s) {
            std::vector<std::pair<uint32_t, uint32_t>> out;
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
        inline void KeepUnique(std::vector<std::pair<uint32_t, uint32_t>>& kmers) {
            std::sort(kmers.begin(), kmers.end());
            std::vector<std::pair<uint32_t, uint32_t>> unique;
            unique.reserve(kmers.size());
            for (size_t i = 0; i < kmers.size();) {
                size_t j = i;
                while (j < kmers.size() && kmers[j].first == kmers[i].first) j++;
                if (j == i + 1) unique.push_back(kmers[i]);
                i = j;
            }
            kmers.swap(unique);
        }
    }

    // The sites of `own` against `other`: nullopt if fewer than kMinCompared of own's bases could be compared (the
    // copies do not pair: another gene, or too different), or own is longer than kMaxLength.
    inline std::optional<Sites> Compare(std::string_view own, std::string_view other) {
        if (own.size() > kMaxLength || own.size() < kKmer || other.size() < kKmer) return std::nullopt;
        auto a = detail::Kmers(own);
        auto b = detail::Kmers(other);
        detail::KeepUnique(a);
        detail::KeepUnique(b);
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
        std::vector<char> compared(own.size(), 0);
        uint32_t previous = 0;
        bool first = true;
        for (auto const& [d, p] : pairs) {
            if (d != diagonal) continue;
            if (!first) std::fill(compared.begin() + previous, compared.begin() + std::min<size_t>(own.size(), p + kKmer), 1);
            previous = p;
            first = false;
        }
        Sites sites;
        size_t n = 0;
        for (size_t i = 0; i < own.size(); i++) {
            if (!compared[i]) continue;
            int64_t const j = static_cast<int64_t>(i) - diagonal;
            if (j < 0 || j >= static_cast<int64_t>(other.size())) continue;
            uint8_t const x = Code(own[i]), y = Code(other[static_cast<size_t>(j)]);
            if (x > 3 || y > 3) continue;
            n++;
            if (x != y) {
                sites.positions.push_back(static_cast<uint16_t>(i));
                sites.bases.push_back(y);
            }
        }
        if (static_cast<double>(n) < kMinCompared * static_cast<double>(own.size())) return std::nullopt;
        sites.compared = n;
        sites.identity = static_cast<float>(1.0 - static_cast<double>(sites.positions.size()) / static_cast<double>(n));
        return sites;
    }

    // What a record covers of the sites: those it covers (M, = and X ops; a deletion skips its sites), those where
    // the read has the species' base (M and =, and X where the read's base is neither), and those where it has the
    // congener's. A record without a sequence ("*") counts nothing.
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

    // The sites of every (species, gene) a run touches, computed once: against the nearest congener
    // (species_neighbours.tsv, nearest first) that has the gene and whose copy pairs with the species' (Compare), of
    // the kCongenersTried nearest. Thread-safe; a shared empty Sites for a species without one. `genomes` is the
    // GenomeLoader (HasGene, GetGeneOMP), a template so that this header needs none of it.
    class Cache {
    public:
        template<typename Genomes>
        std::shared_ptr<Sites const> Get(uint32_t taxid, uint32_t geneid, Genomes& genomes,
                                         species_neighbours::Table const& neighbours) {
            static std::shared_ptr<Sites const> const none = std::make_shared<Sites const>();
            if (neighbours.Empty()) return none;
            uint64_t const key = (static_cast<uint64_t>(taxid) << 32) | geneid;
            {
                std::lock_guard<std::mutex> lock(m_mutex);
                if (auto const it = m_sites.find(key); it != m_sites.end()) return it->second;
            }
            std::shared_ptr<Sites const> result = none;
            if (genomes.HasGene(taxid, geneid)) {
                auto const own = genomes.GetGeneOMP(taxid, geneid).Sequence();
                size_t tried = 0;
                for (auto const& n : neighbours.Of(taxid)) {
                    if (tried >= kCongenersTried) break;
                    if (!genomes.HasGene(n.taxid, geneid)) continue;
                    tried++;
                    auto const other = genomes.GetGeneOMP(n.taxid, geneid).Sequence();
                    if (auto sites = Compare(own.View(), other.View())) {
                        sites->congener = n.taxid;
                        result = std::make_shared<Sites const>(std::move(*sites));
                        break;
                    }
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
