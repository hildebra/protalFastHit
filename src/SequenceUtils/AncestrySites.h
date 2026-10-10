// AncestrySites.h - where a species' copy of a marker gene differs from its congeners' copies, and which side a read
// takes there. The sites are the species' own derived states: the positions where the congeners compared agree on one
// base (kConsensus of those compared at the position) and the species carries another, and the indels the congeners
// share against the species' copy (its own gain or loss of a few bases). A strain of the species carries the species'
// base, and runs through the species' extra bases, there; a species the database lacks that branched off the lineage
// below a site carries the congeners' base, and deletes or inserts what they do. The reads of a false positive are
// mostly a novel congener's, at the same identity as a missed strain's, and this is the one place in the read where
// the two differ (docs/claude/2026-10-07-error-read-signatures). The index's unique k-mers see part of it
// (lu_per_kb); counting the sites with the reads' bases does not lose a site to every sequencing error around it.
//
// 0.7.9 compared the nearest congener alone. Half of the differences between two species are the congener's own
// derived states, where any third species agrees with the species, so a novel congener agreed at 0.68 instead of near
// 0 at r226 v17, and the residual errors were where that signal inverts (docs/claude/2026-10-08-r226-v17). The
// consensus over the congeners (one vote per species, kConsensus of those compared at the position, kMinCongeners of
// them at least) leaves the species' own derived sites; a large genus, whose nearest congener shares most of the
// species' history, gets as many as a small one; with fewer congeners the sites are the nearest's differences, as
// before. With kTolerantFrom congeners or more compared, one of them may carry a third base, neither the species' nor
// the others' (its own change at the position, which says nothing about the species' state): a world grown on known
// trees showed the nine-in-ten rule, which needs all of six congeners, losing a fifth of the species' derived sites
// to such a congener (docs/claude/2026-10-09-ancestry-true-positive-test). A congener that carries the species' base
// still blocks the site: it may share the state by descent, and then a species that branched off below the site
// carries it too.
//
// Two copies are compared along their shared 12-mers, without an alignment: the 12-mers unique to each copy pair them
// on diagonals (the position on one copy minus the position on the other), the pairs are chained in order along the
// copies, a change of diagonal between two chained pairs of at most kMaxIndel is an indel (located by extending the
// exact matches inwards from both pairs; the stretch left between them is not compared), and the bases are compared on
// the stretches between two chained pairs on one diagonal. A position no congener was compared at is no site. Per copy
// the sites are a sorted list of positions with the congeners' base and a short list of indels: a few hundred bytes.
// They are computed once per run for each (species, gene) a record touches (Cache: the gene store holds every
// reference's copy; the congeners are the gene's nearest by alignment in congener_gaps.tsv and the species' nearest in
// species_neighbours.tsv, up to kCongeners) and counted per record from its CIGAR (M: the species' base at every
// site; X: the read's base looked up; D and I: against the indel sites). The gaps features (CongenerGaps.h) take their
// nearest congener from the same table.
#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
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
                                                   // one base other than the species' (or one indel) for a site (three of
                                                   // three, ..., nine of ten)
    inline constexpr size_t kMinCongeners = 3;     // congeners compared at a position for the consensus to apply; with
                                                   // fewer the site is the nearest congener's difference, as in 0.7.9
    inline constexpr size_t kTolerantFrom = 6;     // congeners compared at a position from which one of them may carry
                                                   // a third base (or another indel) at a site: all but one of six
                                                   // carry the congeners' state, the one its own
    inline constexpr double kMinCompared = 0.5;    // of the species' copy compared, else the congener's copy is no use
    inline constexpr int64_t kMaxIndel = 60;       // bases a chained pair may shift the diagonal by: a longer one ends the chain
    inline constexpr size_t kMinIndelLength = 3;   // an indel site is at least this long: shorter ones are frameshifts in
                                                   // the references (assembly errors) or sequencing errors in the reads
    inline constexpr size_t kIndelTolerance = 4;   // bases a record's gap may lie from an indel site (the aligner slides gaps)
    inline constexpr size_t kIndelMargin = 5;      // a record counts an indel site only when aligned this far beyond it
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

    // An indel of the species' copy against the congeners': length > 0, the copy's bases [position, position + length)
    // that the congeners lack (a strain's read runs through them, a congener's read deletes them); length < 0, the
    // -length bases the congeners have before `position` that the copy lacks (a congener's read inserts them there).
    struct Indel {
        uint16_t position = 0;
        int16_t length = 0;
        bool operator==(Indel const&) const = default;
    };

    // The sites of one copy: against one congener's copy (CompareCopies: where the congener differs, with its base, and
    // its indels), or the consensus over its congeners (Consensus: where kConsensus of those compared carry one other
    // base, with that base, and the indels they share).
    struct Sites {
        std::vector<uint16_t> positions;  // on the species' copy, 0-based, ascending
        std::vector<uint8_t> bases;       // the congeners' base at each (Code)
        std::vector<Indel> indels;        // by position
        uint32_t congener = 0;            // the nearest congener compared (the highest identity); 0: none
        uint16_t congeners = 0;           // congeners compared
        size_t compared = 0;              // positions compared with any congener
        float identity = 0;               // of the species' copy and the nearest congener's over the bases compared

        bool Empty() const { return positions.empty() && indels.empty(); }

        // The base sites in [begin, end) of the copy.
        size_t Count(size_t begin, size_t end) const {
            auto const lo = std::lower_bound(positions.begin(), positions.end(), begin);
            auto const hi = std::lower_bound(positions.begin(), positions.end(), end);
            return static_cast<size_t>(hi - lo);
        }

        // The index of the base site at position p, or SIZE_MAX.
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

        // A paired 12-mer: its position on the species' copy and its diagonal (that position minus the one on the other).
        struct Anchor {
            uint32_t position;
            int64_t diagonal;
        };

        // Chains the anchors (sorted by position) from the main diagonal outwards: along it, and across a change of
        // diagonal of at most kMaxIndel that keeps the order on both copies and is followed by another anchor on the new
        // diagonal (a single stray pair is skipped). Returns the chained anchors in order of position.
        inline std::vector<Anchor> Chain(std::vector<Anchor> const& anchors, int64_t main) {
            std::vector<Anchor> chain;
            if (anchors.empty()) return chain;
            size_t first = anchors.size();
            for (size_t i = 0; i < anchors.size(); i++) {
                if (anchors[i].diagonal == main) {
                    first = i;
                    break;
                }
            }
            if (first == anchors.size()) return chain;
            // Whether anchor i may follow `last` on the chain (forward), or precede it (backward): the next anchor on
            // the copy and on the other copy, by at most kMaxIndel of shift, with a second anchor on a new diagonal.
            auto accepts = [&anchors](Anchor const& last, size_t i, bool forward) {
                auto const& a = anchors[i];
                int64_t const shift = a.diagonal - last.diagonal;
                if (shift < -kMaxIndel || shift > kMaxIndel) return false;
                int64_t const own = static_cast<int64_t>(a.position) - static_cast<int64_t>(last.position);
                int64_t const other = own - shift;  // the step on the other copy
                if (forward ? (own <= 0 || other <= 0) : (own >= 0 || other >= 0)) return false;
                if (shift == 0) return true;
                size_t const next = forward ? i + 1 : i - 1;
                return (forward ? i + 1 < anchors.size() : i > 0) && anchors[next].diagonal == a.diagonal;
            };
            std::vector<Anchor> before;  // the chain to the left of `first`, found backwards
            Anchor last = anchors[first];
            for (size_t i = first; i-- > 0;) {
                if (!accepts(last, i, false)) continue;
                before.push_back(anchors[i]);
                last = anchors[i];
            }
            chain.assign(before.rbegin(), before.rend());
            chain.push_back(anchors[first]);
            last = anchors[first];
            for (size_t i = first + 1; i < anchors.size(); i++) {
                if (!accepts(last, i, true)) continue;
                chain.push_back(anchors[i]);
                last = anchors[i];
            }
            return chain;
        }
    }

    // The comparison of a copy with one congener's: where the congener differs (sites, with its indels), and which
    // positions of the copy were compared at all (covered: one char per position, 1 when compared).
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
        // The shared unique k-mers and their diagonals.
        std::vector<detail::Anchor> anchors;
        for (size_t i = 0, j = 0; i < a.size() && j < b.size();) {
            if (a[i].first < b[j].first) i++;
            else if (b[j].first < a[i].first) j++;
            else {
                anchors.push_back({ a[i].second, static_cast<int64_t>(a[i].second) - static_cast<int64_t>(b[j].second) });
                i++;
                j++;
            }
        }
        if (anchors.size() < 2) return std::nullopt;
        // The main diagonal: the one most anchors lie on.
        std::vector<int64_t> diagonals;
        diagonals.reserve(anchors.size());
        for (auto const& x : anchors) diagonals.push_back(x.diagonal);
        std::sort(diagonals.begin(), diagonals.end());
        int64_t main = 0;
        size_t best = 0;
        for (size_t i = 0; i < diagonals.size();) {
            size_t j = i;
            while (j < diagonals.size() && diagonals[j] == diagonals[i]) j++;
            if (j - i > best) {
                best = j - i;
                main = diagonals[i];
            }
            i = j;
        }
        std::sort(anchors.begin(), anchors.end(), [](detail::Anchor const& x, detail::Anchor const& y) {
            return x.position != y.position ? x.position < y.position : x.diagonal < y.diagonal;
        });
        auto const chain = detail::Chain(anchors, main);
        if (chain.size() < 2) return std::nullopt;
        Comparison c;
        c.covered.assign(own.size(), 0);
        // The base at own[i] on the other copy along diagonal d, or 4 beyond it or for another letter.
        auto other_at = [&own, &other](size_t i, int64_t d) -> uint8_t {
            int64_t const j = static_cast<int64_t>(i) - d;
            if (j < 0 || j >= static_cast<int64_t>(other.size())) return 4;
            return Code(other[static_cast<size_t>(j)]);
        };
        // The stretches between consecutive chained anchors: compared base by base on one diagonal; across a change of
        // diagonal the exact matches are extended inwards from both anchors, and what is left between them holds the indel.
        auto cover = [&](size_t from, size_t to, int64_t d) {  // [from, to) compared on diagonal d
            for (size_t i = from; i < to && i < own.size(); i++) {
                uint8_t const x = Code(own[i]), y = other_at(i, d);
                if (x > 3 || y > 3) continue;
                c.covered[i] = 1;
                if (x != y) {
                    c.sites.positions.push_back(static_cast<uint16_t>(i));
                    c.sites.bases.push_back(y);
                }
            }
        };
        for (size_t k = 0; k + 1 < chain.size(); k++) {
            auto const& left = chain[k];
            auto const& right = chain[k + 1];
            if (left.diagonal == right.diagonal) {
                cover(left.position, std::min<size_t>(own.size(), right.position + kKmer), left.diagonal);
                continue;
            }
            cover(left.position, left.position + kKmer, left.diagonal);
            cover(right.position, std::min<size_t>(own.size(), right.position + kKmer), right.diagonal);
            size_t i = left.position + kKmer;
            while (i < right.position && Code(own[i]) < 4 && Code(own[i]) == other_at(i, left.diagonal)) c.covered[i++] = 1;
            size_t j = right.position;
            while (j > i && Code(own[j - 1]) < 4 && Code(own[j - 1]) == other_at(j - 1, right.diagonal)) c.covered[--j] = 1;
            int64_t const shift = right.diagonal - left.diagonal;  // > 0: own has `shift` bases more here, < 0 fewer
            if (static_cast<size_t>(std::llabs(shift)) >= kMinIndelLength) {
                c.sites.indels.push_back({ static_cast<uint16_t>(std::min<size_t>(i, kMaxLength)), static_cast<int16_t>(shift) });
            }
        }
        size_t const n = static_cast<size_t>(std::count(c.covered.begin(), c.covered.end(), 1));
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

    // The sites of `own` against one congener's copy `other`: where other differs, with its base, and its indels
    // (nullopt as CompareCopies).
    inline std::optional<Sites> Compare(std::string_view own, std::string_view other) {
        auto c = CompareCopies(own, other);
        if (!c) return std::nullopt;
        return std::move(c->sites);
    }

    // Whether the congeners compared at a position agree on a state other than the species': `best` of `compared` carry
    // one such state, `other` another one (neither the species' nor the best). kConsensus of them (three of three, nine
    // of ten), or, from kTolerantFrom compared, all but one, when that one carries a state of its own: a congener with
    // the species' state is never tolerated (above).
    inline bool Agree(size_t best, size_t other, size_t compared) {
        if (best == 0) return false;
        if (static_cast<double>(best) >= kConsensus * static_cast<double>(compared)) return true;
        return compared >= kTolerantFrom && best + 1 == compared && best + other == compared;
    }

    // The consensus sites of a copy of `length` over its comparisons with its congeners (one vote each): at the positions
    // compared with kMinCongeners congeners or more, those where the congeners Agree on one base other than the
    // species' (kConsensus of them; from kTolerantFrom all but one that carries a third base); at the positions compared
    // with fewer (a small genus, or a stretch past an indel in the others' copies), those where the nearest congener
    // (the highest identity) differs, with its base, as 0.7.9 had everywhere. The indels alike: a congener compared at
    // the base before or after an indel votes for it (the same position and length) or against it; one that carries
    // another indel at the same position is the tolerated dissenter. congener and identity are the nearest congener's,
    // compared counts the positions any congener was compared at, congeners how many were. Empty without a comparison.
    // `outgroup`, if given (the family's consensus base per position of the copy, column_weights::Columns::consensus,
    // 4 for none), polarises the positions compared with fewer than kMinCongeners congeners: there the nearest
    // congener's difference is a site only where the congener's base is the family's (the species' base is derived),
    // not where the congener's own base is the new one, which the 0.7.9 fallback could not tell (half of a congener
    // pair's differences are the congener's own derived states, and a novel species agrees with the species at them).
    inline Sites Consensus(size_t length, std::vector<Comparison> const& comparisons, std::vector<uint8_t> const* outgroup = nullptr) {
        Sites out;
        if (comparisons.empty() || length > kMaxLength) return out;
        auto const nearest = std::max_element(comparisons.begin(), comparisons.end(), [](Comparison const& a, Comparison const& b) {
            return a.sites.identity < b.sites.identity;
        });
        std::vector<uint8_t> compared(length, 0);
        std::vector<std::array<uint8_t, 4>> votes(length, std::array<uint8_t, 4>{ 0, 0, 0, 0 });
        std::vector<int8_t> fallback(length, -1);  // the nearest congener's base where it differs
        std::vector<std::pair<Indel, uint16_t>> indel_votes;
        for (auto const& c : comparisons) {
            for (size_t i = 0; i < length && i < c.covered.size(); i++) compared[i] += c.covered[i] != 0;
            for (size_t s = 0; s < c.sites.positions.size(); s++) {
                if (c.sites.positions[s] >= length) continue;
                votes[c.sites.positions[s]][c.sites.bases[s] & 3]++;
                if (&c == &*nearest) fallback[c.sites.positions[s]] = static_cast<int8_t>(c.sites.bases[s] & 3);
            }
            for (auto const& indel : c.sites.indels) {
                auto it = std::find_if(indel_votes.begin(), indel_votes.end(), [&indel](auto const& v) { return v.first == indel; });
                if (it == indel_votes.end()) indel_votes.emplace_back(indel, 1);
                else it->second++;
            }
        }
        for (size_t i = 0; i < length; i++) {
            if (compared[i] == 0) continue;
            out.compared++;
            if (compared[i] < kMinCongeners) {
                bool polarised = outgroup != nullptr && i < outgroup->size() && (*outgroup)[i] < 4;
                if (fallback[i] >= 0 && (!polarised || (*outgroup)[i] == static_cast<uint8_t>(fallback[i]))) {
                    out.positions.push_back(static_cast<uint16_t>(i));
                    out.bases.push_back(static_cast<uint8_t>(fallback[i]));
                }
                continue;
            }
            uint8_t best = 0;
            for (uint8_t b = 1; b < 4; b++) {
                if (votes[i][b] > votes[i][best]) best = b;
            }
            size_t const other = static_cast<size_t>(votes[i][0]) + votes[i][1] + votes[i][2] + votes[i][3] - votes[i][best];
            if (Agree(votes[i][best], other, compared[i])) {
                out.positions.push_back(static_cast<uint16_t>(i));
                out.bases.push_back(best);
            }
        }
        // The congeners compared at an indel: those compared at the base before it or at the first base after it.
        auto compared_at = [&](Indel const& indel) {
            size_t const end = static_cast<size_t>(indel.position) + (indel.length > 0 ? static_cast<size_t>(indel.length) : 0);
            size_t n = 0;
            for (auto const& c : comparisons) {
                bool const before = indel.position > 0 && indel.position - 1 < c.covered.size() && c.covered[indel.position - 1];
                bool const after = end < c.covered.size() && c.covered[end];
                n += before || after;
            }
            return n;
        };
        for (auto const& [indel, n] : indel_votes) {
            size_t const at = compared_at(indel);
            bool site = false;
            if (at < kMinCongeners) {
                site = std::find(nearest->sites.indels.begin(), nearest->sites.indels.end(), indel) != nearest->sites.indels.end();
            } else {
                size_t other = 0;  // votes for another indel at the same position (a congener has one indel there at most)
                for (auto const& [o, m] : indel_votes) {
                    if (o.position == indel.position && o.length != indel.length) other += m;
                }
                site = Agree(n, std::min(other, at > n ? at - n : 0), at);
            }
            if (site) out.indels.push_back(indel);
        }
        std::sort(out.indels.begin(), out.indels.end(), [](Indel const& a, Indel const& b) {
            return a.position != b.position ? a.position < b.position : a.length < b.length;
        });
        out.congener = nearest->sites.congener;
        out.identity = nearest->sites.identity;
        out.congeners = static_cast<uint16_t>(std::min<size_t>(comparisons.size(), UINT16_MAX));
        return out;
    }

    // What a record covers of the sites: the base sites it covers (M, = and X ops; a deletion skips its sites), those
    // where the read has the species' base (M and =, and X where the read's base is neither), and those where it has
    // the congeners'; and the indel sites it is aligned kIndelMargin beyond on both sides: those where it has the
    // congeners' indel (a deletion of the species' extra bases, or an insertion of the congeners', of the site's
    // length within kIndelTolerance of it) and those where it runs through with no gap near (the species' state); a
    // record with another gap near an indel site counts it neither way. A record without a sequence ("*") counts no
    // base site.
    struct Counts {
        uint64_t sites = 0, agree = 0, congener = 0;
        uint64_t indel_sites = 0, indel_agree = 0, indel_congener = 0;
    };

    inline Counts Count(Sites const& sites, std::string const& cigar, size_t pos, std::string const& seq) {
        Counts c;
        if (sites.Empty()) return c;
        bool const bases = !sites.positions.empty() && !seq.empty() && seq != "*";
        struct Gap {
            size_t ref, length;
            bool deletion;
        };
        std::vector<Gap> gaps;
        size_t const start = pos == 0 ? 0 : pos - 1;
        size_t ref = start, query = 0, run = 0;
        for (char const op : cigar) {
            if (op >= '0' && op <= '9') {
                run = run * 10 + static_cast<size_t>(op - '0');
                continue;
            }
            if (op == 'M' || op == '=') {
                if (bases) {
                    size_t const n = sites.Count(ref, ref + run);
                    c.sites += n;
                    c.agree += n;
                }
                ref += run;
                query += run;
            } else if (op == 'X') {
                if (bases) {
                    auto it = std::lower_bound(sites.positions.begin(), sites.positions.end(), ref);
                    for (; it != sites.positions.end() && *it < ref + run; ++it) {
                        size_t const q = query + (*it - ref);
                        if (q >= seq.size()) break;
                        c.sites++;
                        if (Code(seq[q]) == sites.bases[static_cast<size_t>(it - sites.positions.begin())]) c.congener++;
                    }
                }
                ref += run;
                query += run;
            } else if (op == 'D' || op == 'N') {
                if (!sites.indels.empty()) gaps.push_back({ ref, run, true });
                ref += run;
            } else if (op == 'I') {
                if (!sites.indels.empty()) gaps.push_back({ ref, run, false });
                query += run;
            } else if (op == 'S') {
                query += run;
            }
            run = 0;
        }
        for (auto const& site : sites.indels) {
            size_t const n = static_cast<size_t>(std::abs(static_cast<int>(site.length)));
            if (n < kMinIndelLength) continue;
            size_t const lo = site.position, hi = site.position + (site.length > 0 ? n : 0);
            if (lo < start + kIndelMargin || hi + kIndelMargin > ref) continue;  // not aligned through the site
            bool matched = false, other = false;
            for (auto const& g : gaps) {
                size_t const g_hi = g.ref + (g.deletion ? g.length : 0);
                if (g_hi + kIndelTolerance < lo || g.ref > hi + kIndelTolerance) continue;
                if (g.deletion == (site.length > 0) && g.length == n) matched = true;
                else other = true;
            }
            if (matched) {
                c.indel_sites++;
                c.indel_congener++;
            } else if (!other) {
                c.indel_sites++;
                c.indel_agree++;
            }
        }
        return c;
    }

    // Each base site a record covers, as Count walks it: fn(position, outcome), the outcome kSpeciesBase where the read
    // has the species' base, kCongenerBase the congeners', kOtherBase another. Nothing without a sequence.
    inline constexpr int kSpeciesBase = 0, kCongenerBase = 1, kOtherBase = 2;

    template<typename Fn>
    void ForEachSite(Sites const& sites, std::string const& cigar, size_t pos, std::string const& seq, Fn&& fn) {
        if (sites.positions.empty() || seq.empty() || seq == "*") return;
        size_t ref = pos == 0 ? 0 : pos - 1, query = 0, run = 0;
        for (char const op : cigar) {
            if (op >= '0' && op <= '9') {
                run = run * 10 + static_cast<size_t>(op - '0');
                continue;
            }
            if (op == 'M' || op == '=' || op == 'X') {
                auto it = std::lower_bound(sites.positions.begin(), sites.positions.end(), ref);
                for (; it != sites.positions.end() && *it < ref + run; ++it) {
                    if (op != 'X') {
                        fn(static_cast<size_t>(*it), kSpeciesBase);
                        continue;
                    }
                    size_t const q = query + (*it - ref);
                    if (q >= seq.size()) break;
                    bool const congener = Code(seq[q]) == sites.bases[static_cast<size_t>(it - sites.positions.begin())];
                    fn(static_cast<size_t>(*it), congener ? kCongenerBase : kOtherBase);
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
        // `outgroup`: the family's consensus base per position of the copy (ColumnWeights.h), for Consensus; may be null.
        template<typename Genomes>
        std::shared_ptr<Sites const> Get(uint32_t taxid, uint32_t geneid, Genomes& genomes,
                                         species_neighbours::Table const& neighbours,
                                         congener_gaps::Table const* gaps = nullptr,
                                         std::vector<uint8_t> const* outgroup = nullptr) {
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
                    if (!comparisons.empty()) result = std::make_shared<Sites const>(Consensus(view.size(), comparisons, outgroup));
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
