#pragma once

// The species a read is from, across its parts: the two mates of a pair whose mates aligned to different genes, or
// the genes of a long read. Each part has candidate alignments to several taxa (a gene and its homologs in the
// relatives). A read is one organism, so its parts are given one taxon: the read's consensus.

#include <algorithm>
#include <climits>
#include <cstdint>
#include <optional>
#include <vector>
#include "AlignmentUtils.h"

namespace protal {
    // A part's candidate alignment: its taxon and score (bitscore, Score(2, 3, 1, 2)).
    struct PartCandidate {
        uint32_t taxid = 0;
        int score = 0;
    };

    struct ReadConsensus {
        uint32_t taxid = 0;
        size_t parts = 0;  // the read's parts with a candidate of the taxon
        // MAPQ of the read's assignment: on the taxon's parts, its scores against the closest other taxon's (0 where
        // that taxon has no candidate); -1 if no other taxon has a candidate on them.
        int mapq = -1;
    };

    // The consensus taxon of a read whose parts have these candidates. Two taxa are compared on the parts where both
    // have a candidate, by the sum of their best scores there; a taxon that scores less than another on their common
    // parts loses to it. A part without a candidate of a taxon is no evidence against it: the database may lack that
    // gene of the read's species, whose reads then align to a relative only. The consensus is a taxon that loses to
    // none (if every taxon loses to one, any), with candidates on the most parts, then the highest total score.
    // nullopt for a read without candidates, or if the consensus has candidates on fewer than `min_parts` parts.
    inline std::optional<ReadConsensus> ConsensusOfRead(std::vector<std::vector<PartCandidate>> const& parts, size_t min_parts = 1) {
        constexpr long kAbsent = LONG_MIN;
        std::vector<uint32_t> taxa;
        for (auto const& part : parts) {
            for (auto const& c : part) {
                if (std::find(taxa.begin(), taxa.end(), c.taxid) == taxa.end()) taxa.push_back(c.taxid);
            }
        }
        if (taxa.empty()) return std::nullopt;
        size_t const n = taxa.size();
        std::vector<std::vector<long>> best(n, std::vector<long>(parts.size(), kAbsent));
        for (size_t p = 0; p < parts.size(); p++) {
            for (auto const& c : parts[p]) {
                auto const t = static_cast<size_t>(std::find(taxa.begin(), taxa.end(), c.taxid) - taxa.begin());
                best[t][p] = std::max(best[t][p], static_cast<long>(c.score));
            }
        }
        // The sums of a and b over the parts both have; false if they have none in common.
        auto common = [&](size_t a, size_t b, long& sum_a, long& sum_b) {
            bool any = false;
            sum_a = sum_b = 0;
            for (size_t p = 0; p < parts.size(); p++) {
                if (best[a][p] == kAbsent || best[b][p] == kAbsent) continue;
                sum_a += best[a][p];
                sum_b += best[b][p];
                any = true;
            }
            return any;
        };
        std::vector<bool> undefeated(n, true);
        for (size_t a = 0; a < n; a++) {
            for (size_t b = 0; b < n && undefeated[a]; b++) {
                long sa, sb;
                if (a != b && common(a, b, sa, sb) && sa < sb) undefeated[a] = false;
            }
        }
        bool const any_undefeated = std::find(undefeated.begin(), undefeated.end(), true) != undefeated.end();
        auto parts_of = [&](size_t t) { return static_cast<size_t>(std::count_if(best[t].begin(), best[t].end(), [&](long s) { return s != kAbsent; })); };
        auto total_of = [&](size_t t) {
            long total = 0;
            for (long s : best[t]) total += s == kAbsent ? 0 : s;
            return total;
        };
        size_t pick = n;
        for (size_t t = 0; t < n; t++) {
            if (any_undefeated && !undefeated[t]) continue;
            if (pick == n || parts_of(t) > parts_of(pick) || (parts_of(t) == parts_of(pick) && total_of(t) > total_of(pick))) pick = t;
        }
        ReadConsensus consensus{ taxa[pick], parts_of(pick), -1 };
        if (consensus.parts < min_parts) return std::nullopt;
        // Confidence: on the consensus taxon's parts, against each other taxon with a candidate there.
        for (size_t u = 0; u < n; u++) {
            if (u == pick) continue;
            long own = 0, other = 0;
            bool shares = false;
            for (size_t p = 0; p < parts.size(); p++) {
                if (best[pick][p] == kAbsent) continue;
                own += best[pick][p];
                if (best[u][p] != kAbsent) {
                    other += best[u][p];
                    shares = true;
                }
            }
            if (!shares) continue;
            int const mapq = own > 0 && own > other ? std::max(0, MAPQv2(static_cast<int>(own), static_cast<int>(std::max(0L, other)))) : 0;
            consensus.mapq = consensus.mapq < 0 ? mapq : std::min(consensus.mapq, mapq);
        }
        return consensus;
    }
}
