#pragma once

// SCRATCH VARIANT (not for commit): ConsensusOfRead with the rule chosen by PROTAL_CONSENSUS:
//   old      taxa compared by the sum of their scores on common parts; all beaten: the most parts (before 2026-10-13)
//   minimax  the same comparison; all beaten: the smallest largest defeat among those on min_parts (the patch)
//   vote     taxa compared by the common parts each wins (the sum breaks a tie); all beaten: the most parts
//   vote_minimax  the vote comparison; all beaten: the smallest largest defeat (in parts won, then score)

#include <algorithm>
#include <climits>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <optional>
#include <vector>
#include "AlignmentUtils.h"

namespace protal {
    struct PartCandidate {
        uint32_t taxid = 0;
        int score = 0;
    };

    struct ReadConsensus {
        uint32_t taxid = 0;
        size_t parts = 0;
        int mapq = -1;
    };

    inline int ConsensusVariant() {
        static int const variant = [] {
            char const* v = std::getenv("PROTAL_CONSENSUS");
            if (v == nullptr || std::strcmp(v, "minimax") == 0) return 1;
            if (std::strcmp(v, "old") == 0) return 0;
            if (std::strcmp(v, "vote") == 0) return 2;
            if (std::strcmp(v, "vote_minimax") == 0) return 3;
            return 1;
        }();
        return variant;
    }

    inline std::optional<ReadConsensus> ConsensusOfRead(std::vector<std::vector<PartCandidate>> const& parts, size_t min_parts = 1) {
        int const variant = ConsensusVariant();
        bool const vote = variant >= 2;
        bool const minimax = variant == 1 || variant == 3;
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
        // How far a falls behind b on their common parts (> 0: a loses), false without common parts. Sum: by the score;
        // vote: by the parts won (scaled), the score's sum breaking a tie.
        auto behind = [&](size_t a, size_t b, long& margin) {
            bool any = false;
            long sum_a = 0, sum_b = 0, won_a = 0, won_b = 0;
            for (size_t p = 0; p < parts.size(); p++) {
                if (best[a][p] == kAbsent || best[b][p] == kAbsent) continue;
                sum_a += best[a][p];
                sum_b += best[b][p];
                won_a += best[a][p] > best[b][p];
                won_b += best[b][p] > best[a][p];
                any = true;
            }
            if (!vote) margin = sum_b - sum_a;
            else margin = won_b != won_a ? (won_b - won_a) * 1000000L : (sum_b - sum_a);
            return any;
        };
        std::vector<long> worst(n, 0);
        for (size_t a = 0; a < n; a++) {
            for (size_t b = 0; b < n; b++) {
                long m;
                if (a != b && behind(a, b, m) && m > 0) worst[a] = std::max(worst[a], m);
            }
        }
        auto parts_of = [&](size_t t) { return static_cast<size_t>(std::count_if(best[t].begin(), best[t].end(), [&](long s) { return s != kAbsent; })); };
        auto total_of = [&](size_t t) {
            long total = 0;
            for (long s : best[t]) total += s == kAbsent ? 0 : s;
            return total;
        };
        bool const any_undefeated = std::find(worst.begin(), worst.end(), 0L) != worst.end();
        std::vector<bool> considered(n);
        for (size_t t = 0; t < n; t++) considered[t] = worst[t] == 0;
        if (!any_undefeated) {
            if (minimax) {
                long least = LONG_MAX;
                for (size_t t = 0; t < n; t++) if (parts_of(t) >= min_parts) least = std::min(least, worst[t]);
                if (least == LONG_MAX) return std::nullopt;
                for (size_t t = 0; t < n; t++) considered[t] = parts_of(t) >= min_parts && worst[t] == least;
            } else {
                std::fill(considered.begin(), considered.end(), true);
            }
        }
        size_t pick = n;
        for (size_t t = 0; t < n; t++) {
            if (!considered[t]) continue;
            if (pick == n || parts_of(t) > parts_of(pick) || (parts_of(t) == parts_of(pick) && total_of(t) > total_of(pick))) pick = t;
        }
        ReadConsensus consensus{ taxa[pick], parts_of(pick), -1 };
        if (consensus.parts < min_parts) return std::nullopt;
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
