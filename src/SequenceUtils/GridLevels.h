// GridLevels.h - the gene grid of a sparse index: which of a gene's k-mer cores the index keeps for a run gated at C
// (--sparse_first_pick grid, docs/claude/2026-10-13-gene-grid-index.md).
//
// The stride pick spaces a read's sampled cores evenly by choosing them on the read, so any key may be one it picks and
// the whole index must be loaded. The grid chooses them on the reference instead: of each gene's cores in gene order (its
// syncmers), one per window of C consecutive cores is kept at C. A read cannot tell which of its cores those are (its
// offset in the gene is unknown), so it looks up every core and only the grid's hit: as many as the stride would pick,
// about every C-th core of the stretch it covers, and a run gated at C needs only the grid's entries in memory.
//
// The windows nest: at C = 2^l a gene's cores fall into windows of 2^l consecutive ordinals, each window's two halves
// are the windows of C/2, and the core a window keeps is the better of its halves' (a tournament). A core's level is the
// last round it won, 0 to 5 (GRID_LEVELS - 1), so the cores kept at C are those of level log2(C) or more, a subset of
// those kept at C/2, and one level per entry serves every gate.
//
// Which core a window keeps: the one a read of the species most likely carries intact and that hits the fewest other
// genes, near the window's middle:
//  - conservation within the species: of the strain alleles (StrainAlleles.h) that cover the core's k-mer, the share
//    that keep it a hit: no edit in the 15-base core, no insertion or deletion in the 31-base k-mer, and at most
//    kMaxFlankSubstitutions in its 16 flank bases (the gate's flank score, 14 of 16). A species without alleles there
//    (one genome, or alleles that do not reach) takes the prior: how conserved the core's columns are within the genera
//    of its family (ColumnWeights.h `within`, the rate a strain's own mutations hit a column); kUnknownPrior without
//    weights. The share and the prior are pooled, the prior as kPriorWeight alleles: (kept + w prior) / (covering + w);
//  - ubiquity: the key's entries in the whole index (every gene and taxon sharing the core), log2 of them over log2 of
//    the most a key keeps (2048), weighted kUbiquityWeight: a core shared by many genes hits them all and costs a long
//    block to scan;
//  - the window's middle: kCentreWeight per window length of distance, so that the kept cores stay spread out where
//    the scores tie;
//  - a core whose key the index leaves out (too many values, Seedmap::BuildValuePointers) is never kept.
#pragma once

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <span>
#include <vector>

#include "Constants.h"
#include "StrainAlleles.h"

namespace protal::grid_levels {
    inline constexpr unsigned kLevels = GRID_LEVELS;
    inline constexpr unsigned kCore = 15;                 // a k-mer's core
    inline constexpr unsigned kFlank = 8;                 // its flank bases on either side
    inline constexpr unsigned kMaxFlankSubstitutions = 2;  // the gate's flank score: 14 of 16
    inline constexpr double kPriorWeight = 2.0;           // the prior counts as this many alleles
    inline constexpr double kUnknownPrior = 0.5;          // a core without alleles or column weights
    inline constexpr uint8_t kPriorCode = 12;             // a within code at or above it counts as fully conserved
    inline constexpr double kUbiquityWeight = 0.5;
    inline constexpr double kMaxEntriesLog2 = 11.0;       // log2(2048): Seedmap::max_key_multiplicity
    inline constexpr double kCentreWeight = 0.05;

    // A level's gate: C = 2^level.
    inline constexpr size_t GateOf(unsigned level) { return size_t{1} << level; }
    // A gate's level: log2(C), C a power of two from 1 to 32.
    inline constexpr unsigned LevelOf(size_t gate) { return gate <= 1 ? 0u : static_cast<unsigned>(std::bit_width(gate) - 1); }

    // A gene's core as the grid scores it.
    struct Core {
        uint32_t position = 0;  // the core's first base in the gene (an entry's position)
        uint32_t entries = 0;   // its key's entries in the index (0: the index leaves the key out)
    };

    // Whether an allele (strain_alleles::Table::View: the representative's bases [begin, end) and its edits there,
    // sorted) keeps a core at `position` a hit: -1 if it does not cover the core's k-mer, else 1 or 0. `max_deletion`: its
    // longest deletion, so that a deletion starting before the k-mer and reaching into it is seen.
    inline int AlleleKeeps(strain_alleles::Table::View const& allele, uint32_t position, uint32_t max_deletion) {
        if (position < kFlank) return -1;
        uint32_t const first = position - kFlank, end = position + kCore + kFlank;  // the k-mer: [first, end)
        if (allele.begin > first || allele.end < end) return -1;
        uint32_t const from = first > max_deletion ? first - max_deletion : 0;
        auto it = std::lower_bound(allele.edits.begin(), allele.edits.end(), from,
                                   [](strain_alleles::Edit const& e, uint32_t p) { return e.pos < p; });
        unsigned flank_substitutions = 0;
        for (; it != allele.edits.end() && it->pos < end; ++it) {
            switch (it->GetKind()) {
                case strain_alleles::kSubstitution:
                    if (it->pos < first) break;
                    if (it->pos >= position && it->pos < position + kCore) return 0;
                    if (++flank_substitutions > kMaxFlankSubstitutions) return 0;
                    break;
                case strain_alleles::kInsertion:  // bases inserted before pos: inside the k-mer if pos is past its first base
                    if (it->pos > first) return 0;
                    break;
                default:  // deletion of [pos, pos + length)
                    if (static_cast<uint32_t>(it->pos) + it->Length() > first) return 0;
                    break;
            }
        }
        return 1;
    }

    // The prior of a core without alleles: the mean within code of its 15 columns (kNoCode, 0, left out), at most
    // kPriorCode, over kPriorCode; kUnknownPrior if none has a code. `within`: the copy's codes by position, or empty.
    inline double Prior(std::span<uint8_t const> within, uint32_t position) {
        double sum = 0;
        unsigned n = 0;
        for (uint32_t p = position; p < position + kCore && p < within.size(); p++) {
            if (within[p] == 0) continue;
            sum += std::min(within[p], kPriorCode);
            n++;
        }
        return n == 0 ? kUnknownPrior : sum / (n * static_cast<double>(kPriorCode));
    }

    // A core's conservation: `kept` of `covering` alleles keep it a hit, pooled with the prior.
    inline double Conservation(unsigned kept, unsigned covering, double prior) {
        return (kept + kPriorWeight * prior) / (covering + kPriorWeight);
    }

    // A core's score before its place in a window: conservation less the ubiquity of its key.
    inline double Score(double conservation, uint32_t entries) {
        double const ubiquity = entries > 1 ? std::log2(static_cast<double>(entries)) / kMaxEntriesLog2 : 0.0;
        return conservation - kUbiquityWeight * ubiquity;
    }

    // Each core's level (levels, resized), from the scores of a gene's cores in gene order; `eligible[i]` 0 for a core
    // the index leaves out, which wins nothing (level 0, as every core that loses its first round). Between two cores
    // as good, the earlier one.
    inline void Levels(std::span<double const> scores, std::span<uint8_t const> eligible, std::vector<uint8_t>& levels,
                       std::vector<int64_t>& winners) {
        size_t const n = scores.size();
        levels.assign(n, 0);
        winners.resize(n);
        for (size_t i = 0; i < n; i++) winners[i] = eligible[i] ? static_cast<int64_t>(i) : -1;
        // Round l: the windows of 2^l ordinals, each the better of its halves' winners (winners[w] for the window w of
        // the round before, compacted in place). A gene of fewer than 2^l cores is one window, whose winner goes on to
        // the last level: a run gated at any C keeps its best core.
        size_t windows = n;
        for (unsigned l = 1; l < kLevels && n > 0; l++) {
            size_t const size = size_t{1} << l;
            size_t const next = (windows + 1) / 2;
            for (size_t w = 0; w < next; w++) {
                int64_t const a = winners[2 * w], b = 2 * w + 1 < windows ? winners[2 * w + 1] : -1;
                double const centre = static_cast<double>(w * size) + (static_cast<double>(size) - 1) / 2;
                auto const value = [&](int64_t i) {
                    return scores[static_cast<size_t>(i)] - kCentreWeight * std::abs(static_cast<double>(i) - centre) / size;
                };
                int64_t best = a < 0 ? b : b < 0 ? a : value(b) > value(a) ? b : a;
                winners[w] = best;
                if (best >= 0) levels[static_cast<size_t>(best)] = static_cast<uint8_t>(l);
            }
            windows = next;
        }
    }

    // What a build's grid made of the genes, for its log: cores, those alleles cover, and per level the cores of exactly
    // that level and their conservation summed (the cores kept at C are those of level log2(C) or more).
    struct Stats {
        uint64_t cores = 0, covered = 0;
        uint64_t by_level[kLevels] = {};
        double conservation[kLevels] = {};
        void Join(Stats const& o) {
            cores += o.cores;
            covered += o.covered;
            for (unsigned l = 0; l < kLevels; l++) {
                by_level[l] += o.by_level[l];
                conservation[l] += o.conservation[l];
            }
        }
        // The cores kept at level l (l or more), and their mean conservation.
        uint64_t Kept(unsigned l) const {
            uint64_t n = 0;
            for (unsigned k = l; k < kLevels; k++) n += by_level[k];
            return n;
        }
        double MeanConservation(unsigned l) const {
            double sum = 0;
            for (unsigned k = l; k < kLevels; k++) sum += conservation[k];
            uint64_t const n = Kept(l);
            return n ? sum / static_cast<double>(n) : 0.0;
        }
    };

    // Every C-th core of a gene at each C = 2^l (the cores whose ordinal is p_l mod 2^l), the phase p_l chosen per level by
    // the summed scores of the cores it keeps (those the index leaves out count nothing), between the two that keep the
    // cores of level l - 1 nested (p_{l-1} and p_{l-1} + 2^{l-1}); a phase that keeps no core of the gene loses, so every
    // gene keeps a core at every C, and between two as good, the first. The stride's spacing, C cores, with conservation
    // choosing where it falls (Choice::Phase).
    inline void PhaseLevels(std::span<double const> scores, std::span<uint8_t const> eligible, std::vector<uint8_t>& levels) {
        size_t const n = scores.size();
        levels.assign(n, 0);
        size_t phase = 0;
        for (unsigned l = 1; l < kLevels && n > 0; l++) {
            size_t const step = size_t{1} << l;
            size_t best = phase;
            double best_sum = 0;
            bool found = false;
            for (size_t const p : { phase, phase + step / 2 }) {
                if (p >= n) continue;
                double sum = 0;
                for (size_t i = p; i < n; i += step) sum += eligible[i] ? scores[i] : 0.0;
                if (!found || sum > best_sum) {
                    best = p;
                    best_sum = sum;
                    found = true;
                }
            }
            phase = best;
            for (size_t i = phase; i < n; i += step) levels[i] = static_cast<uint8_t>(l);
        }
    }

    // How each window's core is chosen (--grid_choice, a developer option of --build): Scored, the best core of each
    // window (Levels: one per window of C cores, anywhere in it); Even, every C-th core of the gene whatever the scores
    // (Levels with equal scores keeps each window's middle core: a stride on the gene, its cores C apart); Phase, every
    // C-th core with conservation choosing the phase (PhaseLevels).
    enum class Choice : uint8_t { Scored, Even, Phase };

    inline char const* ChoiceName(Choice c) { return c == Choice::Even ? "even" : c == Choice::Phase ? "phase" : "scored"; }

    // Reused buffers of GeneLevels, one per thread.
    struct Scratch {
        std::vector<double> scores, conservation, equal;
        std::vector<uint8_t> eligible;
        std::vector<int64_t> winners;
        std::vector<uint32_t> max_deletion;
    };

    // The levels of a gene's cores (in gene order) from its alleles (the copy's, strain_alleles::Table::Of and Get; none:
    // the prior alone) and its within codes by position (empty: none), counted into `stats`, chosen as `choice` says.
    inline void GeneLevels(std::span<Core const> cores, std::span<strain_alleles::Table::View const> alleles,
                           std::span<uint8_t const> within, std::vector<uint8_t>& levels, Scratch& scratch, Stats& stats,
                           Choice choice = Choice::Scored) {
        size_t const n = cores.size();
        scratch.scores.resize(n);
        scratch.conservation.resize(n);
        scratch.eligible.resize(n);
        scratch.max_deletion.assign(alleles.size(), 0);
        for (size_t a = 0; a < alleles.size(); a++) {
            for (auto const& e : alleles[a].edits) {
                if (e.GetKind() == strain_alleles::kDeletion) scratch.max_deletion[a] = std::max<uint32_t>(scratch.max_deletion[a], e.Length());
            }
        }
        for (size_t i = 0; i < n; i++) {
            unsigned kept = 0, covering = 0;
            for (size_t a = 0; a < alleles.size(); a++) {
                int const k = AlleleKeeps(alleles[a], cores[i].position, scratch.max_deletion[a]);
                if (k < 0) continue;
                covering++;
                kept += static_cast<unsigned>(k);
            }
            double const c = Conservation(kept, covering, Prior(within, cores[i].position));
            scratch.conservation[i] = c;
            scratch.scores[i] = Score(c, cores[i].entries);
            scratch.eligible[i] = cores[i].entries > 0;
            stats.covered += covering > 0;
        }
        if (choice == Choice::Phase) {
            PhaseLevels(scratch.scores, scratch.eligible, levels);
        } else if (choice == Choice::Even) {
            scratch.equal.assign(n, 0.0);
            Levels(scratch.equal, scratch.eligible, levels, scratch.winners);
        } else {
            Levels(scratch.scores, scratch.eligible, levels, scratch.winners);
        }
        stats.cores += n;
        for (size_t i = 0; i < n; i++) {
            stats.by_level[levels[i]]++;
            stats.conservation[levels[i]] += scratch.conservation[i];
        }
    }
}
