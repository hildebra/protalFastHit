// AlignmentScreen.h - a k-mer screen before WFA2. A candidate alignment of a read into its gene
// window is refused without aligning when the read and the window share too few k-mers for any
// alignment within the score budget to exist: the q-gram lemma. The candidate fails exactly as
// WFA2 would have failed it (no alignment under the budget exists, so none can be found), and its
// taxon stays one of the read's failed candidates.
//
// The bound: an alignment of score at most B touches (changes a base of, or inserts into) at most
// B / c of the read's k-mers, where c is the smallest score an operation pays per k-mer it touches:
// a mismatch (score `mismatch`) touches at most k k-mers; an insertion run of g read bases
// (`gap_opening` + g * `gap_extension`) at most k - 1 + g; a deletion run at most k - 1. Every
// other k-mer of the read's aligned part is intact: k read bases aligned to k equal, consecutive
// window bases, so it occurs in the window. With free ends, the alignment may leave up to
// `begin_free` read bases before it and `end_free` after it unaligned; the read bases between
// are aligned by every alignment, and only their k-mers count.
//
// On GTDB-scale databases most candidates of real reads fail (docs/claude/2026-10-04-performance-gtdb-scale):
// a failing WFA2 alignment explores its whole budget, where the screen costs a pass over the read
// and the window.
#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

#include "SequenceUtils/KmerUtils.h"

namespace protal {
    class AlignmentScreen {
    public:
        static constexpr size_t kMaxK = 10;    // the stamps are indexed by 2-bit codes of up to 20 bits
        static constexpr size_t kDefaultK = 6;  // for windows of up to kShortWindow bases; longer windows take kLongK
        static constexpr size_t kShortWindow = 512;
        static constexpr size_t kLongK = 7;

        // The k for a window: measured on random reads of 5-25% divergence (test AlignmentScreen.RefusesOnlyWhatWFA2Fails),
        // k = 6 refused most of the failures of 150 bp reads (47%; k = 8: 19%) and k = 7 most of those of 500-1500
        // base windows (39%): a smaller k leaves more k-mers intact per edit but matches more by chance in a long window.
        static size_t KFor(size_t window_length) { return window_length <= kShortWindow ? kDefaultK : kLongK; }

        // The most k-mers an alignment of score at most `max_score - 1` (an alignment whose score reaches max_score
        // is abandoned, WFA2Wrapper2::Alignment) can touch, with the gap-affine penalties given.
        static int64_t MaxTouched(int max_score, int mismatch, int gap_opening, int gap_extension, size_t k_mer = kDefaultK) {
            int64_t const budget = static_cast<int64_t>(max_score) - 1;
            if (budget < 0) return -1;  // nothing aligns
            int64_t const k = static_cast<int64_t>(k_mer);
            int64_t touched = 0;
            if (mismatch > 0) touched = std::max(touched, budget * k / mismatch);
            int64_t const open = static_cast<int64_t>(gap_opening) + static_cast<int64_t>(gap_extension);
            if (open > 0) touched = std::max(touched, budget * k / open);
            if (gap_extension > 0) touched = std::max(touched, budget / gap_extension);
            return touched;
        }

        // The fewest k-mers of `aligned` read bases that occur in the window if an alignment within the budget
        // exists; 0 or less: the screen cannot refuse.
        static int64_t Required(size_t aligned, int max_score, int mismatch, int gap_opening, int gap_extension, size_t k = kDefaultK) {
            if (aligned < k) return 0;
            int64_t const touched = MaxTouched(max_score, mismatch, gap_opening, gap_extension, k);
            if (touched < 0) return static_cast<int64_t>(aligned);  // more than there are: refused
            return static_cast<int64_t>(aligned - k + 1) - touched;
        }

        // Whether an alignment of `read` into `window` of score below max_score may exist, with up to begin_free
        // read bases before and end_free after it left unaligned: true unless the read's aligned part shares
        // fewer k-mers with the window than any such alignment leaves intact. A read k-mer with a base other
        // than ACGT counts as shared (WFA2 matches equal characters, whatever they are); window k-mers with such
        // a base match nothing.
        bool MayAlign(std::string_view read, size_t begin_free, size_t end_free, std::string_view window,
                      int max_score, int mismatch, int gap_opening, int gap_extension, size_t k = 0) {
            if (begin_free + end_free >= read.size()) return true;
            if (k == 0) k = KFor(window.size());
            k = std::min(k, kMaxK);
            size_t const aligned = read.size() - begin_free - end_free;
            int64_t const required = Required(aligned, max_score, mismatch, gap_opening, gap_extension, k);
            if (required <= 0) return true;
            if (required > static_cast<int64_t>(aligned - k + 1)) return false;
            uint32_t const mask = (1u << (2 * k)) - 1;

            Stamp(window, k, mask);
            int64_t shared = 0;
            int64_t remaining = static_cast<int64_t>(aligned - k + 1);  // k-mers not yet counted
            uint32_t code = 0;
            size_t valid = 0;  // ACGT bases in a row ending at the current base
            for (size_t i = begin_free, end = read.size() - end_free; i < end; i++) {
                uint32_t const base = static_cast<uint32_t>(KmerUtils::BaseToInt(read[i]));
                if (base > 3) valid = 0; else { valid++; code = ((code << 2) | base) & mask; }
                if (i + 1 < begin_free + k) continue;
                remaining--;
                // valid < k: a k-mer with a non-ACGT base, counted as shared.
                if (valid < k || m_stamps[code] == m_generation) shared++;
                if (shared >= required) return true;
                if (shared + remaining < required) return false;
            }
            return shared >= required;
        }

        // Candidates refused so far.
        size_t Refused() const { return m_refused; }
        void CountRefused() { m_refused++; }

    private:
        // Marks the window's k-mers with a fresh generation stamp.
        void Stamp(std::string_view window, size_t k, uint32_t mask) {
            if (m_stamps.size() < size_t{1} << (2 * k)) { m_stamps.assign(size_t{1} << (2 * k), 0); m_generation = 0; }
            if (++m_generation == 0) {  // wrapped: every old stamp reads as 0, so start over at 1
                std::fill(m_stamps.begin(), m_stamps.end(), 0);
                m_generation = 1;
            }
            uint32_t code = 0;
            size_t valid = 0;
            for (char const c : window) {
                uint32_t const base = static_cast<uint32_t>(KmerUtils::BaseToInt(c));
                if (base > 3) { valid = 0; continue; }
                valid++;
                code = ((code << 2) | base) & mask;
                if (valid >= k) m_stamps[code] = m_generation;
            }
        }

        std::vector<uint16_t> m_stamps;  // per k-mer code, the generation that last saw it in a window; sized on first use
        uint16_t m_generation = 0;
        size_t m_refused = 0;
    };
}
