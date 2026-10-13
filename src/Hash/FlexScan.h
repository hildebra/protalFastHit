// FlexScan.h - scoring a k-mer core's flex cells against a read's: the inner loop of every seed lookup
// (KmerLookupSM::GetFromLookup). A lookup finds the block of entries of the read k-mer's core and scores each entry's
// flex cell (the 16 bases around the core, 2 bits each) by the bases it shares with the read's (Seedmap::Similarity);
// the entries of the best score are the lookup's seeds. At GTDB r226 a paired-end sample made 2.2G lookups of 70 cells
// each, in a scalar loop that was a large part of the seeding (docs/claude/2026-10-04-performance-gtdb-scale).
//
// ScoreScalar writes every cell's score (0-16) to a byte array and returns the best score and how many cells have it.
// With AVX2, BestAvx2 scores 8 cells per step (the same scores) and returns the best, and TiesAvx2 finds the cells with
// it, 32 scores per step, as bit masks: a scalar walk over the scores for the best cells was a third of a lookup's
// instructions at r226's block sizes (docs/claude/2026-10-06-performance-profiling). With AVX-512, ScanAvx512 does both,
// 16 cells a step, with the instructions AVX2 lacks (VPOPCNTD, ternary logic, funnel shifts, masks): ~9 instructions
// where BestAvx2 takes ~26 for 8 cells (docs/claude/2026-10-08-avx512-assessment.md). Lookups use the CPU's highest
// level (Kernel(); the cap PROTAL_SIMD sets, Utilities/SimdLevel.h); PROTAL_FLEX_SCAN=scalar keeps the scalar scan.
//
// With --flex_indels (kIndels; KmerLookupSM::SetFlexIndels) every kernel scores a cell allowing one insertion or deletion
// in each flank (IndelSimilarity): ONT reads' errors are mostly indels, and an indel shifts the flank's bases beyond it,
// which then match by chance, so one next to the core cost ~6 of the flank's 8 bases. The same kernels, a few times the
// instructions per cell (AVX-512 ~36 for 16 cells against ~9); the ties are found as without it.
#pragma once

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <string_view>
#include <utility>
#include <immintrin.h>

#include "Seedmap.h"
#include "SimdLevel.h"

namespace protal::flex_scan {
    // A flex key or cell (Seedmap::FlexKey): the left flank's 8 bases in bits 16-31 (the k-mer's first base highest, the
    // base next to the core in bits 16-17), the right flank's in bits 0-15 (the base next to the core in bits 14-15, the
    // k-mer's last in bits 0-1). Out from the core is up in the left half and down in the right. kBasePairs: the low bit
    // of every base's two.
    constexpr uint32_t kBasePairs = 0x55555555u, kLeftPairs = 0x55550000u, kRightPairs = 0x00005555u;

    // The bases of a and b that are equal, as the low bit of each base's pair (Seedmap::Similarity counts them).
    constexpr uint32_t EqualBases(uint32_t a, uint32_t b) {
        uint32_t const same = ~(a ^ b);
        return same & (same >> 1) & kBasePairs;
    }

    // The read's flex key and its flanks shifted by one base, for IndelSimilarity. Out from the core, deletion holds at
    // each place the read's base one place nearer the core (a cell base the read lacks puts the read's later bases one
    // place in), insertion the base one place further out (a base the read has extra puts them one place out). A shift
    // brings in at each flank's end a base the key does not hold: deletion's next to the core is never compared (a
    // deletion is past the first differing base), insertion's outermost (the read's base past the k-mer) is not known,
    // and kInsertionKnown leaves it out.
    struct IndelKeys {
        uint32_t same = 0, deletion = 0, insertion = 0;
    };
    constexpr uint32_t kInsertionKnown = 0x15555554u;

    constexpr IndelKeys MakeIndelKeys(uint32_t key) {
        uint32_t const left = key >> 16, right = key & 0xFFFFu;
        return { key, (((left << 2) & 0xFFFFu) << 16) | (right >> 2), ((left >> 2) << 16) | ((right << 2) & 0xFFFFu) };
    }

    // The bases (0-16) a cell shares with the read's flex key allowing one insertion or deletion in each flank: per flank
    // the most of (a) its bases equal to the cell's, as Similarity counts them, and (b) from its first base that differs
    // (out from the core) on, the read's bases one place in (a deletion: that cell base has none) or one place out (an
    // insertion), the bases before as they are. A single indel then costs one base, as a mismatch does: a flank with a
    // difference scores at most 7 (a deletion leaves a cell base without a read base, an insertion the outermost one
    // uncompared), so only a cell equal to the key scores 16. Two errors in a flank score as the shift of the first fits.
    inline uint32_t IndelSimilarity(uint32_t cell, IndelKeys const& keys) {
        uint32_t const equal = EqualBases(cell, keys.same);
        uint32_t const deletion = EqualBases(cell, keys.deletion);
        uint32_t const insertion = EqualBases(cell, keys.insertion) & kInsertionKnown;
        uint32_t const left = ~equal & kLeftPairs, right = ~equal & kRightPairs;  // each flank's differing bases
        // Each flank's first differing base, and the bases from it out (none where the flank has none): the left
        // flank's lowest bit, the right flank's highest.
        uint32_t const first_left = left & (0u - left);
        uint32_t const from_left = 0u - first_left;
        uint32_t const from_right = right == 0 ? 0u : 0xFFFFFFFFu >> __builtin_clz(right);
        uint32_t const from = from_left | from_right;
        uint32_t const beyond = (from_left ^ first_left) | (from_right >> 1);  // from it out, without it
        uint32_t const inner = equal & ~from;
        uint32_t const with_deletion = inner | (deletion & beyond);
        uint32_t const with_insertion = inner | (insertion & from);
        auto flank = [&](int shift) {
            auto bases = [shift](uint32_t x) { return static_cast<uint32_t>(__builtin_popcount((x >> shift) & 0xFFFFu)); };
            return std::max({ bases(equal), bases(with_deletion), bases(with_insertion) });
        };
        return flank(16) + flank(0);
    }

    // (best score, cells with it) of the block's flex cells against key, each cell's score in scores[0, size). The
    // best score starts at 0: a block whose cells share no base with key has best 0, every cell with it. kIndels: the
    // scores of IndelSimilarity.
    template <bool kIndels = false>
    inline std::pair<uint32_t, uint32_t> ScoreScalar(Seedmap::PackedBlock const& block, uint32_t key, uint8_t* scores) {
        uint32_t best = 0, count = 0;
        [[maybe_unused]] IndelKeys const keys = MakeIndelKeys(key);
        for (uint32_t i = 0; i < block.size; i++) {
            uint32_t const score = kIndels ? IndelSimilarity(Seedmap::FlexCell(block, i), keys)
                                           : Seedmap::Similarity(Seedmap::FlexCell(block, i), key);
            scores[i] = static_cast<uint8_t>(score);
            if (score > best) {
                best = score;
                count = 0;
            }
            count += score == best;
        }
        return { best, count };
    }

    // The bytes BestAvx2 may write past a block's last score, and the bytes it reads past the block's last flex cell.
    constexpr size_t kScorePadding = 8;
    constexpr size_t kReadPastCells = 32;
    static_assert(Seedmap::kPackedPadding >= kReadPastCells, "the index's padding must hold BestAvx2's reads past the last block");

    // BestAvx2's steps, as functions of their own (a lambda in a function compiled for AVX2 is not compiled for it).
    // EqualBases on 8 cells: the bits of `known` where cell and key have the same base (ones: all bits set).
    __attribute__((target("avx2"))) inline __m256i EqualBasesAvx2(__m256i cells, __m256i key, __m256i ones, __m256i known) {
        __m256i const same = _mm256_xor_si256(_mm256_xor_si256(cells, key), ones);
        return _mm256_and_si256(_mm256_and_si256(_mm256_srli_epi32(same, 1), same), known);
    }

    // The set bits of each byte, by a 4-bit lookup table (popcount4) of each nibble.
    __attribute__((target("avx2"))) inline __m256i BasesPerByteAvx2(__m256i bits, __m256i popcount4, __m256i nibble) {
        __m256i const low = _mm256_shuffle_epi8(popcount4, _mm256_and_si256(bits, nibble));
        __m256i const high = _mm256_shuffle_epi8(popcount4, _mm256_and_si256(_mm256_srli_epi16(bits, 4), nibble));
        return _mm256_add_epi8(low, high);
    }

    // The best score of the block's flex cells against key (0 if no cell shares a base with it, or the block is empty),
    // every cell's score in scores[0, size) as ScoreScalar gives it. 8 cells per step, the last step masked: cell i is
    // the 32 bits at bit flex_shift (0-7) from block.flex + 4 i, so 8 cells are words i..i+7 shifted right with the low
    // bits of words i+1..i+8 shifted in, and the last step reads up to kReadPastCells bytes past the block's cells
    // (its entries follow them, and the index ends in Seedmap::kPackedPadding bytes). Lanes past the block score 0
    // and do not raise the best; scores must hold size + kScorePadding bytes (the lanes past the block write zeros).
    // kIndels: IndelSimilarity's scores, its steps on 8 cells at once (each flank's bases summed in its 16-bit half).
    template <bool kIndels = false>
    __attribute__((target("avx2"))) inline uint32_t BestAvx2(Seedmap::PackedBlock const& block, uint32_t key, uint8_t* scores) {
        uint32_t const size = block.size;
        __m256i const k = _mm256_set1_epi32(static_cast<int>(key));
        __m256i const ones = _mm256_set1_epi32(-1);
        __m256i const pairs = _mm256_set1_epi32(0x55555555);
        __m256i const nibble = _mm256_set1_epi8(0x0f);
        __m256i const popcount4 = _mm256_setr_epi8(0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
                                                   0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4);
        __m256i const byte_ones = _mm256_set1_epi8(1);
        __m256i const word_ones = _mm256_set1_epi16(1);
        // Byte 0 of each 32-bit lane to the lane's low 4 bytes (per 128-bit half).
        __m256i const gather = _mm256_setr_epi8(0, 4, 8, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                                0, 4, 8, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
        __m256i const lane = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
        __m128i const right = _mm_cvtsi32_si128(static_cast<int>(block.flex_shift));
        __m128i const left = _mm_cvtsi32_si128(static_cast<int>(32 - block.flex_shift));
        [[maybe_unused]] IndelKeys const keys = MakeIndelKeys(key);
        [[maybe_unused]] __m256i const k_deletion = _mm256_set1_epi32(static_cast<int>(keys.deletion));
        [[maybe_unused]] __m256i const k_insertion = _mm256_set1_epi32(static_cast<int>(keys.insertion));
        [[maybe_unused]] __m256i const insertion_known = _mm256_set1_epi32(static_cast<int>(kInsertionKnown));
        [[maybe_unused]] __m256i const left_pairs = _mm256_set1_epi32(static_cast<int>(kLeftPairs));
        [[maybe_unused]] __m256i const right_pairs = _mm256_set1_epi32(static_cast<int>(kRightPairs));
        [[maybe_unused]] __m256i const zero = _mm256_setzero_si256();
        __m256i best8 = _mm256_setzero_si256();
        for (uint32_t i = 0; i < size; i += 8) {
            uint8_t const* p = block.flex + 4 * static_cast<size_t>(i);
            __m256i const a = _mm256_loadu_si256(reinterpret_cast<__m256i const*>(p));
            __m256i const b = _mm256_loadu_si256(reinterpret_cast<__m256i const*>(p + 4));
            __m256i const cells = _mm256_or_si256(_mm256_srl_epi32(a, right), _mm256_sll_epi32(b, left));  // sll by 32: 0
            __m256i const inside = _mm256_cmpgt_epi32(_mm256_set1_epi32(static_cast<int>(size - i)), lane);  // lanes in the block
            __m256i per_cell;
            if constexpr (kIndels) {
                // IndelSimilarity's steps. The right flank's bases from its first differing one on: its highest bit
                // smeared down (the bits are those of base pairs, 2 apart, at most 14 below it).
                __m256i const equal = EqualBasesAvx2(cells, k, ones, pairs);
                __m256i const deletion = EqualBasesAvx2(cells, k_deletion, ones, pairs);
                __m256i const insertion = EqualBasesAvx2(cells, k_insertion, ones, insertion_known);
                __m256i const left_differ = _mm256_andnot_si256(equal, left_pairs), right_differ = _mm256_andnot_si256(equal, right_pairs);
                __m256i const first_left = _mm256_and_si256(left_differ, _mm256_sub_epi32(zero, left_differ));
                __m256i const from_left = _mm256_sub_epi32(zero, first_left);
                __m256i from_right = _mm256_or_si256(right_differ, _mm256_srli_epi32(right_differ, 2));
                from_right = _mm256_or_si256(from_right, _mm256_srli_epi32(from_right, 4));
                from_right = _mm256_or_si256(from_right, _mm256_srli_epi32(from_right, 8));
                __m256i const from = _mm256_or_si256(from_left, from_right);
                __m256i const beyond = _mm256_or_si256(_mm256_xor_si256(from_left, first_left), _mm256_srli_epi32(from_right, 2));
                __m256i const inner = _mm256_andnot_si256(from, equal);
                __m256i const with_deletion = _mm256_or_si256(inner, _mm256_and_si256(deletion, beyond));
                __m256i const with_insertion = _mm256_or_si256(inner, _mm256_and_si256(insertion, from));
                // Each flank's bases in its 16-bit half, the most of the three, the two flanks summed.
                __m256i const flanks_equal = _mm256_maddubs_epi16(BasesPerByteAvx2(equal, popcount4, nibble), byte_ones);
                __m256i const flanks_deletion = _mm256_maddubs_epi16(BasesPerByteAvx2(with_deletion, popcount4, nibble), byte_ones);
                __m256i const flanks_insertion = _mm256_maddubs_epi16(BasesPerByteAvx2(with_insertion, popcount4, nibble), byte_ones);
                __m256i const best_flanks = _mm256_max_epu16(flanks_equal, _mm256_max_epu16(flanks_deletion, flanks_insertion));
                per_cell = _mm256_and_si256(_mm256_madd_epi16(best_flanks, word_ones), inside);
            } else {
                // Equal bases: both bits of a 2-bit pair 0 in cell ^ key (Similarity).
                __m256i const per_byte = BasesPerByteAvx2(EqualBasesAvx2(cells, k, ones, pairs), popcount4, nibble);
                per_cell = _mm256_and_si256(_mm256_madd_epi16(_mm256_maddubs_epi16(per_byte, byte_ones), word_ones), inside);
            }
            best8 = _mm256_max_epu32(best8, per_cell);
            __m256i const packed = _mm256_shuffle_epi8(per_cell, gather);
            uint32_t const first = static_cast<uint32_t>(_mm256_extract_epi32(packed, 0));
            uint32_t const second = static_cast<uint32_t>(_mm256_extract_epi32(packed, 4));
            std::memcpy(scores + i, &first, 4);
            std::memcpy(scores + i + 4, &second, 4);
        }
        __m128i m = _mm_max_epu32(_mm256_castsi256_si128(best8), _mm256_extracti128_si256(best8, 1));
        m = _mm_max_epu32(m, _mm_shuffle_epi32(m, 0x4e));
        m = _mm_max_epu32(m, _mm_shuffle_epi32(m, 0xb1));
        return static_cast<uint32_t>(_mm_cvtsi128_si32(m));
    }

    // Words of tie masks for a block of size cells, and the score bytes TiesAvx2 reads (whole 32-byte steps).
    constexpr size_t TieWords(size_t size) { return (size + 31) / 32; }
    constexpr size_t TieScoreBytes(size_t size) { return 32 * TieWords(size); }

    // The cells of scores[0, size) with score best, as bit masks of 32 cells each (bit j of masks[w]: cell 32 w + j;
    // the bits past size are 0), in masks[0, TieWords(size)); returns how many there are. Reads TieScoreBytes(size)
    // bytes of scores (the bytes past size may hold anything).
    __attribute__((target("avx2"))) inline uint32_t TiesAvx2(uint8_t const* scores, uint32_t size, uint32_t best, uint32_t* masks) {
        __m256i const target = _mm256_set1_epi8(static_cast<char>(best));
        uint32_t const words = static_cast<uint32_t>(TieWords(size));
        uint32_t count = 0;
        for (uint32_t w = 0; w < words; w++) {
            __m256i const s = _mm256_loadu_si256(reinterpret_cast<__m256i const*>(scores + 32 * static_cast<size_t>(w)));
            uint32_t m = static_cast<uint32_t>(_mm256_movemask_epi8(_mm256_cmpeq_epi8(s, target)));
            uint32_t const left = size - 32 * w;
            if (left < 32) m &= (1u << left) - 1;
            masks[w] = m;
            count += static_cast<uint32_t>(__builtin_popcount(m));
        }
        return count;
    }

    // EqualBases on 16 cells (the bits of `known` where cell and key have the same base), for ScanAvx512's kIndels.
    PROTAL_TARGET_AVX512 inline __m512i EqualBasesAvx512(__m512i cells, __m512i key, __m512i known) {
        __m512i const same = _mm512_ternarylogic_epi32(cells, key, key, 0xC3);
        return _mm512_ternarylogic_epi32(_mm512_srli_epi32(same, 1), same, known, 0x80);
    }

    // ScoreScalar, BestAvx2 and TiesAvx2 in one, with AVX-512 (the level of simd::Level::avx512): every cell's score
    // in scores[0, size) (nothing written past it), the cells with the best score as bit masks in masks[0,
    // TieWords(size)) (the bits past size 0), and (best score, cells with it). 16 cells a step: cell i is the funnel
    // shift of words i and i+1 by flex_shift (VBMI2), the equal base pairs by two ternary-logic steps, their count by
    // VPOPCNTD; the last step's lanes past the block are masked off, in the loads too, so nothing past the block's
    // cells (and the 4 bytes after them that FlexCell also reads) is read. The ties then 64 scores a step, by mask
    // compares. Not inlined into its caller (which is not compiled for AVX-512): one call per lookup.
    // kIndels: IndelSimilarity's scores, each flank's bases counted in its 16-bit half (VPOPCNTW), the right flank's
    // first differing base found by VPLZCNTD.
    template <bool kIndels = false>
    PROTAL_TARGET_AVX512 inline std::pair<uint32_t, uint32_t> ScanAvx512(Seedmap::PackedBlock const& block, uint32_t key,
                                                                       uint8_t* scores, uint32_t* masks) {
        uint32_t const size = block.size;
        __m512i const k = _mm512_set1_epi32(static_cast<int>(key));
        __m512i const shift = _mm512_set1_epi32(static_cast<int>(block.flex_shift));
        __m512i const pairs = _mm512_set1_epi32(0x55555555);
        [[maybe_unused]] IndelKeys const keys = MakeIndelKeys(key);
        [[maybe_unused]] __m512i const k_deletion = _mm512_set1_epi32(static_cast<int>(keys.deletion));
        [[maybe_unused]] __m512i const k_insertion = _mm512_set1_epi32(static_cast<int>(keys.insertion));
        [[maybe_unused]] __m512i const insertion_known = _mm512_set1_epi32(static_cast<int>(kInsertionKnown));
        [[maybe_unused]] __m512i const left_pairs = _mm512_set1_epi32(static_cast<int>(kLeftPairs));
        [[maybe_unused]] __m512i const right_pairs = _mm512_set1_epi32(static_cast<int>(kRightPairs));
        [[maybe_unused]] __m512i const all_ones = _mm512_set1_epi32(-1);
        [[maybe_unused]] __m512i const word_ones = _mm512_set1_epi16(1);
        [[maybe_unused]] __m512i const zero = _mm512_setzero_si512();
        __m512i best16 = _mm512_setzero_si512();
        for (uint32_t i = 0; i < size; i += 16) {
            uint32_t const left = size - i;
            __mmask16 const lanes = left >= 16 ? __mmask16{ 0xFFFF } : static_cast<__mmask16>((1u << left) - 1);
            uint8_t const* p = block.flex + 4 * static_cast<size_t>(i);
            __m512i const a = _mm512_maskz_loadu_epi32(lanes, p);
            __m512i const b = _mm512_maskz_loadu_epi32(lanes, p + 4);
            __m512i const cells = _mm512_shrdv_epi32(a, b, shift);  // (b:a) >> shift, the low 32 bits
            __m512i per_cell;
            if constexpr (kIndels) {
                // IndelSimilarity's steps (BestAvx2's, 16 cells at once).
                __m512i const equal = EqualBasesAvx512(cells, k, pairs);
                __m512i const deletion = EqualBasesAvx512(cells, k_deletion, pairs);
                __m512i const insertion = EqualBasesAvx512(cells, k_insertion, insertion_known);
                __m512i const left_differ = _mm512_andnot_si512(equal, left_pairs);
                __m512i const right_differ = _mm512_andnot_si512(equal, right_pairs);
                __m512i const first_left = _mm512_and_si512(left_differ, _mm512_sub_epi32(zero, left_differ));
                __m512i const from_left = _mm512_sub_epi32(zero, first_left);
                // Every bit up to the right flank's highest differing one: all ones shifted right by its leading
                // zeros; by 32 (VPSRLVD's 0) where the flank has none.
                __m512i const from_right = _mm512_srlv_epi32(all_ones, _mm512_lzcnt_epi32(right_differ));
                __m512i const from = _mm512_or_si512(from_left, from_right);
                __m512i const beyond = _mm512_or_si512(_mm512_xor_si512(from_left, first_left), _mm512_srli_epi32(from_right, 1));
                __m512i const inner = _mm512_andnot_si512(from, equal);
                __m512i const with_deletion = _mm512_or_si512(inner, _mm512_and_si512(deletion, beyond));
                __m512i const with_insertion = _mm512_or_si512(inner, _mm512_and_si512(insertion, from));
                __m512i const best_flanks = _mm512_max_epu16(_mm512_popcnt_epi16(equal),
                                                             _mm512_max_epu16(_mm512_popcnt_epi16(with_deletion), _mm512_popcnt_epi16(with_insertion)));
                per_cell = _mm512_maskz_madd_epi16(lanes, best_flanks, word_ones);  // the two flanks summed
            } else {
                // Similarity: equal bases are 2-bit pairs with both bits 0 in cell ^ key. same = ~(cell ^ key) (0xC3: 1
                // where the first two inputs agree); equal = (same >> 1) & same & pairs (0x80: the AND of all three).
                __m512i const same = _mm512_ternarylogic_epi32(cells, k, k, 0xC3);
                __m512i const equal = _mm512_ternarylogic_epi32(_mm512_srli_epi32(same, 1), same, pairs, 0x80);
                per_cell = _mm512_maskz_popcnt_epi32(lanes, equal);
            }
            best16 = _mm512_max_epu32(best16, per_cell);
            _mm512_mask_cvtepi32_storeu_epi8(scores + i, lanes, per_cell);
        }
        uint32_t const best = _mm512_reduce_max_epu32(best16);
        __m512i const target = _mm512_set1_epi8(static_cast<char>(best));
        uint32_t const words = static_cast<uint32_t>(TieWords(size));
        uint32_t count = 0;
        for (uint32_t i = 0; i < size; i += 64) {
            uint32_t const left = size - i;
            __mmask64 const lanes = left >= 64 ? ~__mmask64{ 0 } : static_cast<__mmask64>(_bzhi_u64(~uint64_t{ 0 }, left));
            uint64_t const bits = _cvtmask64_u64(_mm512_mask_cmpeq_epi8_mask(lanes, _mm512_maskz_loadu_epi8(lanes, scores + i), target));
            masks[i / 32] = static_cast<uint32_t>(bits);
            if (i / 32 + 1 < words) masks[i / 32 + 1] = static_cast<uint32_t>(bits >> 32);
            count += static_cast<uint32_t>(__builtin_popcountll(bits));
        }
        return { best, count };
    }

    // The version lookups use: simd::Default() (the CPU's highest level, capped by PROTAL_SIMD), scalar if
    // PROTAL_FLEX_SCAN=scalar.
    inline std::atomic<simd::Level>& Kernel() {
        static std::atomic<simd::Level> level{ [] {
            char const* choice = std::getenv("PROTAL_FLEX_SCAN");
            if (choice && std::string_view(choice) == "scalar") return simd::Level::scalar;
            return simd::Default();
        }() };
        return level;
    }

    // For tests and benchmarks: lookups use `wanted`, or the highest level below it the CPU supports, which is returned.
    inline simd::Level Use(simd::Level wanted) {
        simd::Level const level = simd::Supported(wanted);
        Kernel().store(level, std::memory_order_relaxed);
        return level;
    }
}
