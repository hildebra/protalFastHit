// FlexScan.h - scoring a k-mer core's flex cells against a read's: the inner loop of every seed lookup
// (KmerLookupSM::GetFromLookup). A lookup finds the block of entries of the read k-mer's core and scores each entry's
// flex cell (the 16 bases around the core, 2 bits each) by the bases it shares with the read's (Seedmap::Similarity);
// the entries of the best score are the lookup's seeds. At GTDB r226 a paired-end sample made 2.2G lookups of 70 cells
// each, in a scalar loop that was a large part of the seeding (docs/claude/2026-10-04-performance-gtdb-scale).
//
// Score() writes every cell's score (0-16) to a byte array and returns the best score and how many cells have it; the
// AVX2 version scores 8 cells per step, with the same scores. It is chosen at run time where the CPU has AVX2;
// PROTAL_FLEX_SCAN=scalar in the environment keeps the scalar one (for comparing runs).
#pragma once

#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <string_view>
#include <utility>
#include <immintrin.h>

#include "Seedmap.h"

namespace protal::flex_scan {
    // (best score, cells with it) of the block's flex cells against key, each cell's score in scores[0, size). The
    // best score starts at 0: a block whose cells share no base with key has best 0, every cell with it.
    inline std::pair<uint32_t, uint32_t> ScoreScalar(Seedmap::PackedBlock const& block, uint32_t key, uint8_t* scores) {
        uint32_t best = 0, count = 0;
        for (uint32_t i = 0; i < block.size; i++) {
            uint32_t const score = Seedmap::Similarity(Seedmap::FlexCell(block, i), key);
            scores[i] = static_cast<uint8_t>(score);
            if (score > best) {
                best = score;
                count = 0;
            }
            count += score == best;
        }
        return { best, count };
    }

    // As ScoreScalar, 8 cells at a time. Cell i is the 32 bits at bit flex_shift (0-7) from block.flex + 4 i, so 8 cells
    // are words i..i+7 shifted right with the low bits of words i+1..i+8 shifted in; a block's cells are followed by
    // its entries (and the index by 16 bytes of padding), so the extra word is always inside the index, as the scalar
    // FlexCell's 8-byte read of the last cell is.
    __attribute__((target("avx2"))) inline std::pair<uint32_t, uint32_t> ScoreAvx2(Seedmap::PackedBlock const& block, uint32_t key,
                                                                                    uint8_t* scores) {
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
        __m128i const right = _mm_cvtsi32_si128(static_cast<int>(block.flex_shift));
        __m128i const left = _mm_cvtsi32_si128(static_cast<int>(32 - block.flex_shift));
        __m256i best8 = _mm256_setzero_si256();
        uint32_t i = 0;
        for (; i + 8 <= size; i += 8) {
            uint8_t const* p = block.flex + 4 * static_cast<size_t>(i);
            __m256i const a = _mm256_loadu_si256(reinterpret_cast<__m256i const*>(p));
            __m256i const b = _mm256_loadu_si256(reinterpret_cast<__m256i const*>(p + 4));
            __m256i const cells = _mm256_or_si256(_mm256_srl_epi32(a, right), _mm256_sll_epi32(b, left));  // sll by 32: 0
            // Equal bases: both bits of a 2-bit pair 0 in cell ^ key (Similarity).
            __m256i const same = _mm256_xor_si256(_mm256_xor_si256(cells, k), ones);
            __m256i const equal = _mm256_and_si256(_mm256_and_si256(_mm256_srli_epi32(same, 1), same), pairs);
            __m256i const low = _mm256_shuffle_epi8(popcount4, _mm256_and_si256(equal, nibble));
            __m256i const high = _mm256_shuffle_epi8(popcount4, _mm256_and_si256(_mm256_srli_epi16(equal, 4), nibble));
            __m256i const per_byte = _mm256_add_epi8(low, high);
            __m256i const per_cell = _mm256_madd_epi16(_mm256_maddubs_epi16(per_byte, byte_ones), word_ones);  // 0-16 per lane
            best8 = _mm256_max_epu32(best8, per_cell);
            __m256i const packed = _mm256_shuffle_epi8(per_cell, gather);
            uint32_t const first = static_cast<uint32_t>(_mm256_extract_epi32(packed, 0));
            uint32_t const second = static_cast<uint32_t>(_mm256_extract_epi32(packed, 4));
            std::memcpy(scores + i, &first, 4);
            std::memcpy(scores + i + 4, &second, 4);
        }
        alignas(32) uint32_t lanes[8];
        _mm256_store_si256(reinterpret_cast<__m256i*>(lanes), best8);
        uint32_t best = 0;
        for (uint32_t lane : lanes) best = lane > best ? lane : best;
        for (; i < size; i++) {
            uint32_t const score = Seedmap::Similarity(Seedmap::FlexCell(block, i), key);
            scores[i] = static_cast<uint8_t>(score);
            best = score > best ? score : best;
        }
        // The cells with the best score, 32 scores at a time.
        __m256i const target = _mm256_set1_epi8(static_cast<char>(best));
        uint32_t count = 0, j = 0;
        for (; j + 32 <= size; j += 32) {
            __m256i const s = _mm256_loadu_si256(reinterpret_cast<__m256i const*>(scores + j));
            count += static_cast<uint32_t>(__builtin_popcount(static_cast<uint32_t>(_mm256_movemask_epi8(_mm256_cmpeq_epi8(s, target)))));
        }
        for (; j < size; j++) count += scores[j] == best;
        return { best, count };
    }

    inline bool CpuHasAvx2() {
        __builtin_cpu_init();
        return __builtin_cpu_supports("avx2");
    }

    // Whether Score uses ScoreAvx2: where the CPU has AVX2, unless PROTAL_FLEX_SCAN=scalar.
    inline std::atomic<bool>& Avx2Enabled() {
        static std::atomic<bool> enabled{ [] {
            char const* choice = std::getenv("PROTAL_FLEX_SCAN");
            return CpuHasAvx2() && !(choice && std::string_view(choice) == "scalar");
        }() };
        return enabled;
    }

    // For tests: the AVX2 version on (where the CPU has it) or off.
    inline void UseAvx2(bool use) { Avx2Enabled().store(use && CpuHasAvx2(), std::memory_order_relaxed); }

    inline std::pair<uint32_t, uint32_t> Score(Seedmap::PackedBlock const& block, uint32_t key, uint8_t* scores) {
        return Avx2Enabled().load(std::memory_order_relaxed) ? ScoreAvx2(block, key, scores) : ScoreScalar(block, key, scores);
    }
}
