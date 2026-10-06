// lookup_bench.cpp - seed lookups at GTDB r226's block sizes (the eleventh cluster run: 70.4 flex cells per block; by
// lookups 1.2% blocks of 1 entry, 30.1% of 2-15, 63.2% of 16-255, 5.5% of 256-4095; ~2 seeds per lookup), in an
// arena of packed blocks as Seedmap's packed layout holds them (flex cells of 32 bits, then entries of 41 bits, at bit
// offsets). Each "mate" makes 22 lookups of random blocks; a pair's 44 lookups have their first lines prefetched,
// then each mate's lookups are taken smallest first with the values 16 lookups ahead prefetched, as
// ChainAnchorFinder::PrepareLookups and FindSeeds do.
//
// Variants of KmerLookupSM::GetFromLookup:
//   A  flex_scan::ScoreAvx2 (as in the tree), then a scalar pass over the scores for the best ones (the tree's code)
//   B  ScoreAvx2, then the best ones by 32-byte compares and their bit masks
//   C  one AVX2 scoring pass whose last step is masked (no scalar tail) and no counting pass, then B's tie pass,
//      which also counts
//   D  C, with every line of a block's flex cells prefetched (up to 16 lines), not only the first
// Every variant must give the same seeds (a checksum of every field).
//
// Build (in WSL): g++ -O3 -std=c++20 -march=x86-64 -mtune=generic -I src -I lib lookup_bench.cpp -o lookup_bench
// Run: lookup_bench <arena GB> <pairs> <rounds> [variants]
#include <iostream>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <random>
#include <string>
#include <vector>
#include <sys/mman.h>
#include <immintrin.h>

#include "Hash/FlexScan.h"

using protal::Seedmap;
namespace flex_scan = protal::flex_scan;

namespace {
    constexpr unsigned kPosBits = 14, kTaxonGeneBits = 25, kEntryBits = kTaxonGeneBits + kPosBits + 2;  // r226
    constexpr uint64_t kGenes = 169;  // gene ids 0..168
    constexpr uint64_t kEntryMask = (uint64_t{1} << kEntryBits) - 1, kPosMask = (uint64_t{1} << kPosBits) - 1,
                       kTaxonGeneMask = (uint64_t{1} << kTaxonGeneBits) - 1;
    constexpr uint32_t kMaxUbiquity = 256, kFlexK = 16, kLookupsPerMate = 22;
    constexpr size_t kPrefetchLookups = 16;
    uint64_t const kReciprocal = static_cast<uint64_t>((static_cast<__uint128_t>(1) << 64) / kGenes + 1);

    struct Seed {
        uint32_t taxid, geneid, genepos;
        uint16_t readpos;
        bool unique, unique_two;
    };

    struct Lookup : Seedmap::PackedBlock {
        uint32_t flex_key = 0, read_pos = 0;
    };

    // Seedmap::EntryValue's arithmetic.
    inline uint64_t EntryValue(Seedmap::PackedBlock const& block, uint32_t i) {
        uint64_t const bit = block.entry_shift + static_cast<uint64_t>(i) * kEntryBits;
        __uint128_t x;
        std::memcpy(&x, block.entries + (bit >> 3), 16);
        uint64_t const v = static_cast<uint64_t>(x >> (bit & 7)) & kEntryMask;
        uint64_t const taxon_gene = (v >> kPosBits) & kTaxonGeneMask;
        uint64_t const taxid = static_cast<uint64_t>((static_cast<__uint128_t>(taxon_gene) * kReciprocal) >> 64);
        uint64_t const gene = taxon_gene - taxid * kGenes;
        uint64_t const flags = v >> (kPosBits + kTaxonGeneBits);
        return taxid << 40 | gene << 20 | (v & kPosMask) | flags << 60;
    }

    inline void Emit(std::vector<Seed>& out, Lookup const& l, uint32_t i, bool exact) {
        uint64_t const v = EntryValue(l, i);
        out.push_back({ static_cast<uint32_t>((v >> 40) & 0xfffff), static_cast<uint32_t>((v >> 20) & 0xfffff),
                        static_cast<uint32_t>(v & 0xfffff), static_cast<uint16_t>(l.read_pos + 8),
                        static_cast<bool>((v >> 60) & 1) && exact, static_cast<bool>((v >> 61) & 1) && exact });
    }

    std::vector<uint8_t> g_scores(8192 + 64);
    std::vector<uint32_t> g_masks(256 + 2);

    // A: the tree's GetFromLookup.
    void LookupA(std::vector<Seed>& out, Lookup const& l) {
        if (l.flex != nullptr) {
            auto const [best, count] = flex_scan::ScoreAvx2(l, l.flex_key, g_scores.data());
            if (count > kMaxUbiquity) return;
            for (uint32_t i = 0; i < l.size; i++) {
                if (g_scores[i] == best) Emit(out, l, i, best == kFlexK);
            }
        } else {
            for (uint32_t i = 0; i < l.size; i++) Emit(out, l, i, false);
        }
    }

    // The best cells by 32-byte compares; `count` known (A's) or counted here (count = UINT32_MAX).
    __attribute__((target("avx2"))) void TiesByMasks(std::vector<Seed>& out, Lookup const& l, uint32_t best, uint32_t count) {
        __m256i const target = _mm256_set1_epi8(static_cast<char>(best));
        uint32_t const size = l.size, words = (size + 31) / 32;
        uint32_t total = 0;
        for (uint32_t w = 0; w < words; w++) {
            __m256i const s = _mm256_loadu_si256(reinterpret_cast<__m256i const*>(g_scores.data() + 32 * w));
            uint32_t m = static_cast<uint32_t>(_mm256_movemask_epi8(_mm256_cmpeq_epi8(s, target)));
            uint32_t const left = size - 32 * w;
            if (left < 32) m &= (1u << left) - 1;
            g_masks[w] = m;
            total += static_cast<uint32_t>(__builtin_popcount(m));
        }
        if (count == UINT32_MAX) count = total;
        if (count > kMaxUbiquity) return;
        bool const exact = best == kFlexK;
        for (uint32_t w = 0; w < words; w++) {
            for (uint32_t m = g_masks[w]; m; m &= m - 1) Emit(out, l, 32 * w + static_cast<uint32_t>(__builtin_ctz(m)), exact);
        }
    }

    void LookupB(std::vector<Seed>& out, Lookup const& l) {
        if (l.flex != nullptr) {
            auto const [best, count] = flex_scan::ScoreAvx2(l, l.flex_key, g_scores.data());
            if (count > kMaxUbiquity) return;
            TiesByMasks(out, l, best, count);
        } else {
            for (uint32_t i = 0; i < l.size; i++) Emit(out, l, i, false);
        }
    }

    // ScoreAvx2 with the last (partial) step masked: lanes past the block score 0 and do not raise the best. Reads up to
    // 36 bytes past the block's last cell (the arena is padded). Returns the best score; the scores of lanes past the
    // block are written as 0 and never read as ties (TiesByMasks cuts the masks to the block).
    __attribute__((target("avx2"))) uint32_t ScoreMasked(Seedmap::PackedBlock const& block, uint32_t key, uint8_t* scores) {
        uint32_t const size = block.size;
        __m256i const k = _mm256_set1_epi32(static_cast<int>(key));
        __m256i const ones = _mm256_set1_epi32(-1);
        __m256i const pairs = _mm256_set1_epi32(0x55555555);
        __m256i const nibble = _mm256_set1_epi8(0x0f);
        __m256i const popcount4 = _mm256_setr_epi8(0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
                                                   0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4);
        __m256i const byte_ones = _mm256_set1_epi8(1);
        __m256i const word_ones = _mm256_set1_epi16(1);
        __m256i const gather = _mm256_setr_epi8(0, 4, 8, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                                0, 4, 8, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
        __m256i const lane_index = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
        __m128i const right = _mm_cvtsi32_si128(static_cast<int>(block.flex_shift));
        __m128i const left = _mm_cvtsi32_si128(static_cast<int>(32 - block.flex_shift));
        __m256i best8 = _mm256_setzero_si256();
        for (uint32_t i = 0; i < size; i += 8) {
            uint8_t const* p = block.flex + 4 * static_cast<size_t>(i);
            __m256i const a = _mm256_loadu_si256(reinterpret_cast<__m256i const*>(p));
            __m256i const b = _mm256_loadu_si256(reinterpret_cast<__m256i const*>(p + 4));
            __m256i const cells = _mm256_or_si256(_mm256_srl_epi32(a, right), _mm256_sll_epi32(b, left));
            __m256i const same = _mm256_xor_si256(_mm256_xor_si256(cells, k), ones);
            __m256i const equal = _mm256_and_si256(_mm256_and_si256(_mm256_srli_epi32(same, 1), same), pairs);
            __m256i const low = _mm256_shuffle_epi8(popcount4, _mm256_and_si256(equal, nibble));
            __m256i const high = _mm256_shuffle_epi8(popcount4, _mm256_and_si256(_mm256_srli_epi16(equal, 4), nibble));
            __m256i const per_byte = _mm256_add_epi8(low, high);
            __m256i per_cell = _mm256_madd_epi16(_mm256_maddubs_epi16(per_byte, byte_ones), word_ones);
            // Lanes at or past the block's end: 0.
            __m256i const inside = _mm256_cmpgt_epi32(_mm256_set1_epi32(static_cast<int>(size - i)), lane_index);
            per_cell = _mm256_and_si256(per_cell, inside);
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

    void LookupC(std::vector<Seed>& out, Lookup const& l) {
        if (l.flex != nullptr) {
            uint32_t const best = ScoreMasked(l, l.flex_key, g_scores.data());
            TiesByMasks(out, l, best, UINT32_MAX);
        } else {
            for (uint32_t i = 0; i < l.size; i++) Emit(out, l, i, false);
        }
    }

    inline void PrefetchFirst(Lookup const& l) {
        if (l.flex != nullptr) __builtin_prefetch(l.flex);
        __builtin_prefetch(l.entries);
    }

    inline void PrefetchAll(Lookup const& l) {
        if (l.flex != nullptr) {
            uintptr_t const begin = reinterpret_cast<uintptr_t>(l.flex) & ~uintptr_t{63};
            uintptr_t const end = std::min<uintptr_t>(reinterpret_cast<uintptr_t>(l.flex) + 4 * static_cast<uintptr_t>(l.size), begin + 16 * 64);
            for (uintptr_t a = begin; a < end; a += 64) __builtin_prefetch(reinterpret_cast<void const*>(a));
        }
        __builtin_prefetch(l.entries);
    }

    // ORs w bits of x into the zeroed bit array at bit.
    void PutBits(uint8_t* base, uint64_t bit, uint64_t x, unsigned w) {
        __uint128_t cur;
        std::memcpy(&cur, base + (bit >> 3), 16);
        cur |= static_cast<__uint128_t>(x & ((w == 64) ? ~uint64_t{0} : ((uint64_t{1} << w) - 1))) << (bit & 7);
        std::memcpy(base + (bit >> 3), &cur, 16);
    }

    // A block size by the r226 shares of lookups per group and the groups' mean sizes (8.4, 61.8, 523).
    struct SizeDist {
        double p_in_16_255 = 0, p_in_256_4095 = 0;  // exponents of the power-law-like draws, fitted to the means
        static uint32_t Draw(double u, double lo, double hi, double p) { return static_cast<uint32_t>(lo * std::pow(hi / lo, std::pow(u, p))); }
        static double Mean(double lo, double hi, double p) {
            double s = 0;
            for (int i = 0; i < 20000; i++) s += Draw((i + 0.5) / 20000, lo, hi, p);
            return s / 20000;
        }
        static double Fit(double lo, double hi, double mean) {
            double a = 0.2, b = 20;
            for (int it = 0; it < 60; it++) {
                double const m = std::sqrt(a * b);
                (Mean(lo, hi, m) > mean ? a : b) = m;  // a larger p gives smaller sizes
            }
            return std::sqrt(a * b);
        }
        SizeDist() { p_in_16_255 = Fit(16, 255.99, 61.8); p_in_256_4095 = Fit(256, 4095.99, 523); }
        uint32_t operator()(std::mt19937_64& rng) const {
            std::uniform_real_distribution<double> u(0, 1);
            double const g = u(rng);
            if (g < 0.012) return 1;
            if (g < 0.012 + 0.301) return 2 + static_cast<uint32_t>(rng() % 14);
            if (g < 0.012 + 0.301 + 0.632) return Draw(u(rng), 16, 255.99, p_in_16_255);
            return Draw(u(rng), 256, 4095.99, p_in_256_4095);
        }
    };

    uint32_t Similar(uint32_t key, int mismatches, std::mt19937_64& rng) {
        uint32_t cell = key;
        for (int j = 0; j < mismatches; j++) {
            unsigned const base = static_cast<unsigned>(rng() % 16);
            cell ^= (1u + static_cast<uint32_t>(rng() % 3)) << (2 * base);  // another base there
        }
        return cell;
    }
}

int main(int argc, char** argv) {
    double const arena_gb = argc > 1 ? std::atof(argv[1]) : 1.0;
    size_t const pairs = argc > 2 ? std::strtoull(argv[2], nullptr, 10) : 100000;
    int const rounds = argc > 3 ? std::atoi(argv[3]) : 3;
    std::string const variants = argc > 4 ? argv[4] : "ABCD";
    size_t const bytes = static_cast<size_t>(arena_gb * 1e9);

    // The arena: blocks one after another, each its flex cells (if it has 2+ entries) then its entries.
    uint8_t* arena = static_cast<uint8_t*>(mmap(nullptr, bytes + 4096, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0));
    if (arena == MAP_FAILED) { std::perror("mmap"); return 1; }
    madvise(arena, bytes + 4096, MADV_HUGEPAGE);
    std::mt19937_64 rng(11);
    SizeDist const dist;
    std::vector<Lookup> blocks;
    uint64_t bit = 0;
    double cells = 0, seeds_expected = 0;
    while (true) {
        uint32_t const size = dist(rng);
        uint64_t const need = (size >= 2 ? 32ull * size : 0) + uint64_t{kEntryBits} * size;
        if ((bit + need) / 8 + 64 >= bytes) break;
        Lookup l;
        l.size = size;
        l.flex_key = static_cast<uint32_t>(rng());
        l.read_pos = static_cast<uint32_t>(rng() % 120);
        if (size >= 2) {
            l.flex = arena + (bit >> 3);
            l.flex_shift = static_cast<uint32_t>(bit & 7);
            // 1-3 cells near the key (0-4 bases changed, the same number for each, so they tie), the rest random.
            uint32_t const ties = std::min<uint32_t>(size, 1 + static_cast<uint32_t>(rng() % 3));
            int const changed = static_cast<int>(rng() % 5);
            std::vector<uint32_t> cell(size);
            for (auto& c : cell) c = static_cast<uint32_t>(rng());
            for (uint32_t t = 0; t < ties; t++) cell[rng() % size] = Similar(l.flex_key, changed, rng);
            for (uint32_t i = 0; i < size; i++) PutBits(arena, bit + 32ull * i, cell[i], 32);
            bit += 32ull * size;
        }
        l.entries = arena + (bit >> 3);
        l.entry_shift = static_cast<uint32_t>(bit & 7);
        for (uint32_t i = 0; i < size; i++) {
            uint64_t const taxid = 1 + rng() % 143614, gene = rng() % kGenes, pos = rng() % 16384, flags = rng() % 4;
            PutBits(arena, bit, ((taxid * kGenes + gene) << kPosBits) | pos | flags << (kPosBits + kTaxonGeneBits), kEntryBits);
            bit += kEntryBits;
        }
        cells += size;
        blocks.push_back(l);
    }
    std::printf("arena %.2f GB: %zu blocks, %.1f cells per block\n", bytes / 1e9, blocks.size(), cells / blocks.size());

    // The lookups: per mate 22 random blocks.
    std::vector<uint32_t> order(pairs * 2 * kLookupsPerMate);
    for (auto& o : order) o = static_cast<uint32_t>(rng() % blocks.size());

    std::vector<Seed> seeds;
    std::vector<Lookup> mate;
    auto run = [&](char variant) {
        auto lookup = variant == 'A' ? LookupA : variant == 'B' ? LookupB : LookupC;
        bool const all_lines = variant == 'D';
        uint64_t checksum = 0, total_seeds = 0;
        auto const start = std::chrono::steady_clock::now();
        for (size_t p = 0; p < pairs; p++) {
            uint32_t const* pair = &order[p * 2 * kLookupsPerMate];
            for (uint32_t j = 0; j < 2 * kLookupsPerMate; j++) all_lines ? PrefetchAll(blocks[pair[j]]) : PrefetchFirst(blocks[pair[j]]);
            for (int m = 0; m < 2; m++) {
                mate.clear();
                for (uint32_t j = 0; j < kLookupsPerMate; j++) mate.push_back(blocks[pair[m * kLookupsPerMate + j]]);
                std::sort(mate.begin(), mate.end(), [](Lookup const& a, Lookup const& b) { return a.size < b.size; });
                seeds.clear();
                for (size_t i = 0; i < mate.size(); i++) {
                    if (i + kPrefetchLookups < mate.size()) all_lines ? PrefetchAll(mate[i + kPrefetchLookups]) : PrefetchFirst(mate[i + kPrefetchLookups]);
                    lookup(seeds, mate[i]);
                }
                total_seeds += seeds.size();
                for (auto const& s : seeds) {
                    checksum = checksum * 1000003 + (uint64_t{s.taxid} << 40 ^ uint64_t{s.geneid} << 20 ^ s.genepos ^ uint64_t{s.readpos} << 50 ^
                                                     uint64_t{s.unique} << 62 ^ uint64_t{s.unique_two} << 63);
                }
            }
        }
        double const ns = std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - start).count();
        double const lookups = static_cast<double>(pairs) * 2 * kLookupsPerMate;
        std::printf("%c: %.1f ns per lookup, %.2f seeds per lookup, checksum %016llx\n", variant, ns / lookups,
                    total_seeds / lookups, static_cast<unsigned long long>(checksum));
        std::fflush(stdout);
    };
    for (int r = 0; r < rounds; r++) {
        for (char v : variants) run(v);
    }
    return 0;
}
