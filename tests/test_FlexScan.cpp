// Unit tests for the flex-cell scan of seed lookups (FlexScan.h): the AVX2 and the AVX-512 functions give every cell
// the scalar score, the same best score, and exactly the cells with it (as bit masks) and their count, for blocks of
// any size and bit shift and keys that match cells exactly, partly or not at all. Each runs where the CPU has it. (The bench of both is
// docs/claude/2026-10-06-performance-profiling/scripts/flex_scan_bench.cpp.)
#include <gtest/gtest.h>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <random>
#include <vector>
#include "Hash/FlexScan.h"

using namespace protal;

namespace {
    // A block of `size` random flex cells at bit shift `shift`, followed by random bytes (as a block's entries follow
    // its cells) and 64 bytes of padding (as the index's, Seedmap::kPackedPadding); every fourth cell (by chance) a
    // copy of the key or the key with a few bases changed, so that ties and exact matches occur.
    struct RandomBlock {
        std::vector<uint8_t> bytes;
        Seedmap::PackedBlock block;
        uint32_t key = 0;
        RandomBlock(std::mt19937& rng, uint32_t size, uint32_t shift) {
            key = static_cast<uint32_t>(rng());
            std::vector<uint32_t> cells(size);
            for (auto& c : cells) {
                uint32_t const kind = rng() % 4;
                if (kind == 0) c = key;
                else if (kind == 1) c = key ^ (3u << (2 * (rng() % 16))) ^ (1u << (2 * (rng() % 16)));
                else c = static_cast<uint32_t>(rng());
            }
            size_t const cell_bytes = (static_cast<size_t>(size) * 32 + shift + 7) / 8;
            bytes.assign(cell_bytes + 16 + Seedmap::kPackedPadding, 0);
            for (auto& b : bytes) b = static_cast<uint8_t>(rng());
            // Write the cells at bit `shift` + 32 i, keeping the other bits as they are.
            for (uint32_t i = 0; i < size; i++) {
                uint64_t const bit = shift + 32ull * i;
                uint64_t word;
                std::memcpy(&word, bytes.data() + bit / 8, 8);
                uint64_t const mask = 0xffffffffull << (bit % 8);
                word = (word & ~mask) | (static_cast<uint64_t>(cells[i]) << (bit % 8));
                std::memcpy(bytes.data() + bit / 8, &word, 8);
            }
            block.flex = bytes.data();
            block.flex_shift = shift;
            block.size = size;
            block.entries = bytes.data() + cell_bytes;
            for (uint32_t i = 0; i < size; i++) EXPECT_EQ(Seedmap::FlexCell(block, i), cells[i]) << "the test's own packing";
        }
    };

    // Room for BestAvx2's scores and TiesAvx2's reads of a block of `size` cells.
    size_t ScoreBytes(size_t size) { return std::max(size + flex_scan::kScorePadding, flex_scan::TieScoreBytes(size)); }
}

TEST(FlexScan, Avx2ScoresAndTiesAsTheScalarScan) {
    if (simd::Supported(simd::Level::avx2) != simd::Level::avx2) GTEST_SKIP() << "no AVX2 here";
    std::mt19937 rng(7);
    size_t blocks = 0, exact = 0;
    std::vector<uint32_t> sizes;
    for (uint32_t s = 1; s <= 80; s++) sizes.push_back(s);
    for (uint32_t s : { 95u, 96u, 97u, 127u, 128u, 129u, 255u, 256u, 257u, 1000u, 4095u, 4096u, 5000u }) sizes.push_back(s);
    for (uint32_t size : sizes) {
        for (uint32_t shift = 0; shift < 8; shift++) {
            for (int round = 0; round < 4; round++) {
                RandomBlock b(rng, size, shift);
                if (round == 3) b.key = static_cast<uint32_t>(rng());  // mostly without an exact match
                std::vector<uint8_t> scalar(size), avx2(ScoreBytes(size), 0xff);
                auto const [best, count] = flex_scan::ScoreScalar(b.block, b.key, scalar.data());
                uint32_t const avx2_best = flex_scan::BestAvx2(b.block, b.key, avx2.data());
                ASSERT_EQ(avx2_best, best) << "size " << size << ", shift " << shift;
                ASSERT_TRUE(std::equal(scalar.begin(), scalar.end(), avx2.begin())) << "size " << size << ", shift " << shift;
                // The bytes past the block may hold anything, the best score too: the masks stop at the block's end.
                std::fill(avx2.begin() + size, avx2.end(), static_cast<uint8_t>(best));
                std::vector<uint32_t> masks(flex_scan::TieWords(size) + 1, 0xdeadbeef);
                uint32_t const ties = flex_scan::TiesAvx2(avx2.data(), size, best, masks.data());
                ASSERT_EQ(ties, count) << "size " << size << ", shift " << shift;
                for (uint32_t i = 0; i < 32 * flex_scan::TieWords(size); i++) {
                    bool const bit = (masks[i / 32] >> (i % 32)) & 1;
                    ASSERT_EQ(bit, i < size && scalar[i] == best) << "cell " << i << ", size " << size << ", shift " << shift;
                }
                EXPECT_EQ(masks[flex_scan::TieWords(size)], 0xdeadbeefu) << "no mask word past the block's";
                blocks++;
                exact += best == 16;
            }
        }
    }
    std::cout << blocks << " blocks, " << exact << " with an exact match" << std::endl;
    EXPECT_GT(exact, blocks / 2);
    EXPECT_LT(exact, blocks);
}

TEST(FlexScan, ABlockWithoutASharedBaseHasEveryCellBest) {
    // Key 0 (AAAA...) against cells of only T (11 pairs) shares nothing: best 0, every cell with it, in both versions.
    std::vector<uint8_t> bytes(4 * 40 + flex_scan::kReadPastCells, 0xff);
    Seedmap::PackedBlock block;
    block.flex = bytes.data();
    block.flex_shift = 0;
    block.size = 40;
    std::vector<uint8_t> scores(ScoreBytes(40));
    EXPECT_EQ(flex_scan::ScoreScalar(block, 0, scores.data()), (std::pair<uint32_t, uint32_t>{ 0, 40 }));
    if (simd::Supported(simd::Level::avx2) == simd::Level::avx2) {
        EXPECT_EQ(flex_scan::BestAvx2(block, 0, scores.data()), 0u);
        std::vector<uint32_t> masks(flex_scan::TieWords(40));
        EXPECT_EQ(flex_scan::TiesAvx2(scores.data(), 40, 0, masks.data()), 40u);
        EXPECT_EQ(masks[0], 0xffffffffu);
        EXPECT_EQ(masks[1], 0xffu);
    }
    if (simd::Supported(simd::Level::avx512) == simd::Level::avx512) {
        std::vector<uint32_t> masks(flex_scan::TieWords(40));
        EXPECT_EQ(flex_scan::ScanAvx512(block, 0, scores.data(), masks.data()), (std::pair<uint32_t, uint32_t>{ 0, 40 }));
        EXPECT_EQ(masks[0], 0xffffffffu);
        EXPECT_EQ(masks[1], 0xffu);
    }
}

// ScanAvx512 gives every cell the scalar score and writes no score past the block, the same best score, exactly the
// cells with it as bit masks (none past the block, no mask word past the block's) and their count; for blocks of any
// size and bit shift whose cells end where readable memory ends (the masked loads read nothing past the cells and the
// 4 bytes after them that FlexCell reads too).
TEST(FlexScan, Avx512ScoresAndTiesAsTheScalarScan) {
    if (simd::Supported(simd::Level::avx512) != simd::Level::avx512) GTEST_SKIP() << "no AVX-512 (Ice Lake / Zen 4 level) here";
    std::mt19937 rng(11);
    size_t blocks = 0, exact = 0;
    std::vector<uint32_t> sizes;
    for (uint32_t s = 1; s <= 140; s++) sizes.push_back(s);
    for (uint32_t s : { 191u, 192u, 193u, 255u, 256u, 257u, 1000u, 4095u, 4096u, 5000u }) sizes.push_back(s);
    for (uint32_t size : sizes) {
        for (uint32_t shift = 0; shift < 8; shift++) {
            for (int round = 0; round < 4; round++) {
                RandomBlock b(rng, size, shift);
                if (round == 3) b.key = static_cast<uint32_t>(rng());  // mostly without an exact match
                // The cells copied to the end of a buffer whose last byte is the one FlexCell reads last, so that the
                // address checker (or a page end) would notice a read past them.
                size_t const read_bytes = 4 * static_cast<size_t>(size) + 4 + (shift > 0);
                std::vector<uint8_t> tight(b.bytes.begin(), b.bytes.begin() + static_cast<std::ptrdiff_t>(read_bytes));
                Seedmap::PackedBlock block = b.block;
                block.flex = tight.data();
                std::vector<uint8_t> scalar(size);
                auto const [best, count] = flex_scan::ScoreScalar(block, b.key, scalar.data());
                std::vector<uint8_t> scores(size + 64, 0xee);
                std::vector<uint32_t> masks(flex_scan::TieWords(size) + 1, 0xdeadbeef);
                auto const [avx512_best, avx512_count] = flex_scan::ScanAvx512(block, b.key, scores.data(), masks.data());
                ASSERT_EQ(avx512_best, best) << "size " << size << ", shift " << shift;
                ASSERT_EQ(avx512_count, count) << "size " << size << ", shift " << shift;
                ASSERT_TRUE(std::equal(scalar.begin(), scalar.end(), scores.begin())) << "size " << size << ", shift " << shift;
                ASSERT_TRUE(std::all_of(scores.begin() + size, scores.end(), [](uint8_t x) { return x == 0xee; }))
                    << "no score written past the block, size " << size;
                for (uint32_t i = 0; i < 32 * flex_scan::TieWords(size); i++) {
                    bool const bit = (masks[i / 32] >> (i % 32)) & 1;
                    ASSERT_EQ(bit, i < size && scalar[i] == best) << "cell " << i << ", size " << size << ", shift " << shift;
                }
                EXPECT_EQ(masks[flex_scan::TieWords(size)], 0xdeadbeefu) << "no mask word past the block's";
                blocks++;
                exact += best == 16;
            }
        }
    }
    std::cout << blocks << " blocks, " << exact << " with an exact match" << std::endl;
    EXPECT_GT(exact, blocks / 2);
    EXPECT_LT(exact, blocks);
}

// The level lookups use is the CPU's (capped by PROTAL_SIMD), and Use() lowers a wanted level to what the CPU has.
TEST(FlexScan, LevelsFollowTheCpu) {
    simd::Level const was = flex_scan::Kernel().load();
    EXPECT_LE(static_cast<int>(simd::Default()), static_cast<int>(simd::CpuLevel()));
    EXPECT_EQ(flex_scan::Use(simd::Level::avx512), simd::CpuLevel());
    EXPECT_EQ(flex_scan::Use(simd::Level::scalar), simd::Level::scalar);
    EXPECT_EQ(flex_scan::Kernel().load(), simd::Level::scalar);
    flex_scan::Kernel().store(was);
    std::cout << "this CPU: " << simd::Name(simd::CpuLevel()) << std::endl;
}
