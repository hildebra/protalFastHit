// Unit tests for the flex-cell scan of seed lookups (FlexScan.h): the AVX2 version gives every cell the scalar score,
// and the same best score and count, for blocks of any size and bit shift and keys that match cells exactly, partly or
// not at all; and a bench of both (PROTAL_FLEX_BENCH=1).
#include <gtest/gtest.h>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <random>
#include <vector>
#include "Hash/FlexScan.h"

using namespace protal;

namespace {
    // A block of `size` random flex cells at bit shift `shift`, followed by random bytes (as a block's entries follow
    // its cells) and 16 bytes of padding (as the index's); every fourth cell (by chance) a copy of the key or the key
    // with a few bases changed, so that ties and exact matches occur.
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
            bytes.assign(cell_bytes + 64 + 16, 0);
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
}

TEST(FlexScan, Avx2ScoresAsTheScalarScan) {
    if (!flex_scan::CpuHasAvx2()) GTEST_SKIP() << "no AVX2 here";
    std::mt19937 rng(7);
    size_t blocks = 0, exact = 0;
    std::vector<uint32_t> sizes;
    for (uint32_t s = 1; s <= 80; s++) sizes.push_back(s);
    for (uint32_t s : { 127u, 128u, 129u, 255u, 256u, 257u, 1000u, 4095u, 4096u, 5000u }) sizes.push_back(s);
    for (uint32_t size : sizes) {
        for (uint32_t shift = 0; shift < 8; shift++) {
            for (int round = 0; round < 4; round++) {
                RandomBlock b(rng, size, shift);
                if (round == 3) b.key = static_cast<uint32_t>(rng());  // mostly without an exact match
                std::vector<uint8_t> scalar(size), avx2(size, 0xff);
                auto const s = flex_scan::ScoreScalar(b.block, b.key, scalar.data());
                auto const a = flex_scan::ScoreAvx2(b.block, b.key, avx2.data());
                ASSERT_EQ(a, s) << "size " << size << ", shift " << shift;
                ASSERT_EQ(avx2, scalar) << "size " << size << ", shift " << shift;
                blocks++;
                exact += s.first == 16;
            }
        }
    }
    std::cout << blocks << " blocks, " << exact << " with an exact match" << std::endl;
    EXPECT_GT(exact, blocks / 2);
    EXPECT_LT(exact, blocks);
}

TEST(FlexScan, ABlockWithoutASharedBaseHasEveryCellBest) {
    // Key 0 (AAAA...) against cells of only T (11 pairs) shares nothing: best 0, every cell with it, in both versions.
    std::vector<uint8_t> bytes(4 * 40 + 32, 0xff);
    Seedmap::PackedBlock block;
    block.flex = bytes.data();
    block.flex_shift = 0;
    block.size = 40;
    std::vector<uint8_t> scores(40);
    EXPECT_EQ(flex_scan::ScoreScalar(block, 0, scores.data()), (std::pair<uint32_t, uint32_t>{ 0, 40 }));
    if (flex_scan::CpuHasAvx2()) EXPECT_EQ(flex_scan::ScoreAvx2(block, 0, scores.data()), (std::pair<uint32_t, uint32_t>{ 0, 40 }));
}

// A bench, not a test: PROTAL_FLEX_BENCH=1 scores blocks of 70 cells (the mean at GTDB r226) and of 1,000 with both versions.
TEST(FlexScan, Bench) {
    if (!std::getenv("PROTAL_FLEX_BENCH")) GTEST_SKIP() << "set PROTAL_FLEX_BENCH=1 to run";
    if (!flex_scan::CpuHasAvx2()) GTEST_SKIP() << "no AVX2 here";
    std::mt19937 rng(3);
    for (uint32_t size : { 16u, 70u, 1000u }) {
        std::vector<RandomBlock> blocks;
        for (int i = 0; i < 256; i++) blocks.emplace_back(rng, size, static_cast<uint32_t>(i % 8));
        std::vector<uint8_t> scores(size);
        uint64_t sink = 0;
        size_t const rounds = size < 100 ? 20000 : 2000;
        auto time = [&](bool avx2) {
            auto const start = std::chrono::steady_clock::now();
            for (size_t r = 0; r < rounds; r++) {
                for (auto& b : blocks) {
                    auto const s = avx2 ? flex_scan::ScoreAvx2(b.block, b.key, scores.data()) : flex_scan::ScoreScalar(b.block, b.key, scores.data());
                    sink += s.first + s.second;
                }
            }
            double const ns = std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - start).count();
            return ns / static_cast<double>(rounds * blocks.size() * size);
        };
        double const scalar = time(false), avx2 = time(true);
        std::cout << "blocks of " << size << " cells: scalar " << scalar << " ns per cell, AVX2 " << avx2 << " ns per cell ("
                  << scalar / avx2 << "x)" << (sink == 42 ? "" : "") << std::endl;
    }
}
