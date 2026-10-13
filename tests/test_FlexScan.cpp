// Unit tests for the flex-cell scan of seed lookups (FlexScan.h): the AVX2 and the AVX-512 functions give every cell
// the scalar score, the same best score, and exactly the cells with it (as bit masks) and their count, for blocks of
// any size and bit shift and keys that match cells exactly, partly or not at all. Each runs where the CPU has it. (The bench of both is
// docs/claude/2026-10-06-performance-profiling/scripts/flex_scan_bench.cpp.) The same for the scores of --flex_indels
// (IndelSimilarity), which are first checked against their definition base by base and on reads with one error.
#include <gtest/gtest.h>
#include <algorithm>
#include <array>
#include <cstring>
#include <iostream>
#include <random>
#include <vector>
#include "Hash/FlexScan.h"

using namespace protal;

namespace {
    // A block of `size` random flex cells at bit shift `shift`, followed by random bytes (as a block's entries follow
    // its cells) and 64 bytes of padding (as the index's, Seedmap::kPackedPadding); every fourth cell (by chance) a
    // copy of the key, the key with a few bases changed, or (every eighth) the key with a flank shifted by one base
    // and a base changed, so that ties, exact matches and indels occur.
    struct RandomBlock {
        std::vector<uint8_t> bytes;
        Seedmap::PackedBlock block;
        uint32_t key = 0;
        RandomBlock(std::mt19937& rng, uint32_t size, uint32_t shift) {
            key = static_cast<uint32_t>(rng());
            std::vector<uint32_t> cells(size);
            for (auto& c : cells) {
                uint32_t const kind = rng() % 8;
                if (kind < 2) c = key;
                else if (kind < 4) c = key ^ (3u << (2 * (rng() % 16))) ^ (1u << (2 * (rng() % 16)));
                else if (kind == 4) {
                    flex_scan::IndelKeys const shifted = flex_scan::MakeIndelKeys(key);
                    uint32_t const flank = rng() % 2 ? 0xFFFF0000u : 0x0000FFFFu;
                    c = ((rng() % 2 ? shifted.deletion : shifted.insertion) & flank) | (key & ~flank);
                    if (rng() % 2) c ^= (1u + rng() % 3) << (2 * (rng() % 16));
                }
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

    // The kernels with --flex_indels' scores or without, chosen at run time.
    std::pair<uint32_t, uint32_t> Scalar(bool indels, Seedmap::PackedBlock const& block, uint32_t key, uint8_t* scores) {
        return indels ? flex_scan::ScoreScalar<true>(block, key, scores) : flex_scan::ScoreScalar<false>(block, key, scores);
    }
    uint32_t Avx2(bool indels, Seedmap::PackedBlock const& block, uint32_t key, uint8_t* scores) {
        return indels ? flex_scan::BestAvx2<true>(block, key, scores) : flex_scan::BestAvx2<false>(block, key, scores);
    }
    std::pair<uint32_t, uint32_t> Avx512(bool indels, Seedmap::PackedBlock const& block, uint32_t key, uint8_t* scores, uint32_t* masks) {
        return indels ? flex_scan::ScanAvx512<true>(block, key, scores, masks) : flex_scan::ScanAvx512<false>(block, key, scores, masks);
    }

    // IndelSimilarity as its comment defines it, base by base: base i (1-8, out from the core) of a flex key's left
    // (side 0) or right flank (side 1); per flank the most of the bases equal in place, and from the first that differs
    // on the read's bases one place in or out.
    int Base(uint32_t x, int side, int i) {
        int const pair = side == 0 ? 7 + i : 8 - i;
        return static_cast<int>((x >> (2 * pair)) & 3u);
    }
    uint32_t DefinedFlank(uint32_t cell, uint32_t key, int side) {
        int equal = 0, first = 0;
        for (int i = 1; i <= 8; i++) {
            if (Base(cell, side, i) == Base(key, side, i)) equal++;
            else if (first == 0) first = i;
        }
        if (first == 0) return 8;
        int deletion = first - 1, insertion = first - 1;
        for (int i = first + 1; i <= 8; i++) deletion += Base(cell, side, i) == Base(key, side, i - 1);  // the read lacks base `first`
        for (int i = first; i <= 7; i++) insertion += Base(cell, side, i) == Base(key, side, i + 1);     // the read has one extra
        return static_cast<uint32_t>(std::max({ equal, deletion, insertion }));
    }
    uint32_t DefinedIndelSimilarity(uint32_t cell, uint32_t key) {
        return DefinedFlank(cell, key, 0) + DefinedFlank(cell, key, 1);
    }

    // The flex key of a k-mer from its flanks' bases (0-3) in sequence order: the left flank's 8 (its first base the
    // k-mer's first), the right flank's 8 (its first next to the core).
    uint32_t FlexKeyOf(std::vector<int> const& left, std::vector<int> const& right) {
        uint32_t key = 0;
        for (int b : left) key = (key << 2) | static_cast<uint32_t>(b);
        for (int b : right) key = (key << 2) | static_cast<uint32_t>(b);
        return key;
    }
}

TEST(FlexScan, IndelSimilarityIsItsDefinition) {
    std::mt19937 rng(3);
    size_t cells = 0, shifted = 0;
    for (int round = 0; round < 2000000; round++) {
        uint32_t const key = static_cast<uint32_t>(rng());
        uint32_t cell = static_cast<uint32_t>(rng());
        int const kind = round % 4;
        if (kind == 1) cell = key ^ ((1u + rng() % 3) << (2 * (rng() % 16)));  // one base changed
        if (kind >= 2) {  // a flank of the key shifted one base in or out, maybe a base changed too
            flex_scan::IndelKeys const keys = flex_scan::MakeIndelKeys(key);
            uint32_t const flank = rng() % 2 ? 0xFFFF0000u : 0x0000FFFFu;
            cell = ((kind == 2 ? keys.deletion : keys.insertion) & flank) | (key & ~flank);
            if (rng() % 2) cell ^= (1u + rng() % 3) << (2 * (rng() % 16));
            shifted++;
        }
        uint32_t const score = flex_scan::IndelSimilarity(cell, flex_scan::MakeIndelKeys(key));
        ASSERT_EQ(score, DefinedIndelSimilarity(cell, key)) << std::hex << "cell " << cell << ", key " << key;
        ASSERT_GE(score, Seedmap::Similarity(cell, key)) << "never below the score without indels";
        ASSERT_EQ(score == 16, cell == key) << "16 only for an exact match";
        cells++;
    }
    EXPECT_GT(shifted, cells / 3);
}

// A read with one error in a flank of its k-mer, against the gene's flex cell: a substitution, an insertion or a deletion
// at each place of either flank scores 15 (16 where the read's flank is the gene's all the same, as an insertion of
// the base next to it can leave it), while the score without indels falls by up to 6 for an indel.
TEST(FlexScan, OneErrorInAFlankCostsOneBase) {
    std::mt19937 rng(5);
    std::array<std::array<double, 9>, 3> lost_without{};  // by error and place, the bases lost without indels
    std::array<std::array<int, 9>, 3> counted{};
    for (int round = 0; round < 20000; round++) {
        // The gene: its flanks with two more bases beyond each (a deletion takes the read's flank one base further).
        std::vector<int> outer_left(10), outer_right(10);  // out from the core: [0] next to it
        for (auto& b : outer_left) b = static_cast<int>(rng() % 4);
        for (auto& b : outer_right) b = static_cast<int>(rng() % 4);
        if (round % 5 == 0) for (int i = 1; i < 4; i++) outer_left[i] = outer_left[0];  // homopolymers next to the core
        auto key_of = [](std::vector<int> const& left_out, std::vector<int> const& right_out) {
            std::vector<int> left(left_out.rbegin() + static_cast<std::ptrdiff_t>(left_out.size() - 8), left_out.rend());
            std::vector<int> right(right_out.begin(), right_out.begin() + 8);
            return FlexKeyOf(left, right);
        };
        uint32_t const cell = key_of(outer_left, outer_right);
        int const error = round % 3;  // 0 substitution, 1 insertion, 2 deletion in the read
        int const place = 1 + static_cast<int>(rng() % 8);
        bool const on_left = rng() % 2;
        std::vector<int> read = on_left ? outer_left : outer_right;
        if (error == 0) read[place - 1] = (read[place - 1] + 1 + static_cast<int>(rng() % 3)) % 4;
        if (error == 1) read.insert(read.begin() + (place - 1), static_cast<int>(rng() % 4));
        if (error == 2) read.erase(read.begin() + (place - 1));
        uint32_t const key = on_left ? key_of(read, outer_right) : key_of(outer_left, read);
        uint32_t const score = flex_scan::IndelSimilarity(cell, flex_scan::MakeIndelKeys(key));
        ASSERT_EQ(score, key == cell ? 16u : 15u) << "error " << error << " at place " << place << (on_left ? " left" : " right");
        lost_without[error][place] += 16.0 - Seedmap::Similarity(cell, key);
        counted[error][place]++;
    }
    for (int error = 0; error < 3; error++) {
        std::cout << (error == 0 ? "substitution" : error == 1 ? "insertion   " : "deletion    ") << " at places 1-8, bases lost without indels:";
        for (int place = 1; place <= 8; place++) std::cout << ' ' << lost_without[error][place] / std::max(1, counted[error][place]);
        std::cout << std::endl;
    }
    EXPECT_GT(lost_without[2][1] / counted[2][1], 4.0) << "a deletion next to the core loses most of the flank without indels";
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
            for (int round = 0; round < 8; round++) {
                bool const indels = round >= 4;
                RandomBlock b(rng, size, shift);
                if (round % 4 == 3) b.key = static_cast<uint32_t>(rng());  // mostly without an exact match
                std::vector<uint8_t> scalar(size), avx2(ScoreBytes(size), 0xff);
                auto const [best, count] = Scalar(indels, b.block, b.key, scalar.data());
                uint32_t const avx2_best = Avx2(indels, b.block, b.key, avx2.data());
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
    for (bool indels : { false, true }) {  // shifted, key 0 is still all A
        EXPECT_EQ(Scalar(indels, block, 0, scores.data()), (std::pair<uint32_t, uint32_t>{ 0, 40 }));
        if (simd::Supported(simd::Level::avx2) == simd::Level::avx2) {
            EXPECT_EQ(Avx2(indels, block, 0, scores.data()), 0u);
            std::vector<uint32_t> masks(flex_scan::TieWords(40));
            EXPECT_EQ(flex_scan::TiesAvx2(scores.data(), 40, 0, masks.data()), 40u);
            EXPECT_EQ(masks[0], 0xffffffffu);
            EXPECT_EQ(masks[1], 0xffu);
        }
        if (simd::Supported(simd::Level::avx512) == simd::Level::avx512) {
            std::vector<uint32_t> masks(flex_scan::TieWords(40));
            EXPECT_EQ(Avx512(indels, block, 0, scores.data(), masks.data()), (std::pair<uint32_t, uint32_t>{ 0, 40 }));
            EXPECT_EQ(masks[0], 0xffffffffu);
            EXPECT_EQ(masks[1], 0xffu);
        }
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
            for (int round = 0; round < 8; round++) {
                bool const indels = round >= 4;
                RandomBlock b(rng, size, shift);
                if (round % 4 == 3) b.key = static_cast<uint32_t>(rng());  // mostly without an exact match
                // The cells copied to the end of a buffer whose last byte is the one FlexCell reads last, so that the
                // address checker (or a page end) would notice a read past them.
                size_t const read_bytes = 4 * static_cast<size_t>(size) + 4 + (shift > 0);
                std::vector<uint8_t> tight(b.bytes.begin(), b.bytes.begin() + static_cast<std::ptrdiff_t>(read_bytes));
                Seedmap::PackedBlock block = b.block;
                block.flex = tight.data();
                std::vector<uint8_t> scalar(size);
                auto const [best, count] = Scalar(indels, block, b.key, scalar.data());
                std::vector<uint8_t> scores(size + 64, 0xee);
                std::vector<uint32_t> masks(flex_scan::TieWords(size) + 1, 0xdeadbeef);
                auto const [avx512_best, avx512_count] = Avx512(indels, block, b.key, scores.data(), masks.data());
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
