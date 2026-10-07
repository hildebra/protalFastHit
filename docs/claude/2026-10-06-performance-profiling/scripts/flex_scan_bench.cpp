// flex_scan_bench.cpp - the best cells of a seed lookup's block of flex cells (Hash/FlexScan.h), as a lookup takes them:
// the scalar scan and a walk over its scores, against BestAvx2 with TiesAvx2 and a walk over the masks, for blocks of
// 16, 70 (the mean at GTDB r226) and 1,000 cells at every bit shift. Until 2026-10-06 this was the unit test
// FlexScan.Bench (run with PROTAL_FLEX_BENCH=1); the unit tests check that both give the same cells.
//
// Build: H=<a protal checkout with a build> cc8.sh flex_scan_bench.cpp flex_scan_bench; run: flex_scan_bench
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <random>
#include <vector>

#include "Hash/FlexScan.h"

using protal::Seedmap;
namespace flex_scan = protal::flex_scan;

namespace {
    // A block of `size` random flex cells at bit shift `shift`, followed by random bytes (as a block's entries follow
    // its cells) and the index's padding (Seedmap::kPackedPadding); every fourth cell (by chance) a copy of the key or
    // the key with a few bases changed, so that ties and exact matches occur.
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
            for (uint32_t i = 0; i < size; i++) {  // the cells at bit `shift` + 32 i, the other bits as they are
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
        }
    };

    // Room for BestAvx2's scores and TiesAvx2's reads of a block of `size` cells.
    size_t ScoreBytes(size_t size) { return std::max(size + flex_scan::kScorePadding, flex_scan::TieScoreBytes(size)); }
}

int main() {
    if (!flex_scan::CpuHasAvx2()) {
        std::cerr << "no AVX2 here" << std::endl;
        return 1;
    }
    std::mt19937 rng(3);
    for (uint32_t size : { 16u, 70u, 1000u }) {
        std::vector<RandomBlock> blocks;
        for (int i = 0; i < 256; i++) blocks.emplace_back(rng, size, static_cast<uint32_t>(i % 8));
        std::vector<uint8_t> scores(ScoreBytes(size));
        std::vector<uint32_t> masks(flex_scan::TieWords(size));
        uint64_t sink = 0;
        size_t const rounds = size < 100 ? 20000 : 2000;
        auto time = [&](bool avx2) {
            auto const start = std::chrono::steady_clock::now();
            for (size_t r = 0; r < rounds; r++) {
                for (auto& b : blocks) {
                    if (avx2) {
                        uint32_t const best = flex_scan::BestAvx2(b.block, b.key, scores.data());
                        flex_scan::TiesAvx2(scores.data(), size, best, masks.data());
                        for (uint32_t w = 0; w < flex_scan::TieWords(size); w++) {
                            for (uint32_t m = masks[w]; m; m &= m - 1) sink += 32 * w + static_cast<uint32_t>(__builtin_ctz(m));
                        }
                    } else {
                        auto const [best, count] = flex_scan::ScoreScalar(b.block, b.key, scores.data());
                        for (uint32_t i = 0; i < size; i++) sink += scores[i] == best ? i : 0;
                    }
                }
            }
            double const ns = std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - start).count();
            return ns / static_cast<double>(rounds * blocks.size() * size);
        };
        double const scalar = time(false), avx2 = time(true);
        std::cout << "blocks of " << size << " cells: scalar " << scalar << " ns per cell, AVX2 " << avx2 << " ns per cell ("
                  << scalar / avx2 << "x)" << (sink == 42 ? " " : "") << std::endl;
    }
    return 0;
}
