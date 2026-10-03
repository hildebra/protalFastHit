// SequenceRange::CoverageVector (a difference array) gives the counts of adding every read base by base.
#include <gtest/gtest.h>
#include <random>
#include <vector>
#include "SequenceRange.h"
#include "SequenceRangeHandler.h"

namespace {
    CoverageVec Reference(SequenceRange const& range, std::vector<ReadInfo> const& reads, int strand, uint8_t max_divergence) {
        CoverageVec coverage(range.Length(), 0);
        for (auto const& read : reads) {
            if (strand != SequenceRange::kBothStrands && read.forward != (strand == SequenceRange::kForward)) continue;
            if (read.divergence > max_divergence) continue;
            for (uint32_t i = 0; i < read.length; i++) coverage[i + (read.start - range.Start())]++;
        }
        return coverage;
    }
}

TEST(SequenceRange, CoverageVectorEqualsAddingEveryBase) {
    std::mt19937 rng(7);
    for (int round = 0; round < 200; round++) {
        size_t const start = rng() % 500, length = 1 + rng() % 3000;
        SequenceRange range(start, start + length);
        std::vector<ReadInfo> reads;
        int const n = static_cast<int>(rng() % 120);
        for (int i = 0; i < n; i++) {
            uint32_t const len = 1 + rng() % std::min<size_t>(length, 150);
            uint32_t const from = static_cast<uint32_t>(start + rng() % (length - len + 1));  // reads end at the range's end at most
            ReadInfo read{static_cast<size_t>(i), from, len, (rng() & 1) != 0, static_cast<uint8_t>(rng() % 6)};
            reads.push_back(read);
            range.AddReadInfo(read);
        }
        for (int strand : {SequenceRange::kBothStrands, SequenceRange::kForward, SequenceRange::kReverse})
            for (uint8_t cap : {uint8_t{255}, uint8_t{2}, uint8_t{0}})
                ASSERT_EQ(range.CoverageVector(strand, cap), Reference(range, reads, strand, cap));
    }
}

TEST(SequenceRange, CoveredPortionCountsBasesWithReads) {
    SequenceRangeHandler handler;
    SequenceRange range(10, 60);
    range.AddReadInfo(ReadInfo{0, 10, 20, true});
    range.AddReadInfo(ReadInfo{1, 25, 10, false});   // 25..35 overlaps nothing before 30; covers 10..30 and 25..35
    handler.Add(range);
    EXPECT_EQ(handler.CoveredPortion(), 25u);
    EXPECT_EQ(handler.CoveredPortion(2), 5u);
}
