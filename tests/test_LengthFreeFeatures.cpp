// What makes the presence features independent of the read length and of the share of a read's cores looked up (since
// 2026-10-14, docs/claude/2026-10-14-sparse-index-calls): the seeds of a record's anchor (ZS), the failed candidates without
// the single-seed anchors' (ZF), the identity bound of a read by its aligned bases, and the sample's seed expectation in its
// SAM's header.
#include <gtest/gtest.h>
#include <sstream>
#include <string>
#include <vector>
#include "Alignment/AlignmentUtils.h"
#include "IO/SamHandler.h"
#include "IO/SeedExpectation.h"
#include "Profiling/Strain.h"

using namespace protal;

// A taxon tried only on anchors of a single grid seed (lone) is no failed candidate; one tried on a two-seed anchor too is,
// with the gene of its first two-seed anchor.
TEST(LengthFreeFeatures, FailedCandidatesLeaveOutTheSingleSeedAnchors) {
    using C = FailedCandidate;
    std::vector<C> const attempted{ C(5, 3, false, true), C(1, 2), C(5, 8), C(7, 4, false, true), C(2, 9) };
    EXPECT_EQ(FailedCandidates(attempted, { 2 }), (std::vector<C>{ C(1, 2), C(5, 8) }));
    // Aligned to: none either way (taxon 2, tried on a two-seed anchor, still failed).
    EXPECT_EQ(FailedCandidates(attempted, { 1, 5 }), (std::vector<C>{ C(2, 9) }));
    // Without lone attempts as before.
    EXPECT_EQ(FailedCandidates({ C(3, 7), C(1, 2), C(3, 9) }, {}), (std::vector<C>{ C(1, 2), C(3, 7) }));
}

// The lower end of the 95% Wilson interval: a short read without a difference says less surely that it is identical.
TEST(LengthFreeFeatures, IdentityLowerBoundByAlignedBases) {
    EXPECT_NEAR(IdentityLowerBound(1.0, 150), 0.9750, 1e-4);
    EXPECT_NEAR(IdentityLowerBound(1.0, 75), 0.9513, 1e-4);
    EXPECT_NEAR(IdentityLowerBound(149.0 / 150.0, 150), 0.9632, 1e-4);
    EXPECT_GT(IdentityLowerBound(149.0 / 150.0, 150), IdentityLowerBound(1.0, 75));
    EXPECT_NEAR(IdentityLowerBound(0.999, 15000), 0.9984, 1e-4);  // a long read: near its identity
    EXPECT_EQ(IdentityLowerBound(1.0, 0), 0.0);
    EXPECT_LT(IdentityLowerBound(0.9, 100), 0.9);
}

TEST(LengthFreeFeatures, SeedExpectationLineRoundTrips) {
    SeedExpectation e;
    e.level = 2;
    e.unique = 0.7891234567;
    e.unique_two = 1.0625;
    std::string const line = SeedExpectationLine(e);
    ASSERT_EQ(line.back(), '\n');
    auto const parsed = ParseSeedExpectationLine(line.substr(0, line.size() - 1));
    ASSERT_TRUE(parsed.has_value());
    EXPECT_EQ(parsed->level, 2u);
    EXPECT_EQ(parsed->unique, e.unique);
    EXPECT_EQ(parsed->unique_two, e.unique_two);
    EXPECT_TRUE(ParseSeedExpectationLine(kSamSeedExpectationComment + "two=1 later=x level=1 unique=0.5").has_value());
    for (std::string bad : { "level=1 unique=0.5", "level=1 unique=0 two=1", "level=x unique=1 two=1", "level=1 unique=1 two=" }) {
        EXPECT_FALSE(ParseSeedExpectationLine(kSamSeedExpectationComment + bad).has_value()) << bad;
    }
}

// ZS is written where an anchor's seeds were counted, read back, and absent (0) otherwise; the SAM reader finds the
// header's seed expectation.
TEST(LengthFreeFeatures, SeedsTagAndHeaderLineAreRead) {
    SamEntry sam;
    sam.m_qname = "r1";
    sam.m_flag = 0;
    sam.m_rname = "12_3";
    sam.m_pos = 5;
    sam.m_mapq = 60;
    sam.m_cigar = "4M";
    sam.m_rnext = "*";
    sam.m_pnext = 0;
    sam.m_tlen = 0;
    sam.m_seq = "ACGT";
    sam.m_qual = "IIII";
    sam.m_uniques = 3;
    sam.m_uniques_two = 2;
    sam.m_seeds = 7;
    std::string const with = sam.ToString();
    EXPECT_NE(with.find("\tZU:i:3\tZT:i:2\tZS:i:7"), std::string::npos) << with;
    sam.m_qname = "r2";
    sam.m_seeds = 0;
    std::string const without = sam.ToString();
    EXPECT_EQ(without.find("ZS:i:"), std::string::npos) << without;

    SeedExpectation e;
    e.level = 3;
    e.unique = 0.8;
    e.unique_two = 1.1;
    std::string const text = "@HD\tVN:1.6\n" + SeedExpectationLine(e) + with + "\n" + without + "\n";
    SamReader reader(text, {}, 0);
    SamEntry a, b;
    bool has_a = false, has_b = false;
    ASSERT_TRUE(reader.Next(a, b, has_a, has_b));
    ASSERT_TRUE(has_a);
    EXPECT_EQ(a.m_seeds, 7);
    EXPECT_EQ(a.m_uniques, 3);
    ASSERT_TRUE(reader.Next(a, b, has_a, has_b));
    EXPECT_EQ(a.m_seeds, 0);
    ASSERT_TRUE(reader.SeedExpectationOf().has_value());
    EXPECT_EQ(reader.SeedExpectationOf()->level, 3u);
    EXPECT_EQ(reader.SeedExpectationOf()->unique, 0.8);

    SamReader bad("@CO\tprotal seed expectation: level=1\n" + with + "\n", {}, 0);
    EXPECT_THROW(bad.Next(a, b, has_a, has_b), SamFormatError);
}
