// A sample's spike genes (SpikeGenes.h, since 2026-10-15): genes whose reads pile far beyond what the rest of their
// taxon's genes explain, found per taxon by SpikesOf, and the records per gene a SAM's text gives (CountRecords).
#include <gtest/gtest.h>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include "Profiling/SpikeGenes.h"

namespace sg = protal::spike_genes;

namespace {
    using Records = std::vector<std::pair<uint32_t, uint64_t>>;

    // `genes` hittable genes 1..genes of length `length` each.
    Records Lengths(uint32_t genes, uint64_t length = 1000) {
        Records lengths;
        for (uint32_t g = 1; g <= genes; g++) lengths.emplace_back(g, length);
        return lengths;
    }

    // Genes 1..hit with `each` records each.
    Records Even(uint32_t hit, uint64_t each) {
        Records records;
        for (uint32_t g = 1; g <= hit; g++) records.emplace_back(g, each);
        return records;
    }
}

// The D6331 sample's Akkermansia: 2,108 records on one gene, 1-2 on each of 37 others of 113.
TEST(SpikeGenes, OneGeneHoldingAlmostAllRecordsIsASpike) {
    auto records = Even(37, 2);
    records.emplace_back(71, 2108);
    auto const spikes = sg::SpikesOf(758, records, Lengths(113));
    ASSERT_EQ(spikes.size(), 1u);
    EXPECT_EQ(spikes[0].taxid, 758u);
    EXPECT_EQ(spikes[0].geneid, 71u);
    EXPECT_EQ(spikes[0].records, 2108u);
    // The background leaves out the three genes with the most records per base (gene 71 and two of the others): 70
    // records over 110 genes of equal length.
    EXPECT_NEAR(spikes[0].expected, 70.0 / 110.0, 1e-9);
    EXPECT_GT(spikes[0].Fold(), 1000);
}

// Bacteroides fragilis there: two spikes of about equal size (48% and 47% of its records), neither half on its own.
TEST(SpikeGenes, TwoSpikesTogetherHoldingMostRecordsAreBothTaken) {
    auto records = Even(90, 7);
    records.emplace_back(101, 7761);
    records.emplace_back(105, 7915);
    auto const spikes = sg::SpikesOf(522, records, Lengths(120));
    ASSERT_EQ(spikes.size(), 2u);
    EXPECT_EQ(spikes[0].geneid, 105u);  // the most records per base first
    EXPECT_EQ(spikes[1].geneid, 101u);
}

// A gene well above the rest, but far below half of the taxon's records: not a spike (a relative's reads on a conserved
// gene, or a simulated present species' noisy gene).
TEST(SpikeGenes, AnOutlierHoldingFewRecordsIsNoSpike) {
    auto records = Even(60, 10);
    records.emplace_back(100, 400);  // about 80 times the others' depth, but 40% of the records
    EXPECT_TRUE(sg::SpikesOf(1, records, Lengths(120)).empty());
    auto few = Even(6, 1);
    few.emplace_back(100, sg::kMinRecords - 1);  // most of the records, far deeper, but under kMinRecords
    EXPECT_TRUE(sg::SpikesOf(1, few, Lengths(120)).empty());
    few.back().second = sg::kMinRecords;
    EXPECT_EQ(sg::SpikesOf(1, few, Lengths(120)).size(), 1u);
}

// A taxon with reads on too few other genes keeps its pile: an absent relative whose reads gather on one gene stays as
// it was.
TEST(SpikeGenes, TooFewOtherGenesLeaveThePile) {
    auto records = Even(sg::kMinRestGenes, 1);
    records.emplace_back(50, 5000);
    EXPECT_TRUE(sg::SpikesOf(1, records, Lengths(120)).empty());
}

// Even coverage, deep or shallow, and Poisson-like noise give no spike.
TEST(SpikeGenes, EvenCoverageHasNoSpikes) {
    EXPECT_TRUE(sg::SpikesOf(1, Even(120, 2000), Lengths(120)).empty());
    EXPECT_TRUE(sg::SpikesOf(1, Even(40, 1), Lengths(120)).empty());
    auto records = Even(120, 1000);
    records[5].second = 1400;
    EXPECT_TRUE(sg::SpikesOf(1, records, Lengths(120)).empty());
}

// Genes of other lengths: the expectation is per base, so a long gene with more records is no spike.
TEST(SpikeGenes, ExpectationsArePerBase) {
    Records lengths = Lengths(100, 500);
    lengths.emplace_back(200, 50000);  // a gene 100 times as long
    auto records = Even(100, 5);
    records.emplace_back(200, 500);    // 100 times the records: the same depth
    EXPECT_TRUE(sg::SpikesOf(1, records, lengths).empty());
}

// At most kMaxGenes spikes per taxon; records on genes without a known length are left out.
TEST(SpikeGenes, AtMostMaxGenesAndOnlyKnownGenes) {
    auto records = Even(30, 1);
    for (uint32_t g = 101; g <= 105; g++) records.emplace_back(g, 1000);
    auto const spikes = sg::SpikesOf(1, records, Lengths(120));
    EXPECT_EQ(spikes.size(), sg::kMaxGenes);
    records.emplace_back(999, 100000);  // not a hittable gene of the taxon
    auto const again = sg::SpikesOf(1, records, Lengths(120));
    ASSERT_EQ(again.size(), sg::kMaxGenes);
    for (auto const& spike : again) EXPECT_NE(spike.geneid, 999u);
}

namespace {
    // A sample's counts: taxon `taxid` with `each` records on genes 1..genes, and `pile` more on gene 1 if given.
    void AddTaxon(sg::Counts& counts, uint32_t taxid, uint32_t genes, uint64_t each, uint64_t pile = 0) {
        for (uint32_t g = 1; g <= genes; g++) counts[sg::Key(taxid, g)] += each;
        if (pile) counts[sg::Key(taxid, 1)] += pile;
    }

    std::vector<sg::Spike> SpikesOfAll(sg::Counts const& counts, uint32_t genes) {
        std::map<uint32_t, Records> by_taxon;
        for (auto const& [key, n] : counts) by_taxon[static_cast<uint32_t>(key >> 32)].emplace_back(static_cast<uint32_t>(key), n);
        std::vector<sg::Spike> spikes;
        for (auto& [taxid, records] : by_taxon) {
            auto found = sg::SpikesOf(taxid, records, Lengths(genes));
            spikes.insert(spikes.end(), found.begin(), found.end());
        }
        return spikes;
    }
}

// The Dohlman sample's pattern: the abundant species all have a pile -> spike mode.
TEST(SpikeGenes, SampleWithPilesOnItsAbundantTaxaIsSpiky) {
    sg::Counts counts;
    for (uint32_t t = 1; t <= 8; t++) AddTaxon(counts, t, 60, 2, 5000);
    for (uint32_t t = 9; t <= 20; t++) AddTaxon(counts, t, 40, 3);
    auto const spikes = SpikesOfAll(counts, 120);
    EXPECT_EQ(spikes.size(), 8u);
    auto const check = sg::Judge(counts, spikes);
    EXPECT_EQ(check.spiked_taxa, 8u);
    EXPECT_EQ(check.top_taxa, 10u);
    EXPECT_EQ(check.top_spiked, 8u);
    EXPECT_GT(check.Share(), 0.5);
    EXPECT_TRUE(check.spiky);
}

// A clean sample whose absent relatives pile on a few genes: the piles are thin taxa, not among the most abundant, and a
// small share of the reads -> no spike mode, though the taxa have spike genes by the per-taxon rule.
TEST(SpikeGenes, SampleWhoseOnlyPilesAreThinRelativesIsNotSpiky) {
    sg::Counts counts;
    for (uint32_t t = 1; t <= 12; t++) AddTaxon(counts, t, 110, 40);   // present, even coverage
    AddTaxon(counts, 50, 8, 2, 120);                                   // two relatives' piles
    AddTaxon(counts, 51, 8, 2, 90);
    auto const spikes = SpikesOfAll(counts, 120);
    EXPECT_EQ(spikes.size(), 2u);
    auto const check = sg::Judge(counts, spikes);
    EXPECT_EQ(check.spiked_taxa, 2u);
    EXPECT_EQ(check.top_spiked, 0u);
    EXPECT_LT(check.Share(), sg::kSampleMinShare);
    EXPECT_FALSE(check.spiky);
}

// One abundant taxon with a pile is not the pattern of a library with carried-over products; nor is a share of the reads
// under kSampleMinShare however many top taxa have piles.
TEST(SpikeGenes, OneSpikyTopTaxonOrASmallShareIsNotSpiky) {
    sg::Counts one;
    AddTaxon(one, 1, 60, 2, 5000);
    for (uint32_t t = 2; t <= 12; t++) AddTaxon(one, t, 100, 30);
    auto const check = sg::Judge(one, SpikesOfAll(one, 120));
    EXPECT_EQ(check.top_spiked, 1u);
    EXPECT_FALSE(check.spiky);

    sg::Counts small;
    for (uint32_t t = 1; t <= 3; t++) AddTaxon(small, t, 60, 1, 60);       // piles of 61 records, about x125 the background
    for (uint32_t t = 4; t <= 6; t++) AddTaxon(small, t, 110, 300);       // abundant taxa, even: the piles are ~1% of reads
    auto const spikes = SpikesOfAll(small, 120);
    auto const small_check = sg::Judge(small, spikes);
    EXPECT_GE(small_check.top_spiked, sg::kSampleMinTopSpiked);
    EXPECT_LT(small_check.Share(), sg::kSampleMinShare);
    EXPECT_FALSE(small_check.spiky);
}

TEST(SpikeGenes, ParsesModes) {
    EXPECT_EQ(sg::ParseMode("auto"), sg::Mode::kAuto);
    EXPECT_EQ(sg::ParseMode("on"), sg::Mode::kOn);
    EXPECT_EQ(sg::ParseMode("off"), sg::Mode::kOff);
    EXPECT_FALSE(sg::ParseMode("yes").has_value());
    EXPECT_EQ(sg::ModeName(sg::Mode::kAuto), "auto");
}

// The records per gene of SAM text: header and unmapped (0x4) and secondary (0x100) records left out, supplementary
// records (a long read's other genes) counted, lines it cannot read skipped.
TEST(SpikeGenes, CountsRecordsPerGeneFromSamText) {
    std::string const text =
        "@HD\tVN:1.6\n"
        "@SQ\tSN:758_71\tLN:2700\n"
        "r1\t99\t758_71\t10\t60\t151M\t=\t200\t341\tACGT\t*\n"
        "r1\t147\t758_71\t200\t60\t151M\t=\t10\t-341\tACGT\t*\n"
        "r2\t0\t758_3\t1\t60\t151M\t*\t0\t0\tACGT\t*\n"
        "r3\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\n"
        "r4\t256\t758_3\t1\t0\t151M\t*\t0\t0\tACGT\t*\n"
        "r5\t2048\t12_7\t5\t40\t900M\t*\t0\t0\tACGT\t*\n"
        "\n"
        "broken line without tabs\n"
        "r6\tx\t758_3\t1\n"
        "r7\t0\t758_71\t1\t60\t151M\t*\t0\t0\tACGT\t*";  // no final newline
    sg::Counts counts;
    sg::CountRecords(text, counts);
    EXPECT_EQ(counts[sg::Key(758, 71)], 3u);
    EXPECT_EQ(counts[sg::Key(758, 3)], 1u);
    EXPECT_EQ(counts[sg::Key(12, 7)], 1u);
    EXPECT_EQ(counts.size(), 3u);
    sg::Counts more;
    sg::Add(more, counts);
    sg::Add(more, counts);
    EXPECT_EQ(more[sg::Key(758, 71)], 6u);
}
