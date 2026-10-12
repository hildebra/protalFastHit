// Unit tests for the tallies of --stage1_probe (Stage1Probe.h): the cores' levels, the rules a read's lookups pass, the
// taxa table, the mates' first lookups, and joining the threads' tallies.
#include <gtest/gtest.h>
#include <algorithm>
#include <array>
#include <cstdint>
#include <sstream>
#include <string>
#include <vector>
#include "Core/Stage1Probe.h"

using namespace protal::stage1_probe;

namespace {
    size_t ScoreIndex(uint8_t min_score) {
        for (size_t t = 0; t < kMinScores.size(); t++) {
            if (kMinScores[t] == min_score) return t;
        }
        return SIZE_MAX;
    }

    size_t HitsIndex(uint8_t min_hits) {
        for (size_t h = 0; h < kMinHits.size(); h++) {
            if (kMinHits[h] == min_hits) return h;
        }
        return SIZE_MAX;
    }

    std::vector<std::string> LinesOf(std::string const& text) {
        std::vector<std::string> lines;
        std::istringstream is(text);
        for (std::string line; std::getline(is, line);) lines.push_back(line);
        return lines;
    }
}

// A core is in the sparse index of 1/2^j with a share of about 2^-j, the same core always at the same level.
TEST(Stage1Probe, CoreLevelsHalveWithEachLevel) {
    std::array<uint64_t, kLevels> at_least{};
    uint64_t const n = 1u << 20;
    for (uint64_t core = 0; core < n; core++) {
        uint8_t const level = Level(core * 2654435761u);
        ASSERT_LT(level, kLevels);
        for (size_t j = 0; j <= level; j++) at_least[j]++;
    }
    EXPECT_EQ(at_least[0], n);
    for (size_t j = 1; j < kLevels; j++) {
        double const share = static_cast<double>(at_least[j]) / static_cast<double>(n);
        EXPECT_NEAR(share, 1.0 / static_cast<double>(1u << j), 0.01 / static_cast<double>(1u << j) + 0.002) << "level " << j;
    }
    EXPECT_EQ(Level(12345), Level(12345));
}

// A pair's lookups over both mates against each rule: the share of reads kept and of those with records (recall), the
// taxa table, and what joining two threads' tallies gives.
TEST(Stage1Probe, RulesCountTheLookupsOfBothMates) {
    Tally tally;
    // Aligned pair: mate 1 has a lookup of level 3 scoring 16 and one of level 0 scoring 12, mate 2 one of level 5
    // scoring 14.
    tally.GetMate(0).Add(3, 16, 10);
    tally.GetMate(0).Add(0, 12, 40);
    tally.GetMate(1).Add(5, 14, 7);
    tally.EndRead(2, true, { 77 }, { true, true }, true);
    // Unaligned pair: one chance match each (level 1, score 9; level 4, score 11).
    tally.GetMate(0).Add(1, 9, 3);
    tally.GetMate(1).Add(4, 11, 5);
    tally.EndRead(2, false, {}, { false, false }, true);

    EXPECT_EQ(tally.Reads(), 2u);
    EXPECT_EQ(tally.Aligned(), 1u);
    auto rule = [&](size_t c_level, uint8_t min_score, uint8_t min_hits) {
        return tally.Rule(c_level, ScoreIndex(min_score), HitsIndex(min_hits));
    };
    // c = 1, T = 10: the aligned pair's three lookups, and the unaligned pair's 11.
    EXPECT_EQ(rule(0, 10, 1), (std::pair<double, double>{ 1.0, 1.0 }));
    EXPECT_EQ(rule(0, 10, 2), (std::pair<double, double>{ 0.5, 1.0 }));
    EXPECT_EQ(rule(0, 10, 3), (std::pair<double, double>{ 0.5, 1.0 }));
    EXPECT_EQ(rule(0, 10, 5), (std::pair<double, double>{ 0.0, 0.0 }));
    EXPECT_EQ(rule(0, 12, 1), (std::pair<double, double>{ 0.5, 1.0 }));
    // c = 8 (level 3 and up): 16 and 14, of both mates.
    EXPECT_EQ(rule(3, 14, 2), (std::pair<double, double>{ 0.5, 1.0 }));
    EXPECT_EQ(rule(3, 16, 2), (std::pair<double, double>{ 0.0, 0.0 }));
    EXPECT_EQ(rule(3, 16, 1), (std::pair<double, double>{ 0.5, 1.0 }));
    // c = 32: only mate 2's 14; the unaligned pair's level-4 lookup is not in it.
    EXPECT_EQ(rule(5, 14, 1), (std::pair<double, double>{ 0.5, 1.0 }));
    EXPECT_EQ(rule(5, 15, 1), (std::pair<double, double>{ 0.0, 0.0 }));
    EXPECT_EQ(rule(4, 11, 1), (std::pair<double, double>{ 1.0, 1.0 }));

    std::ostringstream rules;
    tally.WriteRules(rules);
    auto const lines = LinesOf(rules.str());
    ASSERT_EQ(lines.size(), kRules + 1);
    EXPECT_EQ(lines[0], "c\tmin_score\tmin_hits\treads\tkept\tshare_kept\taligned\taligned_kept\trecall");
    EXPECT_EQ(lines[1], "1\t10\t1\t2\t2\t1.000000\t1\t1\t1.000000");

    std::ostringstream taxa;
    tally.WriteTaxa(taxa, [](uint32_t taxid) { return "species " + std::to_string(taxid); });
    auto const taxa_lines = LinesOf(taxa.str());
    ASSERT_EQ(taxa_lines.size(), 2u);
    EXPECT_EQ(taxa_lines[0].rfind("taxid\tname\taligned\tc1_t10_m1\tc1_t10_m2\tc1_t11_m1", 0), 0u) << taxa_lines[0];
    EXPECT_EQ(taxa_lines[1].rfind("77\tspecies 77\t1\t1\t1\t", 0), 0u) << taxa_lines[1];
    // The last column, c = 32, T = 16, m = 2: 0.
    EXPECT_EQ(taxa_lines[1].substr(taxa_lines[1].size() - 2), "\t0");

    // Another thread's tally: an aligned pair without a lookup in the index (kept by no rule).
    Tally other;
    other.EndRead(2, true, { 77 }, { true, false }, true);
    tally.Join(other);
    EXPECT_EQ(tally.Reads(), 3u);
    EXPECT_EQ(tally.Aligned(), 2u);
    EXPECT_EQ(rule(0, 10, 1).second, 0.5);
    std::ostringstream joined;
    tally.WriteTaxa(joined, [](uint32_t) { return std::string("-"); });
    EXPECT_EQ(LinesOf(joined.str())[1].rfind("77\t-\t2\t1\t1\t", 0), 0u);
}

// A mate's first lookups are those of its first read positions, whatever order they were scored in; the table counts
// mates by the best score among their first k and by state.
TEST(Stage1Probe, FirstLookupsAreInReadOrder) {
    Tally tally;
    tally.GetMate(0).Add(0, 9, 30);   // scored in block-size order, not read order
    tally.GetMate(0).Add(0, 16, 10);
    tally.GetMate(0).Add(0, 11, 20);
    tally.GetMate(1).Add(0, 8, 5);
    tally.GetMate(1).Add(0, 16, 40);
    tally.EndRead(2, true, { 5 }, { true, false }, true);
    // A single-end read that does not align, without lookups.
    tally.EndRead(1, false, {}, { false, false }, true);
    // A long read is not counted there.
    tally.GetMate(0).Add(0, 16, 1);
    tally.EndRead(1, true, { 5, 6 }, { true, false }, false);

    EXPECT_EQ(tally.FirstLookupMates(), 3u);
    std::ostringstream os;
    tally.WriteFirstLookups(os);
    auto const lines = LinesOf(os.str());
    // first 1: mate 1's position 10 (16, aligned), mate 2's position 5 (8, its fragment aligned), the single read (0).
    EXPECT_NE(std::find(lines.begin(), lines.end(), "1\t16\t0\t0\t1"), lines.end()) << os.str();
    EXPECT_NE(std::find(lines.begin(), lines.end(), "1\t8\t0\t1\t0"), lines.end()) << os.str();
    EXPECT_NE(std::find(lines.begin(), lines.end(), "1\t0\t1\t0\t0"), lines.end()) << os.str();
    // first 2: mate 2 reaches its 16.
    EXPECT_NE(std::find(lines.begin(), lines.end(), "2\t16\t0\t1\t1"), lines.end()) << os.str();
    // Giving up mates whose first lookup scores 11 or less: the single read and mate 2 (of 3 mates); of the aligned
    // mates none; of the mates of aligned fragments one of two.
    auto const give_up = tally.GiveUp(0, 11);
    EXPECT_DOUBLE_EQ(give_up[0], 2.0 / 3.0);
    EXPECT_DOUBLE_EQ(give_up[1], 0.0);
    EXPECT_DOUBLE_EQ(give_up[2], 0.5);
    // The long read counts in the rules and the taxa (both of its taxa).
    EXPECT_EQ(tally.Reads(), 3u);
    std::ostringstream taxa;
    tally.WriteTaxa(taxa, [](uint32_t) { return std::string(); });
    EXPECT_EQ(LinesOf(taxa.str()).size(), 3u);
}
