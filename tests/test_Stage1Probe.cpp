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
    EXPECT_EQ(lines[0], "c\tmin_score\tmin_hits\treads\tkept\tshare_kept\taligned\taligned_kept\trecall\tlookups_left");
    EXPECT_EQ(lines[1], "1\t10\t1\t2\t2\t1.000000\t1\t1\t1.000000\t1.000000");
    // Lookups left: c = 8, T = 14, m = 2 keeps the aligned pair (3 lookups) and makes the unaligned pair's sampled one
    // (level 4), of 5; c = 32, T = 16 keeps neither, and their samples hold one lookup (level 5).
    EXPECT_DOUBLE_EQ(tally.LookupsLeft(3, ScoreIndex(14), HitsIndex(2)), 4.0 / 5.0);
    EXPECT_DOUBLE_EQ(tally.LookupsLeft(5, ScoreIndex(16), HitsIndex(1)), 1.0 / 5.0);

    std::ostringstream taxa;
    tally.WriteTaxa(taxa, [](uint32_t taxid) { return "species " + std::to_string(taxid); });
    auto const taxa_lines = LinesOf(taxa.str());
    ASSERT_EQ(taxa_lines.size(), 2u);
    EXPECT_EQ(taxa_lines[0].rfind("taxid\tname\taligned\tc1_t10_m1\tc1_t10_m2\tc1_t11_m1", 0), 0u) << taxa_lines[0];
    EXPECT_EQ(taxa_lines[1].rfind("77\tspecies 77\t1\t1\t1\t", 0), 0u) << taxa_lines[1];
    // The last column, the stride's c = 32, T = 16, m = 2: 0 (no lookup here has a core index).
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

// The stride pick of 1/c of a read's cores: ceil(cores / c) of them, c apart, the cores left over split between the
// read's ends (at most one more at its end); every read with a core has a pick.
TEST(Stage1Probe, StridePicksEveryCthCoreCentred) {
    EXPECT_EQ(StrideOffset(22, 8), 2u);  // picks 2, 10, 18: two cores before, three after
    EXPECT_TRUE(StridePicked(2, 22, 8));
    EXPECT_TRUE(StridePicked(10, 22, 8));
    EXPECT_TRUE(StridePicked(18, 22, 8));
    EXPECT_FALSE(StridePicked(0, 22, 8));
    EXPECT_FALSE(StridePicked(22, 22, 8));  // past the read
    EXPECT_FALSE(StridePicked(0, 0, 1));
    for (size_t cores = 1; cores <= 300; cores++) {
        for (size_t j = 0; j < kLevels; j++) {
            size_t const c = size_t{1} << j;
            std::vector<size_t> picks;
            for (size_t i = 0; i < cores + c; i++) {
                if (StridePicked(i, cores, c)) picks.push_back(i);
            }
            ASSERT_EQ(picks.size(), (cores + c - 1) / c) << cores << " cores, c " << c;
            for (size_t p = 1; p < picks.size(); p++) ASSERT_EQ(picks[p] - picks[p - 1], c);
            size_t const before = picks.front(), after = cores - 1 - picks.back();
            ASSERT_TRUE(after == before || after == before + 1) << cores << " cores, c " << c;
        }
    }
}

// The stride a read of few cores is sampled at (StrideFor, --sparse_first_min_picks; F1 audit B2): the gate's c while
// the read has min_picks picks at it, else the largest power of two below that gives it them (every core at worst);
// 0 leaves c. StridePicks counts StridePicked's picks. A 75 bp mate (12 cores) at C 4: 3 picks, at 2 six; 150 bp (28
// cores): 7 at C 4.
TEST(Stage1Probe, ShortReadsGetADenserStride) {
    EXPECT_EQ(StrideFor(12, 4, 6), 2u);
    EXPECT_EQ(StrideFor(28, 4, 6), 4u);
    EXPECT_EQ(StrideFor(5, 4, 6), 1u);
    EXPECT_EQ(StrideFor(12, 4, 0), 4u);
    for (size_t cores = 0; cores <= 120; cores++) {
        for (size_t j = 0; j < kLevels; j++) {
            size_t const c = size_t{1} << j;
            size_t picked = 0;
            for (size_t i = 0; i < cores; i++) picked += StridePicked(i, cores, c);
            ASSERT_EQ(StridePicks(cores, c), picked) << cores << " cores, c " << c;
            ASSERT_EQ(StrideFor(cores, c, 0), c);
            for (size_t const min_picks : { size_t{1}, size_t{3}, size_t{6} }) {
                size_t const s = StrideFor(cores, c, min_picks);
                ASSERT_LE(s, c);
                ASSERT_EQ(s & (s - 1), 0u);
                if (s > 1) ASSERT_GE(StridePicks(cores, s), min_picks) << cores << " cores, c " << c;
                if (s < c) ASSERT_LT(StridePicks(cores, 2 * s), min_picks) << cores << " cores, c " << c;
            }
        }
    }
}

// The stride's rules count the lookups of the picked cores of each mate (its own cores, Begin), side by side with the
// hash sample's; a lookup without a core index is in no stride pick.
TEST(Stage1Probe, StrideRulesCountEachMatesPickedCores) {
    Tally tally;
    // Mate 1, 8 cores: c = 2 picks 0, 2, 4, 6; c = 4 picks 1, 5; c = 8, 16, 32 pick 3.
    tally.GetMate(0).Begin(8);
    tally.GetMate(0).Add(0, 16, 10, 5);
    tally.GetMate(0).Add(5, 14, 3, 3);
    // Mate 2, 3 cores: c = 2 picks 0, 2; c = 4 to 32 pick 1.
    tally.GetMate(1).Begin(3);
    tally.GetMate(1).Add(0, 15, 0, 1);
    tally.GetMate(1).Add(0, 16, 30);  // no core index: hash rules only
    tally.EndRead(2, true, { 77 }, { true, true }, true);

    auto stride = [&](size_t c_level, uint8_t min_score, uint8_t min_hits) {
        return tally.Rule(c_level, ScoreIndex(min_score), HitsIndex(min_hits), Pick::stride).second;
    };
    auto hash = [&](size_t c_level, uint8_t min_score, uint8_t min_hits) {
        return tally.Rule(c_level, ScoreIndex(min_score), HitsIndex(min_hits)).second;
    };
    EXPECT_EQ(stride(0, 14, 3), 1.0);  // c = 1: the three lookups with an index
    EXPECT_EQ(stride(0, 14, 5), 0.0);
    EXPECT_EQ(hash(0, 14, 3), 1.0);    // the hash's c = 1 has all four
    EXPECT_EQ(stride(1, 10, 1), 0.0);  // c = 2 picks none of them
    EXPECT_EQ(stride(2, 15, 2), 1.0);  // c = 4: 16 (mate 1's index 5) and 15 (mate 2's index 1)
    EXPECT_EQ(stride(2, 16, 2), 0.0);
    EXPECT_EQ(stride(3, 14, 2), 1.0);  // c = 8: 14 (index 3) and 15
    EXPECT_EQ(stride(3, 15, 2), 0.0);
    EXPECT_EQ(stride(5, 14, 2), 1.0);  // c = 32 alike
    EXPECT_EQ(hash(5, 14, 2), 0.0);    // the hash's c = 32 has the level-5 lookup alone
    EXPECT_EQ(hash(5, 14, 1), 1.0);

    std::ostringstream rules;
    tally.WriteRules(rules, Pick::stride);
    auto const lines = LinesOf(rules.str());
    ASSERT_EQ(lines.size(), kRules + 1);
    EXPECT_EQ(lines[0], "c\tmin_score\tmin_hits\treads\tkept\tshare_kept\taligned\taligned_kept\trecall\tlookups_left");
    // c = 2, T = 10, m = 1: the first rule of c = 2, kept by neither; its pick holds none of the lookups.
    EXPECT_EQ(lines[1 + RuleIndex(1, 0, 0)], "2\t10\t1\t1\t0\t0.000000\t1\t0\t0.000000\t0.000000");
    // Lookups left of the pair's 4: c = 4, T = 16, m = 2 gives it up after its pick's two (indexes 5 and 1); kept, all.
    EXPECT_DOUBLE_EQ(tally.LookupsLeft(2, ScoreIndex(16), HitsIndex(2), Pick::stride), 0.5);
    EXPECT_DOUBLE_EQ(tally.LookupsLeft(2, ScoreIndex(15), HitsIndex(2), Pick::stride), 1.0);
    EXPECT_DOUBLE_EQ(tally.LookupsLeft(5, ScoreIndex(14), HitsIndex(2)), 0.25);

    std::ostringstream taxa;
    tally.WriteTaxa(taxa, [](uint32_t) { return std::string("-"); });
    auto const taxa_lines = LinesOf(taxa.str());
    ASSERT_EQ(taxa_lines.size(), 2u);
    auto split = [](std::string const& line) {
        std::vector<std::string> fields;
        std::istringstream is(line);
        for (std::string field; std::getline(is, field, '\t');) fields.push_back(field);
        return fields;
    };
    auto const header = split(taxa_lines[0]);
    auto const row = split(taxa_lines[1]);
    ASSERT_EQ(header.size(), 3 + 2 * kTaxonRules);
    ASSERT_EQ(row.size(), header.size());
    auto column = [&](std::string const& name) {
        auto const at = std::find(header.begin(), header.end(), name);
        return at == header.end() ? std::string("missing") : row[static_cast<size_t>(at - header.begin())];
    };
    EXPECT_EQ(column("s8_t14_m2"), "1");
    EXPECT_EQ(column("s8_t15_m2"), "0");
    EXPECT_EQ(column("c32_t14_m2"), "0");
    EXPECT_EQ(column("s32_t14_m2"), "1");
    EXPECT_EQ(column("s2_t10_m1"), "0");

    // Joined into another thread's tally.
    Tally other;
    other.Join(tally);
    EXPECT_EQ(other.Rule(3, ScoreIndex(14), HitsIndex(2), Pick::stride).second, 1.0);
    EXPECT_EQ(RuleText(other, 3, ScoreIndex(14), HitsIndex(2), Pick::stride), "stride c=8 T=14 m=2: 100.00% kept, recall 100.00%");
}
