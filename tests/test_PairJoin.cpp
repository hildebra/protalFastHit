// Unit tests for joining the candidate alignments of two mates into pairs (classify::JoinAlignmentPairs):
// every pair of candidates on one gene in the right orientation, and each candidate that pairs with none on its
// own (the reference below), the same whether each alignment is moved into its last pair or copied.
#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <memory>
#include <random>
#include <string>
#include <tuple>
#include <vector>
#include <unistd.h>
#include "Classify.h"
#include "TestReference.h"

using namespace protal;

namespace {
    // The join takes a genome loader and does not use it: one of a one-gene reference.
    struct OneGeneLoader : test::LoadedReference {
        OneGeneLoader() : LoadedReference({ { 1, { std::string(100, 'A') } } }, "pair join") {}
    };

    // The join written plainly, as the definition: a mate-2 candidate that pairs is not also on its own (until
    // 2026-10-01 the join marked the candidate before it after one in the wrong orientation).
    void ReferenceJoin(PairedAlignmentResultList& pairs, AlignmentResultList& read1, AlignmentResultList& read2) {
        std::vector<bool> selected2(read2.size(), false);
        for (auto& alignment1 : read1) {
            bool paired = false;
            for (size_t j = 0; j < read2.size(); j++) {
                auto& alignment2 = read2[j];
                if (alignment1.Taxid() == alignment2.Taxid() && alignment1.GeneId() == alignment2.GeneId()) {
                    if (!CorrectOrientation(alignment1, alignment2)) continue;
                    pairs.emplace_back(PairedAlignment{ alignment1, alignment2 });
                    selected2[j] = true;
                    paired = true;
                }
            }
            if (!paired) pairs.emplace_back(PairedAlignment{ alignment1, AlignmentResult() });
        }
        for (size_t i = 0; i < selected2.size(); i++) {
            if (!selected2[i]) pairs.emplace_back(PairedAlignment{ AlignmentResult(), read2[i] });
        }
    }

    AlignmentResultList Candidates(std::mt19937& rng, size_t n) {
        AlignmentResultList list;
        for (size_t i = 0; i < n; i++) {
            AlignmentResult r;
            r.Set(1 + rng() % 3, 1 + rng() % 3, static_cast<int32_t>(rng() % 1000), rng() % 2 == 0);
            r.GetAlignmentInfo().cigar = std::string(100 + rng() % 60, "MXID"[rng() % 4]) + std::to_string(i);  // long: not inline
            r.GetAlignmentInfo().compressed_cigar = std::to_string(i) + "M";
            list.push_back(r);
        }
        return list;
    }

    using Key = std::tuple<bool, uint32_t, uint32_t, int32_t, bool, std::string, std::string>;
    Key KeyOf(AlignmentResult const& r) {
        return { r.IsSet(), r.Taxid(), r.GeneId(), r.GenePos(), r.Forward(), r.GetAlignmentInfo().cigar, r.GetAlignmentInfo().compressed_cigar };
    }
}

TEST(PairJoin, MovingGivesThePairsCopyingGave) {
    std::mt19937 rng(3);
    OneGeneLoader genomes;
    GenomeLoader& loader = *genomes.loader;
    size_t pairs_seen = 0, skipped_seen = 0;
    for (int round = 0; round < 2000; round++) {
        auto const read1 = Candidates(rng, rng() % 7), read2 = Candidates(rng, rng() % 7);
        AlignmentResultList r1 = read1, r2 = read2;
        PairedAlignmentResultList expected;
        ReferenceJoin(expected, r1, r2);
        for (bool consume : { false, true }) {
            AlignmentResultList c1 = read1, c2 = read2;
            PairedAlignmentResultList pairs;
            classify::JoinAlignmentPairs(pairs, c1, c2, loader, consume);
            ASSERT_EQ(pairs.size(), expected.size()) << "round " << round << " consume " << consume;
            for (size_t p = 0; p < pairs.size(); p++) {
                EXPECT_EQ(KeyOf(pairs[p].first), KeyOf(expected[p].first)) << "round " << round << " pair " << p;
                EXPECT_EQ(KeyOf(pairs[p].second), KeyOf(expected[p].second)) << "round " << round << " pair " << p;
            }
            if (!consume) {  // without consuming, the candidates are left as they were
                for (size_t i = 0; i < c1.size(); i++) EXPECT_EQ(KeyOf(c1[i]), KeyOf(read1[i]));
                for (size_t i = 0; i < c2.size(); i++) EXPECT_EQ(KeyOf(c2[i]), KeyOf(read2[i]));
            }
        }
        pairs_seen += expected.size();
        // a wrong-orientation candidate before a right one: the case where the paired mark lagged
        for (auto const& a : read1) {
            bool wrong = false;
            for (auto const& b : read2) {
                if (a.Taxid() != b.Taxid() || a.GeneId() != b.GeneId()) continue;
                if (!CorrectOrientation(a, b)) wrong = true;
                else if (wrong) skipped_seen++;
            }
        }
    }
    EXPECT_GT(pairs_seen, 5000u);
    EXPECT_GT(skipped_seen, 50u);
}

// Mate 2 has a candidate on mate 1's gene in the wrong orientation before one in the right orientation: the right one
// pairs and is not on its own as well; the wrong one, and one on another gene, are on their own.
TEST(PairJoin, ACandidateThatPairsIsNotAlsoOnItsOwn) {
    OneGeneLoader genomes;
    auto candidate = [](uint32_t gene, bool forward, std::string const& name) {
        AlignmentResult r;
        r.Set(1, gene, 100, forward);
        r.GetAlignmentInfo().compressed_cigar = name;
        return r;
    };
    AlignmentResultList const read1 = { candidate(1, true, "a") };
    AlignmentResultList const read2 = { candidate(1, true, "wrong"), candidate(1, false, "right"), candidate(2, false, "other") };
    for (bool consume : { false, true }) {
        AlignmentResultList c1 = read1, c2 = read2;
        PairedAlignmentResultList pairs;
        classify::JoinAlignmentPairs(pairs, c1, c2, *genomes.loader, consume);
        std::vector<std::pair<std::string, std::string>> names;
        for (auto const& [first, second] : pairs) {
            names.emplace_back(first.IsSet() ? first.GetAlignmentInfo().compressed_cigar : "-",
                               second.IsSet() ? second.GetAlignmentInfo().compressed_cigar : "-");
        }
        std::vector<std::pair<std::string, std::string>> const expected = { { "a", "right" }, { "-", "wrong" }, { "-", "other" } };
        EXPECT_EQ(names, expected) << "consume " << consume;
    }
}
