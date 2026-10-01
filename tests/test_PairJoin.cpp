// Unit tests for joining the candidate alignments of two mates into pairs (classify::JoinAlignmentPairs):
// moving each alignment into its last pair gives the pairs that copying gave, in the same order, as the join
// did before it moved anything (the reference below), including its marking of mate-2 candidates as paired.
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

using namespace protal;

namespace {
    // The join takes a genome loader and does not use it: one of a one-gene reference.
    struct OneGeneLoader {
        std::filesystem::path dir = std::filesystem::temp_directory_path() / ("protal pair join test " + std::to_string(::getpid()));
        std::unique_ptr<GenomeLoader> loader;
        OneGeneLoader() {
            std::filesystem::create_directories(dir);
            std::string const header = ">1_1\n", gene(100, 'A');
            std::ofstream(dir / "reference.fna", std::ios::binary) << header << gene << '\n';
            std::ofstream(dir / "reference.map", std::ios::binary) << "1\t1\t" << header.size() << '\t' << header.size() + gene.size() << '\n';
            loader = std::make_unique<GenomeLoader>((dir / "reference.fna").string(), (dir / "reference.map").string());
        }
        ~OneGeneLoader() {
            loader.reset();
            std::filesystem::remove_all(dir);
        }
    };

    // The join as it was before candidates were moved, kept here as the definition.
    void ReferenceJoin(PairedAlignmentResultList& pairs, AlignmentResultList& read1, AlignmentResultList& read2) {
        std::vector<bool> selected2(read2.size(), false);
        for (auto& alignment1 : read1) {
            bool paired = false;
            size_t read2_index = 0;
            for (auto& alignment2 : read2) {
                if (alignment1.Taxid() == alignment2.Taxid() && alignment1.GeneId() == alignment2.GeneId()) {
                    if (!CorrectOrientation(alignment1, alignment2)) continue;
                    pairs.emplace_back(PairedAlignment{ alignment1, alignment2 });
                    selected2[read2_index] = true;
                    paired = true;
                }
                read2_index++;
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
        // a wrong-orientation candidate before a right one: the case where the paired mark lags
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
