// Unit tests for the genes' conservation factors (SequenceUtils/GeneConservation.h, gene_conservation.tsv): the factors
// that scale the depth identity margin per gene, as --build estimates and writes them and queries read them.
#include <gtest/gtest.h>
#include <cmath>
#include <cstdint>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include "SequenceUtils/GeneConservation.h"
#include "TestUtil.h"

using namespace protal;
using protal::test::Mutated;
using protal::test::RandomSequence;

TEST(GeneConservation, TheMarginIsFixedForReadErrorsAndScaledForTheRest) {
    using gene_conservation::GeneMargin;
    EXPECT_NEAR(GeneMargin(0.08, 1.0), 0.08, 1e-12);  // a gene of typical conservation: the margin
    EXPECT_NEAR(GeneMargin(0.08, 0.4), 0.05, 1e-12);  // 0.03 + 0.05 x 0.4
    EXPECT_NEAR(GeneMargin(0.08, 1.6), 0.11, 1e-12);
    EXPECT_NEAR(GeneMargin(0.02, 3.0), 0.02, 1e-12);  // below the fixed part: not scaled
    EXPECT_EQ(GeneMargin(1.0, 0.4), 1.0);             // every read counts
}

TEST(GeneConservation, TheTableReadsWhatItWrites) {
    gene_conservation::Table table;
    table.Set(3, 0.4, 25);
    table.Set(120, 1.62, 40);
    EXPECT_EQ(table.Genes(), 2u);
    std::stringstream ss;
    table.Write(ss);
    EXPECT_EQ(ss.str(), "geneid\tfactor\tspecies\n3\t0.4000\t25\n120\t1.6200\t40\n");
    gene_conservation::Table read;
    ASSERT_EQ(read.Read(ss), "");
    EXPECT_EQ(read.Genes(), 2u);
    EXPECT_NEAR(read.Factor(3), 0.4, 1e-6);
    EXPECT_NEAR(read.Factor(120), 1.62, 1e-6);
    EXPECT_EQ(read.Factor(4), 1.0);        // a gene without a factor
    EXPECT_EQ(read.Factor(100000), 1.0);
    auto const [low, high] = read.Range();
    EXPECT_NEAR(low, 0.4, 1e-6);
    EXPECT_NEAR(high, 1.62, 1e-6);
    EXPECT_EQ(gene_conservation::Table().Range(), (std::pair<double, double>{1, 1}));

    std::istringstream two_columns("# comment\n7\t0.9\n");  // the species column is optional
    gene_conservation::Table short_table;
    EXPECT_EQ(short_table.Read(two_columns), "");
    EXPECT_NEAR(short_table.Factor(7), 0.9, 1e-6);
}

TEST(GeneConservation, ReadingStopsAtTheFirstBadLine) {
    auto problem = [](std::string const& content) {
        std::istringstream is(content);
        gene_conservation::Table table;
        return table.Read(is);
    };
    EXPECT_EQ(problem("geneid\tfactor\tspecies\n1\t1.0\t3\nx\t1.0\t3\n"), "line 3: the gene id is not a number");
    EXPECT_EQ(problem("1 1.0\n"), "line 1: expected a gene id and a factor, separated by a tab");
    EXPECT_EQ(problem("1\tfast\n"), "line 1: the factor is not a number");
    EXPECT_EQ(problem("1\t0\n"), "line 1: the factor must be above 0 and at most 100");
    EXPECT_EQ(problem("1\t-0.5\n"), "line 1: the factor must be above 0 and at most 100");
    EXPECT_EQ(problem("1\tnan\n"), "line 1: the factor is not a number");
    EXPECT_EQ(problem("1\t1.0\n1\t1.2\n"), "line 2: gene 1 is listed twice");
    EXPECT_EQ(problem("1048576\t1.0\n"), "line 1: gene id 1048576 is too large");
}

TEST(GeneConservation, TheMashDistanceEstimatesTheShareOfDifferentBases) {
    std::mt19937 rng(7);
    std::string const gene = RandomSequence(3000, rng);
    auto const kmers = gene_conservation::Kmers(gene);
    EXPECT_EQ(gene_conservation::MashDistance(kmers, kmers), 0.0);
    for (double rate : { 0.005, 0.02, 0.05 }) {
        double const d = gene_conservation::MashDistance(kmers, gene_conservation::Kmers(Mutated(gene, rate, rng)));
        EXPECT_NEAR(d, rate, 0.25 * rate + 0.002) << "rate " << rate;
    }
    EXPECT_GT(gene_conservation::MashDistance(kmers, gene_conservation::Kmers(RandomSequence(3000, rng))), 0.25);
    EXPECT_EQ(gene_conservation::Kmers("ACGTNACGTACGTACGTAC").size(), 3u);  // the 14 bases after the N
    EXPECT_EQ(gene_conservation::Kmers("ACGTNACGTACG").size(), 0u);         // no 12 bases without an N
}

TEST(GeneConservation, TheEstimateFollowsHowFastEachGeneDiverges) {
    // 30 species of 21 genes: genes 1, 4, ... diverge at 0.4 times a species' rate, genes 2, 5, ... at 1,
    // genes 3, 6, ... at 1.6; each species' 3 other genomes 0.5-3% from its representative at a gene
    // of rate 1. The full reference lists the representative's own copy too.
    std::mt19937 rng(11);
    std::vector<double> const rate = { 0.4, 1.0, 1.6 };
    std::vector<uint64_t> keys;
    std::map<std::pair<uint64_t, uint64_t>, std::string> reference;
    for (uint64_t taxid = 1; taxid <= 30; taxid++) {
        for (uint64_t gene = 1; gene <= 21; gene++) {
            keys.push_back(gene_conservation::Estimator::Key(taxid, gene));
            reference[{taxid, gene}] = RandomSequence(900, rng);
        }
    }
    // A species with identical genomes and one with only 5 genes inform nothing.
    for (uint64_t gene = 1; gene <= 21; gene++) keys.push_back(gene_conservation::Estimator::Key(31, gene));
    for (uint64_t gene = 1; gene <= 5; gene++) keys.push_back(gene_conservation::Estimator::Key(32, gene));
    gene_conservation::Estimator estimator(keys);
    for (uint64_t taxid = 1; taxid <= 30; taxid++) {
        double const divergence = 0.005 + 0.025 * (taxid - 1) / 29.0;
        for (uint64_t gene = 1; gene <= 21; gene++) {
            auto const& rep = reference[{taxid, gene}];
            auto slot = estimator.Take(taxid, gene);
            ASSERT_TRUE(slot);
            estimator.Add(*slot, rep, rep);  // the representative's own copy
            for (int genome = 0; genome < 3; genome++) {
                slot = estimator.Take(taxid, gene);
                ASSERT_TRUE(slot);
                estimator.Add(*slot, rep, Mutated(rep, divergence * rate[(gene - 1) % 3], rng));
            }
        }
    }
    std::string const same = RandomSequence(900, rng);
    for (uint64_t gene = 1; gene <= 21; gene++) {
        for (int genome = 0; genome < 3; genome++) estimator.Add(*estimator.Take(31, gene), same, same);
    }
    for (uint64_t gene = 1; gene <= 5; gene++) estimator.Add(*estimator.Take(32, gene), same, Mutated(same, 0.02, rng));
    EXPECT_FALSE(estimator.Take(33, 1)) << "not a reference gene";

    auto const estimate = estimator.Finish();
    EXPECT_EQ(estimate.species, 30u);
    EXPECT_EQ(estimate.species_with_copies, 32u);
    EXPECT_EQ(estimate.table.Genes(), 21u);
    for (uint64_t gene = 1; gene <= 21; gene++) {
        double const truth = rate[(gene - 1) % 3];
        double const shrunk = (30 * truth + gene_conservation::kPrior) / (30 + gene_conservation::kPrior);
        EXPECT_NEAR(estimate.table.Factor(gene), shrunk, 0.12) << "gene " << gene << ", rate " << truth;
    }
}

TEST(GeneConservation, CongenersDifferMostOnTheFastGenes) {
    // Genus 1: 6 species from one ancestor, each 4% from it at a gene of rate 1 (8% between two species): genes 1-6 at
    // rate 0.3, 7-12 at 1.7, and gene 13 the same in every species. Genus 2 has one species: nothing to compare.
    std::mt19937 rng(11);
    std::vector<double> const rates = { 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 0 };
    std::vector<std::string> ancestor;
    for (size_t g = 0; g < rates.size(); g++) ancestor.push_back(RandomSequence(1200, rng));
    std::map<uint32_t, std::vector<std::pair<uint64_t, std::string>>> species;
    for (uint32_t taxid = 1; taxid <= 6; taxid++) {
        for (size_t g = 0; g < rates.size(); g++) species[taxid].emplace_back(g + 1, Mutated(ancestor[g], 0.04 * rates[g], rng));
    }
    species[7] = species[1];
    gene_conservation::Table within;
    for (uint64_t g = 1; g <= 12; g++) within.Set(g, g <= 6 ? 0.5 : 1.5, 6);
    within.Set(13, 0.5, 6);
    size_t calls = 0;
    auto const estimate = gene_conservation::CompareCongeners({ { 1, 2, 3, 4, 5, 6 }, { 7 } }, [&](uint32_t taxid) {
        calls++;
        return species.at(taxid);
    }, within, 1);
    EXPECT_EQ(calls, 6u);  // the species of genus 1, once each
    EXPECT_EQ(estimate.genera, 1u);
    EXPECT_EQ(estimate.species, 6u);
    EXPECT_EQ(estimate.pairs, 15u);  // each of 6 species against the next 4, cyclically: all 15 pairs
    ASSERT_EQ(estimate.genes.size(), 13u);
    // A pair's median gene is a slow one (6 slow, 6 fast, 1 the same): slow genes about 1, fast ones about 5.7.
    for (auto const& g : estimate.genes) {
        if (g.geneid <= 6) EXPECT_NEAR(g.between, 1.0, 0.4) << g.geneid;
        else if (g.geneid <= 12) EXPECT_GT(g.between, 3.0) << g.geneid;
        else EXPECT_EQ(g.between, 0.0);
        EXPECT_EQ(g.pairs, 15u) << g.geneid;
        EXPECT_EQ(g.species, 6u) << g.geneid;
        EXPECT_EQ(g.identical, g.geneid == 13 ? 6u : 0u) << g.geneid;
    }
    EXPECT_GT(estimate.spearman, 0.8);
    EXPECT_EQ(estimate.correlated, 13u);
    EXPECT_NEAR(estimate.conserved_between, 1.0, 0.4);
    EXPECT_GT(estimate.fast_between, 3.0);
    EXPECT_NEAR(estimate.conserved_identical, 1.0 / 7, 1e-12);  // gene 13 of the 7 conserved genes
    EXPECT_EQ(estimate.fast_identical, 0.0);
    std::ostringstream os;
    estimate.Write(os);
    EXPECT_EQ(os.str().substr(0, os.str().find('\n')),
              "geneid\twithin_factor\tbetween_factor\tpairs\tspecies\tidentical_share\tnear_identical_share");

    EXPECT_NEAR(gene_conservation::Spearman({ 1, 2, 3, 4 }, { 10, 20, 30, 40 }), 1.0, 1e-12);
    EXPECT_NEAR(gene_conservation::Spearman({ 1, 2, 3, 4 }, { 4, 3, 2, 1 }), -1.0, 1e-12);
    EXPECT_TRUE(std::isnan(gene_conservation::Spearman({ 1, 2 }, { 1, 2 })));
}

TEST(GeneConservation, CopiesBeyondTheCapAreNotCompared) {
    gene_conservation::Estimator estimator({ gene_conservation::Estimator::Key(1, 1) });
    for (uint32_t i = 0; i < gene_conservation::kMaxCopies; i++) EXPECT_TRUE(estimator.Take(1, 1));
    EXPECT_FALSE(estimator.Take(1, 1));
    EXPECT_TRUE(estimator.Finish().table.Empty());  // one species, one gene: no factor
}
