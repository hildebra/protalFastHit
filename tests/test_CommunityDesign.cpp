// Unit tests for simulate_metagenomes' community design (CommunityProfileDesigner): a sample of many
// species with several strains each and few read pairs, as training samples at shallow depths are; congener
// groups; and the sigmas and depths a run's samples take in turn (MetagenomeTypes.h).
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include "RandomForest/CommunityProfileDesigner.h"
#include "RandomForest/MetagenomeSimulator.h"
#include "RandomForest/MetagenomeTypes.h"
#include "RandomForest/PortableRandom.h"

using namespace protal::sim;

namespace {
    std::vector<GenomeRecord> Genomes(int species, int strains) {
        std::vector<GenomeRecord> genomes;
        for (int s = 0; s < species; s++) {
            for (int g = 0; g < strains; g++) {
                std::string const name = "sp" + std::to_string(s);
                genomes.push_back({ "G" + std::to_string(s) + "_" + std::to_string(g),
                                    "d__Bacteria;p__P;c__C;o__O;f__F;g__G" + std::to_string(s) + ";s__G" +
                                    std::to_string(s) + " " + name, "/nonexistent.fna" });
            }
        }
        return genomes;
    }
}

// Every species has three strains in the sample, but most get fewer than three read pairs: each keeps as
// many strains as it has read pairs, instead of the design failing.
TEST(CommunityDesign, RareSpeciesKeepAsManyStrainsAsReadPairs) {
    CommunityProfileDesigner designer(Genomes(20, 3));
    ProfileDesignOptions options;
    options.total_read_pairs = 30;
    options.species_per_sample = 20;
    options.strain_probabilities = { 1.0, 1.0 };
    for (std::uint64_t seed = 1; seed <= 20; seed++) {
        std::mt19937_64 rng(seed);
        auto const assignments = designer.design_profile(options, rng);
        std::map<std::string, std::pair<std::uint64_t, std::size_t>> by_species;  // read pairs, strains
        std::uint64_t total = 0;
        for (auto const& a : assignments) {
            EXPECT_GE(a.read_pairs, 1u) << a.genome.name;
            by_species[a.species].first += a.read_pairs;
            by_species[a.species].second++;
            total += a.read_pairs;
        }
        EXPECT_EQ(total, options.total_read_pairs);
        EXPECT_EQ(by_species.size(), 20u);
        for (auto const& [species, counts] : by_species) {
            EXPECT_LE(counts.second, 3u) << species;
            EXPECT_LE(counts.second, counts.first) << species;
        }
    }
}

namespace {
    // `genera` genera of `per_genus` species each, one genome per species.
    std::vector<GenomeRecord> GeneraOfSpecies(int genera, int per_genus) {
        std::vector<GenomeRecord> genomes;
        for (int g = 0; g < genera; g++) {
            for (int s = 0; s < per_genus; s++) {
                std::string const genus = "G" + std::to_string(g);
                genomes.push_back({ genus + "_" + std::to_string(s),
                                    "d__Bacteria;p__P;c__C;o__O;f__F;g__" + genus + ";s__" + genus + " sp" + std::to_string(s),
                                    "/nonexistent.fna" });
            }
        }
        return genomes;
    }

    // The species of a sample by genus.
    std::map<std::string, std::size_t> SpeciesByGenus(std::vector<GenomeAssignment> const& assignments) {
        std::map<std::string, std::set<std::string>> species;
        for (auto const& a : assignments) species[a.species.substr(0, a.species.find(' '))].insert(a.species);
        std::map<std::string, std::size_t> counts;
        for (auto const& [genus, members] : species) counts[genus] = members.size();
        return counts;
    }
}

// Without congener groups, species are drawn uniformly: among 200 genera of 5, a sample of 20 species has a genus with
// two or more of them by chance (0.73 such genera per sample, 14.7 over 20 samples; about half the samples have one).
// With 0.5:2-4, about half of them come in groups of 2 to 4 congeners, the genera drawn per sample: some three such
// genera in every sample.
TEST(CommunityDesign, CongenerGroupsPutRelativesInEverySample) {
    CommunityProfileDesigner designer(GeneraOfSpecies(200, 5));
    ProfileDesignOptions options;
    options.total_read_pairs = 10000;
    options.species_per_sample = 20;
    std::size_t uniform_pairs = 0;
    std::set<std::string> group_genera;
    for (std::uint64_t seed = 1; seed <= 20; seed++) {
        std::mt19937_64 rng(seed);
        auto const counts = SpeciesByGenus(designer.design_profile(options, rng));
        for (auto const& [genus, n] : counts) uniform_pairs += n >= 2;
    }
    EXPECT_LE(uniform_pairs, 30u) << "uniform draws put few congeners together";  // twice the expected 14.7

    options.congener_share = 0.5;
    options.congener_min = 2;
    options.congener_max = 4;
    std::size_t grouped_pairs = 0;
    for (std::uint64_t seed = 1; seed <= 20; seed++) {
        std::mt19937_64 rng(seed);
        auto const assignments = designer.design_profile(options, rng);
        std::set<std::string> species;
        for (auto const& a : assignments) species.insert(a.species);
        EXPECT_EQ(species.size(), 20u);
        std::size_t grouped = 0;
        for (auto const& [genus, n] : SpeciesByGenus(assignments)) {
            EXPECT_LE(n, 5u);
            if (n >= 2) {
                grouped += n;
                grouped_pairs++;
                group_genera.insert(genus);
            }
        }
        EXPECT_GE(grouped, 8u) << "seed " << seed;   // the target is 10, groups of 2-4, chance pairs aside
        EXPECT_LE(grouped, 14u) << "seed " << seed;
    }
    EXPECT_GT(grouped_pairs, 2 * uniform_pairs) << "the groups, not chance, put the congeners together";
    EXPECT_GE(group_genera.size(), 30u) << "the genera are drawn per sample";
}

// A genus with fewer species than MIN is never a group, and a sample is never more than its species.
TEST(CommunityDesign, CongenerGroupsNeedMinSpeciesAndFitTheSample) {
    CommunityProfileDesigner designer(GeneraOfSpecies(10, 2));
    ProfileDesignOptions options;
    options.total_read_pairs = 1000;
    options.species_per_sample = 6;
    options.congener_share = 1.0;
    options.congener_min = 3;  // no genus has 3
    options.congener_max = 5;
    std::mt19937_64 rng(7);
    auto const assignments = designer.design_profile(options, rng);
    std::set<std::string> species;
    for (auto const& a : assignments) species.insert(a.species);
    EXPECT_EQ(species.size(), 6u);

    options.congener_min = 2;
    options.congener_max = 2;
    std::mt19937_64 rng2(7);
    for (auto const& [genus, n] : SpeciesByGenus(designer.design_profile(options, rng2))) EXPECT_EQ(n, 2u) << genus;
}

TEST(CommunityDesign, CongenerGroupsParse) {
    ProfileDesignOptions options;
    parse_congener_groups("", options);
    EXPECT_EQ(options.congener_share, 0.0);
    parse_congener_groups("0.25:2-5", options);
    EXPECT_EQ(options.congener_share, 0.25);
    EXPECT_EQ(options.congener_min, 2u);
    EXPECT_EQ(options.congener_max, 5u);
    for (std::string bad : { "0.25", "0.25:5-2", "1.5:2-5", "0.25:1-3", "x:2-5", "0.25:2-", "0.25:2-5x", "-0.1:2-3" }) {
        EXPECT_THROW(parse_congener_groups(bad, options), std::runtime_error) << bad;
    }
}

namespace {
    // The read pairs of each species of a sample.
    std::vector<std::uint64_t> PairsBySpecies(std::vector<GenomeAssignment> const& assignments) {
        std::map<std::string, std::uint64_t> pairs;
        for (auto const& a : assignments) pairs[a.species] += a.read_pairs;
        std::vector<std::uint64_t> out;
        for (auto const& [_, n] : pairs) out.push_back(n);
        std::sort(out.begin(), out.end());
        return out;
    }
}

// The lognormal's tail is continuous: few species share the lowest read count (the Poisson-lognormal put ~40-45% of
// them at its lowest weight, 1), and none falls below abundance_floor of the median (here 1/1000; deep samples, so
// that the read counts show the weights).
TEST(CommunityDesign, LognormalTailIsContinuousAndFloored) {
    CommunityProfileDesigner designer(GeneraOfSpecies(100, 5));
    ProfileDesignOptions options;
    options.total_read_pairs = 20'000'000;
    options.species_per_sample = 300;
    options.pln_sigma = 2.5;
    // the share of species at the lowest count, give or take the nudges that make the counts add up
    auto at_lowest = [](std::vector<std::uint64_t> const& pairs) {
        return static_cast<double>(std::count_if(pairs.begin(), pairs.end(), [&](std::uint64_t n) {
                   return n <= pairs.front() + 2; })) / static_cast<double>(pairs.size());
    };
    double tied_lognormal = 0, tied_poisson = 0;
    for (std::uint64_t seed = 1; seed <= 5; seed++) {
        options.distribution = AbundanceDistribution::Lognormal;
        std::mt19937_64 rng(seed);
        auto pairs = PairsBySpecies(designer.design_profile(options, rng));
        ASSERT_EQ(pairs.size(), 300u);
        tied_lognormal += at_lowest(pairs);
        double const median = static_cast<double>(pairs[150]);
        EXPECT_GE(static_cast<double>(pairs.front()) / median, 0.0005) << "seed " << seed;  // the floor, give or take the sample's median
        options.distribution = AbundanceDistribution::PoissonLognormal;
        std::mt19937_64 rng2(seed);
        pairs = PairsBySpecies(designer.design_profile(options, rng2));
        tied_poisson += at_lowest(pairs);
    }
    EXPECT_LT(tied_lognormal / 5, 0.03);
    EXPECT_GT(tied_poisson / 5, 0.30);
    // the floor itself: at a floor of 0.1 a sixth or so of the species (z < ln 0.1 / 2.5 = -0.92) share the lowest weight
    options.distribution = AbundanceDistribution::Lognormal;
    options.abundance_floor = 0.1;
    std::mt19937_64 rng(9);
    auto const pairs = PairsBySpecies(designer.design_profile(options, rng));
    double const floored = static_cast<double>(std::count_if(pairs.begin(), pairs.end(), [&](std::uint64_t n) {
        return n <= pairs.front() + 3; })) / 300.0;  // the floor's count, give or take the rounding's nudges
    EXPECT_NEAR(floored, 0.18, 0.06);
}

// The draws the designs take (PortableRandom.h), the same with every standard library: their distributions.
TEST(PortableRandom, Distributions) {
    std::mt19937_64 rng(11);
    int const n = 200000;
    auto moments = [&](auto draw) {
        double sum = 0, squares = 0;
        for (int i = 0; i < n; ++i) {
            double const x = static_cast<double>(draw());
            sum += x;
            squares += x * x;
        }
        double const mean = sum / n;
        return std::pair<double, double>{mean, squares / n - mean * mean};
    };
    for (double lambda : {0.3, 3.0, 30.0, 1000.0}) {  // multiplication below 10, PTRS above
        auto const [mean, variance] = moments([&] { return portable::Poisson(rng, lambda); });
        EXPECT_NEAR(mean, lambda, 4 * std::sqrt(lambda / n)) << lambda;
        EXPECT_NEAR(variance / lambda, 1.0, 0.03) << lambda;
    }
    for (double shape : {0.5, 2.0, 10.0}) {
        auto const [mean, variance] = moments([&] { return portable::Gamma(rng, shape, 2.0); });
        EXPECT_NEAR(mean / (2 * shape), 1.0, 0.02) << shape;
        EXPECT_NEAR(variance / (4 * shape), 1.0, 0.04) << shape;
    }
    auto const [nb_mean, nb_variance] = moments([&] { return portable::NegativeBinomial(rng, 5, 0.5); });
    EXPECT_NEAR(nb_mean, 5.0, 0.05);
    EXPECT_NEAR(nb_variance, 10.0, 0.3);
    auto const [z_mean, z_variance] = moments([&] { return portable::Normal(rng); });
    EXPECT_NEAR(z_mean, 0.0, 0.01);
    EXPECT_NEAR(z_variance, 1.0, 0.015);
    std::map<std::vector<int>, int> permutations;  // all 6 orders of 3, alike
    for (int i = 0; i < 60000; ++i) {
        std::vector<int> v{1, 2, 3};
        portable::Shuffle(v, rng);
        ++permutations[v];
    }
    ASSERT_EQ(permutations.size(), 6u);
    for (auto const& [order, count] : permutations) EXPECT_NEAR(count, 10000, 400);
    for (int i = 0; i < 1000; ++i) {
        auto const x = portable::Between(rng, 3, 7);
        EXPECT_TRUE(x >= 3 && x <= 7);
    }
    // the same numbers from the same seed
    std::mt19937_64 a(5), b(5);
    for (int i = 0; i < 100; ++i) EXPECT_EQ(portable::Poisson(a, 50.0), portable::Poisson(b, 50.0));
}

// The samples of a run take the abundance sigmas they are given in turn, so that a model does not learn one sigma's
// prior.
TEST(SimulatedSamples, TakeTheirSigmasInTurn) {
    ProfileDesignOptions options;
    options.pln_sigma = 1.3;
    EXPECT_EQ(SigmaForSample(options, 0), 1.3);
    EXPECT_EQ(SigmaForSample(options, 5), 1.3);
    options.pln_sigmas = { 1.3, 2.0 };
    EXPECT_EQ(SigmaForSample(options, 0), 1.3);
    EXPECT_EQ(SigmaForSample(options, 1), 2.0);
    EXPECT_EQ(SigmaForSample(options, 2), 1.3);
    EXPECT_EQ(SigmaForSample(options, 7), 2.0);
}

// And their depths: collect_training_data.py gives a scenario's samples depths of their own (--total_read_pairs a,b,c).
TEST(SimulatedSamples, TakeTheirDepthsInTurn) {
    ProfileDesignOptions options;
    options.total_read_pairs = 1000;
    EXPECT_EQ(ReadPairsForSample(options, 0), 1000u);
    EXPECT_EQ(ReadPairsForSample(options, 3), 1000u);
    options.total_read_pairs_per_sample = { 500, 2000, 1200 };
    EXPECT_EQ(ReadPairsForSample(options, 0), 500u);
    EXPECT_EQ(ReadPairsForSample(options, 1), 2000u);
    EXPECT_EQ(ReadPairsForSample(options, 2), 1200u);
    EXPECT_EQ(ReadPairsForSample(options, 4), 2000u);
}
