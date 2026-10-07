// Unit tests for simulate_metagenomes' community design (CommunityProfileDesigner): a sample of many
// species with several strains each and few read pairs, as training samples at shallow depths are; congener
// groups; and the sigmas and depths a run's samples take in turn (MetagenomeTypes.h).
#include <gtest/gtest.h>
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
