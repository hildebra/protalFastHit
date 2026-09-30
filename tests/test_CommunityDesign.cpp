// Unit tests for simulate_metagenomes' community design (CommunityProfileDesigner): a sample of many
// species with several strains each and few read pairs, as training samples at shallow depths are.
#include <gtest/gtest.h>
#include <map>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include "RandomForest/CommunityProfileDesigner.h"

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
