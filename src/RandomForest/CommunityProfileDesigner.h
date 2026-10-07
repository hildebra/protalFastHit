// SPDX-License-Identifier: GPL-2.0-only
#pragma once

#include <cstdint>
#include <memory>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "MetagenomeTypes.h"

namespace protal::sim {

std::string extract_species(const std::string& taxonomy);

// Per-sample output of the strain pre-assignment phase.
struct SampleStrainAssignment {
    std::unordered_map<std::string, std::vector<GenomeRecord>> forced_strains;
    std::unordered_map<std::string, double> species_min_abundance;
};

class CommunityProfileDesigner {
public:
    explicit CommunityProfileDesigner(std::vector<GenomeRecord> genomes);

    std::vector<GenomeAssignment> design_profile(
        const ProfileDesignOptions& options,
        std::mt19937_64& rng) const;

    // Pre-assign strains across all samples according to strain_sharing specs.
    // Returns one SampleStrainAssignment per sample (indexed 0..sample_count-1).
    std::vector<SampleStrainAssignment> assign_strains_across_samples(
        const std::vector<StrainSharingSpec>& specs,
        std::size_t sample_count,
        std::uint64_t total_read_pairs,
        std::uint64_t paired_read_length,
        const std::unordered_map<std::string, std::uint64_t>& genome_lengths,
        std::mt19937_64& rng) const;

private:
    std::vector<GenomeRecord> genomes_;

    // The genomes by species, and the species by genus and by taxon (each lineage token, with and without its rank
    // prefix): the same for every sample of a run, so made once, at the first design (that throws on a lineage of
    // mixed prefixes, as before).
    struct Groups {
        std::unordered_map<std::string, std::vector<GenomeRecord>> grouped;
        std::vector<std::string const*> species;  // grouped's keys, in its iteration order
        std::unordered_map<std::string, std::vector<std::string>> genus_to_species;
        std::unordered_map<std::string, std::string> species_to_genus;
        std::unordered_map<std::string, std::vector<std::string>> taxon_to_species;
        std::unordered_map<std::string, std::unordered_set<std::string>> species_to_taxa;
    };
    mutable std::unique_ptr<Groups> groups_;
    Groups const& groups() const;

    std::unordered_map<std::string, std::vector<GenomeRecord>> group_by_species() const;
    std::vector<double> draw_weights(std::size_t count, const ProfileDesignOptions& options, std::mt19937_64& rng) const;
};

}  // namespace protal::sim
