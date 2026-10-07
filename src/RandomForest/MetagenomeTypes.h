// SPDX-License-Identifier: GPL-2.0-only
#pragma once

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace protal::sim {

struct GenomeRecord {
    std::string name;
    std::string taxonomy;
    std::filesystem::path fasta_path;
    std::optional<std::uint64_t> genome_length;
};

struct GenomeAssignment {
    GenomeRecord genome;
    std::string species;
    std::uint64_t read_pairs{};
    double relative_abundance{0.0};
    double vertical_coverage{0.0};
    std::uint64_t genome_length{0};
    // ART's --rndSeed for this genome. Recorded in the manifest so that a replay
    // (--from_manifest) reproduces the reads themselves, not just the composition.
    std::optional<std::uint64_t> art_seed;
};

enum class AbundanceDistribution {
    PowerLaw,
    NegativeBinomial,
    PoissonLognormal
};

// Cross-sample strain sharing config for a single species.
struct StrainSharingSpec {
    std::string species;            // species name (resolved against genome table)
    double sample_fraction{1.0};    // fraction of all samples that include this species [0,1]
    std::size_t n_strains{1};       // total distinct strains drawn for this species across all samples
    std::size_t min_occurrence{1};  // each drawn strain must appear in >= this many samples
    double min_vcov{0.0};           // minimum vertical coverage per strain per sample (0 = no floor)
    // Probabilities of adding a 2nd, 3rd, ... conspecific strain to a sample from the forced pool.
    // Same semantics as --strains_per_species but scoped to the n_strains drawn for this species.
    std::vector<double> conspecific_strain_probabilities;
};

struct ProfileDesignOptions {
    std::uint64_t total_read_pairs{100'000};
    std::vector<std::uint64_t> total_read_pairs_per_sample;  // if not empty: the samples' read pairs in turn (sample i
                                                             // gets the i-th, cyclically), so that one run's samples
                                                             // differ in depth (ReadPairsForSample)
    std::size_t species_per_sample{10};
    std::size_t species_per_sample_min{0};  // if > 0, count is drawn uniformly from [min, species_per_sample] per sample
    std::vector<std::size_t> species_per_sample_list;  // if not empty: the samples' species counts in turn (sample i
                                                       // gets the i-th, cyclically), in place of the range: a caller
                                                       // that spreads them itself (collect_training_data.py's scenarios)
    AbundanceDistribution distribution{AbundanceDistribution::PoissonLognormal};
    double powerlaw_alpha{2.0};
    int negative_binomial_r{5};
    double negative_binomial_p{0.5};
    double pln_mu{0.0};
    double pln_sigma{1.3};
    std::vector<double> pln_sigmas;            // if not empty: the samples' sigmas in turn (sample i gets the i-th,
                                               // cyclically), so that one design mixes abundance distributions and a
                                               // model does not learn one sigma's prior (SigmaForSample)
    std::vector<std::string> include_species;  // species that must be present in each sample
    std::unordered_map<std::string, std::size_t> genus_species_counts;  // requested species counts per genus
    std::vector<double> strain_probabilities;  // probabilities for adding 2nd, 3rd, ... strain of a species
    std::unordered_map<std::string, std::size_t> taxon_species_counts;  // requested species counts per taxon token
    bool pick_random_demand_if_fail{false};  // if true, cap genus/taxon demands to species_per_sample instead of failing
    // Congener groups (--congener_groups SHARE:MIN-MAX): about SHARE of each sample's species come in groups of MIN to
    // MAX species of one genus, the genera drawn per sample among those with MIN species or more. Without them
    // (SHARE 0) the species are drawn uniformly, and among thousands of genera congeners hardly ever share a sample,
    // while in real samples they often do, at very different abundances.
    double congener_share{0.0};
    std::size_t congener_min{2};
    std::size_t congener_max{5};

    // Cross-sample strain sharing: consumed by MetagenomeSimulator before the per-sample loop.
    std::vector<StrainSharingSpec> strain_sharing;

    // Per-sample fields populated by MetagenomeSimulator from strain_sharing pre-assignment:
    // Maps species key -> specific GenomeRecord list to use (bypasses normal strain picking).
    std::unordered_map<std::string, std::vector<GenomeRecord>> forced_strains;
    // Maps species key -> minimum relative abundance floor derived from min_vcov.
    std::unordered_map<std::string, double> species_min_abundance;
};

// The Poisson-lognormal sigma of sample `index` (0-based) of a design: the index-th of pln_sigmas, cyclically, or
// pln_sigma when none are given.
inline double SigmaForSample(ProfileDesignOptions const& options, std::size_t index) {
    if (options.pln_sigmas.empty()) return options.pln_sigma;
    return options.pln_sigmas[index % options.pln_sigmas.size()];
}

// The read pairs of sample `index` (0-based) of a design: the index-th of total_read_pairs_per_sample, cyclically, or
// total_read_pairs when none are given.
inline std::uint64_t ReadPairsForSample(ProfileDesignOptions const& options, std::size_t index) {
    if (options.total_read_pairs_per_sample.empty()) return options.total_read_pairs;
    return options.total_read_pairs_per_sample[index % options.total_read_pairs_per_sample.size()];
}

// How a sample's read files are written: BGZF (ISA-L at level 1, _R1.fq.gz) or zstd (one frame at level 3, _R1.fq.zst;
// smaller, and read by protal alike).
enum class ReadsCompression { Bgzf, Zstd };

// How the Illumina paired-end reads are made (IlluminaSimulator.h).
struct IlluminaOptions {
    int read_length{150};
    int fragment_mean{350};
    int fragment_stdev{50};
    std::string sequencer{"HS25"};       // an IlluminaProfile: HS20, HS25, HSXt, NovaSeq, MSv3
    std::optional<double> mean_quality;  // each read's mean base quality (Phred), else the profile's
    int threads{1};
    std::filesystem::path host_folder;   // a host genome (scenarios.prepare_host), for host_pairs
    std::vector<std::uint64_t> host_pairs;  // each sample's host read pairs, in turn (after its community's)
    bool first_reads_only{false};        // only the _R1 files (the same reads as with both)
};

struct SampleOutput {
    std::string sample_name;
    std::filesystem::path read1_path;
    std::filesystem::path read2_path;
    std::vector<GenomeAssignment> assignments;
};

}  // namespace protal::sim
