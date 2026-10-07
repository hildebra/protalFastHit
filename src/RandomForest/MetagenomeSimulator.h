// SPDX-License-Identifier: GPL-2.0-only
#pragma once

#include <filesystem>
#include <optional>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#include "../Utilities/Benchmark.h"
#include "IlluminaSimulator.h"
#include "CommunityProfileDesigner.h"
#include "MetagenomeTypes.h"

namespace protal::sim {

std::vector<GenomeRecord> read_genome_table(const std::filesystem::path& tsv_path);
std::vector<double> parse_strain_probabilities(const std::string& text);
std::vector<std::string> parse_species_list(const std::string& text);
std::unordered_map<std::string, std::size_t> parse_genus_selection(const std::string& text);
std::unordered_map<std::string, std::size_t> parse_taxon_selection(const std::string& text);
// --congener_groups SHARE:MIN-MAX (e.g. "0.25:2-5") into options' congener_share, congener_min and congener_max; ""
// leaves them (no groups). Throws on a share outside 0-1, a MIN below 2 or above MAX, or another form.
void parse_congener_groups(const std::string& text, ProfileDesignOptions& options);
void write_combined_manifest(const std::vector<SampleOutput>& samples, const std::filesystem::path& manifest_path);
void write_sample_manifest(const SampleOutput& sample, const std::filesystem::path& manifest_path);
void write_abundance_matrix(const std::vector<SampleOutput>& samples, const std::filesystem::path& matrix_path);

// Reads a manifest written by write_combined_manifest or write_sample_manifest back
// into a design: which genome at which read depth in which sample, plus the reads' seed (art_seed)
// and fasta path when the manifest carries them. Sample and row order are preserved.
std::vector<SampleOutput> read_manifest(const std::filesystem::path& manifest_path);

// Fills in fasta paths for manifests written before the fasta_path column existed,
// looking each genome up in a genome table by name. Throws naming the genomes it
// cannot resolve.
void resolve_manifest_fasta_paths(std::vector<SampleOutput>& samples, const std::vector<GenomeRecord>& genomes);
class MetagenomeSimulator {
public:
    MetagenomeSimulator(
        std::vector<GenomeRecord> genomes,
        IlluminaOptions illumina_options = {},
        std::uint64_t seed = std::random_device{}());

    std::vector<SampleOutput> simulate_samples(
        const ProfileDesignOptions& profile_options,
        std::size_t sample_count,
        const std::string& sample_prefix,
        const std::filesystem::path& output_dir,
        bool skip_reads = false,
        bool keep_tmp = false);

    // Re-simulates a design read back from a manifest, bypassing the profile designer
    // entirely. Rows carrying an art_seed reproduce their reads exactly; rows without
    // one keep the recorded composition and get fresh read realizations.
    std::vector<SampleOutput> replay_samples(
        std::vector<SampleOutput> design,
        const std::filesystem::path& output_dir,
        bool skip_reads = false,
        bool keep_tmp = false);

    // How the samples' read files are written (default BGZF).
    void set_reads_compression(ReadsCompression compression) { reads_compression_ = compression; }

private:
    std::vector<GenomeRecord> genomes_;
    IlluminaOptions illumina_;
    std::uint64_t seed_;
    CommunityProfileDesigner designer_;
    std::mt19937_64 rng_;
    ReadsCompression reads_compression_ = ReadsCompression::Bgzf;

    // A sample's design (its assignments), prepared (prepare_sample); its reads are written later.
    SampleOutput simulate_single(
        const ProfileDesignOptions& profile_options,
        const std::string& sample_name,
        const std::unordered_map<std::string, std::uint64_t>& genome_lengths,
        std::uint64_t paired_read_length);

    // Draws each assignment's seed of its reads (art_seed in the manifest, a name from the ART days) from rng_, as many
    // as a run draws in the same order whether or not it writes reads, and its coverage: everything of a sample that
    // depends on rng_, so that the reads can be made on threads.
    void prepare_sample(SampleOutput& sample, std::uint64_t paired_read_length);

    // The samples' _R1 and _R2 files (only _R1 with first_reads_only): each assignment's pairs from its seed, then the
    // sample's host pairs (IlluminaSimulator: SimulatePairs), compressed as they are made (BGZF, .fq.gz, or zstd,
    // .fq.zst: reads_compression_), the same bytes for any number of threads. Files that are named pipes are written
    // in place, a sample at a time (a reader such as protal takes them as they are made). With skip_reads, empty
    // placeholder files (a named pipe is left alone).
    void write_all_reads(std::vector<SampleOutput>& samples, const std::filesystem::path& output_dir, bool skip_reads) const;
};

}  // namespace protal::sim
