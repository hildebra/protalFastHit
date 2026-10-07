// SPDX-License-Identifier: GPL-2.0-only
#pragma once

// Illumina paired-end reads made in process (simulate_metagenomes; ART is no longer used): fragments drawn from a
// sample's genomes (and a host's), each made into a pair of reads by a model of Illumina sequencing built from
// published characteristics of the instruments (IlluminaProfile; numbers and sources in
// docs/claude/2026-10-07-illumina-model/README.md):
//
// - a fragment: its length normal (mean, SD; at least a read long, at most the longest contig), its place uniform over
//   the genome's contigs where it fits, either strand; read 1 is its first bases, read 2 its reverse complement's;
// - a base's hidden quality: in the high state, the read's mean by cycle (calibrated so that the qualities written
//   average the instrument's published curve, per read) plus a run's, a cluster's (shared by the pair) and the read's
//   offsets and noise correlated along the read; a read falls into a low state (qualities ~Q10-25) now and then, more
//   often towards its end, for a few cycles, and a few reads collapse to Q2 for the rest of the read;
// - errors: a base is miscalled with its hidden quality's probability (times a calibration factor, more after GG), as
//   a substitution by the instrument's spectrum; Q2-tail bases with a few percent; an N (Q2) now and then (more at
//   the first cycle) and wherever the genome has one; insertions and deletions at their own rates, more in
//   homopolymers of 5 or more; a read is always its full length;
// - the quality written is the hidden one, binned as the instrument bins it.
//
// The rare events of a cycle (the low state's entry, an indel, an N) are drawn as geometric gaps between them rather
// than with a uniform number per cycle, the rates varying by cycle by thinning (a candidate at the highest rate kept at
// the cycle's); normal numbers come from a ziggurat (LongRng::Gaussian). Read 2 draws from a stream of its own.
//
// Reads are named <contig>-<n>/1 and /2 (trace_relatives.py and error_reads.py take the contig from the name), a
// host's h_<n>/1 and /2, n = 1, 2, ... over the sample.

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "LongReadSimulator.h"
#include "ReadPipeline.h"

namespace protal::sim {

struct IlluminaProfile {
    std::string name, description;
    // The published mean quality (written, binned) by cycle of read 1 and read 2: (cycle, Phred) knots, linear between
    // and on at the last two's slope beyond.
    std::vector<std::pair<double, double>> mean[2];
    double run_sd = 1.0;       // SD of a sample's (run's) quality offset
    double cluster_sd = 1.5;   // SD of a pair's offset (both reads)
    double read_sd = 1.0;      // SD of a read's own offset
    double noise_sd = 2.0, noise_sd_end = 2.0, noise_rho = 0.8;  // AR(1) noise along a read (SD at its start, end)
    // The low state: entered at cycle c with low_rate (read 2: low_rate_r2) x exp(low_growth x (c / length - 1)), left
    // with low_exit per cycle; its qualities normal (low_mean, low_sd).
    double low_rate = 0.004, low_rate_r2 = 0.006, low_growth = 2.0, low_exit = 0.25, low_mean = 15.0, low_sd = 5.0;
    double collapse = 0.02;        // share of reads that fall to Q2 from a cycle in their last 30% to their end
    double collapse_error = 0.05;  // the error probability of those Q2 bases
    std::vector<std::pair<int, int>> bins;  // (lowest quality of a bin, the value written); empty: unbinned
    int q_min = 2, q_max = 41;
    double calibration = 1.0;      // a base's error probability is its hidden quality's times this
    double gg_factor = 1.0;        // ... and this after GG (read direction)
    double insertion = 5e-6, deletion = 5e-6, homopolymer_indels = 10.0;  // per base; x this in runs of 5 or more
    double n_rate = 1e-4, n_first = 1e-3;  // an N at Q2, per base (at the first cycle)
    double r2_insert_penalty = 0.005;      // read 2's qualities lower by this per base of fragment above 500
    double substitution[4][4] = {{0, 1, 1, 1}, {1, 0, 1, 1}, {1, 1, 0, 1}, {1, 1, 1, 0}};  // weights, from ACGT to ACGT

    // HS20, HS25, HSXt, NovaSeq, MSv3; throws std::invalid_argument for another name.
    static IlluminaProfile Named(std::string const& name);
    static std::vector<std::string> Names();
};

struct IlluminaSetup {
    IlluminaProfile profile = IlluminaProfile::Named("NovaSeq");
    int read_length = 150;
    double fragment_mean = 350, fragment_sd = 50;
    std::optional<double> mean_quality;  // the qualities written average this (Phred) per read, else the profile's
};

// What a read's errors were (for checking a profile).
struct IlluminaEvents {
    std::uint64_t substitutions = 0, insertions = 0, deletions = 0, ns = 0;
};

class IlluminaModel {
public:
    explicit IlluminaModel(IlluminaSetup setup);

    IlluminaSetup const& Setup() const { return m_setup; }
    // A fragment's length (normal; at least a read long, at most `longest` if it is that long).
    std::uint32_t FragmentLength(LongRng& rng, std::uint64_t longest) const;
    // The two reads of a fragment (its strand): read 1 from its start, read 2 from its reverse complement's; `run`: the
    // sample's quality offset (RunOffset). Read 2 draws from a stream of its own (rng2), so that read 1s do not depend
    // on whether read 2 is made: with r2 null it is not (first reads only).
    void Pair(std::string_view fragment, double run, LongRng& rng, LongRng& rng2, MadeRead& r1, MadeRead* r2,
              std::uint64_t* errors = nullptr) const;
    double RunOffset(LongRng& rng) const { return m_setup.profile.run_sd * rng.Gaussian(); }
    // One read (0: read 1, 1: read 2) of a template, with a quality offset (the run's, the cluster's and the read's).
    void Read(std::string_view templ, int read, double offset, LongRng& rng, MadeRead& out,
              std::uint64_t* errors = nullptr, IlluminaEvents* events = nullptr) const;
    // The high state's hidden mean quality by cycle (0-based), as calibrated.
    std::vector<double> const& HiddenMean(int read) const { return m_mu[read]; }
    // The published (target) mean written quality by cycle, with the shift to mean_quality.
    std::vector<double> const& TargetMean(int read) const { return m_target[read]; }

private:
    void Calibrate();

    IlluminaSetup m_setup;
    std::vector<double> m_mu[2], m_target[2], m_sd[2], m_low_entry[2];
    // The rare events of a cycle, drawn as geometric gaps between candidates at a bound's rate (Read): the low state's
    // entry (its highest rate by cycle), an indel (in homopolymers' rate), an N; and log(1 - rate) of each.
    double m_entry_max[2] = {0, 0}, m_entry_log[2] = {0, 0}, m_indel_max = 0, m_indel_log = 0, m_n_log = 0;
    double m_error[94];
    char m_reported[94];
    double m_sub_cumulative[4][3];
    char m_sub_base[4][3];
};

// The paired-end reads of a run's samples.
struct PairedSample {
    std::string name;
    std::filesystem::path r1, r2;  // .fq.zst or .fq.gz; r2 empty: first reads only (the same as with both)
    struct Genome {
        std::string name;
        std::filesystem::path fasta;
        std::uint64_t pairs = 0, seed = 0;
    };
    std::vector<Genome> genomes;
    std::uint64_t host_pairs = 0, host_seed = 0;  // after the genomes', from PairedOptions::host
    std::uint64_t run_seed = 0;  // the sample's quality offset (IlluminaModel::RunOffset)
};

struct PairedOptions {
    IlluminaSetup setup;
    std::filesystem::path host;  // a folder of scenarios.prepare_host, for samples with host pairs
    int threads = 1;
    std::filesystem::path genome_store;  // genomes from a genome store (GenomeStore.h), or none
    bool plain_pipes = false;            // outputs that are named pipes get plain FASTQ
};

// The samples' reads (pipeline::Run): each genome's pairs, then the host's, in work items of ~600 kB of FASTQ per
// read file, each with a random stream of the genome's seed and its place. -> each sample's totals (reads: pairs).
std::vector<pipeline::Totals> SimulatePairs(std::vector<PairedSample> const& samples, PairedOptions const& options);

// Statistics of a profile's reads, from `pairs` pairs of a random genome: per read, the written quality's mean by
// cycle, the share of bases at Q30 or more, of each written quality, the errors per base (substitutions, insertions,
// deletions, Ns). For checking a profile against its sources (simulate_metagenomes --illumina_report).
struct IlluminaReport {
    std::vector<double> mean_by_cycle[2];
    double q30[2] = {0, 0}, substitutions[2] = {0, 0}, insertions[2] = {0, 0}, deletions[2] = {0, 0}, ns[2] = {0, 0};
    std::vector<double> quality_share[2];  // by written quality, 0-93
};
IlluminaReport ReportProfile(IlluminaSetup const& setup, std::uint64_t pairs, std::uint64_t seed = 1);

}  // namespace protal::sim
