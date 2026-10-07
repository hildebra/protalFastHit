// SPDX-License-Identifier: GPL-2.0-only
#pragma once

// Long reads (PacBio HiFi, Nanopore) and Ultima Genomics single-end reads of a sample's community, made in process
// (simulate_metagenomes --long_samples): templates drawn from the genomes as collect_training_data.py drew them, each
// made into one read by a read model, written compressed as they are made. No template files, no pbsim3 run, no
// Python.
//
// - Templates (as collect_training_data.draw_templates): a read's genome by its weight (relative abundance times
//   length), its length from a gamma distribution of the setup's mean and SD between 100 bases and 1 Mb, its start
//   uniform over the genome's contigs of 100 bases or more (a read ends where its contig does), either strand; drawn
//   in rounds until the sample's bases are reached (the cuts at contigs' ends leave a few bases for another round).
//   A host genome (a folder prepared by scenarios.prepare_host) gives its templates by memory map at random places.
// - Read models: hifi_reads.py's HiFi model and its flow model (Ultima), and pbsim3's qshmm model in its template
//   mode (simulate_by_qshmm_templ of pbsim3 3.0.x, GPL-2.0, Yukiteru Ono), which reads pbsim3's QSHMM-*.model files.
// - Reads are named g<i>x_<n>: i the genome's place in the sample's list, n = 1, 2, ... in the file's order
//   (error_reads.py reads the genome from the name).
//
// The reads do not depend on the threads: a sample's rounds are planned with its own random stream, and each work
// item (some reads of one genome in one round) has a stream of its own; the items are written in their order.

#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "MetagenomeTypes.h"

namespace protal::sim {

// xoshiro256** seeded through splitmix64: small, fast, and the same everywhere.
class LongRng {
public:
    using result_type = std::uint64_t;
    explicit LongRng(std::uint64_t seed);
    static constexpr result_type min() { return 0; }
    static constexpr result_type max() { return ~result_type{0}; }
    result_type operator()();
    double Uniform();                      // [0, 1)
    std::uint64_t Below(std::uint64_t n);  // [0, n), n > 0
    double Normal();                       // N(0, 1)

private:
    std::uint64_t m_s[4];
    bool m_has_spare = false;
    double m_spare = 0.0;
};

// A seed made of two numbers (splitmix64 of their mix), for the random streams of work items.
std::uint64_t MixSeed(std::uint64_t a, std::uint64_t b);

// A read setup as collect_training_data.parse_long_setup reads it: hifi:LENGTH_MEAN:LENGTH_SD:Q_SD,
// ultima:LENGTH_MEAN:LENGTH_SD:Q_MEAN:Q_SD, or qshmm:MODEL:LENGTH_MEAN:LENGTH_SD:ACCURACY_MEAN[:SUB/INS/DEL].
struct LongReadSetup {
    enum class Method { Hifi, Ultima, Qshmm };
    Method method = Method::Hifi;
    std::string model;  // qshmm: the model's name in the setup (the file is given apart: --long_model)
    double length_mean = 15000, length_sd = 3000;
    double q_sd = 3.0, q_mean = 0.0;           // hifi, ultima
    double accuracy = 0.85;                    // qshmm: pbsim3's --accuracy-mean
    long sub_ratio = 6, ins_ratio = 55, del_ratio = 39;  // qshmm: --difference-ratio, pbsim3's defaults

    static LongReadSetup Parse(std::string const& text);  // throws std::invalid_argument
};

struct MadeRead {
    std::string seq, qual;  // bases and Phred+33 qualities
};

// hifi_reads.py's models, read by read (MODEL and FLOW_MODEL there).
namespace hifi {
    inline constexpr double kQMin = 10.0, kQMax = 60.0, kPMax = 0.5, kConfidenceSigma = 1.0;
    inline constexpr int kHomopolymer = 3, kFlowHomopolymer = 2;
    inline constexpr double kSubstitution = 0.30, kInsertion = 0.35;  // the rest deletions, outside homopolymers

    // A read's mean quality (Phred) by its length: Q50 to 5 kb, Q30 at 25 kb, Q20 at 50 kb, linear between and on at
    // the last slope beyond (QUALITY_BY_LENGTH).
    double LengthQuality(double length);

    // A read of `templ` (upper case): with setup.method Ultima the flow model (a mean base quality of q_mean, SD
    // q_sd), else HiFi (a read quality by length, SD q_sd). Counts its errors into `errors` if not null.
    void Mutate(std::string_view templ, LongRng& rng, LongReadSetup const& setup, MadeRead& out,
                std::uint64_t* errors = nullptr);
}

// pbsim3's qshmm model in template mode: a read's accuracy level drawn in [0.75, 1.05] x the mean accuracy (weights
// exp(0.22 x level), at most 100), its qualities by the level's hidden Markov model, and each base's substitution,
// insertion or following deletions by its quality and the difference ratio. Tables quantized as pbsim3 quantizes them
// (to 1/100 for the HMM, 1/1000 for the uniform ones). One change: only levels with an HMM are drawn (pbsim3 gives
// the others a uniform quality: level 100 is Q93 throughout and error-free, ~20% of its reads at a mean of 0.97).
class QshmmModel {
public:
    QshmmModel(std::filesystem::path const& file, double accuracy_mean, long sub_ratio, long ins_ratio,
               long del_ratio);

    void Mutate(std::string_view templ, LongRng& rng, MadeRead& out, std::uint64_t* errors = nullptr) const;

    int DrawAccuracy(LongRng& rng) const;  // the accuracy level (percent) of a read
    int AccuracyMin() const { return m_accuracy_min; }
    int AccuracyMax() const { return m_accuracy_max; }
    bool HasHmm(int accuracy) const { return m_exists[accuracy]; }

private:
    static constexpr int kAccuracyMax = 100, kStates = 50, kQualities = 94;
    template <typename T> using PerAccuracy = std::vector<T>;

    int m_accuracy_min = 0, m_accuracy_max = 0;
    long m_accuracy_draw = 0;                 // accuracy_rand_value
    std::vector<int> m_prob2accuracy;         // [1 .. 100000]
    std::vector<char> m_exists;               // [0 .. 100]
    std::vector<long> m_init_draw;            // qc_rand_value_init[accuracy]
    std::vector<int> m_init2state;            // [accuracy][1 .. 100]
    std::vector<long> m_emis_draw, m_tran_draw;  // [accuracy][state]
    std::vector<int> m_emis2qc, m_tran2state;    // [accuracy][state][1 .. 100]
    std::vector<long> m_freq_draw;            // [accuracy]
    std::vector<int> m_freq2qc;               // [accuracy][1 .. 1000]
    long m_sub_thre[kQualities] = {}, m_ins_thre[kQualities] = {}, m_del_thre[kQualities] = {};
};

// A long read's length: gamma of the setup's mean and SD, redrawn until within [100, 1,000,000]; the mean (clamped)
// if the SD is 0 (collect_training_data.read_length).
std::uint32_t DrawReadLength(LongRng& rng, double mean, double sd);

// Reverse complement as the collector's COMPLEMENT: ACGTN complemented, other letters kept.
void ReverseComplement(std::string& seq);

struct LongGenome {
    std::string name;
    std::filesystem::path fasta;  // the genome, or the host's folder for a host
    double weight = 0.0;          // relative abundance x length (the host: its share of all)
    bool host = false;
};

struct LongSample {
    std::string name;
    std::filesystem::path out;   // .fq.zst (zstd) or .fq.gz (BGZF), written as out + ".partial", then renamed
    std::uint64_t bases = 0;     // template bases to reach
    std::uint64_t seed = 1;
    std::vector<LongGenome> genomes;  // a read of genome i is named g<i>x_<n>
};

struct LongSampleResult {
    std::string name;
    std::uint64_t reads = 0, template_bases = 0, read_bases = 0, errors = 0;
    std::uint32_t rounds = 0;
};

struct LongReadOptions {
    LongReadSetup setup;
    std::filesystem::path model;  // qshmm: the model file
    int threads = 1;
};

// The samples' reads, on `threads` threads; the results in the samples' order. Throws on any failure (a genome that
// cannot be read, a file that cannot be written); the samples written by then stay, the others leave no file.
std::vector<LongSampleResult> SimulateLongReads(std::vector<LongSample> const& samples, LongReadOptions const& options);

// One read of each template of a FASTA (plain, gzip or zstd), named after it, by the setup's model, into `out` (.fq.zst
// or .fq.gz), template i with its own random stream of `seed` and i: what hifi_reads.py and pbsim3 --strategy templ do,
// for comparing the models on the same templates. One thread.
LongSampleResult MutateTemplates(std::filesystem::path const& templates, std::filesystem::path const& out,
                                 LongReadOptions const& options, std::uint64_t seed);

// The samples and genomes of a --long_samples / --long_genomes pair of TSVs (LongReadSimulator.cpp says their columns).
std::vector<LongSample> ReadLongSamples(std::filesystem::path const& samples_tsv,
                                        std::filesystem::path const& genomes_tsv);

}  // namespace protal::sim
