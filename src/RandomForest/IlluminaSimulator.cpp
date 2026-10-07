// SPDX-License-Identifier: GPL-2.0-only
#include "IlluminaSimulator.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>

namespace fs = std::filesystem;

namespace protal::sim {

// ---- the instruments ---------------------------------------------------------------------------------------------
// The numbers and their sources: docs/claude/2026-10-07-illumina-model/README.md. The mean curves are InSilicoSeq's
// (Gourle et al. 2019) HiSeq 2x126, NovaSeq 2x151 and MiSeq 2x301 models as its documentation plots them; the bins
// Illumina's technical notes; the error rates Schirmer et al. 2015/2016, Stoler & Nekrutenko 2021 and Illumina's
// NovaSeq notes; the low-state and offset parameters are fitted so that the written qualities and errors match them
// (simulate_metagenomes --illumina_report).

namespace {

// From A, C, G, T to A, C, G, T. Four-colour chemistry (HiSeq, MiSeq): A<->C and G<->T dominate (Schirmer 2015: A>C
// 66% of A's errors, C>A 58% of C's); two-colour (NovaSeq): transitions, and G (no signal) late.
constexpr double kFourColour[4][4] = {{0, 0.66, 0.20, 0.14}, {0.58, 0, 0.12, 0.30}, {0.30, 0.12, 0, 0.58}, {0.14, 0.20, 0.66, 0}};
constexpr double kTwoColour[4][4] = {{0, 0.25, 0.55, 0.20}, {0.20, 0, 0.25, 0.55}, {0.55, 0.25, 0, 0.20}, {0.20, 0.55, 0.25, 0}};

void Spectrum(IlluminaProfile& p, double const (&s)[4][4]) {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) p.substitution[i][j] = s[i][j];
    }
}

std::vector<IlluminaProfile> const& Profiles() {
    static std::vector<IlluminaProfile> const profiles = [] {
        std::vector<IlluminaProfile> all;
        // HiSeq 2000/2500: InSilicoSeq's HiSeq 2x126 curve; R2 0.5 below R1.
        IlluminaProfile hiseq;
        hiseq.mean[0] = {{1, 32.5}, {5, 32.5}, {8, 37.0}, {70, 36.0}, {110, 35.0}, {126, 30.0}};
        hiseq.mean[1] = {{1, 32.0}, {5, 32.0}, {8, 36.5}, {70, 35.5}, {110, 34.5}, {126, 29.5}};
        hiseq.q_max = 41;
        hiseq.gg_factor = 6.0;
        hiseq.insertion = 2.8e-6;
        hiseq.deletion = 5.1e-6;
        hiseq.noise_sd = 2.0;
        hiseq.noise_sd_end = 3.5;
        Spectrum(hiseq, kFourColour);

        IlluminaProfile p = hiseq;
        p.name = "HS20";
        p.description = "HiSeq 2000, 2x100 bp: qualities unbinned (2-41), a tenth of the reads with a Q2 tail";
        p.collapse = 0.10;
        p.low_rate = 0.020;
        p.low_rate_r2 = 0.032;
        all.push_back(p);

        p = hiseq;
        p.name = "HS25";
        p.description = "HiSeq 2500, 2x125 bp: qualities in Illumina's 8 bins";
        p.bins = {{3, 6}, {10, 15}, {20, 22}, {25, 27}, {30, 33}, {35, 37}, {40, 40}};
        p.collapse = 0.05;
        p.low_rate = 0.020;
        p.low_rate_r2 = 0.032;
        all.push_back(p);

        // HiSeq X Ten / 4000 (four-colour, patterned): InSilicoSeq's NovaSeq curve, the bins seen in its files.
        p = IlluminaProfile{};
        p.name = "HSXt";
        p.description = "HiSeq X Ten / HiSeq 4000, 2x150 bp: qualities in 7 bins (2, 12, 22, 27, 32, 37, 41)";
        p.mean[0] = {{1, 36.0}, {60, 37.0}, {110, 36.5}, {150, 35.5}};
        p.mean[1] = {{1, 36.5}, {60, 36.0}, {110, 34.5}, {151, 31.0}};
        p.bins = {{3, 12}, {17, 22}, {25, 27}, {30, 32}, {35, 37}, {40, 41}};
        p.q_max = 41;
        p.gg_factor = 1.5;
        p.collapse = 0.02;
        p.low_rate = 0.022;
        p.low_rate_r2 = 0.045;
        p.noise_sd_end = 3.0;
        Spectrum(p, kFourColour);
        all.push_back(p);

        p.name = "NovaSeq";
        p.description = "NovaSeq 6000, 2x150 bp: qualities in 4 bins (2, 12, 23, 37), two-colour chemistry";
        p.bins = {{10, 12}, {20, 23}, {30, 37}};
        p.q_max = 40;
        p.low_rate = 0.016;
        p.low_rate_r2 = 0.040;
        Spectrum(p, kTwoColour);
        all.push_back(p);

        // MiSeq v3: InSilicoSeq's MiSeq 2x301 curve; read 2 falls off at its end; high qualities overstated.
        p = IlluminaProfile{};
        p.name = "MSv3";
        p.description = "MiSeq v3, 2x250 or 2x300 bp: qualities unbinned (2-41), read 2 low at its end";
        p.mean[0] = {{1, 33.5}, {9, 33.5}, {10, 37.5}, {100, 37.5}, {160, 35.5}, {220, 32.0}, {264, 28.0}, {300, 23.0}};
        p.mean[1] = {{1, 33.0}, {60, 37.0}, {160, 33.0}, {220, 27.0}, {264, 21.0}, {300, 17.0}};
        p.q_max = 41;
        p.calibration = 1.5;
        p.gg_factor = 6.0;
        p.insertion = 4e-5;
        p.deletion = 2e-5;
        p.noise_sd = 2.0;
        p.noise_sd_end = 7.0;
        p.collapse = 0.04;
        p.low_rate = 0.012;
        p.low_rate_r2 = 0.022;
        Spectrum(p, kFourColour);
        all.push_back(p);
        return all;
    }();
    return profiles;
}

}  // namespace

IlluminaProfile IlluminaProfile::Named(std::string const& name) {
    for (auto const& p : Profiles()) {
        if (p.name == name) return p;
    }
    std::string names;
    for (auto const& n : Names()) names += (names.empty() ? "" : ", ") + n;
    throw std::invalid_argument("Illumina profile " + name + ": expected one of " + names);
}

std::vector<std::string> IlluminaProfile::Names() {
    std::vector<std::string> names;
    for (auto const& p : Profiles()) names.push_back(p.name);
    return names;
}

// ---- the model ---------------------------------------------------------------------------------------------------

namespace {

double Interpolate(std::vector<std::pair<double, double>> const& knots, double x) {
    if (knots.empty()) return 30.0;
    if (knots.size() == 1 || x <= knots.front().first) return knots.front().second;
    std::size_t k = 1;
    while (k + 1 < knots.size() && x > knots[k].first) ++k;  // the segment [k-1, k], the last one beyond its end
    auto const& [x0, y0] = knots[k - 1];
    auto const& [x1, y1] = knots[k];
    return x1 == x0 ? y1 : y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}

inline int Code(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return 4;
    }
}

}  // namespace

IlluminaModel::IlluminaModel(IlluminaSetup setup) : m_setup(std::move(setup)) {
    auto const& p = m_setup.profile;
    int const length = m_setup.read_length;
    if (length < 1) throw std::invalid_argument("the read length must be positive");
    if (!(m_setup.fragment_mean > 0) || m_setup.fragment_sd < 0) {
        throw std::invalid_argument("the fragment length mean must be positive and its SD not negative");
    }
    for (int r = 0; r < 2; ++r) {
        m_target[r].resize(length);
        m_sd[r].resize(length);
        m_low_entry[r].resize(length);
        for (int c = 0; c < length; ++c) {
            m_target[r][c] = Interpolate(p.mean[r], c + 1);
            double const x = length > 1 ? static_cast<double>(c) / (length - 1) : 0.0;
            m_sd[r][c] = p.noise_sd + (p.noise_sd_end - p.noise_sd) * x * x;
            double const rate = r == 0 ? p.low_rate : p.low_rate_r2;
            m_low_entry[r][c] = std::min(0.5, rate * std::exp(p.low_growth * (static_cast<double>(c + 1) / length - 1)));
        }
        if (m_setup.mean_quality) {
            double sum = 0;
            for (double q : m_target[r]) sum += q;
            double const shift = *m_setup.mean_quality - sum / length;
            for (double& q : m_target[r]) q += shift;
        }
        m_mu[r] = m_target[r];
    }
    for (int q = 0; q < 94; ++q) {
        m_error[q] = std::min(0.75, p.calibration * std::pow(10.0, -q / 10.0));
        int written = q;
        if (q > 2) {
            for (auto const& [low, value] : p.bins) {
                if (q >= low) written = value;
            }
            if (!p.bins.empty() && q < p.bins.front().first) written = 2;  // below the lowest bin: Q2, as Illumina writes it
        }
        m_reported[q] = static_cast<char>(33 + std::clamp(written, 0, 93));
    }
    for (int from = 0; from < 4; ++from) {
        double total = 0;
        int k = 0;
        for (int to = 0; to < 4; ++to) {
            if (to == from) continue;
            total += std::max(0.0, p.substitution[from][to]);
            m_sub_cumulative[from][k] = total;
            m_sub_base[from][k] = "ACGT"[to];
            ++k;
        }
        if (!(total > 0)) throw std::invalid_argument("profile " + p.name + ": a base without substitutions");
        for (auto& c : m_sub_cumulative[from]) c /= total;
    }
    for (int r = 0; r < 2; ++r) {
        m_entry_max[r] = *std::max_element(m_low_entry[r].begin(), m_low_entry[r].end());
        m_entry_log[r] = std::log1p(-std::min(m_entry_max[r], 0.999999));
    }
    m_indel_max = (p.insertion + p.deletion) * std::max(1.0, p.homopolymer_indels);
    m_indel_log = std::log1p(-std::min(m_indel_max, 0.999999));
    m_n_log = std::log1p(-std::min(p.n_rate, 0.999999));
    Calibrate();
}

// The hidden mean by cycle that makes the written qualities average the target curve: a few rounds of reads of a
// random template (a fixed seed: the same for every run), each moving the mean by the difference at each cycle, at
// most to 4 above the target (binned instruments write a few values only).
void IlluminaModel::Calibrate() {
    int const length = m_setup.read_length;
    LongRng rng(0x1111a);
    std::string templ(static_cast<std::size_t>(length) + 64, 'A');
    MadeRead read;
    for (int round = 0; round < 6; ++round) {
        for (int r = 0; r < 2; ++r) {
            std::vector<double> sum(length, 0.0);
            int const reads = 2000;
            for (int i = 0; i < reads; ++i) {
                for (auto& c : templ) c = "ACGT"[rng.Below(4)];
                double const offset = m_setup.profile.cluster_sd * rng.Gaussian() + m_setup.profile.read_sd * rng.Gaussian();
                Read(templ, r, offset, rng, read, nullptr);
                for (int c = 0; c < length; ++c) sum[c] += read.qual[c] - 33;
            }
            for (int c = 0; c < length; ++c) {
                double const written = sum[c] / reads;
                m_mu[r][c] = std::min(m_mu[r][c] + std::clamp(m_target[r][c] - written, -5.0, 5.0), m_target[r][c] + 4.0);
            }
        }
    }
}

std::uint32_t IlluminaModel::FragmentLength(LongRng& rng, std::uint64_t longest) const {
    auto const length = static_cast<std::uint64_t>(m_setup.read_length);
    std::uint64_t const most = std::max(length, longest);
    for (int attempt = 0; attempt < 100; ++attempt) {
        double const x = std::round(m_setup.fragment_mean + m_setup.fragment_sd * rng.Gaussian());
        if (x >= static_cast<double>(length) && x <= static_cast<double>(most)) return static_cast<std::uint32_t>(x);
    }
    return static_cast<std::uint32_t>(std::clamp<std::uint64_t>(static_cast<std::uint64_t>(m_setup.fragment_mean), length, most));
}

namespace {

constexpr int kNever = 1 << 30;

// The failures before the first success of Bernoulli trials of probability p (log_q = log(1 - p)): the cycles to skip
// to the next candidate of a rare event, one draw for all of them.
inline int Gap(LongRng& rng, double p, double log_q) {
    if (!(p > 0)) return kNever;
    if (p >= 1) return 0;
    double const g = std::log(1.0 - rng.Uniform()) / log_q;  // log of (0, 1]: >= 0
    return g >= kNever ? kNever : static_cast<int>(g);
}

}  // namespace

void IlluminaModel::Read(std::string_view templ, int read, double offset, LongRng& rng, MadeRead& out,
                         std::uint64_t* errors, IlluminaEvents* events_out) const {
    IlluminaEvents unused;
    IlluminaEvents& events = events_out ? *events_out : unused;
    auto const& p = m_setup.profile;
    int const length = m_setup.read_length;
    out.seq.assign(length, 'N');
    out.qual.assign(length, static_cast<char>(33 + 2));
    int collapse = length;  // the cycle the read falls to Q2 from
    if (rng.Uniform() < p.collapse) {
        int const from = length - std::max(1, (3 * length) / 10);
        collapse = from + static_cast<int>(rng.Below(static_cast<std::uint64_t>(length - from)));
    }
    auto const& mu = m_mu[read];
    auto const& sd = m_sd[read];
    auto const& entry = m_low_entry[read];
    double const rho = p.noise_rho, innovation = std::sqrt(std::max(0.0, 1 - rho * rho));
    double const entry_max = m_entry_max[read], entry_log = m_entry_log[read];
    double noise = rng.Gaussian();  // in SDs
    bool low = false;
    // The cycles of the next candidates of the rare events: the low state's entry (at the highest entry rate, kept at
    // the cycle's), an indel (at the homopolymers' rate, kept at the place's), an N (after the first cycle's own).
    int next_entry = Gap(rng, entry_max, entry_log);
    int next_indel = Gap(rng, m_indel_max, m_indel_log);
    bool const n_first = rng.Uniform() < p.n_first;
    int next_n = 1 + Gap(rng, p.n_rate, m_n_log);
    std::size_t at = 0;  // in the template
    std::uint64_t made = 0;
    for (int c = 0; c < length; ++c) {
        if (c > 0) noise = rho * noise + innovation * rng.Gaussian();
        if (low) {  // the low state ends with low_exit a cycle; the next entry can come from the next cycle on
            if (rng.Uniform() < p.low_exit) {
                low = false;
                next_entry = c + 1 + Gap(rng, entry_max, entry_log);
            }
        } else if (c == next_entry) {
            if (rng.Uniform() * entry_max < entry[c]) {
                low = true;
            } else {
                next_entry = c + 1 + Gap(rng, entry_max, entry_log);
            }
        }
        double const hq = low ? p.low_mean + p.low_sd * rng.Gaussian() : mu[c] + offset + sd[c] * noise;
        int q = std::clamp(static_cast<int>(std::lround(hq)), p.q_min, p.q_max);
        bool const tail = c >= collapse;
        if (tail) q = 2;
        bool const n_call = c == 0 ? n_first : c == next_n;  // no call at this cycle (if it reads a template base)
        if (c > 0 && c == next_n) next_n = c + 1 + Gap(rng, p.n_rate, m_n_log);
        // an insertion (a base not in the template) or a deletion (a template base skipped) before this cycle's base
        if (c == next_indel) {
            next_indel = c + 1 + Gap(rng, m_indel_max, m_indel_log);
            if (at < templ.size()) {
                double const u = rng.Uniform() * m_indel_max;
                std::size_t run = 1;  // the template's homopolymer here, as far as 5
                while (run < 5 && at + run < templ.size() && templ[at + run] == templ[at]) ++run;
                for (std::size_t back = 1; run < 5 && back <= at && templ[at - back] == templ[at]; ++back) ++run;
                double const factor = run >= 5 ? p.homopolymer_indels : 1.0;
                if (u < p.insertion * factor) {
                    out.seq[c] = "ACGT"[rng.Below(4)];
                    out.qual[c] = m_reported[q];
                    ++made;
                    ++events.insertions;
                    continue;
                }
                if (u < (p.insertion + p.deletion) * factor) {
                    ++at;
                    ++made;
                    ++events.deletions;
                }
            }
        }
        if (at >= templ.size()) {  // past the fragment's end (a short fragment and deletions): adapter-like bases
            out.seq[c] = "ACGT"[rng.Below(4)];
            out.qual[c] = m_reported[std::min(q, 10)];
            continue;
        }
        char const base = templ[at];
        int const code = Code(base);
        bool const gg = at >= 2 && templ[at - 1] == 'G' && templ[at - 2] == 'G';
        ++at;
        if (code > 3 || n_call) {  // an N: of the genome, or no call
            out.seq[c] = 'N';
            out.qual[c] = static_cast<char>(33 + 2);
            ++events.ns;
            continue;
        }
        out.qual[c] = m_reported[q];
        double const error = tail ? p.collapse_error : std::min(0.75, m_error[q] * (gg ? p.gg_factor : 1.0));
        if (rng.Uniform() < error) {
            double const v = rng.Uniform();
            int k = 0;
            while (k < 2 && v >= m_sub_cumulative[code][k]) ++k;
            out.seq[c] = m_sub_base[code][k];
            ++made;
            ++events.substitutions;
        } else {
            out.seq[c] = base;
        }
    }
    if (errors) *errors += made;
}

void IlluminaModel::Pair(std::string_view fragment, double run, LongRng& rng, LongRng& rng2, MadeRead& r1, MadeRead* r2,
                         std::uint64_t* errors) const {
    auto const& p = m_setup.profile;
    double const cluster = run + p.cluster_sd * rng.Gaussian();
    Read(fragment, 0, cluster + p.read_sd * rng.Gaussian(), rng, r1, errors);
    if (!r2) return;
    thread_local std::string reverse;
    reverse.assign(fragment);
    ReverseComplement(reverse);
    double const insert = p.r2_insert_penalty * std::max(0.0, static_cast<double>(fragment.size()) - 500.0);
    Read(reverse, 1, cluster + p.read_sd * rng2.Gaussian() - insert, rng2, *r2, errors);
}

IlluminaReport ReportProfile(IlluminaSetup const& setup, std::uint64_t pairs, std::uint64_t seed) {
    IlluminaModel model(setup);
    IlluminaReport report;
    int const length = setup.read_length;
    LongRng rng(seed);
    std::string genome(1'000'000, 'A');
    for (auto& c : genome) c = "ACGT"[rng.Below(4)];
    std::vector<double> sum[2];
    std::uint64_t q30[2] = {0, 0}, bases[2] = {0, 0};
    IlluminaEvents events[2];
    std::vector<std::uint64_t> share[2];
    for (int r = 0; r < 2; ++r) {
        sum[r].assign(length, 0.0);
        share[r].assign(94, 0);
    }
    MadeRead read;
    std::string fragment;
    double const run = 0.0;
    LongRng rng2(MixSeed(seed, 2));  // read 2's stream, as Pair draws read 2
    for (std::uint64_t i = 0; i < pairs; ++i) {
        std::uint32_t const want = model.FragmentLength(rng, genome.size());
        fragment.assign(genome, rng.Below(genome.size() - want + 1), want);
        double const cluster = run + setup.profile.cluster_sd * rng.Gaussian();
        for (int r = 0; r < 2; ++r) {
            LongRng& stream = r == 0 ? rng : rng2;
            std::string templ = fragment;
            if (r == 1) ReverseComplement(templ);
            double offset = cluster + setup.profile.read_sd * stream.Gaussian();
            if (r == 1) offset -= setup.profile.r2_insert_penalty * std::max(0.0, static_cast<double>(want) - 500.0);
            model.Read(templ, r, offset, stream, read, nullptr, &events[r]);
            for (int c = 0; c < length; ++c) {
                int const q = read.qual[c] - 33;
                sum[r][c] += q;
                q30[r] += q >= 30;
                ++share[r][q];
            }
            bases[r] += length;
        }
    }
    for (int r = 0; r < 2; ++r) {
        report.mean_by_cycle[r].resize(length);
        for (int c = 0; c < length; ++c) report.mean_by_cycle[r][c] = sum[r][c] / static_cast<double>(pairs);
        double const n = static_cast<double>(bases[r]);
        report.q30[r] = q30[r] / n;
        report.substitutions[r] = events[r].substitutions / n;
        report.insertions[r] = events[r].insertions / n;
        report.deletions[r] = events[r].deletions / n;
        report.ns[r] = events[r].ns / n;
        report.quality_share[r].resize(94);
        for (int q = 0; q < 94; ++q) report.quality_share[r][q] = share[r][q] / n;
    }
    return report;
}

// ---- the samples' reads --------------------------------------------------------------------------------------------

namespace {

class PairedJob : public pipeline::Job {
public:
    PairedJob(std::vector<PairedSample> const& samples, PairedOptions const& options)
        : m_samples(samples), m_model(options.setup) {
        for (auto const& sample : samples) {
            if (sample.host_pairs > 0 && !m_host) {
                if (options.host.empty()) throw std::runtime_error(sample.name + ": host pairs without a host (--host_folder)");
                m_host = std::make_unique<Host>(options.host);
            }
            LongRng rng(sample.run_seed);
            m_run.push_back(m_model.RunOffset(rng));
        }
        // a work item's pairs: about 600 kB of FASTQ per read file (so that a streamed sample's R1 and R2 stay in step)
        m_pairs_per_item = std::max<std::uint64_t>(64, 600000 / (2 * static_cast<std::uint64_t>(options.setup.read_length) + 60));
    }

    std::size_t Samples() const override { return m_samples.size(); }
    std::vector<fs::path> Outputs(std::size_t s) const override {
        auto const& sample = m_samples[s];
        if (sample.r2.empty()) return {sample.r1};
        return {sample.r1, sample.r2};
    }
    fs::path Fasta(std::size_t s, int g) const override {
        auto const& sample = m_samples[s];
        return g < static_cast<int>(sample.genomes.size()) ? sample.genomes[g].fasta : fs::path();
    }
    std::string GenomeName(std::size_t s, int g) const override {
        auto const& sample = m_samples[s];
        return g < static_cast<int>(sample.genomes.size()) ? sample.genomes[g].name : "host";
    }
    std::uint32_t MinContigLength() const override { return static_cast<std::uint32_t>(m_model.Setup().read_length); }

    // One round: each genome's pairs, then the host's, in items of m_pairs_per_item.
    std::vector<pipeline::Item> Plan(std::size_t s, pipeline::Totals const& so_far) override {
        if (so_far.rounds > 0) return {};
        auto const& sample = m_samples[s];
        std::vector<pipeline::Item> items;
        std::uint64_t first = 0;
        auto add = [&](int g, std::uint64_t pairs) {
            auto const parts = static_cast<std::uint32_t>((pairs + m_pairs_per_item - 1) / m_pairs_per_item);
            for (std::uint32_t part = 0; part < parts; ++part) {
                pipeline::Item item;
                item.genome = g;
                item.part = part;
                item.parts = parts;
                item.first_read = first;
                item.reads = std::min(m_pairs_per_item, pairs - part * m_pairs_per_item);
                first += item.reads;
                items.push_back(std::move(item));
            }
        };
        for (std::size_t g = 0; g < sample.genomes.size(); ++g) {
            if (sample.genomes[g].pairs > 0) add(static_cast<int>(g), sample.genomes[g].pairs);
        }
        if (sample.host_pairs > 0) add(static_cast<int>(sample.genomes.size()), sample.host_pairs);
        return items;
    }

    pipeline::Piece Make(std::size_t s, pipeline::Item const& item, Contigs const* contigs) const override {
        auto const& sample = m_samples[s];
        bool const host = contigs == nullptr;
        std::uint64_t const seed = host ? sample.host_seed : sample.genomes[item.genome].seed;
        LongRng rng(MixSeed(seed, item.part));
        LongRng rng2(MixSeed(MixSeed(seed, item.part), 0x5232));  // read 2's own stream ("R2")
        std::uint64_t const longest = host ? 0 : contigs->longest;
        int const length = m_model.Setup().read_length;
        pipeline::Piece piece;
        std::string fq[2], fragment, name;
        MadeRead r1, r2;
        bool const both = !sample.r2.empty();
        for (std::uint64_t i = 0; i < item.reads; ++i) {
            if (host) {
                std::uint32_t const want = m_model.FragmentLength(rng, ~std::uint64_t{0} >> 1);
                fragment = m_host->Draw(rng, want, static_cast<std::uint32_t>(length));
                name = "h_" + std::to_string(item.first_read + i + 1);
            } else {
                std::uint32_t const want = m_model.FragmentLength(rng, longest);
                std::size_t k = 0;
                std::uint64_t start = 0;
                bool placed = false;
                for (int attempt = 0; attempt < 1000 && !placed; ++attempt) {  // uniform over the places it fits
                    std::uint64_t const at = rng.Below(contigs->starts.back());
                    k = std::upper_bound(contigs->starts.begin(), contigs->starts.end(), at) - contigs->starts.begin();
                    start = at - (k ? contigs->starts[k - 1] : 0);
                    placed = start + want <= contigs->lengths[k];
                }
                if (!placed) {  // only very short contigs: the longest, from its start
                    k = std::max_element(contigs->lengths.begin(), contigs->lengths.end()) - contigs->lengths.begin();
                    start = 0;
                }
                contigs->Extract(k, start, want, fragment);
                name = contigs->names[k] + "-" + std::to_string(item.first_read + i + 1);
            }
            if (rng.Uniform() < 0.5) ReverseComplement(fragment);
            m_model.Pair(fragment, m_run[s], rng, rng2, r1, both ? &r2 : nullptr, &piece.errors);  // read 1s the same either way
            piece.template_bases += fragment.size();
            for (int r = 0; r < (both ? 2 : 1); ++r) {
                MadeRead const& read = r == 0 ? r1 : r2;
                std::string& out = fq[r];
                out += '@';
                out += name;
                out += r == 0 ? "/1\n" : "/2\n";
                out += read.seq;
                out += "\n+\n";
                out += read.qual;
                out += '\n';
                piece.read_bases += read.seq.size();
            }
        }
        piece.reads = item.reads;
        for (int r = 0; r < (both ? 2 : 1); ++r) piece.bytes[r] = std::move(fq[r]);  // packed by the pipeline
        return piece;
    }

private:
    std::vector<PairedSample> const& m_samples;
    IlluminaModel m_model;
    std::unique_ptr<Host> m_host;
    std::vector<double> m_run;
    std::uint64_t m_pairs_per_item = 2000;
};

}  // namespace

std::vector<pipeline::Totals> SimulatePairs(std::vector<PairedSample> const& samples, PairedOptions const& options) {
    if (samples.empty()) return {};
    PairedJob job(samples, options);
    return pipeline::Run(job, {options.threads, options.genome_store, options.plain_pipes});
}

}  // namespace protal::sim
