// SPDX-License-Identifier: GPL-2.0-only
#include "LongReadSimulator.h"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <condition_variable>
#include <cstdio>
#include <cstring>
#include <exception>
#include <fstream>
#include <future>
#include <iterator>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <utility>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "../IO/Bgzf.h"
#include "../Utilities/Zstd.h"
#include "ThreadedGzStream.h"

namespace fs = std::filesystem;

namespace protal::sim {

// ---- random numbers ---------------------------------------------------------------------------------------------

static std::uint64_t SplitMix64(std::uint64_t& x) {
    std::uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static std::uint64_t Rotl(std::uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

LongRng::LongRng(std::uint64_t seed) {
    for (auto& s : m_s) s = SplitMix64(seed);
}

LongRng::result_type LongRng::operator()() {
    std::uint64_t const result = Rotl(m_s[1] * 5, 7) * 9;
    std::uint64_t const t = m_s[1] << 17;
    m_s[2] ^= m_s[0];
    m_s[3] ^= m_s[1];
    m_s[1] ^= m_s[2];
    m_s[0] ^= m_s[3];
    m_s[2] ^= t;
    m_s[3] = Rotl(m_s[3], 45);
    return result;
}

double LongRng::Uniform() { return static_cast<double>((*this)() >> 11) * 0x1.0p-53; }

std::uint64_t LongRng::Below(std::uint64_t n) {  // Lemire's multiply-and-reject: unbiased
    __uint128_t m = static_cast<__uint128_t>((*this)()) * n;
    auto low = static_cast<std::uint64_t>(m);
    if (low < n) {
        std::uint64_t const threshold = -n % n;
        while (low < threshold) {
            m = static_cast<__uint128_t>((*this)()) * n;
            low = static_cast<std::uint64_t>(m);
        }
    }
    return static_cast<std::uint64_t>(m >> 64);
}

double LongRng::Normal() {  // Marsaglia's polar method
    if (m_has_spare) {
        m_has_spare = false;
        return m_spare;
    }
    double u, v, s;
    do {
        u = 2.0 * Uniform() - 1.0;
        v = 2.0 * Uniform() - 1.0;
        s = u * u + v * v;
    } while (s >= 1.0 || s == 0.0);
    double const f = std::sqrt(-2.0 * std::log(s) / s);
    m_spare = v * f;
    m_has_spare = true;
    return u * f;
}

std::uint64_t MixSeed(std::uint64_t a, std::uint64_t b) {
    std::uint64_t x = a ^ Rotl(b, 32) ^ 0x2545f4914f6cdd1dULL;
    SplitMix64(x);
    return SplitMix64(x);
}

// Gamma(shape, scale) by Marsaglia and Tsang, on LongRng alone (no library distribution: the same draws everywhere).
static double Gamma(LongRng& rng, double shape, double scale) {
    if (shape < 1.0) {
        double const u = rng.Uniform();
        return Gamma(rng, shape + 1.0, scale) * std::pow(u > 0 ? u : 0x1.0p-53, 1.0 / shape);
    }
    double const d = shape - 1.0 / 3.0, c = 1.0 / std::sqrt(9.0 * d);
    for (;;) {
        double x, v;
        do {
            x = rng.Normal();
            v = 1.0 + c * x;
        } while (v <= 0.0);
        v = v * v * v;
        double const u = rng.Uniform();
        if (u < 1.0 - 0.0331 * x * x * x * x) return d * v * scale;
        if (u > 0 && std::log(u) < 0.5 * x * x + d * (1.0 - v + std::log(v))) return d * v * scale;
    }
}

// ---- the setup ---------------------------------------------------------------------------------------------------

static std::vector<std::string> Split(std::string const& text, char by) {
    std::vector<std::string> parts;
    std::string part;
    std::istringstream in(text);
    while (std::getline(in, part, by)) parts.push_back(part);
    if (!text.empty() && text.back() == by) parts.emplace_back();
    return parts;
}

static double Number(std::string const& value, std::string const& text) {
    try {
        std::size_t used = 0;
        double const v = std::stod(value, &used);
        if (used == value.size()) return v;
    } catch (std::exception const&) {
    }
    throw std::invalid_argument("read setup " + text + ": " + value + " is not a number");
}

LongReadSetup LongReadSetup::Parse(std::string const& text) {
    auto const parts = Split(text, ':');
    LongReadSetup setup;
    if (parts.empty()) throw std::invalid_argument("an empty read setup");
    if (parts[0] == "ultima") {
        if (parts.size() != 5) throw std::invalid_argument("read setup " + text + ": expected ultima:LENGTH_MEAN:LENGTH_SD:Q_MEAN:Q_SD");
        setup.method = Method::Ultima;
        setup.length_mean = static_cast<int>(Number(parts[1], text));
        setup.length_sd = static_cast<int>(Number(parts[2], text));
        setup.q_mean = Number(parts[3], text);
        setup.q_sd = Number(parts[4], text);
    } else if (parts[0] == "hifi") {
        if (parts.size() != 4) throw std::invalid_argument("read setup " + text + ": expected hifi:LENGTH_MEAN:LENGTH_SD:Q_SD");
        setup.method = Method::Hifi;
        setup.length_mean = static_cast<int>(Number(parts[1], text));
        setup.length_sd = static_cast<int>(Number(parts[2], text));
        setup.q_sd = Number(parts[3], text);
    } else if (parts[0] == "qshmm") {
        if (parts.size() != 5 && parts.size() != 6) {
            throw std::invalid_argument("read setup " + text +
                                        ": expected qshmm:MODEL:LENGTH_MEAN:LENGTH_SD:ACCURACY_MEAN[:SUB/INS/DEL]");
        }
        setup.method = Method::Qshmm;
        setup.model = parts[1];
        setup.length_mean = static_cast<int>(Number(parts[2], text));
        setup.length_sd = static_cast<int>(Number(parts[3], text));
        setup.accuracy = Number(parts[4], text);
        if (parts.size() == 6) {
            auto const ratio = Split(parts[5], '/');
            if (ratio.size() != 3) throw std::invalid_argument("read setup " + text + ": the difference ratio is SUB/INS/DEL");
            setup.sub_ratio = static_cast<long>(Number(ratio[0], text));
            setup.ins_ratio = static_cast<long>(Number(ratio[1], text));
            setup.del_ratio = static_cast<long>(Number(ratio[2], text));
            if (setup.sub_ratio < 0 || setup.ins_ratio < 0 || setup.del_ratio < 0 ||
                setup.sub_ratio + setup.ins_ratio + setup.del_ratio <= 0) {
                throw std::invalid_argument("read setup " + text + ": the difference ratio needs positive parts");
            }
        }
        if (!(setup.accuracy > 0 && setup.accuracy <= 1)) {
            throw std::invalid_argument("read setup " + text + ": the accuracy mean is a fraction, e.g. 0.97");
        }
    } else if (parts[0] == "errhmm") {
        throw std::invalid_argument("read setup " + text + ": pbsim3's errhmm model is not simulated in process (its "
                                    "reads have every quality at Q0); use hifi: for PacBio, qshmm: for Nanopore");
    } else {
        throw std::invalid_argument("read setup " + text + ": expected hifi:, ultima: or qshmm:");
    }
    if (setup.length_mean <= 0 || setup.length_sd < 0) {
        throw std::invalid_argument("read setup " + text + ": the length mean must be positive and its SD not negative");
    }
    return setup;
}

// ---- hifi_reads.py's models --------------------------------------------------------------------------------------

namespace {
    constexpr char kBases[4] = {'A', 'C', 'G', 'T'};

    inline int BaseCode(char c) {
        switch (c) {
            case 'A': return 0;
            case 'C': return 1;
            case 'G': return 2;
            case 'T': return 3;
            default: return 4;
        }
    }
}

namespace hifi {
    double LengthQuality(double length) {
        if (length <= 5000) return 50.0;
        if (length <= 25000) return 50.0 + (length - 5000) * (30.0 - 50.0) / 20000.0;
        // 25 kb to 50 kb, and on at that slope beyond
        return 30.0 + (length - 25000) * (20.0 - 30.0) / 25000.0;
    }

    void Mutate(std::string_view templ, LongRng& rng, LongReadSetup const& setup, MadeRead& out,
                std::uint64_t* errors) {
        bool const flow = setup.method == LongReadSetup::Method::Ultima;
        std::size_t const n = templ.size();
        out.seq.clear();
        out.qual.clear();
        if (n == 0) return;
        thread_local std::vector<std::uint32_t> run;
        thread_local std::vector<double> weight;
        run.resize(n);
        weight.resize(n);
        for (std::size_t i = 0; i < n;) {  // homopolymer runs
            std::size_t j = i + 1;
            while (j < n && templ[j] == templ[i]) ++j;
            for (std::size_t k = i; k < j; ++k) run[k] = static_cast<std::uint32_t>(j - i);
            i = j;
        }
        std::uint32_t const homopolymer = flow ? kFlowHomopolymer : kHomopolymer;
        double sum = 0.0, log_sum = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            double const hp = run[i] >= homopolymer ? (run[i] / 2.0) * (run[i] / 2.0) : 1.0;
            weight[i] = hp * std::exp(kConfidenceSigma * rng.Normal());
            sum += weight[i];
            if (flow) log_sum += std::log10(weight[i]);
        }
        // the bases' error probabilities, in weight (reused)
        if (flow) {
            double const q_read = std::clamp(setup.q_mean + setup.q_sd * rng.Normal(), kQMin, kQMax);
            double const mean_log = log_sum / static_cast<double>(n);
            for (std::size_t i = 0; i < n; ++i) {
                double const phred = q_read - 10.0 * (std::log10(weight[i]) - mean_log);
                weight[i] = std::min(std::pow(10.0, -phred / 10.0), kPMax);
            }
        } else {
            double const q_read = std::clamp(LengthQuality(static_cast<double>(n)) + setup.q_sd * rng.Normal(), kQMin, kQMax);
            double const scale = std::pow(10.0, -q_read / 10.0) * static_cast<double>(n) / sum;
            for (std::size_t i = 0; i < n; ++i) weight[i] = std::min(weight[i] * scale, kPMax);
        }
        out.seq.reserve(n + n / 8 + 16);
        out.qual.reserve(n + n / 8 + 16);
        std::uint64_t made = 0;
        for (std::size_t i = 0; i < n; ++i) {
            double const p = weight[i];
            char const a = templ[i];
            char const q = static_cast<char>(33 + std::clamp(std::nearbyint(-10.0 * std::log10(p)), 1.0, 93.0));
            if (!(rng.Uniform() < p)) {
                out.seq += a;
                out.qual += q;
                continue;
            }
            ++made;
            bool const hp = run[i] >= homopolymer;
            double const kind = rng.Uniform();
            int const code = BaseCode(a);
            if (!hp && kind < kSubstitution && code < 4) {  // a substitution
                out.seq += kBases[(code + 1 + static_cast<int>(rng.Below(3))) % 4];
                out.qual += q;
            } else if (hp ? kind < 0.5 : (kind >= kSubstitution && kind < kSubstitution + kInsertion)) {
                out.seq += a;  // the base, then the inserted one (in a homopolymer, one more of the run)
                out.qual += q;
                out.seq += hp ? a : kBases[rng.Below(4)];
                out.qual += q;
            }  // else a deletion
        }
        if (errors) *errors += made;
    }
}

// ---- pbsim3's qshmm model ----------------------------------------------------------------------------------------
// After simulate_by_qshmm_templ, set_qshmm, set_mut and the quality tables of main in pbsim3's src/pbsim.cpp (pbsim3
// 3.0.x, Yukiteru Ono, GPL-2.0): the same tables, quantized alike; LongRng's Below(n) in place of rand() % n.

QshmmModel::QshmmModel(fs::path const& file, double accuracy_mean, long sub_ratio, long ins_ratio, long del_ratio) {
    constexpr int A = kAccuracyMax + 1, S = kStates + 1, Q = kQualities;
    std::vector<double> ip(A * S, 0.0), ep(static_cast<std::size_t>(A) * S * Q, 0.0), tp(static_cast<std::size_t>(A) * S * S, 0.0);
    m_exists.assign(A, 0);
    std::ifstream in(file);
    if (!in) throw std::runtime_error("cannot read the qshmm model " + file.string());
    std::string line;
    std::size_t line_no = 0;
    while (std::getline(in, line)) {
        ++line_no;
        if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue;
        std::istringstream fields(line);
        int accuracy = -1, state = -1;
        std::string kind;
        fields >> accuracy >> kind >> state;
        auto bad = [&](std::string const& why) {
            return std::runtime_error("qshmm model " + file.string() + ", line " + std::to_string(line_no) + ": " + why);
        };
        if (!fields || accuracy < 0 || accuracy > kAccuracyMax) throw bad("expected ACCURACY IP|EP|TP STATE values");
        if (state < 1) throw bad("state " + std::to_string(state) + " below 1");
        m_exists[accuracy] = 1;
        // pbsim3's arrays are [101][51](...), flat: a model with more than 50 states (QSHMM-ONT-HQ has 56 at accuracy
        // 92) writes its states 51+ into the next accuracy's first rows (and transitions to them into the next state's
        // first entries), and pbsim3's reads follow from that. Written alike here, within the arrays; then only the
        // states 1-50 are read, as there.
        auto at = [&](std::size_t index, std::size_t size) {
            if (index >= size) throw bad("state or value beyond pbsim3's tables");
            return index;
        };
        double v = 0;
        if (kind == "IP") {
            if (!(fields >> v)) throw bad("IP without a value");
            ip[at(static_cast<std::size_t>(accuracy) * S + state, ip.size())] = v;
        } else if (kind == "EP") {
            std::size_t num = 0;
            while (fields >> v) {
                if (num >= Q) throw bad("more than 94 emission probabilities");
                ep[at((static_cast<std::size_t>(accuracy) * S + state) * Q + num++, ep.size())] = v;
            }
        } else if (kind == "TP") {
            std::size_t num = 0;
            while (fields >> v) tp[at((static_cast<std::size_t>(accuracy) * S + state) * S + ++num, tp.size())] = v;
        } else {
            throw bad("unknown record " + kind);
        }
    }

    double qprob[Q];
    for (int i = 0; i < Q; ++i) qprob[i] = std::pow(10, static_cast<double>(i) / -10);
    std::vector<double> uni(A * Q, 0.0);  // uni_ep: a quality (or two) that gives the accuracy
    for (int i = 0; i <= kAccuracyMax; ++i) {
        if (i == kAccuracyMax) {
            uni[i * Q + 93] = 1.0;
            continue;
        }
        double const prob = 1.0 - i / 100.0;
        for (int j = 0; j < Q; ++j) {
            if (prob == qprob[j]) {
                uni[i * Q + j] = 1.0;
                break;
            } else if (prob > qprob[j]) {
                double const rate = (prob - qprob[j]) / (qprob[j - 1] - qprob[j]);
                uni[i * Q + j - 1] = rate;
                uni[i * Q + j] = 1 - rate;
                break;
            }
        }
    }

    // pbsim3 keeps --accuracy-mean to two decimals, and draws a read's level from [0.75, 1.05] x the mean
    double const mean = static_cast<int>(accuracy_mean * 100) * 0.01 * 100;
    long accuracy_max = static_cast<long>(std::floor(mean * 1.05));
    long const accuracy_min = static_cast<long>(std::floor(mean * 0.75));
    if (accuracy_max > kAccuracyMax) accuracy_max = kAccuracyMax;
    if (accuracy_min < 0 || accuracy_min > accuracy_max) throw std::runtime_error("qshmm: the accuracy mean is out of range");
    m_accuracy_min = static_cast<int>(accuracy_min);
    m_accuracy_max = static_cast<int>(accuracy_max);

    long start_wk = 1, end_wk = 0;  // pbsim3 carries end_wk from one table to the next, as here
    m_prob2accuracy.assign(100001, 0);
    double freq_total = 0.0;
    for (long i = accuracy_min; i <= accuracy_max; ++i) freq_total += std::exp(0.22 * i);
    double total = 0.0;
    for (long i = accuracy_min; i <= accuracy_max; ++i) {
        total += std::exp(0.22 * i) / freq_total;
        end_wk = std::min<long>(static_cast<long>(total * 100000 + 0.5), 100000);
        for (long j = start_wk; j <= end_wk; ++j) m_prob2accuracy[j] = static_cast<int>(i);
        if (end_wk >= 100000) break;
        start_wk = end_wk + 1;
    }
    m_accuracy_draw = end_wk;
    if (m_accuracy_draw < 1) throw std::runtime_error("qshmm: the accuracy parameters are not appropriate");

    m_init_draw.assign(A, 0);
    m_init2state.assign(A * 101, 0);
    m_emis_draw.assign(A * S, 0);
    m_tran_draw.assign(A * S, 0);
    m_emis2qc.assign(static_cast<std::size_t>(A) * S * 101, 0);
    m_tran2state.assign(static_cast<std::size_t>(A) * S * 101, 0);
    m_freq_draw.assign(A, 0);
    m_freq2qc.assign(A * 1001, 0);
    // A cumulative table of `steps` entries ([1..steps]) of the nonzero probabilities p(k), k in [first, last].
    auto table = [&](auto prob, int first, int last, long steps, int* out) {
        start_wk = 1;
        double sum = 0.0;
        for (int k = first; k <= last; ++k) {
            double const p = prob(k);
            if (p == 0) continue;
            sum += p;
            end_wk = std::min<long>(static_cast<long>(sum * steps + 0.5), steps);
            for (long l = start_wk; l <= end_wk; ++l) out[l] = k;
            if (end_wk >= steps) break;
            start_wk = end_wk + 1;
        }
        return end_wk;
    };
    for (long i = accuracy_min; i <= accuracy_max; ++i) {
        if (m_exists[i]) {
            m_init_draw[i] = table([&](int j) { return ip[i * S + j]; }, 1, kStates, 100, &m_init2state[i * 101]);
            for (int j = 1; j <= kStates; ++j) {
                m_emis_draw[i * S + j] = table([&](int k) { return ep[(static_cast<std::size_t>(i) * S + j) * Q + k]; },
                                               0, 93, 100, &m_emis2qc[(static_cast<std::size_t>(i) * S + j) * 101]);
            }
            for (int j = 1; j <= kStates; ++j) {
                m_tran_draw[i * S + j] = table([&](int k) { return tp[(static_cast<std::size_t>(i) * S + j) * S + k]; },
                                               1, kStates, 100, &m_tran2state[(static_cast<std::size_t>(i) * S + j) * 101]);
            }
        } else {
            m_freq_draw[i] = table([&](int j) { return uni[i * Q + j]; }, 0, 93, 1000, &m_freq2qc[i * 1001]);
        }
    }

    long const ratio_sum = sub_ratio + ins_ratio + del_ratio;
    double const sub_rate = static_cast<double>(sub_ratio) / ratio_sum, ins_rate = static_cast<double>(ins_ratio) / ratio_sum,
                 del_rate = static_cast<double>(del_ratio) / ratio_sum;
    for (int i = 0; i < Q; ++i) {
        m_sub_thre[i] = static_cast<long>((qprob[i] * sub_rate) * 1000000 + 0.5);
        m_ins_thre[i] = static_cast<long>((qprob[i] * (sub_rate + ins_rate)) * 1000000 + 0.5);
        m_del_thre[i] = static_cast<long>((qprob[i] * del_rate) / (1 + qprob[i] * del_rate) * 1000000 + 0.5);
    }
}

int QshmmModel::DrawAccuracy(LongRng& rng) const {
    return m_prob2accuracy[rng.Below(static_cast<std::uint64_t>(m_accuracy_draw)) + 1];
}

void QshmmModel::Mutate(std::string_view templ, LongRng& rng, MadeRead& out, std::uint64_t* errors) const {
    constexpr int S = kStates + 1;
    int const acc = DrawAccuracy(rng);
    bool const hmm = m_exists[acc];
    std::size_t const n = templ.size();
    out.seq.clear();
    out.qual.clear();
    out.seq.reserve(n + n / 8 + 16);
    out.qual.reserve(n + n / 8 + 16);
    auto draw = [&](long values) {
        if (values < 1) throw std::runtime_error("qshmm: the model has no probabilities for a state it reaches");
        return rng.Below(static_cast<std::uint64_t>(values)) + 1;
    };
    std::uint64_t made = 0;
    int state = 0;
    std::size_t ref = 0;
    while (ref < n) {
        int q;
        if (hmm) {
            state = out.qual.empty() ? m_init2state[acc * 101 + draw(m_init_draw[acc])]
                                     : m_tran2state[(static_cast<std::size_t>(acc) * S + state) * 101 + draw(m_tran_draw[acc * S + state])];
            q = m_emis2qc[(static_cast<std::size_t>(acc) * S + state) * 101 + draw(m_emis_draw[acc * S + state])];
        } else {
            q = m_freq2qc[acc * 1001 + draw(m_freq_draw[acc])];
        }
        char const nt = templ[ref];
        auto const r = static_cast<long>(rng.Below(1000000));
        char base;
        if (r < m_sub_thre[q]) {
            ++made;
            switch (nt) {
                case 'A': base = "TGC"[rng.Below(3)]; break;
                case 'T': base = "AGC"[rng.Below(3)]; break;
                case 'G': base = "ATC"[rng.Below(3)]; break;
                case 'C': base = "ATG"[rng.Below(3)]; break;
                default: base = "ATGC"[rng.Below(4)];
            }
            ++ref;
        } else if (r < m_ins_thre[q]) {
            ++made;
            auto const i = rng.Below(8);
            base = i >= 4 ? nt : "ATGC"[i];
        } else {
            base = nt;
            ++ref;
        }
        out.seq += base;
        out.qual += static_cast<char>(q + 33);
        while (ref < n && static_cast<long>(rng.Below(1000000)) < m_del_thre[q]) {  // deletions after the base
            ++made;
            ++ref;
        }
    }
    if (errors) *errors += made;
}

// ---- templates ---------------------------------------------------------------------------------------------------

constexpr std::uint32_t kMinLength = 100, kMaxLength = 1000000;  // pbsim3's --length-min and --length-max

std::uint32_t DrawReadLength(LongRng& rng, double mean, double sd) {
    if (sd <= 0) {
        return static_cast<std::uint32_t>(std::clamp<long long>(static_cast<long long>(mean), kMinLength, kMaxLength));
    }
    double const shape = (mean / sd) * (mean / sd), scale = sd * sd / mean;
    for (;;) {
        long long const length = std::llround(Gamma(rng, shape, scale));
        if (length >= kMinLength && length <= kMaxLength) return static_cast<std::uint32_t>(length);
    }
}

void ReverseComplement(std::string& seq) {
    std::reverse(seq.begin(), seq.end());
    for (char& c : seq) {
        switch (c) {
            case 'A': c = 'T'; break;
            case 'C': c = 'G'; break;
            case 'G': c = 'C'; break;
            case 'T': c = 'A'; break;
            default: break;  // N and others as they are
        }
    }
}

namespace {

// A genome's contigs of kMinLength bases or more, upper case, and their cumulative start ranges.
struct Contigs {
    std::vector<std::string> seqs;
    std::vector<std::uint64_t> starts;  // cumulative (length - kMinLength + 1)
};

std::shared_ptr<Contigs const> LoadContigs(LongGenome const& genome) {
    protal::ThreadedGzIstream in(genome.fasta.c_str());
    if (!in.rdbuf()->is_open()) throw std::runtime_error(genome.name + ": cannot read " + genome.fasta.string());
    auto contigs = std::make_shared<Contigs>();
    std::string line, seq;
    bool in_record = false;
    auto close = [&] {
        if (in_record && seq.size() >= kMinLength) {
            contigs->starts.push_back((contigs->starts.empty() ? 0 : contigs->starts.back()) + seq.size() - kMinLength + 1);
            contigs->seqs.push_back(std::move(seq));
        }
        seq.clear();
    };
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '>') {
            close();
            in_record = true;
            continue;
        }
        if (!in_record) continue;
        for (char c : line) {
            if (!std::isspace(static_cast<unsigned char>(c))) seq += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        }
    }
    if (in.rdbuf()->read_failed()) {
        throw std::runtime_error(genome.name + ": " + genome.fasta.string() + " is truncated or corrupt (" +
                                 in.rdbuf()->read_error_message() + ")");
    }
    close();
    if (contigs->seqs.empty()) {
        throw std::runtime_error(genome.name + ": no sequence of " + std::to_string(kMinLength) + " bases or more in " +
                                 genome.fasta.string());
    }
    return contigs;
}

// A host genome prepared by scenarios.prepare_host: host.seq (every contig one after the other) by memory map, and
// host.json's contigs ([name, offset, length], ...).
class Host {
public:
    explicit Host(fs::path const& folder) {
        std::ifstream in(folder / "host.json");
        if (!in) throw std::runtime_error("cannot read the host index " + (folder / "host.json").string());
        std::string const json((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        ParseContigs(json, folder);
        fs::path const seq = folder / "host.seq";
        m_fd = ::open(seq.c_str(), O_RDONLY);
        if (m_fd < 0) throw std::runtime_error("cannot read " + seq.string() + ": " + std::strerror(errno));
        struct stat st {};
        if (::fstat(m_fd, &st) != 0) throw std::runtime_error("cannot read " + seq.string() + ": " + std::strerror(errno));
        m_size = static_cast<std::size_t>(st.st_size);
        for (std::size_t k = 0; k < m_offsets.size(); ++k) {
            if (m_offsets[k] + m_lengths[k] > m_size) throw std::runtime_error(seq.string() + " is shorter than its index says");
        }
        void* map = ::mmap(nullptr, m_size, PROT_READ, MAP_SHARED, m_fd, 0);
        if (map == MAP_FAILED) throw std::runtime_error("cannot map " + seq.string() + ": " + std::strerror(errno));
        m_data = static_cast<char const*>(map);
    }
    ~Host() {
        if (m_data) ::munmap(const_cast<char*>(m_data), m_size);
        if (m_fd >= 0) ::close(m_fd);
    }
    Host(Host const&) = delete;
    Host& operator=(Host const&) = delete;

    // As scenarios.Host.draw: `length` bases from a random place (shorter at a contig's end, but at least
    // min(length, 100)), on the forward strand; one too short or with more than 10% N drawn again, up to 50 times.
    std::string Draw(LongRng& rng, std::uint32_t length) const {
        std::uint64_t const least = std::min<std::uint64_t>(length, kMinLength);
        std::string seq;
        for (int attempt = 0; attempt < 50; ++attempt) {
            std::uint64_t const at = rng.Below(m_bases);
            std::size_t const k = std::upper_bound(m_ends.begin(), m_ends.end(), at) - m_ends.begin();
            std::uint64_t const start = at - (m_ends[k] - m_lengths[k]);
            std::uint64_t const end = std::min<std::uint64_t>(m_lengths[k], start + length);
            seq.assign(m_data + m_offsets[k] + start, end - start);
            if (seq.size() >= least && std::count(seq.begin(), seq.end(), 'N') <= 0.1 * static_cast<double>(seq.size())) break;
        }
        return seq;
    }

private:
    void ParseContigs(std::string const& json, fs::path const& folder) {
        auto fail = [&](std::string const& why) {
            return std::runtime_error("host index " + (folder / "host.json").string() + ": " + why);
        };
        std::size_t at = json.find("\"contigs\"");
        if (at == std::string::npos) throw fail("no contigs");
        at = json.find('[', at);
        if (at == std::string::npos) throw fail("no contigs");
        ++at;
        auto skip = [&] {
            while (at < json.size() && std::isspace(static_cast<unsigned char>(json[at]))) ++at;
        };
        auto number = [&]() -> std::uint64_t {
            skip();
            std::size_t used = 0;
            std::uint64_t const v = std::stoull(json.substr(at, 24), &used);
            at += used;
            return v;
        };
        auto expect = [&](char c) {
            skip();
            if (at >= json.size() || json[at] != c) throw fail(std::string("expected '") + c + "'");
            ++at;
        };
        for (;;) {
            skip();
            if (at < json.size() && json[at] == ']') break;
            expect('[');
            expect('"');
            while (at < json.size() && json[at] != '"') at += json[at] == '\\' ? 2 : 1;  // the name, not needed
            expect('"');
            expect(',');
            std::uint64_t const offset = number();
            expect(',');
            std::uint64_t const length = number();
            expect(']');
            m_offsets.push_back(offset);
            m_lengths.push_back(length);
            m_bases += length;
            m_ends.push_back(m_bases);
            skip();
            if (at < json.size() && json[at] == ',') ++at;
        }
        if (m_ends.empty()) throw fail("no contigs");
    }

    std::vector<std::uint64_t> m_offsets, m_lengths, m_ends;
    std::uint64_t m_bases = 0;
    int m_fd = -1;
    std::size_t m_size = 0;
    char const* m_data = nullptr;
};

// ---- compressed output -------------------------------------------------------------------------------------------

enum class Packing { Zstd, Bgzf };

Packing PackingOf(fs::path const& out) {
    std::string const name = out.filename().string();
    auto ends = [&](std::string const& suffix) {
        return name.size() > suffix.size() && name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0;
    };
    if (ends(".zst")) return Packing::Zstd;
    if (ends(".gz")) return Packing::Bgzf;
    throw std::runtime_error(out.string() + ": a sample's reads go to .fq.zst (zstd) or .fq.gz (BGZF)");
}

// A piece of a sample's file: one zstd frame (level 3, with a checksum, as the paired-end reads are written), or
// BGZF blocks (the end-of-file block comes when the file is closed). Concatenated, they are the file.
std::string Pack(std::string const& data, Packing packing) {
    std::string out;
    if (packing == Packing::Bgzf) {
        if (!protal::bgzf::Compress(data.data(), data.size(), out)) throw std::runtime_error("BGZF compression failed");
        return out;
    }
    struct Free {
        void operator()(ZSTD_CCtx* c) const { ZSTD_freeCCtx(c); }
    };
    thread_local std::unique_ptr<ZSTD_CCtx, Free> cctx;
    if (!cctx) {
        std::string error;
        cctx.reset(protal::zstd::MakeCCtx(protal::zstd::Params{3, 0, 1, 0}, error));
        if (!cctx) throw std::runtime_error(error);
    }
    out.resize(ZSTD_compressBound(data.size()));
    std::size_t const size = ZSTD_compress2(cctx.get(), out.data(), out.size(), data.data(), data.size());
    if (ZSTD_isError(size)) throw std::runtime_error(std::string("zstd compression failed: ") + ZSTD_getErrorName(size));
    out.resize(size);
    return out;
}

// ---- the samples' reads on threads ---------------------------------------------------------------------------------

constexpr std::uint64_t kPartBases = 2'000'000;  // template bases of a work item at most (one read more if longer)

struct Item {
    int genome = 0;
    std::uint32_t round = 0, part = 0, parts = 1;
    std::uint64_t first_read = 0;  // reads of the sample before this item's
    std::vector<std::uint32_t> lengths;
};

struct Piece {
    std::string bytes;
    std::uint64_t reads = 0, template_bases = 0, read_bases = 0, errors = 0;
};

struct GenomeSlot {
    std::shared_future<std::shared_ptr<Contigs const>> loaded;
    std::uint32_t parts_left = 0;
};

struct Stream {
    LongSample const* sample = nullptr;
    std::size_t index = 0;
    Packing packing = Packing::Zstd;
    LongRng plan{0};
    std::vector<double> cumulative;
    std::vector<Item> items;
    std::vector<std::optional<Piece>> pieces;
    std::size_t next_issue = 0, next_write = 0, processed = 0, round_end = 0;
    bool planned_all = false, writing = false, closed = false;
    std::uint64_t reads_planned = 0;
    std::map<std::pair<std::uint32_t, int>, GenomeSlot> genomes;
    std::FILE* out = nullptr;
    fs::path partial;
    LongSampleResult result;
};

class Engine {
public:
    Engine(std::vector<LongSample> const& samples, LongReadOptions const& options)
        : m_samples(samples), m_options(options), m_streams(samples.size()) {
        if (options.setup.method == LongReadSetup::Method::Qshmm) {
            if (options.model.empty()) throw std::runtime_error("a qshmm setup needs its model file (--long_model)");
            m_qshmm = std::make_unique<QshmmModel>(options.model, options.setup.accuracy, options.setup.sub_ratio,
                                                   options.setup.ins_ratio, options.setup.del_ratio);
        }
        for (auto const& sample : samples) {
            PackingOf(sample.out);  // fails before any work
            for (auto const& genome : sample.genomes) {
                if (genome.host && !m_hosts.count(genome.fasta.string())) {
                    m_hosts.emplace(genome.fasta.string(), std::make_unique<Host>(genome.fasta));
                }
            }
        }
        m_cap = std::max<std::size_t>(4, 3 * static_cast<std::size_t>(std::max(1, options.threads)));
    }

    ~Engine() {
        for (auto& s : m_streams) {
            if (s.out) {
                std::fclose(s.out);
                std::error_code ec;
                fs::remove(s.partial, ec);
            }
        }
    }

    std::vector<LongSampleResult> Run() {
        int const threads = std::max(1, m_options.threads);
        std::vector<std::thread> pool;
        for (int t = 1; t < threads; ++t) pool.emplace_back([this] { Work(); });
        Work();
        for (auto& thread : pool) thread.join();
        if (m_error) std::rethrow_exception(m_error);
        std::vector<LongSampleResult> results;
        for (auto const& s : m_streams) results.push_back(s.result);
        return results;
    }

private:
    void Work() {
        std::unique_lock<std::mutex> lock(m_mutex);
        for (;;) {
            if (m_error || m_closed == m_samples.size()) break;
            Stream* stream = nullptr;
            if (m_in_flight < m_cap) {
                for (Stream* s : m_active) {
                    if (s->next_issue < s->round_end) {
                        stream = s;
                        break;
                    }
                }
                if (!stream && m_next_sample < m_samples.size()) {  // all open samples wait: open the next
                    try {
                        Open(m_next_sample++, lock);
                    } catch (...) {
                        Fail(std::current_exception());
                    }
                    continue;
                }
            }
            if (!stream) {
                m_cv.wait(lock);
                continue;
            }
            std::size_t const k = stream->next_issue++;
            ++m_in_flight;
            Item const item = stream->items[k];
            LongGenome const& genome = stream->sample->genomes[item.genome];
            std::shared_ptr<std::promise<std::shared_ptr<Contigs const>>> load;
            std::shared_future<std::shared_ptr<Contigs const>> loaded;
            if (!genome.host) {
                auto [slot, fresh] = stream->genomes.try_emplace({item.round, item.genome});
                if (fresh) {
                    load = std::make_shared<std::promise<std::shared_ptr<Contigs const>>>();
                    slot->second.loaded = load->get_future().share();
                    slot->second.parts_left = item.parts;
                }
                loaded = slot->second.loaded;
            }
            lock.unlock();
            std::optional<Piece> piece;
            try {
                if (load) {
                    try {
                        load->set_value(LoadContigs(genome));
                    } catch (...) {
                        load->set_exception(std::current_exception());
                    }
                }
                piece = Make(*stream, item, genome.host ? nullptr : loaded.get().get());
            } catch (...) {
                lock.lock();
                Fail(std::current_exception());
                continue;
            }
            lock.lock();
            Complete(*stream, k, item, std::move(*piece), lock);
        }
        m_cv.notify_all();
    }

    void Fail(std::exception_ptr error) {
        if (!m_error) m_error = error;
        m_cv.notify_all();
    }

    // Opens sample i's file and plans its first round (under the lock).
    void Open(std::size_t i, std::unique_lock<std::mutex>& lock) {
        Stream& s = m_streams[i];
        LongSample const& sample = m_samples[i];
        s.sample = &sample;
        s.index = i;
        s.packing = PackingOf(sample.out);
        s.plan = LongRng(sample.seed);
        s.result.name = sample.name;
        double total = 0.0;
        for (auto const& genome : sample.genomes) {
            total += std::max(0.0, genome.weight);
            s.cumulative.push_back(total);
        }
        if (sample.bases > 0 && !(total > 0)) {
            throw std::runtime_error(sample.name + ": no genome with reads to simulate (relative abundances and lengths are 0)");
        }
        if (!sample.out.parent_path().empty()) fs::create_directories(sample.out.parent_path());
        s.partial = sample.out.string() + ".partial";
        s.out = std::fopen(s.partial.c_str(), "wb");
        if (!s.out) throw std::runtime_error("cannot write " + s.partial.string() + ": " + std::strerror(errno));
        m_active.push_back(&s);
        PlanRound(s);
        Drain(s, lock);  // a sample of no bases is done here
        m_cv.notify_all();
    }

    // The reads of the bases the sample still lacks, genome by genome (collect_training_data.draw_templates's round).
    void PlanRound(Stream& s) {
        LongSample const& sample = *s.sample;
        if (s.result.template_bases >= sample.bases) {
            s.planned_all = true;
            return;
        }
        auto need = static_cast<long long>(sample.bases - s.result.template_bases);
        std::map<int, std::vector<std::uint32_t>> planned;
        double const total = s.cumulative.back();
        while (need > 0) {
            auto g = static_cast<int>(std::upper_bound(s.cumulative.begin(), s.cumulative.end(), s.plan.Uniform() * total) -
                                      s.cumulative.begin());
            g = std::min(g, static_cast<int>(s.cumulative.size()) - 1);
            std::uint32_t const length = DrawReadLength(s.plan, m_options.setup.length_mean, m_options.setup.length_sd);
            planned[g].push_back(length);
            need -= length;
        }
        std::uint32_t const round = s.result.rounds++;
        for (auto& [g, lengths] : planned) {
            std::vector<Item> parts;
            std::uint64_t bases = kPartBases;
            for (std::uint32_t length : lengths) {
                if (bases >= kPartBases) {
                    parts.push_back(Item{g, round, static_cast<std::uint32_t>(parts.size()), 1, 0, {}});
                    bases = 0;
                }
                parts.back().lengths.push_back(length);
                bases += length;
            }
            for (auto& item : parts) {
                item.parts = static_cast<std::uint32_t>(parts.size());
                item.first_read = s.reads_planned;
                s.reads_planned += item.lengths.size();
                s.items.push_back(std::move(item));
            }
        }
        s.round_end = s.items.size();
        s.pieces.resize(s.items.size());
    }

    // An item's reads, as FASTQ compressed (no lock held).
    Piece Make(Stream const& s, Item const& item, Contigs const* contigs) const {
        LongSample const& sample = *s.sample;
        LongGenome const& genome = sample.genomes[item.genome];
        Host const* host = genome.host ? m_hosts.at(genome.fasta.string()).get() : nullptr;
        LongRng rng(MixSeed(MixSeed(sample.seed, item.round), MixSeed(static_cast<std::uint64_t>(item.genome), item.part)));
        Piece piece;
        std::string fastq, templ;
        MadeRead read;
        std::string const prefix = "@g" + std::to_string(item.genome) + "x_";
        for (std::size_t i = 0; i < item.lengths.size(); ++i) {
            std::uint32_t const length = item.lengths[i];
            if (host) {
                templ = host->Draw(rng, length);
            } else {
                auto const& starts = contigs->starts;
                std::uint64_t const at = rng.Below(starts.back());
                std::size_t const k = std::upper_bound(starts.begin(), starts.end(), at) - starts.begin();
                std::uint64_t const start = at - (k ? starts[k - 1] : 0);
                templ.assign(contigs->seqs[k], start, length);
            }
            if (rng.Uniform() < 0.5) ReverseComplement(templ);
            if (m_qshmm) m_qshmm->Mutate(templ, rng, read, &piece.errors);
            else hifi::Mutate(templ, rng, m_options.setup, read, &piece.errors);
            piece.template_bases += templ.size();
            piece.read_bases += read.seq.size();
            fastq += prefix;
            fastq += std::to_string(item.first_read + i + 1);
            fastq += '\n';
            fastq += read.seq;
            fastq += "\n+\n";
            fastq += read.qual;
            fastq += '\n';
        }
        piece.reads = item.lengths.size();
        piece.bytes = Pack(fastq, s.packing);
        return piece;
    }

    void Complete(Stream& s, std::size_t k, Item const& item, Piece piece, std::unique_lock<std::mutex>& lock) {
        s.result.reads += piece.reads;
        s.result.template_bases += piece.template_bases;
        s.result.read_bases += piece.read_bases;
        s.result.errors += piece.errors;
        s.pieces[k] = std::move(piece);
        ++s.processed;
        if (!s.sample->genomes[item.genome].host) {
            auto slot = s.genomes.find({item.round, item.genome});
            if (slot != s.genomes.end() && --slot->second.parts_left == 0) s.genomes.erase(slot);
        }
        if (s.processed == s.round_end && !s.planned_all) PlanRound(s);
        Drain(s, lock);
        m_cv.notify_all();
    }

    // Writes the sample's pieces that are next in order; one thread at a time per sample, the lock released while it
    // writes. Closes the file once the last is in.
    void Drain(Stream& s, std::unique_lock<std::mutex>& lock) {
        if (s.writing || s.closed) return;
        s.writing = true;
        while (!m_error && s.next_write < s.pieces.size() && s.pieces[s.next_write]) {
            std::string bytes = std::move(s.pieces[s.next_write]->bytes);
            s.pieces[s.next_write].reset();
            ++s.next_write;
            --m_in_flight;
            lock.unlock();
            bool const ok = std::fwrite(bytes.data(), 1, bytes.size(), s.out) == bytes.size();
            lock.lock();
            if (!ok) Fail(std::make_exception_ptr(std::runtime_error("writing " + s.partial.string() + " failed: " +
                                                                     std::strerror(errno))));
        }
        if (!m_error && s.planned_all && s.next_write == s.items.size()) {
            lock.unlock();
            std::exception_ptr error;
            try {
                Close(s);
            } catch (...) {
                error = std::current_exception();
            }
            lock.lock();
            if (error) {
                Fail(error);
            } else {
                s.closed = true;
                ++m_closed;
                m_active.erase(std::find(m_active.begin(), m_active.end(), &s));
            }
        }
        s.writing = false;
    }

    void Close(Stream& s) {
        bool ok = true;
        if (s.packing == Packing::Bgzf) {
            ok = std::fwrite(protal::bgzf::kEof, 1, sizeof(protal::bgzf::kEof), s.out) == sizeof(protal::bgzf::kEof);
        } else if (s.items.empty()) {  // no reads: an empty frame, so that the file is valid zstd
            std::string const empty = Pack(std::string(), s.packing);
            ok = std::fwrite(empty.data(), 1, empty.size(), s.out) == empty.size();
        }
        ok = (std::fclose(s.out) == 0) && ok;
        s.out = nullptr;
        if (!ok) throw std::runtime_error("writing " + s.partial.string() + " failed: " + std::strerror(errno));
        fs::rename(s.partial, s.sample->out);
    }

    std::vector<LongSample> const& m_samples;
    LongReadOptions const& m_options;
    std::unique_ptr<QshmmModel> m_qshmm;
    std::unordered_map<std::string, std::unique_ptr<Host>> m_hosts;
    std::vector<Stream> m_streams;
    std::vector<Stream*> m_active;
    std::size_t m_next_sample = 0, m_closed = 0, m_in_flight = 0, m_cap = 4;
    std::exception_ptr m_error;
    std::mutex m_mutex;
    std::condition_variable m_cv;
};

}  // namespace

std::vector<LongSampleResult> SimulateLongReads(std::vector<LongSample> const& samples, LongReadOptions const& options) {
    if (samples.empty()) return {};
    Engine engine(samples, options);
    return engine.Run();
}

LongSampleResult MutateTemplates(fs::path const& templates, fs::path const& out, LongReadOptions const& options,
                                 std::uint64_t seed) {
    Packing const packing = PackingOf(out);
    std::unique_ptr<QshmmModel> qshmm;
    if (options.setup.method == LongReadSetup::Method::Qshmm) {
        if (options.model.empty()) throw std::runtime_error("a qshmm setup needs its model file (--long_model)");
        qshmm = std::make_unique<QshmmModel>(options.model, options.setup.accuracy, options.setup.sub_ratio,
                                             options.setup.ins_ratio, options.setup.del_ratio);
    }
    protal::ThreadedGzIstream in(templates.c_str());
    if (!in.rdbuf()->is_open()) throw std::runtime_error("cannot read " + templates.string());
    fs::path const partial = out.string() + ".partial";
    std::unique_ptr<std::FILE, int (*)(std::FILE*)> file(std::fopen(partial.c_str(), "wb"), &std::fclose);
    if (!file) throw std::runtime_error("cannot write " + partial.string() + ": " + std::strerror(errno));
    LongSampleResult result;
    result.name = templates.filename().string();
    std::string fastq, name, seq, line;
    MadeRead read;
    auto flush = [&](bool last) {
        if (fastq.empty() && !(last && packing == Packing::Zstd && result.reads == 0)) return;
        std::string const bytes = Pack(fastq, packing);
        if (std::fwrite(bytes.data(), 1, bytes.size(), file.get()) != bytes.size()) {
            throw std::runtime_error("writing " + partial.string() + " failed");
        }
        fastq.clear();
    };
    auto make = [&] {
        if (seq.empty()) return;
        LongRng rng(MixSeed(seed, result.reads));
        if (qshmm) qshmm->Mutate(seq, rng, read, &result.errors);
        else hifi::Mutate(seq, rng, options.setup, read, &result.errors);
        result.template_bases += seq.size();
        result.read_bases += read.seq.size();
        ++result.reads;
        fastq += '@' + name + '\n' + read.seq + "\n+\n" + read.qual + '\n';
        if (fastq.size() >= (4u << 20)) flush(false);
    };
    bool in_record = false;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '>') {
            make();
            auto const end = line.find_first_of(" \t\r", 1);
            name = line.substr(1, end == std::string::npos ? std::string::npos : end - 1);
            seq.clear();
            in_record = true;
        } else if (in_record) {
            for (char c : line) {
                if (!std::isspace(static_cast<unsigned char>(c))) seq += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            }
        }
    }
    if (in.rdbuf()->read_failed()) throw std::runtime_error(templates.string() + " is truncated or corrupt");
    make();
    flush(true);
    if (packing == Packing::Bgzf &&
        std::fwrite(protal::bgzf::kEof, 1, sizeof(protal::bgzf::kEof), file.get()) != sizeof(protal::bgzf::kEof)) {
        throw std::runtime_error("writing " + partial.string() + " failed");
    }
    if (std::fclose(file.release()) != 0) throw std::runtime_error("writing " + partial.string() + " failed");
    fs::rename(partial, out);
    result.rounds = 1;
    return result;
}

// ---- the task files ------------------------------------------------------------------------------------------------

static std::vector<std::unordered_map<std::string, std::string>> ReadTsv(fs::path const& path,
                                                                         std::vector<std::string> const& required) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot read " + path.string());
    std::string line;
    if (!std::getline(in, line)) throw std::runtime_error(path.string() + " is empty");
    if (!line.empty() && line.back() == '\r') line.pop_back();
    auto const header = Split(line, '\t');
    for (auto const& column : required) {
        if (std::find(header.begin(), header.end(), column) == header.end()) {
            throw std::runtime_error(path.string() + ": no column " + column);
        }
    }
    std::vector<std::unordered_map<std::string, std::string>> rows;
    std::size_t line_no = 1;
    while (std::getline(in, line)) {
        ++line_no;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        auto const fields = Split(line, '\t');
        if (fields.size() != header.size()) {
            throw std::runtime_error(path.string() + ", line " + std::to_string(line_no) + ": " +
                                     std::to_string(fields.size()) + " columns, the header has " + std::to_string(header.size()));
        }
        auto& row = rows.emplace_back();
        for (std::size_t i = 0; i < header.size(); ++i) row[header[i]] = fields[i];
    }
    return rows;
}

// samples_tsv: sample, out (its reads: .fq.zst or .fq.gz), bases, seed. genomes_tsv: sample, genome, fasta (the
// host's folder for a host), weight, host (1 for a host); a sample's genomes in their order (g<i>x_ in read names).
std::vector<LongSample> ReadLongSamples(fs::path const& samples_tsv, fs::path const& genomes_tsv) {
    std::vector<LongSample> samples;
    std::unordered_map<std::string, std::size_t> index;
    for (auto& row : ReadTsv(samples_tsv, {"sample", "out", "bases", "seed"})) {
        LongSample sample;
        sample.name = row["sample"];
        sample.out = row["out"];
        try {
            sample.bases = std::stoull(row["bases"]);
            sample.seed = std::stoull(row["seed"]);
        } catch (std::exception const&) {
            throw std::runtime_error(samples_tsv.string() + ": sample " + sample.name + ": bases and seed are whole numbers");
        }
        if (!index.emplace(sample.name, samples.size()).second) {
            throw std::runtime_error(samples_tsv.string() + ": sample " + sample.name + " is listed twice");
        }
        samples.push_back(std::move(sample));
    }
    for (auto& row : ReadTsv(genomes_tsv, {"sample", "genome", "fasta", "weight", "host"})) {
        auto it = index.find(row["sample"]);
        if (it == index.end()) throw std::runtime_error(genomes_tsv.string() + ": sample " + row["sample"] + " is not in " + samples_tsv.string());
        LongGenome genome;
        genome.name = row["genome"];
        genome.fasta = row["fasta"];
        try {
            genome.weight = std::stod(row["weight"]);
        } catch (std::exception const&) {
            throw std::runtime_error(genomes_tsv.string() + ": genome " + genome.name + ": the weight is not a number");
        }
        genome.host = row["host"] == "1";
        samples[it->second].genomes.push_back(std::move(genome));
    }
    return samples;
}

}  // namespace protal::sim
