// wfabench - experiment only: WFA2's extend kernels, AVX2 against portable, in one process. Aligns the same
// cases with protal's WFA2 wrapper (gap-affine 4/6/2, ends-free, as protal aligns reads), switching the kernels
// (protal_wfa_extend_use_avx2) between rounds, ROUNDS rounds each, alternated, so that both run in the same
// state of the machine. Cases: 150 bp reads against their gene window (0-15% divergence, a third running past
// the gene's end), and long reads (2-10 kb, 2-15% divergence) end to end against their gene. Prints the
// seconds per round of each and the ratio of the medians.
#include <iostream>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>
#include "Alignment/WFA2Wrapper2.h"
#include "wfa2-extend/extend_dispatch.h"

using namespace protal;

namespace {
    std::string Random(size_t n, std::mt19937& rng) {
        std::string s(n, 'A');
        for (auto& c : s) c = "ACGT"[rng() % 4];
        return s;
    }

    std::string Mutate(std::string const& s, double rate, std::mt19937& rng) {
        std::uniform_real_distribution<double> u(0, 1);
        std::string out;
        for (char c : s) {
            double const r = u(rng);
            if (r < rate * 0.8) out += c == 'A' ? 'C' : 'A';
            else if (r < rate * 0.9) {}
            else if (r < rate) { out += c; out += "ACGT"[rng() % 4]; }
            else out += c;
        }
        return out;
    }

    struct Case { std::string read, ref; int q_end_free, r_begin_free, r_end_free, max_score; };
}

int main(int argc, char** argv) {
    int const rounds = argc > 1 ? std::atoi(argv[1]) : 9;
    // With a divergence (second argument), every read has it, else 0-15% (short) and 2-15% (long).
    double const fixed = argc > 2 ? std::atof(argv[2]) : -1;
    std::mt19937 rng(23);
    std::vector<Case> short_reads, long_reads;
    std::uniform_real_distribution<double> rate_any(0.0, 0.15);
    auto rate = [&](std::mt19937& g) { return fixed >= 0 ? fixed : rate_any(g); };
    for (int i = 0; i < 20000; i++) {
        size_t const len = 150, overhang = i % 3 == 0 ? rng() % 40 : 0;
        std::string const gene = Random(len + 18, rng);
        std::string const read = Mutate(gene.substr(9, len - overhang), rate(rng), rng) + Random(overhang, rng);
        std::string const ref = overhang ? gene.substr(0, 9 + len - overhang) : gene;
        short_reads.push_back({ read, ref, overhang ? static_cast<int>(overhang) + 9 : 0, 18, overhang ? 0 : 18, 61 });
    }
    for (int i = 0; i < 300; i++) {
        std::string const gene = Random(2000 + rng() % 8000, rng);
        long_reads.push_back({ Mutate(gene, fixed >= 0 ? fixed : 0.02 + 0.13 * (rng() % 1000) / 1000.0, rng), gene, 0, 0, 0, INT32_MAX });
    }
    WFA2Wrapper2 aligner(4, 6, 2, 0);
    // The thread's own CPU time: what the other processes on the machine take from it is left out.
    auto cpu = [] {
        timespec t{};
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t);
        return static_cast<double>(t.tv_sec) + t.tv_nsec * 1e-9;
    };
    // Cases [from, to) with one kernel: CPU seconds, and the sum of the scores.
    auto run = [&](std::vector<Case> const& cases, size_t from, size_t to) {
        long checksum = 0;
        double const t0 = cpu();
        for (size_t i = from; i < to; i++) {
            auto const& c = cases[i];
            aligner.Reset();
            aligner.Alignment(c.read, c.ref, 0, c.q_end_free, c.r_begin_free, c.r_end_free, c.max_score);
            checksum += aligner.Success() ? aligner.GetAlignmentScore() : 1;
        }
        return std::pair<double, long>{ cpu() - t0, checksum };
    };
    protal_wfa_extend_use_avx2(true);
    std::printf("AVX2 kernels available: %s\n", protal_wfa_extend_uses_avx2() ? "yes" : "no");
    for (auto const* name : { "short", "long" }) {
        bool const is_short = std::string(name) == "short";
        auto const& cases = is_short ? short_reads : long_reads;
        size_t const batch = is_short ? 1000 : 20;
        run(cases, 0, cases.size());  // warm-up: WFA2's buffers
        // Each round: every batch with both kernels in turn, the first kernel alternating; per round the sums.
        std::vector<double> t[2];
        long sums[2] = { 0, 0 };
        for (int r = 0; r < rounds; r++) {
            double total[2] = { 0, 0 };
            for (size_t from = 0, b = 0; from < cases.size(); from += batch, b++) {
                size_t const to = std::min(cases.size(), from + batch);
                for (int j = 0; j < 2; j++) {
                    int const avx2 = (b + r + j) % 2;
                    protal_wfa_extend_use_avx2(avx2 == 1);
                    auto const [s, sum] = run(cases, from, to);
                    total[avx2] += s;
                    if (r == 0) sums[avx2] += sum;
                }
            }
            t[0].push_back(total[0]);
            t[1].push_back(total[1]);
        }
        for (auto& v : t) std::sort(v.begin(), v.end());
        double const portable = t[0][t[0].size() / 2], avx2 = t[1][t[1].size() / 2];
        std::printf("%s reads (%zu, batches of %zu): portable %.3f s, AVX2 %.3f s CPU per round (medians of %d; ranges "
                    "%.3f-%.3f, %.3f-%.3f); AVX2/portable %.3f; scores %s\n", name, cases.size(), batch, portable, avx2, rounds,
                    t[0].front(), t[0].back(), t[1].front(), t[1].back(), avx2 / portable, sums[0] == sums[1] ? "the same" : "DIFFER");
    }
    return 0;
}
