// wfadiff.cpp: the same protal-like alignments through whichever WFA2-lib it is built against, one line per alignment
// (kind, case, status, score, whether the operations are an alignment of the sequences, the operations), for a diff of
// two libraries' classic tracebacks (and, with "singletrack" as argv[1], the fork's Singletrack against classic).
// Gap-affine 4/6/2 with WFA2's default wf-adaptive heuristic, as protal aligns; budgets as protal sets them.
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <random>
#include <string>
#include "wfa2-lib/bindings/cpp/WFAligner.hpp"

static std::string Random(size_t n, std::mt19937& rng) {
    std::string s(n, 'A');
    for (auto& c : s) c = "ACGT"[rng() % 4];
    return s;
}

static std::string Mutate(std::string const& s, double rate, std::mt19937& rng) {
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

static bool Valid(std::string const& ops, std::string const& ref, std::string const& read) {
    size_t r = 0, q = 0;
    for (char const op : ops) {
        if (op == 'M' || op == 'X') {
            if (r >= ref.size() || q >= read.size() || (ref[r] == read[q]) != (op == 'M')) return false;
            r++, q++;
        } else if (op == 'D') r++;
        else if (op == 'I') q++;
        else return false;
    }
    return r == ref.size() && q == read.size();
}

int main(int argc, char** argv) {
    bool const single = argc > 1 && std::strcmp(argv[1], "singletrack") == 0;
#ifdef SINGLETRACK
    auto const memory = single ? wfa::WFAligner::MemorySingletrack : wfa::WFAligner::MemoryHigh;
#else
    auto const memory = wfa::WFAligner::MemoryHigh;
    if (single) return 2;
#endif
    std::mt19937 rng(20261011);
    std::uniform_real_distribution<double> rate(0.0, 0.25);
    for (int kind = 0; kind < 4; kind++) {
        size_t const n = kind == 0 ? 20000 : kind == 1 ? 1500 : 6000;
        for (size_t i = 0; i < n; i++) {
            std::string ref, read;
            int qe = 0, rb = 0, re = 0, budget = 0;
            size_t x_drop = 0;
            bool e2e = false;
            if (kind == 0) {  // short reads in gene windows with dovetails (Cases in test_WFA2Wrapper.cpp), X-drop 1000
                size_t const overhang = i % 3 == 0 ? rng() % 40 : 0;
                std::string const gene = Random(168, rng);
                read = Mutate(gene.substr(9, 150 - overhang), rate(rng) * 0.6, rng) + Random(overhang, rng);
                ref = overhang ? gene.substr(0, 159 - overhang) : gene;
                qe = overhang ? static_cast<int>(overhang) + 9 : 0, rb = 18, re = overhang ? 0 : 18, budget = 61, x_drop = 1000;
            } else if (kind == 1) {  // long reads in gene windows, a fifth unrelated
                size_t const len = 1000 + rng() % 2000;
                ref = Random(len + 320, rng);
                read = i % 5 == 4 ? Random(len, rng) : Mutate(ref.substr(160, len), rate(rng) * 0.6, rng);
                rb = re = 200, budget = static_cast<int>(std::ceil(0.15 * static_cast<double>(len))) * 4 + 1;
            } else if (kind == 2) {  // long reads' flanks: the read's end free, the reference's start not, a quarter unrelated
                ref = Random(150 + rng() % 450, rng);
                read = i % 4 == 3 ? Random(200 + rng() % 400, rng) : Mutate(ref, rate(rng), rng) + Random(rng() % 300, rng);
                qe = static_cast<int>(read.size()), re = rng() % 2 ? 0 : static_cast<int>(ref.size() / 4);
                budget = i % 2 ? 2000 : static_cast<int>(std::ceil(0.15 * static_cast<double>(read.size()))) * 4 + 1;
            } else {  // pieces between a long read's links, end to end
                ref = Random(20 + rng() % 400, rng);
                read = Mutate(ref, rate(rng) * 0.6, rng);
                budget = 400, e2e = true;
            }
            wfa::WFAlignerGapAffine aligner(4, 6, 2, wfa::WFAligner::Alignment, memory);
            if (x_drop) aligner.setHeuristicXDrop(static_cast<int>(x_drop), 1);
            aligner.setMaxAlignmentSteps(budget);
            int const status = e2e ? aligner.alignEnd2End(ref, read) : aligner.alignEndsFree(ref, rb, re, read, 0, qe);
            std::string const ops = status == wfa::WFAligner::StatusAlgCompleted ? aligner.getAlignment() : "";
            std::printf("%d\t%zu\t%d\t%d\t%d\t%s\n", kind, i, status, status == wfa::WFAligner::StatusAlgCompleted ? -aligner.getAlignmentScore() : 0,
                        ops.empty() ? -1 : static_cast<int>(Valid(ops, ref, read)), ops.c_str());
        }
    }
}
