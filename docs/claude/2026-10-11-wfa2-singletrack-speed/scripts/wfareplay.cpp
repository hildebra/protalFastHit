// wfareplay.cpp: replays the WFA2 calls of protal runs (dumped by the scratch WFA2Wrapper2: kind A = ends-free, E = end to
// end; x_drop, free ends, max_score, read, reference) through protal's WFA2Wrapper2, as protal makes them: one aligner per
// kind and X-drop, reused. Variants, chosen at run time and alternated round by round so that each sees the same machine:
//   classic | singletrack        the traceback
//   classic/2 | singletrack/2    penalties halved (2,3,1 for 4,6,2): every score is even, so WFA2's odd score steps are
//                                empty; the budget becomes ceil(max/2) and the score is doubled back
// Prints each variant's thread CPU seconds per round (median, min) and writes per call "status score cigar-hash" to
// <out>.<variant> for the identity checks.
//   wfareplay CALLS ROUNDS OUT VARIANT...
#include <algorithm>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include "WFA2Wrapper2.h"

using namespace protal;

struct Call {
    char kind;
    size_t x_drop;
    int qb, qe, rb, re, max_score;
    std::string query, ref;
};

struct Variant {
    std::string name;
    WfaTraceback traceback;
    bool halve;
};

static double ThreadSeconds() {
    timespec t;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t);
    return static_cast<double>(t.tv_sec) + 1e-9 * static_cast<double>(t.tv_nsec);
}

static uint64_t Fnv(std::string const& s) {
    uint64_t h = 1469598103934665603ULL;
    for (unsigned char c : s) { h ^= c; h *= 1099511628211ULL; }
    return h;
}

int main(int argc, char** argv) {
    if (argc < 5) { std::cerr << "usage: wfareplay CALLS ROUNDS OUT VARIANT...\n"; return 2; }
    std::vector<Call> calls;
    {
        std::ifstream in(argv[1]);
        std::string line;
        while (std::getline(in, line)) {
            std::istringstream is(line);
            Call c;
            std::string kind;
            is >> kind >> c.x_drop >> c.qb >> c.qe >> c.rb >> c.re >> c.max_score >> c.query >> c.ref;
            if (kind.empty()) continue;
            c.kind = kind[0];
            calls.push_back(std::move(c));
        }
    }
    int const rounds = std::atoi(argv[2]);
    std::string const out = argv[3];
    std::vector<Variant> variants;
    for (int i = 4; i < argc; i++) {
        std::string const v = argv[i];
        bool const halve = v.size() > 2 && v.substr(v.size() - 2) == "/2";
        std::string const base = halve ? v.substr(0, v.size() - 2) : v;
        variants.push_back({ v, base == "singletrack" ? WfaTraceback::Singletrack : WfaTraceback::Classic, halve });
    }
    std::cerr << calls.size() << " calls\n";
    std::vector<std::vector<double>> seconds(variants.size());
    for (int r = 0; r < rounds; r++) {
        for (size_t vi = 0; vi < variants.size(); vi++) {
            Variant const& v = variants[vi];
            int const x = v.halve ? 2 : 4, o = v.halve ? 3 : 6, e = v.halve ? 1 : 2;
            std::map<std::pair<char, size_t>, std::unique_ptr<WFA2Wrapper2>> aligners;
            std::vector<std::string> results;
            if (r == 0) results.reserve(calls.size());
            double const start = ThreadSeconds();
            for (auto const& c : calls) {
                auto& aligner = aligners[{ c.kind, c.x_drop }];
                if (!aligner) aligner = std::make_unique<WFA2Wrapper2>(x, o, e, c.x_drop, v.traceback);
                int const budget = v.halve ? (c.max_score == INT32_MAX ? INT32_MAX : (c.max_score + 1) / 2) : c.max_score;
                aligner->Reset();
                if (c.kind == 'E') aligner->EndToEnd(c.query, c.ref, budget);
                else aligner->Alignment(c.query, c.ref, c.qb, c.qe, c.rb, c.re, budget);
                if (r == 0) {
                    bool const ok = aligner->Success();
                    int const score = ok ? -aligner->GetAlignmentScore() * (v.halve ? 2 : 1) : 0;
                    char buf[64];
                    std::snprintf(buf, sizeof buf, "%d %d %016llx", ok ? 1 : 0, score,
                                  ok ? static_cast<unsigned long long>(Fnv(aligner->Cigar())) : 0ULL);
                    results.emplace_back(buf);
                }
            }
            seconds[vi].push_back(ThreadSeconds() - start);
            if (r == 0) {
                std::string name = v.halve ? v.name.substr(0, v.name.size() - 2) + "_half" : v.name;
                std::ofstream f(out + "." + name);
                for (auto const& s : results) f << s << '\n';
            }
        }
    }
    for (size_t vi = 0; vi < variants.size(); vi++) {
        auto s = seconds[vi];
        std::sort(s.begin(), s.end());
        std::printf("%-16s median %.3f s, min %.3f s (%d rounds)\n", variants[vi].name.c_str(), s[s.size() / 2], s[0], rounds);
    }
    std::printf("fallbacks to classic: %zu\n", WFA2Wrapper2::SingletrackFallbacks());
}
