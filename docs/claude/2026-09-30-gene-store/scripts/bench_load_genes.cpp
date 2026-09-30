// Loading a GTDB-sized reference: bench_load_genes DIR GENES THREADS [generate]
// Builds DIR/reference.fna and DIR/reference.map with GENES genes of 900-1199 random bases (once,
// with "generate"), then times GenomeLoader's reference.map parse and LoadAllGenomes with THREADS
// threads, reads every gene once (Sequence(), one base of it), and prints the resident memory. Compiled
// against a checkout's src/, so the same file measures the byte-per-base store and the 2-bit store:
//   g++ -O2 -std=c++20 -fopenmp -I src -I lib bench_load_genes.cpp -o bench -pthread -lzstd ...
#include <chrono>
#include <iostream>
#include <algorithm>
#include <vector>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <random>
#include <string>
#include "SequenceUtils/GenomeLoader.h"
using namespace protal;
static size_t Rss() { std::ifstream s("/proc/self/status"); std::string l; while (std::getline(s, l)) if (l.rfind("VmRSS", 0) == 0) return std::stoull(l.substr(6)) * 1024; return 0; }
static size_t Hwm() { std::ifstream s("/proc/self/status"); std::string l; while (std::getline(s, l)) if (l.rfind("VmHWM", 0) == 0) return std::stoull(l.substr(6)) * 1024; return 0; }
int main(int argc, char** argv) {
    std::string dir = argv[1]; size_t genes = std::stoull(argv[2]); int threads = std::stoi(argv[3]);
    std::string fna = dir + "/reference.fna", map = dir + "/reference.map";
    if (argc > 4) {
        std::mt19937_64 rng(1); std::ofstream f(fna, std::ios::binary), m(map); std::string s; size_t off = 0;
        const char acgt[4] = {'A', 'C', 'G', 'T'};
        for (size_t g = 0; g < genes; g++) {
            size_t taxid = g / 120 + 1, gene = g % 120 + 1, len = 900 + rng() % 300;
            std::string h = ">" + std::to_string(taxid) + "_" + std::to_string(gene) + "\n"; f << h; off += h.size();
            s.resize(len); for (size_t i = 0; i < len; i += 8) { uint64_t r = rng(); for (size_t j = 0; j < 8 && i + j < len; j++) s[i + j] = acgt[(r >> (2 * j)) & 3]; }
            f << s << '\n'; m << taxid << '\t' << gene << '\t' << off << '\t' << off + len << '\n'; off += len + 1;
        }
        return 0;
    }
    auto now = [] { return std::chrono::steady_clock::now(); };
    auto sec = [](auto a, auto b) { return std::chrono::duration<double>(b - a).count(); };
    auto t0 = now();
    GenomeLoader loader(db::DbFile::OnDisk(fna), db::DbFile::OnDisk(map), threads);
    auto t1 = now();
    size_t after_map = Rss();
    loader.LoadAllGenomes(threads);
    auto t2 = now();
    size_t after_load = Rss();
    uint64_t sum = 0, bases = 0;
    for (auto& [id, genome] : loader.GetGenomeMap()) for (auto const& g : genome.GetGeneList()) if (g.IsSet()) { auto const s = g.Sequence(); sum += s[s.size() / 2]; bases += s.size(); }
    auto t3 = now();
    printf("%zu genes, %.2f Gbase, %d threads: reference.map %.2f s (RSS %.2f GB), LoadAllGenomes %.2f s (RSS %.2f GB), read all genes %.2f s (%.0f ns each), peak RSS %.2f GB  [%llu]\n",
           genes, bases / 1e9, threads, sec(t0, t1), after_map / 1e9, sec(t1, t2), after_load / 1e9, sec(t2, t3), sec(t2, t3) * 1e9 / genes, Hwm() / 1e9, (unsigned long long)sum);
}
