// The syncmer scan of protal's read loop alone (SimpleKmerHandler<ClosedSyncmer>, k = 31, m = 15, s = 7, t = 2, as
// RunProtal sets it up), timed on the reads of a FASTQ: build it against two trees and run them alternately.
//   g++ -O3 -std=c++20 -I<tree>/src -I<tree>/src/SequenceUtils ... bench_syncmers.cpp -o bench_syncmers
//   bench_syncmers READS.fq [READS] [ROUNDS]   -> ns per read, the minimum and median over the rounds; a checksum
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "SequenceUtils/KmerIterator.h"

using namespace protal;

int main(int argc, char** argv) {
    if (argc < 2) { std::cerr << "bench_syncmers READS.fq [READS] [ROUNDS]" << std::endl; return 1; }
    size_t const max_reads = argc > 2 ? std::stoul(argv[2]) : 200000;
    int const rounds = argc > 3 ? std::stoi(argv[3]) : 15;
    std::vector<std::string> reads;
    std::ifstream in(argv[1]);
    std::string line;
    for (size_t n = 0; reads.size() < max_reads && std::getline(in, line); n++) {
        if (n % 4 == 1) reads.push_back(line);
    }
    ClosedSyncmer syncmer{15, 7, 2, true};
    SimpleKmerHandler<ClosedSyncmer> handler{31, 15, syncmer};
    KmerList kmers;
    std::vector<double> ns;
    uint64_t checksum = 0;
    for (int r = 0; r < rounds; r++) {
        checksum = 0;
        auto t0 = std::chrono::steady_clock::now();
        for (auto const& read : reads) {
            handler(std::string_view(read), kmers);
            checksum += kmers.size();
            if (!kmers.empty()) checksum ^= kmers.back().first;
        }
        auto t1 = std::chrono::steady_clock::now();
        ns.push_back(std::chrono::duration<double, std::nano>(t1 - t0).count() / static_cast<double>(reads.size()));
    }
    std::sort(ns.begin(), ns.end());
    std::printf("%zu reads, AVX2 %d: min %.1f ns/read, median %.1f; checksum %llu\n", reads.size(), handler.UsesAvx2() ? 1 : 0,
                ns.front(), ns[ns.size() / 2], static_cast<unsigned long long>(checksum));
    return 0;
}
