// readbench - experiment only (multithreading audit, 2026-10-01). Reads a FASTQ pair through protal's
// input path (ThreadedGzIstream per file, SeqReaderPE copies sharing critical(reader)) on T threads that
// do nothing else with the reads, or spin `work_us` microseconds per pair: the most pairs per second the
// input path hands to alignment threads. Built against the instrumented tree (MTAUDIT counters on stderr).
#include "SequenceUtils/SeqReader.h"
#include "IO/ThreadedGzStream.h"
#include <omp.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>

int main(int argc, char** argv) {
    if (argc < 4) { std::fprintf(stderr, "usage: readbench R1 R2 threads [work_us]\n"); return 2; }
    int const threads = std::atoi(argv[3]);
    double const work_us = argc > 4 ? std::atof(argv[4]) : 0;
    protal::ThreadedGzIstream is1(argv[1]), is2(argv[2]);
    protal::SeqReaderPE global(is1, is2);
    size_t pairs = 0, bases = 0;
    auto const t0 = std::chrono::steady_clock::now();
#pragma omp parallel num_threads(threads) reduction(+:pairs,bases)
    {
        FastxRecord r1, r2;
        protal::SeqReaderPE reader(global);
        while (reader(r1, r2)) {
            pairs++;
            bases += r1.sequence.size() + r2.sequence.size();
            if (work_us > 0) {
                auto const until = std::chrono::steady_clock::now() + std::chrono::nanoseconds(static_cast<long>(work_us * 1000));
                while (std::chrono::steady_clock::now() < until) {}
            }
        }
    }
    double const s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("threads\t%d\twork_us\t%g\tpairs\t%zu\tseconds\t%.3f\tpairs_per_s\t%.0f\tMB_per_s\t%.1f\n",
                threads, work_us, pairs, s, pairs / s, bases / s / 1e6);
    return 0;
}
