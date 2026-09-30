// The read pairs protal's reader delivers when threads do nothing else with them: SeqReaderPE on
// two ThreadedGzIstreams (protal's read input), THREADS threads taking batches under the reader
// lock. With enough threads the lock (or inflating) is the limit.
//   bench_reader R1 R2 THREADS
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iostream>  // SeqReader.h up to 94f6a12 uses std::cerr without including it
#include <omp.h>
#include "IO/ThreadedGzStream.h"
#include "SequenceUtils/SeqReader.h"

using namespace protal;

int main(int argc, char** argv) {
    if (argc != 4) { std::fprintf(stderr, "usage: bench_reader R1 R2 THREADS\n"); return 2; }
    int const threads = std::atoi(argv[3]);
    auto const start = std::chrono::steady_clock::now();
    ThreadedGzIstream is1(argv[1]), is2(argv[2]);
    SeqReaderPE global{ is1, is2 };
    size_t pairs = 0, bases = 0;
#pragma omp parallel num_threads(threads) reduction(+:pairs, bases)
    {
        SeqReaderPE reader{ global };
        FastxRecord a, b;
        while (reader(a, b)) {
            pairs++;
            bases += a.sequence.size() + b.sequence.size() + a.quality.size() + b.id.size();
        }
    }
    double const s = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::printf("threads %d\tpairs %zu\tbases %zu\t%.2f s\t%.0fk pairs/s\n", threads, pairs, bases, s, pairs / s / 1000);
    return 0;
}
