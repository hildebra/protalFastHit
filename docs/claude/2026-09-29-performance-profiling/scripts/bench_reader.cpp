// Read throughput of SeqReaderPE alone: THREADS threads take pairs and do nothing else.
//   bench_reader {gz|threaded} R1 R2 THREADS
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <omp.h>
#include "gzstream.h"
#include "IO/ThreadedGzStream.h"
#include "SequenceUtils/SeqReader.h"

using namespace protal;

template<typename Stream>
static void Run(char const* r1, char const* r2, int threads, char const* label) {
    auto const start = std::chrono::steady_clock::now();
    Stream is1(r1), is2(r2);
    SeqReaderPE global{ is1, is2 };
    size_t pairs = 0, bases = 0;
#pragma omp parallel num_threads(threads) reduction(+:pairs, bases)
    {
        SeqReaderPE reader{ global };
        FastxRecord a, b;
        while (reader(a, b)) {
            pairs++;
            bases += a.sequence.size() + b.sequence.size();
        }
    }
    double const s = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::printf("%s\tthreads=%d\tpairs=%zu\tbases=%zu\t%.2f s\t%.0f pairs/s\n", label, threads, pairs, bases, s, pairs / s);
}

int main(int argc, char** argv) {
    if (argc != 5) { std::fprintf(stderr, "bench_reader {gz|threaded} R1 R2 THREADS\n"); return 2; }
    int const threads = std::atoi(argv[4]);
    if (std::strcmp(argv[1], "gz") == 0) Run<igzstream>(argv[2], argv[3], threads, "igzstream");
    else Run<ThreadedGzIstream>(argv[2], argv[3], threads, "threaded");
    return 0;
}
