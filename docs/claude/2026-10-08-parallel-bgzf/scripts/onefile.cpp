// onefile - experiment only: one file through ThreadedGzIstream, taken in 64 KB reads by one thread that does
// nothing else; prints MB/s of inflated bytes and the inflating threads used (PROTAL_INFLATE_THREADS).
#include "IO/ThreadedGzStream.h"
#include <chrono>
#include <cstdio>
int main(int argc, char** argv) {
    if (argc < 2) return 2;
    protal::ThreadedGzIstream is(argv[1]);
    static char buffer[1 << 16];
    size_t bytes = 0;
    auto const t0 = std::chrono::steady_clock::now();
    while (is.read(buffer, sizeof(buffer)) || is.gcount() > 0) bytes += static_cast<size_t>(is.gcount());
    double const s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("threads\t%u\tMB\t%.0f\tseconds\t%.3f\tMB_per_s\t%.0f\tfailed\t%d\n", is.rdbuf()->inflate_threads(), bytes / 1e6, s, bytes / s / 1e6,
                is.rdbuf()->read_failed() ? 1 : 0);
}
