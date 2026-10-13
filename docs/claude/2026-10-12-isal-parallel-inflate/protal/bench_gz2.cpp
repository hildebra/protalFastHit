// Reads a file through ThreadedGzIstream with N inflate threads; prints wall time, CPU time and GB/s of output.
// usage: bench_gz FILE THREADS
#include <sys/resource.h>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "IO/ThreadedGzStream.h"

static double CpuSeconds() {
    rusage u {};
    getrusage(RUSAGE_SELF, &u);
    return static_cast<double>(u.ru_utime.tv_sec + u.ru_stime.tv_sec) + 1e-6 * static_cast<double>(u.ru_utime.tv_usec + u.ru_stime.tv_usec);
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: bench_gz FILE THREADS\n");
        return 2;
    }
    std::string const path = argv[1];
    unsigned const threads = static_cast<unsigned>(std::atoi(argv[2]));
    protal::ThreadedGzStreambuf::SetInflateThreads(threads);
    protal::ThreadedGzStreambuf::SetGzipThreads(threads);
    if (char const* p = std::getenv("PGZ_PIECE")) protal::ParallelGzip::SetPieceBytes(static_cast<size_t>(std::atol(p)));
    auto const t0 = std::chrono::steady_clock::now();
    double const c0 = CpuSeconds();
    protal::ThreadedGzIstream is(path.c_str());
    std::vector<char> buffer(1 << 20);
    size_t total = 0;
    unsigned long long check = 0;
    while (true) {
        std::streamsize const n = is.rdbuf()->sgetn(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        if (n <= 0) break;
        total += static_cast<size_t>(n);
        for (std::streamsize i = 0; i < n; i += 4096) check += static_cast<unsigned char>(buffer[static_cast<size_t>(i)]);
    }
    double const wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    double const cpu = CpuSeconds() - c0;
    long rss = 0;
    {
        rusage u {};
        getrusage(RUSAGE_SELF, &u);
        rss = u.ru_maxrss;
    }
    std::printf("%s\tthreads %u\tparallel %d\tbytes %zu\twall %.3f\tcpu %.3f\tGB/s %.3f\tmaxrss_mb %.0f\tfailed %d\tcheck %llu\t%s\n",
                path.c_str(), is.rdbuf()->inflate_threads(), is.rdbuf()->parallel_gzip() ? 1 : 0, total, wall, cpu,
                static_cast<double>(total) / wall / 1e9, static_cast<double>(rss) / 1024.0, is.rdbuf()->read_failed() ? 1 : 0,
                check, is.rdbuf()->read_error_message().c_str());
    return is.rdbuf()->read_failed() ? 1 : 0;
}
