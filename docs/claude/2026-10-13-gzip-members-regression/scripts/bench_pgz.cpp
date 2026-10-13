// Reads a file through ThreadedGzIstream with THREADS gzip threads (1: ISA-L on one thread, more: ParallelGzip's
// pieces, or with PGZ_NEW a gzip file of several members by members); prints wall and CPU seconds, GB/s of output, max
// RSS and how the file was read (ParallelGzip::Counts; with PGZ_NEW whether by members and the restarts on one thread).
// usage: bench_pgz FILE THREADS
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
        std::fprintf(stderr, "usage: bench_pgz FILE THREADS\n");
        return 2;
    }
    unsigned const threads = static_cast<unsigned>(std::atoi(argv[2]));
    protal::ThreadedGzStreambuf::SetInflateThreads(threads);
    protal::ThreadedGzStreambuf::SetGzipThreads(threads);
    auto const t0 = std::chrono::steady_clock::now();
    double const c0 = CpuSeconds();
    protal::ThreadedGzIstream is(argv[1]);
    std::vector<char> buffer(1 << 20);
    size_t total = 0;
    while (true) {
        std::streamsize const n = is.rdbuf()->sgetn(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        if (n <= 0) break;
        total += static_cast<size_t>(n);
    }
    double const wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    double const cpu = CpuSeconds() - c0;
    rusage u {};
    getrusage(RUSAGE_SELF, &u);
    auto const c = is.rdbuf()->parallel_gzip_counts();
#ifdef PGZ_NEW
    int const members = is.rdbuf()->opened_as_gzip_members() ? 1 : 0;
    size_t const restarts = is.rdbuf()->gzip_member_restarts();
#else
    int const members = 0;
    size_t const restarts = 0;
#endif
    std::printf("threads %u\tparallel %d\tmembers %d\trestarts %zu\tbytes %zu\twall %.3f\tcpu %.3f\tGB/s %.3f\tmaxrss_mb %.0f\t"
                "ahead %zu\tin_order %zu\tstreamed %zu\tfailed %d %s\n",
                threads, is.rdbuf()->parallel_gzip() ? 1 : 0, members, restarts, total, wall, cpu,
                static_cast<double>(total) / wall / 1e9, static_cast<double>(u.ru_maxrss) / 1024.0, c.ahead, c.in_order,
                c.streamed, is.rdbuf()->read_failed() ? 1 : 0, is.rdbuf()->read_error_message().c_str());
    return is.rdbuf()->read_failed() ? 1 : 0;
}
