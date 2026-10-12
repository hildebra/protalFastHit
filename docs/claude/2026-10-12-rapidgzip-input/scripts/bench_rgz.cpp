// rapidgzip's ParallelGzipReader on its own: FILE THREADS CHUNK_MIB READ_MIB STATS(0/1)
#include <sys/resource.h>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <vector>

#include <filereader/Standard.hpp>
#include <rapidgzip/ParallelGzipReader.hpp>

static double CpuSeconds() {
    rusage u {};
    getrusage(RUSAGE_SELF, &u);
    return static_cast<double>(u.ru_utime.tv_sec + u.ru_stime.tv_sec) + 1e-6 * static_cast<double>(u.ru_utime.tv_usec + u.ru_stime.tv_usec);
}

int main(int argc, char** argv) {
    if (argc < 6) return 2;
    unsigned const threads = static_cast<unsigned>(std::atoi(argv[2]));
    size_t const chunk = static_cast<size_t>(std::atof(argv[3]) * (1 << 20));
    size_t const readsize = static_cast<size_t>(std::atof(argv[4]) * (1 << 20));
    bool const stats = std::atoi(argv[5]) != 0;
    auto const t0 = std::chrono::steady_clock::now();
    double const c0 = CpuSeconds();
    size_t total = 0;
    {
        rapidgzip::ParallelGzipReader<> reader(std::make_unique<rapidgzip::StandardFileReader>(std::string(argv[1])), threads, chunk);
        reader.setCRC32Enabled(true);
        reader.setKeepIndex(false);
        if (stats) {
            reader.setStatisticsEnabled(true);
            reader.setShowProfileOnDestruction(true);
        }
        std::vector<char> buffer(readsize);
        while (true) {
            size_t const n = reader.read(buffer.data(), buffer.size());
            if (n == 0) break;
            total += n;
        }
    }
    double const wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    rusage u {};
    getrusage(RUSAGE_SELF, &u);
    std::printf("threads %u chunk %.1f MiB read %.1f MiB: %zu bytes, wall %.3f s, cpu %.3f s, %.3f GB/s, maxrss %.0f MB\n", threads,
                static_cast<double>(chunk) / (1 << 20), static_cast<double>(readsize) / (1 << 20), total, wall, CpuSeconds() - c0,
                static_cast<double>(total) / wall / 1e9, static_cast<double>(u.ru_maxrss) / 1024.0);
    return 0;
}
