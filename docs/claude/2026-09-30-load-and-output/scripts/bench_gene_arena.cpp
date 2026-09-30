// GenomeLoader::LoadAllGenomes at a given scale and thread count: wall time and the CPU it takes
// (the serial sizing of the gene sequences, then the parallel fill), and the peak RSS.
//   bench_gene_arena reference.fna reference.map THREADS
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <sys/resource.h>
#include "GenomeLoader.h"

static double Seconds(timeval const& t) { return static_cast<double>(t.tv_sec) + t.tv_usec / 1e6; }

int main(int argc, char** argv) {
    if (argc != 4) { std::cerr << "usage: bench_gene_arena fna map threads\n"; return 2; }
    int const threads = std::atoi(argv[3]);
    protal::GenomeLoader loader(protal::db::DbFile::OnDisk(argv[1]), protal::db::DbFile::OnDisk(argv[2]), threads);
    rusage before{}, after{};
    getrusage(RUSAGE_SELF, &before);
    auto const t = std::chrono::steady_clock::now();
    loader.LoadAllGenomes(threads);
    double const wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count();
    getrusage(RUSAGE_SELF, &after);
    std::cout << "genes " << loader.GeneCount() << "\tthreads " << threads << "\tLoadAllGenomes " << wall << " s wall, "
              << Seconds(after.ru_utime) - Seconds(before.ru_utime) << " s user, "
              << Seconds(after.ru_stime) - Seconds(before.ru_stime) << " s sys, peak RSS " << after.ru_maxrss / 1024
              << " MB" << std::endl;
    return 0;
}
