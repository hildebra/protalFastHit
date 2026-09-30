// protal's gene tables at a given scale and thread count: GenomeLoader reading reference.map and
// LoadUniqueKmers, each timed (the parallel parser of GenomeLoader.h).
//   bench_gene_tables reference.fna reference.map unique_kmers.tsv THREADS
#include <chrono>
#include <cstdlib>
#include <iostream>
#include "GenomeLoader.h"

int main(int argc, char** argv) {
    if (argc != 5) { std::cerr << "usage: bench_gene_tables fna map unique_kmers threads\n"; return 2; }
    int const threads = std::atoi(argv[4]);
    using clock = std::chrono::steady_clock;
    auto secs = [](clock::time_point a) { return std::chrono::duration<double>(clock::now() - a).count(); };
    auto t = clock::now();
    protal::GenomeLoader loader(protal::db::DbFile::OnDisk(argv[1]), protal::db::DbFile::OnDisk(argv[2]), threads);
    double const map_s = secs(t);
    t = clock::now();
    loader.LoadUniqueKmers(std::string(argv[3]), threads);
    double const uk_s = secs(t);
    std::cout << "genes " << loader.GeneCount() << "\tthreads " << threads << "\treference.map " << map_s
              << " s\tunique_kmers.tsv " << uk_s << " s" << std::endl;
    return 0;
}
