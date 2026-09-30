// protal's per-gene start-up steps at a given scale: GenomeLoader reading reference.map,
// LoadUniqueKmers, and WriteSamHeader (to a file), each timed.
//   bench_genes reference.fna reference.map unique_kmers.tsv header_out.sam
#include <chrono>
#include <fstream>
#include <iostream>
#include "GenomeLoader.h"

int main(int argc, char** argv) {
    if (argc != 5) { std::cerr << "usage: bench_genes fna map unique_kmers header_out\n"; return 2; }
    using clock = std::chrono::steady_clock;
    auto secs = [](clock::time_point a) { return std::chrono::duration<double>(clock::now() - a).count(); };
    auto t = clock::now();
    protal::GenomeLoader loader(protal::db::DbFile::OnDisk(argv[1]), protal::db::DbFile::OnDisk(argv[2]));
    double const map_s = secs(t);
    t = clock::now();
    loader.LoadUniqueKmers(std::string(argv[3]));
    double const uk_s = secs(t);
    t = clock::now();
    {
        std::ofstream os(argv[4]);
        loader.WriteSamHeader(os);
    }
    double const header_s = secs(t);
    std::cout << "genes " << loader.GeneCount() << "\treference.map " << map_s << " s\tunique_kmers.tsv " << uk_s
              << " s\tSAM header " << header_s << " s" << std::endl;
    return 0;
}
