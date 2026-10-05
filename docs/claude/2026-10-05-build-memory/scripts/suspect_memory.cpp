// suspect_memory.cpp - memory of the suspect-copy scan of protal --build (Build.h WriteSuspectCopies:
// every species' copy of each gene sketched, the sketches gathered by gene, gene_incongruence::Scan) at the
// scale of GTDB r226's training database (13.63M gene copies of 168 genes, v10 log).
//
//   suspect_memory THREADS [SPECIES=81130] [GENES=168] [copy|move]
//
// The sketches are random (kSketchSize sorted hashes, no two copies alike), so Scan finds no near pairs and
// is quick; what it holds per gene (the sorted hash table of all sketches) is the same as for real genes.
// The collection loop is Build.h's, with the sketch made from random hashes instead of the gene; `copy` merges
// the threads' sketches as 18dceda does (copies), `move` as the change does. Prints the resident memory at each
// step and after the copies are freed (what malloc keeps), then after malloc_trim.
#include "GeneIncongruence.h"
#include "rss.h"

#include <malloc.h>
#include <omp.h>

#include <chrono>
#include <cstdio>
#include <random>
#include <vector>

using namespace protal;

int main(int argc, char** argv) {
    int const threads = argc > 1 ? std::atoi(argv[1]) : 6;
    size_t const species = argc > 2 ? std::stoul(argv[2]) : 81130;
    size_t const genes = argc > 3 ? std::stoul(argv[3]) : 168;
    bool const move = argc > 4 && std::string(argv[4]) == "move";
    // Taxonomy: species 1..S, five per genus, genera below one family each, one order, class, phylum, domain.
    std::unordered_map<uint32_t, gene_incongruence::Lineage> lineages;
    for (uint32_t t = 1; t <= species; t++) {
        gene_incongruence::Lineage l;
        uint32_t const genus = 1000000 + (t - 1) / 5;
        l.ancestor = { genus, 2000000 + (genus - 1000000) / 20, 3000000, 3000001, 3000002, 3000003 };
        lineages[t] = l;
    }
    std::vector<uint32_t> taxids(species);
    for (size_t t = 0; t < species; t++) taxids[t] = static_cast<uint32_t>(t + 1);
    double const base = RssMb();
    std::printf("%zu species x %zu genes = %zu copies, %d threads, sketches %s; resident %.0f MB\n", species, genes, species * genes,
                threads, move ? "moved" : "copied", base);

    auto t0 = std::chrono::steady_clock::now();
    PeakSampler collect_peak;
    std::vector<std::vector<gene_incongruence::Copy>> copies;
    {
        #pragma omp parallel num_threads(threads)
        {
            std::vector<std::vector<gene_incongruence::Copy>> mine;
            std::mt19937 rng(omp_get_thread_num() + 1);
            std::vector<uint32_t> sketch;
            #pragma omp for schedule(dynamic, 16)
            for (size_t t = 0; t < taxids.size(); t++) {
                for (size_t i = 0; i < genes; i++) {
                    sketch.resize(gene_incongruence::kSketchSize);
                    for (auto& h : sketch) h = rng();
                    std::sort(sketch.begin(), sketch.end());
                    if (mine.size() <= i + 1) mine.resize(i + 2);
                    mine[i + 1].push_back({ taxids[t], std::vector<uint32_t>(sketch) });
                }
            }
            #pragma omp critical(suspect_copies_merge)
            {
                if (copies.size() < mine.size()) copies.resize(mine.size());
                for (size_t g = 0; g < mine.size(); g++) {
                    if (move) copies[g].insert(copies[g].end(), std::make_move_iterator(mine[g].begin()), std::make_move_iterator(mine[g].end()));
                    else copies[g].insert(copies[g].end(), mine[g].begin(), mine[g].end());
                }
            }
        }
        for (auto& gene : copies) {
            std::sort(gene.begin(), gene.end(), [](auto const& a, auto const& b) { return a.taxid < b.taxid; });
        }
    }
    double const c = collect_peak.Stop();
    double const held = RssMb();
    auto t1 = std::chrono::steady_clock::now();
    PeakSampler scan_peak;
    auto const result = gene_incongruence::Scan(copies, lineages, 0.02, threads);
    double const s = scan_peak.Stop();
    auto t2 = std::chrono::steady_clock::now();
    copies = {};
    double const freed = RssMb();
    malloc_trim(0);
    double const trimmed = RssMb();
    auto secs = [](auto a, auto b) { return std::chrono::duration<double>(b - a).count(); };
    std::printf("collect %6.1f s  peak +%6.0f MB, the sketches then hold +%.0f MB\n", secs(t0, t1), c - base, held - base);
    std::printf("scan    %6.1f s  peak +%6.0f MB (%zu candidates, %zu pairs)\n", secs(t1, t2), s - base, result.candidates, result.pairs.size());
    std::printf("freed: +%.0f MB still resident; after malloc_trim +%.0f MB; VmHWM %.0f MB\n", freed - base, trimmed - base, HwmMb());
    return 0;
}
