// genome_bench.cpp - what an anchor's gene costs to reach at GTDB r226 size: 143,614 genomes of 101 genes in a
// tsl::sparse_map<taxid, Genome> as GenomeLoader holds them (the real Genome and Gene classes), looked up as
// GenomeLoader::GetGenome does (contains, then at: two probes) and Genome::GetGeneOMP (the loaded flag, the gene
// vector), against one find(), a flat taxid -> Genome* table, and a flat taxid -> Gene* table. Each lookup reads the
// gene's length and its packed bytes' first word (the sequence an extension reads), from a 2 GB arena standing in for
// the 4 GB of packed genes. Independent lookups (as the anchors of a read are) and a dependent chain (each lookup's
// result picks the next), both in ns per lookup.
//
// Build: cc.sh genome_bench.cpp genome_bench; run: genome_bench <lookups> <rounds>
#include <iostream>
#include <chrono>
#include <cstdio>
#include <random>
#include <vector>
#include <sys/mman.h>

#include "SequenceUtils/GenomeLoader.h"

using protal::Gene;
using protal::Genome;

int main(int argc, char** argv) {
    size_t const lookups = argc > 1 ? std::strtoull(argv[1], nullptr, 10) : 5000000;
    int const rounds = argc > 2 ? std::atoi(argv[2]) : 3;
    constexpr size_t kGenomes = 143614, kGenes = 101, kGeneLength = 1100;
    size_t const arena_bytes = size_t{2} << 30;
    uint8_t* arena = static_cast<uint8_t*>(mmap(nullptr, arena_bytes, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0));
    madvise(arena, arena_bytes, MADV_HUGEPAGE);
    std::mt19937_64 rng(3);

    // Taxids as GTDB's internal ones: up to 2^18, not all used.
    std::vector<size_t> taxids;
    for (size_t t = 1; taxids.size() < kGenomes; t++) if (rng() % 100 < 60) taxids.push_back(t);
    tsl::sparse_map<size_t, Genome> map;
    for (size_t t : taxids) map.insert({ t, Genome(t) });
    size_t const max_taxid = taxids.back();
    std::vector<Genome*> by_taxid(max_taxid + 1, nullptr);
    std::vector<Gene*> genes_by_taxid(max_taxid + 1, nullptr);
    size_t offset = 0;
    for (size_t t : taxids) {
        auto& genome = map.find(t).value();
        for (size_t g = 1; g <= kGenes; g++) genome.AddGene(g, g, 0, kGeneLength, nullptr);
        for (auto& gene : genome.Genes()) {
            gene.SetPacked(arena + (offset % (arena_bytes - 4096)));
            offset += protal::packed::Bytes(kGeneLength) + 7;
        }
        genome.MarkLoaded();
        by_taxid[t] = &genome;
        genes_by_taxid[t] = genome.Genes().data();
    }
    std::memset(arena, 1, arena_bytes);  // touched: no first-touch faults in the timing
    std::printf("%zu genomes (sizeof(Genome) %zu, sizeof(Gene) %zu), max taxid %zu\n", map.size(), sizeof(Genome), sizeof(Gene), max_taxid);

    std::vector<std::pair<uint32_t, uint32_t>> pairs(lookups);
    for (auto& p : pairs) p = { static_cast<uint32_t>(taxids[rng() % taxids.size()]), static_cast<uint32_t>(1 + rng() % kGenes) };

    auto first_word = [](Gene& gene) {
        uint64_t w;
        std::memcpy(&w, gene.MutablePacked(), 8);
        return w + gene.GetLength();
    };
    auto as_loader = [&](uint32_t t, uint32_t g) -> Gene& {  // GenomeLoader::GetGenome + Genome::GetGeneOMP
        if (!map.contains(t)) std::exit(10);
        return map.at(t).GetGeneOMP(g);
    };
    auto one_find = [&](uint32_t t, uint32_t g) -> Gene& {
        auto it = map.find(t);
        if (it == map.end()) std::exit(10);
        return it.value().GetGeneOMP(g);
    };
    auto flat_genome = [&](uint32_t t, uint32_t g) -> Gene& { return by_taxid[t]->GetGeneOMP(g); };
    auto flat_gene = [&](uint32_t t, uint32_t g) -> Gene& { return genes_by_taxid[t][g - 1]; };

    auto time = [&](char const* name, auto&& get) {
        uint64_t sink = 0;
        auto start = std::chrono::steady_clock::now();
        for (auto const& [t, g] : pairs) sink += first_word(get(t, g));
        double const independent = std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - start).count() / lookups;
        // Dependent: the next pair's index depends on the last word read.
        start = std::chrono::steady_clock::now();
        size_t i = 0;
        for (size_t n = 0; n < lookups; n++) {
            auto const& [t, g] = pairs[i];
            uint64_t const w = first_word(get(t, g));
            sink += w;
            i = (i + 1 + (w & 1)) % lookups;
        }
        double const dependent = std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - start).count() / lookups;
        std::printf("%-34s independent %6.1f ns, dependent %6.1f ns (%llu)\n", name, independent, dependent, static_cast<unsigned long long>(sink & 1));
    };
    for (int r = 0; r < rounds; r++) {
        time("GetGenome + GetGeneOMP (2 probes)", as_loader);
        time("one find + GetGeneOMP", one_find);
        time("flat taxid -> Genome*", flat_genome);
        time("flat taxid -> Gene*", flat_gene);
    }
    return 0;
}
