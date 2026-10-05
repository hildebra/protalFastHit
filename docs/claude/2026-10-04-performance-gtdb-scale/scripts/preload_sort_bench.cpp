// The preload's sort at r226 size: 143,614 genomes x ~101 of 120 gene ids, reference written gene by gene (taxids
// ascending in each gene's block), genes listed genome by genome. std::sort of (start, gene) pairs, of gene pointers
// (the old way), and __gnu_parallel::sort of the pairs.
#include <parallel/algorithm>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <random>
#include <vector>
#include <omp.h>
struct Gene { uint64_t where; uint32_t id, length, a, b, c, d; };  // 32 bytes, as protal's
int main() {
    const int genomes = 143614, ids = 120;
    std::mt19937 rng(7);
    std::vector<std::vector<Gene>> by_genome(genomes, std::vector<Gene>(ids));
    std::vector<uint64_t> block(ids + 1, 0);
    std::vector<std::vector<char>> has(genomes, std::vector<char>(ids));
    for (int g = 0; g < genomes; g++) for (int i = 0; i < ids; i++) has[g][i] = rng() % 100 < 84;
    uint64_t offset = 0;
    for (int i = 0; i < ids; i++) for (int g = 0; g < genomes; g++) if (has[g][i]) {
        auto& gene = by_genome[g][i]; gene.where = offset; gene.length = 900 + rng() % 400; gene.id = i + 1; offset += gene.length + 12; }
    std::vector<std::pair<uint64_t, Gene*>> pairs; std::vector<Gene*> pointers;
    for (int g = 0; g < genomes; g++) for (int i = 0; i < ids; i++) if (has[g][i]) { pairs.emplace_back(by_genome[g][i].where, &by_genome[g][i]); pointers.push_back(&by_genome[g][i]); }
    printf("%zu genes, %d threads\n", pairs.size(), omp_get_max_threads());
    auto time = [](const char* what, auto&& f) { auto t = std::chrono::steady_clock::now(); f();
        printf("%-34s %.3f s\n", what, std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count()); };
    auto by = [](auto const& a, auto const& b) { return a.first < b.first; };
    { auto p = pointers; time("gene pointers, std::sort", [&] { std::sort(p.begin(), p.end(), [](Gene* a, Gene* b) { return a->where < b->where; }); }); }
    { auto p = pairs; time("pairs, std::sort", [&] { std::sort(p.begin(), p.end(), by); }); }
    { auto p = pairs; time("pairs, __gnu_parallel::sort", [&] { __gnu_parallel::sort(p.begin(), p.end(), by); }); }
    { auto p = pairs; time("pairs, is_sorted (no sort)", [&] { volatile bool s = std::is_sorted(p.begin(), p.end(), by); (void)s; }); }
}
