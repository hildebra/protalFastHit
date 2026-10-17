// What spike_genes::Tally costs per MB of SAM text: the text from stdin, cut into blocks of about 1 MB at line ends (as
// the output handlers hand them over), each block added once on one thread; best of 5 passes.
#include <chrono>
#include <cstdio>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include "Profiling/SpikeGenes.h"

int main() {
    std::string text((std::istreambuf_iterator<char>(std::cin)), std::istreambuf_iterator<char>());
    std::vector<std::string_view> blocks;
    size_t start = 0;
    while (start < text.size()) {
        size_t end = std::min(text.size(), start + (size_t{1} << 20));
        if (end < text.size()) end = text.find('\n', end) + 1;
        if (end == 0) end = text.size();
        blocks.emplace_back(text.data() + start, end - start);
        start = end;
    }
    double best = 1e9;
    size_t records = 0, genes = 0;
    for (int pass = 0; pass < 5; pass++) {
        protal::spike_genes::Tally tally;
        auto const t0 = std::chrono::steady_clock::now();
        for (auto const& block : blocks) tally.Add(block);
        auto const counts = tally.Take();
        double const s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        best = std::min(best, s);
        records = 0;
        for (auto const& [_, n] : counts) records += n;
        genes = counts.size();
    }
    std::printf("%.1f MB in %zu blocks, %zu records counted on %zu genes: %.4f s, %.1f ms/GB, %.1f ns/record\n",
                text.size() / 1e6, blocks.size(), records, genes, best, best * 1e3 / (text.size() / 1e9), best * 1e9 / records);
}
