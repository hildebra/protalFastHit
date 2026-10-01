// What decoding a gene costs against decoding a window of it (SequenceUtils/PackedSequence.h):
//   g++ -O3 -march=x86-64-v3 -I<src>/SequenceUtils bench_window.cpp -o bench_window && ./bench_window
// The genes are packed and kept apart in memory (many genes, so that each decode starts cold, as in an index of
// thousands of genes) and decoded in random order.
#include <chrono>
#include <cstdio>
#include <random>
#include <string>
#include <vector>
#include "PackedSequence.h"

using namespace protal;

int main() {
    std::mt19937_64 rng(1);
    for (size_t length : { 1000u, 2000u, 4000u, 6000u, 20000u }) {
        size_t const genes = std::max<size_t>(64, (256u << 20) / packed::Bytes(length));  // 256 MB of packed genes
        std::vector<uint8_t> store(genes * packed::Bytes(length) + 64);
        std::string s(length, 'A');
        for (auto& c : s) c = "ACGT"[rng() % 4];
        for (size_t g = 0; g < genes; g++) packed::Pack(s.data(), length, store.data() + g * packed::Bytes(length));
        size_t const calls = 200000;
        std::vector<size_t> order(calls), begins(calls);
        for (size_t i = 0; i < calls; i++) { order[i] = rng() % genes; begins[i] = rng() % (length - 200); }
        auto run = [&](size_t width) {
            size_t sink = 0;
            auto t0 = std::chrono::steady_clock::now();
            for (size_t i = 0; i < calls; i++) {
                GeneSequence seq(store.data() + order[i] * packed::Bytes(length), length, begins[i], begins[i] + width);
                sink += static_cast<unsigned char>(seq[begins[i] + width / 2]);
            }
            double ns = std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - t0).count() / calls;
            return std::pair<double, size_t>(ns, sink);
        };
        auto whole = run(length), window = run(200);
        std::printf("gene of %5zu bases: whole %7.1f ns, window of 200 bases %6.1f ns (%.1fx)%s\n", length, whole.first, window.first,
                    whole.first / window.first, whole.second == window.second + 1 ? "" : "");
    }
}
