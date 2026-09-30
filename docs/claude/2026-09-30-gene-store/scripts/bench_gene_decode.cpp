// Cost of decoding a gene from the 2-bit store (src/SequenceUtils/PackedSequence.h), the table and the
// AVX2 path, against copying the same bases from a byte-per-base arena (what Sequence() cost before):
//   g++ -O2 -I src docs/claude/2026-09-30-gene-store/scripts/bench_gene_decode.cpp -o bench_gene_decode -pthread
// 16.6M-gene-like store scaled to 1 Gbase; genes of the given length at random byte-aligned places;
// ns per gene, and GB/s of bases produced. Also packing speed (the load path): 1 Gbase table vs AVX2.
#include <chrono>
#include <cstdio>
#include <random>
#include <vector>
#include "SequenceUtils/PackedSequence.h"
using namespace protal;
int main() {
    constexpr size_t N = size_t{1} << 30;
    std::mt19937_64 rng(1);
    std::vector<char> bytes(N);
    const char acgt[4] = {'A', 'C', 'G', 'T'};
    for (size_t i = 0; i < N; i += 8) { uint64_t r = rng(); for (int j = 0; j < 8; j++) bytes[i + j] = acgt[(r >> (2 * j)) & 3]; }
    std::vector<uint8_t> packed(N / 4);
    auto clock = [] { return std::chrono::steady_clock::now(); };
    for (bool avx2 : {false, true}) {
        if (avx2 && !packed::CpuHasAvx2()) continue;
        packed::UseAvx2(avx2);
        auto t0 = clock();
        packed::Pack(bytes.data(), N, packed.data());
        double s = std::chrono::duration<double>(clock() - t0).count();
        printf("pack %s: %.2f s for %.0f Mbase (%.2f GB/s of bases)\n", avx2 ? "avx2 " : "table", s, N / 1e6, N / s / 1e9);
    }
    for (size_t L : {300, 1000, 3000}) {
        const size_t reps = 1'000'000;
        std::vector<size_t> pos(reps); for (auto& p : pos) p = (rng() % (N - 2 * L - 8)) & ~size_t{3};
        std::vector<char> out(L + 256);
        auto time = [&](auto f) { auto t0 = clock(); uint64_t sink = 0; for (size_t p : pos) sink += f(p);
            double ns = std::chrono::duration<double, std::nano>(clock() - t0).count() / reps; if (sink == 42) puts(""); return ns; };
        double t_copy = time([&](size_t p) { memcpy(out.data(), &bytes[p], L); return out[L / 2]; });
        packed::UseAvx2(false);
        double t_tab = time([&](size_t p) { packed::Unpack(&packed[p / 4], L, out.data()); return out[L / 2]; });
        double t_avx = -1;
        if (packed::CpuHasAvx2()) { packed::UseAvx2(true); t_avx = time([&](size_t p) { packed::Unpack(&packed[p / 4], L, out.data()); return out[L / 2]; }); }
        double t_seq = -1;
        if (packed::CpuHasAvx2()) t_seq = time([&](size_t p) { GeneSequence g(&packed[p / 4], L); return g[L / 2]; });
        printf("gene of %4zu bases: byte copy %6.1f ns | unpack table %6.1f ns | unpack avx2 %6.1f ns | GeneSequence avx2 %6.1f ns\n", L, t_copy, t_tab, t_avx, t_seq);
    }
}
