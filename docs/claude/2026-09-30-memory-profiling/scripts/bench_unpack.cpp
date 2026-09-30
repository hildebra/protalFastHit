// Cost of reading a gene window from a 2-bit store instead of a byte-per-base arena:
//   g++ -O2 -march=x86-64 bench_unpack.cpp -o bench_unpack && ./bench_unpack
// 1 Gbase store (250 MB packed, 1 GB as bytes); random windows of W bases, as an alignment or seed
// extension reads them. Reports ns per window: touching the bytes arena (sum of the window, the
// baseline: what the code does today) and unpacking with a 256-entry table (scalar, baseline x86-64).
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <random>
#include <vector>
int main() {
    constexpr size_t N = size_t{1} << 30;
    std::mt19937_64 rng(1);
    std::vector<uint8_t> bytes(N);
    std::vector<uint8_t> packed(N / 4);
    const char acgt[4] = {'A', 'C', 'G', 'T'};
    for (size_t i = 0; i < N; i += 8) { uint64_t r = rng(); for (int j = 0; j < 8; j++) bytes[i + j] = acgt[(r >> (2 * j)) & 3]; }
    for (size_t i = 0; i < N; i += 4) packed[i / 4] = (bytes[i] >> 1 & 3) | ((bytes[i + 1] >> 1 & 3) << 2) | ((bytes[i + 2] >> 1 & 3) << 4) | ((bytes[i + 3] >> 1 & 3) << 6);
    // 'A'=0x41 'C'=0x43 'G'=0x47 'T'=0x54: (c>>1)&3 = 0,1,3,2 -> a code, decoded below with its own table
    uint8_t code2char[4]; for (int c = 0; c < 4; c++) for (int k = 0; k < 4; k++) if ((acgt[k] >> 1 & 3) == c) code2char[c] = acgt[k];
    uint32_t lut[256];
    for (int b = 0; b < 256; b++) { uint8_t o[4]; for (int j = 0; j < 4; j++) o[j] = code2char[(b >> (2 * j)) & 3]; memcpy(&lut[b], o, 4); }
    for (size_t W : {64, 150, 350, 1000}) {
        const size_t reps = 2'000'000;
        std::vector<size_t> pos(reps); for (auto& p : pos) p = rng() % (N - 2 * W - 8) & ~size_t{3};
        std::vector<uint8_t> out(W + 16);
        auto time = [&](auto f) { auto t0 = std::chrono::steady_clock::now(); uint64_t sink = 0; for (size_t p : pos) sink += f(p);
            double ns = std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - t0).count() / reps; if (sink == 42) puts(""); return ns; };
        double t_bytes = time([&](size_t p) { memcpy(out.data(), &bytes[p], W); return out[W / 2]; });
        double t_lut = time([&](size_t p) { const uint8_t* s = &packed[p / 4]; uint32_t* o = (uint32_t*)out.data(); for (size_t i = 0; i < (W + 3) / 4; i++) o[i] = lut[s[i]]; return out[W / 2]; });
        printf("window %5zu bases: bytes memcpy %7.1f ns   2-bit table %7.1f ns   (+%.1f ns)%s\n", W, t_bytes, t_lut, t_lut - t_bytes, t_pdep >= 0 ? "" : "");
        (void)t_pdep;
    }
}
