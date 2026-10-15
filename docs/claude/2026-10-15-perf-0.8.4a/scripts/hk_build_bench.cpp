// hk_build_bench - experiment only (docs/claude/2026-10-15-perf-0.8.4a, section 11): the time HashedKeys::Build takes
// for a compact key map of r226's size (2^30 cores, 9.2% used, slot counts geometric with mean 44), on 1 to N threads.
// Compiled against either version of HashedKeys.h (the include path picks it).
//   hk_build_bench THREADS...
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>
#include "Hash/CompactKeys.h"
#include "Hash/HashedKeys.h"

int main(int argc, char** argv) {
    unsigned const bits = 30;
    uint64_t const cores = uint64_t{1} << bits;
    std::mt19937_64 rng(42);
    std::bernoulli_distribution used(0.092);
    std::geometric_distribution<uint64_t> extra(1.0 / 44.0);
    std::vector<uint16_t> slots(cores, 0);
    for (uint64_t c = 0; c < cores; c++) {
        if (used(rng)) slots[c] = static_cast<uint16_t>(std::min<uint64_t>(1 + extra(rng), 3071));
    }
    protal::CompactKeys compact;
    compact.Init(cores);
    for (uint64_t c = 0; c < cores; c++) {
        if (slots[c]) compact.MarkUsed(c);
    }
    uint64_t const n_used = compact.CountRanks();
    compact.AllocateStarts();
    uint64_t next = 0, rank = 0;
    for (uint64_t c = 0; c < cores; c++) {
        if (!slots[c]) continue;
        if (rank % protal::CompactKeys::kGroup == 0) compact.SetBase(rank, next);
        next += slots[c];
        rank++;
    }
    if (n_used % protal::CompactKeys::kGroup == 0) compact.SetBase(n_used, next);
    next = 0;
    rank = 0;
    for (uint64_t c = 0; c < cores; c++) {
        if (!slots[c]) continue;
        if (!compact.SetOffset(rank++, next)) {
            std::fprintf(stderr, "offset does not fit\n");
            return 2;
        }
        next += slots[c];
    }
    if (!compact.SetOffset(n_used, next)) return 2;
    std::vector<uint16_t>().swap(slots);
    std::printf("compact: %llu used cores, %llu slots\n", static_cast<unsigned long long>(n_used), static_cast<unsigned long long>(next));
    for (int a = 1; a < argc; a++) {
        int const threads = std::atoi(argv[a]);
        for (int rep = 0; rep < 2; rep++) {
            protal::HashedKeys hashed;
            auto const t0 = std::chrono::steady_clock::now();
            bool const ok = hashed.Build(compact, bits, threads);
            double const s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            std::printf("threads %d\trep %d\tok %d\tbuild %.3f s\tpast home %llu\tmax %llu\t%.3f GB\n", threads, rep, ok ? 1 : 0, s,
                        static_cast<unsigned long long>(hashed.PastHome()), static_cast<unsigned long long>(hashed.MaxDisplacement()),
                        static_cast<double>(hashed.Bytes()) / 1e9);
        }
    }
    return 0;
}
