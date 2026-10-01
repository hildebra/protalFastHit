// What a random read of a 3 GB table costs on this machine, which bounds the seeding stage (22 key map
// lookups per read). Dependent loads give the latency, independent ones the memory-level parallelism,
// prefetching a batch ahead what software pipelining could win. The table is 24-byte blocks like the
// key map (134M of them), backed by transparent huge pages where the system allows.
//   g++ -O3 -march=x86-64-v3 bench_mem.cpp -o bench_mem && ./bench_mem [huge|small]
#include <sys/mman.h>

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>

using Clock = std::chrono::steady_clock;
static double Ns(Clock::time_point a, Clock::time_point b) { return std::chrono::duration<double, std::nano>(b - a).count(); }

int main(int argc, char** argv) {
    bool const huge = !(argc > 1 && std::strcmp(argv[1], "small") == 0);
    size_t const blocks = size_t{134} << 20, bytes = blocks * 24;
    void* mem = mmap(nullptr, bytes, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (mem == MAP_FAILED) return 1;
    madvise(mem, bytes, huge ? MADV_HUGEPAGE : MADV_NOHUGEPAGE);
    auto* t = static_cast<uint8_t*>(mem);
    std::memset(t, 1, bytes);  // touch every page
    std::mt19937_64 rng(1);
    size_t const n = 1 << 22;
    std::vector<uint64_t> idx(n);
    for (auto& i : idx) i = (rng() % blocks) * 24;
    uint64_t sink = 0;
    auto run = [&](char const* what, auto&& body) {
        auto t0 = Clock::now();
        body();
        std::printf("%-46s %7.1f ns per access\n", what, Ns(t0, Clock::now()) / static_cast<double>(n));
    };
    std::printf("table %.1f GB, %s pages\n", bytes / 1e9, huge ? "huge (madvise)" : "4 KB");
    run("independent loads (a lookup per key, sum)", [&] { for (size_t i = 0; i < n; i++) sink += t[idx[i]]; });
    run("dependent loads (pointer chase through the table)", [&] {
        uint64_t p = idx[0];
        for (size_t i = 0; i < n; i++) p = idx[(p + t[p]) % n];
        sink += p;
    });
    for (size_t batch : {8, 22, 44}) {
        // a read: `batch` lookups, then a dependent decision (branch on loaded data), then the next read
        char what[80];
        std::snprintf(what, sizeof what, "reads of %zu lookups, branch on each, no prefetch", batch);
        run(what, [&] {
            for (size_t r = 0; r + batch <= n; r += batch)
                for (size_t j = 0; j < batch; j++) { uint8_t v = t[idx[r + j]]; if (v == 7) sink += 3; else sink += v; }
        });
        std::snprintf(what, sizeof what, "reads of %zu lookups, all prefetched first", batch);
        run(what, [&] {
            for (size_t r = 0; r + batch <= n; r += batch) {
                for (size_t j = 0; j < batch; j++) __builtin_prefetch(t + idx[r + j]);
                for (size_t j = 0; j < batch; j++) { uint8_t v = t[idx[r + j]]; if (v == 7) sink += 3; else sink += v; }
            }
        });
        std::snprintf(what, sizeof what, "reads of %zu lookups, next read's prefetched", batch);
        run(what, [&] {
            for (size_t r = 0; r + 2 * batch <= n; r += batch) {
                for (size_t j = 0; j < batch; j++) __builtin_prefetch(t + idx[r + batch + j]);
                for (size_t j = 0; j < batch; j++) { uint8_t v = t[idx[r + j]]; if (v == 7) sink += 3; else sink += v; }
            }
        });
    }
    std::printf("(%llu)\n", static_cast<unsigned long long>(sink));
    return 0;
}
