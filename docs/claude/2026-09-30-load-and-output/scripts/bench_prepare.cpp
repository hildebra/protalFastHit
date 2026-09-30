// LoadAllGenomes' serial step: every gene's std::string sized and zero-filled (Gene::PrepareSequence),
// against the same done by 8 threads, and against one arena allocation. n genes of 900-1199 bytes.
//   bench_prepare n_genes
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <thread>
#include <vector>

int main(int argc, char** argv) {
    size_t const n = argc > 1 ? std::strtoull(argv[1], nullptr, 10) : 4000000;
    using clock = std::chrono::steady_clock;
    auto secs = [](clock::time_point a) { return std::chrono::duration<double>(clock::now() - a).count(); };
    std::vector<size_t> len(n);
    size_t total = 0;
    for (size_t i = 0; i < n; i++) total += len[i] = 900 + (i * 7919) % 300;
    {
        std::vector<std::string> genes(n);
        auto t = clock::now();
        for (size_t i = 0; i < n; i++) genes[i].assign(len[i], '\0');
        std::printf("%zu genes, %.2f GB: serial assign %.2f s\n", n, total / 1e9, secs(t));
    }
    {
        std::vector<std::string> genes(n);
        auto t = clock::now();
        std::vector<std::thread> pool;
        for (int w = 0; w < 8; w++) pool.emplace_back([&, w] { for (size_t i = w; i < n; i += 8) genes[i].assign(len[i], '\0'); });
        for (auto& th : pool) th.join();
        std::printf("  8 threads assign %.2f s\n", secs(t));
    }
    {
        auto t = clock::now();
        char* arena = static_cast<char*>(std::malloc(total));
        std::vector<std::thread> pool;
        size_t const part = (total + 7) / 8;
        for (int w = 0; w < 8; w++) pool.emplace_back([&, w] { size_t b = w * part, e = std::min(total, b + part); if (b < e) std::memset(arena + b, 0, e - b); });
        for (auto& th : pool) th.join();
        std::printf("  one arena, touched by 8 threads %.2f s\n", secs(t));
        std::free(arena);
    }
    return 0;
}
