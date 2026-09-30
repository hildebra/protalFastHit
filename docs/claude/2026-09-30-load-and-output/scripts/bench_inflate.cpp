// Inflating one gzip file (from the page cache) the way ThreadedGzStreambuf does: gzread in reads of
// 1 MB with a 128 KB input buffer, with zlib or zlib-ng; and libdeflate on the whole file in memory
// (one member only), as the upper bound. Built once per library:
//   g++ -O2 -DWITH_ZLIB bench_inflate.cpp -lz
//   g++ -O2 -DWITH_ZLIB_NG -DWITH_GZFILEOP -I<zlib-ng build dir> bench_inflate.cpp <build dir>/libz-ng.a
//   g++ -O2 -DWITH_LIBDEFLATE bench_inflate.cpp -ldeflate
//   bench_inflate FILE.gz [REPEATS]
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>
#if defined(WITH_ZLIB)
#include <zlib.h>
#define LIBRARY "zlib " ZLIB_VERSION
#elif defined(WITH_ZLIB_NG)
#include <zlib-ng.h>
#define LIBRARY "zlib-ng " ZLIBNG_VERSION
#elif defined(WITH_LIBDEFLATE)
#include <libdeflate.h>
#define LIBRARY "libdeflate " LIBDEFLATE_VERSION_STRING
#endif

#if defined(WITH_LIBDEFLATE)
// The file and a buffer for its inflated bytes, read and touched once, outside the timing (a fresh
// buffer's page faults would cost more than the inflating).
struct Input {
    std::string in;
    std::vector<char> out;
    explicit Input(char const* path) {
        std::ifstream is(path, std::ios::binary);
        in.assign(std::istreambuf_iterator<char>(is), {});
        // ISIZE, the last 4 bytes: the inflated size modulo 2^32 (the files here are smaller).
        auto const* end = reinterpret_cast<unsigned char const*>(in.data() + in.size());
        out.assign(end[-4] | end[-3] << 8 | end[-2] << 16 | static_cast<size_t>(end[-1]) << 24, 1);
    }
    size_t Inflate() {
        size_t total = 0;
        libdeflate_decompressor* d = libdeflate_alloc_decompressor();
        if (libdeflate_gzip_decompress(d, in.data(), in.size(), out.data(), out.size(), &total) != LIBDEFLATE_SUCCESS) total = 0;
        libdeflate_free_decompressor(d);
        return total;
    }
};
#else
// The file read with gzread into one reused 1 MB block, as ThreadedGzStreambuf reads it.
struct Input {
    char const* path;
    std::vector<char> block = std::vector<char>(size_t{1} << 20);
    explicit Input(char const* p) : path(p) {}
    size_t Inflate() {
        size_t total = 0;
#if defined(WITH_ZLIB)
        gzFile f = gzopen(path, "rb");
        gzbuffer(f, 1u << 17);
        int n;
        while ((n = gzread(f, block.data(), static_cast<unsigned>(block.size()))) > 0) total += static_cast<size_t>(n);
        gzclose(f);
#else
        gzFile f = zng_gzopen(path, "rb");
        zng_gzbuffer(f, 1u << 17);
        int32_t n;
        while ((n = zng_gzread(f, block.data(), static_cast<uint32_t>(block.size()))) > 0) total += static_cast<size_t>(n);
        zng_gzclose(f);
#endif
        return total;
    }
};
#endif

int main(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: bench_inflate FILE.gz [REPEATS]\n"); return 2; }
    int const repeats = argc > 2 ? std::atoi(argv[2]) : 3;
    Input input(argv[1]);
    for (int r = 0; r < repeats; r++) {
        auto const start = std::chrono::steady_clock::now();
        size_t const bytes = input.Inflate();
        double const s = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        std::printf("%s\t%zu bytes\t%.3f s\t%.0f MB/s\n", LIBRARY, bytes, s, bytes / s / 1e6);
    }
    return 0;
}
