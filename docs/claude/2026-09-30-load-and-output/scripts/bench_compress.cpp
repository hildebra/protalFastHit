// Compression libraries on protal's inputs and outputs, one thread, data in memory.
//   bench_compress inflate FILE.gz   zlib streaming inflate (as ThreadedGzStreambuf, via inflate())
//                                    against libdeflate's whole-member gunzip
//   bench_compress deflate FILE      FILE (a SAM) cut into blocks, each compressed as its own gzip
//                                    member (zlib, libdeflate) or zstd frame; ratio, compression and
//                                    decompression speed
// g++ -O2 -march=x86-64-v3 bench_compress.cpp -o bench_compress -lz -ldeflate -lzstd
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iterator>
#include <string>
#include <vector>
#include <libdeflate.h>
#include <zlib.h>
#include <zstd.h>

static std::string Slurp(char const* path) {
    std::ifstream is(path, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(is), {});
}

static double Seconds(std::function<void()> const& f) {
    auto const t0 = std::chrono::steady_clock::now();
    f();
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

static double Best(int reps, std::function<void()> const& f) {
    double best = 1e30;
    for (int i = 0; i < reps; i++) best = std::min(best, Seconds(f));
    return best;
}

static size_t ZlibInflate(std::string const& gz, std::vector<char>& chunk) {
    z_stream s{};
    inflateInit2(&s, 16 + MAX_WBITS);
    s.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(gz.data()));
    s.avail_in = static_cast<uInt>(gz.size());
    size_t total = 0;
    for (;;) {
        s.next_out = reinterpret_cast<Bytef*>(chunk.data());
        s.avail_out = static_cast<uInt>(chunk.size());
        int const r = inflate(&s, Z_NO_FLUSH);
        total += chunk.size() - s.avail_out;
        if (r == Z_STREAM_END) {
            if (s.avail_in == 0) break;
            inflateReset(&s);  // next member
            continue;
        }
        if (r != Z_OK) { std::fprintf(stderr, "inflate error %d\n", r); std::exit(1); }
    }
    inflateEnd(&s);
    return total;
}

static int Inflate(char const* path) {
    std::string const gz = Slurp(path);
    uint32_t isize;
    std::memcpy(&isize, gz.data() + gz.size() - 4, 4);
    std::vector<char> chunk(1 << 20), out(isize);
    size_t n = 0;
    double const tz = Best(3, [&] { n = ZlibInflate(gz, chunk); });
    libdeflate_decompressor* d = libdeflate_alloc_decompressor();
    size_t got = 0;
    double const tl = Best(3, [&] {
        if (libdeflate_gzip_decompress(d, gz.data(), gz.size(), out.data(), out.size(), &got) != LIBDEFLATE_SUCCESS) {
            std::fprintf(stderr, "libdeflate failed (multi-member or > 4 GB?)\n");
            std::exit(1);
        }
    });
    libdeflate_free_decompressor(d);
    double const tc = Best(3, [&] { volatile uint32_t c = crc32(0, reinterpret_cast<Bytef const*>(out.data()), out.size()); (void)c; });
    double const mb = n / 1e6;
    std::printf("%s: %.0f MB inflated (%zu MB gz)\n", path, mb, gz.size() / 1000000);
    std::printf("  zlib inflate()          %.3f s  %6.0f MB/s\n", tz, mb / tz);
    std::printf("  libdeflate gunzip       %.3f s  %6.0f MB/s  (%.2fx)\n", tl, mb / tl, tz / tl);
    std::printf("  zlib crc32 alone        %.3f s  %6.0f MB/s\n", tc, mb / tc);
    return got == n ? 0 : 1;
}

struct Codec {
    std::string name;
    std::function<size_t(char const*, size_t, std::vector<char>&)> compress;       // returns bytes written
    std::function<size_t(char const*, size_t, char*, size_t)> decompress;          // returns bytes out
};

static int Deflate(char const* path) {
    std::string const data = Slurp(path);
    std::vector<Codec> codecs;
    for (int level : {1, 6}) {
        codecs.push_back({"zlib " + std::to_string(level),
            [level](char const* in, size_t n, std::vector<char>& out) {
                z_stream s{};
                deflateInit2(&s, level, Z_DEFLATED, 16 + MAX_WBITS, 8, Z_DEFAULT_STRATEGY);
                out.resize(deflateBound(&s, n) + 64);
                s.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(in));
                s.avail_in = static_cast<uInt>(n);
                s.next_out = reinterpret_cast<Bytef*>(out.data());
                s.avail_out = static_cast<uInt>(out.size());
                deflate(&s, Z_FINISH);
                size_t const w = s.total_out;
                deflateEnd(&s);
                return w;
            },
            [](char const* in, size_t n, char* out, size_t cap) {
                z_stream s{};
                inflateInit2(&s, 16 + MAX_WBITS);
                s.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(in));
                s.avail_in = static_cast<uInt>(n);
                s.next_out = reinterpret_cast<Bytef*>(out);
                s.avail_out = static_cast<uInt>(cap);
                inflate(&s, Z_FINISH);
                size_t const w = s.total_out;
                inflateEnd(&s);
                return w;
            }});
    }
    for (int level : {1, 4, 6}) {
        codecs.push_back({"libdeflate " + std::to_string(level),
            [level](char const* in, size_t n, std::vector<char>& out) {
                thread_local libdeflate_compressor* c = nullptr;
                thread_local int c_level = -1;
                if (c_level != level) { if (c) libdeflate_free_compressor(c); c = libdeflate_alloc_compressor(level); c_level = level; }
                out.resize(libdeflate_gzip_compress_bound(c, n));
                return libdeflate_gzip_compress(c, in, n, out.data(), out.size());
            },
            [](char const* in, size_t n, char* out, size_t cap) {
                thread_local libdeflate_decompressor* d = libdeflate_alloc_decompressor();
                size_t got = 0;
                libdeflate_gzip_decompress(d, in, n, out, cap, &got);
                return got;
            }});
    }
    for (int level : {1, 3}) {
        codecs.push_back({"zstd " + std::to_string(level),
            [level](char const* in, size_t n, std::vector<char>& out) {
                out.resize(ZSTD_compressBound(n));
                return ZSTD_compress(out.data(), out.size(), in, n, level);
            },
            [](char const* in, size_t n, char* out, size_t cap) { return ZSTD_decompress(out, cap, in, n); }});
    }
    std::printf("%s: %.0f MB\n", path, data.size() / 1e6);
    for (size_t block : {size_t{64} << 10, size_t{1} << 20}) {
        for (auto const& codec : codecs) {
            std::vector<std::vector<char>> frames((data.size() + block - 1) / block);
            size_t total = 0;
            double const tc = Best(2, [&] {
                total = 0;
                for (size_t i = 0; i < frames.size(); i++) {
                    size_t const off = i * block, n = std::min(block, data.size() - off);
                    size_t const w = codec.compress(data.data() + off, n, frames[i]);
                    frames[i].resize(w);
                    total += w;
                }
            });
            std::vector<char> out(block);
            size_t back = 0;
            double const td = Best(2, [&] {
                back = 0;
                for (auto const& f : frames) back += codec.decompress(f.data(), f.size(), out.data(), out.size());
            });
            std::printf("  %4zu KB blocks  %-13s ratio %5.2f  compress %6.0f MB/s  decompress %6.0f MB/s%s\n", block >> 10,
                        codec.name.c_str(), double(data.size()) / total, data.size() / 1e6 / tc, data.size() / 1e6 / td,
                        back == data.size() ? "" : "  ROUND TRIP FAILED");
        }
    }
    return 0;
}

int main(int argc, char** argv) {
    if (argc == 3 && std::string(argv[1]) == "inflate") return Inflate(argv[2]);
    if (argc == 3 && std::string(argv[1]) == "deflate") return Deflate(argv[2]);
    std::fprintf(stderr, "usage: bench_compress inflate FILE.gz | deflate FILE\n");
    return 2;
}
