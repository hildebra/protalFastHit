// Bgzf.h - gzip as protal writes it, and reads it where it can, with libdeflate.
//
// protal writes gzip (.sam.gz, the simulator's .fq.gz) as BGZF (SAM specification, section 4.1):
// gzip members of at most 64 KB, each with its size in a "BC" extra field, then an empty member
// as end-of-file marker. zcat, gzip and zlib read it as any multi-member gzip file; htslib reads
// its blocks independently. Every block is compressed and decompressed whole, which is what
// libdeflate does (inflating about 3x as fast as zlib and 2x as fast as zlib-ng, but it cannot
// stream). Other gzip files, e.g. the single-member .fq.gz of sequencers, are read with zlib-ng's
// streaming inflate (ThreadedGzStream.h).
#pragma once

#include <libdeflate.h>

#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <thread>
#include <vector>

namespace protal::bgzf {
    inline constexpr size_t kBlockInput = 0xff00;      // input bytes per block, as htslib
    inline constexpr size_t kMaxBlock = 0x10000;       // a block's size must fit its 16-bit field
    inline constexpr size_t kHeaderBytes = 18, kFooterBytes = 8;
    inline constexpr int kLevel = 6;                   // gzip's default, as pigz
    inline constexpr unsigned char kHeader[16] = { 0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, 'B', 'C', 2, 0 };
    inline constexpr unsigned char kEof[28] = { 0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, 'B', 'C', 2, 0,
                                                0x1b, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    inline void PutLE16(unsigned char* p, uint32_t v) {
        p[0] = static_cast<unsigned char>(v & 0xff);
        p[1] = static_cast<unsigned char>(v >> 8 & 0xff);
    }

    inline void PutLE32(unsigned char* p, uint32_t v) {
        for (int i = 0; i < 4; i++) p[i] = static_cast<unsigned char>(v >> (8 * i) & 0xff);
    }

    inline uint32_t LE32(unsigned char const* p) {
        return uint32_t(p[0]) | uint32_t(p[1]) << 8 | uint32_t(p[2]) << 16 | uint32_t(p[3]) << 24;
    }

    // A libdeflate compressor at one level, for one thread.
    class Deflater {
    public:
        explicit Deflater(int level) : m_c(libdeflate_alloc_compressor(level)) {}
        ~Deflater() {
            if (m_c) libdeflate_free_compressor(m_c);
        }
        Deflater(Deflater const&) = delete;
        Deflater& operator=(Deflater const&) = delete;

        // Raw deflate of [in, in + size) into out; the compressed size, or 0 if it does not fit.
        size_t Compress(char const* in, size_t size, unsigned char* out, size_t capacity) {
            return m_c ? libdeflate_deflate_compress(m_c, in, size, out, capacity) : 0;
        }

    private:
        libdeflate_compressor* m_c;
    };

    // Appends [data, data + size) to out as BGZF blocks, deflated at kLevel. A block that does not
    // shrink below the 64 KB limit (random data) is stored instead. False only if libdeflate fails
    // (out of memory).
    inline bool Compress(char const* data, size_t size, std::string& out) {
        thread_local Deflater deflater(kLevel);
        thread_local Deflater store(0);
        constexpr size_t kPayload = kMaxBlock - kHeaderBytes - kFooterBytes;
        for (size_t offset = 0; offset < size; offset += kBlockInput) {
            size_t const n = std::min(kBlockInput, size - offset);
            size_t const start = out.size();
            out.resize(start + kMaxBlock);
            auto* block = reinterpret_cast<unsigned char*>(out.data() + start);
            size_t compressed = deflater.Compress(data + offset, n, block + kHeaderBytes, kPayload);
            if (compressed == 0) compressed = store.Compress(data + offset, n, block + kHeaderBytes, kPayload);
            if (compressed == 0) return false;
            size_t const total = kHeaderBytes + compressed + kFooterBytes;
            std::memcpy(block, kHeader, sizeof(kHeader));
            PutLE16(block + 16, static_cast<uint32_t>(total - 1));
            PutLE32(block + kHeaderBytes + compressed, libdeflate_crc32(0, data + offset, n));
            PutLE32(block + kHeaderBytes + compressed + 4, static_cast<uint32_t>(n));
            out.resize(start + total);
        }
        return true;
    }

    // Whether the 18 bytes at h begin a BGZF block (gzip, FEXTRA with the BC subfield first).
    inline bool IsBlockHeader(unsigned char const* h) {
        return h[0] == 0x1f && h[1] == 0x8b && h[2] == 8 && (h[3] & 4) && h[10] == 6 && h[11] == 0 &&
               h[12] == 'B' && h[13] == 'C' && h[14] == 2 && h[15] == 0;
    }

    // The size of the block whose header is at h, header and footer included.
    inline size_t BlockSize(unsigned char const* h) {
        return (size_t(h[16]) | size_t(h[17]) << 8) + 1;
    }

    // Decompresses the whole block [block, block + size) into out (capacity bytes, at least its
    // content): its content size in n, or false and error for a corrupt block.
    inline bool DecompressBlock(unsigned char const* block, size_t size, char* out, size_t capacity, size_t& n, std::string& error) {
        struct Inflater {
            libdeflate_decompressor* d = libdeflate_alloc_decompressor();
            ~Inflater() {
                if (d) libdeflate_free_decompressor(d);
            }
        };
        thread_local Inflater inflater;
        if (size < kHeaderBytes + kFooterBytes || !IsBlockHeader(block) || BlockSize(block) != size) {
            error = "not a BGZF block";
            return false;
        }
        uint32_t const crc = LE32(block + size - 8), isize = LE32(block + size - 4);
        if (isize > capacity) {
            error = "a BGZF block holds more than 64 KB";
            return false;
        }
        if (!inflater.d) {
            error = "out of memory";
            return false;
        }
        if (libdeflate_deflate_decompress(inflater.d, block + kHeaderBytes, size - kHeaderBytes - kFooterBytes, out, isize,
                                          nullptr) != LIBDEFLATE_SUCCESS) {
            error = "a corrupt BGZF block";
            return false;
        }
        if (libdeflate_crc32(0, out, isize) != crc) {
            error = "a BGZF block fails its CRC check";
            return false;
        }
        n = isize;
        return true;
    }

    // Whether the file begins with a BGZF block header.
    inline bool StartsAsBgzf(std::string const& path) {
        std::ifstream is(path, std::ios::binary);
        unsigned char b[kHeaderBytes] = {};
        is.read(reinterpret_cast<char*>(b), sizeof(b));
        return is.gcount() == static_cast<std::streamsize>(sizeof(b)) && IsBlockHeader(b);
    }

    // Whether the file ends with the BGZF end-of-file block.
    inline bool EndsWithEof(std::string const& path) {
        std::error_code ec;
        auto const size = std::filesystem::file_size(path, ec);
        if (ec || size < sizeof(kEof)) return false;
        std::ifstream is(path, std::ios::binary);
        is.seekg(static_cast<std::streamoff>(size - sizeof(kEof)));
        unsigned char b[sizeof(kEof)] = {};
        is.read(reinterpret_cast<char*>(b), sizeof(b));
        return is.gcount() == static_cast<std::streamsize>(sizeof(b)) && std::memcmp(b, kEof, sizeof(b)) == 0;
    }

    // Writes src as BGZF to dst with `threads` threads. Chunks of whole blocks are compressed in
    // parallel and written in order, so the output is the same for any thread count. Returns an
    // error message, empty on success.
    inline std::string CompressFile(std::string const& src, std::string const& dst, int threads) {
        constexpr size_t kChunk = 64 * kBlockInput;  // ~4 MB, whole blocks
        std::FILE* in = std::fopen(src.c_str(), "rb");
        if (!in) return "cannot read " + src + ": " + std::strerror(errno);
        std::FILE* out = std::fopen(dst.c_str(), "wb");
        if (!out) {
            std::string const error = "cannot write " + dst + ": " + std::strerror(errno);
            std::fclose(in);
            return error;
        }
        size_t const workers = static_cast<size_t>(std::max(1, threads));
        std::vector<std::string> chunks(workers), packed(workers);
        std::vector<char> failed(workers);
        std::string error;
        bool done = false;
        while (!done && error.empty()) {
            size_t filled = 0;
            for (; filled < workers; filled++) {
                chunks[filled].resize(kChunk);
                size_t const n = std::fread(chunks[filled].data(), 1, kChunk, in);
                chunks[filled].resize(n);
                if (n < kChunk) {
                    if (std::ferror(in)) error = "reading " + src + " failed";
                    done = true;
                    if (n > 0) filled++;
                    break;
                }
            }
            auto compress = [&](size_t i) {
                packed[i].clear();
                failed[i] = !Compress(chunks[i].data(), chunks[i].size(), packed[i]);
            };
            std::vector<std::thread> pool;
            for (size_t i = 1; i < filled; i++) pool.emplace_back(compress, i);
            if (filled > 0) compress(0);
            for (auto& t : pool) t.join();
            for (size_t i = 0; i < filled && error.empty(); i++) {
                if (failed[i]) error = "compressing " + src + " failed";
                else if (std::fwrite(packed[i].data(), 1, packed[i].size(), out) != packed[i].size()) error = "writing " + dst + " failed: " + std::strerror(errno);
            }
        }
        if (error.empty() && std::fwrite(kEof, 1, sizeof(kEof), out) != sizeof(kEof)) error = "writing " + dst + " failed";
        std::fclose(in);
        if (std::fclose(out) != 0 && error.empty()) error = "writing " + dst + " failed: " + std::strerror(errno);
        return error;
    }
}
