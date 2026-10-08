// Probe of the fix in Bgzf.h: Crc32cWord against the crc32 instruction, the cost of placing a Deflater's stream, and
// whether small FASTQ blocks get the same bytes from Deflaters at different addresses: placed (the fix) and on the heap
// (as before the fix). Also block 6 of the failing run, in this process (another address than the run's).
#include "IO/Bgzf.h"

#include <nmmintrin.h>

#include <chrono>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <memory>
#include <random>
#include <string>
#include <thread>
#include <vector>

// The Deflater before the fix: its stream on the heap.
class HeapDeflater {
public:
    HeapDeflater() : m_stream(std::make_unique<isal_zstream>()), m_level_buffer(protal::bgzf::kLevelBuffer) {}
    size_t Compress(char const* in, size_t size, unsigned char* out, size_t capacity) {
        isal_zstream& stream = *m_stream;
        isal_deflate_stateless_init(&stream);
        stream.level = protal::bgzf::kLevel;
        stream.level_buf = m_level_buffer.data();
        stream.level_buf_size = static_cast<uint32_t>(m_level_buffer.size());
        stream.end_of_stream = 1;
        stream.flush = NO_FLUSH;
        stream.next_in = reinterpret_cast<uint8_t*>(const_cast<char*>(in));
        stream.avail_in = static_cast<uint32_t>(size);
        stream.next_out = out;
        stream.avail_out = static_cast<uint32_t>(capacity);
        return isal_deflate_stateless(&stream) == COMP_OK ? stream.total_out : 0;
    }
    void const* Address() const { return m_stream.get(); }

private:
    std::unique_ptr<isal_zstream> m_stream;
    std::vector<uint8_t> m_level_buffer;
};

// Blocks of two FASTQ records, the second named as the first, the names random (the unit test's corpus).
static std::vector<std::string> NamedPairs(size_t count) {
    std::mt19937 rng(11);
    std::vector<std::string> blocks;
    for (size_t b = 0; b < count; ++b) {
        std::string name = "@";
        for (int i = 0; i < 3; i++) name += static_cast<char>('A' + rng() % 26);
        std::string block;
        for (int r = 1; r <= 2; ++r) {
            std::string seq(50, 'A');
            for (auto& c : seq) c = "ACGT"[rng() % 4];
            block += name + "_" + std::to_string(r) + "\n" + seq + "\n+\n" + std::string(50, 'I') + "\n";
        }
        blocks.push_back(std::move(block));
    }
    return blocks;
}

template <typename D>
static std::vector<std::string> CompressAll(D& d, std::vector<std::string> const& blocks) {
    std::vector<std::string> out;
    std::vector<unsigned char> buffer(protal::bgzf::kMaxBlock);
    for (auto const& block : blocks) {
        size_t const n = d.Compress(block.data(), block.size(), buffer.data(), buffer.size());
        out.emplace_back(reinterpret_cast<char const*>(buffer.data()), n);
    }
    return out;
}

int main(int argc, char** argv) {
    namespace bgzf = protal::bgzf;
    std::mt19937_64 rng(1);
    int bad = 0;
    for (int i = 0; i < 1000000; ++i) {
        auto const w = static_cast<uint32_t>(rng());
        bad += bgzf::detail::Crc32cWord(w) != _mm_crc32_u32(0, w);
    }
    std::printf("Crc32cWord against _mm_crc32_u32: %d of 1000000 differ\n", bad);

    for (int count : {20, 100}) {
        auto const t0 = std::chrono::steady_clock::now();
        std::vector<std::unique_ptr<bgzf::Deflater>> placed;
        for (int i = 0; i < count; ++i) placed.push_back(std::make_unique<bgzf::Deflater>());
        double const ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        int ok = 0;
        for (auto const& d : placed) ok += d->Placed();
        std::printf("%d placed Deflaters at once in %.1f ms (%.3f ms each), %d placed\n", count, ms, ms / count, ok);
    }
    {  // created and freed one after the other: the place is used again
        auto const t0 = std::chrono::steady_clock::now();
        for (int i = 0; i < 100; ++i) bgzf::Deflater d;
        double const ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        std::printf("100 Deflaters one after the other in %.1f ms\n", ms);
    }

    auto const blocks = NamedPairs(3000);
    std::vector<std::unique_ptr<HeapDeflater>> heap;
    std::vector<std::unique_ptr<bgzf::Deflater>> placed;
    for (int i = 0; i < 8; ++i) {
        heap.push_back(std::make_unique<HeapDeflater>());
        placed.push_back(std::make_unique<bgzf::Deflater>());
    }
    auto const heap0 = CompressAll(*heap[0], blocks), placed0 = CompressAll(*placed[0], blocks);
    for (size_t i = 1; i < heap.size(); ++i) {
        auto const h = CompressAll(*heap[i], blocks), p = CompressAll(*placed[i], blocks);
        int dh = 0, dp = 0;
        for (size_t b = 0; b < blocks.size(); ++b) {
            dh += h[b] != heap0[b];
            dp += p[b] != placed0[b];
        }
        std::printf("deflater %zu: heap (G %4u vs %4u under 255) %3d of %zu blocks differ from deflater 0; placed %d\n", i,
                    bgzf::detail::Crc32cWord(static_cast<uint32_t>(reinterpret_cast<uintptr_t>(heap[i]->Address()) >> 16)) & 255,
                    bgzf::detail::Crc32cWord(static_cast<uint32_t>(reinterpret_cast<uintptr_t>(heap[0]->Address()) >> 16)) & 255,
                    dh, blocks.size(), dp);
    }
    // bgzf::Compress on other threads (their own Deflaters) against this one
    auto compress_blocks = [&] {
        std::string out;
        for (auto const& block : blocks) bgzf::Compress(block.data(), block.size(), out);
        return out;
    };
    std::string const here = compress_blocks();
    std::vector<std::string> there(4);
    std::vector<std::thread> threads;
    for (auto& t : there) threads.emplace_back([&t, &compress_blocks] { t = compress_blocks(); });
    for (auto& t : threads) t.join();
    int same = 0;
    for (auto const& t : there) same += t == here;
    std::printf("bgzf::Compress on 4 other threads: %d of 4 the same as on the main thread\n", same);

    if (argc > 3) {  // block 6 of the failing run: the placed stream's bytes
        std::ifstream f(argv[1], std::ios::binary);
        std::string const all((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
        std::string const block = all.substr(std::strtoull(argv[2], nullptr, 10), std::strtoull(argv[3], nullptr, 10));
        std::vector<unsigned char> out(0x10000);
        size_t const n = placed[0]->Compress(block.data(), block.size(), out.data(), out.size());
        std::printf("block 6 with the placed stream: %zu bytes of deflate\n", n);
    }
    return 0;
}
