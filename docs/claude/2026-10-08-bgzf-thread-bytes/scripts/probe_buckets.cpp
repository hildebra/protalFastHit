// Which value of the third byte's bucket (G = CRC-32C(stream address >> 16) & hash_mask) changes a block's output, and
// why MapStream misses. usage: probe_buckets CONTENT OFFSET LENGTH
#include "IO/Bgzf.h"

#include <nmmintrin.h>
#include <sys/mman.h>

#include <cerrno>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <map>
#include <new>
#include <string>
#include <vector>

static std::uint64_t Fnv(unsigned char const* p, size_t n) {
    std::uint64_t h = 1469598103934665603ULL;
    for (size_t i = 0; i < n; ++i) h = (h ^ p[i]) * 1099511628211ULL;
    return h;
}

// The output of compressing `block` with the stream at `place`.
static std::uint64_t CompressAt(void* place, std::vector<uint8_t>& level, std::string const& block, size_t* size = nullptr) {
    std::vector<unsigned char> out(0x10000);
    auto* stream = new (place) isal_zstream();
    isal_deflate_stateless_init(stream);
    stream->level = 1;
    stream->level_buf = level.data();
    stream->level_buf_size = static_cast<uint32_t>(level.size());
    stream->end_of_stream = 1;
    stream->flush = NO_FLUSH;
    stream->next_in = reinterpret_cast<uint8_t*>(const_cast<char*>(block.data()));
    stream->avail_in = static_cast<uint32_t>(block.size());
    stream->next_out = out.data();
    stream->avail_out = 0x10000 - 26;
    isal_deflate_stateless(stream);
    if (size) *size = stream->total_out;
    return Fnv(out.data(), stream->total_out);
}

static uint32_t Hash4(std::string const& s, size_t at) {
    uint32_t w;
    std::memcpy(&w, s.data() + at, 4);
    return _mm_crc32_u32(0, w);
}

int main(int argc, char** argv) {
    namespace bgzf = protal::bgzf;
    std::ifstream f(argv[1], std::ios::binary);
    std::string const all((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    std::string const block = all.substr(std::strtoull(argv[2], nullptr, 10), std::strtoull(argv[3], nullptr, 10));
    uint32_t const mask = 4095;  // 3580 bytes
    std::printf("block hashes & %u: pos0 %u pos1 %u pos2 %u pos3 %u; at 1859: %u %u %u %u\n", mask, Hash4(block, 0) & mask,
                Hash4(block, 1) & mask, Hash4(block, 2) & mask, Hash4(block, 3) & mask, Hash4(block, 1859) & mask,
                Hash4(block, 1860) & mask, Hash4(block, 1861) & mask, Hash4(block, 1862) & mask);

    // One stream place per bucket value G: find a 64 KB-aligned address with that G in a big reserved region.
    size_t const region_size = size_t{1} << 34;  // 16 GB of address space, PROT_NONE until used
    auto* region = static_cast<unsigned char*>(mmap(nullptr, region_size, PROT_NONE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0));
    if (region == MAP_FAILED) return 3;
    uintptr_t const base = (reinterpret_cast<uintptr_t>(region) + 0xffff) & ~uintptr_t(0xffff);
    std::map<uint32_t, uintptr_t> place_of;  // G -> address
    for (uintptr_t a = base; a + 0x20000 < reinterpret_cast<uintptr_t>(region) + region_size && place_of.size() <= mask; a += 0x10000) {
        place_of.emplace(bgzf::detail::Crc32cWord(static_cast<uint32_t>(a >> 16)) & mask, a);
    }
    std::printf("%zu of %u bucket values have a place\n", place_of.size(), mask + 1);
    std::vector<uint8_t> level(ISAL_DEF_LVL1_DEFAULT);
    std::map<std::uint64_t, std::vector<uint32_t>> by_output;
    std::map<std::uint64_t, size_t> sizes;
    for (auto const& [g, a] : place_of) {
        void* const p = reinterpret_cast<void*>(a);
        mprotect(p, 0x20000, PROT_READ | PROT_WRITE);
        size_t n = 0;
        auto const h = CompressAt(p, level, block, &n);
        by_output[h].push_back(g);
        sizes[h] = n;
        mprotect(p, 0x20000, PROT_NONE);
        madvise(p, 0x20000, MADV_DONTNEED);
    }
    for (auto const& [h, gs] : by_output) {
        std::printf("output %016llx (%zu bytes): %zu bucket values", static_cast<unsigned long long>(h), sizes[h], gs.size());
        if (gs.size() < 10) {
            for (auto g : gs) std::printf(" %u", g);
        }
        std::printf("\n");
    }
    munmap(region, region_size);

    // MapStream's misses: where its scan stops and why.
    for (int k = 0; k < 6; ++k) {
        size_t const size = sizeof(isal_zstream);
        void* const probe = mmap(nullptr, size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        uintptr_t at = reinterpret_cast<uintptr_t>(probe) & ~(bgzf::detail::kStreamStep - 1);
        munmap(probe, size);
        long scanned = 0;
        auto const t0 = std::chrono::steady_clock::now();
        for (int tried = 0; tried < 4 && at > (uintptr_t{1} << 32); at -= bgzf::detail::kStreamStep) {
            ++scanned;
            if (!bgzf::detail::IsStreamPlace(at)) continue;
            tried++;
            void* const p = mmap(reinterpret_cast<void*>(at), size, PROT_READ | PROT_WRITE,
                                 MAP_PRIVATE | MAP_ANONYMOUS | MAP_FIXED_NOREPLACE, -1, 0);
            int const e = errno;
            double const ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
            std::printf("probe %p: candidate %#lx after %ld steps (%.1f ms): %s\n", probe, at, scanned, ms,
                        p == reinterpret_cast<void*>(at) ? "mapped" : p == MAP_FAILED ? std::strerror(e) : "elsewhere");
            if (p == reinterpret_cast<void*>(at)) break;  // kept mapped, so the next probe goes elsewhere
            if (p != MAP_FAILED) munmap(p, size);
        }
    }
    std::FILE* maps = std::fopen("/proc/self/maps", "r");
    char line[512];
    int lines = 0;
    while (std::fgets(line, sizeof(line), maps)) ++lines;
    std::printf("%d mappings\n", lines);
    return 0;
}
