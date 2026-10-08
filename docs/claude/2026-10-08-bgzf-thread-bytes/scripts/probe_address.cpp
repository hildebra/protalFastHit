// Probe: ISA-L's level-1 stateless deflate of one input, with its isal_zstream at different addresses (64 KB apart).
// Prints how many distinct outputs there are, and whether the outputs of the failing run's two files appear.
// usage: probe_address CONTENT OFFSET LENGTH BGZF_ONE BGZF_FOUR BLOCK_OFFSET [COUNT]
#include <isa-l/igzip_lib.h>
#include <sys/mman.h>

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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

static std::string Slurp(char const* path) {
    std::ifstream f(path, std::ios::binary);
    return std::string((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}

// The deflate payload of the BGZF block at `at`.
static std::string Payload(std::string const& bgzf, size_t at) {
    size_t const size = (static_cast<unsigned char>(bgzf[at + 16]) | static_cast<unsigned char>(bgzf[at + 17]) << 8) + 1;
    return bgzf.substr(at + 18, size - 26);
}

int main(int argc, char** argv) {
    if (argc < 7) return 2;
    std::string const all = Slurp(argv[1]);
    size_t const offset = std::strtoull(argv[2], nullptr, 10), length = std::strtoull(argv[3], nullptr, 10);
    size_t const at = std::strtoull(argv[6], nullptr, 10);
    int const count = argc > 7 ? std::atoi(argv[7]) : 8192;
    std::string const one = Payload(Slurp(argv[4]), at), four = Payload(Slurp(argv[5]), at);
    std::uint64_t const h_one = Fnv(reinterpret_cast<unsigned char const*>(one.data()), one.size());
    std::uint64_t const h_four = Fnv(reinterpret_cast<unsigned char const*>(four.data()), four.size());
    std::printf("one-thread payload %zu bytes %016llx, four-thread %zu bytes %016llx\n", one.size(),
                static_cast<unsigned long long>(h_one), four.size(), static_cast<unsigned long long>(h_four));

    std::string const block = all.substr(offset, length);
    std::vector<uint8_t> level_buffer(ISAL_DEF_LVL1_DEFAULT);
    std::vector<unsigned char> out(0x10000);
    size_t const region_size = static_cast<size_t>(count) * 0x10000 + 0x40000;
    auto* region = static_cast<unsigned char*>(
        mmap(nullptr, region_size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0));
    if (region == MAP_FAILED) return 3;
    std::map<std::uint64_t, std::pair<size_t, int>> outputs;  // checksum -> size, addresses
    std::vector<std::uintptr_t> one_at;
    for (int k = 0; k < count; ++k) {
        void* place = region + static_cast<size_t>(k) * 0x10000;
        auto* stream = new (place) isal_zstream();  // value-initialized, as std::make_unique<isal_zstream>()
        isal_deflate_stateless_init(stream);
        stream->level = 1;
        stream->level_buf = level_buffer.data();
        stream->level_buf_size = static_cast<uint32_t>(level_buffer.size());
        stream->end_of_stream = 1;
        stream->flush = NO_FLUSH;
        stream->next_in = reinterpret_cast<uint8_t*>(const_cast<char*>(block.data()));
        stream->avail_in = static_cast<uint32_t>(block.size());
        stream->next_out = out.data();
        stream->avail_out = 0x10000 - 26;
        if (isal_deflate_stateless(stream) != COMP_OK) return 4;
        std::uint64_t const h = Fnv(out.data(), stream->total_out);
        auto& o = outputs[h];
        o.first = stream->total_out;
        ++o.second;
        if (h == h_one) one_at.push_back(reinterpret_cast<std::uintptr_t>(place));
        madvise(place, 0x20000, MADV_DONTNEED);
    }
    std::printf("%d addresses, %zu distinct outputs:\n", count, outputs.size());
    for (auto const& [h, o] : outputs) {
        std::printf("  %016llx %zu bytes x%d%s\n", static_cast<unsigned long long>(h), o.first, o.second,
                    h == h_one ? "  <- the one-thread run's" : h == h_four ? "  <- the four-thread runs'" : "");
    }
    for (size_t i = 0; i < one_at.size() && i < 5; ++i) std::printf("  one-thread output at %#lx\n", one_at[i]);
    return 0;
}
