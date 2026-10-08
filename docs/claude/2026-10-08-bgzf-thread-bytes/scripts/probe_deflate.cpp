// Probe: does ISA-L's stateless level-1 deflate (protal::bgzf::Deflater) of the same input depend on anything but the
// input? Compresses [offset, offset + length) of a file under several conditions and prints the compressed size and a
// checksum of the output.
#include "IO/Bgzf.h"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iterator>
#include <map>
#include <random>
#include <string>
#include <vector>

static std::uint64_t Fnv(unsigned char const* p, size_t n) {
    std::uint64_t h = 1469598103934665603ULL;
    for (size_t i = 0; i < n; ++i) h = (h ^ p[i]) * 1099511628211ULL;
    return h;
}

// Compresses `size` bytes at `in` with deflater d; prints the size and checksum.
static size_t Run(protal::bgzf::Deflater& d, char const* in, size_t size, char const* label) {
    std::vector<unsigned char> out(0x10000);
    size_t const n = d.Compress(in, size, out.data(), out.size() - 26);
    std::printf("%-58s %6zu bytes  %016llx\n", label, n, static_cast<unsigned long long>(Fnv(out.data(), n)));
    return n;
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::fprintf(stderr, "usage: probe FILE OFFSET LENGTH\n");
        return 2;
    }
    std::ifstream f(argv[1], std::ios::binary);
    std::string const all((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    size_t const offset = std::strtoull(argv[2], nullptr, 10), length = std::strtoull(argv[3], nullptr, 10);
    std::string const block = all.substr(offset, length);

    // 1. a fresh deflater, the input at the end of a buffer of exactly its size (valgrind sees any read beyond it)
    {
        protal::bgzf::Deflater d;
        char* exact = static_cast<char*>(std::malloc(length));
        std::memcpy(exact, block.data(), length);
        Run(d, exact, length, "fresh, exact-size buffer");
        std::free(exact);
    }
    // 2. what follows the input in its buffer
    std::vector<std::pair<std::string, std::string>> tails = {
        {"zeros", std::string(4096, '\0')},
        {"0xff", std::string(4096, '\xff')},
        {"the file's next bytes", all.substr(offset + length, 4096)},
        {"the input again", block.substr(0, std::min<size_t>(4096, length))},
        {"its last 300 bytes again", block.substr(length - 300) + std::string(3796, 'A')},
    };
    std::mt19937 gen(1);
    for (int r = 0; r < 3; ++r) {
        std::string t(4096, 0);
        for (auto& c : t) c = "ACGT\n+@"[gen() % 7];
        tails.push_back({"random FASTQ-like letters " + std::to_string(r), t});
    }
    for (auto const& [name, tail] : tails) {
        protal::bgzf::Deflater d;
        std::string buf = block + tail;
        Run(d, buf.data(), length, ("fresh, followed by " + name).c_str());
    }
    // 3. the same deflater after other blocks (the file's blocks before it, at 65280 each)
    {
        protal::bgzf::Deflater d;
        for (size_t at = 0; at + 0xff00 <= offset; at += 0xff00) {
            std::vector<unsigned char> out(0x10000);
            d.Compress(all.data() + at, 0xff00, out.data(), out.size() - 26);
        }
        std::string buf = block + std::string(4096, '\0');
        Run(d, buf.data(), length, "after the file's blocks, followed by zeros");
    }
    // 4. a search: tails of the input's own bytes repeated from each of its last 300 offsets
    std::map<size_t, int> sizes;
    for (size_t back = 1; back <= 600 && back <= length; ++back) {
        protal::bgzf::Deflater d;
        std::string buf = block + block.substr(length - back) + std::string(4096, 'N');
        std::vector<unsigned char> out(0x10000);
        size_t const n = d.Compress(buf.data(), length, out.data(), out.size() - 26);
        ++sizes[n];
    }
    std::printf("input followed by its own last 1-600 bytes:");
    for (auto const& [n, count] : sizes) std::printf("  %zu bytes x%d", n, count);
    std::printf("\n");
    return 0;
}
