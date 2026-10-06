// make_bigdb.cpp - a large single-file database for timing protal --add_model: the members of an existing
// database.protal (their frames copied), a padding member padding.bin of PAD_GB GB of incompressible data in 64 MB
// frames after them, and the models of the existing database last.
//
//   make_bigdb SOURCE_DB OUT_DB PAD_GB THREADS
#include "Utilities/Database.h"
#include "ReadType.h"

#include <iostream>

using namespace protal;

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "usage: make_bigdb SOURCE_DB OUT_DB PAD_GB THREADS" << std::endl;
        return 1;
    }
    std::string const source = argv[1], out = argv[2];
    uint64_t const pad_gb = std::stoull(argv[3]);
    int const threads = std::stoi(argv[4]);
    std::string error;
    auto const bundle = db::Bundle::Open(source, error);
    if (!bundle) {
        std::cerr << source << ": " << error << std::endl;
        return 1;
    }
    std::string const pad = out + ".padding.zst";
    uint64_t constexpr kFrame = uint64_t{64} << 20;
    zstd::Params const params{1, 0, threads, kFrame};
    auto const written = zstd::WriteSeekable(pad, pad_gb * 16, params, [](size_t index, std::vector<char>& frame) {
        frame.resize(kFrame);
        uint64_t x = 0x9e3779b97f4a7c15ull * (index + 1);
        for (size_t i = 0; i + 8 <= frame.size(); i += 8) {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            std::memcpy(frame.data() + i, &x, 8);
        }
        return std::string();
    }, error);
    if (!written) {
        std::cerr << "padding: " << error << std::endl;
        return 1;
    }
    auto const models = AllModelFiles();
    std::vector<db::Source> sources, last;
    for (auto const& member : bundle->Members()) {
        bool const model = std::find(models.begin(), models.end(), member.name) != models.end();
        (model ? last : sources).push_back({member.name, source, member.frames});
    }
    sources.push_back({"padding.bin", pad});
    sources.insert(sources.end(), last.begin(), last.end());
    auto const size = db::Write(out, sources, params, error);
    std::filesystem::remove(pad);
    if (!size) {
        std::cerr << out << ": " << error << std::endl;
        return 1;
    }
    std::cout << out << ": " << *size << " bytes, " << sources.size() << " members" << std::endl;
    return 0;
}
