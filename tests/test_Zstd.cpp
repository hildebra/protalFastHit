// Unit tests for zstd-compressed database files: the compressing and decompressing streams,
// failure modes (truncated, corrupt, trailing data), rewinding, CompressFile, and loading a
// compressed reference.fna.zst with GenomeLoader.
#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <sysexits.h>
#include <unistd.h>
#include "Utilities/Zstd.h"
#include "Hash/Seedmap.h"
#include "SequenceUtils/GenomeLoader.h"

namespace fs = std::filesystem;
using namespace protal;

namespace {
    struct TempDir {
        fs::path dir;
        TempDir() {
            dir = fs::temp_directory_path() / ("protal_zstdtest_" + std::to_string(::getpid()));
            fs::create_directories(dir);
        }
        ~TempDir() { fs::remove_all(dir); }
        std::string operator/(std::string const& name) const { return (dir / name).string(); }
    };

    // Compressible, but not trivially: random DNA with repeats.
    std::string TestData(size_t size, unsigned seed = 1) {
        std::mt19937 rng(seed);
        std::string s;
        s.reserve(size);
        while (s.size() < size) {
            if (s.size() > 1000 && rng() % 4 == 0) {
                size_t const from = rng() % (s.size() - 500);
                s.append(s, from, std::min<size_t>(300, size - s.size()));
            } else {
                s.push_back("ACGT"[rng() % 4]);
            }
        }
        return s;
    }

    void WriteCompressed(std::string const& path, std::string const& data, zstd::Params params = {3, 0, 1},
                         bool pledge = true) {
        zstd::OStream os(path, params, pledge ? std::optional<uint64_t>(data.size()) : std::nullopt);
        // Mixed write sizes: small writes go through the buffer, large ones straight to zstd.
        size_t pos = 0, step = 1;
        while (pos < data.size()) {
            size_t const n = std::min(step, data.size() - pos);
            os.write(data.data() + pos, static_cast<std::streamsize>(n));
            pos += n;
            step = step < (size_t{4} << 20) ? step * 7 : 1;
        }
        ASSERT_TRUE(os.Close()) << os.Buffer().Error();
    }

    std::string ReadAll(zstd::InputFile& in, size_t chunk) {
        std::string out, buffer(chunk, '\0');
        while (in.Stream().read(buffer.data(), static_cast<std::streamsize>(chunk)) || in.Stream().gcount() > 0) {
            out.append(buffer.data(), static_cast<size_t>(in.Stream().gcount()));
            if (!in.Stream()) break;
        }
        return out;
    }

    std::string Slurp(std::string const& path) {
        std::ifstream is(path, std::ios::binary);
        return {std::istreambuf_iterator<char>(is), std::istreambuf_iterator<char>()};
    }

    void Spit(std::string const& path, std::string const& data) {
        std::ofstream os(path, std::ios::binary);
        os << data;
    }
}

TEST(Zstd, RoundTripWithSmallAndLargeReads) {
    TempDir tmp;
    for (size_t size : {size_t{1}, size_t{1000}, size_t{3} << 20, size_t{20} << 20}) {
        std::string const data = TestData(size, static_cast<unsigned>(size));
        WriteCompressed(tmp / "a.zst", data);
        EXPECT_TRUE(zstd::IsCompressed(tmp / "a.zst"));
        for (size_t chunk : {size_t{1}, size_t{4096}, size_t{64} << 20}) {
            if (chunk == 1 && size > 100000) continue;  // byte by byte only for the small ones
            zstd::InputFile in(tmp / "a.zst");
            ASSERT_TRUE(in.IsOpen());
            EXPECT_TRUE(in.Compressed());
            EXPECT_EQ(ReadAll(in, chunk), data) << "size " << size << ", chunk " << chunk;
            EXPECT_FALSE(in.Stream().bad());
        }
    }
}

TEST(Zstd, MultithreadedLongWindowAndUnpledgedSize) {
    TempDir tmp;
    std::string const data = TestData(size_t{24} << 20, 7);
    for (bool pledge : {true, false}) {
        WriteCompressed(tmp / "b.zst", data, {9, 27, 4}, pledge);
        EXPECT_LT(fs::file_size(tmp / "b.zst"), data.size() / 3);
        zstd::InputFile in(tmp / "b.zst");
        EXPECT_EQ(ReadAll(in, size_t{1} << 20), data);
    }
}

TEST(Zstd, TellAndGetline) {
    TempDir tmp;
    std::string const data = ">1_1\nACGT\n>1_2\nGGGCCC\n";
    WriteCompressed(tmp / "c.zst", data);
    zstd::InputFile in(tmp / "c.zst");
    std::string line;
    std::vector<long> positions;
    while (std::getline(in.Stream(), line)) positions.push_back(static_cast<long>(in.Stream().tellg()));
    EXPECT_EQ(positions, (std::vector<long>{5, 10, 15, 22}));
}

TEST(Zstd, RewindRawAndCompressed) {
    TempDir tmp;
    std::string const data = TestData(size_t{10} << 20);
    Spit(tmp / "raw", data);
    WriteCompressed(tmp / "raw.zst", data);
    for (auto name : {"raw", "raw.zst"}) {
        zstd::InputFile in(tmp / name);
        EXPECT_EQ(ReadAll(in, 1 << 16), data);
        ASSERT_TRUE(in.Rewind());
        EXPECT_EQ(ReadAll(in, 1 << 16), data) << name;
    }
}

TEST(Zstd, ConcatenatedFramesAreOneStream) {
    TempDir tmp;
    WriteCompressed(tmp / "x.zst", "first part|");
    WriteCompressed(tmp / "y.zst", "second part");
    Spit(tmp / "xy.zst", Slurp(tmp / "x.zst") + Slurp(tmp / "y.zst"));
    zstd::InputFile in(tmp / "xy.zst");
    EXPECT_EQ(ReadAll(in, 7), "first part|second part");
}

TEST(Zstd, TruncatedCorruptAndTrailingDataFail) {
    TempDir tmp;
    std::string const data = TestData(size_t{2} << 20);
    WriteCompressed(tmp / "good.zst", data);
    std::string const good = Slurp(tmp / "good.zst");

    Spit(tmp / "truncated.zst", good.substr(0, good.size() / 2));
    std::string corrupt = good;
    corrupt[corrupt.size() / 2] ^= 0x5a;
    Spit(tmp / "corrupt.zst", corrupt);
    Spit(tmp / "trailing.zst", good + "garbage");

    for (auto name : {"truncated.zst", "corrupt.zst", "trailing.zst"}) {
        zstd::InputFile in(tmp / name);
        std::string const read = ReadAll(in, 1 << 16);
        EXPECT_TRUE(in.Stream().bad()) << name << " was read without an error";
        EXPECT_NE(read, data + "garbage") << name;
    }
}

TEST(Zstd, PledgedSizeMismatchFailsOnClose) {
    TempDir tmp;
    zstd::OStream os(tmp / "short.zst", {3, 0, 1}, 100);
    os << "only a few bytes";
    EXPECT_FALSE(os.Close());
}

TEST(Zstd, ResolvePrefersTheRawFile) {
    TempDir tmp;
    EXPECT_EQ(zstd::Resolve(tmp / "f"), tmp / "f");  // neither: the raw name
    Spit(tmp / "f.zst", "z");
    EXPECT_EQ(zstd::Resolve(tmp / "f"), tmp / "f.zst");
    Spit(tmp / "f", "r");
    EXPECT_EQ(zstd::Resolve(tmp / "f"), tmp / "f");
}

TEST(Zstd, CompressFileVerifiesAndRenames) {
    TempDir tmp;
    std::string const data = TestData(size_t{5} << 20);
    Spit(tmp / "reference.fna", data);
    std::string error;
    ASSERT_TRUE(zstd::CompressFile(tmp / "reference.fna", tmp / "reference.fna.zst", {3, 27, 2}, true, error)) << error;
    EXPECT_FALSE(fs::exists(tmp / "reference.fna.zst.partial"));
    zstd::InputFile in(tmp / "reference.fna.zst");
    EXPECT_EQ(ReadAll(in, 1 << 20), data);

    EXPECT_FALSE(zstd::CompressFile(tmp / "missing.fna", tmp / "missing.fna.zst", {3, 0, 1}, true, error));
    EXPECT_FALSE(fs::exists(tmp / "missing.fna.zst"));
}

// The index loader stops with exit 8 on a compressed index that holds too little data, and on
// one whose zstd frame is cut off.
TEST(Zstd, TruncatedCompressedIndexExits8) {
    TempDir tmp;
    std::ostringstream header;
    auto put = [&](auto v) { header.write(reinterpret_cast<char const*>(&v), sizeof(v)); };
    constexpr size_t keys = size_t{1} << 30, total = keys + ((keys / 8) + 1) * 4;
    put(Seedmap::kFileMagic);
    put(Seedmap::kIndexVersionMajor);
    put(Seedmap::kIndexVersionMinor);
    put(Seedmap::kIndexVersionPatch);
    for (size_t v : {keys, total, size_t{8}, size_t{8}, size_t{3}, size_t{0}}) put(v);
    WriteCompressed(tmp / "short.prx.zst", header.str() + TestData(size_t{1} << 20));
    EXPECT_EXIT({ Seedmap map; map.Load(tmp / "short.prx.zst"); }, testing::ExitedWithCode(8), "truncated or corrupt");

    std::string const frame = Slurp(tmp / "short.prx.zst");
    Spit(tmp / "cut.prx.zst", frame.substr(0, frame.size() / 2));
    EXPECT_EXIT({ Seedmap map; map.Load(tmp / "cut.prx.zst"); }, testing::ExitedWithCode(8), "truncated or corrupt");
}

// reference.fna.zst: preloading gives exactly the genes of reference.fna; loading a single gene
// on demand is refused.
TEST(Zstd, GenomeLoaderReadsACompressedReference) {
    TempDir tmp;
    std::vector<std::pair<std::string, std::string>> records = {
            {"1_1", TestData(300, 1)}, {"1_2", TestData(1200, 2)}, {"2_1", TestData(80, 3)}, {"10_3", TestData(5000, 4)}};
    {
        std::ofstream fna(tmp / "reference.fna", std::ios::binary);
        std::ofstream map(tmp / "reference.map");
        size_t offset = 0;
        for (auto const& [name, seq] : records) {
            std::string const header = ">" + name + "\n";
            fna << header << seq << '\n';
            offset += header.size();
            map << name.substr(0, name.find('_')) << '\t' << name.substr(name.find('_') + 1) << '\t'
                << offset << '\t' << offset + seq.size() << '\n';
            offset += seq.size() + 1;
        }
    }
    std::string error;
    ASSERT_TRUE(zstd::CompressFile(tmp / "reference.fna", tmp / "reference.fna.zst", {19, 27, 1}, true, error)) << error;

    GenomeLoader raw(tmp / "reference.fna", tmp / "reference.map");
    GenomeLoader compressed(tmp / "reference.fna.zst", tmp / "reference.map");
    EXPECT_FALSE(raw.IsCompressed());
    EXPECT_TRUE(compressed.IsCompressed());
    raw.LoadAllGenomes();
    compressed.LoadAllGenomes();
    EXPECT_TRUE(compressed.AllGenomesLoaded());
    for (auto const& [name, seq] : records) {
        size_t const taxid = std::stoul(name.substr(0, name.find('_')));
        size_t const gene = std::stoul(name.substr(name.find('_') + 1));
        EXPECT_EQ(raw.GetGenome(taxid).GetGene(gene).Sequence(), seq) << name;
        EXPECT_EQ(compressed.GetGenome(taxid).GetGene(gene).Sequence(), seq) << name;
    }

    GenomeLoader lazy(tmp / "reference.fna.zst", tmp / "reference.map");
    EXPECT_EXIT(lazy.GetGenome(1).GetGeneOMP(1), testing::ExitedWithCode(EX_SOFTWARE), "must be preloaded");
}
