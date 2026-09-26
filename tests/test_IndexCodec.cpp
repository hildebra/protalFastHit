// Unit tests for the column format of the compressed index (Hash/IndexCodec.h), on small indexes
// laid out as Seedmap::BuildValuePointers lays them out.
#include <gtest/gtest.h>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <vector>
#include <unistd.h>
#include "Hash/IndexCodec.h"

namespace fs = std::filesystem;
using namespace protal;

namespace {
    struct TempDir {
        fs::path dir;
        TempDir() {
            dir = fs::temp_directory_path() / ("protal_codectest_" + std::to_string(::getpid()));
            fs::create_directories(dir);
        }
        ~TempDir() { fs::remove_all(dir); }
        std::string operator/(std::string const& name) const { return (dir / name).string(); }
    };

    // Keys with c entries take c cells, or c + ceil(c/2) with flex keys (c >= 2); block offsets and
    // per-key offsets are prefix sums, the last control block holds the number of values.
    struct SmallIndex {
        index_codec::Layout layout;
        std::vector<uint16_t> keymap;
        std::vector<uint64_t> values;
        std::string header = std::string("an index header\0with bytes", 26);

        SmallIndex(uint64_t blocks, unsigned seed, double empty_share) {
            std::mt19937_64 rng(seed);
            layout.blocks = blocks;
            keymap.assign(layout.KeymapCells(), 0);
            uint64_t v = 0;
            for (uint64_t b = 0; b < blocks; b++) {
                uint16_t* cells = keymap.data() + b * layout.CellsPerBlock();
                std::memcpy(cells, &v, 8);
                uint64_t local = 0;
                for (uint64_t j = 0; j < layout.keys_per_block; j++) {
                    cells[4 + j] = static_cast<uint16_t>(local);
                    std::uniform_real_distribution<double> u(0, 1);
                    uint64_t const c = u(rng) < empty_share ? 0 : 1 + rng() % (rng() % 4 == 0 ? 40 : 3);
                    uint64_t const n = c >= 2 ? c + (c + 1) / 2 : c;
                    uint64_t const taxid = 1 + rng() % 5000;
                    for (uint64_t i = 0; i < n - c; i++) values.push_back(rng());  // flex cells
                    for (uint64_t i = 0; i < c; i++) {
                        uint64_t const t = taxid + rng() % 20, g = 1 + rng() % 168, p = rng() % 4000, f = rng() % 4;
                        values.push_back(t << 40 | g << 20 | p | f << 60);
                    }
                    local += n;
                }
                v += local;
            }
            std::memcpy(keymap.data() + blocks * layout.CellsPerBlock(), &v, 8);
            layout.values = values.size();
        }
    };

    std::string Slurp(std::string const& path) {
        std::ifstream is(path, std::ios::binary);
        return {std::istreambuf_iterator<char>(is), std::istreambuf_iterator<char>()};
    }

    void Spit(std::string const& path, std::string const& data) {
        std::ofstream os(path, std::ios::binary);
        os << data;
    }

    // Decodes path into fresh arrays; returns the error message.
    std::string DecodeFile(std::string const& path, int threads, std::vector<uint16_t>& km, std::vector<uint64_t>& vals,
                           index_codec::Container* out = nullptr) {
        std::string error;
        auto table = zstd::ReadSeekTable(path, error);
        if (!table) return "no seek table " + error;
        auto container = index_codec::ReadContainer(path, *table, error);
        if (!container) return "no container " + error;
        km.assign(container->layout.KeymapCells(), 0xabcd);
        vals.assign(container->layout.values, 0xabcdef);
        if (out) *out = *container;
        return index_codec::Decode(path, *table, *container, km.data(), vals.data(), threads);
    }
}

TEST(IndexCodec, RoundTripWithManyChunksAndAnyThreadCount) {
    TempDir tmp;
    SmallIndex index(3000, 1, 0.5);
    for (uint64_t chunk : {uint64_t{4096}, uint64_t{50000}, uint64_t{64} << 20}) {
        size_t raw_chunks = 99;
        std::string error;
        auto written = index_codec::Write(tmp / "index.zst", index.header, index.layout, index.keymap.data(),
                                          index.values.data(), {3, 0, 4, chunk}, error, &raw_chunks);
        ASSERT_TRUE(written) << error;
        EXPECT_EQ(raw_chunks, 0u);
        EXPECT_TRUE(index_codec::IsSplitIndex(tmp / "index.zst"));
        EXPECT_EQ(index_codec::Verify(tmp / "index.zst", index.header, index.layout, index.keymap.data(), index.values.data(), 3), "");
        for (int threads : {1, 4}) {
            std::vector<uint16_t> km;
            std::vector<uint64_t> vals;
            index_codec::Container container;
            ASSERT_EQ(DecodeFile(tmp / "index.zst", threads, km, vals, &container), "") << chunk;
            EXPECT_EQ(container.index_header, index.header);
            EXPECT_EQ(container.layout.blocks, index.layout.blocks);
            if (chunk == 4096) EXPECT_GT(container.chunks.size(), 10u);
            EXPECT_EQ(km, index.keymap) << "chunks of " << chunk << ", " << threads << " threads";
            EXPECT_EQ(vals, index.values) << "chunks of " << chunk << ", " << threads << " threads";
        }
    }
}

TEST(IndexCodec, ColumnsCompressBetterThanRawCells) {
    TempDir tmp;
    SmallIndex index(20000, 2, 0.3);
    std::string error;
    auto columns = index_codec::Write(tmp / "index.zst", index.header, index.layout, index.keymap.data(),
                                      index.values.data(), {19, 0, 4, uint64_t{1} << 20}, error);
    ASSERT_TRUE(columns) << error;
    std::string raw(reinterpret_cast<char const*>(index.keymap.data()), index.keymap.size() * 2);
    raw.append(reinterpret_cast<char const*>(index.values.data()), index.values.size() * 8);
    std::vector<char> compressed(ZSTD_compressBound(raw.size()));
    size_t const plain = ZSTD_compress(compressed.data(), compressed.size(), raw.data(), raw.size(), 19);
    EXPECT_LT(*columns, plain * 9 / 10) << "columns " << *columns << " bytes, raw cells " << plain;
}

TEST(IndexCodec, NonCanonicalBlocksAreKeptAsRawCells) {
    TempDir tmp;
    SmallIndex index(500, 3, 0.2);
    index.keymap[100 * index.layout.CellsPerBlock() + 4] = 1;  // key 0 of block 100 does not start at the block
    size_t raw_chunks = 0;
    std::string error;
    ASSERT_TRUE(index_codec::Write(tmp / "index.zst", index.header, index.layout, index.keymap.data(), index.values.data(),
                                   {3, 0, 2, 2048}, error, &raw_chunks)) << error;
    EXPECT_EQ(raw_chunks, 1u);
    std::vector<uint16_t> km;
    std::vector<uint64_t> vals;
    ASSERT_EQ(DecodeFile(tmp / "index.zst", 3, km, vals), "");
    EXPECT_EQ(km, index.keymap);
    EXPECT_EQ(vals, index.values);
}

TEST(IndexCodec, EmptyIndex) {
    TempDir tmp;
    SmallIndex index(64, 4, 1.0);
    ASSERT_TRUE(index.values.empty());
    std::string error;
    ASSERT_TRUE(index_codec::Write(tmp / "index.zst", index.header, index.layout, index.keymap.data(), index.values.data(),
                                   {3, 0, 2, 256}, error)) << error;
    std::vector<uint16_t> km;
    std::vector<uint64_t> vals;
    ASSERT_EQ(DecodeFile(tmp / "index.zst", 2, km, vals), "");
    EXPECT_EQ(km, index.keymap);
}

TEST(IndexCodec, CorruptAndTruncatedFilesFail) {
    TempDir tmp;
    SmallIndex index(2000, 5, 0.4);
    std::string error;
    ASSERT_TRUE(index_codec::Write(tmp / "index.zst", index.header, index.layout, index.keymap.data(), index.values.data(),
                                   {3, 0, 2, 8192}, error)) << error;
    std::string const good = Slurp(tmp / "index.zst");
    auto table = zstd::ReadSeekTable(tmp / "index.zst", error);
    ASSERT_TRUE(table);

    std::string corrupt = good;  // inside chunk 3
    corrupt[table->frames[3].compressed_offset + table->frames[3].compressed_size / 2] ^= 0x5a;
    Spit(tmp / "corrupt.zst", corrupt);
    std::vector<uint16_t> km;
    std::vector<uint64_t> vals;
    EXPECT_NE(DecodeFile(tmp / "corrupt.zst", 2, km, vals), "");
    EXPECT_NE(index_codec::Verify(tmp / "corrupt.zst", index.header, index.layout, index.keymap.data(), index.values.data(), 2), "");

    Spit(tmp / "truncated.zst", good.substr(0, good.size() - 30));  // the seek table is gone
    EXPECT_FALSE(index_codec::IsSplitIndex(tmp / "truncated.zst"));
    EXPECT_TRUE(index_codec::StartsWithMagic(tmp / "truncated.zst"));

    // A different index does not verify.
    SmallIndex other(2000, 6, 0.4);
    EXPECT_NE(index_codec::Verify(tmp / "index.zst", index.header, other.layout, other.keymap.data(), other.values.data(), 2), "");
}
