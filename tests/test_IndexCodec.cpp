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
    // per-key offsets are prefix sums, the last control block holds the number of values. The keys
    // of the first single_cell_blocks blocks have 0 or 1 entries, so those blocks have no flex cells.
    // With empty_runs, the blocks of every second run of that many (from block empty_runs on) are empty.
    struct SmallIndex {
        index_codec::Layout layout;
        std::vector<uint16_t> keymap;
        std::vector<uint64_t> values;
        std::string header = std::string("an index header\0with bytes", 26);

        SmallIndex(uint64_t blocks, unsigned seed, double empty_share, uint64_t single_cell_blocks = 0, uint64_t empty_runs = 0) {
            std::mt19937_64 rng(seed);
            layout.blocks = blocks;
            keymap.assign(layout.KeymapCells(), 0);
            uint64_t v = 0;
            for (uint64_t b = 0; b < blocks; b++) {
                uint16_t* cells = keymap.data() + b * layout.CellsPerBlock();
                bool const empty_block = empty_runs > 0 && (b / empty_runs) % 2 == 1;
                std::memcpy(cells, &v, 8);
                uint64_t local = 0;
                for (uint64_t j = 0; j < layout.keys_per_block; j++) {
                    cells[4 + j] = static_cast<uint16_t>(local);
                    std::uniform_real_distribution<double> u(0, 1);
                    uint64_t const c = empty_block || u(rng) < empty_share ? 0 : b < single_cell_blocks ? 1 : 1 + rng() % (rng() % 4 == 0 ? 40 : 3);
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

    // What the header of a chunk (see IndexCodec.h) says about it.
    struct ChunkInfo {
        uint8_t mode = 0;
        uint64_t entries = 0, flex = 0;
    };

    std::vector<ChunkInfo> ChunkInfos(std::string const& path) {
        std::string error;
        auto table = zstd::ReadSeekTable(path, error);
        std::vector<ChunkInfo> infos;
        if (!table) return infos;
        infos.resize(table->frames.size() - 1);
        zstd::ForEachFrame(path, *table, 1, 1, [&](size_t frame, char const* data, size_t, size_t) -> std::string {
            ChunkInfo& info = infos[frame - 1];
            info.mode = static_cast<uint8_t>(data[0]);
            std::memcpy(&info.entries, data + 8, 8);
            std::memcpy(&info.flex, data + 16, 8);
            return "";
        });
        return infos;
    }

    // The split form of chunk c as the format describes it (IndexCodec.h), written plainly from copies of the
    // chunk's cells: header, bitmap, then byte planes of the counts (2), flex cells (8), taxids, gene ids and
    // positions (the widths their maxima need) and flags (1).
    std::string ReferenceSplit(SmallIndex const& index, index_codec::Chunk const& c) {
        auto const& l = index.layout;
        std::string bitmap((c.blocks + 7) / 8, '\0');
        std::vector<uint64_t> counts, flex, entries;
        for (uint64_t b = c.first_block; b < c.first_block + c.blocks; b++) {
            uint16_t const* cells = index.keymap.data() + b * l.CellsPerBlock();
            uint64_t start = 0, end = 0;
            std::memcpy(&start, cells, 8);
            std::memcpy(&end, cells + l.CellsPerBlock(), 8);
            if (end == start) continue;
            bitmap[(b - c.first_block) >> 3] = static_cast<char>(bitmap[(b - c.first_block) >> 3] | 1 << ((b - c.first_block) & 7));
            for (uint64_t j = 0; j < l.keys_per_block; j++) {
                uint64_t const offset = cells[l.ctrl_cells + j], next = j + 1 < l.keys_per_block ? cells[l.ctrl_cells + j + 1] : end - start;
                uint64_t const n = next - offset, f = n >= 2 ? (n + 2) / 3 : 0;
                counts.push_back(n);
                auto const* v = index.values.data() + start + offset;
                flex.insert(flex.end(), v, v + f);
                entries.insert(entries.end(), v + f, v + n);
            }
        }
        auto field = [](uint64_t e, int shift) { return (e >> shift) & 0xFFFFF; };
        auto width = [&](int shift) {
            uint64_t max = 0;
            for (uint64_t e : entries) max = std::max(max, field(e, shift));
            return max == 0 ? 0 : max < 256 ? 1 : max < 65536 ? 2 : 3;
        };
        int const tw = width(40), gw = width(20), pw = width(0);
        std::string out(index_codec::kChunkHeaderBytes, '\0');
        out[0] = static_cast<char>(index_codec::kModeSplit);
        out[1] = static_cast<char>(tw);
        out[2] = static_cast<char>(gw);
        out[3] = static_cast<char>(pw);
        uint64_t const ne = entries.size(), nf = flex.size(), filled = counts.size() / l.keys_per_block;
        std::memcpy(out.data() + 8, &ne, 8);
        std::memcpy(out.data() + 16, &nf, 8);
        std::memcpy(out.data() + 24, &filled, 8);
        out += bitmap;
        auto planes = [&out](std::vector<uint64_t> const& cells, int bytes, auto&& value) {
            for (int b = 0; b < bytes; b++) {
                for (uint64_t cell : cells) out.push_back(static_cast<char>(value(cell) >> (8 * b)));
            }
        };
        auto same = [](uint64_t x) { return x; };
        planes(counts, 2, same);
        planes(flex, 8, same);
        planes(entries, tw, [&](uint64_t e) { return field(e, 40); });
        planes(entries, gw, [&](uint64_t e) { return field(e, 20); });
        planes(entries, pw, [&](uint64_t e) { return field(e, 0); });
        planes(entries, 1, [](uint64_t e) { return e >> 60; });
        return out;
    }

    // Writes index with chunks of about chunk_bytes, then checks that every chunk is stored split,
    // that it reads back as it was written with any thread count, and that Verify accepts it.
    // Returns the chunks' headers.
    std::vector<ChunkInfo> RoundTrip(SmallIndex const& index, std::string const& path, uint64_t chunk_bytes) {
        size_t raw_chunks = 99;
        std::string error;
        auto written = index_codec::Write(path, index.header, index.layout, index.keymap.data(), index.values.data(),
                                          {3, 0, 4, chunk_bytes}, error, &raw_chunks);
        EXPECT_TRUE(written) << error;
        EXPECT_EQ(raw_chunks, 0u) << "chunks of " << chunk_bytes;
        EXPECT_EQ(index_codec::Verify(path, index.header, index.layout, index.keymap.data(), index.values.data(), 3), "");
        for (int threads : {1, 4}) {
            std::vector<uint16_t> km;
            std::vector<uint64_t> vals;
            EXPECT_EQ(DecodeFile(path, threads, km, vals), "") << "chunks of " << chunk_bytes << ", " << threads << " threads";
            EXPECT_EQ(km, index.keymap) << "chunks of " << chunk_bytes << ", " << threads << " threads";
            EXPECT_EQ(vals, index.values) << "chunks of " << chunk_bytes << ", " << threads << " threads";
        }
        return ChunkInfos(path);
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

// Runs of empty blocks shorter and longer than the 64 the decoder fills at once (from a chunk's every 64th block),
// at many offsets to the chunks' starts, between filled blocks.
TEST(IndexCodec, LongRunsOfEmptyBlocks) {
    TempDir tmp;
    for (uint64_t run : {63, 64, 65, 129, 300}) {
        SmallIndex index(5000, static_cast<unsigned>(run), 0.3, 0, run);
        for (uint64_t chunk : {uint64_t{3000}, uint64_t{40000}, uint64_t{4} << 20}) {
            RoundTrip(index, tmp / ("runs" + std::to_string(run) + "_" + std::to_string(chunk) + ".zst"), chunk);
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
    // No value cells at all: Verify decodes into empty value buffers, which have no data pointer.
    EXPECT_EQ(index_codec::Verify(tmp / "index.zst", index.header, index.layout, index.keymap.data(), index.values.data(), 2), "");
}

// Keys with 0 or 1 cells have no flex cells; a chunk of such keys has an empty flex column, and the
// decoder must not copy from the (null) data pointer of an empty buffer.
TEST(IndexCodec, ChunksWithoutFlexCells) {
    TempDir tmp;
    SmallIndex index(3000, 7, 0.3, 3000);
    ASSERT_FALSE(index.values.empty());
    for (uint64_t chunk : {uint64_t{4096}, uint64_t{50000}, uint64_t{64} << 20}) {
        auto const infos = RoundTrip(index, tmp / "index.zst", chunk);
        ASSERT_FALSE(infos.empty());
        if (chunk == 4096) EXPECT_GT(infos.size(), 10u);
        for (auto const& info : infos) {
            EXPECT_EQ(info.mode, index_codec::kModeSplit);
            EXPECT_EQ(info.flex, 0u) << "chunks of " << chunk;
            EXPECT_GT(info.entries, 0u) << "chunks of " << chunk;
        }
    }
}

// A chunk that starts with keys of one cell and has flex keys only after its first batch of cells
// (DecodeChunk composes the cells in batches of 65536): the flex cells still land at the right keys.
TEST(IndexCodec, FlexCellsAfterAFirstBatchWithoutAny) {
    TempDir tmp;
    SmallIndex index(22000, 8, 0.2, 20000);
    auto const infos = RoundTrip(index, tmp / "index.zst", uint64_t{64} << 20);
    ASSERT_EQ(infos.size(), 1u);
    EXPECT_EQ(infos[0].mode, index_codec::kModeSplit);
    EXPECT_GT(infos[0].flex, 0u);
    EXPECT_GT(infos[0].entries, uint64_t{65536});
}

// The writer composes the planes straight from the index's cells (no copies of them): its chunks are, byte for
// byte, the format's plain description of them, as the copying writer wrote them before.
TEST(IndexCodec, ChunksAreTheFormatsBytes) {
    TempDir tmp;
    for (auto const& index : { SmallIndex(3000, 11, 0.5), SmallIndex(22000, 8, 0.2, 20000), SmallIndex(5000, 65, 0.3, 0, 65) }) {
        for (uint64_t chunk : { uint64_t{4096}, uint64_t{50000}, uint64_t{64} << 20 }) {
            std::string error;
            ASSERT_TRUE(index_codec::Write(tmp / "index.zst", index.header, index.layout, index.keymap.data(), index.values.data(),
                                           {3, 0, 4, chunk}, error)) << error;
            auto const table = zstd::ReadSeekTable(tmp / "index.zst", error);
            ASSERT_TRUE(table) << error;
            auto const container = index_codec::ReadContainer(tmp / "index.zst", *table, error);
            ASSERT_TRUE(container) << error;
            size_t compared = 0;
            std::string const e = zstd::ForEachFrame(tmp / "index.zst", *table, 1, 1, [&](size_t frame, char const* data, size_t size, size_t) {
                std::string const expected = ReferenceSplit(index, container->chunks[frame - 1]);
                compared++;
                return std::string(data, size) == expected ? std::string() : "chunk " + std::to_string(frame) + " differs";
            });
            EXPECT_EQ(e, "") << "chunks of " << chunk;
            EXPECT_EQ(compared, container->chunks.size());
        }
    }
}

// Verify compares the values as it decodes them: one changed flex cell, entry or raw cell is found.
TEST(IndexCodec, VerifyFindsOneChangedValue) {
    TempDir tmp;
    SmallIndex index(2000, 9, 0.3);
    std::string error;
    ASSERT_TRUE(index_codec::Write(tmp / "index.zst", index.header, index.layout, index.keymap.data(), index.values.data(),
                                   {3, 0, 2, 8192}, error)) << error;
    ASSERT_EQ(index_codec::Verify(tmp / "index.zst", index.header, index.layout, index.keymap.data(), index.values.data(), 2), "");
    for (size_t cell : { size_t{0}, index.values.size() / 2, index.values.size() - 1 }) {
        auto changed = index.values;
        changed[cell] ^= uint64_t{1} << 33;
        std::string const e = index_codec::Verify(tmp / "index.zst", index.header, index.layout, index.keymap.data(), changed.data(), 2);
        EXPECT_NE(e.find("differs from the index"), std::string::npos) << "cell " << cell << ": " << e;
    }
    // A chunk kept as raw cells (its first key does not start at its block).
    SmallIndex odd(500, 3, 0.2);
    odd.keymap[4] = 1;
    size_t raw_chunks = 0;
    ASSERT_TRUE(index_codec::Write(tmp / "odd.zst", odd.header, odd.layout, odd.keymap.data(), odd.values.data(), {3, 0, 2, 2048}, error,
                                   &raw_chunks)) << error;
    ASSERT_EQ(raw_chunks, 1u);
    ASSERT_EQ(index_codec::Verify(tmp / "odd.zst", odd.header, odd.layout, odd.keymap.data(), odd.values.data(), 2), "");
    auto changed = odd.values;
    changed[1] ^= 1;
    EXPECT_NE(index_codec::Verify(tmp / "odd.zst", odd.header, odd.layout, odd.keymap.data(), changed.data(), 2).find("differs from the index"),
              std::string::npos);
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
