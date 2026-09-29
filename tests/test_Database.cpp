// Unit tests for the single-file database (Utilities/Database.h): writing it from seekable and other
// files, reading members sequentially and in parallel, the index's column chunks and the reference
// read through member frames, where --db points (Locate), and failures on truncated or corrupt files.
#include <gtest/gtest.h>
#include <cctype>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <vector>
#include <unistd.h>
#include "Utilities/Database.h"
#include "Utilities/ReferenceFingerprint.h"
#include "Hash/IndexCodec.h"
#include "SequenceUtils/GenomeLoader.h"
#include "Taxonomy/Taxonomy.h"

namespace fs = std::filesystem;
using namespace protal;

namespace {
    struct TempDir {
        fs::path dir;
        TempDir() {
            dir = fs::temp_directory_path() / ("protal_dbtest_" + std::to_string(::getpid()));
            fs::remove_all(dir);
            fs::create_directories(dir);
        }
        ~TempDir() { fs::remove_all(dir); }
        std::string operator/(std::string const& name) const { return (dir / name).string(); }
    };

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

    std::string Slurp(std::string const& path) {
        std::ifstream is(path, std::ios::binary);
        return {std::istreambuf_iterator<char>(is), std::istreambuf_iterator<char>()};
    }

    std::string Spit(std::string const& path, std::string const& data) {
        std::ofstream os(path, std::ios::binary);
        os << data;
        return path;
    }

    // Collects what ParallelRead delivers.
    struct StringSink : zstd::Sink {
        std::string data;
        explicit StringSink(size_t size) : data(size, '\0') {}
        void Copy(uint64_t offset, char const* bytes, size_t size) override {
            ASSERT_LE(offset + size, data.size());
            std::memcpy(data.data() + offset, bytes, size);
        }
    };

    std::string Content(db::DbFile const& file) {
        std::string error;
        auto content = file.ReadAll(error);
        EXPECT_TRUE(content.has_value()) << file.Name() << ": " << error;
        return content.value_or("");
    }

    zstd::Params SmallFrames(uint64_t frame_size, int threads = 2) {
        return {3, 0, threads, frame_size};
    }
}

// Seekable sources are copied frame by frame, the others compressed into frames; every member reads
// back as its source, sequentially and with any number of threads.
TEST(Database, MembersReadBackAsTheirSources) {
    TempDir tmp;
    std::string const seekable_data = TestData(20000, 1), raw_data = TestData(9000, 2), single_data = TestData(3000, 3);
    std::string error;
    Spit(tmp / "seekable.txt", seekable_data);
    ASSERT_TRUE(zstd::CompressFile(tmp / "seekable.txt", tmp / "seekable.zst", SmallFrames(1500), false, error)) << error;
    Spit(tmp / "single.txt", single_data);
    ASSERT_TRUE(zstd::CompressFile(tmp / "single.txt", tmp / "single.zst", {3, 0, 1, 0}, false, error)) << error;
    ASSERT_FALSE(zstd::IsSeekable(tmp / "single.zst"));
    Spit(tmp / "raw.txt", raw_data);
    Spit(tmp / "empty.txt", "");

    std::vector<db::Source> const sources = {{"a.zst-member", tmp / "seekable.zst"},
                                             {"raw.txt", tmp / "raw.txt"},
                                             {"empty.txt", tmp / "empty.txt"},
                                             {"single.txt", tmp / "single.zst"}};
    auto const written = db::Write(tmp / "database.protal", sources, SmallFrames(4096), error);
    ASSERT_TRUE(written) << error;
    EXPECT_EQ(*written, fs::file_size(tmp / "database.protal"));
    EXPECT_FALSE(fs::exists(tmp / "database.protal.partial"));

    // The whole file is plain seekable zstd whose content starts with the directory.
    zstd::InputFile whole(tmp / "database.protal");
    char magic[8] = {};
    whole.Stream().read(magic, 8);
    EXPECT_EQ(std::string(magic, 8), "PROTALDB");

    auto const bundle = db::Bundle::Open(tmp / "database.protal", error);
    ASSERT_TRUE(bundle) << error;
    ASSERT_EQ(bundle->Members().size(), 4u);
    EXPECT_EQ(bundle->Members()[0].name, "a.zst-member");
    EXPECT_EQ(bundle->Members()[0].frames.frames.size(), zstd::ReadSeekTable(tmp / "seekable.zst", error)->frames.size());
    EXPECT_EQ(bundle->Members()[1].frames.frames.size(), 3u);  // 9000 bytes in 4096-byte frames
    EXPECT_EQ(bundle->Members()[2].frames.frames.size(), 0u);
    EXPECT_EQ(bundle->Find("missing"), nullptr);

    std::pair<std::string, std::string> const expected[] = {
            {"a.zst-member", seekable_data}, {"raw.txt", raw_data}, {"empty.txt", ""}, {"single.txt", single_data}};
    for (auto const& [name, data] : expected) {
        auto const file = db::DbFile::InBundle(*bundle, name);
        EXPECT_TRUE(file.Exists());
        EXPECT_TRUE(file.InBundle());
        EXPECT_TRUE(file.Compressed());
        EXPECT_EQ(file.Size().value_or(0), data.size()) << name;
        EXPECT_EQ(Content(file), data) << name;
        for (int threads : {1, 3, 8}) {
            StringSink sink(data.size());
            EXPECT_EQ(file.ParallelRead(threads, sink, error), data.size()) << name;
            EXPECT_TRUE(error.empty()) << error;
            EXPECT_EQ(sink.data, data) << name << ", " << threads << " threads";
        }
        // tellg follows the content position across frames.
        auto input = file.Open();
        std::string part(2000, '\0');
        input->Stream().read(part.data(), 2000);
        if (data.size() >= 2000) EXPECT_EQ(input->Stream().tellg(), std::streampos(2000)) << name;
    }
    EXPECT_FALSE(db::DbFile::InBundle(*bundle, "missing").Exists());
    EXPECT_FALSE(db::DbFile::InBundle(*bundle, "missing").ReadAll(error));
}

TEST(Database, OtherFilesAreNotSingleFileDatabases) {
    TempDir tmp;
    std::string error;
    Spit(tmp / "raw.txt", TestData(5000));
    ASSERT_TRUE(zstd::CompressFile(tmp / "raw.txt", tmp / "seekable.zst", SmallFrames(1000), false, error)) << error;
    for (auto const& path : {tmp / "raw.txt", tmp / "seekable.zst", tmp / "missing"}) {
        EXPECT_FALSE(db::Bundle::Open(path, error)) << path;
        EXPECT_TRUE(error.empty()) << path << ": " << error;
        EXPECT_FALSE(db::IsBundle(path));
    }
}

// Truncated or corrupt files fail with a message, not with wrong content.
TEST(Database, TruncatedAndCorruptFilesFail) {
    TempDir tmp;
    std::string error;
    std::string const data = TestData(30000, 7);
    Spit(tmp / "data.txt", data);
    ASSERT_TRUE(db::Write(tmp / "db.protal", {{"data.txt", tmp / "data.txt"}}, SmallFrames(8192), error)) << error;
    std::string const bytes = Slurp(tmp / "db.protal");

    Spit(tmp / "cut.protal", bytes.substr(0, bytes.size() - 20));
    EXPECT_FALSE(db::Bundle::Open(tmp / "cut.protal", error));
    EXPECT_NE(error.find("seek table"), std::string::npos) << error;

    // A flipped byte in the middle of a member frame: its checksum no longer matches.
    std::string corrupt = bytes;
    corrupt[bytes.size() / 2] ^= 0x5a;
    Spit(tmp / "corrupt.protal", corrupt);
    error.clear();
    auto const bundle = db::Bundle::Open(tmp / "corrupt.protal", error);
    ASSERT_TRUE(bundle) << error;
    auto const file = db::DbFile::InBundle(*bundle, "data.txt");
    EXPECT_FALSE(file.ReadAll(error));
    StringSink sink(data.size());
    error.clear();
    file.ParallelRead(2, sink, error);
    EXPECT_FALSE(error.empty());
}

// A database rewritten from its own members' frames with one member replaced (as --add_model does):
// the others are byte-identical frames, the new one reads back.
TEST(Database, RewriteWithOneMemberReplaced) {
    TempDir tmp;
    std::string error;
    std::string const big = TestData(40000, 3), model = "<PMML>old</PMML>", replaced = "<PMML>new model</PMML>";
    Spit(tmp / "big.txt", big);
    Spit(tmp / "model.xml", model);
    Spit(tmp / "new.xml", replaced);
    ASSERT_TRUE(db::Write(tmp / "database.protal", {{"big.txt", tmp / "big.txt"}, {"model_pe.xml", tmp / "model.xml"}},
                          SmallFrames(4096), error)) << error;
    auto const old = db::Bundle::Open(tmp / "database.protal", error);
    ASSERT_TRUE(old) << error;
    std::string const old_bytes = Slurp(tmp / "database.protal");
    auto const& old_big = old->Members()[0].frames.frames;

    std::vector<db::Source> sources = {{"big.txt", old->Path(), old->Members()[0].frames},
                                       {"model_pe.xml", tmp / "new.xml"},
                                       {"model_se.xml", old->Path(), old->Members()[1].frames}};
    ASSERT_TRUE(db::Write(tmp / "database.protal", sources, SmallFrames(4096), error)) << error;
    auto const rewritten = db::Bundle::Open(tmp / "database.protal", error);
    ASSERT_TRUE(rewritten) << error;
    ASSERT_EQ(rewritten->Members().size(), 3u);
    EXPECT_EQ(Content(db::DbFile::InBundle(*rewritten, "big.txt")), big);
    EXPECT_EQ(Content(db::DbFile::InBundle(*rewritten, "model_pe.xml")), replaced);
    EXPECT_EQ(Content(db::DbFile::InBundle(*rewritten, "model_se.xml")), model);
    std::string const new_bytes = Slurp(tmp / "database.protal");
    auto const& new_big = rewritten->Members()[0].frames.frames;
    ASSERT_EQ(new_big.size(), old_big.size());
    for (size_t f = 0; f < new_big.size(); f++) {
        EXPECT_EQ(new_bytes.substr(new_big[f].compressed_offset, new_big[f].compressed_size),
                  old_bytes.substr(old_big[f].compressed_offset, old_big[f].compressed_size)) << "frame " << f;
    }
}

TEST(Database, ModelsPerReadType) {
    EXPECT_EQ(db::ModelCandidates("pe"), (std::vector<std::string>{"model_pe.xml", "model.xml", "random_forest.xml"}));
    EXPECT_EQ(db::ModelCandidates("se"), std::vector<std::string>{"model_se.xml"});
    EXPECT_EQ(db::ModelCandidates("pb"), std::vector<std::string>{"model_PB.xml"});
    EXPECT_EQ(db::ModelCandidates("ont"), std::vector<std::string>{"model_ONT.xml"});
    EXPECT_EQ(db::FindReadType("illumina"), nullptr);
    EXPECT_EQ(db::AllModelFiles().size(), 6u);
}

TEST(Database, WriteRejectsBadMembers) {
    TempDir tmp;
    std::string error;
    Spit(tmp / "a.txt", "a");
    EXPECT_FALSE(db::Write(tmp / "x.protal", {{"../a.txt", tmp / "a.txt"}}, SmallFrames(100), error));
    EXPECT_FALSE(db::Write(tmp / "x.protal", {{"a.txt", tmp / "a.txt"}, {"a.txt", tmp / "a.txt"}}, SmallFrames(100), error));
    EXPECT_FALSE(db::Write(tmp / "x.protal", {{"a.txt", tmp / "missing"}}, SmallFrames(100), error));
    EXPECT_FALSE(db::Write(tmp / "x.protal", {{"a.txt", tmp / "a.txt"}}, SmallFrames(0), error));
    EXPECT_FALSE(fs::exists(tmp / "x.protal"));
    EXPECT_FALSE(fs::exists(tmp / "x.protal.partial"));
}

// Separate files next to database.protal take precedence; a file given as --db is the database.
TEST(Database, LocatePrefersSeparateFiles) {
    TempDir tmp;
    fs::create_directories(tmp / "only");
    fs::create_directories(tmp / "both");
    fs::create_directories(tmp / "none");
    Spit(tmp / "only/database.protal", "x");
    Spit(tmp / "both/database.protal", "x");
    Spit(tmp / "both/index.prx.zst", "x");

    auto only = db::Locate(tmp / "only");
    EXPECT_EQ(only.bundle, tmp / "only/database.protal");
    EXPECT_EQ(only.dir, tmp / "only");
    auto both = db::Locate(tmp / "both");
    EXPECT_TRUE(both.bundle.empty());
    EXPECT_EQ(both.unused_bundle, tmp / "both/database.protal");
    auto none = db::Locate(tmp / "none");
    EXPECT_TRUE(none.bundle.empty());
    EXPECT_TRUE(none.unused_bundle.empty());
    auto file = db::Locate(tmp / "both/database.protal");
    EXPECT_EQ(file.bundle, tmp / "both/database.protal");
    EXPECT_EQ(file.dir, tmp / "both");
}

// The index's column chunks decode from the member's frames as from index.prx.zst.
TEST(Database, IndexChunksDecodeFromTheMember) {
    TempDir tmp;
    index_codec::Layout layout;
    layout.blocks = 3000;
    std::vector<uint16_t> keymap(layout.KeymapCells(), 0);
    std::vector<uint64_t> values;
    std::mt19937_64 rng(5);
    uint64_t v = 0;
    for (uint64_t b = 0; b < layout.blocks; b++) {
        uint16_t* cells = keymap.data() + b * layout.CellsPerBlock();
        std::memcpy(cells, &v, 8);
        uint64_t local = 0;
        for (uint64_t j = 0; j < layout.keys_per_block; j++) {
            cells[4 + j] = static_cast<uint16_t>(local);
            uint64_t const c = rng() % 3 == 0 ? 1 + rng() % 4 : 0;
            uint64_t const n = c >= 2 ? c + (c + 1) / 2 : c;
            for (uint64_t i = 0; i < n - c; i++) values.push_back(rng());
            for (uint64_t i = 0; i < c; i++) values.push_back((1 + rng() % 900) << 40 | (1 + rng() % 120) << 20 | rng() % 3000);
            local += n;
        }
        v += local;
    }
    std::memcpy(keymap.data() + layout.blocks * layout.CellsPerBlock(), &v, 8);
    layout.values = values.size();

    std::string error;
    ASSERT_TRUE(index_codec::Write(tmp / "index.prx.zst", "header", layout, keymap.data(), values.data(), SmallFrames(20000), error)) << error;
    Spit(tmp / "reference.map", "1\t1\t0\t10\n");
    ASSERT_TRUE(db::Write(tmp / "database.protal", {{"reference.map", tmp / "reference.map"}, {"index.prx", tmp / "index.prx.zst"}},
                          SmallFrames(1 << 20), error)) << error;
    auto const bundle = db::Bundle::Open(tmp / "database.protal", error);
    ASSERT_TRUE(bundle) << error;
    auto const index = db::DbFile::InBundle(*bundle, "index.prx");
    auto const container = index_codec::ReadContainer(index.Path(), index.Frames(), error);
    ASSERT_TRUE(container) << error;
    EXPECT_GT(container->chunks.size(), 2u);
    for (int threads : {1, 4}) {
        std::vector<uint16_t> km(layout.KeymapCells(), 1);
        std::vector<uint64_t> vals(layout.values, 1);
        EXPECT_EQ(index_codec::Decode(index.Path(), index.Frames(), *container, km.data(), vals.data(), threads), "");
        EXPECT_EQ(km, keymap);
        EXPECT_EQ(vals, values);
    }
}

// The reference, its map and the taxonomy load from a single-file database as from the files, with
// the same fingerprint.
TEST(Database, ReferenceAndTaxonomyLoadFromMembers) {
    TempDir tmp;
    std::vector<std::pair<std::string, std::string>> records;
    for (int i = 1; i <= 30; i++) records.emplace_back("1_" + std::to_string(i), TestData(300 + 41 * i, 10 + i));
    {
        std::ofstream fna(tmp / "reference.fna", std::ios::binary);
        std::ofstream map(tmp / "reference.map");
        size_t offset = 0;
        for (auto const& [name, seq] : records) {
            std::string const header = ">" + name + "\n";
            fna << header << seq << '\n';
            offset += header.size();
            map << "1\t" << name.substr(2) << '\t' << offset << '\t' << offset + seq.size() << '\n';
            offset += seq.size() + 1;
        }
    }
    Spit(tmp / "internal_taxonomy.dmp",
         "id\tparent_id\texternal_id\tname\trank\tlevel\trep_genome\n"
         "1\t3\t0\ts__Mockella alpha\tspecies\t7\tGCF_1\n"
         "2\t2\t0\troot\tno rank\t0\t\n"
         "3\t2\t0\tg__Mockella\tgenus\t6\t\n");
    std::string error;
    ASSERT_TRUE(db::Write(tmp / "database.protal",
                          {{"reference.fna", tmp / "reference.fna"}, {"reference.map", tmp / "reference.map"},
                           {"internal_taxonomy.dmp", tmp / "internal_taxonomy.dmp"}},
                          SmallFrames(1000), error)) << error;
    auto const bundle = db::Bundle::Open(tmp / "database.protal", error);
    ASSERT_TRUE(bundle) << error;
    auto const fna = db::DbFile::InBundle(*bundle, "reference.fna");
    auto const map = db::DbFile::InBundle(*bundle, "reference.map");
    EXPECT_GT(fna.Frames().frames.size(), 5u);

    for (int threads : {1, 4}) {
        GenomeLoader loader(fna, map);
        EXPECT_TRUE(loader.IsCompressed());
        loader.LoadAllGenomes(threads);
        for (int i = 1; i <= 30; i++) {
            EXPECT_EQ(loader.GetGenome(1).GetGene(i).Sequence(), records[i - 1].second) << "gene " << i << ", " << threads << " threads";
        }
    }
    EXPECT_EQ(ReferenceFingerprint::Of(map, fna), ReferenceFingerprint::Of(tmp / "reference.map", tmp / "reference.fna"));

    auto const taxonomy_file = db::DbFile::InBundle(*bundle, "internal_taxonomy.dmp");
    auto input = taxonomy_file.Open();
    taxonomy::IntTaxonomy taxonomy(input->Stream(), taxonomy_file.Name());
    EXPECT_EQ(taxonomy.Get("s__Mockella alpha"), 1u);
    EXPECT_EQ(taxonomy.root_id, 2u);
}
