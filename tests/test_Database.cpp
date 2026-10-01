// Unit tests for the single-file database (Utilities/Database.h): writing it from seekable and other
// files, reading members sequentially and in parallel, the index's column chunks and the reference
// read through member frames, where --db points (Locate), and failures on truncated or corrupt files;
// the genes' conservation factors (gene_conservation.tsv): reading, writing and estimating them.
#include <gtest/gtest.h>
#include <cctype>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>
#include "ReadType.h"
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
    EXPECT_EQ(ModelCandidates(ReadType::Paired), (std::vector<std::string>{"model_pe.xml", "model.xml", "random_forest.xml"}));
    EXPECT_EQ(ModelCandidates(ReadType::Single), std::vector<std::string>{"model_se.xml"});
    EXPECT_EQ(ModelCandidates(ReadType::PacBio), std::vector<std::string>{"model_PB.xml"});
    EXPECT_EQ(ModelCandidates(ReadType::ONT), std::vector<std::string>{"model_ONT.xml"});
    EXPECT_EQ(ReadTypeFromToken("illumina"), std::nullopt);
    EXPECT_EQ(ReadTypeTokens(), "pe, se, pb, ont");
    EXPECT_EQ(AllModelFiles().size(), 6u);
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
    for (auto const& found : {only, both, none, file}) EXPECT_TRUE(found.missing.empty()) << found.missing;
}

// A path with nothing at it is missing, not a folder of separate files; one that cannot be looked at
// says why.
TEST(Database, LocateReportsAMissingPath) {
    TempDir tmp;
    Spit(tmp / "file", "x");
    for (auto const& path : {tmp / "nothing", tmp / "nothing/database.protal", tmp / "file/database.protal", std::string()}) {
        auto const missing = db::Locate(path);
        EXPECT_EQ(missing.missing, "does not exist") << path;
        EXPECT_TRUE(missing.bundle.empty()) << path;
        EXPECT_TRUE(missing.unused_bundle.empty()) << path;
    }
    fs::create_directories(tmp / "closed/db");
    fs::permissions(tmp / "closed", fs::perms::none);
    auto const closed = db::Locate(tmp / "closed/db");
    fs::permissions(tmp / "closed", fs::perms::owner_all);
    if (::geteuid() != 0) EXPECT_EQ(closed.missing.rfind("cannot be accessed: ", 0), 0u) << closed.missing;  // root may look
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

// gene_conservation.tsv (SequenceUtils/GeneConservation.h): the factors that scale the depth identity
// margin per gene, as --build estimates and writes them and queries read them.
namespace {
    std::string RandomBases(std::mt19937& rng, size_t n) {
        static constexpr char kBases[] = "ACGT";
        std::string s(n, 'A');
        for (auto& c : s) c = kBases[rng() % 4];
        return s;
    }

    // seq with a share `rate` of its bases substituted.
    std::string Mutate(std::mt19937& rng, std::string seq, double rate) {
        static constexpr char kBases[] = "ACGT";
        std::bernoulli_distribution change(rate);
        for (auto& c : seq) {
            if (!change(rng)) continue;
            char other = c;
            while (other == c) other = kBases[rng() % 4];
            c = other;
        }
        return seq;
    }
}

TEST(GeneConservation, TheMarginIsFixedForReadErrorsAndScaledForTheRest) {
    using gene_conservation::GeneMargin;
    EXPECT_NEAR(GeneMargin(0.08, 1.0), 0.08, 1e-12);  // a gene of typical conservation: the margin
    EXPECT_NEAR(GeneMargin(0.08, 0.4), 0.05, 1e-12);  // 0.03 + 0.05 x 0.4
    EXPECT_NEAR(GeneMargin(0.08, 1.6), 0.11, 1e-12);
    EXPECT_NEAR(GeneMargin(0.02, 3.0), 0.02, 1e-12);  // below the fixed part: not scaled
    EXPECT_EQ(GeneMargin(1.0, 0.4), 1.0);             // every read counts
}

TEST(GeneConservation, TheTableReadsWhatItWrites) {
    gene_conservation::Table table;
    table.Set(3, 0.4, 25);
    table.Set(120, 1.62, 40);
    EXPECT_EQ(table.Genes(), 2u);
    std::stringstream ss;
    table.Write(ss);
    EXPECT_EQ(ss.str(), "geneid\tfactor\tspecies\n3\t0.4000\t25\n120\t1.6200\t40\n");
    gene_conservation::Table read;
    ASSERT_EQ(read.Read(ss), "");
    EXPECT_EQ(read.Genes(), 2u);
    EXPECT_NEAR(read.Factor(3), 0.4, 1e-6);
    EXPECT_NEAR(read.Factor(120), 1.62, 1e-6);
    EXPECT_EQ(read.Factor(4), 1.0);        // a gene without a factor
    EXPECT_EQ(read.Factor(100000), 1.0);
    auto const [low, high] = read.Range();
    EXPECT_NEAR(low, 0.4, 1e-6);
    EXPECT_NEAR(high, 1.62, 1e-6);
    EXPECT_EQ(gene_conservation::Table().Range(), (std::pair<double, double>{1, 1}));

    std::istringstream two_columns("# comment\n7\t0.9\n");  // the species column is optional
    gene_conservation::Table short_table;
    EXPECT_EQ(short_table.Read(two_columns), "");
    EXPECT_NEAR(short_table.Factor(7), 0.9, 1e-6);
}

TEST(GeneConservation, ReadingStopsAtTheFirstBadLine) {
    auto problem = [](std::string const& content) {
        std::istringstream is(content);
        gene_conservation::Table table;
        return table.Read(is);
    };
    EXPECT_EQ(problem("geneid\tfactor\tspecies\n1\t1.0\t3\nx\t1.0\t3\n"), "line 3: the gene id is not a number");
    EXPECT_EQ(problem("1 1.0\n"), "line 1: expected a gene id and a factor, separated by a tab");
    EXPECT_EQ(problem("1\tfast\n"), "line 1: the factor is not a number");
    EXPECT_EQ(problem("1\t0\n"), "line 1: the factor must be above 0 and at most 100");
    EXPECT_EQ(problem("1\t-0.5\n"), "line 1: the factor must be above 0 and at most 100");
    EXPECT_EQ(problem("1\tnan\n"), "line 1: the factor is not a number");
    EXPECT_EQ(problem("1\t1.0\n1\t1.2\n"), "line 2: gene 1 is listed twice");
    EXPECT_EQ(problem("1048576\t1.0\n"), "line 1: gene id 1048576 is too large");
}

TEST(GeneConservation, TheMashDistanceEstimatesTheShareOfDifferentBases) {
    std::mt19937 rng(7);
    std::string const gene = RandomBases(rng, 3000);
    auto const kmers = gene_conservation::Kmers(gene);
    EXPECT_EQ(gene_conservation::MashDistance(kmers, kmers), 0.0);
    for (double rate : { 0.005, 0.02, 0.05 }) {
        double const d = gene_conservation::MashDistance(kmers, gene_conservation::Kmers(Mutate(rng, gene, rate)));
        EXPECT_NEAR(d, rate, 0.25 * rate + 0.002) << "rate " << rate;
    }
    EXPECT_GT(gene_conservation::MashDistance(kmers, gene_conservation::Kmers(RandomBases(rng, 3000))), 0.25);
    EXPECT_EQ(gene_conservation::Kmers("ACGTNACGTACGTACGTAC").size(), 3u);  // the 14 bases after the N
    EXPECT_EQ(gene_conservation::Kmers("ACGTNACGTACG").size(), 0u);         // no 12 bases without an N
}

TEST(GeneConservation, TheEstimateFollowsHowFastEachGeneDiverges) {
    // 30 species of 21 genes: genes 1, 4, ... diverge at 0.4 times a species' rate, genes 2, 5, ... at 1,
    // genes 3, 6, ... at 1.6; each species' 3 other genomes 0.5-3% from its representative at a gene
    // of rate 1. The full reference lists the representative's own copy too.
    std::mt19937 rng(11);
    std::vector<double> const rate = { 0.4, 1.0, 1.6 };
    std::vector<uint64_t> keys;
    std::map<std::pair<uint64_t, uint64_t>, std::string> reference;
    for (uint64_t taxid = 1; taxid <= 30; taxid++) {
        for (uint64_t gene = 1; gene <= 21; gene++) {
            keys.push_back(gene_conservation::Estimator::Key(taxid, gene));
            reference[{taxid, gene}] = RandomBases(rng, 900);
        }
    }
    // A species with identical genomes and one with only 5 genes inform nothing.
    for (uint64_t gene = 1; gene <= 21; gene++) keys.push_back(gene_conservation::Estimator::Key(31, gene));
    for (uint64_t gene = 1; gene <= 5; gene++) keys.push_back(gene_conservation::Estimator::Key(32, gene));
    gene_conservation::Estimator estimator(keys);
    for (uint64_t taxid = 1; taxid <= 30; taxid++) {
        double const divergence = 0.005 + 0.025 * (taxid - 1) / 29.0;
        for (uint64_t gene = 1; gene <= 21; gene++) {
            auto const& rep = reference[{taxid, gene}];
            auto slot = estimator.Take(taxid, gene);
            ASSERT_TRUE(slot);
            estimator.Add(*slot, rep, rep);  // the representative's own copy
            for (int genome = 0; genome < 3; genome++) {
                slot = estimator.Take(taxid, gene);
                ASSERT_TRUE(slot);
                estimator.Add(*slot, rep, Mutate(rng, rep, divergence * rate[(gene - 1) % 3]));
            }
        }
    }
    std::string const same = RandomBases(rng, 900);
    for (uint64_t gene = 1; gene <= 21; gene++) {
        for (int genome = 0; genome < 3; genome++) estimator.Add(*estimator.Take(31, gene), same, same);
    }
    for (uint64_t gene = 1; gene <= 5; gene++) estimator.Add(*estimator.Take(32, gene), same, Mutate(rng, same, 0.02));
    EXPECT_FALSE(estimator.Take(33, 1)) << "not a reference gene";

    auto const estimate = estimator.Finish();
    EXPECT_EQ(estimate.species, 30u);
    EXPECT_EQ(estimate.species_with_copies, 32u);
    EXPECT_EQ(estimate.table.Genes(), 21u);
    for (uint64_t gene = 1; gene <= 21; gene++) {
        double const truth = rate[(gene - 1) % 3];
        double const shrunk = (30 * truth + gene_conservation::kPrior) / (30 + gene_conservation::kPrior);
        EXPECT_NEAR(estimate.table.Factor(gene), shrunk, 0.12) << "gene " << gene << ", rate " << truth;
    }
}

TEST(GeneConservation, CongenersDifferMostOnTheFastGenes) {
    // Genus 1: 6 species from one ancestor, each 4% from it at a gene of rate 1 (8% between two species): genes 1-6 at
    // rate 0.3, 7-12 at 1.7, and gene 13 the same in every species. Genus 2 has one species: nothing to compare.
    std::mt19937 rng(11);
    std::vector<double> const rates = { 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 0 };
    std::vector<std::string> ancestor;
    for (size_t g = 0; g < rates.size(); g++) ancestor.push_back(RandomBases(rng, 1200));
    std::map<uint32_t, std::vector<std::pair<uint64_t, std::string>>> species;
    for (uint32_t taxid = 1; taxid <= 6; taxid++) {
        for (size_t g = 0; g < rates.size(); g++) species[taxid].emplace_back(g + 1, Mutate(rng, ancestor[g], 0.04 * rates[g]));
    }
    species[7] = species[1];
    gene_conservation::Table within;
    for (uint64_t g = 1; g <= 12; g++) within.Set(g, g <= 6 ? 0.5 : 1.5, 6);
    within.Set(13, 0.5, 6);
    size_t calls = 0;
    auto const estimate = gene_conservation::CompareCongeners({ { 1, 2, 3, 4, 5, 6 }, { 7 } }, [&](uint32_t taxid) {
        calls++;
        return species.at(taxid);
    }, within, 1);
    EXPECT_EQ(calls, 6u);  // the species of genus 1, once each
    EXPECT_EQ(estimate.genera, 1u);
    EXPECT_EQ(estimate.species, 6u);
    EXPECT_EQ(estimate.pairs, 15u);  // each of 6 species against the next 4, cyclically: all 15 pairs
    ASSERT_EQ(estimate.genes.size(), 13u);
    // A pair's median gene is a slow one (6 slow, 6 fast, 1 the same): slow genes about 1, fast ones about 5.7.
    for (auto const& g : estimate.genes) {
        if (g.geneid <= 6) EXPECT_NEAR(g.between, 1.0, 0.4) << g.geneid;
        else if (g.geneid <= 12) EXPECT_GT(g.between, 3.0) << g.geneid;
        else EXPECT_EQ(g.between, 0.0);
        EXPECT_EQ(g.pairs, 15u) << g.geneid;
        EXPECT_EQ(g.species, 6u) << g.geneid;
        EXPECT_EQ(g.identical, g.geneid == 13 ? 6u : 0u) << g.geneid;
    }
    EXPECT_GT(estimate.spearman, 0.8);
    EXPECT_EQ(estimate.correlated, 13u);
    EXPECT_NEAR(estimate.conserved_between, 1.0, 0.4);
    EXPECT_GT(estimate.fast_between, 3.0);
    EXPECT_NEAR(estimate.conserved_identical, 1.0 / 7, 1e-12);  // gene 13 of the 7 conserved genes
    EXPECT_EQ(estimate.fast_identical, 0.0);
    std::ostringstream os;
    estimate.Write(os);
    EXPECT_EQ(os.str().substr(0, os.str().find('\n')),
              "geneid\twithin_factor\tbetween_factor\tpairs\tspecies\tidentical_share\tnear_identical_share");

    EXPECT_NEAR(gene_conservation::Spearman({ 1, 2, 3, 4 }, { 10, 20, 30, 40 }), 1.0, 1e-12);
    EXPECT_NEAR(gene_conservation::Spearman({ 1, 2, 3, 4 }, { 4, 3, 2, 1 }), -1.0, 1e-12);
    EXPECT_TRUE(std::isnan(gene_conservation::Spearman({ 1, 2 }, { 1, 2 })));
}

TEST(GeneConservation, CopiesBeyondTheCapAreNotCompared) {
    gene_conservation::Estimator estimator({ gene_conservation::Estimator::Key(1, 1) });
    for (uint32_t i = 0; i < gene_conservation::kMaxCopies; i++) EXPECT_TRUE(estimator.Take(1, 1));
    EXPECT_FALSE(estimator.Take(1, 1));
    EXPECT_TRUE(estimator.Finish().table.Empty());  // one species, one gene: no factor
}
