// Unit tests for the single-file database (Utilities/Database.h): writing it from seekable and other
// files, reading members sequentially and in parallel, the index's column chunks and the reference
// read through member frames, where --db points (Locate), and failures on truncated or corrupt files.
#include <gtest/gtest.h>
#include <algorithm>
#include <cctype>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <map>
#include <optional>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <unistd.h>
#include "ReadType.h"
#include "Utilities/Database.h"
#include "Utilities/ReferenceFingerprint.h"
#include "Hash/IndexCodec.h"
#include "SequenceUtils/GenomeLoader.h"
#include "Taxonomy/Taxonomy.h"
#include "TestUtil.h"

namespace fs = std::filesystem;
using namespace protal;
using namespace protal::test;

namespace {

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
    ScratchDir tmp;
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
    ScratchDir tmp;
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
    ScratchDir tmp;
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

namespace {
    using Listed = std::tuple<std::string, uint64_t, uint64_t>;  // a member's name, first frame and number of frames

    // A directory: the magic, the version, the member count (that of `members` unless given), then the members.
    std::string Directory(std::vector<Listed> const& members, uint64_t version = db::kVersion,
                          std::optional<uint64_t> count = std::nullopt) {
        std::string directory(db::kMagic, db::kMagic + 8);
        db::detail::PutU64(directory, version);
        db::detail::PutU64(directory, count.value_or(members.size()));
        for (auto const& [name, first, frames] : members) {
            db::detail::PutU64(directory, name.size());
            directory += name;
            db::detail::PutU64(directory, first);
            db::detail::PutU64(directory, frames);
        }
        return directory;
    }

    // A file in the single-file format: `directory` in a raw frame (as db::Write writes it), `frames` frames of a few
    // bytes ("frame 0", "frame 1", ...) and the seek table, which gives the directory's content size as
    // `directory_size` if set.
    std::string WriteWithDirectory(ScratchDir const& tmp, std::string const& name, std::string const& directory, size_t frames,
                                   std::optional<uint64_t> directory_size = std::nullopt) {
        std::string const path = tmp / name;
        zstd::FrameWriter out(path);
        std::string const head = db::detail::RawFrame(directory);
        EXPECT_TRUE(out.Add(head.data(), head.size(), directory_size.value_or(directory.size())));
        for (size_t i = 0; i < frames; i++) {
            std::string const content = "frame " + std::to_string(i);
            std::string const frame = db::detail::RawFrame(content);
            EXPECT_TRUE(out.Add(frame.data(), frame.size(), content.size()));
        }
        EXPECT_TRUE(out.Finish()) << out.Error();
        return path;
    }
}

// Bundle::Open checks a directory before it uses it: each kind of damage fails with its own message.
TEST(Database, ADamagedDirectoryIsRejected) {
    ScratchDir tmp;
    std::string error;
    std::vector<Listed> const members = {{"a.txt", 1, 1}, {"b.txt", 2, 2}};
    auto const good = db::Bundle::Open(WriteWithDirectory(tmp, "good.protal", Directory(members), 3), error);
    ASSERT_TRUE(good) << error;
    EXPECT_EQ(Content(db::DbFile::InBundle(*good, "b.txt")), "frame 1frame 2");

    std::vector<std::pair<std::string, std::string>> const cases = {  // directory, message
        {Directory(members, db::kVersion + 1), "version " + std::to_string(db::kVersion + 1) + " of the single-file format"},
        {std::string(db::kMagic, db::kMagic + 8) + "abc", "invalid directory: too short"},
        {Directory(members, db::kVersion, 5), "invalid directory: member count"},  // more members than frames
        {Directory({{"a.txt", 1, 1}, {"b.txt", 3, 1}}), "the members do not tile the file"},
        {Directory({{"a.txt", 1, 1}, {"b.txt", 2, 1}}), "the members do not cover the file"},
        {Directory({{"a.txt", 1, 1}, {"a.txt", 2, 2}}), "member a.txt is listed twice"},
        {Directory({{"../a.txt", 1, 1}, {"b.txt", 2, 2}}), "member name '../a.txt' is not a file name"},
        {Directory({{db::kFileName, 1, 1}, {"b.txt", 2, 2}}), "member name '" + db::kFileName + "' is not a file name"},
        {Directory({{"", 1, 1}, {"b.txt", 2, 2}}), "invalid directory: member name (corrupt file?)"},
        {Directory(members) + "x", "unexpected data after the member list"},
    };
    for (size_t i = 0; i < cases.size(); i++) {
        auto const& [directory, message] = cases[i];
        error.clear();
        EXPECT_FALSE(db::Bundle::Open(WriteWithDirectory(tmp, "bad" + std::to_string(i) + ".protal", directory, 3), error)) << message;
        EXPECT_NE(error.find(message), std::string::npos) << "expected: " << message << "; got: " << error;
    }
    // A first frame of more than 16 MB (as its seek table entry says) that starts as a directory.
    error.clear();
    EXPECT_FALSE(db::Bundle::Open(WriteWithDirectory(tmp, "large.protal", Directory(members), 3, (uint64_t{16} << 20) + 1), error));
    EXPECT_NE(error.find("invalid directory: larger than 16 MB"), std::string::npos) << error;
}

// A database rewritten from its own members' frames with one member replaced (as --add_model does):
// the others are byte-identical frames, the new one reads back.
TEST(Database, RewriteWithOneMemberReplaced) {
    ScratchDir tmp;
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

// The directory frame holds its content in raw blocks: any zstd decoder reads it, and its size depends on the
// content's length alone.
TEST(Database, TheDirectoryFrameIsRaw) {
    for (size_t n : {1, 1000, 131072, 131073, 300000}) {
        std::string const content = TestData(n, static_cast<unsigned>(n));
        std::string const frame = db::detail::RawFrame(content);
        EXPECT_EQ(frame.size(), db::detail::RawFrameSize(n)) << n;
        EXPECT_EQ(ZSTD_getFrameContentSize(frame.data(), frame.size()), n);
        std::string out(n, '\0');
        size_t const got = ZSTD_decompress(out.data(), out.size(), frame.data(), frame.size());
        ASSERT_FALSE(ZSTD_isError(got)) << n << ": " << ZSTD_getErrorName(got);
        EXPECT_EQ(got, n);
        EXPECT_EQ(out, content) << n;
    }
    ScratchDir tmp;
    std::string error;
    Spit(tmp / "a.txt", TestData(5000));
    ASSERT_TRUE(db::Write(tmp / "database.protal", {{"a.txt", tmp / "a.txt"}}, SmallFrames(4096), error)) << error;
    auto const bundle = db::Bundle::Open(tmp / "database.protal", error);
    ASSERT_TRUE(bundle) << error;
    EXPECT_EQ(bundle->DirectoryFrame().compressed_size, db::detail::RawFrameSize(bundle->DirectoryFrame().decompressed_size));
}

namespace {
    // A database of `big` (copied from a seekable file), other.txt, and two models, the models last.
    struct ModelDatabase {
        ScratchDir tmp;
        std::string big = TestData(40000, 3), other = TestData(7000, 4), pe = "<PMML>pe</PMML>", se = "<PMML>se</PMML>";
        std::string path = tmp / "database.protal";

        ModelDatabase() {
            std::string error;
            Spit(tmp / "big.txt", big);
            EXPECT_TRUE(zstd::CompressFile(tmp / "big.txt", tmp / "big.zst", SmallFrames(1500), false, error)) << error;
            EXPECT_TRUE(db::Write(path, {{"big.txt", tmp / "big.zst"}, {"other.txt", Spit(tmp / "other.txt", other)},
                                         {"model_pe.xml", Spit(tmp / "pe.xml", pe)}, {"model_se.xml", Spit(tmp / "se.xml", se)}},
                                  SmallFrames(4096), error)) << error;
        }

        db::Bundle Open() const {
            std::string error;
            auto bundle = db::Bundle::Open(path, error);
            EXPECT_TRUE(bundle) << error;
            return *bundle;
        }

        // The bundle's members, those named in `replace` from these files instead.
        static std::vector<db::Source> Sources(db::Bundle const& bundle, std::map<std::string, std::string> const& replace) {
            std::vector<db::Source> sources;
            for (auto const& member : bundle.Members()) {
                auto const it = replace.find(member.name);
                if (it != replace.end()) sources.push_back({member.name, it->second});
                else sources.push_back({member.name, bundle.Path(), member.frames});
            }
            return sources;
        }
    };
}

// Replacing the models (the last members) writes only the end of the file and the directory: the other members'
// bytes stay where they were, the models read back, the file grows or shrinks with them, and no journal is left.
TEST(Database, ModelsAreReplacedInPlace) {
    ModelDatabase d;
    std::string error;
    std::string const before = Slurp(d.path);
    auto const old = d.Open();
    uint64_t const directory = old.DirectoryFrame().compressed_size;
    uint64_t const tail = old.Members()[2].frames.frames.front().compressed_offset;

    std::string const big_pe = TestData(9000, 8);  // three frames instead of one
    auto sources = ModelDatabase::Sources(old, {{"model_pe.xml", Spit(d.tmp / "new_pe.xml", big_pe)}});
    auto const first = db::InPlaceFrom(old, sources);
    ASSERT_EQ(first, std::optional<size_t>(2));
    auto const written = db::ReplaceTail(old, sources, *first, SmallFrames(4096, 3), error);
    ASSERT_TRUE(written) << error;
    EXPECT_EQ(*written, fs::file_size(d.path));
    EXPECT_FALSE(fs::exists(d.path + db::kJournalExtension));
    std::string const after = Slurp(d.path);
    EXPECT_EQ(after.substr(directory, tail - directory), before.substr(directory, tail - directory));
    EXPECT_NE(after.substr(0, directory), before.substr(0, directory));  // the frame counts changed

    auto const now = d.Open();
    ASSERT_EQ(now.Members().size(), 4u);
    EXPECT_TRUE(db::SameFrames(now.Members()[0].frames, old.Members()[0].frames));
    EXPECT_TRUE(db::SameFrames(now.Members()[1].frames, old.Members()[1].frames));
    EXPECT_EQ(now.Members()[2].frames.frames.size(), 3u);
    EXPECT_EQ(Content(db::DbFile::InBundle(now, "big.txt")), d.big);
    EXPECT_EQ(Content(db::DbFile::InBundle(now, "other.txt")), d.other);
    EXPECT_EQ(Content(db::DbFile::InBundle(now, "model_pe.xml")), big_pe);
    EXPECT_EQ(Content(db::DbFile::InBundle(now, "model_se.xml")), d.se);  // copied from the old end of the file

    // Back to one frame: the file shrinks; then the last member alone.
    sources = ModelDatabase::Sources(now, {{"model_pe.xml", d.tmp / "pe.xml"}});
    ASSERT_TRUE(db::ReplaceTail(now, sources, *db::InPlaceFrom(now, sources), SmallFrames(4096), error)) << error;
    EXPECT_EQ(Slurp(d.path), before);
    auto const again = d.Open();
    sources = ModelDatabase::Sources(again, {{"model_se.xml", Spit(d.tmp / "new_se.xml", "<PMML>se 2</PMML>")}});
    ASSERT_EQ(db::InPlaceFrom(again, sources), std::optional<size_t>(3));
    ASSERT_TRUE(db::ReplaceTail(again, sources, 3, SmallFrames(4096), error)) << error;
    EXPECT_EQ(Content(db::DbFile::InBundle(d.Open(), "model_se.xml")), "<PMML>se 2</PMML>");
    EXPECT_EQ(Content(db::DbFile::InBundle(d.Open(), "model_pe.xml")), d.pe);
}

// In place only with the same members in the same order, something to replace, a tail that fits in memory and a
// raw directory frame: a file written before (its directory compressed) is rewritten once by Write, which makes
// its directory raw.
TEST(Database, InPlaceNeedsTheSameMembersAndARawDirectory) {
    ModelDatabase d;
    std::string error;
    auto const bundle = d.Open();
    std::string const pe = Spit(d.tmp / "new_pe.xml", "<PMML>new</PMML>");
    auto sources = ModelDatabase::Sources(bundle, {});
    EXPECT_FALSE(db::InPlaceFrom(bundle, sources)) << "nothing changes";
    sources = ModelDatabase::Sources(bundle, {{"model_pe.xml", pe}});
    EXPECT_TRUE(db::InPlaceFrom(bundle, sources));
    EXPECT_FALSE(db::InPlaceFrom(bundle, sources, 10)) << "more than max_bytes";
    auto added = sources;
    added.push_back({"model_ONT.xml", pe});
    EXPECT_FALSE(db::InPlaceFrom(bundle, added));
    auto reordered = sources;
    std::swap(reordered[2], reordered[3]);
    EXPECT_FALSE(db::InPlaceFrom(bundle, reordered));
    auto renamed = sources;
    renamed[2].name = "model_PB.xml";
    EXPECT_FALSE(db::InPlaceFrom(bundle, renamed));

    // The same file with its directory compressed, as protal wrote it before.
    std::string const bytes = Slurp(d.path);
    auto const table = zstd::ReadSeekTable(d.path, error);
    ASSERT_TRUE(table) << error;
    auto const& frames = table->frames;
    std::string const listing = db::detail::Directory(*db::detail::Plan(ModelDatabase::Sources(bundle, {}), SmallFrames(4096), error));
    std::string compressed(ZSTD_compressBound(listing.size()), '\0');
    compressed.resize(ZSTD_compress(compressed.data(), compressed.size(), listing.data(), listing.size(), 3));
    std::string entries;
    zstd::PutLE32(entries, static_cast<uint32_t>(compressed.size()));
    zstd::PutLE32(entries, static_cast<uint32_t>(listing.size()));
    for (size_t f = 1; f < frames.size(); f++) {
        zstd::PutLE32(entries, static_cast<uint32_t>(frames[f].compressed_size));
        zstd::PutLE32(entries, static_cast<uint32_t>(frames[f].decompressed_size));
    }
    auto const& last = frames.back();
    Spit(d.tmp / "old.protal", compressed + bytes.substr(frames[0].compressed_size, last.compressed_offset + last.compressed_size -
                                                         frames[0].compressed_size) + zstd::SeekTableFrame(entries));
    auto const old = db::Bundle::Open(d.tmp / "old.protal", error);
    ASSERT_TRUE(old) << error;
    EXPECT_EQ(Content(db::DbFile::InBundle(*old, "big.txt")), d.big);
    sources = ModelDatabase::Sources(*old, {{"model_pe.xml", pe}});
    EXPECT_FALSE(db::InPlaceFrom(*old, sources));
    ASSERT_TRUE(db::Write(old->Path(), sources, SmallFrames(4096), error)) << error;
    auto const rewritten = db::Bundle::Open(d.tmp / "old.protal", error);
    ASSERT_TRUE(rewritten) << error;
    EXPECT_TRUE(db::InPlaceFrom(*rewritten, ModelDatabase::Sources(*rewritten, {{"model_se.xml", pe}})));
}

// A replacement that stopped half way (its journal next to the file, which then does not open) is undone from the
// journal, if the journal is of that file.
TEST(Database, AnInterruptedReplacementIsWrittenBack) {
    ModelDatabase d;
    std::string error;
    std::string const before = Slurp(d.path);
    auto const old = d.Open();
    db::detail::Journal journal;
    uint64_t const directory = old.DirectoryFrame().compressed_size;
    journal.tail_offset = old.Members()[2].frames.frames.front().compressed_offset;
    journal.check_offset = std::max<uint64_t>(directory, journal.tail_offset - std::min<uint64_t>(journal.tail_offset, 4096));
    journal.directory = before.substr(0, directory);
    journal.check = before.substr(journal.check_offset, journal.tail_offset - journal.check_offset);
    journal.tail = before.substr(journal.tail_offset);

    auto const sources = ModelDatabase::Sources(old, {{"model_pe.xml", Spit(d.tmp / "new_pe.xml", TestData(9000, 8))}});
    ASSERT_TRUE(db::ReplaceTail(old, sources, 2, SmallFrames(4096), error)) << error;
    std::string const after = Slurp(d.path);
    std::string const journal_path = d.path + db::kJournalExtension;
    auto interrupt = [&]() {
        Spit(d.path, after.substr(0, after.size() - 20));  // the seek table not yet complete
        EXPECT_FALSE(db::Bundle::Open(d.path, error));
    };

    interrupt();
    Spit(journal_path, journal.Serialize().substr(0, 100));
    EXPECT_FALSE(db::RestoreFromJournal(d.path, error));
    EXPECT_NE(error.find("not a complete journal"), std::string::npos) << error;
    auto other = journal;
    other.check[0] ^= 1;
    Spit(journal_path, other.Serialize());
    EXPECT_FALSE(db::RestoreFromJournal(d.path, error));
    EXPECT_NE(error.find("is not of"), std::string::npos) << error;
    EXPECT_EQ(Slurp(d.path), after.substr(0, after.size() - 20)) << "not touched";

    Spit(journal_path, journal.Serialize());
    EXPECT_TRUE(db::RestoreFromJournal(d.path, error)) << error;
    EXPECT_EQ(Slurp(d.path), before);
    EXPECT_FALSE(fs::exists(journal_path));
    EXPECT_EQ(Content(db::DbFile::InBundle(d.Open(), "model_pe.xml")), d.pe);

    // A journal left after the database was written anew is removed with it.
    Spit(journal_path, journal.Serialize());
    ASSERT_TRUE(db::Write(d.path, ModelDatabase::Sources(d.Open(), {{"model_se.xml", d.tmp / "pe.xml"}}), SmallFrames(4096), error)) << error;
    EXPECT_FALSE(fs::exists(journal_path));
}

// Each read type's model files, in order of precedence (paired-end reads: also those of databases from before read
// types), all among the files a database may hold; each token names its read type.
TEST(Database, ModelsPerReadType) {
    EXPECT_EQ(ModelCandidates(ReadType::Paired), (std::vector<std::string>{"model_pe.xml", "model.xml", "random_forest.xml"}));
    EXPECT_EQ(ModelCandidates(ReadType::Single), std::vector<std::string>{"model_se.xml"});
    EXPECT_EQ(ModelCandidates(ReadType::PacBio), std::vector<std::string>{"model_PB.xml"});
    EXPECT_EQ(ModelCandidates(ReadType::ONT), std::vector<std::string>{"model_ONT.xml"});
    auto const all = AllModelFiles();
    EXPECT_EQ(std::set<std::string>(all.begin(), all.end()).size(), all.size()) << "each file once";
    for (auto const type : {ReadType::Paired, ReadType::Single, ReadType::PacBio, ReadType::ONT}) {
        for (auto const& name : ModelCandidates(type)) EXPECT_NE(std::find(all.begin(), all.end(), name), all.end()) << name;
    }
    EXPECT_EQ(ReadTypeFromToken("pe"), std::optional<ReadType>(ReadType::Paired));
    EXPECT_EQ(ReadTypeFromToken("se"), std::optional<ReadType>(ReadType::Single));
    EXPECT_EQ(ReadTypeFromToken("pb"), std::optional<ReadType>(ReadType::PacBio));
    EXPECT_EQ(ReadTypeFromToken("ont"), std::optional<ReadType>(ReadType::ONT));
    EXPECT_EQ(ReadTypeFromToken("illumina"), std::nullopt);
    EXPECT_EQ(ReadTypeFromToken(""), std::nullopt);
    EXPECT_EQ(ReadTypeTokens(), "pe, se, pb, ont");
}

TEST(Database, WriteRejectsBadMembers) {
    ScratchDir tmp;
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
    ScratchDir tmp;
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
    ScratchDir tmp;
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
    ScratchDir tmp;
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
    ScratchDir tmp;
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
