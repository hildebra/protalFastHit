// Unit tests for SamFile.h: SAM files written by several threads, plain, as BGZF (.sam.gz) or as
// seekable zstd (.sam.zst), with the header after the records or first; reading them back, and
// telling a complete file from one cut at a block or frame boundary.
#include <gtest/gtest.h>
#include <zlib-ng.h>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#include <unistd.h>
#include "IO/SamFile.h"

using namespace protal;
namespace fs = std::filesystem;

namespace {
    struct ScratchDir {
        fs::path path;
        ScratchDir() {
            path = fs::temp_directory_path() / ("protal samfile test " + std::to_string(::getpid()));
            fs::create_directories(path);
        }
        ~ScratchDir() { fs::remove_all(path); }
        std::string File(std::string const& name) const { return (path / name).string(); }
    };

    std::string Slurp(std::string const& path) {
        std::ifstream is(path, std::ios::binary);
        return std::string(std::istreambuf_iterator<char>(is), {});
    }

    // Everything a stream holds, read through std::istream (which turns a read error into badbit,
    // as for the SAM reader's getline).
    std::string ReadAll(std::istream& is) {
        std::string text;
        char buffer[1 << 16];
        while (is.read(buffer, sizeof(buffer)) || is.gcount() > 0) text.append(buffer, static_cast<size_t>(is.gcount()));
        return text;
    }

    // All of a SAM file as SamInput reads it.
    std::string ReadBack(std::string const& path) {
        SamInput input(path);
        EXPECT_TRUE(input.IsOpen()) << path;
        EXPECT_EQ(input.Problem(), "") << path;
        std::string text = ReadAll(input.Stream());
        EXPECT_FALSE(input.ReadFailed()) << input.ReadError();
        return text;
    }

    std::string GzipRead(std::string const& path) {
        gzFile f = zng_gzopen(path.c_str(), "rb");
        std::string text;
        char buffer[1 << 16];
        int32_t n;
        while ((n = zng_gzread(f, buffer, sizeof(buffer))) > 0) text.append(buffer, static_cast<size_t>(n));
        zng_gzclose(f);
        return text;
    }

    // Block `b` of thread `t`: a few "records" naming genes (t, b % 7) and (t, 100 + b % 3).
    std::string Block(int t, int b, std::vector<uint64_t>& genes) {
        std::string block;
        for (int r = 0; r < 3; r++) {
            block += "t" + std::to_string(t) + "b" + std::to_string(b) + "r" + std::to_string(r) + "\t0\t" +
                     std::to_string(t) + "_" + std::to_string(b % 7) + "\t1\t60\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIII\n";
        }
        genes.push_back(SamGeneKey(t, b % 7));
        genes.push_back(SamGeneKey(t, 100 + b % 3));
        genes.push_back(SamGeneKey(t, b % 7));  // duplicates are fine
        return block;
    }

    // Writes kThreads x kBlocks blocks from kThreads threads; returns the blocks and the genes named.
    constexpr int kThreads = 4, kBlocks = 400;
    std::pair<std::vector<std::string>, std::set<uint64_t>> WriteFromThreads(SamSink& sink) {
        std::vector<std::vector<std::string>> blocks(kThreads);
        std::set<uint64_t> named;
        std::vector<std::vector<uint64_t>> all_genes(kThreads);
        std::vector<std::thread> threads;
        for (int t = 0; t < kThreads; t++) {
            threads.emplace_back([&, t] {
                std::vector<uint64_t> genes;
                for (int b = 0; b < kBlocks; b++) {
                    std::string const block = Block(t, b, genes);
                    all_genes[t].insert(all_genes[t].end(), genes.begin(), genes.end());
                    blocks[t].push_back(block);
                    sink.Write(block.data(), block.size(), genes);
                    EXPECT_TRUE(genes.empty());
                }
            });
        }
        for (auto& th : threads) th.join();
        std::vector<std::string> flat;
        for (int t = 0; t < kThreads; t++) {
            flat.insert(flat.end(), blocks[t].begin(), blocks[t].end());
            named.insert(all_genes[t].begin(), all_genes[t].end());
        }
        return { flat, named };
    }

    // The text after the header holds every block once, each whole.
    void ExpectBlocks(std::string const& records, std::vector<std::string> blocks) {
        size_t total = 0;
        for (auto const& b : blocks) {
            EXPECT_NE(records.find(b), std::string::npos) << "a block is missing or split";
            total += b.size();
        }
        EXPECT_EQ(records.size(), total);
    }

    std::string const kHeader = "@HD\tVN:1.6\n@SQ\tSN:1_1\tLN:50\n@CO\tprotal read type: pe\n";
}

TEST(SamFile, CompressionFollowsTheName) {
    EXPECT_EQ(SamCompressionOf("a.sam"), SamCompression::None);
    EXPECT_EQ(SamCompressionOf("dir.gz/a.sam"), SamCompression::None);
    EXPECT_EQ(SamCompressionOf("a.sam.gz"), SamCompression::Gzip);
    EXPECT_EQ(SamCompressionOf("a.sam.zst"), SamCompression::Zstd);
    EXPECT_EQ(UncompressedSamName("a.sam.gz"), "a.sam");
    EXPECT_EQ(UncompressedSamName("a.sam.zst"), "a.sam");
    EXPECT_EQ(UncompressedSamName("a.sam"), "a.sam");
    auto const [taxid, gene] = SamGeneOfKey(SamGeneKey(123456, 789));
    EXPECT_EQ(taxid, 123456u);
    EXPECT_EQ(gene, 789u);
}

TEST(SamFile, HeaderAfterTheRecordsInEveryFormat) {
    for (auto const& name : { "out.sam", "out.sam.gz", "out.sam.zst" }) {
        SCOPED_TRACE(name);
        ScratchDir dir;
        auto const path = dir.File(name);
        std::vector<std::string> blocks;
        std::set<uint64_t> named;
        {
            SamOutput out(path, SamCompressionOf(path));
            ASSERT_TRUE(out.Ok()) << out.Error();
            std::tie(blocks, named) = WriteFromThreads(out);
            EXPECT_TRUE(fs::exists(path + SamOutput::kRecordsSuffix));
            EXPECT_FALSE(fs::exists(path));
            auto const genes = out.Genes();
            EXPECT_EQ(genes, std::vector<uint64_t>(named.begin(), named.end()));  // sorted, each once
            ASSERT_TRUE(out.Finish(kHeader)) << out.Error();
        }
        EXPECT_FALSE(fs::exists(path + SamOutput::kRecordsSuffix));
        auto const text = ReadBack(path);
        ASSERT_EQ(text.substr(0, kHeader.size()), kHeader);
        ExpectBlocks(text.substr(kHeader.size()), blocks);
    }
}

TEST(SamFile, AHeaderGivenUpFrontIsWrittenFirst) {
    for (auto const& name : { "out.sam", "out.sam.gz", "out.sam.zst" }) {
        SCOPED_TRACE(name);
        ScratchDir dir;
        auto const path = dir.File(name);
        SamOutput out(path, SamCompressionOf(path), kHeader);
        auto [blocks, named] = WriteFromThreads(out);
        EXPECT_FALSE(fs::exists(path + SamOutput::kRecordsSuffix));  // the records go straight into the SAM
        ASSERT_TRUE(out.Finish()) << out.Error();
        auto const text = ReadBack(path);
        ASSERT_EQ(text.substr(0, kHeader.size()), kHeader);
        ExpectBlocks(text.substr(kHeader.size()), blocks);
    }
}

TEST(SamFile, GzipIsBgzfThatZlibNgReads) {
    ScratchDir dir;
    auto const path = dir.File("out.sam.gz");
    std::vector<std::string> blocks;
    {
        SamOutput out(path, SamCompression::Gzip);
        blocks = WriteFromThreads(out).first;
        ASSERT_TRUE(out.Finish(kHeader));
    }
    // Every member is a BGZF block of at most 64 KB whose size field is right; the last is the EOF block.
    std::string const bytes = Slurp(path);
    size_t offset = 0, members = 0;
    while (offset < bytes.size()) {
        ASSERT_GE(bytes.size() - offset, 28u);
        auto const* b = reinterpret_cast<unsigned char const*>(bytes.data() + offset);
        ASSERT_EQ(std::memcmp(b, bgzf::kHeader, 16), 0) << "member " << members;
        size_t const size = (b[16] | b[17] << 8) + 1;
        ASSERT_LE(size, bgzf::kMaxBlock);
        offset += size;
        members++;
    }
    EXPECT_EQ(offset, bytes.size());
    EXPECT_GT(members, 2u);
    EXPECT_EQ(bytes.compare(bytes.size() - 28, 28, reinterpret_cast<char const*>(bgzf::kEof), 28), 0);
    EXPECT_TRUE(bgzf::StartsAsBgzf(path));
    EXPECT_TRUE(bgzf::EndsWithEof(path));
    // zlib-ng's gzread (a reader other than protal's, which is ISA-L's) reads all members, as zcat does.
    auto const text = GzipRead(path);
    ASSERT_EQ(text.substr(0, kHeader.size()), kHeader);
    ExpectBlocks(text.substr(kHeader.size()), blocks);
}

TEST(SamFile, BgzfStoresBlocksThatDoNotShrink) {
    ScratchDir dir;
    auto const path = dir.File("random.sam.gz");
    std::mt19937_64 rng(7);
    std::string data(300000, '\0');
    for (auto& c : data) c = static_cast<char>(rng());
    {
        SamOutput out(path, SamCompression::Gzip);
        std::vector<uint64_t> genes;
        out.Write(data.data(), data.size(), genes);
        ASSERT_TRUE(out.Finish(""));
    }
    EXPECT_EQ(GzipRead(path), data);
    SamInput input(path);
    EXPECT_EQ(input.Problem(), "");
}

TEST(SamFile, ZstdIsSeekableWithAMarker) {
    ScratchDir dir;
    auto const path = dir.File("out.sam.zst");
    std::vector<std::string> blocks;
    {
        SamOutput out(path, SamCompression::Zstd);
        blocks = WriteFromThreads(out).first;
        ASSERT_TRUE(out.Finish(kHeader));
    }
    EXPECT_TRUE(sam_zstd::StartsWithMarker(path));
    std::string error;
    auto const table = zstd::ReadSeekTable(path, error);
    ASSERT_TRUE(table.has_value()) << error;
    size_t total = kHeader.size();
    for (auto const& b : blocks) total += b.size();
    EXPECT_EQ(table->DecompressedSize(), total);
    EXPECT_GT(table->frames.size(), 2u);  // the marker, the header, the records' frames
    EXPECT_EQ(table->frames.front().decompressed_size, 0u);
}

TEST(SamFile, AFileCutAtABlockOrFrameBoundaryIsIncomplete) {
    ScratchDir dir;
    // BGZF without its end-of-file block, which is where a cut at a block boundary leaves it.
    auto const gz = dir.File("out.sam.gz");
    {
        SamOutput out(gz, SamCompression::Gzip);
        WriteFromThreads(out);
        ASSERT_TRUE(out.Finish(kHeader));
    }
    fs::resize_file(gz, fs::file_size(gz) - 28);
    EXPECT_NE(SamInput(gz).Problem().find("end-of-file block"), std::string::npos);

    // zstd cut after its last frame: the seek table is gone.
    auto const zst = dir.File("out.sam.zst");
    {
        SamOutput out(zst, SamCompression::Zstd);
        WriteFromThreads(out);
        ASSERT_TRUE(out.Finish(kHeader));
    }
    std::string error;
    auto const table = zstd::ReadSeekTable(zst, error);
    ASSERT_TRUE(table.has_value());
    auto const& last = table->frames.back();
    fs::resize_file(zst, last.compressed_offset + last.compressed_size);
    EXPECT_NE(SamInput(zst).Problem(), "");

    // Cut inside a block or frame, the gzip reader and zstd notice while reading. The zstd file is cut in the middle of
    // its middle frame: the threads write many small frames, and half its size fell on a frame's end now and then (a
    // whole frame read, then the end of the file: Problem() above tells that the seek table is missing).
    auto const& middle = table->frames[table->frames.size() / 2];
    ASSERT_GE(middle.compressed_size, 2u);
    for (auto const& path : { gz, zst }) {
        fs::resize_file(path, path == zst ? middle.compressed_offset + middle.compressed_size / 2 : fs::file_size(path) / 2);
        SamInput input(path);
        ReadAll(input.Stream());
        EXPECT_TRUE(input.ReadFailed()) << path;
        EXPECT_NE(input.ReadError(), "");
    }

    // A gzip file from another tool (one member, no BGZF) is complete without the EOF block.
    auto const plain_gz = dir.File("other.sam.gz");
    {
        gzFile f = zng_gzopen(plain_gz.c_str(), "wb");
        zng_gzwrite(f, kHeader.data(), static_cast<uint32_t>(kHeader.size()));
        zng_gzclose(f);
    }
    EXPECT_EQ(SamInput(plain_gz).Problem(), "");
    EXPECT_EQ(ReadBack(plain_gz), kHeader);
}

TEST(SamFile, DiscardAndFailuresLeaveNoFiles) {
    ScratchDir dir;
    auto const path = dir.File("out.sam.zst");
    {
        SamOutput out(path, SamCompression::Zstd);
        std::vector<uint64_t> genes;
        out.Write(kHeader.data(), kHeader.size(), genes);
        // Destroyed without Finish (an interrupted sample): nothing is left behind.
    }
    EXPECT_FALSE(fs::exists(path));
    EXPECT_FALSE(fs::exists(path + SamOutput::kRecordsSuffix));
    {
        SamOutput out(path, SamCompression::None, kHeader);
        out.Discard();
    }
    EXPECT_FALSE(fs::exists(path));

    SamOutput unwritable(dir.File("no such dir/out.sam"), SamCompression::Gzip);
    EXPECT_FALSE(unwritable.Ok());
    EXPECT_NE(unwritable.Error().find("cannot write"), std::string::npos) << unwritable.Error();
    std::vector<uint64_t> genes = { 1 };
    unwritable.Write(kHeader.data(), kHeader.size(), genes);  // skipped, no crash
    EXPECT_FALSE(unwritable.Finish(kHeader));
}

namespace {
    // The blocks of one thread, written in order, so that two SAMs hold the same text.
    std::vector<std::string> WriteInOrder(SamSink& sink, int blocks) {
        std::vector<std::string> written;
        std::vector<uint64_t> genes;
        for (int b = 0; b < blocks; b++) {
            written.push_back(Block(0, b, genes));
            sink.Write(written.back().data(), written.back().size(), genes);
        }
        return written;
    }

    // What both zstd readers make of a SAM: the reading thread's, and frames on threads of their own.
    std::string ReadBothWays(std::string const& path) {
        std::string const text = ReadBack(path);
        SamInput parallel(path, 3);
        EXPECT_EQ(parallel.Problem(), "");
        EXPECT_EQ(ReadAll(parallel.Stream()), text);
        EXPECT_FALSE(parallel.ReadFailed()) << parallel.ReadError();
        return text;
    }

    std::string Incompressible(size_t size, uint64_t seed) {
        std::mt19937_64 rng(seed);
        std::string text(size, '\0');
        for (auto& c : text) c = static_cast<char>(rng());
        return text;
    }
}

TEST(SamFile, ZstdRecordsGoStraightIntoTheSamBehindRoomForTheHeader) {
    ScratchDir dir;
    auto const copied = dir.File("copied.sam.zst"), placed = dir.File("placed.sam.zst");
    size_t const room = sam_zstd::HeaderRoom(5000);
    {
        SamOutput a(copied, SamCompression::Zstd), b(placed, SamCompression::Zstd, std::nullopt, room);
        WriteInOrder(a, 3000);
        WriteInOrder(b, 3000);
        EXPECT_TRUE(fs::exists(copied + SamOutput::kRecordsSuffix));
        EXPECT_FALSE(fs::exists(placed + SamOutput::kRecordsSuffix));  // the records are in the SAM already
        EXPECT_TRUE(fs::exists(placed));
        ASSERT_TRUE(a.Finish(kHeader)) << a.Error();
        ASSERT_TRUE(b.Finish(kHeader)) << b.Error();
    }
    EXPECT_FALSE(fs::exists(placed + SamOutput::kRecordsSuffix));
    EXPECT_EQ(ReadBothWays(placed), ReadBothWays(copied));
    EXPECT_TRUE(sam_zstd::StartsWithMarker(placed));
    // The seek table: the marker, the header, the padding up to the records (all three at the
    // start, the skippable ones without content), then the records' frames as in the copied SAM.
    std::string error;
    auto const ta = zstd::ReadSeekTable(copied, error), tb = zstd::ReadSeekTable(placed, error);
    ASSERT_TRUE(ta.has_value() && tb.has_value()) << error;
    ASSERT_EQ(tb->frames.size(), ta->frames.size() + 1);
    EXPECT_EQ(tb->frames[0].decompressed_size, 0u);
    EXPECT_EQ(tb->frames[2].decompressed_size, 0u);
    uint64_t const records_at = tb->frames[3].compressed_offset;
    EXPECT_EQ(records_at, room);
    for (size_t i = 3; i < tb->frames.size(); i++) {
        EXPECT_EQ(tb->frames[i].compressed_size, ta->frames[i - 1].compressed_size);
        EXPECT_EQ(tb->frames[i].compressed_offset - records_at, ta->frames[i - 1].compressed_offset - ta->frames[2].compressed_offset);
    }
    // The files differ by the padding frame (never written) and its seek table entry.
    EXPECT_EQ(fs::file_size(placed) - fs::file_size(copied), records_at - ta->frames[2].compressed_offset + 8);
}

TEST(SamFile, AZstdHeaderTooLargeForItsRoomMovesTheRecordsBehindIt) {
    ScratchDir dir;
    size_t const room = sam_zstd::HeaderRoom(0);
    size_t const marker = sam_zstd::MarkerFrame().size();
    // Incompressible headers, whose frames are their size plus a few bytes, from well inside the
    // room to past it: the padding frame fits, fits without content, or cannot fit (1-7 bytes).
    for (size_t size = room - marker - 64; size <= room - marker + 8; size++) {
        SCOPED_TRACE(size);
        auto const path = dir.File("out" + std::to_string(size) + ".sam.zst");
        std::string const header = Incompressible(size, size);
        std::vector<std::string> blocks;
        {
            SamOutput out(path, SamCompression::Zstd, std::nullopt, room);
            blocks = WriteInOrder(out, 50);
            ASSERT_TRUE(out.Finish(header)) << out.Error();
        }
        EXPECT_FALSE(fs::exists(path + SamOutput::kRecordsSuffix));
        auto const text = ReadBothWays(path);
        ASSERT_EQ(text.substr(0, header.size()), header);
        ExpectBlocks(text.substr(header.size()), blocks);
        std::string error;
        EXPECT_TRUE(zstd::ReadSeekTable(path, error).has_value()) << error;
    }
}

TEST(SamFile, AZstdSamWithRoomButNoRecords) {
    ScratchDir dir;
    auto const path = dir.File("empty.sam.zst");
    {
        SamOutput out(path, SamCompression::Zstd, std::nullopt, sam_zstd::HeaderRoom(10));
        ASSERT_TRUE(out.Finish(kHeader)) << out.Error();
    }
    EXPECT_EQ(ReadBothWays(path), kHeader);
    EXPECT_LT(fs::file_size(path), 1000u);  // no room left over
    {
        SamOutput out(path, SamCompression::Zstd, std::nullopt, sam_zstd::HeaderRoom(10));
        std::vector<uint64_t> genes;
        out.Write(kHeader.data(), kHeader.size(), genes);
        // Destroyed without Finish: nothing is left behind.
    }
    EXPECT_FALSE(fs::exists(path));
    EXPECT_FALSE(fs::exists(path + SamOutput::kRecordsSuffix));
}

TEST(SamFile, HeaderRoomIsEightBytesAGeneWithinItsBounds) {
    EXPECT_EQ(sam_zstd::HeaderRoom(0), size_t{16} << 10);
    EXPECT_EQ(sam_zstd::HeaderRoom(8192), size_t{64} << 10);
    EXPECT_EQ(sam_zstd::HeaderRoom(10000), size_t{80} << 10);  // 80,000 bytes, in whole 4 KB
    EXPECT_EQ(sam_zstd::HeaderRoom(1000000000), size_t{1} << 20);
    // No more than 1/16 of the read files, when their size is known.
    EXPECT_EQ(sam_zstd::HeaderRoom(1000000000, 1600000), size_t{100} << 10);  // 100,000 bytes, in whole 4 KB
    EXPECT_EQ(sam_zstd::HeaderRoom(1000000000, 1000), size_t{16} << 10);
    EXPECT_EQ(sam_zstd::HeaderRoom(10000, uint64_t{1} << 40), size_t{80} << 10);
    // 4 bytes more a taxon for the unaligned reads' failed candidates, no more than 1/16 of the input, at most 1 MB more.
    EXPECT_EQ(sam_zstd::HeaderRoom(10000, 0, 1000), size_t{84} << 10);  // 84,000 bytes, in whole 4 KB
    EXPECT_EQ(sam_zstd::HeaderRoom(1000000000, 0, 1000000000), size_t{2} << 20);
    EXPECT_EQ(sam_zstd::HeaderRoom(1000000000, 1600000, 1000000), size_t{196} << 10);  // 200,000 bytes, in whole 4 KB
    // The room is for zstd only: other formats collect the records in the temporary file.
    ScratchDir dir;
    auto const path = dir.File("out.sam.gz");
    SamOutput out(path, SamCompression::Gzip, std::nullopt, sam_zstd::HeaderRoom(0));
    EXPECT_TRUE(fs::exists(path + SamOutput::kRecordsSuffix));
    out.Discard();
}

TEST(SamFile, StreamSinkCollectsTheGenes) {
    std::ostringstream os;
    SamStreamSink sink(os);
    auto [blocks, named] = WriteFromThreads(sink);
    ExpectBlocks(os.str(), blocks);
    EXPECT_EQ(sink.Genes(), std::vector<uint64_t>(named.begin(), named.end()));
}
