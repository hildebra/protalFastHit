// Unit tests for run-level plumbing: the failure collector behind the exit code, and the BGZF writer
// the simulator writes its reads with (reproducible output).
#include <gtest/gtest.h>
#include <zlib-ng.h>
#include <algorithm>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <random>
#include <sstream>
#include <string>
#include <unistd.h>
#include "IO/Bgzf.h"
#include "Utilities/RunStatus.h"
#include "TestUtil.h"

namespace fs = std::filesystem;
using protal::test::ScratchDir;
using protal::test::Slurp;

namespace {
    // A gzip file's content as a reader other than protal's (ISA-L) reads it (zlib-ng; zcat reads the same).
    std::string Gunzip(fs::path const& p) {
        gzFile f = zng_gzopen(p.c_str(), "rb");
        std::string text;
        char buffer[1 << 16];
        int32_t n;
        while ((n = zng_gzread(f, buffer, sizeof(buffer))) > 0) text.append(buffer, static_cast<size_t>(n));
        zng_gzclose(f);
        return text;
    }

    // FASTQ-like text of at least `size` bytes.
    std::string Reads(size_t size) {
        std::mt19937 rng(3);
        std::string text;
        while (text.size() < size) {
            std::string seq(150, 'A');
            for (auto& c : seq) c = "ACGT"[rng() % 4];
            text += "@read" + std::to_string(text.size()) + "\n" + seq + "\n+\n" + std::string(150, 'I') + "\n";
        }
        return text;
    }
}

TEST(RunStatus, ExitCodeReflectsFailures) {
    protal::RunStatus status;
    std::ostringstream report;
    EXPECT_TRUE(status.Ok());
    EXPECT_EQ(status.Finish(report), 0);
    EXPECT_TRUE(report.str().empty());

    status.Fail("sample A: no SAM file");
    status.Fail("qcmsa failed for B");
    EXPECT_FALSE(status.Ok());
    EXPECT_EQ(status.Finish(report), 1);
    EXPECT_NE(report.str().find("2 error(s)"), std::string::npos);
    EXPECT_NE(report.str().find("qcmsa failed for B"), std::string::npos);
}

// bgzf::Writer, which the simulator appends each genome's reads to: a block starts every kBlockInput bytes of the
// content however it arrives, so pieces of any size (across the writer's ~4 MB flushes) give the same file. That file
// is BGZF (blocks of at most 64 KB with their size, then the end-of-file block) and gzip readers read the content back.
TEST(BgzfWriter, PiecesOfAnySizeGiveTheSameFile) {
    namespace bgzf = protal::bgzf;
    ScratchDir dir;
    std::string const text = Reads(65 * bgzf::kBlockInput + 12345);  // one flush of 64 blocks, then more
    std::string first;
    for (size_t piece : { text.size(), size_t{7}, bgzf::kBlockInput, size_t{100003} }) {
        auto const out = dir.path / ("streamed_" + std::to_string(piece) + ".fq.gz");
        bgzf::Writer writer(out.string());
        for (size_t at = 0; at < text.size(); at += piece) writer.Write(text.data() + at, std::min(piece, text.size() - at));
        ASSERT_TRUE(writer.Close()) << writer.Error();
        std::string const bytes = Slurp(out);
        if (first.empty()) first = bytes;
        EXPECT_EQ(bytes, first) << "pieces of " << piece;
    }

    // Its blocks: each a BGZF member whose size field is right, kBlockInput bytes of content but the last; then the
    // end-of-file block.
    size_t offset = 0, blocks = 0, content = 0;
    while (offset + sizeof(bgzf::kEof) < first.size()) {
        auto const* b = reinterpret_cast<unsigned char const*>(first.data() + offset);
        ASSERT_EQ(std::memcmp(b, bgzf::kHeader, sizeof(bgzf::kHeader)), 0) << "block " << blocks;
        size_t const size = bgzf::BlockSize(b);
        ASSERT_LE(size, bgzf::kMaxBlock);
        ASSERT_LE(offset + size, first.size());
        size_t const isize = bgzf::LE32(b + size - 4);
        EXPECT_EQ(isize, std::min(bgzf::kBlockInput, text.size() - content)) << "block " << blocks;
        content += isize;
        offset += size;
        blocks++;
    }
    EXPECT_EQ(content, text.size());
    EXPECT_EQ(blocks, (text.size() + bgzf::kBlockInput - 1) / bgzf::kBlockInput);
    EXPECT_EQ(first.compare(offset, std::string::npos, reinterpret_cast<char const*>(bgzf::kEof), sizeof(bgzf::kEof)), 0);
    EXPECT_EQ(Gunzip(dir.path / ("streamed_" + std::to_string(text.size()) + ".fq.gz")), text);

    // Nothing written: only the end-of-file block.
    auto const empty = dir.path / "empty.fq.gz";
    bgzf::Writer writer(empty.string());
    ASSERT_TRUE(writer.Close());
    EXPECT_EQ(Slurp(empty), std::string(reinterpret_cast<char const*>(bgzf::kEof), sizeof(bgzf::kEof)));
    EXPECT_EQ(Gunzip(empty), "");
    // A file that cannot be written says so.
    bgzf::Writer bad((dir.path / "no_such_dir" / "x.gz").string());
    EXPECT_FALSE(bad.Error().empty());
    EXPECT_FALSE(bad.Close());
}
