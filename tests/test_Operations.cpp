// Unit tests for run-level plumbing: the failure collector behind the exit code, and in-place
// compression of the simulator's reads (BGZF, reproducible output).
#include <gtest/gtest.h>
#include <zlib-ng.h>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <random>
#include <sstream>
#include <string>
#include <unistd.h>
#include "Utilities/Compressor.h"
#include "Utilities/RunStatus.h"

namespace fs = std::filesystem;

namespace {
    // A scratch directory whose name contains a space, removed afterwards.
    struct ScratchDir {
        fs::path path;
        ScratchDir() {
            path = fs::temp_directory_path() / ("protal ops test " + std::to_string(::getpid()));
            fs::create_directories(path);
        }
        ~ScratchDir() { fs::remove_all(path); }
    };

    std::string ReadAll(fs::path const& p) {
        std::ifstream in(p, std::ios::binary);
        return std::string(std::istreambuf_iterator<char>(in), {});
    }

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

    // ~10 MB of FASTQ-like text: several of the compressor's 4 MB chunks.
    std::string Reads() {
        std::mt19937 rng(3);
        std::string text;
        while (text.size() < (size_t{10} << 20)) {
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

TEST(Compressor, CompressesInPlaceToBgzf) {
    ScratchDir dir;
    auto reads = dir.path / "sample $(touch INJECTED); x_R1.fq";  // no shell is involved
    std::string const text = Reads();
    std::ofstream(reads, std::ios::binary) << text;

    Compressor::compressInPlace(reads, 3);

    fs::path const gz = reads.string() + ".gz";
    EXPECT_FALSE(fs::exists(reads));
    ASSERT_TRUE(fs::exists(gz));
    EXPECT_FALSE(fs::exists(dir.path / "INJECTED"));
    EXPECT_EQ(Gunzip(gz), text);
    EXPECT_TRUE(protal::bgzf::StartsAsBgzf(gz.string()));
    EXPECT_TRUE(protal::bgzf::EndsWithEof(gz.string()));
    EXPECT_LT(fs::file_size(gz), text.size() / 2);
}

TEST(Compressor, SameContentGivesIdenticalBytesWithAnyThreads) {
    // No name or time stamp in the gzip headers, and the blocks do not depend on the threads.
    ScratchDir dir;
    std::string const text = Reads();
    std::vector<std::string> packed;
    for (int threads : { 1, 2, 5 }) {
        auto file = dir.path / ("reads_" + std::to_string(threads) + ".fq");
        std::ofstream(file, std::ios::binary) << text;
        fs::last_write_time(file, fs::file_time_type::clock::now() - std::chrono::hours(threads));
        Compressor::compressInPlace(file, threads);
        packed.push_back(ReadAll(file.string() + ".gz"));
    }
    EXPECT_EQ(packed[0], packed[1]);
    EXPECT_EQ(packed[0], packed[2]);
}

TEST(Compressor, AnEmptyFileAndAMissingOne) {
    ScratchDir dir;
    auto empty = dir.path / "empty.fq";
    std::ofstream(empty).close();
    Compressor::compressInPlace(empty, 2);
    EXPECT_EQ(Gunzip(empty.string() + ".gz"), "");
    EXPECT_TRUE(protal::bgzf::EndsWithEof(empty.string() + ".gz"));

    EXPECT_THROW(Compressor::compressInPlace(dir.path / "missing.fq", 1), std::invalid_argument);
    EXPECT_FALSE(fs::exists(dir.path / "missing.fq.gz"));
}

TEST(Compressor, TheStreamingWriterGivesTheSameBytes) {
    // bgzf::Writer, which the simulator appends each genome's reads to: the blocks start every kBlockInput bytes of
    // the content however it arrives, so the file is what compressing the whole content gives.
    ScratchDir dir;
    std::string text = Reads();
    while (text.size() < 3 * protal::bgzf::kBlockInput * 64) text += text;  // several flushes
    auto whole = dir.path / "whole.fq";
    std::ofstream(whole, std::ios::binary) << text;
    Compressor::compressInPlace(whole, 2);
    std::string const expected = ReadAll(whole.string() + ".gz");
    for (size_t piece : { size_t{1} << 20, size_t{7}, protal::bgzf::kBlockInput, size_t{100003} }) {
        auto out = dir.path / ("streamed_" + std::to_string(piece) + ".fq.gz");
        protal::bgzf::Writer writer(out.string());
        for (size_t at = 0; at < text.size(); at += piece) writer.Write(text.data() + at, std::min(piece, text.size() - at));
        ASSERT_TRUE(writer.Close()) << writer.Error();
        EXPECT_EQ(ReadAll(out.string()), expected) << "pieces of " << piece;
    }
    // Nothing written: only the end-of-file block, as for an empty file.
    auto empty = dir.path / "empty.fq.gz";
    protal::bgzf::Writer writer(empty.string());
    ASSERT_TRUE(writer.Close());
    EXPECT_EQ(Gunzip(empty), "");
    EXPECT_TRUE(protal::bgzf::EndsWithEof(empty.string()));
    // A file that cannot be written says so.
    protal::bgzf::Writer bad((dir.path / "no_such_dir" / "x.gz").string());
    EXPECT_FALSE(bad.Error().empty());
    EXPECT_FALSE(bad.Close());
}
