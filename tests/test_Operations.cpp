// Unit tests for run-level plumbing: the failure collector behind the exit code, and in-place
// compression of the simulator's reads (BGZF, reproducible output).
#include <gtest/gtest.h>
#include <zlib.h>
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

    // A gzip file's content as zlib (zcat, gzip) reads it.
    std::string Gunzip(fs::path const& p) {
        gzFile f = gzopen(p.c_str(), "rb");
        std::string text;
        char buffer[1 << 16];
        int n;
        while ((n = gzread(f, buffer, sizeof(buffer))) > 0) text.append(buffer, static_cast<size_t>(n));
        gzclose(f);
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
