// Unit tests for run-level plumbing: the failure collector behind the exit code, and in-place SAM
// compression (no shell, reproducible output).
#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <sys/stat.h>
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

    // Stand-in for pigz that needs only gzip: drops "-p N" and passes everything else on.
    std::string WritePigzStandIn(fs::path const& dir) {
        auto script = dir / "pigz-standin";
        std::ofstream(script) << "#!/usr/bin/env bash\n"
                                 "args=()\n"
                                 "while [ $# -gt 0 ]; do\n"
                                 "  if [ \"$1\" = \"-p\" ]; then shift 2; continue; fi\n"
                                 "  args+=(\"$1\"); shift\n"
                                 "done\n"
                                 "exec gzip \"${args[@]}\"\n";
        ::chmod(script.c_str(), 0755);
        return script.string();
    }

    std::string ReadAll(fs::path const& p) {
        std::ifstream in(p, std::ios::binary);
        return std::string(std::istreambuf_iterator<char>(in), {});
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

TEST(Compressor, HandlesShellMetacharactersWithoutAShell) {
    ScratchDir dir;
    auto pigz = WritePigzStandIn(dir.path);
    auto sam = dir.path / "sample $(touch INJECTED); x.sam";
    std::ofstream(sam) << "@HD\tVN:1.6\n";

    Compressor::compressInPlace(sam, 2, pigz);

    EXPECT_FALSE(fs::exists(sam));
    EXPECT_TRUE(fs::exists(sam.string() + ".gz"));
    EXPECT_FALSE(fs::exists(dir.path / "INJECTED"));
    EXPECT_FALSE(fs::exists("INJECTED"));
}

TEST(Compressor, SameContentGivesIdenticalBytes) {
    // -n keeps name and time stamp out of the gzip header.
    ScratchDir dir;
    auto pigz = WritePigzStandIn(dir.path);
    auto a = dir.path / "a.sam";
    auto b = dir.path / "b_other_name.sam";
    std::ofstream(a) << "r1\t0\t1_1\t1\t60\t4M\t*\t0\t4\tACGT\tIIII\n";
    std::ofstream(b) << "r1\t0\t1_1\t1\t60\t4M\t*\t0\t4\tACGT\tIIII\n";
    fs::last_write_time(b, fs::last_write_time(a) - std::chrono::hours(5));

    Compressor::compressInPlace(a, 1, pigz);
    Compressor::compressInPlace(b, 1, pigz);
    EXPECT_EQ(ReadAll(a.string() + ".gz"), ReadAll(b.string() + ".gz"));
}

TEST(Compressor, ReportsAMissingProgram) {
    ScratchDir dir;
    auto sam = dir.path / "a.sam";
    std::ofstream(sam) << "x\n";
    EXPECT_THROW(Compressor::compressInPlace(sam, 1, "no-such-pigz-binary"), std::runtime_error);
    EXPECT_TRUE(fs::exists(sam));  // the input is left alone
}
