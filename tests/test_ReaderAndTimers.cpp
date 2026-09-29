// Unit tests for the input stream that inflates in a thread of its own (ThreadedGzStream), the
// read pairs threads take from it (SeqReaderPE), and the stage timers (Benchmark).
#include <gtest/gtest.h>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <omp.h>
#include <unistd.h>
#include "gzstream.h"
#include "IO/ThreadedGzStream.h"
#include "SequenceUtils/SeqReader.h"
#include "Utilities/Benchmark.h"

using namespace protal;
namespace fs = std::filesystem;

namespace {
    struct ScratchDir {
        fs::path path;
        ScratchDir() {
            path = fs::temp_directory_path() / ("protal reader test " + std::to_string(::getpid()));
            fs::create_directories(path);
        }
        ~ScratchDir() { fs::remove_all(path); }

        std::string Plain(std::string const& name, std::string const& content) const {
            auto file = path / name;
            std::ofstream(file, std::ios::binary) << content;
            return file.string();
        }

        std::string Gzip(std::string const& name, std::string const& content) const {
            auto file = (path / name).string();
            ogzstream os(file.c_str());
            os << content;
            os.close();
            return file;
        }
    };

    // A FASTQ file of n reads named <prefix><i>, several MB for n in the ten thousands.
    std::string Fastq(size_t n, std::string const& prefix) {
        std::string const bases = "ACGTTGCAAGCTAGCTTACGGATCCGATTACAGGCATCGATCGGCTAGCATCGACTAGCATTACGACTACGATCAGCTACGA";
        std::ostringstream os;
        for (size_t i = 0; i < n; i++) {
            std::string seq;
            for (size_t j = 0; j < 150; j++) seq += bases[(i * 7 + j * 13) % bases.size()];
            os << '@' << prefix << i << "/1\n" << seq << "\n+\n" << std::string(150, 'I') << '\n';
        }
        return os.str();
    }

    std::vector<std::string> Lines(std::istream& is) {
        std::vector<std::string> lines;
        std::string line;
        while (std::getline(is, line)) lines.push_back(line);
        return lines;
    }
}

TEST(ThreadedGzStream, ReadsAGzipFileAsIgzstreamDoes) {
    ScratchDir dir;
    auto const content = Fastq(40000, "r");  // ~12 MB: each block is filled and reused several times
    ASSERT_GT(content.size(), 2 * ThreadedGzStreambuf::kBlocks * ThreadedGzStreambuf::kBlockSize);
    auto const path = dir.Gzip("reads.fq.gz", content);

    ThreadedGzIstream threaded(path.c_str());
    igzstream reference(path.c_str());
    EXPECT_EQ(Lines(threaded), Lines(reference));
    EXPECT_FALSE(threaded.rdbuf()->read_failed());
}

TEST(ThreadedGzStream, ReadsAPlainFile) {
    ScratchDir dir;
    auto const content = Fastq(5000, "p");
    ThreadedGzIstream is(dir.Plain("reads.fq", content).c_str());
    std::stringstream expected(content);
    EXPECT_EQ(Lines(is), Lines(expected));
    EXPECT_FALSE(is.rdbuf()->read_failed());
}

TEST(ThreadedGzStream, ATruncatedFileReadsAsAPrefixAndIsReported) {
    ScratchDir dir;
    auto const content = Fastq(40000, "t");
    auto const path = dir.Gzip("reads.fq.gz", content);
    fs::resize_file(path, fs::file_size(path) / 2);

    ThreadedGzIstream is(path.c_str());
    std::string read((std::istreambuf_iterator<char>(is)), std::istreambuf_iterator<char>());
    EXPECT_TRUE(is.rdbuf()->read_failed());
    EXPECT_FALSE(is.rdbuf()->read_error_message().empty());
    EXPECT_GT(read.size(), 0u);
    EXPECT_LT(read.size(), content.size());
    EXPECT_EQ(content.compare(0, read.size(), read), 0);

    igzstream reference(path.c_str());
    std::string ignored((std::istreambuf_iterator<char>(reference)), std::istreambuf_iterator<char>());
    EXPECT_TRUE(reference.rdbuf()->read_failed());  // as igzstream reports it
}

TEST(ThreadedGzStream, AMissingFileFailsToOpen) {
    ScratchDir dir;
    ThreadedGzIstream is((dir.path / "none.fq.gz").string().c_str());
    EXPECT_FALSE(is.good());
    EXPECT_FALSE(is.rdbuf()->is_open());
    std::string line;
    EXPECT_FALSE(std::getline(is, line));
    EXPECT_FALSE(is.rdbuf()->read_failed());
}

TEST(ThreadedGzStream, ClosingBeforeTheEndStopsTheInflatingThread) {
    ScratchDir dir;
    auto const path = dir.Gzip("reads.fq.gz", Fastq(40000, "c"));
    for (int i = 0; i < 20; i++) {  // the inflating thread is blocked on full blocks or still inflating
        ThreadedGzIstream is(path.c_str());
        std::string line;
        ASSERT_TRUE(std::getline(is, line));
        EXPECT_EQ(line, "@c0/1");
        if (i % 2) is.close();  // else the destructor closes it
    }
}

TEST(ThreadedGzStream, ThreadsTakeEveryPairOnceInStep) {
    ScratchDir dir;
    size_t const n = 30000;
    auto const r1 = dir.Gzip("r1.fq.gz", Fastq(n, "pair"));
    auto const r2 = dir.Gzip("r2.fq.gz", Fastq(n, "pair"));
    ThreadedGzIstream is1(r1.c_str()), is2(r2.c_str());
    SeqReaderPE reader_global{ is1, is2 };

    std::vector<std::string> ids;
    size_t out_of_step = 0;
#pragma omp parallel num_threads(4) shared(reader_global, ids, out_of_step) default(none)
    {
        SeqReaderPE reader{ reader_global };
        FastxRecord record1, record2;
        std::vector<std::string> mine;
        size_t mismatches = 0;
        while (reader(record1, record2)) {
            mismatches += record1.id != record2.id || record1.sequence != record2.sequence;
            mine.push_back(record1.id);
        }
#pragma omp critical(test_ids)
        {
            ids.insert(ids.end(), mine.begin(), mine.end());
            out_of_step += mismatches;
            reader_global.UpdateSuccess(reader);
        }
    }
    EXPECT_EQ(out_of_step, 0u);
    EXPECT_TRUE(reader_global.Success());
    EXPECT_EQ(ids.size(), n);
    EXPECT_EQ(std::set<std::string>(ids.begin(), ids.end()).size(), n);
    EXPECT_FALSE(is1.rdbuf()->read_failed() || is2.rdbuf()->read_failed());
}

TEST(Benchmark, SumsIntervalsShorterThanAMicrosecond) {
    Benchmark bm{"short"};
    for (int i = 0; i < 1000; i++) {
        bm.Start();
        bm.Stop();
    }
    EXPECT_GT(bm.Seconds(), 0.0);  // flooring each interval to whole microseconds gave 0
    EXPECT_EQ(bm.Threads(), 1u);
}

TEST(Benchmark, PrintingKeepsTheSum) {
    Benchmark global{"stage", 0};
    for (int t = 0; t < 4; t++) {
        Benchmark local{"stage"};
        local.Start();
        usleep(2000);
        local.Stop();
        global.Join(local);
    }
    EXPECT_EQ(global.Threads(), 4u);
    double const seconds = global.Seconds();
    EXPECT_GE(seconds, 0.008);
    EXPECT_DOUBLE_EQ(global.MeanSeconds(), seconds / 4);
    testing::internal::CaptureStdout();
    global.PrintResults();
    global.PrintResults();
    auto const printed = testing::internal::GetCapturedStdout();
    EXPECT_DOUBLE_EQ(global.Seconds(), seconds);  // printing divided the sum by the threads, each time
    EXPECT_NE(printed.find("stage took"), std::string::npos) << printed;
    EXPECT_NE(printed.find("mean over 4 threads"), std::string::npos) << printed;
    EXPECT_NEAR(static_cast<double>(global.GetDuration(Time::milliseconds)), seconds * 1000, 1.0);
}
