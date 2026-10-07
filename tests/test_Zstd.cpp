// Unit tests for zstd-compressed database files: the compressing and decompressing streams,
// failure modes (truncated, corrupt, trailing data), rewinding, CompressFile, the seekable format
// with parallel reading (index scatter), and loading a reference.fna.zst without a seek table with
// GenomeLoader (GenomeLoaderPacked tests the seekable one).
#include <gtest/gtest.h>
#include <atomic>
#include <chrono>
#include <cstring>
#include <mutex>
#include <thread>
#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <sysexits.h>
#include <unistd.h>
#include "Utilities/Zstd.h"
#include "Hash/Seedmap.h"
#include "SequenceUtils/GenomeLoader.h"
#include "TestUtil.h"

namespace fs = std::filesystem;
using namespace protal;
using namespace protal::test;

namespace {
    void WriteCompressed(std::string const& path, std::string const& data, zstd::Params params = {3, 0, 1},
                         bool pledge = true) {
        zstd::OStream os(path, params, pledge ? std::optional<uint64_t>(data.size()) : std::nullopt);
        // Mixed write sizes: small writes go through the buffer, large ones straight to zstd.
        size_t pos = 0, step = 1;
        while (pos < data.size()) {
            size_t const n = std::min(step, data.size() - pos);
            os.write(data.data() + pos, static_cast<std::streamsize>(n));
            pos += n;
            step = step < (size_t{4} << 20) ? step * 7 : 1;
        }
        ASSERT_TRUE(os.Close()) << os.Buffer().Error();
    }

    std::string ReadAll(zstd::InputFile& in, size_t chunk) {
        std::string out, buffer(chunk, '\0');
        while (in.Stream().read(buffer.data(), static_cast<std::streamsize>(chunk)) || in.Stream().gcount() > 0) {
            out.append(buffer.data(), static_cast<size_t>(in.Stream().gcount()));
            if (!in.Stream()) break;
        }
        return out;
    }

}

// 3 MB: WriteCompressed's writes grow past the stream's 1 MB buffer, so the last goes to zstd unbuffered.
TEST(Zstd, RoundTripWithSmallAndLargeReads) {
    ScratchDir tmp;
    for (size_t size : {size_t{1}, size_t{1000}, size_t{3} << 20}) {
        std::string const data = TestData(size, static_cast<unsigned>(size));
        WriteCompressed(tmp / "a.zst", data);
        EXPECT_TRUE(zstd::IsCompressed(tmp / "a.zst"));
        for (size_t chunk : {size_t{1}, size_t{4096}, size_t{64} << 20}) {
            if (chunk == 1 && size > 100000) continue;  // byte by byte only for the small ones
            zstd::InputFile in(tmp / "a.zst");
            ASSERT_TRUE(in.IsOpen());
            EXPECT_TRUE(in.Compressed());
            EXPECT_EQ(ReadAll(in, chunk), data) << "size " << size << ", chunk " << chunk;
            EXPECT_FALSE(in.Stream().bad());
        }
    }
}

TEST(Zstd, MultithreadedLongWindowAndUnpledgedSize) {
    // Without multithreading, libzstd compresses on one thread and the stream only warns: nothing here would test it.
    if (ZSTD_cParam_getBounds(ZSTD_c_nbWorkers).upperBound == 0) GTEST_SKIP() << "this libzstd has no multithreading";
    ScratchDir tmp;
    std::string const data = TestData(size_t{4} << 20, 7);
    for (bool pledge : {true, false}) {
        WriteCompressed(tmp / "b.zst", data, {9, 27, 4}, pledge);
        EXPECT_LT(fs::file_size(tmp / "b.zst"), data.size() / 3);
        zstd::InputFile in(tmp / "b.zst");
        EXPECT_EQ(ReadAll(in, size_t{1} << 20), data);
    }
}

TEST(Zstd, TellAndGetline) {
    ScratchDir tmp;
    std::string const data = ">1_1\nACGT\n>1_2\nGGGCCC\n";
    WriteCompressed(tmp / "c.zst", data);
    zstd::InputFile in(tmp / "c.zst");
    std::string line;
    std::vector<long> positions;
    while (std::getline(in.Stream(), line)) positions.push_back(static_cast<long>(in.Stream().tellg()));
    EXPECT_EQ(positions, (std::vector<long>{5, 10, 15, 22}));
}

TEST(Zstd, RewindRawAndCompressed) {
    ScratchDir tmp;
    std::string const data = TestData(size_t{10} << 20);
    Spit(tmp / "raw", data);
    WriteCompressed(tmp / "raw.zst", data);
    for (auto name : {"raw", "raw.zst"}) {
        zstd::InputFile in(tmp / name);
        EXPECT_EQ(ReadAll(in, 1 << 16), data);
        ASSERT_TRUE(in.Rewind());
        EXPECT_EQ(ReadAll(in, 1 << 16), data) << name;
    }
}

TEST(Zstd, ConcatenatedFramesAreOneStream) {
    ScratchDir tmp;
    WriteCompressed(tmp / "x.zst", "first part|");
    WriteCompressed(tmp / "y.zst", "second part");
    Spit(tmp / "xy.zst", Slurp(tmp / "x.zst") + Slurp(tmp / "y.zst"));
    zstd::InputFile in(tmp / "xy.zst");
    EXPECT_EQ(ReadAll(in, 7), "first part|second part");
}

// The content size comes from the frame header, or is counted when the file starts with a
// skippable frame (as pzstd writes them) or the header does not record it.
TEST(Zstd, UncompressedSize) {
    ScratchDir tmp;
    std::string const data = TestData(300000);
    Spit(tmp / "raw", data);
    WriteCompressed(tmp / "pledged.zst", data);
    WriteCompressed(tmp / "unpledged.zst", data, {3, 0, 1}, false);
    std::string skippable = "\x50\x2a\x4d\x18";
    skippable += std::string("\x05\x00\x00\x00", 4) + "hello";
    Spit(tmp / "skippable.zst", skippable + Slurp(tmp / "pledged.zst"));
    for (auto name : {"raw", "pledged.zst", "unpledged.zst", "skippable.zst"}) {
        EXPECT_EQ(zstd::UncompressedSize(tmp / name), std::optional<uint64_t>(data.size())) << name;
    }
    EXPECT_TRUE(zstd::IsCompressed(tmp / "skippable.zst"));
    zstd::InputFile in(tmp / "skippable.zst");
    EXPECT_EQ(ReadAll(in, 1 << 16), data);
    EXPECT_EQ(zstd::UncompressedSize(tmp / "missing"), std::nullopt);
}

TEST(Zstd, TruncatedCorruptAndTrailingDataFail) {
    ScratchDir tmp;
    std::string const data = TestData(size_t{2} << 20);
    WriteCompressed(tmp / "good.zst", data);
    std::string const good = Slurp(tmp / "good.zst");

    Spit(tmp / "truncated.zst", good.substr(0, good.size() / 2));
    std::string corrupt = good;
    corrupt[corrupt.size() / 2] ^= 0x5a;
    Spit(tmp / "corrupt.zst", corrupt);
    Spit(tmp / "trailing.zst", good + "garbage");

    // What a truncated or corrupt file gives before the error is less than its content, and its start.
    for (auto name : {"truncated.zst", "corrupt.zst"}) {
        zstd::InputFile in(tmp / name);
        std::string const read = ReadAll(in, 1 << 16);
        EXPECT_TRUE(in.Stream().bad()) << name << " was read without an error";
        EXPECT_LT(read.size(), data.size()) << name;
        EXPECT_EQ(data.compare(0, read.size(), read), 0) << name << ": what was read is not the content's start";
    }
    // Bytes after the frame that are no frame: the whole content, then the error.
    zstd::InputFile trailing(tmp / "trailing.zst");
    EXPECT_EQ(ReadAll(trailing, 1 << 16), data);
    EXPECT_TRUE(trailing.Stream().bad()) << "trailing.zst was read without an error";
}

TEST(Zstd, PledgedSizeMismatchFailsOnClose) {
    ScratchDir tmp;
    zstd::OStream os(tmp / "short.zst", {3, 0, 1}, 100);
    os << "only a few bytes";
    EXPECT_FALSE(os.Close());
}

TEST(Zstd, ResolvePrefersTheRawFile) {
    ScratchDir tmp;
    EXPECT_EQ(zstd::Resolve(tmp / "f"), tmp / "f");  // neither: the raw name
    Spit(tmp / "f.zst", "z");
    EXPECT_EQ(zstd::Resolve(tmp / "f"), tmp / "f.zst");
    Spit(tmp / "f", "r");
    EXPECT_EQ(zstd::Resolve(tmp / "f"), tmp / "f");
}

TEST(Zstd, CompressFileVerifiesAndRenames) {
    ScratchDir tmp;
    std::string const data = TestData(size_t{5} << 20);
    Spit(tmp / "reference.fna", data);
    std::string error;
    ASSERT_TRUE(zstd::CompressFile(tmp / "reference.fna", tmp / "reference.fna.zst", {3, 27, 2}, true, error)) << error;
    EXPECT_FALSE(fs::exists(tmp / "reference.fna.zst.partial"));
    zstd::InputFile in(tmp / "reference.fna.zst");
    EXPECT_EQ(ReadAll(in, 1 << 20), data);

    EXPECT_FALSE(zstd::CompressFile(tmp / "missing.fna", tmp / "missing.fna.zst", {3, 0, 1}, true, error));
    EXPECT_FALSE(fs::exists(tmp / "missing.fna.zst"));
}

// CompressFile's check of what it wrote: a source of two frames without a seek table (cat a.zst b.zst) has the first
// frame's size in its header, which UncompressedSize takes for the whole content; the check finds more, CompressFile
// fails, and neither the file nor its partial copy is left.
TEST(Zstd, CompressFileFailsWhenWhatItWroteIsNotTheSource) {
    ScratchDir tmp;
    WriteCompressed(tmp / "a.zst", "first part|");
    WriteCompressed(tmp / "b.zst", "second part");
    Spit(tmp / "ab.zst", Slurp(tmp / "a.zst") + Slurp(tmp / "b.zst"));
    ASSERT_EQ(zstd::UncompressedSize(tmp / "ab.zst"), std::optional<uint64_t>(11));
    std::string error;
    EXPECT_FALSE(zstd::CompressFile(tmp / "ab.zst", tmp / "out.zst", {3, 0, 1, 1 << 20}, true, error));
    EXPECT_NE(error.find("verifying"), std::string::npos) << error;
    EXPECT_NE(error.find("22 of 11 bytes compared"), std::string::npos) << error;
    EXPECT_FALSE(fs::exists(tmp / "out.zst"));
    EXPECT_FALSE(fs::exists(tmp / "out.zst.partial"));
}

namespace {
    // Collects ParallelRead's output in one buffer (optionally offering it as Direct memory).
    struct BufferSink : zstd::Sink {
        std::string data;
        bool direct;
        explicit BufferSink(size_t size, bool direct = false) : data(size, '\0'), direct(direct) {}
        char* Direct(uint64_t offset, size_t size) override {
            return direct && offset + size <= data.size() ? data.data() + offset : nullptr;
        }
        void Copy(uint64_t offset, char const* src, size_t size) override {
            ASSERT_LE(offset + size, data.size());
            std::memcpy(data.data() + offset, src, size);
        }
    };

    // The parts one after another, written in the seekable format in frames of frame_size bytes (CompressFrames);
    // returns them.
    std::string WriteSeekable(std::string const& path, std::vector<std::string> const& parts, uint64_t frame_size,
                              int threads = 4, int level = 3) {
        std::string all;
        for (auto const& p : parts) all += p;
        std::istringstream is(all);
        zstd::StreamReader reader(is);
        std::string error;
        auto const written = zstd::CompressFrames(reader, path, {level, 0, threads, frame_size}, error);
        EXPECT_TRUE(written.has_value()) << error;
        EXPECT_EQ(written.value_or(0), fs::file_size(path));
        return all;
    }
}

TEST(ZstdSeekable, RoundTripWithAnyThreadCount) {
    ScratchDir tmp;
    // Parts of odd sizes: frames do not line up with them.
    std::string const data = WriteSeekable(tmp / "s.zst", {TestData(12345, 1), TestData(size_t{5} << 20, 2), TestData(777, 3),
                                                           TestData(size_t{4} << 20, 4)}, size_t{1} << 20);
    std::string error;
    auto const table = zstd::ReadSeekTable(tmp / "s.zst", error);
    ASSERT_TRUE(table) << error;
    EXPECT_EQ(table->frames.size(), (data.size() + (1 << 20) - 1) >> 20);
    EXPECT_EQ(table->DecompressedSize(), data.size());
    for (size_t i = 0; i + 1 < table->frames.size(); i++) EXPECT_EQ(table->frames[i].decompressed_size, size_t{1} << 20);
    EXPECT_TRUE(zstd::IsSeekable(tmp / "s.zst"));
    EXPECT_EQ(zstd::UncompressedSize(tmp / "s.zst"), std::optional<uint64_t>(data.size()));

    for (int threads : {1, 3, 8, 64}) {
        for (bool direct : {false, true}) {
            BufferSink sink(data.size(), direct);
            EXPECT_EQ(zstd::ParallelRead(tmp / "s.zst", threads, sink, error), data.size());
            EXPECT_TRUE(error.empty()) << error;
            EXPECT_EQ(sink.data, data) << threads << " threads, direct " << direct;
        }
    }
    // Plain sequential readers (and the zstd CLI) see the same content: frames, then a skippable frame.
    zstd::InputFile in(tmp / "s.zst");
    EXPECT_EQ(ReadAll(in, 1 << 16), data);
}

TEST(ZstdSeekable, ParallelReadOfRawAndSingleFrameFiles) {
    ScratchDir tmp;
    std::string const data = TestData(size_t{3} << 20);
    Spit(tmp / "raw", data);
    WriteCompressed(tmp / "single.zst", data);
    EXPECT_FALSE(zstd::IsSeekable(tmp / "single.zst"));
    for (auto name : {"raw", "single.zst"}) {
        std::string error;
        BufferSink sink(data.size(), true);
        EXPECT_EQ(zstd::ParallelRead(tmp / name, 4, sink, error), data.size()) << name;
        EXPECT_TRUE(error.empty()) << error;
        EXPECT_EQ(sink.data, data) << name;
    }
}

TEST(ZstdSeekable, CorruptFramesAndTablesFail) {
    ScratchDir tmp;
    std::string const data = WriteSeekable(tmp / "good.zst", {TestData(size_t{4} << 20)}, size_t{1} << 20);
    std::string const good = Slurp(tmp / "good.zst");
    std::string error;
    auto const table = zstd::ReadSeekTable(tmp / "good.zst", error);
    ASSERT_TRUE(table);

    std::string corrupt = good;  // a byte in the middle of frame 3
    corrupt[table->frames[2].compressed_offset + table->frames[2].compressed_size / 2] ^= 0x5a;
    Spit(tmp / "corrupt.zst", corrupt);
    std::string wrong_size = good;  // frame 2's decompressed size in the seek table
    size_t const entry = good.size() - 9 - (table->frames.size() - 1) * 8 + 4;
    wrong_size[entry] ^= 0x01;
    Spit(tmp / "wrong_size.zst", wrong_size);
    Spit(tmp / "truncated.zst", good.substr(0, good.size() - 20));  // no seek table: read as one stream

    for (auto name : {"corrupt.zst", "wrong_size.zst", "truncated.zst"}) {
        error.clear();
        BufferSink sink(data.size() + 16);
        zstd::ParallelRead(tmp / name, 4, sink, error);
        EXPECT_FALSE(error.empty()) << name << " was read without an error";
    }
}

namespace {
    // A stream read to its end with read() (which turns an exception of the buffer into badbit).
    std::string ReadToEnd(std::istream& in) {
        std::vector<char> buffer(size_t{1} << 16);
        std::string all;
        while (in.read(buffer.data(), static_cast<std::streamsize>(buffer.size())) || in.gcount() > 0) {
            all.append(buffer.data(), static_cast<size_t>(in.gcount()));
        }
        return all;
    }
}

// ForEachFrame with a budget (a load into memory that becomes resident as it is written): a worker takes another frame
// only while the output of the frames not yet started holds the buffers of the workers decoding, so near the end the
// workers stop one by one; every frame is still handled once, with its content. With an ample budget only the last frame
// waits for the others (nothing is left to write after it); without one, no worker stops early.
TEST(ZstdSeekable, ForEachFrameStopsWorkersNearTheEndOfALoad) {
    ScratchDir tmp;
    std::string const data = WriteSeekable(tmp / "load.zst", {TestData(size_t{3} << 20, 7)}, size_t{1} << 16);  // 48 frames
    std::string error;
    auto const table = zstd::ReadSeekTable(tmp / "load.zst", error);
    ASSERT_TRUE(table) << error;
    size_t const frames = table->frames.size();
    ASSERT_GE(frames, 40u);
    uint64_t buffer = 0;
    for (auto const& f : table->frames) buffer = std::max<uint64_t>(buffer, f.compressed_size + f.decompressed_size);
    for (uint64_t output_per_frame : {buffer, buffer / 3, uint64_t{1} << 40}) {
        SCOPED_TRACE(output_per_frame);
        zstd::LoadBudget budget;
        budget.output.assign(frames, output_per_frame);
        std::mutex mutex;
        size_t running = 0;
        std::vector<size_t> handled(frames, 0), running_at_start(frames, 0);
        std::string content(data.size(), '\0');
        std::string const e = zstd::ForEachFrame(tmp / "load.zst", *table, 0, 8, [&](size_t frame, char const* d, size_t n, size_t) {
            {
                std::lock_guard<std::mutex> lock(mutex);
                running_at_start[frame] = ++running;
                handled[frame]++;
            }
            std::memcpy(content.data() + table->frames[frame].decompressed_offset, d, n);
            std::this_thread::sleep_for(std::chrono::milliseconds(3));  // the workers overlap
            std::lock_guard<std::mutex> lock(mutex);
            running--;
            return std::string();
        }, &budget);
        EXPECT_EQ(e, "");
        EXPECT_EQ(content, data);
        EXPECT_EQ(handled, std::vector<size_t>(frames, 1));
        for (size_t i = 0; i < frames; i++) {
            uint64_t const rest = (frames - 1 - i) * output_per_frame;  // the frames after i, started after it
            EXPECT_TRUE(running_at_start[i] == 1 || running_at_start[i] * buffer <= rest)
                << "frame " << i << " started beside " << running_at_start[i] - 1 << " others, " << rest << " bytes to come";
        }
        if (output_per_frame == (uint64_t{1} << 40)) {
            // Only the last frame waits for the others: once it starts, nothing is left to write but its own output.
            EXPECT_LE(budget.frames_on_fewer, 1u);
        } else {
            // With a frame's output of one buffer, 8 workers go on until 8 frames are left; of a third, until 24.
            EXPECT_GE(budget.frames_on_fewer, output_per_frame == buffer ? 4u : 15u);
            EXPECT_LT(budget.frames_on_fewer, frames);
        }
    }
    // With an allowance of four buffers beyond the output (a query run's alignment memory): the frames started beside others
    // stay within what is left to write plus the allowance, and fewer frames start on fewer workers than without it.
    {
        size_t fewer_without = 0;
        for (uint64_t allowance : {uint64_t{0}, 4 * buffer}) {
            SCOPED_TRACE(allowance);
            zstd::LoadBudget budget;
            budget.output.assign(frames, buffer);
            budget.allowance = allowance;
            std::mutex mutex;
            size_t running = 0;
            std::vector<size_t> running_at_start(frames, 0);
            EXPECT_EQ(zstd::ForEachFrame(tmp / "load.zst", *table, 0, 8, [&](size_t frame, char const*, size_t, size_t) {
                {
                    std::lock_guard<std::mutex> lock(mutex);
                    running_at_start[frame] = ++running;
                }
                std::this_thread::sleep_for(std::chrono::milliseconds(3));
                std::lock_guard<std::mutex> lock(mutex);
                running--;
                return std::string();
            }, &budget), "");
            for (size_t i = 0; i < frames; i++) {
                uint64_t const rest = (frames - 1 - i) * buffer + allowance;
                EXPECT_TRUE(running_at_start[i] == 1 || running_at_start[i] * buffer <= rest) << "frame " << i;
            }
            if (allowance == 0) fewer_without = budget.frames_on_fewer;
            else EXPECT_LT(budget.frames_on_fewer, fewer_without);
        }
    }
    // Without a budget: the same frames, each once.
    std::atomic<size_t> count{0};
    EXPECT_EQ(zstd::ForEachFrame(tmp / "load.zst", *table, 0, 8, [&](size_t, char const*, size_t, size_t) {
        count++;
        return std::string();
    }), "");
    EXPECT_EQ(count.load(), frames);
}

TEST(ZstdSeekable, ParallelFramesAreReadInOrderOnAnyThreadCount) {
    ScratchDir tmp;
    // Frames of 64 KB: many, more than the threads may hold ahead.
    std::string const data = WriteSeekable(tmp / "s.zst", {TestData(12345, 1), TestData(size_t{3} << 20, 2), TestData(777, 3)}, size_t{1} << 16);
    std::string error;
    auto const table = zstd::ReadSeekTable(tmp / "s.zst", error);
    ASSERT_TRUE(table) << error;
    for (size_t threads : {1, 2, 4, 9}) {
        for (size_t ahead : {1, 2, 5, 64}) {
            zstd::ParallelFrameStreambuf buffer(tmp / "s.zst", *table, threads, ahead);
            ASSERT_TRUE(buffer.IsOpen());
            std::istream in(&buffer);
            EXPECT_EQ(ReadToEnd(in), data) << threads << " threads, " << ahead << " ahead";
            EXPECT_FALSE(in.bad());
            EXPECT_TRUE(buffer.Error().empty());
        }
    }
    // A reader that stops early ends the threads, also those waiting for room.
    for (size_t threads : {1, 4}) {
        zstd::ParallelFrameStreambuf buffer(tmp / "s.zst", *table, threads, 3);
        std::istream in(&buffer);
        std::string first(100, '\0');
        in.read(first.data(), 100);
        EXPECT_EQ(first, data.substr(0, 100));
    }
    zstd::ParallelFrameStreambuf missing(tmp / "missing.zst", *table, 2, 2);
    EXPECT_FALSE(missing.IsOpen());
}

TEST(ZstdSeekable, ParallelFramesStopAtACorruptFrame) {
    ScratchDir tmp;
    std::string const data = WriteSeekable(tmp / "good.zst", {TestData(size_t{4} << 20)}, size_t{1} << 20);
    std::string corrupt = Slurp(tmp / "good.zst");
    std::string error;
    auto const table = zstd::ReadSeekTable(tmp / "good.zst", error);
    ASSERT_TRUE(table);
    corrupt[table->frames[2].compressed_offset + table->frames[2].compressed_size / 2] ^= 0x5a;  // in frame 3
    Spit(tmp / "corrupt.zst", corrupt);
    for (size_t threads : {1, 3}) {
        testing::internal::CaptureStderr();
        zstd::ParallelFrameStreambuf buffer(tmp / "corrupt.zst", *table, threads, 2);
        std::istream in(&buffer);
        std::string const read = ReadToEnd(in);
        std::string const log = testing::internal::GetCapturedStderr();
        // The frames before it, then the error, printed once.
        EXPECT_EQ(read, data.substr(0, size_t{2} << 20)) << threads << " threads";
        EXPECT_TRUE(in.bad());
        EXPECT_FALSE(buffer.Error().empty());
        EXPECT_EQ(log.find("Error reading"), log.rfind("Error reading")) << log;
        EXPECT_NE(log.find("Error reading"), std::string::npos) << log;
    }
}

TEST(ZstdSeekable, CompressFileWritesAndRecompressesFrames) {
    ScratchDir tmp;
    std::string const data = TestData(size_t{3} << 20);
    Spit(tmp / "reference.fna", data);
    std::string error;
    ASSERT_TRUE(zstd::CompressFile(tmp / "reference.fna", tmp / "reference.fna.zst", {3, 27, 4, 256 << 10}, true, error)) << error;
    auto table = zstd::ReadSeekTable(tmp / "reference.fna.zst", error);
    ASSERT_TRUE(table);
    EXPECT_EQ(table->frames.size(), 12u);

    // A single-frame file is recompressed in place into frames.
    WriteCompressed(tmp / "single.zst", data);
    ASSERT_TRUE(zstd::CompressFile(tmp / "single.zst", tmp / "single.zst", {3, 0, 2, 1 << 20}, true, error)) << error;
    EXPECT_TRUE(zstd::IsSeekable(tmp / "single.zst"));
    zstd::InputFile in(tmp / "single.zst");
    EXPECT_EQ(ReadAll(in, 1 << 16), data);
}

// The index sink scatters a file into key map and values, whatever the frame boundaries.
TEST(ZstdSeekable, IndexSinkScattersIntoKeymapAndValues) {
    ScratchDir tmp;
    std::string const header = TestData(72, 5), keymap = TestData(50000, 6), values = TestData(80000, 7);
    for (uint64_t frame : {uint64_t{1000}, uint64_t{50072}, uint64_t{1} << 20}) {
        WriteSeekable(tmp / "index.zst", {header, keymap, values}, frame);
        std::string km(keymap.size(), '\0'), vals(values.size(), '\0');
        IndexSink sink(header.size(), km.data(), km.size(), vals.data(), vals.size());
        std::string error;
        EXPECT_EQ(zstd::ParallelRead(tmp / "index.zst", 4, sink, error), header.size() + keymap.size() + values.size());
        EXPECT_TRUE(error.empty()) << error;
        EXPECT_EQ(km, keymap) << "frames of " << frame;
        EXPECT_EQ(vals, values) << "frames of " << frame;
    }
}

// The index loader stops with exit 8 on a compressed index that holds too little data, and on
// one whose zstd frame is cut off.
TEST(Zstd, TruncatedCompressedIndexExits8) {
    ScratchDir tmp;
    std::ostringstream header;
    auto put = [&](auto v) { header.write(reinterpret_cast<char const*>(&v), sizeof(v)); };
    constexpr size_t keys = size_t{1} << 30, total = keys + ((keys / 8) + 1) * 4;
    put(Seedmap::kFileMagic);
    put(Seedmap::kIndexVersionMajor);
    put(Seedmap::kIndexVersionMinor);
    put(Seedmap::kIndexVersionPatch);
    for (size_t v : {keys, total, size_t{8}, size_t{8}, size_t{3}, size_t{0}}) put(v);
    WriteCompressed(tmp / "short.prx.zst", header.str() + TestData(size_t{1} << 20));
    EXPECT_EXIT({ Seedmap map; map.Load(tmp / "short.prx.zst"); }, testing::ExitedWithCode(8), "truncated or corrupt");

    std::string const frame = Slurp(tmp / "short.prx.zst");
    Spit(tmp / "cut.prx.zst", frame.substr(0, frame.size() / 2));
    EXPECT_EXIT({ Seedmap map; map.Load(tmp / "cut.prx.zst"); }, testing::ExitedWithCode(8), "truncated or corrupt");
}

// A reference.fna.zst of one frame without a seek table (as the zstd command writes it) is preloaded as one stream, on
// any number of threads: its genes are those of reference.fna. A single gene is not loaded on demand from it.
TEST(Zstd, GenomeLoaderReadsASingleFrameReference) {
    ScratchDir tmp;
    std::vector<std::pair<std::string, std::string>> records = {
            {"1_1", TestData(300, 1)}, {"1_2", TestData(1200, 2)}, {"2_1", TestData(80, 3)}, {"10_3", TestData(5000, 4)}};
    std::string fna, map;
    for (auto const& [name, seq] : records) {
        std::string const header = ">" + name + "\n";
        fna += header;
        map += name.substr(0, name.find('_')) + '\t' + name.substr(name.find('_') + 1) + '\t' + std::to_string(fna.size()) + '\t' +
               std::to_string(fna.size() + seq.size()) + '\n';
        fna += seq + '\n';
    }
    Spit(tmp / "reference.map", map);
    WriteCompressed(tmp / "reference.fna.zst", fna);
    ASSERT_FALSE(zstd::IsSeekable(tmp / "reference.fna.zst"));

    for (int threads : {1, 3}) {
        GenomeLoader compressed(tmp / "reference.fna.zst", tmp / "reference.map");
        EXPECT_TRUE(compressed.IsCompressed());
        compressed.LoadAllGenomes(threads);
        EXPECT_TRUE(compressed.AllGenomesLoaded());
        for (auto const& [name, seq] : records) {
            size_t const taxid = std::stoul(name.substr(0, name.find('_')));
            size_t const gene = std::stoul(name.substr(name.find('_') + 1));
            EXPECT_EQ(compressed.GetGenome(taxid).GetGene(gene).Sequence(), seq) << name << ", " << threads << " threads";
        }
    }

    GenomeLoader lazy(tmp / "reference.fna.zst", tmp / "reference.map");
    EXPECT_EXIT(lazy.GetGenome(1).GetGeneOMP(1), testing::ExitedWithCode(EX_SOFTWARE), "must be preloaded");
}
