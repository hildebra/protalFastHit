// Unit tests for the input stream that inflates in a thread of its own (ThreadedGzStream), the
// read pairs threads take from it (SeqReaderPE), and the stage timers (Benchmark).
#include <gtest/gtest.h>
#include <algorithm>
#include <cerrno>
#include <csignal>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>
#include <omp.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zstd.h>
#include <zstd_errors.h>
#include <zlib-ng.h>
#include "gzstream.h"
#include "IO/ThreadedGzStream.h"
#include "SequenceUtils/FastaBatches.h"
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

namespace {
    // content written as BGZF (bgzf::CompressFile, as the simulator and bgzip write reads).
    std::string Bgzf(ScratchDir const& dir, std::string const& name, std::string const& content) {
        auto const plain = dir.Plain(name + ".plain", content);
        auto const path = (dir.path / name).string();
        EXPECT_EQ(bgzf::CompressFile(plain, path, 3), "");
        return path;
    }

    std::string ReadAllOf(ThreadedGzIstream& is) {
        std::string text;
        char buffer[1 << 16];
        while (is.read(buffer, sizeof(buffer)) || is.gcount() > 0) text.append(buffer, static_cast<size_t>(is.gcount()));
        return text;
    }

    // content as gzip of another writer than protal (zlib-ng, by default at gzip's default level), one member. With
    // `full_header`, its header has every optional field: an extra field, a name, a comment and a header CRC.
    std::string GzipBytes(std::string_view text, bool full_header = false, int level = 6) {
        zng_stream zs {};
        EXPECT_EQ(zng_deflateInit2(&zs, level, Z_DEFLATED, 15 + 16, 8, Z_DEFAULT_STRATEGY), Z_OK);  // 15 + 16: a gzip wrapper
        static unsigned char extra[] = { 'X', 'Y', 4, 0, 1, 2, 3, 4 };
        static unsigned char name[] = "reads.fq", comment[] = "a comment";
        zng_gz_header header {};
        header.extra = extra;
        header.extra_len = sizeof(extra);
        header.name = name;
        header.comment = comment;
        header.hcrc = 1;
        if (full_header) EXPECT_EQ(zng_deflateSetHeader(&zs, &header), Z_OK);
        std::string out(zng_deflateBound(&zs, text.size()) + 64, '\0');
        zs.next_in = reinterpret_cast<uint8_t const*>(text.data());
        zs.avail_in = static_cast<uint32_t>(text.size());
        zs.next_out = reinterpret_cast<uint8_t*>(out.data());
        zs.avail_out = static_cast<uint32_t>(out.size());
        EXPECT_EQ(zng_deflate(&zs, Z_FINISH), Z_STREAM_END);
        out.resize(zs.total_out);
        zng_deflateEnd(&zs);
        return out;
    }

    std::string FileBytes(std::string const& path) {
        std::ifstream is(path, std::ios::binary);
        return std::string((std::istreambuf_iterator<char>(is)), std::istreambuf_iterator<char>());
    }

    // text as one zstd frame with a checksum (as the zstd command writes it). window_log > 0: streamed
    // without a size, so that the frame keeps a window of 2^window_log (zstd --long=N writes such frames).
    std::string ZstdBytes(std::string_view text, int window_log = 0) {
        ZSTD_CCtx* cctx = ZSTD_createCCtx();
        ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, 3);
        ZSTD_CCtx_setParameter(cctx, ZSTD_c_checksumFlag, 1);
        std::string out(ZSTD_compressBound(text.size()) + 1024, '\0');
        size_t written = 0;
        if (window_log == 0) {
            written = ZSTD_compress2(cctx, out.data(), out.size(), text.data(), text.size());
            EXPECT_FALSE(ZSTD_isError(written)) << ZSTD_getErrorName(written);
        } else {
            ZSTD_CCtx_setParameter(cctx, ZSTD_c_windowLog, window_log);
            ZSTD_CCtx_setParameter(cctx, ZSTD_c_enableLongDistanceMatching, 1);
            ZSTD_outBuffer o{ out.data(), out.size(), 0 };
            for (size_t from = 0; from < text.size(); from += 100000) {  // in pieces: the size stays unknown
                ZSTD_inBuffer i{ text.data() + from, std::min<size_t>(100000, text.size() - from), 0 };
                while (i.pos < i.size) EXPECT_FALSE(ZSTD_isError(ZSTD_compressStream2(cctx, &o, &i, ZSTD_e_continue)));
            }
            ZSTD_inBuffer none{ nullptr, 0, 0 };
            while (ZSTD_compressStream2(cctx, &o, &none, ZSTD_e_end) > 0) {}
            written = o.pos;
        }
        ZSTD_freeCCtx(cctx);
        out.resize(written);
        return out;
    }

    // A skippable zstd frame of n bytes (pzstd and the seekable format write such frames between and after others).
    std::string SkippableFrame(size_t n) {
        std::string frame = "\x5a\x2a\x4d\x18";  // magic 0x184D2A5A, little-endian
        for (int i = 0; i < 4; i++) frame += static_cast<char>((n >> (8 * i)) & 0xff);
        return frame + std::string(n, 's');
    }
}

// A zstd file (by its first bytes, whatever its name) is decompressed with libzstd in the stream's own thread:
// one frame, several (cat a.zst b.zst), skippable frames among them, an empty frame, and a frame of a window
// larger than libzstd decodes by default (zstd --long=28 writes it).
TEST(ThreadedGzStream, ReadsAZstdFile) {
    ScratchDir dir;
    auto const content = Fastq(40000, "z");  // ~12 MB, more than all blocks
    std::string several;
    for (size_t from = 0; from < content.size(); from += 1000003) several += ZstdBytes(std::string_view(content).substr(from, 1000003));
    std::string const with_skippable = SkippableFrame(100) + ZstdBytes(content.substr(0, 5000000)) + SkippableFrame(70000) +
                                       ZstdBytes(content.substr(5000000)) + SkippableFrame(16);
    std::string const long_window = ZstdBytes(content, 28);
    {  // libzstd's defaults do not decode it (a window over 2^27)
        ZSTD_DCtx* dctx = ZSTD_createDCtx();
        std::string out(1 << 20, '\0');
        ZSTD_inBuffer i{ long_window.data(), long_window.size(), 0 };
        ZSTD_outBuffer o{ out.data(), out.size(), 0 };
        size_t const ret = ZSTD_decompressStream(dctx, &o, &i);
        ZSTD_freeDCtx(dctx);
        ASSERT_TRUE(ZSTD_isError(ret));
        EXPECT_EQ(ZSTD_getErrorCode(ret), ZSTD_error_frameParameter_windowTooLarge) << ZSTD_getErrorName(ret);
    }
    std::vector<std::pair<std::string, std::string>> const inputs{
        { "one.fq.zst", ZstdBytes(content) }, { "several.fq.zst", several }, { "skippable.fq.zst", with_skippable },
        { "long.fq.zst", long_window }, { "named_as_plain.fq", ZstdBytes(content) } };
    for (auto const& [name, bytes] : inputs) {
        SCOPED_TRACE(name);
        ThreadedGzIstream is(dir.Plain(name, bytes).c_str());
        EXPECT_EQ(ReadAllOf(is), content);
        EXPECT_FALSE(is.rdbuf()->read_failed()) << is.rdbuf()->read_error_message();
    }
    ThreadedGzIstream empty(dir.Plain("empty.fq.zst", ZstdBytes("")).c_str());
    EXPECT_EQ(ReadAllOf(empty), "");
    EXPECT_FALSE(empty.rdbuf()->read_failed()) << empty.rdbuf()->read_error_message();
}

// A zstd file cut inside a frame, with a corrupt frame, or with data after its frames that is no frame, reads as
// a prefix of its content and says why; a cut between frames cannot be told from the end, as with gzip members.
TEST(ThreadedGzStream, ACutOrCorruptZstdFileIsReported) {
    ScratchDir dir;
    auto const first = Fastq(10000, "c"), second = Fastq(10000, "d");
    std::string const one = ZstdBytes(first), two = ZstdBytes(second);
    int files = 0;
    auto read = [&](std::string const& bytes, std::string& error) {
        ThreadedGzIstream is(dir.Plain("reads" + std::to_string(files++) + ".fq.zst", bytes).c_str());
        auto const text = ReadAllOf(is);
        error = is.rdbuf()->read_failed() ? is.rdbuf()->read_error_message() : "";
        return text;
    };
    auto const npos = std::string::npos;
    std::string error;
    auto text = read(one + two.substr(0, two.size() / 2), error);
    EXPECT_NE(error.find("the file ends inside zstd frame 2 (truncated file?)"), npos) << error;
    EXPECT_EQ(text.compare(0, first.size(), first), 0);
    EXPECT_LT(text.size(), first.size() + second.size());
    EXPECT_EQ((first + second).compare(0, text.size(), text), 0);

    text = read(one.substr(0, one.size() - 3), error);  // only the checksum's last bytes missing
    EXPECT_NE(error.find("the file ends inside zstd frame 1 (truncated file?)"), npos) << error;

    std::string corrupt = two;
    for (size_t i = corrupt.size() / 2; i < corrupt.size() / 2 + 16; i++) corrupt[i] = static_cast<char>(corrupt[i] ^ 0x5a);
    text = read(one + corrupt, error);
    EXPECT_NE(error.find("in zstd frame 2 (corrupt file?)"), npos) << error;
    EXPECT_EQ(text.compare(0, first.size(), first), 0);

    EXPECT_EQ(read(one + "not a zstd frame", error), first);
    EXPECT_NE(error.find("in zstd frame 2 (corrupt file?)"), npos) << error;
    EXPECT_EQ(read(one + two, error), first + second);  // whole frames: a file cut between them reads as complete
    EXPECT_EQ(error, "");
}

// gzip that is not BGZF (as sequencers write it) is inflated as a stream (ISA-L), also when protal did
// not write it: here zlib-ng writes it, as one member and as several (all members are read).
TEST(ThreadedGzStream, ReadsGzipOfAnotherWriter) {
    ScratchDir dir;
    auto const content = Fastq(20000, "g");  // ~6 MB, several blocks
    std::string const one = GzipBytes(content);
    std::string several;
    for (size_t from = 0; from < content.size(); from += 1000003) several += GzipBytes(std::string_view(content).substr(from, 1000003));
    for (auto const& [name, bytes] : { std::pair{ "one.fq.gz", one }, std::pair{ "several.fq.gz", several } }) {
        SCOPED_TRACE(name);
        auto const path = dir.Plain(name, bytes);
        ASSERT_FALSE(bgzf::StartsAsBgzf(path));
        ThreadedGzIstream is(path.c_str());
        EXPECT_EQ(ReadAllOf(is), content);
        EXPECT_FALSE(is.rdbuf()->read_failed()) << is.rdbuf()->read_error_message();
    }
}

// A gzip header with every optional field (an extra field, a name, a comment, a header CRC) split between two reads
// of the file, at every place: ISA-L's inflate reads such a header wrongly when it comes in two pieces (ThreadedGzStream
// reads headers with a state of its own). The first member is stored (level 0) and sized so that the second one's
// header begins s bytes before the end of the reader's first read; then members of random sizes with such headers
// through a pipe written a few bytes at a time.
TEST(ThreadedGzStream, AHeaderSplitBetweenReadsIsRead) {
    ScratchDir dir;
    auto const second = Fastq(300, "h");
    std::string const two = GzipBytes(second, true);
    auto const text = Fastq(1000, "s");  // ~330 KB, more than a read
    size_t constexpr kHeader = 10 + (2 + 8) + 9 + 10 + 2;  // fixed part, extra field, name, comment, header CRC
    for (size_t s = 1; s <= kHeader + 2; s++) {
        SCOPED_TRACE(s);
        size_t const target = ThreadedGzStreambuf::kInputBuffer - s;
        std::string first = text.substr(0, target), one;
        for (int i = 0; i < 10 && one.size() != target; i++) {
            one = GzipBytes(first, true, 0);
            first = text.substr(0, first.size() + target - one.size());
        }
        ASSERT_EQ(one.size(), target);
        ThreadedGzIstream is(dir.Plain("split" + std::to_string(s) + ".fq.gz", one + two).c_str());
        EXPECT_EQ(ReadAllOf(is), first + second);
        EXPECT_FALSE(is.rdbuf()->read_failed()) << is.rdbuf()->read_error_message();
    }

    std::signal(SIGPIPE, SIG_IGN);
    std::mt19937 rng(5);
    std::string content, bytes;
    for (int m = 0; m < 300; m++) {
        auto const piece = Fastq(1 + rng() % 40, "p" + std::to_string(m));
        content += piece;
        bytes += GzipBytes(piece, true);
    }
    auto const fifo = (dir.path / "members.fifo").string();
    ASSERT_EQ(::mkfifo(fifo.c_str(), 0600), 0);
    std::thread writer([&fifo, &bytes] {
        std::ofstream os(fifo, std::ios::binary);
        std::mt19937 chunks(9);
        for (size_t at = 0; at < bytes.size();) {
            size_t const n = std::min<size_t>(1 + chunks() % 97, bytes.size() - at);
            os.write(bytes.data() + at, static_cast<std::streamsize>(n));
            os.flush();
            at += n;
        }
    });
    ThreadedGzIstream is(fifo.c_str());
    EXPECT_EQ(ReadAllOf(is), content);
    EXPECT_FALSE(is.rdbuf()->read_failed()) << is.rdbuf()->read_error_message();
    writer.join();
}

TEST(ThreadedGzStream, ReadsABgzfFileBlockByBlock) {
    ScratchDir dir;
    auto const content = Fastq(40000, "b");  // ~12 MB: ~200 BGZF blocks, several output blocks
    auto const path = Bgzf(dir, "reads.fq.gz", content);
    ASSERT_TRUE(bgzf::StartsAsBgzf(path));
    ThreadedGzIstream is(path.c_str());
    EXPECT_EQ(ReadAllOf(is), content);
    EXPECT_FALSE(is.rdbuf()->read_failed()) << is.rdbuf()->read_error_message();

    auto const empty = Bgzf(dir, "empty.fq.gz", "");  // only the end-of-file block
    ThreadedGzIstream none(empty.c_str());
    EXPECT_EQ(ReadAllOf(none), "");
    EXPECT_FALSE(none.rdbuf()->read_failed());
}

TEST(ThreadedGzStream, ACutOrCorruptBgzfFileIsReported) {
    ScratchDir dir;
    auto const content = Fastq(20000, "x");
    // Cut at a block boundary (the end-of-file block gone), inside a block, or with a flipped byte.
    for (int variant = 0; variant < 3; variant++) {
        SCOPED_TRACE(variant);
        auto const path = Bgzf(dir, "reads" + std::to_string(variant) + ".fq.gz", content);
        auto const size = fs::file_size(path);
        if (variant == 0) fs::resize_file(path, size - sizeof(bgzf::kEof));
        if (variant == 1) fs::resize_file(path, size / 2);
        if (variant == 2) {
            std::fstream f(path, std::ios::in | std::ios::out | std::ios::binary);
            f.seekp(static_cast<std::streamoff>(size / 3));
            char c = 0;
            f.read(&c, 1);
            f.seekp(static_cast<std::streamoff>(size / 3));
            c = static_cast<char>(c ^ 0x5a);
            f.write(&c, 1);
        }
        ThreadedGzIstream is(path.c_str());
        auto const read = ReadAllOf(is);
        EXPECT_TRUE(is.rdbuf()->read_failed());
        EXPECT_FALSE(is.rdbuf()->read_error_message().empty());
        EXPECT_LT(read.size(), variant == 0 ? content.size() + 1 : content.size());
        EXPECT_EQ(content.compare(0, read.size(), read), 0);  // what was read is a prefix
        if (variant == 0) EXPECT_NE(is.rdbuf()->read_error_message().find("end-of-file block"), std::string::npos);
    }
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

// Why a file cannot be read, for protal's message: missing, a directory, or not readable.
TEST(ThreadedGzStream, AFileThatCannotBeReadSaysWhy) {
    ScratchDir dir;
    ThreadedGzIstream missing((dir.path / "none.fq.gz").string().c_str());
    EXPECT_EQ(missing.rdbuf()->open_error(), std::strerror(ENOENT));
    ThreadedGzIstream directory(dir.path.string().c_str());
    EXPECT_FALSE(directory.rdbuf()->is_open());
    EXPECT_EQ(directory.rdbuf()->open_error(), "it is a directory");
    if (::geteuid() != 0) {  // root reads any file
        auto const path = dir.Plain("locked.fq", Fastq(1, "l"));
        fs::permissions(path, fs::perms::none);
        ThreadedGzIstream locked(path.c_str());
        EXPECT_FALSE(locked.rdbuf()->is_open());
        EXPECT_EQ(locked.rdbuf()->open_error(), std::strerror(EACCES));
    }
}

// A pipe (a FIFO, process substitution) is read from its one descriptor: plain, gzip, BGZF (which a
// pipe cannot be peeked at for, so it is read as gzip members), or zstd.
TEST(ThreadedGzStream, ReadsAPipe) {
    ScratchDir dir;
    std::signal(SIGPIPE, SIG_IGN);  // a reader that stops early fails the test, not the process
    auto const content = Fastq(20000, "f");  // ~6 MB, more than a pipe holds
    std::vector<std::pair<std::string, std::string>> const inputs{
        { "plain", content }, { "gzip", GzipBytes(content) }, { "bgzf", FileBytes(Bgzf(dir, "reads.fq.gz", content)) },
        { "zstd", ZstdBytes(content) } };
    for (auto const& [name, bytes] : inputs) {
        SCOPED_TRACE(name);
        auto const fifo = (dir.path / (name + ".fifo")).string();
        ASSERT_EQ(::mkfifo(fifo.c_str(), 0600), 0);
        std::thread writer([&fifo, &bytes] { std::ofstream(fifo, std::ios::binary) << bytes; });
        ThreadedGzIstream is(fifo.c_str());
        EXPECT_EQ(ReadAllOf(is), content);
        EXPECT_FALSE(is.rdbuf()->read_failed()) << is.rdbuf()->read_error_message();
        writer.join();
    }
}

// After a gzip member comes another member, zero bytes of padding, or the end of the file.
// Anything else is a damaged member header, reported rather than taken for the end of the file
// (zlib's gzread ignores it); a corrupt or cut member is reported by its number.
TEST(ThreadedGzStream, WhatFollowsAGzipMemberMustBeAMember) {
    ScratchDir dir;
    auto const first = Fastq(3000, "m"), second = Fastq(3000, "n");
    std::string const one = GzipBytes(first), two = GzipBytes(second);
    int files = 0;
    auto read = [&](std::string const& bytes, std::string& error) {
        ThreadedGzIstream is(dir.Plain("reads" + std::to_string(files++) + ".fq.gz", bytes).c_str());
        auto const text = ReadAllOf(is);
        error = is.rdbuf()->read_failed() ? is.rdbuf()->read_error_message() : "";
        return text;
    };
    auto const npos = std::string::npos;
    std::string error;
    EXPECT_EQ(read(one + two + std::string(70000, '\0'), error), first + second);
    EXPECT_EQ(error, "");

    std::string damaged = two;
    damaged[1] = 'X';
    EXPECT_EQ(read(one + damaged, error), first);
    EXPECT_NE(error.find("the data at byte " + std::to_string(one.size()) + ", after gzip member 1, is no gzip member"), npos) << error;
    EXPECT_EQ(read(one + std::string(200000, '\0') + "x", error), first);  // padding that is not all zero
    EXPECT_NE(error.find("after gzip member 1, is no gzip member"), npos) << error;

    std::string corrupt = two;
    for (size_t i = corrupt.size() / 2; i < corrupt.size() / 2 + 16; i++) corrupt[i] = static_cast<char>(corrupt[i] ^ 0x5a);
    auto const text = read(one + corrupt, error);
    EXPECT_NE(error.find("in gzip member 2 (corrupt file?)"), npos) << error;
    EXPECT_EQ(text.compare(0, first.size(), first), 0);

    EXPECT_EQ(read(one + two.substr(0, two.size() - 4), error).substr(0, first.size()), first);  // no length field
    EXPECT_NE(error.find("the file ends inside gzip member 2 (truncated file?)"), npos) << error;

    // A BGZF file goes on with members of other kinds (cat a.bgzf.gz b.gz), which are held to the same.
    auto const bgzf = FileBytes(Bgzf(dir, "first.fq.gz", first));
    EXPECT_EQ(read(bgzf + two + one, error), first + second + first);
    EXPECT_EQ(error, "");
    EXPECT_EQ(read(bgzf + two + damaged, error), first + second);
    EXPECT_NE(error.find("the data at byte " + std::to_string(bgzf.size() + two.size())), npos) << error;
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

// The single-end readers, in batches of records (SeqReaderSE) and in blocks of bytes (SeqReader),
// share the reader lock as the pairs do: each read once, on more threads than cores too.
TEST(ThreadedGzStream, ThreadsTakeEverySingleEndReadOnce) {
    ScratchDir dir;
    size_t const n = 30000;
    auto const path = dir.Gzip("se.fq.gz", Fastq(n, "se"));
    for (bool const batches : { true, false }) {
        SCOPED_TRACE(batches ? "SeqReaderSE" : "SeqReader");
        ThreadedGzIstream is(path.c_str());
        SeqReaderSE batch_global{ is };
        SeqReader block_global{ is };
        std::vector<std::string> ids;
#pragma omp parallel num_threads(8) shared(batch_global, block_global, ids, batches) default(none)
        {
            SeqReaderSE batch_reader{ batch_global };
            SeqReader block_reader{ block_global };
            FastxRecord record;
            std::vector<std::string> mine;
            while (batches ? batch_reader(record) : block_reader(record)) mine.push_back(record.id);
#pragma omp critical(test_ids)
            {
                ids.insert(ids.end(), mine.begin(), mine.end());
                batch_global.UpdateSuccess(batch_reader);
            }
        }
        EXPECT_TRUE(batch_global.Success());
        EXPECT_EQ(ids.size(), n);
        EXPECT_EQ(std::set<std::string>(ids.begin(), ids.end()).size(), n);
        EXPECT_FALSE(is.rdbuf()->read_failed());
    }
}

TEST(ThreadedGzStream, TakeLinesCutsWholeLinesAcrossBlocks) {
    ScratchDir dir;
    // Short, empty and CRLF lines, one longer than a block, and a last line without its '\n'.
    std::string content;
    for (int i = 0; i < 30000; i++) content += (i % 5 == 0 ? std::string() : std::string(i % 97, 'a' + i % 26)) + (i % 7 ? "\n" : "\r\n");
    content += std::string(ThreadedGzStreambuf::kBlockSize + 12345, 'L') + "\n" + Fastq(20000, "t") + "last";
    std::stringstream reference(content);
    auto const lines = Lines(reference);
    for (auto const& path : { dir.Gzip("lines.gz", content), dir.Plain("lines.txt", content) }) {
        SCOPED_TRACE(path);
        ThreadedGzIstream is(path.c_str());
        std::string taken;
        size_t count = 0, n = 0;
        while ((n = is.rdbuf()->TakeLines(7, taken)) > 0) {
            count += n;
            ASSERT_EQ(taken.back(), '\n');
            if (n < 7) break;
        }
        EXPECT_EQ(is.rdbuf()->TakeLines(7, taken), 0u);
        EXPECT_EQ(count, lines.size());
        EXPECT_EQ(taken, content + "\n");
    }
}

namespace {
    std::vector<std::string> Records(std::function<bool(FastxRecord&)> const& next) {
        std::vector<std::string> records;
        FastxRecord record;
        while (next(record)) records.push_back(record.header + "|" + record.id + "|" + record.sequence + "|" + record.quality);
        return records;
    }
}

// SeqReaderSE parses FASTQ batches outside the reader lock (BufferedFastxReader::NextFastq), taken
// from protal's input stream with TakeLines or from any other stream line by line; SeqReader still
// parses them as before (ReadNextSequence).
TEST(SeqReader, FastqRecordsAreTheSameFromAnyStream) {
    ScratchDir dir;
    std::string content = "@r1 some description\nACGT\n+\nIIII\n@r2\tx\r\nACGTT\r\n+r2\r\nIIIII\r\n";
    content += Fastq(100, "b");  // several batches of 32
    std::string const long_read(ThreadedGzStreambuf::kBlockSize * 3 / 2, 'C');
    content += "@long\n" + long_read + "\n+\n" + std::string(long_read.size(), '5') + "\n@end/1\nAC\n+\nII";

    std::istringstream oracle_stream(content);
    SeqReader oracle(oracle_stream);
    auto const expected = Records([&](FastxRecord& r) { return oracle(r); });
    ASSERT_EQ(expected.size(), 104u);
    EXPECT_EQ(expected[0], "@r1 some description|r1|ACGT|IIII");
    EXPECT_EQ(expected[1], "@r2\tx|r2|ACGTT|IIIII");
    EXPECT_EQ(expected.back(), "@end/1|end/1|AC|II");

    std::istringstream plain_stream(content);
    ThreadedGzIstream plain(dir.Plain("reads.fq", content).c_str()), gz(dir.Gzip("reads.fq.gz", content).c_str());
    ThreadedGzIstream zst(dir.Plain("reads.fq.zst", ZstdBytes(content)).c_str());
    for (std::istream* is : std::initializer_list<std::istream*>{ &plain_stream, &plain, &gz, &zst }) {
        SeqReaderSE reader(*is);
        EXPECT_EQ(Records([&](FastxRecord& r) { return reader(r); }), expected);
        EXPECT_TRUE(reader.Success());
    }
}

namespace {
    using FastaRecords = std::vector<std::pair<std::string, std::string>>;  // header, sequence

    FastaRecords SeqReaderRecords(std::string const& text) {
        std::istringstream is(text);
        SeqReader reader(is);
        FastxRecord record;
        FastaRecords records;
        while (reader(record)) records.emplace_back(record.header, record.sequence);
        return records;
    }

    FastaRecords BatchRecords(std::string const& text, size_t bytes, size_t& batches) {
        std::istringstream is(text);
        FastaBatches reader(is);
        std::string batch, scratch;
        FastaRecords records;
        batches = 0;
        while (reader.Next(batch, bytes)) {
            batches++;
            EXPECT_EQ(batch.front(), '>') << "batch " << batches;
            ForEachFastaRecord(batch, scratch, [&](std::string_view header, std::string_view sequence) {
                records.emplace_back(header, sequence);
            });
        }
        EXPECT_EQ(reader.Error(), "");
        return records;
    }
}

// The index build's parallel passes read the reference in batches of whole records (FastaBatches)
// and parse each in a thread (ForEachFastaRecord): the records SeqReader gives, for any batch size.
TEST(FastaBatches, GiveTheRecordsOfSeqReaderInBatchesOfAnySize) {
    std::mt19937 rng(3);
    std::string text;
    for (int i = 0; i < 300; i++) {
        size_t const length = i == 100 ? 300000 : 50 + rng() % 2000;  // one record longer than most batches
        std::string sequence;
        for (size_t j = 0; j < length; j++) sequence += "ACGT"[rng() % 4];
        std::string const eol = i % 7 == 0 ? "\r\n" : "\n";
        std::string const blanks = i % 11 == 0 ? " \t" : "";
        text += ">" + std::to_string(i % 17 + 1) + "_" + std::to_string(i + 1) + (i % 5 == 0 ? " a description" : "") + eol;
        if (i % 3 == 0) {  // several lines of 60 bases
            for (size_t j = 0; j < length; j += 60) text += sequence.substr(j, 60) + blanks + eol;
        } else {
            text += sequence + blanks + eol;
        }
    }
    text.pop_back();  // the last record ends without a line end
    auto const expected = SeqReaderRecords(text);
    ASSERT_EQ(expected.size(), 300u);
    for (size_t bytes : { size_t{1}, size_t{100}, size_t{4096}, size_t{65536}, size_t{1} << 20, size_t{1} << 26 }) {
        size_t batches = 0;
        EXPECT_EQ(BatchRecords(text, bytes, batches), expected) << bytes << " bytes per batch";
        // A batch ends before the first record that starts `bytes` or more into it.
        EXPECT_LE(batches, expected.size());
        if (bytes == 1) EXPECT_EQ(batches, expected.size()) << "one record per batch";
        // At most `bytes` plus one record (2,100 bytes here, but for the long one) each
        if (bytes <= 4096) EXPECT_GE(batches, (text.size() - 300000) / (bytes + 2100)) << bytes << " bytes per batch";
        if (bytes >= size_t{1} << 20) EXPECT_EQ(batches, 1u);
    }
}

TEST(FastaBatches, AFileThatIsNotFastaIsAnError) {
    std::string batch;
    std::istringstream fastq("@r1\nACGT\n+\nIIII\n");
    FastaBatches reader(fastq);
    EXPECT_FALSE(reader.Next(batch, 1024));
    EXPECT_NE(reader.Error().find("starts with '>'"), std::string::npos) << reader.Error();
    std::istringstream empty("");
    FastaBatches none(empty);
    EXPECT_FALSE(none.Next(batch, 1024));
    EXPECT_EQ(none.Error(), "");
}

TEST(SeqReader, AMalformedFastqRecordIsAnError) {
    ScratchDir dir;
    std::string const content = "@r1\nACGT\n+\nIIII\nr2\nACGT\n+\nIIII\n";
    std::istringstream plain_stream(content);
    ThreadedGzIstream gz(dir.Gzip("bad.fq.gz", content).c_str());
    for (std::istream* is : std::initializer_list<std::istream*>{ &plain_stream, &gz }) {
        SeqReaderSE reader(*is);
        FastxRecord record;
        testing::internal::CaptureStderr();
        ASSERT_TRUE(reader(record));
        EXPECT_EQ(record.id, "r1");
        EXPECT_FALSE(reader(record));
        EXPECT_NE(testing::internal::GetCapturedStderr().find("malformed FASTQ file (exp. '@', saw r2)"), std::string::npos);
        EXPECT_FALSE(reader.Success());
    }
}

// A record cut short (a truncated file), without its '+' line, with more or fewer qualities than
// bases, or without a name is an error after the records before it, in SeqReaderSE's batches and
// SeqReader's lines alike.
TEST(SeqReader, AnIncompleteOrUnevenFastqRecordIsAnError) {
    std::string const good = "@r1\nACGT\n+\nIIII\n";
    std::vector<std::pair<std::string, std::string>> const cases{
        { good + "@r2\nACGT\n", "malformed FASTQ file: read r2 is incomplete (truncated file?)" },
        { good + "@r2\nACGT\n+\n", "malformed FASTQ file: read r2 is incomplete (truncated file?)" },
        { good + "@r2\nACGT\nIIII\n@r3\nACGT\n+\nIIII\n", "malformed FASTQ file: read r2 has no '+' line after its sequence" },
        { good + "@r2\nACGT\n+\nIII\n", "malformed FASTQ file: read r2 has 4 bases but 3 qualities" },
        { good + "@r2\nACGT\n+\nII", "malformed FASTQ file: read r2 has 4 bases but 2 qualities" },
        { good + "@\nACGT\n+\nIIII\n", "a read without a name" },
    };
    for (auto const& [content, message] : cases) {
        SCOPED_TRACE(content);
        std::istringstream batch_stream(content), line_stream(content);
        SeqReaderSE batches(batch_stream);
        SeqReader lines(line_stream);
        for (auto const& [next, success] : std::initializer_list<std::pair<std::function<bool(FastxRecord&)>, std::function<bool()>>>{
                 { [&](FastxRecord& r) { return batches(r); }, [&] { return batches.Success(); } },
                 { [&](FastxRecord& r) { return lines(r); }, [&] { return lines.Success(); } } }) {
            FastxRecord record;
            testing::internal::CaptureStderr();
            ASSERT_TRUE(next(record));
            EXPECT_EQ(record.id, "r1");
            EXPECT_FALSE(next(record));
            EXPECT_NE(testing::internal::GetCapturedStderr().find(message), std::string::npos);
            EXPECT_FALSE(success());
        }
    }
}

// Lowercase bases (soft-masked reads) are read as uppercase, from FASTQ and FASTA alike.
TEST(SeqReader, BasesAreReadInUppercase) {
    for (std::string const content : { "@r1\nacgtn\n+\nIIIII\n@r2\nAcGt\n+\nIIII\n", ">r1\nacgtn\n>r2\nAc\ngt\n" }) {
        SCOPED_TRACE(content);
        std::istringstream batch_stream(content), line_stream(content);
        SeqReaderSE batches(batch_stream);
        SeqReader lines(line_stream);
        for (auto const& next : std::initializer_list<std::function<bool(FastxRecord&)>>{
                 [&](FastxRecord& r) { return batches(r); }, [&](FastxRecord& r) { return lines(r); } }) {
            std::vector<std::string> sequences;
            FastxRecord record;
            while (next(record)) sequences.push_back(record.sequence);
            EXPECT_EQ(sequences, (std::vector<std::string>{ "ACGTN", "ACGT" }));
        }
    }
}

// With a maximum length (short reads: MAX_SHORT_READ_LENGTH), the first longer read stops the
// reader and fails it, also after the first 100 reads protal checks before it starts.
TEST(SeqReader, AReadOverTheMaximumLengthStopsTheReader) {
    std::string const content = Fastq(150, "s") + "@long\n" + std::string(1001, 'A') + "\n+\n" + std::string(1001, 'I') + "\n" + Fastq(40, "t");
    std::istringstream is(content);
    SeqReaderSE reader(is, 0, 1000);
    EXPECT_EQ(Records([&](FastxRecord& r) { return reader(r); }).size(), 150u);
    EXPECT_FALSE(reader.Success());
    EXPECT_EQ(reader.TooLong(), 1001u);
    EXPECT_EQ(reader.TooLongId(), "long");

    std::istringstream unused("");
    SeqReaderSE joined(unused);  // the reader the threads' copies are joined into
    joined.UpdateSuccess(SeqReaderSE(joined));
    EXPECT_TRUE(joined.Success());
    joined.UpdateSuccess(reader);
    EXPECT_FALSE(joined.Success());
    EXPECT_EQ(joined.TooLong(), 1001u);
    EXPECT_EQ(joined.TooLongId(), "long");

    std::istringstream again(content);
    SeqReaderSE unlimited(again);
    EXPECT_EQ(Records([&](FastxRecord& r) { return unlimited(r); }).size(), 191u);
    EXPECT_TRUE(unlimited.Success());
}

// Mates are named alike when their names are the same, or the same but for a last 1 and 2.
TEST(SeqReaderPE, MatesAreNamedAlikeButForTheirNumber) {
    EXPECT_TRUE(SeqReaderPE::MatesNamedAlike("r7", "r7"));
    EXPECT_TRUE(SeqReaderPE::MatesNamedAlike("r7/1", "r7/2"));
    EXPECT_TRUE(SeqReaderPE::MatesNamedAlike("r7.1", "r7.2"));
    EXPECT_FALSE(SeqReaderPE::MatesNamedAlike("r7/2", "r7/1"));
    EXPECT_FALSE(SeqReaderPE::MatesNamedAlike("r7/1", "r8/2"));
    EXPECT_FALSE(SeqReaderPE::MatesNamedAlike("r7", "r71"));
    EXPECT_FALSE(SeqReaderPE::MatesNamedAlike("", "2"));

    // Pairs whose mates' names differ are counted (for protal's warning), and the first is kept.
    std::string const r1 = Fastq(50, "a"), r2 = Fastq(10, "a") + Fastq(40, "b");
    std::istringstream is1(r1), is2(r2);
    SeqReaderPE reader(is1, is2);
    FastxRecord a, b;
    size_t pairs = 0;
    while (reader(a, b)) pairs++;
    EXPECT_EQ(pairs, 50u);
    EXPECT_TRUE(reader.Success());
    EXPECT_EQ(reader.NameMismatches(), 40u);
    EXPECT_EQ(reader.FirstNameMismatch(), "a10/1, b0/1");
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

namespace {
    // Busy-waits, so that an interval lasts at least as long as asked whatever the scheduler does.
    void Spin(std::chrono::microseconds length) {
        auto const end = std::chrono::steady_clock::now() + length;
        while (std::chrono::steady_clock::now() < end) {}
    }
}

// The tests below check the sampling by which calls are timed, which needs no clock, and the estimate
// only from below: a busy machine can make a timed interval longer, never shorter than its spin.
TEST(Benchmark, ASampledStageTimesTheFirstCallAndThenEveryNth) {
    Benchmark bm{"per read", 0, Benchmark::kPerRead};
    std::vector<int> timed;
    for (int i = 0; i < 400; i++) {
        bm.Start();
        if (bm.Timing()) timed.push_back(i);
        bm.Stop();
    }
    EXPECT_EQ(timed, (std::vector<int>{ 0, 61, 122, 183, 244, 305, 366 }));
    EXPECT_EQ(bm.TimedCalls(), timed.size());
    EXPECT_FALSE(bm.Timing());
}

TEST(Benchmark, ASampledStageWithAnOddPeriodTimesBothMatesOfAPair) {
    // A stage called once per mate: an even period would time the same mate every time.
    Benchmark bm{"per mate", 0, Benchmark::kPerRead};
    size_t first_mate = 0, second_mate = 0;
    for (int i = 0; i < 6100; i++) {
        bm.Start();
        if (bm.Timing()) (i % 2 ? second_mate : first_mate)++;
        bm.Stop();
    }
    EXPECT_EQ(first_mate + second_mate, 100u);
    EXPECT_NEAR(static_cast<double>(first_mate), 50.0, 1.0);
}

TEST(Benchmark, ASampledStageEstimatesTheTimeOfAllItsCalls) {
    Benchmark bm{"per read", 0, Benchmark::kPerRead};
    constexpr int kCalls = 6100;  // 100 timed calls
    for (int i = 0; i < kCalls; i++) {
        bm.Start();
        Spin(std::chrono::microseconds(20));
        bm.Stop();
    }
    EXPECT_GE(bm.Seconds(), 0.99 * kCalls * 20e-6);
    EXPECT_EQ(bm.Threads(), 1u);
}

TEST(Benchmark, ASampledStageScalesTheFirstCallToTheCallsMade) {
    Benchmark bm{"per read", 0, Benchmark::kPerRead};
    bm.Start();
    Spin(std::chrono::microseconds(2000));
    bm.Stop();
    EXPECT_GE(bm.Seconds(), 0.002);  // a short run has a value from its first call
    for (int i = 1; i < 61; i++) {   // 60 more calls, none timed: nothing is added, the mean times 61 is
        bm.Start();
        EXPECT_FALSE(bm.Timing());
        bm.Stop();
    }
    EXPECT_EQ(bm.TimedCalls(), 1u);
    EXPECT_GE(bm.Seconds(), 61 * 0.002);
}

TEST(Benchmark, ASampledStageJoinsAndPrintsItsEstimate) {
    Benchmark global{"stage", 0};
    for (int t = 0; t < 3; t++) {
        Benchmark local{"stage", 0, Benchmark::kPerRead};
        for (int i = 0; i < 122; i++) {  // 2 timed calls
            local.Start();
            Spin(std::chrono::microseconds(500));
            local.Stop();
        }
        global.Join(local);
    }
    EXPECT_EQ(global.Threads(), 3u);
    EXPECT_GE(global.Seconds(), 3 * 122 * 500e-6 * 0.99);
    testing::internal::CaptureStdout();
    global.PrintResults();
    auto const printed = testing::internal::GetCapturedStdout();
    EXPECT_NE(printed.find("stage took"), std::string::npos) << printed;
    EXPECT_NE(printed.find("mean over 3 threads"), std::string::npos) << printed;
}

TEST(Benchmark, AStartedStageIsStoppedWhenPrinted) {
    Benchmark bm{"open"};
    bm.Start();
    Spin(std::chrono::microseconds(3000));
    testing::internal::CaptureStdout();
    bm.PrintResults();
    testing::internal::GetCapturedStdout();
    EXPECT_GE(bm.Seconds(), 0.003);
    double const once = bm.Seconds();
    testing::internal::CaptureStdout();
    bm.PrintResults();
    testing::internal::GetCapturedStdout();
    EXPECT_DOUBLE_EQ(bm.Seconds(), once);  // not stopped, and not added, a second time
}
