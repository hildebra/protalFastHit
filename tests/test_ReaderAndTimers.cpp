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
#include "ReadType.h"
#include "SequenceUtils/FastaBatches.h"
#include "SequenceUtils/SeqReader.h"
#include "Utilities/Benchmark.h"
#include "TestUtil.h"

using namespace protal;
namespace fs = std::filesystem;
using protal::test::ScratchDir;

namespace {
    // content written gzip-compressed (gzstream) into dir; returns the file's path.
    std::string Gzip(ScratchDir const& dir, std::string const& name, std::string const& content) {
        auto const file = dir / name;
        ogzstream os(file.c_str());
        os << content;
        os.close();
        return file;
    }

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

namespace {
    // content written as BGZF (bgzf::Writer, as the simulator writes reads; bgzip writes the same format).
    std::string Bgzf(ScratchDir const& dir, std::string const& name, std::string const& content) {
        auto const path = (dir.path / name).string();
        bgzf::Writer writer(path);
        writer.Write(content.data(), content.size());
        EXPECT_TRUE(writer.Close()) << writer.Error();
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
        ThreadedGzIstream is(dir.Write(name, bytes).c_str());
        EXPECT_EQ(ReadAllOf(is), content);
        EXPECT_FALSE(is.rdbuf()->read_failed()) << is.rdbuf()->read_error_message();
    }
    ThreadedGzIstream empty(dir.Write("empty.fq.zst", ZstdBytes("")).c_str());
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
        ThreadedGzIstream is(dir.Write("reads" + std::to_string(files++) + ".fq.zst", bytes).c_str());
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
        auto const path = dir.Write(name, bytes);
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
        ThreadedGzIstream is(dir.Write("split" + std::to_string(s) + ".fq.gz", one + two).c_str());
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
    {
        ThreadedGzIstream is(fifo.c_str());
        EXPECT_EQ(ReadAllOf(is), content);
        EXPECT_FALSE(is.rdbuf()->read_failed()) << is.rdbuf()->read_error_message();
    }  // closed before the join: a writer left blocked on a full pipe by a reader that stopped early fails, not hangs
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
    auto const path = Gzip(dir, "reads.fq.gz", content);
    fs::resize_file(path, fs::file_size(path) / 2);

    ThreadedGzIstream is(path.c_str());
    std::string read((std::istreambuf_iterator<char>(is)), std::istreambuf_iterator<char>());
    EXPECT_TRUE(is.rdbuf()->read_failed());
    EXPECT_FALSE(is.rdbuf()->read_error_message().empty());
    EXPECT_GT(read.size(), 0u);
    EXPECT_LT(read.size(), content.size());
    EXPECT_EQ(content.compare(0, read.size(), read), 0);
}

// A file that cannot be read fails to open, reads nothing (not as a read error) and says why, for protal's message:
// missing, a directory, or not readable.
TEST(ThreadedGzStream, AFileThatCannotBeReadSaysWhy) {
    ScratchDir dir;
    ThreadedGzIstream missing((dir.path / "none.fq.gz").string().c_str());
    EXPECT_FALSE(missing.good());
    EXPECT_FALSE(missing.rdbuf()->is_open());
    std::string line;
    EXPECT_FALSE(std::getline(missing, line));
    EXPECT_FALSE(missing.rdbuf()->read_failed());
    EXPECT_EQ(missing.rdbuf()->open_error(), std::strerror(ENOENT));
    ThreadedGzIstream directory(dir.path.string().c_str());
    EXPECT_FALSE(directory.rdbuf()->is_open());
    EXPECT_EQ(directory.rdbuf()->open_error(), "it is a directory");
    if (::geteuid() != 0) {  // root reads any file
        auto const path = dir.Write("locked.fq", Fastq(1, "l"));
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
        {
            ThreadedGzIstream is(fifo.c_str());
            EXPECT_EQ(ReadAllOf(is), content);
            EXPECT_FALSE(is.rdbuf()->read_failed()) << is.rdbuf()->read_error_message();
        }  // closed before the join (see AHeaderSplitBetweenReadsIsRead)
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
        ThreadedGzIstream is(dir.Write("reads" + std::to_string(files++) + ".fq.gz", bytes).c_str());
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

namespace {
    // Sets the inflating threads for gzip files opened in its scope, BGZF (protal's) and not (ParallelGzip's); back to
    // one after.
    struct InflateThreads {
        explicit InflateThreads(unsigned n) {
            ThreadedGzStreambuf::SetInflateThreads(n);
            ThreadedGzStreambuf::SetGzipThreads(n);
        }
        ~InflateThreads() {
            ThreadedGzStreambuf::SetInflateThreads(1);
            ThreadedGzStreambuf::SetGzipThreads(1);
        }
    };

    struct ReadResult {
        std::string text, error;
        bool failed = false;
        unsigned threads = 0;
        bool parallel = false;  // a gzip file inflated in pieces on several threads
        bool members = false;   // a gzip file of several members inflated by members on several threads
        bool switched = false;  // read in pieces, it went on by members
        size_t restarts = 0;    // how often it went back to one thread for the rest of the file
    };

    ReadResult ReadWithThreads(std::string const& path, unsigned threads) {
        InflateThreads const setting(threads);
        ThreadedGzIstream is(path.c_str());
        ReadResult r;
        r.threads = is.rdbuf()->inflate_threads();
        r.parallel = is.rdbuf()->parallel_gzip();
        r.members = is.rdbuf()->opened_as_gzip_members();
        r.text = ReadAllOf(is);
        r.failed = is.rdbuf()->read_failed();
        r.error = is.rdbuf()->read_error_message();
        r.restarts = is.rdbuf()->gzip_member_restarts();
        r.switched = is.rdbuf()->switched_to_gzip_members();
        return r;
    }

    // The longest gzip member read whole for the gzip files opened in its scope (0: none read by members, so that a file
    // of several members is inflated in pieces); the default after.
    struct MaxGzipMember {
        explicit MaxGzipMember(size_t bytes) { ThreadedGzStreambuf::SetMaxGzipMemberBytes(bytes); }
        ~MaxGzipMember() { ThreadedGzStreambuf::SetMaxGzipMemberBytes(ThreadedGzStreambuf::kMaxGzipMember); }
    };
}

// BGZF inflated on several threads gives the bytes one thread gives: files of no, one and many blocks, BGZF files
// concatenated (end-of-file blocks inside), and BGZF followed by gzip members of another kind.
TEST(ThreadedGzStream, BgzfOnSeveralThreadsReadsAsOnOne) {
    ScratchDir dir;
    auto const big = Fastq(40000, "p");  // ~12 MB: ~200 blocks, a dozen output blocks
    auto const small = Fastq(30, "s");
    std::string many_text, many_bytes;
    for (int i = 0; i < 60; i++) {  // 60 files of 1-3 blocks each, one after the other
        auto const part = Fastq(50 + 97 * static_cast<size_t>(i), "c" + std::to_string(i) + "_");
        many_text += part;
        many_bytes += FileBytes(Bgzf(dir, "part" + std::to_string(i) + ".fq.gz", part));
    }
    std::vector<std::pair<std::string, std::string>> files = {  // (path, content)
        { Bgzf(dir, "empty.fq.gz", ""), "" },
        { Bgzf(dir, "small.fq.gz", small), small },
        { Bgzf(dir, "big.fq.gz", big), big },
        { dir.Write("many.fq.gz", many_bytes), many_text },
        { dir.Write("mixed.fq.gz", FileBytes(Bgzf(dir, "first.fq.gz", big)) + GzipBytes(small) + GzipBytes(big)), big + small + big },
    };
    for (auto const& [path, content] : files) {
        for (unsigned threads : { 1u, 2u, 3u, 4u, 7u }) {
            SCOPED_TRACE(path + ", " + std::to_string(threads) + " threads");
            auto const r = ReadWithThreads(path, threads);
            EXPECT_EQ(r.threads, threads);
            EXPECT_EQ(r.text.size(), content.size());
            EXPECT_TRUE(r.text == content);
            EXPECT_FALSE(r.failed) << r.error;
        }
    }
    // A file that is not compressed is read on its own thread whatever the setting; a zstd file's frames are
    // decompressed on as many as a BGZF file's blocks (gzip of another kind: in pieces, below).
    for (auto const& [path, threads] : { std::pair{ dir.Write("plain.fq", small), 1u },
                                         std::pair{ dir.Write("small.fq.zst", ZstdBytes(small)), 4u } }) {
        EXPECT_EQ(ReadWithThreads(path, 4).threads, threads);
        EXPECT_EQ(ReadWithThreads(path, 4).text, small);
    }
}

namespace {
    // A FASTQ file of n reads of random bases and qualities: it compresses about 3x (Fastq's repeats compress
    // ~100x), so that a gzip file of a few MB of it spans several of ParallelGzip's pieces.
    std::string RandomFastq(size_t n, std::string const& prefix, unsigned seed) {
        std::mt19937 rng(seed);
        std::string text;
        for (size_t i = 0; i < n; i++) {
            size_t const length = 100 + rng() % 101;
            text += '@' + prefix + std::to_string(i) + "/1\n";
            for (size_t j = 0; j < length; j++) text += "ACGT"[rng() % 4];
            text += "\n+\n";
            for (size_t j = 0; j < length; j++) text += static_cast<char>('#' + rng() % 40);
            text += '\n';
        }
        return text;
    }
}

// A gzip file that is not BGZF, inflated on several threads, gives the bytes ISA-L gives on one: in pieces
// (ParallelGzip) one member of several pieces (at gzip's default level, the fastest, stored blocks: none found, read in
// order), a header with every optional field, a tiny file; by members, as its second member begins near its start,
// members of several sizes one after the other, gzip followed by BGZF blocks, zero padding after the last member, an
// empty last member.
TEST(ThreadedGzStream, GzipOnSeveralThreadsReadsAsOnOne) {
    ScratchDir dir;
    auto const big = RandomFastq(40000, "r", 1);  // ~12 MB, ~4.5 MB as gzip: ~9 pieces of 512 KB
    auto const small = Fastq(30, "s");
    std::string several_text, several_bytes;
    std::mt19937 rng(3);
    for (int m = 0; m < 25; m++) {
        auto const part = RandomFastq(1 + rng() % 3000, "m" + std::to_string(m) + "_", static_cast<unsigned>(m));
        several_text += part;
        several_bytes += GzipBytes(part);
    }
    std::vector<std::pair<std::string, std::string>> const files = {  // (path, content)
        { dir.Write("big.fq.gz", GzipBytes(big)), big },
        { dir.Write("fast.fq.gz", GzipBytes(big, false, 1)), big },
        { dir.Write("stored.fq.gz", GzipBytes(big.substr(0, 3000000), false, 0)), big.substr(0, 3000000) },
        { dir.Write("header.fq.gz", GzipBytes(big, true)), big },
        { dir.Write("several.fq.gz", several_bytes), several_text },
        { dir.Write("then_bgzf.fq.gz", GzipBytes(small) + FileBytes(Bgzf(dir, "b.fq.gz", big))), small + big },
        { dir.Write("tiny.fq.gz", GzipBytes(small)), small },
        { dir.Write("padded.fq.gz", GzipBytes(big) + GzipBytes(small) + std::string(70000, '\0')), big + small },
        { dir.Write("empty.fq.gz", GzipBytes("")), "" },
        { dir.Write("empty_last.fq.gz", GzipBytes(big) + GzipBytes("")), big },
    };
    std::set<std::string> const by_members = { "several.fq.gz", "then_bgzf.fq.gz", "padded.fq.gz", "empty_last.fq.gz" };
    for (auto const& [path, content] : files) {
        ASSERT_FALSE(bgzf::StartsAsBgzf(path));
        auto const one = ReadWithThreads(path, 1);
        EXPECT_EQ(one.threads, 1u);
        EXPECT_FALSE(one.parallel);
        EXPECT_FALSE(one.members);
        EXPECT_TRUE(one.text == content) << path;
        EXPECT_FALSE(one.failed) << path << ": " << one.error;
        bool const members = by_members.count(fs::path(path).filename().string()) > 0;
        for (unsigned threads : { 2u, 3u, 4u, 7u }) {
            SCOPED_TRACE(path + ", " + std::to_string(threads) + " threads");
            auto const r = ReadWithThreads(path, threads);
            EXPECT_EQ(r.parallel, !members);
            EXPECT_EQ(r.members, members);
            EXPECT_EQ(r.restarts, 0u);
            EXPECT_EQ(r.threads, threads);
            EXPECT_EQ(r.text.size(), content.size());
            EXPECT_TRUE(r.text == content);
            EXPECT_FALSE(r.failed) << r.error;
        }
    }
}

namespace {
    // Pieces of `bytes` of the compressed file for the gzip readers opened in its scope, the workers `gap` pieces ahead
    // of the reader (ParallelGzip::SetPieceBytes), so that a file of a few MB has hundreds; the defaults after.
    struct PieceBytes {
        explicit PieceBytes(size_t bytes, size_t gap = 0) { ParallelGzip::SetPieceBytes(bytes, gap); }
        ~PieceBytes() { ParallelGzip::SetPieceBytes(0); }
    };
}

// With pieces of 16 KiB, a few MB of gzip are hundreds of pieces: the workers decode pieces ahead (their 32 KB windows
// labelled, filled in once known; here from the piece after the reader's on, and as far ahead as by default), the
// reader's thread the others in order, and a slot without a block start belongs to the piece before. The bytes are the file's on any number of threads: one member at gzip's levels 1, 6 and 9,
// FASTQ of repeats (long copies from the window; few blocks), members of random sizes at every level with zero
// padding after them (in pieces too: none read by members here).
TEST(ThreadedGzStream, GzipPiecesDecodedAheadReadAsTheFile) {
    ScratchDir dir;
    MaxGzipMember const in_pieces(0);
    auto const random = RandomFastq(20000, "q", 21);  // ~6 MB, ~2.3 MB as gzip: ~140 pieces
    auto const repeats = Fastq(60000, "f");
    std::string several_text, several_bytes;
    std::mt19937 rng(5);
    for (int m = 0; m < 12; m++) {
        auto const part = RandomFastq(500 + rng() % 4000, "s" + std::to_string(m) + "_", 100 + static_cast<unsigned>(m));
        several_text += part;
        several_bytes += GzipBytes(part, false, 1 + m % 9);
    }
    struct File {
        std::string path, content;
        bool ahead;  // enough pieces that the workers decode some ahead
    };
    std::vector<File> const files = {
        { dir.Write("l1.fq.gz", GzipBytes(random, false, 1)), random, false },  // fixed codes only: one piece
        { dir.Write("l2.fq.gz", GzipBytes(random, false, 2)), random, true },
        { dir.Write("l6.fq.gz", GzipBytes(random)), random, true },
        { dir.Write("l9.fq.gz", GzipBytes(random, false, 9)), random, true },
        { dir.Write("repeats.fq.gz", GzipBytes(repeats)), repeats, false },
        { dir.Write("several.fq.gz", several_bytes + std::string(100, '\0')), several_text, true },
    };
    for (auto const& file : files) for (size_t const gap : { size_t{ 1 }, size_t{ 0 } }) {
        size_t ahead = 0, in_order = 0;
        PieceBytes const pieces(16 << 10, gap);
        for (unsigned threads : { 2u, 3u, 5u }) {
            SCOPED_TRACE(file.path + ", " + std::to_string(threads) + " threads, gap " + std::to_string(gap));
            InflateThreads const setting(threads);
            ThreadedGzIstream is(file.path.c_str());
            ASSERT_TRUE(is.rdbuf()->parallel_gzip());
            auto const text = ReadAllOf(is);
            EXPECT_FALSE(is.rdbuf()->read_failed()) << is.rdbuf()->read_error_message();
            EXPECT_EQ(text.size(), file.content.size());
            EXPECT_TRUE(text == file.content);
            auto const counts = is.rdbuf()->parallel_gzip_counts();
            ahead += counts.ahead;
            in_order += counts.in_order + counts.streamed;
        }
        EXPECT_GT(in_order, 0u) << file.path;
        if (file.ahead && gap == 1) EXPECT_GT(ahead, 0u) << file.path;
    }
}

// A block header found where there is none (here a real one, copied into data that gzip stores rather than
// compresses, at gzip's level 0): the piece from it is never used, the range before it is decoded on to the next start
// (or the member's end), and the bytes are the file's.
TEST(ThreadedGzStream, AFalseBlockStartIsNeverUsed) {
    ScratchDir dir;
    PieceBytes const pieces(16 << 10);
    // A dynamic block's header: the first block of a large gzip file (after its 10-byte gzip header), not its last.
    std::string const header = GzipBytes(RandomFastq(20000, "h", 9)).substr(10, 400);
    std::mt19937 rng(17);
    std::string content;
    while (content.size() < 1500000) {
        for (int i = 0; i < 90000; i++) content += static_cast<char>('A' + rng() % 26);
        content += header;
        content += '\n';
    }
    auto const path = dir.Write("stored.gz", GzipBytes(content, false, 0));
    for (unsigned threads : { 2u, 4u }) {
        SCOPED_TRACE(std::to_string(threads) + " threads");
        InflateThreads const setting(threads);
        ThreadedGzIstream is(path.c_str());
        ASSERT_TRUE(is.rdbuf()->parallel_gzip());
        auto const text = ReadAllOf(is);
        EXPECT_FALSE(is.rdbuf()->read_failed()) << is.rdbuf()->read_error_message();
        EXPECT_EQ(text.size(), content.size());
        EXPECT_TRUE(text == content);
        EXPECT_GT(is.rdbuf()->parallel_gzip_counts().dropped_starts, 0u);
    }
}

// A damaged gzip file fails in pieces on several threads as on one: cut in a member or in its trailer, a corrupt
// member, a wrong CRC-32 or length (ISIZE), a damaged second member header, data after the last member (here not read
// by members: ADamagedGzipMembersFileReadsAlikeOnAnyThreads). What is read before a cut is the file's; a corrupt member
// may give wrong bytes before the damage shows, on one thread too (ISA-L checks the CRC-32 at the member's end).
TEST(ThreadedGzStream, ADamagedGzipFileFailsOnSeveralThreads) {
    ScratchDir dir;
    MaxGzipMember const in_pieces(0);
    auto const content = RandomFastq(30000, "d", 7);
    std::string const good = GzipBytes(content), small = GzipBytes(Fastq(3000, "n"));
    struct Variant {
        std::string name, bytes;
        bool prefix = false;  // what is read is a prefix of the content
    };
    std::vector<Variant> variants{
        { "cut in half", good.substr(0, good.size() / 2), true },
        { "cut in the trailer", good.substr(0, good.size() - 4), true },
        { "no trailer", good.substr(0, good.size() - 8), true },
    };
    for (size_t at : { good.size() / 5, good.size() / 2, good.size() - 3000 }) {
        std::string corrupt = good;
        for (size_t i = at; i < at + 16; i++) corrupt[i] = static_cast<char>(corrupt[i] ^ 0x5a);
        variants.push_back({ "corrupt at " + std::to_string(at), corrupt });
    }
    std::string crc = good;
    crc[crc.size() - 6] = static_cast<char>(crc[crc.size() - 6] ^ 1);
    variants.push_back({ "wrong CRC-32", crc });
    std::string length = good;
    length[length.size() - 2] = static_cast<char>(length[length.size() - 2] ^ 1);
    variants.push_back({ "wrong length", length });
    std::string damaged = small;
    damaged[1] = 'X';
    variants.push_back({ "damaged second header", good + damaged, true });
    variants.push_back({ "data after the last member", good + small + "not gzip", true });
    variants.push_back({ "a byte after the last member", good + "x", true });
    for (size_t v = 0; v < variants.size(); v++) {
        auto const path = dir.Write("damaged" + std::to_string(v) + ".fq.gz", variants[v].bytes);
        auto const one = ReadWithThreads(path, 1);
        EXPECT_TRUE(one.failed) << variants[v].name;
        for (unsigned const threads : { 2u, 4u, 5u }) {  // 5: pieces of 16 KiB, the workers from the next one on
            SCOPED_TRACE(variants[v].name + ", " + std::to_string(threads) + " threads");
            PieceBytes const pieces(threads == 5 ? 16 << 10 : 0, threads == 5 ? 1 : 0);
            auto const r = ReadWithThreads(path, threads);
            EXPECT_TRUE(r.parallel);
            EXPECT_EQ(r.threads, threads);
            EXPECT_TRUE(r.failed);
            EXPECT_NE(r.error.find("corrupt gzip data"), std::string::npos) << r.error;
            if (variants[v].prefix) {
                size_t const n = std::min(r.text.size(), content.size());
                EXPECT_EQ(content.compare(0, n, r.text, 0, n), 0);
            }
        }
    }
}

// A gzip file cut anywhere fails on several threads, and what is read before is the file's (rapidgzip, which read gzip
// files in pieces until 2026-10-13, read many cut files as shorter ones, some as empty ones). Files of 100 KB to 5 MB of FASTQ at
// gzip levels 1, 6 and 9, cut at several places.
TEST(ThreadedGzStream, ACutGzipFileAlwaysFailsOnSeveralThreads) {
    ScratchDir dir;
    auto const text = RandomFastq(16000, "c", 13);  // ~5 MB
    int file = 0, parallel = 0, reads = 0;
    for (size_t const bytes : { size_t{100000}, size_t{500000}, size_t{2000000}, text.size() }) {
        auto const content = text.substr(0, bytes);
        for (int const level : { 1, 6, 9 }) {
            std::string const whole = GzipBytes(content, false, level);
            for (double const at : { 0.0, 0.25, 0.5, 0.9, 0.999 }) {
                size_t const cut = std::max<size_t>(20, static_cast<size_t>(static_cast<double>(whole.size()) * at));
                auto const path = dir.Write("cut" + std::to_string(file++) + ".fq.gz", whole.substr(0, cut));
                for (unsigned const threads : { 2u, 4u, 5u }) {  // 5: pieces of 16 KiB, the workers from the next one on
                    SCOPED_TRACE(std::to_string(bytes) + " bytes, level " + std::to_string(level) + ", cut at " +
                                 std::to_string(cut) + " of " + std::to_string(whole.size()) + ", " + std::to_string(threads) + " threads");
                    PieceBytes const pieces(threads == 5 ? 16 << 10 : 0, threads == 5 ? 1 : 0);
                    auto const r = ReadWithThreads(path, threads);
                    EXPECT_TRUE(r.failed) << r.text.size() << " bytes read, " << (r.parallel ? "in pieces" : "in order");
                    size_t const n = std::min(r.text.size(), content.size());
                    EXPECT_EQ(content.compare(0, n, r.text, 0, n), 0);
                    parallel += r.parallel;
                    reads++;
                }
            }
        }
    }
    // A file cut into its header is read on one thread; the others in pieces.
    EXPECT_GT(parallel, reads * 3 / 4);
}

// A pipe is read on one thread whatever the setting: its pieces would have to be kept.
TEST(ThreadedGzStream, AGzipPipeIsInflatedOnOneThread) {
    ScratchDir dir;
    std::signal(SIGPIPE, SIG_IGN);
    auto const content = RandomFastq(5000, "p", 2);
    std::string const bytes = GzipBytes(content);
    auto const fifo = (dir.path / "reads.fifo").string();
    ASSERT_EQ(::mkfifo(fifo.c_str(), 0600), 0);
    std::thread writer([&fifo, &bytes] { std::ofstream(fifo, std::ios::binary) << bytes; });
    {
        InflateThreads const setting(4);
        ThreadedGzIstream is(fifo.c_str());
        EXPECT_EQ(is.rdbuf()->inflate_threads(), 1u);
        EXPECT_FALSE(is.rdbuf()->parallel_gzip());
        EXPECT_EQ(ReadAllOf(is), content);
        EXPECT_FALSE(is.rdbuf()->read_failed()) << is.rdbuf()->read_error_message();
    }
    writer.join();
}

namespace {
    // text as gzip members of min to max bytes of it each (sizes by rng), at levels 0 to 9 in turn (stored members among
    // them), every fifth with every optional header field and an empty member now and then: as tools that compress a
    // file in chunks write it.
    std::string GzipMembers(std::string_view text, size_t min, size_t max, unsigned seed) {
        std::mt19937 rng(seed);
        std::string out;
        int m = 0;
        for (size_t from = 0; from < text.size(); m++) {
            size_t const n = std::min<size_t>(text.size() - from, min + rng() % (max - min + 1));
            out += GzipBytes(text.substr(from, n), m % 5 == 4, m % 10);
            if (m % 13 == 7) out += GzipBytes("");
            from += n;
        }
        return out;
    }
}

// A gzip file of many small members (as tools that compress a file in chunks write it) is read by members on several
// threads, as BGZF is, and gives the bytes one thread gives: members of 1-64 KB at every level, headers with every
// optional field, empty members, BGZF blocks after them, zero padding after the last. With members read whole up to
// 256 KiB, a member longer than that among them, from which on the file is read on one thread. A file whose second
// member begins further on than that is read in pieces.
TEST(ThreadedGzStream, GzipMembersOnSeveralThreadsReadAsOnOne) {
    ScratchDir dir;
    auto const text = RandomFastq(20000, "u", 41);  // ~6 MB
    auto const tail = Fastq(3000, "t");
    auto const stored = RandomFastq(3000, "l", 47);  // ~0.9 MB, as much stored
    std::string const members = GzipMembers(text, 1000, 64000, 43);
    struct File {
        std::string name, bytes, content;
        size_t max_member;  // the longest member read whole
        bool by_members;
    };
    std::vector<File> const files = {
        { "members.fq.gz", members, text, ThreadedGzStreambuf::kMaxGzipMember, true },
        { "then_bgzf.fq.gz", members + FileBytes(Bgzf(dir, "b.fq.gz", tail)), text + tail, ThreadedGzStreambuf::kMaxGzipMember, true },
        { "padded.fq.gz", members + std::string(70000, '\0'), text, ThreadedGzStreambuf::kMaxGzipMember, true },
        { "long.fq.gz", GzipMembers(text.substr(0, 1000000), 1000, 30000, 53) + GzipBytes(stored, false, 0) + GzipMembers(tail, 500, 5000, 59),
          text.substr(0, 1000000) + stored + tail, 256 << 10, true },
        { "late.fq.gz", GzipBytes(text) + GzipMembers(tail, 500, 5000, 61), text + tail, 256 << 10, false },
    };
    for (auto const& file : files) {
        auto const path = dir.Write(file.name, file.bytes);
        MaxGzipMember const limit(file.max_member);
        auto const one = ReadWithThreads(path, 1);
        EXPECT_FALSE(one.members);
        EXPECT_TRUE(one.text == file.content) << file.name;
        for (unsigned threads : { 2u, 3u, 4u, 7u }) {
            SCOPED_TRACE(file.name + ", " + std::to_string(threads) + " threads");
            auto const r = ReadWithThreads(path, threads);
            EXPECT_EQ(r.members, file.by_members);
            EXPECT_EQ(r.parallel, !file.by_members);
            EXPECT_EQ(r.threads, threads);
            EXPECT_EQ(r.restarts, 0u);
            EXPECT_EQ(r.text.size(), file.content.size());
            EXPECT_TRUE(r.text == file.content);
            EXPECT_FALSE(r.failed) << r.error;
        }
    }
}

// A gzip member header found by chance in a member's data (here written into data that gzip stores rather than
// compresses, at level 0) cuts the file wrongly. The batch that ends there fails, the file is read on one thread from
// that batch on, and the bytes are the file's. One such header follows a length of just over an output block (1 MiB),
// so that a batch ends at it; the others follow small ones, and the members they cut stay within their batch.
TEST(ThreadedGzStream, AGzipMemberHeaderFoundByChanceIsNeverUsed) {
    ScratchDir dir;
    auto const before = RandomFastq(3000, "a", 61), after = RandomFastq(3000, "z", 67);
    // A trailer's length, then a member header (no flags, extra flags 0, Unix) and a stored block of 5 bytes.
    auto const fake = [](uint32_t isize) {
        std::string f;
        for (int i = 0; i < 4; i++) f += static_cast<char>((isize >> (8 * i)) & 0xff);
        return f + std::string("\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\x03\x01\x05\x00\xfa\xff", 15) + "hello";
    };
    std::mt19937 rng(71);
    std::string stored;
    for (int i = 0; i < 40; i++) {
        for (int j = 0; j < 3000; j++) stored += static_cast<char>('A' + rng() % 26);
        stored += fake(i == 20 ? (1u << 20) + 1 : 16);
    }
    std::string const content = before + stored + after;
    auto const path = dir.Write("chance.gz", GzipMembers(before, 1000, 20000, 73) + GzipBytes(stored, false, 0) +
                                                 GzipMembers(after, 1000, 20000, 79));
    auto const one = ReadWithThreads(path, 1);
    ASSERT_TRUE(one.text == content);
    for (unsigned threads : { 2u, 4u }) {
        SCOPED_TRACE(std::to_string(threads) + " threads");
        auto const r = ReadWithThreads(path, threads);
        EXPECT_TRUE(r.members);
        EXPECT_FALSE(r.failed) << r.error;
        EXPECT_EQ(r.text.size(), content.size());
        EXPECT_TRUE(r.text == content);
        EXPECT_GE(r.restarts, 1u);
    }
}

// A damaged gzip file of many members reads the same by members on several threads as on one: the same bytes up to the
// damage and the same error (its member's number, its byte), wherever the damage is: a flipped byte in a member's data,
// a damaged member header, a wrong CRC-32 or length, a cut inside a member, data after the last member.
TEST(ThreadedGzStream, ADamagedGzipMembersFileReadsAlikeOnAnyThreads) {
    ScratchDir dir;
    auto const text = RandomFastq(12000, "d", 83);  // ~3.6 MB
    std::vector<std::string> parts;               // the members
    std::mt19937 rng(89);
    for (size_t from = 0; from < text.size();) {
        size_t const n = std::min<size_t>(text.size() - from, 2000 + rng() % 30000);
        parts.push_back(GzipBytes(std::string_view(text).substr(from, n)));
        from += n;
    }
    auto const join = [&parts](size_t damaged, std::string const& with) {
        std::string all;
        for (size_t i = 0; i < parts.size(); i++) all += i == damaged ? with : parts[i];
        return all;
    };
    auto const flipped = [](std::string bytes, size_t at, char by) {
        bytes[at] = static_cast<char>(bytes[at] ^ by);
        return bytes;
    };
    size_t const mid = parts.size() / 2;
    std::string const& m = parts[mid];
    std::string const whole = join(SIZE_MAX, "");
    std::vector<std::pair<std::string, std::string>> const variants = {
        { "flipped in the first member", join(0, flipped(parts[0], parts[0].size() / 2, 0x5a)) },
        { "flipped at a third", join(mid, flipped(m, m.size() / 3, 0x5a)) },
        { "flipped at a half", join(mid, flipped(m, m.size() / 2, 0x5a)) },
        { "damaged header", join(mid, flipped(m, 1, 0x33)) },
        { "wrong CRC-32", join(mid, flipped(m, m.size() - 6, 1)) },
        { "wrong length", join(mid, flipped(m, m.size() - 2, 1)) },
        { "cut", whole.substr(0, whole.size() * 2 / 3) },
        { "data after the last member", whole + "not gzip" },
    };
    for (size_t v = 0; v < variants.size(); v++) {
        auto const path = dir.Write("damaged" + std::to_string(v) + ".fq.gz", variants[v].second);
        auto const one = ReadWithThreads(path, 1);
        EXPECT_TRUE(one.failed) << variants[v].first;
        for (unsigned threads : { 2u, 4u, 7u }) {
            SCOPED_TRACE(variants[v].first + ", " + std::to_string(threads) + " threads");
            auto const r = ReadWithThreads(path, threads);
            EXPECT_TRUE(r.members);
            EXPECT_EQ(r.failed, one.failed);
            EXPECT_EQ(r.error, one.error);
            EXPECT_EQ(r.text.size(), one.text.size());
            EXPECT_TRUE(r.text == one.text);
        }
    }
}

// A gzip file read in pieces (its first member longer than the longest read whole, here 64 KiB) goes on by members
// where small members follow a member, and gives the bytes one thread gives: small members after a large one, after two
// large ones, and with a large member again after them (the rest then on one thread). Damage after the switch reads
// as on one thread: the same bytes and error.
TEST(ThreadedGzStream, AGzipFileInPiecesGoesOnByMembers) {
    ScratchDir dir;
    MaxGzipMember const limit(64 << 10);
    auto const big = RandomFastq(3000, "b", 91);   // ~0.9 MB, ~0.3 MB as gzip
    auto const next = RandomFastq(2000, "n", 97);  // ~0.2 MB as gzip
    auto const text = RandomFastq(6000, "s", 93);  // in members of 1-20 KB
    std::string const members = GzipMembers(text, 1000, 20000, 99);
    std::vector<std::pair<std::string, std::string>> const files = {  // (bytes, content)
        { GzipBytes(big) + members, big + text },
        { GzipBytes(big) + GzipBytes(next) + members, big + next + text },
        { GzipBytes(big) + members + GzipBytes(next) + GzipMembers(text, 1000, 20000, 101), big + text + next + text },
    };
    for (size_t f = 0; f < files.size(); f++) {
        auto const path = dir.Write("mixed" + std::to_string(f) + ".fq.gz", files[f].first);
        for (unsigned threads : { 2u, 4u }) {
            for (size_t const piece : { size_t{ 0 }, size_t{ 16 } << 10 }) {
                SCOPED_TRACE("file " + std::to_string(f) + ", " + std::to_string(threads) + " threads, pieces of " + std::to_string(piece));
                PieceBytes const pieces(piece);
                auto const r = ReadWithThreads(path, threads);
                EXPECT_TRUE(r.parallel);
                EXPECT_FALSE(r.members);
                EXPECT_TRUE(r.switched);
                EXPECT_EQ(r.restarts, 0u);
                EXPECT_EQ(r.text.size(), files[f].second.size());
                EXPECT_TRUE(r.text == files[f].second);
                EXPECT_FALSE(r.failed) << r.error;
            }
        }
    }
    std::string const& good = files[0].first;
    size_t const after = good.size() - members.size();  // where the members begin
    std::string flipped = good;
    flipped[after + members.size() / 2] = static_cast<char>(flipped[after + members.size() / 2] ^ 0x5a);
    std::string header = good;
    size_t const second = after + gzip_member::Next(reinterpret_cast<unsigned char const*>(members.data()), members.size(), 0, 64 << 10);
    header[second + 1] = static_cast<char>(header[second + 1] ^ 0x33);
    std::vector<std::pair<std::string, std::string>> const damaged = {
        { "flipped", flipped }, { "damaged header", header }, { "cut", good.substr(0, after + members.size() * 2 / 3) },
        { "data after the last member", good + "not gzip" },
    };
    for (size_t v = 0; v < damaged.size(); v++) {
        auto const path = dir.Write("damaged" + std::to_string(v) + ".fq.gz", damaged[v].second);
        auto const one = ReadWithThreads(path, 1);
        EXPECT_TRUE(one.failed) << damaged[v].first;
        for (unsigned threads : { 2u, 4u }) {
            SCOPED_TRACE(damaged[v].first + ", " + std::to_string(threads) + " threads");
            auto const r = ReadWithThreads(path, threads);
            EXPECT_TRUE(r.switched);
            EXPECT_EQ(r.failed, one.failed);
            EXPECT_EQ(r.error, one.error);
            EXPECT_EQ(r.text.size(), one.text.size());
            EXPECT_TRUE(r.text == one.text);
        }
    }
}

// Read pairs taken by several threads from two gzip files inflated on several threads each, in pieces (one member
// each) or by members: every pair once.
TEST(ThreadedGzStream, ThreadsTakeEveryPairOnceFromParallelGzip) {
    ScratchDir dir;
    size_t const n = 30000;
    std::mt19937 rng(11);
    std::string r1, r2;
    for (size_t i = 0; i < n; i++) {
        std::string seq(100 + i % 50, 'A');
        for (auto& c : seq) c = "ACGT"[rng() % 4];
        r1 += "@r" + std::to_string(i) + "/1\n" + seq + "\n+\n" + std::string(seq.size(), 'I') + "\n";
        r2 += "@r" + std::to_string(i) + "/2\n" + seq + "\n+\n" + std::string(seq.size(), 'J') + "\n";
    }
    for (bool const members : { false, true }) {
        SCOPED_TRACE(members ? "by members" : "in pieces");
        auto const p1 = dir.Write(members ? "R1m.fq.gz" : "R1.fq.gz", members ? GzipMembers(r1, 1000, 40000, 3) : GzipBytes(r1));
        auto const p2 = dir.Write(members ? "R2m.fq.gz" : "R2.fq.gz", members ? GzipMembers(r2, 1000, 40000, 5) : GzipBytes(r2));
        InflateThreads const setting(3);
        ThreadedGzIstream is1(p1.c_str()), is2(p2.c_str());
        ASSERT_EQ(is1.rdbuf()->parallel_gzip(), !members);
        ASSERT_EQ(is1.rdbuf()->opened_as_gzip_members(), members);
        ASSERT_EQ(is1.rdbuf()->inflate_threads(), 3u);
        SeqReaderPE global(is1, is2);
        std::vector<std::vector<size_t>> taken(4);
#pragma omp parallel num_threads(4)
        {
            SeqReaderPE reader(global);
            FastxRecord a, b;
            auto& mine = taken[static_cast<size_t>(omp_get_thread_num())];
            while (reader(a, b)) {
                if (a.id.substr(0, a.id.size() - 2) == b.id.substr(0, b.id.size() - 2)) mine.push_back(std::stoul(a.id.substr(1)));
            }
        }
        std::vector<size_t> all;
        for (auto const& t : taken) all.insert(all.end(), t.begin(), t.end());
        std::sort(all.begin(), all.end());
        ASSERT_EQ(all.size(), n);
        for (size_t i = 0; i < n; i++) ASSERT_EQ(all[i], i);
        EXPECT_FALSE(is1.rdbuf()->read_failed()) << is1.rdbuf()->read_error_message();
        EXPECT_FALSE(is2.rdbuf()->read_failed()) << is2.rdbuf()->read_error_message();
    }
}

// A cut or corrupt BGZF file reads the same on several threads as on one: the same bytes up to the damage, the same
// error, wherever the damage is (a flipped byte in each of many blocks, a damaged header, a cut at or inside a block,
// a missing end-of-file block, a gzip member after the blocks that is cut).
TEST(ThreadedGzStream, ADamagedBgzfFileReadsAlikeOnAnyThreads) {
    ScratchDir dir;
    auto const content = Fastq(40000, "d");
    auto const good = FileBytes(Bgzf(dir, "good.fq.gz", content));
    std::vector<std::string> variants;
    for (size_t at : { size_t{100}, good.size() / 7, good.size() / 3, good.size() / 2, good.size() - 400 }) {
        std::string flipped = good;
        flipped[at] = static_cast<char>(flipped[at] ^ 0x5a);
        variants.push_back(flipped);
    }
    std::string header = good;
    size_t const second_block = bgzf::BlockSize(reinterpret_cast<unsigned char const*>(good.data()));
    header[second_block] = 0;  // the second block's gzip magic (a changed BC subfield would still be a valid gzip member)
    variants.push_back(header);
    variants.push_back(good.substr(0, good.size() - sizeof(bgzf::kEof)));  // no end-of-file block
    variants.push_back(good.substr(0, good.size() / 2));                   // cut inside a block
    variants.push_back(good.substr(0, second_block * 5));                  // cut at a block boundary
    variants.push_back(good.substr(0, second_block * 5 + 10));             // cut inside a header
    std::string const member = GzipBytes(content);
    variants.push_back(good + member.substr(0, member.size() / 2));        // a cut gzip member after the blocks
    for (size_t v = 0; v < variants.size(); v++) {
        auto const path = dir.Write("damaged" + std::to_string(v) + ".fq.gz", variants[v]);
        auto const one = ReadWithThreads(path, 1);
        EXPECT_TRUE(one.failed) << "variant " << v;
        EXPECT_EQ(content.compare(0, std::min(one.text.size(), content.size()), one.text, 0, std::min(one.text.size(), content.size())), 0);
        for (unsigned threads : { 2u, 4u, 7u }) {
            SCOPED_TRACE("variant " + std::to_string(v) + ", " + std::to_string(threads) + " threads");
            auto const r = ReadWithThreads(path, threads);
            EXPECT_EQ(r.failed, one.failed);
            EXPECT_EQ(r.error, one.error);
            EXPECT_EQ(r.text.size(), one.text.size());
            EXPECT_TRUE(r.text == one.text);
        }
    }
}

// Read pairs taken by several threads from two BGZF files inflated on several threads each: every pair once.
TEST(ThreadedGzStream, ThreadsTakeEveryPairOnceFromParallelBgzf) {
    ScratchDir dir;
    size_t const n = 30000;
    std::string r1, r2;
    for (size_t i = 0; i < n; i++) {
        std::string const seq(100 + i % 50, "ACGT"[i % 4]);
        r1 += "@r" + std::to_string(i) + "/1\n" + seq + "\n+\n" + std::string(seq.size(), 'I') + "\n";
        r2 += "@r" + std::to_string(i) + "/2\n" + seq + "\n+\n" + std::string(seq.size(), 'J') + "\n";
    }
    auto const p1 = Bgzf(dir, "R1.fq.gz", r1), p2 = Bgzf(dir, "R2.fq.gz", r2);
    InflateThreads const setting(3);
    ThreadedGzIstream is1(p1.c_str()), is2(p2.c_str());
    ASSERT_EQ(is1.rdbuf()->inflate_threads(), 3u);
    SeqReaderPE global(is1, is2);
    std::vector<std::vector<size_t>> taken(4);
#pragma omp parallel num_threads(4)
    {
        SeqReaderPE reader(global);
        FastxRecord a, b;
        auto& mine = taken[static_cast<size_t>(omp_get_thread_num())];
        while (reader(a, b)) {
            if (a.id.substr(0, a.id.size() - 2) == b.id.substr(0, b.id.size() - 2)) mine.push_back(std::stoul(a.id.substr(1)));
        }
    }
    std::vector<size_t> all;
    for (auto const& t : taken) all.insert(all.end(), t.begin(), t.end());
    std::sort(all.begin(), all.end());
    ASSERT_EQ(all.size(), n);
    for (size_t i = 0; i < n; i++) ASSERT_EQ(all[i], i);
}

// A zstd file's frames decompressed on several threads (as a BGZF file's blocks) give the bytes one thread gives: many
// frames of random sizes (as protal's simulator writes a frame per batch of reads) with skippable and empty frames among
// them, a frame larger than an output block, a frame whose header gives no size (a stream's, zstd --long here) after
// others and first (the rest of the file then decompressed as a stream), and one frame alone.
TEST(ThreadedGzStream, ZstdFramesOnSeveralThreadsReadAsOnOne) {
    ScratchDir dir;
    auto const text = RandomFastq(30000, "z", 23);  // ~9 MB
    std::mt19937 rng(29);
    std::string frames;
    for (size_t from = 0; from < text.size();) {
        size_t const n = std::min<size_t>(text.size() - from, 1000 + rng() % 400000);
        frames += ZstdBytes(std::string_view(text).substr(from, n));
        if (rng() % 5 == 0) frames += SkippableFrame(rng() % 300);
        if (rng() % 7 == 0) frames += ZstdBytes("");
        from += n;
    }
    std::vector<std::pair<std::string, std::string>> const files = {  // (name, bytes), each of `text`
        { "frames.fq.zst", frames },
        { "large.fq.zst", ZstdBytes(text.substr(0, 100000)) + ZstdBytes(text.substr(100000, 3000000)) +
                              ZstdBytes(text.substr(3100000)) },
        { "stream_after.fq.zst", ZstdBytes(text.substr(0, 2000000)) + ZstdBytes(text.substr(2000000, 1000000), 27) +
                                     ZstdBytes(text.substr(3000000)) },
        { "stream_first.fq.zst", ZstdBytes(text.substr(0, 2000000), 27) + ZstdBytes(text.substr(2000000)) },
        { "one.fq.zst", ZstdBytes(text) },
    };
    for (auto const& [name, bytes] : files) {
        auto const path = dir.Write(name, bytes);
        for (unsigned threads : { 1u, 2u, 3u, 4u, 7u }) {
            SCOPED_TRACE(name + ", " + std::to_string(threads) + " threads");
            auto const r = ReadWithThreads(path, threads);
            // A file whose first frame gives no size is read as a stream from the start, on its own thread.
            EXPECT_EQ(r.threads, name == "stream_first.fq.zst" ? 1u : threads);
            EXPECT_EQ(r.text.size(), text.size());
            EXPECT_TRUE(r.text == text);
            EXPECT_FALSE(r.failed) << r.error;
        }
    }
}

// A cut or corrupt zstd file fails on several threads as on one, what is read before the damage being the file's: cut
// inside a frame, the last bytes of a frame's checksum missing, a corrupt frame, data after the last frame.
TEST(ThreadedGzStream, ADamagedZstdFileFailsOnSeveralThreads) {
    ScratchDir dir;
    auto const first = RandomFastq(5000, "a", 31), second = RandomFastq(5000, "b", 37), content = first + second;
    std::string const one = ZstdBytes(first), two = ZstdBytes(second);
    std::string corrupt = two;
    for (size_t i = corrupt.size() / 2; i < corrupt.size() / 2 + 16; i++) corrupt[i] = static_cast<char>(corrupt[i] ^ 0x5a);
    std::vector<std::pair<std::string, std::string>> const variants = {  // (bytes, what the error says)
        { one + two.substr(0, two.size() / 2), "the file ends inside zstd frame 2 (truncated file?)" },
        { one.substr(0, one.size() - 3), "the file ends inside zstd frame 1 (truncated file?)" },
        { one + corrupt, "in zstd frame 2 (corrupt file?)" },
        { one + "not a zstd frame", "zstd frame" },
    };
    for (size_t v = 0; v < variants.size(); v++) {
        auto const path = dir.Write("damaged" + std::to_string(v) + ".fq.zst", variants[v].first);
        for (unsigned threads : { 1u, 3u }) {
            SCOPED_TRACE(variants[v].second + ", " + std::to_string(threads) + " threads");
            auto const r = ReadWithThreads(path, threads);
            EXPECT_TRUE(r.failed);
            EXPECT_NE(r.error.find(variants[v].second), std::string::npos) << r.error;
            size_t const n = std::min(r.text.size(), content.size());
            EXPECT_EQ(content.compare(0, n, r.text, 0, n), 0);
        }
    }
}

TEST(ThreadedGzStream, InflateThreadsFollowTheRunsThreads) {
    EXPECT_EQ(ThreadedGzStreambuf::InflateThreadsFor(1, 2), 1u);
    EXPECT_EQ(ThreadedGzStreambuf::InflateThreadsFor(8, 2), 1u);
    EXPECT_EQ(ThreadedGzStreambuf::InflateThreadsFor(16, 2), 2u);
    EXPECT_EQ(ThreadedGzStreambuf::InflateThreadsFor(32, 2), 4u);
    EXPECT_EQ(ThreadedGzStreambuf::InflateThreadsFor(128, 2), 4u);
    EXPECT_EQ(ThreadedGzStreambuf::InflateThreadsFor(12, 1), 3u);
    // Other gzip in pieces: a thread per 4 and file, 1 to 8.
    EXPECT_EQ(ThreadedGzStreambuf::GzipThreadsFor(8, 2), 1u);
    EXPECT_EQ(ThreadedGzStreambuf::GzipThreadsFor(16, 2), 2u);
    EXPECT_EQ(ThreadedGzStreambuf::GzipThreadsFor(32, 2), 4u);
    EXPECT_EQ(ThreadedGzStreambuf::GzipThreadsFor(48, 2), 6u);
    EXPECT_EQ(ThreadedGzStreambuf::GzipThreadsFor(128, 2), 8u);
    EXPECT_EQ(ThreadedGzStreambuf::GzipThreadsFor(8, 1), 2u);
    EXPECT_EQ(ThreadedGzStreambuf::GzipThreadsFor(16, 1), 4u);
}

TEST(ThreadedGzStream, ClosingBeforeTheEndStopsTheInflatingThread) {
    ScratchDir dir;
    auto const path = Gzip(dir, "reads.fq.gz", Fastq(40000, "c"));
    for (int i = 0; i < 20; i++) {  // the inflating thread is blocked on full blocks or still inflating
        ThreadedGzIstream is(path.c_str());
        std::string line;
        ASSERT_TRUE(std::getline(is, line));
        EXPECT_EQ(line, "@c0/1");
        if (i % 2) is.close();  // else the destructor closes it
    }
    // A BGZF file on several threads: the dispatcher and the workers blocked or busy when it is closed; a gzip file
    // in pieces alike.
    auto const bgzf = Bgzf(dir, "reads.bgzf.gz", Fastq(40000, "b"));
    auto const parallel = dir.Write("random.fq.gz", GzipBytes(RandomFastq(40000, "g", 4)));
    InflateThreads const setting(4);
    for (int i = 0; i < 20; i++) {
        ThreadedGzIstream is(bgzf.c_str());
        std::string line;
        ASSERT_TRUE(std::getline(is, line));
        EXPECT_EQ(line, "@b0/1");
        if (i % 2) is.close();
    }
    for (int i = 0; i < 20; i++) {
        ThreadedGzIstream is(parallel.c_str());
        ASSERT_TRUE(is.rdbuf()->parallel_gzip());
        std::string line;
        ASSERT_TRUE(std::getline(is, line));
        EXPECT_EQ(line, "@g0/1");
        if (i % 2) is.close();
    }
    // A gzip file of many members by members, alike; and one read in pieces that goes on by members, closed before,
    // while and after it does.
    auto const members = dir.Write("members.fq.gz", GzipMembers(RandomFastq(40000, "m", 5), 2000, 50000, 7));
    for (int i = 0; i < 20; i++) {
        ThreadedGzIstream is(members.c_str());
        ASSERT_TRUE(is.rdbuf()->opened_as_gzip_members());
        std::string line;
        ASSERT_TRUE(std::getline(is, line));
        EXPECT_EQ(line, "@m0/1");
        if (i % 2) is.close();
    }
    MaxGzipMember const limit(64 << 10);
    auto const mixed = dir.Write("mixed.fq.gz", GzipBytes(RandomFastq(3000, "x", 9)) + GzipMembers(RandomFastq(20000, "y", 11), 2000, 30000, 13));
    for (int i = 0; i < 30; i++) {
        ThreadedGzIstream is(mixed.c_str());
        ASSERT_TRUE(is.rdbuf()->parallel_gzip());
        std::string line;
        for (int l = 0; l < 4 * 200 * i && std::getline(is, line); l++) {}  // 200 reads more each time (3,000 before the members)
        if (i % 2) is.close();
    }
}

TEST(ThreadedGzStream, ThreadsTakeEveryPairOnceInStep) {
    ScratchDir dir;
    size_t const n = 30000;
    auto const r1 = Gzip(dir, "r1.fq.gz", Fastq(n, "pair"));
    auto const r2 = Gzip(dir, "r2.fq.gz", Fastq(n, "pair"));
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
    auto const path = Gzip(dir, "se.fq.gz", Fastq(n, "se"));
    for (bool const batches : { true, false }) {
        SCOPED_TRACE(batches ? "SeqReaderSE" : "SeqReader");
        ThreadedGzIstream is(path.c_str());
        SeqReaderSE batch_global{ is };
        SeqReader block_global{ is };
        std::vector<std::string> ids;
        bool blocks_read = true;  // every thread's SeqReader (which has no UpdateSuccess) ended without an error
#pragma omp parallel num_threads(8) shared(batch_global, block_global, ids, batches, blocks_read) default(none)
        {
            SeqReaderSE batch_reader{ batch_global };
            SeqReader block_reader{ block_global };
            FastxRecord record;
            std::vector<std::string> mine;
            while (batches ? batch_reader(record) : block_reader(record)) mine.push_back(record.id);
#pragma omp critical(test_ids)
            {
                ids.insert(ids.end(), mine.begin(), mine.end());
                if (batches) batch_global.UpdateSuccess(batch_reader);
                else blocks_read = blocks_read && block_reader.Success();
            }
        }
        EXPECT_TRUE(batch_global.Success());
        EXPECT_TRUE(blocks_read);
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
    for (auto const& path : { Gzip(dir, "lines.gz", content), dir.Write("lines.txt", content) }) {
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

// A sample's lanes: several files read one after the other as one stream, each in its own format (plain, gzip, BGZF
// on one thread or several, zstd), an empty one among them, and a '\n' after a file whose last line lacks it.
TEST(ThreadedGzStream, SeveralFilesReadAsOneStream) {
    ScratchDir dir;
    auto const a = Fastq(20000, "a"), b = Fastq(3000, "b"), c = Fastq(15000, "c"), d = Fastq(500, "d");
    std::string const b_open = b.substr(0, b.size() - 1);  // its last line without '\n'
    for (unsigned const threads : { 1u, 3u }) {
        SCOPED_TRACE("inflate threads " + std::to_string(threads));
        ThreadedGzStreambuf::SetInflateThreads(threads);
        ThreadedGzStreambuf::SetGzipThreads(threads);
        std::vector<std::string> const files{ dir.Write("a.fq", a), Gzip(dir, "b.fq.gz", b_open), dir.Write("empty.fq", ""),
                                              Bgzf(dir, "c.fq.gz", c), dir.Write("d.fq.zst", ZstdBytes(d)) };
        std::vector<std::string> opened;
        ThreadedGzIstream is(files, [&opened](std::string const& path, ThreadedGzStreambuf const&) { opened.push_back(path); });
        EXPECT_EQ(ReadAllOf(is), a + b + c + d);
        EXPECT_FALSE(is.rdbuf()->read_failed()) << is.rdbuf()->read_error_message();
        EXPECT_EQ(opened, files);
        EXPECT_EQ(is.rdbuf()->files(), files.size());
        EXPECT_EQ(is.rdbuf()->current_path(), files.back());

        ThreadedGzIstream by_line(files);  // getline (FASTA's way) through the same stream
        std::stringstream reference(a + b + c + d);
        EXPECT_EQ(Lines(by_line), Lines(reference));
    }
    ThreadedGzStreambuf::SetInflateThreads(1);
    ThreadedGzStreambuf::SetGzipThreads(1);
}

// Of several files, a batch ends with its file once it holds a line, so that the mate files of lanes stay in step
// (SeqReaderPE); one that begins at a file's end goes on into the next file, and a file's open last line ends its
// batch with the '\n' it lacked.
TEST(ThreadedGzStream, ABatchEndsWithItsFile) {
    ScratchDir dir;
    auto const a = Fastq(10, "a"), b = Fastq(7, "b"), c = Fastq(5, "c");
    std::vector<std::string> const files{ dir.Write("a.fq", a), Gzip(dir, "b.fq.gz", b.substr(0, b.size() - 1)),
                                          dir.Write("empty.fq", ""), dir.Write("c.fq", c) };
    ThreadedGzIstream is(files);
    std::vector<std::string> batches;
    std::vector<size_t> counts;
    while (true) {
        std::string batch;
        size_t const n = is.rdbuf()->TakeLines(4 * 32, batch);
        if (n == 0) break;
        counts.push_back(n);
        batches.push_back(batch);
    }
    EXPECT_EQ(counts, (std::vector<size_t>{ 40, 28, 20 }));
    EXPECT_EQ(batches, (std::vector<std::string>{ a, b, c }));

    // Batches of 3 lines: the file's last batch is short, and the next file's begins with its first line.
    ThreadedGzIstream again(files);
    std::string all;
    counts.clear();
    for (size_t n; (n = again.rdbuf()->TakeLines(12, all)) > 0;) counts.push_back(n);
    EXPECT_EQ(counts, (std::vector<size_t>{ 12, 12, 12, 4, 12, 12, 4, 12, 8 }));
    EXPECT_EQ(all, a + b + c);
}

// Of several files, one that is cut or corrupt, or cannot be read when its turn comes, ends the stream there: what
// came before it is read, nothing after it, and the error names the file. A first file that cannot be read fails
// the open, as one file does.
TEST(ThreadedGzStream, AFileOfSeveralThatFailsEndsTheStream) {
    ScratchDir dir;
    auto const a = Fastq(2000, "a"), b = Fastq(20000, "b"), c = Fastq(2000, "c");
    auto const good_a = dir.Write("a.fq", a), good_c = dir.Write("c.fq", c);
    auto const gz = FileBytes(Gzip(dir, "full.fq.gz", b));
    auto const cut = dir.Write("cut.fq.gz", gz.substr(0, gz.size() / 2));
    auto const missing = (dir.path / "missing.fq").string();
    auto const npos = std::string::npos;
    {
        ThreadedGzIstream is(std::vector<std::string>{ good_a, cut, good_c });
        auto const text = ReadAllOf(is);
        EXPECT_EQ(text.compare(0, a.size(), a), 0);
        EXPECT_LT(text.size(), a.size() + b.size());
        EXPECT_EQ(text.find("@c0/1"), npos);
        ASSERT_TRUE(is.rdbuf()->read_failed());
        EXPECT_EQ(is.rdbuf()->read_error_message().rfind(cut + ": ", 0), 0u) << is.rdbuf()->read_error_message();
    }
    {
        ThreadedGzIstream is(std::vector<std::string>{ good_a, missing, good_c });
        EXPECT_EQ(ReadAllOf(is), a);
        ASSERT_TRUE(is.rdbuf()->read_failed());
        EXPECT_EQ(is.rdbuf()->read_error_message(), missing + ": it cannot be read (" + std::strerror(ENOENT) + ")");
    }
    ThreadedGzIstream first_missing(std::vector<std::string>{ missing, good_a });
    EXPECT_FALSE(first_missing.rdbuf()->is_open());
    EXPECT_EQ(first_missing.rdbuf()->current_path(), missing);
    EXPECT_EQ(first_missing.rdbuf()->open_error(), std::strerror(ENOENT));
}

// Read pairs of a sample's lanes (two files a mate, matched by position) on several threads: every pair once and
// in step; a lane whose mates differ in number fails the reader, even where the totals agree, instead of pairing
// the next lane's reads with the wrong mates.
TEST(SeqReader, PairsOfLanesStayInStep) {
    ScratchDir dir;
    auto r1 = std::vector<std::string>{ Gzip(dir, "l1_r1.fq.gz", Fastq(30000, "x")), dir.Write("l2_r1.fq", Fastq(20000, "y")) };
    auto r2 = std::vector<std::string>{ dir.Write("l1_r2.fq", Fastq(30000, "x")), Gzip(dir, "l2_r2.fq.gz", Fastq(20000, "y")) };
    ThreadedGzIstream is1(r1), is2(r2);
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
    EXPECT_FALSE(reader_global.MatesOutOfStep());
    EXPECT_EQ(ids.size(), 50000u);
    EXPECT_EQ(std::set<std::string>(ids.begin(), ids.end()).size(), 50000u);

    // Lane 1's first file has a read more, lane 2's second file one: 53 reads a side.
    ThreadedGzIstream uneven1(std::vector<std::string>{ dir.Write("u1_r1.fq", Fastq(33, "p")), dir.Write("u2_r1.fq", Fastq(20, "q")) });
    ThreadedGzIstream uneven2(std::vector<std::string>{ dir.Write("u1_r2.fq", Fastq(32, "p")), dir.Write("u2_r2.fq", Fastq(21, "q")) });
    testing::internal::CaptureStderr();
    SeqReaderPE reader{ uneven1, uneven2 };
    FastxRecord record1, record2;
    size_t pairs = 0;
    while (reader(record1, record2)) pairs++;
    testing::internal::GetCapturedStderr();
    EXPECT_FALSE(reader.Success());
    EXPECT_TRUE(reader.MatesOutOfStep());
    EXPECT_LT(pairs, 53u);
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
    ThreadedGzIstream plain(dir.Write("reads.fq", content).c_str()), gz(Gzip(dir, "reads.fq.gz", content).c_str());
    ThreadedGzIstream zst(dir.Write("reads.fq.zst", ZstdBytes(content)).c_str());
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
    ThreadedGzIstream gz(Gzip(dir, "bad.fq.gz", content).c_str());
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
    constexpr int kPeriod = Benchmark::kPerRead;
    Benchmark bm{"per read", 0, Benchmark::kPerRead};
    std::vector<int> timed;
    for (int i = 0; i < 6 * kPeriod + 7; i++) {
        bm.Start();
        if (bm.Timing()) timed.push_back(i);
        bm.Stop();
    }
    EXPECT_EQ(timed, (std::vector<int>{ 0, kPeriod, 2 * kPeriod, 3 * kPeriod, 4 * kPeriod, 5 * kPeriod, 6 * kPeriod }));
    EXPECT_EQ(bm.TimedCalls(), timed.size());
    EXPECT_FALSE(bm.Timing());
}

TEST(Benchmark, ASampledStageWithAnOddPeriodTimesBothMatesOfAPair) {
    // A stage called once per mate: an even period would time the same mate every time.
    ASSERT_EQ(Benchmark::kPerRead % 2, 1u);
    Benchmark bm{"per mate", 0, Benchmark::kPerRead};
    size_t first_mate = 0, second_mate = 0;
    for (uint32_t i = 0; i < 100 * Benchmark::kPerRead; i++) {
        bm.Start();
        if (bm.Timing()) (i % 2 ? second_mate : first_mate)++;
        bm.Stop();
    }
    EXPECT_EQ(first_mate + second_mate, 100u);
    EXPECT_NEAR(static_cast<double>(first_mate), 50.0, 1.0);
}

TEST(Benchmark, ASampledStageEstimatesTheTimeOfAllItsCalls) {
    Benchmark bm{"per read", 0, Benchmark::kPerRead};
    uint32_t const calls = 20 * Benchmark::kPerRead;  // 20 timed calls
    for (uint32_t i = 0; i < calls; i++) {
        bm.Start();
        Spin(std::chrono::microseconds(20));
        bm.Stop();
    }
    EXPECT_GE(bm.Seconds(), 0.99 * calls * 20e-6);
    EXPECT_EQ(bm.Threads(), 1u);
}

TEST(Benchmark, ASampledStageScalesTheFirstCallToTheCallsMade) {
    Benchmark bm{"per read", 0, Benchmark::kPerRead};
    bm.Start();
    Spin(std::chrono::microseconds(2000));
    bm.Stop();
    EXPECT_GE(bm.Seconds(), 0.002);  // a short run has a value from its first call
    for (uint32_t i = 1; i < Benchmark::kPerRead; i++) {  // the other calls of the period, none timed: nothing is added,
        bm.Start();                                       // the mean times their number is
        EXPECT_FALSE(bm.Timing());
        bm.Stop();
    }
    EXPECT_EQ(bm.TimedCalls(), 1u);
    EXPECT_GE(bm.Seconds(), Benchmark::kPerRead * 0.002);
}

TEST(Benchmark, ASampledStageJoinsAndPrintsItsEstimate) {
    Benchmark global{"stage", 0};
    uint32_t const calls = 2 * Benchmark::kPerRead;  // 2 timed calls
    for (int t = 0; t < 3; t++) {
        Benchmark local{"stage", 0, Benchmark::kPerRead};
        for (uint32_t i = 0; i < calls; i++) {
            local.Start();
            Spin(std::chrono::microseconds(100));
            local.Stop();
        }
        global.Join(local);
    }
    EXPECT_EQ(global.Threads(), 3u);
    EXPECT_GE(global.Seconds(), 3 * calls * 100e-6 * 0.99);
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

TEST(SeqReader, ReadsWithoutQualitiesGetTheReadTypesQuality) {
    std::istringstream fasta(">r1\nACGT\n>r2\nAC\nGT\n");
    SeqReaderSE reader(fasta, FastaQualityChar(ReadType::PacBio));
    FastxRecord record;
    ASSERT_TRUE(reader(record));
    EXPECT_EQ(record.quality, "????");  // Q30
    ASSERT_TRUE(reader(record));
    EXPECT_EQ(record.sequence, "ACGT");
    EXPECT_EQ(record.quality, "????");
    EXPECT_FALSE(reader(record));

    std::istringstream fastq("@r1\nACGT\n+\nIII#\n");
    SeqReaderSE keeps(fastq, FastaQualityChar(ReadType::PacBio));
    ASSERT_TRUE(keeps(record));
    EXPECT_EQ(record.quality, "III#") << "FASTQ keeps its qualities";

    std::istringstream mate1(">p/1\nAAAA\n"), mate2(">p/2\nCCC\n");
    SeqReaderPE pairs(mate1, mate2, FastaQualityChar(ReadType::Paired));
    FastxRecord record2;
    ASSERT_TRUE(pairs(record, record2));
    EXPECT_EQ(record.quality, "????");
    EXPECT_EQ(record2.quality, "???");
}
