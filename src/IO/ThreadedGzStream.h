#pragma once

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zstd.h>

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <functional>
#include <istream>
#include <memory>
#include <mutex>
#include <new>
#include <streambuf>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "Bgzf.h"
#include "ParallelGzip.h"

namespace protal {

    // Reads a gzip-compressed (or plain) file as igzstream does, but inflates it in a thread of its
    // own, up to kBlocks blocks ahead of the reader. protal's alignment threads take reads in
    // batches under one lock (SeqReaderPE, SeqReaderSE); with igzstream the lock holder also
    // inflated the input, which capped a run at one core's inflate speed (~200k read pairs/s)
    // whatever -t was. Here the lock holder only copies bytes that are already inflated, and the
    // two files of a pair inflate in parallel. A BGZF file (Bgzf.h: bgzip's, protal's) is
    // inflated block by block, any other gzip file streamed, both with ISA-L (igzip, Bgzf.h); a
    // single-member .fq.gz inflates at about 1 GB/s of FASTQ on one thread, 1.6-1.8x zlib-ng's speed,
    // which inflated it before and capped paired-end runs at 32 threads near 1.8M pairs/s
    // (docs/claude/2026-10-06-performance-profiling). A zstd file (its frames, one after the other,
    // skippable ones too, any window zstd --long writes) is read with libzstd's streaming
    // decompression. The format is told from the first bytes, not from the name; a file of none of
    // them is read as it is.
    //
    // As with igzstream, a truncated or corrupt file reads as one that ends early; read_failed()
    // tells the two apart once reading has stopped. A BGZF file must end with its end-of-file
    // block, so a cut at a block boundary is found too (it may go on with gzip members of other
    // kinds, as cat a.bgzf.gz b.gz writes); in other gzip files, what follows a member must be
    // another member (or zero bytes of padding), so a damaged member header is an error rather
    // than the end; a zstd file must not end inside a frame, and holds zstd frames only. The file
    // is opened once and read from that descriptor only, so pipes and process substitution
    // (<(zcat ...)) work; only a regular file is checked for BGZF.
    // The buffer is read by one thread at a time (the reader lock serialises it); its inflating
    // thread is stopped by close().
    //
    // A BGZF file's blocks are independent deflate streams, so they can inflate on several threads
    // (SetInflateThreads, protal: InflateThreadsFor its -t): one thread then reads the compressed blocks in
    // order and hands them out in batches of up to 1 MB of output, worker threads inflate the
    // batches with ISA-L, and the reader takes them in file order, as if one thread had inflated
    // them: the same bytes, and an error (a corrupt block, a cut file) where the file has it, with
    // everything before it read. One gzip member cannot be split this way (its deflate blocks are
    // found only by decoding it from the start, and refer back up to 32 KB): another gzip file, if
    // it is a regular file and gets more than one thread, is inflated by rapidgzip on that many
    // threads of its own (ParallelGzip.h), and this file's thread copies its output into the blocks;
    // with one thread, from a pipe (which rapidgzip would have to buffer), or ending in zero bytes
    // (padding, which rapidgzip fails on), ISA-L inflates it here on one thread. The bytes are the
    // same either way; a cut file, a corrupt member, a wrong CRC-32 or length, a damaged member header
    // or data after the last member fail either way, in other words (what was read before may differ:
    // a corrupt member's bytes up to where the damage shows; ParallelGzip checks what rapidgzip does not).
    //
    // Several files (open with a list: a sample's lanes) read as one stream, one after the other, each in
    // its own format: a file whose last line lacks its '\n' gets one, TakeLines ends a batch with its file
    // (so that the mate files of lanes stay in step, see TakeLines), and a truncated or corrupt file, or one
    // that cannot be opened when its turn comes, ends the stream there, named in read_error_message().
    class ThreadedGzStreambuf : public std::streambuf {
    public:
        // Called as each file of the stream is opened (the first by open, the others as the stream reaches them).
        using OnOpen = std::function<void(std::string const& path, ThreadedGzStreambuf const& buf)>;

        static constexpr size_t kBlockSize = size_t{1} << 20;  // inflated bytes per block
        static constexpr size_t kBlocks = 4;                   // with one inflating thread; 2 per worker + 2 with more
        static constexpr uint32_t kInputBuffer = 1u << 17;     // compressed bytes per read()
        static constexpr unsigned kMaxInflateThreads = 16;
        // The fewest threads rapidgzip gets for a gzip file that is not BGZF: on FASTQ it inflates ~0.3 GB/s a
        // thread, as its first-pass decoder does nearly all of it (a FASTQ file's copies of names and qualities
        // keep pointing back to the unknown 32 KB before a piece, so ISA-L seldom takes over), against ISA-L's
        // 1.2-1.6 GB/s on one (docs/claude/2026-10-12-rapidgzip-input): fewer than ~6 are slower than one.
        static constexpr unsigned kMinParallelGzipThreads = 8;

        // Threads that inflate a BGZF file opened from now on (1: the file's own thread does it; at most
        // kMaxInflateThreads). PROTAL_INFLATE_THREADS in the environment overrides it, and GzipThreads (for
        // comparing runs).
        static void SetInflateThreads(unsigned threads) {
            InflateThreadsSetting().store(std::clamp(threads, 1u, kMaxInflateThreads), std::memory_order_relaxed);
        }
        static unsigned InflateThreads() {
            if (unsigned const n = InflateThreadsOfTheEnvironment()) return n;
            return InflateThreadsSetting().load(std::memory_order_relaxed);
        }
        // protal's choice for a run on `threads` threads reading `files` files at once: a worker per 4
        // alignment threads and file, 1 to 4 per file (a paired-end run on 32 threads: 4 per file, ~4 GB/s
        // of FASTQ each; one ISA-L thread inflates ~1 GB/s).
        static unsigned InflateThreadsFor(size_t threads, size_t files) {
            size_t const per_file = threads / (4 * std::max<size_t>(files, 1));
            return static_cast<unsigned>(std::clamp<size_t>(per_file, 1, 4));
        }

        // Threads that inflate a gzip file that is not BGZF (a regular file) opened from now on: rapidgzip's, with
        // more than one (ParallelGzip); 1: ISA-L on the file's own thread.
        static void SetGzipThreads(unsigned threads) {
            GzipThreadsSetting().store(std::clamp(threads, 1u, kMaxInflateThreads), std::memory_order_relaxed);
        }
        static unsigned GzipThreads() {
            if (unsigned const n = InflateThreadsOfTheEnvironment()) return n;
            return GzipThreadsSetting().load(std::memory_order_relaxed);
        }
        // protal's choice: a thread per 2 alignment threads and file, if that is kMinParallelGzipThreads or more (a
        // paired-end run on 32 threads or more, one file on 16 or more), at most kMaxInflateThreads; else 1. The
        // threads come on top of -t; they wait for room, not spin, while the alignment takes the reads slower than
        // they inflate them.
        static unsigned GzipThreadsFor(size_t threads, size_t files) {
            size_t const per_file = threads / (2 * std::max<size_t>(files, 1));
            if (per_file < kMinParallelGzipThreads) return 1;
            return static_cast<unsigned>(std::min<size_t>(per_file, kMaxInflateThreads));
        }

        ThreadedGzStreambuf() = default;
        ThreadedGzStreambuf(ThreadedGzStreambuf const&) = delete;
        ThreadedGzStreambuf& operator=(ThreadedGzStreambuf const&) = delete;
        ~ThreadedGzStreambuf() override { close(); }

        // False if the file cannot be read; open_error() then says why.
        bool open(char const* path) {
            return open(std::vector<std::string>{ path });
        }

        // Several files, read one after the other as one stream; each is opened once the one before it has been read.
        // False if the first cannot be read (open_error() says why); a later one that cannot be read ends the stream
        // there (read_failed()).
        bool open(std::vector<std::string> paths) {
            if (is_open() || paths.empty()) return false;
            m_paths = std::move(paths);
            m_path = 0;
            m_files_begun = 0;
            m_last_byte = '\n';
            return OpenFile(m_paths.front().c_str());
        }

        // set_on_open(f): f is called as each file is opened from the next open() on.
        void set_on_open(OnOpen on_open) { m_on_open = std::move(on_open); }

        // The file being read (the last one opened, or the one that could not be opened).
        std::string const& current_path() const {
            static std::string const none;
            return m_path < m_paths.size() ? m_paths[m_path] : none;
        }

        // The files the stream reads (one, unless opened with a list).
        size_t files() const { return m_paths.size(); }

    private:
        bool OpenFile(char const* path) {
            if (is_open()) return false;
            m_open_error.clear();
            m_fd = ::open(path, O_RDONLY | O_CLOEXEC);
            if (m_fd < 0) {
                m_open_error = std::strerror(errno);
                return false;
            }
            struct stat st {};
            if (::fstat(m_fd, &st) != 0 || S_ISDIR(st.st_mode)) {
                m_open_error = S_ISDIR(st.st_mode) ? "it is a directory" : std::strerror(errno);
                ::close(m_fd);
                m_fd = -1;
                return false;
            }
            // pread leaves the offset where it is; a pipe cannot be peeked at, and is read as gzip.
            unsigned char head[bgzf::kHeaderBytes] = {};
            ssize_t const peeked = S_ISREG(st.st_mode) ? ::pread(m_fd, head, sizeof(head), 0) : 0;
            m_bgzf = peeked == static_cast<ssize_t>(sizeof(head)) && bgzf::IsBlockHeader(head);
            m_offset = 0;
            m_members = 0;
            m_parallel.reset();
            m_parallel_error.clear();
            unsigned const threads = m_bgzf ? InflateThreads() : GzipThreads();
            if (m_bgzf) {
                m_block_in.resize(bgzf::kMaxBlock);
                m_last_empty = false;
            } else if (peeked >= 2 && head[0] == 0x1f && head[1] == 0x8b && threads > 1) {
                // A gzip file that is not BGZF, on several threads: rapidgzip's (ParallelGzip); on one if it cannot
                // be had (m_parallel_error says why), or if the file ends in four zero bytes: rapidgzip fails on zero
                // padding after the last member, which gzip and FillStream skip. (An empty last member ends so too,
                // and one of 4 GiB times n bytes; read on one thread, those are read all the same.)
                unsigned char tail[4] = { 1, 1, 1, 1 };
                bool const zero_tail = st.st_size >= 4 && ::pread(m_fd, tail, sizeof(tail), st.st_size - 4) == 4 &&
                                       tail[0] == 0 && tail[1] == 0 && tail[2] == 0 && tail[3] == 0;
                if (zero_tail) m_parallel_error = "the file ends in zero bytes (padding?), which it does not read";
                else m_parallel = ParallelGzip::Open(m_fd, threads, m_parallel_error);
            }
            if (!m_bgzf && !m_parallel && !StartGzip(0, m_open_error)) {
                ::close(m_fd);
                m_fd = -1;
                return false;
            }
            m_opened_bgzf = m_bgzf;
            m_workers = m_bgzf ? threads : 1;
            m_parallel_threads = m_parallel ? threads : 0;
            size_t const blocks = m_workers > 1 ? 2 * size_t{m_workers} + 2 : kBlocks;
            m_blocks.assign(blocks, std::vector<char>(kBlockSize));
            m_jobs.assign(blocks, Job{});
            m_free.clear();
            m_ready.clear();
            m_todo.clear();
            for (size_t i = 0; i < blocks; i++) m_free.push_back(i);
            m_current = kNone;
            m_stop = false;
            m_done = false;
            m_failed = false;
            m_error.clear();
            m_dispatch_error.clear();
            setg(nullptr, nullptr, nullptr);
            if (m_workers > 1) {
                m_thread = std::thread(&ThreadedGzStreambuf::Dispatch, this);
                for (unsigned w = 0; w < m_workers; w++) m_worker_threads.emplace_back(&ThreadedGzStreambuf::Work, this);
            } else {
                m_thread = std::thread(&ThreadedGzStreambuf::Inflate, this);
            }
            m_files_begun++;
            m_last_byte = '\n';
            if (m_on_open) m_on_open(path, *this);
            return true;
        }

        // After the end of a file read whole, the next file of the stream: false after the last, after a file that ended
        // in an error (nothing more is read), or if the next cannot be opened (an error, as for a damaged file).
        bool NextFile() {
            {
                std::lock_guard<std::mutex> lock(m_mutex);
                if (!m_error.empty() || m_failed) return false;
            }
            if (m_path + 1 >= m_paths.size()) return false;
            close();
            m_path++;
            if (OpenFile(m_paths[m_path].c_str())) return true;
            std::lock_guard<std::mutex> lock(m_mutex);
            m_error = "it cannot be read (" + m_open_error + ")";
            return false;
        }

    public:
        // The threads inflating the open file: BGZF blocks, or rapidgzip's for other gzip (1: the file's own thread).
        unsigned inflate_threads() const { return m_parallel_threads > 0 ? m_parallel_threads : m_workers; }

        // Whether the open file began as BGZF (else gzip, zstd or uncompressed).
        bool opened_as_bgzf() const { return m_opened_bgzf; }

        // Whether the open file, gzip but not BGZF, is inflated by rapidgzip on several threads (ParallelGzip).
        bool parallel_gzip() const { return m_parallel_threads > 0; }

        // Why a gzip file that would have been inflated on several threads is inflated on one (empty if none).
        std::string const& parallel_gzip_error() const { return m_parallel_error; }

        bool is_open() const { return m_fd >= 0; }

        // Why the last open() failed.
        std::string const& open_error() const { return m_open_error; }

        void close() {
            if (!is_open()) return;
            {
                std::lock_guard<std::mutex> lock(m_mutex);
                m_stop = true;
            }
            m_cv.notify_all();
            if (m_thread.joinable()) m_thread.join();
            for (auto& worker : m_worker_threads) worker.join();
            m_worker_threads.clear();
            m_jobs.clear();
            m_parallel.reset();  // its threads end; its copy of the descriptor is closed
            m_parallel_threads = 0;
            m_inflate.reset();
            if (m_zstd) ZSTD_freeDCtx(m_zstd);
            m_zstd = nullptr;
            ::close(m_fd);
            m_fd = -1;
            setg(nullptr, nullptr, nullptr);
            m_blocks.clear();
        }

        // A truncated or corrupt file reads as one that ends early; this tells them apart.
        bool read_failed() const {
            std::lock_guard<std::mutex> lock(m_mutex);
            return !m_error.empty();
        }

        // Why reading failed; with several files, after the name of the file it failed in.
        std::string read_error_message() const {
            std::lock_guard<std::mutex> lock(m_mutex);
            if (m_error.empty() || m_paths.size() < 2) return m_error;
            return current_path() + ": " + m_error;
        }

        // Appends up to `lines` whole lines to out, each with its '\n' (a last line without one gets
        // it, as getline reads such a line too), and returns how many. The lines are found with
        // memchr in the inflated block and copied at once, so the reader lock is held for little
        // more than a copy (BufferedFastxReader::LoadBatch).
        // Of several files, a batch ends with its file once it has a line: the batches of the mate files
        // of a sample's lanes end with the same pair if the lanes' mates match, and otherwise a lane's
        // extra reads fail the sample (SeqReaderPE) instead of pairing the next lane's reads wrongly.
        size_t TakeLines(size_t lines, std::string& out) {
            size_t taken = 0;
            bool partial = false;  // out ends inside a line
            bool ending = false;   // the file ended inside out's last line: the batch ends with that line
            while (taken < lines) {
                if (gptr() == egptr()) {
                    size_t const files = m_files_begun;
                    if (traits_type::eq_int_type(underflow(), traits_type::eof())) break;
                    if (m_files_begun != files) {  // the next file
                        if (partial) ending = true;    // its line gets the '\n' the last file lacked
                        else if (taken > 0) break;
                    }
                }
                char* const begin = gptr();
                char* const end = egptr();
                char* p = begin;
                while (taken < lines) {
                    auto* const newline = static_cast<char*>(std::memchr(p, '\n', static_cast<size_t>(end - p)));
                    if (!newline) {
                        p = end;
                        break;
                    }
                    p = newline + 1;
                    taken++;
                }
                out.append(begin, static_cast<size_t>(p - begin));
                partial = p[-1] != '\n';
                gbump(static_cast<int>(p - begin));
                if (ending && !partial) break;
            }
            if (partial && taken < lines) {
                out.push_back('\n');
                taken++;
            }
            return taken;
        }

    protected:
        // The current file's next block; at its end the next file's (NextFile), after a '\n' if the file's
        // last line lacked one.
        int_type underflow() override {
            if (gptr() < egptr()) return traits_type::to_int_type(*gptr());
            if (egptr() > eback()) m_last_byte = egptr()[-1];
            while (is_open()) {
                int_type const next = UnderflowFile();
                if (!traits_type::eq_int_type(next, traits_type::eof())) return next;
                bool const open_line = m_last_byte != '\n';
                if (!NextFile()) break;
                if (open_line) {
                    setg(m_newline, m_newline, m_newline + 1);
                    return traits_type::to_int_type('\n');
                }
            }
            return traits_type::eof();
        }

    private:
        int_type UnderflowFile() {
            if (!is_open()) return traits_type::eof();
            std::unique_lock<std::mutex> lock(m_mutex);
            if (m_current != kNone) {
                m_free.push_back(m_current);
                m_current = kNone;
                m_cv.notify_all();
            }
            // The blocks in file order: the next one once it is inflated (a worker may finish a later one first).
            while (!m_failed) {
                m_cv.wait(lock, [this] { return (!m_ready.empty() && m_jobs[m_ready.front()].complete) || (m_ready.empty() && m_done); });
                if (m_ready.empty()) {
                    // The end, or what stopped the dispatcher, after everything before it was read.
                    if (!m_dispatch_error.empty() && m_error.empty()) m_error = m_dispatch_error;
                    break;
                }
                size_t const index = m_ready.front();
                m_ready.pop_front();
                Job const& job = m_jobs[index];
                if (!job.error.empty()) {
                    // A corrupt block: what precedes it in the batch is read, then the file ends here.
                    if (m_error.empty()) m_error = job.error;
                    m_failed = true;
                    m_stop = true;
                    m_cv.notify_all();
                }
                if (job.size == 0) {
                    m_free.push_back(index);
                    m_cv.notify_all();
                    continue;
                }
                m_current = index;
                char* data = m_blocks[index].data();
                setg(data, data, data + job.size);
                return traits_type::to_int_type(*gptr());
            }
            setg(nullptr, nullptr, nullptr);
            return traits_type::eof();
        }

    private:
        static constexpr size_t kNone = SIZE_MAX;

        // An output block's batch: with several inflating threads, the compressed BGZF blocks to inflate into it.
        struct Job {
            std::vector<unsigned char> in;               // whole BGZF blocks, one after the other
            std::vector<size_t> ends;                    // where each ends in `in`
            uint64_t offset = 0;                         // the file offset of the first
            size_t size = 0;                             // inflated bytes in the output block
            bool complete = false;                       // inflated (or filled): the reader may take it
            std::string error;                           // the block that could not be inflated, if one
        };


        void Inflate() {
            bool end = false;
            while (!end) {
                size_t index;
                {
                    std::unique_lock<std::mutex> lock(m_mutex);
                    m_cv.wait(lock, [this] { return !m_free.empty() || m_stop; });
                    if (m_stop) break;
                    index = m_free.front();
                    m_free.pop_front();
                }
                char* data = m_blocks[index].data();
                size_t size = 0;
                std::string error;
                end = m_parallel ? FillParallel(data, size, error)
                    : m_bgzf     ? FillBgzf(data, size, error)
                                 : FillStream(data, size, error);
                {
                    std::lock_guard<std::mutex> lock(m_mutex);
                    m_jobs[index].size = size;
                    m_jobs[index].error.clear();
                    m_jobs[index].complete = true;
                    if (size > 0) m_ready.push_back(index);
                    else m_free.push_back(index);
                    if (!error.empty() && m_error.empty()) m_error = error;
                }
                m_cv.notify_all();
            }
            {
                std::lock_guard<std::mutex> lock(m_mutex);
                m_done = true;
            }
            m_cv.notify_all();
        }

        // With several inflating threads (a BGZF file): reads the blocks in order into batches (ReadBgzfJob), one per
        // free output block, and queues them for the workers (Work) and, in the same order, for the reader. If a gzip
        // member of another kind follows the BGZF blocks, it is inflated here, as Inflate does, into the next blocks.
        // What stops the reading (a cut file, a damaged header, a missing end-of-file block) is the reader's error once
        // it has read every block before it.
        void Dispatch() {
            bool end = false;
            while (!end) {
                size_t index;
                {
                    std::unique_lock<std::mutex> lock(m_mutex);
                    m_cv.wait(lock, [this] { return !m_free.empty() || m_stop; });
                    if (m_stop) break;
                    index = m_free.front();
                    m_free.pop_front();
                }
                Job& job = m_jobs[index];  // this thread's until it is queued
                job.complete = false;
                job.size = 0;
                job.error.clear();
                std::string error;
                bool const bgzf = m_bgzf;
                if (bgzf) {
                    end = ReadBgzfJob(job, error);
                } else {
                    size_t size = 0;
                    end = FillStream(m_blocks[index].data(), size, error);
                    job.size = size;
                    job.complete = true;
                }
                {
                    std::lock_guard<std::mutex> lock(m_mutex);
                    if (bgzf && !job.ends.empty()) {
                        m_ready.push_back(index);
                        m_todo.push_back(index);
                    } else if (!bgzf && job.size > 0) {
                        m_ready.push_back(index);
                    } else {
                        m_free.push_back(index);
                    }
                    if (!error.empty() && m_dispatch_error.empty()) m_dispatch_error = error;
                }
                m_cv.notify_all();
            }
            {
                std::lock_guard<std::mutex> lock(m_mutex);
                m_done = true;
            }
            m_cv.notify_all();
        }

        // A worker: inflates queued batches block by block into their output blocks (bgzf::DecompressBlock, ISA-L, a
        // state per thread), until the dispatcher is done and nothing is queued, or the file is closed. A corrupt
        // block ends its batch there, with the error and the byte where the block starts.
        void Work() {
            while (true) {
                size_t index;
                {
                    std::unique_lock<std::mutex> lock(m_mutex);
                    m_cv.wait(lock, [this] { return !m_todo.empty() || m_done || m_stop; });
                    if (m_stop || m_todo.empty()) return;
                    index = m_todo.front();
                    m_todo.pop_front();
                }
                Job& job = m_jobs[index];
                char* out = m_blocks[index].data();
                size_t size = 0, begin = 0;
                uint64_t offset = job.offset;
                std::string error;
                for (size_t const end : job.ends) {
                    size_t n = 0;
                    if (!bgzf::DecompressBlock(job.in.data() + begin, end - begin, out + size, kBlockSize - size, n, error)) {
                        error += " at byte " + std::to_string(offset);
                        break;
                    }
                    size += n;
                    offset += end - begin;
                    begin = end;
                }
                {
                    std::lock_guard<std::mutex> lock(m_mutex);
                    job.size = size;
                    job.error = std::move(error);
                    job.complete = true;
                }
                m_cv.notify_all();
            }
        }

        // FillBgzf without the inflating: reads whole BGZF blocks into job (job.in, each block's end in job.ends)
        // while their contents, as their footers give them, leave room for one more block of 64 KB in an output block;
        // true at the end of the file or on an error (then set). A gzip member of another kind after the blocks
        // switches the file to FillStream (StartGzip), and false is returned.
        bool ReadBgzfJob(Job& job, std::string& error) {
            job.in.clear();
            job.ends.clear();
            job.offset = m_offset;
            size_t out = 0;
            while (out + bgzf::kMaxBlock <= kBlockSize) {
                unsigned char* head = m_block_in.data();
                size_t const got = ReadUpTo(head, bgzf::kHeaderBytes, error);
                if (!error.empty()) return true;
                if (got == 0) {
                    if (!m_last_empty) error = "the BGZF end-of-file block is missing (truncated file?)";
                    return true;
                }
                if (got < bgzf::kHeaderBytes || !bgzf::IsBlockHeader(head)) {
                    if (got >= 2 && head[0] == 0x1f && head[1] == 0x8b) return !StartGzip(got, error);
                    error = got < bgzf::kHeaderBytes ? "the file ends inside a BGZF block header (truncated file?)"
                                                     : "no BGZF block at byte " + std::to_string(m_offset) +
                                                       " (a damaged block header?)";
                    return true;
                }
                size_t const block_size = bgzf::BlockSize(head);
                if (block_size < bgzf::kHeaderBytes + bgzf::kFooterBytes) {
                    error = "an invalid BGZF block size at byte " + std::to_string(m_offset);
                    return true;
                }
                size_t const start = job.in.size();
                job.in.resize(start + block_size);
                std::memcpy(job.in.data() + start, head, bgzf::kHeaderBytes);
                size_t const rest = block_size - bgzf::kHeaderBytes;
                if (ReadUpTo(job.in.data() + start + bgzf::kHeaderBytes, rest, error) != rest) {
                    job.in.resize(start);
                    if (error.empty()) error = "the file ends inside a BGZF block (truncated file?)";
                    return true;
                }
                uint32_t const isize = bgzf::LE32(job.in.data() + start + block_size - 4);
                job.ends.push_back(job.in.size());
                m_last_empty = isize == 0;
                m_members++;
                m_offset += block_size;
                out += isize;
            }
            return false;
        }

        // Reads up to size bytes at the file's current offset; the bytes read (fewer only at its end).
        size_t ReadUpTo(unsigned char* dst, size_t size, std::string& error) {
            size_t done = 0;
            while (done < size) {
                ssize_t const n = ::read(m_fd, dst + done, size - done);
                if (n < 0 && errno == EINTR) continue;
                if (n < 0) {
                    error = std::string("read error: ") + std::strerror(errno);
                    break;
                }
                if (n == 0) break;
                done += static_cast<size_t>(n);
            }
            return done;
        }

        // Inflates BGZF blocks into data (from size) while a whole block still fits; true at the end
        // of the file or on an error (then set).
        bool FillBgzf(char* data, size_t& size, std::string& error) {
            while (size + bgzf::kMaxBlock <= kBlockSize) {
                unsigned char* block = m_block_in.data();
                size_t const got = ReadUpTo(block, bgzf::kHeaderBytes, error);
                if (!error.empty()) return true;
                if (got == 0) {
                    if (!m_last_empty) error = "the BGZF end-of-file block is missing (truncated file?)";
                    return true;
                }
                if (got < bgzf::kHeaderBytes || !bgzf::IsBlockHeader(block)) {
                    // A gzip member of another kind after the BGZF blocks (cat a.bgzf.gz b.gz): on as gzip.
                    if (got >= 2 && block[0] == 0x1f && block[1] == 0x8b) {
                        return !StartGzip(got, error) || FillStream(data, size, error);
                    }
                    error = got < bgzf::kHeaderBytes ? "the file ends inside a BGZF block header (truncated file?)"
                                                     : "no BGZF block at byte " + std::to_string(m_offset) +
                                                       " (a damaged block header?)";
                    return true;
                }
                size_t const block_size = bgzf::BlockSize(block);
                if (block_size < bgzf::kHeaderBytes + bgzf::kFooterBytes) {
                    error = "an invalid BGZF block size at byte " + std::to_string(m_offset);
                    return true;
                }
                size_t const rest = block_size - bgzf::kHeaderBytes;
                if (ReadUpTo(block + bgzf::kHeaderBytes, rest, error) != rest) {
                    if (error.empty()) error = "the file ends inside a BGZF block (truncated file?)";
                    return true;
                }
                size_t n = 0;
                if (!bgzf::DecompressBlock(block, block_size, data + size, kBlockSize - size, n, error)) {
                    error += " at byte " + std::to_string(m_offset);
                    return true;
                }
                m_last_empty = n == 0;
                m_members++;
                size += n;
                m_offset += block_size;
            }
            return false;
        }

        // Copies the file's next inflated bytes from rapidgzip (m_parallel) into data (from size) until the block is
        // full; true at the end of the file or on an error (then set).
        bool FillParallel(char* data, size_t& size, std::string& error) {
            while (size < kBlockSize) {
                size_t const n = m_parallel->Read(data + size, kBlockSize - size);
                if (n == 0) {
                    if (!m_parallel->Error().empty()) error = "corrupt gzip data (corrupt or truncated file?): " + m_parallel->Error();
                    return true;
                }
                size += n;
            }
            return false;
        }

        // Reads the rest of the file with FillStream, from its first `got` bytes in m_block_in (read
        // there by FillBgzf, else none); false if the inflate state cannot be had (error then set).
        bool StartGzip(size_t got, std::string& error) {
            m_bgzf = false;
            m_block_in.resize(std::max<size_t>(m_block_in.size(), kInputBuffer));
            m_in_pos = 0;
            m_in_end = got;
            m_in_eof = got > 0 && got < bgzf::kHeaderBytes;  // FillBgzf read fewer only at the end
            m_mode = Mode::Look;
            try {
                m_inflate = std::make_unique<inflate_state>();  // ~85 KB
            } catch (std::bad_alloc const&) {
                error = "out of memory for the gzip inflate state";
                return false;
            }
            isal_inflate_init(m_inflate.get());
            return true;
        }

        // A gzip member begins at the unread input: its header is read next, then its deflate data and trailer,
        // which ISA-L checks (CRC-32 and length). The header is read with a state of its own (m_gz_header):
        // isal_inflate's own header reading (crc_flag ISAL_GZIP) fails on a header with a header CRC, or with two
        // of name, comment and extra field, that is split between two reads of the file (ISA-L 2.32.1).
        void StartMember() {
            isal_inflate_reset(m_inflate.get());
            m_inflate->crc_flag = ISAL_GZIP_NO_HDR_VER;
            isal_gzip_header_init(&m_gz_header);  // no buffers: the name, comment and extra field are skipped
            m_header_read = false;
        }

        // What an ISA-L error code says, as zlib words it.
        static char const* InflateError(int ret) {
            switch (ret) {
                case ISAL_INVALID_BLOCK: return "invalid block type";
                case ISAL_INVALID_SYMBOL: return "invalid code";
                case ISAL_INVALID_LOOKBACK: return "invalid distance too far back";
                case ISAL_INVALID_WRAPPER: return "incorrect header check";
                case ISAL_UNSUPPORTED_METHOD: return "unknown compression method";
                case ISAL_INCORRECT_CHECKSUM: return "incorrect data or length check";
                default: return "corrupt gzip data";
            }
        }

        // A zstd frame's magic number, or a skippable frame's (pzstd and the seekable format write them).
        static bool IsZstdFrame(unsigned char const* in) {
            return (in[0] == 0x28 && in[1] == 0xb5 && in[2] == 0x2f && in[3] == 0xfd) ||
                   ((in[0] & 0xf0) == 0x50 && in[1] == 0x2a && in[2] == 0x4d && in[3] == 0x18);
        }

        // FillZstd's decompression context; false if it cannot be made (error then set).
        bool StartZstd(std::string& error) {
            m_zstd = ZSTD_createDCtx();
            if (!m_zstd) {
                error = "cannot allocate a zstd decompression context";
                return false;
            }
            // The largest windows (zstd --long up to 2 GB); the memory follows the frame's own.
            ZSTD_DCtx_setParameter(m_zstd, ZSTD_d_windowLogMax, ZSTD_dParam_getBounds(ZSTD_d_windowLogMax).upperBound);
            m_zstd_in_frame = false;
            return true;
        }

        // Makes at least `want` unread bytes available in m_block_in unless the file ends first; false
        // on a read error (then set). Unread bytes move to the front first.
        bool Available(size_t want, std::string& error) {
            if (m_in_end - m_in_pos >= want || m_in_eof) return true;
            std::memmove(m_block_in.data(), m_block_in.data() + m_in_pos, m_in_end - m_in_pos);
            m_in_end -= m_in_pos;
            m_in_pos = 0;
            while (m_in_end < want && !m_in_eof) {
                size_t const got = ReadUpTo(m_block_in.data() + m_in_end, m_block_in.size() - m_in_end, error);
                if (!error.empty()) return false;
                if (got == 0) m_in_eof = true;
                m_in_end += got;
            }
            return true;
        }

        // The rest of the file is zero bytes (padding, which gzip ignores too).
        bool RestIsZero(std::string& error) {
            while (true) {
                for (size_t i = m_in_pos; i < m_in_end; i++) {
                    if (m_block_in[i] != 0) return false;
                }
                m_offset += m_in_end - m_in_pos;
                m_in_pos = m_in_end;
                if (!Available(1, error) || m_in_end == m_in_pos) return error.empty();
            }
        }

        // Inflates gzip members (decompresses a zstd file, or copies a file that is neither) into data
        // (from size) until the block is full; true at the end of the file or on an error (then set).
        // After a gzip member, the file must end, hold zero padding, or go on with another member:
        // anything else is taken for a damaged member, not ignored as zlib's gzread does.
        bool FillStream(char* data, size_t& size, std::string& error) {
            while (size < kBlockSize) {
                if (m_mode == Mode::Look) {
                    if (!Available(4, error)) return true;
                    size_t const left = m_in_end - m_in_pos;
                    if (left == 0) return true;  // the end of the file
                    unsigned char const* in = m_block_in.data() + m_in_pos;
                    if (left >= 2 && in[0] == 0x1f && in[1] == 0x8b) {
                        StartMember();
                        m_members++;
                        m_mode = Mode::Gzip;
                    } else if (m_members == 0 && left >= 4 && IsZstdFrame(in)) {
                        if (!StartZstd(error)) return true;
                        m_mode = Mode::Zstd;
                    } else if (m_members == 0) {
                        m_mode = Mode::Copy;  // neither gzip nor zstd: read as it is
                    } else {
                        uint64_t const at = m_offset;
                        if (!RestIsZero(error) && error.empty()) {
                            error = "the data at byte " + std::to_string(at) + ", after gzip member " + std::to_string(m_members) +
                                    ", is no gzip member (a damaged member header?)";
                        }
                        return true;
                    }
                }
                if (m_mode == Mode::Zstd) return FillZstd(data, size, error);
                if (m_mode == Mode::Copy) {
                    if (!Available(1, error)) return true;
                    size_t const n = std::min(m_in_end - m_in_pos, kBlockSize - size);
                    if (n == 0) return true;
                    std::memcpy(data + size, m_block_in.data() + m_in_pos, n);
                    m_in_pos += n;
                    m_offset += n;
                    size += n;
                    continue;
                }
                // Mode::Gzip: the member's header, then its deflate data and trailer. At the member's end ISA-L
                // gives back the bytes it read past the trailer, so the unread input starts where the member ends.
                if (!Available(1, error)) return true;
                if (m_in_end == m_in_pos) {
                    error = "the file ends inside gzip member " + std::to_string(m_members) + " (truncated file?)";
                    return true;
                }
                inflate_state& state = *m_inflate;
                size_t const had = m_in_end - m_in_pos;
                size_t const size_before = size;
                state.next_in = m_block_in.data() + m_in_pos;
                state.avail_in = static_cast<uint32_t>(had);
                int ret;
                if (!m_header_read) {
                    ret = isal_read_gzip_header(&state, &m_gz_header);
                    m_header_read = ret == ISAL_DECOMP_OK;
                    if (ret > 0) ret = ISAL_DECOMP_OK;  // ISAL_END_INPUT: more of the header to come
                } else {
                    state.next_out = reinterpret_cast<uint8_t*>(data + size);
                    state.avail_out = static_cast<uint32_t>(kBlockSize - size);
                    ret = isal_inflate(&state);
                    size = kBlockSize - state.avail_out;
                }
                size_t const used = had - state.avail_in;
                m_in_pos += used;
                m_offset += used;
                if (ret != ISAL_DECOMP_OK) {
                    error = std::string(InflateError(ret)) + " in gzip member " + std::to_string(m_members) + " (corrupt file?)";
                    return true;
                }
                if (m_header_read && state.block_state == ISAL_BLOCK_FINISH) {
                    m_mode = Mode::Look;
                } else if (used == 0 && size == size_before && size < kBlockSize) {
                    // Nothing consumed or made although there is room and input: more input is needed.
                    std::string more;
                    if (!Available(had + 1, more)) { error = more; return true; }
                    if (m_in_end - m_in_pos == had) {
                        error = "the file ends inside gzip member " + std::to_string(m_members) + " (truncated file?)";
                        return true;
                    }
                }
            }
            return false;
        }

        // Decompresses zstd frames into data (from size) until the block is full; true at the end of the
        // file or on an error (then set). ZSTD_decompressStream reads frame after frame, skippable ones
        // too, and returns 0 where one ends: the file must end there, not inside a frame.
        bool FillZstd(char* data, size_t& size, std::string& error) {
            while (size < kBlockSize) {
                if (!Available(1, error)) return true;
                size_t const had = m_in_end - m_in_pos;
                if (had == 0 && !m_zstd_in_frame) return true;  // the end of the file, after a whole frame
                ZSTD_inBuffer in {m_block_in.data() + m_in_pos, had, 0};
                ZSTD_outBuffer out {data + size, kBlockSize - size, 0};
                size_t const ret = ZSTD_decompressStream(m_zstd, &out, &in);
                m_in_pos += in.pos;
                m_offset += in.pos;
                size += out.pos;
                if (ZSTD_isError(ret)) {
                    error = std::string(ZSTD_getErrorName(ret)) + " in zstd frame " + std::to_string(m_members + 1) +
                            " (corrupt file?)";
                    return true;
                }
                if (in.pos == 0 && out.pos == 0) {
                    // No progress: the frame needs more input than there is unread (or, with none, it is cut).
                    std::string more;
                    if (!Available(had + 1, more)) { error = more; return true; }
                    if (m_in_end - m_in_pos == had) {
                        error = "the file ends inside zstd frame " + std::to_string(m_members + 1) + " (truncated file?)";
                        return true;
                    }
                    continue;
                }
                m_zstd_in_frame = ret != 0;
                if (ret == 0) m_members++;  // a frame ended
            }
            return false;
        }

        enum class Mode { Look, Gzip, Zstd, Copy };
        std::vector<std::string> m_paths;                // the stream's files, read one after the other
        size_t m_path = 0;                               // the one open (or that could not be opened)
        size_t m_files_begun = 0;                        // files opened so far (TakeLines ends a batch at a new one)
        char m_last_byte = '\n';                         // the last byte the reader took of the open file
        char m_newline[1] = { '\n' };                    // the get area of the '\n' a file's open last line gets
        OnOpen m_on_open;
        bool m_bgzf = false;                             // read with FillBgzf, else FillStream
        int m_fd = -1;
        std::string m_open_error;
        std::vector<unsigned char> m_block_in;           // one compressed BGZF block, or FillStream's input
        uint64_t m_offset = 0;                           // of the next unread byte of the file
        size_t m_in_pos = 0, m_in_end = 0;               // FillStream: unread input in m_block_in
        bool m_in_eof = false;                           // FillStream: the file has no more bytes
        size_t m_members = 0;                            // gzip members begun (BGZF blocks read; zstd frames ended)
        Mode m_mode = Mode::Look;
        std::unique_ptr<inflate_state> m_inflate;        // a gzip member's inflate state (ISA-L, ~85 KB)
        isal_gzip_header m_gz_header {};                 // its header's, while it is read
        bool m_header_read = false;                      // the member's header has been read
        ZSTD_DCtx* m_zstd = nullptr;                     // a zstd file's decompression
        bool m_zstd_in_frame = false;                    // a zstd frame has begun and not ended
        bool m_last_empty = false;                       // the last block read was empty (the EOF marker)
        static std::atomic<unsigned>& InflateThreadsSetting() {
            static std::atomic<unsigned> threads{ 1 };
            return threads;
        }
        static std::atomic<unsigned>& GzipThreadsSetting() {
            static std::atomic<unsigned> threads{ 1 };
            return threads;
        }
        // PROTAL_INFLATE_THREADS's value (at most kMaxInflateThreads), 0 if it is not set.
        static unsigned InflateThreadsOfTheEnvironment() {
            if (char const* env = std::getenv("PROTAL_INFLATE_THREADS")) {
                long const n = std::strtol(env, nullptr, 10);
                if (n >= 1) return static_cast<unsigned>(std::min<long>(n, kMaxInflateThreads));
            }
            return 0;
        }

        std::thread m_thread;                            // Inflate, or with several inflating threads Dispatch
        unsigned m_workers = 1;                          // threads inflating BGZF blocks (1: m_thread does)
        std::unique_ptr<ParallelGzip> m_parallel;        // a gzip file's inflating on rapidgzip's threads (Inflate copies)
        unsigned m_parallel_threads = 0;                 // its threads (0: none)
        std::string m_parallel_error;                    // why a gzip file is not inflated on several threads after all
        bool m_opened_bgzf = false;                      // the file began as BGZF
        std::vector<std::thread> m_worker_threads;       // Work, with several
        std::vector<std::vector<char>> m_blocks;
        std::vector<Job> m_jobs;                         // per output block
        std::deque<size_t> m_free;                       // blocks the inflating thread (dispatcher) may fill
        std::deque<size_t> m_ready;                      // blocks in file order; the reader takes the first once complete
        std::deque<size_t> m_todo;                       // batches for the workers
        size_t m_current = kNone;                        // the block the get area points into
        bool m_stop = false;                             // close() asks the inflating threads to end
        bool m_done = false;                             // no more blocks will be queued
        bool m_failed = false;                           // the reader reached a corrupt block: nothing more is read
        std::string m_dispatch_error;                    // what stopped the dispatcher, the reader's once it gets there
        std::string m_error;
        mutable std::mutex m_mutex;
        std::condition_variable m_cv;
    };

    // The istream for ThreadedGzStreambuf, a drop-in for igzstream when reading.
    class ThreadedGzIstream : public std::istream {
    public:
        explicit ThreadedGzIstream(char const* path) : ThreadedGzIstream(std::vector<std::string>{ path }) {}

        // Several files read one after the other as one stream (ThreadedGzStreambuf::open); on_open, if given, is
        // called as each is opened.
        explicit ThreadedGzIstream(std::vector<std::string> paths, ThreadedGzStreambuf::OnOpen on_open = {}) : std::istream(nullptr) {
            std::istream::rdbuf(&m_buf);  // clears the badbit of the null buffer
            m_buf.set_on_open(std::move(on_open));
            if (!m_buf.open(std::move(paths))) setstate(std::ios::badbit);
        }

        ThreadedGzStreambuf* rdbuf() { return &m_buf; }

        void close() { m_buf.close(); }

    private:
        ThreadedGzStreambuf m_buf;
    };
}
