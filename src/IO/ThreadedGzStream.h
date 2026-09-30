#pragma once

#include <fcntl.h>
#include <unistd.h>
#include <zlib.h>

#include <cerrno>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <deque>
#include <istream>
#include <mutex>
#include <streambuf>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "Bgzf.h"

namespace protal {

    // Reads a gzip-compressed (or plain) file as igzstream does, but inflates it in a thread of its
    // own, up to kBlocks blocks ahead of the reader. protal's alignment threads take reads in
    // batches under one lock (SeqReaderPE, SeqReaderSE); with igzstream the lock holder also
    // inflated the input, which capped a run at one core's inflate speed (~200k read pairs/s)
    // whatever -t was. Here the lock holder only copies bytes that are already inflated, and the
    // two files of a pair inflate in parallel. A BGZF file (Bgzf.h: bgzip's, protal's) is
    // inflated block by block with libdeflate, about 3x faster; any other file with zlib, which
    // streams (libdeflate cannot).
    //
    // As with igzstream, a truncated or corrupt file reads as one that ends early; read_failed()
    // tells the two apart once reading has stopped. A BGZF file must end with its end-of-file
    // block, so a cut at a block boundary is found too. The buffer is read by one thread at a time
    // (the reader lock serialises it); its inflating thread is stopped by close().
    class ThreadedGzStreambuf : public std::streambuf {
    public:
        static constexpr size_t kBlockSize = size_t{1} << 20;  // inflated bytes per block
        static constexpr size_t kBlocks = 4;
        static constexpr unsigned kZlibBuffer = 1u << 17;      // compressed bytes per read()

        ThreadedGzStreambuf() = default;
        ThreadedGzStreambuf(ThreadedGzStreambuf const&) = delete;
        ThreadedGzStreambuf& operator=(ThreadedGzStreambuf const&) = delete;
        ~ThreadedGzStreambuf() override { close(); }

        bool open(char const* path) {
            if (is_open()) return false;
            m_bgzf = bgzf::StartsAsBgzf(path);
            if (m_bgzf) {
                m_fd = ::open(path, O_RDONLY | O_CLOEXEC);
                if (m_fd < 0) return false;
                m_block_in.resize(bgzf::kMaxBlock);
                m_offset = 0;
                m_last_empty = false;
            } else {
                m_file = gzopen(path, "rb");
                if (!m_file) return false;
                gzbuffer(m_file, kZlibBuffer);
            }
            m_blocks.assign(kBlocks, std::vector<char>(kBlockSize));
            m_free.clear();
            m_ready.clear();
            for (size_t i = 0; i < kBlocks; i++) m_free.push_back(i);
            m_current = kNone;
            m_stop = false;
            m_done = false;
            m_error.clear();
            setg(nullptr, nullptr, nullptr);
            m_thread = std::thread(&ThreadedGzStreambuf::Inflate, this);
            return true;
        }

        bool is_open() const { return m_file != nullptr || m_fd >= 0; }

        void close() {
            if (!is_open()) return;
            {
                std::lock_guard<std::mutex> lock(m_mutex);
                m_stop = true;
            }
            m_cv.notify_all();
            if (m_thread.joinable()) m_thread.join();
            if (m_file) gzclose(m_file);
            if (m_fd >= 0) ::close(m_fd);
            m_file = nullptr;
            m_fd = -1;
            setg(nullptr, nullptr, nullptr);
            m_blocks.clear();
        }

        // zlib reports a truncated or corrupt file as the end of the file; this tells them apart.
        bool read_failed() const {
            std::lock_guard<std::mutex> lock(m_mutex);
            return !m_error.empty();
        }

        std::string read_error_message() const {
            std::lock_guard<std::mutex> lock(m_mutex);
            return m_error;
        }

        // Appends up to `lines` whole lines to out, each with its '\n' (a last line without one gets
        // it, as getline reads such a line too), and returns how many. The lines are found with
        // memchr in the inflated block and copied at once, so the reader lock is held for little
        // more than a copy (BufferedFastxReader::LoadBatch).
        size_t TakeLines(size_t lines, std::string& out) {
            size_t taken = 0;
            bool partial = false;  // out ends inside a line
            while (taken < lines) {
                if (gptr() == egptr() && traits_type::eq_int_type(underflow(), traits_type::eof())) break;
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
            }
            if (partial && taken < lines) {
                out.push_back('\n');
                taken++;
            }
            return taken;
        }

    protected:
        int_type underflow() override {
            if (gptr() < egptr()) return traits_type::to_int_type(*gptr());
            if (!is_open()) return traits_type::eof();
            std::unique_lock<std::mutex> lock(m_mutex);
            if (m_current != kNone) {
                m_free.push_back(m_current);
                m_current = kNone;
                m_cv.notify_all();
            }
            m_cv.wait(lock, [this] { return !m_ready.empty() || m_done; });
            if (m_ready.empty()) {
                setg(nullptr, nullptr, nullptr);
                return traits_type::eof();
            }
            auto const [index, size] = m_ready.front();
            m_ready.pop_front();
            m_current = index;
            char* data = m_blocks[index].data();
            setg(data, data, data + size);
            return traits_type::to_int_type(*gptr());
        }

    private:
        static constexpr size_t kNone = SIZE_MAX;

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
                if (m_bgzf) {
                    end = FillBgzf(data, size, error);
                } else {
                    while (size < kBlockSize) {
                        int const n = gzread(m_file, data + size, static_cast<unsigned>(kBlockSize - size));
                        if (n > 0) {
                            size += static_cast<size_t>(n);
                            continue;
                        }
                        // A gzip file that ends early reads as a normal end of file, but leaves Z_BUF_ERROR.
                        int errnum = Z_OK;
                        char const* message = gzerror(m_file, &errnum);
                        if (n < 0 || errnum != Z_OK) error = (message && *message) ? message : "read error";
                        end = true;
                        break;
                    }
                }
                {
                    std::lock_guard<std::mutex> lock(m_mutex);
                    if (size > 0) m_ready.emplace_back(index, size);
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
                    error = got < bgzf::kHeaderBytes ? "the file ends inside a BGZF block header (truncated file?)"
                                                     : "no BGZF block at byte " + std::to_string(m_offset) +
                                                       " (gzip members of other kinds after BGZF ones?)";
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
                size += n;
                m_offset += block_size;
            }
            return false;
        }

        gzFile m_file = nullptr;
        bool m_bgzf = false;                             // read with FillBgzf, from m_fd
        int m_fd = -1;
        std::vector<unsigned char> m_block_in;           // one compressed BGZF block
        uint64_t m_offset = 0;                           // of the next BGZF block in the file
        bool m_last_empty = false;                       // the last block read was empty (the EOF marker)
        std::thread m_thread;
        std::vector<std::vector<char>> m_blocks;
        std::deque<size_t> m_free;                       // blocks the inflating thread may fill
        std::deque<std::pair<size_t, size_t>> m_ready;   // (block, bytes) in file order
        size_t m_current = kNone;                        // the block the get area points into
        bool m_stop = false;                             // close() asks the inflating thread to end
        bool m_done = false;                             // no more blocks will come
        std::string m_error;
        mutable std::mutex m_mutex;
        std::condition_variable m_cv;
    };

    // The istream for ThreadedGzStreambuf, a drop-in for igzstream when reading.
    class ThreadedGzIstream : public std::istream {
    public:
        explicit ThreadedGzIstream(char const* path) : std::istream(nullptr) {
            std::istream::rdbuf(&m_buf);  // clears the badbit of the null buffer
            if (!m_buf.open(path)) setstate(std::ios::badbit);
        }

        ThreadedGzStreambuf* rdbuf() { return &m_buf; }

        void close() { m_buf.close(); }

    private:
        ThreadedGzStreambuf m_buf;
    };
}
