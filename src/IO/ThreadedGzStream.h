#pragma once

#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <istream>
#include <mutex>
#include <streambuf>
#include <string>
#include <thread>
#include <utility>
#include <vector>
#include <zlib.h>

namespace protal {

    // Reads a gzip-compressed (or plain) file as igzstream does, but inflates it in a thread of its
    // own, up to kBlocks blocks ahead of the reader. protal's alignment threads take reads in
    // batches under one lock (SeqReaderPE, SeqReaderSE); with igzstream the lock holder also
    // inflated the input, which capped a run at one core's inflate speed (~200k read pairs/s)
    // whatever -t was. Here the lock holder only copies bytes that are already inflated, and the
    // two files of a pair inflate in parallel.
    //
    // As with igzstream, a truncated or corrupt file reads as one that ends early; read_failed()
    // tells the two apart once reading has stopped. The buffer is read by one thread at a time
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
            if (m_file) return false;
            m_file = gzopen(path, "rb");
            if (!m_file) return false;
            gzbuffer(m_file, kZlibBuffer);
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

        bool is_open() const { return m_file != nullptr; }

        void close() {
            if (!m_file) return;
            {
                std::lock_guard<std::mutex> lock(m_mutex);
                m_stop = true;
            }
            m_cv.notify_all();
            if (m_thread.joinable()) m_thread.join();
            gzclose(m_file);
            m_file = nullptr;
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

    protected:
        int_type underflow() override {
            if (gptr() < egptr()) return traits_type::to_int_type(*gptr());
            if (!m_file) return traits_type::eof();
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

        gzFile m_file = nullptr;
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
