// SamChunks.h - a SAM stream cut into chunks of whole reads, so that one sample can be profiled on
// several threads (Profiler::ProfileSam): a thread of its own reads (and decompresses) the stream
// ahead, and each chunk is parsed on its own. A chunk ends before a record line whose QNAME differs
// from that of the record line before it, so a read's records (both mates, all its candidates, a long
// read's genes) stay in one chunk; header and blank lines stay in the chunk they are in.
#pragma once

#include <algorithm>
#include <chrono>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <deque>
#include <exception>
#include <istream>
#include <mutex>
#include <optional>
#include <streambuf>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

namespace protal::sam_chunks {

    // The QNAME of the SAM line [begin, end) (without its newline), or nullopt for a header or blank line.
    inline std::optional<std::string_view> RecordName(char const* begin, char const* end) {
        if (end > begin && end[-1] == '\r') end--;
        if (begin == end || *begin == '@') return std::nullopt;
        auto const* tab = static_cast<char const*>(std::memchr(begin, '\t', static_cast<size_t>(end - begin)));
        return std::string_view(begin, static_cast<size_t>((tab ? tab : end) - begin));
    }

    // Where the first chunk of text[0, size) ends: at the start of the last complete record line whose QNAME
    // differs from that of the record line before it. 0 if there is none: no complete line, or one read only.
    inline size_t ChunkEnd(char const* text, size_t size) {
        size_t line_end = size;  // one past the newline that ends the current line
        while (line_end > 0 && text[line_end - 1] != '\n') line_end--;
        std::optional<std::string_view> next;  // the QNAME of the record line after the current one
        size_t next_start = 0;
        while (line_end > 0) {
            size_t line_start = line_end - 1;
            while (line_start > 0 && text[line_start - 1] != '\n') line_start--;
            if (auto const name = RecordName(text + line_start, text + line_end - 1)) {
                if (next && *name != *next) return next_start;
                next = name;
                next_start = line_start;
            }
            line_end = line_start;
        }
        return 0;
    }

    struct Chunk {
        std::string text;
        size_t index = 0;       // its place in the stream, from 0
        size_t first_line = 0;  // lines before it
        bool last = false;      // the stream's last chunk (it may end without a newline)
    };

    // Reads a stream in a thread of its own into chunks of about `bytes` bytes, at most `ahead` of them
    // waiting to be taken; a read longer than that goes into one larger chunk. The stream is the reader's
    // until the last chunk is taken or the reader is stopped.
    class ChunkReader {
    public:
        ChunkReader(std::istream& is, size_t bytes, size_t ahead) :
                m_is(is), m_bytes(std::max<size_t>(bytes, 1)), m_ahead(std::max<size_t>(ahead, 1)) {
            m_thread = std::thread(&ChunkReader::Run, this);
        }

        ~ChunkReader() { Halt(); }

        ChunkReader(ChunkReader const&) = delete;
        ChunkReader& operator=(ChunkReader const&) = delete;

        // The next chunk in stream order, or false after the last one (or after Stop).
        bool Next(Chunk& chunk) {
            std::unique_lock<std::mutex> lock(m_mutex);
            m_cv.wait(lock, [this] { return !m_ready.empty() || m_done; });
            if (m_ready.empty()) return false;
            chunk = std::move(m_ready.front());
            m_ready.pop_front();
            m_cv.notify_all();
            return true;
        }

        // Ends reading (also before the end of the stream) and waits for the thread; throws what reading
        // threw (a failed allocation, say).
        void Stop() {
            Halt();
            if (m_failure) std::rethrow_exception(std::exchange(m_failure, nullptr));
        }

        // Where the reading thread spent its time, and the text it read; complete once the last chunk is taken (or after Stop).
        struct Times {
            double read = 0;  // in the stream's read: decompressing, or waiting for the decompressing threads
            double cut = 0;   // cutting the text into chunks of whole reads and counting their lines
            double wait = 0;  // waiting for room: the chunks read ahead not yet taken
            uint64_t bytes = 0;
            size_t chunks = 0;
        };
        Times GetTimes() const {
            std::lock_guard<std::mutex> lock(m_mutex);
            return m_times;
        }

        // Whether the stream went bad (a read error); known once the last chunk is taken.
        bool Bad() const {
            std::lock_guard<std::mutex> lock(m_mutex);
            return m_bad;
        }

    private:
        void Halt() {
            {
                std::lock_guard<std::mutex> lock(m_mutex);
                m_stop = true;
            }
            m_cv.notify_all();
            if (m_thread.joinable()) m_thread.join();
        }

        void Run() {
            try {
                Read();
            } catch (...) {
                std::lock_guard<std::mutex> lock(m_mutex);
                m_failure = std::current_exception();
            }
            {
                std::lock_guard<std::mutex> lock(m_mutex);
                m_done = true;
            }
            m_cv.notify_all();
        }

        using Clock = std::chrono::steady_clock;
        static double Since(Clock::time_point& start) {
            auto const now = Clock::now();
            double const seconds = std::chrono::duration<double>(now - start).count();
            start = now;
            return seconds;
        }

        void Read() {
            std::string buffer;
            size_t index = 0, lines = 0;
            Times times;
            while (true) {
                auto clock = Clock::now();
                size_t const had = buffer.size();
                buffer.resize(had + m_bytes);
                m_is.read(buffer.data() + had, static_cast<std::streamsize>(m_bytes));
                buffer.resize(had + static_cast<size_t>(m_is.gcount()));
                times.bytes += static_cast<uint64_t>(m_is.gcount());
                times.read += Since(clock);
                bool const end = !m_is;
                size_t cut = end ? buffer.size() : ChunkEnd(buffer.data(), buffer.size());
                if (cut == 0 && !end) {  // one read so far: read on
                    times.cut += Since(clock);
                    continue;
                }
                Chunk chunk;
                chunk.index = index++;
                chunk.first_line = lines;
                chunk.last = end;
                std::string rest(buffer, cut);
                buffer.resize(cut);
                chunk.text = std::move(buffer);
                buffer = std::move(rest);
                lines += static_cast<size_t>(std::count(chunk.text.begin(), chunk.text.end(), '\n'));
                times.cut += Since(clock);
                times.chunks++;
                if (!Put(std::move(chunk), end && m_is.bad(), times, clock) || end) return;
            }
        }

        // False if the reader was stopped. `bad`: the stream went bad (set with the last chunk). The wait for room is
        // added to `times`, which are then the reader's.
        bool Put(Chunk&& chunk, bool bad, Times& times, Clock::time_point& clock) {
            std::unique_lock<std::mutex> lock(m_mutex);
            m_cv.wait(lock, [this] { return m_ready.size() < m_ahead || m_stop; });
            times.wait += Since(clock);
            m_times = times;
            if (m_stop) return false;
            m_bad = m_bad || bad;
            m_ready.push_back(std::move(chunk));
            m_cv.notify_all();
            return true;
        }

        std::istream& m_is;
        size_t const m_bytes;
        size_t const m_ahead;
        std::deque<Chunk> m_ready;
        bool m_done = false;
        bool m_stop = false;
        bool m_bad = false;
        Times m_times;
        std::exception_ptr m_failure;
        mutable std::mutex m_mutex;
        std::condition_variable m_cv;
        std::thread m_thread;
    };

    // Runs work(i) for every i in [0, n) on up to `threads` threads (the calling one among them), each
    // taking the next i that is left. The first exception a call throws is thrown once all have ended.
    template<typename Work>
    void ParallelFor(size_t n, size_t threads, Work&& work) {
        size_t const workers = std::max<size_t>(1, std::min(threads, n));
        std::mutex mutex;
        size_t next = 0;
        std::exception_ptr failure;
        auto run = [&]() {
            while (true) {
                size_t i;
                {
                    std::lock_guard<std::mutex> lock(mutex);
                    if (next >= n || failure) return;
                    i = next++;
                }
                try {
                    work(i);
                } catch (...) {
                    std::lock_guard<std::mutex> lock(mutex);
                    if (!failure) failure = std::current_exception();
                }
            }
        };
        std::vector<std::thread> pool;
        pool.reserve(workers - 1);
        for (size_t t = 1; t < workers; t++) pool.emplace_back(run);
        run();
        for (auto& thread : pool) thread.join();
        if (failure) std::rethrow_exception(failure);
    }
}
