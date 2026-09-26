// Zstd.h - zstd-compressed database files (index.prx.zst, reference.fna.zst).
//
// Readers open database files through InputFile, which takes a raw or a zstd-compressed file
// (detected by its magic number) and exposes a std::istream either way. Decompression reads the
// file in large chunks on a background thread, so reading from (slow) network storage overlaps
// with decompressing. Writers use OStream, which compresses with several threads and records a
// content checksum that is verified when the file is read back.
#pragma once

#include <zstd.h>

#include <fcntl.h>
#include <unistd.h>

#include <algorithm>
#include <cerrno>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <thread>
#include <vector>

namespace protal::zstd {
    inline const std::string kExtension = ".zst";

    struct Params {
        int level = 19;       // 1..ZSTD_maxCLevel(); decompression speed barely depends on it
        int window_log = 27;  // long-distance matching window, log2 bytes; 0 turns it off
        int threads = 1;      // compression workers
    };

    // True if the file starts with the zstd frame magic number (bytes 28 b5 2f fd).
    inline bool IsCompressed(std::string const& path) {
        std::ifstream is(path, std::ios::binary);
        unsigned char magic[4] = {0, 0, 0, 0};
        is.read(reinterpret_cast<char*>(magic), 4);
        return is.gcount() == 4 && magic[0] == 0x28 && magic[1] == 0xb5 && magic[2] == 0x2f && magic[3] == 0xfd;
    }

    // The file to read for a database file: the raw file if it exists, else its .zst sibling if
    // that exists, else the raw name (which then fails to open with the usual message).
    inline std::string Resolve(std::string const& raw_path) {
        if (std::filesystem::exists(raw_path)) return raw_path;
        if (std::filesystem::exists(raw_path + kExtension)) return raw_path + kExtension;
        return raw_path;
    }

    // Reads a file sequentially in large chunks on a background thread.
    class ReadAhead {
    public:
        struct Chunk {
            std::vector<char> data;
            size_t size = 0;
        };

        ReadAhead(int fd, size_t chunk_size, size_t chunks) : m_fd(fd) {
            m_chunks.resize(chunks);
            for (auto& chunk : m_chunks) {
                chunk.data.resize(chunk_size);
                m_free.push_back(&chunk);
            }
            m_thread = std::thread([this] { Run(); });
        }

        ~ReadAhead() {
            {
                std::lock_guard<std::mutex> lock(m_mutex);
                m_stop = true;
            }
            m_cv.notify_all();
            m_thread.join();
        }

        ReadAhead(ReadAhead const&) = delete;
        ReadAhead& operator=(ReadAhead const&) = delete;

        // The next chunk in file order, or nullptr at the end of the file or after a read error
        // (see Error()). The previously returned chunk is recycled.
        Chunk* Next() {
            std::unique_lock<std::mutex> lock(m_mutex);
            if (m_current) {
                m_free.push_back(m_current);
                m_current = nullptr;
                m_cv.notify_all();
            }
            m_cv.wait(lock, [this] { return !m_filled.empty() || m_done; });
            if (m_filled.empty()) return nullptr;
            m_current = m_filled.front();
            m_filled.pop_front();
            return m_current;
        }

        std::string Error() {
            std::lock_guard<std::mutex> lock(m_mutex);
            return m_error;
        }

    private:
        void Run() {
            while (true) {
                Chunk* chunk = nullptr;
                {
                    std::unique_lock<std::mutex> lock(m_mutex);
                    m_cv.wait(lock, [this] { return m_stop || !m_free.empty(); });
                    if (m_stop) return;
                    chunk = m_free.front();
                    m_free.pop_front();
                }
                size_t got = 0;
                std::string error;
                while (got < chunk->data.size()) {
                    ssize_t const r = ::read(m_fd, chunk->data.data() + got, chunk->data.size() - got);
                    if (r < 0 && errno == EINTR) continue;
                    if (r < 0) {
                        error = std::strerror(errno);
                        break;
                    }
                    if (r == 0) break;
                    got += static_cast<size_t>(r);
                }
                chunk->size = got;
                bool const last = !error.empty() || got < chunk->data.size();
                {
                    std::lock_guard<std::mutex> lock(m_mutex);
                    if (got > 0 && error.empty()) m_filled.push_back(chunk);
                    else m_free.push_back(chunk);
                    if (!error.empty()) m_error = error;
                    if (last) m_done = true;
                }
                m_cv.notify_all();
                if (last) return;
            }
        }

        int m_fd;
        std::vector<Chunk> m_chunks;
        std::deque<Chunk*> m_free;
        std::deque<Chunk*> m_filled;
        Chunk* m_current = nullptr;
        bool m_done = false;
        bool m_stop = false;
        std::string m_error;
        std::mutex m_mutex;
        std::condition_variable m_cv;
        std::thread m_thread;
    };

    // Decompressing input stream buffer over a zstd file (one or more frames). Errors (read
    // errors, corrupt data, a checksum mismatch, a truncated file) are printed once and thrown
    // from the read call; std::istream turns that into badbit.
    class IStreambuf : public std::streambuf {
    public:
        static constexpr size_t kChunkSize = size_t{8} << 20;   // read size, large for network storage
        static constexpr size_t kChunks = 3;
        static constexpr size_t kBufferSize = size_t{1} << 20;  // get area for small reads

        explicit IStreambuf(std::string path) : m_path(std::move(path)), m_buffer(kBufferSize) {
            m_fd = ::open(m_path.c_str(), O_RDONLY | O_CLOEXEC);
            if (m_fd < 0) {
                m_open_error = std::strerror(errno);
                return;
            }
            posix_fadvise(m_fd, 0, 0, POSIX_FADV_SEQUENTIAL);
            m_dctx = ZSTD_createDCtx();
            if (!m_dctx) {
                m_open_error = "cannot allocate a zstd decompression context";
                return;
            }
            // Accept the largest windows (--long up to 2 GB); memory follows the actual frame.
            ZSTD_DCtx_setParameter(m_dctx, ZSTD_d_windowLogMax, ZSTD_dParam_getBounds(ZSTD_d_windowLogMax).upperBound);
            Start();
        }

        ~IStreambuf() override {
            m_reader.reset();
            if (m_dctx) ZSTD_freeDCtx(m_dctx);
            if (m_fd >= 0) ::close(m_fd);
        }

        IStreambuf(IStreambuf const&) = delete;
        IStreambuf& operator=(IStreambuf const&) = delete;

        bool IsOpen() const { return m_open_error.empty(); }
        std::string const& OpenError() const { return m_open_error; }
        std::string const& Error() const { return m_error; }

        // Back to the start of the file.
        bool Rewind() {
            if (!IsOpen()) return false;
            m_reader.reset();
            if (::lseek(m_fd, 0, SEEK_SET) != 0) return false;
            ZSTD_DCtx_reset(m_dctx, ZSTD_reset_session_only);
            m_error.clear();
            Start();
            return true;
        }

    protected:
        int_type underflow() override {
            if (gptr() < egptr()) return traits_type::to_int_type(*gptr());
            size_t const n = Decompress(m_buffer.data(), m_buffer.size());
            setg(m_buffer.data(), m_buffer.data(), m_buffer.data() + n);
            return n ? traits_type::to_int_type(m_buffer[0]) : traits_type::eof();
        }

        // Large reads (the index arrays) are decompressed straight into the caller's memory.
        std::streamsize xsgetn(char* s, std::streamsize n) override {
            std::streamsize done = 0;
            while (done < n) {
                std::streamsize const available = egptr() - gptr();
                if (available > 0) {
                    std::streamsize const k = std::min(available, n - done);
                    std::memcpy(s + done, gptr(), static_cast<size_t>(k));
                    gbump(static_cast<int>(k));
                    done += k;
                } else if (n - done >= static_cast<std::streamsize>(m_buffer.size())) {
                    size_t const got = Decompress(s + done, static_cast<size_t>(n - done));
                    if (got == 0) break;
                    done += static_cast<std::streamsize>(got);
                } else if (traits_type::eq_int_type(underflow(), traits_type::eof())) {
                    break;
                }
            }
            return done;
        }

        // Only reports the position (tellg); the stream cannot seek.
        pos_type seekoff(off_type off, std::ios_base::seekdir dir, std::ios_base::openmode) override {
            if (off != 0 || dir != std::ios_base::cur) return pos_type(off_type(-1));
            return pos_type(static_cast<off_type>(m_delivered - static_cast<uint64_t>(egptr() - gptr())));
        }

        pos_type seekpos(pos_type pos, std::ios_base::openmode which) override {
            pos_type const here = seekoff(0, std::ios_base::cur, which);
            return pos == here ? here : pos_type(off_type(-1));
        }

    private:
        void Start() {
            m_reader = std::make_unique<ReadAhead>(m_fd, kChunkSize, kChunks);
            m_in = {nullptr, 0, 0};
            m_input_done = false;
            m_seen_input = false;
            m_frame_bytes_left = 1;
            m_delivered = 0;
            setg(m_buffer.data(), m_buffer.data(), m_buffer.data());
        }

        [[noreturn]] void Fail(std::string const& what) {
            if (m_error.empty()) {
                m_error = what;
                std::cerr << "Error reading " << m_path << ": " << what << std::endl;
            }
            throw std::runtime_error(m_path + ": " + what);
        }

        // Decompresses up to capacity bytes into dst; 0 only at the clean end of the data.
        size_t Decompress(char* dst, size_t capacity) {
            if (!IsOpen()) return 0;
            if (!m_error.empty()) Fail(m_error);
            ZSTD_outBuffer out{dst, capacity, 0};
            while (out.pos < out.size) {
                if (m_in.pos == m_in.size && !m_input_done) {
                    ReadAhead::Chunk* chunk = m_reader->Next();
                    if (chunk) {
                        m_in = {chunk->data.data(), chunk->size, 0};
                        m_seen_input = true;
                    } else {
                        m_input_done = true;
                        m_in = {nullptr, 0, 0};
                        std::string const error = m_reader->Error();
                        if (!error.empty()) Fail("read error: " + error);
                        if (!m_seen_input) Fail("the file is empty");
                    }
                }
                size_t const out_before = out.pos, in_before = m_in.pos;
                size_t const ret = ZSTD_decompressStream(m_dctx, &out, &m_in);
                if (ZSTD_isError(ret)) Fail(std::string(ZSTD_getErrorName(ret)) + " (corrupt or not a zstd file?)");
                bool const progress = out.pos != out_before || m_in.pos != in_before;
                // Without progress, ret only hints at the header of a possible next frame.
                if (progress) m_frame_bytes_left = ret;
                if (m_input_done && m_in.pos == m_in.size && !progress) {
                    if (m_frame_bytes_left != 0) Fail("the file ends inside a zstd frame (truncated or corrupt file?)");
                    break;
                }
            }
            m_delivered += out.pos;
            return out.pos;
        }

        std::string m_path;
        std::string m_open_error;
        std::string m_error;
        int m_fd = -1;
        ZSTD_DCtx* m_dctx = nullptr;
        std::unique_ptr<ReadAhead> m_reader;
        ZSTD_inBuffer m_in{nullptr, 0, 0};
        bool m_input_done = false;
        bool m_seen_input = false;
        size_t m_frame_bytes_left = 1;  // hint of the last decompression call that made progress; 0: frame complete
        uint64_t m_delivered = 0;       // decompressed bytes handed out so far
        std::vector<char> m_buffer;
    };

    // Compressing output stream buffer: one zstd frame with a content checksum, written with
    // Params::threads workers. Close() finishes the frame; it returns false on any error.
    class OStreambuf : public std::streambuf {
    public:
        static constexpr size_t kInBufferSize = size_t{1} << 20;
        static constexpr size_t kOutBufferSize = size_t{4} << 20;

        OStreambuf(std::string path, Params const& params, std::optional<uint64_t> content_size)
                : m_path(std::move(path)), m_in_buffer(kInBufferSize), m_out_buffer(kOutBufferSize) {
            setp(m_in_buffer.data(), m_in_buffer.data() + m_in_buffer.size());
            m_file = std::fopen(m_path.c_str(), "wb");
            if (!m_file) {
                SetError(std::string("cannot open for writing: ") + std::strerror(errno));
                return;
            }
            m_cctx = ZSTD_createCCtx();
            if (!m_cctx) {
                SetError("cannot allocate a zstd compression context");
                return;
            }
            Set(ZSTD_c_compressionLevel, params.level, "compression level");
            Set(ZSTD_c_checksumFlag, 1, "checksum flag");
            if (params.window_log > 0) {
                Set(ZSTD_c_enableLongDistanceMatching, 1, "long-distance matching");
                Set(ZSTD_c_windowLog, params.window_log, "window log");
            }
            if (params.threads > 1 && ZSTD_isError(ZSTD_CCtx_setParameter(m_cctx, ZSTD_c_nbWorkers, params.threads))) {
                std::cerr << "Warning: this libzstd has no multithreading support, compressing " << m_path
                          << " with one thread" << std::endl;
            }
            if (content_size && ZSTD_isError(ZSTD_CCtx_setPledgedSrcSize(m_cctx, *content_size))) {
                SetError("cannot set the content size");
            }
        }

        ~OStreambuf() override {
            Close();
            if (m_cctx) ZSTD_freeCCtx(m_cctx);
        }

        OStreambuf(OStreambuf const&) = delete;
        OStreambuf& operator=(OStreambuf const&) = delete;

        bool Close() {
            if (m_closed) return m_error.empty();
            m_closed = true;
            if (m_error.empty() && FlushPutArea()) Feed(nullptr, 0, ZSTD_e_end);
            if (m_file) {
                if (std::fclose(m_file) != 0 && m_error.empty()) SetError(std::string("closing failed: ") + std::strerror(errno));
                m_file = nullptr;
            }
            return m_error.empty();
        }

        std::string const& Error() const { return m_error; }
        uint64_t BytesIn() const { return m_bytes_in + static_cast<uint64_t>(pptr() - pbase()); }
        uint64_t BytesOut() const { return m_bytes_out; }

    protected:
        int_type overflow(int_type ch) override {
            if (!FlushPutArea()) return traits_type::eof();
            if (!traits_type::eq_int_type(ch, traits_type::eof())) {
                *pptr() = traits_type::to_char_type(ch);
                pbump(1);
            }
            return traits_type::not_eof(ch);
        }

        std::streamsize xsputn(char const* s, std::streamsize n) override {
            if (n <= epptr() - pptr()) {
                std::memcpy(pptr(), s, static_cast<size_t>(n));
                pbump(static_cast<int>(n));
                return n;
            }
            if (!FlushPutArea()) return 0;
            if (n < static_cast<std::streamsize>(m_in_buffer.size())) {
                std::memcpy(pptr(), s, static_cast<size_t>(n));
                pbump(static_cast<int>(n));
                return n;
            }
            // Large writes (the index arrays) are compressed straight from the caller's memory.
            m_bytes_in += static_cast<uint64_t>(n);
            return Feed(s, static_cast<size_t>(n), ZSTD_e_continue) ? n : 0;
        }

        // Flushing the stream does not end a zstd block: that would only cost compression.
        int sync() override { return m_error.empty() ? 0 : -1; }

    private:
        void Set(ZSTD_cParameter parameter, int value, char const* what) {
            if (m_error.empty() && ZSTD_isError(ZSTD_CCtx_setParameter(m_cctx, parameter, value))) {
                SetError(std::string("invalid zstd ") + what + " " + std::to_string(value));
            }
        }

        void SetError(std::string const& what) {
            if (m_error.empty()) {
                m_error = what;
                std::cerr << "Error writing " << m_path << ": " << what << std::endl;
            }
        }

        bool FlushPutArea() {
            size_t const n = static_cast<size_t>(pptr() - pbase());
            setp(m_in_buffer.data(), m_in_buffer.data() + m_in_buffer.size());
            if (n == 0) return m_error.empty();
            m_bytes_in += n;
            return Feed(m_in_buffer.data(), n, ZSTD_e_continue);
        }

        bool Feed(char const* data, size_t size, ZSTD_EndDirective mode) {
            if (!m_error.empty()) return false;
            ZSTD_inBuffer in{data, size, 0};
            while (true) {
                ZSTD_outBuffer out{m_out_buffer.data(), m_out_buffer.size(), 0};
                size_t const remaining = ZSTD_compressStream2(m_cctx, &out, &in, mode);
                if (ZSTD_isError(remaining)) {
                    SetError(ZSTD_getErrorName(remaining));
                    return false;
                }
                if (out.pos > 0) {
                    if (std::fwrite(m_out_buffer.data(), 1, out.pos, m_file) != out.pos) {
                        SetError(std::string("write failed: ") + std::strerror(errno));
                        return false;
                    }
                    m_bytes_out += out.pos;
                }
                bool const finished = mode == ZSTD_e_end ? remaining == 0 : in.pos == in.size;
                if (finished) return true;
            }
        }

        std::string m_path;
        std::string m_error;
        std::FILE* m_file = nullptr;
        ZSTD_CCtx* m_cctx = nullptr;
        std::vector<char> m_in_buffer;
        std::vector<char> m_out_buffer;
        uint64_t m_bytes_in = 0;
        uint64_t m_bytes_out = 0;
        bool m_closed = false;
    };

    // std::ostream writing one zstd frame to a file. Close() before checking for success.
    class OStream : public std::ostream {
    public:
        OStream(std::string const& path, Params const& params, std::optional<uint64_t> content_size = std::nullopt)
                : std::ostream(nullptr), m_buf(path, params, content_size) {
            rdbuf(&m_buf);
            if (!m_buf.Error().empty()) setstate(std::ios_base::badbit);
        }

        bool Close() {
            flush();
            bool const ok = m_buf.Close() && !fail();
            if (!ok) setstate(std::ios_base::badbit);
            return ok;
        }

        OStreambuf const& Buffer() const { return m_buf; }

    private:
        OStreambuf m_buf;
    };

    // A database file opened for sequential reading, raw or zstd-compressed.
    class InputFile {
    public:
        static constexpr size_t kRawBufferSize = size_t{4} << 20;

        explicit InputFile(std::string path) : m_path(std::move(path)) {
            m_compressed = IsCompressed(m_path);
            if (m_compressed) {
                m_zbuf = std::make_unique<IStreambuf>(m_path);
                m_zin = std::make_unique<std::istream>(m_zbuf.get());
                if (!m_zbuf->IsOpen()) m_zin->setstate(std::ios_base::badbit);
                m_stream = m_zin.get();
            } else {
                m_raw_buffer.resize(kRawBufferSize);
                m_raw.rdbuf()->pubsetbuf(m_raw_buffer.data(), static_cast<std::streamsize>(m_raw_buffer.size()));
                m_raw.open(m_path, std::ios::binary);
                m_stream = &m_raw;
            }
        }

        InputFile(InputFile const&) = delete;
        InputFile& operator=(InputFile const&) = delete;

        bool IsOpen() const { return m_compressed ? m_zbuf->IsOpen() : m_raw.is_open(); }
        bool Compressed() const { return m_compressed; }
        std::string const& Path() const { return m_path; }
        std::istream& Stream() { return *m_stream; }

        // Back to the start (the same std::istream object stays valid).
        bool Rewind() {
            m_stream->clear();
            if (m_compressed) return m_zbuf->Rewind();
            m_raw.seekg(0, std::ios::beg);
            return static_cast<bool>(m_raw);
        }

    private:
        std::string m_path;
        bool m_compressed = false;
        std::vector<char> m_raw_buffer;
        std::ifstream m_raw;
        std::unique_ptr<IStreambuf> m_zbuf;
        std::unique_ptr<std::istream> m_zin;
        std::istream* m_stream = nullptr;
    };

    // Compresses src into dst (via dst.partial, renamed at the end). With verify, the new file is
    // decompressed and compared with src before the rename. Returns false with a message in error.
    inline bool CompressFile(std::string const& src, std::string const& dst, Params const& params, bool verify,
                             std::string& error) {
        std::error_code ec;
        uint64_t const size = std::filesystem::file_size(src, ec);
        if (ec) {
            error = "cannot read " + src + ": " + ec.message();
            return false;
        }
        std::string const partial = dst + ".partial";
        auto fail = [&](std::string const& what) {
            error = what;
            std::filesystem::remove(partial, ec);
            return false;
        };
        std::vector<char> buffer(size_t{8} << 20);
        {
            std::ifstream in(src, std::ios::binary);
            OStream out(partial, params, size);
            while (in && out) {
                in.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
                out.write(buffer.data(), in.gcount());
            }
            if (!out.Close()) return fail("writing " + partial + " failed: " + out.Buffer().Error());
            if (in.bad() || !in.eof()) return fail("reading " + src + " failed");
        }
        if (verify) {
            std::ifstream original(src, std::ios::binary);
            InputFile copy(partial);
            std::vector<char> decompressed(buffer.size());
            uint64_t compared = 0;
            while (true) {
                original.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
                copy.Stream().read(decompressed.data(), static_cast<std::streamsize>(decompressed.size()));
                if (copy.Stream().bad()) return fail("verifying " + partial + " failed: it cannot be decompressed");
                std::streamsize const a = original.gcount(), b = copy.Stream().gcount();
                if (a != b || std::memcmp(buffer.data(), decompressed.data(), static_cast<size_t>(a)) != 0) {
                    return fail("verifying " + partial + " failed: content differs after byte " + std::to_string(compared));
                }
                compared += static_cast<uint64_t>(a);
                if (a == 0) break;
            }
            if (compared != size) return fail("verifying " + partial + " failed: " + std::to_string(compared) + " of " +
                                              std::to_string(size) + " bytes compared");
        }
        std::filesystem::rename(partial, dst, ec);
        if (ec) return fail("cannot rename " + partial + " to " + dst + ": " + ec.message());
        return true;
    }
}
