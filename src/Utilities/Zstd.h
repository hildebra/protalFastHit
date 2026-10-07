// Zstd.h - zstd-compressed database files (index.prx.zst, reference.fna.zst).
//
// Database files are written in zstd's seekable format (CompressFrames): independent frames of
// Params::frame_size bytes, each with a content checksum, followed by a seek table. ParallelRead
// decompresses such a file with several threads, frame by frame, straight into the caller's
// memory; it reads raw files in chunks the same way, and any other zstd file sequentially.
// InputFile gives a sequential std::istream over a raw or zstd file (reading ahead on a
// background thread); plain zstd tools read the seekable files too. ParallelReadFrames and
// FramesInput read some of a file's frames only: a member of the single-file database
// (Database.h).
#pragma once

#include <zstd.h>

#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
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
        int window_log = 27;  // long-distance matching window, log2 bytes (capped at the frame size); 0: off
        int threads = 1;      // compression workers
        uint64_t frame_size = uint64_t{64} << 20;  // bytes per independent frame (seekable); 0: one frame
    };

    // True if the file starts with a zstd frame (magic bytes 28 b5 2f fd) or a skippable frame
    // (50..5f 2a 4d 18, e.g. written by pzstd).
    inline bool IsCompressed(std::string const& path) {
        std::ifstream is(path, std::ios::binary);
        unsigned char magic[4] = {0, 0, 0, 0};
        is.read(reinterpret_cast<char*>(magic), 4);
        if (is.gcount() != 4) return false;
        bool const frame = magic[0] == 0x28 && magic[1] == 0xb5 && magic[2] == 0x2f && magic[3] == 0xfd;
        bool const skippable = (magic[0] & 0xf0) == 0x50 && magic[1] == 0x2a && magic[2] == 0x4d && magic[3] == 0x18;
        return frame || skippable;
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

        // Decompresses up to capacity bytes into dst; 0 only at the clean end of the data. Bytes
        // decoded before an error are returned first; the error is thrown by the next call.
        size_t Decompress(char* dst, size_t capacity) {
            if (!IsOpen()) return 0;
            if (!m_error.empty()) Fail(m_error);
            ZSTD_outBuffer out{dst, capacity, 0};
            std::string error;
            while (out.pos < out.size && error.empty()) {
                if (m_in.pos == m_in.size && !m_input_done) {
                    ReadAhead::Chunk* chunk = m_reader->Next();
                    if (chunk) {
                        m_in = {chunk->data.data(), chunk->size, 0};
                        m_seen_input = true;
                    } else {
                        m_input_done = true;
                        m_in = {nullptr, 0, 0};
                        std::string const read_error = m_reader->Error();
                        if (!read_error.empty()) error = "read error: " + read_error;
                        else if (!m_seen_input) error = "the file is empty";
                        if (!error.empty()) break;
                    }
                }
                size_t const out_before = out.pos, in_before = m_in.pos;
                size_t const ret = ZSTD_decompressStream(m_dctx, &out, &m_in);
                if (ZSTD_isError(ret)) {
                    error = std::string(ZSTD_getErrorName(ret)) + " (corrupt or not a zstd file?)";
                    break;
                }
                bool const progress = out.pos != out_before || m_in.pos != in_before;
                // Without progress, ret only hints at the header of a possible next frame.
                if (progress) m_frame_bytes_left = ret;
                if (m_input_done && m_in.pos == m_in.size && !progress) {
                    if (m_frame_bytes_left != 0) error = "the file ends inside a zstd frame (truncated or corrupt file?)";
                    break;
                }
            }
            if (!error.empty()) {
                if (out.pos == 0) Fail(error);
                if (m_error.empty()) {
                    m_error = error;
                    std::cerr << "Error reading " << m_path << ": " << error << std::endl;
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

    // A database file's content, read sequentially (InputFile, FramesInput).
    class Input {
    public:
        virtual ~Input() = default;
        virtual bool IsOpen() const = 0;
        virtual std::istream& Stream() = 0;
    };

    // A database file opened for sequential reading, raw or zstd-compressed.
    class InputFile : public Input {
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

        bool IsOpen() const override { return m_compressed ? m_zbuf->IsOpen() : m_raw.is_open(); }
        bool Compressed() const { return m_compressed; }
        std::string const& Path() const { return m_path; }
        std::istream& Stream() override { return *m_stream; }

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

    // ---- Seekable format ------------------------------------------------------------------------
    // zstd's contrib/seekable_format: independent frames, then a skippable frame (magic 0x184D2A5E)
    // holding one entry per frame (compressed size, decompressed size; 4 bytes each, optionally a
    // 4-byte checksum) and a footer (frame count, descriptor byte, magic 0x8F92EAB1).

    inline constexpr uint32_t kSeekTableMagic = 0x8F92EAB1;
    inline constexpr uint32_t kSeekTableFrameMagic = 0x184D2A5E;

    struct SeekTable {
        struct Frame {
            uint64_t compressed_offset, compressed_size, decompressed_offset, decompressed_size;
        };
        std::vector<Frame> frames;

        uint64_t DecompressedSize() const {
            return frames.empty() ? 0 : frames.back().decompressed_offset + frames.back().decompressed_size;
        }
    };

    inline uint32_t ReadLE32(unsigned char const* p) {
        return uint32_t(p[0]) | uint32_t(p[1]) << 8 | uint32_t(p[2]) << 16 | uint32_t(p[3]) << 24;
    }

    inline void PutLE32(std::string& s, uint32_t v) {
        for (int i = 0; i < 4; i++) s.push_back(static_cast<char>((v >> (8 * i)) & 0xff));
    }

    // The seek table if the file ends with one. A seek table that does not fit the file sets error.
    inline std::optional<SeekTable> ReadSeekTable(std::string const& path, std::string& error) {
        std::error_code ec;
        uint64_t const size = std::filesystem::file_size(path, ec);
        if (ec || size < 17) return std::nullopt;
        std::ifstream is(path, std::ios::binary);
        unsigned char footer[9];
        is.seekg(static_cast<std::streamoff>(size - 9));
        is.read(reinterpret_cast<char*>(footer), 9);
        if (!is || ReadLE32(footer + 5) != kSeekTableMagic) return std::nullopt;
        uint32_t const n = ReadLE32(footer);
        unsigned char const descriptor = footer[4];
        if (descriptor & 0x7c) {
            error = "invalid seek table (reserved bits set)";
            return std::nullopt;
        }
        uint64_t const entry = (descriptor & 0x80) ? 12 : 8;
        uint64_t const table = uint64_t(n) * entry + 9;
        if (table + 8 > size) {
            error = "invalid seek table (larger than the file)";
            return std::nullopt;
        }
        uint64_t const table_start = size - table - 8;
        std::vector<unsigned char> bytes(table + 8);
        is.seekg(static_cast<std::streamoff>(table_start));
        is.read(reinterpret_cast<char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
        if (!is || ReadLE32(bytes.data()) != kSeekTableFrameMagic || ReadLE32(bytes.data() + 4) != table) {
            error = "invalid seek table (bad frame header)";
            return std::nullopt;
        }
        SeekTable result;
        result.frames.reserve(n);
        uint64_t compressed = 0, decompressed = 0;
        for (uint32_t i = 0; i < n; i++) {
            unsigned char const* e = bytes.data() + 8 + i * entry;
            uint64_t const c = ReadLE32(e), d = ReadLE32(e + 4);
            result.frames.push_back({compressed, c, decompressed, d});
            compressed += c;
            decompressed += d;
        }
        if (compressed != table_start) {
            error = "invalid seek table (its frames take " + std::to_string(compressed) + " bytes, the file has " +
                    std::to_string(table_start) + ")";
            return std::nullopt;
        }
        return result;
    }

    inline bool IsSeekable(std::string const& path) {
        std::string error;
        return ReadSeekTable(path, error).has_value();
    }

    // Receives a file's content from ParallelRead. Calls come from several threads at once, but
    // always for disjoint byte ranges.
    class Sink {
    public:
        virtual ~Sink() = default;
        // Memory for [offset, offset + size) if the whole range goes to one place, else nullptr.
        virtual char* Direct(uint64_t offset, size_t size) { return nullptr; }
        virtual void Copy(uint64_t offset, char const* data, size_t size) = 0;
    };

    inline bool PreadAll(int fd, char* dst, size_t size, uint64_t offset) {
        while (size > 0) {
            ssize_t const r = ::pread(fd, dst, size, static_cast<off_t>(offset));
            if (r < 0 && errno == EINTR) continue;
            if (r <= 0) return false;
            dst += r;
            size -= static_cast<size_t>(r);
            offset += static_cast<uint64_t>(r);
        }
        return true;
    }

    inline size_t WorkerCount(size_t units, int threads) {
        return std::max<size_t>(1, std::min<size_t>(static_cast<size_t>(std::max(threads, 1)), units));
    }

    // An allocator of anonymous mappings: memory freed goes back to the system at once (munmap), not into malloc's
    // arenas, where glibc keeps freed blocks under its mmap threshold resident (in each thread's arena). For the frame
    // buffers of ForEachFrame's workers, freed when a worker stops before the end of a load.
    template <typename T>
    struct MappedAllocator {
        using value_type = T;
        MappedAllocator() = default;
        template <typename U>
        MappedAllocator(MappedAllocator<U> const&) {}
        T* allocate(size_t n) {
            void* p = ::mmap(nullptr, std::max<size_t>(n, 1) * sizeof(T), PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
            if (p == MAP_FAILED) throw std::bad_alloc();
            return static_cast<T*>(p);
        }
        void deallocate(T* p, size_t n) { ::munmap(p, std::max<size_t>(n, 1) * sizeof(T)); }
        template <typename U>
        bool operator==(MappedAllocator<U> const&) const { return true; }
        template <typename U>
        bool operator!=(MappedAllocator<U> const&) const { return false; }
    };
    using FrameBuffer = std::vector<char, MappedAllocator<char>>;

    // Runs work(unit, worker) for every unit in [0, units) on WorkerCount(units, threads) threads
    // (worker in [0, WorkerCount)); stops at the first failure. work returns an error message,
    // empty on success.
    template<typename Work>
    std::string ParallelFor(size_t units, int threads, Work&& work) {
        std::atomic<size_t> next{0};
        std::atomic<bool> failed{false};
        std::mutex mutex;
        std::string first_error;
        auto run = [&](size_t worker) {
            for (size_t i; !failed && (i = next++) < units;) {
                std::string error;
                try {
                    error = work(i, worker);
                } catch (std::exception const& e) {
                    error = e.what();
                }
                if (!error.empty()) {
                    std::lock_guard<std::mutex> lock(mutex);
                    if (first_error.empty()) first_error = error;
                    failed = true;
                }
            }
        };
        size_t const n = WorkerCount(units, threads);
        std::vector<std::thread> pool;
        for (size_t t = 1; t < n; t++) pool.emplace_back(run, t);
        run(0);
        for (auto& t : pool) t.join();
        return first_error;
    }

    // Reads the content of the frames `table` lists in a file (all of a seekable file, or one member
    // of a single-file database) into sink with up to `threads` threads, each frame decompressed
    // straight into sink.Direct memory where possible. Returns the number of content bytes; on
    // failure sets error.
    inline uint64_t ParallelReadFrames(std::string const& path, SeekTable const& table, int threads, Sink& sink,
                                       std::string& error) {
        int const fd = ::open(path.c_str(), O_RDONLY | O_CLOEXEC);
        if (fd < 0) {
            error = std::string("cannot open the file: ") + std::strerror(errno);
            return 0;
        }
        size_t const units = table.frames.size();
        int const window_max = ZSTD_dParam_getBounds(ZSTD_d_windowLogMax).upperBound;
        struct State {
            ZSTD_DCtx* dctx = nullptr;
            std::vector<char> input, output;
            State() = default;
            State(State const&) = delete;
            ~State() { if (dctx) ZSTD_freeDCtx(dctx); }
        };
        std::vector<State> states(WorkerCount(units, threads));
        std::string const failure = ParallelFor(units, threads, [&](size_t i, size_t worker) -> std::string {
            State* s = &states[worker];
            auto const& frame = table.frames[i];
            std::string const where = "frame " + std::to_string(i + 1) + " of " + std::to_string(units);
            s->input.resize(frame.compressed_size);
            if (!PreadAll(fd, s->input.data(), frame.compressed_size, frame.compressed_offset)) {
                return where + ": read error";
            }
            unsigned long long const content = ZSTD_getFrameContentSize(s->input.data(), s->input.size());
            if (content == ZSTD_CONTENTSIZE_ERROR) return where + ": not a zstd frame (corrupt file?)";
            if (content != ZSTD_CONTENTSIZE_UNKNOWN && content != frame.decompressed_size) {
                return where + ": holds " + std::to_string(content) + " bytes, the seek table says " +
                       std::to_string(frame.decompressed_size);
            }
            size_t const n = static_cast<size_t>(frame.decompressed_size);
            char* dst = sink.Direct(frame.decompressed_offset, n);
            if (!dst) {
                s->output.resize(n);
                dst = s->output.data();
            }
            if (!s->dctx) {
                s->dctx = ZSTD_createDCtx();
                if (!s->dctx) return std::string("cannot allocate a zstd decompression context");
                ZSTD_DCtx_setParameter(s->dctx, ZSTD_d_windowLogMax, window_max);
            }
            size_t const got = ZSTD_decompressDCtx(s->dctx, dst, n, s->input.data(), s->input.size());
            if (ZSTD_isError(got)) return where + ": " + ZSTD_getErrorName(got) + " (corrupt file?)";
            if (got != n) return where + ": " + std::to_string(got) + " bytes, the seek table says " + std::to_string(n);
            if (dst == s->output.data()) sink.Copy(frame.decompressed_offset, dst, n);
            return "";
        });
        ::close(fd);
        if (!failure.empty()) {
            error = failure;
            return 0;
        }
        return table.DecompressedSize();
    }

    // Reads the whole content of a raw or zstd file into sink with up to `threads` threads: a raw
    // file in 64 MB chunks, a seekable zstd file frame by frame (ParallelReadFrames), any other zstd
    // file sequentially. Returns the number of content bytes; on failure sets error.
    inline uint64_t ParallelRead(std::string const& path, int threads, Sink& sink, std::string& error) {
        constexpr uint64_t kRawChunk = uint64_t{64} << 20;
        bool const compressed = IsCompressed(path);
        std::optional<SeekTable> table;
        if (compressed) {
            table = ReadSeekTable(path, error);
            if (!error.empty()) return 0;
            if (table) return ParallelReadFrames(path, *table, threads, sink, error);
        }
        if (compressed) {
            // Not seekable (e.g. from the zstd CLI): one stream.
            InputFile in(path);
            if (!in.IsOpen()) {
                error = "cannot open the file";
                return 0;
            }
            std::vector<char> buffer(size_t{8} << 20);
            uint64_t offset = 0;
            while (true) {
                in.Stream().read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
                size_t const got = static_cast<size_t>(in.Stream().gcount());
                if (got > 0) sink.Copy(offset, buffer.data(), got);
                offset += got;
                if (!in.Stream()) break;
            }
            if (in.Stream().bad()) error = "the file cannot be decompressed (truncated or corrupt file?)";
            return offset;
        }

        int const fd = ::open(path.c_str(), O_RDONLY | O_CLOEXEC);
        if (fd < 0) {
            error = std::string("cannot open the file: ") + std::strerror(errno);
            return 0;
        }
        std::error_code ec;
        uint64_t const total = std::filesystem::file_size(path, ec);
        size_t const units = static_cast<size_t>((total + kRawChunk - 1) / kRawChunk);
        std::vector<std::vector<char>> buffers(WorkerCount(units, threads));
        std::string const failure = ParallelFor(units, threads, [&](size_t i, size_t worker) -> std::string {
            uint64_t const offset = i * kRawChunk;
            size_t const n = static_cast<size_t>(std::min(kRawChunk, total - offset));
            char* dst = sink.Direct(offset, n);
            if (!dst) {
                buffers[worker].resize(n);
                dst = buffers[worker].data();
            }
            if (!PreadAll(fd, dst, n, offset)) return "read error at byte " + std::to_string(offset);
            if (dst == buffers[worker].data()) sink.Copy(offset, dst, n);
            return "";
        });
        ::close(fd);
        if (!failure.empty()) {
            error = failure;
            return 0;
        }
        return total;
    }

    // Decompresses frame `index` of a seekable file (opened as fd) into out (vectors of char).
    template <typename In, typename Out>
    inline std::string ReadFrame(int fd, SeekTable const& table, size_t index, In& input, Out& out, ZSTD_DCtx* dctx) {
        auto const& frame = table.frames[index];
        std::string const where = "frame " + std::to_string(index + 1) + " of " + std::to_string(table.frames.size());
        input.resize(frame.compressed_size);
        if (!PreadAll(fd, input.data(), frame.compressed_size, frame.compressed_offset)) return where + ": read error";
        unsigned long long const content = ZSTD_getFrameContentSize(input.data(), input.size());
        if (content == ZSTD_CONTENTSIZE_ERROR) return where + ": not a zstd frame";
        if (content != ZSTD_CONTENTSIZE_UNKNOWN && content != frame.decompressed_size) {
            return where + ": holds " + std::to_string(content) + " bytes, the seek table says " +
                   std::to_string(frame.decompressed_size);
        }
        out.resize(frame.decompressed_size);
        size_t const got = ZSTD_decompressDCtx(dctx, out.data(), out.size(), input.data(), input.size());
        if (ZSTD_isError(got)) return where + ": " + ZSTD_getErrorName(got);
        if (got != out.size()) return where + ": " + std::to_string(got) + " bytes, the seek table says " + std::to_string(out.size());
        return "";
    }

    // What the frames of a load write into memory that becomes resident only as it is written (calloc'd arrays: the
    // index's key map and values), for ForEachFrame to stop workers near the end: output[i] is the bytes frame
    // first + i writes. Without it, the workers' buffers sit on top of an output that is nearly all resident when the
    // load ends: its peak.
    struct LoadBudget {
        std::vector<uint64_t> output;
        // What the workers' buffers may hold beyond that: memory the run takes later anyway (a query run's alignment), so
        // that buffers this size do not raise its peak. More workers decode to the end.
        uint64_t allowance = 0;
        size_t frames_on_fewer = 0;  // set by ForEachFrame: the frames started after a worker had stopped for memory
    };

    // Calls handle(index, data, size, worker) for each frame of a seekable file from `first` on, decompressed, on up
    // to `threads` threads (worker in [0, WorkerCount)); handle returns an error message, empty on success; the first
    // failure stops the workers. A worker's buffers (its compressed and decompressed frame) are freed, back to the
    // system, when it stops. With a budget, a worker takes another frame only while the output of the frames not yet
    // started, plus the budget's allowance, is at least the buffers of the workers decoding (each counted at the
    // largest frame's, compressed and decompressed), else it stops: near the end of the load the workers stop one by
    // one, so that what is written and what is buffered stay within the output's full size plus the allowance and one
    // worker's buffers (the last worker always goes on; the last frame, after which nothing is left to write, starts
    // once the others are done). The frames still are all handled, each once; only the number decoded at a time changes.
    template<typename Handle>
    std::string ForEachFrame(std::string const& path, SeekTable const& table, size_t first, int threads, Handle&& handle,
                             LoadBudget* budget = nullptr) {
        int const fd = ::open(path.c_str(), O_RDONLY | O_CLOEXEC);
        if (fd < 0) return std::string("cannot open the file: ") + std::strerror(errno);
        size_t const units = table.frames.size() > first ? table.frames.size() - first : 0;
        size_t const workers = WorkerCount(units, threads);
        // rest[i]: the output of the frames from first + i on; buffer: what a decoding worker holds at most.
        std::vector<uint64_t> rest(units + 1, 0);
        uint64_t buffer = 0;
        if (budget) {
            budget->frames_on_fewer = 0;
            for (size_t i = units; i-- > 0;) rest[i] = rest[i + 1] + (i < budget->output.size() ? budget->output[i] : 0);
            for (size_t i = 0; i < units; i++) {
                buffer = std::max(buffer, table.frames[first + i].compressed_size + table.frames[first + i].decompressed_size);
            }
        }
        int const window_max = ZSTD_dParam_getBounds(ZSTD_d_windowLogMax).upperBound;
        std::mutex mutex;
        size_t next = 0, active = workers, fewer = 0;
        bool failed = false;
        std::string first_error;
        auto run = [&](size_t worker) {
            ZSTD_DCtx* dctx = nullptr;
            FrameBuffer input, output;  // unmapped when this worker stops
            while (true) {
                size_t i;
                {
                    std::lock_guard<std::mutex> lock(mutex);
                    bool const over = budget && active > 1 && active * buffer > rest[next + 1 <= units ? next + 1 : units] + budget->allowance;
                    if (failed || next >= units || over) {
                        active--;
                        break;
                    }
                    if (active < workers) fewer++;
                    i = next++;
                }
                std::string error;
                try {
                    if (!dctx) {
                        dctx = ZSTD_createDCtx();
                        if (!dctx) error = "cannot allocate a zstd decompression context";
                        else ZSTD_DCtx_setParameter(dctx, ZSTD_d_windowLogMax, window_max);
                    }
                    if (error.empty()) error = ReadFrame(fd, table, first + i, input, output, dctx);
                    if (error.empty()) error = handle(first + i, output.data(), output.size(), worker);
                } catch (std::exception const& e) {
                    error = e.what();
                }
                if (!error.empty()) {
                    std::lock_guard<std::mutex> lock(mutex);
                    if (first_error.empty()) first_error = error;
                    failed = true;
                }
            }
            if (dctx) ZSTD_freeDCtx(dctx);
        };
        std::vector<std::thread> pool;
        for (size_t t = 1; t < workers; t++) pool.emplace_back(run, t);
        run(0);
        for (auto& t : pool) t.join();
        ::close(fd);
        if (budget) budget->frames_on_fewer = fewer;
        return first_error;
    }

    // Sequential stream buffer over the frames `table` lists in a file (e.g. one member of a
    // single-file database, Database.h), decompressed one frame at a time. Errors are printed once
    // and thrown from the read call, as in IStreambuf.
    class FrameStreambuf : public std::streambuf {
    public:
        FrameStreambuf(std::string path, SeekTable table, std::string name) :
                m_path(std::move(path)), m_name(std::move(name)), m_table(std::move(table)) {
            m_fd = ::open(m_path.c_str(), O_RDONLY | O_CLOEXEC);
            if (m_fd < 0) {
                m_open_error = std::strerror(errno);
                return;
            }
            m_dctx = ZSTD_createDCtx();
            if (!m_dctx) {
                m_open_error = "cannot allocate a zstd decompression context";
                return;
            }
            ZSTD_DCtx_setParameter(m_dctx, ZSTD_d_windowLogMax, ZSTD_dParam_getBounds(ZSTD_d_windowLogMax).upperBound);
        }

        ~FrameStreambuf() override {
            if (m_dctx) ZSTD_freeDCtx(m_dctx);
            if (m_fd >= 0) ::close(m_fd);
        }

        FrameStreambuf(FrameStreambuf const&) = delete;
        FrameStreambuf& operator=(FrameStreambuf const&) = delete;

        bool IsOpen() const { return m_open_error.empty(); }
        std::string const& OpenError() const { return m_open_error; }

    protected:
        int_type underflow() override {
            if (gptr() < egptr()) return traits_type::to_int_type(*gptr());
            while (IsOpen() && m_next < m_table.frames.size()) {
                size_t const index = m_next++;
                std::string const error = ReadFrame(m_fd, m_table, index, m_input, m_frame, m_dctx);
                if (!error.empty()) {
                    std::cerr << "Error reading " << m_name << ": " << error << std::endl;
                    throw std::runtime_error(m_name + ": " + error);
                }
                m_start = m_table.frames[index].decompressed_offset;
                setg(m_frame.data(), m_frame.data(), m_frame.data() + m_frame.size());
                if (!m_frame.empty()) return traits_type::to_int_type(m_frame[0]);
            }
            return traits_type::eof();
        }

        // Only reports the position (tellg); the stream cannot seek.
        pos_type seekoff(off_type off, std::ios_base::seekdir dir, std::ios_base::openmode) override {
            if (off != 0 || dir != std::ios_base::cur) return pos_type(off_type(-1));
            return pos_type(static_cast<off_type>(m_start + static_cast<uint64_t>(gptr() - eback())));
        }

        pos_type seekpos(pos_type pos, std::ios_base::openmode which) override {
            pos_type const here = seekoff(0, std::ios_base::cur, which);
            return pos == here ? here : pos_type(off_type(-1));
        }

    private:
        std::string m_path, m_name, m_open_error;
        SeekTable m_table;
        int m_fd = -1;
        ZSTD_DCtx* m_dctx = nullptr;
        size_t m_next = 0;
        uint64_t m_start = 0;  // content offset of the frame in the get area
        std::vector<char> m_input, m_frame;
    };

    // The content of some frames of a file (FrameStreambuf) as a std::istream.
    class FramesInput : public Input {
    public:
        FramesInput(std::string const& path, SeekTable table, std::string name) :
                m_buffer(path, std::move(table), std::move(name)), m_stream(&m_buffer) {
            if (!m_buffer.IsOpen()) m_stream.setstate(std::ios_base::badbit);
        }
        bool IsOpen() const override { return m_buffer.IsOpen(); }
        std::istream& Stream() override { return m_stream; }

    private:
        FrameStreambuf m_buffer;
        std::istream m_stream;
    };

    // A seekable zstd file's content read in order, its frames decompressed by `threads` threads of their
    // own ahead of the reader: at most `ahead` frames are decompressed or waiting to be read, each in a
    // buffer of its own. A frame that cannot be read or decompressed stops the threads; the reader gets the
    // frames before it, then the error, as from IStreambuf: printed once (Error() keeps it) and thrown from
    // the read call, which std::istream turns into badbit.
    class ParallelFrameStreambuf : public std::streambuf {
    public:
        ParallelFrameStreambuf(std::string path, SeekTable table, size_t threads, size_t ahead) :
                m_path(std::move(path)), m_table(std::move(table)), m_slots(std::max<size_t>(ahead, 1)) {
            m_fd = ::open(m_path.c_str(), O_RDONLY | O_CLOEXEC);
            if (m_fd < 0) {
                m_open_error = std::strerror(errno);
                return;
            }
            posix_fadvise(m_fd, 0, 0, POSIX_FADV_SEQUENTIAL);
            for (size_t t = 0; t < std::max<size_t>(threads, 1); t++) m_workers.emplace_back(&ParallelFrameStreambuf::Work, this);
        }

        ~ParallelFrameStreambuf() override {
            {
                std::lock_guard<std::mutex> lock(m_mutex);
                m_stop = true;
            }
            m_cv.notify_all();
            for (auto& worker : m_workers) worker.join();
            if (m_fd >= 0) ::close(m_fd);
        }

        ParallelFrameStreambuf(ParallelFrameStreambuf const&) = delete;
        ParallelFrameStreambuf& operator=(ParallelFrameStreambuf const&) = delete;

        bool IsOpen() const { return m_open_error.empty(); }
        std::string const& OpenError() const { return m_open_error; }
        // Why reading stopped early (as IStreambuf::Error); for the reading thread.
        std::string const& Error() const { return m_error; }

    protected:
        int_type underflow() override {
            if (gptr() < egptr()) return traits_type::to_int_type(*gptr());
            if (!IsOpen()) return traits_type::eof();
            std::unique_lock<std::mutex> lock(m_mutex);
            if (m_reading) {
                m_slots[m_released % m_slots.size()].ready = false;
                m_released++;
                m_reading = false;
                m_cv.notify_all();
            }
            setg(nullptr, nullptr, nullptr);
            while (m_released < m_table.frames.size()) {
                size_t const frame = m_released;
                Slot& slot = m_slots[frame % m_slots.size()];
                m_cv.wait(lock, [&] { return slot.ready && slot.frame == frame; });
                if (!slot.error.empty()) {
                    std::string const error = slot.error;
                    lock.unlock();
                    if (m_error.empty()) {
                        m_error = error;
                        std::cerr << "Error reading " << m_path << ": " << error << std::endl;
                    }
                    throw std::runtime_error(m_path + ": " + error);
                }
                if (slot.data.empty()) {  // a frame without content: the marker frame of protal's SAMs
                    slot.ready = false;
                    m_released++;
                    m_cv.notify_all();
                    continue;
                }
                m_reading = true;
                setg(slot.data.data(), slot.data.data(), slot.data.data() + slot.data.size());
                return traits_type::to_int_type(*gptr());
            }
            return traits_type::eof();
        }

    private:
        struct Slot {
            std::vector<char> data;
            size_t frame = SIZE_MAX;  // the frame it holds once ready
            bool ready = false;
            std::string error;
        };

        // Takes the next frame while it is fewer than m_slots frames ahead of the reader, and decompresses it
        // into its slot (the slot of the frame m_slots before it, which the reader has released).
        void Work() {
            ZSTD_DCtx* dctx = ZSTD_createDCtx();
            if (dctx) ZSTD_DCtx_setParameter(dctx, ZSTD_d_windowLogMax, ZSTD_dParam_getBounds(ZSTD_d_windowLogMax).upperBound);
            std::vector<char> input, output;
            while (true) {
                size_t frame;
                {
                    std::unique_lock<std::mutex> lock(m_mutex);
                    m_cv.wait(lock, [this] {
                        return m_stop || m_failed || m_claimed >= m_table.frames.size() || m_claimed < m_released + m_slots.size();
                    });
                    if (m_stop || m_failed || m_claimed >= m_table.frames.size()) break;
                    frame = m_claimed++;
                    std::swap(output, m_slots[frame % m_slots.size()].data);  // the slot's buffer, to reuse
                }
                std::string error = dctx ? Decompress(frame, dctx, input, output) : "cannot allocate a zstd decompression context";
                {
                    std::lock_guard<std::mutex> lock(m_mutex);
                    Slot& slot = m_slots[frame % m_slots.size()];
                    std::swap(slot.data, output);
                    slot.frame = frame;
                    slot.error = std::move(error);
                    slot.ready = true;
                    if (!slot.error.empty()) m_failed = true;
                }
                m_cv.notify_all();
            }
            if (dctx) ZSTD_freeDCtx(dctx);
        }

        std::string Decompress(size_t index, ZSTD_DCtx* dctx, std::vector<char>& input, std::vector<char>& output) const {
            auto const& frame = m_table.frames[index];
            input.resize(frame.compressed_size);
            if (!PreadAll(m_fd, input.data(), frame.compressed_size, frame.compressed_offset)) {
                return "read error at byte " + std::to_string(frame.compressed_offset);
            }
            output.resize(frame.decompressed_size);
            size_t const got = ZSTD_decompressDCtx(dctx, output.data(), output.size(), input.data(), input.size());
            if (ZSTD_isError(got)) return std::string(ZSTD_getErrorName(got)) + " (corrupt or not a zstd file?)";
            if (got != output.size()) {
                return "frame " + std::to_string(index + 1) + " holds " + std::to_string(got) + " bytes, the seek table says " +
                       std::to_string(output.size()) + " (corrupt file?)";
            }
            return {};
        }

        std::string m_path;
        std::string m_open_error;
        std::string m_error;
        SeekTable m_table;
        int m_fd = -1;
        std::vector<Slot> m_slots;
        size_t m_claimed = 0;   // frames taken by the threads
        size_t m_released = 0;  // frames the reader is done with
        bool m_reading = false; // the reader reads frame m_released
        bool m_stop = false;
        bool m_failed = false;
        std::mutex m_mutex;
        std::condition_variable m_cv;
        std::vector<std::thread> m_workers;
    };

    // The seek table (a skippable frame) of the frames whose sizes `entries` lists: per frame its compressed and its
    // content size, 4 bytes each, little-endian.
    inline std::string SeekTableFrame(std::string const& entries) {
        std::string seek;
        PutLE32(seek, kSeekTableFrameMagic);
        PutLE32(seek, static_cast<uint32_t>(entries.size() + 9));
        seek += entries;
        PutLE32(seek, static_cast<uint32_t>(entries.size() / 8));
        seek.push_back('\0');  // descriptor: no per-frame checksums in the table (each frame has one)
        PutLE32(seek, kSeekTableMagic);
        return seek;
    }

    // Writes a file in the seekable format: zstd frames (compressed already) in order, then the
    // seek table (Finish). Failures make Add and Finish return false, with the reason in Error().
    class FrameWriter {
    public:
        explicit FrameWriter(std::string const& path) {
            m_out = std::fopen(path.c_str(), "wb");
            if (!m_out) m_error = std::string("cannot open for writing: ") + std::strerror(errno);
        }
        ~FrameWriter() {
            if (m_out) std::fclose(m_out);
        }
        FrameWriter(FrameWriter const&) = delete;
        FrameWriter& operator=(FrameWriter const&) = delete;

        bool Ok() const { return m_error.empty(); }
        std::string const& Error() const { return m_error; }
        uint64_t Frames() const { return m_frames; }
        uint64_t Written() const { return m_written; }

        // Appends a zstd frame of `size` bytes whose content has `content_size` bytes.
        bool Add(char const* data, size_t size, uint64_t content_size) {
            if (!Ok()) return false;
            if (size > 0xffffffffu || content_size > 0xffffffffu) return Fail("a frame is larger than 4 GB");
            if (m_frames == 0xffffffffu) return Fail("more than 2^32 frames");
            if (std::fwrite(data, 1, size, m_out) != size) return Fail(std::string("write failed: ") + std::strerror(errno));
            PutLE32(m_table, static_cast<uint32_t>(size));
            PutLE32(m_table, static_cast<uint32_t>(content_size));
            m_written += size;
            m_frames++;
            return true;
        }

        // Writes the seek table and closes the file.
        bool Finish() {
            if (m_out && Ok()) {
                std::string const seek = SeekTableFrame(m_table);
                if (std::fwrite(seek.data(), 1, seek.size(), m_out) != seek.size()) Fail(std::string("write failed: ") + std::strerror(errno));
                m_written += seek.size();
            }
            if (m_out && std::fclose(m_out) != 0) Fail(std::string("closing failed: ") + std::strerror(errno));
            m_out = nullptr;
            return Ok();
        }

        bool Fail(std::string const& what) {
            if (m_error.empty()) m_error = what;
            return false;
        }

    private:
        std::FILE* m_out = nullptr;
        std::string m_error;
        std::string m_table;
        uint64_t m_frames = 0, m_written = 0;
    };

    // A frame's compressed bytes, room for ZSTD_compressBound of its content: not zeroed, so only the
    // bytes zstd writes become resident (a vector's resize writes every byte of the bound, of which a
    // frame needs half or less: 64 threads held ~1 GB of zeros writing an r226 index, ~3 GB packing its reference).
    struct CompressedBuffer {
        std::unique_ptr<char[]> data;
        size_t capacity = 0;
        size_t size = 0;

        char* Reserve(size_t bytes) {
            if (bytes > capacity) {
                data.reset(new char[bytes]);
                capacity = bytes;
            }
            return data.get();
        }
    };

    // A compression context for frames with content checksums; nullptr (and error set) if params
    // are invalid.
    inline ZSTD_CCtx* MakeCCtx(Params const& params, std::string& error) {
        ZSTD_CCtx* cctx = ZSTD_createCCtx();
        bool ok = cctx && !ZSTD_isError(ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, params.level)) &&
                  !ZSTD_isError(ZSTD_CCtx_setParameter(cctx, ZSTD_c_checksumFlag, 1));
        if (ok && params.window_log > 0) {
            ok = !ZSTD_isError(ZSTD_CCtx_setParameter(cctx, ZSTD_c_enableLongDistanceMatching, 1)) &&
                 !ZSTD_isError(ZSTD_CCtx_setParameter(cctx, ZSTD_c_windowLog, params.window_log));
        }
        if (!ok) {
            error = "cannot set up zstd compression (level " + std::to_string(params.level) + ", window log " +
                    std::to_string(params.window_log) + ")";
            if (cctx) ZSTD_freeCCtx(cctx);
            return nullptr;
        }
        return cctx;
    }

    // Writes `count` frames to path in the seekable format: make(index, out) fills frame index's
    // content (an error message, empty on success); up to params.threads frames are made and
    // compressed at once, and written in order, followed by the seek table. Returns the bytes
    // written, or nullopt with a message in error.
    template<typename Make>
    std::optional<uint64_t> WriteSeekable(std::string const& path, size_t count, Params const& params, Make&& make,
                                          std::string& error) {
        int const threads = std::max(1, params.threads);
        struct Slot {
            ZSTD_CCtx* cctx = nullptr;
            std::vector<char> in;
            CompressedBuffer out;
            Slot() = default;
            Slot(Slot const&) = delete;
            ~Slot() { if (cctx) ZSTD_freeCCtx(cctx); }
        };
        std::vector<Slot> slots(static_cast<size_t>(threads));
        for (auto& slot : slots) {
            if (!(slot.cctx = MakeCCtx(params, error))) return std::nullopt;
        }
        FrameWriter out(path);
        if (!out.Ok()) {
            error = out.Error();
            return std::nullopt;
        }
        for (size_t batch_start = 0; batch_start < count && error.empty(); batch_start += slots.size()) {
            size_t const batch = std::min(slots.size(), count - batch_start);
            error = ParallelFor(batch, threads, [&](size_t b, size_t) -> std::string {
                Slot& slot = slots[b];
                slot.in.clear();
                std::string e = make(batch_start + b, slot.in);
                if (!e.empty()) return e;
                if (slot.in.size() > 0xffffffffu) return "frame " + std::to_string(batch_start + b) + " is larger than 4 GB";
                size_t const bound = ZSTD_compressBound(slot.in.size());
                size_t const r = ZSTD_compress2(slot.cctx, slot.out.Reserve(bound), bound, slot.in.data(), slot.in.size());
                if (ZSTD_isError(r)) return ZSTD_getErrorName(r);
                slot.out.size = r;
                return "";
            });
            for (size_t b = 0; b < batch && error.empty(); b++) {
                if (!out.Add(slots[b].out.data.get(), slots[b].out.size, slots[b].in.size())) error = out.Error();
            }
        }
        if (error.empty() && !out.Finish()) error = out.Error();
        if (!error.empty()) return std::nullopt;
        return out.Written();
    }

    // A source of data to compress, read in order.
    class Reader {
    public:
        virtual ~Reader() = default;
        // Up to size bytes into dst; returns the number read, 0 at the end.
        virtual size_t Read(char* dst, size_t size) = 0;
        virtual bool Ok() const { return true; }
    };

    class StreamReader : public Reader {
    public:
        explicit StreamReader(std::istream& is) : m_is(is) {}
        size_t Read(char* dst, size_t size) override {
            m_is.read(dst, static_cast<std::streamsize>(size));
            return static_cast<size_t>(m_is.gcount());
        }
        bool Ok() const override { return !m_is.bad(); }

    private:
        std::istream& m_is;
    };

    // Compresses reader's data into out: frames of params.frame_size bytes, each compressed on its
    // own with a content checksum, params.threads frames at a time. Returns false with a message in
    // error.
    inline bool CompressFramesTo(Reader& reader, FrameWriter& out, Params const& params, std::string& error) {
        size_t const frame_size = static_cast<size_t>(params.frame_size);
        if (frame_size == 0 || frame_size > 0xffffffffu) {
            error = "frame size must be between 1 byte and 4 GB";
            return false;
        }
        int const threads = std::max(1, params.threads);
        struct Slot {
            ZSTD_CCtx* cctx = nullptr;
            std::vector<char> in;
            CompressedBuffer out;
            size_t in_size = 0;
            ~Slot() { if (cctx) ZSTD_freeCCtx(cctx); }
        };
        std::vector<Slot> slots(static_cast<size_t>(threads));
        for (auto& slot : slots) {
            if (!(slot.cctx = MakeCCtx(params, error))) return false;
        }
        bool end = false;
        while (!end) {
            // Read up to `threads` frames in order, compress them in parallel, write them in order.
            size_t batch = 0;
            while (batch < slots.size() && !end) {
                Slot& slot = slots[batch];
                slot.in.resize(frame_size);
                size_t n = 0;
                while (n < frame_size) {
                    size_t const got = reader.Read(slot.in.data() + n, frame_size - n);
                    if (got == 0) break;
                    n += got;
                }
                if (n > 0) {
                    slot.in_size = n;
                    batch++;
                }
                if (n < frame_size) end = true;
            }
            if (!reader.Ok()) {
                error = "reading the data to compress failed";
                break;
            }
            error = ParallelFor(batch, threads, [&](size_t b, size_t) -> std::string {
                Slot& slot = slots[b];
                size_t const bound = ZSTD_compressBound(slot.in_size);
                size_t const r = ZSTD_compress2(slot.cctx, slot.out.Reserve(bound), bound, slot.in.data(), slot.in_size);
                if (ZSTD_isError(r)) return ZSTD_getErrorName(r);
                slot.out.size = r;
                return "";
            });
            if (!error.empty()) break;
            for (size_t b = 0; b < batch; b++) {
                if (!out.Add(slots[b].out.data.get(), slots[b].out.size, slots[b].in_size)) {
                    error = out.Error();
                    break;
                }
            }
            if (!error.empty()) break;
        }
        return error.empty();
    }

    // Writes reader's data to path in the seekable format (CompressFramesTo), then the seek table.
    // Returns the number of bytes written, or nullopt with a message in error.
    inline std::optional<uint64_t> CompressFrames(Reader& reader, std::string const& path, Params const& params,
                                                  std::string& error) {
        FrameWriter out(path);
        if (!out.Ok()) {
            error = out.Error();
            return std::nullopt;
        }
        if (!CompressFramesTo(reader, out, params, error)) return std::nullopt;
        if (!out.Finish()) {
            error = out.Error();
            return std::nullopt;
        }
        return out.Written();
    }

    // Size of the file's content: the file size for a raw file; for a zstd file the total of its
    // seek table, else the content size in the frame header (single-frame files from the zstd CLI),
    // else counted by decompressing. nullopt if the file cannot be read.
    inline std::optional<uint64_t> UncompressedSize(std::string const& path) {
        std::error_code ec;
        if (!IsCompressed(path)) {
            uint64_t const size = std::filesystem::file_size(path, ec);
            return ec ? std::nullopt : std::optional<uint64_t>(size);
        }
        {
            std::string error;
            if (auto table = ReadSeekTable(path, error)) return table->DecompressedSize();
        }
        {
            unsigned char header[18];  // the largest zstd frame header
            std::ifstream is(path, std::ios::binary);
            is.read(reinterpret_cast<char*>(header), sizeof(header));
            // Only a regular first frame tells the size (a skippable one reports 0).
            bool const regular = is.gcount() >= 4 && header[0] == 0x28 && header[1] == 0xb5 && header[2] == 0x2f &&
                                 header[3] == 0xfd;
            unsigned long long const size = ZSTD_getFrameContentSize(header, static_cast<size_t>(is.gcount()));
            if (regular && size != ZSTD_CONTENTSIZE_UNKNOWN && size != ZSTD_CONTENTSIZE_ERROR) return size;
        }
        InputFile in(path);
        if (!in.IsOpen()) return std::nullopt;
        std::vector<char> buffer(size_t{8} << 20);
        uint64_t total = 0;
        while (in.Stream().read(buffer.data(), static_cast<std::streamsize>(buffer.size())) || in.Stream().gcount() > 0) {
            total += static_cast<uint64_t>(in.Stream().gcount());
            if (!in.Stream()) break;
        }
        return in.Stream().bad() ? std::nullopt : std::optional<uint64_t>(total);
    }

    // Compresses src (raw, or zstd to recompress it) into dst via dst.partial, renamed at the end:
    // in the seekable format if params.frame_size > 0, else as one frame. dst may be src. With
    // verify, the new file is decompressed and compared with src first. Returns false with a
    // message in error.
    inline bool CompressFile(std::string const& src, std::string const& dst, Params const& params, bool verify,
                             std::string& error) {
        std::error_code ec;
        auto const content_size = UncompressedSize(src);
        if (!content_size) {
            error = "cannot read " + src;
            return false;
        }
        uint64_t const size = *content_size;
        std::string const partial = dst + ".partial";
        auto fail = [&](std::string const& what) {
            error = what;
            std::filesystem::remove(partial, ec);
            return false;
        };
        std::vector<char> buffer(size_t{8} << 20);
        {
            InputFile in(src);
            if (!in.IsOpen()) return fail("cannot open " + src);
            if (params.frame_size > 0) {
                StreamReader reader(in.Stream());
                std::string frame_error;
                if (!CompressFrames(reader, partial, params, frame_error)) return fail("writing " + partial + " failed: " + frame_error);
            } else {
                OStream out(partial, params, size);
                while (in.Stream() && out) {
                    in.Stream().read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
                    out.write(buffer.data(), in.Stream().gcount());
                }
                if (!out.Close()) return fail("writing " + partial + " failed: " + out.Buffer().Error());
            }
            if (in.Stream().bad() || !in.Stream().eof()) return fail("reading " + src + " failed");
        }
        if (verify) {
            InputFile original(src);
            InputFile copy(partial);
            std::vector<char> decompressed(buffer.size());
            uint64_t compared = 0;
            while (true) {
                original.Stream().read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
                copy.Stream().read(decompressed.data(), static_cast<std::streamsize>(decompressed.size()));
                if (original.Stream().bad()) return fail("verifying " + partial + " failed: cannot re-read " + src);
                if (copy.Stream().bad()) return fail("verifying " + partial + " failed: it cannot be decompressed");
                std::streamsize const a = original.Stream().gcount(), b = copy.Stream().gcount();
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
