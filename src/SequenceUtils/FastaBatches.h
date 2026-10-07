#pragma once

// A FASTA file in batches of whole records, for the index build's parallel passes (Build.h): one
// thread reads raw bytes and cuts them at record starts, and the thread that takes a batch parses
// its records (ForEachFastaRecord), with what SeqReader gives for FASTA. A seekable zstd file whose
// frames each begin at a record (FastaFrames: full_reference.fna.zst as gtdb_to_protal_db.py writes it)
// is read on several threads at once, a frame each.

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cstring>
#include <istream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <string_view>
#include <vector>
#include <fcntl.h>
#include <unistd.h>
#include <zstd.h>
#include "Zstd.h"

namespace protal {

    class FastaBatches {
    public:
        explicit FastaBatches(std::istream& is) : m_is(is) {}

        // The next whole records into batch: those that start within its first `bytes` (so about
        // `bytes`, more by the rest of the last record); false at the end of the input or on an
        // error (Error()).
        bool Next(std::string& batch, size_t bytes) {
            batch.swap(m_rest);
            m_rest.clear();
            if (!m_error.empty()) return false;
            bytes = std::max<size_t>(bytes, 1);
            size_t want = bytes;
            while (true) {
                if (!m_end && batch.size() < want) Read(batch, want - batch.size());
                if (!m_error.empty()) return false;
                if (batch.empty()) return false;
                if (!m_started) {
                    m_started = true;
                    if (batch[0] != '>') {
                        m_error = "a FASTA file starts with '>', this one with '" + std::string(1, batch[0]) + "'";
                        return false;
                    }
                }
                // Cut before the first record that starts at or after `bytes` (a line end followed by
                // '>'): the batch holds `bytes` and the rest of the record there.
                if (batch.size() > bytes) {
                    size_t const cut = batch.find("\n>", bytes - 1);
                    if (cut != std::string::npos) {
                        m_rest.assign(batch, cut + 1, std::string::npos);
                        batch.resize(cut + 1);
                        return true;
                    }
                }
                if (m_end) return true;  // the rest of the file, whole records
                want = std::max(want, batch.size()) * 2;  // that record goes on: read on
            }
        }

        std::string const& Error() const { return m_error; }

    private:
        void Read(std::string& batch, size_t size) {
            size_t const old = batch.size();
            batch.resize(old + size);
            m_is.read(batch.data() + old, static_cast<std::streamsize>(size));
            size_t const got = static_cast<size_t>(m_is.gcount());
            batch.resize(old + got);
            if (got < size) m_end = true;
            // A read or decompression error is not the end of the file: the pass would go on without the rest of it.
            if (m_is.bad()) m_error = "a read error (a truncated or corrupt file?)";
        }

        std::istream& m_is;
        std::string m_rest;  // read after the last whole record of the previous batch
        bool m_end = false;
        bool m_started = false;
        std::string m_error;
    };

    // The line without its trailing white space (as BufferedFastxReader's StripString).
    inline std::string_view StrippedLine(std::string_view line) {
        while (!line.empty() && std::isspace(static_cast<unsigned char>(line.back()))) line.remove_suffix(1);
        return line;
    }

    // Calls f(header, sequence) for each record of batch (whole records, as FastaBatches gives them):
    // the header line, and the sequence's lines joined, each without trailing white space, as
    // BufferedFastxReader reads FASTA. A one-line sequence is a view into batch; several lines are
    // joined in scratch.
    template<typename F>
    void ForEachFastaRecord(std::string_view batch, std::string& scratch, F&& f) {
        size_t pos = 0;
        size_t const n = batch.size();
        auto line_end = [&](size_t from) {
            size_t const end = batch.find('\n', from);
            return end == std::string_view::npos ? n : end;
        };
        while (pos < n) {
            size_t end = line_end(pos);
            std::string_view const header = StrippedLine(batch.substr(pos, end - pos));
            pos = end < n ? end + 1 : n;
            std::string_view single;
            size_t lines = 0;
            while (pos < n && batch[pos] != '>') {
                end = line_end(pos);
                std::string_view const line = StrippedLine(batch.substr(pos, end - pos));
                if (lines == 0) {
                    single = line;
                } else {
                    if (lines == 1) scratch.assign(single);
                    scratch.append(line);
                }
                lines++;
                pos = end < n ? end + 1 : n;
            }
            f(header, lines <= 1 ? single : std::string_view(scratch));
        }
    }

    // One zstd frame of a file (at `offset`, `size` bytes), decompressed as a stream: the frame's window (128 MB for
    // gtdb_to_protal_db.py's --long=27) and two small buffers, never the whole content. A decompression or read
    // error is thrown from the read call, which an std::istream turns into badbit (FastaBatches then reports it).
    class FrameStreambuf final : public std::streambuf {
    public:
        static constexpr size_t kInput = size_t{1} << 20, kOutput = size_t{1} << 20;

        // output, input: the bytes decompressed and read at a time (a few and a block's worth for a look at the frame's first
        // byte: then only its first block is read and decompressed).
        FrameStreambuf(int fd, ZSTD_DCtx* dctx, uint64_t offset, uint64_t size, uint64_t content, size_t output = kOutput,
                       size_t input = kInput)
                : m_fd(fd), m_dctx(dctx), m_next(offset), m_left(size), m_content(content), m_in(std::max<size_t>(input, 1)),
                  m_out(std::max<size_t>(output, 1)) {
            ZSTD_DCtx_reset(m_dctx, ZSTD_reset_session_only);
            m_input = { m_in.data(), 0, 0 };
            setg(m_out.data(), m_out.data(), m_out.data());
        }

        // What went wrong, if the read call threw (the istream keeps only badbit).
        std::string const& Failure() const { return m_failure; }

    protected:
        int_type underflow() override {
            if (gptr() < egptr()) return traits_type::to_int_type(*gptr());
            while (true) {
                if (m_input.pos == m_input.size && m_left > 0) {
                    size_t const n = static_cast<size_t>(std::min<uint64_t>(m_left, m_in.size()));
                    if (!zstd::PreadAll(m_fd, m_in.data(), n, m_next)) Fail("read error");
                    m_input = { m_in.data(), n, 0 };
                    m_next += n;
                    m_left -= n;
                }
                if (m_done) {
                    if (m_input.pos < m_input.size || m_left > 0) Fail("data after the end of the zstd frame");
                    if (m_produced != m_content) {
                        Fail("the frame holds " + std::to_string(m_produced) + " bytes, the seek table says " + std::to_string(m_content));
                    }
                    return traits_type::eof();
                }
                ZSTD_outBuffer output{ m_out.data(), m_out.size(), 0 };
                size_t const r = ZSTD_decompressStream(m_dctx, &output, &m_input);
                if (ZSTD_isError(r)) Fail(std::string("zstd: ") + ZSTD_getErrorName(r));
                if (r == 0) m_done = true;  // the frame is complete
                if (output.pos > 0) {
                    m_produced += output.pos;
                    setg(m_out.data(), m_out.data(), m_out.data() + output.pos);
                    return traits_type::to_int_type(*gptr());
                }
                if (!m_done && m_input.pos == m_input.size && m_left == 0) Fail("the zstd frame is cut short");
            }
        }

    private:
        [[noreturn]] void Fail(std::string what) {
            m_failure = what;
            throw std::runtime_error(what);
        }

        int m_fd;
        ZSTD_DCtx* m_dctx;
        uint64_t m_next, m_left, m_content, m_produced = 0;
        bool m_done = false;
        std::vector<char> m_in, m_out;
        ZSTD_inBuffer m_input{};
        std::string m_failure;
    };

    // A FASTA file in frames that each begin at a record: a seekable zstd file (its seek table, zstd's seekable format)
    // as gtdb_to_protal_db.py writes full_reference.fna.zst, a frame per marker file (or gene). Each frame is a FASTA
    // file of its own, so several threads read the file at once, each its frames (Reader), instead of one thread
    // decompressing for all: the uniqueness check and gene conservation of --build read 86 GB at GTDB r226.
    class FastaFrames {
    public:
        // The frames of the file at path; nullopt if it is not a seekable zstd file of two frames or more whose every
        // non-empty frame begins with '>' (read it as one stream then), with error set if it cannot be read at all.
        static std::optional<FastaFrames> Open(std::string const& path, std::string& error) {
            if (!zstd::IsCompressed(path)) return std::nullopt;
            auto table = zstd::ReadSeekTable(path, error);
            if (!table || table->frames.size() < 2) return std::nullopt;
            FastaFrames frames;
            frames.m_path = path;
            frames.m_table = std::move(*table);
            // Each frame must begin a record: its first byte (the first block decompressed, no more).
            int const fd = ::open(path.c_str(), O_RDONLY | O_CLOEXEC);
            if (fd < 0) {
                error = std::string("cannot open ") + path + ": " + std::strerror(errno);
                return std::nullopt;
            }
            ZSTD_DCtx* dctx = MakeDCtx();
            bool aligned = dctx != nullptr;
            try {
                for (auto const& frame : frames.m_table.frames) {
                    if (!aligned || frame.decompressed_size == 0) continue;
                    FrameStreambuf buffer(fd, dctx, frame.compressed_offset, frame.compressed_size, frame.decompressed_size, 16,
                                          size_t{256} << 10);
                    aligned = buffer.sgetc() == '>';
                }
            } catch (std::exception const& e) {
                error = path + ": " + e.what();
                aligned = false;
            }
            if (dctx) ZSTD_freeDCtx(dctx);
            ::close(fd);
            if (!aligned) return std::nullopt;
            return frames;
        }

        size_t Frames() const { return m_table.frames.size(); }
        std::string const& Path() const { return m_path; }

        static ZSTD_DCtx* MakeDCtx() {
            ZSTD_DCtx* dctx = ZSTD_createDCtx();
            if (dctx) ZSTD_DCtx_setParameter(dctx, ZSTD_d_windowLogMax, ZSTD_dParam_getBounds(ZSTD_d_windowLogMax).upperBound);
            return dctx;
        }

        // One thread's reader: frames taken one after another (Take), each in batches of whole records (Next).
        class Reader {
        public:
            explicit Reader(FastaFrames const& frames) : m_frames(&frames) {
                m_fd = ::open(frames.m_path.c_str(), O_RDONLY | O_CLOEXEC);
                m_dctx = MakeDCtx();
                if (m_fd < 0 || !m_dctx) m_error = "cannot open " + frames.m_path;
            }
            ~Reader() {
                if (m_fd >= 0) ::close(m_fd);
                if (m_dctx) ZSTD_freeDCtx(m_dctx);
            }
            Reader(Reader const&) = delete;
            Reader& operator=(Reader const&) = delete;

            // Starts frame f (the previous one, if any, is dropped).
            bool Take(size_t f) {
                m_batches.reset();
                m_stream.reset();
                m_buffer.reset();
                if (!m_error.empty()) return false;
                auto const& frame = m_frames->m_table.frames[f];
                m_buffer = std::make_unique<FrameStreambuf>(m_fd, m_dctx, frame.compressed_offset, frame.compressed_size, frame.decompressed_size);
                m_stream = std::make_unique<std::istream>(m_buffer.get());
                m_batches = std::make_unique<FastaBatches>(*m_stream);
                m_frame = f;
                return true;
            }

            // The next whole records of the current frame into batch (about `bytes`, FastaBatches::Next); false at its
            // end or on an error (Error()).
            bool Next(std::string& batch, size_t bytes) {
                if (!m_batches) return false;
                if (m_batches->Next(batch, bytes)) return true;
                if (!m_batches->Error().empty() && m_error.empty()) {
                    m_error = m_frames->m_path + ", frame " + std::to_string(m_frame + 1) + " of " + std::to_string(m_frames->Frames()) +
                              ": " + (m_buffer->Failure().empty() ? m_batches->Error() : m_buffer->Failure());
                }
                m_batches.reset();
                return false;
            }

            std::string const& Error() const { return m_error; }

        private:
            FastaFrames const* m_frames;
            int m_fd = -1;
            ZSTD_DCtx* m_dctx = nullptr;
            size_t m_frame = 0;
            std::unique_ptr<FrameStreambuf> m_buffer;
            std::unique_ptr<std::istream> m_stream;
            std::unique_ptr<FastaBatches> m_batches;
            std::string m_error;
        };

    private:
        std::string m_path;
        zstd::SeekTable m_table;
    };
}
