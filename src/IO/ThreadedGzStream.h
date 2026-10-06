#pragma once

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zstd.h>

#include <algorithm>
#include <cerrno>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <deque>
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
    class ThreadedGzStreambuf : public std::streambuf {
    public:
        static constexpr size_t kBlockSize = size_t{1} << 20;  // inflated bytes per block
        static constexpr size_t kBlocks = 4;
        static constexpr uint32_t kInputBuffer = 1u << 17;     // compressed bytes per read()

        ThreadedGzStreambuf() = default;
        ThreadedGzStreambuf(ThreadedGzStreambuf const&) = delete;
        ThreadedGzStreambuf& operator=(ThreadedGzStreambuf const&) = delete;
        ~ThreadedGzStreambuf() override { close(); }

        // False if the file cannot be read; open_error() then says why.
        bool open(char const* path) {
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
            m_bgzf = S_ISREG(st.st_mode) && ::pread(m_fd, head, sizeof(head), 0) == static_cast<ssize_t>(sizeof(head)) &&
                     bgzf::IsBlockHeader(head);
            m_offset = 0;
            m_members = 0;
            if (m_bgzf) {
                m_block_in.resize(bgzf::kMaxBlock);
                m_last_empty = false;
            } else if (!StartGzip(0, m_open_error)) {
                ::close(m_fd);
                m_fd = -1;
                return false;
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
                end = m_bgzf ? FillBgzf(data, size, error) : FillStream(data, size, error);
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
