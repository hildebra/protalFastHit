// SamFile.h - SAM files as protal writes and reads them: plain (.sam), gzip (.sam.gz, written as
// BGZF blocks with ISA-L, Bgzf.h) or zstd (.sam.zst, written in zstd's seekable format),
// chosen by the file name.
//
// Writing (SamOutput): the output handlers of all alignment threads hand over blocks of whole
// reads' records together with the genes those records name. The thread that hands a block over
// compresses it, and only appending it to the file is serialised, so compression runs alongside
// alignment instead of in a pass of its own afterwards. The records collect in a temporary file;
// once they are complete, Finish writes the SAM: its header, which lists only the genes that the
// records name (SAM needs an @SQ line only for the references that records use, and a database has
// millions of genes), then the records. A header given up front (every gene, --full_sam_header) is
// written first and the records follow it directly. A .sam.zst can instead take its records
// straight into the SAM behind room left for the header, which a skippable frame pads (header_room),
// so that they are not copied: the room the header does not need is never written.
//
// Reading (SamInput): a plain, gzip (BGZF or not) or zstd file as a std::istream, with the checks
// that it is complete: the gzip readers' and zstd's own, and the end markers of the formats
// written here, which a file cut exactly at a block or frame boundary lacks: BGZF's end-of-file
// block, and the seek table of a .sam.zst that starts with protal's marker frame.
#pragma once

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zstd.h>

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <istream>
#include <memory>
#include <mutex>
#include <optional>
#include <ostream>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

#include "Bgzf.h"
#include "ThreadedGzStream.h"
#include "Zstd.h"

namespace protal {
    enum class SamCompression { None, Gzip, Zstd };

    // By the file name: ".gz" is gzip, ".zst" zstd, any other name plain SAM.
    inline SamCompression SamCompressionOf(std::string_view path) {
        if (path.ends_with(".gz")) return SamCompression::Gzip;
        if (path.ends_with(".zst")) return SamCompression::Zstd;
        return SamCompression::None;
    }

    // The name without its compression extension (".gz", ".zst").
    inline std::string UncompressedSamName(std::string const& path) {
        if (path.ends_with(".gz")) return path.substr(0, path.size() - 3);
        if (path.ends_with(".zst")) return path.substr(0, path.size() - 4);
        return path;
    }

    // A gene of the database as one number, as SAM records name it (<taxid>_<gene id>).
    inline uint64_t SamGeneKey(uint64_t taxid, uint64_t geneid) {
        return taxid << 32 | geneid;
    }

    inline std::pair<uint64_t, uint64_t> SamGeneOfKey(uint64_t key) {
        return { key >> 32, key & 0xffffffffu };
    }

    // ---- zstd ------------------------------------------------------------------------------------
    // zstd's seekable format (Zstd.h): independent frames with content checksums, then the seek
    // table. The file starts with a skippable frame (which zstd tools skip) that marks it as
    // protal's, so that a reader knows a seek table must end it. With room left for the header,
    // another skippable frame fills what the header does not need, up to the records.
    namespace sam_zstd {
        inline constexpr int kLevel = 3;                                   // zstd's default
        inline constexpr size_t kFrameInput = size_t{1} << 20;
        inline constexpr uint32_t kMarkerMagic = 0x184D2A50;               // skippable frame, variant 0
        inline constexpr uint32_t kPaddingMagic = 0x184D2A51;              // skippable frame, variant 1
        inline constexpr size_t kSkippableHeader = 8;                      // magic and size
        inline constexpr std::string_view kMarkerText = "protal SAM, zstd seekable format: a seek table ends the file";

        inline std::string MarkerFrame() {
            std::string frame;
            zstd::PutLE32(frame, kMarkerMagic);
            zstd::PutLE32(frame, static_cast<uint32_t>(kMarkerText.size()));
            frame += kMarkerText;
            return frame;
        }

        // Room for the marker and header of a SAM whose records name at most `genes` genes, aligned
        // from `input_bytes` of read files (0: unknown): 8 bytes a gene (a header compresses to about
        // 3.5 bytes per @SQ line) and no more than 1/16 of the input (headers measured at most 3% of
        // it), at least 16 KB and at most 1 MB, in whole 4 KB. With the unaligned reads' failed candidates
        // of up to `taxa` taxa in it (FailedCandidatesLine), 4 bytes more a taxon (about 3 compressed),
        // no more than 1/16 of the input either, at most 1 MB more. What the header does not need stays
        // unwritten; a header that needs more is written as without room (SamOutput).
        inline size_t HeaderRoom(size_t genes, uint64_t input_bytes = 0, size_t taxa = 0) {
            constexpr size_t kUnit = size_t{4} << 10;
            size_t room = genes * 8;
            if (input_bytes > 0) room = std::min<uint64_t>(room, input_bytes / 16);
            room = std::clamp<size_t>(room, size_t{16} << 10, size_t{1} << 20);
            size_t more = taxa * 4;
            if (input_bytes > 0) more = std::min<uint64_t>(more, input_bytes / 16);
            room += std::min<size_t>(more, size_t{1} << 20);
            return (room + kUnit - 1) / kUnit * kUnit;
        }

        using Frames = std::vector<std::pair<uint32_t, uint32_t>>;  // (compressed, content) bytes per frame

        // Appends [data, data + size) to out as zstd frames of at most kFrameInput content bytes each,
        // and their sizes to frames. False only if zstd fails.
        inline bool Compress(char const* data, size_t size, std::string& out, Frames& frames) {
            struct Context {
                ZSTD_CCtx* cctx = nullptr;
                Context() {
                    std::string error;
                    cctx = zstd::MakeCCtx(zstd::Params{ kLevel, 0, 1, 0 }, error);
                }
                ~Context() {
                    if (cctx) ZSTD_freeCCtx(cctx);
                }
            };
            thread_local Context context;
            if (!context.cctx) return false;
            for (size_t offset = 0; offset < size; offset += kFrameInput) {
                size_t const n = std::min(kFrameInput, size - offset);
                size_t const start = out.size();
                out.resize(start + ZSTD_compressBound(n));
                size_t const written = ZSTD_compress2(context.cctx, out.data() + start, out.size() - start, data + offset, n);
                if (ZSTD_isError(written)) return false;
                out.resize(start + written);
                frames.emplace_back(static_cast<uint32_t>(written), static_cast<uint32_t>(n));
            }
            return true;
        }

        inline bool StartsWithMarker(std::string const& path) {
            std::string const marker = MarkerFrame();
            std::ifstream is(path, std::ios::binary);
            std::string head(marker.size(), '\0');
            is.read(head.data(), static_cast<std::streamsize>(head.size()));
            return is.gcount() == static_cast<std::streamsize>(head.size()) && head == marker;
        }

        // The seek table frame for frames, in file order.
        inline std::string SeekTable(Frames const& frames) {
            std::string entries;
            for (auto const& [compressed, content] : frames) {
                zstd::PutLE32(entries, compressed);
                zstd::PutLE32(entries, content);
            }
            std::string table;
            zstd::PutLE32(table, zstd::kSeekTableFrameMagic);
            zstd::PutLE32(table, static_cast<uint32_t>(entries.size() + 9));
            table += entries;
            zstd::PutLE32(table, static_cast<uint32_t>(frames.size()));
            table.push_back('\0');  // descriptor: no per-frame checksums in the table (each frame has one)
            zstd::PutLE32(table, zstd::kSeekTableMagic);
            return table;
        }
    }

    // ---- Writing ---------------------------------------------------------------------------------

    // Where output handlers put SAM records: blocks of whole reads' records, each line ending in a
    // newline, from any thread at any time, with the genes they name (SamGeneKey of RNAME and
    // RNEXT; duplicates allowed). Write empties `genes`.
    // The reads that seeded on taxa but aligned nowhere go in as one unmapped record each (FLAG 4, ZF tag), or, unless
    // SetUnmappedRecords(true), as counts per taxon (AddFailedCandidates) for the header (FailedCandidatesLine).
    class SamSink {
    public:
        virtual ~SamSink() = default;
        virtual void Write(char const* data, size_t size, std::vector<uint64_t>& genes) = 0;

        // Whether output handlers write unaligned reads' unmapped records (--write_unmapped_reads) rather than count them.
        bool WritesUnmappedRecords() const { return m_unmapped_records; }
        void SetUnmappedRecords(bool write) { m_unmapped_records = write; }

        // Adds `counts` (reads per taxid) to the sink's and sets them to 0.
        void AddFailedCandidates(std::vector<uint32_t>& counts) {
            std::lock_guard<std::mutex> lock(m_failed_mutex);
            if (m_failed.size() < counts.size()) m_failed.resize(counts.size(), 0);
            for (size_t taxid = 0; taxid < counts.size(); taxid++) {
                m_failed[taxid] += counts[taxid];
                counts[taxid] = 0;
            }
        }

        // The unaligned reads per taxid that AddFailedCandidates counted.
        std::vector<uint64_t> FailedCandidates() const {
            std::lock_guard<std::mutex> lock(m_failed_mutex);
            return m_failed;
        }

    private:
        bool m_unmapped_records = true;
        mutable std::mutex m_failed_mutex;
        std::vector<uint64_t> m_failed;
    };

    // A SAM file written by several threads (see the top of this file). Errors (a file that cannot
    // be written, a full disk) are kept in Error(); later writes are then skipped.
    class SamOutput : public SamSink {
    public:
        static inline const std::string kRecordsSuffix = ".records.partial";

        // Writes the SAM `path`, compressed as `compression`. With a header, the header goes first
        // and the records follow; without one, the records collect in path + kRecordsSuffix until
        // Finish(header). For zstd, `header_room` bytes (sam_zstd::HeaderRoom) instead take the
        // records into `path` itself, after that much room for the marker and the header.
        SamOutput(std::string path, SamCompression compression, std::optional<std::string> const& header = std::nullopt,
                  size_t header_room = 0) :
                m_path(std::move(path)), m_records_path(m_path + kRecordsSuffix), m_compression(compression),
                m_header_first(header.has_value()) {
            if (m_header_first) {
                m_fd = Create(m_path);
                if (m_fd >= 0) WriteHead(*header);
            } else if (compression == SamCompression::Zstd && header_room >= sam_zstd::MarkerFrame().size() + sam_zstd::kSkippableHeader) {
                m_room = header_room;
                m_fd = Create(m_path);
                if (m_fd >= 0 && ::lseek(m_fd, static_cast<off_t>(m_room), SEEK_SET) < 0) {
                    Fail("cannot write " + m_path + ": " + std::strerror(errno));
                }
            } else {
                m_fd = Create(m_records_path);
            }
        }

        ~SamOutput() override {
            if (!m_finished) Discard();
        }

        SamOutput(SamOutput const&) = delete;
        SamOutput& operator=(SamOutput const&) = delete;

        bool Ok() const {
            std::lock_guard<std::mutex> lock(m_mutex);
            return m_error.empty();
        }

        std::string Error() const {
            std::lock_guard<std::mutex> lock(m_mutex);
            return m_error;
        }

        std::string const& Path() const { return m_path; }

        void Write(char const* data, size_t size, std::vector<uint64_t>& genes) override {
            // Compressed by the calling thread, outside the lock.
            thread_local std::string packed;
            thread_local sam_zstd::Frames frames;
            packed.clear();
            frames.clear();
            bool const packed_ok = size == 0 || Pack(data, size, packed, frames);
            std::sort(genes.begin(), genes.end());
            genes.erase(std::unique(genes.begin(), genes.end()), genes.end());

            std::lock_guard<std::mutex> lock(m_mutex);
            m_genes.insert(genes.begin(), genes.end());
            genes.clear();
            if (!m_error.empty() || size == 0) return;
            if (!packed_ok) {
                Fail("compressing the records failed");
                return;
            }
            bool const plain = m_compression == SamCompression::None;
            if (!WriteAll(plain ? data : packed.data(), plain ? size : packed.size())) return;
            m_record_frames.insert(m_record_frames.end(), frames.begin(), frames.end());
        }

        // The genes that the records name, sorted.
        std::vector<uint64_t> Genes() const {
            std::lock_guard<std::mutex> lock(m_mutex);
            std::vector<uint64_t> genes(m_genes.begin(), m_genes.end());
            std::sort(genes.begin(), genes.end());
            return genes;
        }

        // How Finish placed the header, for the run's log: its size as written (compressed, with zstd's marker), the room
        // the records were written behind (zstd; 0: none, they collected in the records file), and whether the records
        // were copied behind the header (no room, or a header too large for it), how many bytes, and in how long.
        struct HeadPlacement {
            uint64_t header_bytes = 0, room = 0, copied_bytes = 0;
            bool copied = false;
            double copy_seconds = 0;
        };
        HeadPlacement const& Placement() const { return m_placement; }

        // Completes the SAM after the last Write: the header (unless one was given up front), the
        // records and the format's end (BGZF's end-of-file block, zstd's seek table). The temporary
        // records file is removed. Returns Ok().
        bool Finish(std::string const& header = {}) {
            if (m_finished) return Ok();
            m_finished = true;
            if (!Ok()) {
                Discard();
                return false;
            }
            if (!m_header_first && m_room > 0) {
                if (!PlaceHead(header)) return Abort();
            } else if (!m_header_first) {
                if (!CloseFd()) return Abort();
                m_fd = Create(m_path);
                if (m_fd < 0 || !WriteHead(header) || !AppendRecords()) return Abort();
                std::remove(m_records_path.c_str());
            }
            if (m_compression == SamCompression::Gzip) {
                WriteAll(reinterpret_cast<char const*>(bgzf::kEof), sizeof(bgzf::kEof));
            } else if (m_compression == SamCompression::Zstd) {
                sam_zstd::Frames all = m_head_frames;
                all.insert(all.end(), m_record_frames.begin(), m_record_frames.end());
                if (all.size() > 0xffffffffu) {
                    Fail("more than 2^32 zstd frames");
                } else {
                    std::string const table = sam_zstd::SeekTable(all);
                    WriteAll(table.data(), table.size());
                }
            }
            if (!CloseFd() || !Ok()) return Abort();
            return true;
        }

        // Removes what was written (the SAM and the temporary records file).
        void Discard() {
            m_finished = true;
            CloseFd();
            std::remove(m_records_path.c_str());
            std::remove(m_path.c_str());
        }

    private:
        std::string m_path;
        std::string m_records_path;
        SamCompression m_compression;
        bool m_header_first;
        uint64_t m_room = 0;               // zstd: where the records start in m_path (0: in m_records_path)
        bool m_finished = false;
        int m_fd = -1;
        mutable std::mutex m_mutex;
        std::string m_error;
        std::unordered_set<uint64_t> m_genes;
        sam_zstd::Frames m_head_frames;    // zstd: the marker and the header's frames
        sam_zstd::Frames m_record_frames;  // zstd: the records' frames, in file order
        HeadPlacement m_placement;

        void Fail(std::string const& what) {
            if (m_error.empty()) m_error = what;
        }

        bool Abort() {
            Discard();
            return false;
        }

        int Create(std::string const& path) {
            int const fd = ::open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC, 0666);
            if (fd < 0) Fail("cannot write " + path + ": " + std::strerror(errno));
            return fd;
        }

        bool CloseFd() {
            if (m_fd < 0) return true;
            int const fd = m_fd;
            m_fd = -1;
            if (::close(fd) != 0) {
                Fail("writing " + m_path + " failed: " + std::strerror(errno));
                return false;
            }
            return true;
        }

        // Appends to the open file; false (and Error() set) on failure.
        bool WriteAll(char const* data, size_t size) {
            while (size > 0) {
                ssize_t const n = ::write(m_fd, data, size);
                if (n < 0) {
                    if (errno == EINTR) continue;
                    Fail("writing " + m_path + " failed: " + std::strerror(errno));
                    return false;
                }
                data += n;
                size -= static_cast<size_t>(n);
            }
            return true;
        }

        // Writes at `offset` of the open file; false (and Error() set) on failure.
        bool PWriteAll(char const* data, size_t size, uint64_t offset) {
            while (size > 0) {
                ssize_t const n = ::pwrite(m_fd, data, size, static_cast<off_t>(offset));
                if (n < 0) {
                    if (errno == EINTR) continue;
                    Fail("writing " + m_path + " failed: " + std::strerror(errno));
                    return false;
                }
                data += n;
                size -= static_cast<size_t>(n);
                offset += static_cast<uint64_t>(n);
            }
            return true;
        }

        bool Pack(char const* data, size_t size, std::string& out, sam_zstd::Frames& frames) const {
            switch (m_compression) {
                case SamCompression::None: return true;
                case SamCompression::Gzip: return bgzf::Compress(data, size, out);
                case SamCompression::Zstd: return sam_zstd::Compress(data, size, out, frames);
            }
            return false;
        }

        // The start of the file: zstd's marker frame, then the header.
        bool WriteHead(std::string const& header) {
            std::string packed;
            sam_zstd::Frames frames;
            if (m_compression == SamCompression::Zstd) {
                packed = sam_zstd::MarkerFrame();
                frames.emplace_back(static_cast<uint32_t>(packed.size()), 0);
            }
            if (!header.empty() && !Pack(header.data(), header.size(), packed, frames)) {
                Fail("compressing the SAM header failed");
                return false;
            }
            bool const plain = m_compression == SamCompression::None;
            m_head_frames = std::move(frames);
            m_placement.header_bytes = plain ? header.size() : packed.size();
            return WriteAll(plain ? header.data() : packed.data(), plain ? header.size() : packed.size());
        }

        // With room for the header (zstd): the marker and the header's frames go to the start of
        // the SAM, and a skippable frame from there to the records, whose content is never written
        // (a hole where the file system has them). Cutting that out of the file
        // (FALLOC_FL_COLLAPSE_RANGE) would make the file system write the records to disk first:
        // 0.5-1.6 s for 400 MB just written. A header too large for the room is written as
        // without room: the records are copied behind it.
        bool PlaceHead(std::string const& header) {
            std::string head = sam_zstd::MarkerFrame();
            sam_zstd::Frames frames{ { static_cast<uint32_t>(head.size()), 0 } };
            if (!header.empty() && !sam_zstd::Compress(header.data(), header.size(), head, frames)) {
                Fail("compressing the SAM header failed");
                return false;
            }
            off_t const end = ::lseek(m_fd, 0, SEEK_CUR);  // where the records end
            if (end < 0) {
                Fail("writing " + m_path + " failed: " + std::strerror(errno));
                return false;
            }
            m_placement.header_bytes = head.size();
            m_placement.room = m_room;
            bool const records = static_cast<uint64_t>(end) > m_room;
            if (records && head.size() + sam_zstd::kSkippableHeader > m_room) return MoveRecordsBehind(header);
            // Without records, the header is all there is before the seek table.
            uint64_t const records_at = records ? m_room : head.size();
            if (records_at > head.size()) {
                uint64_t const padding = records_at - head.size();
                zstd::PutLE32(head, sam_zstd::kPaddingMagic);
                zstd::PutLE32(head, static_cast<uint32_t>(padding - sam_zstd::kSkippableHeader));
                frames.emplace_back(static_cast<uint32_t>(padding), 0);
            }
            m_head_frames = std::move(frames);
            if (!PWriteAll(head.data(), head.size(), 0)) return false;
            if (::lseek(m_fd, 0, SEEK_END) < 0) {
                Fail("writing " + m_path + " failed: " + std::strerror(errno));
                return false;
            }
            return true;
        }

        // The header did not fit its room: the records, from the room's end on, are copied behind
        // it, as when they collect in the temporary records file.
        bool MoveRecordsBehind(std::string const& header) {
            if (!CloseFd()) return false;
            std::error_code ec;
            std::filesystem::rename(m_path, m_records_path, ec);
            if (ec) {
                Fail("cannot move " + m_path + " to " + m_records_path + ": " + ec.message());
                return false;
            }
            m_fd = Create(m_path);
            if (m_fd < 0 || !WriteHead(header) || !AppendRecords(m_room)) return false;
            std::remove(m_records_path.c_str());
            return true;
        }

        // Copies the records file, from byte `from` on, to the end of the SAM (in the kernel where it can).
        bool AppendRecords(uint64_t from = 0) {
            int const in = ::open(m_records_path.c_str(), O_RDONLY | O_CLOEXEC);
            if (in < 0) {
                Fail("cannot read " + m_records_path + ": " + std::strerror(errno));
                return false;
            }
            struct stat st{};
            if (::fstat(in, &st) != 0 || ::lseek(in, static_cast<off_t>(from), SEEK_SET) < 0) {
                Fail("cannot read " + m_records_path + ": " + std::strerror(errno));
                ::close(in);
                return false;
            }
            if (static_cast<uint64_t>(st.st_size) < from) {
                Fail("reading " + m_records_path + " failed (it ended early)");
                ::close(in);
                return false;
            }
            auto left = static_cast<uint64_t>(st.st_size) - from;
            m_placement.copied = true;
            m_placement.copied_bytes = left;
            auto const copy_start = std::chrono::steady_clock::now();
            struct CopyTime {
                HeadPlacement& placement;
                std::chrono::steady_clock::time_point start;
                ~CopyTime() { placement.copy_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count(); }
            } copy_time{ m_placement, copy_start };
#ifdef __linux__
            while (left > 0) {
                ssize_t const n = ::copy_file_range(in, nullptr, m_fd, nullptr, left, 0);
                if (n > 0) {
                    left -= static_cast<uint64_t>(n);
                    continue;
                }
                if (n < 0 && errno == EINTR) continue;
                break;  // not supported here, or an error: the loop below copies the rest or reports it
            }
#endif
            std::vector<char> buffer(size_t{4} << 20);
            while (left > 0) {
                ssize_t const n = ::read(in, buffer.data(), std::min<uint64_t>(buffer.size(), left));
                if (n < 0 && errno == EINTR) continue;
                if (n <= 0) {
                    Fail("reading " + m_records_path + " failed" + (n < 0 ? std::string(": ") + std::strerror(errno) : std::string(" (it ended early)")));
                    ::close(in);
                    return false;
                }
                if (!WriteAll(buffer.data(), static_cast<size_t>(n))) {
                    ::close(in);
                    return false;
                }
                left -= static_cast<uint64_t>(n);
            }
            ::close(in);
            return true;
        }
    };

    // ---- Reading ---------------------------------------------------------------------------------

    // A SAM file for reading: plain, gzip or zstd, told apart by its content. Problem() says why a
    // file is incomplete before it is read (a missing end marker, see the top of this file);
    // ReadFailed() and ReadError() tell, once reading stopped, whether it ended early because the
    // data are truncated or corrupt. `decompress_threads`: 0 decompresses a zstd file in the thread that
    // reads it; more decompress a seekable one (protal's .sam.zst) frame by frame on that many threads of
    // their own, ahead of the reader (zstd::ParallelFrameStreambuf).
    class SamInput {
    public:
        explicit SamInput(std::string const& path, size_t decompress_threads = 0) {
            if (zstd::IsCompressed(path)) {
                std::optional<zstd::SeekTable> table;
                bool const marker = sam_zstd::StartsWithMarker(path);
                if (marker || decompress_threads > 0) {
                    std::string error;
                    table = zstd::ReadSeekTable(path, error);
                    if (marker && !table) m_problem = error.empty() ? "the seek table that ends it is missing" : error;
                }
                if (table && decompress_threads > 0) {
                    m_pbuf = std::make_unique<zstd::ParallelFrameStreambuf>(path, std::move(*table), decompress_threads,
                                                                             4 * decompress_threads + 4);
                    m_zin = std::make_unique<std::istream>(m_pbuf.get());
                    if (!m_pbuf->IsOpen()) m_zin->setstate(std::ios_base::badbit);
                } else {
                    m_zbuf = std::make_unique<zstd::IStreambuf>(path);
                    m_zin = std::make_unique<std::istream>(m_zbuf.get());
                    if (!m_zbuf->IsOpen()) m_zin->setstate(std::ios_base::badbit);
                }
                m_stream = m_zin.get();
            } else {
                // Plain or gzip (BGZF or not, inflated with ISA-L), in a thread of its own.
                if (bgzf::StartsAsBgzf(path) && !bgzf::EndsWithEof(path)) m_problem = "the BGZF end-of-file block is missing";
                m_gz = std::make_unique<ThreadedGzIstream>(path.c_str());
                m_stream = m_gz.get();
            }
        }

        SamInput(SamInput const&) = delete;
        SamInput& operator=(SamInput const&) = delete;

        bool IsOpen() const {
            if (m_zbuf) return m_zbuf->IsOpen();
            if (m_pbuf) return m_pbuf->IsOpen();
            return m_gz->rdbuf()->is_open();
        }
        std::istream& Stream() { return *m_stream; }
        std::string const& Problem() const { return m_problem; }

        bool ReadFailed() const {
            if (m_zbuf) return !m_zbuf->Error().empty() || m_zin->bad();
            if (m_pbuf) return !m_pbuf->Error().empty() || m_zin->bad();
            return m_gz->rdbuf()->read_failed();
        }

        std::string ReadError() const {
            std::string const& error = m_zbuf ? m_zbuf->Error() : m_pbuf ? m_pbuf->Error() : m_gz->rdbuf()->read_error_message();
            return error.empty() && (m_zbuf || m_pbuf) ? "read error" : error;
        }

    private:
        std::unique_ptr<ThreadedGzIstream> m_gz;
        std::unique_ptr<zstd::IStreambuf> m_zbuf;
        std::unique_ptr<zstd::ParallelFrameStreambuf> m_pbuf;
        std::unique_ptr<std::istream> m_zin;
        std::istream* m_stream = nullptr;
        std::string m_problem;
    };
}
