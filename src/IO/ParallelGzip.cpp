// ParallelGzip.cpp - ParallelGzip on rapidgzip's ParallelGzipReader (lib/rapidgzip.cmake).

#include "ParallelGzip.h"

#include <sys/stat.h>

#include <cstdint>
#include <cstring>
#include <exception>
#include <utility>

#include <filereader/Standard.hpp>
#include <rapidgzip/ParallelGzipReader.hpp>

namespace protal {

    // rapidgzip ends a file it cannot decode to its end as if it ended there: a gzip file cut inside a member reads as
    // one that is shorter, and nothing is said. It checks a member's CRC-32 only while the pieces it hands over follow
    // each other in the file, and its length (ISIZE) only where it read the member's header in the same piece. So the
    // pieces are checked here as they are read: each must begin where the one before ended in the file, each member's
    // bytes must number its ISIZE (modulo 2^32), and at the end the bytes must end with a member, whose trailer ends
    // the file. A cut file, data after the last member and a piece left out all fail.
    struct ParallelGzip::Impl {
        using Reader = rapidgzip::ParallelGzipReader<>;
        using Chunk = rapidgzip::ChunkData;

        Reader reader;
        uint64_t file_bits = 0;               // the file's size in bits
        uint64_t total = 0;                   // bytes read so far
        uint64_t member_start = 0;            // where the current member's bytes begin
        uint64_t last_end = 0;                // where the last whole member's bytes end
        uint64_t last_trailer_bits = 0;       // and its trailer in the file (in bits)
        size_t members = 0;                   // whole members read
        Chunk const* chunk = nullptr;         // the piece being read
        size_t next_footer = 0;               // its next member end (footer) not yet passed
        uint64_t chunk_end_bits = UINT64_MAX; // where the last piece ended in the file (in bits); none yet
        std::string error;

        Impl(int fd, unsigned threads, uint64_t file_bytes) :
            reader(std::make_unique<rapidgzip::StandardFileReader>(fd), threads, kChunkBytes), file_bits(8 * file_bytes) {
            reader.setCRC32Enabled(true);
            // The file is read once, front to back: the windows of the pieces read need not be kept for seeking.
            reader.setKeepIndex(false);
        }

        // Copies `size` bytes of `piece` from `offset` to out (the piece's bytes from there are the file's next ones),
        // checking the piece's place and the members ending in it.
        void Take(std::shared_ptr<Chunk> const& piece, size_t offset, size_t size, char* out) {
            if (piece.get() != chunk || offset == 0) {
                if (offset == 0 && chunk_end_bits != UINT64_MAX && piece->encodedOffsetInBits != chunk_end_bits &&
                    error.empty()) {
                    error = "a piece of the file was left out (bits " + std::to_string(chunk_end_bits) + " to " +
                            std::to_string(piece->encodedOffsetInBits) + ")";
                }
                chunk = piece.get();
                next_footer = 0;
                chunk_end_bits = piece->encodedOffsetInBits + piece->encodedSizeInBits;
            }
            uint64_t const piece_start = total - offset;  // where the piece's bytes begin among the file's
            auto const& footers = piece->footers;
            while (next_footer < footers.size() && footers[next_footer].blockBoundary.decodedOffset <= offset + size) {
                auto const& footer = footers[next_footer++];
                uint64_t const end = piece_start + footer.blockBoundary.decodedOffset;
                uint32_t const length = static_cast<uint32_t>(end - member_start);
                if (length != footer.gzipFooter.uncompressedSize && error.empty()) {
                    error = "gzip member " + std::to_string(members + 1) + " holds " + std::to_string(end - member_start) +
                            " bytes, its trailer says " + std::to_string(footer.gzipFooter.uncompressedSize) +
                            " (modulo 2^32)";
                }
                member_start = end;
                last_end = end;
                last_trailer_bits = footer.blockBoundary.encodedOffset;
                members++;
            }
            size_t copied = 0;
            for (auto it = rapidgzip::deflate::DecodedData::Iterator(*piece, offset, size); static_cast<bool>(it); ++it) {
                auto const& [buffer, bytes] = *it;
                std::memcpy(out + copied, buffer, bytes);
                copied += bytes;
            }
            total += size;
        }

        // At the end of the file: the bytes end with a whole member, whose trailer ends the file.
        void CheckTheEnd() {
            if (!error.empty()) return;
            if (members == 0) {
                error = "the file ends inside gzip member 1 (truncated file?)";
            } else if (last_end != total) {
                error = "the file ends inside gzip member " + std::to_string(members + 1) + " (truncated file?)";
            } else if (last_trailer_bits != file_bits) {
                error = "the data at byte " + std::to_string(last_trailer_bits / 8) + ", after gzip member " +
                        std::to_string(members) + ", is no gzip member that rapidgzip reads (a damaged member header?)";
            }
        }
    };

    ParallelGzip::ParallelGzip() = default;
    ParallelGzip::~ParallelGzip() = default;

    std::unique_ptr<ParallelGzip> ParallelGzip::Open(int fd, unsigned threads, std::string& error) {
        struct stat st {};
        if (::fstat(fd, &st) != 0) {
            error = std::strerror(errno);
            return nullptr;
        }
        std::unique_ptr<ParallelGzip> gz(new ParallelGzip());
        try {
            gz->m_impl = std::make_unique<Impl>(fd, threads, static_cast<uint64_t>(st.st_size));
        } catch (std::exception const& e) {
            error = e.what();
            return nullptr;
        }
        return gz;
    }

    size_t ParallelGzip::Read(char* out, size_t size) {
        if (!m_error.empty() || !m_impl) return 0;
        Impl& impl = *m_impl;
        size_t copied = 0;
        try {
            impl.reader.read([&impl, out, &copied](std::shared_ptr<Impl::Chunk> const& piece, size_t offset, size_t bytes) {
                impl.Take(piece, offset, bytes, out + copied);
                copied += bytes;
            }, size);
            if (copied < size || impl.reader.eof()) impl.CheckTheEnd();
        } catch (std::exception const& e) {
            m_error = e.what();
            if (m_error.empty()) m_error = "rapidgzip stopped with an error";
            return 0;
        }
        if (!impl.error.empty()) {
            m_error = impl.error;
            return 0;
        }
        return copied;
    }

}  // namespace protal
