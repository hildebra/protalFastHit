// ParallelGzip.h - one gzip file inflated on several threads (rapidgzip, lib/rapidgzip), for
// ThreadedGzStreambuf.
//
// A gzip member (as gzip and most tools write .fastq.gz files: the whole file one member) cannot be cut
// at block boundaries as BGZF can: its deflate blocks are found only by decoding it from the start, and
// refer back up to 32 KB. rapidgzip decodes pieces of the file (1 MiB of it each) on its threads at once,
// each from a deflate block it finds by trial near the piece's start, with the 32 KB before it unknown at
// first (its bytes marked, and filled in once the piece before is done); once 32 KB of a piece are known,
// ISA-L inflates the rest of it. Read gives the file's bytes in order, its members one after the other;
// each member's CRC-32 and length are checked, and that the file ends with a whole member (rapidgzip
// alone reads a cut file as a shorter one: ParallelGzip.cpp). ParallelGzip.cpp is the only file that
// includes rapidgzip, whose headers are large.
#pragma once

#include <cstddef>
#include <memory>
#include <string>

namespace protal {

    class ParallelGzip {
    public:
        // Bytes of the compressed file each piece decoded on a thread holds. rapidgzip's default of 4 MiB held ~780 MB
        // on 4 threads for a FASTQ file, 1 MiB ~200 MB, ~5% slower (docs/claude/2026-10-12-rapidgzip-input).
        static constexpr size_t kChunkBytes = size_t{1} << 20;

        // A reader of the regular file open as fd, from its start, on `threads` threads (rapidgzip's own):
        // fd is duplicated, so its offset and the descriptor stay the caller's. Null if the reader cannot
        // be made (error then says why).
        static std::unique_ptr<ParallelGzip> Open(int fd, unsigned threads, std::string& error);

        ~ParallelGzip();
        ParallelGzip(ParallelGzip const&) = delete;
        ParallelGzip& operator=(ParallelGzip const&) = delete;

        // Up to size inflated bytes into out, the next ones of the file; 0 at its end, or on an error
        // (Error() then says what: a corrupt member, a CRC-32 that does not match, a truncated file).
        size_t Read(char* out, size_t size);

        std::string const& Error() const { return m_error; }

    private:
        ParallelGzip();

        struct Impl;
        std::unique_ptr<Impl> m_impl;
        std::string m_error;
    };

}  // namespace protal
