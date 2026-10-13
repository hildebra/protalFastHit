// Reads a gzip file through ThreadedGzIstream in pieces of PIECE_BYTES on THREADS threads; prints how the pieces
// were decoded (ParallelGzip::Counts) and the bytes read.
// usage: pgz_counts FILE THREADS PIECE_BYTES
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "IO/ThreadedGzStream.h"

int main(int argc, char** argv) {
    if (argc < 4) return 2;
    unsigned const threads = static_cast<unsigned>(std::atoi(argv[2]));
    protal::ParallelGzip::SetPieceBytes(static_cast<size_t>(std::atol(argv[3])), argc > 4 ? static_cast<size_t>(std::atol(argv[4])) : 0);
    protal::ThreadedGzStreambuf::SetGzipThreads(threads);
    protal::ThreadedGzIstream is(argv[1]);
    std::vector<char> buffer(1 << 20);
    size_t total = 0;
    while (true) {
        std::streamsize const n = is.rdbuf()->sgetn(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        if (n <= 0) break;
        total += static_cast<size_t>(n);
    }
    auto const c = is.rdbuf()->parallel_gzip_counts();
    std::printf("%s\tthreads %u\tparallel %d\tbytes %zu\tahead %zu\tin_order %zu\tstreamed %zu\tdropped %zu\tunused %zu\tfailed %d %s\n",
                argv[1], threads, is.rdbuf()->parallel_gzip() ? 1 : 0, total, c.ahead, c.in_order, c.streamed,
                c.dropped_starts, c.unused, is.rdbuf()->read_failed() ? 1 : 0, is.rdbuf()->read_error_message().c_str());
    return 0;
}
