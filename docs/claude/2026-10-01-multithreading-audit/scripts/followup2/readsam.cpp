// readsam - experiment only: reads a SAM through SamInput (with `decompress` threads) and the profiler's chunk reader
// (SamChunks.h), doing nothing with the chunks: how fast the reader path hands text to the parsers.
#include "IO/SamFile.h"
#include "Profiling/SamChunks.h"
#include <chrono>
#include <cstdio>
#include <cstdlib>

int main(int argc, char** argv) {
    if (argc < 3) { std::fprintf(stderr, "usage: readsam SAM DECOMPRESS_THREADS\n"); return 2; }
    size_t const threads = std::strtoul(argv[2], nullptr, 10);
    auto const t0 = std::chrono::steady_clock::now();
    protal::SamInput input(argv[1], threads);
    protal::sam_chunks::ChunkReader reader(input.Stream(), size_t{1} << 20, 48);
    protal::sam_chunks::Chunk chunk;
    size_t bytes = 0, chunks = 0;
    while (reader.Next(chunk)) { bytes += chunk.text.size(); chunks++; }
    reader.Stop();
    double const s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("decompress_threads\t%zu\tbytes\t%zu\tchunks\t%zu\tseconds\t%.3f\tMB_per_s\t%.0f\tfailed\t%d\n",
                threads, bytes, chunks, s, bytes / s / 1e6, static_cast<int>(input.ReadFailed()));
    return 0;
}
