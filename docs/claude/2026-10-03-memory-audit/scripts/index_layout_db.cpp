// index_layout_db database.protal [threads]: layout_stats.h on the index inside a single-file database,
// decoded chunk by chunk (IndexCodec.h) without holding the index: memory is a few chunks (~64 MB raw
// each) per thread. For the r226 database on a cluster node: ~18 GB read, minutes. Also takes a
// seekable index.prx.zst in protal's column format.
// Build in the checkout's root (zstd headers and library needed; in a conda environment add -I $CONDA_PREFIX/include -L $CONDA_PREFIX/lib):
//   g++ -O2 -std=c++20 -I src -I src/Utilities -I src/Hash docs/claude/2026-10-03-memory-audit/scripts/index_layout_db.cpp -o index_layout_db -lzstd -pthread
#include "layout_stats.h"
#include "IndexCodec.h"
#include "Database.h"
#include <mutex>
#include <string>
using protal::zstd::SeekTable;
namespace zstd = protal::zstd;

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "usage: index_layout_db database.protal|index.prx.zst [threads]\n"); return 1; }
    std::string const path = argv[1];
    int const threads = argc > 2 ? std::atoi(argv[2]) : 4;
    std::string error;
    zstd::SeekTable frames;
    if (auto bundle = protal::db::Bundle::Open(path, error)) {
        auto const* member = bundle->Find("index.prx");
        if (!member) { fprintf(stderr, "the database has no index.prx member\n"); return 1; }
        frames = member->frames;
    } else if (!error.empty()) {
        fprintf(stderr, "%s: %s\n", path.c_str(), error.c_str()); return 1;
    } else {
        auto table = zstd::ReadSeekTable(path, error);
        if (!table) { fprintf(stderr, "%s: not a single-file database or seekable index (%s)\n", path.c_str(), error.c_str()); return 1; }
        frames = *table;
    }
    auto const container = protal::index_codec::ReadContainer(path, frames, error);
    if (!container) { fprintf(stderr, "no index in protal's column format: %s\n", error.c_str()); return 1; }
    auto const& l = container->layout;
    fprintf(stderr, "%zu chunks, %lu blocks of %lu keys, %lu value slots\n", container->chunks.size(), l.blocks, l.keys_per_block, l.values);
    size_t const workers = zstd::WorkerCount(container->chunks.size(), threads);
    std::vector<layout::Stats> stats(workers);
    std::vector<std::vector<uint16_t>> km(workers);
    std::vector<std::vector<uint64_t>> vals(workers);
    uint64_t const cpb = l.CellsPerBlock();
    uint64_t const kpb = l.keys_per_block;
    uint64_t const cc = l.ctrl_cells;
    error = zstd::ForEachFrame(path, frames, 1, threads, [&](size_t frame, char const* data, size_t size, size_t w) -> std::string {
        auto const& ch = container->chunks[frame - 1];
        km[w].resize(ch.blocks * cpb);
        vals[w].resize(ch.values);
        std::string const e = protal::index_codec::detail::DecodeChunk(data, size, l, ch, km[w].data(), vals[w].data());
        if (!e.empty()) return "chunk " + std::to_string(frame) + ": " + e;
        for (uint64_t b = 0; b < ch.blocks; b++) {
            uint16_t const* cells = km[w].data() + b * cpb;
            uint64_t start; std::memcpy(&start, cells, 8);
            uint64_t next = ch.first_value + ch.values;
            if (b + 1 < ch.blocks) std::memcpy(&next, cells + cpb, 8);
            if (next == start) continue;
            for (uint64_t j = 0; j < kpb; j++) {
                uint64_t const a = cells[cc + j], end = j + 1 < kpb ? cells[cc + j + 1] : next - start;
                stats[w].AddKey(vals[w].data() + (start - ch.first_value) + a, end - a);
            }
        }
        return "";
    });
    if (!error.empty()) { fprintf(stderr, "%s\n", error.c_str()); return 1; }
    layout::Stats all;
    for (auto const& s : stats) all.Merge(s);
    if (all.slots != l.values) fprintf(stderr, "warning: %lu slots counted, layout says %lu\n", all.slots, l.values);
    all.Print(l.blocks * kpb);
    return 0;
}
