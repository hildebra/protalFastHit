// codec_memory.cpp - memory of the index write (index_codec::Write) and its read-back (index_codec::Verify), the
// last phase of protal --build, on an index of GTDB r226's composition but a smaller key space.
//
//   codec_memory OUT_FILE THREADS [KEY_BITS=26] [FRAME_MB=64]
//
// Compiled against a source tree's headers (build_harnesses.sh), so the same file measures HEAD and a change.
//
// The index: 2^KEY_BITS keys in control blocks of 8, 9.2% of the keys filled (r226 v7, memory audit section 9),
// half of them with 1-8 values, half with a geometric number of mean 50 (about 29 per filled key; r226: 29),
// at most 2047; each key's cells as Seedmap lays them out (n = e + ceil(e/2) cells, the first (n + 2) / 3 flex
// cells, then the entries: taxid 18 bits, gene id < 168, position 14 bits, two flag bits). A chunk thus holds
// about as many value cells as an r226 chunk (~6M). Prints the resident memory above the index while each
// stage runs, its time and the file size.
#include "IndexCodec.h"
#include "rss.h"

#include <chrono>
#include <cstdio>
#include <random>
#include <vector>

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: %s OUT_FILE THREADS [KEY_BITS] [FRAME_MB]\n", argv[0]);
        return 2;
    }
    std::string const path = argv[1];
    int const threads = std::atoi(argv[2]);
    int const key_bits = argc > 3 ? std::atoi(argv[3]) : 26;
    int const frame_mb = argc > 4 ? std::atoi(argv[4]) : 64;
    using namespace protal;
    using namespace protal::index_codec;

    Layout l;
    l.blocks = uint64_t{1} << (key_bits - 3);
    std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> unit(0, 1);
    std::geometric_distribution<int> many(1.0 / 50);
    std::vector<uint16_t> counts(l.blocks * 8, 0);  // cells per key
    uint64_t total = 0;
    for (uint64_t k = 0; k < counts.size(); k++) {
        if (unit(rng) >= 0.092) continue;
        int const e = std::min(2047, unit(rng) < 0.5 ? 1 + static_cast<int>(rng() % 8) : 1 + many(rng));
        uint64_t const n = e == 1 ? 1 : e + (e + 1) / 2;
        counts[k] = static_cast<uint16_t>(n);
        total += n;
    }
    l.values = total;
    std::vector<uint16_t> keymap(l.KeymapCells(), 0);
    std::vector<uint64_t> values(l.values);
    uint64_t v = 0;
    for (uint64_t b = 0; b < l.blocks; b++) {
        uint16_t* cells = keymap.data() + b * l.CellsPerBlock();
        std::memcpy(cells, &v, 8);
        uint64_t local = 0;
        for (uint64_t j = 0; j < 8; j++) {
            cells[l.ctrl_cells + j] = static_cast<uint16_t>(local);
            uint64_t const n = counts[b * 8 + j];
            uint64_t const f = detail::FlexCells(n);
            for (uint64_t i = 0; i < n; i++) {
                values[v + local + i] = i < f ? rng()
                    : (rng() & 3) << 60 | (rng() % (1u << 18)) << 40 | (rng() % 168) << 20 | (rng() % (1u << 14));
            }
            local += n;
        }
        v += local;
    }
    std::memcpy(keymap.data() + l.blocks * l.CellsPerBlock(), &v, 8);
    counts = {};
    std::string const header = "harness header";
    zstd::Params params{3, 27, threads, uint64_t(frame_mb) << 20};

    double const base = RssMb();
    std::printf("index: %llu blocks, %llu value cells (%.2f GB values, %.2f GB key map); resident %.0f MB\n",
                (unsigned long long)l.blocks, (unsigned long long)l.values, l.values * 8 / 1e9, keymap.size() * 2 / 1e9, base);
    std::string error;
    size_t raw = 0;
    auto t0 = std::chrono::steady_clock::now();
    PeakSampler write_peak;
    auto const written = Write(path, header, l, keymap.data(), values.data(), params, error, &raw);
    double const w = write_peak.Stop();
    auto t1 = std::chrono::steady_clock::now();
    if (!written) {
        std::fprintf(stderr, "write: %s\n", error.c_str());
        return 1;
    }
    double const after_write = RssMb();
    PeakSampler verify_peak;
    error = Verify(path, header, l, keymap.data(), values.data(), threads);
    double const r = verify_peak.Stop();
    auto t2 = std::chrono::steady_clock::now();
    if (!error.empty()) {
        std::fprintf(stderr, "verify: %s\n", error.c_str());
        return 1;
    }
    std::string table_error;
    auto const table = zstd::ReadSeekTable(path, table_error);
    size_t const chunks = table ? table->frames.size() - 1 : 0;
    auto secs = [](auto a, auto b) { return std::chrono::duration<double>(b - a).count(); };
    std::printf("threads %d, frame %d MB: %zu chunks (%zu raw), %.1f MB per chunk of values; file %.2f GB\n", threads, frame_mb,
                chunks, raw, l.values * 8.0 / 1048576 / std::max<size_t>(chunks, 1), *written / 1e9);
    std::printf("write  %6.1f s  peak +%6.0f MB above the index (%.1f MB per worker); after it +%.0f MB\n", secs(t0, t1), w - base,
                (w - base) / std::min<size_t>(threads, chunks), after_write - base);
    std::printf("verify %6.1f s  peak +%6.0f MB above the index (%.1f MB per worker); VmHWM %.0f MB\n", secs(t1, t2), r - base,
                (r - base) / std::min<size_t>(threads, chunks), HwmMb());
    // Decoding as a query run's load does (Decode into fresh arrays; LoadColumns packs from the same DecodeChunk):
    // best of three, and a check that it gives the index back.
    {
        auto const container = ReadContainer(path, *table, table_error);
        std::vector<uint16_t> km(l.KeymapCells());
        std::vector<uint64_t> vals(l.values);
        double best = 1e9;
        for (int round = 0; round < 3; round++) {
            auto d0 = std::chrono::steady_clock::now();
            error = Decode(path, *table, *container, km.data(), vals.data(), threads);
            best = std::min(best, secs(d0, std::chrono::steady_clock::now()));
            if (!error.empty()) {
                std::fprintf(stderr, "decode: %s\n", error.c_str());
                return 1;
            }
        }
        bool const same = km == keymap && vals == values;
        std::printf("decode %6.2f s  (best of 3), %s\n", best, same ? "the index back" : "NOT THE INDEX");
    }
    std::remove(path.c_str());
    return 0;
}
