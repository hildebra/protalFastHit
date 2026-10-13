// Speed and accuracy of rapidgzip's dynamic-Huffman block finder on one gzip member.
//   bfinder <in.gz> <bounds.tsv> <pieceMiB>
// 1. full scan: every candidate in the whole deflate stream (time, MB/s, true/false positives vs bounds.tsv)
// 2. piece starts: from every pieceMiB of compressed data, the first candidate that is a true block start
#include <rapidgzip/blockfinder/DynamicHuffman.hpp>
#include <filereader/BufferView.hpp>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <memory>
#include <unordered_set>
#include <vector>

static double now()
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

int main(int argc, char** argv)
{
    if (argc < 4) return 1;
    int fd = open(argv[1], O_RDONLY);
    struct stat sb{};
    fstat(fd, &sb);
    const size_t n = size_t(sb.st_size);
    const void* m = mmap(nullptr, n, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, 0);
    close(fd);

    std::vector<uint64_t> bits;
    FILE* f = fopen(argv[2], "r");
    char line[512];
    if (!fgets(line, sizeof(line), f)) return 1;
    unsigned long long x, y;
    while (fscanf(f, "%llu\t%llu\n", &x, &y) == 2) bits.push_back(x);
    fclose(f);
    const std::unordered_set<uint64_t> truth(bits.begin(), bits.end());
    const uint64_t startBit = bits.front(), endBit = uint64_t(n - 8) * 8;

    // 1. Full scan.
    {
        rapidgzip::gzip::BitReader br(std::make_unique<rapidgzip::BufferViewFileReader>(m, n));
        size_t found = 0, truePos = 0, falsePos = 0;
        const double t0 = now();
        uint64_t pos = startBit + 1;
        for (;;) {
            br.seekTo(pos);
            const size_t c = rapidgzip::blockfinder::seekToNonFinalDynamicDeflateBlock(br, endBit);
            if (c == std::numeric_limits<size_t>::max() || c >= endBit) break;
            ++found;
            if (truth.count(c)) ++truePos; else ++falsePos;
            pos = c + 1;
        }
        const double t1 = now();
        printf("finder_full_scan\tcompressed_MB\t%.1f\ts\t%.3f\tMBps\t%.1f\tcandidates\t%zu\ttrue_block_starts\t%zu\tof\t%zu\tfalse_positives\t%zu\n",
               (endBit - startBit) / 8e6, t1 - t0, (endBit - startBit) / 8e6 / (t1 - t0), found, truePos, bits.size() - 1, falsePos);
    }
    // 2. Piece starts as a real implementation would find them.
    {
        rapidgzip::gzip::BitReader br(std::make_unique<rapidgzip::BufferViewFileReader>(m, n));
        const uint64_t step = uint64_t(atof(argv[3]) * 1048576.0) * 8;
        size_t pieces = 0, falseSkipped = 0;
        uint64_t scannedBits = 0;
        const double t0 = now();
        for (uint64_t target = startBit + step; target < endBit; target += step) {
            uint64_t pos = target;
            for (;;) {
                br.seekTo(pos);
                const size_t c = rapidgzip::blockfinder::seekToNonFinalDynamicDeflateBlock(br, endBit);
                if (c == std::numeric_limits<size_t>::max()) { pos = endBit; break; }
                if (truth.count(c)) { pos = c; break; }  // a real implementation would verify by decoding
                ++falseSkipped;
                pos = c + 1;
            }
            if (pos >= endBit) break;
            scannedBits += pos - target;
            ++pieces;
        }
        const double t1 = now();
        printf("finder_piece_starts\tpiece_MiB\t%s\tpieces\t%zu\ts\t%.4f\tms_per_piece\t%.3f\tmean_scan_KiB\t%.1f\tfalse_candidates_skipped\t%zu\n",
               argv[3], pieces, t1 - t0, pieces ? 1000 * (t1 - t0) / pieces : 0.0, pieces ? scannedBits / 8.0 / 1024 / pieces : 0.0,
               falseSkipped);
    }
    return 0;
}
