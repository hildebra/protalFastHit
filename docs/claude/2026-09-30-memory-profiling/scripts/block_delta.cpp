// How compressible are the index's value blocks in memory? block_delta index.prx
// For every key with several values (a flex block, see Seedmap.h) the values are coded one after
// the other as the previous value's differences: taxid (zigzag varint), gene (one byte, 0 = same),
// position (zigzag varint), and the two flags in the gene byte; a block of n values costs a length
// byte plus that. Compared with the 8 bytes per value the index uses now. Keys with one value stay 8 bytes.
// Stored order and, for comparison, the block sorted by (taxid, gene, position).
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
static inline int varint(uint64_t v) { int n = 1; while (v >= 128) { v >>= 7; n++; } return n; }
static inline uint64_t zz(int64_t d) { return (uint64_t)((d << 1) ^ (d >> 63)); }
int main(int, char** argv) {
    int fd = open(argv[1], O_RDONLY);
    struct stat st; fstat(fd, &st);
    auto* base = (const uint8_t*)mmap(nullptr, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    uint32_t feat = 0; size_t off = 8 + 12; memcpy(&feat, base + off, 4); off += 4;
    if (feat & 4) off += 16;
    uint64_t keymap_size, total, cbs, kpb, shift, values;
    memcpy(&keymap_size, base + off, 8); memcpy(&total, base + off + 8, 8); memcpy(&cbs, base + off + 16, 8);
    memcpy(&kpb, base + off + 24, 8); memcpy(&shift, base + off + 32, 8); memcpy(&values, base + off + 40, 8);
    off += 48;
    auto* km = (const uint16_t*)(base + off);
    auto* slots = (const uint64_t*)(base + off + total * 2);
    const size_t K = 8, C = 4;
    uint64_t single = 0, multi_vals = 0, multi_flex_bytes = 0, bytes_stored = 0, bytes_sorted = 0, blocks = 0;
    std::vector<uint64_t> tmp;
    for (uint64_t g = 0; g < keymap_size / K; g++) {
        uint64_t b = g * (K + C), s0, s1; memcpy(&s0, km + b, 8); memcpy(&s1, km + b + K + C, 8);
        if (s1 == s0) continue;
        for (size_t j = 0; j < K; j++) {
            uint64_t a = km[b + C + j], e = (j + 1 == K) ? s1 - s0 : km[b + C + j + 1];
            uint64_t S = e - a; if (!S) continue;
            if (S == 1) { single++; continue; }
            uint64_t f = (S + 2) / 3, n = S - f;
            const uint64_t* v = slots + s0 + a + f;
            tmp.assign(v, v + n);
            multi_vals += n; multi_flex_bytes += 4 * n; blocks++;
            for (int pass = 0; pass < 2; pass++) {
                if (pass == 1) std::sort(tmp.begin(), tmp.end(), [](uint64_t x, uint64_t y) { return (x & ((1ull << 60) - 1)) < (y & ((1ull << 60) - 1)); });
                uint64_t cost = 1 + 8, pt = tmp[0] >> 40 & 0xfffff, pg = tmp[0] >> 20 & 0xfffff, pp = tmp[0] & 0xfffff;  // first value whole
                for (size_t i = 1; i < n; i++) {
                    uint64_t t = tmp[i] >> 40 & 0xfffff, ge = tmp[i] >> 20 & 0xfffff, p = tmp[i] & 0xfffff;
                    cost += varint(zz((int64_t)t - (int64_t)pt)) + 1 + (ge != pg ? varint(ge) : 0) + varint(zz((int64_t)p - (int64_t)pp));
                    pt = t; pg = ge; pp = p;
                }
                (pass ? bytes_sorted : bytes_stored) += cost;
            }
        }
    }
    double raw = multi_vals * 8.0;
    printf("single-value keys %lu (8 bytes each); keys with several values %lu holding %lu values\n", single, blocks, multi_vals);
    printf("those values: now %.1f MB (8 B each) + flex %.1f MB\n", raw / 1e6, multi_flex_bytes / 1e6);
    printf("delta coded, stored order: %.1f MB (%.2f B/value, %.1f%% of now); sorted by taxid: %.1f MB (%.2f B/value, %.1f%%)\n",
           bytes_stored / 1e6, bytes_stored / (double)multi_vals, 100.0 * bytes_stored / raw, bytes_sorted / 1e6, bytes_sorted / (double)multi_vals, 100.0 * bytes_sorted / raw);
}
