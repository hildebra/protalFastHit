// index_layout index.prx: layout_stats.h on a raw index (format 1 "PRXSEEDM" or 2 "PRXSEED2"; protal
// --decompress_db gives one). g++ -O2 -std=c++17 -o index_layout index_layout.cpp
#include "layout_stats.h"
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "usage: index_layout index.prx\n"); return 1; }
    int fd = open(argv[1], O_RDONLY);
    if (fd < 0) { perror("open"); return 1; }
    struct stat st; fstat(fd, &st);
    auto* base = (const uint8_t*)mmap(nullptr, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    uint64_t magic; memcpy(&magic, base, 8);
    size_t off = 8 + 12;  // magic, 3 version words
    if (magic == 0x5052585345454432ull) {  // PRXSEED2: features, maybe a fingerprint
        uint32_t feat = 0; memcpy(&feat, base + off, 4); off += 4;
        if (feat & 4) off += 16;
    } else if (magic != 0x505258534545444dull) { fprintf(stderr, "not a raw protal index\n"); return 1; }
    uint64_t keymap_size, total, cbs, kpb, shift, values;
    memcpy(&keymap_size, base + off, 8); memcpy(&total, base + off + 8, 8); memcpy(&cbs, base + off + 16, 8);
    memcpy(&kpb, base + off + 24, 8); memcpy(&shift, base + off + 32, 8); memcpy(&values, base + off + 40, 8);
    off += 48;
    if (kpb != 8 || cbs != 8) { fprintf(stderr, "unexpected layout\n"); return 1; }
    auto* km = (const uint16_t*)(base + off);
    auto* slots = (const uint64_t*)(base + off + total * 2);
    const size_t K = 8, C = 4;
    layout::Stats stats;
    for (uint64_t g = 0; g < keymap_size / K; g++) {
        uint64_t b = g * (K + C), s0, s1; memcpy(&s0, km + b, 8); memcpy(&s1, km + b + K + C, 8);
        if (s1 == s0) continue;
        for (size_t j = 0; j < K; j++) {
            uint64_t a = km[b + C + j], end = (j + 1 == K) ? s1 - s0 : km[b + C + j + 1];
            stats.AddKey(slots + s0 + a, end - a);
        }
    }
    if (stats.slots != values) fprintf(stderr, "warning: %lu slots counted, header says %lu\n", stats.slots, values);
    stats.Print(keymap_size);
    return 0;
}
