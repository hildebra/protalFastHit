// Where the index's memory goes: index_stats index.prx   (a raw index, e.g. from protal --unpack_db)
// Reads the layout documented in Seedmap.h: header, key map (16-bit cells: per 8 keys one 64-bit
// control block holding the first value slot of the group, then 8 offsets), then 8-byte value slots.
// A key with n >= 2 slots... see Seedmap::Get: block of S slots = S - ceil(S/3) values + ceil(S/3) flex slots.
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#include <algorithm>
int main(int argc, char** argv) {
    int fd = open(argv[1], O_RDONLY);
    struct stat st; fstat(fd, &st);
    auto* base = (const uint8_t*)mmap(nullptr, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    uint64_t magic; uint32_t feat = 0; memcpy(&magic, base, 8);
    size_t off = 8 + 12;  // magic, 3 version words
    memcpy(&feat, base + off, 4); off += 4;
    if (feat & 4) off += 16;  // reference fingerprint (hash, size)
    uint64_t keymap_size, total, cbs, kpb, shift, values;
    memcpy(&keymap_size, base + off, 8); memcpy(&total, base + off + 8, 8); memcpy(&cbs, base + off + 16, 8);
    memcpy(&kpb, base + off + 24, 8); memcpy(&shift, base + off + 32, 8); memcpy(&values, base + off + 40, 8);
    off += 48;
    auto* km = (const uint16_t*)(base + off);
    fprintf(stderr, "keys %lu cells %lu slots %lu (%.2f GB keymap, %.2f GB values)\n", keymap_size, total, values, total * 2 / 1e9, values * 8 / 1e9);
    const size_t K = 8, C = 4;  // keys per control block, cells per control block
    // bins by slots per key: 0,1,2,3-4,5-8,... ; count keys, slots
    const int NB = 20; uint64_t keys[NB] = {}, slots[NB] = {}; uint64_t nonempty = 0, flexslots = 0, vals = 0;
    for (uint64_t g = 0; g < keymap_size / K; g++) {
        uint64_t b = g * (K + C);
        uint64_t s0, s1; memcpy(&s0, km + b, 8); memcpy(&s1, km + b + K + C, 8);
        if (s1 == s0) continue;
        for (size_t j = 0; j < K; j++) {
            uint64_t a = km[b + C + j], e = (j + 1 == K) ? s1 - s0 : km[b + C + j + 1];
            uint64_t S = e - a; if (!S) continue;
            int bin = S <= 2 ? (int)S : 1 + 64 - __builtin_clzll(S - 1);
            keys[bin]++; slots[bin] += S; nonempty++;
            uint64_t f = S >= 2 ? (S + 2) / 3 : 0; flexslots += f; vals += S - f;
        }
    }
    printf("non-empty keys %lu of %lu (%.1f%%); values %lu, flex slots %lu (%.1f%% of slots)\n", nonempty, keymap_size, 100.0 * nonempty / keymap_size, vals, flexslots, 100.0 * flexslots / values);
    printf("%-14s %14s %10s %14s %10s\n", "slots/key", "keys", "% keys", "slots", "% slots");
    for (int i = 1; i < NB; i++) if (keys[i]) {
        char lab[32]; if (i <= 2) snprintf(lab, 32, "%d", i); else snprintf(lab, 32, "%lu-%lu", (1ul << (i - 2)) + 1, 1ul << (i - 1));
        printf("%-14s %14lu %9.2f%% %14lu %9.2f%%\n", lab, keys[i], 100.0 * keys[i] / nonempty, slots[i], 100.0 * slots[i] / values);
    }
}
