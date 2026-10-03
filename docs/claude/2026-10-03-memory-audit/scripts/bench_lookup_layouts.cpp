// What a k-mer lookup's block work costs under other in-memory layouts of the index values:
//   bench_lookup_layouts index.prx [lookups]
// The index's values are copied into four layouts, and for a sample of keys (drawn in proportion to
// their entries, i.e. as reads from database species hit them) each layout does what
// KmerLookupSM::GetFromLookup does: compare every flex cell of the key with the read's (Similarity),
// take the entries at the best similarity (none if more than 256), and read their fields.
//   now         8-byte slots: flex cells first, then 8-byte entries (random access by index)
//   6-byte      entries 6 bytes each in one array, flex cells 4 bytes each in another (random access)
//   bit-packed  per key: minima and bit widths of taxid, gene, position; every entry in w bits; flex
//               cells 32 bits (random access by bit offset)
//   delta       per key: entries as varint deltas (taxid zigzag, gene if changed, position zigzag,
//               flags), flex cells as XOR varints: the whole key is decoded sequentially, then read
// The read's flex cell is the key's first flex cell with one base changed, so that similar but not
// identical cells are the common case, as for a read with a mismatch in the flanks.
// g++ -O2 -std=c++17 -mpopcnt -o bench_lookup_layouts bench_lookup_layouts.cpp
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <numeric>
#include <random>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

static inline uint32_t Similarity(uint32_t a, uint32_t b) {
    return __builtin_popcount(((~(a ^ b) >> 1) & ~(a ^ b)) & 0b01010101010101010101010101010101);
}
static inline int bits(uint64_t v) { return v == 0 ? 0 : 64 - __builtin_clzll(v); }
static inline uint64_t zz(int64_t d) { return ((uint64_t)d << 1) ^ (uint64_t)(d >> 63); }
static inline int64_t unzz(uint64_t v) { return (int64_t)(v >> 1) ^ -(int64_t)(v & 1); }
static inline void put_varint(std::vector<uint8_t>& out, uint64_t v) { while (v >= 128) { out.push_back(uint8_t(v | 128)); v >>= 7; } out.push_back(uint8_t(v)); }
static inline uint64_t get_varint(const uint8_t*& p) { uint64_t v = 0; int s = 0; while (*p & 128) { v |= uint64_t(*p & 127) << s; s += 7; p++; } v |= uint64_t(*p) << s; p++; return v; }

struct Key { uint64_t slot; uint32_t e; uint64_t off6, offf, bitoff, deltaoff, dictoff, bitoff2; uint32_t mn_t, mn_g, mn_p, nd; uint8_t w_t, w_g, w_p, w_d; };

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "usage: bench_lookup_layouts index.prx [lookups]\n"); return 1; }
    size_t const lookups = argc > 2 ? std::strtoull(argv[2], nullptr, 10) : 4000000;
    int fd = open(argv[1], O_RDONLY);
    if (fd < 0) { perror("open"); return 1; }
    struct stat st; fstat(fd, &st);
    auto* base = (const uint8_t*)mmap(nullptr, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    uint64_t magic; memcpy(&magic, base, 8);
    size_t off = 8 + 12;
    if (magic == 0x5052585345454432ull) { uint32_t feat = 0; memcpy(&feat, base + off, 4); off += 4; if (feat & 4) off += 16; }
    else if (magic != 0x505258534545444dull) { fprintf(stderr, "not a raw protal index\n"); return 1; }
    uint64_t keymap_size, total, cbs, kpb, shift, values;
    memcpy(&keymap_size, base + off, 8); memcpy(&total, base + off + 8, 8); memcpy(&cbs, base + off + 16, 8);
    memcpy(&kpb, base + off + 24, 8); memcpy(&shift, base + off + 32, 8); memcpy(&values, base + off + 40, 8);
    off += 48;
    auto* km = (const uint16_t*)(base + off);
    std::vector<uint64_t> slots((const uint64_t*)(base + off + total * 2), (const uint64_t*)(base + off + total * 2) + values);

    // Keys and the other layouts.
    std::vector<Key> keys;
    std::vector<uint8_t> ent6; std::vector<uint32_t> flex4; std::vector<uint64_t> packed; std::vector<uint8_t> delta;
    uint64_t bitpos = 0;
    std::vector<uint32_t> dictflex; std::vector<uint64_t> packed2; uint64_t bitpos2 = 0; std::vector<uint32_t> sorted_flex;
    auto put_bits = [](std::vector<uint64_t>& buf, uint64_t& pos, uint64_t x, int w) {
        if ((pos + w + 7) / 8 + 8 > buf.size() * 8) buf.resize(buf.size() * 2 + 1024, 0);
        auto* bytes = (uint8_t*)buf.data() + pos / 8; uint64_t cur; memcpy(&cur, bytes, 8); cur |= x << (pos % 8); memcpy(bytes, &cur, 8);
        pos += w;
    };
    const size_t K = 8, C = 4;
    for (uint64_t g = 0; g < keymap_size / K; g++) {
        uint64_t b = g * (K + C), s0, s1; memcpy(&s0, km + b, 8); memcpy(&s1, km + b + K + C, 8);
        if (s1 == s0) continue;
        for (size_t j = 0; j < K; j++) {
            uint64_t a = km[b + C + j], end = (j + 1 == K) ? s1 - s0 : km[b + C + j + 1];
            uint64_t S = end - a; if (!S) continue;
            uint64_t f = S >= 2 ? (S + 2) / 3 : 0, e = S - f;
            Key k{}; k.slot = s0 + a; k.e = (uint32_t)e;
            const uint64_t* v = slots.data() + k.slot + f;
            auto* fc = (const uint32_t*)(slots.data() + k.slot);
            // 6-byte
            k.off6 = ent6.size(); k.offf = flex4.size();
            for (uint64_t i = 0; i < e; i++) { uint8_t bts[8]; memcpy(bts, &v[i], 8); uint64_t x = (v[i] & ((1ull << 60) - 1)) | ((v[i] >> 60) << 42); memcpy(bts, &x, 8); ent6.insert(ent6.end(), bts, bts + 6); }
            if (e > 1) flex4.insert(flex4.end(), fc, fc + e);
            // bit-packed
            uint64_t mn_t = ~0ull, mx_t = 0, mn_g = ~0ull, mx_g = 0, mn_p = ~0ull, mx_p = 0;
            for (uint64_t i = 0; i < e; i++) {
                uint64_t t = v[i] >> 40 & 0xfffff, ge = v[i] >> 20 & 0xfffff, p = v[i] & 0xfffff;
                mn_t = std::min(mn_t, t); mx_t = std::max(mx_t, t); mn_g = std::min(mn_g, ge); mx_g = std::max(mx_g, ge); mn_p = std::min(mn_p, p); mx_p = std::max(mx_p, p);
            }
            k.mn_t = mn_t; k.mn_g = mn_g; k.mn_p = mn_p; k.w_t = bits(mx_t - mn_t); k.w_g = bits(mx_g - mn_g); k.w_p = bits(mx_p - mn_p);
            int w = k.w_t + k.w_g + k.w_p + 2;
            k.bitoff = bitpos;
            for (uint64_t i = 0; i < e; i++) {
                uint64_t t = (v[i] >> 40 & 0xfffff) - mn_t, ge = (v[i] >> 20 & 0xfffff) - mn_g, p = (v[i] & 0xfffff) - mn_p, fl = v[i] >> 60;
                uint64_t x = t | ge << k.w_t | p << (k.w_t + k.w_g) | fl << (k.w_t + k.w_g + k.w_p);
                if ((bitpos + w + 7) / 8 + 8 > packed.size() * 8) packed.resize(packed.size() * 2 + 1024, 0);
                auto* bytes = (uint8_t*)packed.data() + bitpos / 8; uint64_t cur; memcpy(&cur, bytes, 8); cur |= x << (bitpos % 8); memcpy(bytes, &cur, 8);
                bitpos += w;
            }
            // bit-packed with a flex dictionary: the key's distinct flex cells, each entry with the index of its cell
            k.dictoff = dictflex.size(); k.nd = 0; k.w_d = 0; k.bitoff2 = bitpos2;
            if (e > 1) {
                sorted_flex.assign(fc, fc + e); std::sort(sorted_flex.begin(), sorted_flex.end());
                sorted_flex.erase(std::unique(sorted_flex.begin(), sorted_flex.end()), sorted_flex.end());
                k.nd = sorted_flex.size(); k.w_d = bits(k.nd - 1);
                dictflex.insert(dictflex.end(), sorted_flex.begin(), sorted_flex.end());
                for (uint64_t i = 0; i < e; i++) {
                    uint64_t t = (v[i] >> 40 & 0xfffff) - mn_t, ge = (v[i] >> 20 & 0xfffff) - mn_g, p = (v[i] & 0xfffff) - mn_p, fl = v[i] >> 60;
                    uint64_t di = std::lower_bound(sorted_flex.begin(), sorted_flex.end(), fc[i]) - sorted_flex.begin();
                    uint64_t x = t | ge << k.w_t | p << (k.w_t + k.w_g) | fl << (k.w_t + k.w_g + k.w_p) | di << (k.w_t + k.w_g + k.w_p + 2);
                    put_bits(packed2, bitpos2, x, w + k.w_d);
                }
            }
            // delta
            k.deltaoff = delta.size();
            if (e > 1) {
                uint64_t pt = v[0] >> 40 & 0xfffff, pg = v[0] >> 20 & 0xfffff, pp = v[0] & 0xfffff; uint32_t pf = fc[0];
                uint8_t first[8]; uint64_t x = (v[0] & ((1ull << 60) - 1)) | ((v[0] >> 60) << 42); memcpy(first, &x, 8); delta.insert(delta.end(), first, first + 6);
                delta.insert(delta.end(), (uint8_t*)&pf, (uint8_t*)&pf + 4);
                for (uint64_t i = 1; i < e; i++) {
                    uint64_t t = v[i] >> 40 & 0xfffff, ge = v[i] >> 20 & 0xfffff, p = v[i] & 0xfffff;
                    put_varint(delta, zz((int64_t)t - (int64_t)pt));
                    delta.push_back(uint8_t((v[i] >> 60) << 1 | (ge != pg)));
                    if (ge != pg) delta.push_back(uint8_t(ge));
                    put_varint(delta, zz((int64_t)p - (int64_t)pp));
                    put_varint(delta, fc[i] ^ pf);
                    pt = t; pg = ge; pp = p; pf = fc[i];
                }
            }
            keys.push_back(k);
        }
    }
    munmap((void*)base, st.st_size);  // km pointed into the file until here
    // Reads of 8 bytes at the last entry of a layout must stay inside its buffer.
    delta.resize(delta.size() + 16, 0);
    ent6.resize(ent6.size() + 16, 0);
    flex4.resize(flex4.size() + 4, 0);
    packed.resize(packed.size() + 4, 0);
    packed2.resize(packed2.size() + 4, 0);
    dictflex.resize(dictflex.size() + 4, 0);
    uint64_t entries = 0; for (auto& k : keys) entries += k.e;
    printf("keys %zu, entries %lu; bytes: now %.1f MB, 6-byte %.1f MB, bit-packed %.1f MB (flex %.1f MB), delta %.1f MB\n", keys.size(), entries,
           slots.size() * 8 / 1e6, (ent6.size() + flex4.size() * 4) / 1e6, bitpos / 8 / 1e6 + keys.size() * 0.0, flex4.size() * 4 / 1e6, delta.size() / 1e6);

    // Sample: keys in proportion to their entries, a query flex cell per lookup (one base changed), 4 bins by e.
    std::vector<uint32_t> by_entry; by_entry.reserve(entries);
    for (uint32_t i = 0; i < keys.size(); i++) for (uint32_t j = 0; j < keys[i].e; j++) by_entry.push_back(i);
    std::mt19937_64 rng(42);
    std::vector<uint32_t> sample(lookups); std::vector<uint32_t> query(lookups);
    for (size_t i = 0; i < lookups; i++) {
        sample[i] = by_entry[rng() % by_entry.size()];
        Key const& k = keys[sample[i]];
        uint32_t f = k.e > 1 ? *(const uint32_t*)(slots.data() + k.slot) : (uint32_t)rng();
        f ^= (1u + (uint32_t)(rng() % 3)) << (2 * (rng() % 16));
        query[i] = f;
    }
    auto bin = [](uint32_t e) { return e == 1 ? 0 : e <= 8 ? 1 : e <= 64 ? 2 : 3; };
    const char* bin_names[4] = { "1", "2-8", "9-64", "65+" };
    std::vector<uint32_t> idx_by_bin[4];
    for (uint32_t i = 0; i < lookups; i++) idx_by_bin[bin(keys[sample[i]].e)].push_back(i);
    std::vector<uint16_t> sims; sims.reserve(4096); std::vector<uint64_t> dec; dec.reserve(4096); std::vector<uint8_t> dsims; dsims.reserve(4096);
    uint64_t sink = 0;

    auto best = [&](uint32_t q, auto&& flex_at, uint32_t e, uint32_t& max, uint32_t& count) {
        sims.clear(); max = 0; count = 0;
        for (uint32_t i = 0; i < e; i++) { uint32_t s = Similarity(flex_at(i), q); sims.push_back(s); if (s > max) { max = s; count = 0; } count += s == max; }
    };
    auto run = [&](const char* name, auto&& lookup) {
        printf("%-11s", name);
        double all = 0;
        for (int b = 0; b < 4; b++) {
            auto const& idx = idx_by_bin[b];
            if (idx.empty()) { printf(" %8s", "-"); continue; }
            auto t0 = std::chrono::steady_clock::now();
            for (uint32_t i : idx) lookup(i);
            double ns = std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - t0).count();
            all += ns;
            printf(" %8.1f", ns / idx.size());
        }
        printf(" %10.1f  (%lu)\n", all / lookups, sink & 0xff);
    };
    printf("%-11s %8s %8s %8s %8s %10s   ns per lookup, by entries of the key\n", "layout", bin_names[0], bin_names[1], bin_names[2], bin_names[3], "all");
    for (int rep = 0; rep < 2; rep++) {
        run("now", [&](uint32_t i) {
            Key const& k = keys[sample[i]]; const uint64_t* cells = slots.data() + k.slot;
            if (k.e == 1) { sink += cells[0] >> 40 & 0xfffff; sink += cells[0] & 0xfffff; return; }
            uint64_t f = (k.e * 3 / 2 + 2) / 3; const uint64_t* v = cells + f; auto* fc = (const uint32_t*)cells;
            uint32_t max, count; best(query[i], [&](uint32_t j) { return fc[j]; }, k.e, max, count);
            if (count > 256) return;
            for (uint32_t j = 0; j < k.e; j++) if (sims[j] == max) { sink += v[j] >> 40 & 0xfffff; sink += v[j] >> 20 & 0xfffff; sink += v[j] & 0xfffff; sink += v[j] >> 60; }
        });
        run("6-byte", [&](uint32_t i) {
            Key const& k = keys[sample[i]]; const uint8_t* e6 = ent6.data() + k.off6;
            auto entry = [&](uint32_t j) { uint64_t x; memcpy(&x, e6 + 6 * j, 8); return x & ((1ull << 48) - 1); };
            if (k.e == 1) { uint64_t x = entry(0); sink += x >> 40 & 3; sink += x & 0xfffff; return; }
            const uint32_t* fc = flex4.data() + k.offf;
            uint32_t max, count; best(query[i], [&](uint32_t j) { return fc[j]; }, k.e, max, count);
            if (count > 256) return;
            for (uint32_t j = 0; j < k.e; j++) if (sims[j] == max) { uint64_t x = entry(j); sink += x >> 40 & 0xfffff; sink += x >> 20 & 0xfffff; sink += x & 0xfffff; sink += x >> 42; }
        });
        run("bit-packed", [&](uint32_t i) {
            Key const& k = keys[sample[i]]; int w = k.w_t + k.w_g + k.w_p + 2;
            auto entry = [&](uint32_t j) { uint64_t bp = k.bitoff + (uint64_t)j * w, x; memcpy(&x, (const uint8_t*)packed.data() + bp / 8, 8); x >>= bp % 8; x &= (1ull << w) - 1;
                                           return std::tuple<uint64_t, uint64_t, uint64_t, uint64_t>{ (x & ((1ull << k.w_t) - 1)) + k.mn_t, (x >> k.w_t & ((1ull << k.w_g) - 1)) + k.mn_g, (x >> (k.w_t + k.w_g) & ((1ull << k.w_p) - 1)) + k.mn_p, x >> (w - 2) }; };
            if (k.e == 1) { uint64_t x; memcpy(&x, ent6.data() + k.off6, 8); sink += x >> 40 & 3; sink += x & 0xfffff; return; }  // a single entry stays a plain 6-byte entry
            const uint32_t* fc = flex4.data() + k.offf;
            uint32_t max, count; best(query[i], [&](uint32_t j) { return fc[j]; }, k.e, max, count);
            if (count > 256) return;
            for (uint32_t j = 0; j < k.e; j++) if (sims[j] == max) { auto [t, g, p, fl] = entry(j); sink += t + g + p + fl; }
        });
        run("packed+dict", [&](uint32_t i) {
            Key const& k = keys[sample[i]];
            if (k.e == 1) { uint64_t x; memcpy(&x, ent6.data() + k.off6, 8); sink += x >> 40 & 3; sink += x & 0xfffff; return; }
            int const w = k.w_t + k.w_g + k.w_p + 2 + k.w_d;
            const uint32_t* dict = dictflex.data() + k.dictoff;
            dsims.clear(); for (uint32_t d = 0; d < k.nd; d++) dsims.push_back((uint8_t)Similarity(dict[d], query[i]));
            auto raw = [&](uint32_t j) { uint64_t bp = k.bitoff2 + (uint64_t)j * w, x; memcpy(&x, (const uint8_t*)packed2.data() + bp / 8, 8); x >>= bp % 8; return x & ((1ull << w) - 1); };
            sims.clear(); uint32_t max = 0, count = 0;
            for (uint32_t j = 0; j < k.e; j++) { uint32_t s = dsims[raw(j) >> (w - k.w_d)]; sims.push_back(s); if (s > max) { max = s; count = 0; } count += s == max; }
            if (count > 256) return;
            for (uint32_t j = 0; j < k.e; j++) if (sims[j] == max) { uint64_t x = raw(j); sink += (x & ((1ull << k.w_t) - 1)) + k.mn_t; sink += (x >> k.w_t & ((1ull << k.w_g) - 1)) + k.mn_g; sink += (x >> (k.w_t + k.w_g) & ((1ull << k.w_p) - 1)) + k.mn_p; sink += x >> (w - k.w_d - 2) & 3; }
        });
        run("delta", [&](uint32_t i) {
            Key const& k = keys[sample[i]];
            if (k.e == 1) { uint64_t x; memcpy(&x, ent6.data() + k.off6, 8); sink += x >> 40 & 3; sink += x & 0xfffff; return; }
            const uint8_t* p = delta.data() + k.deltaoff;
            uint64_t x; memcpy(&x, p, 8); x &= (1ull << 48) - 1; p += 6; uint32_t f; memcpy(&f, p, 4); p += 4;
            uint64_t t = x >> 40 & 0xfffff, g = x >> 20 & 0xfffff, pos = x & 0xfffff, fl = x >> 42;
            dec.clear(); sims.clear(); uint32_t max = 0, count = 0;
            auto take = [&]() { uint32_t s = Similarity(f, query[i]); sims.push_back(s); if (s > max) { max = s; count = 0; } count += s == max;
                                dec.push_back(t << 40 | g << 20 | pos | fl << 60); };
            take();
            for (uint32_t j = 1; j < k.e; j++) {
                t = (uint64_t)((int64_t)t + unzz(get_varint(p)));
                uint8_t gb = *p++; fl = gb >> 1; if (gb & 1) g = *p++;
                pos = (uint64_t)((int64_t)pos + unzz(get_varint(p)));
                f ^= (uint32_t)get_varint(p);
                take();
            }
            if (count > 256) return;
            for (uint32_t j = 0; j < k.e; j++) if (sims[j] == max) { sink += dec[j] >> 40 & 0xfffff; sink += dec[j] >> 20 & 0xfffff; sink += dec[j] & 0xfffff; sink += dec[j] >> 60; }
        });
    }
    return 0;
}
