// What the index's value blocks would cost under other in-memory layouts. Shared by index_layout.cpp
// (a raw index.prx) and index_layout_db.cpp (the index in database.protal, chunk by chunk).
//
// Layout today (Seedmap.h): per key with e values, S = e + ceil(e/2) slots of 8 bytes: ceil(S/3) slots
// of 32-bit flex cells (one per value), then e 8-byte entries (taxid 20 bits, gene 20, position 20,
// flags 2). A key with one value has one 8-byte entry and no flex cell.
//
// Measured, bytes per entry and in total:
//   now            8 B per slot
//   6-byte entries entries in 6 bytes (48 bits hold 42), flex cells 4 bytes, no half-slot padding
//   5-byte entries entries in 5 bytes (if the fields fit: see the widths printed)
//   delta          per key: first entry whole (6 B), then taxid delta (zigzag varint), gene (1 B if
//                  changed, flagged), position delta (zigzag varint), flags in the gene byte; flex
//                  cells as XOR with the previous one (varint, 1 byte if equal)
//   bit-packed     per key: 10 header bytes (widths, minima), then every entry in w_t + w_g + w_p + 2
//                  bits (the widths the key's values need above the key's minima) and its flex cell in
//                  32 bits, or in ceil(log2(distinct flex cells)) bits plus a dictionary of them
// and the share of entries in "ubiquitous groups": entries of a key whose flex cell occurs more than
// 256 (--max_key_ubiquity's default) or 2048 times in that key. Such entries can never be returned by
// a lookup (KmerLookupSM::GetFromLookup: at the maximal similarity their count exceeds the limit), so
// one representative per group would preserve every result.
#pragma once
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>

namespace layout {
    static inline int varint(uint64_t v) { int n = 1; while (v >= 128) { v >>= 7; n++; } return n; }
    static inline uint64_t zz(int64_t d) { return ((uint64_t)d << 1) ^ (uint64_t)(d >> 63); }
    static inline int bits(uint64_t v) { return v == 0 ? 0 : 64 - __builtin_clzll(v); }

    struct Stats {
        static constexpr int NB = 16;
        uint64_t nonempty = 0, singles = 0, multi_keys = 0, entries = 0, multi_entries = 0, flex_slots = 0, odd_keys = 0, slots = 0;
        uint64_t b_delta = 0, b_packed = 0, b_packed_dict = 0;
        uint64_t ubiq256_entries = 0, ubiq256_groups = 0, ubiq2048_entries = 0, ubiq2048_groups = 0;
        uint64_t max_tax = 0, max_gene = 0, max_pos = 0;
        uint64_t distinct_flex = 0;  // over the keys with several values: their distinct flex cells
        uint64_t bin_keys[NB] = {}, bin_entries[NB] = {};
        std::vector<uint64_t> ent; std::vector<uint32_t> flex, sorted;

        // One key's S slots at `cells` (flex cells first, then the entries).
        void AddKey(const uint64_t* cells, uint64_t S) {
            if (!S) return;
            nonempty++; slots += S;
            uint64_t f = S >= 2 ? (S + 2) / 3 : 0, e = S - f;
            entries += e; flex_slots += f;
            int bin = e == 1 ? 0 : std::min(NB - 1, bits(e - 1));  // 1 | 2 | 3-4 | 5-8 | ...
            bin_keys[bin]++; bin_entries[bin] += e;
            const uint64_t* v = cells + f;
            for (uint64_t i = 0; i < e; i++) {
                max_tax = std::max(max_tax, v[i] >> 40 & 0xfffff); max_gene = std::max(max_gene, v[i] >> 20 & 0xfffff);
                max_pos = std::max(max_pos, v[i] & 0xfffff);
            }
            if (S == 1) { singles++; return; }
            multi_keys++; multi_entries += e; odd_keys += e & 1;
            ent.assign(v, v + e);
            auto* fc = (const uint32_t*)cells;
            flex.assign(fc, fc + e);
            // delta
            uint64_t cost = 6 + 4, pt = ent[0] >> 40 & 0xfffff, pg = ent[0] >> 20 & 0xfffff, pp = ent[0] & 0xfffff;
            uint32_t pf = flex[0];
            for (uint64_t i = 1; i < e; i++) {
                uint64_t t = ent[i] >> 40 & 0xfffff, ge = ent[i] >> 20 & 0xfffff, p = ent[i] & 0xfffff;
                cost += varint(zz((int64_t)t - (int64_t)pt)) + 1 + (ge != pg ? 1 : 0) + varint(zz((int64_t)p - (int64_t)pp));
                cost += varint(flex[i] ^ pf);
                pt = t; pg = ge; pp = p; pf = flex[i];
            }
            b_delta += cost;
            // bit-packed
            uint64_t mn_t = ~0ull, mx_t = 0, mn_g = ~0ull, mx_g = 0, mn_p = ~0ull, mx_p = 0;
            for (uint64_t x : ent) {
                uint64_t t = x >> 40 & 0xfffff, ge = x >> 20 & 0xfffff, p = x & 0xfffff;
                mn_t = std::min(mn_t, t); mx_t = std::max(mx_t, t); mn_g = std::min(mn_g, ge); mx_g = std::max(mx_g, ge);
                mn_p = std::min(mn_p, p); mx_p = std::max(mx_p, p);
            }
            int w = bits(mx_t - mn_t) + bits(mx_g - mn_g) + bits(mx_p - mn_p) + 2;
            sorted = flex; std::sort(sorted.begin(), sorted.end());
            uint64_t distinct = std::unique(sorted.begin(), sorted.end()) - sorted.begin();
            distinct_flex += distinct;
            b_packed += 10 + (e * (w + 32) + 7) / 8;
            b_packed_dict += 10 + distinct * 4 + (e * (w + bits(distinct - 1)) + 7) / 8;
            // ubiquitous groups: runs of equal flex cells in the sorted copy
            std::sort(flex.begin(), flex.end());
            for (uint64_t i = 0; i < e;) {
                uint64_t k = i; while (k < e && flex[k] == flex[i]) k++;
                uint64_t n = k - i;
                if (n > 256) { ubiq256_entries += n; ubiq256_groups++; }
                if (n > 2048) { ubiq2048_entries += n; ubiq2048_groups++; }
                i = k;
            }
        }

        void Merge(Stats const& o) {
            nonempty += o.nonempty; singles += o.singles; multi_keys += o.multi_keys; entries += o.entries;
            multi_entries += o.multi_entries; flex_slots += o.flex_slots; odd_keys += o.odd_keys; slots += o.slots;
            b_delta += o.b_delta; b_packed += o.b_packed; b_packed_dict += o.b_packed_dict;
            ubiq256_entries += o.ubiq256_entries; ubiq256_groups += o.ubiq256_groups;
            ubiq2048_entries += o.ubiq2048_entries; ubiq2048_groups += o.ubiq2048_groups;
            distinct_flex += o.distinct_flex;
            max_tax = std::max(max_tax, o.max_tax); max_gene = std::max(max_gene, o.max_gene); max_pos = std::max(max_pos, o.max_pos);
            for (int i = 0; i < NB; i++) { bin_keys[i] += o.bin_keys[i]; bin_entries[i] += o.bin_entries[i]; }
        }

        void Print(uint64_t keymap_size) const {
            double now = slots * 8.0, single_b = singles * 8.0;
            printf("keys %lu, non-empty %lu (%.1f%%); slots %lu = %.3f GB; entries %lu (%.1f per non-empty key), flex slots %lu\n",
                   keymap_size, nonempty, 100.0 * nonempty / keymap_size, slots, now / 1e9, entries, (double)entries / nonempty, flex_slots);
            printf("single-value keys %lu (%.1f%% of entries, %.1f%% of bytes); keys with several values %lu holding %lu entries (%.1f per key; %lu keys with an odd count pad half a slot)\n",
                   singles, 100.0 * singles / entries, 100.0 * single_b / now, multi_keys, multi_entries, (double)multi_entries / multi_keys, odd_keys);
            printf("distinct flex cells in the keys with several values: %lu (%.1f%% of their %lu entries; the rest repeat a cell of their key)\n",
                   distinct_flex, 100.0 * distinct_flex / multi_entries, multi_entries);
            printf("field maxima: taxid %lu (%d bits), gene %lu (%d bits), position %lu (%d bits); 2 flag bits\n",
                   max_tax, bits(max_tax), max_gene, bits(max_gene), max_pos, bits(max_pos));
            printf("%-22s %12s %14s %9s\n", "entries per key", "keys", "entries", "% entries");
            for (int i = 0; i < NB; i++) if (bin_keys[i]) {
                char lab[32]; if (i == 0) snprintf(lab, 32, "1"); else if (i == 1) snprintf(lab, 32, "2"); else snprintf(lab, 32, "%lu-%lu", (1ul << (i - 1)) + 1, 1ul << i);
                printf("%-22s %12lu %14lu %8.2f%%\n", lab, bin_keys[i], bin_entries[i], 100.0 * bin_entries[i] / entries);
            }
            auto line = [&](const char* name, double bytes) {
                printf("%-42s %10.1f MB  %5.2f B/entry  %6.1f%% of now\n", name, bytes / 1e6, bytes / entries, 100.0 * bytes / now);
            };
            printf("\nlayout (all keys)\n");
            line("now (8 B slots)", now);
            line("6-byte entries + 4-byte flex", 6.0 * entries + 4.0 * multi_entries);
            line("5-byte entries + 4-byte flex", 5.0 * entries + 4.0 * multi_entries);
            line("delta (singles 6 B)", single_b * 6 / 8 + b_delta);
            line("bit-packed, flex 32 bits (singles 6 B)", single_b * 6 / 8 + b_packed);
            line("bit-packed, flex dictionary (singles 6 B)", single_b * 6 / 8 + b_packed_dict);
            printf("\nubiquitous groups (entries of a key with the same flex cell): > 256: %lu groups, %lu entries (%.2f%% of entries, %.1f MB of 12 B cells); > 2048: %lu groups, %lu entries (%.2f%%)\n",
                   ubiq256_groups, ubiq256_entries, 100.0 * ubiq256_entries / entries, ubiq256_entries * 12.0 / 1e6, ubiq2048_groups, ubiq2048_entries, 100.0 * ubiq2048_entries / entries);
        }
    };
}
