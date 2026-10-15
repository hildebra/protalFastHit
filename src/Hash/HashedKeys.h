// HashedKeys.h - the key map of a loaded index as a hash table of its used cores: one 64-byte line a lookup, as the
// direct key map reads, in a third of its memory (docs/claude/2026-10-15-perf-0.8.4a, section 10).
//
// The direct key map (3.2 GB whatever the reference) answers a lookup from one line; the compact one (CompactKeys.h,
// ~0.45 GB at r226) from two dependent ones, the core's presence and rank, then its start, which cost r226 runs ~25 ns
// a lookup (pe seeding +18%). Here each used core has a 64-bit entry in a table of 64-byte lines, 8 entries a line. A
// core's home line is chosen by a bijective mix of the core, so that the cores of related genomes, which cluster in
// core order (a 128-core bucket per line overflowed for 17.5% of the reduced r226's cores), spread over the lines. An
// entry holds the mixed core's bits below the line index (its tag: with the line, the core), how many lines past its
// home it lies (0-31), its first slot and its slot count (12 bits, never 0, so that an empty entry is 0). A core whose
// home line filled up lies in a following line (linear probing), so a lookup reads on only while the line it read is
// full: a line with an empty entry holds every core of that home that lies in or past it. 2^b lines for the used cores
// at a load of 0.75 at most: 1.07 GB at r226 (98.7 M cores, b = 24, load 0.735).
#pragma once

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>
#include <sys/mman.h>
#include <unistd.h>
#include "CompactKeys.h"

namespace protal {
    class HashedKeys {
    public:
        static constexpr uint64_t kEntries = 8;          // 64-bit entries a line
        static constexpr unsigned kDispBits = 5;          // lines past its home a core may lie: 0-31
        static constexpr uint64_t kMaxDisp = (uint64_t{1} << kDispBits) - 1;
        static constexpr unsigned kCountShift = 52;       // the slot count's 12 bits (a core has at most 3,071 slots)
        static constexpr double kMaxLoad = 0.75;

        HashedKeys() = default;
        HashedKeys(HashedKeys const&) = delete;
        HashedKeys& operator=(HashedKeys const&) = delete;
        ~HashedKeys() { std::free(m_lines); }

        // A bijection of the cores [0, 2^core_bits): odd multipliers and xor-shifts, so that neighbouring cores get lines far
        // apart.
        [[gnu::always_inline]] inline uint64_t Mix(uint64_t core) const {
            uint64_t x = core & m_core_mask;
            x = (x * 0x2545F491ull) & m_core_mask;
            x ^= x >> 15;
            x = (x * 0x1B873593ull) & m_core_mask;
            x ^= x >> 13;
            return x;
        }

        // The table of `compact`'s used cores (cores of `core_bits` bits) on `threads` threads, at most `max_load` used cores
        // per entry of the home lines (the lines a power of two; tests crowd it). The cores are put in partitions of
        // lines by their home line and each partition filled on its own thread in core order; a core that would lie past
        // its partition's last line is put after all partitions, in core order too, so that the table is the same on any
        // number of threads. False (nothing kept) if a core would lie more than kMaxDisp lines past its home, or its start
        // does not fit its entry.
        bool Build(CompactKeys const& compact, unsigned core_bits, int threads, double max_load = kMaxLoad) {
            Free();
            m_core_bits = core_bits;
            m_core_mask = (uint64_t{1} << core_bits) - 1;
            m_used = compact.UsedCores();
            m_bits = 6;
            while (m_bits < core_bits && static_cast<double>(m_used) > max_load * kEntries * static_cast<double>(uint64_t{1} << m_bits)) m_bits++;
            m_tag_bits = core_bits - m_bits;
            m_tag_mask = (uint64_t{1} << m_tag_bits) - 1;
            m_start_shift = m_tag_bits + kDispBits;
            if (m_start_shift >= kCountShift) return Fail();
            m_start_mask = (uint64_t{1} << (kCountShift - m_start_shift)) - 1;
            m_low_mask = (uint64_t{1} << m_start_shift) - 1;
            m_home_lines = uint64_t{1} << m_bits;
            m_lines_total = m_home_lines + kMaxDisp;
            m_lines = static_cast<uint64_t*>(AllocateAdvised(m_lines_total * 64));
            if (!m_lines) return Fail();
            threads = std::max(threads, 1);
            int64_t const zero_chunks = static_cast<int64_t>((m_lines_total + (uint64_t{1} << 16) - 1) >> 16);
#pragma omp parallel for num_threads(threads) schedule(static)
            for (int64_t z = 0; z < zero_chunks; z++) {
                uint64_t const first = static_cast<uint64_t>(z) << 16, last = std::min(m_lines_total, first + (uint64_t{1} << 16));
                std::memset(m_lines + first * kEntries, 0, (last - first) * 64);
            }

            // The used cores by partition of home lines, in core order within each: counted per range of presence lines,
            // then placed. A range's used cores are read from the presence words, a set bit at a time.
            uint64_t const cores = uint64_t{1} << core_bits;
            unsigned const part_bits = std::min<unsigned>(m_bits, 12);
            uint64_t const parts = uint64_t{1} << part_bits, part_shift = m_bits - part_bits;
            uint64_t const presence_lines = (cores + CompactKeys::kCoresPerLine - 1) / CompactKeys::kCoresPerLine;
            uint64_t const ranges = std::min<uint64_t>(256, presence_lines);
            uint64_t const range_lines = (presence_lines + ranges - 1) / ranges;
            uint64_t const* const presence = compact.PresenceData();
            std::vector<uint32_t> counts(ranges * parts, 0);
            auto each_used = [&](uint64_t r, auto&& on_core) {
                uint64_t const first = r * range_lines, last = std::min(presence_lines, first + range_lines);
                for (uint64_t line = first; line < last; line++) {
                    for (uint64_t w = 0; w < 7; w++) {
                        for (uint64_t bits = presence[line * 8 + w]; bits != 0; bits &= bits - 1) {
                            uint64_t const core = line * CompactKeys::kCoresPerLine + w * 64 + static_cast<uint64_t>(__builtin_ctzll(bits));
                            if (core < cores) on_core(core);
                        }
                    }
                }
            };
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
            for (int64_t r = 0; r < static_cast<int64_t>(ranges); r++) {
                uint32_t* const c = counts.data() + static_cast<uint64_t>(r) * parts;
                each_used(static_cast<uint64_t>(r), [&](uint64_t core) { c[(Mix(core) >> m_tag_bits) >> part_shift]++; });
            }
            std::vector<uint64_t> part_first(parts + 1, 0);
            std::vector<uint64_t> offset(ranges * parts, 0);
            for (uint64_t p = 0, at = 0; p < parts; p++) {
                part_first[p] = at;
                for (uint64_t r = 0; r < ranges; r++) {
                    offset[r * parts + p] = at;
                    at += counts[r * parts + p];
                }
                part_first[p + 1] = at;
            }
            std::vector<uint32_t>().swap(counts);
            std::vector<uint32_t> by_part(m_used);
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
            for (int64_t r = 0; r < static_cast<int64_t>(ranges); r++) {
                uint64_t* const o = offset.data() + static_cast<uint64_t>(r) * parts;
                each_used(static_cast<uint64_t>(r), [&](uint64_t core) { by_part[o[(Mix(core) >> m_tag_bits) >> part_shift]++] = static_cast<uint32_t>(core); });
            }
            std::vector<uint64_t>().swap(offset);

            // Each partition filled on its own; the cores that would lie past it kept for after.
            std::vector<std::vector<uint32_t>> spilled(parts);
            std::atomic<bool> fits{ true };
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
            for (int64_t p = 0; p < static_cast<int64_t>(parts); p++) {
                uint64_t const end_line = (static_cast<uint64_t>(p) + 1) << part_shift;
                for (uint64_t i = part_first[static_cast<uint64_t>(p)]; i < part_first[static_cast<uint64_t>(p) + 1]; i++) {
                    int const put = Put(compact, by_part[i], end_line);
                    if (put < 0) fits = false;
                    else if (put == 0) spilled[static_cast<uint64_t>(p)].push_back(by_part[i]);
                }
            }
            for (auto const& cores_spilled : spilled) {
                for (uint32_t core : cores_spilled) {
                    if (Put(compact, core, m_lines_total) <= 0) fits = false;
                }
            }
            if (!fits) return Fail();
            m_max_disp = m_max_disp_seen.load();
            for (uint64_t d = 1; d <= kMaxDisp; d++) m_past += m_by_disp[d].load();
            return true;
        }

        // The slots of a core's values, as Seedmap::Locate gives them: the first and how many; false if it has none. A line's
        // 8 entries are compared without a branch (at most one is the core's: its entry is or-ed in; an empty one, 0, adds
        // nothing), so that the entry's place in the line costs no misprediction, and the loop compiles to vector
        // compares where the caller is cloned for x86-64-v3.
        [[gnu::always_inline]] inline bool Locate(uint64_t core, uint64_t& start, uint64_t& slots) const {
            uint64_t const p = Mix(core), tag = p & m_tag_mask;
            uint64_t const* w = m_lines + (p >> m_tag_bits) * kEntries;
            for (uint64_t d = 0;; d++, w += kEntries) {
                uint64_t const want = tag | d << m_tag_bits;
                uint64_t hit = 0, empty = 0;
                for (uint64_t e = 0; e < kEntries; e++) {
                    uint64_t const x = w[e];
                    hit |= (x & m_low_mask) == want ? x : 0;
                    empty |= x == 0;
                }
                if (hit != 0) {
                    start = (hit >> m_start_shift) & m_start_mask;
                    slots = hit >> kCountShift;
                    return true;
                }
                if (empty != 0 || d >= m_max_disp) return false;
            }
        }

        // The line of a core's home (a lookup's one read, mostly).
        [[gnu::always_inline]] inline void Prefetch(uint64_t core) const { __builtin_prefetch(m_lines + (Mix(core) >> m_tag_bits) * kEntries); }

        uint64_t UsedCores() const { return m_used; }
        uint64_t Bytes() const { return m_lines_total * 64; }
        unsigned LineBits() const { return m_bits; }
        double Load() const { return m_home_lines ? static_cast<double>(m_used) / static_cast<double>(kEntries * m_home_lines) : 0.0; }
        uint64_t PastHome() const { return m_past; }        // cores in a line past their home's (a lookup of them reads two)
        uint64_t MaxDisplacement() const { return m_max_disp; }
        uint64_t const* LineData() const { return m_lines; }

    private:
        // Core `core` (used in `compact`) into the first line from its home with an empty entry, before line `end_line`:
        // 1 put, 0 not before end_line (put later), -1 never (more than kMaxDisp lines past its home, or its start too
        // large for its entry). Lines of a partition are written by its thread alone; the cores put after them by one.
        int Put(CompactKeys const& compact, uint64_t core, uint64_t end_line) {
            uint64_t start = 0, slots = 0;
            if (!compact.Locate(core, start, slots) || slots == 0 || slots >= (uint64_t{1} << (64 - kCountShift)) || start > m_start_mask) return -1;
            uint64_t const p = Mix(core), home = p >> m_tag_bits, tag = p & m_tag_mask;
            for (uint64_t d = 0; d <= kMaxDisp; d++) {
                if (home + d >= end_line) return 0;
                uint64_t* const w = m_lines + (home + d) * kEntries;
                for (uint64_t e = 0; e < kEntries; e++) {
                    if (w[e] != 0) continue;
                    w[e] = tag | d << m_tag_bits | start << m_start_shift | slots << kCountShift;
                    m_by_disp[d].fetch_add(1, std::memory_order_relaxed);
                    uint64_t seen = m_max_disp_seen.load(std::memory_order_relaxed);
                    while (d > seen && !m_max_disp_seen.compare_exchange_weak(seen, d, std::memory_order_relaxed)) {}
                    return 1;
                }
            }
            return -1;
        }

        bool Fail() {
            Free();
            return false;
        }

        // A block on huge pages where the system gives them, as CompactKeys allocates (CompactKeys::AdvisedRange).
        static void* AllocateAdvised(uint64_t bytes) {
            constexpr uint64_t kHugePage = CompactKeys::kHugePage;
            if (bytes < kHugePage) return std::aligned_alloc(64, (std::max<uint64_t>(bytes, 1) + 63) / 64 * 64);
            uint64_t const rounded = (bytes + kHugePage - 1) / kHugePage * kHugePage;
            void* const data = std::aligned_alloc(kHugePage, rounded);
#ifdef MADV_HUGEPAGE
            if (data) madvise(data, rounded, MADV_HUGEPAGE);
#endif
            return data;
        }

        void Free() {
            std::free(m_lines);
            m_lines = nullptr;
            m_used = m_home_lines = m_lines_total = m_past = m_max_disp = 0;
            m_max_disp_seen = 0;
            for (auto& n : m_by_disp) n = 0;
        }

        uint64_t* m_lines = nullptr;
        unsigned m_core_bits = 0, m_bits = 0, m_tag_bits = 0, m_start_shift = 0;
        uint64_t m_core_mask = 0, m_tag_mask = 0, m_start_mask = 0, m_low_mask = 0;
        uint64_t m_used = 0, m_home_lines = 0, m_lines_total = 0, m_past = 0, m_max_disp = 0;
        std::atomic<uint64_t> m_max_disp_seen{ 0 };
        std::atomic<uint64_t> m_by_disp[kMaxDisp + 1]{};
    };
}
