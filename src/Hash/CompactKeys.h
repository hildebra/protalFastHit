// CompactKeys.h - the key map of a loaded index in the room its used cores need (docs/claude/2026-10-13-compact-key-map.md).
//
// Seedmap's key map is addressed directly: a 16-bit offset for each of the 2^30 cores and an 8-byte control cell per 8,
// 3.2 GB whatever the reference holds, of which r226 uses 9.2% (98.7 M cores). This holds the same answer, a core's
// first slot and its number of slots, for the used cores only:
//  - presence: a bit per core, kCoresPerLine (448) to a 64-byte line whose last word is the number of used cores before
//    the line, so that a core's presence and its rank (its index among the used cores) take one cache line;
//  - starts: the first slot of each used core in rank order, and one past the last one's (the values' end), as a 64-bit
//    base per kGroup (64) ranks and a 24-bit offset from it per rank. A core has at most 3,071 slots (2,047 entries and
//    their flex cells; a core of max_key_multiplicity values or more is left out of the index), so a group spans less
//    than 2^24; the builder checks it.
// A core's slots are its start to the next rank's. About 153 MB of presence and 3 bytes a used core: ~0.45 GB at r226.
//
// A lookup counts bits (std::popcount) for its rank. protal is built for plain x86-64, where that is a call to libgcc's
// software count: the lookup functions are inlined always, so that the loops compiled for x86-64-v3 as well
// (PROTAL_CLONE_V3: ChainAnchorFinder::LookupKmers, SparsePrepare, SparseHits) count with popcnt
// (docs/claude/2026-10-14-sparse-load-key-map).
#pragma once

#include <algorithm>
#include <atomic>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>
#include <sys/mman.h>
#include <unistd.h>
#include "TargetClones.h"

namespace protal {
    class CompactKeys {
    public:
        static constexpr uint64_t kCoresPerLine = 448;  // 7 words of presence bits, then the rank count
        static constexpr uint64_t kGroup = 64;          // ranks per 64-bit base of the starts
        static constexpr uint64_t kMaxGroupSpan = uint64_t{1} << 24;
        static constexpr uint64_t kHugePage = uint64_t{1} << 21;

        CompactKeys() = default;
        CompactKeys(CompactKeys const&) = delete;
        CompactKeys& operator=(CompactKeys const&) = delete;
        ~CompactKeys() { Free(); }

        // Room for `cores` cores (2^30), none used.
        void Init(uint64_t cores) {
            Free();
            m_cores = cores;
            m_lines = (cores + kCoresPerLine - 1) / kCoresPerLine;
            m_presence = static_cast<uint64_t*>(AllocateAdvised(m_lines * 64));
            if (!m_presence) Fail("presence lines", m_lines * 64);
            std::memset(m_presence, 0, m_lines * 64);
        }

        // Core `core` has values. Threads may mark cores of the same line (atomic per word).
        void MarkUsed(uint64_t core) {
            uint64_t const line = core / kCoresPerLine, i = core % kCoresPerLine;
            std::atomic_ref<uint64_t>(m_presence[line * 8 + (i >> 6)]).fetch_or(uint64_t{1} << (i & 63), std::memory_order_relaxed);
        }

        // After every MarkUsed: each line's count of the used cores before it; returns the used cores.
        PROTAL_CLONE_V3 uint64_t CountRanks() {
            uint64_t used = 0;
            for (uint64_t line = 0; line < m_lines; line++) {
                uint64_t* const w = m_presence + line * 8;
                w[7] = used;
                for (int k = 0; k < 7; k++) used += static_cast<uint64_t>(std::popcount(w[k]));
            }
            m_used = used;
            return used;
        }

        // Room for the starts of `m_used` cores and the end.
        void AllocateStarts() {
            std::free(m_bases);
            std::free(m_offsets);
            m_groups = (m_used + 1 + kGroup - 1) / kGroup;
            m_bases = static_cast<uint64_t*>(AllocateAdvised(m_groups * sizeof(uint64_t)));
            m_offsets_bytes = 3 * (m_used + 1) + 8;  // a 4-byte load of the last offset stays inside
            m_offsets = static_cast<uint8_t*>(AllocateAdvised(m_offsets_bytes));
            if (!m_bases || !m_offsets) Fail("starts", m_groups * sizeof(uint64_t) + m_offsets_bytes);
            std::memset(m_offsets, 0, m_offsets_bytes);
        }

        // The start of rank `rank` (0 to m_used, the last the values' end): first the base of each group (SetBase for the
        // ranks that are multiples of kGroup), then every rank's offset (SetOffset). False if a rank is more than
        // kMaxGroupSpan past its group's base.
        void SetBase(uint64_t rank, uint64_t start) { m_bases[rank / kGroup] = start; }
        bool SetOffset(uint64_t rank, uint64_t start) {
            uint64_t const offset = start - m_bases[rank / kGroup];
            if (start < m_bases[rank / kGroup] || offset >= kMaxGroupSpan) return false;
            uint8_t* const p = m_offsets + 3 * rank;
            p[0] = static_cast<uint8_t>(offset);
            p[1] = static_cast<uint8_t>(offset >> 8);
            p[2] = static_cast<uint8_t>(offset >> 16);
            return true;
        }

        // A used core's rank; `core` must be used (Used()).
        [[gnu::always_inline]] inline uint64_t Rank(uint64_t core) const {
            uint64_t const* const w = m_presence + (core / kCoresPerLine) * 8;
            uint64_t const i = core % kCoresPerLine, word = i >> 6;
            uint64_t rank = w[7] + static_cast<uint64_t>(std::popcount(w[word] & ((uint64_t{1} << (i & 63)) - 1)));
            for (uint64_t k = 0; k < word; k++) rank += static_cast<uint64_t>(std::popcount(w[k]));
            return rank;
        }

        [[gnu::always_inline]] inline bool Used(uint64_t core) const {
            uint64_t const i = core % kCoresPerLine;
            return (m_presence[(core / kCoresPerLine) * 8 + (i >> 6)] >> (i & 63)) & 1;
        }

        // The first slot of rank `rank` (0 to the used cores).
        [[gnu::always_inline]] inline uint64_t Start(uint64_t rank) const {
            uint32_t v;
            std::memcpy(&v, m_offsets + 3 * rank, 4);
            return m_bases[rank / kGroup] + (v & 0xFFFFFF);
        }

        // The slots of a core's values, as Seedmap::Locate gives them: the first and how many; false if it has none.
        [[gnu::always_inline]] inline bool Locate(uint64_t core, uint64_t& start, uint64_t& slots) const {
            uint64_t const* const w = m_presence + (core / kCoresPerLine) * 8;
            uint64_t const i = core % kCoresPerLine, word = i >> 6, bit = uint64_t{1} << (i & 63);
            if (!(w[word] & bit)) return false;
            uint64_t rank = w[7] + static_cast<uint64_t>(std::popcount(w[word] & (bit - 1)));
            for (uint64_t k = 0; k < word; k++) rank += static_cast<uint64_t>(std::popcount(w[k]));
            start = Start(rank);
            slots = Start(rank + 1) - start;
            return true;
        }

        // The cache line of a core's presence and rank (a lookup's first read).
        [[gnu::always_inline]] inline void Prefetch(uint64_t core) const { __builtin_prefetch(m_presence + (core / kCoresPerLine) * 8); }

        // Its start and the next (the second read; the presence line should be there by now).
        [[gnu::always_inline]] inline void PrefetchStart(uint64_t core) const {
            if (!Used(core)) return;
            uint64_t const rank = Rank(core);
            __builtin_prefetch(m_offsets + 3 * rank);
            __builtin_prefetch(m_bases + rank / kGroup);
        }

        uint64_t UsedCores() const { return m_used; }
        uint64_t Cores() const { return m_cores; }
        uint64_t Bytes() const { return m_lines * 64 + m_groups * sizeof(uint64_t) + m_offsets_bytes; }
        uint64_t* PresenceData() const { return m_presence; }
        uint64_t PresenceBytes() const { return m_lines * 64; }
        uint8_t* OffsetData() const { return m_offsets; }
        uint64_t OffsetBytes() const { return m_offsets_bytes; }

        // The part of a block [data, data + bytes) that huge-page advice covers: its whole pages (madvise takes a
        // page-aligned start, else it refuses the advice), none for a block under a huge page.
        static std::pair<uintptr_t, uintptr_t> AdvisedRange(uintptr_t data, uint64_t bytes, uintptr_t page) {
            if (bytes < kHugePage) return { data, data };
            uintptr_t const begin = (data + page - 1) & ~(page - 1), end = (data + bytes) & ~(page - 1);
            return end > begin ? std::pair{ begin, end } : std::pair{ data, data };
        }

    private:
        [[noreturn]] static void Fail(char const* what, uint64_t bytes) {
            std::cerr << "Cannot allocate the compact key map's " << what << " (" << bytes << " bytes)" << std::endl;
            std::exit(8);
        }

        // A block for lookups at random, on huge pages where the system gives them (transparent huge pages in "madvise"
        // mode; see Seedmap::AdviseHugePages): one of a huge page or more starts on a huge page, so that all of it can,
        // and is advised before it is first touched. Until 2026-10-14 the advice took the 64-byte-aligned start as it
        // was, which madvise refuses: the compact key map was on 4 KB pages (docs/claude/2026-10-14-sparse-load-key-map).
        static void* AllocateAdvised(uint64_t bytes) {
            if (bytes < kHugePage) return std::aligned_alloc(64, (std::max<uint64_t>(bytes, 1) + 63) / 64 * 64);
            uint64_t const rounded = (bytes + kHugePage - 1) / kHugePage * kHugePage;
            void* const data = std::aligned_alloc(kHugePage, rounded);
#ifdef MADV_HUGEPAGE
            if (data) {
                auto const [begin, end] = AdvisedRange(reinterpret_cast<uintptr_t>(data), rounded, static_cast<uintptr_t>(sysconf(_SC_PAGESIZE)));
                if (end > begin) madvise(reinterpret_cast<void*>(begin), end - begin, MADV_HUGEPAGE);
            }
#endif
            return data;
        }

        void Free() {
            std::free(m_presence);
            std::free(m_bases);
            std::free(m_offsets);
            m_presence = nullptr;
            m_bases = nullptr;
            m_offsets = nullptr;
            m_cores = m_lines = m_used = m_groups = m_offsets_bytes = 0;
        }

        uint64_t* m_presence = nullptr;  // m_lines lines of 8 words: 7 of presence, the rank count before the line
        uint64_t* m_bases = nullptr;     // per kGroup ranks, the first one's start
        uint8_t* m_offsets = nullptr;    // per rank (and the end), 3 bytes: its start less its group's base
        uint64_t m_cores = 0, m_lines = 0, m_used = 0, m_groups = 0, m_offsets_bytes = 0;
    };
}
