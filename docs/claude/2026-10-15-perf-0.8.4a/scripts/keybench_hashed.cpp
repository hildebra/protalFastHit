// keybench - experiment only (sparse load and key map profiling, 2026-10-14): protal's index loads and the key map step
// of a lookup (a k-mer's core -> the slots of its values), in protal's layouts and in proposed ones.
//
//   keybench load DATABASE.protal LEVEL COMPACT THREADS [GRIDBITS]
//       one index load as protal makes it (LEVEL > 0: a sparse load of that grid level; COMPACT 1: --compact_keys),
//       timed; PROTAL_LOAD_PHASES=1 prints the phases (an instrumented Seedmap.h).
//   keybench keys DATABASE.protal READS[.gz] MAX_READS PAIRED(0/1) ROUNDS [LEVEL]
//       the reads' k-mers (as protal takes them) located in each layout, one thread, the layouts alternating per batch
//       of reads: D direct key map (protal's default), C compact (CompactKeys.h, --compact_keys), M compact with each
//       group's base in its offsets' line, B one 64-byte line per 128 cores (presence, base, end offsets; overflow
//       to a side table). Each layout's lookups give the same slots (checked). PAIRED: the next read's keys fetched
//       while this one is looked up (Classify's short-read loop), else one stream per read (long reads).
//   keybench synth DENSITY MEAN_SLOTS LOOKUPS ROUNDS [READ_KMERS PAIRED]
//       the same on a synthetic key map: DENSITY of the 2^30 cores used (r226: 0.092), slot counts geometric with
//       MEAN_SLOTS; the lookups 99% of used cores (r226: 99.2% of a read's cores have values), in "reads" of READ_KMERS.
#include <bits/stdc++.h>
#include <omp.h>
#include <sys/resource.h>

#include "Hash/Seedmap.h"
#include "Hash/CompactKeys.h"
#include "SequenceUtils/KmerIterator.h"
#include "SequenceUtils/Minimizer.h"
#include "SequenceUtils/SeqReader.h"
#include "SequenceUtils/GenomeLoader.h"
#include "ThreadedGzStream.h"
#include "Utilities/Database.h"

using namespace protal;

namespace {
    double Now() {
        timespec ts;
        clock_gettime(CLOCK_MONOTONIC, &ts);
        return static_cast<double>(ts.tv_sec) + 1e-9 * static_cast<double>(ts.tv_nsec);
    }

    std::string Status(char const* field) {
        std::ifstream is("/proc/self/status");
        std::string line;
        while (std::getline(is, line)) {
            if (line.rfind(field, 0) == 0) return line;
        }
        return "";
    }

    // A block for a layout's lines on huge pages (2 MB-aligned, advised before first touched); with KEYBENCH_OLD_ADVICE as
    // CompactKeys allocated before the fix (64-byte aligned, the advice refused).
    void* HugeAlloc(uint64_t bytes) {
        uint64_t const huge = uint64_t{1} << 21;
        if (std::getenv("KEYBENCH_OLD_ADVICE") || bytes < huge) return std::aligned_alloc(64, (bytes + 63) / 64 * 64);
        uint64_t const rounded = (bytes + huge - 1) / huge * huge;
        void* const p = std::aligned_alloc(huge, rounded);
        madvise(p, rounded, MADV_HUGEPAGE);
        return p;
    }

    // The direct key map's cells on huge pages, as protal holds them (Seedmap::AllocateKeymap and AdviseHugePages).
    struct HugeCells {
        uint16_t* p = nullptr;
        uint64_t n = 0;
        explicit HugeCells(uint64_t cells) : n(cells) {
            uint64_t const huge = uint64_t{1} << 21, bytes = (cells * 2 + huge - 1) / huge * huge;
            p = static_cast<uint16_t*>(std::aligned_alloc(huge, bytes));
            madvise(p, bytes, MADV_HUGEPAGE);
            std::memset(p, 0, cells * 2);
        }
        ~HugeCells() { std::free(p); }
        uint16_t* data() { return p; }
        uint64_t size() const { return n; }
    };

    // Anonymous memory on huge pages, from /proc/self/smaps_rollup.
    std::string HugePages() {
        std::ifstream is("/proc/self/smaps_rollup");
        std::string line;
        while (std::getline(is, line)) {
            if (line.rfind("AnonHugePages", 0) == 0) return line;
        }
        return "";
    }

    // The key map interface every layout gives: a core's slots (first, count), and what to fetch ahead.
    // Direct: protal's key map (Seedmap::Locate on a direct map; Seedmap::PrefetchKey).
    struct DirectKeys {
        uint16_t const* km = nullptr;  // 2^27 blocks of 12 cells, then the final control block
        static uint64_t BlockIndex(uint64_t core) { return (core >> 3) * 12; }
        void Prefetch(uint64_t core) const {
            char const* const block = reinterpret_cast<char const*>(km + BlockIndex(core));
            __builtin_prefetch(block);
            __builtin_prefetch(block + 31);
        }
        void PrefetchStart(uint64_t) const {}
        bool Locate(uint64_t core, uint64_t& start, uint64_t& slots) const {
            uint64_t const b = BlockIndex(core), e = b + 12;
            uint64_t s0, s1;
            std::memcpy(&s0, km + b, 8);
            std::memcpy(&s1, km + e, 8);
            if (s0 == s1) return false;
            uint64_t const k = b + 4 + (core & 7);
            start = s0 + km[k];
            uint64_t const end = (k + 1) == e ? s1 : s0 + km[k + 1];
            if (end <= start) return false;
            slots = end - start;
            return true;
        }
        static constexpr bool kTwoStep = false;
        static constexpr size_t kAhead = 8, kStartsAhead = 0;
    };

    // D2: the direct key map, a core's end read from its own block for the first 7 cores of a block (only the 8th needs
    // the next block's base), and the line of the bytes it reads prefetched: fewer lines a lookup.
    struct DirectKeys2 {
        uint16_t const* km = nullptr;
        void Prefetch(uint64_t core) const {
            char const* const block = reinterpret_cast<char const*>(km + (core >> 3) * 12);
            uint64_t const j = core & 7;
            __builtin_prefetch(block);
            __builtin_prefetch(block + (j == 7 ? 31 : 8 + 2 * j + 3));
        }
        void PrefetchStart(uint64_t) const {}
        bool Locate(uint64_t core, uint64_t& start, uint64_t& slots) const {
            uint64_t const b = (core >> 3) * 12, j = core & 7;
            uint64_t s0, end;
            std::memcpy(&s0, km + b, 8);
            start = s0 + km[b + 4 + j];
            if (j < 7) {
                end = s0 + km[b + 5 + j];
            } else {
                std::memcpy(&end, km + b + 12, 8);
            }
            if (end <= start) return false;
            slots = end - start;
            return true;
        }
        static constexpr bool kTwoStep = false;
        static constexpr size_t kAhead = 8, kStartsAhead = 0;
    };

    // C: CompactKeys as protal has it.
    struct CompactNow {
        CompactKeys const* c = nullptr;
        void Prefetch(uint64_t core) const { c->Prefetch(core); }
        void PrefetchStart(uint64_t core) const { c->PrefetchStart(core); }
        bool Locate(uint64_t core, uint64_t& start, uint64_t& slots) const { return c->Locate(core, start, slots); }
        static constexpr bool kTwoStep = true;
        static constexpr size_t kAhead = 16, kStartsAhead = 8;
    };

    // Cl: the same with twice the lead (presence lines 32 k-mers ahead, starts 16).
    struct CompactLead : CompactNow {
        static constexpr size_t kAhead = 32, kStartsAhead = 16;
    };

    // M: CompactKeys' presence lines, and the starts in 64-byte lines of kRanks ranks: a 40-bit base (the start of the
    // group's first rank) and kRanks 24-bit end offsets (each rank's end less the base), so that a core's start and end
    // come from one line, and no second array (the bases) is read.
    struct MergedStarts {
        static constexpr uint64_t kRanks = 19;  // 5 + 19 * 3 = 62 bytes
        CompactKeys const* c = nullptr;
        uint8_t* lines = nullptr;
        uint64_t n_lines = 0;
        bool ok = true;
        ~MergedStarts() { std::free(lines); }
        void Build(CompactKeys const& compact) {
            c = &compact;
            uint64_t const used = compact.UsedCores();
            n_lines = (used + kRanks - 1) / kRanks;
            lines = static_cast<uint8_t*>(HugeAlloc(std::max<uint64_t>(n_lines, 1) * 64));
            std::memset(lines, 0, std::max<uint64_t>(n_lines, 1) * 64);
            for (uint64_t r = 0; r < used; r++) {
                uint8_t* const line = lines + (r / kRanks) * 64;
                uint64_t const base = compact.Start((r / kRanks) * kRanks);
                if (r % kRanks == 0) std::memcpy(line, &base, 5);
                uint64_t const off = compact.Start(r + 1) - base;
                if (off >= (uint64_t{1} << 24)) ok = false;
                std::memcpy(line + 5 + 3 * (r % kRanks), &off, 3);
            }
        }
        void Prefetch(uint64_t core) const { c->Prefetch(core); }
        void PrefetchStart(uint64_t core) const {
            if (!c->Used(core)) return;
            __builtin_prefetch(lines + (c->Rank(core) / kRanks) * 64);
        }
        bool Locate(uint64_t core, uint64_t& start, uint64_t& slots) const {
            uint64_t const* const w = c->PresenceData() + (core / CompactKeys::kCoresPerLine) * 8;
            uint64_t const i = core % CompactKeys::kCoresPerLine, word = i >> 6, bit = uint64_t{1} << (i & 63);
            if (!(w[word] & bit)) return false;
            uint64_t rank = w[7] + static_cast<uint64_t>(std::popcount(w[word] & (bit - 1)));
            for (uint64_t k = 0; k < word; k++) rank += static_cast<uint64_t>(std::popcount(w[k]));
            uint8_t const* const line = lines + (rank / kRanks) * 64;
            uint64_t const j = rank % kRanks;
            uint64_t base = 0, end = 0, begin = 0;
            std::memcpy(&base, line, 8);
            base &= (uint64_t{1} << 40) - 1;
            uint32_t v;
            std::memcpy(&v, line + 5 + 3 * j, 4);
            end = v & 0xFFFFFF;
            if (j > 0) {
                std::memcpy(&v, line + 5 + 3 * (j - 1), 4);
                begin = v & 0xFFFFFF;
            }
            start = base + begin;
            slots = end - begin;
            return true;
        }
        uint64_t Bytes() const { return c->PresenceBytes() + n_lines * 64; }
        static constexpr bool kTwoStep = true;
        static constexpr size_t kAhead = 16, kStartsAhead = 8;
    };

    // B: one 64-byte line per 128 cores: their presence bits (2 words), a header (the start of the line's first used
    // core in 40 bits; bit 40 set: overflow, its id above), and the end offsets of up to 20 used cores (16 bits, each
    // used core's end less the base). A line of more than 20 used cores, or whose ends span 2^16 slots, keeps its ends
    // in a side table (32-bit, from the header's overflow id). One line a lookup, no rank count across lines.
    struct BucketLines {
        static constexpr uint64_t kCores = 128, kSlots = 20;
        uint64_t* lines = nullptr;  // 8 words a line
        uint64_t n_lines = 0;
        std::vector<uint64_t> overflow_first;  // per overflow line, its first end in overflow_ends
        std::vector<uint32_t> overflow_ends;
        uint64_t used_cores = 0, overflow_cores = 0;  // the used cores, and those in overflow lines
        std::array<uint64_t, 65> histogram{};         // lines by used cores (64: 64 or more), the empty ones not counted
        ~BucketLines() { std::free(lines); }
        // The lines by their used cores, in bands, and the share of used cores in overflow lines (a lookup's second read).
        std::string Describe() const {
            uint64_t lines_used = 0, bands[5] = {};
            for (size_t i = 1; i < histogram.size(); i++) {
                lines_used += histogram[i];
                bands[i <= 4 ? 0 : i <= 10 ? 1 : i <= 20 ? 2 : i <= 32 ? 3 : 4] += histogram[i];
            }
            char s[512];
            std::snprintf(s, sizeof(s), "%" PRIu64 " of %" PRIu64 " lines used; by used cores 1-4 %" PRIu64 ", 5-10 %" PRIu64 ", 11-20 %" PRIu64
                          ", 21-32 %" PRIu64 ", 33+ %" PRIu64 "; %.3f%% of used cores in overflow lines",
                          lines_used, n_lines, bands[0], bands[1], bands[2], bands[3], bands[4],
                          used_cores ? 100.0 * static_cast<double>(overflow_cores) / static_cast<double>(used_cores) : 0.0);
            return s;
        }
        template<typename LocateFn>
        void Build(uint64_t cores, LocateFn&& locate) {
            n_lines = cores / kCores;
            lines = static_cast<uint64_t*>(HugeAlloc(n_lines * 64));
            std::memset(lines, 0, n_lines * 64);
            std::vector<uint64_t> starts, ends;
            for (uint64_t l = 0; l < n_lines; l++) {
                uint64_t* const w = lines + l * 8;
                starts.clear();
                ends.clear();
                for (uint64_t i = 0; i < kCores; i++) {
                    uint64_t start, slots;
                    if (!locate(l * kCores + i, start, slots)) continue;
                    w[i >> 6] |= uint64_t{1} << (i & 63);
                    starts.push_back(start);
                    ends.push_back(start + slots);
                }
                if (starts.empty()) continue;
                used_cores += starts.size();
                histogram[std::min<size_t>(starts.size(), histogram.size() - 1)]++;
                uint64_t const base = starts[0];
                bool fits = starts.size() <= kSlots && ends.back() - base < (uint64_t{1} << 16);
                // The used cores' values follow each other (a core's end is the next one's start).
                for (size_t k = 1; k < starts.size(); k++) fits = fits && starts[k] == ends[k - 1];
                if (fits) {
                    w[2] = base;
                    auto* const e = reinterpret_cast<uint16_t*>(w + 3);
                    for (size_t k = 0; k < ends.size(); k++) e[k] = static_cast<uint16_t>(ends[k] - base);
                } else {
                    overflow_cores += starts.size();
                    w[2] = base | uint64_t{1} << 40 | static_cast<uint64_t>(overflow_first.size()) << 41;
                    overflow_first.push_back(overflow_ends.size());
                    for (size_t k = 0; k < ends.size(); k++) {
                        // Ends relative to the base, and each start where it is not the previous end (never, in protal's
                        // key map) would need its own field: the bench checks the slots against the direct map.
                        overflow_ends.push_back(static_cast<uint32_t>(ends[k] - base));
                    }
                }
            }
        }
        void Prefetch(uint64_t core) const { __builtin_prefetch(lines + (core / kCores) * 8); }
        void PrefetchStart(uint64_t) const {}
        bool Locate(uint64_t core, uint64_t& start, uint64_t& slots) const {
            uint64_t const* const w = lines + (core / kCores) * 8;
            uint64_t const i = core % kCores, word = i >> 6, bit = uint64_t{1} << (i & 63);
            if (!(w[word] & bit)) return false;
            uint64_t const j = static_cast<uint64_t>(std::popcount(w[word] & (bit - 1))) + (word ? static_cast<uint64_t>(std::popcount(w[0])) : 0);
            uint64_t const h = w[2], base = h & ((uint64_t{1} << 40) - 1);
            uint64_t begin, end;
            if (!(h >> 40 & 1)) {
                auto const* const e = reinterpret_cast<uint16_t const*>(w + 3);
                end = e[j];
                begin = j ? e[j - 1] : 0;
            } else {
                uint32_t const* const e = overflow_ends.data() + overflow_first[h >> 41];
                end = e[j];
                begin = j ? e[j - 1] : 0;
            }
            start = base + begin;
            slots = end - begin;
            return true;
        }
        uint64_t Bytes() const { return n_lines * 64 + overflow_first.size() * 8 + overflow_ends.size() * 4; }
        static constexpr bool kTwoStep = false;
        static constexpr size_t kAhead = 8, kStartsAhead = 0;
    };

    // H (2026-10-15): the used cores in a hash table of 64-byte lines, 8 entries a line, each core's home line chosen by
    // a bijective mix of the core, so that the cores of related genomes (which cluster in core order and sank B) spread
    // over the lines. An entry is 64 bits: the mixed core's bits below the line index (its tag), how many lines past its
    // home it lies (0-31), its first slot, and its slot count (12 bits, never 0, so an empty entry is 0). One line a
    // lookup; a core whose home line filled up lies in a following line (linear probing), so a lookup reads on only while
    // the line it read is full (a line with an empty entry: no core of that home lies further).
    struct HashedLines {
        static constexpr uint64_t kCoreBits = 30, kCoreMask = (uint64_t{1} << kCoreBits) - 1;
        static constexpr unsigned kDispBits = 5, kCountShift = 52;
        uint64_t* lines = nullptr;
        unsigned bits = 0, tag_bits = 0, start_shift = 0;
        uint64_t n_home = 0, n_lines = 0, max_disp = 0, displaced = 0, used = 0, start_mask = 0, low_mask = 0, tag_mask = 0;
        std::array<uint64_t, 32> by_disp{};
        bool ok = true;
        bool branchless = false;  // Hb: the line's entries compared without a branch (HashedKeys.h's Locate)
        bool prefetch_next = false;  // Hp: Prefetch fetches the home line's next line too
        ~HashedLines() { std::free(lines); }
        // A bijection of [0, 2^30): odd multipliers mod 2^30 and xor-shifts.
        static uint64_t Mix(uint64_t core) {
            uint64_t x = core & kCoreMask;
            x = (x * 0x2545F491ull) & kCoreMask;
            x ^= x >> 15;
            x = (x * 0x1B873593ull) & kCoreMask;
            x ^= x >> 13;
            return x;
        }
        // `load`: the used cores per entry of the home lines at most (the line count a power of two).
        template<typename LocateFn>
        void Build(uint64_t cores, LocateFn&& locate, double load) {
            for (uint64_t core = 0; core < cores; core++) {
                uint64_t s, n;
                if (locate(core, s, n)) used++;
            }
            bits = 10;
            while (bits < 26 && static_cast<double>(used) > load * 8.0 * static_cast<double>(uint64_t{1} << bits)) bits++;
            tag_bits = static_cast<unsigned>(kCoreBits) - bits;
            tag_mask = (uint64_t{1} << tag_bits) - 1;
            start_shift = tag_bits + kDispBits;
            start_mask = (uint64_t{1} << (kCountShift - start_shift)) - 1;
            low_mask = (uint64_t{1} << start_shift) - 1;
            n_home = uint64_t{1} << bits;
            n_lines = n_home + 32;
            lines = static_cast<uint64_t*>(HugeAlloc(n_lines * 64));
            std::memset(lines, 0, n_lines * 64);
            for (uint64_t core = 0; core < cores; core++) {
                uint64_t s, n;
                if (!locate(core, s, n)) continue;
                if (s > start_mask || n >= 4096 || n == 0) {
                    ok = false;
                    continue;
                }
                uint64_t const p = Mix(core), home = p >> tag_bits, tag = p & tag_mask;
                bool placed = false;
                for (uint64_t d = 0; d < 32 && !placed; d++) {
                    uint64_t* const w = lines + (home + d) * 8;
                    for (int e = 0; e < 8; e++) {
                        if (w[e] != 0) continue;
                        w[e] = tag | d << tag_bits | s << start_shift | n << kCountShift;
                        placed = true;
                        max_disp = std::max(max_disp, d);
                        by_disp[d]++;
                        displaced += d > 0;
                        break;
                    }
                }
                if (!placed) ok = false;
            }
        }
        void Prefetch(uint64_t core) const {
            if (!prefetch_next) {
                __builtin_prefetch(lines + (Mix(core) >> tag_bits) * 8);
                return;
            }
            uint64_t const* const w = lines + (Mix(core) >> tag_bits) * 8;
            __builtin_prefetch(w);
            __builtin_prefetch(w + 8);  // Hp: the next line too (a core past its home, a full home line)
        }
        void PrefetchStart(uint64_t) const {}
        bool Locate(uint64_t core, uint64_t& start, uint64_t& slots) const {
            uint64_t const p = Mix(core), home = p >> tag_bits, tag = p & tag_mask;
            for (uint64_t d = 0;; d++) {
                uint64_t const* const w = lines + (home + d) * 8;
                uint64_t const want = tag | d << tag_bits;
                if (branchless) {
                    // The 8 entries without a branch: the core's entry (at most one) or-ed in, an empty one adds nothing.
                    uint64_t hit = 0, empty = 0;
                    for (int e = 0; e < 8; e++) {
                        uint64_t const x = w[e];
                        hit |= (x & low_mask) == want ? x : 0;
                        empty |= x == 0;
                    }
                    if (hit != 0) {
                        start = (hit >> start_shift) & start_mask;
                        slots = hit >> kCountShift;
                        return true;
                    }
                    if (empty != 0 || d >= max_disp) return false;
                    continue;
                }
                bool full = true;
                for (int e = 0; e < 8; e++) {
                    uint64_t const x = w[e];
                    if (x != 0 && (x & low_mask) == want) {
                        start = (x >> start_shift) & start_mask;
                        slots = x >> kCountShift;
                        return true;
                    }
                    full &= x != 0;
                }
                if (!full || d >= max_disp) return false;
            }
        }
        uint64_t Bytes() const { return n_lines * 64; }
        std::string Describe() const {
            char s[512];
            std::snprintf(s, sizeof(s), "2^%u lines (%.3f GB), %" PRIu64 " used cores, load %.3f, %" PRIu64 " (%.2f%%) past their home line (1: %" PRIu64 ", 2: %" PRIu64
                          ", 3+: %" PRIu64 "), max %" PRIu64 ", tag %u bits, start %u bits, ok %d",
                          bits, Bytes() / 1e9, used, static_cast<double>(used) / (8.0 * static_cast<double>(n_home)), displaced,
                          used ? 100.0 * static_cast<double>(displaced) / static_cast<double>(used) : 0.0, by_disp[1], by_disp[2],
                          displaced - by_disp[1] - by_disp[2], max_disp, tag_bits, kCountShift - start_shift, ok ? 1 : 0);
            return s;
        }
        static constexpr bool kTwoStep = false;
        static constexpr size_t kAhead = 8, kStartsAhead = 0;
    };

    // KEYBENCH_VALUES_GB: a values array of that many GB, whose lines each lookup then reads as FindSeeds reads a block
    // (the read's lookups' blocks fetched kValuesAhead lookups ahead, then each block's first lines read), so that the key
    // map's misses compete with the values' for the memory, as in a real run. Off (null) by default.
    uint8_t* g_values = nullptr;
    uint64_t g_values_lines = 0;
    constexpr size_t kValuesAhead = 16, kValueLines = 12;
    inline uint64_t ValueLine(uint64_t start) { return (start * 0x9E3779B97F4A7C15ull >> 20) % g_values_lines; }
    inline uint64_t ValueLines(uint64_t slots) { return std::min<uint64_t>(kValueLines, 1 + slots * 6 / 64); }
    thread_local std::vector<std::pair<uint64_t, uint64_t>> t_found;
    inline uint64_t ReadValues() {
        uint64_t sum = 0;
        size_t const n = t_found.size();
        auto fetch = [&](size_t i) {
            uint8_t const* const p = g_values + ValueLine(t_found[i].first) * 64;
            uint64_t const lines = ValueLines(t_found[i].second);
            for (uint64_t l = 0; l < lines && ValueLine(t_found[i].first) + l < g_values_lines; l++) __builtin_prefetch(p + 64 * l);
        };
        for (size_t i = 0; i < std::min(kValuesAhead, n); i++) fetch(i);
        for (size_t i = 0; i < n; i++) {
            if (i + kValuesAhead < n) fetch(i + kValuesAhead);
            uint64_t const line = ValueLine(t_found[i].first), lines = ValueLines(t_found[i].second);
            for (uint64_t l = 0; l < lines && line + l < g_values_lines; l++) sum += g_values[(line + l) * 64 + (l & 63)];
        }
        t_found.clear();
        return sum;
    }

    // A read's keys looked up as ChainAnchorFinder::LookupKmers looks them up: keys `ahead` k-mers ahead (twice as far
    // for a two-step layout), a two-step layout's start half as far. `prime_starts`: the first starts fetched before
    // the loop too (keys fetched a read earlier). Returns a checksum of the slots found.
    template<typename Keys>
    inline uint64_t LookupRead(Keys const& keys, uint64_t const* cores, size_t n, int prime_starts, uint64_t& found) {
        size_t const ahead = Keys::kAhead, starts_ahead = Keys::kStartsAhead;
        for (size_t j = 0; j < std::min(ahead, n); j++) keys.Prefetch(cores[j]);
        if (Keys::kTwoStep && (prime_starts == 1 || prime_starts == 2)) {
            for (size_t j = 0; j < std::min(prime_starts == 2 ? n : starts_ahead, n); j++) keys.PrefetchStart(cores[j]);
        }
        uint64_t sum = 0;
        for (size_t j = 0; j < n; j++) {
            if (j + ahead < n) keys.Prefetch(cores[j + ahead]);
            if (Keys::kTwoStep && j + starts_ahead < n) keys.PrefetchStart(cores[j + starts_ahead]);
            uint64_t start, slots;
            if (keys.Locate(cores[j], start, slots)) {
                sum += start * 31 + slots;
                found++;
                if (g_values) t_found.emplace_back(start, slots);
            }
        }
        if (g_values) sum += ReadValues();
        return sum;
    }

    struct Reads {
        std::vector<uint64_t> cores;   // every read's cores, in order
        std::vector<uint64_t> begin;   // read i: cores [begin[i], begin[i + 1])
    };

    struct Tally {
        double seconds = 0;
        uint64_t lookups = 0, found = 0, sum = 0;
    };

    // Spin about `ns` nanoseconds without touching memory (the rest of a read's work between its lookups and the next's).
    void Spin(uint64_t ns) {
        if (ns == 0) return;
        double const end = Now() + 1e-9 * static_cast<double>(ns);
        uint64_t x = 1;
        while (Now() < end) {
            for (int i = 0; i < 64; i++) x = x * 6364136223846793005ull + 1442695040888963407ull;
        }
        asm volatile("" : : "r"(x));
    }

    // One layout as the bench runs it. `prime` (short reads, whose keys were fetched a read earlier), for a two-step
    // layout: 0 nothing more (protal, 2026-10-14); 1 the first starts fetched before the lookup loop; 2 every start of
    // the read fetched before the loop; 3 the next read's starts fetched halfway through the gap (a read ahead, after
    // its keys).
    template<typename Keys>
    struct Run {
        std::string name;
        Keys const* keys;
        int prime;
    };

    // Each layout in turn per batch of reads (the order rotating), `rounds` times over the reads. Paired (short reads):
    // read r's lookups timed, then the next read's keys fetched (Classify's PrefetchKeys of the next pair) and `gap` ns of
    // other work, as protal seeds a read (pair) ahead. Else (long reads) each read's lookups streamed, nothing ahead.
    template<typename... Keys>
    void Bench(Reads const& reads, bool paired, int rounds, uint64_t gap, Run<Keys> const&... runs) {
        constexpr size_t kLayouts = sizeof...(Keys);
        std::array<Tally, kLayouts> tally{};
        size_t const n = reads.begin.size() - 2;
        constexpr size_t kBatch = 200;
        auto run = [&](auto const& r, size_t first, size_t last, Tally& t) {
            auto const& keys = *r.keys;
            // A read's keys fetched, then the gap (with prime 3 its starts fetched halfway).
            auto ahead = [&](size_t i) {
                for (size_t j = reads.begin[i]; j < reads.begin[i + 1]; j++) keys.Prefetch(reads.cores[j]);
                if (r.prime == 3) {
                    Spin(gap / 2);
                    for (size_t j = reads.begin[i]; j < reads.begin[i + 1]; j++) keys.PrefetchStart(reads.cores[j]);
                    Spin(gap - gap / 2);
                } else {
                    Spin(gap);
                }
            };
            if (paired) ahead(first);
            for (size_t i = first; i < last; i++) {
                uint64_t const* const c = reads.cores.data() + reads.begin[i];
                size_t const k = reads.begin[i + 1] - reads.begin[i];
                double const t0 = Now();
                t.sum += LookupRead(keys, c, k, paired ? r.prime : 0, t.found);
                t.seconds += Now() - t0;
                t.lookups += k;
                if (paired) ahead(i + 1);
            }
        };
        for (int round = 0; round < rounds; round++) {
            for (size_t b = 0, batch = 0; b < n; b += kBatch, batch++) {
                size_t const last = std::min(n, b + kBatch);
                size_t const rot = (batch + static_cast<size_t>(round)) % kLayouts;
                std::array<std::function<void()>, kLayouts> jobs;
                size_t idx = 0;
                ((jobs[idx] = [&, idx]() { run(runs, b, last, tally[idx]); }, idx++), ...);
                for (size_t li = 0; li < kLayouts; li++) jobs[(li + rot) % kLayouts]();
            }
        }
        std::array<std::string, kLayouts> names;
        size_t idx = 0;
        ((names[idx++] = runs.name), ...);
        std::printf("layout\tlookups\tfound\tns_per_lookup\tchecksum\n");
        for (size_t i = 0; i < kLayouts; i++) {
            std::printf("%s\t%" PRIu64 "\t%" PRIu64 "\t%.2f\t%" PRIu64 "\n", names[i].c_str(), tally[i].lookups, tally[i].found,
                        1e9 * tally[i].seconds / static_cast<double>(tally[i].lookups), tally[i].sum);
        }
    }

    // Every layout over the same reads.
    // 2026-10-15: D2 (protal's direct map), Cf (protal's compact map: every start fetched before the loop), H (hashed lines)
    // at two loads, and B; the other variants of 2026-10-14 left out.
    void BenchAll(Reads const& reads, bool paired, int rounds, uint64_t gap, DirectKeys const& d, CompactNow const& c,
                  HashedLines const& h, HashedLines const& h2, BucketLines const& b) {
        DirectKeys2 d2{ d.km };
        // KEYBENCH_H2_BRANCHLESS: the second table is H's, at its load, compared without a branch (Hb); else H at load 0.5.
        HashedLines const& hb = h2;
        if (paired) {
            Bench(reads, paired, rounds, gap, Run<DirectKeys2>{ "D2_direct_one_line", &d2, 0 },
                  Run<CompactNow>{ "Cf_compact_all_starts", &c, 2 }, Run<HashedLines>{ h.prefetch_next ? "Hp_hashed_branchless_next_line" : "H_hashed_lines", &h, 0 },
                  Run<HashedLines>{ hb.branchless ? "Hb_hashed_lines_branchless" : "H2_hashed_lines_lower_load", &h2, 0 },
                  Run<BucketLines>{ "B_bucket128", &b, 0 });
        } else {
            Bench(reads, paired, rounds, gap, Run<DirectKeys2>{ "D2_direct_one_line", &d2, 0 },
                  Run<CompactNow>{ "C_compact", &c, 0 }, Run<HashedLines>{ h.prefetch_next ? "Hp_hashed_branchless_next_line" : "H_hashed_lines", &h, 0 },
                  Run<HashedLines>{ hb.branchless ? "Hb_hashed_lines_branchless" : "H2_hashed_lines_lower_load", &h2, 0 },
                  Run<BucketLines>{ "B_bucket128", &b, 0 });
        }
    }

    // The values array (KEYBENCH_VALUES_GB), on huge pages and written, so that its pages exist before the timing.
    void MakeValues() {
        char const* const env = std::getenv("KEYBENCH_VALUES_GB");
        if (!env) return;
        uint64_t const bytes = static_cast<uint64_t>(std::stod(env) * 1e9) / 64 * 64;
        g_values = static_cast<uint8_t*>(HugeAlloc(bytes));
        std::memset(g_values, 1, bytes);
        g_values_lines = bytes / 64;
        std::fprintf(stderr, "values array: %.2f GB, each found lookup reads up to %zu of its lines, %zu lookups ahead\n", bytes / 1e9,
                     kValueLines, kValuesAhead);
    }

    // CompactKeys of the direct key map's used cores (as Seedmap::CompactKeymap builds it).
    void BuildCompact(CompactKeys& ck, DirectKeys const& d, uint64_t cores) {
        ck.Init(cores);
        uint64_t s, n, used = 0, end = 0;
        for (uint64_t core = 0; core < cores; core++) {
            if (d.Locate(core, s, n)) ck.MarkUsed(core);
        }
        ck.CountRanks();
        ck.AllocateStarts();
        for (uint64_t core = 0; core < cores; core++) {
            if (!d.Locate(core, s, n)) continue;
            if (used % CompactKeys::kGroup == 0) ck.SetBase(used, s);
            used++;
            end = s + n;
        }
        if (used % CompactKeys::kGroup == 0) ck.SetBase(used, end);
        uint64_t rank = 0;
        for (uint64_t core = 0; core < cores; core++) {
            if (d.Locate(core, s, n)) ck.SetOffset(rank++, s);
        }
        ck.SetOffset(rank, end);
    }

    template<typename A, typename B>
    bool SameSlots(A const& a, B const& b, uint64_t cores, uint64_t step) {
        for (uint64_t core = 0; core < cores; core += step) {
            uint64_t s1 = 0, n1 = 0, s2 = 0, n2 = 0;
            bool const f1 = a.Locate(core, s1, n1), f2 = b.Locate(core, s2, n2);
            if (f1 != f2 || (f1 && (s1 != s2 || n1 != n2))) {
                std::fprintf(stderr, "core %" PRIu64 ": %d %" PRIu64 " %" PRIu64 " vs %d %" PRIu64 " %" PRIu64 "\n", core, f1, s1, n1, f2, s2, n2);
                return false;
            }
        }
        return true;
    }

    std::optional<Seedmap::PackedLayout> LayoutOf(db::Bundle const& bundle, int threads, bool grid_bits) {
        auto gene_table = db::DbFile::InBundle(bundle, gene_table_file::kFileName);
        GenomeLoader genomes(db::DbFile::InBundle(bundle, "reference.fna"), db::DbFile::InBundle(bundle, "reference.map"), threads,
                             gene_table.Exists() ? std::optional<db::DbFile>(gene_table) : std::nullopt,
                             db::DbFile::InBundle(bundle, "unique_kmers.tsv").Size());
        auto const [max_taxid, max_gene, max_length] = genomes.IndexFieldMaxima();
        Seedmap probe;
        auto layout = Seedmap::PackedLayout::For(max_taxid, max_gene, max_length + probe.m_flex_k);
        if (grid_bits) layout.level_bits = GRID_LEVEL_BITS;
        return layout;
    }

    void LoadMap(Seedmap& map, db::Bundle const& bundle, Seedmap::PackedLayout const& layout, int level, bool compact, int threads) {
        map.SetLoadLevel(static_cast<uint8_t>(level));
        map.SetLoadAllowance(static_cast<uint64_t>(threads) * (uint64_t{24} << 20));
        map.SetCompactLoad(compact);
        map.Load(db::DbFile::InBundle(bundle, "index.prx"), threads, &layout);
    }

    Reads ReadCores(std::string const& path, size_t max_reads, bool full_mask) {
        ClosedSyncmer minimizer{ 15, 7, 2, full_mask };
        SimpleKmerHandler<ClosedSyncmer> kmer_handler{ 31, 15, minimizer };
        Seedmap probe;
        Reads reads;
        reads.begin.push_back(0);
        ThreadedGzIstream is(path.c_str());
        SeqReaderSE reader(is);
        FastxRecord r;
        KmerList k;
        while (reads.begin.size() <= max_reads && reader(r)) {
            k.clear();
            kmer_handler(std::string_view(r.sequence), k);
            for (auto const& [kmer, pos] : k) reads.cores.push_back(probe.MainKey(kmer));
            reads.begin.push_back(reads.cores.size());
        }
        reads.begin.push_back(reads.cores.size());  // a sentinel for the read-ahead
        return reads;
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: keybench load|keys|synth ...\n");
        return 1;
    }
    std::string const mode = argv[1];
    if (mode == "load") {
        if (argc < 6) return 1;
        int const level = std::stoi(argv[3]);
        bool const compact = std::stoi(argv[4]) != 0;
        int const threads = std::stoi(argv[5]);
        bool const grid_bits = argc > 6 ? std::stoi(argv[6]) != 0 : level > 0;
        std::string error;
        auto bundle = db::Bundle::Open(argv[2], error);
        if (!bundle) { std::fprintf(stderr, "%s\n", error.c_str()); return 1; }
        auto const layout = LayoutOf(*bundle, threads, grid_bits);
        Seedmap map;
        rusage ru0{}, ru1{};
        getrusage(RUSAGE_SELF, &ru0);
        double const t0 = Now();
        LoadMap(map, *bundle, *layout, level, compact, threads);
        double const seconds = Now() - t0;
        getrusage(RUSAGE_SELF, &ru1);
        auto tv = [](timeval const& a, timeval const& b) { return static_cast<double>(b.tv_sec - a.tv_sec) + 1e-6 * static_cast<double>(b.tv_usec - a.tv_usec); };
        std::printf("load\tlevel %d\tcompact %d\tthreads %d\twall %.3f s\tuser %.3f s\tsys %.3f s\tminor faults %ld\t%s\t%s\n", level, compact,
                    threads, seconds, tv(ru0.ru_utime, ru1.ru_utime), tv(ru0.ru_stime, ru1.ru_stime), ru1.ru_minflt - ru0.ru_minflt,
                    Status("VmHWM").c_str(), map.MemoryDescription().c_str());
        if (std::getenv("KEYBENCH_CHECKSUM")) {
            // Every core's slots and its cells in the file's layout (UnpackKey), hashed in core order.
            uint64_t h = 1469598103934665603ull, used = 0, cells = 0;
            auto mix = [&h](uint64_t x) { h = (h ^ x) * 1099511628211ull; };
            std::vector<uint64_t> out;
            for (uint64_t core = 0; core < (uint64_t{1} << 30); core++) {
                uint64_t start, slots;
                if (!map.Locate(core, start, slots)) continue;
                used++;
                mix(core);
                mix(start);
                mix(slots);
                out.resize(slots);
                map.UnpackKey(start, slots, out.data());
                for (uint64_t x : out) mix(x);
                cells += slots;
            }
            std::printf("checksum\tlevel %d\tcompact %d\tthreads %d\t%" PRIu64 " cores\t%" PRIu64 " slots\t%016" PRIx64 "\n", level, compact, threads, used,
                        cells, h);
        }
        return 0;
    }
    if (mode == "keys") {
        if (argc < 7) return 1;
        size_t const max_reads = std::stoul(argv[4]);
        bool const paired = std::stoi(argv[5]) != 0;
        int const rounds = std::stoi(argv[6]);
        int const level = argc > 7 ? std::stoi(argv[7]) : 0;
        uint64_t const gap = argc > 8 ? std::stoull(argv[8]) : 2000;
        std::string error;
        auto bundle = db::Bundle::Open(argv[2], error);
        if (!bundle) { std::fprintf(stderr, "%s\n", error.c_str()); return 1; }
        auto const layout = LayoutOf(*bundle, 2, level > 0);
        uint64_t const cores = uint64_t{1} << 30, blocks = cores >> 3;
        // The loaded direct key map, rebuilt in the bench's own array from Seedmap::Locate (Seedmap keeps its array
        // private), and the loaded index freed.
        HugeCells km(blocks * 12 + 4);
        bool full_mask = true;
        {
            auto direct = std::make_unique<Seedmap>();
            LoadMap(*direct, *bundle, *layout, level, false, 2);
            std::fprintf(stderr, "direct: %s\n", direct->MemoryDescription().c_str());
            full_mask = direct->UsesFullSyncmerMask();
            uint64_t base = 0;
            for (uint64_t b = 0; b < blocks; b++) {
                uint16_t* const cells = km.data() + b * 12;
                std::memcpy(cells, &base, 8);
                uint64_t local = 0;
                for (uint64_t j = 0; j < 8; j++) {
                    cells[4 + j] = static_cast<uint16_t>(local);
                    uint64_t start, slots;
                    if (direct->Locate(b * 8 + j, start, slots)) {
                        if (start != base + local) { std::fprintf(stderr, "direct map not contiguous at core %" PRIu64 "\n", b * 8 + j); return 2; }
                        local += slots;
                    }
                }
                base += local;
            }
            std::memcpy(km.data() + blocks * 12, &base, 8);
        }
        DirectKeys d{ km.data() };
        static CompactKeys ck;
        BuildCompact(ck, d, cores);
        CompactNow c{ &ck };
        auto located = [&](uint64_t core, uint64_t& s, uint64_t& n) { return d.Locate(core, s, n); };
        HashedLines h, h2;
        h.Build(cores, located, 0.75);
        if (std::getenv("KEYBENCH_H_PAIR")) h.branchless = h.prefetch_next = true;  // Hp: H branchless, the next line fetched too
        h2.branchless = std::getenv("KEYBENCH_H2_BRANCHLESS") != nullptr;  // Hb: H's load, compared without a branch
        h2.Build(cores, located, h2.branchless ? 0.75 : 0.5);
        BucketLines bl;
        bl.Build(cores, located);
        std::fprintf(stderr, "key maps: direct %.3f GB, compact %.3f GB (%" PRIu64 " used cores), hashed %.3f GB, hashed at a lower load %.3f GB, buckets %.3f GB (%zu overflow lines of %" PRIu64 ")\n",
                     km.size() * 2 / 1e9, ck.Bytes() / 1e9, ck.UsedCores(), h.Bytes() / 1e9, h2.Bytes() / 1e9, bl.Bytes() / 1e9, bl.overflow_first.size(), bl.n_lines);
        std::fprintf(stderr, "hashed lines: %s\nhashed lines, lower load: %s\n", h.Describe().c_str(), h2.Describe().c_str());
        std::fprintf(stderr, "bucket lines: %s\n", bl.Describe().c_str());
        bool const same = h.ok && h2.ok && SameSlots(d, c, cores, 1) && SameSlots(d, h, cores, 1) && SameSlots(d, h2, cores, 1) && SameSlots(d, bl, cores, 1);
        std::fprintf(stderr, "every core the same slots in every layout: %s\n", same ? "yes" : "NO");
        if (!same) return 3;
        Reads const reads = ReadCores(argv[3], max_reads, full_mask);
        std::fprintf(stderr, "%zu reads, %zu cores (%.1f a read), %s, gap %" PRIu64 " ns\n", reads.begin.size() - 2, reads.cores.size(),
                     static_cast<double>(reads.cores.size()) / static_cast<double>(reads.begin.size() - 2), paired ? "read ahead" : "streamed", gap);
        MakeValues();
        std::fprintf(stderr, "huge pages: %s\n", HugePages().c_str());
        BenchAll(reads, paired, rounds, gap, d, c, h, h2, bl);
        return 0;
    }
    if (mode == "synth") {
        if (argc < 6) return 1;
        double const density = std::stod(argv[2]);
        double const mean_slots = std::stod(argv[3]);
        size_t const lookups = std::stoul(argv[4]);
        int const rounds = std::stoi(argv[5]);
        size_t const read_kmers = argc > 6 ? std::stoul(argv[6]) : 27;
        bool const paired = argc > 7 ? std::stoi(argv[7]) != 0 : true;
        uint64_t const gap = argc > 8 ? std::stoull(argv[8]) : 2000;
        uint64_t const cores = uint64_t{1} << 30, blocks = cores >> 3;
        std::mt19937_64 rng(42);
        std::bernoulli_distribution used(density);
        std::geometric_distribution<uint64_t> extra(1.0 / mean_slots);
        HugeCells km(blocks * 12 + 4);
        uint64_t base = 0, n_used = 0;
        std::vector<uint32_t> used_cores;
        for (uint64_t b = 0; b < blocks; b++) {
            uint16_t* const cells = km.data() + b * 12;
            std::memcpy(cells, &base, 8);
            uint64_t local = 0;
            for (uint64_t j = 0; j < 8; j++) {
                cells[4 + j] = static_cast<uint16_t>(local);
                if (used(rng)) {
                    local += std::min<uint64_t>(1 + extra(rng), 3071);
                    n_used++;
                    used_cores.push_back(static_cast<uint32_t>(b * 8 + j));
                }
            }
            if (local > 60000) { std::fprintf(stderr, "block too full\n"); return 2; }
            base += local;
        }
        std::memcpy(km.data() + blocks * 12, &base, 8);
        DirectKeys d{ km.data() };
        static CompactKeys ck;
        BuildCompact(ck, d, cores);
        CompactNow c{ &ck };
        auto located = [&](uint64_t core, uint64_t& s, uint64_t& n) { return d.Locate(core, s, n); };
        HashedLines h, h2;
        h.Build(cores, located, 0.75);
        if (std::getenv("KEYBENCH_H_PAIR")) h.branchless = h.prefetch_next = true;  // Hp: H branchless, the next line fetched too
        h2.branchless = std::getenv("KEYBENCH_H2_BRANCHLESS") != nullptr;  // Hb: H's load, compared without a branch
        h2.Build(cores, located, h2.branchless ? 0.75 : 0.5);
        BucketLines bl;
        bl.Build(cores, located);
        std::fprintf(stderr, "synthetic: %" PRIu64 " used cores (%.4f), %" PRIu64 " slots; key maps: direct %.3f GB, compact %.3f GB, hashed %.3f GB, hashed at a lower load %.3f GB, buckets %.3f GB (%zu overflow lines of %" PRIu64 ")\n",
                     n_used, static_cast<double>(n_used) / static_cast<double>(cores), base, km.size() * 2 / 1e9, ck.Bytes() / 1e9, h.Bytes() / 1e9,
                     h2.Bytes() / 1e9, bl.Bytes() / 1e9, bl.overflow_first.size(), bl.n_lines);
        std::fprintf(stderr, "hashed lines: %s\nhashed lines, lower load: %s\n", h.Describe().c_str(), h2.Describe().c_str());
        std::fprintf(stderr, "bucket lines: %s\n", bl.Describe().c_str());
        bool const same = h.ok && h2.ok && SameSlots(d, c, cores, 7) && SameSlots(d, h, cores, 7) && SameSlots(d, h2, cores, 7) && SameSlots(d, bl, cores, 7);
        std::fprintf(stderr, "every 7th core the same slots in every layout: %s\n", same ? "yes" : "NO");
        if (!same) return 3;
        Reads reads;
        reads.begin.push_back(0);
        std::uniform_int_distribution<size_t> pick(0, used_cores.size() - 1);
        std::uniform_int_distribution<uint64_t> any(0, cores - 1);
        std::bernoulli_distribution hit(0.99);
        for (size_t i = 0; i < lookups; i++) {
            reads.cores.push_back(hit(rng) ? used_cores[pick(rng)] : any(rng));
            if (reads.cores.size() % read_kmers == 0) reads.begin.push_back(reads.cores.size());
        }
        if (reads.begin.back() != reads.cores.size()) reads.begin.push_back(reads.cores.size());
        reads.begin.push_back(reads.cores.size());
        std::vector<uint32_t>().swap(used_cores);
        MakeValues();
        std::fprintf(stderr, "huge pages: %s\n", HugePages().c_str());
        BenchAll(reads, paired, rounds, gap, d, c, h, h2, bl);
        return 0;
    }
    return 1;
}
