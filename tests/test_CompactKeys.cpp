// The compact key map (CompactKeys.h, --compact_keys): the key map of the used cores alone gives every core the slots the
// direct key map gives it (none to a core without values), so every lookup reads the same values.
#include <gtest/gtest.h>
#include <algorithm>
#include <bit>
#include <cstdint>
#include <random>
#include <string>
#include <vector>
#include "Hash/CompactKeys.h"
#include "Hash/Seedmap.h"
#include "TestUtil.h"

using protal::CompactKeys;
using protal::Seedmap;
using protal::ValueEntry;
using protal::test::ScratchDir;

// Cores used at random, the edges of presence lines (448 cores) and of start groups (64 ranks) among them, each of 1 to
// 3,071 slots: each used core gets its rank, its start and its slots; the others none. A start 2^24 or more past its
// group's base is refused.
TEST(CompactKeys, LocatesTheUsedCoresAlone) {
    uint64_t const cores = uint64_t{1} << 16;
    std::mt19937_64 rng(3);
    std::vector<uint64_t> slots(cores, 0);
    for (uint64_t c = 0; c < cores; c++) {
        if (rng() % 5 == 0) slots[c] = 1 + rng() % 3071;
    }
    for (uint64_t c : { uint64_t{0}, uint64_t{447}, uint64_t{448}, uint64_t{449}, uint64_t{895}, cores - 1 }) slots[c] = 1 + c % 7;
    CompactKeys compact;
    compact.Init(cores);
    for (uint64_t c = 0; c < cores; c++) {
        if (slots[c]) compact.MarkUsed(c);
    }
    uint64_t const used = compact.CountRanks();
    compact.AllocateStarts();
    std::vector<uint64_t> start(cores, 0);
    uint64_t next = 1000;  // the values need not start at slot 0
    for (uint64_t c = 0; c < cores; c++) {
        if (!slots[c]) continue;
        start[c] = next;
        next += slots[c];
    }
    uint64_t rank = 0;
    for (uint64_t c = 0; c < cores; c++) {
        if (!slots[c]) continue;
        if (rank % CompactKeys::kGroup == 0) compact.SetBase(rank, start[c]);
        rank++;
    }
    ASSERT_EQ(rank, used);
    if (used % CompactKeys::kGroup == 0) compact.SetBase(used, next);
    rank = 0;
    for (uint64_t c = 0; c < cores; c++) {
        if (!slots[c]) continue;
        ASSERT_TRUE(compact.SetOffset(rank, start[c]));
        EXPECT_EQ(compact.Rank(c), rank) << "core " << c;
        rank++;
    }
    ASSERT_TRUE(compact.SetOffset(used, next));
    for (uint64_t c = 0; c < cores; c++) {
        uint64_t s = 0, n = 0;
        bool const found = compact.Locate(c, s, n);
        ASSERT_EQ(found, slots[c] > 0) << "core " << c;
        EXPECT_EQ(compact.Used(c), slots[c] > 0) << "core " << c;
        if (!found) continue;
        EXPECT_EQ(s, start[c]) << "core " << c;
        EXPECT_EQ(n, slots[c]) << "core " << c;
    }
    EXPECT_EQ(compact.UsedCores(), used);
    EXPECT_LT(compact.Bytes(), cores * 3);  // the direct key map's 3 bytes a core

    CompactKeys wide;
    wide.Init(128);
    wide.MarkUsed(0);
    wide.MarkUsed(1);
    ASSERT_EQ(wide.CountRanks(), 2u);
    wide.AllocateStarts();
    wide.SetBase(0, 0);
    EXPECT_TRUE(wide.SetOffset(0, 0));
    EXPECT_TRUE(wide.SetOffset(1, CompactKeys::kMaxGroupSpan - 1));
    EXPECT_FALSE(wide.SetOffset(2, CompactKeys::kMaxGroupSpan));
}

namespace {
    // Maps of 10-base cores (Seedmap(size_t)), as test_PackedIndex's and test_GridLevels': the core from bit 16, the
    // flex part's halves below and above it.
    constexpr size_t kCoreBases = 10;
    constexpr unsigned kFlexHigh = 16 + 2 * kCoreBases;
    constexpr uint64_t kCoreMask = (uint64_t{1} << (2 * kCoreBases)) - 1;

    uint64_t Key(uint64_t core, uint64_t flex) {
        return (flex >> 16) << kFlexHigh | core << 16 | (flex & 0xffff);
    }
}

// The direct key map gives each core the slots its values take, the cores in order from slot 0 (a core of e >= 2 values
// e + ceil(e/2): its flex cells beside them), and nothing to a core without values, those of empty blocks included.
// Seedmap::Locate reads a core's end from its own block and the next block's first slot only for a block's last core
// (since 2026-10-14): whole blocks used, a block's last core alone, and the last core of the key map (whose end is the
// final control cell) are among the cores.
TEST(CompactKeys, TheDirectKeyMapGivesEveryCoreItsSlots) {
    std::mt19937_64 rng(5);
    Seedmap map(kCoreBases);
    uint64_t const cores = uint64_t{1} << (2 * kCoreBases);
    std::vector<uint64_t> values(cores, 0);
    std::vector<uint64_t> keys;
    std::vector<uint64_t> picked{ 0, 1, 2, 3, 4, 5, 6, 7, 15, 16, 23, 4095, cores - 9, cores - 1 };
    for (size_t k = 0; k < 30000; k++) picked.push_back(rng() & kCoreMask);
    for (size_t k = 0; k < picked.size(); k++) {
        size_t const n = k % 4 == 0 ? 1 : 2 + rng() % 30;
        for (size_t i = 0; i < n; i++) keys.push_back(Key(picked[k], rng() & 0xffffffffu));
        values[picked[k]] += n;
    }
    for (uint64_t key : keys) map.CountUpKey(map.MainKey(key));
    auto const layout = Seedmap::PackedLayout::For(1000, 100, 1000 + 16);
    map.BuildValuePointers(2, &layout, false);
    uint64_t next = 0;
    for (uint64_t c = 0; c < cores; c++) {
        uint64_t start = 0, slots = 0;
        bool const found = map.Locate(c, start, slots);
        ASSERT_EQ(found, values[c] > 0) << "core " << c;
        if (!found) continue;
        ASSERT_EQ(start, next) << "core " << c;
        ASSERT_EQ(slots, values[c] >= 2 ? values[c] + (values[c] + 1) / 2 : values[c]) << "core " << c;
        next += slots;
    }
    EXPECT_GT(next, keys.size());
}

// The compact key map's blocks are advised as huge pages, which madvise takes from a page-aligned start alone (until
// 2026-10-14 the 64-byte-aligned start was passed as it was, and the advice refused): a block of a huge page or more
// starts on a huge page, and the advice covers its whole pages; a smaller block gets none.
TEST(CompactKeys, BlocksAreAdvisedAsHugePagesFromAPageBoundary) {
    uint64_t const page = 4096, huge = CompactKeys::kHugePage;
    auto const [b0, e0] = CompactKeys::AdvisedRange(0x10040, 8 * huge, page);
    EXPECT_EQ(b0, 0x11000u);
    EXPECT_EQ(e0, (0x10040 + 8 * huge) & ~(page - 1));
    auto const [b1, e1] = CompactKeys::AdvisedRange(0x200000, 4 * huge, page);
    EXPECT_EQ(b1, 0x200000u);
    EXPECT_EQ(e1, 0x200000u + 4 * huge);
    auto const [b2, e2] = CompactKeys::AdvisedRange(0x10040, huge - 1, page);
    EXPECT_EQ(b2, e2);
    CompactKeys keys;
    keys.Init(CompactKeys::kCoresPerLine * (huge / 64));  // presence lines of 2 MB
    EXPECT_EQ(keys.PresenceBytes(), huge);
    EXPECT_EQ(reinterpret_cast<uintptr_t>(keys.PresenceData()) % huge, 0u);
}

// An index built as --build builds it (packed), then its key map compacted (Seedmap::CompactKeymap): every core of the
// 2^20 locates the same slots, so every lookup reads the same packed values; the direct key map is freed, the index line
// says so, and a second compaction does nothing.
TEST(CompactKeys, AnIndexLooksUpTheSameAfterItsKeyMapIsCompacted) {
    std::mt19937_64 rng(11);
    Seedmap map(kCoreBases);
    struct Value { uint64_t key, taxid, gene, pos; };
    std::vector<Value> values;
    for (size_t k = 0; k < 20000; k++) {
        uint64_t const core = rng() & kCoreMask;
        size_t const n = k % 3 == 0 ? 1 : k % 101 == 0 ? 200 + rng() % 400 : 2 + rng() % 20;
        for (size_t i = 0; i < n; i++) values.push_back({ Key(core, rng() & 0xffffffffu), 1 + rng() % 143614, 1 + rng() % 168, rng() % 12884 });
    }
    for (auto const& v : values) map.CountUpKey(map.MainKey(v.key));
    auto const layout = Seedmap::PackedLayout::For(143614, 168, 12883 + 16);
    map.BuildValuePointers(2, &layout, false);
    for (auto const& v : values) {
        ValueEntry entry;
        entry.Put(v.taxid, v.gene, v.pos);
        ASSERT_TRUE(map.PutOwned(v.key, entry.value));
    }
    uint64_t const cores = uint64_t{1} << (2 * kCoreBases);
    struct Slots { bool found; uint64_t start, slots; };
    std::vector<Slots> before(cores);
    uint64_t used = 0;
    for (uint64_t c = 0; c < cores; c++) {
        before[c].found = map.Locate(c, before[c].start, before[c].slots);
        used += before[c].found;
    }
    ASSERT_GT(used, 10000u);
    uint64_t const direct = map.KeymapBytes();
    ASSERT_TRUE(map.CompactKeymap(2));
    EXPECT_TRUE(map.HasCompactKeys());
    EXPECT_LT(map.KeymapBytes(), direct);
    EXPECT_NE(map.MemoryDescription().find("compact: " + std::to_string(used) + " cores with values of " + std::to_string(cores)),
              std::string::npos) << map.MemoryDescription();
    for (uint64_t c = 0; c < cores; c++) {
        uint64_t start = 0, slots = 0;
        bool const found = map.Locate(c, start, slots);
        ASSERT_EQ(found, before[c].found) << "core " << c;
        if (!found) continue;
        EXPECT_EQ(start, before[c].start) << "core " << c;
        EXPECT_EQ(slots, before[c].slots) << "core " << c;
    }
    EXPECT_FALSE(map.CompactKeymap(2)) << "compacted twice";
}

namespace {
    struct GridValue { uint64_t key, taxid, gene, pos, level; };

    // An index with grid levels (as test_GridLevels' FillGrid): keys of 1 to ~600 values, each with a level of 0 to 5.
    std::vector<GridValue> FillGridIndex(Seedmap& map, unsigned seed) {
        std::mt19937_64 rng(seed);
        std::vector<GridValue> values;
        for (size_t k = 0; k < 6000; k++) {
            uint64_t const core = rng() & kCoreMask;
            size_t const n = k % 3 == 0 ? 1 : k % 101 == 0 ? 200 + rng() % 400 : 2 + rng() % 20;
            for (size_t i = 0; i < n; i++) {
                uint64_t const level = std::min<uint64_t>(static_cast<uint64_t>(std::countr_zero(rng() | (uint64_t{1} << 5))), 5);
                values.push_back({ Key(core, rng() & 0xffffffffu), 1 + rng() % 143614, 1 + rng() % 168, rng() % 12884, level });
            }
        }
        map.SetGridLevels();
        for (auto const& v : values) map.CountUpKey(map.MainKey(v.key));
        auto layout = Seedmap::PackedLayout::For(143614, 168, 12883 + 16);
        layout.level_bits = protal::GRID_LEVEL_BITS;
        map.BuildValuePointers(2, &layout, false);
        for (auto const& v : values) {
            ValueEntry entry;
            entry.Put(v.taxid, v.gene, v.pos);
            entry.value |= v.level << protal::GRID_LEVEL_SHIFT;
            EXPECT_TRUE(map.PutOwned(v.key, entry.value));
        }
        return values;
    }

    // Every core's slots and entries (values and flex cells) of two loaded indexes are the same.
    void ExpectSameIndex(Seedmap const& a, Seedmap const& b, std::vector<GridValue> const& values) {
        uint64_t const cores = uint64_t{1} << (2 * kCoreBases);
        for (uint64_t c = 0; c < cores; c++) {
            uint64_t sa = 0, na = 0, sb = 0, nb = 0;
            bool const fa = a.Locate(c, sa, na), fb = b.Locate(c, sb, nb);
            ASSERT_EQ(fa, fb) << "core " << c;
            if (fa) {
                EXPECT_EQ(sa, sb) << "core " << c;
                EXPECT_EQ(na, nb) << "core " << c;
            }
        }
        for (auto const& v : values) {
            Seedmap::PackedBlock ba, bb;
            bool const fa = a.GetPacked(v.key, ba), fb = b.GetPacked(v.key, bb);
            ASSERT_EQ(fa, fb);
            if (!fa) continue;
            ASSERT_EQ(ba.size, bb.size);
            for (uint32_t i = 0; i < ba.size; i++) {
                EXPECT_EQ(a.EntryValue(ba, i), b.EntryValue(bb, i));
                if (ba.flex) EXPECT_EQ(Seedmap::FlexCell(ba, i), Seedmap::FlexCell(bb, i));
            }
        }
    }
}

// A load that builds the compact key map as it decodes (SetCompactLoad, --compact_keys), whole and sparse (grid levels
// 2 and 5), on one thread and three, holds the same index as the load of the direct key map: every core the same
// slots, every entry the same value and flex cell; the compact one never holds the direct key map, and its index line
// says so. (The whole direct load is the reference: its key map is decoded as before; the sparse loads' entries are
// checked against the file's in test_GridLevels.)
TEST(CompactKeys, ALoadBuildsTheCompactKeyMapAsItDecodes) {
    Seedmap map(kCoreBases);
    auto const values = FillGridIndex(map, 17);
    ScratchDir dir("compactload");
    protal::zstd::Params params;
    params.level = 3;
    params.window_log = 0;
    params.threads = 2;
    params.frame_size = uint64_t{1} << 16;  // many chunks: their edges among the cores
    uint64_t written = 0;
    size_t raw_chunks = 0;
    ASSERT_EQ(map.SaveCompressed(dir / "index.prx.zst", params, written, raw_chunks), "");
    auto layout = Seedmap::PackedLayout::For(143614, 168, 12883 + 16);
    layout.level_bits = protal::GRID_LEVEL_BITS;
    for (uint8_t level : { 0, 2, 5 }) {
        for (int threads : { 1, 3 }) {
            Seedmap direct(kCoreBases), compact(kCoreBases);
            direct.SetLoadLevel(level);
            compact.SetLoadLevel(level);
            compact.SetCompactLoad(true);
            direct.Load(dir / "index.prx.zst", threads, &layout);
            compact.Load(dir / "index.prx.zst", threads, &layout);
            EXPECT_FALSE(direct.HasCompactKeys());
            ASSERT_TRUE(compact.HasCompactKeys()) << "level " << int(level);
            EXPECT_LT(compact.KeymapBytes(), direct.KeymapBytes());
            EXPECT_NE(compact.MemoryDescription().find("compact: "), std::string::npos);
            EXPECT_EQ(compact.PackedEntries(), direct.PackedEntries());
            ExpectSameIndex(direct, compact, values);
        }
    }
}
