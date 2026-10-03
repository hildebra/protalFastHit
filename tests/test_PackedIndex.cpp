// The packed layout of a query run's index (Seedmap::PackedLayout, Pack, GetPacked): every flex cell
// and entry of the 8-byte layout must come back, from Pack in several threads and from a load of the
// column format that packs chunk by chunk; a value outside the layout must stop the load.
#include <gtest/gtest.h>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <random>
#include <string>
#include <vector>
#include <unistd.h>
#include "Hash/Seedmap.h"
#include "Hash/KmerLookup.h"

using protal::Seedmap;
using protal::ValueEntry;

namespace {
    struct SmallValue { uint64_t key, taxid, gene, pos; };

    // Keys as the build sees them (62 bits: the 15-mer core in bits 16-45, the flex part around it),
    // with 1 to ~300 values each; taxids, genes and positions in GTDB r226's ranges.
    std::vector<SmallValue> SmallValues(unsigned seed, size_t keys) {
        std::mt19937_64 rng(seed);
        std::vector<SmallValue> values;
        for (size_t k = 0; k < keys; k++) {
            uint64_t const core = rng() & ((uint64_t{1} << 30) - 1);
            size_t const n = k % 7 == 0 ? 1 : k % 11 == 0 ? 2 : k % 23 == 0 ? 150 + rng() % 150 : 3 + rng() % 30;
            for (size_t i = 0; i < n; i++) {
                uint64_t const flex = rng() & 0xffffffffu;
                uint64_t const key = (flex >> 16) << 46 | core << 16 | (flex & 0xffff);
                values.push_back({ key, 1 + rng() % 143614, 1 + rng() % 168, rng() % 12884 });
            }
        }
        return values;
    }

    // Builds the index as --build does (count, lay out, place), and sets some flags.
    void Fill(Seedmap& map, std::vector<SmallValue> const& values) {
        for (auto const& v : values) map.CountUpKey(map.MainKey(v.key));
        map.BuildValuePointers(2);
        for (auto const& v : values) {
            ValueEntry entry;
            entry.Put(v.taxid, v.gene, v.pos);
            ASSERT_TRUE(map.PutOwned(v.key, entry.value));
        }
        for (size_t i = 0; i < values.size(); i += 3) {
            ValueEntry *start, *end;
            uint32_t *fb, *fe;
            map.Get(values[i].key, start, end, fb, fe);
            ASSERT_NE(start, nullptr);
            for (auto* e = start; e < end; e++) {
                if ((e - start) % 2 == 1) e->SetFlagNonUnique();
                if ((e - start) % 3 == 0) e->SetFlagUniqueDistanceMinTwo();
            }
        }
    }

    struct Stored { std::vector<uint32_t> flex; std::vector<uint64_t> entries; };

    // What the 8-byte layout holds for each value's key.
    std::vector<Stored> Record(Seedmap& map, std::vector<SmallValue> const& values) {
        std::vector<Stored> stored;
        for (auto const& v : values) {
            ValueEntry *start, *end;
            uint32_t *fb = nullptr, *fe = nullptr;
            map.Get(v.key, start, end, fb, fe);
            Stored s;
            for (auto* e = start; e < end; e++) s.entries.push_back(e->value);
            if (fb) for (size_t i = 0; i < s.entries.size(); i++) s.flex.push_back(fb[i]);
            stored.push_back(std::move(s));
        }
        return stored;
    }

    void ExpectPackedEquals(Seedmap const& map, std::vector<SmallValue> const& values, std::vector<Stored> const& stored) {
        size_t checked = 0;
        for (size_t k = 0; k < values.size(); k++) {
            Seedmap::PackedBlock block;
            ASSERT_TRUE(map.GetPacked(values[k].key, block)) << "key " << k;
            ASSERT_EQ(block.size, stored[k].entries.size()) << "key " << k;
            EXPECT_EQ(block.flex != nullptr, !stored[k].flex.empty()) << "key " << k;
            for (uint32_t i = 0; i < block.size; i++) {
                EXPECT_EQ(map.EntryValue(block, i), stored[k].entries[i]) << "key " << k << " entry " << i;
                if (block.flex) EXPECT_EQ(Seedmap::FlexCell(block, i), stored[k].flex[i]) << "key " << k << " cell " << i;
                checked++;
            }
        }
        EXPECT_GT(checked, values.size());
    }
}

TEST(PackedIndex, HoldsEveryFlexCellAndEntryOfTheStoredLayout) {
    auto const values = SmallValues(5, 3000);
    Seedmap map;
    Fill(map, values);
    auto const stored = Record(map, values);
    auto const layout = Seedmap::PackedLayout::For(143614, 168, 12883 + 16);
    EXPECT_EQ(layout.EntryBits(), 42u);
    EXPECT_EQ(layout.SlotBits(), 50u);
    map.Pack(layout, 4);
    EXPECT_TRUE(map.IsPacked());
    EXPECT_EQ(map.PackedEntries(), values.size());
    ExpectPackedEquals(map, values, stored);

    // The lookups a read makes see the same values: GetFromLookup against the stored cells.
    protal::KmerLookupSM lookup(map, 256);
    size_t compared = 0;
    for (size_t k = 0; k < values.size(); k += 17) {
        protal::LookupList results;
        size_t kmer = values[k].key;
        lookup.Get(results, kmer, 7);
        auto const& s = stored[k];
        std::vector<uint64_t> expected;
        if (s.flex.empty()) {
            expected = s.entries;
        } else {
            uint32_t const query = static_cast<uint32_t>(map.FlexKey(kmer));
            uint32_t max = 0;
            for (auto cell : s.flex) max = std::max(max, Seedmap::Similarity(cell, query));
            for (size_t i = 0; i < s.flex.size(); i++) {
                if (Seedmap::Similarity(s.flex[i], query) == max) expected.push_back(s.entries[i]);
            }
        }
        ASSERT_EQ(results.size(), expected.size()) << "key " << k;
        for (size_t i = 0; i < results.size(); i++) {
            ValueEntry e;
            e.value = expected[i];
            auto const [taxid, gene, pos] = e.Get();
            EXPECT_EQ(results[i].taxid, taxid);
            EXPECT_EQ(results[i].geneid, gene);
            EXPECT_EQ(results[i].genepos, pos);
            compared++;
        }
    }
    EXPECT_GT(compared, 100u);
}

TEST(PackedIndex, TheColumnFormatPacksChunkByChunkToTheSameValues) {
    auto const values = SmallValues(6, 2000);
    Seedmap map;
    Fill(map, values);
    auto const stored = Record(map, values);
    auto const dir = std::filesystem::temp_directory_path() / ("protal_packtest_" + std::to_string(::getpid()));
    std::filesystem::create_directories(dir);
    std::string const path = (dir / "index.prx.zst").string();
    uint64_t written = 0;
    size_t raw_chunks = 0;
    protal::zstd::Params params;
    params.level = 3;
    params.window_log = 0;
    params.threads = 2;
    params.frame_size = uint64_t{1} << 20;  // many chunks, so that chunk boundaries fall inside shared bytes
    ASSERT_EQ(map.SaveCompressed(path, params, written, raw_chunks), "");
    auto const layout = Seedmap::PackedLayout::For(143614, 168, 12883 + 16);
    Seedmap loaded;
    loaded.Load(path, 3, &layout);
    EXPECT_TRUE(loaded.IsPacked());
    EXPECT_EQ(loaded.PackedEntries(), values.size());
    ExpectPackedEquals(loaded, values, stored);
    EXPECT_TRUE(loaded.Layout() == layout);
    std::filesystem::remove_all(dir);
}

TEST(PackedIndex, AValueOutsideTheLayoutStopsTheLoad) {
    auto const values = SmallValues(7, 200);
    Seedmap map;
    Fill(map, values);
    auto const narrow = Seedmap::PackedLayout::For(1000, 168, 12883 + 16);  // taxids go to 143614
    EXPECT_EXIT(map.Pack(narrow, 2), testing::ExitedWithCode(8), "outside the reference");
}

TEST(PackedIndex, PutBitsSharesBytesBetweenRanges) {
    // Two ranges meet inside byte 6 (bit 53): each writes its fields with the other's bytes unsafe.
    std::vector<uint8_t> bits(32, 0);
    Seedmap::PutBits(bits.data(), 11, 0x2aaaaaaaaaull, 42, 0, 6);   // bits 11-52: the first range (bytes 0-5 safe)
    Seedmap::PutBits(bits.data(), 53, 0x155555555ull, 33, 7, 32);   // bits 53-85: the second (bytes 7-31 safe)
    Seedmap::PutBits(bits.data(), 86, 0xfffff, 20, 7, 32);
    auto read = [&](uint64_t bit, unsigned w) {
        __uint128_t x;
        std::memcpy(&x, bits.data() + (bit >> 3), 16);
        return static_cast<uint64_t>(x >> (bit & 7)) & ((uint64_t{1} << w) - 1);
    };
    EXPECT_EQ(read(0, 11), 0u);
    EXPECT_EQ(read(11, 42), 0x2aaaaaaaaaull);
    EXPECT_EQ(read(53, 33), 0x155555555ull);
    EXPECT_EQ(read(86, 20), 0xfffffu);
    EXPECT_EQ(read(106, 20), 0u);
}

TEST(PackedIndex, LayoutWidthsFollowTheReference) {
    auto const l = Seedmap::PackedLayout::For(143614, 168, 12883);
    EXPECT_EQ(l.taxid_bits, 18u);
    EXPECT_EQ(l.gene_bits, 8u);
    EXPECT_EQ(l.pos_bits, 14u);
    EXPECT_EQ(Seedmap::PackedLayout::For(0, 0, 0).EntryBits(), 5u);
    EXPECT_EQ(Seedmap::PackedLayout::For(~0ull, ~0ull, ~0ull).EntryBits(), 62u);  // the file's own widths
    // A key of S slots holds e = S - ceil(S/3) entries and (S >= 2) e flex cells: they must fit S * SlotBits().
    for (unsigned w : { 5u, 30u, 42u, 50u, 62u }) {
        Seedmap::PackedLayout layout{ w - 2, 0, 0 };
        for (uint64_t S = 1; S < 4000; S++) {
            uint64_t const e = S - (S >= 2 ? (S + 2) / 3 : 0);
            EXPECT_LE((S >= 2 ? 32 * e : 0) + layout.EntryBits() * e, S * layout.SlotBits()) << "W " << w << " S " << S;
        }
    }
}
