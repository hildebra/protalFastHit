// The packed layout of a query run's index (Seedmap::PackedLayout, Pack, GetPacked): every flex cell
// and entry of the 8-byte layout must come back, from Pack in several threads and from a load of the
// column format that packs chunk by chunk; a value outside the layout must stop the load.
#include <gtest/gtest.h>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <random>
#include <string>
#include <vector>
#include <unistd.h>
#include "Hash/Seedmap.h"
#include "Hash/KmerLookup.h"
#include "TestUtil.h"

using protal::Seedmap;
using protal::ValueEntry;
using protal::test::ScratchDir;
using protal::test::Slurp;

namespace {
    // The maps here hold cores of 10 bases (Seedmap(size_t)), a key map of 3 MB instead of 3 GB. A key is laid out as
    // a 15-base map's: the core from bit 16, the flex part's halves below it (bits 0-15) and above it (from kFlexHigh).
    constexpr size_t kCoreBases = 10;
    constexpr unsigned kFlexHigh = 16 + 2 * kCoreBases;
    constexpr uint64_t kCoreMask = (uint64_t{1} << (2 * kCoreBases)) - 1;

    uint64_t Key(uint64_t core, uint64_t flex) {
        return (flex >> 16) << kFlexHigh | core << 16 | (flex & 0xffff);
    }

    // Restores flex_scan's level when it goes out of scope, also when an assertion ends the test early.
    struct FlexScanSetting {
        protal::simd::Level const level = protal::flex_scan::Kernel().load();
        ~FlexScanSetting() { protal::flex_scan::Kernel().store(level); }
    };

    struct SmallValue { uint64_t key, taxid, gene, pos; };

    // Keys with 1 to ~300 values each; taxids, genes and positions in GTDB r226's ranges.
    std::vector<SmallValue> SmallValues(unsigned seed, size_t keys) {
        std::mt19937_64 rng(seed);
        std::vector<SmallValue> values;
        for (size_t k = 0; k < keys; k++) {
            uint64_t const core = rng() & kCoreMask;
            size_t const n = k % 7 == 0 ? 1 : k % 11 == 0 ? 2 : k % 23 == 0 ? 150 + rng() % 150 : 3 + rng() % 30;
            for (size_t i = 0; i < n; i++) {
                uint64_t const flex = rng() & 0xffffffffu;
                values.push_back({ Key(core, flex), 1 + rng() % 143614, 1 + rng() % 168, rng() % 12884 });
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

    // As Fill, but into the packed layout from the start (as --build since 2026-10-05), the same flags set through
    // the packed entries.
    void FillPacked(Seedmap& map, std::vector<SmallValue> const& values, Seedmap::PackedLayout const& layout) {
        for (auto const& v : values) map.CountUpKey(map.MainKey(v.key));
        map.BuildValuePointers(2, &layout);
        ASSERT_TRUE(map.IsPacked());
        for (auto const& v : values) {
            ValueEntry entry;
            entry.Put(v.taxid, v.gene, v.pos);
            ASSERT_TRUE(map.PutOwned(v.key, entry.value));
        }
        for (size_t i = 0; i < values.size(); i += 3) {
            Seedmap::PackedBlock block;
            ASSERT_TRUE(map.GetPacked(values[i].key, block));
            for (uint32_t e = 0; e < block.size; e++) {
                if (e % 2 == 1) map.ClearUniqueFlag(block, e);
                if (e % 3 == 0) map.SetUniqueDistanceTwoFlag(block, e);
            }
        }
    }

    // The column format in chunks of about frame_size bytes (many, so that chunk boundaries fall inside shared bytes).
    protal::zstd::Params ColumnParams(int threads = 2, uint64_t frame_size = uint64_t{1} << 16) {
        protal::zstd::Params params;
        params.level = 3;
        params.window_log = 0;
        params.threads = threads;
        params.frame_size = frame_size;
        return params;
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
    Seedmap map(kCoreBases);
    Fill(map, values);
    auto const stored = Record(map, values);
    auto const layout = Seedmap::PackedLayout::For(143614, 168, 12883 + 16);
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

// A lookup's seeds come out the same with the AVX2 scan (BestAvx2, TiesAvx2), the AVX-512 one (ScanAvx512) and the
// scalar one: the same entries in the same order with the same flags, and the same lookups dropped as too ubiquitous,
// for keys whose flex part matches a stored one exactly, in a few bases or not at all. Each level the CPU has.
TEST(PackedIndex, LookupsGiveTheSameSeedsAtEveryVectorLevel) {
    using protal::simd::Level;
    std::vector<Level> levels;
    for (Level level : { Level::avx2, Level::avx512 }) {
        if (protal::simd::Supported(level) == level) levels.push_back(level);
    }
    if (levels.empty()) GTEST_SKIP() << "no AVX2 here";
    if (levels.size() < 2) std::cout << "no AVX-512 here: AVX2 against scalar only" << std::endl;
    FlexScanSetting const restore;
    auto const values = SmallValues(9, 3000);
    Seedmap map(kCoreBases);
    Fill(map, values);
    map.Pack(Seedmap::PackedLayout::For(143614, 168, 12883 + 16), 4);
    std::mt19937_64 rng(4);
    // Base j of a key's flex part: at bit 2j (j < 8) or kFlexHigh + 2 (j - 8).
    auto change_base = [&](uint64_t key) {
        unsigned const j = static_cast<unsigned>(rng() % 16);
        unsigned const bit = j < 8 ? 2 * j : kFlexHigh + 2 * (j - 8);
        return key ^ ((1 + rng() % 3) << bit);
    };
    size_t compared = 0, seeds = 0, dropped = 0;
    for (size_t ubiquity : { size_t{256}, size_t{3} }) {
        protal::KmerLookupSM lookup(map, ubiquity);
        for (size_t k = 0; k < values.size(); k += 3) {
            for (int variant = 0; variant < 4; variant++) {
                uint64_t key = values[k].key;
                if (variant >= 1) key = change_base(key);
                if (variant >= 2) key = change_base(change_base(key));
                if (variant == 3) key = (key & (kCoreMask << 16)) | Key(0, rng() & 0xffffffffu);
                size_t kmer = key;
                std::vector<protal::LookupPointer> pointers;
                lookup.Get(pointers, kmer, static_cast<uint32_t>(rng() % 120));
                ASSERT_EQ(pointers.size(), 1u) << "the core is stored";
                protal::LookupList without;
                protal::flex_scan::Use(Level::scalar);
                bool const taken_without = lookup.GetFromLookup(without, pointers[0]);
                for (Level level : levels) {
                    protal::LookupList with;
                    ASSERT_EQ(protal::flex_scan::Use(level), level);
                    bool const taken_with = lookup.GetFromLookup(with, pointers[0]);
                    ASSERT_EQ(taken_with, taken_without) << protal::simd::Name(level) << ", key " << k << " variant " << variant;
                    ASSERT_EQ(with.size(), without.size()) << protal::simd::Name(level) << ", key " << k << " variant " << variant;
                    for (size_t i = 0; i < with.size(); i++) {
                        EXPECT_EQ(with[i].taxid, without[i].taxid);
                        EXPECT_EQ(with[i].geneid, without[i].geneid);
                        EXPECT_EQ(with[i].genepos, without[i].genepos);
                        EXPECT_EQ(with[i].readpos, without[i].readpos);
                        EXPECT_EQ(with[i].unique, without[i].unique);
                        EXPECT_EQ(with[i].unique_dist_two, without[i].unique_dist_two);
                    }
                }
                compared++;
                seeds += without.size();
                dropped += !taken_without;
            }
        }
    }
    EXPECT_GT(compared, 7000u);
    EXPECT_GT(seeds, compared);
    EXPECT_GT(dropped, 100u) << "lookups of more than 3 tied cells are dropped";
}

TEST(PackedIndex, TheColumnFormatPacksChunkByChunkToTheSameValues) {
    auto const values = SmallValues(6, 2000);
    Seedmap map(kCoreBases);
    Fill(map, values);
    auto const stored = Record(map, values);
    ScratchDir dir("packtest");
    uint64_t written = 0;
    size_t raw_chunks = 0;
    ASSERT_EQ(map.SaveCompressed(dir / "index.prx.zst", ColumnParams(), written, raw_chunks), "");
    auto const layout = Seedmap::PackedLayout::For(143614, 168, 12883 + 16);
    Seedmap loaded(kCoreBases);
    loaded.Load(dir / "index.prx.zst", 3, &layout);
    EXPECT_TRUE(loaded.IsPacked());
    EXPECT_EQ(loaded.PackedEntries(), values.size());
    ExpectPackedEquals(loaded, values, stored);
    EXPECT_TRUE(loaded.Layout() == layout);
}

// --build fills the packed layout directly: it holds what the 8-byte build holds once packed, the files it
// writes (the column format, and the raw index unpacked) are byte for byte the 8-byte build's, and the build's
// lookups (GetExact, GetSingleEntry) see the same entries.
TEST(PackedIndex, ABuildIntoThePackedLayoutWritesTheSameIndex) {
    auto const values = SmallValues(8, 3000);
    auto const layout = Seedmap::PackedLayout::For(143614, 168, 12883 + 16);
    Seedmap wide(kCoreBases);
    Fill(wide, values);
    auto const stored = Record(wide, values);
    Seedmap packed(kCoreBases);
    FillPacked(packed, values, layout);
    EXPECT_EQ(packed.PackedEntries(), values.size());
    ExpectPackedEquals(packed, values, stored);

    ScratchDir dir("packbuild");
    for (auto* map : { &wide, &packed }) {
        uint64_t written = 0;
        size_t raw_chunks = 0;
        std::string const name = map == &wide ? "wide" : "packed";
        ASSERT_EQ(map->SaveCompressed(dir / (name + ".prx.zst"), ColumnParams(3), written, raw_chunks), "") << name;
        EXPECT_EQ(raw_chunks, 0u);
        std::ofstream os(dir / (name + ".prx"), std::ios::binary);
        map->Save(os);
    }
    EXPECT_EQ(Slurp(dir / "packed.prx.zst"), Slurp(dir / "wide.prx.zst"));
    std::string const raw = Slurp(dir / "packed.prx");
    EXPECT_EQ(raw, Slurp(dir / "wide.prx"));
    EXPECT_EQ(raw.size(), wide.SerializedSize());

    protal::KmerLookupSM lookup(packed, 256);
    size_t exact_found = 0, singles = 0;
    for (size_t k = 0; k < values.size(); k += 5) {
        size_t kmer = values[k].key;
        Seedmap::PackedBlock block;
        std::vector<uint32_t> exact;
        size_t const compared = lookup.GetExact(kmer, block, exact);
        auto const& s = stored[k];
        if (s.flex.empty()) {
            EXPECT_EQ(compared, 0u);
            EXPECT_TRUE(exact.empty());
            ASSERT_TRUE(lookup.GetSingleEntry(kmer, block)) << "key " << k;
            EXPECT_EQ(lookup.Entry(block, 0).value, s.entries[0]);
            singles++;
            continue;
        }
        EXPECT_FALSE(lookup.GetSingleEntry(kmer, block));
        EXPECT_EQ(compared, s.entries.size());
        std::vector<uint32_t> expected;
        for (uint32_t i = 0; i < s.flex.size(); i++) {
            if (s.flex[i] == static_cast<uint32_t>(packed.FlexKey(kmer))) expected.push_back(i);
        }
        EXPECT_EQ(exact, expected) << "key " << k;
        exact_found += exact.size();
    }
    EXPECT_GT(exact_found, 100u);
    EXPECT_GT(singles, 10u);
}

TEST(PackedIndex, ABuildStopsAtAValueWiderThanItsLayout) {
    auto const values = SmallValues(9, 50);
    auto const narrow = Seedmap::PackedLayout::For(1000, 168, 12883 + 16);  // taxids go to 143614
    Seedmap map(kCoreBases);
    for (auto const& v : values) map.CountUpKey(map.MainKey(v.key));
    map.BuildValuePointers(1, &narrow);
    EXPECT_EXIT({
        for (auto const& v : values) {
            ValueEntry entry;
            entry.Put(v.taxid, v.gene, v.pos);
            map.PutOwned(v.key, entry.value);
        }
    }, testing::ExitedWithCode(8), "wider than the reference");
}

TEST(PackedIndex, AValueOutsideTheLayoutStopsTheLoad) {
    auto const values = SmallValues(7, 200);
    Seedmap map(kCoreBases);
    Fill(map, values);
    auto const narrow = Seedmap::PackedLayout::For(1000, 168, 12883 + 16);  // taxids go to 143614
    EXPECT_EXIT(map.Pack(narrow, 2), testing::ExitedWithCode(8), "outside the reference");
}

// The column format is packed key by key as its chunks are decoded (index_codec::ValueSink): a value outside the
// layout stops that load too, with the same message.
TEST(PackedIndex, AValueOutsideTheLayoutStopsTheColumnLoad) {
    auto const values = SmallValues(7, 200);
    Seedmap map(kCoreBases);
    Fill(map, values);
    ScratchDir dir("packtest_narrow");
    std::string const path = dir / "index.prx.zst";
    uint64_t written = 0;
    size_t raw_chunks = 0;
    ASSERT_EQ(map.SaveCompressed(path, ColumnParams(1), written, raw_chunks), "");
    auto const narrow = Seedmap::PackedLayout::For(1000, 168, 12883 + 16);
    EXPECT_EXIT({ Seedmap loaded(kCoreBases); loaded.Load(path, 2, &narrow); }, testing::ExitedWithCode(8), "outside the reference");
}

// A column-format index whose container describes another layout than its header (a corrupt file) stops the load,
// packed or not, before anything is decoded.
TEST(PackedIndex, AColumnIndexWhoseLayoutIsNotItsHeadersStopsTheLoad) {
    Seedmap empty(kCoreBases);  // the header: 4^10 keys, no values
    auto layout = empty.CodecLayout();
    layout.blocks /= 2;         // the container: half the blocks
    std::vector<uint16_t> const keymap(layout.KeymapCells(), 0);
    std::vector<uint64_t> const values;
    ScratchDir dir("packtest_layout");
    std::string const path = dir / "index.prx.zst";
    std::string error;
    ASSERT_TRUE(protal::index_codec::Write(path, empty.HeaderBytes(), layout, keymap.data(), values.data(), ColumnParams(1), error))
        << error;
    auto const pack = Seedmap::PackedLayout::For(1000, 168, 12883 + 16);
    EXPECT_EXIT({ Seedmap loaded(kCoreBases); loaded.Load(path); }, testing::ExitedWithCode(8), "does not match its header");
    EXPECT_EXIT({ Seedmap loaded(kCoreBases); loaded.Load(path, 2, &pack); }, testing::ExitedWithCode(8), "does not match its header");
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
    EXPECT_EQ(l.genes, 169u);
    EXPECT_EQ(l.taxon_gene_bits, 25u);  // 143614 * 169 + 168 < 2^25; apart, 18 + 8
    EXPECT_EQ(l.pos_bits, 14u);
    EXPECT_EQ(Seedmap::PackedLayout::For(0, 0, 0).EntryBits(), 4u);
    EXPECT_EQ(Seedmap::PackedLayout::For(1000, 255, 0).taxon_gene_bits, 18u);  // 256 genes: as taxid and gene apart
    EXPECT_EQ(Seedmap::PackedLayout::For(~0ull, ~0ull, ~0ull).EntryBits(), 62u);  // the file's own widths
    // A key of S slots holds e = S - ceil(S/3) entries and (S >= 2) e flex cells: its region, from bit
    // floor(a * F / 32) to floor((a + S) * F / 32) for its first slot a and F = SlotBits32(), must hold them whatever
    // a is (the pattern repeats every 32 slots), and F is the least that does.
    auto fits = [](uint64_t F, unsigned w) {
        for (uint64_t S = 1; S < 4000; S++) {
            uint64_t const e = S - (S >= 2 ? (S + 2) / 3 : 0);
            uint64_t const need = (S >= 2 ? 32 * e : 0) + w * e;
            for (uint64_t a = 0; a < 32; a++) {
                if (((a + S) * F >> 5) - (a * F >> 5) < need) return false;
            }
        }
        return true;
    };
    for (unsigned w : { 4u, 5u, 30u, 41u, 42u, 50u, 62u }) {
        Seedmap::PackedLayout layout{ 1, w - 3, 1 };
        ASSERT_EQ(layout.EntryBits(), w);
        EXPECT_TRUE(fits(layout.SlotBits32(), w)) << "W " << w;
        EXPECT_FALSE(fits(layout.SlotBits32() - 1, w)) << "W " << w;
    }
}

// An entry's taxid and gene are one number, taxid * genes + gene, split again by a multiplication
// (Seedmap::DivideByGenes): every value comes back for 2 genes (also the layout of a reference whose only gene
// id is 0), a power of two, r226's 169 and the file's 2^20, with the largest and smallest taxids and genes; and the
// division is exact over the whole range.
TEST(PackedIndex, TaxidAndGeneAsOneNumberComeBackForEveryGeneCount) {
    struct Case { uint64_t max_taxid, max_gene; };
    for (Case const c : { Case{ 143614, 168 }, Case{ 5, 0 }, Case{ 70000, 1 }, Case{ 1000, 255 }, Case{ 999, 999 },
                          Case{ (1u << 20) - 1, (1u << 20) - 1 } }) {
        std::mt19937_64 rng(c.max_taxid ^ c.max_gene);
        std::vector<SmallValue> values;
        for (size_t k = 0; k < 400; k++) {
            uint64_t const core = rng() & kCoreMask;
            size_t const n = k % 5 == 0 ? 1 : 2 + rng() % 40;
            for (size_t i = 0; i < n; i++) {
                uint64_t const flex = rng() & 0xffffffffu;
                uint64_t const taxid = i % 4 == 0 ? c.max_taxid : i % 4 == 1 ? 0 : rng() % (c.max_taxid + 1);
                uint64_t const gene = i % 3 == 0 ? c.max_gene : i % 3 == 1 ? 0 : rng() % (c.max_gene + 1);
                values.push_back({ Key(core, flex), taxid, gene, rng() % 12884 });
            }
        }
        Seedmap map(kCoreBases);
        Fill(map, values);
        auto const stored = Record(map, values);
        auto const layout = Seedmap::PackedLayout::For(c.max_taxid, c.max_gene, 12883);
        map.Pack(layout, 2);
        ExpectPackedEquals(map, values, stored);

        EXPECT_EQ(layout.genes, std::max<uint64_t>(c.max_gene, 1) + 1);
        uint64_t const genes = layout.genes, end = uint64_t{1} << layout.taxon_gene_bits;
        size_t wrong = 0;
        auto check = [&](uint64_t x) { wrong += map.DivideByGenes(x) != x / genes; };
        for (uint64_t x = 0; x < std::min<uint64_t>(end, 1u << 20); x++) check(x);
        for (uint64_t x = end - std::min<uint64_t>(end, 1u << 20); x < end; x++) check(x);
        for (int i = 0; i < 1000000; i++) check(rng() % end);
        EXPECT_EQ(wrong, 0u) << "genes " << genes;
    }
}
