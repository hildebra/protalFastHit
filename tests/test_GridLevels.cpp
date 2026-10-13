// The gene grid of a sparse index (GridLevels.h, --sparse_first_pick grid): the tournament that gives each of a gene's
// cores its level, the scores it compares (strain alleles, column weights, ubiquity), the levels kept by the index's
// values (packed and in the file), a sparse load of one level, and lookups on the whole index that count only the
// entries of a level seeding exactly as lookups on a sparse index of that level.
#include <gtest/gtest.h>
#include <algorithm>
#include <cstdint>
#include <map>
#include <random>
#include <string>
#include <vector>
#include "Hash/Seedmap.h"
#include "Hash/KmerLookup.h"
#include "SequenceUtils/GridLevels.h"
#include "TestUtil.h"

using protal::Seedmap;
using protal::ValueEntry;
using protal::test::ScratchDir;
namespace grid = protal::grid_levels;
namespace sa = protal::strain_alleles;

namespace {
    std::vector<uint8_t> LevelsOf(std::vector<double> const& scores, std::vector<uint8_t> const& eligible) {
        std::vector<uint8_t> levels;
        std::vector<int64_t> winners;
        grid::Levels(scores, eligible, levels, winners);
        return levels;
    }

    // An allele over [begin, end) with these edits (sorted).
    sa::Table::View View(uint16_t begin, uint16_t end, std::vector<sa::Edit> const& edits) {
        return { begin, end, std::span<sa::Edit const>(edits) };
    }
}

// Each window of 2^l consecutive cores keeps exactly one core at level l or more, for every l and any scores (a window
// the gene ends inside too, and a gene shorter than the window as one): the cores a run gated at C = 2^l counts are 1/C
// of the gene's, about every C-th, and those of C are among those of C/2.
TEST(GridLevels, EachWindowKeepsOneCoreAtItsLevel) {
    std::mt19937_64 rng(5);
    for (size_t n : { 1u, 2u, 3u, 7u, 64u, 100u, 257u }) {
        std::vector<double> scores(n);
        for (auto& s : scores) s = static_cast<double>(rng() % 1000) / 1000.0;
        auto const levels = LevelsOf(scores, std::vector<uint8_t>(n, 1));
        for (unsigned l = 1; l < grid::kLevels; l++) {
            size_t const size = size_t{1} << l;
            for (size_t w = 0; w * size < n; w++) {
                size_t kept = 0;
                for (size_t i = w * size; i < std::min(n, (w + 1) * size); i++) kept += levels[i] >= l;
                EXPECT_EQ(kept, 1u) << "n " << n << " level " << l << " window " << w;
            }
        }
    }
}

// The better score wins a round; between equal ones, the core nearer the window's middle, then the earlier one. With
// equal scores the kept cores are spread out: of 8 cores, C 2 keeps 0, 2, 4, 6, C 4 keeps 2 and 6 (the middles), C 8 and
// beyond keep 2. A core the index leaves out never wins.
TEST(GridLevels, TheBetterCoreWinsAndTiesGoToTheMiddle) {
    EXPECT_EQ(LevelsOf({ 0, 0, 1, 0 }, { 1, 1, 1, 1 }), (std::vector<uint8_t>{ 1, 0, 5, 0 }));
    EXPECT_EQ(LevelsOf(std::vector<double>(8, 0.5), std::vector<uint8_t>(8, 1)), (std::vector<uint8_t>{ 1, 0, 5, 0, 1, 0, 2, 0 }));
    EXPECT_EQ(LevelsOf({ 0.9, 0.1, 0.2, 0.3 }, { 0, 1, 1, 1 }), (std::vector<uint8_t>{ 0, 1, 0, 5 }));
    EXPECT_EQ(LevelsOf({ 0.5, 0.5 }, { 0, 0 }), (std::vector<uint8_t>{ 0, 0 }));
    EXPECT_EQ(LevelsOf({ 0.5 }, { 1 }), (std::vector<uint8_t>{ 5 })) << "a gene of one core keeps it at every C";
}

// --grid_choice phase: every C-th core at each C (exactly C apart, the windows' one core each), the phase per level the
// one whose kept cores score more in sum, nested in the level below's; a gene keeps a core at every C. even (Levels with
// equal scores): every C-th core whatever the scores.
TEST(GridLevels, ThePhaseChoiceKeepsEveryCthCoreWhereTheScoresSumHighest) {
    auto phase = [](std::vector<double> const& scores, std::vector<uint8_t> const& eligible) {
        std::vector<uint8_t> levels;
        grid::PhaseLevels(scores, eligible, levels);
        return levels;
    };
    // Odd cores score more: phase 1 at C 2; at C 4 core 3 and 7 (0.9 each) beat 1 and 5 (0.6); at C 8, 3 beats 7 (tie:
    // the first); a gene of 8 cores keeps 3 from C 8 on.
    std::vector<double> const scores{ 0.1, 0.6, 0.1, 0.9, 0.1, 0.6, 0.1, 0.9 };
    EXPECT_EQ(phase(scores, std::vector<uint8_t>(8, 1)), (std::vector<uint8_t>{ 0, 1, 0, 5, 0, 1, 0, 2 }));
    // A core the index leaves out counts nothing: with 3 and 7 out, C 4 keeps 1 and 5.
    EXPECT_EQ(phase(scores, { 1, 1, 1, 0, 1, 1, 1, 0 }), (std::vector<uint8_t>{ 0, 5, 0, 1, 0, 2, 0, 1 }));
    EXPECT_EQ(phase({ 0.5 }, { 1 }), (std::vector<uint8_t>{ 5 })) << "a gene of one core keeps it at every C";
    EXPECT_EQ(phase({ 0.1, 0.9, 0.1 }, { 1, 1, 1 }), (std::vector<uint8_t>{ 0, 5, 0 }));
    std::mt19937_64 rng(9);
    for (size_t n : { 2u, 5u, 64u, 100u, 257u }) {
        std::vector<double> random(n);
        for (auto& s : random) s = static_cast<double>(rng() % 1000) / 1000.0;
        auto const levels = phase(random, std::vector<uint8_t>(n, 1));
        for (unsigned l = 1; l < grid::kLevels; l++) {
            std::vector<size_t> kept;
            for (size_t i = 0; i < n; i++) if (levels[i] >= l) kept.push_back(i);
            ASSERT_FALSE(kept.empty()) << "n " << n << " level " << l;
            for (size_t k = 1; k < kept.size(); k++) EXPECT_EQ(kept[k] - kept[k - 1], size_t{1} << l) << "n " << n << " level " << l;
        }
    }
    // even: the scores ignored, the cores C apart (8 cores: C 2 keeps 0, 2, 4, 6; C 4 keeps 2 and 6).
    std::vector<sa::Table::View> const no_alleles;
    std::vector<grid::Core> cores(8);
    for (uint32_t i = 0; i < 8; i++) cores[i] = { 20 + 3 * i, 1 };
    grid::Scratch scratch;
    grid::Stats stats;
    std::vector<uint8_t> levels;
    grid::GeneLevels(cores, no_alleles, {}, levels, scratch, stats, grid::Choice::Even);
    EXPECT_EQ(levels, (std::vector<uint8_t>{ 1, 0, 5, 0, 1, 0, 2, 0 }));
}

// An allele keeps a core a hit if it covers the core's 31-base k-mer (8 flank bases either side of the 15-base core)
// without an edit in the core, an indel in the k-mer, or more than two substitutions in the flanks.
TEST(GridLevels, AnAlleleKeepsACoreWithoutEditsInItsCoreOrIndels) {
    uint32_t const core = 100;  // the k-mer: [92, 123)
    auto keeps = [&](uint16_t begin, uint16_t end, std::vector<sa::Edit> const& edits, uint32_t max_deletion = 0) {
        return grid::AlleleKeeps(View(begin, end, edits), core, max_deletion);
    };
    EXPECT_EQ(keeps(0, 500, {}), 1);
    EXPECT_EQ(keeps(93, 500, {}), -1) << "does not cover the k-mer";
    EXPECT_EQ(keeps(0, 122, {}), -1);
    EXPECT_EQ(keeps(92, 123, {}), 1);
    EXPECT_EQ(keeps(0, 500, { sa::Edit::Substitution(105, 1) }), 0) << "in the core";
    EXPECT_EQ(keeps(0, 500, { sa::Edit::Substitution(100, 1) }), 0);
    EXPECT_EQ(keeps(0, 500, { sa::Edit::Substitution(114, 1) }), 0);
    EXPECT_EQ(keeps(0, 500, { sa::Edit::Substitution(93, 1), sa::Edit::Substitution(120, 2) }), 1) << "two in the flanks";
    EXPECT_EQ(keeps(0, 500, { sa::Edit::Substitution(93, 1), sa::Edit::Substitution(95, 1), sa::Edit::Substitution(120, 2) }), 0);
    EXPECT_EQ(keeps(0, 500, { sa::Edit::Substitution(91, 1), sa::Edit::Substitution(123, 1) }), 1) << "outside the k-mer";
    EXPECT_EQ(keeps(0, 500, { sa::Edit::Indel(sa::kInsertion, 110, 2) }), 0);
    EXPECT_EQ(keeps(0, 500, { sa::Edit::Indel(sa::kInsertion, 92, 2) }), 1) << "inserted before the k-mer's first base";
    EXPECT_EQ(keeps(0, 500, { sa::Edit::Indel(sa::kInsertion, 122, 2) }), 0);
    EXPECT_EQ(keeps(0, 500, { sa::Edit::Indel(sa::kDeletion, 80, 15) }, 15), 0) << "a deletion reaching into the k-mer";
    EXPECT_EQ(keeps(0, 500, { sa::Edit::Indel(sa::kDeletion, 80, 5) }, 5), 1);
    EXPECT_EQ(keeps(0, 500, { sa::Edit::Indel(sa::kDeletion, 123, 3) }, 3), 1);
}

// Without alleles a core's conservation is the prior: its columns' within codes (0: none) up to kPriorCode, averaged.
TEST(GridLevels, ThePriorIsTheCoreColumnsConservation) {
    std::vector<uint8_t> within(40, 0);
    EXPECT_DOUBLE_EQ(grid::Prior(within, 10), grid::kUnknownPrior);
    EXPECT_DOUBLE_EQ(grid::Prior({}, 10), grid::kUnknownPrior);
    std::fill(within.begin() + 10, within.begin() + 25, 15);
    EXPECT_DOUBLE_EQ(grid::Prior(within, 10), 1.0);
    std::fill(within.begin() + 10, within.begin() + 25, 6);
    EXPECT_DOUBLE_EQ(grid::Prior(within, 10), 0.5);
    within[12] = 0;
    EXPECT_DOUBLE_EQ(grid::Prior(within, 10), 0.5) << "a column without a code is left out";
    EXPECT_DOUBLE_EQ(grid::Conservation(4, 4, 0.5), (4 + 2 * 0.5) / 6);
    EXPECT_DOUBLE_EQ(grid::Conservation(0, 0, 0.7), 0.7);
}

// A gene's levels from its alleles and its cores' keys: a core an allele breaks loses to its neighbour, and of two cores
// as conserved the one whose key fewer genes share wins; the counts go to the stats.
TEST(GridLevels, AGeneKeepsTheCoresItsAllelesKeepAndFewerGenesShare) {
    grid::Scratch scratch;
    grid::Stats stats;
    std::vector<uint8_t> levels;
    std::vector<grid::Core> cores = { { 20, 1 }, { 30, 1 }, { 40, 1 }, { 50, 1 } };
    std::vector<sa::Edit> edits = { sa::Edit::Substitution(24, 2) };  // in the first core
    std::vector<sa::Table::View> alleles = { View(0, 300, edits), View(0, 300, {}) };
    grid::GeneLevels(cores, alleles, {}, levels, scratch, stats);
    EXPECT_EQ(levels[0], 0);
    EXPECT_GE(levels[1], 1);
    EXPECT_EQ(stats.cores, 4u);
    EXPECT_EQ(stats.covered, 4u);
    EXPECT_EQ(stats.Kept(0), 4u);
    EXPECT_EQ(stats.Kept(1), 2u);
    EXPECT_EQ(stats.Kept(2), 1u);
    EXPECT_GT(stats.MeanConservation(1), stats.MeanConservation(0));
    // No alleles: the cores tie on conservation; the second's key has 1000 entries, so the first wins.
    cores = { { 20, 1 }, { 30, 1000 } };
    grid::GeneLevels(cores, {}, {}, levels, scratch, stats);
    EXPECT_EQ(levels, (std::vector<uint8_t>{ 5, 0 }));
    cores = { { 20, 1000 }, { 30, 1 } };
    grid::GeneLevels(cores, {}, {}, levels, scratch, stats);
    EXPECT_EQ(levels, (std::vector<uint8_t>{ 0, 5 }));
    // A core whose key the index leaves out never wins.
    cores = { { 20, 0 }, { 30, 1000 } };
    grid::GeneLevels(cores, {}, {}, levels, scratch, stats);
    EXPECT_EQ(levels, (std::vector<uint8_t>{ 0, 5 }));
}

namespace {
    // Maps of 10-base cores (Seedmap(size_t)), as test_PackedIndex's: the core from bit 16, the flex part's halves below
    // and above it.
    constexpr size_t kCoreBases = 10;
    constexpr unsigned kFlexHigh = 16 + 2 * kCoreBases;
    constexpr uint64_t kCoreMask = (uint64_t{1} << (2 * kCoreBases)) - 1;

    uint64_t Key(uint64_t core, uint64_t flex) {
        return (flex >> 16) << kFlexHigh | core << 16 | (flex & 0xffff);
    }

    struct GridValue { uint64_t key, taxid, gene, pos, level; };

    // Keys of 1 to ~300 values, each value with a level of 0 to 5 (half of them 0, as a grid's), positions below 2^17.
    std::vector<GridValue> GridValues(unsigned seed, size_t keys) {
        std::mt19937_64 rng(seed);
        std::vector<GridValue> values;
        for (size_t k = 0; k < keys; k++) {
            uint64_t const core = rng() & kCoreMask;
            size_t const n = k % 7 == 0 ? 1 : k % 11 == 0 ? 2 : k % 23 == 0 ? 150 + rng() % 150 : 3 + rng() % 30;
            uint64_t const flex_base = rng() & 0xffffffffu;
            for (size_t i = 0; i < n; i++) {
                // Some flex cells shared within a key, so that lookups tie.
                uint64_t const flex = i % 4 == 3 ? flex_base : rng() & 0xffffffffu;
                uint64_t const level = std::min<uint64_t>(static_cast<uint64_t>(std::countr_zero(rng() | (uint64_t{1} << 5))), 5);
                values.push_back({ Key(core, flex), 1 + rng() % 143614, 1 + rng() % 168, rng() % 12884, level });
            }
        }
        return values;
    }

    Seedmap::PackedLayout GridLayout(bool levels) {
        auto layout = Seedmap::PackedLayout::For(143614, 168, 12883 + 16);
        if (levels) layout.level_bits = protal::GRID_LEVEL_BITS;
        return layout;
    }

    // An index with grid levels, built as --build does into the packed layout.
    void FillGrid(Seedmap& map, std::vector<GridValue> const& values) {
        map.SetGridLevels();
        for (auto const& v : values) map.CountUpKey(map.MainKey(v.key));
        auto const layout = GridLayout(true);
        map.BuildValuePointers(2, &layout, false);
        for (auto const& v : values) {
            ValueEntry entry;
            entry.Put(v.taxid, v.gene, v.pos);
            entry.value |= v.level << protal::GRID_LEVEL_SHIFT;
            ASSERT_TRUE(map.PutOwned(v.key, entry.value));
        }
    }

    protal::zstd::Params ColumnParams() {
        protal::zstd::Params params;
        params.level = 3;
        params.window_log = 0;
        params.threads = 2;
        params.frame_size = uint64_t{1} << 16;
        return params;
    }

    // Each key's entries (values without levels) and their levels, in order.
    struct KeyEntries { std::vector<uint64_t> values; std::vector<uint32_t> levels, flex; };
    std::map<uint64_t, KeyEntries> EntriesByKey(Seedmap const& map, std::vector<GridValue> const& values) {
        std::map<uint64_t, KeyEntries> out;
        for (auto const& v : values) {
            uint64_t const core = map.MainKey(v.key);
            if (out.contains(core)) continue;
            Seedmap::PackedBlock block;
            if (!map.GetPacked(v.key, block)) {
                out[core];
                continue;
            }
            auto& e = out[core];
            for (uint32_t i = 0; i < block.size; i++) {
                e.values.push_back(map.EntryValue(block, i));
                e.levels.push_back(map.EntryLevel(block, i));
                if (block.flex) e.flex.push_back(Seedmap::FlexCell(block, i));
            }
        }
        return out;
    }
}

// An index with grid levels keeps each entry's level beside it (packed: above its flags; in the file: in its position
// field's top bits), and gives its position without them; a load without room for the levels drops them.
TEST(GridIndex, EntriesKeepTheirLevelsPackedAndInTheFile) {
    auto const values = GridValues(3, 800);
    Seedmap map(kCoreBases);
    FillGrid(map, values);
    EXPECT_TRUE(map.HasGridLevels());
    EXPECT_TRUE(map.KeepsGridLevels());
    std::map<uint64_t, std::vector<uint64_t>> levels_of;  // by core, in placing order
    for (auto const& v : values) levels_of[map.MainKey(v.key)].push_back(v.level);
    auto const built = EntriesByKey(map, values);
    for (auto const& [core, e] : built) {
        EXPECT_EQ(e.levels, std::vector<uint32_t>(levels_of[core].begin(), levels_of[core].end())) << "core " << core;
        for (auto const v : e.values) EXPECT_LT(v & 0xfffff, uint64_t{12884}) << "a position without its level";
    }
    ScratchDir dir("gridtest");
    uint64_t written = 0;
    size_t raw_chunks = 0;
    ASSERT_EQ(map.SaveCompressed(dir / "index.prx.zst", ColumnParams(), written, raw_chunks), "");
    auto const features = Seedmap::IndexFeatures(protal::db::DbFile::OnDisk(dir / "index.prx.zst"));
    ASSERT_TRUE(features);
    EXPECT_TRUE(*features & Seedmap::kFeatureGridLevels);
    for (bool keep : { true, false }) {
        Seedmap loaded(kCoreBases);
        auto const layout = GridLayout(keep);
        loaded.Load(dir / "index.prx.zst", 3, &layout);
        EXPECT_TRUE(loaded.HasGridLevels());
        EXPECT_EQ(loaded.KeepsGridLevels(), keep);
        auto const back = EntriesByKey(loaded, values);
        for (auto const& [core, e] : built) {
            EXPECT_EQ(back.at(core).values, e.values) << "core " << core;
            EXPECT_EQ(back.at(core).levels, keep ? e.levels : std::vector<uint32_t>(e.levels.size(), 0)) << "core " << core;
        }
    }
    // The raw file too.
    map.Save(dir / "index.prx");
    Seedmap raw(kCoreBases);
    auto const layout = GridLayout(true);
    raw.Load(dir / "index.prx", 2, &layout);
    auto const raw_back = EntriesByKey(raw, values);
    for (auto const& [core, e] : built) EXPECT_EQ(raw_back.at(core).levels, e.levels) << "core " << core;
}

// A sparse load of level 2 holds, of each key, its entries of level 2 or more and their flex cells, in order, and no
// key whose entries are all below; its entries are the file's of those levels.
TEST(GridIndex, ASparseLoadHoldsTheEntriesOfItsLevel) {
    auto const values = GridValues(4, 1500);
    Seedmap map(kCoreBases);
    FillGrid(map, values);
    auto const whole = EntriesByKey(map, values);
    ScratchDir dir("gridsparse");
    uint64_t written = 0;
    size_t raw_chunks = 0;
    ASSERT_EQ(map.SaveCompressed(dir / "index.prx.zst", ColumnParams(), written, raw_chunks), "");
    for (uint8_t level : { 1, 2, 4 }) {
        Seedmap sparse(kCoreBases);
        sparse.SetLoadLevel(level);
        auto const layout = GridLayout(true);
        sparse.Load(dir / "index.prx.zst", 3, &layout);
        EXPECT_EQ(sparse.LoadedLevel(), level);
        EXPECT_EQ(sparse.FileEntries(), values.size());
        auto const back = EntriesByKey(sparse, values);
        uint64_t kept = 0, kept_cores = 0;
        for (auto const& [core, e] : whole) {
            KeyEntries want;
            for (size_t i = 0; i < e.values.size(); i++) {
                if (e.levels[i] < level) continue;
                want.values.push_back(e.values[i]);
                want.levels.push_back(e.levels[i]);
                if (!e.flex.empty()) want.flex.push_back(e.flex[i]);
            }
            if (want.values.size() < 2) want.flex.clear();  // a key of one value has no flex cell
            kept += want.values.size();
            kept_cores += !want.values.empty();
            EXPECT_EQ(back.at(core).values, want.values) << "level " << int(level) << " core " << core;
            EXPECT_EQ(back.at(core).levels, want.levels) << "level " << int(level) << " core " << core;
            EXPECT_EQ(back.at(core).flex, want.flex) << "level " << int(level) << " core " << core;
        }
        EXPECT_EQ(sparse.PackedEntries(), kept);
        EXPECT_LT(kept, values.size());
        EXPECT_NE(sparse.MemoryDescription().find("a sparse index: the entries of grid level " + std::to_string(level)), std::string::npos);
        // The cores that keep an entry, of those with entries in the file (the key map's used slots, whole and sparse).
        EXPECT_NE(sparse.MemoryDescription().find("on " + std::to_string(kept_cores) + " of its " + std::to_string(whole.size()) +
                                                  " cores with entries"), std::string::npos) << sparse.MemoryDescription();
    }
}

// Lookups that count a level's entries alone give the same lookups (none for a core without entries of the level) and
// the same seeds, flags and best scores on the whole index as on a sparse index of that level: a run gated with the grid
// pick seeds alike whether it loads the whole index or the sparse one (--sparse_index off or auto).
TEST(GridIndex, LookupsOnTheWholeIndexSeedAsOnTheSparseIndex) {
    auto const values = GridValues(9, 1200);
    Seedmap map(kCoreBases);
    FillGrid(map, values);
    ScratchDir dir("gridlookups");
    uint64_t written = 0;
    size_t raw_chunks = 0;
    ASSERT_EQ(map.SaveCompressed(dir / "index.prx.zst", ColumnParams(), written, raw_chunks), "");
    auto const layout = GridLayout(true);
    Seedmap whole(kCoreBases);
    whole.Load(dir / "index.prx.zst", 2, &layout);
    std::mt19937_64 rng(11);
    for (size_t max_ubiquity : { size_t{3}, size_t{256} }) {
        for (uint8_t level : { 1, 3 }) {
            Seedmap sparse(kCoreBases);
            sparse.SetLoadLevel(level);
            sparse.Load(dir / "index.prx.zst", 2, &layout);
            protal::KmerLookupSM on_whole(whole, max_ubiquity), on_sparse(sparse, max_ubiquity);
            size_t compared = 0, seeded = 0;
            for (size_t q = 0; q < values.size(); q += 2) {
                // A query of the value's core, its flex part as stored or with some bases changed.
                uint64_t flex = (values[q].key >> kFlexHigh) << 16 | (values[q].key & 0xffff);
                if (q % 3 == 1) flex ^= uint64_t{1} << (rng() % 32);
                if (q % 5 == 2) flex = rng() & 0xffffffffu;
                size_t query = Key(map.MainKey(values[q].key), flex);
                std::vector<protal::LookupPointer> a, b;
                on_whole.Get(a, query, 7, level);
                on_sparse.Get(b, query, 7, level);
                ASSERT_EQ(a.size(), b.size()) << "query " << q;
                if (a.empty()) continue;
                EXPECT_EQ(a[0].level_entries, b[0].size);
                EXPECT_EQ(b[0].level_entries, b[0].size);
                // A key of one value in the sparse index has no flex cells: the gate scores it from its gene
                // (ChainAnchorFinder::SingleValueScore), as it would the one cell the whole index has for it.
                if (b[0].flex != nullptr) EXPECT_EQ(on_whole.BestScore(a[0]), on_sparse.BestScore(b[0])) << "query " << q;
                protal::LookupList sa_seeds, sb_seeds;
                bool const ra = on_whole.GetFromLookup(sa_seeds, a[0]), rb = on_sparse.GetFromLookup(sb_seeds, b[0]);
                EXPECT_EQ(ra, rb) << "query " << q;
                ASSERT_EQ(sa_seeds.size(), sb_seeds.size()) << "query " << q;
                for (size_t s = 0; s < sa_seeds.size(); s++) {
                    EXPECT_EQ(sa_seeds[s].ToString(), sb_seeds[s].ToString()) << "query " << q;
                    EXPECT_EQ(sa_seeds[s].unique, sb_seeds[s].unique);
                    EXPECT_EQ(sa_seeds[s].unique_dist_two, sb_seeds[s].unique_dist_two);
                    EXPECT_EQ(sa_seeds[s].single, sb_seeds[s].single);
                }
                compared++;
                seeded += !sa_seeds.empty();
            }
            EXPECT_GT(compared, 100u);
            EXPECT_GT(seeded, 50u);
        }
    }
}

// A sparse load of an index without grid levels stops with a message (exit 8): it has no levels to keep.
TEST(GridIndex, ASparseLoadNeedsGridLevels) {
    Seedmap map(kCoreBases);
    std::vector<GridValue> values = GridValues(2, 50);
    for (auto const& v : values) map.CountUpKey(map.MainKey(v.key));
    auto const layout = GridLayout(false);
    map.BuildValuePointers(1, &layout, false);
    for (auto const& v : values) {
        ValueEntry entry;
        entry.Put(v.taxid, v.gene, v.pos);
        ASSERT_TRUE(map.PutOwned(v.key, entry.value));
    }
    ScratchDir dir("gridnone");
    uint64_t written = 0;
    size_t raw_chunks = 0;
    ASSERT_EQ(map.SaveCompressed(dir / "index.prx.zst", ColumnParams(), written, raw_chunks), "");
    auto const features = Seedmap::IndexFeatures(protal::db::DbFile::OnDisk(dir / "index.prx.zst"));
    ASSERT_TRUE(features);
    EXPECT_FALSE(*features & Seedmap::kFeatureGridLevels);
    EXPECT_EXIT({
        Seedmap sparse(kCoreBases);
        sparse.SetLoadLevel(2);
        auto const with_levels = GridLayout(true);
        sparse.Load(dir / "index.prx.zst", 1, &with_levels);
    }, testing::ExitedWithCode(8), "has no grid levels");
}
