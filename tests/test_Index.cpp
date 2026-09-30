// Unit tests for loading index.prx: malformed files must stop with a clear message (exit 8)
// instead of loading garbage sizes and reading out of bounds later.
#include <gtest/gtest.h>
#include <cstdint>
#include <random>
#include <sstream>
#include <string>
#include "Hash/Seedmap.h"
#include "SequenceUtils/Minimizer.h"

using protal::Seedmap;

namespace {
    template<typename T>
    void Put(std::ostream& os, T value) { os.write(reinterpret_cast<char const*>(&value), sizeof(T)); }

    // A valid file header followed by the six layout fields of an index.
    std::string Header(size_t keymap_size, size_t keymap_size_total, size_t ctrl_block_bytes,
                       size_t keys_per_block, size_t bitshift, size_t values_size) {
        std::ostringstream os;
        Put<uint64_t>(os, Seedmap::kFileMagic);
        Put<uint32_t>(os, Seedmap::kIndexVersionMajor);
        Put<uint32_t>(os, Seedmap::kIndexVersionMinor);
        Put<uint32_t>(os, Seedmap::kIndexVersionPatch);
        for (size_t v : { keymap_size, keymap_size_total, ctrl_block_bytes, keys_per_block, bitshift, values_size }) Put<size_t>(os, v);
        return os.str();
    }

    constexpr size_t kKeys = size_t{1} << 30;                      // 15-mer core
    constexpr size_t kTotal = kKeys + ((kKeys / 8) + 1) * 4;         // 8 keys per block, 4 cells each

    void LoadFrom(std::string bytes) {
        Seedmap map;
        std::istringstream in(std::move(bytes));
        map.Load(in, "test.prx");
    }
}

TEST(IndexLoad, RejectsATruncatedHeader) {
    auto header = Header(kKeys, kTotal, 8, 8, 3, 0);
    EXPECT_EXIT(LoadFrom(header.substr(0, header.size() - 4)), testing::ExitedWithCode(8), "ends inside its header");
}

TEST(IndexLoad, RejectsAnotherKeySize) {
    EXPECT_EXIT(LoadFrom(Header(size_t{1} << 28, kTotal, 8, 8, 3, 0)), testing::ExitedWithCode(8), "this protal expects");
}

TEST(IndexLoad, RejectsAnInconsistentLayout) {
    EXPECT_EXIT(LoadFrom(Header(kKeys, kTotal, 8, 8, 4, 0)), testing::ExitedWithCode(8), "control block layout");
    EXPECT_EXIT(LoadFrom(Header(kKeys, kTotal + 1, 8, 8, 3, 0)), testing::ExitedWithCode(8), "the layout implies");
}

TEST(IndexLoad, RejectsMissingData) {
    // Header claims a full key map and 10 values, but no data follows.
    EXPECT_EXIT(LoadFrom(Header(kKeys, kTotal, 8, 8, 3, 10)), testing::ExitedWithCode(8), "truncated or corrupt");
}

namespace {
    // Closed-syncmer selection written out plainly, for one s-mer mask.
    bool ReferenceSyncmer(uint64_t key, uint32_t mask, int k = 15, int s = 7, int t = 2) {
        int n = k - s + 1;
        uint64_t min = UINT64_MAX;
        int min_index = 0;
        for (int i = 0; i < n; i++) {
            uint64_t smer = (key >> (n * 2 - (i + 1) * 2)) & mask;
            if (smer < min) { min = smer; min_index = i; }
        }
        return min_index == t || min_index == n - 1 - t;
    }
}

TEST(ClosedSyncmer, LegacyAndFullMasksMatchTheirDefinition) {
    protal::ClosedSyncmer legacy{15, 7, 2, false};
    protal::ClosedSyncmer full{15, 7, 2, true};
    EXPECT_FALSE(legacy.UsesFullSmerMask());
    EXPECT_TRUE(full.UsesFullSmerMask());

    std::mt19937_64 rng(3);
    int legacy_mismatches = 0, full_mismatches = 0, differ = 0;
    for (int i = 0; i < 200000; i++) {
        uint64_t key = rng() & ((uint64_t{1} << 30) - 1);
        bool l = legacy(key), f = full(key);
        legacy_mismatches += l != ReferenceSyncmer(key, 0xFF);   // format 1: last 4 bases of each 7-mer
        full_mismatches += f != ReferenceSyncmer(key, 0x3FFF);   // format 2: whole 7-mers
        differ += l != f;
    }
    EXPECT_EQ(legacy_mismatches, 0);
    EXPECT_EQ(full_mismatches, 0);
    EXPECT_GT(differ, 0);  // the two formats really sample differently
}

TEST(IndexHeader, FeaturesRoundTripAndFormat1IsLegacy) {
    std::stringstream v2;
    {
        Seedmap built;
        built.SaveHeader(v2);
    }
    Seedmap loaded;
    loaded.LoadHeader(v2);
    EXPECT_TRUE(loaded.UsesFullSyncmerMask());
    EXPECT_TRUE(loaded.ChecksSingleEntryUniques());

    std::stringstream v1;
    Put<uint64_t>(v1, Seedmap::kFileMagic);
    Put<uint32_t>(v1, Seedmap::kIndexVersionMajor);
    Put<uint32_t>(v1, Seedmap::kIndexVersionMinor);
    Put<uint32_t>(v1, Seedmap::kIndexVersionPatch);
    Seedmap legacy;
    legacy.LoadHeader(v1);
    EXPECT_FALSE(legacy.UsesFullSyncmerMask());
    EXPECT_FALSE(legacy.ChecksSingleEntryUniques());
    EXPECT_NE(legacy.FeatureDescription().find("legacy"), std::string::npos);
}

// The unique k-mer statistics find, for large cores, the values with another within flex distance 1
// by sorting with each position masked (FlexNeighbours); it must agree with comparing all pairs.
TEST(UniqueKmers, FlexNeighboursMatchAllPairs) {
    std::mt19937 rng(7);
    for (size_t n : { 2, 3, 17, 300, 1000 }) {
        std::vector<uint32_t> flex(n);
        for (auto& cell : flex) cell = rng();
        // Some parts one position away from another, some equal, some two positions away.
        for (size_t i = 1; i < n; i += 5) flex[i] = flex[i - 1] ^ (1u << (2 * (rng() % 16)));
        for (size_t i = 3; i < n; i += 11) flex[i] = flex[i - 2];
        for (size_t i = 4; i < n; i += 13) flex[i] = flex[i - 1] ^ (3u << 2) ^ (2u << 20);
        std::vector<uint8_t> close;
        std::vector<std::pair<uint32_t, uint32_t>> keyed;
        Seedmap::FlexNeighbours(flex.data(), n, 16, close, keyed);
        size_t near = 0;
        for (size_t e = 0; e < n; e++) {
            bool expected = false;
            for (size_t o = 0; o < n; o++) expected |= o != e && Seedmap::Similarity(flex[e], flex[o]) + 1 >= 16;
            EXPECT_EQ(bool(close[e]), expected) << "n " << n << ", part " << e;
            near += expected;
        }
        if (n >= 17) EXPECT_GT(near, 0u);
        if (n >= 17) EXPECT_LT(near, n);
    }
}
