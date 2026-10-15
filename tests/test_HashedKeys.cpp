// The hashed key map (HashedKeys.h, --hashed_keys): the compact key map's used cores in a hash table of lines give every
// core the slots the compact map gives it (none to a core without values), so every lookup reads the same values; also
// where lines fill up and cores lie past their home line, and the table is the same on any number of threads.
#include <gtest/gtest.h>
#include <cstdint>
#include <cstring>
#include <random>
#include <vector>
#include "Hash/CompactKeys.h"
#include "Hash/HashedKeys.h"

using protal::CompactKeys;
using protal::HashedKeys;

namespace {
    // A compact key map of `slots[c]` slots for each used core c (0: unused), the values from slot `first` on.
    void MakeCompact(CompactKeys& compact, std::vector<uint64_t> const& slots, uint64_t first) {
        uint64_t const cores = slots.size();
        compact.Init(cores);
        for (uint64_t c = 0; c < cores; c++) {
            if (slots[c]) compact.MarkUsed(c);
        }
        uint64_t const used = compact.CountRanks();
        compact.AllocateStarts();
        uint64_t next = first, rank = 0;
        std::vector<uint64_t> start(cores, 0);
        for (uint64_t c = 0; c < cores; c++) {
            if (!slots[c]) continue;
            start[c] = next;
            if (rank % CompactKeys::kGroup == 0) compact.SetBase(rank, next);
            next += slots[c];
            rank++;
        }
        if (used % CompactKeys::kGroup == 0) compact.SetBase(used, next);
        rank = 0;
        for (uint64_t c = 0; c < cores; c++) {
            if (slots[c]) ASSERT_TRUE(compact.SetOffset(rank++, start[c]));
        }
        ASSERT_TRUE(compact.SetOffset(used, next));
    }

    // Used cores at about `density`, in runs as related genomes' cores cluster, each of 1 to 3,071 slots.
    std::vector<uint64_t> ClusteredSlots(uint64_t cores, double density, uint64_t seed) {
        std::mt19937_64 rng(seed);
        std::vector<uint64_t> slots(cores, 0);
        std::uniform_real_distribution<double> u(0.0, 1.0);
        for (uint64_t c = 0; c < cores; c++) {
            if (u(rng) < density / 4) {
                for (uint64_t k = 0; k < 4 + rng() % 5 && c + k < cores; k++) slots[c + k] = 1 + rng() % 3071;
                c += 8;
            }
        }
        slots[0] = 3071;
        slots[cores - 1] = 1;
        return slots;
    }

    void ExpectSameSlots(CompactKeys const& compact, HashedKeys const& hashed, uint64_t cores) {
        uint64_t found = 0;
        for (uint64_t c = 0; c < cores; c++) {
            uint64_t s1 = 0, n1 = 0, s2 = 0, n2 = 0;
            bool const f1 = compact.Locate(c, s1, n1), f2 = hashed.Locate(c, s2, n2);
            ASSERT_EQ(f1, f2) << "core " << c;
            if (f1) {
                ASSERT_EQ(s1, s2) << "core " << c;
                ASSERT_EQ(n1, n2) << "core " << c;
                found++;
            }
        }
        EXPECT_EQ(found, hashed.UsedCores());
    }
}

// 2^20 cores, ~9% used in clusters (r226's density), the values from slot 2^33 on (r226's index has billions of slots):
// every core the compact map's slots, at most 0.75 used cores per entry.
TEST(HashedKeys, EveryCoreGetsTheCompactMapsSlots) {
    unsigned const bits = 20;
    uint64_t const cores = uint64_t{1} << bits;
    auto const slots = ClusteredSlots(cores, 0.092, 7);
    CompactKeys compact;
    MakeCompact(compact, slots, uint64_t{1} << 33);
    HashedKeys hashed;
    ASSERT_TRUE(hashed.Build(compact, bits, 3));
    EXPECT_EQ(hashed.UsedCores(), compact.UsedCores());
    EXPECT_LE(hashed.Load(), HashedKeys::kMaxLoad);
    EXPECT_GT(hashed.Load(), HashedKeys::kMaxLoad / 2);
    ExpectSameSlots(compact, hashed, cores);
}

// A crowded table (0.9 used cores per entry): many cores lie past their home line, and a lookup of a core without values
// reads on past full lines; every core still gets the compact map's slots.
TEST(HashedKeys, ACrowdedTableFindsEveryCoreAndNoOther) {
    unsigned const bits = 16;
    uint64_t const cores = uint64_t{1} << bits;
    std::mt19937_64 rng(11);
    std::vector<uint64_t> slots(cores, 0);
    uint64_t const want = static_cast<uint64_t>(0.9 * HashedKeys::kEntries * 512);  // 2^9 lines at load 0.9
    for (uint64_t placed = 0; placed < want;) {
        uint64_t const c = rng() % cores;
        if (slots[c]) continue;
        slots[c] = 1 + rng() % 3071;
        placed++;
    }
    CompactKeys compact;
    MakeCompact(compact, slots, 5);
    HashedKeys hashed;
    ASSERT_TRUE(hashed.Build(compact, bits, 2, 0.9));
    EXPECT_EQ(hashed.LineBits(), 9u);
    EXPECT_GT(hashed.PastHome(), 0u);
    EXPECT_GT(hashed.MaxDisplacement(), 0u);
    ExpectSameSlots(compact, hashed, cores);
}

// The partitions are filled in core order and the cores spilled past them after all, so the table is the same on one
// thread and on four.
TEST(HashedKeys, TheSameTableOnAnyNumberOfThreads) {
    unsigned const bits = 18;
    uint64_t const cores = uint64_t{1} << bits;
    auto const slots = ClusteredSlots(cores, 0.3, 5);
    CompactKeys compact;
    MakeCompact(compact, slots, 0);
    HashedKeys one, four;
    ASSERT_TRUE(one.Build(compact, bits, 1));
    ASSERT_TRUE(four.Build(compact, bits, 4));
    ASSERT_EQ(one.Bytes(), four.Bytes());
    EXPECT_EQ(std::memcmp(one.LineData(), four.LineData(), one.Bytes()), 0);
    ExpectSameSlots(compact, four, cores);
}

// The mix of the cores is a bijection: no two cores of 2^16 share a home line and tag.
TEST(HashedKeys, TheMixIsABijection) {
    unsigned const bits = 16;
    uint64_t const cores = uint64_t{1} << bits;
    std::vector<uint64_t> slots(cores, 0);
    slots[3] = 1;
    CompactKeys compact;
    MakeCompact(compact, slots, 0);
    HashedKeys hashed;
    ASSERT_TRUE(hashed.Build(compact, bits, 1));
    std::vector<bool> seen(cores, false);
    for (uint64_t c = 0; c < cores; c++) {
        uint64_t const m = hashed.Mix(c);
        ASSERT_LT(m, cores);
        ASSERT_FALSE(seen[m]) << "core " << c;
        seen[m] = true;
    }
}
