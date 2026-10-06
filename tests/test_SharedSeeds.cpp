// Unit tests for the seeds the anchor finder sorts (SharedSeeds, ChainAnchorFinder.h): those whose taxon and gene another
// seed of the read has, and a few others, sorted, so that FindPairs finds the same groups in the same order as in all
// seeds sorted. Against the reference: all seeds sorted by their 128-bit key, the runs of one taxon and gene of two or more.
#include <gtest/gtest.h>
#include <algorithm>
#include <random>
#include <vector>
#include "SequenceUtils/GenomeLoader.h"
#include "Hash/KmerLookup.h"
#include "Core/AlignmentStrategy.h"
#include "Core/ChainAnchorFinder.h"

using namespace protal;

namespace {
    bool BySortKey(LookupResult const& a, LookupResult const& b) { return a.SortKey() < b.SortKey(); }

    // What FindPairs pairs: the runs of two or more seeds of one taxon and gene of a sorted list, in order.
    SeedList Runs(SeedList const& sorted) {
        SeedList runs;
        for (size_t i = 0; i < sorted.size();) {
            size_t j = i + 1;
            while (j < sorted.size() && sorted[j].FromSameSequence(sorted[i])) j++;
            if (j - i > 1) runs.insert(runs.end(), sorted.begin() + static_cast<long>(i), sorted.begin() + static_cast<long>(j));
            i = j;
        }
        return runs;
    }

    bool Same(LookupResult const& a, LookupResult const& b) { return a.SortKey() == b.SortKey(); }
}

// Reads of 0 to 5000 seeds, from few or many taxa and genes (all shared, none shared, mixed), with a k-mer's hits at
// several places of one gene and seeds repeated exactly, the largest taxon and gene ids among them; one SharedSeeds for
// all. The sorted seeds are sorted, seeds of the read, hold every shared seed, and FindPairs finds the same runs in them
// as in all seeds sorted; the seeds of genes of their own among them are few.
TEST(SharedSeeds, FindPairsSeesTheRunsOfAllSeedsSorted) {
    std::mt19937 rng(41);
    SharedSeeds shared_seeds;
    SeedList sorted_seeds;
    size_t seeds_total = 0, shared_total = 0, singletons = 0, singletons_sorted = 0;
    for (int n = 0; n < 4000; n++) {
        size_t const count = n % 97 == 0 ? 1000 + rng() % 4000 : rng() % 120;
        uint32_t const taxa = 1 + rng() % (n % 3 == 0 ? 3 : 2000);
        uint32_t const genes = 1 + rng() % (n % 5 == 0 ? 2 : 120);
        SeedList seeds;
        for (size_t i = 0; i < count; i++) {
            if (!seeds.empty() && rng() % 10 == 0) {  // a repeat: the same seed, or another place of the same gene
                LookupResult again = seeds[rng() % seeds.size()];
                if (rng() % 2) again.genepos = rng() % 3000;
                seeds.push_back(again);
                continue;
            }
            uint32_t taxid = rng() % taxa, geneid = rng() % genes;
            if (rng() % 50 == 0) taxid = UINT32_MAX - rng() % 2;
            if (rng() % 50 == 0) geneid = (1u << 20) - 1;
            seeds.emplace_back(taxid, geneid, rng() % 3000, rng() % 150, rng() % 2, rng() % 2);
        }
        SeedList const before = seeds;
        shared_seeds.Sort(seeds, sorted_seeds);
        ASSERT_EQ(seeds.size(), before.size());
        for (size_t i = 0; i < seeds.size(); i++) ASSERT_TRUE(Same(seeds[i], before[i])) << "read " << n << ": seeds changed";
        ASSERT_TRUE(std::is_sorted(sorted_seeds.begin(), sorted_seeds.end(), BySortKey)) << "read " << n;

        SeedList all = seeds;
        std::sort(all.begin(), all.end(), BySortKey);
        SeedList const expected = Runs(all), got = Runs(sorted_seeds);
        ASSERT_EQ(got.size(), expected.size()) << "read " << n << " of " << count << " seeds";
        for (size_t i = 0; i < expected.size(); i++) ASSERT_TRUE(Same(got[i], expected[i])) << "read " << n << ", seed " << i;
        // Every sorted seed is a seed of the read: the sorted ones are a sub-multiset of all, as sorted lists.
        ASSERT_TRUE(std::includes(all.begin(), all.end(), sorted_seeds.begin(), sorted_seeds.end(), BySortKey)) << "read " << n;

        seeds_total += count;
        shared_total += expected.size();
        singletons += count - expected.size();
        singletons_sorted += sorted_seeds.size() - got.size();
    }
    std::cout << seeds_total << " seeds, " << shared_total << " of them shared; of the " << singletons << " others "
              << singletons_sorted << " sorted along" << std::endl;
    EXPECT_GT(shared_total, seeds_total / 10);
    EXPECT_LT(shared_total, seeds_total * 9 / 10);
    EXPECT_LT(singletons_sorted, singletons / 10);
}

// No seeds or one: nothing to sort; two of one gene: both; two of different genes: no run for FindPairs.
TEST(SharedSeeds, OneSeedOrNoneHasNoneShared) {
    SharedSeeds shared_seeds;
    SeedList sorted_seeds{ LookupResult(1, 2, 3, 4) };
    shared_seeds.Sort({}, sorted_seeds);
    EXPECT_TRUE(sorted_seeds.empty());
    shared_seeds.Sort({ LookupResult(1, 2, 3, 4) }, sorted_seeds);
    EXPECT_TRUE(sorted_seeds.empty());
    shared_seeds.Sort({ LookupResult(1, 2, 3, 4), LookupResult(1, 2, 3, 4) }, sorted_seeds);
    EXPECT_EQ(sorted_seeds.size(), 2u);
    shared_seeds.Sort({ LookupResult(1, 2, 3, 4), LookupResult(1, 3, 3, 4) }, sorted_seeds);
    EXPECT_TRUE(Runs(sorted_seeds).empty());
}
