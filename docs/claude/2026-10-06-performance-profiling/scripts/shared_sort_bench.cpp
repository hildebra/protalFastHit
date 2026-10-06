// shared_sort_bench.cpp - the seed sort at GTDB r226's seed mix: all seeds sorted as 128-bit keys (ChainAnchorFinder::Sort
// up to cdce3c0) against SharedSeeds::Sort (byte counters per hash of taxon and gene, the seeds of counters of two or more
// sorted) and its first version, ExactSharedSeeds below (taxon and gene counted exactly in an open-addressing table).
// The seed lists are sort_bench.cpp's (51 seeds per mate, r226 43, from ~22 lookups at random taxa and genes of r226's 143,614
// genomes x 169 genes, a third of the mates with an anchor's 22 seeds among them): ~14% of the seeds shared, as the
// twelfth cluster run counted
// 12.6% (pe). All three must give FindPairs the same groups: the runs of two or more of the full sort, in order.
//
// Build: cc8.sh shared_sort_bench.cpp shared_sort_bench; run: shared_sort_bench <mates> <rounds> [variant 0-2]
#include <iostream>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <random>
#include <vector>

#include "SequenceUtils/GenomeLoader.h"
#include "Hash/KmerLookup.h"
#include "Core/AlignmentStrategy.h"
#include "Core/ChainAnchorFinder.h"

using protal::LookupResult;
using protal::SeedList;

namespace {
    class ExactSharedSeeds {
    public:
        // The shared seeds of `seeds`, sorted, into `shared`; `seeds` is left as it is.
        void Sort(SeedList const& seeds, SeedList& shared) {
            shared.clear();
            m_sort_keys.clear();
            if (seeds.size() < 2) return;
            Count(seeds);
            for (size_t i = 0; i < seeds.size(); i++) {
                if (m_slots[m_slot_of[i]].count > 1) m_sort_keys.push_back(seeds[i].SortKey());
            }
            std::sort(m_sort_keys.begin(), m_sort_keys.end());
            shared.resize(m_sort_keys.size());
            for (size_t i = 0; i < m_sort_keys.size(); i++) shared[i] = LookupResult::FromSortKey(m_sort_keys[i]);
        }

    private:
        // Each seed's slot in m_slot_of, and the slot's count of the read's seeds of its taxon and gene.
        void Count(SeedList const& seeds) {
            size_t const n = seeds.size();
            size_t bits = 6;
            while ((size_t{ 1 } << bits) < 2 * n) bits++;
            if (m_slots.size() < (size_t{ 1 } << bits)) {
                m_slots.assign(size_t{ 1 } << bits, Slot{});
                m_stamp = 0;
            }
            if (++m_stamp == 0) {  // wrapped: every old stamp reads as 0, so start over at 1
                for (auto& slot : m_slots) slot.stamp = 0;
                m_stamp = 1;
            }
            size_t const mask = (size_t{ 1 } << bits) - 1;
            m_slot_of.resize(n);
            for (size_t i = 0; i < n; i++) {
                uint64_t const key = (static_cast<uint64_t>(seeds[i].taxid) << 32) | seeds[i].geneid;
                size_t slot = static_cast<size_t>((key * 0x9e3779b97f4a7c15ull) >> (64 - bits));
                while (m_slots[slot].stamp == m_stamp && m_slots[slot].key != key) slot = (slot + 1) & mask;
                if (m_slots[slot].stamp != m_stamp) m_slots[slot] = Slot{ key, m_stamp, 0 };
                m_slots[slot].count++;
                m_slot_of[i] = static_cast<uint32_t>(slot);
            }
        }

        struct Slot {
            uint64_t key = 0;
            uint32_t stamp = 0;
            uint32_t count = 0;
        };
        std::vector<Slot> m_slots;
        std::vector<uint32_t> m_slot_of;  // per seed of the read, its slot
        uint32_t m_stamp = 0;
        std::vector<LookupResult::SortKeyType> m_sort_keys;  // kept for their capacity
    };
}

int main(int argc, char** argv) {
    size_t const mates = argc > 1 ? std::strtoull(argv[1], nullptr, 10) : 200000;
    int const rounds = argc > 2 ? std::atoi(argv[2]) : 3;
    int const only = argc > 3 ? std::atoi(argv[3]) : -1;  // one variant (for callgrind), else all three
    constexpr uint64_t kGenes = 169;
    std::mt19937_64 rng(5);
    std::vector<std::vector<LookupResult>> lists(mates);
    size_t total = 0;
    for (auto& list : lists) {
        bool const present = rng() % 3 == 0;
        uint32_t const taxid = 1 + rng() % 143614, gene = rng() % kGenes, start = rng() % 1000;
        for (int l = 0; l < 22; l++) {
            uint16_t const readpos = static_cast<uint16_t>(8 + rng() % 127);
            int const n = static_cast<int>(rng() % 5);
            for (int j = 0; j < n; j++) {
                list.emplace_back(1 + rng() % 143614, rng() % kGenes, rng() % 16000, readpos, rng() & 1, rng() & 1);
            }
            if (present) list.emplace_back(taxid, gene, start + readpos, readpos, true, rng() & 1);
        }
        std::shuffle(list.begin(), list.end(), rng);
        total += list.size();
    }

    std::vector<LookupResult> work, shared;
    std::vector<LookupResult::SortKeyType> k128;
    protal::SharedSeeds shared_seeds;
    ExactSharedSeeds exact_shared_seeds;
    auto checksum = [](LookupResult const& s, uint64_t c) {
        return c * 1000003 + (uint64_t{s.taxid} << 40 ^ uint64_t{s.geneid} << 20 ^ s.genepos ^ uint64_t{s.readpos} << 50 ^
                              uint64_t{s.unique} << 62 ^ uint64_t{s.unique_dist_two} << 63);
    };
    size_t shared_total = 0, sorted_total = 0;
    for (int r = 0; r < rounds; r++) {
        for (int variant = 0; variant < 3; variant++) {
            if (only >= 0 && variant != only) continue;
            uint64_t c = 0;
            double ns = 0;
            shared_total = 0;
            sorted_total = 0;
            for (auto const& list : lists) {
                work = list;
                auto const start = std::chrono::steady_clock::now();
                if (variant == 0) {  // up to cdce3c0: every seed sorted
                    k128.resize(work.size());
                    for (size_t i = 0; i < work.size(); i++) k128[i] = work[i].SortKey();
                    std::sort(k128.begin(), k128.end());
                    for (size_t i = 0; i < work.size(); i++) work[i] = LookupResult::FromSortKey(k128[i]);
                } else if (variant == 1) {
                    shared_seeds.Sort(work, shared);
                } else {
                    exact_shared_seeds.Sort(work, shared);
                }
                ns += std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - start).count();
                // What FindPairs sees: the runs of one taxon and gene of two or more, in order.
                auto const& sorted = variant == 0 ? work : shared;
                sorted_total += sorted.size();
                for (size_t i = 0; i < sorted.size();) {
                    size_t j = i + 1;
                    while (j < sorted.size() && sorted[j].FromSameSequence(sorted[i])) j++;
                    if (j - i > 1) for (size_t x = i; x < j; x++) c = checksum(sorted[x], c);
                    shared_total += j - i > 1 ? j - i : 0;
                    i = j;
                }
            }
            std::printf("%s: %.1f ns per seed, %.1f%% of %zu seeds shared, %.1f%% sorted, checksum %016llx\n",
                        variant == 0 ? "all sorted     " : variant == 1 ? "byte counters  " : "exact counting ",
                        ns / static_cast<double>(total), 100.0 * static_cast<double>(shared_total) / static_cast<double>(total), total,
                        100.0 * static_cast<double>(sorted_total) / static_cast<double>(total), static_cast<unsigned long long>(c));
        }
    }
    return 0;
}
