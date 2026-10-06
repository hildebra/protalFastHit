// sort_bench.cpp - ChainAnchorFinder::Sort at GTDB r226's seed counts (43 seeds per mate, ~2 per lookup from ~22
// lookups): the tree's 128-bit keys (LookupResult::SortKey, std::sort, back to LookupResult) against 64-bit keys
// holding the same fields in the same order where they fit (taxon and gene as taxid * genes + gene, 25 bits at r226;
// read position 16; gene position 14; 2 flags = 57 bits), and the same 64-bit keys built straight from a packed index
// entry (no taxid / gene split before the sort). Every variant must give the same seed order.
//
// Build: cc.sh sort_bench.cpp sort_bench; run: sort_bench <mates> <rounds>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <random>
#include <vector>

#include "Hash/KmerLookup.h"

using protal::LookupResult;

namespace {
    constexpr uint64_t kGenes = 169;
    constexpr unsigned kPosBits = 14, kReadBits = 16;
    uint64_t const kReciprocal = static_cast<uint64_t>((static_cast<__uint128_t>(1) << 64) / kGenes + 1);

    inline uint64_t Key64(LookupResult const& s) {
        uint64_t const tg = uint64_t{s.taxid} * kGenes + s.geneid;
        return tg << (kReadBits + kPosBits + 2) | uint64_t{s.readpos} << (kPosBits + 2) | uint64_t{s.genepos} << 2 |
               uint64_t{s.unique} << 1 | uint64_t{s.unique_dist_two};
    }
    inline LookupResult From64(uint64_t key) {
        uint64_t const tg = key >> (kReadBits + kPosBits + 2);
        uint64_t const taxid = static_cast<uint64_t>((static_cast<__uint128_t>(tg) * kReciprocal) >> 64);
        LookupResult s;
        s.taxid = static_cast<uint32_t>(taxid);
        s.geneid = static_cast<uint32_t>(tg - taxid * kGenes);
        s.readpos = static_cast<uint16_t>(key >> (kPosBits + 2));
        s.genepos = static_cast<uint32_t>((key >> 2) & ((1u << kPosBits) - 1));
        s.unique = (key >> 1) & 1;
        s.unique_dist_two = key & 1;
        return s;
    }
}

int main(int argc, char** argv) {
    size_t const mates = argc > 1 ? std::strtoull(argv[1], nullptr, 10) : 200000;
    int const rounds = argc > 2 ? std::atoi(argv[2]) : 3;
    std::mt19937_64 rng(5);
    // Per mate: 22 lookups at read positions 0-134, each giving 0-4 seeds (mean ~2); a third of the mates from a
    // species in the database: their lookups also hit that taxon's gene on one diagonal (an anchor's seeds).
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
    std::printf("%zu mates, %.1f seeds per mate\n", mates, static_cast<double>(total) / mates);

    std::vector<LookupResult> work;
    std::vector<LookupResult::SortKeyType> k128;
    std::vector<uint64_t> k64;
    auto checksum = [](std::vector<LookupResult> const& v, uint64_t c) {
        for (auto const& s : v) c = c * 1000003 + (uint64_t{s.taxid} << 40 ^ uint64_t{s.geneid} << 20 ^ s.genepos ^ uint64_t{s.readpos} << 50 ^ uint64_t{s.unique} << 62 ^ uint64_t{s.unique_dist_two} << 63);
        return c;
    };
    for (int r = 0; r < rounds; r++) {
        for (int variant = 0; variant < 2; variant++) {
            uint64_t c = 0;
            double ns = 0;
            for (auto const& list : lists) {
                work = list;
                auto const start = std::chrono::steady_clock::now();
                if (variant == 0) {  // the tree's Sort
                    k128.resize(work.size());
                    for (size_t i = 0; i < work.size(); i++) k128[i] = work[i].SortKey();
                    std::sort(k128.begin(), k128.end());
                    for (size_t i = 0; i < work.size(); i++) work[i] = LookupResult::FromSortKey(k128[i]);
                } else {
                    k64.resize(work.size());
                    for (size_t i = 0; i < work.size(); i++) k64[i] = Key64(work[i]);
                    std::sort(k64.begin(), k64.end());
                    for (size_t i = 0; i < work.size(); i++) work[i] = From64(k64[i]);
                }
                ns += std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - start).count();
                c = checksum(work, c);
            }
            std::printf("%s: %.1f ns per seed, checksum %016llx\n", variant == 0 ? "128-bit" : " 64-bit", ns / total, static_cast<unsigned long long>(c));
        }
    }
    return 0;
}
