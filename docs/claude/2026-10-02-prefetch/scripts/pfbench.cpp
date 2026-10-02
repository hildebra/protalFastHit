// pfbench - experiment only: protal's seeding (ChainAnchorFinder::FindSeeds: the key map lookups of a read's
// k-mers, then the values of the smallest lookups) on the k-mers of real reads against the real index, with and
// without software prefetching. The variants run alternately, one batch of reads each, so that they see the same
// state of the machine; time is the thread's CPU time (stalls included).
//
//   pfbench DATABASE.protal READS.fq[.gz] [READS [ROUNDS]]
//
// V0 as protal: Get for every k-mer, sort the lookups by size, GetFromLookup until enough seeds.
// V1 the read's own control blocks prefetched first, then V0 (round 2's prototype, without the values).
// V2 the control blocks of the read two ahead prefetched, then V0.
// V3 two stages: the control blocks of the read two ahead prefetched; the next read's lookups made (Get, its blocks
//    in cache by then) and their values prefetched; this read's precomputed lookups sorted and read.
// V4 V3 with the blocks four reads ahead and the lookups two ahead.
// Everything outside protal first, so that the define opens only protal's classes (Seedmap's key map).
#include <bits/stdc++.h>
#include <emmintrin.h>
#include <fcntl.h>
#include <omp.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zstd.h>
#include <robin/robin_set.h>
#include <robin/robin_map.h>
#include <sparse_map.h>
#include <sparse_set.h>
#define private public
#include "Hash/Seedmap.h"
#include "Hash/KmerLookup.h"
#undef private
#include "SequenceUtils/KmerIterator.h"
#include "SequenceUtils/Minimizer.h"
#include "SequenceUtils/SeqReader.h"
#include "ThreadedGzStream.h"
#include "Utilities/Database.h"

#include <algorithm>
#include <cstdio>
#include <ctime>
#include <string>
#include <vector>

using namespace protal;

namespace {
    double ThreadSeconds() {
        timespec ts;
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts);
        return static_cast<double>(ts.tv_sec) + 1e-9 * static_cast<double>(ts.tv_nsec);
    }

    constexpr size_t kMaxSeedSize = 128, kMinSuccessfulLookups = 4;  // protal's defaults

    struct Seeding {
        Seedmap& map;
        KmerLookupSM lookup;
        std::vector<LookupPointer> lookups;
        SeedList seeds;
        size_t seeds_total = 0, lookups_total = 0;

        Seeding(Seedmap& m) : map(m), lookup(m, 256) {}

        inline void PrefetchBlock(uint64_t kmer) {
            uint64_t const idx = map.ControlBlockIndex(map.MainKey(kmer));
            char const* p = reinterpret_cast<char const*>(map.m_keymap + idx);
            __builtin_prefetch(p);
            __builtin_prefetch(p + 31);  // the block and the next block's start: 32 bytes, maybe a second line
        }

        void Lookups(KmerList& kmers, std::vector<LookupPointer>& out) {
            out.clear();
            for (auto [mmer, pos] : kmers) lookup.Get(out, mmer, pos);
        }

        static void PrefetchValues(std::vector<LookupPointer> const& ls) {
            for (auto const& l : ls) {
                __builtin_prefetch(l.values_begin);
                if (l.flex_begin) __builtin_prefetch(l.flex_begin);
            }
        }

        // FindSeeds after the Get loop, on lookups made already.
        void FromLookups(std::vector<LookupPointer>& ls) {
            std::sort(ls.begin(), ls.end(), [](LookupPointer const& a, LookupPointer const& b) { return a.size < b.size; });
            seeds.clear();
            uint32_t previous_size = 0, successful = 0;
            for (size_t i = 0; i < ls.size(); i++) {
                lookup.GetFromLookup(seeds, ls[i]);
                successful += seeds.size() > previous_size;
                if (seeds.size() > kMaxSeedSize && successful > kMinSuccessfulLookups) break;
                previous_size = static_cast<uint32_t>(seeds.size());
            }
            seeds_total += seeds.size();
            lookups_total += ls.size();
        }

        void FindSeeds(KmerList& kmers) {
            Lookups(kmers, lookups);
            FromLookups(lookups);
        }

        // V5: pipelined within the read: the block of the k-mer D ahead in the Get loop, the values of the lookup E ahead
        // in the values loop.
        void FindSeedsPipelined(KmerList& kmers, size_t d, size_t e) {
            lookups.clear();
            size_t const n = kmers.size();
            for (size_t j = 0; j < std::min(d, n); j++) PrefetchBlock(kmers[j].first);
            for (size_t j = 0; j < n; j++) {
                if (j + d < n) PrefetchBlock(kmers[j + d].first);
                lookup.Get(lookups, kmers[j].first, kmers[j].second);
            }
            std::sort(lookups.begin(), lookups.end(), [](LookupPointer const& a, LookupPointer const& b) { return a.size < b.size; });
            for (size_t i = 0; i < std::min(e, lookups.size()); i++) { __builtin_prefetch(lookups[i].values_begin); if (lookups[i].flex_begin) __builtin_prefetch(lookups[i].flex_begin); }
            seeds.clear();
            uint32_t previous_size = 0, successful = 0;
            for (size_t i = 0; i < lookups.size(); i++) {
                if (i + e < lookups.size()) { __builtin_prefetch(lookups[i + e].values_begin); if (lookups[i + e].flex_begin) __builtin_prefetch(lookups[i + e].flex_begin); }
                lookup.GetFromLookup(seeds, lookups[i]);
                successful += seeds.size() > previous_size;
                if (seeds.size() > kMaxSeedSize && successful > kMinSuccessfulLookups) break;
                previous_size = static_cast<uint32_t>(seeds.size());
            }
            seeds_total += seeds.size();
            lookups_total += lookups.size();
        }
    };
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: pfbench DATABASE.protal READS.fq[.gz] [READS [ROUNDS]]\n");
        return 1;
    }
    size_t const max_reads = argc > 3 ? std::stoul(argv[3]) : 400000;
    int const rounds = argc > 4 ? std::stoi(argv[4]) : 1;
    std::string error;
    auto bundle = db::Bundle::Open(argv[1], error);
    if (!bundle) { std::fprintf(stderr, "%s\n", error.c_str()); return 1; }
    Seedmap map;
    double t0 = ThreadSeconds();
    map.Load(db::DbFile::InBundle(*bundle, "index.prx"), 6);
    std::fprintf(stderr, "index loaded (%.1f s thread CPU)\n", ThreadSeconds() - t0);

    ClosedSyncmer minimizer{15, 7, 2, map.UsesFullSyncmerMask()};
    SimpleKmerHandler kmer_handler{31, 15, minimizer};
    std::vector<KmerList> reads;
    {
        ThreadedGzIstream is(argv[2]);
        SeqReaderSE reader(is);
        FastxRecord record;
        KmerList kmers;
        while (reads.size() < max_reads && reader(record)) {
            kmers.clear();
            kmer_handler(std::string_view(record.sequence), kmers);
            reads.push_back(kmers);
        }
    }
    size_t kmers_total = 0;
    for (auto const& r : reads) kmers_total += r.size();
    std::fprintf(stderr, "%zu reads, %.1f k-mers each\n", reads.size(), static_cast<double>(kmers_total) / static_cast<double>(reads.size()));

    constexpr int kVariants = 7;
    char const* names[kVariants] = { "V0 as protal", "V1 own blocks first", "V2 blocks 2 ahead", "V3 2 stages (2,1)", "V4 2 stages (4,2)", "V5 in-read (4,4)", "V6 in-read (8,16)" };
    double seconds[kVariants] = {};
    size_t done[kVariants] = {}, seeds[kVariants] = {};
    Seeding s(map);
    std::vector<std::vector<LookupPointer>> ring(8);  // V3/V4: lookups made ahead, by read index % 8
    constexpr size_t kBatch = 1000;
    for (int round = 0; round < rounds; round++) {
        for (size_t b = 0; b * kBatch + 8 < reads.size(); b++) {
            int const v = static_cast<int>((b + static_cast<size_t>(round)) % kVariants);
            size_t const first = b * kBatch, last = std::min(reads.size() - 8, first + kBatch);
            size_t const seeds_before = s.seeds_total;
            // Variants that work ahead start their pipeline before the clock (as protal's loop would be in its stride).
            int const block_ahead = v == 2 ? 2 : v == 3 ? 2 : v == 4 ? 4 : 0;
            int const lookup_ahead = v == 3 ? 1 : v == 4 ? 2 : 0;
            if (v == 3 || v == 4) {
                for (int a = 0; a < block_ahead; a++) for (auto [m, p] : reads[first + a]) s.PrefetchBlock(m);
                for (int a = 0; a < lookup_ahead; a++) s.Lookups(reads[first + a], ring[(first + a) % 8]);
            } else if (v == 2) {
                for (int a = 0; a < block_ahead; a++) for (auto [m, p] : reads[first + a]) s.PrefetchBlock(m);
            }
            double const start = ThreadSeconds();
            for (size_t r = first; r < last; r++) {
                switch (v) {
                    case 0: s.FindSeeds(reads[r]); break;
                    case 1:
                        for (auto [m, p] : reads[r]) s.PrefetchBlock(m);
                        s.FindSeeds(reads[r]);
                        break;
                    case 2:
                        for (auto [m, p] : reads[r + 2]) s.PrefetchBlock(m);
                        s.FindSeeds(reads[r]);
                        break;
                    case 5: s.FindSeedsPipelined(reads[r], 4, 4); break;
                    case 6: s.FindSeedsPipelined(reads[r], 8, 16); break;
                    default: {
                        for (auto [m, p] : reads[r + block_ahead]) s.PrefetchBlock(m);
                        auto& next = ring[(r + lookup_ahead) % 8];
                        s.Lookups(reads[r + lookup_ahead], next);
                        Seeding::PrefetchValues(next);
                        s.FromLookups(ring[r % 8]);
                    }
                }
            }
            seconds[v] += ThreadSeconds() - start;
            done[v] += last - first;
            seeds[v] += s.seeds_total - seeds_before;
        }
    }
    std::printf("variant\treads\tns_per_read\tseeds_per_read\trelative\n");
    for (int v = 0; v < kVariants; v++) {
        std::printf("%s\t%zu\t%.0f\t%.2f\t%.3f\n", names[v], done[v], 1e9 * seconds[v] / static_cast<double>(done[v]),
                    static_cast<double>(seeds[v]) / static_cast<double>(done[v]),
                    (seconds[v] / static_cast<double>(done[v])) / (seconds[0] / static_cast<double>(done[0])));
    }
    return 0;
}
