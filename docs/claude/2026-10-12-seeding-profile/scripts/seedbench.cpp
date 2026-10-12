// seedbench - experiment only: protal's seeding of paired-end reads (ChainAnchorFinder::LookupKmers / PrepareLookups /
// FindSeeds, KmerLookupSM::GetFromLookup) against a real database's index, the reads' k-mers taken as protal takes them,
// in protal's pipeline and in variants of it. Each variant makes the same lookups in the same order; the exact ones
// emit the same seeds (a hash over every seed, in order, is printed per variant). Variants alternate per batch of
// pairs on every thread, so that they see the same state of the machine; their time is the wall time
// inside the seeding (as protal's "Seeding" timer; the threads one to a core), not the gaps that stand for the rest of a pair's work.
//
//   seedbench DATABASE.protal R1.fq[.gz] R2.fq[.gz] [PAIRS [THREADS [ROUNDS [GAP1_NS [GAP2_NS]]]]]
//
// protal's paired-end loop, per pair p: PrepareLookups(p) (the lookups made, the key map lines prefetched 8 k-mers
// ahead, then the first line of every lookup's flex cells and entries prefetched); the next pair read and its k-mers
// taken, PrefetchKeys(p + 1) (GAP1: ~1.5 us at r226); FindSeeds(p) (smallest block first, the values prefetched 16
// lookups ahead, GetFromLookup until enough seeds); then the rest of the pair (GAP2: ~6 us at r226). The gaps are
// spins that touch no memory.
//
//   P0 as protal before 2026-10-12 (the first line of each block's flex cells and entries prefetched)
//   P1 the lookups made one pair earlier: values prefetched a whole pair ahead, keys two pairs ahead      (exact)
//   P2 as protal, every line of a block prefetched (its flex cells and entries, up to kWholeBlockLines)   (exact)
//   P3 P1 and P2                                                                                          (exact)
//   P4 four pairs' lookups made and prefetched together, then their seeds (more misses in flight)          (exact)
//   B0 the key map only: the lookups made, nothing of the values read                                    (part)
//   B1 B0 and every block scanned (the flex cells scored), no entries read, no seeds                     (part)
//   B2 as protal, but each pair seeded twice and the second time measured: the caches warm (compute)    (part)
#include <bits/stdc++.h>  // protal's headers take the standard headers as included before them
#include <atomic>
#include <chrono>
#include <cinttypes>
#include <cstdio>
#include <ctime>
#include <memory>
#include <string>
#include <vector>
#include <omp.h>

#include "Hash/Seedmap.h"
#include "Hash/KmerLookup.h"
#include "Hash/FlexScan.h"
#include "SequenceUtils/KmerIterator.h"
#include "SequenceUtils/Minimizer.h"
#include "SequenceUtils/SeqReader.h"
#include "SequenceUtils/GenomeLoader.h"
#include "ThreadedGzStream.h"
#include "Utilities/Database.h"

using namespace protal;

namespace {
    // Wall time from the monotonic clock (vDSO, ~20 ns a call): the thread CPU time is a system call, ~0.1-0.2 us, four
    // times a pair. Threads pinned one to a core (taskset, or a whole node) see little else.
    double ThreadSeconds() {
        timespec ts;
        clock_gettime(CLOCK_MONOTONIC, &ts);
        return static_cast<double>(ts.tv_sec) + 1e-9 * static_cast<double>(ts.tv_nsec);
    }

    // Spin about `ns` nanoseconds of wall time without touching memory (a stand-in for the rest of a pair's work).
    void Spin(uint64_t ns) {
        if (ns == 0) return;
        auto const end = std::chrono::steady_clock::now() + std::chrono::nanoseconds(ns);
        uint64_t x = 1;
        while (std::chrono::steady_clock::now() < end) {
            for (int i = 0; i < 64; i++) x = x * 6364136223846793005ull + 1442695040888963407ull;
        }
        asm volatile("" : : "r"(x));
    }

    constexpr size_t kMaxSeedSize = 128, kMinSuccessfulLookups = 4, kMaxUbiquity = 256;  // protal's defaults
    constexpr size_t kPrefetchKmers = 8, kPrefetchLookups = 16;                          // ChainAnchorFinder's
    constexpr size_t kWholeBlockLines = 32;

    struct Counts {
        uint64_t pairs = 0, kmers = 0, lookups = 0, cells = 0, seeds = 0, hash = 1469598103934665603ull;
        double seconds = 0;
    };

    // One mate's seeding state, as a ChainAnchorFinder holds it.
    struct Mate {
        std::vector<LookupPointer> lookups;
        SeedList seeds;
    };

    struct Seeder {
        Seedmap& map;
        KmerLookupSM lookup;
        uint64_t entry_bits;
        bool whole_blocks = false;

        Seeder(Seedmap& m, uint64_t bits) : map(m), lookup(m, kMaxUbiquity), entry_bits(bits) {}

        void PrefetchKeys(KmerList const& kmers) {
            for (auto const& [kmer, pos] : kmers) lookup.PrefetchKey(kmer);
        }

        void PrefetchValues(LookupPointer const& l) const {
            if (!whole_blocks) {
                // protal before 2026-10-12: the first line of the flex cells and of the entries (P0, B1); since then
                // KmerLookupSM::PrefetchValues fetches whole blocks, as P2 does.
                if (l.flex != nullptr) __builtin_prefetch(l.flex);
                __builtin_prefetch(l.entries);
                return;
            }
            // The block's flex cells and entries are one stretch: from the flex cells' start (or the entries') to the
            // entries' end, at entry_bits each.
            uint8_t const* begin = l.flex ? l.flex : l.entries;
            uint64_t const entry_bytes = (static_cast<uint64_t>(l.size) * entry_bits + l.entry_shift + 7) / 8;
            uint8_t const* end = l.entries + entry_bytes;
            size_t lines = 0;
            for (uint8_t const* p = begin; p < end + 63 && lines < kWholeBlockLines; p += 64, lines++) __builtin_prefetch(p);
        }

        // ChainAnchorFinder::LookupKmers.
        void LookupKmers(KmerList const& kmers, std::vector<LookupPointer>& out) {
            out.clear();
            size_t const n = kmers.size();
            for (size_t j = 0; j < std::min(kPrefetchKmers, n); j++) lookup.PrefetchKey(kmers[j].first);
            for (size_t j = 0; j < n; j++) {
                if (j + kPrefetchKmers < n) lookup.PrefetchKey(kmers[j + kPrefetchKmers].first);
                size_t kmer = kmers[j].first;
                lookup.Get(out, kmer, static_cast<uint32_t>(kmers[j].second));
            }
        }

        // ChainAnchorFinder::PrepareLookups.
        void Prepare(KmerList const& kmers, Mate& mate) {
            LookupKmers(kmers, mate.lookups);
            for (auto const& l : mate.lookups) PrefetchValues(l);
        }

        // ChainAnchorFinder::FindSeeds on prepared lookups.
        void FindSeeds(Mate& mate, Counts& c) {
            auto& lookups = mate.lookups;
            mate.seeds.clear();
            std::sort(lookups.begin(), lookups.end(), [](LookupPointer const& a, LookupPointer const& b) { return a.size < b.size; });
            for (size_t i = 0; i < std::min(kPrefetchLookups, lookups.size()); i++) PrefetchValues(lookups[i]);
            uint32_t previous_size = 0, successful = 0;
            for (size_t i = 0; i < lookups.size(); i++) {
                if (i + kPrefetchLookups < lookups.size()) PrefetchValues(lookups[i + kPrefetchLookups]);
                c.lookups++;
                c.cells += lookups[i].size;
                lookup.GetFromLookup(mate.seeds, lookups[i]);
                successful += mate.seeds.size() > previous_size;
                if (mate.seeds.size() > kMaxSeedSize && successful > kMinSuccessfulLookups) break;
                previous_size = static_cast<uint32_t>(mate.seeds.size());
            }
            c.seeds += mate.seeds.size();
            for (auto const& s : mate.seeds) {
                uint64_t const v[3] = { (uint64_t{ s.taxid } << 32) | s.geneid, (uint64_t{ s.genepos } << 16) | s.readpos,
                                        (uint64_t{ s.unique } << 2) | (uint64_t{ s.unique_dist_two } << 1) | s.single };
                for (uint64_t x : v) c.hash = (c.hash ^ x) * 1099511628211ull;
            }
        }

        // B1: every block scanned (the flex cells scored with the kernel GetFromLookup uses), no entries, no seeds; the best
        // score of each block counted by its size class (best[class][score]: lookups, cells).
        std::vector<uint8_t> scores;
        std::vector<uint32_t> masks;
        std::array<std::array<std::array<uint64_t, 2>, 17>, 4> best{};
        static size_t SizeClass(uint32_t size) { return size < 16 ? 0 : size < 256 ? 1 : size < 4096 ? 2 : 3; }
        void ScanOnly(Mate& mate, Counts& c) {
            auto& lookups = mate.lookups;
            std::sort(lookups.begin(), lookups.end(), [](LookupPointer const& a, LookupPointer const& b) { return a.size < b.size; });
            uint32_t sink = 0;
            for (size_t i = 0; i < lookups.size(); i++) {
                if (i + kPrefetchLookups < lookups.size()) PrefetchValues(lookups[i + kPrefetchLookups]);
                auto const& l = lookups[i];
                c.lookups++;
                c.cells += l.size;
                if (!l.flex) continue;
                size_t const bytes = std::max<size_t>(l.size + flex_scan::kScorePadding, flex_scan::TieScoreBytes(l.size)) + 64;
                if (scores.size() < bytes) scores.resize(bytes);
                if (masks.size() < flex_scan::TieWords(l.size) + 2) masks.resize(flex_scan::TieWords(l.size) + 2);
                simd::Level const level = flex_scan::Kernel().load(std::memory_order_relaxed);
                uint32_t b = 0;
                if (level == simd::Level::avx512) b = flex_scan::ScanAvx512(l, l.flex_key, scores.data(), masks.data()).first;
                else if (level == simd::Level::avx2) b = flex_scan::BestAvx2(l, l.flex_key, scores.data());
                else b = flex_scan::ScoreScalar(l, l.flex_key, scores.data()).first;
                sink += b;
                auto& cell = best[SizeClass(l.size)][std::min<uint32_t>(b, 16)];
                cell[0]++;
                cell[1] += l.size;
            }
            c.hash += sink;
        }
    };

    char const* kNames[] = { "P0 as protal", "P1 values a pair ahead", "P2 whole blocks prefetched", "P3 P1+P2",
                             "P4 four pairs together", "B0 key map only", "B1 key map + scans", "B2 warm caches (compute)" };
    constexpr int kVariants = 8;
    constexpr bool kExact[kVariants] = { true, true, true, true, true, false, false, false };
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::fprintf(stderr, "usage: seedbench DATABASE.protal R1.fq[.gz] R2.fq[.gz] [PAIRS [THREADS [ROUNDS [GAP1_NS [GAP2_NS]]]]]\n");
        return 1;
    }
    size_t const max_pairs = argc > 4 ? std::stoul(argv[4]) : 500000;
    int const threads = argc > 5 ? std::stoi(argv[5]) : 8;
    int const rounds = argc > 6 ? std::stoi(argv[6]) : 3;
    uint64_t const gap1 = argc > 7 ? std::stoull(argv[7]) : 1500;
    uint64_t const gap2 = argc > 8 ? std::stoull(argv[8]) : 6000;

    std::string error;
    auto bundle = db::Bundle::Open(argv[1], error);
    if (!bundle) { std::fprintf(stderr, "%s\n", error.c_str()); return 1; }
    auto const t_start = std::chrono::steady_clock::now();
    auto gene_table = db::DbFile::InBundle(*bundle, gene_table_file::kFileName);
    GenomeLoader genomes(db::DbFile::InBundle(*bundle, "reference.fna"), db::DbFile::InBundle(*bundle, "reference.map"), threads,
                         gene_table.Exists() ? std::optional<db::DbFile>(gene_table) : std::nullopt,
                         db::DbFile::InBundle(*bundle, "unique_kmers.tsv").Size());
    auto const [max_taxid, max_gene, max_length] = genomes.IndexFieldMaxima();
    Seedmap map;
    Seedmap::PackedLayout const layout = Seedmap::PackedLayout::For(max_taxid, max_gene, max_length + map.m_flex_k);
    map.Load(db::DbFile::InBundle(*bundle, "index.prx"), threads, &layout);
    std::fprintf(stderr, "index loaded in %.1f s: %s\n",
                 std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count(), map.MemoryDescription().c_str());

    // The pairs' k-mers, as protal takes them (RunProtal: ClosedSyncmer{15, 7, 2}, SimpleKmerHandler{31, 15}).
    ClosedSyncmer minimizer{ 15, 7, 2, map.UsesFullSyncmerMask() };
    SimpleKmerHandler<ClosedSyncmer> kmer_handler{ 31, 15, minimizer };
    std::vector<KmerList> mates1, mates2;
    {
        ThreadedGzIstream is1(argv[2]), is2(argv[3]);
        SeqReaderSE reader1(is1), reader2(is2);
        FastxRecord r1, r2;
        KmerList k;
        while (mates1.size() < max_pairs && reader1(r1) && reader2(r2)) {
            k.clear();
            kmer_handler(std::string_view(r1.sequence), k);
            mates1.push_back(k);
            k.clear();
            kmer_handler(std::string_view(r2.sequence), k);
            mates2.push_back(k);
        }
    }
    size_t const n = mates1.size();
    size_t kmers = 0;
    for (size_t i = 0; i < n; i++) kmers += mates1[i].size() + mates2[i].size();
    std::fprintf(stderr, "%zu pairs, %.1f k-mers a mate; %d threads, %d rounds, gaps %" PRIu64 " and %" PRIu64 " ns; flex scan %s\n", n,
                 static_cast<double>(kmers) / (2.0 * static_cast<double>(n)), threads, rounds, gap1, gap2,
                 flex_scan::Kernel().load() == simd::Level::avx512 ? "AVX-512" : flex_scan::Kernel().load() == simd::Level::avx2 ? "AVX2" : "scalar");

    constexpr size_t kBatch = 500;
    std::vector<std::array<Counts, kVariants>> per_thread(static_cast<size_t>(threads));
    std::array<std::array<std::array<uint64_t, 2>, 17>, 4> best_scores{};  // B1's blocks by size class and best score
#pragma omp parallel num_threads(threads)
    {
        int const t = omp_get_thread_num();
        size_t const begin = n * static_cast<size_t>(t) / static_cast<size_t>(threads), end = n * static_cast<size_t>(t + 1) / static_cast<size_t>(threads);
        Seeder s(map, layout.EntryBits());
        std::array<Mate, 2> cur, next;            // P0-P3: this pair's mates and the next pair's (P1, P3)
        std::array<std::array<Mate, 2>, 4> four;  // P4
        auto& counts = per_thread[static_cast<size_t>(t)];
        for (int round = 0; round < rounds; round++) {
            for (size_t b = begin, batch = 0; b + 4 < end; b += kBatch, batch++) {
                int const v = static_cast<int>((batch + static_cast<size_t>(round) + static_cast<size_t>(t)) % kVariants);
                size_t const last = std::min(end - 4, b + kBatch);
                Counts& c = counts[static_cast<size_t>(v)];
                s.whole_blocks = v == 2 || v == 3;
                double busy = 0;
                if (v == 0 || v == 2 || v == 7) {
                    s.PrefetchKeys(mates1[b]);
                    s.PrefetchKeys(mates2[b]);
                    for (size_t p = b; p < last; p++) {
                        double t0 = ThreadSeconds();
                        s.Prepare(mates1[p], cur[0]);
                        s.Prepare(mates2[p], cur[1]);
                        if (v != 7) busy += ThreadSeconds() - t0;
                        Spin(gap1);
                        t0 = ThreadSeconds();
                        s.PrefetchKeys(mates1[p + 1]);
                        s.PrefetchKeys(mates2[p + 1]);
                        if (v == 7) {  // the pair seeded once untimed, then again (its lookups made again) timed
                            Counts scratch;
                            s.FindSeeds(cur[0], scratch);
                            s.FindSeeds(cur[1], scratch);
                            t0 = ThreadSeconds();
                            s.Prepare(mates1[p], cur[0]);
                            s.Prepare(mates2[p], cur[1]);
                        }
                        s.FindSeeds(cur[0], c);
                        s.FindSeeds(cur[1], c);
                        busy += ThreadSeconds() - t0;
                        Spin(gap2);
                        c.pairs++;
                        c.kmers += mates1[p].size() + mates2[p].size();
                    }
                } else if (v == 1 || v == 3) {
                    // The lookups of pair p + 1 made (and their values prefetched) while pair p is seeded and the rest
                    // of it runs; the keys of pair p + 2 prefetched then.
                    s.PrefetchKeys(mates1[b]);
                    s.PrefetchKeys(mates2[b]);
                    s.Prepare(mates1[b], next[0]);
                    s.Prepare(mates2[b], next[1]);
                    s.PrefetchKeys(mates1[b + 1]);
                    s.PrefetchKeys(mates2[b + 1]);
                    for (size_t p = b; p < last; p++) {
                        std::swap(cur, next);
                        double t0 = ThreadSeconds();
                        s.Prepare(mates1[p + 1], next[0]);
                        s.Prepare(mates2[p + 1], next[1]);
                        busy += ThreadSeconds() - t0;
                        Spin(gap1);
                        t0 = ThreadSeconds();
                        s.PrefetchKeys(mates1[p + 2]);
                        s.PrefetchKeys(mates2[p + 2]);
                        s.FindSeeds(cur[0], c);
                        s.FindSeeds(cur[1], c);
                        busy += ThreadSeconds() - t0;
                        Spin(gap2);
                        c.pairs++;
                        c.kmers += mates1[p].size() + mates2[p].size();
                    }
                } else if (v == 4) {
                    // Four pairs at a time: all their lookups made, all values prefetched, then each pair's seeds;
                    // the next four's keys prefetched meanwhile. The gaps as for four pairs, in the same places.
                    for (size_t q = 0; q < 4; q++) { s.PrefetchKeys(mates1[b + q]); s.PrefetchKeys(mates2[b + q]); }
                    for (size_t p = b; p + 4 <= last; p += 4) {
                        double t0 = ThreadSeconds();
                        for (size_t q = 0; q < 4; q++) { s.Prepare(mates1[p + q], four[q][0]); s.Prepare(mates2[p + q], four[q][1]); }
                        busy += ThreadSeconds() - t0;
                        Spin(4 * gap1);
                        t0 = ThreadSeconds();
                        for (size_t q = 0; q < 4 && p + 4 + q < end; q++) { s.PrefetchKeys(mates1[p + 4 + q]); s.PrefetchKeys(mates2[p + 4 + q]); }
                        for (size_t q = 0; q < 4; q++) {
                            s.FindSeeds(four[q][0], c);
                            s.FindSeeds(four[q][1], c);
                            if (q < 3) { busy += ThreadSeconds() - t0; Spin(gap2); t0 = ThreadSeconds(); }
                        }
                        busy += ThreadSeconds() - t0;
                        Spin(gap2);
                        c.pairs += 4;
                        for (size_t q = 0; q < 4; q++) c.kmers += mates1[p + q].size() + mates2[p + q].size();
                    }
                } else {
                    // B0: the lookups only; B1: and the scans.
                    s.PrefetchKeys(mates1[b]);
                    s.PrefetchKeys(mates2[b]);
                    for (size_t p = b; p < last; p++) {
                        double t0 = ThreadSeconds();
                        s.LookupKmers(mates1[p], cur[0].lookups);
                        s.LookupKmers(mates2[p], cur[1].lookups);
                        if (v == 6) for (auto& m : cur) for (auto const& l : m.lookups) s.PrefetchValues(l);
                        busy += ThreadSeconds() - t0;
                        Spin(gap1);
                        t0 = ThreadSeconds();
                        s.PrefetchKeys(mates1[p + 1]);
                        s.PrefetchKeys(mates2[p + 1]);
                        if (v == 6) {
                            s.ScanOnly(cur[0], c);
                            s.ScanOnly(cur[1], c);
                        } else {
                            c.lookups += cur[0].lookups.size() + cur[1].lookups.size();
                            for (auto& m : cur) for (auto const& l : m.lookups) c.cells += l.size;
                        }
                        busy += ThreadSeconds() - t0;
                        Spin(gap2);
                        c.pairs++;
                        c.kmers += mates1[p].size() + mates2[p].size();
                    }
                }
                c.seconds += busy;
            }
        }
#pragma omp critical(best_scores)
        for (size_t k = 0; k < best_scores.size(); k++) {
            for (size_t b = 0; b < 17; b++) {
                best_scores[k][b][0] += s.best[k][b][0];
                best_scores[k][b][1] += s.best[k][b][1];
            }
        }
    }
    std::array<Counts, kVariants> total{};
    for (auto const& tc : per_thread) {
        for (int v = 0; v < kVariants; v++) {
            auto& c = total[static_cast<size_t>(v)];
            c.pairs += tc[static_cast<size_t>(v)].pairs;
            c.kmers += tc[static_cast<size_t>(v)].kmers;
            c.lookups += tc[static_cast<size_t>(v)].lookups;
            c.cells += tc[static_cast<size_t>(v)].cells;
            c.seeds += tc[static_cast<size_t>(v)].seeds;
            c.seconds += tc[static_cast<size_t>(v)].seconds;
        }
    }
    // The exact variants' seeds against P0's, pair for pair: seeded again on one thread, untimed, and hashed in order.
    {
        Seeder s(map, layout.EntryBits());
        std::array<uint64_t, kVariants> hashes{};
        for (int v = 0; v < kVariants; v++) {
            if (!kExact[v]) continue;
            s.whole_blocks = v == 2 || v == 3;
            Counts c;
            Mate m;
            size_t const check = std::min<size_t>(n, 20000);
            for (size_t p = 0; p < check; p++) {
                s.Prepare(mates1[p], m); s.FindSeeds(m, c);
                s.Prepare(mates2[p], m); s.FindSeeds(m, c);
            }
            hashes[static_cast<size_t>(v)] = c.hash;
        }
        std::printf("variant\tpairs\tkmers\tlookups\tcells_per_lookup\tseeds_per_lookup\tns_per_kmer\tns_per_lookup\tus_per_pair\trelative\tseeds_hash\n");
        double const base = total[0].lookups ? total[0].seconds / static_cast<double>(total[0].lookups) : 0;
        for (int v = 0; v < kVariants; v++) {
            auto const& c = total[static_cast<size_t>(v)];
            double const per = c.lookups ? c.seconds / static_cast<double>(c.lookups) : 0;
            std::printf("%s\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%.1f\t%.2f\t%.1f\t%.1f\t%.3f\t%.3f\t%s\n", kNames[v], c.pairs, c.kmers, c.lookups,
                        c.lookups ? static_cast<double>(c.cells) / static_cast<double>(c.lookups) : 0.0,
                        c.lookups ? static_cast<double>(c.seeds) / static_cast<double>(c.lookups) : 0.0,
                        c.kmers ? 1e9 * c.seconds / static_cast<double>(c.kmers) : 0.0, 1e9 * per,
                        c.pairs ? 1e6 * c.seconds / static_cast<double>(c.pairs) : 0.0, base > 0 ? per / base : 0.0,
                        kExact[v] ? (hashes[static_cast<size_t>(v)] == hashes[0] ? "same" : "DIFFERENT") : "-");
        }
    }
    // B1's blocks with flex cells, by size class: the share of lookups (and of cells) whose best score is each value.
    {
        char const* classes[4] = { "2-15", "16-255", "256-4095", ">=4096" };
        std::printf("\nblock_size\tlookups\tcells\tshare_of_lookups_by_best_score(16 15 14 13 12 <=11)\tshare_of_cells_by_best_score(16 15 14 13 12 <=11)\n");
        for (size_t k = 0; k < 4; k++) {
            uint64_t lookups = 0, cells = 0;
            for (size_t b = 0; b < 17; b++) { lookups += best_scores[k][b][0]; cells += best_scores[k][b][1]; }
            if (lookups == 0) continue;
            std::printf("%s\t%" PRIu64 "\t%" PRIu64 "\t", classes[k], lookups, cells);
            for (size_t which = 0; which < 2; which++) {
                uint64_t const whole = which ? cells : lookups;
                uint64_t low = 0;
                for (size_t b = 0; b <= 11; b++) low += best_scores[k][b][which];
                for (int b = 16; b >= 12; b--) std::printf("%.3f ", static_cast<double>(best_scores[k][static_cast<size_t>(b)][which]) / static_cast<double>(whole));
                std::printf("%.3f%s", static_cast<double>(low) / static_cast<double>(whole), which ? "\n" : "\t");
            }
        }
    }
    return 0;
}
