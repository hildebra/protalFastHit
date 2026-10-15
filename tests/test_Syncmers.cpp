// Unit tests for the syncmer extraction of reads and genes (SimpleKmerHandler<ClosedSyncmer>): the
// whole-sequence scan must give exactly the k-mers of the window-by-window definition, which the
// index was built with, for both s-mer masks (index formats 1 and 2), evaluating windows one by one
// and, where the CPU has them, 8 at a time with AVX2 and 16 at a time with AVX-512, for several k-mer, core and
// s-mer lengths; and ClosedSyncmer's selection of one core.
#include <gtest/gtest.h>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include "SequenceUtils/KmerIterator.h"

using namespace protal;

namespace {
    constexpr size_t kK = 31, kM = 15;

    // Each window on its own, from the definition: the forward and the reverse (complemented) k-mer,
    // their middle 15-mer cores, the canonical one (the reverse on a tie) tested by ClosedSyncmer.
    KmerList BruteForce(std::string const& seq, ClosedSyncmer& syncmer, size_t k = kK, size_t m = kM) {
        KmerList list;
        if (seq.size() < k) return list;
        for (size_t p = 0; p + k <= seq.size(); p++) {
            uint64_t fwd = 0, rev = 0;
            for (size_t i = 0; i < k; i++) {
                fwd |= KmerUtils::BaseToInt(seq[p + i], 0) << (2 * (k - 1 - i));
                rev |= KmerUtils::BaseToIntC(seq[p + i], 0) << (2 * i);
            }
            uint64_t const mask = (uint64_t{1} << (2 * m)) - 1;
            uint64_t core_f = (fwd >> (k - m)) & mask, core_r = (rev >> (k - m)) & mask;
            uint64_t core = core_f < core_r ? core_f : core_r;
            if (syncmer(core)) list.emplace_back(core_f < core_r ? fwd : rev, p);
        }
        return list;
    }

    std::vector<std::string> Sequences() {
        std::mt19937 rng(17);
        std::vector<std::string> seqs = { "", "ACGT", std::string(30, 'A'), std::string(31, 'A'), std::string(32, 'C'),
                                          std::string(200, 'A'), std::string(200, 'N') };
        for (std::string unit : { "AC", "AAT", "ACGT", "AAAACCCC", "ATATATATGC" }) {  // repeats: ties among s-mers
            std::string s;
            while (s.size() < 300) s += unit;
            seqs.push_back(s);
        }
        for (size_t len : { 31u, 32u, 45u, 100u, 150u, 151u, 250u, 1000u, 5000u }) {
            for (int copy = 0; copy < 20; copy++) {
                std::string s(len, 'A');
                for (auto& c : s) c = "ACGT"[rng() % 4];
                if (copy % 4 == 1) for (size_t i = rng() % 37; i < len; i += 37) s[i] = 'N';                // Ns
                if (copy % 4 == 2) for (size_t i = 0; i < len; i += 5) s[i] = static_cast<char>(s[i] + 32);  // lower case
                if (copy % 4 == 3) s[rng() % len] = "RYKM-."[rng() % 6];                                      // other symbols
                seqs.push_back(s);
                seqs.push_back(KmerUtils::ReverseComplement(s));
            }
        }
        return seqs;
    }

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

// The selection of a core itself, for the s-mer mask of each index format.
TEST(ClosedSyncmer, LegacyAndFullMasksMatchTheirDefinition) {
    ClosedSyncmer legacy{15, 7, 2, false};
    ClosedSyncmer full{15, 7, 2, true};
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

// The levels the scan is tested at here: scalar, and AVX2 and AVX-512 where the CPU has them (a note for the others).
static std::vector<simd::Level> Levels() {
    std::vector<simd::Level> levels{ simd::Level::scalar };
    for (simd::Level level : { simd::Level::avx2, simd::Level::avx512 }) {
        if (simd::Supported(level) == level) levels.push_back(level);
        else std::cout << "no " << simd::Name(level) << " here: not tested" << std::endl;
    }
    return levels;
}

TEST(Syncmers, ScanGivesTheKmersOfTheDefinition) {
    auto const seqs = Sequences();
    for (simd::Level level : Levels()) {
        for (bool full_mask : { true, false }) {
            ClosedSyncmer syncmer{kM, 7, 2, full_mask};
            SimpleKmerHandler<ClosedSyncmer> handler{kK, kM, syncmer};
            ASSERT_TRUE(handler.Scans());
            ASSERT_EQ(handler.Use(level), level);
            SimpleKmerHandler<ClosedSyncmer> copy{handler};  // each thread works on a copy
            ASSERT_TRUE(copy.Scans());
            ASSERT_EQ(copy.Level(), level);
            KmerList scanned, windows;
            size_t total = 0;
            for (auto const& seq : seqs) {
                auto const expected = BruteForce(seq, syncmer);
                copy(std::string_view(seq), scanned);
                EXPECT_EQ(scanned, expected) << simd::Name(level) << ", full mask " << full_mask << ", sequence " << seq.substr(0, 60);
                EXPECT_EQ(copy.TotalKmers(), seq.size() >= kK ? seq.size() - kK + 1 : 0);
                EXPECT_EQ(copy.TotalMinimizers(), expected.size());
                copy.WindowByWindow(std::string_view(seq), windows);
                EXPECT_EQ(windows, expected);
                total += expected.size();
            }
            EXPECT_GT(total, 10000u);
        }
    }
}

// Other k-mer, core and s-mer lengths the scan applies to (the core's place in the k-mer, the s-mers per core and t
// change): the k-mers of the definition at every level.
TEST(Syncmers, ScanGivesTheDefinitionsKmersForOtherShapes) {
    struct Shape { size_t k, m; uint8_t s, t; };
    auto const seqs = Sequences();
    size_t shapes = 0;
    for (Shape shape : { Shape{ 31, 13, 5, 1 }, Shape{ 25, 13, 6, 2 }, Shape{ 21, 11, 5, 1 }, Shape{ 31, 11, 4, 0 },
                         Shape{ 15, 15, 7, 2 }, Shape{ 29, 15, 3, 4 }, Shape{ 17, 9, 2, 3 }, Shape{ 31, 15, 14, 0 } }) {
        for (simd::Level level : Levels()) {
            ClosedSyncmer syncmer{static_cast<uint8_t>(shape.m), shape.s, shape.t, true};
            SimpleKmerHandler<ClosedSyncmer> handler{shape.k, shape.m, syncmer};
            ASSERT_TRUE(handler.Scans()) << "k " << shape.k << ", m " << shape.m << ", s " << int{shape.s};
            ASSERT_EQ(handler.Use(level), level);
            KmerList scanned;
            for (auto const& seq : seqs) {
                handler(std::string_view(seq), scanned);
                ASSERT_EQ(scanned, BruteForce(seq, syncmer, shape.k, shape.m))
                    << simd::Name(level) << ", k " << shape.k << ", m " << shape.m << ", s " << int{shape.s} << ", length " << seq.size();
            }
        }
        shapes++;
    }
    EXPECT_EQ(shapes, 8u);
}

TEST(Syncmers, TheLevelIsTheCpusUnlessCapped) {
    ClosedSyncmer syncmer{kM, 7, 2, true};
    SimpleKmerHandler<ClosedSyncmer> handler{kK, kM, syncmer};
    EXPECT_EQ(handler.Level(), simd::Default());  // the CPU's, capped by PROTAL_SIMD
    EXPECT_EQ(handler.Use(simd::Level::avx512), simd::CpuLevel());
    EXPECT_EQ(handler.Use(simd::Level::scalar), simd::Level::scalar);
    EXPECT_EQ(handler.Level(), simd::Level::scalar);
}

// The scan codes 64 (AVX-512) or 32 (AVX2) characters at a time, the rest masked or one by one: any byte, any length.
TEST(Syncmers, ScanCodesEveryByteAndEveryLengthAsTheDefinition) {
    std::mt19937 rng(23);
    std::vector<std::string> seqs;
    for (size_t len = kK; len < kK + 160; len++) {  // every length around the 32- and 64-character steps
        std::string s(len, 'A');
        for (auto& c : s) c = "ACGTACGTACGTacgtN"[rng() % 17];
        seqs.push_back(s);
    }
    for (int copy = 0; copy < 200; copy++) {        // bases with a sprinkle of every other byte
        std::string s(150, 'A');
        for (auto& c : s) c = "ACGT"[rng() % 4];
        for (int i = 0; i < 6; i++) s[rng() % s.size()] = static_cast<char>(1 + rng() % 255);
        seqs.push_back(s);
    }
    std::string every(256 + kK, 'A');               // each byte value once, in a k-mer of bases
    for (int b = 0; b < 256; b++) every[kK / 2 + b] = static_cast<char>(b);
    seqs.push_back(every);
    for (simd::Level level : Levels()) {
        ClosedSyncmer syncmer{kM, 7, 2, true};
        SimpleKmerHandler<ClosedSyncmer> handler{kK, kM, syncmer};
        ASSERT_EQ(handler.Use(level), level);
        KmerList scanned;
        for (auto const& seq : seqs) {
            handler(std::string_view(seq), scanned);
            ASSERT_EQ(scanned, BruteForce(seq, syncmer)) << simd::Name(level) << ", length " << seq.size();
        }
    }
}

// A read's k-mers over an ambiguous base are not looked up (DropReadAmbiguousKmers; F1 audit L10, they were looked up as
// if the base were A or T): the list of a read with an N at p is that of the same read with A there, less the k-mers
// whose window covers p (every other window is the same in both), at every level; a read of A, C, G, T only, lowercase
// too, keeps its list. (The list's second is the window's first position: ScanGivesTheKmersOfTheDefinition.)
TEST(Syncmers, ReadKmersOverAnAmbiguousBaseAreLeftOut) {
    std::mt19937 rng(23);
    for (simd::Level level : Levels()) {
        ClosedSyncmer syncmer{kM, 7, 2, true};
        SimpleKmerHandler<ClosedSyncmer> handler{kK, kM, syncmer};
        ASSERT_EQ(handler.Use(level), level);
        ASSERT_EQ(handler.K(), kK);
        size_t dropped = 0;
        for (int round = 0; round < 60; round++) {
            std::string read(150, 'A');
            for (auto& c : read) c = "ACGT"[rng() % 4];
            KmerList clean;
            handler(std::string_view(read), clean);
            KmerList kept = clean;
            DropReadAmbiguousKmers(read, handler, kept);
            EXPECT_EQ(kept, clean);
            std::string lower = read;
            for (auto& c : lower) c = static_cast<char>(c | 0x20);
            KmerList lower_list;
            handler(std::string_view(lower), lower_list);
            kept = lower_list;
            DropReadAmbiguousKmers(lower, handler, kept);
            EXPECT_EQ(kept, lower_list);

            size_t const p = rng() % read.size();
            std::string n_read = read, a_read = read;
            n_read[p] = round % 3 == 0 ? 'n' : round % 3 == 1 ? 'N' : 'R';
            a_read[p] = 'A';
            KmerList with_n, a_list, expected;
            handler(std::string_view(n_read), with_n);
            handler(std::string_view(a_read), a_list);
            DropReadAmbiguousKmers(n_read, handler, with_n);
            for (auto const& kmer : a_list) {
                if (!(kmer.second <= p && p < kmer.second + kK)) expected.push_back(kmer);
            }
            EXPECT_EQ(with_n, expected) << simd::Name(level) << ", N at " << p;
            dropped += a_list.size() - expected.size();
        }
        EXPECT_GT(dropped, 0u);
    }
}
