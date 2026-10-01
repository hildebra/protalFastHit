// Unit tests for the syncmer extraction of reads and genes (SimpleKmerHandler<ClosedSyncmer>): the
// whole-sequence scan must give exactly the k-mers of the window-by-window definition, which the
// index was built with, for both s-mer masks (index formats 1 and 2), evaluating windows one by one
// and, where the CPU has it, 8 at a time with AVX2.
#include <gtest/gtest.h>
#include <random>
#include <string>
#include <vector>
#include "SequenceUtils/KmerIterator.h"

using namespace protal;

namespace {
    constexpr size_t kK = 31, kM = 15;

    // Each window on its own, from the definition: the forward and the reverse (complemented) k-mer,
    // their middle 15-mer cores, the canonical one (the reverse on a tie) tested by ClosedSyncmer.
    KmerList BruteForce(std::string const& seq, ClosedSyncmer& syncmer) {
        KmerList list;
        if (seq.size() < kK) return list;
        for (size_t p = 0; p + kK <= seq.size(); p++) {
            uint64_t fwd = 0, rev = 0;
            for (size_t i = 0; i < kK; i++) {
                fwd |= KmerUtils::BaseToInt(seq[p + i], 0) << (2 * (kK - 1 - i));
                rev |= KmerUtils::BaseToIntC(seq[p + i], 0) << (2 * i);
            }
            uint64_t const mask = (uint64_t{1} << (2 * kM)) - 1;
            uint64_t core_f = (fwd >> (kK - kM)) & mask, core_r = (rev >> (kK - kM)) & mask;
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
}

TEST(Syncmers, ScanGivesTheKmersOfTheDefinition) {
    auto const seqs = Sequences();
    for (bool avx2 : { false, true }) {
        for (bool full_mask : { true, false }) {
            ClosedSyncmer syncmer{kM, 7, 2, full_mask};
            SimpleKmerHandler<ClosedSyncmer> handler{kK, kM, syncmer};
            ASSERT_TRUE(handler.Scans());
            handler.UseAvx2(avx2);
            if (avx2 && !handler.UsesAvx2()) GTEST_SKIP() << "this CPU has no AVX2; the one-by-one evaluation passed";
            SimpleKmerHandler<ClosedSyncmer> copy{handler};  // each thread works on a copy
            ASSERT_TRUE(copy.Scans());
            ASSERT_EQ(copy.UsesAvx2(), avx2);
            KmerList scanned, windows;
            size_t total = 0;
            for (auto const& seq : seqs) {
                auto const expected = BruteForce(seq, syncmer);
                copy(std::string_view(seq), scanned);
                EXPECT_EQ(scanned, expected) << "AVX2 " << avx2 << ", full mask " << full_mask << ", sequence " << seq.substr(0, 60);
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

TEST(Syncmers, AVX2IsUsedWhereTheCpuHasIt) {
    ClosedSyncmer syncmer{kM, 7, 2, true};
    SimpleKmerHandler<ClosedSyncmer> handler{kK, kM, syncmer};
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
    EXPECT_EQ(handler.UsesAvx2(), static_cast<bool>(__builtin_cpu_supports("avx2")));
#endif
    handler.UseAvx2(false);
    EXPECT_FALSE(handler.UsesAvx2());
}

TEST(Syncmers, ScanReusesItsBuffersAcrossLengths) {
    ClosedSyncmer syncmer{kM, 7, 2, true};
    SimpleKmerHandler<ClosedSyncmer> handler{kK, kM, syncmer};
    auto const seqs = Sequences();
    KmerList list;
    for (int round = 0; round < 2; round++) {  // long, short, long: buffers shrink and grow
        handler.UseAvx2(round == 1);
        for (auto it = seqs.rbegin(); it != seqs.rend(); ++it) {
            handler(std::string_view(*it), list);
            EXPECT_EQ(list, BruteForce(*it, syncmer));
        }
    }
}

// The scan codes 32 characters at a time (AVX2) and the rest one by one: any byte, any length.
TEST(Syncmers, ScanCodesEveryByteAndEveryLengthAsTheDefinition) {
    std::mt19937 rng(23);
    std::vector<std::string> seqs;
    for (size_t len = kK; len < kK + 100; len++) {  // every length around the 32-character steps
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
    for (bool avx2 : { false, true }) {
        ClosedSyncmer syncmer{kM, 7, 2, true};
        SimpleKmerHandler<ClosedSyncmer> handler{kK, kM, syncmer};
        handler.UseAvx2(avx2);
        if (avx2 && !handler.UsesAvx2()) GTEST_SKIP() << "this CPU has no AVX2";
        KmerList scanned;
        for (auto const& seq : seqs) {
            handler(std::string_view(seq), scanned);
            ASSERT_EQ(scanned, BruteForce(seq, syncmer)) << "AVX2 " << avx2 << ", length " << seq.size();
        }
    }
}
