// Unit tests for the k-mer screen before WFA2 (AlignmentScreen.h): a candidate it refuses is one WFA2
// fails within the same budget, with free ends and Ns; and through SimpleAlignmentHandler the same
// alignments with the screen on and off.
#include <gtest/gtest.h>
#include <algorithm>
#include <array>
#include <filesystem>
#include <fstream>
#include <memory>
#include <random>
#include <string>
#include <vector>
#include <unistd.h>
#include "Alignment/AlignmentScreen.h"
#include "Core/AlignmentStrategy.h"
#include "TestReference.h"

using namespace protal;

namespace {
    constexpr int kMismatch = 4, kGapOpen = 6, kGapExtend = 2;

    std::string RandomBases(std::mt19937& rng, size_t n) {
        std::string s(n, 'A');
        for (auto& c : s) c = "ACGT"[rng() % 4];
        return s;
    }

    // A read of `length` bases copied from `gene` at `start` (random bases beyond the gene) with edits at the rate
    // `divergence` (90% mismatches, 5% deletions, 5% insertions) and Ns at `n_rate`.
    std::string Derive(std::mt19937& rng, std::string const& gene, long start, size_t length, double divergence, double n_rate) {
        std::uniform_real_distribution<double> u(0, 1);
        std::string read;
        for (long i = start; read.size() < length; i++) {
            char c = i >= 0 && i < static_cast<long>(gene.size()) ? gene[static_cast<size_t>(i)] : "ACGT"[rng() % 4];
            double const r = u(rng);
            if (r < divergence * 0.9) c = "ACGT"[(std::string("ACGT").find(c) + 1 + rng() % 3) % 4];
            else if (r < divergence * 0.95) continue;
            else if (r < divergence) { read += "ACGT"[rng() % 4]; if (read.size() == length) break; }
            if (u(rng) < n_rate) c = 'N';
            read += c;
        }
        return read;
    }
}

TEST(AlignmentScreen, TheBoundForShortReads) {
    // 150 bp inside the gene at -a 0.9: a budget of 61 (MaxScore), so a score of at most 60: 15 mismatches touch
    // 120 8-mers at most, leaving 23 of the 143 intact.
    EXPECT_EQ(AlignmentScreen::MaxTouched(61, kMismatch, kGapOpen, kGapExtend, 8), 120);
    EXPECT_EQ(AlignmentScreen::Required(150, 61, kMismatch, kGapOpen, kGapExtend, 8), 23);
    // With the k the screen takes for a short window, 6: 90 touched, 55 of 145 intact.
    EXPECT_EQ(AlignmentScreen::KFor(150), 6u);
    EXPECT_EQ(AlignmentScreen::Required(150, 61, kMismatch, kGapOpen, kGapExtend), 55);
    // A long read's window of 1500 bases: 600 of score, 1200 touched, 293 of 1493 intact.
    EXPECT_EQ(AlignmentScreen::Required(1500, 601, kMismatch, kGapOpen, kGapExtend, 8), 293);
    // Too short to screen, and nothing to refuse without a budget.
    EXPECT_EQ(AlignmentScreen::Required(AlignmentScreen::kDefaultK - 1, 61, kMismatch, kGapOpen, kGapExtend), 0);
    EXPECT_LT(AlignmentScreen::Required(150, 1000, kMismatch, kGapOpen, kGapExtend), 0);
    // A budget of 0: no alignment at all.
    EXPECT_EQ(AlignmentScreen::Required(150, 0, kMismatch, kGapOpen, kGapExtend), 150);
}

TEST(AlignmentScreen, AReadOfItsOwnGenePasses) {
    std::mt19937 rng(3);
    AlignmentScreen screen;
    std::string const gene = RandomBases(rng, 400);
    for (int n = 0; n < 200; n++) {
        std::string const read = Derive(rng, gene, 100, 150, 0.05, 0);
        EXPECT_TRUE(screen.MayAlign(read, 0, 0, std::string_view(gene).substr(80, 200), 61, kMismatch, kGapOpen, kGapExtend));
    }
    // An unrelated read is refused.
    size_t refused = 0;
    for (int n = 0; n < 200; n++) {
        refused += !screen.MayAlign(RandomBases(rng, 150), 0, 0, std::string_view(gene).substr(80, 200), 61, kMismatch, kGapOpen, kGapExtend);
    }
    EXPECT_EQ(refused, 200u);
}

// What the screen refuses, WFA2 fails: random reads of 2-25% divergence from random genes, running past the gene's
// ends or not, with the free ends and the budget AlignAnchor would give them, aligned ends-free by WFA2 with the
// same budget. The screen's power (refused among the failures) is printed.
TEST(AlignmentScreen, RefusesOnlyWhatWFA2Fails) {
    std::mt19937 rng(11);
    std::uniform_real_distribution<double> u(0, 1);
    AlignmentScreen screen;
    WFA2Wrapper2 wfa(kMismatch, kGapOpen, kGapExtend, 0);
    std::vector<size_t> const ks = { 5, 6, 7, 8, 9, 10 };
    std::vector<size_t> refused_k(ks.size(), 0), refused_short_k(ks.size(), 0), refused_long_k(ks.size(), 0);
    size_t cases = 0, failed = 0, failed_short = 0, failed_long = 0, refused = 0, refused_but_aligned = 0;
    for (int n = 0; n < 20000; n++) {
        size_t const gene_length = 300 + rng() % 1200;
        std::string const gene = RandomBases(rng, gene_length);
        size_t const read_length = n % 4 == 0 ? 500 + rng() % 1000 : 100 + rng() % 151;  // long-read windows now and then
        long const start = static_cast<long>(rng() % (gene_length + 60)) - 40;
        double const divergence = u(rng) < 0.3 ? u(rng) * 0.05 : 0.05 + u(rng) * 0.2;
        std::string const read = Derive(rng, gene, start, read_length, divergence, n % 7 == 0 ? 0.01 : 0);
        // As SimpleAlignmentHandler::AlignAnchor sets the window up: the read's bases past the gene's ends are free,
        // plus a dovetail of 9; the window is the gene's overlapping part with up to 9 bases either side.
        long const abs_pos = start;
        int const begin_free = abs_pos < 0 ? static_cast<int>(-abs_pos) + 9 : 0;
        long const past_end = abs_pos + static_cast<long>(read_length) - static_cast<long>(gene_length);
        int const end_free = past_end > 0 ? static_cast<int>(past_end) + 9 : 0;
        if (begin_free + end_free >= static_cast<int>(read_length)) continue;
        long const ref_start = std::max(0L, abs_pos - 9), ref_end = std::min(static_cast<long>(gene_length), abs_pos + static_cast<long>(read_length) + 9);
        if (ref_end <= ref_start) continue;
        std::string const window = gene.substr(static_cast<size_t>(ref_start), static_cast<size_t>(ref_end - ref_start));
        int const overlap = static_cast<int>(read_length) - begin_free - end_free;
        int const max_score = SimpleAlignmentHandler::MaxScore(0.9, static_cast<uint32_t>(overlap));
        int const ref_begin_free = abs_pos > 0 ? 18 : 0, ref_end_free = past_end < 0 ? 18 : 0;

        bool const may = screen.MayAlign(read, static_cast<size_t>(begin_free), static_cast<size_t>(end_free), window, max_score, kMismatch, kGapOpen, kGapExtend);
        wfa.Reset();
        wfa.Alignment(read, window, begin_free, end_free, ref_begin_free, ref_end_free, max_score);
        bool const aligned = wfa.Success();
        bool const is_long = read_length > 400;
        cases++;
        failed += !aligned;
        (is_long ? failed_long : failed_short) += !aligned;
        refused += !may;
        refused_but_aligned += !may && aligned;
        if (!may) EXPECT_FALSE(aligned) << "case " << n << ": divergence " << divergence << ", read " << read_length << " at " << start;
        // Every k: exact as well, and its power on short reads and long-read windows.
        for (size_t j = 0; j < ks.size(); j++) {
            bool const may_k = screen.MayAlign(read, static_cast<size_t>(begin_free), static_cast<size_t>(end_free), window, max_score, kMismatch, kGapOpen, kGapExtend, ks[j]);
            if (!may_k) EXPECT_FALSE(aligned) << "k " << ks[j] << ", case " << n;
            refused_k[j] += !may_k;
            (is_long ? refused_long_k[j] : refused_short_k[j]) += !may_k;
        }
    }
    std::cout << cases << " candidates: WFA2 failed " << failed << " (short reads " << failed_short << ", long-read windows " << failed_long
              << "), the screen refused " << refused << " (" << (failed ? 100 * refused / failed : 0) << "% of the failures), refused but aligned "
              << refused_but_aligned << std::endl;
    for (size_t j = 0; j < ks.size(); j++) {
        std::cout << "  k = " << ks[j] << ": refused " << refused_k[j] << " (" << 100 * refused_k[j] / std::max<size_t>(failed, 1) << "% of the failures; short reads "
                  << 100 * refused_short_k[j] / std::max<size_t>(failed_short, 1) << "%, long-read windows " << 100 * refused_long_k[j] / std::max<size_t>(failed_long, 1) << "%)" << std::endl;
    }
    ASSERT_GT(cases, 15000u);
    EXPECT_EQ(refused_but_aligned, 0u);
    EXPECT_GT(refused, failed / 10);
}

namespace {
    // The screen's answer from its definition, without its early exits: the read k-mers fully inside the aligned part
    // that occur in the window (whose k-mers count when all their bases are A, C, G or T; a read k-mer with any other
    // base counts as shared), at least Required of them.
    bool DefinitionAnswer(std::string_view read, size_t begin_free, size_t end_free, std::string_view window, int max_score, size_t k) {
        if (begin_free + end_free >= read.size()) return true;
        size_t const aligned = read.size() - begin_free - end_free;
        int64_t const required = AlignmentScreen::Required(aligned, max_score, kMismatch, kGapOpen, kGapExtend, k);
        if (required <= 0) return true;
        auto code = [k](std::string_view s, size_t at, uint32_t& c) {
            c = 0;
            for (size_t j = 0; j < k; j++) {
                uint64_t const b = KmerUtils::BaseToInt(s[at + j]);
                if (b > 3) return false;
                c |= static_cast<uint32_t>(b) << (2 * j);
            }
            return true;
        };
        // Whether each k-mer occurs in the window, for every k up to kMaxK; the k-mers set are cleared again after use.
        static std::vector<char> in_window(size_t{ 1 } << (2 * AlignmentScreen::kMaxK), 0);
        static std::vector<uint32_t> set;
        set.clear();
        uint32_t c = 0;
        for (size_t p = 0; p + k <= window.size(); p++) {
            if (code(window, p, c) && !in_window[c]) {
                in_window[c] = 1;
                set.push_back(c);
            }
        }
        int64_t shared = 0;
        for (size_t p = begin_free; p + k <= read.size() - end_free; p++) shared += !code(read, p, c) || in_window[c];
        for (uint32_t const s : set) in_window[s] = 0;
        return shared >= required;
    }

    // Restores the packed sequences' AVX2 switch (packed::UseAvx2) as it was.
    struct RestoreAvx2 {
        bool const was = packed::Avx2Enabled().load();
        ~RestoreAvx2() { packed::UseAvx2(was); }
    };
}

// The read side from the read's strands packed once (ReadKmers, ReadStretch): for stretches of a read's forward strand and
// of its reverse complement at any offset (four k-mers to a load or not), and long-read windows [s, e) as the stretch
// from n - e of the reverse strand, with Ns, ambiguity codes and lower case in the read (or none: the fast path), every
// k up to kMaxK, free ends and budgets, the answer is the definition's, and that of the screen packing the candidate's
// read for it alone. The window side from the gene's packed bytes (MayAlignPacked) and from its decoded window (MayAlign)
// alike: genes with Ns, ambiguity codes and lower case (packed as bases, as the decoded window has them), some shorter
// than k, windows at any offset (inside a packed byte or not), up to the gene's last base, shorter than k or empty, read
// in exactly the gene's bytes. The reverse strand is given, or made from the read. The reads come one after another into
// the same buffers, ReadKmers and screen, so a read must not see the strands or stamps of the one before.
TEST(AlignmentScreen, TheReadsSharedCodesGiveTheDefinitionsAnswer) {
    RestoreAvx2 restore;
    std::mt19937 rng(37);
    std::uniform_real_distribution<double> u(0, 1);
    AlignmentScreen screen;
    AlignmentScreen::ReadKmers kmers;
    std::string read, rev;
    std::vector<size_t> const ks = { 0, 4, 5, 6, 7, 8, 9, 10 };
    size_t cases = 0, refused = 0, passed = 0;
    for (int n = 0; n < 1500; n++) {
        size_t const gene_length = n % 50 == 0 ? 1 + rng() % 12 : 50 + rng() % 3000;
        std::string gene = RandomBases(rng, gene_length);
        for (auto& c : gene) if (rng() % 50 == 0) c = "NRYKMSWnacgt"[rng() % 12];  // packed as a base
        std::vector<uint8_t> packed(packed::Bytes(gene_length), 0);  // exactly the gene's bytes
        packed::Pack(gene.data(), gene_length, packed.data());
        std::string decoded(gene_length, '\0');
        packed::Unpack(packed.data(), gene_length, decoded.data());
        size_t const read_length = 20 + rng() % (n % 4 == 0 ? 4000 : 300);
        long const start = static_cast<long>(rng() % gene_length) - 30;
        double const divergence = u(rng) < 0.5 ? u(rng) * 0.05 : u(rng) * 0.3;
        read.assign(Derive(rng, decoded, start, read_length, divergence, n % 3 == 0 ? 0.02 : 0));
        if (n % 5 != 0) for (auto& c : read) if (rng() % 100 == 0) c = "RYacgtn-"[rng() % 8];
        KmerUtils::ReverseComplementInto(read, rev);
        if (n % 2 == 0) kmers.Set(read);
        else kmers.Set(read, rev);
        packed::UseAvx2(n % 4 != 3);  // the strands packed and marked without AVX2 now and then
        for (int candidate = 0; candidate < 12; candidate++) {
            bool const forward = rng() % 2;
            size_t const s = rng() % read_length, e = s + 1 + rng() % (read_length - s);  // the window [s, e) of the read
            std::string const stretch = forward ? read.substr(s, e - s) : KmerUtils::ReverseComplement(std::string_view(read).substr(s, e - s));
            size_t const offset = forward ? s : read_length - e;
            ASSERT_EQ(stretch, (forward ? read : rev).substr(offset, e - s));
            size_t const length = stretch.size();
            size_t const begin = rng() % (gene_length + 1);
            size_t const end = candidate % 4 == 0 ? gene_length : std::min(gene_length, begin + rng() % 1600);
            size_t const begin_free = rng() % 3 == 0 ? rng() % (length / 2 + 1) : 0;
            size_t const end_free = rng() % 3 == 0 ? rng() % (length / 2 + 1) : 0;
            int const max_score = 1 + static_cast<int>(rng() % (length + 100));
            std::string_view const window = std::string_view(decoded).substr(begin, end - begin);
            for (size_t k : ks) {
                size_t const k_used = std::min(k == 0 ? AlignmentScreen::KFor(end - begin) : k, AlignmentScreen::kMaxK);
                bool const expected = DefinitionAnswer(stretch, begin_free, end_free, window, max_score, k_used);
                bool const shared = screen.MayAlignPacked(kmers, forward, offset, length, begin_free, end_free, packed.data(), packed.size(),
                                                          begin, end, max_score, kMismatch, kGapOpen, kGapExtend, k);
                bool const own = screen.MayAlignPacked(stretch, begin_free, end_free, packed.data(), packed.size(), begin, end,
                                                       max_score, kMismatch, kGapOpen, kGapExtend, k);
                bool const text = screen.MayAlign(stretch, begin_free, end_free, window, max_score, kMismatch, kGapOpen, kGapExtend, k);
                ASSERT_EQ(shared, expected) << "read " << n << " (" << read_length << " bases), " << (forward ? "forward" : "reverse")
                                            << " window " << s << "-" << e << ", k " << k << ", gene window " << begin << "-" << end
                                            << ", free " << begin_free << "/" << end_free << ", budget " << max_score;
                ASSERT_EQ(own, expected) << "read " << n << ", candidate " << candidate << ", k " << k;
                ASSERT_EQ(text, expected) << "read " << n << ", candidate " << candidate << ", k " << k;
                cases++;
                refused += !expected;
                passed += expected;
            }
        }
    }
    std::cout << cases << " screens: " << refused << " refused, " << passed << " passed, the same from the read's packed strands" << std::endl;
    EXPECT_GT(refused, cases / 10);
    EXPECT_GT(passed, cases / 10);
}

namespace {
    constexpr size_t kGenes = 12, kGeneLength = 1200;

    // Random genes of taxid 1, loaded (as test_AnchoredAlignment.cpp).
    struct RandomReference : test::LoadedReference {
        std::vector<std::string> const& genes = LoadedReference::genes.at(1);

        RandomReference() : LoadedReference({ { 1, Genes() } }, "screen") {}

        static std::vector<std::string> Genes() {
            std::mt19937 rng(29);
            return test::RandomGenes(std::vector<size_t>(kGenes, kGeneLength), rng);
        }
    };

    // The maximal exact runs of at least min_length on the diagonal.
    ChainList ExactRuns(std::string const& read, std::string const& gene, long diagonal, size_t min_length = 15) {
        ChainList runs;
        size_t i = 0;
        while (i < read.size()) {
            long const g = diagonal + static_cast<long>(i);
            if (g < 0 || g >= static_cast<long>(gene.size()) || read[i] != gene[static_cast<size_t>(g)]) { i++; continue; }
            size_t j = i;
            while (j < read.size() && diagonal + static_cast<long>(j) < static_cast<long>(gene.size()) && read[j] == gene[static_cast<size_t>(diagonal + static_cast<long>(j))]) j++;
            if (j - i >= min_length) runs.emplace_back(static_cast<uint32_t>(diagonal + static_cast<long>(i)), static_cast<uint16_t>(i), static_cast<uint16_t>(j - i));
            i = j;
        }
        return runs;
    }
}

// Through the handler, on random anchors of reads of every divergence and of either strand: the screen on and off give
// the same outcome, score and CIGAR for every candidate, from the anchor's exact matches and as whole windows; for one
// anchor (AlignAnchor, the screen packing the candidate's read for it alone) and as protal runs it (operator(), the
// screen taking the k-mers of the read's strands packed once for all its anchors).
TEST(AlignmentScreen, TheHandlerAlignsTheSameWithAndWithoutIt) {
    RandomReference ref;
    WFA2Wrapper2 aligner{kMismatch, kGapOpen, kGapExtend, 1000};
    SimpleAlignmentHandler with(*ref.loader, aligner, 31, 3, 0.9, false), without(*ref.loader, aligner, 31, 3, 0.9, false);
    without.SetAlignmentScreen(false);
    std::mt19937 rng(17);
    std::uniform_real_distribution<double> u(0, 1);
    size_t cases = 0, aligned = 0, reverse_aligned = 0;
    auto expect_same = [](AlignmentResult const& a, AlignmentResult const& b, int n) {
        EXPECT_EQ(a.AlignmentScore(), b.AlignmentScore()) << "case " << n;
        EXPECT_EQ(a.GetAlignmentInfo().cigar, b.GetAlignmentInfo().cigar) << "case " << n;
        EXPECT_EQ(a.GetAlignmentInfo().gene_alignment_start, b.GetAlignmentInfo().gene_alignment_start) << "case " << n;
        EXPECT_EQ(a.Forward(), b.Forward()) << "case " << n;
    };
    for (bool anchored : { true, false }) {
        with.SetAnchoredAlignment(anchored);
        without.SetAnchoredAlignment(anchored);
        for (int n = 0; n < 3000; n++) {
            uint32_t const gene_id = 1 + rng() % kGenes;
            std::string const& gene = ref.genes[gene_id - 1];
            long const start = static_cast<long>(rng() % (kGeneLength + 60)) - 40;
            double const divergence = u(rng) < 0.4 ? u(rng) * 0.05 : u(rng) * 0.25;
            std::string const segment = Derive(rng, gene, start, 150, divergence, 0.002);  // the read on the gene's strand
            ChainList runs = ExactRuns(segment, gene, start);
            if (runs.empty()) continue;
            ChainList chain = rng() % 2 ? ChainList{ runs[rng() % runs.size()] } : runs;
            // A read of the other strand is the segment's reverse complement; its anchor's chain is on the segment.
            bool const forward = rng() % 2;
            std::string const read = forward ? segment : KmerUtils::ReverseComplement(segment);
            std::string const rev = KmerUtils::ReverseComplement(read);
            std::string id = "r";
            ChainAlignmentAnchor anchor(1, gene_id, forward);
            anchor.chain = chain;
            ChainAlignmentAnchor a = anchor, b = anchor;
            AlignmentResult ra, rb;
            size_t const screened_before = with.m_screened_alignments;
            bool const oa = with.AlignAnchor(a, ra, read, rev, false, id);
            bool const ob = without.AlignAnchor(b, rb, read, rev, false, id);
            cases++;
            ASSERT_EQ(oa, ob) << "case " << n << " divergence " << divergence;
            if (oa) {
                aligned++;
                reverse_aligned += !forward;
                expect_same(ra, rb, n);
            }
            size_t const screened_alone = with.m_screened_alignments - screened_before;
            AlignmentAnchorList anchors_with{ anchor }, anchors_without{ anchor };
            AlignmentResultList la, lb;
            with(anchors_with, la, read, rev, 3, id);
            without(anchors_without, lb, read, rev, 3, id);
            ASSERT_EQ(with.m_screened_alignments - screened_before - screened_alone, screened_alone) << "case " << n;
            ASSERT_EQ(la.size(), lb.size()) << "case " << n;
            // operator() keeps the alignments of at least the -a identity, here AlignAnchor's.
            ASSERT_EQ(la.size(), oa && ra.GetAlignmentInfo().GetProxyANI() >= 0.9 ? 1u : 0u) << "case " << n;
            if (!la.empty()) {
                expect_same(la[0], lb[0], n);
                expect_same(la[0], ra, n);
            }
        }
    }
    std::cout << cases << " candidates, " << aligned << " aligned (" << reverse_aligned << " of reverse reads); the screen refused "
              << with.m_screened_alignments << " of " << with.m_attempted_alignments << " (WFA2 ran " << with.m_anchored_alignments +
              with.m_whole_window_alignments << " times with it, " << without.m_anchored_alignments + without.m_whole_window_alignments
              << " without)" << std::endl;
    ASSERT_GT(cases, 4000u);
    EXPECT_GT(reverse_aligned, aligned / 4);
    EXPECT_GT(with.m_screened_alignments, 0u);
    EXPECT_EQ(without.m_screened_alignments, 0u);
}

namespace {
    // The LCS by dynamic programming, with the indel bound's matches: equal bases, and a base other than A, C, G or T
    // (either case; KmerUtils::BaseToInt) against any.
    size_t LcsByDp(std::string_view a, std::string_view b) {
        std::vector<size_t> row(b.size() + 1, 0), next(b.size() + 1, 0);
        std::vector<uint64_t> codes(b.size());
        for (size_t j = 0; j < b.size(); j++) codes[j] = KmerUtils::BaseToInt(b[j]);
        for (char const x : a) {
            uint64_t const cx = KmerUtils::BaseToInt(x);
            for (size_t j = 0; j < b.size(); j++) {
                uint64_t const cy = codes[j];
                bool const match = cx > 3 || cy > 3 || cx == cy;
                next[j + 1] = match ? row[j] + 1 : std::max(row[j + 1], next[j]);
            }
            std::swap(row, next);
        }
        return row[b.size()];
    }

    // The indel bound's answer from its definition: n_req + m_req bases required, the LCS by dynamic programming.
    bool IndelDefinitionAnswer(std::string_view read, size_t begin_free, size_t end_free, std::string_view window, size_t ref_begin_free,
                               size_t ref_end_free, int max_score) {
        size_t const n_req = begin_free + end_free < read.size() ? read.size() - begin_free - end_free : 0;
        size_t const m_req = ref_begin_free + ref_end_free < window.size() ? window.size() - ref_begin_free - ref_end_free : 0;
        // q (n_req + m_req - 2 LCS) > 2 (max_score - 1) refuses, q = min(4, 2 * 2) = 4
        int64_t const lower = 4 * (static_cast<int64_t>(n_req + m_req) - 2 * static_cast<int64_t>(LcsByDp(read, window)));
        return lower <= 2 * (static_cast<int64_t>(max_score) - 1);
    }

    // An ONT-like copy of `source`: substitutions, and insertions and deletions three times as likely in a homopolymer, at
    // the rate `divergence` in all (about half substitutions); Ns at n_rate. With `gene_of`, the source position of each
    // base of the copy (-1 for an inserted one).
    std::string OntCopy(std::mt19937& rng, std::string const& source, double divergence, double n_rate = 0,
                        std::vector<long>* gene_of = nullptr) {
        std::uniform_real_distribution<double> u(0, 1);
        std::string read;
        for (size_t i = 0; i < source.size(); i++) {
            char c = source[i];
            bool const homopolymer = (i > 0 && source[i - 1] == c) || (i + 1 < source.size() && source[i + 1] == c);
            double const indel = divergence * 0.25 * (homopolymer ? 3 : 1);
            double const r = u(rng);
            if (r < indel) continue;  // deletion
            if (r < 2 * indel) {      // insertion
                read += homopolymer ? c : "ACGT"[rng() % 4];
                if (gene_of) gene_of->push_back(-1);
            }
            if (u(rng) < divergence * 0.5) c = "ACGT"[(std::string("ACGT").find(c) + 1 + rng() % 3) % 4];
            if (u(rng) < n_rate) c = 'N';
            read += c;
            if (gene_of) gene_of->push_back(static_cast<long>(i));
        }
        return read;
    }

    // `source` with substitutions at a rate between `low` and `high`: a relative of a gene.
    std::string Relative(std::mt19937& rng, std::string source, double low, double high) {
        std::uniform_real_distribution<double> u(0, 1);
        double const rate = low + u(rng) * (high - low);
        for (auto& c : source) if (u(rng) < rate) c = "ACGT"[(std::string("ACGT").find(c) + 1 + rng() % 3) % 4];
        return source;
    }
}

// The bit-parallel LCS is the dynamic programming one, for lengths around the 64-bit words and with Ns, ambiguity codes
// and lower case on either side; and the bound's answer is its definition's, from the gene's packed bytes (whose Ns and
// codes are packed as bases, as the decoded window has them) and from the decoded window, with free ends at either end of
// either sequence and any budget, the exits included. The candidates come one after another into the same buffers.
TEST(IndelBound, TheLcsAndTheAnswerAreTheDefinitions) {
    std::mt19937 rng(53);
    std::uniform_real_distribution<double> u(0, 1);
    IndelBound bound;
    auto noisy = [&rng](std::string s) {
        for (auto& c : s) if (rng() % 40 == 0) c = "NRYKnacgt"[rng() % 9];
        return s;
    };
    for (size_t const n : { 0, 1, 2, 63, 64, 65, 127, 128, 129, 300 }) {
        for (size_t const m : { 0, 1, 63, 64, 65, 128, 129, 191, 257 }) {
            for (int repeat = 0; repeat < 3; repeat++) {
                std::string const window = noisy(RandomBases(rng, m));
                std::string const read = repeat == 0 ? noisy(RandomBases(rng, n)) :
                                         noisy(OntCopy(rng, window + RandomBases(rng, n), 0.1)).substr(0, n);
                ASSERT_EQ(bound.Lcs(read, window), LcsByDp(read, window)) << "read " << n << ", window " << m << ", repeat " << repeat;
            }
        }
    }
    size_t cases = 0, refused = 0, passed = 0;
    for (int n = 0; n < 1500; n++) {
        size_t const gene_length = n % 50 == 0 ? 1 + rng() % 12 : 50 + rng() % 700;
        std::string const gene = noisy(RandomBases(rng, gene_length));
        std::vector<uint8_t> packed(packed::Bytes(gene_length), 0);  // exactly the gene's bytes
        packed::Pack(gene.data(), gene_length, packed.data());
        std::string decoded(gene_length, '\0');
        packed::Unpack(packed.data(), gene_length, decoded.data());
        size_t const begin = rng() % (gene_length + 1);
        size_t const end = n % 4 == 0 ? gene_length : std::min(gene_length, begin + rng() % 900);
        std::string const window(std::string_view(decoded).substr(begin, end - begin));
        // From the window with errors, with flanks that run past it, or unrelated.
        double const divergence = u(rng) < 0.5 ? u(rng) * 0.1 : u(rng) * 0.4;
        std::string const read = noisy(n % 3 == 0 ? RandomBases(rng, 20 + rng() % 800) :
                                       RandomBases(rng, rng() % 100) + OntCopy(rng, window, divergence) + RandomBases(rng, rng() % 100));
        size_t const begin_free = rng() % 3 == 0 ? rng() % (read.size() / 2 + 1) : rng() % 3 == 0 ? read.size() : 0;
        size_t const end_free = rng() % 3 == 0 ? rng() % (read.size() / 2 + 1) : 0;
        size_t const ref_begin_free = rng() % 3 == 0 ? rng() % 20 : 0;
        size_t const ref_end_free = rng() % 5 == 0 ? rng() % (window.size() + 2) : 0;
        int const max_score = rng() % 10 == 0 ? static_cast<int>(rng() % 3) :
                              1 + static_cast<int>(rng() % static_cast<uint32_t>(2 * (read.size() + window.size()) + 10));
        bool const expected = IndelDefinitionAnswer(read, begin_free, end_free, window, ref_begin_free, ref_end_free, max_score);
        bool const from_packed = bound.MayAlignPacked(read, begin_free, end_free, packed.data(), begin, end, ref_begin_free, ref_end_free,
                                                      max_score, kMismatch, kGapOpen, kGapExtend);
        bool const from_text = bound.MayAlign(read, begin_free, end_free, window, ref_begin_free, ref_end_free, max_score, kMismatch,
                                              kGapOpen, kGapExtend);
        ASSERT_EQ(from_packed, expected) << "case " << n << ": read " << read.size() << " (free " << begin_free << "/" << end_free
                                         << "), window " << begin << "-" << end << " (free " << ref_begin_free << "/" << ref_end_free
                                         << "), budget " << max_score;
        ASSERT_EQ(from_text, expected) << "case " << n;
        cases++;
        refused += !expected;
        passed += expected;
    }
    std::cout << cases << " candidates: " << refused << " refused, " << passed << " passed, as the definition" << std::endl;
    EXPECT_GT(refused, cases / 10);
    EXPECT_GT(passed, cases / 10);
    // Without a budget nothing aligns; a penalty that gives no bound (a mismatch of 0) refuses nothing.
    EXPECT_FALSE(IndelBound().MayAlign("ACGT", 0, 0, "ACGT", 0, 0, 0, kMismatch, kGapOpen, kGapExtend));
    EXPECT_TRUE(IndelBound().MayAlign("AAAA", 0, 0, "CCCC", 0, 0, 1, 0, kGapOpen, kGapExtend));
    EXPECT_EQ(IndelBound::MaxRefusedLcs(1000, 601, kMismatch, kGapOpen, kGapExtend), 349);  // 4 (1000 - 2 L) >= 1201
    EXPECT_EQ(IndelBound::MaxRefusedLcs(100, 601, kMismatch, kGapOpen, kGapExtend), -1);
}

// What the indel bound refuses, WFA2 fails: ONT-like long-read windows (LongReadAligner: the gene with a margin of
// 100 + L / 20 bases either side, the read's bases past the gene free) of the gene itself at 2-15% divergence, of a
// relative 10-25% away, or of unrelated sequence, aligned ends-free by WFA2 without X-drop at ONT's floor of 0.85, as
// RunProtal aligns long reads. Its power (refused among the failures, by kind) is printed: the unrelated windows, most
// of the failing ONT candidates at r226, are refused, which the k-mer screen cannot at 0.85.
TEST(IndelBound, RefusesOnlyWhatWFA2Fails) {
    std::mt19937 rng(61);
    std::uniform_real_distribution<double> u(0, 1);
    IndelBound bound;
    AlignmentScreen screen;
    WFA2Wrapper2 wfa(kMismatch, kGapOpen, kGapExtend, 0);
    char const* const kinds[] = { "own gene", "relative", "unrelated" };
    std::array<size_t, 3> cases{}, failed{}, refused{}, screened{};
    size_t refused_but_aligned = 0;
    for (int n = 0; n < 900; n++) {
        int const kind = n % 3;
        size_t const gene_length = n % 5 == 0 ? 250 + rng() % 250 : 500 + rng() % 1500;
        std::string const gene = RandomBases(rng, gene_length);
        size_t const margin = 100 + gene_length / 20;
        std::string const source = kind == 0 ? gene : kind == 1 ? Relative(rng, gene, 0.1, 0.25) : RandomBases(rng, gene_length);
        double const divergence = kind == 0 ? 0.02 + u(rng) * 0.13 : 0.05 + u(rng) * 0.05;
        // The read's window from `start` (gene coordinates): before the gene by up to the margin, or inside its first third.
        long const start = static_cast<long>(rng() % (margin + gene_length / 3)) - static_cast<long>(margin);
        std::string const before = start < 0 ? RandomBases(rng, static_cast<size_t>(-start)) : "";
        std::string const read = OntCopy(rng, before + source.substr(static_cast<size_t>(std::max(0L, start))) + RandomBases(rng, margin),
                                         divergence, n % 7 == 0 ? 0.01 : 0);
        // As AlignAnchor sets the window up from the chain's diagonals: the read's bases before and after the gene free,
        // plus 9; the gene from 9 before the read's start, with a dovetail of 18 free gene bases where the read starts in it.
        size_t const begin_free = start < 0 ? static_cast<size_t>(-start) + 9 : 0;
        size_t const end_free = margin + 9;
        if (begin_free + end_free >= read.size()) continue;
        size_t const ref_start = static_cast<size_t>(std::max(0L, start - 9));
        size_t const ref_begin_free = start > 0 ? 18 : 0;
        std::string const window = gene.substr(ref_start);
        int const overlap = static_cast<int>(gene_length - static_cast<size_t>(std::max(0L, start)));
        int const max_score = SimpleAlignmentHandler::MaxScore(0.85, static_cast<uint32_t>(overlap));

        bool const may = bound.MayAlign(read, begin_free, end_free, window, ref_begin_free, 0, max_score, kMismatch, kGapOpen, kGapExtend);
        bool const may_screen = screen.MayAlign(read, begin_free, end_free, window, max_score, kMismatch, kGapOpen, kGapExtend);
        wfa.Reset();
        wfa.Alignment(read, window, static_cast<int>(begin_free), static_cast<int>(end_free), static_cast<int>(ref_begin_free), 0, max_score);
        bool const aligned = wfa.Success();
        if (!may) EXPECT_FALSE(aligned) << kinds[kind] << " case " << n << ": divergence " << divergence << ", gene " << gene_length;
        refused_but_aligned += !may && aligned;
        cases[kind]++;
        failed[kind] += !aligned;
        refused[kind] += !may;
        screened[kind] += !may_screen;
    }
    for (int kind = 0; kind < 3; kind++) {
        std::cout << kinds[kind] << ": " << cases[kind] << " candidates, WFA2 failed " << failed[kind] << ", the indel bound refused "
                  << refused[kind] << ", the k-mer screen " << screened[kind] << std::endl;
    }
    EXPECT_EQ(refused_but_aligned, 0u);
    // Nearly every unrelated window is refused, and few of the own gene's.
    EXPECT_GT(cases[2], 250u);
    EXPECT_GE(refused[2], failed[2] * 95 / 100);
    EXPECT_LT(refused[0], cases[0] / 4);
}

namespace {
    // A long read's window over a gene as LongReadAligner aligns it: `margin` random bases either side of `source` (the
    // gene, a relative or unrelated sequence) copied with ONT-like errors; and its chain: the exact runs of at least 20
    // bases on the gene along the true path within 6 diagonals of the first (ChainAnchorFinder), or, for unrelated
    // sequence, a planted exact 24-mer of the gene, as a spurious seed pair.
    std::pair<std::string, ChainList> LongReadWindow(std::mt19937& rng, std::string const& gene, std::string const& source, size_t margin,
                                                     double divergence, bool related) {
        std::vector<long> gene_of(margin, -1);
        std::string read = RandomBases(rng, margin) + OntCopy(rng, source, divergence, 0, &gene_of) + RandomBases(rng, margin);
        gene_of.resize(read.size(), -1);
        ChainList runs;
        if (related) {
            size_t i = 0;
            while (i < read.size()) {
                if (gene_of[i] < 0 || read[i] != gene[static_cast<size_t>(gene_of[i])]) { i++; continue; }
                size_t j = i + 1;
                while (j < read.size() && gene_of[j] == gene_of[j - 1] + 1 && read[j] == gene[static_cast<size_t>(gene_of[j])]) j++;
                if (j - i >= 20) runs.emplace_back(static_cast<uint32_t>(gene_of[i]), static_cast<uint16_t>(i), static_cast<uint16_t>(j - i));
                i = j;
            }
        } else {
            size_t const g = rng() % (gene.size() - 24), at = margin + g * (read.size() - 2 * margin) / gene.size();
            read.replace(at, 24, gene, g, 24);
            runs.emplace_back(static_cast<uint32_t>(g), static_cast<uint16_t>(at), static_cast<uint16_t>(24));
        }
        ChainList chain;
        if (runs.empty()) return { read, chain };
        long const d0 = static_cast<long>(runs.front().genepos) - static_cast<long>(runs.front().readpos);
        for (auto const& run : runs) {
            if (std::abs(static_cast<long>(run.genepos) - static_cast<long>(run.readpos) - d0) <= 6) chain.push_back(run);
        }
        return { read, chain };
    }
}

// Through the handler as RunProtal sets it up for ONT reads (-a 0.85, no X-drop, anchored through the chain with
// indels, the indel bound in place of the k-mer screen), and with the k-mer screen instead (--no_indel_bound): the same
// outcome, score, start and CIGAR for every candidate, of the gene itself, a relative or unrelated sequence, of either
// strand, from the chain or as a whole window.
TEST(IndelBound, TheHandlerAlignsLongReadsTheSameWithAndWithoutIt) {
    RandomReference ref;
    WFA2Wrapper2 aligner{kMismatch, kGapOpen, kGapExtend, 0};
    SimpleAlignmentHandler with(*ref.loader, aligner, 31, 3, 0.85, false), without(*ref.loader, aligner, 31, 3, 0.85, false);
    for (auto* handler : { &with, &without }) handler->SetAnchoredIndels(true);
    with.SetAlignmentScreen(false);
    with.SetIndelBound(true);
    std::mt19937 rng(71);
    std::uniform_real_distribution<double> u(0, 1);
    size_t cases = 0, aligned = 0;
    for (bool anchored : { true, false }) {
        with.SetAnchoredAlignment(anchored);
        without.SetAnchoredAlignment(anchored);
        for (int n = 0; n < 240; n++) {
            int const kind = n % 3;  // the gene, a relative, unrelated
            uint32_t const gene_id = 1 + rng() % kGenes;
            std::string const& gene = ref.genes[gene_id - 1];
            std::string const source = kind == 0 ? gene : kind == 1 ? Relative(rng, gene, 0.1, 0.2) : RandomBases(rng, kGeneLength);
            double const divergence = kind == 0 ? 0.02 + u(rng) * 0.06 : 0.05 + u(rng) * 0.05;  // the own gene mostly aligns
            auto [segment, chain] = LongReadWindow(rng, gene, source, 100 + kGeneLength / 20, divergence, kind != 2);
            if (chain.empty()) continue;
            // A read of the other strand is the segment's reverse complement; its anchor's chain is on the segment.
            bool const forward = rng() % 2;
            std::string const read = forward ? segment : KmerUtils::ReverseComplement(segment);
            std::string const rev = KmerUtils::ReverseComplement(read);
            std::string id = "r";
            ChainAlignmentAnchor a(1, gene_id, forward), b(1, gene_id, forward);
            a.chain = chain;
            b.chain = chain;
            AlignmentResult ra, rb;
            bool const oa = with.AlignAnchor(a, ra, read, rev, false, id);
            bool const ob = without.AlignAnchor(b, rb, read, rev, false, id);
            cases++;
            ASSERT_EQ(oa, ob) << "case " << n << ", kind " << kind << ", divergence " << divergence;
            if (!oa) continue;
            aligned++;
            EXPECT_EQ(ra.AlignmentScore(), rb.AlignmentScore()) << "case " << n;
            EXPECT_EQ(ra.GetAlignmentInfo().cigar, rb.GetAlignmentInfo().cigar) << "case " << n;
            EXPECT_EQ(ra.GetAlignmentInfo().gene_alignment_start, rb.GetAlignmentInfo().gene_alignment_start) << "case " << n;
        }
    }
    std::cout << cases << " long-read candidates, " << aligned << " aligned; the indel bound refused " << with.m_indel_refused
              << ", the k-mer screen in its place " << without.m_screened_alignments << " (WFA2 ran " << with.m_anchored_alignments + with.m_whole_window_alignments
              << " times with it, " << without.m_anchored_alignments + without.m_whole_window_alignments << " without)" << std::endl;
    ASSERT_GT(cases, 400u);
    EXPECT_GT(aligned, cases / 5);
    EXPECT_GT(with.m_indel_refused, cases / 4);
    EXPECT_EQ(without.m_indel_refused, 0u);
    EXPECT_EQ(with.m_screened_alignments, 0u);
    EXPECT_EQ(with.m_attempted_alignments, with.m_screened_alignments + with.m_indel_refused + with.m_anchored_alignments +
                                           with.m_whole_window_alignments);
}
