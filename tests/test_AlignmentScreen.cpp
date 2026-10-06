// Unit tests for the k-mer screen before WFA2 (AlignmentScreen.h): a candidate it refuses is one WFA2
// fails within the same budget, with free ends and Ns; and through SimpleAlignmentHandler the same
// alignments with the screen on and off.
#include <gtest/gtest.h>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <memory>
#include <random>
#include <string>
#include <vector>
#include <unistd.h>
#include "Alignment/AlignmentScreen.h"
#include "Core/AlignmentStrategy.h"

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
    EXPECT_EQ(AlignmentScreen::Required(150, 61, kMismatch, kGapOpen, kGapExtend), 55);
    // A long read's window of 1500 bases: 600 of score, 1200 touched, 293 of 1493 intact.
    EXPECT_EQ(AlignmentScreen::Required(1500, 601, kMismatch, kGapOpen, kGapExtend, 8), 293);
    // Too short to screen, and nothing to refuse without a budget.
    EXPECT_EQ(AlignmentScreen::Required(5, 61, kMismatch, kGapOpen, kGapExtend), 0);
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

// The screen from a gene's packed bytes (MayAlignPacked) answers as the screen of its decoded window (MayAlign): genes
// with Ns and ambiguity codes (packed as bases, as the decoded window has them), windows at any offset (inside a packed
// byte or not), up to the gene's last base, shorter than k or empty, read in exactly the gene's bytes; reads of every
// divergence with Ns, free ends, budgets and k; one screen for both, its stamps reused between them.
TEST(AlignmentScreen, PackedWindowsGiveTheDecodedWindowsAnswer) {
    std::mt19937 rng(23);
    std::uniform_real_distribution<double> u(0, 1);
    AlignmentScreen screen;
    std::vector<size_t> const ks = { 0, 4, 5, 6, 7, 8, 9, 10 };
    size_t cases = 0, refused = 0, passed = 0;
    for (int n = 0; n < 6000; n++) {
        size_t const gene_length = 1 + rng() % 2500;
        std::string gene = RandomBases(rng, gene_length);
        for (auto& c : gene) if (rng() % 50 == 0) c = "NRYKMSWnacgt"[rng() % 12];
        std::vector<uint8_t> packed(packed::Bytes(gene_length), 0);  // exactly the gene's bytes
        packed::Pack(gene.data(), gene_length, packed.data());
        std::string decoded(gene_length, '\0');
        packed::Unpack(packed.data(), gene_length, decoded.data());
        size_t const begin = rng() % (gene_length + 1);
        size_t const end = n % 5 == 0 ? gene_length : std::min(gene_length, begin + rng() % 1600);
        size_t const read_length = 5 + rng() % (n % 3 == 0 ? 1500 : 200);
        long const start = static_cast<long>(begin) + static_cast<long>(rng() % 40) - 20;
        double const divergence = u(rng) < 0.5 ? u(rng) * 0.05 : u(rng) * 0.3;
        std::string const read = Derive(rng, decoded, start, read_length, divergence, n % 4 == 0 ? 0.02 : 0);
        size_t const begin_free = rng() % 3 == 0 ? rng() % (read_length / 2 + 1) : 0;
        size_t const end_free = rng() % 3 == 0 ? rng() % (read_length / 2 + 1) : 0;
        int const max_score = 1 + static_cast<int>(rng() % 400);
        std::string_view const window = std::string_view(decoded).substr(begin, end - begin);
        for (size_t k : ks) {
            bool const from_text = screen.MayAlign(read, begin_free, end_free, window, max_score, kMismatch, kGapOpen, kGapExtend, k);
            bool const from_packed = screen.MayAlignPacked(read, begin_free, end_free, packed.data(), packed.size(), begin, end,
                                                           max_score, kMismatch, kGapOpen, kGapExtend, k);
            ASSERT_EQ(from_packed, from_text) << "case " << n << ", k " << k << ": gene " << gene_length << ", window " << begin << "-"
                                              << end << ", read " << read_length << ", free " << begin_free << "/" << end_free;
            cases++;
            refused += !from_text;
            passed += from_text;
        }
    }
    std::cout << cases << " screens: " << refused << " refused, " << passed << " passed, the same from the packed bytes" << std::endl;
    EXPECT_GT(refused, cases / 10);
    EXPECT_GT(passed, cases / 10);
}

namespace {
    constexpr size_t kGenes = 12, kGeneLength = 1200;

    // Random genes of taxid 1 in a temporary directory, loaded (as test_AnchoredAlignment.cpp).
    struct RandomReference {
        std::vector<std::string> genes;
        std::filesystem::path dir;
        std::unique_ptr<GenomeLoader> loader;

        RandomReference() {
            std::mt19937 rng(29);
            dir = std::filesystem::temp_directory_path() / ("protal_screen_test_" + std::to_string(::getpid()));
            std::filesystem::create_directories(dir);
            std::ofstream fna(dir / "reference.fna"), map(dir / "reference.map");
            size_t offset = 0;
            for (size_t id = 1; id <= kGenes; id++) {
                std::string const seq = RandomBases(rng, kGeneLength);
                std::string header = ">1_" + std::to_string(id) + "\n";
                fna << header << seq << '\n';
                map << "1\t" << id << '\t' << offset + header.size() << '\t' << offset + header.size() + seq.size() << '\n';
                offset += header.size() + seq.size() + 1;
                genes.push_back(seq);
            }
            fna.close();
            map.close();
            loader = std::make_unique<GenomeLoader>((dir / "reference.fna").string(), (dir / "reference.map").string());
            loader->LoadAllGenomes();
        }
        ~RandomReference() { std::filesystem::remove_all(dir); }
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

// Through the handler, on random anchors of reads of every divergence: the screen on and off give the same
// outcome, score and CIGAR for every candidate, from the anchor's exact matches and as whole windows.
TEST(AlignmentScreen, TheHandlerAlignsTheSameWithAndWithoutIt) {
    RandomReference ref;
    WFA2Wrapper2 aligner{kMismatch, kGapOpen, kGapExtend, 1000};
    SimpleAlignmentHandler with(*ref.loader, aligner, 31, 3, 0.9, false), without(*ref.loader, aligner, 31, 3, 0.9, false);
    without.SetAlignmentScreen(false);
    std::mt19937 rng(17);
    std::uniform_real_distribution<double> u(0, 1);
    size_t cases = 0, aligned = 0;
    for (bool anchored : { true, false }) {
        with.SetAnchoredAlignment(anchored);
        without.SetAnchoredAlignment(anchored);
        for (int n = 0; n < 3000; n++) {
            uint32_t const gene_id = 1 + rng() % kGenes;
            std::string const& gene = ref.genes[gene_id - 1];
            long const start = static_cast<long>(rng() % (kGeneLength + 60)) - 40;
            double const divergence = u(rng) < 0.4 ? u(rng) * 0.05 : u(rng) * 0.25;
            std::string const read = Derive(rng, gene, start, 150, divergence, 0.002);
            ChainList runs = ExactRuns(read, gene, start);
            if (runs.empty()) continue;
            ChainList chain = rng() % 2 ? ChainList{ runs[rng() % runs.size()] } : runs;
            std::string const rev = KmerUtils::ReverseComplement(read);
            std::string id = "r";
            ChainAlignmentAnchor a(1, gene_id, true), b(1, gene_id, true);
            a.chain = chain;
            b.chain = chain;
            AlignmentResult ra, rb;
            bool const oa = with.AlignAnchor(a, ra, read, rev, false, id);
            bool const ob = without.AlignAnchor(b, rb, read, rev, false, id);
            cases++;
            ASSERT_EQ(oa, ob) << "case " << n << " divergence " << divergence;
            if (oa) {
                aligned++;
                EXPECT_EQ(ra.AlignmentScore(), rb.AlignmentScore());
                EXPECT_EQ(ra.GetAlignmentInfo().cigar, rb.GetAlignmentInfo().cigar);
                EXPECT_EQ(ra.GetAlignmentInfo().gene_alignment_start, rb.GetAlignmentInfo().gene_alignment_start);
            }
        }
    }
    std::cout << cases << " candidates, " << aligned << " aligned; the screen refused " << with.m_screened_alignments
              << " of " << with.m_attempted_alignments << " (WFA2 ran " << with.m_anchored_alignments + with.m_whole_window_alignments
              << " times with it, " << without.m_anchored_alignments + without.m_whole_window_alignments << " without)" << std::endl;
    ASSERT_GT(cases, 4000u);
    EXPECT_GT(with.m_screened_alignments, 0u);
    EXPECT_EQ(without.m_screened_alignments, 0u);
    EXPECT_EQ(with.m_attempted_alignments, without.m_attempted_alignments);
}
