// Unit tests for aligning a read from its anchor's exact matches (AnchoredAligner, as
// SimpleAlignmentHandler::AlignAnchor uses it) against aligning the whole read into its window:
// the same alignments where the best path runs through the anchor, the whole-window alignment for
// chains the anchored one does not handle, and valid alignments throughout.
#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <memory>
#include <random>
#include <string>
#include <vector>
#include <unistd.h>
#include "Core/AlignmentStrategy.h"

using namespace protal;

namespace {
    constexpr size_t kGenes = 20, kGeneLength = 1500;

    // Random genes of taxid 1 (genes 1..kGenes) in a temporary directory, loaded.
    struct RandomReference {
        std::vector<std::string> genes;
        std::filesystem::path dir;
        std::unique_ptr<GenomeLoader> loader;

        RandomReference() {
            std::mt19937 rng(23);
            dir = std::filesystem::temp_directory_path() / ("protal_anchored_test_" + std::to_string(::getpid()));
            std::filesystem::create_directories(dir);
            std::ofstream fna(dir / "reference.fna"), map(dir / "reference.map");
            size_t offset = 0;
            for (size_t id = 1; id <= kGenes; id++) {
                std::string seq(kGeneLength, 'A');
                for (auto& c : seq) c = "ACGT"[rng() % 4];
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

    // Both methods on the same reference, as protal sets them up (-a 0.9).
    struct Handlers {
        WFA2Wrapper2 aligner{4, 6, 2, 1000};
        SimpleAlignmentHandler whole, anchored;
        explicit Handlers(GenomeLoader& loader) :
                whole(loader, aligner, 31, 3, 0.9, false), anchored(loader, aligner, 31, 3, 0.9, false) {
            whole.SetAnchoredAlignment(false);
            anchored.SetAnchoredAlignment(true);
        }
    };

    struct Outcome {
        bool aligned = false;
        int score = 0;
        int start = 0;
        std::string cigar;
    };

    Outcome Align(SimpleAlignmentHandler& handler, uint32_t gene, std::string read, ChainList chain) {
        ChainAlignmentAnchor anchor(1, gene, true);
        anchor.chain = std::move(chain);
        AlignmentResult result;
        std::string rev = KmerUtils::ReverseComplement(read), id = "r";
        Outcome o;
        o.aligned = handler.AlignAnchor(anchor, result, read, rev, false, id);
        if (o.aligned) {
            o.score = result.AlignmentScore();
            o.start = result.GetAlignmentInfo().gene_alignment_start;
            o.cigar = result.GetAlignmentInfo().cigar;
        }
        return o;
    }

    // The maximal exact runs of at least min_length on the diagonal (gene position - read position).
    ChainList ExactRuns(std::string const& read, std::string const& gene, long diagonal, size_t min_length = 15) {
        ChainList runs;
        size_t i = 0;
        while (i < read.size()) {
            long const g = diagonal + static_cast<long>(i);
            if (g < 0 || g >= static_cast<long>(gene.size()) || read[i] != gene[g]) { i++; continue; }
            size_t j = i;
            while (j < read.size() && diagonal + static_cast<long>(j) < static_cast<long>(gene.size()) && read[j] == gene[diagonal + j]) j++;
            if (j - i >= min_length) runs.emplace_back(static_cast<uint32_t>(diagonal + static_cast<long>(i)), static_cast<uint16_t>(i), static_cast<uint16_t>(j - i));
            i = j;
        }
        return runs;
    }

    void ExpectSame(Outcome const& whole, Outcome const& anchored) {
        EXPECT_EQ(anchored.aligned, whole.aligned);
        EXPECT_EQ(anchored.score, whole.score);
        EXPECT_EQ(anchored.start, whole.start);
        EXPECT_EQ(anchored.cigar, whole.cigar);
    }
}

TEST(AnchoredAlignment, AReadInsideItsGene) {
    RandomReference ref;
    Handlers h(*ref.loader);
    std::string read = ref.genes[0].substr(400, 150);
    read[20] = read[20] == 'A' ? 'C' : 'A';  // mismatches left and right of the link
    read[130] = read[130] == 'A' ? 'C' : 'A';
    ChainList chain = { ChainLink(421, 21, 109) };
    auto whole = Align(h.whole, 1, read, chain), anchored = Align(h.anchored, 1, read, chain);
    ASSERT_TRUE(anchored.aligned);
    ExpectSame(whole, anchored);
    EXPECT_EQ(anchored.start, 400);
    EXPECT_EQ(h.anchored.m_anchored_alignments, 1u);
    EXPECT_EQ(h.whole.m_whole_window_alignments, 1u);
}

TEST(AnchoredAlignment, ReadsRunningPastTheGeneEnds) {
    RandomReference ref;
    Handlers h(*ref.loader);
    std::mt19937 rng(3);
    auto random = [&rng](size_t n) { std::string s(n, 'A'); for (auto& c : s) c = "ACGT"[rng() % 4]; return s; };
    // 20 bases before the gene's start, 30 past its end: soft-clipped by both methods alike.
    std::string before = random(20) + ref.genes[1].substr(0, 130);
    std::string after = ref.genes[1].substr(kGeneLength - 120) + random(30);
    auto w1 = Align(h.whole, 2, before, ExactRuns(before, ref.genes[1], -20));
    auto a1 = Align(h.anchored, 2, before, ExactRuns(before, ref.genes[1], -20));
    auto w2 = Align(h.whole, 2, after, ExactRuns(after, ref.genes[1], kGeneLength - 120));
    auto a2 = Align(h.anchored, 2, after, ExactRuns(after, ref.genes[1], kGeneLength - 120));
    ASSERT_TRUE(a1.aligned && a2.aligned);
    ExpectSame(w1, a1);
    ExpectSame(w2, a2);
    EXPECT_EQ(a1.cigar.find_first_not_of('S'), 20u);
    EXPECT_EQ(a2.cigar.size() - 1 - a2.cigar.find_last_not_of('S'), 30u);
    EXPECT_EQ(h.anchored.m_anchored_alignments, 2u);
}

TEST(AnchoredAlignment, GapsInEitherFlank) {
    RandomReference ref;
    Handlers h(*ref.loader);
    std::string const& gene = ref.genes[2];
    // 3 bases deleted from the read 10 bases in, 2 inserted 20 bases before its end; the link in between.
    std::string read = gene.substr(600, 10) + gene.substr(613, 110) + "GT" + gene.substr(723, 28);
    ASSERT_EQ(read.size(), 150u);
    ChainList chain = ExactRuns(read, gene, 603);  // the diagonal after the deletion
    ASSERT_FALSE(chain.empty());
    ChainList link = { chain.front() };
    auto whole = Align(h.whole, 3, read, link), anchored = Align(h.anchored, 3, read, link);
    ASSERT_TRUE(anchored.aligned);
    EXPECT_EQ(anchored.score, whole.score);
    EXPECT_NE(anchored.cigar.find('D'), std::string::npos);
    EXPECT_NE(anchored.cigar.find('I'), std::string::npos);
}

TEST(AnchoredAlignment, AnNInsideTheLinkCountsAsAMismatch) {
    RandomReference ref;
    Handlers h(*ref.loader);
    std::string read = ref.genes[3].substr(200, 150);
    read[75] = 'N';  // protal's seeds read N as A, so a link may span it
    ChainList chain = { ChainLink(200, 0, 150) };
    auto whole = Align(h.whole, 4, read, chain), anchored = Align(h.anchored, 4, read, chain);
    ASSERT_TRUE(anchored.aligned);
    ExpectSame(whole, anchored);
    EXPECT_EQ(anchored.cigar[75], 'X');
}

TEST(AnchoredAlignment, LinksOnOneDiagonalAndBetweenThem) {
    RandomReference ref;
    Handlers h(*ref.loader);
    std::string read = ref.genes[4].substr(300, 150);
    for (size_t i : { 40u, 41u, 90u, 92u, 94u, 96u, 98u }) read[i] = read[i] == 'A' ? 'C' : 'A';  // 2, then 5 mismatches between links
    ChainList chain = ExactRuns(read, ref.genes[4], 300);
    ASSERT_GE(chain.size(), 3u);
    auto whole = Align(h.whole, 5, read, chain), anchored = Align(h.anchored, 5, read, chain);
    ASSERT_TRUE(anchored.aligned);
    ExpectSame(whole, anchored);
}

TEST(AnchoredAlignment, SeedsThatImplyAnIndelAreAlignedAsAWhole) {
    RandomReference ref;
    Handlers h(*ref.loader);
    std::string const& gene = ref.genes[5];
    std::string read = gene.substr(500, 70) + gene.substr(574, 80);  // 4 bases deleted in the middle
    ChainList chain = { ExactRuns(read, gene, 500).front(), ExactRuns(read, gene, 504).back() };
    auto whole = Align(h.whole, 6, read, chain), anchored = Align(h.anchored, 6, read, chain);
    ExpectSame(whole, anchored);
    EXPECT_EQ(h.anchored.m_anchored_alignments, 0u);
    EXPECT_EQ(h.anchored.m_whole_window_alignments, 1u);
}

// A flank with at most one mismatch on its link's diagonal is aligned without WFA2 (AnchoredAligner::Flank): the
// same status and operations as WFA2 gives, for reads starting and ending in or past the window, free ends of either
// side and budgets down to the edge, with mismatches and Ns anywhere in the flanks.
TEST(AnchoredAlignment, UngappedFlanksAreWhatWFA2Gives) {
    std::mt19937 rng(7);
    auto base = [&rng]() { return "ACGT"[rng() % 4]; };
    WFA2Wrapper2 aligner{4, 6, 2, 1000};
    AnchoredAligner fast, wfa;
    wfa.SetUngappedFlanks(false);
    size_t aligned = 0, failed = 0;
    std::string fast_ops, wfa_ops;
    for (int n = 0; n < 20000; n++) {
        // The window is gene [20, 20 + length); the gene has 20 bases more at either end, which a read may start or
        // end in (bases past the window). The link: 15-30 exact bases in the middle of the read.
        size_t const length = 60 + rng() % 120;
        std::string gene(length + 40, 'A');
        for (auto& c : gene) c = base();
        size_t const start = 20 + rng() % 25 - 12, end = 20 + length + rng() % 25 - 12;
        std::string read = gene.substr(start, end - start);
        size_t const link_length = 15 + rng() % 16, link_read = (read.size() - link_length) / 2;
        for (int i = static_cast<int>(rng() % 5); i > 0; i--) {  // up to 4 changes, in the flanks only
            size_t const p = rng() % read.size();
            if (p >= link_read && p < link_read + link_length) continue;
            read[p] = rng() % 8 == 0 ? 'N' : base();
        }
        AlignmentWindow w;
        w.ref_start = 20;
        w.ref_end = 20 + length;
        w.ref_begin_free = static_cast<int>(rng() % 15);
        w.ref_end_free = static_cast<int>(rng() % 15);
        w.read_begin_free = static_cast<int>(rng() % 15);
        w.read_end_free = static_cast<int>(rng() % 15);
        w.max_score = 1 + static_cast<int>(rng() % 30);
        ChainList chain = { ChainLink(static_cast<uint32_t>(start + link_read), static_cast<uint16_t>(link_read), static_cast<uint16_t>(link_length)) };
        auto const a = fast.Align(read, gene, chain, w, aligner, fast_ops);
        auto const b = wfa.Align(read, gene, chain, w, aligner, wfa_ops);
        ASSERT_EQ(a, b) << "case " << n << ": " << read;
        if (a == AnchoredAligner::Status::Aligned) {
            ASSERT_EQ(fast_ops, wfa_ops) << "case " << n << ": " << read;
            aligned++;
        }
        failed += a == AnchoredAligner::Status::Failed;
    }
    std::cout << aligned << " aligned, " << failed << " failed, " << fast.UngappedFlanks() << " flanks without WFA2" << std::endl;
    EXPECT_GT(aligned, 5000u);
    EXPECT_GT(failed, 1000u);
    EXPECT_GT(fast.UngappedFlanks(), 10000u);
    EXPECT_EQ(wfa.UngappedFlanks(), 0u);
}

// A flank with 2 or 3 mismatches on its link's diagonal is aligned without WFA2 where no alignment with a gap costs as
// much (AnchoredAligner::NoCheaperGap): the same status and operations as WFA2 gives, with up to 8 changes (Ns
// included) and sometimes an indel in the read, free ends of either side, budgets down to the edge, and genes of random
// bases or of short repeats (where a gap can cost the same as the mismatches).
TEST(AnchoredAlignment, FlanksWithTwoOrThreeMismatchesAreWhatWFA2Gives) {
    std::mt19937 rng(11);
    auto base = [&rng]() { return "ACGT"[rng() % 4]; };
    WFA2Wrapper2 aligner{4, 6, 2, 1000};
    AnchoredAligner fast, wfa;
    wfa.SetUngappedFlanks(false);
    size_t aligned = 0, failed = 0;
    std::string fast_ops, wfa_ops;
    for (int n = 0; n < 30000; n++) {
        size_t const length = 60 + rng() % 120;
        std::string gene(length + 40, 'A');
        size_t const unit = rng() % 3 == 0 ? 1 + rng() % 4 : 0;  // a third of the genes repeat a unit of 1-4 bases, with changes
        for (size_t i = 0; i < gene.size(); i++) gene[i] = unit == 0 || rng() % 10 == 0 ? base() : gene[i - std::min(i, unit)];
        if (unit > 0) for (size_t i = 0; i < unit; i++) gene[i] = base();
        std::string read = gene.substr(20, length);  // the window is the gene [20, 20 + length): read and window equal
        size_t const link_length = 15 + rng() % 16, link_read = (read.size() - link_length) / 2;
        for (int i = static_cast<int>(rng() % 9); i > 0; i--) {
            size_t const p = rng() % read.size();
            if (p >= link_read && p < link_read + link_length) continue;
            read[p] = rng() % 8 == 0 ? 'N' : base();
        }
        size_t link_in_read = link_read;
        if (rng() % 4 == 0) {  // a deletion from the read, outside the link: the flanks then differ in length or need a gap
            size_t const p = rng() % read.size();
            if (p < link_read) { read.erase(p, 1); link_in_read--; }
            else if (p >= link_read + link_length) read.erase(p, 1);
        }
        AlignmentWindow w;
        w.ref_start = 20;
        w.ref_end = 20 + length;
        if (rng() % 2) {
            w.ref_begin_free = static_cast<int>(rng() % 12);
            w.ref_end_free = static_cast<int>(rng() % 12);
            w.read_begin_free = static_cast<int>(rng() % 12);
            w.read_end_free = static_cast<int>(rng() % 12);
        }
        w.max_score = 1 + static_cast<int>(rng() % 40);
        if (read.size() < link_read + link_length) continue;
        ChainList chain = { ChainLink(static_cast<uint32_t>(20 + link_read), static_cast<uint16_t>(link_in_read), static_cast<uint16_t>(link_length)) };
        auto const a = fast.Align(read, gene, chain, w, aligner, fast_ops);
        auto const b = wfa.Align(read, gene, chain, w, aligner, wfa_ops);
        ASSERT_EQ(a, b) << "case " << n << ": " << read;
        if (a == AnchoredAligner::Status::Aligned) {
            ASSERT_EQ(fast_ops, wfa_ops) << "case " << n << ": " << read;
            aligned++;
        }
        failed += a == AnchoredAligner::Status::Failed;
    }
    std::cout << aligned << " aligned, " << failed << " failed, " << fast.UngappedFlanks() << " flanks without WFA2" << std::endl;
    EXPECT_GT(aligned, 5000u);
    EXPECT_GT(failed, 1000u);
    EXPECT_GT(fast.UngappedFlanks(), 20000u);
}

namespace {
    // A long read's window as LongReadAligner aligns it: 100 random bases, the whole gene with ONT-like errors
    // (substitutions, deletions and insertions at the given rates), 100 random bases. Returns the read and its
    // exact runs of at least 20 bases along the true path, on whatever diagonal the indels put them.
    std::pair<std::string, ChainList> LongReadOf(std::string const& gene, std::mt19937& rng, double sub, double del, double ins) {
        std::uniform_real_distribution<double> u(0, 1);
        auto base = [&rng]() { return "ACGT"[rng() % 4]; };
        std::string read;
        std::vector<long> gene_of;  // per read base, its gene position, -1 for an inserted or random base
        for (int i = 0; i < 100; i++) { read += base(); gene_of.push_back(-1); }
        for (size_t g = 0; g < gene.size(); g++) {
            double const r = u(rng);
            if (r < del) continue;
            if (r < del + ins) { read += base(); gene_of.push_back(-1); }
            char c = gene[g];
            if (u(rng) < sub) c = "ACGT"[(std::string("ACGT").find(c) + 1 + rng() % 3) % 4];
            read += c;
            gene_of.push_back(static_cast<long>(g));
        }
        for (int i = 0; i < 100; i++) { read += base(); gene_of.push_back(-1); }
        ChainList runs;
        size_t i = 0;
        while (i < read.size()) {
            if (gene_of[i] < 0 || read[i] != gene[gene_of[i]]) { i++; continue; }
            size_t j = i + 1;
            while (j < read.size() && gene_of[j] == gene_of[j - 1] + 1 && read[j] == gene[gene_of[j]]) j++;
            if (j - i >= 20) runs.emplace_back(static_cast<uint32_t>(gene_of[i]), static_cast<uint16_t>(i), static_cast<uint16_t>(j - i));
            i = j;
        }
        return { read, runs };
    }

    // The runs within 6 diagonals of the first, as the anchor finder chains a long read's seeds
    // (ChainAnchorFinder::FindAnchorsSingleRef).
    ChainList AsChained(ChainList const& runs) {
        ChainList chain;
        if (runs.empty()) return chain;
        long const d0 = static_cast<long>(runs.front().genepos) - static_cast<long>(runs.front().readpos);
        for (auto const& run : runs) {
            if (std::abs(static_cast<long>(run.genepos) - static_cast<long>(run.readpos) - d0) <= 6) chain.push_back(run);
        }
        return chain;
    }

    // Both methods as protal sets them up for ONT reads (-a 0.85, no X-drop).
    struct LongReadHandlers {
        WFA2Wrapper2 aligner{4, 6, 2, 0};
        SimpleAlignmentHandler whole, anchored;
        explicit LongReadHandlers(GenomeLoader& loader) :
                whole(loader, aligner, 31, 3, 0.85, false), anchored(loader, aligner, 31, 3, 0.85, false) {
            whole.SetAnchoredAlignment(false);
            anchored.SetAnchoredAlignment(true);
            anchored.SetAnchoredIndels(true);
        }
    };
}

// Long reads through every link of their chain (AnchoredAligner::AllowIndels): seeds on diagonals their indels shift,
// a chain that covers part of the gene, more links found in the rest (Reseed). Alignments as good as aligning the
// whole read into the whole window, and valid.
TEST(AnchoredAlignment, LongReadsThroughTheirChain) {
    RandomReference ref;
    LongReadHandlers h(*ref.loader);
    std::mt19937 rng(17);
    size_t cases = 0, both = 0, worse = 0, much_worse = 0;
    long whole_sum = 0, anchored_sum = 0;
    for (int n = 0; n < 60; n++) {
        uint32_t const gene_id = 1 + rng() % kGenes;
        double const rate = n % 2 ? 0.02 : 0.005;  // ONT-like and HiFi-like reads
        auto [read, runs] = LongReadOf(ref.genes[gene_id - 1], rng, rate, rate, rate * 0.7);
        ChainList const chain = AsChained(runs);
        if (chain.empty()) continue;
        cases++;
        auto whole = Align(h.whole, gene_id, read, chain), anchored = Align(h.anchored, gene_id, read, chain);
        if (anchored.aligned) {
            EXPECT_EQ(std::count_if(anchored.cigar.begin(), anchored.cigar.end(), [](char c) { return c != 'D'; }),
                      static_cast<long>(read.size()));
        }
        if (!whole.aligned || !anchored.aligned) continue;
        both++;
        whole_sum += whole.score;
        anchored_sum += anchored.score;
        worse += anchored.score < whole.score;
        much_worse += anchored.score < whole.score - std::abs(whole.score) / 50;
    }
    std::cout << cases << " long reads: both aligned " << both << ", anchored worse " << worse << " (by over 2%: " << much_worse
              << "), scores " << anchored_sum << " anchored, " << whole_sum << " whole; through the chain "
              << h.anchored.m_anchored_alignments << ", as whole windows " << h.anchored.m_whole_window_alignments << std::endl;
    ASSERT_GE(cases, 50u);
    EXPECT_GE(both, cases - 2);
    EXPECT_LE(much_worse, 1u);
    EXPECT_GE(static_cast<double>(anchored_sum), static_cast<double>(whole_sum) - 0.005 * std::abs(static_cast<double>(whole_sum)));
    // The window's right end comes from the last link's diagonal for chains with indels (SimpleAlignmentHandler::AlignAnchor,
    // AlignmentOrientation::Update): the read's bases past the gene's end are free by the end's own diagonal. Before, the first
    // link's diagonal placed the end, and where the indels had moved it by more than the 9 bases of dovetail the chain's
    // right flank did not fit and the whole-window alignment took the read (9 of 60 here; 4 now, chains the anchored aligner
    // does not handle for other reasons).
    EXPECT_EQ(h.anchored.m_anchored_alignments + h.anchored.m_whole_window_alignments, cases);
    EXPECT_GE(10 * h.anchored.m_anchored_alignments, 9 * cases);
}

// Links that overlap (exact seeds on two diagonals reach into a homopolymer that lost bases) are cut where they
// overlap, and the bases between them on different diagonals become the indel.
TEST(AnchoredAlignment, OverlappingLinksOfALongReadAreCut) {
    std::mt19937 rng(9);
    std::string gene(609, 'A');
    for (auto& c : gene) c = "CGT"[rng() % 3];
    for (size_t i = 300; i < 309; i++) gene[i] = 'A';           // 9 As
    std::string const read = gene.substr(0, 303) + gene.substr(306);  // 3 of them lost: 6 As
    // Exact on diagonal 0 up to the read's 6th A (read 306), and on diagonal 3 from its 1st A (read 300).
    ChainList chain = { ChainLink(0, 0, 306), ChainLink(303, 300, 306) };
    for (auto const& link : chain) ASSERT_EQ(read.substr(link.readpos, link.length), gene.substr(link.genepos, link.length));
    AlignmentWindow w;
    w.ref_start = 0;
    w.ref_end = gene.size();
    w.max_score = 100;
    WFA2Wrapper2 aligner{4, 6, 2, 0};
    AnchoredAligner anchored;
    anchored.AllowIndels(true);
    std::string ops;
    ASSERT_EQ(anchored.Align(read, gene, chain, w, aligner, ops), AnchoredAligner::Status::Aligned);
    EXPECT_EQ(ops, std::string(306, 'M') + "DDD" + std::string(300, 'M'));
    // Without indels allowed, the same chain is left to the whole-window alignment.
    AnchoredAligner short_reads;
    EXPECT_EQ(short_reads.Align(read, gene, chain, w, aligner, ops), AnchoredAligner::Status::NotApplicable);
}

// Simulated reads as protal meets them: from their gene or a relative (up to 15% divergence), with
// indels, Ns, and running past gene ends, anchored at an exact run on their diagonal or at several.
TEST(AnchoredAlignment, AgreesWithTheWholeReadAlignment) {
    RandomReference ref;
    Handlers h(*ref.loader);
    std::mt19937 rng(41);
    std::uniform_real_distribution<double> u(0, 1);
    size_t cases = 0, both = 0, same_score = 0, same_cigar = 0, only_whole = 0, only_anchored = 0, anchored_better = 0, whole_better = 0;
    for (int n = 0; n < 6000; n++) {
        uint32_t const gene_id = 1 + rng() % kGenes;
        std::string const& gene = ref.genes[gene_id - 1];
        long const start = static_cast<long>(rng() % (kGeneLength + 60)) - 40;  // some reads run past either end
        double const divergence = u(rng) < 0.5 ? u(rng) * 0.03 : u(rng) * 0.15;
        std::string read;
        for (long i = start; static_cast<long>(read.size()) < 150; i++) {
            char c = i >= 0 && i < static_cast<long>(gene.size()) ? gene[i] : "ACGT"[rng() % 4];
            double const r = u(rng);
            if (r < divergence * 0.9) c = "ACGT"[(std::string("ACGT").find(c) + 1 + rng() % 3) % 4];
            else if (r < divergence * 0.95) continue;                                                    // deletion
            else if (r < divergence) { read += "ACGT"[rng() % 4]; if (read.size() == 150) break; }       // insertion
            if (u(rng) < 0.002) c = 'N';
            read += c;
        }
        ChainList runs = ExactRuns(read, gene, start);
        if (runs.empty()) continue;
        ChainList chain;
        if (rng() % 2) chain = { runs[rng() % runs.size()] };  // one link
        else chain = runs;                                     // every run on the diagonal
        auto whole = Align(h.whole, gene_id, read, chain), anchored = Align(h.anchored, gene_id, read, chain);
        cases++;
        if (anchored.aligned) {
            EXPECT_EQ(std::count_if(anchored.cigar.begin(), anchored.cigar.end(), [](char c) { return c != 'D'; }), 150);
        }
        if (whole.aligned && anchored.aligned) {
            both++;
            same_score += whole.score == anchored.score;
            same_cigar += whole.cigar == anchored.cigar && whole.start == anchored.start;
            anchored_better += anchored.score > whole.score;
            whole_better += whole.score > anchored.score;
        } else if (whole.aligned) {
            only_whole++;
        } else if (anchored.aligned) {
            only_anchored++;
        }
    }
    RecordProperty("cases", static_cast<int>(cases));
    std::cout << cases << " anchors: both aligned " << both << ", same score " << same_score << ", same alignment " << same_cigar
              << ", anchored better " << anchored_better << ", whole better " << whole_better << ", only whole " << only_whole
              << ", only anchored " << only_anchored << "; anchored " << h.anchored.m_anchored_alignments << ", whole "
              << h.anchored.m_whole_window_alignments << std::endl;
    ASSERT_GT(cases, 4000u);
    EXPECT_GT(both, cases / 3);
    EXPECT_GE(static_cast<double>(same_score), 0.99 * static_cast<double>(both));
    EXPECT_GE(static_cast<double>(same_cigar), 0.98 * static_cast<double>(both));
    EXPECT_LE(only_whole + only_anchored, cases / 200);
    EXPECT_GT(h.anchored.m_anchored_alignments, cases / 2);
}

// The alignment handler takes the reverse complement of the read from its caller (the anchor finder
// has it), or makes it itself: either way the same alignments, for reads of either strand, and the
// caller's string is neither changed nor copied per anchor.
TEST(SimpleAlignmentHandler, ACallersReverseComplementGivesTheSameAlignments) {
    RandomReference ref;
    Handlers own(*ref.loader), given(*ref.loader);
    std::mt19937 rng(5);
    size_t aligned = 0;
    for (bool reverse : { false, true }) {
        for (int round = 0; round < 20; round++) {
            uint32_t const gene = 1 + rng() % kGenes;
            size_t const start = 100 + rng() % 1000;
            std::string segment = ref.genes[gene - 1].substr(start, 150);
            for (int m = 0; m < 3; m++) {
                char& c = segment[10 + rng() % 130];
                c = c == 'A' ? 'C' : 'A';
            }
            std::string const read = reverse ? KmerUtils::ReverseComplement(segment) : segment;
            ChainAlignmentAnchor anchor(1, gene, !reverse);
            // the chain is in the orientation of the alignment, the gene's strand: AlignAnchor aligns `rev` of a
            // reverse anchor, which is the segment
            anchor.chain = ExactRuns(segment, ref.genes[gene - 1], static_cast<long>(start));
            ASSERT_FALSE(anchor.chain.empty());
            std::string const rev = KmerUtils::ReverseComplement(read);
            std::string const read_before = read, rev_before = rev;
            std::string header = "r";
            AlignmentAnchorList a, b;
            a.push_back(anchor);
            b.push_back(anchor);
            AlignmentResultList from_own, from_given;
            own.anchored(a, from_own, read, 3, header);
            given.anchored(b, from_given, read, rev, 3, header);
            ASSERT_EQ(from_own.size(), from_given.size());
            for (size_t i = 0; i < from_own.size(); i++) {
                EXPECT_EQ(from_own[i].AlignmentScore(), from_given[i].AlignmentScore());
                EXPECT_EQ(from_own[i].GetAlignmentInfo().cigar, from_given[i].GetAlignmentInfo().cigar);
                EXPECT_EQ(from_own[i].GetAlignmentInfo().gene_alignment_start, from_given[i].GetAlignmentInfo().gene_alignment_start);
                EXPECT_EQ(from_own[i].Forward(), from_given[i].Forward());
            }
            aligned += from_given.size();
            EXPECT_EQ(read, read_before);
            EXPECT_EQ(rev, rev_before);
        }
    }
    EXPECT_GT(aligned, 30u);
}
