// Unit tests for protal's WFA2 wrapper: --x_drop reaches WFA2 (also in the copies alignment threads
// use), its default changes no alignment, and an alignment X-drop breaks is reported as failed
// rather than as a successful alignment that does not fit the read.
#include <gtest/gtest.h>
#include <random>
#include <string>
#include "Alignment/WFA2Wrapper2.h"

using namespace protal;

namespace {
    std::string Random(size_t n, std::mt19937& rng) {
        std::string s(n, 'A');
        for (auto& c : s) c = "ACGT"[rng() % 4];
        return s;
    }

    // Substitutions, deletions and insertions at `rate`.
    std::string Mutate(std::string const& s, double rate, std::mt19937& rng) {
        std::uniform_real_distribution<double> u(0, 1);
        std::string out;
        for (char c : s) {
            double const r = u(rng);
            if (r < rate * 0.8) out += c == 'A' ? 'C' : 'A';
            else if (r < rate * 0.9) {}
            else if (r < rate) { out += c; out += "ACGT"[rng() % 4]; }
            else out += c;
        }
        return out;
    }

    size_t ReadBases(std::string const& operations) {
        size_t n = 0;
        for (char c : operations) n += c != 'D';
        return n;
    }

    struct Case {
        std::string read, ref;
        int q_end_free, r_begin_free, r_end_free;
    };

    // Reads as protal aligns them: in a gene window with dovetails, some running past the gene's end.
    std::vector<Case> Cases(size_t n, std::mt19937& rng) {
        std::vector<Case> cases;
        std::uniform_real_distribution<double> rate(0.0, 0.15);
        for (size_t i = 0; i < n; i++) {
            size_t const len = 150;
            size_t const overhang = i % 3 == 0 ? rng() % 40 : 0;
            std::string const gene = Random(len + 18, rng);
            std::string read = Mutate(gene.substr(9, len - overhang), rate(rng), rng) + Random(overhang, rng);
            std::string const ref = overhang ? gene.substr(0, 9 + len - overhang) : gene;
            cases.push_back({ read, ref, overhang ? static_cast<int>(overhang) + 9 : 0, 18, overhang ? 0 : 18 });
        }
        return cases;
    }
}

// A read of a simulated sample (w900) in its gene window, as protal aligns it: the window's first 18 bases are free
// (dovetail), the rest must align. Without X-drop it aligns over the whole read. WFA2-lib v2.3.6 with X-drop 50 reports
// it completed, with operations for 151 read bases where the read has 150: the wrapper reports it failed, and so does a
// copy of it (as alignment threads use), which shows that the copy has the X-drop too.
TEST(WFA2Wrapper, AnAlignmentXDropCutsShortFails) {
    std::string const read = "TCTCACCACGTTAGAGCTTGCAGCTGCTTCCCGAGGTCCAGAGCGCCGCAGTGATTTGCAGCGAAGTAAACGCCAATCTGCTACACACC"
                             "GCAATTGTATGTGGGTGGATATTAACACTTATTACGTCTACGCCCTGCGTGGACTGTCACC";
    std::string const ref = "CCTCCGGGCTGTCACCACGCTAGAGCTTGCAGCTGCTTCCCGAGGTCCGGAGCGCCGCAGTGATTTGCAGCGAAGTAAACGCCAATCTGC"
                            "TCACCGCAATTGTATGTGGGTGGATATTAACACTTATTACGTCTACGCCCTGCGTGGACTGTCACCTGT";
    WFA2Wrapper2 uncut(4, 6, 2, 0);
    uncut.Alignment(read, ref, 0, 0, 18, 0, 61);
    ASSERT_TRUE(uncut.Success());
    EXPECT_EQ(ReadBases(uncut.Cigar()), read.size());

    wfa::WFAlignerGapAffine raw(4, 6, 2, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryHigh);
    raw.setHeuristicXDrop(50, 1);
    raw.setMaxAlignmentSteps(61);
    if (raw.alignEndsFree(ref, 18, 0, read, 0, 0) != wfa::WFAligner::StatusAlgCompleted || ReadBases(raw.getAlignment()) == read.size()) {
        GTEST_SKIP() << "this WFA2-lib no longer reports the alignment X-drop cut short as completed: the wrapper's check of "
                        "such alignments (WFA2Wrapper2::CoversText) has nothing to catch here";
    }
    WFA2Wrapper2 cut(4, 6, 2, 50);
    cut.Alignment(read, ref, 0, 0, 18, 0, 61);
    EXPECT_FALSE(cut.Success());
    WFA2Wrapper2 copy(cut);  // without the X-drop it would align as `uncut`; without the wrapper's check, succeed
    copy.Alignment(read, ref, 0, 0, 18, 0, 61);
    EXPECT_FALSE(copy.Success());
}

TEST(WFA2Wrapper, AlignmentsItReportsCoverTheRead) {
    std::mt19937 rng(11);
    auto const cases = Cases(3000, rng);
    for (size_t x_drop : { size_t{0}, size_t{1000}, size_t{50}, size_t{20} }) {
        WFA2Wrapper2 aligner(4, 6, 2, x_drop);
        WFA2Wrapper2 copy(aligner);  // each alignment thread works on a copy
        size_t aligned = 0;
        for (auto const& c : cases) {
            copy.Reset();
            copy.Alignment(c.read, c.ref, 0, c.q_end_free, c.r_begin_free, c.r_end_free, 61);
            if (!copy.Success()) continue;
            aligned++;
            EXPECT_EQ(ReadBases(copy.Cigar()), c.read.size()) << "x_drop " << x_drop;
        }
        EXPECT_GT(aligned, cases.size() / 2) << "x_drop " << x_drop;
    }
}

TEST(WFA2Wrapper, TheDefaultXDropChangesNoAlignment) {
    std::mt19937 rng(5);
    auto const cases = Cases(1000, rng);
    WFA2Wrapper2 off(4, 6, 2, 0), deflt(4, 6, 2, 1000);
    for (auto const& c : cases) {
        off.Reset();
        deflt.Reset();
        off.Alignment(c.read, c.ref, 0, c.q_end_free, c.r_begin_free, c.r_end_free, 61);
        deflt.Alignment(c.read, c.ref, 0, c.q_end_free, c.r_begin_free, c.r_end_free, 61);
        ASSERT_EQ(off.Success(), deflt.Success());
        if (off.Success()) {
            EXPECT_EQ(off.Cigar(), deflt.Cigar());
            EXPECT_EQ(off.GetAlignmentScore(), deflt.GetAlignmentScore());
        }
    }
}
