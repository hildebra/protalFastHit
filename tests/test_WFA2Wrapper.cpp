// Unit tests for protal's WFA2 wrapper: --x_drop reaches WFA2 (also in the copies alignment threads
// use), its default changes no alignment, and an alignment X-drop breaks is reported as failed
// rather than as a successful alignment that does not fit the read. The Singletrack traceback
// (--wfa_traceback) gives the classic traceback's outcomes and scores with operations that are a
// valid alignment of that score, and build_settings.tsv records it. The wrapper's halved penalties,
// and X-drop left out where it cannot prune, give WFA2's alignments with the penalties as given.
#include <gtest/gtest.h>
#include <cmath>
#include <random>
#include <sstream>
#include <string>
#include "Alignment/WFA2Wrapper2.h"
#include "Utilities/BuildSettings.h"

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

    // Long reads as protal aligns them: a 1-3 kb gene in a window with ~160 free bases at either end, the read part at up
    // to 15% divergence (ONT's identity floor), some of it unrelated sequence.
    std::vector<Case> LongCases(size_t n, std::mt19937& rng) {
        std::vector<Case> cases;
        std::uniform_real_distribution<double> rate(0.0, 0.15);
        for (size_t i = 0; i < n; i++) {
            size_t const len = 1000 + rng() % 2000;
            std::string const window = Random(len + 320, rng);
            std::string const read = i % 5 == 4 ? Random(len, rng) : Mutate(window.substr(160, len), rate(rng), rng);
            cases.push_back({ read, window, 0, 200, 200 });
        }
        return cases;
    }

    // The gap-affine score (mismatch 4, a gap of L bases 6 + 2L, a match 0) of `ops` aligning `read` (text) to `ref`
    // (pattern), the reference's first `r_begin_free` and last `r_end_free` bases and the read's last `q_end_free` free
    // (WFA2's ends-free: a free end is a leading or trailing run of D or I, paid beyond its allowance as a gap of the
    // rest). -1 if `ops` is no alignment of the two: an M of different bases, an X of equal ones, or bases left over.
    int AlignmentScore(std::string const& ops, std::string const& read, std::string const& ref, int r_begin_free, int r_end_free,
                       int q_end_free) {
        size_t r = 0, q = 0;
        for (char const op : ops) {
            if (op == 'M' || op == 'X') {
                if (r >= ref.size() || q >= read.size() || (ref[r] == read[q]) != (op == 'M')) return -1;
                r++, q++;
            } else if (op == 'D') {
                if (r++ >= ref.size()) return -1;
            } else if (op == 'I') {
                if (q++ >= read.size()) return -1;
            } else {
                return -1;
            }
        }
        if (r != ref.size() || q != read.size()) return -1;
        auto run = [](std::string const& s, size_t from, size_t to) { size_t n = from; while (n < to && s[n] == s[from]) n++; return n - from; };
        size_t const lead = ops.empty() || (ops[0] != 'D' && ops[0] != 'I') ? 0 : run(ops, 0, ops.size());
        std::string const rev(ops.rbegin(), ops.rend());
        size_t const trail = lead == ops.size() || rev.empty() || (rev[0] != 'D' && rev[0] != 'I') ? 0 : run(rev, 0, rev.size());
        auto paid = [](size_t length, int free) { return length <= static_cast<size_t>(free) ? 0 : 6 + 2 * static_cast<int>(length - free); };
        int score = 0;
        if (lead) score += paid(lead, ops[0] == 'D' ? r_begin_free : 0);
        if (trail) score += paid(trail, rev[0] == 'D' ? r_end_free : q_end_free);
        for (size_t i = lead; i < ops.size() - trail;) {
            if (ops[i] == 'M') { i++; continue; }
            if (ops[i] == 'X') { score += 4; i++; continue; }
            size_t const n = run(ops, i, ops.size() - trail);
            score += 6 + 2 * static_cast<int>(n);
            i += n;
        }
        return score;
    }

    // Classic and Singletrack aligners on `cases` with a budget each (`max_score`, or with 0 the long reads' 0.15 of the read
    // as mismatches): the same outcome and score, Singletrack's operations an alignment of that score. -> how many
    // alignments succeeded and how many of them have other operations.
    std::pair<size_t, size_t> CompareTracebacks(std::vector<Case> const& cases, size_t x_drop, int max_score) {
        WFA2Wrapper2 classic(4, 6, 2, x_drop, WfaTraceback::Classic);
        WFA2Wrapper2 single(4, 6, 2, x_drop, WfaTraceback::Singletrack);
        size_t aligned = 0, other = 0;
        for (auto const& c : cases) {
            int const budget = max_score ? max_score : static_cast<int>(std::ceil(0.15 * static_cast<double>(c.read.size()))) * 4 + 1;
            classic.Reset();
            single.Reset();
            classic.Alignment(c.read, c.ref, 0, c.q_end_free, c.r_begin_free, c.r_end_free, budget);
            single.Alignment(c.read, c.ref, 0, c.q_end_free, c.r_begin_free, c.r_end_free, budget);
            EXPECT_EQ(classic.Success(), single.Success()) << "x_drop " << x_drop << ", read of " << c.read.size();
            if (!classic.Success() || !single.Success()) continue;
            aligned++;
            EXPECT_EQ(classic.GetAlignmentScore(), single.GetAlignmentScore());
            // The operations are an alignment, of at most the score (AlignmentScore reads a gap next to a free end as
            // using the free allowance first: a path that ends on a paid gap before a shorter free end costs WFA2 more).
            int const classic_ops = AlignmentScore(classic.Cigar(), c.read, c.ref, c.r_begin_free, c.r_end_free, c.q_end_free);
            int const single_ops = AlignmentScore(single.Cigar(), c.read, c.ref, c.r_begin_free, c.r_end_free, c.q_end_free);
            EXPECT_GE(classic_ops, 0) << classic.Cigar() << "\n" << c.ref << "\n" << c.read;
            EXPECT_GE(single_ops, 0) << single.Cigar() << "\n" << c.ref << "\n" << c.read;
            EXPECT_LE(classic_ops, -classic.GetAlignmentScore());
            EXPECT_LE(single_ops, -single.GetAlignmentScore());
            other += classic.Cigar() != single.Cigar();
        }
        return { aligned, other };
    }
}

// Short reads with dovetails (budget of 150 bp at 0.85), with and without X-drop: Singletrack aligns what the classic
// traceback aligns, with the same score, and its operations are an alignment of that score.
TEST(WfaTraceback, SingletrackGivesTheClassicScoresOnShortReads) {
    std::mt19937 rng(17);
    auto const cases = Cases(3000, rng);
    for (size_t x_drop : { size_t{0}, size_t{1000} }) {
        auto const [aligned, other] = CompareTracebacks(cases, x_drop, 61);
        EXPECT_GT(aligned, cases.size() / 2);
        RecordProperty("short_x" + std::to_string(x_drop) + "_other_operations", std::to_string(other) + " of " + std::to_string(aligned));
    }
}

// Long reads in gene-long windows at ONT's identity floor (budget 0.15 of the read, as mismatches), a fifth unrelated.
TEST(WfaTraceback, SingletrackGivesTheClassicScoresOnLongReads) {
    std::mt19937 rng(23);
    auto const cases = LongCases(60, rng);
    auto const [aligned, other] = CompareTracebacks(cases, 0, 0);
    EXPECT_GT(aligned, cases.size() / 2);
    EXPECT_LT(aligned, cases.size());  // the unrelated reads fail
    RecordProperty("long_other_operations", std::to_string(other) + " of " + std::to_string(aligned));
}

// A PacBio read's flank (r226-like simulated sample, realworld_s1): the read's end free, the reference's start not, a poor
// alignment that wf-adaptive prunes. The fork's Singletrack backtrace found no path (a 3-base deletion opened from an
// M-wavefront element outside its pruned bounds, which the forward pass still read) and ended the process
// (exit(-1)); with protal's patch it traces the classic score's alignment.
TEST(WfaTraceback, SingletrackFollowsPathsThroughPrunedWavefronts) {
    std::string const ref = "AACCTACTAGAGCACTTACCCGACCTAGACGACTGCGTCAAGTCCGCCGTACCTTCCTGCCAACCGGAACATGAGGCCCATTTCGAGCAGCGCTTTCC"
                            "AGAACGGGAACCCCAGCTTCTGCTACCGGTCCCTCGTCGGCCCTATCAGGAACCGGTCGAACGCGTCGTACCGGTGGTGCGGGGAGTGCCGTGGCCAC"
                            "CCTTTCGGCAGGCACCACCAGTCGCCGCCGTCGCGTACCTTCCGCGCTTTATCCATCCTACTCGGCAGGTGGTCTTATTACCCAAACAGCAAGCGGCG"
                            "GCGCAGACGGTGGTA";
    std::string const read = "CCGCTACCCCCGACGGCCGTCCCGCATCCGCTTGTTGACTGGCGCTGAGTCTACGCCCGCCAGGCCGCCGCGCCTAGCGCAGTACGCAGTCAAGGAGGA"
                             "CCAGCCGATGTTGCCGTCAGTCTGGCCGGAAACCGGGCGGACCGACGGGTGGCCGAGCTCGTGGTGGTGCCCTTAGCAGGCGGGCCGGCGGGCGACGAC"
                             "TATGAGCAGTATCCCCCCTATCTAGCGCCTGTAGACAAAGGTGTCTGCGCGAGCGAGCGACGGAAGCCGCTGCGCCCGAAAGCGGAAGTACACCCACCG"
                             "CGGGTGGTAGTCAAGCGTCAGGCCTAGCCGAAGCTCCAACTCCATGGTGTTGCCTACCAGCCCGTCCGGCGGGAACCTACGCTAGTGGTGCATCACTGG"
                             "CGGCTTGAACGGCATCTAGTCGCCGCGCCCGGCGCCGAGCGGCCAGGAAGCTGGGGCCCAGCACTCGTCGTCGCAGCCGTCTTACGGGCCGGGCGGGTG"
                             "GCGTCACCTCAGGGTGTACACGCGGAAAAGGCCTAGGGCGACTATCTCA";
    ASSERT_EQ(ref.size(), 309u);
    ASSERT_EQ(read.size(), 544u);
    int const q_end_free = static_cast<int>(read.size());
    for (size_t x_drop : { size_t{0}, size_t{1000} }) {
        WFA2Wrapper2 classic(4, 6, 2, x_drop, WfaTraceback::Classic), single(4, 6, 2, x_drop, WfaTraceback::Singletrack);
        classic.Alignment(read, ref, 0, q_end_free, 0, 0, 1000);
        size_t const fallbacks = WFA2Wrapper2::SingletrackFallbacks();
        single.Alignment(read, ref, 0, q_end_free, 0, 0, 1000);
        ASSERT_TRUE(classic.Success());
        ASSERT_TRUE(single.Success());
        EXPECT_EQ(-classic.GetAlignmentScore(), 740);  // wf-adaptive's; 586 without it
        EXPECT_EQ(single.GetAlignmentScore(), classic.GetAlignmentScore());
        EXPECT_EQ(AlignmentScore(single.Cigar(), read, ref, 0, 0, q_end_free), 740);
        EXPECT_EQ(WFA2Wrapper2::SingletrackFallbacks(), fallbacks) << "the patched backtrace finds the path itself";
    }
}

// Flanks like that one at random: a reference piece, the read its copy at up to 25% divergence with unrelated sequence
// after it, the read's end free and the reference's start not; and unrelated reads. Singletrack's outcomes and scores
// are the classic ones, each of its alignments one of that score (fallbacks to the classic traceback included).
TEST(WfaTraceback, SingletrackWithTheReadsEndFree) {
    std::mt19937 rng(31);
    std::uniform_real_distribution<double> rate(0.0, 0.25);
    std::vector<Case> cases;
    for (size_t i = 0; i < 400; i++) {
        std::string const ref = Random(150 + rng() % 450, rng);
        std::string const read = i % 4 == 3 ? Random(200 + rng() % 400, rng) : Mutate(ref, rate(rng), rng) + Random(rng() % 300, rng);
        cases.push_back({ read, ref, static_cast<int>(read.size()), 0, static_cast<int>(rng() % 2 ? 0 : ref.size() / 4) });
    }
    size_t const fallbacks = WFA2Wrapper2::SingletrackFallbacks();
    auto const [aligned, other] = CompareTracebacks(cases, 0, 2000);
    EXPECT_GT(aligned, cases.size() / 2);
    RecordProperty("read_end_free_other_operations", std::to_string(other) + " of " + std::to_string(aligned));
    RecordProperty("read_end_free_fallbacks", std::to_string(WFA2Wrapper2::SingletrackFallbacks() - fallbacks));
}

// End to end (the pieces between a long read's links): the same scores.
TEST(WfaTraceback, SingletrackEndToEnd) {
    std::mt19937 rng(29);
    WFA2Wrapper2 classic(4, 6, 2, 0, WfaTraceback::Classic), single(4, 6, 2, 0, WfaTraceback::Singletrack);
    for (size_t i = 0; i < 500; i++) {
        std::string const ref = Random(20 + rng() % 400, rng);
        std::string const read = Mutate(ref, 0.12, rng);
        classic.Reset();
        single.Reset();
        classic.EndToEnd(read, ref, 400);
        single.EndToEnd(read, ref, 400);
        ASSERT_EQ(classic.Success(), single.Success());
        if (!classic.Success()) continue;
        EXPECT_EQ(classic.GetAlignmentScore(), single.GetAlignmentScore());
        EXPECT_EQ(AlignmentScore(single.Cigar(), read, ref, 0, 0, 0), -single.GetAlignmentScore());
    }
}

// Alignment threads work on copies, and a long read's pieces on an aligner made from the window's: both keep the traceback.
TEST(WfaTraceback, ACopyKeepsItsTraceback) {
    WFA2Wrapper2 single(4, 6, 2, 1000, WfaTraceback::Singletrack);
    WFA2Wrapper2 copy(single);
    EXPECT_EQ(copy.Traceback(), WfaTraceback::Singletrack);
    EXPECT_EQ(WFA2Wrapper2(4, 6, 2, 0).Traceback(), WfaTraceback::Classic);
}

TEST(WfaTraceback, NamesAndResolution) {
    EXPECT_EQ(ParseWfaTraceback("classic"), WfaTraceback::Classic);
    EXPECT_EQ(ParseWfaTraceback("singletrack"), WfaTraceback::Singletrack);
    EXPECT_FALSE(ParseWfaTraceback("auto"));
    EXPECT_FALSE(ParseWfaTraceback("Singletrack"));
    EXPECT_EQ(ParseWfaTraceback(WfaTracebackName(WfaTraceback::Singletrack)), WfaTraceback::Singletrack);
    // given beats the database's; a database without build_settings.tsv was built classic
    EXPECT_EQ(ResolveWfaTraceback(WfaTraceback::Classic, WfaTraceback::Singletrack), WfaTraceback::Classic);
    EXPECT_EQ(ResolveWfaTraceback(std::nullopt, WfaTraceback::Singletrack), WfaTraceback::Singletrack);
    EXPECT_EQ(ResolveWfaTraceback(std::nullopt, std::nullopt), WfaTraceback::Classic);
}

TEST(BuildSettings, WrittenSettingsReadBack) {
    for (auto const traceback : { WfaTraceback::Classic, WfaTraceback::Singletrack }) {
        build_settings::Settings written;
        written.wfa_traceback = traceback;
        std::stringstream file;
        build_settings::Write(file, written);
        build_settings::Settings read;
        EXPECT_EQ(build_settings::Read(file, read), "");
        EXPECT_EQ(read.wfa_traceback, traceback);
    }
}

TEST(BuildSettings, UnknownSettingsSkippedBadValuesRefused) {
    build_settings::Settings settings;
    std::istringstream later("# comment\n\nlater_setting\t7\r\nwfa_traceback\tsingletrack\r\n");
    EXPECT_EQ(build_settings::Read(later, settings), "");
    EXPECT_EQ(settings.wfa_traceback, WfaTraceback::Singletrack);

    build_settings::Settings none;
    std::istringstream empty("# nothing set\n");
    EXPECT_EQ(build_settings::Read(empty, none), "");
    EXPECT_FALSE(none.wfa_traceback);

    std::istringstream bad("wfa_traceback\tauto\n");
    EXPECT_NE(build_settings::Read(bad, settings).find("line 1"), std::string::npos);
    std::istringstream untabbed("wfa_traceback singletrack\n");
    EXPECT_NE(build_settings::Read(untabbed, settings).find("no tab"), std::string::npos);
}

// The fragment count (--fragment_count, 2026-10-12): written and read back beside the traceback; a file without it
// says nothing (mapq); a value other than mapq or soft is refused.
TEST(BuildSettings, FragmentCountWrittenAndReadBack) {
    EXPECT_EQ(ParseFragmentCount("soft"), FragmentCount::Soft);
    EXPECT_EQ(ParseFragmentCount(FragmentCountName(FragmentCount::Mapq)), FragmentCount::Mapq);
    EXPECT_FALSE(ParseFragmentCount("auto"));
    for (auto const count : { FragmentCount::Mapq, FragmentCount::Soft }) {
        build_settings::Settings written;
        written.wfa_traceback = WfaTraceback::Singletrack;
        written.fragment_count = count;
        std::stringstream file;
        build_settings::Write(file, written);
        build_settings::Settings read;
        EXPECT_EQ(build_settings::Read(file, read), "");
        EXPECT_EQ(read.fragment_count, count);
        EXPECT_EQ(read.wfa_traceback, WfaTraceback::Singletrack);
    }
    build_settings::Settings older;
    std::istringstream before("wfa_traceback\tsingletrack\n");
    EXPECT_EQ(build_settings::Read(before, older), "");
    EXPECT_FALSE(older.fragment_count);
    std::istringstream bad("fragment_count\thard\n");
    EXPECT_NE(build_settings::Read(bad, older).find("mapq or soft"), std::string::npos);
}

// The sparse-first gates (--sparse_first, 2026-10-12): auto, one C for every read type, or read types and theirs; a C
// is 0 (off) or a power of two up to 32. Written per read type with their score and hits, read back; a later
// version's read type is skipped, a bad C refused; a file without them says nothing (off).
TEST(BuildSettings, SparseFirstGatesParsedWrittenAndReadBack) {
    auto const index = [](ReadType type) { return static_cast<size_t>(type); };
    auto const automatic = ParseSparseFirstGates("auto");
    ASSERT_TRUE(automatic);
    EXPECT_EQ(*automatic, SparseFirstGates());
    auto const all = ParseSparseFirstGates("4");
    ASSERT_TRUE(all);
    for (auto const& gate : *all) EXPECT_EQ(gate, size_t{4});
    auto const some = ParseSparseFirstGates("pe:1,pb:8,ont:0");
    ASSERT_TRUE(some);
    EXPECT_EQ((*some)[index(ReadType::Paired)], size_t{1});
    EXPECT_FALSE((*some)[index(ReadType::Single)]);
    EXPECT_EQ((*some)[index(ReadType::PacBio)], size_t{8});
    EXPECT_EQ((*some)[index(ReadType::ONT)], size_t{0});
    EXPECT_EQ(SparseFirstGatesText(*some), "pe 1, pb 8, ont off");
    EXPECT_EQ(SparseFirstGatesText(SparseFirstGates()), "none");
    for (char const* bad : { "3", "64", "-1", "", "pe", "pe:", "pe:3", "xx:2", "pe:2,", "pe:2;se:4", "auto,pe:2" }) {
        EXPECT_FALSE(ParseSparseFirstGates(bad)) << bad;
    }

    build_settings::Settings written;
    written.wfa_traceback = WfaTraceback::Singletrack;
    written.sparse_first = *some;
    written.sparse_first_score = 13;
    written.sparse_first_hits = 2;
    std::stringstream file;
    build_settings::Write(file, written);
    EXPECT_NE(file.str().find("sparse_first_pe\t1\nsparse_first_pb\t8\nsparse_first_ont\t0\nsparse_first_score\t13\n"
                              "sparse_first_hits\t2\n"), std::string::npos) << file.str();
    build_settings::Settings read;
    EXPECT_EQ(build_settings::Read(file, read), "");
    EXPECT_EQ(read.sparse_first, *some);
    EXPECT_EQ(read.sparse_first_score, size_t{13});
    EXPECT_EQ(read.sparse_first_hits, size_t{2});

    build_settings::Settings later;
    std::istringstream unknown("sparse_first_xr\t4\nsparse_first_se\t16\n");
    EXPECT_EQ(build_settings::Read(unknown, later), "");
    EXPECT_EQ(later.sparse_first[index(ReadType::Single)], size_t{16});
    build_settings::Settings none;
    std::istringstream older("wfa_traceback\tsingletrack\n");
    EXPECT_EQ(build_settings::Read(older, none), "");
    EXPECT_EQ(none.sparse_first, SparseFirstGates());
    EXPECT_FALSE(none.sparse_first_score);
    for (char const* bad : { "sparse_first_pe\t3\n", "sparse_first_score\t17\n", "sparse_first_hits\t0\n" }) {
        std::istringstream is(bad);
        build_settings::Settings settings;
        EXPECT_NE(build_settings::Read(is, settings).find("line 1"), std::string::npos) << bad;
    }
}

// The gates' picks (--sparse_first_pick, 2026-10-13): auto, hash or stride for every read type, or read types and
// theirs. Written per read type after the gates, read back; the gates' lines are not taken for picks nor the picks'
// for gates; a bad pick refused; a file without them says nothing (hash).
TEST(BuildSettings, SparseFirstPicksParsedWrittenAndReadBack) {
    auto const index = [](ReadType type) { return static_cast<size_t>(type); };
    auto const automatic = ParseSparseFirstPicks("auto");
    ASSERT_TRUE(automatic);
    EXPECT_EQ(*automatic, SparseFirstPicks());
    auto const all = ParseSparseFirstPicks("stride");
    ASSERT_TRUE(all);
    for (auto const& pick : *all) EXPECT_EQ(pick, SparseFirstPick::Stride);
    auto const some = ParseSparseFirstPicks("pe:stride,pb:hash");
    ASSERT_TRUE(some);
    EXPECT_EQ((*some)[index(ReadType::Paired)], SparseFirstPick::Stride);
    EXPECT_FALSE((*some)[index(ReadType::Single)]);
    EXPECT_EQ((*some)[index(ReadType::PacBio)], SparseFirstPick::Hash);
    EXPECT_EQ(SparseFirstPicksText(*some), "pe stride, pb hash");
    EXPECT_EQ(SparseFirstPicksText(SparseFirstPicks()), "none");
    for (char const* bad : { "strided", "8", "", "pe", "pe:", "pe:8", "xx:stride", "pe:stride,", "auto,pe:hash" }) {
        EXPECT_FALSE(ParseSparseFirstPicks(bad)) << bad;
    }

    build_settings::Settings written;
    written.sparse_first = *ParseSparseFirstGates("pe:4,pb:8");
    written.sparse_first_score = 14;
    written.sparse_first_hits = 1;
    written.sparse_first_pick = *some;
    std::stringstream file;
    build_settings::Write(file, written);
    EXPECT_NE(file.str().find("sparse_first_pick_pe\tstride\nsparse_first_pick_pb\thash\n"), std::string::npos) << file.str();
    build_settings::Settings read;
    EXPECT_EQ(build_settings::Read(file, read), "");
    EXPECT_EQ(read.sparse_first_pick, *some);
    EXPECT_EQ(read.sparse_first, written.sparse_first);

    build_settings::Settings later;
    std::istringstream unknown("sparse_first_pick_xr\tstride\nsparse_first_pick_ont\tstride\n");
    EXPECT_EQ(build_settings::Read(unknown, later), "");
    EXPECT_EQ(later.sparse_first_pick[index(ReadType::ONT)], SparseFirstPick::Stride);
    EXPECT_EQ(later.sparse_first, SparseFirstGates());
    build_settings::Settings none;
    std::istringstream older("sparse_first_pe\t2\n");
    EXPECT_EQ(build_settings::Read(older, none), "");
    EXPECT_EQ(none.sparse_first_pick, SparseFirstPicks());
    std::istringstream bad("sparse_first_pick_pe\t8\n");
    build_settings::Settings settings;
    EXPECT_NE(build_settings::Read(bad, settings).find("hash or stride"), std::string::npos);
}

// The pick a database's runs gate with by default (DatabasePick): as recorded; hash for a gate recorded without a pick
// (databases built on 2026-10-12/13, before the pick could be chosen); stride for a read type with neither (an
// ungated database's runs gated by --sparse_first: protal's default since 2026-10-13).
TEST(BuildSettings, DatabasePickDefaultsToStrideButHashForOlderGates) {
    auto const index = [](ReadType type) { return static_cast<size_t>(type); };
    build_settings::Settings settings;
    std::istringstream file("sparse_first_pe\t2\nsparse_first_se\t4\nsparse_first_pick_se\tstride\nsparse_first_pick_pb\thash\n");
    ASSERT_EQ(build_settings::Read(file, settings), "");
    EXPECT_EQ(build_settings::DatabasePick(settings, index(ReadType::Paired)), SparseFirstPick::Hash);
    EXPECT_EQ(build_settings::DatabasePick(settings, index(ReadType::Single)), SparseFirstPick::Stride);
    EXPECT_EQ(build_settings::DatabasePick(settings, index(ReadType::PacBio)), SparseFirstPick::Hash);
    EXPECT_EQ(build_settings::DatabasePick(settings, index(ReadType::ONT)), SparseFirstPick::Stride);
    EXPECT_EQ(build_settings::DatabasePick(build_settings::Settings(), index(ReadType::Paired)), SparseFirstPick::Stride);
}

// A sample is profiled with the model of the gate nearest its own (NearestGate): on the scale of log2 C, off one step
// below C 1, the lower of two as near.
TEST(BuildSettings, NearestGateOnTheLogScaleLowerOnTies) {
    std::vector<size_t> const ladder = { 1, 4, 8 };
    EXPECT_EQ(NearestGate(4, ladder), 4u);
    EXPECT_EQ(NearestGate(0, ladder), 1u);
    EXPECT_EQ(NearestGate(1, ladder), 1u);
    EXPECT_EQ(NearestGate(2, ladder), 1u);  // 1 and 4 one step away each: the lower
    EXPECT_EQ(NearestGate(16, ladder), 8u);
    EXPECT_EQ(NearestGate(32, ladder), 8u);
    EXPECT_EQ(NearestGate(2, std::vector<size_t>{ 8, 16, 32 }), 8u);
    EXPECT_EQ(NearestGate(0, std::vector<size_t>{ 0, 4 }), 0u);
    EXPECT_EQ(NearestGate(1, std::vector<size_t>{ 0, 4 }), 0u);  // off one step below 1, 4 two above
    EXPECT_EQ(NearestGate(8, std::vector<size_t>{ 4 }), 4u);
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

namespace {
    // WFA2 as the wrapper set it up until 2026-10-11: one aligner with the penalties as given (4/6/2) and X-drop added to
    // wf-adaptive, Singletrack's failed backtraces aligned again classic, and with X-drop an alignment that does not
    // cover the read counted as failed.
    struct Unscaled {
        wfa::WFAlignerGapAffine aligner, classic;
        size_t x_drop;
        wfa::WFAligner* last = nullptr;
        bool success = false;

        Unscaled(size_t x_drop, WfaTraceback traceback) :
                aligner(4, 6, 2, wfa::WFAligner::Alignment,
                        traceback == WfaTraceback::Singletrack ? wfa::WFAligner::MemorySingletrack : wfa::WFAligner::MemoryHigh),
                classic(4, 6, 2, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryHigh),
                x_drop(x_drop) {
            if (x_drop) {
                aligner.setHeuristicXDrop(static_cast<int>(x_drop), 1);
                classic.setHeuristicXDrop(static_cast<int>(x_drop), 1);
            }
        }

        void Align(Case const& c, int max_score, bool end_to_end) {
            auto run = [&](wfa::WFAlignerGapAffine& a) {
                a.setMaxAlignmentSteps(max_score);
                last = &a;
                return end_to_end ? a.alignEnd2End(c.ref, c.read)
                                  : a.alignEndsFree(c.ref, c.r_begin_free, c.r_end_free, c.read, 0, c.q_end_free);
            };
            auto status = run(aligner);
            if (static_cast<int>(status) == WF_STATUS_UNATTAINABLE) status = run(classic);
            success = status == wfa::WFAligner::StatusAlgCompleted && (x_drop == 0 || ReadBases(last->getAlignment()) == c.read.size());
        }
    };

    // The wrapper and Unscaled on `cases` at budgets around each alignment's score (score - 1 fails, score + 1 succeeds:
    // the wrapper's halved step limit must cut where WFA2 cut with the penalties as given), at 61 (short reads' budget),
    // 0 and none: the same outcome, score and operations. end_to_end: the whole read against the whole reference, the
    // free ends ignored. -> the successful alignments compared.
    size_t CompareUnscaled(std::vector<Case> const& cases, size_t x_drop, WfaTraceback traceback, bool end_to_end) {
        WFA2Wrapper2 wrapper(4, 6, 2, x_drop, traceback);
        Unscaled unscaled(x_drop, traceback);
        size_t compared = 0;
        for (auto const& c : cases) {
            unscaled.Align(c, INT32_MAX, end_to_end);
            int const score = unscaled.success ? -unscaled.last->getAlignmentScore() : 61;
            for (int const budget : { score - 1, score, score + 1, score + 2, 61, 0, 1, INT32_MAX }) {
                unscaled.Align(c, budget, end_to_end);
                wrapper.Reset();
                if (end_to_end) wrapper.EndToEnd(c.read, c.ref, budget);
                else wrapper.Alignment(c.read, c.ref, 0, c.q_end_free, c.r_begin_free, c.r_end_free, budget);
                EXPECT_EQ(wrapper.Success(), unscaled.success) << "x_drop " << x_drop << ", budget " << budget << " (score "
                                                               << score << "), read of " << c.read.size();
                if (!wrapper.Success() || !unscaled.success) continue;
                compared++;
                EXPECT_EQ(wrapper.GetAlignmentScore(), unscaled.last->getAlignmentScore()) << "budget " << budget;
                EXPECT_EQ(wrapper.Cigar(), unscaled.last->getAlignment()) << "budget " << budget;
            }
        }
        return compared;
    }
}

// The wrapper aligns with the penalties divided by 2 (2/3/1) and without X-drop where X-drop cannot prune (WFA2Wrapper2's
// class comment): WFA2's alignments with 4/6/2 and X-drop, the same outcome, score and operations, at budgets on either
// side of each score. Short reads take the halved aligner at X-drop 0 and 1000, the X-drop one at 50 and 20; long reads
// the halved one at 0 only.
TEST(WFA2Wrapper, HalvedPenaltiesGiveTheSameAlignments) {
    std::mt19937 rng(17);
    auto const cases = Cases(400, rng);
    auto const long_cases = LongCases(24, rng);
    for (auto traceback : { WfaTraceback::Classic, WfaTraceback::Singletrack }) {
        for (size_t x_drop : { size_t{0}, size_t{1000}, size_t{50}, size_t{20} }) {
            EXPECT_GT(CompareUnscaled(cases, x_drop, traceback, false), cases.size()) << "x_drop " << x_drop;
            EXPECT_GT(CompareUnscaled(cases, x_drop, traceback, true), cases.size()) << "x_drop " << x_drop;
            EXPECT_GT(CompareUnscaled(long_cases, x_drop, traceback, false), 0u) << "x_drop " << x_drop;
        }
    }
}

// X-drop just above the wrapper's bound (read, reference and budget adding up to less than twice X-drop) prunes nothing:
// the halved aligner without X-drop gives what WFA2 with X-drop gives. Just below it the wrapper aligns with X-drop.
TEST(WFA2Wrapper, XDropAboveTheBoundPrunesNothing) {
    std::mt19937 rng(23);
    auto cases = Cases(300, rng);
    for (auto const& c : LongCases(12, rng)) cases.push_back(c);
    size_t compared = 0;
    for (auto const& c : cases) {
        for (int const budget : { 61, static_cast<int>(c.read.size()) }) {
            int64_t const sum = static_cast<int64_t>(c.read.size() + c.ref.size()) + budget;
            size_t const x_drop = static_cast<size_t>(sum / 2 + 1);
            WFA2Wrapper2 wrapper(4, 6, 2, x_drop);
            ASSERT_FALSE(wrapper.XDropMayPrune(c.read.size(), c.ref.size(), budget));
            ASSERT_TRUE(WFA2Wrapper2(4, 6, 2, x_drop - 1).XDropMayPrune(c.read.size(), c.ref.size(), budget));
            Unscaled unscaled(x_drop, WfaTraceback::Classic);
            unscaled.Align(c, budget, false);
            wrapper.Alignment(c.read, c.ref, 0, c.q_end_free, c.r_begin_free, c.r_end_free, budget);
            ASSERT_EQ(wrapper.Success(), unscaled.success) << "x_drop " << x_drop << ", read of " << c.read.size();
            if (!wrapper.Success()) continue;
            compared++;
            EXPECT_EQ(wrapper.GetAlignmentScore(), unscaled.last->getAlignmentScore());
            EXPECT_EQ(wrapper.Cigar(), unscaled.last->getAlignment());
        }
    }
    EXPECT_GT(compared, cases.size() / 2);
    // protal's short reads at the default X-drop take the halved aligner; long reads (aligned without X-drop) would not.
    EXPECT_FALSE(WFA2Wrapper2(4, 6, 2, 1000).XDropMayPrune(150, 400, 61));
    EXPECT_TRUE(WFA2Wrapper2(4, 6, 2, 1000).XDropMayPrune(2000, 2400, 1201));
}
