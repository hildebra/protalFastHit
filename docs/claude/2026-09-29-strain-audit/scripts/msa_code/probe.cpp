// Audit probes for protal's MSA construction (Strain.h MSA()). Prints rows; EXPECTs document what
// a correct MSA would hold, so failures are the findings.
#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include <unistd.h>
#include "Profiling/Profiler.h"
#include "Profiling/Strain.h"

using namespace protal;

namespace {
    constexpr char kReference[] = "ACGTTGCAAGGCTTACCGATGACTGAAACCGGTTTACGATCGGTAGCATG";  // 50 bp

    SamEntry MakeSam(std::string seq, std::string cigar, POS_t pos, FLAG_t flag = 0, std::string qual = "") {
        SamEntry sam;
        sam.m_qname = "read";
        sam.m_flag = flag;
        sam.m_rname = "1_1";
        sam.m_pos = pos;
        sam.m_mapq = 60;
        sam.m_qual = qual.empty() ? std::string(seq.size(), 'I') : std::move(qual);
        sam.m_seq = std::move(seq);
        sam.m_cigar = std::move(cigar);
        return sam;
    }

    struct TinyReference {
        std::filesystem::path dir;
        std::unique_ptr<GenomeLoader> loader;
        TinyReference() {
            dir = std::filesystem::temp_directory_path() / ("protal_probe_" + std::to_string(::getpid()));
            std::filesystem::create_directories(dir);
            std::string header = ">1_1\n";
            std::string gene = kReference;
            { std::ofstream fna(dir / "reference.fna"); fna << header << gene << '\n'; }
            { std::ofstream map(dir / "reference.map"); map << "1\t1\t" << header.size() << '\t' << header.size() + gene.size() << '\n'; }
            loader = std::make_unique<GenomeLoader>((dir / "reference.fna").string(), (dir / "reference.map").string());
            loader->LoadAllGenomes();
        }
        ~TinyReference() { std::filesystem::remove_all(dir); }
        Gene& gene() { return loader->GetGenome(1).GetGeneOMP(1); }
    };

    char Alt(char c) { return c == 'A' ? 'C' : 'A'; }

    // protal's defaults: snp_min_cov 2, phred sum 90, mean qual 15, strand filter, 3 IUPAC alleles.
    struct Run {
        std::vector<std::string> rows;
        std::string ref_row;
        bool ok;
    };
    Run RunMSA(std::vector<StrainLevelContainer*> strains, std::string const& reference, std::vector<double> afs = {},
               size_t max_alleles = 3, uint32_t min_cov = 2) {
        if (afs.empty()) afs.assign(strains.size(), 0.0);
        MSASequenceItems items;
        for (auto* s : strains) {
            if (!s) { items.emplace_back(OptionalMSASequenceItem{}); continue; }
            items.emplace_back(OptionalMSASequenceItem{ { SharedAlignmentRegion::GetSNPs(s->GetVariantHandler()), s->GetSequenceRangeHandler() } });
        }
        MSAVector msa(items.size());
        MSARow ref_row;
        Run r;
        r.ok = MSA(items, reference, msa, min_cov, 90, afs, true, 15, nullptr, &ref_row, max_alleles);
        for (auto& row : msa) r.rows.emplace_back(row.begin(), row.end());
        r.ref_row = std::string(ref_row.begin(), ref_row.end());
        return r;
    }
    void Print(std::string const& name, Run const& r) {
        std::cout << "== " << name << " ok=" << r.ok << "\n  ref: " << r.ref_row << '\n';
        for (size_t i = 0; i < r.rows.size(); i++) std::cout << "  s" << i << ":  " << r.rows[i] << '\n';
    }
    std::string BinString(StrainLevelContainer const& s, VariantPos pos) {
        auto const& v = s.GetVariantHandler().GetVariants();
        auto it = v.find(pos);
        return it == v.end() ? "(no bin)" : VariantHandler::VariantBinToString(it->second);
    }
}

// Q1/Q3: a position covered only by reads with N there.
TEST(Probe, NOnlyPositionIsWrittenAsReference) {
    TinyReference ref; std::string const reference = ref.gene().Sequence();
    StrainLevelContainer s(ref.gene());
    auto nread = reference; nread[10] = 'N';
    for (int i = 0; i < 3; i++) ASSERT_TRUE(s.AddSam(MakeSam(nread, "50M", 1, i % 2 ? 0x10 : 0), i, true));
    s.PostProcess(2, 2, 0.0, 15, 90, true);
    auto r = RunMSA({ &s }, reference);
    Print("N-only at 10 (3 reads, all N)", r);
    EXPECT_NE(r.rows[0][10], reference[10]) << "no read shows the reference base at 10";

    StrainLevelContainer s2(ref.gene());
    for (int i = 0; i < 2; i++) ASSERT_TRUE(s2.AddSam(MakeSam(nread, "50M", 1, i % 2 ? 0x10 : 0), i, true));
    ASSERT_TRUE(s2.AddSam(MakeSam(reference, "50M", 1, 0), 2, true));
    s2.PostProcess(2, 2, 0.0, 15, 90, true);
    auto r2 = RunMSA({ &s2 }, reference);
    Print("2 N reads + 1 ref read at 10, min_cov 2", r2);
    EXPECT_NE(r2.rows[0][10], reference[10]) << "one read supports the reference base, min_cov is 2";
}

// Q4: the two mates of one fragment overlap: counted as two reads, on both strands.
TEST(Probe, OverlappingMatesOfOneFragmentPassMinCovAndStrand) {
    TinyReference ref; std::string const reference = ref.gene().Sequence();
    StrainLevelContainer s(ref.gene());
    auto read = reference; read[10] = Alt(reference[10]);
    // One fragment (read_id 7): read 1 forward, read 2 reverse, both over the whole gene.
    ASSERT_TRUE(s.AddSam(MakeSam(read, "10M1X39M", 1, 0x1 | 0x40 | 0x20), 7, true));
    ASSERT_TRUE(s.AddSam(MakeSam(read, "10M1X39M", 1, 0x1 | 0x80 | 0x10), 7, true));
    s.PostProcess(2, 2, 0.0, 15, 90, true);
    std::cout << "bin 10: " << BinString(s, 10) << " cov10=" << s.GetSequenceRangeHandler().CalculateCoverageVector2()[10] << '\n';
    auto r = RunMSA({ &s }, reference);
    Print("one fragment, overlapping mates, SNP at 10", r);
    EXPECT_EQ(r.rows[0], std::string(50, '-')) << "a single fragment is depth 1 < snp_min_cov 2";
}

// Q3: a read that AddVariantsFromSam rejects still adds its range (coverage), and its X ops before
// the failing M are kept as variants.
TEST(Probe, RejectedReadKeepsCoverageAndVariants) {
    TinyReference ref; std::string const reference = ref.gene().Sequence();
    StrainLevelContainer s(ref.gene());
    ASSERT_TRUE(s.AddSam(MakeSam(reference, "50M", 1, 0), 0, true));
    auto bad = reference; bad[10] = Alt(reference[10]); bad[30] = Alt(reference[30]);  // 30 mismatches under M
    EXPECT_FALSE(s.AddSam(MakeSam(bad, "10M1X39M", 1, 0), 1, true));
    EXPECT_FALSE(s.AddSam(MakeSam(bad, "10M1X39M", 1, 0x10), 2, true));
    auto cov = s.GetSequenceRangeHandler().CalculateCoverageVector2();
    std::cout << "cov[10]=" << cov[10] << " cov[40]=" << cov[40] << " bin10: " << BinString(s, 10) << '\n';
    s.PostProcess(2, 2, 0.0, 15, 90, true);
    auto r = RunMSA({ &s }, reference);
    Print("1 accepted ref read + 2 rejected reads", r);
    EXPECT_EQ(cov[40], 1u) << "rejected reads must not count as coverage";
    EXPECT_EQ(r.rows[0], std::string(50, '-')) << "only one accepted read: depth 1 < min_cov 2";
}

// Q2: an insertion and a SNP at the same position compete for one consensus.
TEST(Probe, InsertionAndSnpAtOnePositionCompete) {
    TinyReference ref; std::string const reference = ref.gene().Sequence();
    char const alt = Alt(reference[10]);
    std::string read = reference.substr(0, 10) + "CC" + std::string(1, alt) + reference.substr(11);
    // (a) Every read carries the insertion before 10 AND the SNP at 10; the SNP base has Q20.
    StrainLevelContainer a(ref.gene());
    for (int i = 0; i < 4; i++) {
        std::string q(read.size(), 'I'); q[12] = '5';
        ASSERT_TRUE(a.AddSam(MakeSam(read, "10M2I1X39M", 1, i % 2 ? 0x10 : 0, q), i, true));
    }
    a.PostProcess(2, 2, 0.0, 15, 90, true);
    std::cout << "bin10 (a): " << BinString(a, 10) << '\n';
    // (b) Same reads, the insertion bases have Q20 instead.
    StrainLevelContainer b(ref.gene());
    for (int i = 0; i < 4; i++) {
        std::string q(read.size(), 'I'); q[10] = '5'; q[11] = '5';
        ASSERT_TRUE(b.AddSam(MakeSam(read, "10M2I1X39M", 1, i % 2 ? 0x10 : 0, q), i, true));
    }
    b.PostProcess(2, 2, 0.0, 15, 90, true);
    std::cout << "bin10 (b): " << BinString(b, 10) << '\n';
    auto r = RunMSA({ &a, &b }, reference);
    Print("INS CC before 10 + SNP at 10 in every read: (a) SNP Q20, (b) INS Q20", r);
    std::string const truth = reference.substr(0, 10) + "CC" + std::string(1, alt) + reference.substr(11);
    EXPECT_EQ(r.rows[0], truth);
    EXPECT_EQ(r.rows[1], truth);

    // (c) 4 reads with the insertion (reference base at 10), 6 reads with the SNP (no insertion).
    StrainLevelContainer c(ref.gene());
    std::string ins_read = reference.substr(0, 10) + "CC" + reference.substr(10);
    std::string snp_read = reference; snp_read[10] = alt;
    for (int i = 0; i < 4; i++) ASSERT_TRUE(c.AddSam(MakeSam(ins_read, "10M2I40M", 1, i % 2 ? 0x10 : 0), i, true));
    for (int i = 0; i < 6; i++) ASSERT_TRUE(c.AddSam(MakeSam(snp_read, "10M1X39M", 1, i % 2 ? 0x10 : 0), 10 + i, true));
    c.PostProcess(2, 2, 0.0, 15, 90, true);
    std::cout << "bin10 (c): " << BinString(c, 10) << '\n';
    auto rc = RunMSA({ &c }, reference);
    Print("(c) 4 reads INS+ref base at 10, 6 reads SNP at 10 (AF 0, 3 IUPAC alleles)", rc);
    char const iupac = [&]{ std::vector<char> al{ reference[10], alt }; std::sort(al.begin(), al.end()); return IUPACCode(al); }();
    EXPECT_EQ(rc.rows[0][10], iupac) << "4 of 10 reads show the reference base at 10, 6 the SNP";
}

// Q2: the consensus is picked by quality sum, and the reference allele always gets Q40.
TEST(Probe, TopAlleleByQualitySumNotObservations) {
    TinyReference ref; std::string const reference = ref.gene().Sequence();
    char const alt = Alt(reference[10]);
    auto snp_read = reference; snp_read[10] = alt;
    // 7 of 10 reads carry the SNP; all bases Q15 ('0'), as ONT reads.
    StrainLevelContainer s(ref.gene());
    std::string q15(50, '0');
    for (int i = 0; i < 7; i++) ASSERT_TRUE(s.AddSam(MakeSam(snp_read, "10M1X39M", 1, i % 2 ? 0x10 : 0, q15), i, true));
    for (int i = 0; i < 3; i++) ASSERT_TRUE(s.AddSam(MakeSam(reference, "50M", 1, i % 2 ? 0x10 : 0, q15), 10 + i, true));
    s.PostProcess(2, 2, 0.2, 15, 90, true);
    std::cout << "bin10: " << BinString(s, 10) << '\n';
    auto r1 = RunMSA({ &s }, reference, { 0.2 }, 1);
    Print("7/10 reads SNP at 10, Q15, --snp_max_alleles 1", r1);
    EXPECT_EQ(r1.rows[0][10], alt) << "the top allele by observations is the SNP (70%)";

    // Q30 (FASTA reads of pe/se/pb): 5 of 9 reads carry the SNP.
    StrainLevelContainer t(ref.gene());
    std::string q30(50, '?');
    for (int i = 0; i < 5; i++) ASSERT_TRUE(t.AddSam(MakeSam(snp_read, "10M1X39M", 1, i % 2 ? 0x10 : 0, q30), i, true));
    for (int i = 0; i < 4; i++) ASSERT_TRUE(t.AddSam(MakeSam(reference, "50M", 1, i % 2 ? 0x10 : 0, q30), 10 + i, true));
    t.PostProcess(2, 2, 0.0, 15, 90, true);
    auto r2 = RunMSA({ &t }, reference, { 0.0 }, 1);
    Print("5/9 reads SNP at 10, Q30, --snp_max_alleles 1", r2);
    EXPECT_EQ(r2.rows[0][10], alt);

    // Top allele fails, a second passes: 1 reference read (Q40), 2 SNP reads (Q15, both strands).
    StrainLevelContainer u(ref.gene());
    ASSERT_TRUE(u.AddSam(MakeSam(reference, "50M", 1, 0), 0, true));
    for (int i = 0; i < 2; i++) ASSERT_TRUE(u.AddSam(MakeSam(snp_read, "10M1X39M", 1, i % 2 ? 0x10 : 0, q15), 1 + i, true));
    u.PostProcess(2, 2, 0.0, 15, 90, true);
    std::cout << "bin10 (u): " << BinString(u, 10) << '\n';
    auto r3 = RunMSA({ &u }, reference);
    Print("1 ref read Q40 + 2 SNP reads Q15 (fwd+rev) at 10, defaults", r3);
    EXPECT_EQ(r3.rows[0][10], alt) << "the SNP passes every filter, the reference allele (1 read) does not";
}

// Q1: deletion and insertion columns across samples and the reference row.
TEST(Probe, DeletionAndInsertionColumnsLineUp) {
    TinyReference ref; std::string const reference = ref.gene().Sequence();
    StrainLevelContainer a(ref.gene()), b(ref.gene());
    std::string del = reference.substr(0, 20) + reference.substr(23);
    for (int i = 0; i < 4; i++) ASSERT_TRUE(a.AddSam(MakeSam(del, "20M3D27M", 1, i % 2 ? 0x10 : 0), i, true));
    std::string ins = reference.substr(0, 21) + "GG" + reference.substr(21);
    for (int i = 0; i < 4; i++) ASSERT_TRUE(b.AddSam(MakeSam(ins, "21M2I29M", 1, i % 2 ? 0x10 : 0), i, true));
    a.PostProcess(2, 2, 0.0, 15, 90, true);
    b.PostProcess(2, 2, 0.0, 15, 90, true);
    auto r = RunMSA({ &a, &b, nullptr }, reference);
    Print("A: 3bp DEL at 20-22; B: GG inserted before 21; C absent", r);
    EXPECT_EQ(r.ref_row, reference.substr(0, 21) + "--" + reference.substr(21));
    EXPECT_EQ(r.rows[0], reference.substr(0, 20) + "-" + "--" + "--" + reference.substr(23));
    EXPECT_EQ(r.rows[1], reference.substr(0, 21) + "GG" + reference.substr(21));
    EXPECT_EQ(r.rows[2], std::string(52, '-'));
}

// Q3: a deletion right after a soft clip is not called but covers the deleted base.
TEST(Probe, DeletionAfterSoftClipCountsAsReferenceCoverage) {
    TinyReference ref; std::string const reference = ref.gene().Sequence();
    StrainLevelContainer s(ref.gene());
    std::string read = "TT" + reference.substr(1, 48);
    for (int i = 0; i < 3; i++) ASSERT_TRUE(s.AddSam(MakeSam(read, "2S1D48M", 1, i % 2 ? 0x10 : 0), i, true));
    s.PostProcess(2, 2, 0.0, 15, 90, true);
    auto r = RunMSA({ &s }, reference);
    Print("3 reads 2S1D48M at pos 1 (no read base at 0)", r);
    EXPECT_NE(r.rows[0][0], reference[0]) << "no read has a base at position 0";
}

// Q3/Q2: reads with a deletion are counted as reference support inside it; MultiAllelicPositions
// counts a position MSA writes as a gap.
TEST(Probe, InsideADeletion) {
    TinyReference ref; std::string const reference = ref.gene().Sequence();
    StrainLevelContainer s(ref.gene());
    std::string del = reference.substr(0, 20) + reference.substr(23);
    for (int i = 0; i < 6; i++) ASSERT_TRUE(s.AddSam(MakeSam(del, "20M3D27M", 1, i % 2 ? 0x10 : 0), i, true));
    auto snp = reference; snp[21] = Alt(reference[21]);
    for (int i = 0; i < 2; i++) ASSERT_TRUE(s.AddSam(MakeSam(snp, "21M1X28M", 1, i % 2 ? 0x10 : 0), 10 + i, true));
    for (int i = 0; i < 2; i++) ASSERT_TRUE(s.AddSam(MakeSam(reference, "50M", 1, i % 2 ? 0x10 : 0), 20 + i, true));
    s.PostProcess(2, 2, 0.0, 15, 90, true);
    std::cout << "bin21: " << BinString(s, 21) << '\n';
    auto cov = s.GetSequenceRangeHandler().CalculateCoverageVector2();
    size_t multi = MultiAllelicPositions(s.GetVariantHandler().GetVariants(), cov, 2, 90, 0.0, true, 15, 3);
    auto r = RunMSA({ &s }, reference);
    Print("6 reads DEL 20-22, 2 reads SNP at 21, 2 ref reads", r);
    std::cout << "MultiAllelicPositions=" << multi << " IUPAC cells in row=" <<
        std::count_if(r.rows[0].begin(), r.rows[0].end(), [](char c){ return std::string("ACGTN-").find(c) == std::string::npos; }) << '\n';
}
