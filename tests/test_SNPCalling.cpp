// Unit tests for SNP calling: quality parsing, the reference allele's strand handling,
// variant ownership and matching, IUPAC codes and SAM strand flags.
#include <gtest/gtest.h>
#include <string>
#include <vector>
#include "Profiling/Strain.h"

using namespace protal;

namespace {
    SamEntry MakeSam(std::string seq, std::string qual, std::string cigar, POS_t pos, FLAG_t flag) {
        SamEntry sam;
        sam.m_qname = "read";
        sam.m_flag = flag;
        sam.m_rname = "1_1";
        sam.m_pos = pos;
        sam.m_mapq = 60;
        sam.m_seq = std::move(seq);
        sam.m_qual = std::move(qual);
        sam.m_cigar = std::move(cigar);
        return sam;
    }
}

TEST(PhredScore, DecodesPhred33) {
    EXPECT_EQ(VariantHandler::PhredScore('!'), 0);
    EXPECT_EQ(VariantHandler::PhredScore('#'), 2);
    EXPECT_EQ(VariantHandler::PhredScore('+'), 10);
    EXPECT_EQ(VariantHandler::PhredScore('I'), 40);
}

TEST(PhredScore, ClampsBelowOffsetInsteadOfWrapping) {
    EXPECT_EQ(VariantHandler::PhredScore(' '), 0);
    EXPECT_EQ(VariantHandler::PhredScore('\0'), 0);
}

TEST(ExtractVariants, LowQualitySnpKeepsLowQuality) {
    // The read carries G instead of the reference T at position 3 (0-based), with quality '#' (Q2).
    std::string reference = "ACGTACGTAC";
    VariantHandler handler(reference);
    auto sam = MakeSam("ACGGACGTAC", "III#IIIIII", "3M1X6M", 1, 0);

    ASSERT_TRUE(handler.AddVariantsFromSam(sam));
    ASSERT_TRUE(handler.HasVariantBin(3));
    auto& bin = handler.GetVariantBin(3);
    ASSERT_EQ(bin.size(), 1u);
    EXPECT_EQ(bin.front().GetVariant(), 'G');
    EXPECT_EQ(bin.front().QualitySum(), 2u);
}

TEST(StrandFilter, EveryAlleleIsJudgedByTheStrandsOfTheSite) {
    // A site of 10 forward and 10 reverse reads.
    auto site = [](uint32_t alt_forward, uint32_t alt_reverse) {
        Variant ref(10, 'A', 'A');
        ref.SetObservations(10 - alt_forward, 10 - alt_reverse);
        Variant alt(10, 'G', 'A');
        for (uint32_t i = 0; i < alt_forward; i++) alt.AddObservation(30, true);
        for (uint32_t i = 0; i < alt_reverse; i++) alt.AddObservation(30, false);
        return VariantBin{ ref, alt };
    };
    auto bin = site(2, 0);
    EXPECT_TRUE(PassesStrandFilter(bin[1], bin)) << "2 of 10 forward reads: chance";
    bin = site(6, 0);
    EXPECT_FALSE(PassesStrandFilter(bin[1], bin)) << "6 of the forward reads and none of the reverse: bias";
    EXPECT_TRUE(PassesStrandFilter(bin[0], bin)) << "the reference is on both strands";
    bin = site(6, 5);
    EXPECT_TRUE(PassesStrandFilter(bin[1], bin));

    // All reads on one strand (common at low depth): no evidence of bias, for any allele.
    Variant ref(10, 'A', 'A');
    ref.SetObservations(0, 3);
    Variant alt(10, 'G', 'A');
    alt.AddObservation(30, false);
    alt.AddObservation(30, false);
    VariantBin one_strand{ ref, alt };
    EXPECT_TRUE(PassesStrandFilter(one_strand[0], one_strand));
    EXPECT_TRUE(PassesStrandFilter(one_strand[1], one_strand));

    // The reference allele is tested too: 8 forward reads show it, 8 reverse reads the SNP.
    Variant split_ref(10, 'A', 'A');
    split_ref.SetObservations(8, 0);
    Variant split_alt(10, 'G', 'A');
    for (int i = 0; i < 8; i++) split_alt.AddObservation(30, false);
    VariantBin split{ split_ref, split_alt };
    EXPECT_FALSE(PassesStrandFilter(split[0], split));
    EXPECT_FALSE(PassesStrandFilter(split[1], split));
}

TEST(StrandFilter, ReferenceConsensusPassesMsaFilter) {
    // One forward-only error read at a well-covered position must not turn the REF call into N.
    std::string reference(50, 'A');
    VariantHandler handler(reference);
    Variant alt(10, 'G', 'A');
    alt.AddObservation(30, true);
    VariantBin bin{alt};
    handler.PostProcessSNPBin(bin, 10, 10, 2, 2, 0.0, 0, 0, true);

    auto ref_it = std::find_if(bin.begin(), bin.end(), [](Variant const& v) { return v.IsReference(); });
    ASSERT_NE(ref_it, bin.end());
    EXPECT_EQ(ref_it->Observations(), 19u);
    EXPECT_TRUE(ref_it->GetValid());
    EXPECT_TRUE(VariantPass(*ref_it, bin, 50, 3, 0.0, 20, /*require_strand=*/true, 0));
}

TEST(PostProcessSNPBin, FiltersBinsWithoutReferenceReads) {
    // Every read carries a variant (coverage == observations): no REF allele is added, but the
    // variants are still filtered rather than left valid by default.
    std::string reference(50, 'A');
    VariantHandler handler(reference);
    Variant alt(10, 'G', 'A');
    alt.AddObservation(30, true);
    VariantBin bin{alt};
    handler.PostProcessSNPBin(bin, 1, 0, 5, 5, 0.0, 0, 0, false);

    ASSERT_EQ(bin.size(), 1u);
    EXPECT_FALSE(bin.front().GetValid());
}

TEST(VariantMatch, DeletionsMatchDeletions) {
    std::string deleted = "CG";
    Variant a(VariantType::DEL, 5, 'A', deleted);
    Variant b(VariantType::DEL, 5, 'A', deleted);
    std::string other = "CGT";
    Variant c(VariantType::DEL, 5, 'A', other);
    EXPECT_TRUE(a.Match(b));
    EXPECT_FALSE(a.Match(c));
}

TEST(VariantOwnership, CopiesOwnTheirIndelSequence) {
    std::string inserted = "TTG";
    std::vector<Variant> bin;
    {
        Variant ins(VariantType::INS, 7, 'A', inserted);
        bin.push_back(ins);
        Variant assigned(3, 'C', 'A');
        assigned = ins;
        bin.push_back(assigned);
    }
    for (int i = 0; i < 100; i++) bin.push_back(bin.front());  // force reallocation and copies

    for (auto const& v : bin) {
        EXPECT_TRUE(v.IsINS());
        EXPECT_EQ(v.GetStructural(), "TTG");
        EXPECT_EQ(v.GetStructuralSize(), 3u);
    }
}

TEST(IUPAC, TwoAndThreeAlleleCodes) {
    std::vector<char> ag{'A', 'G'};
    std::vector<char> ct{'C', 'T'};
    std::vector<char> acg{'A', 'C', 'G'};
    std::vector<char> single{'T'};
    EXPECT_EQ(IUPACCode(ag), 'R');
    EXPECT_EQ(IUPACCode(ct), 'Y');
    EXPECT_EQ(IUPACCode(acg), 'V');
    EXPECT_EQ(IUPACCode(single), 'T');
}

TEST(SamFlags, Read2StrandComesFromItsOwnFlagBit) {
    const FLAG_t paired = 0x1, reverse = 0x10, mate_reverse = 0x20, read1 = 0x40, read2 = 0x80;

    auto r2_reverse = MakeSam("ACGT", "IIII", "4M", 1, paired | read2 | reverse);
    auto r2_forward_mate_reverse = MakeSam("ACGT", "IIII", "4M", 1, paired | read2 | mate_reverse);
    auto r1_reverse = MakeSam("ACGT", "IIII", "4M", 1, paired | read1 | reverse);

    EXPECT_TRUE(r2_reverse.IsReversed());
    EXPECT_FALSE(r2_forward_mate_reverse.IsReversed());
    EXPECT_TRUE(r1_reverse.IsReversed());
}

TEST(IsAlignmentValid, AcceptsMatchesAndRejectsMismatchesAndOverruns) {
    std::string gene = "ACGTACGTACGTACGTACGT";
    AlignmentInfo info;
    info.gene_alignment_start = 4;
    info.compressed_cigar = "8M";
    EXPECT_TRUE(IsAlignmentValid(info, "ACGTACGT", gene, 0, true));
    EXPECT_TRUE(IsAlignmentValid(info, "ACGNACGT", gene, 0, true));   // N matches anything
    EXPECT_FALSE(IsAlignmentValid(info, "ACGTTCGT", gene, 0, true));  // mismatch inside M

    // CIGARs running past the read or the gene are invalid, and must not be read out of bounds.
    info.compressed_cigar = "12M";
    EXPECT_FALSE(IsAlignmentValid(info, "ACGTACGT", gene, 0, true));
    info.gene_alignment_start = 16;
    info.compressed_cigar = "8M";
    EXPECT_FALSE(IsAlignmentValid(info, "ACGTACGT", gene, 0, true));

    // Soft clips and insertions consume read bases only.
    info.gene_alignment_start = 4;
    info.compressed_cigar = "2S4M1I3M";
    EXPECT_TRUE(IsAlignmentValid(info, "TTACGTGACG", gene, 0, true));
}
