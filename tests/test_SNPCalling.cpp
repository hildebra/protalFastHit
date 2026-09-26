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

TEST(StrandFilter, ReferenceAlleleIsExempt) {
    Variant ref(10, 'A', 'A');
    ref.SetObservations(20);  // inferred from coverage: no strand information
    EXPECT_TRUE(ref.IsReference());
    EXPECT_FALSE(ref.HasFwdAndRev());
    EXPECT_TRUE(ref.PassesStrandFilter());

    Variant alt(10, 'G', 'A');
    alt.AddObservation(30, true);
    alt.AddObservation(30, true);
    EXPECT_FALSE(alt.PassesStrandFilter());
    alt.AddObservation(30, false);
    EXPECT_TRUE(alt.PassesStrandFilter());
}

TEST(StrandFilter, ReferenceConsensusPassesMsaFilter) {
    // One forward-only error read at a well-covered position must not turn the REF call into N.
    std::string reference(50, 'A');
    VariantHandler handler(reference);
    Variant alt(10, 'G', 'A');
    alt.AddObservation(30, true);
    VariantBin bin{alt};
    handler.PostProcessSNPBin(bin, 20, 2, 2, 0.0, 0, 0, true);

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
    handler.PostProcessSNPBin(bin, 1, 5, 5, 0.0, 0, 0, false);

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
