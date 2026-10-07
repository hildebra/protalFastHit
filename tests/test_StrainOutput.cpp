// Unit tests for the strain output: variant calls from CIGARs (N bases, deletions, indels at the
// ends of an alignment), counters beyond 16 bits, deletion columns in the MSA, and the depth
// estimate behind abundances.
#include <gtest/gtest.h>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <type_traits>
#include <vector>
#include <unistd.h>
#include "Profiling/Profiler.h"
#include "Profiling/Strain.h"
#include "TestReference.h"

using namespace protal;

namespace {
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

    constexpr char kReference[] = "ACGTTGCAAGGCTTACCGATGACTGAAACCGGTTTACGATCGGTAGCATG";  // 50 bp
}

TEST(ExtractVariants, NBasesAreUncalledNotAlleles) {
    std::string reference = kReference;
    VariantHandler handler(reference);
    auto with_n = reference;
    with_n[10] = 'N';
    auto with_snp = reference;
    with_snp[10] = 'T';
    ASSERT_TRUE(handler.AddVariantsFromSam(MakeSam(with_n, "10M1X39M", 1)));
    EXPECT_FALSE(handler.HasVariantBin(10)) << "an N is no allele";
    ASSERT_TRUE(handler.AddVariantsFromSam(MakeSam(with_snp, "10M1X39M", 1, 0x10)));
    ASSERT_TRUE(handler.AddVariantsFromSam(MakeSam(reference, "50M", 1)));

    // Three reads cover position 10: one N and one reference G (forward), one T (reverse).
    handler.PostProcessSNPs(CoverageVec(50, 2), CoverageVec(50, 1), 1, 1, 0.0, 0, 0, false);
    auto& bin = handler.GetVariantBin(10);
    auto ref = std::find_if(bin.begin(), bin.end(), [](Variant const& v) { return v.IsReference(); });
    ASSERT_NE(ref, bin.end());
    EXPECT_EQ(ref->Observations(), 1u) << "the N read is not reference support";
    EXPECT_EQ(std::count_if(bin.begin(), bin.end(), [](Variant const& v) { return v.GetVariant() == 'N'; }), 0);
}

TEST(ExtractVariants, DeletionsGetTheirFlankingQuality) {
    std::string reference = kReference;
    VariantHandler handler(reference);
    auto read = reference.substr(0, 20) + reference.substr(22);
    auto qual = std::string(read.size(), 'I');
    qual[19] = '5';  // Q20 before the deletion, Q40 after it
    ASSERT_TRUE(handler.AddVariantsFromSam(MakeSam(read, "20M2D28M", 1, 0, qual)));
    ASSERT_TRUE(handler.HasVariantBin(20));
    auto& bin = handler.GetVariantBin(20);
    ASSERT_EQ(bin.size(), 1u);
    EXPECT_TRUE(bin.front().IsDEL());
    EXPECT_EQ(bin.front().GetStructural(), reference.substr(20, 2));
    EXPECT_EQ(bin.front().QualitySum(), 20u);
}

TEST(ExtractVariants, IndelsAtAlignmentEndsAreIgnored) {
    std::string reference = kReference;
    VariantHandler handler(reference);
    // Leading insertion, behind a soft clip too, and an insertion after the gene's last base.
    ASSERT_TRUE(handler.AddVariantsFromSam(MakeSam("TTTTT" + reference.substr(0, 15), "5I15M", 1)));
    ASSERT_TRUE(handler.AddVariantsFromSam(MakeSam("AATTT" + reference.substr(0, 15), "2S3I15M", 1)));
    ASSERT_TRUE(handler.AddVariantsFromSam(MakeSam(reference.substr(35) + "GGGGG", "15M5I", 36)));
    ASSERT_TRUE(handler.AddVariantsFromSam(MakeSam(reference.substr(0, 10), "10M2D", 1)));
    EXPECT_TRUE(handler.GetVariants().empty());

    // Between aligned bases, an insertion is called at the base it precedes.
    ASSERT_TRUE(handler.AddVariantsFromSam(MakeSam(reference.substr(0, 5) + "CC" + reference.substr(5, 5), "5M2I5M", 1)));
    ASSERT_TRUE(handler.HasVariantBin(5));
    EXPECT_TRUE(handler.GetVariantBin(5).front().IsINS());

    // A CIGAR running past the gene is rejected, not read out of bounds.
    EXPECT_FALSE(handler.AddVariantsFromSam(MakeSam(reference.substr(40) + "ACGTA", "15M", 41)));
}

TEST(ExtractVariants, LowQualitySnpKeepsLowQuality) {
    // The read carries G instead of the reference T at position 3 (0-based), with quality '#' (Q2).
    std::string reference = "ACGTACGTAC";
    VariantHandler handler(reference);
    auto sam = MakeSam("ACGGACGTAC", "3M1X6M", 1, 0, "III#IIIIII");

    ASSERT_TRUE(handler.AddVariantsFromSam(sam));
    ASSERT_TRUE(handler.HasVariantBin(3));
    auto& bin = handler.GetVariantBin(3);
    ASSERT_EQ(bin.size(), 1u);
    EXPECT_EQ(bin.front().GetVariant(), 'G');
    EXPECT_EQ(bin.front().QualitySum(), 2u);
}

TEST(Counters, CountBeyondSixteenBits) {
    Variant snp(3, 'G', 'A');
    for (int i = 0; i < 70000; i++) snp.AddObservation(30, i % 2 == 0);
    EXPECT_EQ(snp.Observations(), 70000u);

    SequenceRange range(0, 10);
    for (size_t i = 0; i < 70000; i++) range.AddReadInfo(ReadInfo{ i, 0, 10, true });
    auto coverage = range.CoverageVector();
    ASSERT_EQ(coverage.size(), 10u);
    EXPECT_EQ(coverage[5], 70000u);
}

namespace {
    constexpr char kSecondGene[] = "TTGACCGTAGCATGCAAGTCCGATTGACGTAACGGTCATGCAGTTCAGAC";  // 50 bp

    // Taxon 1's genes 1_1, 1_2, ... (kReference alone by default), loaded.
    struct TinyReference : test::LoadedReference {
        explicit TinyReference(std::vector<std::string> const& genes = { kReference })
            : LoadedReference({ { 1, genes } }, "strain reference") {}
    };
}

TEST(Coverage, AReadOverSeveralGapsMergesAllTheirRanges) {
    // Reads at 0-10, 20-30 and 40-50, then one from 5 to 45 that bridges both gaps.
    SequenceRangeHandler ranges;
    auto read = [&](size_t start, size_t end) {
        SequenceRange range(start, end);
        range.AddReadInfo(ReadInfo{ 0, static_cast<uint32_t>(start), static_cast<uint32_t>(end - start), true });
        ranges.Merge(std::move(range));
    };
    read(0, 10);
    read(20, 30);
    read(40, 50);
    ASSERT_EQ(ranges.Size(), 3u);
    read(5, 45);
    ASSERT_EQ(ranges.Size(), 1u);

    auto cov = ranges.CalculateCoverageVector2();
    ASSERT_EQ(cov.size(), 50u);
    std::vector<std::pair<size_t, uint32_t>> expected = { { 0, 1 }, { 5, 2 }, { 15, 1 }, { 25, 2 }, { 35, 1 }, { 44, 2 }, { 45, 1 }, { 49, 1 } };
    for (auto [pos, depth] : expected) EXPECT_EQ(cov[pos], depth) << "position " << pos;
}

TEST(MSA, DeletionsBecomeGaps) {
    TinyReference ref;
    auto& gene = ref.loader->GetGenome(1).GetGeneOMP(1);
    std::string reference(gene.Sequence());
    StrainLevelContainer strain(gene);
    auto read = reference.substr(0, 20) + reference.substr(22);
    for (size_t i = 0; i < 4; i++) {
        ASSERT_TRUE(strain.AddSam(MakeSam(read, "20M2D28M", 1, i % 2 ? 0x10 : 0), i, true));
    }
    strain.PostProcess(2, 2, 0.0, 0, 0, false);

    MSASequenceItems items;
    items.emplace_back(OptionalMSASequenceItem{ { SharedAlignmentRegion::GetSNPs(strain.GetVariantHandler()), strain.InformativeCoverage() } });
    MSAVector msa(1);
    ASSERT_TRUE(MSA(items, reference, msa, 2, 50));
    std::string row(msa[0].begin(), msa[0].end());
    EXPECT_EQ(row, reference.substr(0, 20) + "--" + reference.substr(22));
}

namespace {
    // The row of one sample's gene 1_1 after adding reads (sequence, CIGAR, position, reverse strand,
    // fragment id), with the default SNP filters of protal (2 reads, strand test, 3 alleles).
    struct OneSample {
        TinyReference ref;
        std::string reference;
        StrainLevelContainer strain;
        OneSample() : reference(ref.loader->GetGenome(1).GetGeneOMP(1).Sequence()), strain(ref.loader->GetGenome(1).GetGeneOMP(1)) {}

        bool Add(std::string const& seq, std::string const& cigar, POS_t pos, bool reverse, size_t fragment, std::string qual = "") {
            return strain.AddSam(MakeSam(seq, cigar, pos, reverse ? 0x10 : 0, qual), fragment, true);
        }

        std::string Row(size_t max_alleles = 3, uint32_t min_cov = 2, uint32_t min_depth = 1) {
            strain.PostProcess(min_cov, min_cov, 0.0, 15, 90, true);
            MSASequenceItems items;
            items.emplace_back(OptionalMSASequenceItem{ { SharedAlignmentRegion::GetSNPs(strain.GetVariantHandler()), strain.InformativeCoverage() } });
            MSAVector msa(1);
            if (!MSA(items, reference, msa, min_cov, 90, 0.0, true, 15, nullptr, nullptr, max_alleles, min_depth)) return "";
            return std::string(msa[0].begin(), msa[0].end());
        }

        std::string WithBase(size_t pos, char base) const {
            auto s = reference;
            s[pos] = base;
            return s;
        }
        char Other(size_t pos) const { return reference[pos] == 'A' ? 'C' : 'A'; }
    };
}

TEST(MSA, ARejectedReadLeavesNothing) {
    // A read with M at a mismatch after an X: it does not fit the gene, so neither its SNP nor its
    // coverage counts.
    OneSample s;
    auto read = s.WithBase(10, s.Other(10));
    read[30] = s.Other(30);
    EXPECT_FALSE(s.Add(read, "10M1X39M", 1, false, 1));
    EXPECT_TRUE(s.strain.GetVariantHandler().GetVariants().empty());
    EXPECT_TRUE(s.strain.InformativeCoverage().empty());
}

TEST(MSA, AFragmentCountsOnce) {
    // Both mates of one fragment carry the SNP at 25 (and overlap there); another fragment shows the
    // reference. One molecule against one: neither allele has the 2 reads a mixture needs, so the
    // site is N (counted twice, the SNP would win 2 to 1). Two fragments with the SNP make it the call.
    OneSample s;
    auto snp = s.WithBase(25, s.Other(25));
    ASSERT_TRUE(s.Add(snp.substr(0, 40), "25M1X14M", 1, false, 1));
    ASSERT_TRUE(s.Add(snp.substr(10), "15M1X24M", 11, true, 1));
    ASSERT_TRUE(s.Add(s.reference, "50M", 1, false, 2));
    auto cov = s.strain.InformativeCoverage();
    EXPECT_EQ(cov[25], 2u);
    EXPECT_EQ(cov[5], 2u);
    EXPECT_EQ(cov[45], 2u);
    EXPECT_EQ(s.Row()[25], 'N');

    OneSample t;
    ASSERT_TRUE(t.Add(snp.substr(0, 40), "25M1X14M", 1, false, 1));
    ASSERT_TRUE(t.Add(snp.substr(10), "15M1X24M", 11, true, 3));
    ASSERT_TRUE(t.Add(t.reference, "50M", 1, false, 2));
    EXPECT_EQ(t.Row()[25], t.Other(25));
}

TEST(MSA, OneReadIsEnoughWhereTheReadsAgree) {
    // One read, with a SNP at 25 and a Q10 mismatch at 30: every position is written from it, the
    // SNP too, but the Q10 base fails the quality filter and is N.
    auto add = [](OneSample& s) {
        auto read = s.WithBase(25, s.Other(25));
        read[30] = s.Other(30);
        std::string qual(50, 'I');
        qual[30] = '+';  // Q10
        return s.Add(read, "25M1X4M1X19M", 1, false, 1, qual);
    };
    OneSample s;
    ASSERT_TRUE(add(s));
    auto const row = s.Row();
    ASSERT_EQ(row.size(), 50u);
    EXPECT_EQ(row.substr(0, 25), s.reference.substr(0, 25));
    EXPECT_EQ(row[25], s.Other(25));
    EXPECT_EQ(row[30], 'N');
    EXPECT_EQ(row.substr(31), s.reference.substr(31));

    // With --msa_min_depth 2, one read writes nothing.
    OneSample t;
    ASSERT_TRUE(add(t));
    EXPECT_EQ(t.Row(3, 2, 2), "");
}

TEST(MSA, ReadsWithoutABaseAreNoReferenceSupport) {
    // 4 reads carry a SNP at 25; 2 reads span it with a deletion of 24-26. The deletion's reads have no
    // base at 25, so they are no support for the reference there: the cell is the SNP, not a mixture.
    OneSample s;
    auto snp = s.WithBase(25, s.Other(25));
    for (size_t i = 0; i < 4; i++) ASSERT_TRUE(s.Add(snp, "25M1X24M", 1, i % 2, i));
    auto deleted = s.reference.substr(0, 24) + s.reference.substr(27);
    for (size_t i = 4; i < 6; i++) ASSERT_TRUE(s.Add(deleted, "24M3D23M", 1, i % 2, i));
    EXPECT_EQ(s.Row()[25], s.Other(25));
}

TEST(MSA, ASnpOnOneStrandOfFewReadsPasses) {
    // Two reads, both reverse, carry the SNP: at this depth a one-strand allele is no strand bias.
    OneSample s;
    auto snp = s.WithBase(25, s.Other(25));
    for (size_t i = 0; i < 2; i++) ASSERT_TRUE(s.Add(snp, "25M1X24M", 1, true, i));
    EXPECT_EQ(s.Row()[25], s.Other(25));

    // 10 reads (5 per strand) show the reference, 8 reverse reads the SNP: strand bias, the SNP fails.
    OneSample t;
    for (size_t i = 0; i < 10; i++) ASSERT_TRUE(t.Add(t.reference, "50M", 1, i % 2, i));
    for (size_t i = 10; i < 18; i++) ASSERT_TRUE(t.Add(snp, "25M1X24M", 1, true, i));
    EXPECT_EQ(t.Row()[25], t.reference[25]);
}

TEST(MSA, TheCallIsTheMostObservedAllele) {
    // 7 reads carry the SNP at Q15, 3 show the reference: the SNP is the call, whatever the qualities
    // (the reference allele's quality is taken as 40).
    OneSample s;
    auto snp = s.WithBase(25, s.Other(25));
    std::string qual(50, 'I');
    qual[25] = '0';  // Q15
    for (size_t i = 0; i < 7; i++) ASSERT_TRUE(s.Add(snp, "25M1X24M", 1, i % 2, i, qual));
    for (size_t i = 7; i < 10; i++) ASSERT_TRUE(s.Add(s.reference, "50M", 1, i % 2, i));
    EXPECT_EQ(s.Row(1)[25], s.Other(25));
}

TEST(MSA, MismatchClustersAtReadEndsAreNotCalled) {
    // A read that starts with mismatches (as reads starting at an insertion are written): its first
    // bases, before 5 matching ones, are neither SNPs nor coverage. A single mismatch near an end is
    // kept.
    OneSample s;
    auto read = s.reference;
    read[0] = s.Other(0);
    read[1] = s.Other(1);
    read[3] = s.Other(3);
    for (size_t i = 0; i < 4; i++) ASSERT_TRUE(s.Add(read, "2X1M1X46M", 1, i % 2, i));
    EXPECT_TRUE(s.strain.GetVariantHandler().GetVariants().empty());
    auto cov = s.strain.InformativeCoverage();
    EXPECT_EQ(cov[3], 0u);
    EXPECT_EQ(cov[4], 4u);

    OneSample t;
    auto one = t.WithBase(3, t.Other(3));
    for (size_t i = 0; i < 4; i++) ASSERT_TRUE(t.Add(one, "3M1X46M", 1, i % 2, i));
    EXPECT_EQ(t.Row()[3], t.Other(3));
}

TEST(MSA, AnInsertionAndASnpAtOnePosition) {
    // 6 reads carry a SNP at 25; 4 reads carry the reference base at 25 behind an insertion of "GG".
    // The insertion fills insertion columns, and the base is a mixture of the SNP and the reference.
    OneSample s;
    auto snp = s.WithBase(25, s.Other(25));
    for (size_t i = 0; i < 6; i++) ASSERT_TRUE(s.Add(snp, "25M1X24M", 1, i % 2, i));
    auto inserted = s.reference.substr(0, 25) + "GG" + s.reference.substr(25);
    for (size_t i = 6; i < 10; i++) ASSERT_TRUE(s.Add(inserted, "25M2I25M", 1, i % 2, i));
    auto row = s.Row();
    ASSERT_EQ(row.size(), 52u);
    EXPECT_EQ(row.substr(25, 2), "GG");
    std::vector<char> both = { s.reference[25], s.Other(25) };
    std::sort(both.begin(), both.end());
    EXPECT_EQ(row[27], IUPACCode(both));
}

TEST(Counters, ObservationsKeepTheirReadsDivergence) {
    // Observations of reads at 99% and 90% identity: a copy limited to 95% keeps the first only, with
    // their strands and qualities.
    Variant snp(3, 'G', 'A');
    for (int i = 0; i < 4; i++) snp.AddObservation(30, i % 2 == 0, DivergenceBin(0.99));
    for (int i = 0; i < 3; i++) snp.AddObservation(20, true, DivergenceBin(0.90));
    auto const own = snp.WithMaxDivergence(MaxDivergenceBin(0.95));
    EXPECT_EQ(own.Observations(), 4u);
    EXPECT_EQ(own.ObservationsForward(), 2u);
    EXPECT_EQ(own.ObservationsReverse(), 2u);
    EXPECT_EQ(own.QualitySum(), 120u);
    EXPECT_EQ(snp.WithMaxDivergence(MaxDivergenceBin(0)).Observations(), 7u);
    EXPECT_EQ(DivergenceBin(1.0), 0u);
    EXPECT_EQ(DivergenceBin(0.0), 127u);
    EXPECT_EQ(MaxDivergenceBin(0.95), DivergenceBin(0.95));
}

TEST(MSA, TheRowIsMadeOfTheTaxonsOwnReads) {
    // 6 reads of the strain (99% identity) show the reference at 25; 4 reads of a relative (90%) carry
    // a SNP there and cover 40-49 alone. From every read, 25 is a mixture and 40-49 are called; from
    // the reads of at least 95% identity, 25 is the reference base and 40-49 have no read.
    OneSample s;
    for (size_t i = 0; i < 6; i++) {
        ASSERT_TRUE(s.strain.AddSam(MakeSam(s.reference.substr(0, 40), "40M", 1, i % 2 ? 0x10 : 0), i, true, 0.99));
    }
    auto relative = s.WithBase(25, s.Other(25));
    for (size_t i = 6; i < 10; i++) {
        ASSERT_TRUE(s.strain.AddSam(MakeSam(relative, "25M1X24M", 1, i % 2 ? 0x10 : 0), i, true, 0.90));
    }
    s.strain.PostProcess(2, 2, 0.0, 15, 90, true);
    auto row = [&](double min_identity) {
        MSASequenceItems items;
        items.emplace_back(OptionalMSASequenceItem{ s.strain.MSAItem(min_identity, 2, 0.15, 15, 90, true) });
        MSAVector msa(1);
        if (!MSA(items, s.reference, msa, 2, 90, 0.15, true, 15, nullptr, nullptr, 3)) return std::string();
        return std::string(msa[0].begin(), msa[0].end());
    };
    auto const all = row(0), own = row(0.95);
    ASSERT_EQ(all.size(), s.reference.size());
    ASSERT_EQ(own.size(), s.reference.size());
    std::vector<char> both = { s.reference[25], s.Other(25) };
    std::sort(both.begin(), both.end());
    EXPECT_EQ(all[25], IUPACCode(both));
    EXPECT_EQ(all.substr(40), s.reference.substr(40));
    EXPECT_EQ(own[25], s.reference[25]);
    EXPECT_EQ(own.substr(0, 25), s.reference.substr(0, 25));
    EXPECT_EQ(own.substr(40), std::string(10, '-'));

    // The informative coverage (the MSA's depth, and .meta.tsv's) counts the same reads.
    auto const [bins, coverage] = s.strain.MSAItem(0.95, 2, 0.15, 15, 90, true);
    EXPECT_EQ(coverage[25], 6u);
    EXPECT_TRUE(bins.empty()) << "no allele of the strain's reads differs from the reference";
}

TEST(MSA, EachSampleTakesItsOwnMinimumAlleleFrequency) {
    // Two samples with the same reads: 3 of 10 carry a SNP at position 10. With no minimum allele
    // frequency, the SNP and the reference base are written as an IUPAC code; with 0.5 (the other
    // sample, as for a noisier read type), the SNP does not pass.
    TinyReference ref;
    auto& gene = ref.loader->GetGenome(1).GetGeneOMP(1);
    std::string const reference(gene.Sequence());
    std::string snp_read = reference;
    snp_read[10] = reference[10] == 'A' ? 'C' : 'A';
    std::vector<std::unique_ptr<StrainLevelContainer>> strains;
    MSASequenceItems items;
    for (int sample = 0; sample < 2; sample++) {
        strains.push_back(std::make_unique<StrainLevelContainer>(gene));
        auto& strain = *strains.back();
        for (size_t i = 0; i < 10; i++) {
            bool const snp = i < 3;
            ASSERT_TRUE(strain.AddSam(MakeSam(snp ? snp_read : reference, snp ? "10M1X39M" : "50M", 1, i % 2 ? 0x10 : 0), i, true));
        }
        strain.PostProcess(2, 2, 0.0, 0, 0, false);
        items.emplace_back(OptionalMSASequenceItem{ { SharedAlignmentRegion::GetSNPs(strain.GetVariantHandler()), strain.InformativeCoverage() } });
    }
    MSAVector msa(2);
    ASSERT_TRUE(MSA(items, reference, msa, 2, 0, std::vector<double>{ 0.0, 0.5 }, false, 0, nullptr, nullptr, 2));
    std::string const loose(msa[0].begin(), msa[0].end()), strict(msa[1].begin(), msa[1].end());
    ASSERT_EQ(loose.size(), reference.size());
    ASSERT_EQ(reference[10], 'G');
    EXPECT_EQ(loose[10], 'R') << "A and G";
    EXPECT_EQ(strict[10], 'G') << "the reference's base alone";
    EXPECT_EQ(loose.substr(11), strict.substr(11));
}

TEST(Abundance, BlendedDepthIsUnbiasedAtLowCoverage) {
    // 200 genes of 1 kb, 100 bp reads: the number of reads per gene is Poisson(10 x depth).
    std::mt19937 rng(7);
    for (double depth : { 0.03, 0.08, 0.2, 0.5, 1.0, 3.0 }) {
        double sum = 0, median_sum = 0;
        int const replicates = 40;
        for (int r = 0; r < replicates; r++) {
            std::poisson_distribution<int> reads(10 * depth);
            std::vector<double> gene_depths;
            size_t mapped_bases = 0;
            for (int g = 0; g < 200; g++) {
                int n = reads(rng);
                if (n == 0) continue;
                gene_depths.push_back(n * 100 / 1000.0);
                mapped_bases += n * 100;
            }
            std::sort(gene_depths.begin(), gene_depths.end());
            double median = profiler::Taxon::Median(gene_depths);
            sum += profiler::Taxon::BlendedDepth(median, mapped_bases, 200 * 1000, gene_depths.size(), 200);
            median_sum += median;
        }
        double estimate = sum / replicates;
        EXPECT_NEAR(estimate / depth, 1.0, 0.15) << "depth " << depth << ", median alone gives " << median_sum / replicates;
    }
}

TEST(Abundance, BlendedDepthHasNoStep) {
    // Around both ends of both ramps, a tiny change of the input changes the estimate a little.
    size_t const expected_genes = 1000000;
    auto at = [&](double hit_fraction, double median) {
        return profiler::Taxon::BlendedDepth(median, 400000, 1000000, static_cast<size_t>(hit_fraction * expected_genes), expected_genes);
    };
    for (double f : { 0.80, 0.95 }) {
        EXPECT_NEAR(at(f - 1e-5, 0.3), at(f + 1e-5, 0.3), 1e-4) << "hit fraction " << f;
    }
    for (double m : { 0.5, 1.0 }) {
        EXPECT_NEAR(at(0.5, m - 1e-5), at(0.5, m + 1e-5), 1e-4) << "median " << m;
    }
    for (double f : { 0.25, 0.5 }) {
        EXPECT_NEAR(at(f - 1e-5, 2.0), at(f + 1e-5, 2.0), 1e-3) << "hit fraction " << f << " at a high median";
    }
    EXPECT_DOUBLE_EQ(at(0.5, 0.1), 0.4);   // low coverage: all expected genes
    EXPECT_DOUBLE_EQ(at(1.0, 0.7), 0.7);   // every gene hit: the median
    EXPECT_DOUBLE_EQ(at(0.5, 2.0), 2.0);   // high median depth: the median
    EXPECT_DOUBLE_EQ(at(0.02, 2.0), 0.4);  // ...but not from a couple of reads on a short gene
}

TEST(Abundance, StrongOwnEvidenceNeedsDepthFromTheTaxonsOwnReads) {
    // What unreported_species.tsv lists a taxon for, whatever its score: its own reads give 1x or
    // more, on 90% of its genes (here its one gene), and most of its bases are its own reads'.
    TinyReference ref;
    std::string const reference(ref.loader->GetGenome(1).GetGeneOMP(1).Sequence());
    profiler::MicrobialProfile profile(*ref.loader);
    int read = 0;  // each read its own fragment
    profile.SetDepthIdentityMargin(0.04);
    auto own = MakeSam(reference.substr(0, 20), "20M", 1);
    for (int i = 0; i < 2; i++) ASSERT_TRUE(profile.AddSam(1, 1, own, 1.0, true, ++read));
    auto const& taxon = profile.GetTaxa().at(1);
    EXPECT_FALSE(profiler::StrongOwnEvidence(taxon)) << "0.8x";
    ASSERT_TRUE(profile.AddSam(1, 1, own, 1.0, true, ++read));
    EXPECT_TRUE(profiler::StrongOwnEvidence(taxon)) << "1.2x";
    auto relative = MakeSam(reference.substr(20, 20), "15M5X", 21);  // identity 0.75
    for (int i = 0; i < 4; i++) ASSERT_TRUE(profile.AddSam(1, 1, relative, 1.0, true, ++read));
    EXPECT_FALSE(profiler::StrongOwnEvidence(taxon)) << "most bases from a relative's reads";
}

TEST(Abundance, DepthCountsOnlyTheTaxonsOwnReads) {
    TinyReference ref;
    std::string reference(ref.loader->GetGenome(1).GetGeneOMP(1).Sequence());
    profiler::MicrobialProfile profile(*ref.loader);
    int read = 0;  // each read its own fragment
    profile.SetDepthIdentityMargin(0.04);
    auto own = MakeSam(reference.substr(0, 20), "20M", 1);
    auto relative = MakeSam(reference.substr(20, 20), "15M5X", 21);  // identity 0.75
    for (auto const* sam : { &own, &own, &relative, &relative }) {
        ASSERT_TRUE(profile.AddSam(1, 1, *sam, 1.0, true, ++read));
    }
    auto& taxon = profile.GetTaxa().at(1);
    // One 50 bp gene, every gene hit: the depth is the median gene depth of the own reads.
    EXPECT_NEAR(taxon.VerticalCoverage(true), 40.0 / 50, 1e-9);
    EXPECT_NEAR(taxon.LowIdentityShare(), 0.5, 1e-9);

    profile.SetDepthIdentityMargin(1);  // every read counts
    EXPECT_NEAR(taxon.VerticalCoverage(), 80.0 / 50, 1e-9);
    EXPECT_NEAR(taxon.LowIdentityShare(), 0.0, 1e-9);
}

TEST(Abundance, AFragmentsOverlapCountsOnceInTheDepth) {
    // On one 50 bp gene: a pair whose mates overlap by 10 bases (1-30 and 21-50) adds 50 bases, as the strain MSA
    // counts it, not 60; a pair whose mates do not overlap 40; a single read 50. The depth does not depend on how
    // often a library's mates overlap; the records' own bases (rates per aligned base) still count both mates.
    TinyReference ref;
    std::string reference(ref.loader->GetGenome(1).GetGeneOMP(1).Sequence());
    profiler::MicrobialProfile profile(*ref.loader);
    auto mate = [&](size_t start, size_t length, FLAG_t flag) {
        return MakeSam(reference.substr(start, length), std::to_string(length) + "M", static_cast<POS_t>(start + 1), flag);
    };
    ASSERT_TRUE(profile.AddSam(1, 1, mate(0, 30, 0x1 | 0x40), 1.0, true, 1));
    ASSERT_TRUE(profile.AddSam(1, 1, mate(20, 30, 0x1 | 0x80 | 0x10), 1.0, true, 1));
    ASSERT_TRUE(profile.AddSam(1, 1, mate(0, 20, 0x1 | 0x40), 1.0, true, 2));
    ASSERT_TRUE(profile.AddSam(1, 1, mate(30, 20, 0x1 | 0x80 | 0x10), 1.0, true, 2));
    ASSERT_TRUE(profile.AddSam(1, 1, mate(0, 50, 0), 1.0, true, 3));
    auto const& taxon = profile.GetTaxa().at(1);
    auto const& gene = taxon.GetGenes().at(1);
    EXPECT_EQ(gene.m_fragment_bases, 140u);
    EXPECT_EQ(gene.m_mapped_length, 150u);
    EXPECT_NEAR(gene.VerticalCoverage(), 140.0 / 50, 1e-9);
    EXPECT_NEAR(taxon.VerticalCoverage(), 140.0 / 50, 1e-9);
    EXPECT_EQ(taxon.TotalLength(), 150u);
    EXPECT_EQ(taxon.Fragments(), 3u);
}

TEST(Abundance, TheMsaTakesReadsByAStricterMarginThanTheDepth) {
    // Reads 5% below the best ones (a distant strain's, or a relative's) count towards the depth with
    // the default margin (0.08); the strain MSA takes reads within --msa_identity_margin (0.04) only.
    TinyReference ref;
    std::string reference(ref.loader->GetGenome(1).GetGeneOMP(1).Sequence());
    profiler::MicrobialProfile profile(*ref.loader);
    int read = 0;  // each read its own fragment
    profile.SetDepthIdentityMargin(0.08);
    auto own = MakeSam(reference.substr(0, 20), "20M", 1);
    auto strain = MakeSam(reference.substr(20, 20), "19M1X", 21);  // identity 0.95
    for (auto const* sam : { &own, &own, &strain, &strain }) {
        ASSERT_TRUE(profile.AddSam(1, 1, *sam, 1.0, true, ++read));
    }
    auto& taxon = profile.GetTaxa().at(1);
    EXPECT_DOUBLE_EQ(taxon.TopIdentity(), 1.0);
    EXPECT_NEAR(taxon.OwnIdentityThreshold(1), 0.92, 1e-9);
    EXPECT_NEAR(taxon.VerticalCoverage(true), 80.0 / 50, 1e-9);  // both kinds of reads
    EXPECT_NEAR(taxon.IdentityThreshold(0.04), 0.96, 1e-9);       // the MSA's: the strain's reads out
    EXPECT_EQ(taxon.IdentityThreshold(1), 0);
}

TEST(Abundance, TheDepthMarginScalesWithTheGenesConservation) {
    // Two 50 bp genes of one species: gene 1 conserved (factor 0.4: margin 0.03 + 0.05 x 0.4 = 0.05),
    // gene 2 fast (1.6: 0.11). Reads 10% below the best on the fast gene are the species' own (a
    // strain's) and count; a read 6% below on the conserved gene is a relative's and does not.
    TinyReference ref({ kReference, kSecondGene });
    GenomeLoader& loader = *ref.loader;
    auto own = MakeSam(kReference, "50M", 1);
    auto relative = MakeSam(kReference, "47M3X", 1);  // identity 0.94
    auto strain = MakeSam(kSecondGene, "45M5X", 1);    // identity 0.90
    // A taxon takes the loader's factors, and whether they scale the margin, when its first read comes.
    auto fill = [&](profiler::MicrobialProfile& profile) -> profiler::Taxon& {
        int read = 0;  // each read its own fragment
        profile.SetDepthIdentityMargin(0.08);
        for (auto const* sam : { &own, &own, &relative }) EXPECT_TRUE(profile.AddSam(1, 1, *sam, 1.0, true, ++read));
        for (int i = 0; i < 2; i++) EXPECT_TRUE(profile.AddSam(1, 2, strain, 1.0, true, ++read));
        return profile.GetTaxa().at(1);
    };
    profiler::MicrobialProfile without(loader);
    auto& unscaled = fill(without);
    ASSERT_DOUBLE_EQ(unscaled.TopIdentity(), 1.0);

    // Without factors: 0.08 on both genes, the strain's reads out and the relative's in.
    EXPECT_NEAR(unscaled.OwnIdentityThreshold(1), 0.92, 1e-9);
    EXPECT_NEAR(unscaled.OwnIdentityThreshold(2), 0.92, 1e-9);
    unscaled.VerticalCoverage(true);
    EXPECT_NEAR(unscaled.LowIdentityShare(), 100.0 / 250, 1e-9);
    EXPECT_EQ(unscaled.ConservationPattern(), (std::pair<double, double>{ 0.0, 0.5 }));

    gene_conservation::Table factors;
    factors.Set(1, 0.4);
    factors.Set(2, 1.6);
    loader.SetGeneConservation(factors);
    // The factors alone (--gene_conservation none, the default) give the conservation pattern, not the margin:
    // the conserved gene's depth is 3, the fast one's 2.
    EXPECT_NEAR(unscaled.OwnIdentityThreshold(1), 0.92, 1e-9);
    auto const [ratio, share] = unscaled.ConservationPattern();
    EXPECT_NEAR(ratio, std::log2(3.001 / 2.001), 1e-9);
    EXPECT_DOUBLE_EQ(share, 0.5);

    loader.SetScaleDepthMargin(true);  // --gene_conservation db
    profiler::MicrobialProfile profile(loader);
    auto& taxon = fill(profile);
    EXPECT_NEAR(taxon.OwnIdentityThreshold(1), 0.95, 1e-6);
    EXPECT_NEAR(taxon.OwnIdentityThreshold(2), 0.89, 1e-6);
    EXPECT_NEAR(taxon.VerticalCoverage(true), 100.0 / 50, 1e-9);  // two reads' worth on each gene
    EXPECT_NEAR(taxon.LowIdentityShare(), 50.0 / 250, 1e-9);       // the relative's read only

    profile.SetDepthIdentityMargin(1);  // every read counts, whatever the gene
    EXPECT_EQ(taxon.OwnIdentityThreshold(1), 0);
    EXPECT_NEAR(taxon.LowIdentityShare(), 0.0, 1e-9);
}

TEST(Abundance, TheRecordsDepthOnConservedAgainstFastGenes) {
    // conserved_fast_record_ratio: every best record (NoteRecord, before the MAPQ filter) on gene 1 (factor 0.4)
    // against gene 2 (1.6), both 50 bp: 150 bases on gene 1 (one of the records MAPQ 0), 50 on gene 2.
    TinyReference ref({ kReference, kSecondGene });
    GenomeLoader& loader = *ref.loader;
    auto own = MakeSam(kReference, "50M", 1);
    auto ambiguous = MakeSam(kReference, "50M", 1);
    ambiguous.m_mapq = 0;
    auto fast = MakeSam(kSecondGene, "45M5X", 1);
    auto fill = [&](profiler::MicrobialProfile& profile) -> profiler::Taxon const& {
        for (auto const* sam : { &own, &own, &ambiguous }) profile.NoteRecord(1, 1, *sam);
        profile.NoteRecord(1, 2, fast);
        EXPECT_TRUE(profile.AddSam(1, 1, own, 1.0));
        profile.ApplyRecordEvidence();
        return profile.GetTaxa().at(1);
    };
    std::map<std::string, double> features;

    profiler::MicrobialProfile without(loader);  // no factors: 0
    for (auto const& [name, value] : profiler::TaxonFeatures(fill(without))) features[name] = value;
    EXPECT_EQ(features.at("conserved_fast_record_ratio"), 0.0);
    EXPECT_EQ(features.at("conserved_fast_kept_ratio"), 0.0);

    gene_conservation::Table factors;
    factors.Set(1, 0.4);
    factors.Set(2, 1.6);
    loader.SetGeneConservation(factors);  // for the features; the margin is not scaled
    profiler::MicrobialProfile profile(loader);
    auto const& taxon = fill(profile);
    EXPECT_NEAR(taxon.RecordConservedFastRatio(), std::log2(3.001 / 1.001), 1e-9);
    for (auto const& [name, value] : profiler::TaxonFeatures(taxon)) features[name] = value;
    EXPECT_NEAR(features.at("conserved_fast_record_ratio"), std::log2(3.001 / 1.001), 1e-9);
    // Only the one record the filters keep is a hit: the hit genes' pattern sees gene 1 alone, and the kept records'
    // depth is 1 on gene 1 and 0 on gene 2.
    EXPECT_EQ(taxon.ConservationPattern(), (std::pair<double, double>{ 0.0, 1.0 }));
    EXPECT_NEAR(features.at("conserved_fast_kept_ratio"), std::log2(1.001 / 0.001), 1e-9);
}

// tsl::sparse_map copies the values of a bucket on every insert into it unless they move without
// throwing.
static_assert(std::is_nothrow_move_constructible_v<profiler::Gene>);
static_assert(std::is_nothrow_move_constructible_v<profiler::Taxon>);

TEST(Abundance, ReleasingReadDataKeepsWhatLaterStagesRead) {
    // Once a sample's outputs are written, its reads' records are freed. A taxon that enters the strain MSAs keeps what
    // the strain stage reads of its genes (Gene::KeepForMSAs), the same as the records gave; depth and counters stay.
    TinyReference ref;
    std::string reference(ref.loader->GetGenome(1).GetGeneOMP(1).Sequence());
    // 40 bases with one SNP: 97.5% identical, within the 0.04 identity margin of the taxon's best reads.
    std::string snp = reference.substr(0, 40);
    snp[5] = snp[5] == 'A' ? 'C' : 'A';
    profiler::MSAReleaseParameters msa;
    msa.identity_margin = 0.04;
    msa.item = { 0, 1, 0.15, 15, 0, false };
    msa.count_keys = { { 2, 60 }, { 0, 0 } };
    for (bool enters : { true, false }) {
        profiler::MicrobialProfile profile(*ref.loader);
        profile.SetDepthIdentityMargin(0.04);
        auto first = MakeSam(reference.substr(0, 20), "20M", 1);
        auto second = MakeSam(reference.substr(20, 20), "20M", 21);
        auto variant = MakeSam(snp, "5M1X34M", 1, 0x10);
        int read_id = 0;  // each record a fragment of its own: a fragment's second record skips what its first covered
        for (auto const* sam : { &first, &first, &second, &variant, &variant }) {
            ASSERT_TRUE(profile.AddSam(1, 1, *sam, 1.0, true, read_id++, false));
        }
        auto& taxon = profile.GetTaxa().at(1);
        double const depth = taxon.VerticalCoverage();
        size_t const length = taxon.TotalLength();
        auto& gene = taxon.GetGenes().at(1);
        ASSERT_GT(gene.GetStrainLevel().GetSequenceRangeHandler().Size(), 0u);
        auto item = msa.item;
        item.min_identity = taxon.IdentityThreshold(msa.identity_margin);
        auto const before = gene.MSAItem(item);
        ASSERT_FALSE(before.first.empty());
        auto const counts = gene.AlleleSNPCounts(0, 0);
        size_t const covered = gene.Coverage();

        taxon.ReleaseReadData(enters ? std::optional(msa) : std::nullopt);
        EXPECT_TRUE(gene.m_read_identities.empty());
        EXPECT_EQ(gene.GetStrainLevel().GetSequenceRangeHandler().Size(), 0u);
        EXPECT_EQ(gene.m_msa_evidence != nullptr, enters);
        if (enters) {
            auto const after = gene.MSAItem(item);
            EXPECT_EQ(after.second, before.second);
            ASSERT_EQ(after.first.size(), before.first.size());
            for (size_t b = 0; b < before.first.size(); b++) {
                ASSERT_EQ(after.first[b].size(), before.first[b].size());
                for (size_t v = 0; v < before.first[b].size(); v++) {
                    auto const& x = before.first[b][v];
                    auto const& y = after.first[b][v];
                    EXPECT_TRUE(y.IsCompact());
                    EXPECT_EQ(y.Position(), x.Position());
                    EXPECT_EQ(y.GetVariant(), x.GetVariant());
                    EXPECT_EQ(y.ObservationsForward(), x.ObservationsForward());
                    EXPECT_EQ(y.ObservationsReverse(), x.ObservationsReverse());
                    EXPECT_EQ(y.QualitySum(), x.QualitySum());
                    EXPECT_EQ(y.MeanQuality(), x.MeanQuality());
                    EXPECT_EQ(y.GetValid(), x.GetValid());
                }
            }
            auto const kept = gene.AlleleSNPCounts(0, 0);
            EXPECT_EQ(std::tie(kept.filtered, kept.mono, kept.bi, kept.tri, kept.tetra),
                      std::tie(counts.filtered, counts.mono, counts.bi, counts.tri, counts.tetra));
            EXPECT_EQ(gene.Coverage(), covered);
            // What it did not keep cannot be asked for (it would be read from the records, which are gone).
            EXPECT_THROW(gene.AlleleSNPCounts(5, 5), std::logic_error);
            auto other = item;
            other.min_cov = 7;
            EXPECT_THROW(gene.MSAItem(other), std::logic_error);
        }
        EXPECT_EQ(taxon.VerticalCoverage(), depth);
        EXPECT_EQ(taxon.TotalLength(), length);
        EXPECT_EQ(taxon.TotalHits(), 5u);
    }
}

// A packed MSA item unpacks to the alleles the strain MSA reads: per allele its kind, position, bases, observations per
// strand, quality sum and flags, INDEL bases included; and the coverage, in whatever width its largest value needs.
TEST(PackedMSAItem, UnpacksWhatTheMSAReads) {
    std::string inserted = "GT";
    Variant snp(7, 'T', 'A');
    snp.AddObservation(30, true);
    snp.AddObservation(12, false);
    snp.SetValid(false);
    Variant reference(7, 'A', 'A');
    reference.SetObservations(5, 2);
    reference.SetMajorAllele(true);
    Variant insertion(VariantType::INS, 9, 'C', inserted);
    insertion.AddObservation(35, true);
    std::pair<VariantVec, CoverageVec> item{ { { snp, reference }, { insertion } }, {} };
    for (CoverageVec coverage : { CoverageVec{ 0, 3, 255 }, CoverageVec{ 1, 256, 65535 }, CoverageVec{ 2, 65536, 4000000000u }, CoverageVec{} }) {
        item.second = coverage;
        PackedMSAItem packed(item);
        auto const back = packed.Unpack();
        EXPECT_EQ(back.second, coverage);
        ASSERT_EQ(back.first.size(), item.first.size());
        for (size_t b = 0; b < item.first.size(); b++) {
            ASSERT_EQ(back.first[b].size(), item.first[b].size());
            for (size_t a = 0; a < item.first[b].size(); a++) {
                auto const& x = item.first[b][a];
                auto const& y = back.first[b][a];
                EXPECT_TRUE(y.IsCompact());
                EXPECT_EQ(y.GetType(), x.GetType());
                EXPECT_EQ(y.Position(), x.Position());
                EXPECT_EQ(y.GetReference(), x.GetReference());
                EXPECT_EQ(y.GetVariant(), x.GetVariant());
                EXPECT_EQ(y.GetStructuralSize(), x.GetStructuralSize());
                EXPECT_EQ(y.StructuralBases() != nullptr, x.StructuralBases() != nullptr);
                if (x.StructuralBases()) EXPECT_EQ(*y.StructuralBases(), *x.StructuralBases());
                EXPECT_EQ(y.ObservationsForward(), x.ObservationsForward());
                EXPECT_EQ(y.ObservationsReverse(), x.ObservationsReverse());
                EXPECT_EQ(y.QualitySum(), x.QualitySum());
                EXPECT_EQ(y.MeanQuality(), x.MeanQuality());
                EXPECT_EQ(y.GetValid(), x.GetValid());
                EXPECT_EQ(y.IsMajorAllele(), x.IsMajorAllele());
            }
        }
    }
}

namespace {
    // A model that scores every taxon 0.5: one leaf (on present_genes, which it does not split on).
    profiler::TaxonFilterObj HalfModel(double knob) {
        std::string const xml = R"(<?xml version="1.0" encoding="UTF-8"?>
<PMML version="4.4">
 <Header/>
 <DataDictionary>
  <DataField name="truth" optype="categorical" dataType="string"><Value value="FALSE"/><Value value="TRUE"/></DataField>
  <DataField name="present_genes" optype="continuous" dataType="double"/>
 </DataDictionary>
 <TreeModel functionName="classification" splitCharacteristic="binarySplit">
  <MiningSchema>
   <MiningField name="truth" usageType="predicted"/>
   <MiningField name="present_genes"/>
  </MiningSchema>
  <Node score="TRUE"><True/><ScoreDistribution value="FALSE" recordCount="1"/><ScoreDistribution value="TRUE" recordCount="1"/></Node>
 </TreeModel>
</PMML>
)";
        return profiler::TaxonFilterObj(cpmml::Model::from_string(xml), knob);
    }
}

// A sample's taxa that do not enter the strain MSAs are dropped once its outputs are written, unless every taxon's
// numbers are asked for (--taxon_statistics); without strain MSAs all of them are.
TEST(Abundance, ReleasingReadDataDropsTheTaxaNoMsaTakes) {
    TinyReference ref;
    std::string reference(ref.loader->GetGenome(1).GetGeneOMP(1).Sequence());
    profiler::MSAReleaseParameters msa;
    msa.identity_margin = 0.04;
    msa.item = { 0, 2, 0.15, 15, 90, true };
    for (int knob : { 0, 1 }) {  // the model scores 0.5: at knob 0 every taxon enters, at 1 none does
        for (bool keep_all : { false, true }) {
            for (bool strains : { true, false }) {
                profiler::MicrobialProfile profile(*ref.loader);
                auto sam = MakeSam(reference.substr(0, 20), "20M", 1);
                ASSERT_TRUE(profile.AddSam(1, 1, sam, 1.0, true, 0, false));
                profile.ReleaseReadData(HalfModel(knob), strains ? std::optional(msa) : std::nullopt, keep_all);
                bool const kept = keep_all || (strains && knob == 0);
                EXPECT_EQ(profile.GetTaxa().contains(1), kept) << "knob " << knob << " keep_all " << keep_all << " strains " << strains;
            }
        }
    }
}
