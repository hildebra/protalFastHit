// Unit tests for the strain output: variant calls from CIGARs (N bases, deletions, indels at the
// ends of an alignment), counters beyond 16 bits, deletion columns in the MSA, and the depth
// estimate behind abundances.
#include <gtest/gtest.h>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <type_traits>
#include <vector>
#include <unistd.h>
#include "Profiling/Profiler.h"
#include "Profiling/Strain.h"

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

    // Three reads cover position 10: one N, one T, one reference G.
    handler.PostProcessSNPs(CoverageVec(50, 3), 1, 1, 0.0, 0, 0, false);
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
    struct TinyReference {
        std::filesystem::path dir;
        std::unique_ptr<GenomeLoader> loader;

        TinyReference() {
            dir = std::filesystem::temp_directory_path() / ("protal_strain_test_" + std::to_string(::getpid()));
            std::filesystem::create_directories(dir);
            std::string header = ">1_1\n";
            std::string gene = kReference;
            { std::ofstream fna(dir / "reference.fna"); fna << header << gene << '\n'; }
            { std::ofstream map(dir / "reference.map"); map << "1\t1\t" << header.size() << '\t' << header.size() + gene.size() << '\n'; }
            loader = std::make_unique<GenomeLoader>((dir / "reference.fna").string(), (dir / "reference.map").string());
            loader->LoadAllGenomes();
        }
        ~TinyReference() { std::filesystem::remove_all(dir); }
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

    ranges.CalculateCoverageVector();
    EXPECT_EQ(ranges.CalculateCoverageVector(), cov);  // recomputed, not appended to
}

TEST(MSA, DeletionsBecomeGaps) {
    TinyReference ref;
    auto& gene = ref.loader->GetGenome(1).GetGeneOMP(1);
    std::string reference = gene.Sequence();
    StrainLevelContainer strain(gene);
    auto read = reference.substr(0, 20) + reference.substr(22);
    for (size_t i = 0; i < 4; i++) {
        ASSERT_TRUE(strain.AddSam(MakeSam(read, "20M2D28M", 1, i % 2 ? 0x10 : 0), i, true));
    }
    strain.PostProcess(2, 2, 0.0, 0, 0, false);

    MSASequenceItems items;
    items.emplace_back(OptionalMSASequenceItem{ { SharedAlignmentRegion::GetSNPs(strain.GetVariantHandler()), strain.GetSequenceRangeHandler() } });
    MSAVector msa(1);
    ASSERT_TRUE(MSA(items, reference, msa, 2, 50));
    std::string row(msa[0].begin(), msa[0].end());
    EXPECT_EQ(row, reference.substr(0, 20) + "--" + reference.substr(22));
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

TEST(Abundance, DepthCountsOnlyTheTaxonsOwnReads) {
    TinyReference ref;
    std::string reference = ref.loader->GetGenome(1).GetGeneOMP(1).Sequence();
    profiler::MicrobialProfile profile(*ref.loader);
    profile.SetDepthIdentityMargin(0.04);
    auto own = MakeSam(reference.substr(0, 20), "20M", 1);
    auto relative = MakeSam(reference.substr(20, 20), "15M5X", 21);  // identity 0.75
    for (auto const* sam : { &own, &own, &relative, &relative }) {
        ASSERT_TRUE(profile.AddSam(1, 1, *sam, 1.0));
    }
    auto& taxon = profile.GetTaxa().at(1);
    // One 50 bp gene, every gene hit: the depth is the median gene depth of the own reads.
    EXPECT_NEAR(taxon.VerticalCoverage(true), 40.0 / 50, 1e-9);
    EXPECT_NEAR(taxon.LowIdentityShare(), 0.5, 1e-9);

    profile.SetDepthIdentityMargin(1);  // every read counts
    EXPECT_NEAR(taxon.VerticalCoverage(), 80.0 / 50, 1e-9);
    EXPECT_NEAR(taxon.LowIdentityShare(), 0.0, 1e-9);
}

// tsl::sparse_map copies the values of a bucket on every insert into it unless they move without
// throwing.
static_assert(std::is_nothrow_move_constructible_v<profiler::Gene>);
static_assert(std::is_nothrow_move_constructible_v<profiler::Taxon>);

TEST(Abundance, ReleasingReadDataKeepsWhatLaterStagesRead) {
    // Once a sample's outputs are written, its reads' identities are freed, and its variants and
    // read ranges unless the strain stage needs them; depth and counters stay.
    TinyReference ref;
    std::string reference = ref.loader->GetGenome(1).GetGeneOMP(1).Sequence();
    for (bool keep_strain_data : { true, false }) {
        profiler::MicrobialProfile profile(*ref.loader);
        profile.SetDepthIdentityMargin(0.04);
        auto first = MakeSam(reference.substr(0, 20), "20M", 1);
        auto second = MakeSam(reference.substr(20, 20), "20M", 21);
        for (auto const* sam : { &first, &first, &second }) {
            ASSERT_TRUE(profile.AddSam(1, 1, *sam, 1.0, true, 0, false));
        }
        auto& taxon = profile.GetTaxa().at(1);
        double const depth = taxon.VerticalCoverage();
        size_t const length = taxon.TotalLength();
        ASSERT_GT(taxon.GetGenes().at(1).GetStrainLevel().GetSequenceRangeHandler().Size(), 0u);

        taxon.ReleaseReadData(keep_strain_data);
        auto const& gene = taxon.GetGenes().at(1);
        EXPECT_TRUE(gene.m_read_identities.empty());
        EXPECT_EQ(gene.GetStrainLevel().GetSequenceRangeHandler().Size() > 0, keep_strain_data);
        EXPECT_EQ(taxon.VerticalCoverage(), depth);
        EXPECT_EQ(taxon.TotalLength(), length);
        EXPECT_EQ(taxon.TotalHits(), 3u);
    }
}

TEST(ModelFeatures, NamesAreUniqueAndValuesKeepTheirPrecision) {
    Genome no_genome(0);
    auto features = profiler::TaxonFeatures(profiler::Taxon(no_genome));
    std::set<std::string> names;
    for (auto const& [name, _] : features) EXPECT_TRUE(names.insert(name).second) << name << " twice";
    for (auto const* name : { "RAF0", "RA4", "su_rate_ref", "lu_rate_ref", "lsu_rate_ref", "mean_mapq", "lu_per_read" }) {
        EXPECT_TRUE(names.contains(name)) << name;
    }

    EXPECT_EQ(profiler::FeatureString(0.5), "0.5");
    EXPECT_EQ(profiler::FeatureString(3), "3");
    EXPECT_EQ(std::stod(profiler::FeatureString(1.25e-7)), 1.25e-7);
    EXPECT_EQ(std::stod(profiler::FeatureString(0.1 + 0.2)), 0.1 + 0.2);
}
