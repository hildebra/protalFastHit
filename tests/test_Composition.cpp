// Unit tests for the profile's composition (Composition.h): the scanned reads' SAM header line, the species' genome sizes
// in species_priors.tsv, the shares the called species explain, the average genome size and its quantiles, the unknown
// share and the species it would be, and the profile that ends with the unknown share ("?").
#include <gtest/gtest.h>
#include <cmath>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include "Profiling/Profiler.h"
#include "Profiling/Composition.h"
#include "IO/ScannedReads.h"
#include "SequenceUtils/SpeciesPriors.h"
#include "TestReference.h"

using namespace protal;
using protal::test::RandomSequence;
using Reference = protal::test::LoadedReference;

namespace {
    // A record of the gene's own bases from 0-based `start`.
    std::string Record(Reference const& ref, std::string const& qname, int flag, uint32_t taxid, uint32_t gene, int start, int length) {
        std::string const seq = ref.genes.at(taxid)[gene - 1].substr(static_cast<size_t>(start), static_cast<size_t>(length));
        return qname + '\t' + std::to_string(flag) + '\t' + std::to_string(taxid) + "_" + std::to_string(gene) + '\t' + std::to_string(start + 1) +
               "\t60\t" + std::to_string(length) + "M\t*\t0\t0\t" + seq + '\t' + std::string(seq.size(), 'I') + "\tZA:Z:*\n";
    }

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

    constexpr char kTaxonomy[] =
            "id\tparent_id\texternal_id\tname\trank\tlevel\trep_genome\n"
            "1\t6\t0\ts__Mockella alpha\tspecies\t7\tGCF_1\n"
            "2\t6\t0\ts__Mockella beta\tspecies\t7\tGCF_2\n"
            "4\t4\t0\troot\tno rank\t0\t\n"
            "5\t4\t0\td__Bacteria\tdomain\t1\t\n"
            "6\t5\t0\tg__Mockella\tgenus\t6\t\n";

    constexpr char kPriorsHeader[] = "taxid\trep_genome\tmarkers\tduplicate_markers\tcheckm_completeness\tcheckm_contamination\tani_radius\t"
                                     "mean_intra_ani\tmin_intra_ani\tclustered_genomes\tgenome_size\tsized_genomes\trep_genome_size\t"
                                     "marker_bases\tmarker_share\n";
}

TEST(Composition, TheScannedReadsLineRoundTrips) {
    ScannedReads const scanned{ 1234, 2468, 370200 };
    std::string line = ScannedReadsLine(scanned);
    EXPECT_EQ(line, "@CO\tprotal scanned reads: fragments=1234 reads=2468 bases=370200\n");
    line.pop_back();
    auto const parsed = ParseScannedReadsLine(line);
    ASSERT_TRUE(parsed.has_value());
    EXPECT_EQ(parsed->fragments, 1234u);
    EXPECT_EQ(parsed->reads, 2468u);
    EXPECT_EQ(parsed->bases, 370200u);
    // Fields in any order, unknown ones ignored; a missing or broken one fails.
    EXPECT_TRUE(ParseScannedReadsLine(kSamScannedReadsComment + "bases=3 future=x reads=2 fragments=1").has_value());
    EXPECT_FALSE(ParseScannedReadsLine(kSamScannedReadsComment + "fragments=1 reads=2").has_value());
    EXPECT_FALSE(ParseScannedReadsLine(kSamScannedReadsComment + "fragments=1 reads=2 bases=x").has_value());
    EXPECT_FALSE(ParseScannedReadsLine("@CO\tsomething else").has_value());
}

TEST(Composition, TheSharesAsComputedByHand) {
    // A: depth 10 at 2 Mb, B: depth 5 at 4 Mb, C: depth 1 without a size (counted at the average size). 100 Mb read, 10%
    // of them the mates' overlap.
    std::vector<composition::Species> const species = { { 1, 10, 2e6 }, { 2, 5, 4e6 }, { 3, 1, -1 } };
    auto const r = composition::Compute(species, ScannedReads{ 400000, 800000, 100000000 }, 0.9, true);
    double const ags = (10 * 2e6 + 5 * 4e6) / 15;
    EXPECT_EQ(r.species, 3u);
    EXPECT_EQ(r.sized, 2u);
    EXPECT_NEAR(r.average_genome_size, ags, 1e-6);
    EXPECT_NEAR(r.genome_equivalents, 16, 1e-12);
    EXPECT_NEAR(r.scanned_fragment_bases, 90e6, 1e-6);
    double const attributed = 40e6 + ags;
    EXPECT_NEAR(r.attributed_bases, attributed, 1e-3);
    EXPECT_NEAR(r.attributed_share, attributed / 90e6, 1e-12);
    EXPECT_NEAR(r.attributed_fragments, attributed / 90e6 * 400000, 1e-6);
    double const unknown = (90e6 - attributed) / ags;
    EXPECT_NEAR(r.unknown_genome_equivalents, unknown, 1e-9);
    EXPECT_NEAR(r.unknown_share, unknown / (16 + unknown), 1e-12);
    EXPECT_EQ(r.median_depth, 5);
    EXPECT_EQ(r.lowest_depth, 1);
    EXPECT_NEAR(r.missing_species_at_median_depth, unknown / 5, 1e-9);
    EXPECT_NEAR(r.missing_species_at_lowest_depth, unknown, 1e-9);
    // Depth-weighted: two thirds of the cells have 2 Mb.
    EXPECT_EQ(r.genome_size_quantiles[0], 2e6);
    EXPECT_EQ(r.genome_size_quantiles[2], 2e6);
    EXPECT_EQ(r.genome_size_quantiles[3], 4e6);
    EXPECT_EQ(r.genome_size_quantiles[4], 4e6);
    EXPECT_NEAR(r.genome_size_sd, std::sqrt((10 * (2e6 - ags) * (2e6 - ags) + 5 * (4e6 - ags) * (4e6 - ags)) / 15), 1e-3);
    EXPECT_NEAR(r.FragmentsOf(10, 2e6), 10 * 2e6 / 90e6 * 400000, 1e-6);
    EXPECT_NEAR(r.FragmentsOf(1, -1), ags / 90e6 * 400000, 1e-6);  // without a size: the average one

    std::ostringstream table;
    composition::WriteHeader(table);
    composition::WriteRow(table, "S1", r);
    std::istringstream lines(table.str());
    std::string header, row;
    std::getline(lines, header);
    std::getline(lines, row);
    auto const columns = std::count(header.begin(), header.end(), '\t');
    EXPECT_EQ(std::count(row.begin(), row.end(), '\t'), columns);
    EXPECT_EQ(row.substr(0, row.find('\t', row.find('\t', 3) + 1)), "S1\t400000\t800000");
    EXPECT_NE(composition::Summary(r).find("species at the median depth"), std::string::npos);
}

TEST(Composition, OverExplainedNothingCalledOrNothingKnown) {
    // The called species explain more than was read (their depth or sizes overestimated): nothing unknown.
    auto const over = composition::Compute({ { 1, 10, 2e6 } }, ScannedReads{ 1000, 2000, 10000000 }, 1, true);
    EXPECT_NEAR(over.attributed_share, 2.0, 1e-12);
    EXPECT_EQ(over.unknown_share, 0);
    // Nothing called: every genome sequenced is unknown, whatever the size.
    auto const none = composition::Compute({}, ScannedReads{ 1000, 2000, 300000 }, 1, true);
    EXPECT_EQ(none.attributed_share, 0);
    EXPECT_EQ(none.unknown_share, 1);
    EXPECT_TRUE(std::isnan(none.average_genome_size));
    EXPECT_NE(composition::Summary(none).find("0.0% of 1000 fragments"), std::string::npos);
    // A database without sizes: no unknown share, even for a sample without calls.
    EXPECT_FALSE(composition::Compute({}, ScannedReads{ 1000, 2000, 300000 }, 1, false).HasUnknownShare());
    auto const unsized = composition::Compute({ { 1, 10, -1 } }, ScannedReads{ 1000, 2000, 300000 }, 1, false);
    EXPECT_FALSE(unsized.HasUnknownShare());
    EXPECT_TRUE(std::isnan(unsized.attributed_share));
    EXPECT_NE(composition::Summary(unsized).find("none with a genome size"), std::string::npos);
    // A SAM that does not say what was read: the sizes, but no shares.
    auto const unscanned = composition::Compute({ { 1, 10, 2e6 } }, std::nullopt, 1, true);
    EXPECT_FALSE(unscanned.HasUnknownShare());
    EXPECT_EQ(unscanned.average_genome_size, 2e6);
    EXPECT_TRUE(std::isnan(unscanned.attributed_share));
    EXPECT_TRUE(std::isnan(unscanned.FragmentsOf(10, 2e6)));
    EXPECT_NE(composition::Summary(unscanned).find("does not say how many reads were scanned"), std::string::npos);
    std::ostringstream table;
    composition::WriteRow(table, "S2", unscanned);
    EXPECT_EQ(table.str().substr(0, 12), "S2\tNA\tNA\tNA\t");
}

TEST(Composition, SpeciesPriorsGiveGenomeSizesByTheirHeader) {
    species_priors::Table table;
    std::istringstream is(std::string(kPriorsHeader) +
                          "1\tGCF_1\t120\t0\t99\t1\t95\t98\t96\t10\t4500000\t8\t4400000\t130000\t0.02955\n"
                          "2\tGCF_2\t120\t0\t-1\t-1\t-1\t-1\t-1\t1\t-1\t0\t-1\t129000\t-1\n"
                          "3\tGCF_3\t120\t0\t-1\t-1\t-1\t-1\t-1\t1\t0\t0\t0\t0\t-1\n");
    ASSERT_EQ(table.Read(is), "");
    EXPECT_EQ(table.Get(1).genome_size, 4500000);
    EXPECT_EQ(table.Get(1).rep_genome_size, 4400000);
    EXPECT_EQ(table.Get(1).marker_bases, 130000);
    EXPECT_EQ(table.Get(2).genome_size, species_priors::kUnknown);
    EXPECT_EQ(table.Get(3).genome_size, species_priors::kUnknown);  // 0 says nothing
    EXPECT_EQ(table.WithGenomeSize(), 1u);
    // The columns by name, wherever they are after the tenth; a table without them (before 2026-10-08) has no sizes.
    std::istringstream moved("taxid\trep_genome\tmarkers\tduplicate_markers\tcheckm_completeness\tcheckm_contamination\tani_radius\t"
                             "mean_intra_ani\tmin_intra_ani\tclustered_genomes\tother\tgenome_size\n"
                             "1\tGCF_1\t120\t0\t99\t1\t95\t98\t96\t10\tx\t3000000\n");
    species_priors::Table second;
    ASSERT_EQ(second.Read(moved), "");
    EXPECT_EQ(second.Get(1).genome_size, 3000000);
    EXPECT_EQ(second.Get(1).completeness, 99);
    std::istringstream old("taxid\trep_genome\tmarkers\tduplicate_markers\tcheckm_completeness\tcheckm_contamination\tani_radius\t"
                           "mean_intra_ani\tmin_intra_ani\tclustered_genomes\n1\tGCF_1\t120\t0\t99\t1\t95\t98\t96\t10\n");
    species_priors::Table third;
    ASSERT_EQ(third.Read(old), "");
    EXPECT_EQ(third.WithGenomeSize(), 0u);
    // A row shorter than the header's size columns, or a size that is no number, is an error.
    std::istringstream short_row(std::string(kPriorsHeader) + "1\tGCF_1\t120\t0\t99\t1\t95\t98\t96\t10\n");
    EXPECT_NE(species_priors::Table().Read(short_row), "");
    std::istringstream bad(std::string(kPriorsHeader) + "1\tGCF_1\t120\t0\t99\t1\t95\t98\t96\t10\tbig\t1\t1\t1\t1\n");
    EXPECT_NE(species_priors::Table().Read(bad), "");
}

// Two species with a gene each: a pair whose mates share 50 bases on species 1's gene, a single read on species 2's. The
// SAM header says 10 fragments of 3000 bases were scanned; species 1's genome has 4000 bp, species 2's 2000 bp.
TEST(Composition, TheProfileReadsTheScannedReadsAndEndsWithTheUnknownShare) {
    std::mt19937 rng(41);
    Reference ref({ { 1, { RandomSequence(400, rng) } }, { 2, { RandomSequence(400, rng) } } }, "composition reference");
    species_priors::Table priors;
    std::istringstream is(std::string(kPriorsHeader) + "1\tGCF_1\t1\t0\t-1\t-1\t-1\t-1\t-1\t1\t4000\t1\t4000\t400\t0.1\n"
                                                       "2\tGCF_2\t1\t0\t-1\t-1\t-1\t-1\t-1\t1\t2000\t1\t2000\t400\t0.2\n");
    ASSERT_EQ(priors.Read(is), "");
    ref.loader->SetSpeciesPriors(priors);
    taxonomy::IntTaxonomy taxonomy(ref.Write("taxonomy.dmp", kTaxonomy));
    std::string const sam = ref.Header() + ScannedReadsLine({ 10, 20, 3000 }) +
                            Record(ref, "p", 99, 1, 1, 10, 100) + Record(ref, "p", 147, 1, 1, 60, 100) + Record(ref, "s", 0, 2, 1, 10, 100);

    for (size_t threads : { 1, 3 }) {
        SCOPED_TRACE(threads);
        profiler::Profiler profiler(*ref.loader);
        profiler.SetDepthIdentityMargin(0.08);
        profiler.SetChunkBytes(200);
        profiler::MicrobialProfile profile(*ref.loader);
        auto const path = ref.Write("sample" + std::to_string(threads) + ".sam", sam);
        ASSERT_EQ(profiler.ProfileSam(path, profile, {}, 2, 2, 0.0, 15, 0, false, threads), "");
        ASSERT_TRUE(profiler.Scanned().has_value());
        EXPECT_EQ(profiler.Scanned()->bases, 3000u);
        EXPECT_NEAR(profile.FragmentBaseShare(), 250.0 / 300.0, 1e-12);

        auto filter = HalfModel(0.5);
        auto const species = profile.CalledSpecies(filter);
        ASSERT_EQ(species.size(), 2u);
        EXPECT_NEAR(species[0].depth, 150.0 / 400, 1e-12);
        EXPECT_NEAR(species[1].depth, 100.0 / 400, 1e-12);
        EXPECT_EQ(species[0].genome_size, 4000);
        auto result = composition::Compute(species, profiler.Scanned(), profile.FragmentBaseShare(), true);
        // 2500 fragment bases read, 1500 + 500 explained; 500 left at the average 3200 bp: 0.15625 genomes beside 0.625.
        EXPECT_NEAR(result.attributed_share, 0.8, 1e-12);
        EXPECT_NEAR(result.average_genome_size, 3200, 1e-9);
        EXPECT_NEAR(result.unknown_share, 0.2, 1e-12);
        profile.SetComposition(result, true);

        std::ostringstream os, os_total, os_genes;
        profile.WriteSparseProfile(taxonomy, filter, os, &os_total, &os_genes, threads);
        std::istringstream lines(os.str());
        std::vector<std::string> rows;
        for (std::string line; std::getline(lines, line);) rows.push_back(line);
        ASSERT_EQ(rows.size(), 3u);
        auto share = [](std::string const& row) { return std::stod(row.substr(row.rfind('\t') + 1)); };
        EXPECT_NEAR(share(rows[0]), 0.6 * 0.8, 1e-5);
        EXPECT_NEAR(share(rows[1]), 0.4 * 0.8, 1e-5);
        EXPECT_EQ(rows[2].substr(0, 4), "?\t?\t");
        EXPECT_NEAR(share(rows[2]), 0.2, 1e-5);
        // .profile.log: each taxon's genome size and the fragments its depth implies (species 1: 1500 of 2500 bases).
        std::istringstream log(os_total.str());
        std::string header, first;
        std::getline(log, header);
        std::getline(log, first);
        auto field = [](std::string const& text, size_t index) {
            size_t start = 0;
            for (size_t i = 0; i < index; i++) start = text.find('\t', start) + 1;
            return text.substr(start, text.find('\t', start) - start);
        };
        size_t column = 0;
        for (size_t start = 0; field(header, column) != "GenomeSize"; column++) ASSERT_LT(start++, 100u);
        EXPECT_EQ(field(first, column), "4000");
        EXPECT_EQ(field(first, column + 1), "6");  // 1500 / 2500 of 10 fragments

        // --no_unknown_share: the called species' shares of their depth, no "?" line.
        profile.SetComposition(result, false);
        std::ostringstream plain;
        profile.WriteSparseProfile(taxonomy, filter, plain, nullptr, nullptr, threads);
        EXPECT_EQ(plain.str().find('?'), std::string::npos);
        std::istringstream plain_lines(plain.str());
        std::string line;
        std::getline(plain_lines, line);
        EXPECT_NEAR(share(line), 0.6, 1e-5);
    }

    // A SAM without the line: nothing scanned known.
    profiler::Profiler profiler(*ref.loader);
    profiler::MicrobialProfile profile(*ref.loader);
    std::string const old = ref.Header() + Record(ref, "s", 0, 2, 1, 10, 100);
    ASSERT_EQ(profiler.ProfileSam(ref.Write("old.sam", old), profile), "");
    EXPECT_FALSE(profiler.Scanned().has_value());
    // A broken line stops the SAM.
    profiler::MicrobialProfile broken(*ref.loader);
    std::string const bad = ref.Header() + kSamScannedReadsComment + "fragments=1\n" + Record(ref, "s", 0, 2, 1, 10, 100);
    EXPECT_NE(profiler.ProfileSam(ref.Write("bad.sam", bad), broken), "");
}
