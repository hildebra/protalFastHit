// Unit tests for the species' NCBI taxids (NcbiTaxids.h): ncbi_taxids.tsv read by its header, its errors, and the profile
// of --ncbi_taxids, whose lines end with each species' NCBI taxid and name.
#include <gtest/gtest.h>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include "Profiling/Profiler.h"
#include "Profiling/Composition.h"
#include "IO/ScannedReads.h"
#include "Taxonomy/NcbiTaxids.h"
#include "SequenceUtils/SpeciesPriors.h"
#include "TestReference.h"

using namespace protal;
using protal::test::RandomSequence;
using Reference = protal::test::LoadedReference;

namespace {
    constexpr char kHeader[] = "taxid\trep_genome\tspecies\tncbi_taxid\tncbi_name\tnamed\tncbi_genomes\tgenomes\tshared_by\n";

    std::string Read(ncbi_taxids::Table& table, std::string const& text) {
        std::istringstream is(text);
        return table.Read(is);
    }

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

    // Fields of a tab-separated line.
    std::vector<std::string> Fields(std::string const& line) {
        std::vector<std::string> fields;
        std::istringstream is(line);
        for (std::string field; std::getline(is, field, '\t');) fields.push_back(field);
        return fields;
    }
}

TEST(NcbiTaxids, TheTableIsReadByItsHeader) {
    ncbi_taxids::Table table;
    ASSERT_EQ(Read(table, std::string(kHeader) +
                          "1\tGCF_1\tEscherichia coli\t562\tEscherichia coli\t1\t41484\t44640\t20\n"
                          "# a comment\n"
                          "2\tGCF_2\tPetralouisia sp910576975\t297314\tuncultured Lachnospiraceae bacterium\t0\t136\t136\t1048\n"), "");
    EXPECT_EQ(table.Size(), 2u);
    EXPECT_EQ(table.Named(), 1u);
    ASSERT_NE(table.Get(1), nullptr);
    EXPECT_EQ(table.Get(1)->ncbi_taxid, 562u);
    EXPECT_EQ(table.Get(1)->name, "Escherichia coli");
    EXPECT_TRUE(table.Get(1)->named);
    EXPECT_EQ(table.Get(2)->name, "uncultured Lachnospiraceae bacterium");
    EXPECT_FALSE(table.Get(2)->named);
    EXPECT_EQ(table.Get(3), nullptr);
    // The columns by name, in any order, the others skipped.
    ncbi_taxids::Table moved;
    ASSERT_EQ(Read(moved, "taxid\tnamed\tother\tncbi_name\tncbi_taxid\n7\t1\tx\tMockella alpha\t1000001\r\n"), "");
    EXPECT_EQ(moved.Get(7)->ncbi_taxid, 1000001u);
    EXPECT_EQ(moved.Get(7)->name, "Mockella alpha");
}

TEST(NcbiTaxids, BrokenTablesAreRefused) {
    std::string const row = "1\tGCF_1\tEscherichia coli\t562\tEscherichia coli\t1\t1\t1\t1\n";
    for (std::string const& text : {
            row,                                                                   // no header
            std::string("taxid\trep_genome\tncbi_taxid\tncbi_name\n") + row,       // no named column
            std::string(kHeader) + "1\tGCF_1\tEscherichia coli\t562\n",            // a short line
            std::string(kHeader) + "1\tGCF_1\tEscherichia coli\tx\tEscherichia coli\t1\t1\t1\t1\n",   // no number
            std::string(kHeader) + "1\tGCF_1\tEscherichia coli\t0\tEscherichia coli\t1\t1\t1\t1\n",   // taxid 0
            std::string(kHeader) + "1\tGCF_1\tEscherichia coli\t562\tEscherichia coli\tyes\t1\t1\t1\n",
            std::string(kHeader) + row + row,                                      // a species twice
            std::string("") }) {
        SCOPED_TRACE(text);
        ncbi_taxids::Table table;
        EXPECT_NE(Read(table, text), "");
    }
}

// Two species with a gene each, both called; the table knows the first only. The profile's lines end with the NCBI taxid
// and name (NA, NA for the second), the "?" line with ?, ?; without the table they are as before.
TEST(NcbiTaxids, TheProfileEndsWithTheNcbiColumns) {
    std::mt19937 rng(43);
    Reference ref({ { 1, { RandomSequence(400, rng) } }, { 2, { RandomSequence(400, rng) } } }, "ncbi reference");
    species_priors::Table priors;
    std::istringstream is("taxid\trep_genome\tmarkers\tduplicate_markers\tcheckm_completeness\tcheckm_contamination\tani_radius\t"
                          "mean_intra_ani\tmin_intra_ani\tclustered_genomes\tgenome_size\n"
                          "1\tGCF_1\t1\t0\t-1\t-1\t-1\t-1\t-1\t1\t4000\n2\tGCF_2\t1\t0\t-1\t-1\t-1\t-1\t-1\t1\t2000\n");
    ASSERT_EQ(priors.Read(is), "");
    ref.loader->SetSpeciesPriors(priors);
    taxonomy::IntTaxonomy taxonomy(ref.Write("taxonomy.dmp", kTaxonomy));
    std::string const sam = ref.Header() + ScannedReadsLine({ 10, 20, 3000 }) +
                            Record(ref, "p", 99, 1, 1, 10, 100) + Record(ref, "p", 147, 1, 1, 60, 100) + Record(ref, "s", 0, 2, 1, 10, 100);
    ncbi_taxids::Table ncbi;
    ASSERT_EQ(Read(ncbi, std::string(kHeader) + "1\tGCF_1\tMockella alpha\t1000001\tMockella alpha\t1\t2\t3\t1\n"), "");

    profiler::Profiler profiler(*ref.loader);
    profiler.SetDepthIdentityMargin(0.08);
    profiler::MicrobialProfile profile(*ref.loader);
    ASSERT_EQ(profiler.ProfileSam(ref.Write("sample.sam", sam), profile, {}, 2, 2, 0.0, 15, 0, false, 1), "");
    auto filter = HalfModel(0.5);
    profile.SetComposition(composition::Compute(profile.CalledSpecies(filter), profiler.Scanned(), profile.FragmentBaseShare(), true), true);

    std::ostringstream plain, with;
    profile.WriteSparseProfile(taxonomy, filter, plain, nullptr, nullptr, 1);
    profile.WriteSparseProfile(taxonomy, filter, with, nullptr, nullptr, 1, &ncbi);
    std::istringstream plain_lines(plain.str()), with_lines(with.str());
    std::vector<std::vector<std::string>> rows;
    for (std::string a, b; std::getline(plain_lines, a) && std::getline(with_lines, b);) {
        auto const fields = Fields(b);
        ASSERT_EQ(fields.size(), 5u) << b;
        EXPECT_EQ(b.substr(0, a.size() + 1), a + '\t') << "the columns before are the same";
        rows.push_back(fields);
    }
    ASSERT_EQ(rows.size(), 3u);
    EXPECT_EQ(rows[0][0], "GCF_1");
    EXPECT_EQ(std::vector<std::string>(rows[0].begin() + 3, rows[0].end()), (std::vector<std::string>{ "1000001", "Mockella alpha" }));
    EXPECT_EQ(rows[1][0], "GCF_2");
    EXPECT_EQ(std::vector<std::string>(rows[1].begin() + 3, rows[1].end()), (std::vector<std::string>{ "NA", "NA" }));
    EXPECT_EQ(rows[2][0], "?");
    EXPECT_EQ(std::vector<std::string>(rows[2].begin() + 3, rows[2].end()), (std::vector<std::string>{ "?", "?" }));
}
