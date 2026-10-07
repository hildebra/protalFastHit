// Unit tests for the model's inputs and contract: the features of a taxon (their names, values and precision), the
// check that a model fits protal, and the knobs by the sample's depth that a model's header carries.
#include <gtest/gtest.h>
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include "Profiling/Profiler.h"
#include "TestReference.h"

using namespace protal;

TEST(ModelFeatures, NormalizedFeaturesOfATaxon) {
    // Genome 1 has two hittable genes of 10 and 8 bp. Three mates of two fragments land on gene 1_1,
    // one with an insertion.
    test::LoadedReference ref({ { 1, { "ACGTACGTAA", "ccggttaa" } }, { 2, { "GGGGCCCCAT" } } }, "model features");
    profiler::MicrobialProfile profile(*ref.loader);
    auto sam = [](std::string seq, std::string cigar) {
        SamEntry s;
        s.m_qname = "r";
        s.m_rname = "1_1";
        s.m_pos = 1;
        s.m_mapq = 60;
        s.m_qual = std::string(seq.size(), 'I');
        s.m_seq = std::move(seq);
        s.m_cigar = std::move(cigar);
        return s;
    };
    auto exact = sam("ACGTACGTAA", "10M");
    auto insertion = sam("ACGTAGCGTAA", "5M1I5M");  // identity 10/11
    ASSERT_TRUE(profile.AddSam(1, 1, exact, 1.0, true, 0));
    ASSERT_TRUE(profile.AddSam(1, 1, exact, 1.0, true, 0));  // its mate: the same fragment
    ASSERT_TRUE(profile.AddSam(1, 1, insertion, 1.0, true, 1));
    auto const& taxon = profile.GetTaxa().at(1);

    EXPECT_EQ(taxon.TotalHits(), 3u);
    EXPECT_EQ(taxon.Fragments(), 2u);
    EXPECT_NEAR(taxon.BaseIdentity(), (10 + 10 + 10 * 10.0 / 11) / 30, 1e-9);
    EXPECT_NEAR(taxon.TopIdentity(), 1.0, 1e-6);
    EXPECT_NEAR(taxon.HitGeneFraction(), 0.5, 1e-12);
    // Two fragments over genes of 10 and 8 bp: 1 - (8/18)^2 + 1 - (10/18)^2 genes expected, 1 hit.
    double const expected_genes = 2 - std::pow(8.0 / 18, 2) - std::pow(10.0 / 18, 2);
    EXPECT_NEAR(taxon.GenePresenceRatio(), 1 / expected_genes, 1e-9);
    // Fragments per gene 2 and 0 against 2 * 10/18 and 2 * 8/18.
    double const e1 = 2 * 10.0 / 18, e2 = 2 * 8.0 / 18;
    EXPECT_NEAR(taxon.GeneDispersion(), (2 - e1) * (2 - e1) / e1 + e2, 1e-9);
}

namespace {
    // A one-split tree on `field`, predicting `yes` or `no`, as PMML.
    std::string OneSplitModel(std::string const& field, std::string const& yes = "TRUE", std::string const& no = "FALSE") {
        return R"(<?xml version="1.0" encoding="UTF-8"?>
<PMML version="4.4">
 <Header/>
 <DataDictionary>
  <DataField name="truth" optype="categorical" dataType="string"><Value value=")" + no + R"("/><Value value=")" + yes + R"("/></DataField>
  <DataField name=")" + field + R"(" optype="continuous" dataType="double"/>
 </DataDictionary>
 <TreeModel functionName="classification" splitCharacteristic="binarySplit">
  <MiningSchema>
   <MiningField name="truth" usageType="predicted"/>
   <MiningField name=")" + field + R"("/>
  </MiningSchema>
  <Node score=")" + no + R"(">
   <True/>
   <Node score=")" + yes + R"("><SimplePredicate field=")" + field + R"(" operator="greaterThan" value="10"/><ScoreDistribution value=")" + no + R"(" recordCount="1"/><ScoreDistribution value=")" + yes + R"(" recordCount="9"/></Node>
   <Node score=")" + no + R"("><True/><ScoreDistribution value=")" + no + R"(" recordCount="9"/><ScoreDistribution value=")" + yes + R"(" recordCount="1"/></Node>
  </Node>
 </TreeModel>
</PMML>
)";
    }
}

TEST(ModelFeatures, TheModelMustFitProtal) {
    // As protal checks the model it loads (LoadModel): parsed from its text, which ModelContractProblemInXml reads too.
    auto problem = [](std::string const& xml) {
        return profiler::ModelContractProblemInXml(profiler::TaxonFilterForest(cpmml::Model::from_string(xml), 0.5), xml);
    };
    EXPECT_EQ(problem(OneSplitModel("present_genes")), "");
    EXPECT_EQ(problem(OneSplitModel("gene_presence_ratio")), "");
    EXPECT_EQ(problem(OneSplitModel("coverage_of_moon")), "it needs 1 input(s) protal does not compute: coverage_of_moon");
    EXPECT_NE(problem(OneSplitModel("present_genes", "yes", "no")).find("has no value TRUE"), std::string::npos);
}

namespace {
    // A model of one leaf on `fragments`, with `extension` in its header, where the trainer writes the depth knobs (cPMML
    // loads past it), and `after_header` after it.
    std::string KnobModel(std::string const& extension, std::string const& after_header = "") {
        return "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<PMML version=\"4.4\">\n <Header description=\"m\">\n  " + extension +
               "\n  <Application name=\"test\"/>\n </Header>\n" + after_header + R"( <DataDictionary>
  <DataField name="truth" optype="categorical" dataType="string"><Value value="FALSE"/><Value value="TRUE"/></DataField>
  <DataField name="fragments" optype="continuous" dataType="double"/>
 </DataDictionary>
 <TreeModel functionName="classification" splitCharacteristic="binarySplit">
  <MiningSchema>
   <MiningField name="truth" usageType="predicted"/>
   <MiningField name="fragments"/>
  </MiningSchema>
  <Node score="FALSE"><True/><ScoreDistribution value="FALSE" recordCount="0.4"/><ScoreDistribution value="TRUE" recordCount="0.6"/></Node>
 </TreeModel>
</PMML>
)";
    }

    std::string KnobBins(std::string const& value) { return R"(<Extension name="protal_depth_knobs" value=")" + value + R"("/>)"; }
    std::string KnobCurve(std::string const& value) { return R"(<Extension name="protal_depth_knob_curve" value=")" + value + R"("/>)"; }
}

TEST(ModelFeatures, TheDepthKnobBinIsTheDigitsOfTheSamplesFragmentsLessOne) {
    // 2 to 6, as machine_learning_cmdline.py's depth_bins.
    EXPECT_EQ(profiler::DepthKnobBin(0), 2);
    EXPECT_EQ(profiler::DepthKnobBin(999), 2);
    EXPECT_EQ(profiler::DepthKnobBin(1000), 3);
    EXPECT_EQ(profiler::DepthKnobBin(99999), 4);
    EXPECT_EQ(profiler::DepthKnobBin(1000000), 6);
    EXPECT_EQ(profiler::DepthKnobBin(5000000000ull), 6);
}

TEST(ModelFeatures, DepthKnobBinsComeFromTheirExtensionInTheHeader) {
    std::map<int, double> knobs;
    EXPECT_EQ(profiler::ParseDepthKnobs(KnobModel(KnobBins("2:0.31,4:0.4")), knobs), "");
    EXPECT_EQ(knobs, (std::map<int, double>{ { 2, 0.31 }, { 4, 0.4 } }));
    // None: another extension, one after the header, or a knob curve.
    EXPECT_EQ(profiler::ParseDepthKnobs(KnobModel(R"(<Extension name="other" value="2:0.3"/>)"), knobs), "");
    EXPECT_TRUE(knobs.empty());
    EXPECT_EQ(profiler::ParseDepthKnobs(KnobModel("", KnobBins("2:0.3")), knobs), "");
    EXPECT_TRUE(knobs.empty());
    EXPECT_EQ(profiler::ParseDepthKnobs(KnobModel(KnobCurve("2.000:0.2,4.000:0.8")), knobs), "");
    EXPECT_TRUE(knobs.empty());
}

TEST(ModelFeatures, MalformedDepthKnobBinsAreAProblem) {
    // A bin outside 2-6, a knob outside 0-1, a bin twice, no colon, trailing characters, no bins: no knobs.
    for (std::string const value : { "7:0.3", "2:1.5", "2:0.3,2:0.4", "2-0.3", "2:0.3x", "" }) {
        std::map<int, double> knobs = { { 3, 0.5 } };
        auto const problem = profiler::ParseDepthKnobs(KnobModel(KnobBins(value)), knobs);
        EXPECT_NE(problem.find("its depth knobs are malformed"), std::string::npos) << value << ": " << problem;
        EXPECT_TRUE(knobs.empty()) << value;
    }
}

TEST(ModelFeatures, AModelsKnobComesFromTheBinOfTheSamplesFragments) {
    std::string const xml = KnobModel(KnobBins("2:0.31,4:0.4"));
    std::map<int, double> knobs;
    ASSERT_EQ(profiler::ParseDepthKnobs(xml, knobs), "");
    profiler::TaxonFilterForest model(cpmml::Model::from_string(xml), 0.5);
    EXPECT_EQ(profiler::ModelContractProblemInXml(model, xml), "");
    EXPECT_FALSE(model.DepthKnob(500).has_value());  // none set yet
    model.SetDepthKnobs(knobs);
    EXPECT_EQ(model.DepthKnob(500), 0.31);
    EXPECT_FALSE(model.DepthKnob(5000).has_value());  // bin 3: --knob's default
    EXPECT_EQ(model.DepthKnob(50000), 0.4);
    auto const copy = model.WithKnob(0.2);
    EXPECT_EQ(copy.GetKnob(), 0.2);
    EXPECT_EQ(copy.DepthKnobs(), knobs);
}

TEST(ModelFeatures, ADepthKnobCurveIsLinearInLog10OfTheSamplesFragments) {
    // The trainer's since 0.7.3: log10 of the sample's fragments : knob, linear between the points, the ends' beyond them
    // (machine_learning_cmdline.py knob_at).
    std::string const xml = KnobModel(KnobCurve("2.000:0.2,4.000:0.8"));
    profiler::DepthKnobCurve curve;
    EXPECT_EQ(profiler::ParseDepthKnobCurve(xml, curve), "");
    EXPECT_EQ(curve, (profiler::DepthKnobCurve{ { 2.0, 0.2 }, { 4.0, 0.8 } }));
    profiler::TaxonFilterForest curved(cpmml::Model::from_string(xml), 0.5);
    EXPECT_FALSE(curved.HasDepthKnobs());
    curved.SetDepthKnobCurve(curve);
    EXPECT_TRUE(curved.HasDepthKnobs());
    EXPECT_DOUBLE_EQ(*curved.DepthKnob(0), 0.2);        // below the first point: its knob
    EXPECT_DOUBLE_EQ(*curved.DepthKnob(100), 0.2);
    EXPECT_DOUBLE_EQ(*curved.DepthKnob(1000), 0.5);     // halfway in log10
    EXPECT_NEAR(*curved.DepthKnob(3162), 0.65, 1e-4);   // 10^3.5
    EXPECT_DOUBLE_EQ(*curved.DepthKnob(10000), 0.8);
    EXPECT_DOUBLE_EQ(*curved.DepthKnob(5000000000ull), 0.8);  // deeper than any trained: the deepest knob
    EXPECT_EQ(curved.WithKnob(0.3).GetDepthKnobCurve(), curve);
    EXPECT_DOUBLE_EQ(profiler::DepthKnobAt({ { 1.5, 0.4 } }, 7), 0.4);  // one point: its knob everywhere
    EXPECT_DOUBLE_EQ(profiler::DepthKnobAt({}, 7), 0.5);
    EXPECT_EQ(profiler::ParseDepthKnobCurve(KnobModel(""), curve), "");  // none
    EXPECT_TRUE(curve.empty());
}

TEST(ModelFeatures, MalformedDepthKnobCurvesAreAProblem) {
    // x not increasing or outside 0-12, a knob outside 0-1, no colon, trailing characters, no points: no curve.
    for (std::string const value : { "3:0.3,2:0.4", "2:0.3,2:0.4", "13:0.3", "-1:0.3", "2:1.5", "2-0.3", "2:0.3x", "" }) {
        profiler::DepthKnobCurve curve = { { 1.0, 0.5 } };
        auto const problem = profiler::ParseDepthKnobCurve(KnobModel(KnobCurve(value)), curve);
        EXPECT_NE(problem.find("its depth knob curve is malformed"), std::string::npos) << value << ": " << problem;
        EXPECT_TRUE(curve.empty()) << value;
    }
}

TEST(ModelFeatures, NamesAreUniqueAndValuesKeepTheirPrecision) {
    Genome no_genome(0);
    auto features = profiler::TaxonFeatures(profiler::Taxon(no_genome));
    std::set<std::string> names;
    for (auto const& [name, _] : features) EXPECT_TRUE(names.insert(name).second) << name << " twice";
    for (auto const* name : { "RAF0", "RA4", "su_rate_ref", "lu_rate_ref", "lsu_rate_ref", "mean_mapq", "lu_per_read",
                              "fragments", "gene_presence_ratio", "identity", "variant_sites_per_kb" }) {
        EXPECT_TRUE(names.contains(name)) << name;
    }

    EXPECT_EQ(profiler::FeatureString(0.5), "0.5");
    EXPECT_EQ(profiler::FeatureString(3), "3");
    EXPECT_EQ(profiler::FeatureString(-1), "-1");
    EXPECT_EQ(std::stod(profiler::FeatureString(1.25e-7)), 1.25e-7);
    EXPECT_EQ(std::stod(profiler::FeatureString(0.1 + 0.2)), 0.1 + 0.2);
    // Subnormal values are written as 0: cPMML reads the values back with stod, which throws on underflow (an EM share of
    // 1e-311 made a run crash). The smallest normal value is kept.
    EXPECT_EQ(profiler::FeatureString(9.0946131981813e-311), "0");
    EXPECT_EQ(profiler::FeatureString(-4e-320), "0");
    EXPECT_EQ(profiler::FeatureString(2.2250738585072014e-308), "2.2250738585072014e-308");
}
