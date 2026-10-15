// Unit tests for what ties a run to the database and models it is profiled with (since 2026-10-15): the run settings
// build_settings.tsv records and runs follow (BuildSettings.h RunSettings, Options' auto), the SAM header's settings line
// (SamSettings.h), the versions of the features' definitions (FeatureVersions.h) and the features' training ranges
// (FeatureRanges.h).
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include "Options.h"
#include "Profiling/Profiler.h"
#include "FeatureRanges.h"
#include "FeatureVersions.h"
#include "SamSettings.h"

using namespace protal;

namespace {
    size_t index(ReadType type) { return static_cast<size_t>(type); }

    build_settings::RunSettings EveryRunSetting() {
        build_settings::RunSettings s;
        s.weighted_site_shift = true;
        s.allele_scores = false;
        s.align_top = 5;
        s.adaptive_candidates = 0;
        s.mate_guidance = "off";
        s.max_score_ani = { 0.9, 0.9, 0.9, 0.85 };
        s.unique_kmers = "cores";
        s.keep_foreign_genes = true;
        s.keep_suspect_copies = false;
        s.depth_identity_margin = 0.08;
        return s;
    }
}

// Every run setting is written as its name and value and read back; the decimals as their shortest text.
TEST(RunSettings, WrittenAndReadBack) {
    build_settings::Settings written;
    written.wfa_traceback = WfaTraceback::Singletrack;
    written.run = EveryRunSetting();
    std::stringstream file;
    build_settings::Write(file, written);
    std::string const text = file.str();
    for (auto const* line : { "weighted_site_shift\ton\n", "allele_scores\toff\n", "align_top\t5\n", "adaptive_candidates\t0\n",
                              "mate_guidance\toff\n", "max_score_ani_pe\t0.9\n", "max_score_ani_ont\t0.85\n",
                              "unique_kmers\tcores\n", "keep_foreign_genes\ton\n", "keep_suspect_copies\toff\n",
                              "depth_identity_margin\t0.08\n" }) {
        EXPECT_NE(text.find(line), std::string::npos) << line;
    }
    build_settings::Settings read;
    ASSERT_EQ(build_settings::Read(file, read), "");
    EXPECT_EQ(read.run, written.run);
    EXPECT_EQ(read.wfa_traceback, WfaTraceback::Singletrack);
    EXPECT_EQ(build_settings::RunFields(read.run).size(), build_settings::RunSettingNames().size());
    for (size_t i = 0; i < build_settings::RunSettingNames().size(); i++) {
        EXPECT_EQ(build_settings::RunFields(read.run)[i].first, build_settings::RunSettingNames()[i]);
    }
}

// A file of an older protal records none (each nothing); one of a later protal may hold settings and read types this one
// does not know (skipped) and modes it does not know (kept as text); a value that cannot be is refused, naming its line.
// The new lines are none of the older keys, so an older protal skips them as unknown (no sparse_first_ prefix).
TEST(RunSettings, OlderFilesSayNothingLaterOnesAreTolerated) {
    build_settings::Settings older;
    std::istringstream before("wfa_traceback\tsingletrack\nfragment_count\tmapq\nsparse_first_pe\t4\n");
    ASSERT_EQ(build_settings::Read(before, older), "");
    EXPECT_EQ(older.run, build_settings::RunSettings());
    EXPECT_TRUE(build_settings::RunFields(older.run).empty());

    build_settings::Settings later;
    std::istringstream after("max_score_ani_hifi\t0.95\nmate_guidance\tvote\nunique_kmers\tgrid\nnew_setting\t1\n"
                             "allele_scores\ton\n");
    ASSERT_EQ(build_settings::Read(after, later), "");
    EXPECT_EQ(later.run.mate_guidance, std::optional<std::string>("vote"));
    EXPECT_FALSE(build_settings::KnownMateGuidance(*later.run.mate_guidance));
    for (auto const* mode : { "consensus", "guide", "off" }) EXPECT_TRUE(build_settings::KnownMateGuidance(mode)) << mode;
    EXPECT_FALSE(build_settings::KnownMateGuidance("on"));
    EXPECT_EQ(later.run.unique_kmers, std::optional<std::string>("grid"));
    EXPECT_EQ(later.run.allele_scores, std::optional<bool>(true));
    for (auto const& floor : later.run.max_score_ani) EXPECT_FALSE(floor);

    for (std::string const name : build_settings::RunSettingNames()) EXPECT_NE(name.rfind("sparse_first_", 0), 0u) << name;

    for (auto const* bad : { "allele_scores\tyes\n", "align_top\tthree\n", "adaptive_candidates\t-1\n",
                             "depth_identity_margin\t-0.1\n", "max_score_ani_ont\tnan\n", "keep_foreign_genes\t1\n",
                             "mate_guidance\t\n" }) {
        build_settings::Settings settings;
        std::istringstream is(std::string("# a comment\n") + bad);
        EXPECT_NE(build_settings::Read(is, settings).find("line 2"), std::string::npos) << bad;
    }
}

TEST(RunSettings, DecimalsReadBackExactly) {
    for (double value : { 0.9, 0.85, 0.08, 1.0, 0.0, 0.123456789 }) {
        auto const text = build_settings::Decimal(value);
        ASSERT_TRUE(build_settings::ParseDecimal(text)) << text;
        EXPECT_EQ(*build_settings::ParseDecimal(text), value) << text;
    }
    EXPECT_EQ(build_settings::Decimal(0.9), "0.9");
    EXPECT_EQ(build_settings::Decimal(0.85), "0.85");
    EXPECT_FALSE(build_settings::ParseDecimal("0.9x"));
    EXPECT_FALSE(build_settings::ParseDecimal(""));
    EXPECT_FALSE(build_settings::ParseDecimal("inf"));
}

// A run is compared with the settings the database records only, and only for the settings the run is about.
TEST(RunSettings, DifferencesOfTheRecordedSettingsOnly) {
    build_settings::RunSettings recorded;
    recorded.allele_scores = true;
    recorded.align_top = 3;
    recorded.max_score_ani[index(ReadType::ONT)] = 0.85;
    auto used = EveryRunSetting();  // allele_scores off, align_top 5, ONT 0.85, ...
    auto const all = [](std::string const&) { return true; };
    auto differences = build_settings::Differences(build_settings::RunFields(recorded), build_settings::RunFields(used), all);
    ASSERT_EQ(differences.size(), 2u);
    EXPECT_EQ(differences[0], "allele_scores off (the database's: on)");
    EXPECT_EQ(differences[1], "align_top 5 (the database's: 3)");
    auto const no_short = [](std::string const& name) { return name != build_settings::kAlleleScores; };
    differences = build_settings::Differences(build_settings::RunFields(recorded), build_settings::RunFields(used), no_short);
    ASSERT_EQ(differences.size(), 1u);
    used.max_score_ani[index(ReadType::ONT)] = 0.8;
    differences = build_settings::Differences(build_settings::RunFields(recorded), build_settings::RunFields(used), no_short);
    ASSERT_EQ(differences.size(), 2u);
    EXPECT_EQ(differences[1], "max_score_ani_ont 0.8 (the database's: 0.85)");
}

// A run's options left at auto take the database's settings; those given stay; -a given applies to every read type, and
// without it each read type takes the database's floor, else its own default (ONT 0.85).
TEST(RunSettings, AutoTakesTheDatabasesGivenStays) {
    OptionsData d;
    d.run_settings_given.align_top = 7;
    d.align_top = 7;
    d.run_settings_given.keep_foreign_genes = false;
    Options options(d);
    EXPECT_EQ(options.GetMaxScoreAni(ReadType::ONT), 0.85);
    EXPECT_EQ(options.GetMaxScoreAni(ReadType::Paired), DEFAULT_MAX_SCORE_ANI);
    auto recorded = EveryRunSetting();
    recorded.max_score_ani = { 0.92, std::nullopt, 0.88, 0.8 };
    options.FollowDatabaseRunSettings(recorded);
    EXPECT_EQ(options.GetAlignTop(), 7u);                   // given
    EXPECT_TRUE(options.DropForeignGenes());                // given off: left out, though the database keeps them
    EXPECT_TRUE(options.WeightedSiteShift());               // the database's
    EXPECT_TRUE(options.NoAlleleScores());
    EXPECT_EQ(options.GetAdaptiveCandidates(), 0u);
    EXPECT_EQ(options.GetMateGuidance(), MateGuidanceMode::Off);
    EXPECT_EQ(options.GetMaxScoreAni(ReadType::Paired), 0.92);
    EXPECT_EQ(options.GetMaxScoreAni(ReadType::Single), DEFAULT_MAX_SCORE_ANI);  // not recorded: protal's default
    EXPECT_EQ(options.GetMaxScoreAni(ReadType::ONT), 0.8);
    auto const used = options.RunSettingsInUse();
    EXPECT_EQ(used.align_top, std::optional<size_t>(7));
    EXPECT_EQ(used.mate_guidance, std::optional<std::string>("off"));
    EXPECT_EQ(used.max_score_ani[index(ReadType::PacBio)], std::optional<double>(0.88));
    EXPECT_EQ(used.unique_kmers, std::optional<std::string>("cores"));

    // A mode this protal does not know leaves its default.
    Options plain{ OptionsData() };
    recorded.mate_guidance = "vote";
    plain.FollowDatabaseRunSettings(recorded);
    EXPECT_EQ(plain.GetMateGuidance(), kDefaultMateGuidance);
    // A database that records no mode (built before 2026-10-15) trained its models with guide; a mode given stays.
    Options older{ OptionsData() };
    older.FollowDatabaseRunSettings(build_settings::RunSettings());
    EXPECT_EQ(older.GetMateGuidance(), MateGuidanceMode::Guide);
    EXPECT_EQ(older.RunSettingsInUse().mate_guidance, std::optional<std::string>("guide"));
    OptionsData consensus_given;
    consensus_given.run_settings_given.mate_guidance = "consensus";
    consensus_given.mate_guidance = MateGuidanceMode::Consensus;
    Options given{ consensus_given };
    given.FollowDatabaseRunSettings(build_settings::RunSettings());
    EXPECT_EQ(given.GetMateGuidance(), MateGuidanceMode::Consensus);
    recorded.mate_guidance = "guide";
    given.FollowDatabaseRunSettings(recorded);
    EXPECT_EQ(given.GetMateGuidance(), MateGuidanceMode::Consensus);

    // -a given: every read type's, whatever the database records.
    OptionsData a;
    a.max_score_ani = 0.95;
    a.max_score_ani_given = true;
    a.run_settings_given.max_score_ani.fill(0.95);
    Options floors(a);
    floors.FollowDatabaseRunSettings(recorded);
    for (auto const& info : kReadTypes) EXPECT_EQ(floors.GetMaxScoreAni(info.type), 0.95);
}

// --profile_only: a SAM's gate becomes its sample's, unless --sparse_first of its read type gives one.
TEST(RunSettings, ASamsGateIsItsSamples) {
    OptionsData d;
    d.prefix_list = { "a", "b" };
    Options options(d);
    options.SetSampleGate(1, 4, SparseFirstPick::Grid);
    EXPECT_EQ(options.GetSparseFirst(size_t{0}), 0u);
    EXPECT_EQ(options.GetSparseFirst(size_t{1}), 4u);
    EXPECT_EQ(options.GetSparseFirstPick(size_t{1}), SparseFirstPick::Grid);
    options.SetSampleGate(5, 8, std::nullopt);  // no such sample
    OptionsData given = d;
    given.sparse_first[index(ReadType::Paired)] = 2;
    Options told(given);
    told.SetSampleGate(0, 4, std::nullopt);
    EXPECT_EQ(told.GetSparseFirst(size_t{0}), 2u);
}

// The SAM header's settings line: written as key=value, read back with the keys of a later protal kept; a header without
// it gives nothing, and only the header is read.
TEST(SamSettings, LineWrittenAndReadBack) {
    sam_settings::Fields const fields = { { "read_type", "pe" }, { "gate", "4" }, { "pick", "stride" }, { "unique_kmers", "cores" },
                                          { "traceback", "singletrack" }, { "allele_scores", "on" }, { "max_score_ani", "0.9" } };
    std::string const line = sam_settings::SamSettingsLine(fields);
    EXPECT_EQ(line, "@CO\tprotal settings: read_type=pe gate=4 pick=stride unique_kmers=cores traceback=singletrack "
                    "allele_scores=on max_score_ani=0.9\n");
    auto const parsed = sam_settings::ParseSamSettingsLine(line.substr(0, line.size() - 1));
    ASSERT_TRUE(parsed);
    EXPECT_EQ(*parsed, fields);
    EXPECT_EQ(sam_settings::Get(*parsed, "gate"), std::optional<std::string>("4"));
    EXPECT_FALSE(sam_settings::Get(*parsed, "mate_guidance"));

    auto const later = sam_settings::ParseSamSettingsLine("@CO\tprotal settings: gate=2 novel=x junk =y gate=8 empty=\r");
    ASSERT_TRUE(later);
    EXPECT_EQ(sam_settings::Get(*later, "gate"), std::optional<std::string>("8"));
    EXPECT_EQ(sam_settings::Get(*later, "novel"), std::optional<std::string>("x"));
    EXPECT_EQ(sam_settings::Get(*later, "empty"), std::optional<std::string>(""));
    EXPECT_EQ(later->size(), 3u);
    EXPECT_FALSE(sam_settings::ParseSamSettingsLine("@CO\tprotal read type: pe"));

    std::istringstream sam("@HD\tVN:1.6\n@SQ\tSN:1_1\tLN:100\n@CO\tprotal read type: pe\n" + line +
                           "r1\t0\t1_1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\n");
    auto const read = sam_settings::ReadSamSettings(sam);
    ASSERT_TRUE(read);
    EXPECT_EQ(*read, fields);
    std::istringstream older("@HD\tVN:1.6\n@CO\tprotal read type: pe\nr1\t0\t1_1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\n"
                             "@CO\tprotal settings: gate=4\n");
    EXPECT_FALSE(sam_settings::ReadSamSettings(older));
}

// The table holds every feature the model may be given (TaxonFeatures) once, and nothing else, each with a version, the
// date of its definition and a note; --feature_versions prints it.
TEST(FeatureVersions, TheTableHoldsEveryFeatureOnce) {
    Genome no_genome(0);
    std::vector<std::string> features;
    for (auto const& [name, _] : profiler::TaxonFeatures(profiler::Taxon(no_genome))) features.push_back(name);
    std::set<std::string> listed;
    size_t i = 0;
    for (auto const& f : feature_versions::kFeatureVersions) {
        EXPECT_TRUE(listed.insert(std::string(f.name)).second) << f.name << " twice";
        EXPECT_GE(f.version, 1) << f.name;
        EXPECT_EQ(f.changed.size(), 10u) << f.name;
        EXPECT_TRUE(f.changed.substr(0, 3) == "202" && f.changed[4] == '-' && f.changed[7] == '-') << f.name;
        EXPECT_FALSE(f.note.empty()) << f.name;
        EXPECT_EQ(f.note.find('\t'), std::string_view::npos) << f.name;
        if (i < features.size()) EXPECT_EQ(f.name, features[i]) << "the table's order is TaxonFeatures'";
        i++;
    }
    EXPECT_EQ(i, features.size());
    for (auto const& name : features) EXPECT_TRUE(listed.contains(name)) << name << " has no version";

    std::ostringstream os;
    feature_versions::WriteTable(os);
    std::istringstream lines(os.str());
    std::string line;
    std::getline(lines, line);
    EXPECT_EQ(line, "feature\tversion\tchanged\tnote");
    size_t rows = 0;
    while (std::getline(lines, line)) {
        rows++;
        EXPECT_EQ(std::count(line.begin(), line.end(), '\t'), 3) << line;
    }
    EXPECT_EQ(rows, features.size());
}

// A model's recorded versions: parsed; the model's features of another version (or none recorded) make one warning
// naming them; a model of this protal's versions none.
TEST(FeatureVersions, ParsedAndCompared) {
    feature_versions::Recorded recorded;
    ASSERT_EQ(feature_versions::Parse(feature_versions::ExtensionValue(), recorded), "");
    EXPECT_EQ(recorded.size(), std::size(feature_versions::kFeatureVersions));
    std::vector<std::string> const features = { "depth", "lu_per_kb", "failed_candidate_rate" };
    EXPECT_EQ(feature_versions::Mismatch("model_pe.xml", recorded, features), "");

    auto const depth = "depth:" + std::to_string(feature_versions::Find("depth")->version);  // as this protal's
    ASSERT_EQ(feature_versions::Parse(depth + ",lu_per_kb:1,failed_candidate_rate:1,other:3", recorded), "");
    auto const lu = feature_versions::Find("lu_per_kb");
    auto const failed = feature_versions::Find("failed_candidate_rate");
    ASSERT_TRUE(lu && failed);
    ASSERT_GE(lu->version, 2);
    ASSERT_GE(failed->version, 2);
    auto const warning = feature_versions::Mismatch("model_pe.xml", recorded, features);
    EXPECT_EQ(warning, "the model model_pe.xml was trained on other definitions of feature(s) lu_per_kb, failed_candidate_rate "
                       "(versions in the model 1, 1; in this protal " + std::to_string(lu->version) + ", " +
                       std::to_string(failed->version) + ", protal --feature_versions): retrain it, or expect calls that differ "
                       "from its training");
    ASSERT_EQ(feature_versions::Parse(depth, recorded), "");
    EXPECT_NE(feature_versions::Mismatch("m.xml", recorded, { "depth", "breadth_ratio" }).find("records no version of feature(s) breadth_ratio"),
              std::string::npos);
    auto const comparison = feature_versions::Compare(recorded, { "depth", "made_up" });
    EXPECT_TRUE(comparison.differ.empty());
    ASSERT_TRUE(feature_versions::Parse("made_up:2", recorded).empty());
    EXPECT_EQ(feature_versions::Compare(recorded, { "made_up" }).protal_versions, std::vector<int>{ 0 });

    for (auto const* bad : { "depth", "depth:0", "depth:x", ":1", "depth:1,depth:2", "depth:1,", "depth:1.5" }) {
        EXPECT_NE(feature_versions::Parse(bad, recorded), "") << bad;
        EXPECT_TRUE(recorded.empty()) << bad;
    }
    EXPECT_EQ(feature_versions::Parse("", recorded), "");
}

// The features' training ranges: parsed (five fields each, in order), kept for the model's inputs; a taxon's values
// outside [min, max] are hard, outside [p01, p99] soft; a subnormal value counts as 0 and a NaN is not compared.
TEST(FeatureRanges, ParsedAndChecked) {
    feature_ranges::Ranges ranges;
    ASSERT_EQ(feature_ranges::Parse("depth:0:0.5:90:200,excess_scaled_median:-0.05:-0.01:0.03:0.2,identity:0.8:0.9:1:1", ranges), "");
    ASSERT_EQ(ranges.size(), 3u);
    EXPECT_EQ(ranges.at("excess_scaled_median").p01, -0.01);
    for (auto const* bad : { "depth:0:1:2", "depth:0:1:2:3:4", "depth:2:1:3:4", "depth:0:1:2:3,depth:0:1:2:3", ":0:1:2:3",
                             "depth:0:1:x:3", "depth:0:1:nan:3", "depth:0:1:2:3," }) {
        EXPECT_NE(feature_ranges::Parse(bad, ranges), "") << bad;
        EXPECT_TRUE(ranges.empty()) << bad;
    }
    ASSERT_EQ(feature_ranges::Parse("depth:0:0.5:90:200,excess_scaled_median:-0.05:-0.01:0.03:0.2,identity:0.8:0.9:1:1", ranges), "");
    auto const kept = feature_ranges::Restricted(ranges, { "depth", "excess_scaled_median", "unranged" });
    EXPECT_EQ(kept.size(), 2u);

    std::vector<std::pair<std::string, double>> features = {
        { "depth", 250 }, { "excess_scaled_median", -0.99 }, { "identity", 0.85 }, { "unranged", 1e9 } };
    auto outside = feature_ranges::Check(ranges, features);
    EXPECT_EQ(outside.hard, (std::vector<std::string>{ "depth", "excess_scaled_median" }));
    EXPECT_EQ(outside.soft, (std::vector<std::string>{ "depth", "excess_scaled_median", "identity" }));
    EXPECT_EQ(feature_ranges::Joined(outside.hard), "depth,excess_scaled_median");
    features = { { "depth", std::numeric_limits<double>::denorm_min() }, { "identity", std::nan("") } };
    outside = feature_ranges::Check(ranges, features);
    EXPECT_TRUE(outside.hard.empty());
    EXPECT_EQ(outside.soft, std::vector<std::string>{ "depth" });  // 0, below p01 0.5
    EXPECT_EQ(feature_ranges::Joined({}), "");
    // A value at a bound of training is inside, though the trainer wrote the bound to 9 digits (%.9g): 2/3 as
    // 0.666666667 (above it), 1/3 as 0.333333333 (below it); a value beyond the rounding is outside.
    feature_ranges::Ranges rounded;
    ASSERT_EQ(feature_ranges::Parse("share:0.666666667:0.7:0.8:0.9,rate:0:0.1:0.2:0.333333333", rounded), "");
    features = { { "share", 2.0 / 3.0 }, { "rate", 1.0 / 3.0 } };
    EXPECT_TRUE(feature_ranges::Check(rounded, features).hard.empty());
    features = { { "share", 0.6666666 }, { "rate", 0.3333334 } };
    EXPECT_EQ(feature_ranges::Check(rounded, features).hard, (std::vector<std::string>{ "share", "rate" }));

    EXPECT_TRUE(feature_ranges::NearKnob(0.31, 0.5, 0));
    EXPECT_FALSE(feature_ranges::NearKnob(0.29, 0.5, 99));
    EXPECT_TRUE(feature_ranges::NearKnob(0.0, 0.5, 100));
}

// A sample's line: only with a hard one, the features by how many taxa they are outside of, the first five.
TEST(FeatureRanges, SampleLine) {
    feature_ranges::SampleSummary summary;
    summary.Add({});
    EXPECT_EQ(summary.Line("s1"), "");
    summary.Add({ { "a", "b" }, { "a", "b", "c" } });
    summary.Add({ { "b" }, { "b" } });
    for (auto const* name : { "c", "d", "e", "f" }) summary.Add({ { name }, { name } });
    auto const line = summary.Line("s1");
    EXPECT_EQ(line.rfind("Sample s1: 6 of 7 taxa near or above the knob have features outside the model's training range (b: 2 "
                         "taxa, a: 1 taxa, c: 1 taxa, d: 1 taxa, e: 1 taxa, and 1 more feature(s)); their probabilities are "
                         "extrapolated (6 have features outside", 0), 0u) << line;
}

// The model's input features: its active mining fields, each once, not its target.
TEST(FeatureRanges, ModelActiveFeatures) {
    std::string const xml = "<PMML><Header><Extension name=\"protal_feature_ranges\" value=\"depth:0:0:1:1\"/></Header>"
                            "<MiningSchema><MiningField name=\"label\" usageType=\"target\"/><MiningField name=\"depth\"/>"
                            "<MiningField name=\"identity\" usageType=\"active\"/></MiningSchema><Segment><MiningSchema>"
                            "<MiningField name=\"depth\"/><MiningField name=\"x\" usageType=\"supplementary\"/></MiningSchema>"
                            "</Segment></PMML>";
    EXPECT_EQ(feature_ranges::ModelActiveFeatures(xml), (std::vector<std::string>{ "depth", "identity" }));
    EXPECT_EQ(profiler::HeaderExtension(xml, feature_ranges::kExtension), std::optional<std::string>("depth:0:0:1:1"));
}
