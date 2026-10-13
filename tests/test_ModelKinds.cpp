// Unit tests for the kinds of presence models a database holds (since 2026-10-13): gradient-boosted trees and a random
// forest of each read type at each gate, the kind told from a model's PMML, --model_kind as a run and --add_model take
// it, the model a sample takes of the kind asked for; and the notes database.protal carries for its users
// (build_metadata.tsv last, after the models).
#include <gtest/gtest.h>
#include <algorithm>
#include <filesystem>
#include <map>
#include <set>
#include <string>
#include <vector>
#include "Build.h"
#include "Profiling/Profiler.h"
#include "TestUtil.h"

using namespace protal;
using protal::test::ScratchDir;
namespace fs = std::filesystem;

namespace {
    // A model's PMML down to its first segmentation, combining its trees by `method` (nothing: a lone tree).
    std::string Pmml(std::string const& method, std::string const& header = "<Header/>") {
        std::string const segmentation = method.empty() ? "" : "<Segmentation multipleModelMethod=\"" + method + "\">";
        return "<?xml version=\"1.0\"?>\n<PMML version=\"4.4\">\n " + header + "\n <MiningModel functionName=\"classification\">" +
               segmentation + "<Segment id=\"1\"><True/><TreeModel functionName=\"regression\"><Node score=\"0\"><True/></Node>"
               "</TreeModel></Segment>" + (method.empty() ? "" : "</Segmentation>") + "</MiningModel>\n</PMML>\n";
    }

    std::vector<std::string> Names(std::vector<db::Source> const& sources) {
        std::vector<std::string> names;
        for (auto const& source : sources) names.push_back(source.name);
        return names;
    }
}

// A model of another kind sits beside the file of its read type and gate (model_pe_c1.xml -> model_pe_c1_forest.xml),
// and every such file is a model file of the database (packed last, replaced in place), each name once.
TEST(ModelKinds, EachModelFileHasOneOfEveryKind) {
    EXPECT_EQ(kModelKinds.size(), 2u);
    EXPECT_TRUE(IsModelKind("gbm"));
    EXPECT_TRUE(IsModelKind("forest"));
    EXPECT_FALSE(IsModelKind("auto"));
    EXPECT_FALSE(IsModelKind("rf"));
    EXPECT_EQ(ModelKindNames(), "gbm or forest");
    EXPECT_EQ(KindModelFile("model_pe.xml", "forest"), "model_pe_forest.xml");
    EXPECT_EQ(KindModelFile(GateModelFile(ReadType::PacBio, 32), "gbm"), "model_PB_c32_gbm.xml");
    auto const all = AllModelFiles();
    EXPECT_EQ(std::set<std::string>(all.begin(), all.end()).size(), all.size()) << "each file once";
    for (auto const& info : kReadTypes) {
        for (auto const& kind : kModelKinds) {
            EXPECT_NE(std::find(all.begin(), all.end(), KindModelFile(info.model_file, kind)), all.end());
            for (size_t gate : kModelGates) {
                auto const name = KindModelFile(GateModelFile(info.type, gate), kind);
                EXPECT_NE(std::find(all.begin(), all.end(), name), all.end()) << name;
            }
        }
    }
}

// The kind of a model, as its first segmentation combines its trees: chained into a logit (machine_learning_cmdline.py's
// boosting), averaged or voted (its forests, R's); a placeholder by its annotation; any other model has none.
TEST(ModelKinds, TheKindIsHowTheModelCombinesItsTrees) {
    EXPECT_EQ(profiler::ModelKindOf(Pmml("modelChain")), "gbm");
    EXPECT_EQ(profiler::ModelKindOf(Pmml("average")), "forest");
    EXPECT_EQ(profiler::ModelKindOf(Pmml("majorityVote")), "forest");
    EXPECT_EQ(profiler::ModelKindOf(Pmml("weightedAverage")), "forest");
    EXPECT_EQ(profiler::ModelKindOf(Pmml("sum")), "");
    EXPECT_EQ(profiler::ModelKindOf(Pmml("")), "");
    EXPECT_EQ(profiler::ModelKindOf(Pmml("average", "<Header><Annotation>protal:placeholder</Annotation></Header>")), "placeholder");
    // Only the first segmentation counts, not a later, nested one's.
    std::string nested = Pmml("average");
    nested.insert(nested.find("<TreeModel"), "<MiningModel functionName=\"regression\"><Segmentation multipleModelMethod=\"modelChain\">"
                                             "</Segmentation></MiningModel>");
    EXPECT_EQ(profiler::ModelKindOf(nested), "forest");
}

// --model_kind of a run: one kind for every read type, read types and theirs (the others auto), or auto: the database's
// default files (no kind).
TEST(ModelKinds, ARunsKindsByReadType) {
    std::string error;
    auto kinds = Options::ModelKindsOf("auto", error);
    ASSERT_TRUE(kinds) << error;
    for (auto const& kind : *kinds) EXPECT_EQ(kind, "");
    kinds = Options::ModelKindsOf("forest", error);
    ASSERT_TRUE(kinds) << error;
    for (auto const& kind : *kinds) EXPECT_EQ(kind, "forest");
    kinds = Options::ModelKindsOf("pe:forest,ont:gbm,se:auto", error);
    ASSERT_TRUE(kinds) << error;
    EXPECT_EQ(*kinds, (std::array<std::string, kReadTypeCount>{ "forest", "", "", "gbm" }));
    for (std::string const bad : { "rf", "pe:rf", "illumina:gbm", "pe", "pe:forest,,se:gbm", "" }) {
        error.clear();
        EXPECT_FALSE(Options::ModelKindsOf(bad, error)) << bad;
        EXPECT_NE(error.find("--model_kind is auto, gbm or forest"), std::string::npos) << error;
    }
}

// --add_model's models with their kinds at the same places ('-': the read type's default file), and without
// --model_kind (auto) all of them the default files.
TEST(ModelKinds, AddedModelsWithTheirKinds) {
    OptionsData d;
    d.add_model = "a.xml,b.xml,c.xml";
    d.read_type = "pe,pe,se";
    d.model_gate = "-,1,-";
    d.model_kind = "-,forest,gbm";
    auto const models = Options(d).AddModels();
    ASSERT_EQ(models.size(), 3u);
    EXPECT_EQ(models[0].kind, std::nullopt);
    EXPECT_EQ(models[1].kind, std::optional<std::string>("forest"));
    EXPECT_EQ(models[1].gate, std::optional<size_t>(1));
    EXPECT_EQ(models[2].kind, std::optional<std::string>("gbm"));
    EXPECT_EQ(models[2].type, ReadType::Single);
    d.model_kind = "auto";
    for (auto const& model : Options(d).AddModels()) EXPECT_EQ(model.kind, std::nullopt);
}

// A sample takes the model of the kind asked for at the gate nearest its own (SetModelGates' members), else the
// database's default file there; a model named (--model) stays.
TEST(ModelKinds, ASampleTakesTheModelOfItsKindAtItsGate) {
    ScratchDir dir;
    OptionsData d;
    d.database_path = dir.path.string();
    Options options(d);
    std::array<std::vector<size_t>, kReadTypeCount> gates;
    std::array<std::map<size_t, std::string>, kReadTypeCount> kinds;
    gates[0] = { 1, 4 };
    kinds[0] = { { 1, "model_pe_c1_forest.xml" }, { 4, "model_pe.xml" } };  // the database's default is a forest
    options.SetModelGates(gates, { 4, 0, 0, 0 }, kinds);
    auto const file = [&options]() { return fs::path(options.SampleModelDbFile(0).Path()).filename().string(); };
    options.SetSparseFirst({ 4, 0, 0, 0 }, 14, 1);
    EXPECT_EQ(file(), "model_pe.xml");
    options.SetSparseFirst({ 1, 0, 0, 0 }, 14, 1);
    EXPECT_EQ(file(), "model_pe_c1_forest.xml");
    options.SetSparseFirst({ 0, 0, 0, 0 }, 14, 1);  // off: nearest 1
    EXPECT_EQ(file(), "model_pe_c1_forest.xml");
    // Without a kind: the default files at the gates.
    options.SetModelGates(gates, { 4, 0, 0, 0 });
    EXPECT_EQ(file(), "model_pe_c1.xml");
    d.model = "mine";
    Options named(d);
    named.SetModelGates(gates, { 4, 0, 0, 0 }, kinds);
    named.SetSparseFirst({ 1, 0, 0, 0 }, 14, 1);
    EXPECT_EQ(fs::path(named.SampleModelDbFile(0).Path()).filename().string(), "mine.xml");
}

// The notes (build_metadata.tsv, gene2geneid.tsv, genome2tiid.tsv) are packed from the folder and stay there;
// --add_tables stores them; build_metadata.tsv goes last, after the models, which follow the other members.
TEST(ModelKinds, TheNotesArePackedAndBuildMetadataGoesLast) {
    auto const tables = Options::AddableTables();
    for (auto const& note : Options::NoteFiles()) {
        EXPECT_TRUE(Options::IsNoteFile(note)) << note;
        EXPECT_NE(std::find(tables.begin(), tables.end(), note), tables.end()) << note;
        EXPECT_TRUE(build::IsKeptSource({ note, "/x/" + note }));
    }
    EXPECT_FALSE(Options::IsNoteFile("reference.map"));
    EXPECT_FALSE(build::IsKeptSource({ "reference.map", "/x/reference.map" }));

    std::vector<db::Source> sources = { { "build_metadata.tsv", "m" }, { "model_pe_forest.xml", "f" }, { "index.prx", "i" },
                                        { "model_pe.xml", "p" }, { "gene2geneid.tsv", "g" }, { "reference.map", "r" },
                                        { "model_pe_c1_forest.xml", "c" } };
    build::ModelsLast(sources);
    EXPECT_EQ(Names(sources), (std::vector<std::string>{ "index.prx", "gene2geneid.tsv", "reference.map", "model_pe_forest.xml",
                                                         "model_pe.xml", "model_pe_c1_forest.xml", "build_metadata.tsv" }));
}

// What a database's models cover, for the logs of --build and --add_model: each read type's other gates, then its other
// kinds' (at its own gate and others); a read type without its default model is listed with what it has.
TEST(ModelKinds, TheModelsOfEachKindAreListed) {
    auto const [with, without] = build::ModelCoverage({ "model_pe.xml", "model_pe_c1.xml", "model_pe_forest.xml",
                                                        "model_pe_c1_forest.xml", "model_se.xml", "model_PB_c8_forest.xml" });
    EXPECT_EQ(with, "pe (also gates 1; forest at its gate, 1), se");
    EXPECT_EQ(without, "pb (forest at 8), ont");
    EXPECT_EQ(build::ModelCoverage({ "model_pe.xml", "model_pe_c4.xml" }).first, "pe (also gates 4)");
}
