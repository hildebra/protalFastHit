// The model rules protal shares with the trainer (scripts/machine_learning_cmdline.py) on the golden vectors of
// tests/data/golden_model_rules.tsv, which scripts/test_model_pmml.py checks the trainer on: the prior adjusted to a
// sample (context::SampleAdjusted), the calls at a target share of false calls (context::FalseCallKnob) and the knob
// curve over the sample's depth (profiler::DepthKnobAt).
#include <gtest/gtest.h>
#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include "Profiling/Profiler.h"
#include "Profiling/SampleContext.h"

using namespace protal;
namespace ctx = protal::profiler::context;

namespace {
    // tests/data, from the build's compile definition or else beside this file.
    std::filesystem::path DataDir() {
#ifdef PROTAL_TESTS_DIR
        return std::filesystem::path(PROTAL_TESTS_DIR) / "data";
#else
        return std::filesystem::path(__FILE__).parent_path() / "data";
#endif
    }

    // A case of the file: its rule, its name and its name=value fields.
    struct Case {
        std::string rule;
        std::string name;
        std::map<std::string, std::string> fields;

        std::string const& Field(std::string const& key) const {
            auto const it = fields.find(key);
            if (it == fields.end()) throw std::runtime_error(rule + " " + name + ": no field " + key);
            return it->second;
        }
        double Number(std::string const& key) const { return std::stod(Field(key)); }
    };

    std::vector<Case> Cases(std::string const& rule) {
        auto const path = DataDir() / "golden_model_rules.tsv";
        std::ifstream in(path);
        if (!in) throw std::runtime_error("cannot read " + path.string());
        std::vector<Case> cases;
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::stringstream ss(line);
            Case c;
            std::string field;
            std::getline(ss, c.rule, '\t');
            std::getline(ss, c.name, '\t');
            while (std::getline(ss, field, '\t')) {
                auto const eq = field.find('=');
                if (eq == std::string::npos) throw std::runtime_error("no name=value field: " + field);
                c.fields[field.substr(0, eq)] = field.substr(eq + 1);
            }
            if (c.rule == rule) cases.push_back(std::move(c));
        }
        return cases;
    }

    // "0.1,0.2*3" -> { 0.1, 0.2, 0.2, 0.2 }
    std::vector<double> Numbers(std::string const& text) {
        std::vector<double> values;
        std::stringstream ss(text);
        std::string item;
        while (std::getline(ss, item, ',')) {
            if (item.empty()) continue;
            auto const star = item.find('*');
            size_t const copies = star == std::string::npos ? 1 : std::stoul(item.substr(star + 1));
            values.insert(values.end(), copies, std::stod(item.substr(0, star)));
        }
        return values;
    }

    // "2:0.2,4:0.8" -> { {2, 0.2}, {4, 0.8} }
    std::vector<std::pair<double, double>> Curve(std::string const& text) {
        std::vector<std::pair<double, double>> curve;
        std::stringstream ss(text);
        std::string item;
        while (std::getline(ss, item, ',')) {
            if (item.empty()) continue;
            auto const colon = item.find(':');
            curve.emplace_back(std::stod(item.substr(0, colon)), std::stod(item.substr(colon + 1)));
        }
        return curve;
    }
}

TEST(GoldenModelRules, TheFileHoldsEveryRule) {
    EXPECT_GE(Cases("SampleAdjusted").size(), 5u);
    EXPECT_GE(Cases("FalseCallKnob").size(), 8u);
    EXPECT_GE(Cases("DepthKnobAt").size(), 8u);
}

TEST(GoldenModelRules, SampleAdjusted) {
    for (auto const& c : Cases("SampleAdjusted")) {
        double rate = -1;
        auto const adjusted = ctx::SampleAdjusted(Numbers(c.Field("q")), c.Number("prior"), &rate);
        auto const expected = Numbers(c.Field("adjusted"));
        EXPECT_NEAR(rate, c.Number("rate"), 1e-9) << c.name;
        ASSERT_EQ(adjusted.size(), expected.size()) << c.name;
        for (size_t i = 0; i < expected.size(); i++) EXPECT_NEAR(adjusted[i], expected[i], 1e-9) << c.name << " " << i;
    }
}

TEST(GoldenModelRules, FalseCallKnob) {
    for (auto const& c : Cases("FalseCallKnob")) {
        auto const scores = Numbers(c.Field("scores"));
        auto const calls = ctx::FalseCallKnob(scores, Curve(c.Field("curve")), c.Number("prior"), c.Number("fdr"));
        // protal calls every taxon whose score reaches the knob: those tied with the last one counted too.
        auto const called = std::count_if(scores.begin(), scores.end(), [&](double s) { return s >= calls.knob; });
        EXPECT_EQ(static_cast<size_t>(called), std::stoul(c.Field("called"))) << c.name;
        if (c.Field("knob") == "none") {
            EXPECT_EQ(calls.called, 0u) << c.name;
            EXPECT_GT(calls.knob, 1.0) << c.name;
        } else {
            EXPECT_EQ(calls.knob, c.Number("knob")) << c.name;
        }
        EXPECT_NEAR(calls.sample_prior, c.Number("rate"), 1e-9) << c.name;
    }
}

TEST(GoldenModelRules, DepthKnobAt) {
    for (auto const& c : Cases("DepthKnobAt")) {
        auto const fragments = static_cast<size_t>(std::stoull(c.Field("fragments")));
        EXPECT_NEAR(profiler::DepthKnobAt(Curve(c.Field("curve")), fragments), c.Number("knob"), 1e-12) << c.name;
    }
}
