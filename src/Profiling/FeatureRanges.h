// FeatureRanges.h - the range of each feature a presence model was trained on, and the taxa whose features fall outside
// it (since 2026-10-15). The trainer writes them into the model's header:
//   <Extension name="protal_feature_ranges" value="name:min:p01:p99:max,name:min:p01:p99:max,..."/>
// per feature of the model its least value over the training rows, the 1st and the 99th percentile and the largest, as
// C-locale decimals (%.9g). A tree scores a value beyond the training's like the training's most extreme rows, so the
// probability of a taxon with such a value is an extrapolation: placeholder base qualities put every taxon of a sample
// beyond the long-read models' strongest feature (docs/claude/2026-10-15-f1-audit-since-0.6.0a.md, C3 and C4). When a
// sample is profiled, each taxon that is called, within 0.2 below the knob, or of at least 100 fragments is checked
// (NearKnob): a feature outside [min, max] is outside the training range (hard), outside [p01, p99] beyond its central 98%
// (soft). The sample gets one line when a taxon has a hard one (SampleSummary::Line), and each taxon its hard ones in the
// .profile.log column outside_training_range. A model without the extension is not checked.
// Also the model's input features (ModelActiveFeatures), which FeatureVersions.h checks too.
#pragma once

#include <algorithm>
#include <charconv>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace protal::feature_ranges {
    inline constexpr std::string_view kExtension = "protal_feature_ranges";

    // A taxon is checked when its score is within this of the knob or above it, or it has kCheckedFragments or more.
    inline constexpr double kKnobMargin = 0.2;
    inline constexpr size_t kCheckedFragments = 100;

    struct Range {
        double min = 0, p01 = 0, p99 = 0, max = 0;
    };
    using Ranges = std::map<std::string, Range, std::less<>>;

    // The ranges of the extension's value into `ranges` (empty first); an error text if it is malformed: an item without
    // its five fields, a number that is none (or not finite), the four out of order, a feature twice.
    inline std::string Parse(std::string_view text, Ranges& ranges) {
        ranges.clear();
        auto number = [](std::string_view field, double& value) {
            auto const [end, error] = std::from_chars(field.data(), field.data() + field.size(), value);
            return !field.empty() && error == std::errc() && end == field.data() + field.size() && std::isfinite(value);
        };
        std::string const expected = "name:min:p01:p99:max, min <= p01 <= p99 <= max, each feature once";
        while (!text.empty()) {
            size_t const comma = text.find(',');
            std::string_view const item = text.substr(0, comma);
            text = comma == std::string_view::npos ? std::string_view() : text.substr(comma + 1);
            std::vector<std::string_view> fields;
            for (std::string_view rest = item;;) {
                size_t const colon = rest.find(':');
                fields.push_back(rest.substr(0, colon));
                if (colon == std::string_view::npos) break;
                rest = rest.substr(colon + 1);
            }
            Range r;
            bool const trailing = comma != std::string_view::npos && text.empty();  // "a:0:0:1:1,"
            bool const ok = !trailing && fields.size() == 5 && !fields[0].empty() && number(fields[1], r.min) && number(fields[2], r.p01) &&
                            number(fields[3], r.p99) && number(fields[4], r.max) && r.min <= r.p01 && r.p01 <= r.p99 &&
                            r.p99 <= r.max && !ranges.contains(fields[0]);
            if (!ok) {
                ranges.clear();
                return "its feature ranges are malformed ('" + std::string(item) + "': " + expected + ")";
            }
            ranges.emplace(std::string(fields[0]), r);
        }
        return {};
    }

    // Only the ranges of `features` (the model's inputs).
    inline Ranges Restricted(Ranges const& ranges, std::vector<std::string> const& features) {
        Ranges kept;
        for (auto const& name : features) {
            auto const it = ranges.find(name);
            if (it != ranges.end()) kept.emplace(it->first, it->second);
        }
        return kept;
    }

    // Whether a taxon of model score `score` at knob `knob` with `fragments` is checked.
    inline bool NearKnob(double score, double knob, size_t fragments) {
        return score >= knob - kKnobMargin || fragments >= kCheckedFragments;
    }

    // The features of a taxon outside the training range: hard, outside [min, max]; soft, outside [p01, p99] (the hard
    // ones too). Each value as the model is given it (a subnormal value is 0, profiler::FeatureString); a value that is
    // no number is not compared. In the order of `features`.
    struct Outside {
        std::vector<std::string> hard, soft;
    };

    template<typename FeatureList>
    Outside Check(Ranges const& ranges, FeatureList const& features) {
        Outside outside;
        for (auto const& [name, raw] : features) {
            auto const it = ranges.find(name);
            if (it == ranges.end()) continue;
            double value = raw;
            if (std::isnan(value)) continue;
            if (std::fabs(value) < std::numeric_limits<double>::min()) value = 0;
            auto const& r = it->second;
            if (value < r.min || value > r.max) outside.hard.emplace_back(name);
            if (value < r.p01 || value > r.p99) outside.soft.emplace_back(name);
        }
        return outside;
    }

    // "a,b,c" (empty for none), the .profile.log column.
    inline std::string Joined(std::vector<std::string> const& names) {
        std::string text;
        for (auto const& name : names) text += (text.empty() ? "" : ",") + name;
        return text;
    }

    // What a sample's checked taxa (NearKnob) show.
    struct SampleSummary {
        size_t checked = 0;    // taxa checked
        size_t hard_taxa = 0;  // of them with a feature outside [min, max]
        size_t soft_taxa = 0;  // and with one outside [p01, p99]
        std::map<std::string, size_t> hard_counts;  // per feature, the taxa it is outside [min, max] of

        void Add(Outside const& outside) {
            checked++;
            hard_taxa += !outside.hard.empty();
            soft_taxa += !outside.soft.empty();
            for (auto const& name : outside.hard) hard_counts[name]++;
        }

        // "Sample S: K of M taxa near or above the knob have features outside the model's training range (feature: n
        // taxa, ... the 5 most frequent); their probabilities are extrapolated", and the taxa beyond the central 98% of
        // the training's values; "" when no taxon has a feature outside the range.
        std::string Line(std::string const& sample) const {
            if (hard_taxa == 0) return {};
            std::vector<std::pair<std::string, size_t>> features(hard_counts.begin(), hard_counts.end());
            std::stable_sort(features.begin(), features.end(), [](auto const& a, auto const& b) { return a.second > b.second; });
            std::string list;
            for (size_t i = 0; i < features.size() && i < 5; i++) {
                list += (i ? ", " : "") + features[i].first + ": " + std::to_string(features[i].second) + " taxa";
            }
            if (features.size() > 5) list += ", and " + std::to_string(features.size() - 5) + " more feature(s)";
            return "Sample " + sample + ": " + std::to_string(hard_taxa) + " of " + std::to_string(checked) + " taxa near or above "
                   "the knob have features outside the model's training range (" + list + "); their probabilities are "
                   "extrapolated (" + std::to_string(soft_taxa) + " have features outside the central 98% of the training's "
                   "values; the taxa's features in outside_training_range of the .profile.log)";
        }
    };

    // The active input features of the PMML model `xml` (its MiningFields that are no target), in its order.
    inline std::vector<std::string> ModelActiveFeatures(std::string const& xml) {
        auto attribute = [](std::string_view tag, std::string_view name) -> std::string {
            std::string const key = " " + std::string(name) + "=\"";
            auto start = tag.find(key);
            if (start == std::string_view::npos) return {};
            start += key.size();
            auto const end = tag.find('"', start);
            return end == std::string_view::npos ? std::string() : std::string(tag.substr(start, end - start));
        };
        std::vector<std::string> features;
        std::set<std::string, std::less<>> seen;
        for (size_t pos = xml.find("<MiningField "); pos != std::string::npos; pos = xml.find("<MiningField ", pos + 1)) {
            auto const end = xml.find('>', pos);
            if (end == std::string::npos) break;
            std::string_view const tag(xml.data() + pos, end - pos + 1);
            auto const name = attribute(tag, "name");
            auto usage = attribute(tag, "usageType");
            if (usage.empty()) usage = "active";
            if (usage != "active" || name.empty() || seen.contains(name)) continue;
            seen.insert(name);
            features.push_back(name);
        }
        return features;
    }
}
