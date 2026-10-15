// SamSettings.h - how a SAM was aligned, in a line of its header (since 2026-10-15): "@CO\tprotal settings: " and the
// settings the run aligned the sample with, key=value separated by spaces:
//   read_type=pe gate=4 pick=stride score=14 hits=1 unique_kmers=cores traceback=singletrack allele_scores=on
//   weighted_site_shift=off mate_guidance=on align_top=3 adaptive_candidates=7 max_score_ani=0.9
// (allele_scores, weighted_site_shift and adaptive_candidates for short reads only, mate_guidance for paired-end reads
// only: the others do not use them). --profile_only reads it: the sample is profiled with the database's model nearest
// the SAM's gate, and is warned when the SAM was aligned otherwise than the database's models were trained
// (RunProtal: ApplyRunSettings). A SAM without the line (an older protal's) is profiled as before. A later protal's keys
// and values are kept as text: the reader neither refuses nor interprets them.
#pragma once

#include <istream>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace protal::sam_settings {
    inline const std::string kSamSettingsComment = "@CO\tprotal settings: ";

    // The keys, in the line's order.
    inline constexpr std::string_view kReadType = "read_type";
    inline constexpr std::string_view kGate = "gate";
    inline constexpr std::string_view kPick = "pick";
    inline constexpr std::string_view kScore = "score";
    inline constexpr std::string_view kHits = "hits";
    inline constexpr std::string_view kMinPicks = "min_picks";  // the stride pick's least picks of a read (since 2026-10-15)
    inline constexpr std::string_view kUniqueKmers = "unique_kmers";
    inline constexpr std::string_view kTraceback = "traceback";
    inline constexpr std::string_view kAlleleScores = "allele_scores";
    inline constexpr std::string_view kWeightedSiteShift = "weighted_site_shift";
    inline constexpr std::string_view kMateGuidance = "mate_guidance";
    inline constexpr std::string_view kAlignTop = "align_top";
    inline constexpr std::string_view kAdaptiveCandidates = "adaptive_candidates";
    inline constexpr std::string_view kMaxScoreAni = "max_score_ani";

    using Fields = std::vector<std::pair<std::string, std::string>>;

    // The header line of `fields`, with its newline. A value is written as it is: protal's values hold no space.
    inline std::string SamSettingsLine(Fields const& fields) {
        std::string line = kSamSettingsComment;
        for (size_t i = 0; i < fields.size(); i++) line += (i ? " " : "") + fields[i].first + "=" + fields[i].second;
        return line + '\n';
    }

    // The fields of a header line (without its newline) that starts with kSamSettingsComment; nullopt for another line.
    // A field without '=' or without a key is skipped; a key given twice keeps its last value.
    inline std::optional<Fields> ParseSamSettingsLine(std::string_view line) {
        if (!line.empty() && line.back() == '\r') line.remove_suffix(1);
        if (line.substr(0, kSamSettingsComment.size()) != kSamSettingsComment) return std::nullopt;
        std::string_view rest = line.substr(kSamSettingsComment.size());
        Fields fields;
        while (!rest.empty()) {
            size_t const space = rest.find(' ');
            std::string_view const field = rest.substr(0, space);
            rest = space == std::string_view::npos ? std::string_view() : rest.substr(space + 1);
            size_t const equals = field.find('=');
            if (equals == std::string_view::npos || equals == 0) continue;
            std::string key(field.substr(0, equals)), value(field.substr(equals + 1));
            bool replaced = false;
            for (auto& [k, v] : fields) {
                if (k != key) continue;
                v = value;
                replaced = true;
            }
            if (!replaced) fields.emplace_back(std::move(key), std::move(value));
        }
        return fields;
    }

    // The value of `key`, or nothing.
    inline std::optional<std::string> Get(Fields const& fields, std::string_view key) {
        for (auto const& [k, v] : fields) {
            if (k == key) return v;
        }
        return std::nullopt;
    }

    // The settings line of a SAM stream's header (read up to its first record), or nothing without one.
    inline std::optional<Fields> ReadSamSettings(std::istream& is) {
        std::string line;
        while (std::getline(is, line)) {
            if (line.empty() || line == "\r") continue;
            if (line[0] != '@') break;
            if (auto fields = ParseSamSettingsLine(line)) return fields;
        }
        return std::nullopt;
    }
}
