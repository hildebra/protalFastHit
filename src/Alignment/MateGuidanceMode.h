// MateGuidanceMode.h - how a pair whose mates did not align together as its best candidate is settled
// (--mate_guidance, MateGuidance.h). consensus: both mates vote, each taxon of either mate's candidates scored by the
// sum of the mates' best scores on it, a mate without one looked for there (VoteMates); guide: a mate sure of its
// alignment (MAPQ 20 or more against its own others) guides the other one to its taxon, and the two mates of a pair on
// two genes take their candidates' consensus (GuideMate, ConsensusOfRead; protal 0.7.1 to 0.8.4); off: neither guide
// nor rescue, the consensus of the candidates only. The models' features depend on it, so a database records the
// setting its models were trained with (build_settings.tsv, key mate_guidance, these names).
#pragma once

#include <optional>
#include <string_view>

namespace protal {
    enum class MateGuidanceMode { Off, Guide, Consensus };

    // A run's mode without --mate_guidance: consensus since 2026-10-15 (guide before; the calls of the reduced r226
    // database's test sets and mocks did not change with models trained on guide, and fewer fragments go to the wrong
    // species: docs/claude/2026-10-15-mate-consensus).
    inline constexpr MateGuidanceMode kDefaultMateGuidance = MateGuidanceMode::Consensus;

    // The name of `mode` in --mate_guidance and build_settings.tsv.
    inline std::string_view MateGuidanceName(MateGuidanceMode mode) {
        switch (mode) {
            case MateGuidanceMode::Off: return "off";
            case MateGuidanceMode::Guide: return "guide";
            default: return "consensus";
        }
    }

    // The mode named `name` (consensus, guide, off), or nothing if it names none.
    inline std::optional<MateGuidanceMode> ParseMateGuidance(std::string_view name) {
        if (name == "consensus") return MateGuidanceMode::Consensus;
        if (name == "guide") return MateGuidanceMode::Guide;
        if (name == "off") return MateGuidanceMode::Off;
        return std::nullopt;
    }
}
