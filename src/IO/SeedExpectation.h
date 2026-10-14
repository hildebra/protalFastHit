// SeedExpectation.h - how a sample's unique seeds compare with its genes' unique k-mers, carried in the SAM header to the
// profile. A read's ZU and ZT count the unique values among the seeds of its anchor (ZS); the share of a gene's index values
// that are unique (unique_kmers.tsv) is what a read of the species finds where every core is looked up. A read seeded from
// a grid level alone (a sparse index, --sparse_index emulate) finds that level's values, whose unique share the grid's
// choice of cores changes: the factors scale the genes' shares to it (RunProtal::SeedExpectations; the features
// unique_seed_enrichment and unique2_seed_enrichment, docs/claude/2026-10-14-sparse-index-calls). Without the line (a run
// seeded from every core, a SAM of an older protal) both are 1.
#pragma once

#include <charconv>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>

namespace protal {
    struct SeedExpectation {
        uint32_t level = 0;       // the grid level the sample's reads were seeded from (0: every value)
        double unique = 1.0;      // the level's unique share of values over the whole index's (ZU)
        double unique_two = 1.0;  // and of those unique at distance two (ZT)
    };

    inline const std::string kSamSeedExpectationComment = "@CO\tprotal seed expectation: ";

    // The header line of `e`, with its newline: "@CO\tprotal seed expectation: level=L unique=U two=T", the factors as the
    // shortest text that reads back to the same double.
    inline std::string SeedExpectationLine(SeedExpectation const& e) {
        auto text = [](double value) {
            char buf[32];
            auto const [end, ec] = std::to_chars(buf, buf + sizeof(buf), value);
            return ec == std::errc() ? std::string(buf, end) : std::to_string(value);
        };
        return kSamSeedExpectationComment + "level=" + std::to_string(e.level) + " unique=" + text(e.unique) + " two=" +
               text(e.unique_two) + '\n';
    }

    // The values of a header line (without its newline) that starts with kSamSeedExpectationComment; nullopt if a field is
    // missing, no number, or a factor not positive. Other fields are ignored.
    inline std::optional<SeedExpectation> ParseSeedExpectationLine(std::string_view line) {
        if (line.substr(0, kSamSeedExpectationComment.size()) != kSamSeedExpectationComment) return std::nullopt;
        std::string_view rest = line.substr(kSamSeedExpectationComment.size());
        SeedExpectation e;
        bool level = false, unique = false, two = false;
        while (!rest.empty()) {
            size_t const space = rest.find(' ');
            std::string_view const field = rest.substr(0, space);
            rest = space == std::string_view::npos ? std::string_view() : rest.substr(space + 1);
            size_t const equals = field.find('=');
            if (equals == std::string_view::npos) continue;
            std::string_view const key = field.substr(0, equals), text = field.substr(equals + 1);
            if (text.empty()) return std::nullopt;
            if (key == "level") {
                auto const [p, ec] = std::from_chars(text.data(), text.data() + text.size(), e.level);
                if (ec != std::errc() || p != text.data() + text.size()) return std::nullopt;
                level = true;
            } else if (key == "unique" || key == "two") {
                double& target = key == "unique" ? e.unique : e.unique_two;
                auto const [p, ec] = std::from_chars(text.data(), text.data() + text.size(), target);
                if (ec != std::errc() || p != text.data() + text.size() || !(target > 0)) return std::nullopt;
                (key == "unique" ? unique : two) = true;
            }
        }
        if (!level || !unique || !two) return std::nullopt;
        return e;
    }
}
