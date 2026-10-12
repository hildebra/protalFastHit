// BuildSettings.h - build_settings.tsv: the settings a database was built with that a run should match. --build writes
// it into the database folder and packs it into database.protal; a run reads it before aligning. It holds the WFA2
// traceback (WfaTraceback.h): the strain alleles, congener gaps and column weights are aligned with it, and the models
// are trained on runs that align with it; and how the models' `fragments` feature counts (FragmentCount, since
// 2026-10-12; mapq when not said). A database built before 2026-10-11 has no build_settings.tsv; it was built with the
// classic traceback.
#pragma once

#include <istream>
#include <optional>
#include <ostream>
#include <string>
#include <string_view>

#include "WfaTraceback.h"

namespace protal {
    // How the model's `fragments` feature counts a taxon's fragments (--fragment_count): mapq, those the MAPQ and length
    // filters keep (as before 2026-10-12); soft, those plus the filtered ones by the share the read EM leaves the taxon
    // (a divergent strain's reads that tie with a congener's reference fall to MAPQ 0 and counted for nothing).
    enum class FragmentCount { Mapq, Soft };

    inline std::optional<FragmentCount> ParseFragmentCount(std::string_view name) {
        if (name == "mapq") return FragmentCount::Mapq;
        if (name == "soft") return FragmentCount::Soft;
        return std::nullopt;
    }

    inline std::string_view FragmentCountName(FragmentCount count) {
        return count == FragmentCount::Soft ? "soft" : "mapq";
    }
}

namespace protal::build_settings {
    inline const std::string kFileName = "build_settings.tsv";
    inline constexpr std::string_view kWfaTraceback = "wfa_traceback";
    inline constexpr std::string_view kFragmentCount = "fragment_count";

    struct Settings {
        std::optional<WfaTraceback> wfa_traceback;    // nothing if the file does not say
        std::optional<FragmentCount> fragment_count;  // nothing if the file does not say (mapq)
    };

    // build_settings.tsv: a line per setting, its name and value tab-separated; lines starting with '#' and empty lines
    // are skipped, as are settings of later protal versions this one does not know. Returns an error text, or "" when
    // every line was read.
    inline std::string Read(std::istream& is, Settings& settings) {
        std::string line;
        size_t number = 0;
        while (std::getline(is, line)) {
            number++;
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty() || line[0] == '#') continue;
            size_t const tab = line.find('\t');
            if (tab == std::string::npos) return "line " + std::to_string(number) + ": no tab between name and value";
            std::string_view const name = std::string_view(line).substr(0, tab);
            std::string_view const value = std::string_view(line).substr(tab + 1);
            if (name == kWfaTraceback) {
                auto const traceback = ParseWfaTraceback(value);
                if (!traceback) {
                    return "line " + std::to_string(number) + ": " + std::string(kWfaTraceback) + " is classic or singletrack, not '" +
                           std::string(value) + "'";
                }
                settings.wfa_traceback = traceback;
            } else if (name == kFragmentCount) {
                auto const count = ParseFragmentCount(value);
                if (!count) {
                    return "line " + std::to_string(number) + ": " + std::string(kFragmentCount) + " is mapq or soft, not '" +
                           std::string(value) + "'";
                }
                settings.fragment_count = count;
            }
        }
        return "";
    }

    inline void Write(std::ostream& os, Settings const& settings) {
        os << "# The settings this protal database was built with that a run should match (protal --build)\n";
        if (settings.wfa_traceback) os << kWfaTraceback << '\t' << WfaTracebackName(*settings.wfa_traceback) << '\n';
        if (settings.fragment_count) os << kFragmentCount << '\t' << FragmentCountName(*settings.fragment_count) << '\n';
    }
}
