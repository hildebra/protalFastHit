// ScannedReads.h - how many reads the aligner read for a sample, carried in the SAM header to the profile, whose
// composition (Composition.h) compares them with the reads the called species explain. protal writes the line after
// the reads are aligned, with the header (a SAM with --full_sam_header, whose header goes first, has none).
#pragma once

#include <charconv>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>

namespace protal {
    struct ScannedReads {
        uint64_t fragments = 0;  // read pairs, or single reads
        uint64_t reads = 0;      // reads (two per pair)
        uint64_t bases = 0;      // their bases
    };

    inline const std::string kSamScannedReadsComment = "@CO\tprotal scanned reads: ";

    // The header line of `scanned`, with its newline: "@CO\tprotal scanned reads: fragments=F reads=R bases=B".
    inline std::string ScannedReadsLine(ScannedReads const& scanned) {
        return kSamScannedReadsComment + "fragments=" + std::to_string(scanned.fragments) + " reads=" + std::to_string(scanned.reads) +
               " bases=" + std::to_string(scanned.bases) + '\n';
    }

    // The counts of a header line (without its newline) that starts with kSamScannedReadsComment; nullopt if a field
    // is missing or no number. Fields other than these three are ignored.
    inline std::optional<ScannedReads> ParseScannedReadsLine(std::string_view line) {
        if (line.substr(0, kSamScannedReadsComment.size()) != kSamScannedReadsComment) return std::nullopt;
        std::string_view rest = line.substr(kSamScannedReadsComment.size());
        ScannedReads scanned;
        bool fragments = false, reads = false, bases = false;
        while (!rest.empty()) {
            size_t const space = rest.find(' ');
            std::string_view const field = rest.substr(0, space);
            rest = space == std::string_view::npos ? std::string_view() : rest.substr(space + 1);
            size_t const equals = field.find('=');
            if (equals == std::string_view::npos) continue;
            std::string_view const key = field.substr(0, equals), text = field.substr(equals + 1);
            uint64_t* const target = key == "fragments" ? &scanned.fragments : key == "reads" ? &scanned.reads :
                                     key == "bases" ? &scanned.bases : nullptr;
            if (!target) continue;
            auto const [p, ec] = std::from_chars(text.data(), text.data() + text.size(), *target);
            if (ec != std::errc() || p != text.data() + text.size() || text.empty()) return std::nullopt;
            (target == &scanned.fragments ? fragments : target == &scanned.reads ? reads : bases) = true;
        }
        if (!fragments || !reads || !bases) return std::nullopt;
        return scanned;
    }
}
