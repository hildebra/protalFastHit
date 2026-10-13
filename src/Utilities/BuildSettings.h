// BuildSettings.h - build_settings.tsv: the settings a database was built with that a run should match. --build writes
// it into the database folder and packs it into database.protal; a run reads it before aligning. It holds the WFA2
// traceback (WfaTraceback.h): the strain alleles, congener gaps and column weights are aligned with it, and the models
// are trained on runs that align with it; how the models' `fragments` feature counts (FragmentCount, since
// 2026-10-12; mapq when not said); and the sparse-first gate of each read type (--sparse_first, since 2026-10-12; off
// when not said) and its pick (--sparse_first_pick, since 2026-10-13; when not said hash if the gate is said, else
// stride: DatabasePick): the reads it gives up count no failed candidates, so the models are trained on runs gated
// alike. A
// database built before 2026-10-11 has no build_settings.tsv; it was built with the classic traceback.
#pragma once

#include <array>
#include <bit>
#include <charconv>
#include <cstddef>
#include <cstdlib>
#include <istream>
#include <optional>
#include <ostream>
#include <string>
#include <string_view>
#include <vector>

#include "ReadType.h"
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

    // The sparse-first gate (--sparse_first, ChainAnchorFinder::SetSparseFirst) of each read type, by ReadType's index:
    // C (1/C of a read's cores looked up first), 0 off; nothing where not said.
    using SparseFirstGates = std::array<std::optional<size_t>, kReadTypeCount>;

    // A gate's C: 0 (off), 1, 2, 4, 8, 16 or 32.
    inline std::optional<size_t> ParseSparseFirstLevel(std::string_view text) {
        size_t value = 0;
        auto const [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
        if (error != std::errc() || end != text.data() + text.size() || value > 32 || (value & (value - 1)) != 0) {
            return std::nullopt;
        }
        return value;
    }

    // A value per read type, by ReadType's index, from "auto" (nothing for every read type), one value for every read
    // type, or read types and their values, comma-separated ("pe:1,se:1,pb:8,ont:1"; a read type not listed stays
    // nothing). `parse` reads a value (nothing if invalid). Nothing if the text is none of these.
    template<typename T, typename Parse>
    std::optional<std::array<std::optional<T>, kReadTypeCount>> ParsePerReadType(std::string_view text, Parse const& parse) {
        std::array<std::optional<T>, kReadTypeCount> values;
        if (text == "auto") return values;
        if (text.find(':') == std::string_view::npos) {
            std::optional<T> const value = parse(text);
            if (!value) return std::nullopt;
            values.fill(value);
            return values;
        }
        while (!text.empty()) {
            size_t const comma = text.find(',');
            std::string_view const item = text.substr(0, comma);
            text = comma == std::string_view::npos ? std::string_view() : text.substr(comma + 1);
            if (comma != std::string_view::npos && text.empty()) return std::nullopt;  // a trailing comma
            size_t const colon = item.find(':');
            if (colon == std::string_view::npos) return std::nullopt;
            auto const type = ReadTypeFromToken(std::string(item.substr(0, colon)));
            std::optional<T> const value = parse(item.substr(colon + 1));
            if (!type || !value) return std::nullopt;
            values[static_cast<size_t>(*type)] = value;
        }
        return values;
    }

    // --sparse_first's value (and a map's SPARSE_FIRST): "auto" (nothing given: every read type the database's), one C
    // for every read type, or read types and their C, comma-separated ("pe:1,se:1,pb:8,ont:1"; a read type not listed
    // stays as the database has it). Nothing if the text is none of these.
    inline std::optional<SparseFirstGates> ParseSparseFirstGates(std::string_view text) {
        return ParsePerReadType<size_t>(text, ParseSparseFirstLevel);
    }

    // Which 1/C of a read's cores a sparse-first gate looks up first (--sparse_first_pick, since 2026-10-13): those whose
    // hash falls in 1/C of its range (stage1_probe::Level; the only pick before), every C-th in read order
    // (stage1_probe::StridePicked), spread evenly over the read, or (grid) the 1/C of each gene's cores the database keeps
    // at C (GridLevels.h): every core looked up, only those entries counting, and a kept read seeded from them alone, as a
    // sparse index of 1/C of the values (--sparse_index) holds them.
    enum class SparseFirstPick { Hash, Stride, Grid };
    using SparseFirstPicks = std::array<std::optional<SparseFirstPick>, kReadTypeCount>;

    inline std::optional<SparseFirstPick> ParseSparseFirstPick(std::string_view name) {
        if (name == "hash") return SparseFirstPick::Hash;
        if (name == "stride") return SparseFirstPick::Stride;
        if (name == "grid") return SparseFirstPick::Grid;
        return std::nullopt;
    }

    inline std::string_view SparseFirstPickName(SparseFirstPick pick) {
        return pick == SparseFirstPick::Grid ? "grid" : pick == SparseFirstPick::Stride ? "stride" : "hash";
    }

    // --sparse_first_pick's value: "auto", hash, stride or grid for every read type, or read types and theirs
    // ("pe:stride,pb:hash"), as --sparse_first's.
    inline std::optional<SparseFirstPicks> ParseSparseFirstPicks(std::string_view text) {
        return ParsePerReadType<SparseFirstPick>(text, ParseSparseFirstPick);
    }

    // "pe stride, pb hash" (the read types given), for messages; "none" if none is given.
    inline std::string SparseFirstPicksText(SparseFirstPicks const& picks) {
        std::string text;
        for (auto const& info : kReadTypes) {
            auto const& pick = picks[static_cast<size_t>(info.type)];
            if (pick) text += (text.empty() ? "" : ", ") + info.token + " " + std::string(SparseFirstPickName(*pick));
        }
        return text.empty() ? "none" : text;
    }

    // Of the gates a read type has models of (`trained`, not empty), the one nearest `gate` on the scale of C's log2,
    // off (0) one step below C 1; of two as near, the lower (the more sensitive) one. A run gated at `gate` is profiled
    // with that gate's model (Options::ModelGate).
    inline size_t NearestGate(size_t gate, std::vector<size_t> const& trained) {
        auto const step = [](size_t c) { return c == 0 ? -1 : static_cast<int>(std::bit_width(c)) - 1; };
        size_t best = trained.front();
        for (size_t c : trained) {
            int const d = std::abs(step(c) - step(gate)), best_d = std::abs(step(best) - step(gate));
            if (d < best_d || (d == best_d && c < best)) best = c;
        }
        return best;
    }

    // "pe 1, se 1, pb 8, ont 1" (the read types given; "off" for 0), for messages; "none" if none is given.
    inline std::string SparseFirstGatesText(SparseFirstGates const& gates) {
        std::string text;
        for (auto const& info : kReadTypes) {
            auto const& gate = gates[static_cast<size_t>(info.type)];
            if (!gate) continue;
            text += (text.empty() ? "" : ", ") + info.token + " " + (*gate ? std::to_string(*gate) : std::string("off"));
        }
        return text.empty() ? "none" : text;
    }
}

namespace protal::build_settings {
    inline const std::string kFileName = "build_settings.tsv";
    inline constexpr std::string_view kWfaTraceback = "wfa_traceback";
    inline constexpr std::string_view kFragmentCount = "fragment_count";
    // sparse_first_pe, sparse_first_se, sparse_first_pb, sparse_first_ont: the read type's gate (C, 0 off); and the
    // gates' flank score and hits.
    inline constexpr std::string_view kSparseFirst = "sparse_first_";
    inline constexpr std::string_view kSparseFirstScore = "sparse_first_score";
    inline constexpr std::string_view kSparseFirstHits = "sparse_first_hits";
    // sparse_first_pick_pe .. sparse_first_pick_ont (since 2026-10-13): the read type's gate's pick, hash or stride. A
    // protal from before skips these lines (no read type "pick_pe") and gates by hash.
    inline constexpr std::string_view kSparseFirstPick = "sparse_first_pick_";

    struct Settings {
        std::optional<WfaTraceback> wfa_traceback;    // nothing if the file does not say
        std::optional<FragmentCount> fragment_count;  // nothing if the file does not say (mapq)
        SparseFirstGates sparse_first;                // per read type; nothing if the file does not say (off)
        std::optional<size_t> sparse_first_score;     // nothing if the file does not say (14)
        std::optional<size_t> sparse_first_hits;      // nothing if the file does not say (1)
        SparseFirstPicks sparse_first_pick;           // per read type; nothing if the file does not say (DatabasePick)
    };

    // The pick a database's runs of reads of type index `t` gate with by default, and its models were trained with: as
    // build_settings.tsv records it. A database that records the type's gate but not its pick was built in the two
    // days before the pick could be chosen (2026-10-12/13) and gated by hash; one that records neither was trained
    // ungated, and a run gated on it picks by stride, protal's default since 2026-10-13.
    inline SparseFirstPick DatabasePick(Settings const& settings, size_t t) {
        if (settings.sparse_first_pick[t]) return *settings.sparse_first_pick[t];
        return settings.sparse_first[t] ? SparseFirstPick::Hash : SparseFirstPick::Stride;
    }

    // A whole number for a setting's value, or nothing.
    inline std::optional<size_t> ParseCount(std::string_view text) {
        size_t value = 0;
        auto const [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
        if (error != std::errc() || end != text.data() + text.size() || text.empty()) return std::nullopt;
        return value;
    }

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
            } else if (name == kSparseFirstScore || name == kSparseFirstHits) {
                auto const count = ParseCount(value);
                if (!count || (name == kSparseFirstScore ? *count > 16 : *count == 0)) {
                    return "line " + std::to_string(number) + ": " + std::string(name) + " is " +
                           (name == kSparseFirstScore ? "0 to 16" : "1 or more") + ", not '" + std::string(value) + "'";
                }
                (name == kSparseFirstScore ? settings.sparse_first_score : settings.sparse_first_hits) = count;
            } else if (name.substr(0, kSparseFirstPick.size()) == kSparseFirstPick) {
                auto const type = ReadTypeFromToken(std::string(name.substr(kSparseFirstPick.size())));
                if (!type) continue;  // a later version's read type
                auto const pick = ParseSparseFirstPick(value);
                if (!pick) {
                    return "line " + std::to_string(number) + ": " + std::string(name) + " is hash, stride or grid, not '" +
                           std::string(value) + "'";
                }
                settings.sparse_first_pick[static_cast<size_t>(*type)] = pick;
            } else if (name.substr(0, kSparseFirst.size()) == kSparseFirst) {
                // Another read type than protal knows (a later version's) is skipped, as other unknown settings are.
                auto const type = ReadTypeFromToken(std::string(name.substr(kSparseFirst.size())));
                if (!type) continue;
                auto const level = ParseSparseFirstLevel(value);
                if (!level) {
                    return "line " + std::to_string(number) + ": " + std::string(name) + " is 0, 1, 2, 4, 8, 16 or 32, not '" +
                           std::string(value) + "'";
                }
                settings.sparse_first[static_cast<size_t>(*type)] = level;
            }
        }
        return "";
    }

    inline void Write(std::ostream& os, Settings const& settings) {
        os << "# The settings this protal database was built with that a run should match (protal --build)\n";
        if (settings.wfa_traceback) os << kWfaTraceback << '\t' << WfaTracebackName(*settings.wfa_traceback) << '\n';
        if (settings.fragment_count) os << kFragmentCount << '\t' << FragmentCountName(*settings.fragment_count) << '\n';
        for (auto const& info : kReadTypes) {
            auto const& gate = settings.sparse_first[static_cast<size_t>(info.type)];
            if (gate) os << kSparseFirst << info.token << '\t' << *gate << '\n';
        }
        if (settings.sparse_first_score) os << kSparseFirstScore << '\t' << *settings.sparse_first_score << '\n';
        if (settings.sparse_first_hits) os << kSparseFirstHits << '\t' << *settings.sparse_first_hits << '\n';
        for (auto const& info : kReadTypes) {
            auto const& pick = settings.sparse_first_pick[static_cast<size_t>(info.type)];
            if (pick) os << kSparseFirstPick << info.token << '\t' << SparseFirstPickName(*pick) << '\n';
        }
    }
}
