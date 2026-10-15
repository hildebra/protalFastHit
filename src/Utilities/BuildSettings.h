// BuildSettings.h - build_settings.tsv: the settings a database was built with that a run should match. --build writes
// it into the database folder and packs it into database.protal; a run reads it before aligning. It holds the WFA2
// traceback (WfaTraceback.h): the strain alleles, congener gaps and column weights are aligned with it, and the models
// are trained on runs that align with it; how the models' `fragments` feature counts (FragmentCount, since
// 2026-10-12; mapq when not said); and the sparse-first gate of each read type (--sparse_first, since 2026-10-12; off
// when not said) and its pick (--sparse_first_pick, since 2026-10-13; when not said hash if the gate is said, else
// stride: DatabasePick): the reads it gives up count no failed candidates, so the models are trained on runs gated
// alike. Since 2026-10-15 also the run settings the models' features depend on (RunSettings: the site shift, the allele
// scores, --align_top, --adaptive_candidates, mate guidance, the identity floor of each read type, the unique k-mers'
// scheme, foreign genes, suspect copies, the depth identity margin), which runs follow when their options are left at
// auto and are warned about otherwise (Run.cpp ApplyRunSettings). A
// database built before 2026-10-11 has no build_settings.tsv; it was built with the classic traceback.
#pragma once

#include <array>
#include <bit>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <istream>
#include <optional>
#include <ostream>
#include <string>
#include <string_view>
#include <utility>
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

    // The run settings (since 2026-10-15), one line each: on/off, a whole number, a decimal (the shortest text that reads
    // back to the same double), or a mode's name. A protal from before skips them as settings it does not know.
    inline constexpr std::string_view kWeightedSiteShift = "weighted_site_shift";      // on, off
    inline constexpr std::string_view kAlleleScores = "allele_scores";                 // on, off
    inline constexpr std::string_view kAlignTop = "align_top";                         // N
    inline constexpr std::string_view kAdaptiveCandidates = "adaptive_candidates";     // N
    inline constexpr std::string_view kMateGuidance = "mate_guidance";                 // on, off (a later protal's modes too)
    inline constexpr std::string_view kMaxScoreAni = "max_score_ani_";                 // + pe, se, pb, ont: the identity floor
    inline constexpr std::string_view kUniqueKmers = "unique_kmers";                   // cores, anchors
    inline constexpr std::string_view kKeepForeignGenes = "keep_foreign_genes";        // on, off
    inline constexpr std::string_view kKeepSuspectCopies = "keep_suspect_copies";      // on, off
    inline constexpr std::string_view kDepthIdentityMargin = "depth_identity_margin";  // decimal
    // unique_kmers: ZS, ZU and ZT from the gene's cores an alignment covers (unique_cores.bin), or from its anchor's seeds.
    inline constexpr std::string_view kUniqueCores = "cores";
    inline constexpr std::string_view kUniqueAnchors = "anchors";
    // The mate guidance modes this protal knows; a mode of a later protal is kept as text (KnownMateGuidance).
    inline constexpr std::string_view kMateGuidanceOn = "on";
    inline constexpr std::string_view kMateGuidanceOff = "off";

    // The settings of a run that change its alignments or the features its models read: those --build ran with
    // (build_settings.tsv), those a run was given (nothing: auto, the database's), or those it uses. Each is nothing where
    // not said. mate_guidance and unique_kmers are kept as text, so that a later protal's modes read without an error.
    struct RunSettings {
        std::optional<bool> weighted_site_shift;
        std::optional<bool> allele_scores;
        std::optional<size_t> align_top;
        std::optional<size_t> adaptive_candidates;
        std::optional<std::string> mate_guidance;
        std::array<std::optional<double>, kReadTypeCount> max_score_ani{};  // the identity floor (-a) of each read type
        std::optional<std::string> unique_kmers;
        std::optional<bool> keep_foreign_genes;
        std::optional<bool> keep_suspect_copies;
        std::optional<double> depth_identity_margin;
        bool operator==(RunSettings const&) const = default;
    };

    struct Settings {
        std::optional<WfaTraceback> wfa_traceback;    // nothing if the file does not say
        std::optional<FragmentCount> fragment_count;  // nothing if the file does not say (mapq)
        SparseFirstGates sparse_first;                // per read type; nothing if the file does not say (off)
        std::optional<size_t> sparse_first_score;     // nothing if the file does not say (14)
        std::optional<size_t> sparse_first_hits;      // nothing if the file does not say (1)
        SparseFirstPicks sparse_first_pick;           // per read type; nothing if the file does not say (DatabasePick)
        RunSettings run;                              // each nothing if the file does not say (built before 2026-10-15)
    };

    // A decimal as the shortest text that reads back to the same double, in the C locale ("0.9", "0.85").
    inline std::string Decimal(double value) {
        char buffer[64];
        auto const [end, error] = std::to_chars(buffer, buffer + sizeof(buffer), value);
        return error == std::errc() ? std::string(buffer, end) : std::string("nan");
    }

    // A finite decimal, or nothing.
    inline std::optional<double> ParseDecimal(std::string_view text) {
        double value = 0;
        auto const [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
        if (text.empty() || error != std::errc() || end != text.data() + text.size() || !std::isfinite(value)) return std::nullopt;
        return value;
    }

    inline std::optional<bool> ParseOnOff(std::string_view text) {
        if (text == "on") return true;
        if (text == "off") return false;
        return std::nullopt;
    }

    inline std::string OnOff(bool on) {
        return on ? "on" : "off";
    }

    // Whether this protal knows mate guidance mode `mode` (on, off).
    inline bool KnownMateGuidance(std::string_view mode) {
        return mode == kMateGuidanceOn || mode == kMateGuidanceOff;
    }

    // The run settings that are said, as build_settings.tsv names them and their values as text, in the file's order:
    // what Write writes, and what runs and databases are compared by (Differences).
    using Fields = std::vector<std::pair<std::string, std::string>>;

    inline Fields RunFields(RunSettings const& s) {
        Fields fields;
        auto add = [&fields](std::string_view key, std::string value) { fields.emplace_back(std::string(key), std::move(value)); };
        if (s.weighted_site_shift) add(kWeightedSiteShift, OnOff(*s.weighted_site_shift));
        if (s.allele_scores) add(kAlleleScores, OnOff(*s.allele_scores));
        if (s.align_top) add(kAlignTop, std::to_string(*s.align_top));
        if (s.adaptive_candidates) add(kAdaptiveCandidates, std::to_string(*s.adaptive_candidates));
        if (s.mate_guidance) add(kMateGuidance, *s.mate_guidance);
        for (auto const& info : kReadTypes) {
            auto const& floor = s.max_score_ani[static_cast<size_t>(info.type)];
            if (floor) add(std::string(kMaxScoreAni) + info.token, Decimal(*floor));
        }
        if (s.unique_kmers) add(kUniqueKmers, *s.unique_kmers);
        if (s.keep_foreign_genes) add(kKeepForeignGenes, OnOff(*s.keep_foreign_genes));
        if (s.keep_suspect_copies) add(kKeepSuspectCopies, OnOff(*s.keep_suspect_copies));
        if (s.depth_identity_margin) add(kDepthIdentityMargin, Decimal(*s.depth_identity_margin));
        return fields;
    }

    // Every run setting's name, in the file's order (max_score_ani_ with each read type).
    inline std::vector<std::string> RunSettingNames() {
        std::vector<std::string> names = { std::string(kWeightedSiteShift), std::string(kAlleleScores), std::string(kAlignTop),
                                           std::string(kAdaptiveCandidates), std::string(kMateGuidance) };
        for (auto const& info : kReadTypes) names.push_back(std::string(kMaxScoreAni) + info.token);
        for (auto name : { kUniqueKmers, kKeepForeignGenes, kKeepSuspectCopies, kDepthIdentityMargin }) names.emplace_back(name);
        return names;
    }

    // The value of `key` in `fields`, or nothing.
    inline std::optional<std::string> FieldValue(Fields const& fields, std::string_view key) {
        for (auto const& [name, value] : fields) {
            if (name == key) return value;
        }
        return std::nullopt;
    }

    // The settings of `used` whose values differ from those `recorded` says (a setting not recorded is not compared), as
    // "name value (the database's: value)", for the keys `compared` admits.
    template<typename Compared>
    std::vector<std::string> Differences(Fields const& recorded, Fields const& used, Compared const& compared) {
        std::vector<std::string> differences;
        for (auto const& [name, value] : used) {
            if (!compared(name)) continue;
            auto const database = FieldValue(recorded, name);
            if (database && *database != value) differences.push_back(name + " " + value + " (the database's: " + *database + ")");
        }
        return differences;
    }

    // A run setting's line into `run`: the error text of a value it cannot be, or "" (also for a name that is no run
    // setting: then `known` is false).
    inline std::string ReadRunSetting(std::string_view name, std::string_view value, RunSettings& run, bool& known) {
        known = true;
        auto on_off = [&](std::optional<bool>& target) -> std::string {
            target = ParseOnOff(value);
            return target ? "" : std::string(name) + " is on or off, not '" + std::string(value) + "'";
        };
        auto count = [&](std::optional<size_t>& target) -> std::string {
            size_t number = 0;
            auto const [end, error] = std::from_chars(value.data(), value.data() + value.size(), number);
            if (value.empty() || error != std::errc() || end != value.data() + value.size()) {
                return std::string(name) + " is a whole number, not '" + std::string(value) + "'";
            }
            target = number;
            return "";
        };
        auto decimal = [&](std::optional<double>& target) -> std::string {
            target = ParseDecimal(value);
            if (!target || *target < 0) {
                target.reset();
                return std::string(name) + " is a decimal of 0 or more, not '" + std::string(value) + "'";
            }
            return "";
        };
        auto text = [&](std::optional<std::string>& target) -> std::string {
            if (value.empty()) return std::string(name) + " has no value";
            target = std::string(value);
            return "";
        };
        if (name == kWeightedSiteShift) return on_off(run.weighted_site_shift);
        if (name == kAlleleScores) return on_off(run.allele_scores);
        if (name == kAlignTop) return count(run.align_top);
        if (name == kAdaptiveCandidates) return count(run.adaptive_candidates);
        if (name == kMateGuidance) return text(run.mate_guidance);
        if (name == kUniqueKmers) return text(run.unique_kmers);
        if (name == kKeepForeignGenes) return on_off(run.keep_foreign_genes);
        if (name == kKeepSuspectCopies) return on_off(run.keep_suspect_copies);
        if (name == kDepthIdentityMargin) return decimal(run.depth_identity_margin);
        if (name.substr(0, kMaxScoreAni.size()) == kMaxScoreAni) {
            auto const type = ReadTypeFromToken(std::string(name.substr(kMaxScoreAni.size())));
            if (!type) {  // a later version's read type
                known = false;
                return "";
            }
            return decimal(run.max_score_ani[static_cast<size_t>(*type)]);
        }
        known = false;
        return "";
    }

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
            bool run_setting = false;
            if (std::string const problem = ReadRunSetting(name, value, settings.run, run_setting); run_setting) {
                if (!problem.empty()) return "line " + std::to_string(number) + ": " + problem;
                continue;
            }
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
        for (auto const& [name, value] : RunFields(settings.run)) os << name << '\t' << value << '\n';
    }
}
