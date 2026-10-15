// protal::Run: the run's start-up (options, the database's gene tables, the tables and models loaded side by side), then
// its stages, each in a unit of its own (RunStages.h).

#include <algorithm>
#include <atomic>
#include <chrono>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>
#include "RunProtal.h"
#include "ProfileStage.h"
#include "Benchmark.h"
#include "BuildInfo.h"
#include "CpuLog.h"
#include "FeatureRanges.h"
#include "FeatureVersions.h"
#include "MemoryUse.h"
#include "RunStatus.h"
#include "ScratchSpace.h"
#include "SeqReader.h"
#include "TaxonStatisticsOutput.h"
#include "ThreadedGzStream.h"
#include "protal_config.h"

namespace protal {

    static const std::string PROTAL_LOGO =
            "                                             ,,  \n"
            "`7MM\"\"\"Mq.                   mm            `7MM  \n"
            "  MM   `MM.                  MM              MM  \n"
            "  MM   ,M9 `7Mb,od8 ,pW\"Wq.mmMMmm  ,6\"Yb.    MM  \n"
            "  MMmmdM9    MM' \"'6W'   `Wb MM   8)   MM    MM  \n"
            "  MM         MM    8M     M8 MM    ,pm9MM    MM  \n"
            "  MM         MM    YA.   ,A9 MM   8M   MM    MM  \n"
            ".JMML.     .JMML.   `Ybmd9'  `Mbmo`Moo9^Yo..JMML.\n\n";



    static void PrintLogo(std::ostream& os = std::cout) {
        os << PROTAL_LOGO << std::endl;
    }

    static void PrintProtalInformation(std::ostream& os = std::cout) {
        // os << "\nOption --output_names has recently been changed to --prefix."<< std::endl;
    }

    // The database's models of kind `kind` (--model_kind) of reads of `type`, by the gate they were trained at, each by its
    // member name: its KindModelFile of the type's model file at the database's gate `own` and of GateModelFile at the
    // others; where it has none, the default file there if that is of the kind (profiler::ModelKindOf: a database whose
    // default models are of that kind, or one built before kinds), which is read to tell (a few MB at most).
    static std::map<size_t, std::string> KindModels(Options const& options, ReadType type, size_t own, std::string const& kind) {
        std::map<size_t, std::string> members;
        std::vector<size_t> gates = { own };
        for (size_t gate : kModelGates) {
            if (gate != own) gates.push_back(gate);
        }
        for (size_t gate : gates) {
            std::string plain = gate == own ? Info(type).model_file : GateModelFile(type, gate);
            std::string const kinded = KindModelFile(plain, kind);
            if (options.DatabaseFileNamed(kinded).Exists()) {
                members[gate] = kinded;
                continue;
            }
            if (gate == own) {  // the type's model, or a paired-end model of a database from before read types
                for (auto const& candidate : ModelCandidates(type)) {
                    if (!options.DatabaseFileNamed(candidate).Exists()) continue;
                    plain = candidate;
                    break;
                }
            }
            auto const file = options.DatabaseFileNamed(plain);
            if (!file.Exists()) continue;
            std::string error;
            auto const xml = file.ReadAll(error);
            if (xml && profiler::ModelKindOf(*xml) == kind) members[gate] = plain;
        }
        return members;
    }

    // Each read type's sparse-first gate (--sparse_first): given, else (auto) as the database records it (off where it
    // does not say), and the gates' score and hits likewise; a map's SPARSE_FIRST sets a sample's own. The reads a gate
    // gives up count no failed candidates, which the models' failed_candidate_rate uses: the database's models were
    // trained on runs gated as it records, and a run gated otherwise (a read type's gate, or a sample's) is warned.
    // --stage1_probe measures what a gate would give up: the database's gates are not used with it. Each gate's pick
    // (--sparse_first_pick: hash or stride) likewise: given, else the database's (build_settings::DatabasePick: stride
    // where it does not say, hash for a database that records its gates but no pick), a map's SPARSE_FIRST_PICK a
    // sample's own; a sample gated with another pick than its database's is gated otherwise. A database may hold models
    // trained at other gates than its own (GateModelFile): each sample is profiled with the model of the gate nearest
    // its own (Options::ModelGate), and is gated otherwise only if that gate is not its own.
    // --sparse_seeding stands for an index of the hash pick's cores: a sample it would seed with the stride pick stops
    // the run (exit 2). With a kind of models asked for (--model_kind), the gates are those the read type has models of
    // that kind at (KindModels); a read type of the run's samples with none stops the run (exit 2).
    // With --profile_only, `sam_gates` holds the gate and pick each sample's SAM was aligned with, where its header says
    // (SamSettings.h; CheckSamSettings): those samples are warned as a run's would be, the others (an older protal's SAMs)
    // are not, as before.
    using SamGates = std::vector<std::optional<std::pair<size_t, SparseFirstPick>>>;
    static void ApplySparseFirst(Options& options, build_settings::Settings const& settings, SamGates const& sam_gates,
                                 std::ostream& out = std::cout) {
        bool const sams_given = ProfilesExistingSams(options);  // their gates, as their headers say (sam_gates)
        std::array<size_t, kReadTypeCount> gates{};
        std::array<SparseFirstPick, kReadTypeCount> picks{};
        std::array<SparseFirstPick, kReadTypeCount> database_picks{};
        std::array<size_t, kReadTypeCount> database_gates{};
        std::array<std::vector<size_t>, kReadTypeCount> model_gates;
        std::array<std::map<size_t, std::string>, kReadTypeCount> kind_models;
        for (size_t t = 0; t < kReadTypeCount; t++) {
            auto const type = static_cast<ReadType>(t);
            auto const& given = options.SparseFirstGiven()[t];
            database_gates[t] = settings.sparse_first[t].value_or(0);
            gates[t] = given ? *given : options.Stage1Probe() ? 0 : database_gates[t];
            database_picks[t] = build_settings::DatabasePick(settings, t);
            auto const& pick_given = options.SparseFirstPickGiven()[t];
            picks[t] = pick_given ? *pick_given : database_picks[t];
            auto const& kind = options.ModelKind(type);
            if (!kind.empty() && options.ModelNamed(type).empty()) {
                // The gates it has models of the kind asked for at.
                kind_models[t] = KindModels(options, type, database_gates[t], kind);
                for (auto const& [gate, member] : kind_models[t]) model_gates[t].push_back(gate);
                if (kind_models[t].empty() && options.AnySample(type) && !options.NoProfile()) {
                    std::cerr << "Error: --model_kind " << kind << ": the database has no " << kind << " model of " << ReadTypeName(type)
                              << " reads (" << KindModelFile(Info(type).model_file, kind) << ", or a " << Info(type).model_file
                              << " of that kind; other gates likewise). build_gtdb_database.py trains both kinds (--model "
                              << "gbm,forest, its default since 2026-10-13); protal --add_model FILE --read_type " << Info(type).token
                              << " --model_kind " << kind << " stores one. Leave out --model_kind for the database's default models."
                              << std::endl;
                    exit(2);
                }
                continue;
            }
            // The gates the read type has models of: the database's own, and those of its models at other gates.
            model_gates[t].push_back(database_gates[t]);
            for (size_t gate : kModelGates) {
                if (gate != database_gates[t] && options.GateModelDbFile(type, gate).Exists()) {
                    model_gates[t].push_back(gate);
                }
            }
            std::sort(model_gates[t].begin(), model_gates[t].end());
        }
        options.SetModelGates(model_gates, database_gates, kind_models);
        size_t const recorded_score = settings.sparse_first_score.value_or(14);
        size_t const recorded_hits = settings.sparse_first_hits.value_or(1);
        options.SetSparseFirst(gates, options.SparseFirstScoreGiven().value_or(recorded_score),
                               options.SparseFirstHitsGiven().value_or(recorded_hits));
        options.SetSparseFirstPick(picks);
        size_t const recorded_min_picks = settings.sparse_first_min_picks.value_or(0);
        options.SetSparseFirstMinPicks(options.SparseFirstMinPicksGiven().value_or(recorded_min_picks));
        if (options.SparseSeeding()) {
            for (auto index : options.GetRange()) {
                if (options.GetSparseFirst(index) <= 1 || options.GetSparseFirstPick(index) != SparseFirstPick::Stride) continue;
                std::cerr << "Error: --sparse_seeding seeds as an index of the hash pick's cores would; sample "
                          << options.GetSampleId(index) << " is gated with the stride pick (--sparse_first_pick, the map's "
                          << "SPARSE_FIRST_PICK or the database's): give --sparse_first_pick hash" << std::endl;
                exit(2);
            }
        }
        bool const recorded_any = std::any_of(settings.sparse_first.begin(), settings.sparse_first.end(),
                                              [](auto const& g) { return g && *g > 0; });
        if (!recorded_any && !options.AnySparseFirst()) return;  // ungated database and run: nothing to say
        // The read types of this run's samples, as the run gates them and as the database does (a gate's pick named
        // where it is stride; the gates of its models at other gates in brackets); the samples gated otherwise than the
        // model they are profiled with was trained (a read type's gate or pick, or the map's SPARSE_FIRST or
        // SPARSE_FIRST_PICK, without a model of that gate).
        std::string run, database;
        size_t others = 0;
        std::string example, example_gates;
        auto const level = [](size_t c, SparseFirstPick pick) {
            return c ? std::to_string(c) + (c > 1 && pick != SparseFirstPick::Hash ? " " + std::string(SparseFirstPickName(pick)) : "")
                     : std::string("off");
        };
        for (auto const& info : kReadTypes) {
            if (!options.AnySample(info.type)) continue;
            size_t const t = static_cast<size_t>(info.type);
            run += (run.empty() ? "" : ", ") + info.token + " " + level(gates[t], picks[t]);
            database += (database.empty() ? "" : ", ") + info.token + " " + level(database_gates[t], database_picks[t]);
            std::string also;
            for (size_t gate : model_gates[t]) {
                if (gate != database_gates[t]) also += (also.empty() ? "" : ", ") + (gate ? std::to_string(gate) : std::string("off"));
            }
            if (!also.empty()) database += " (models also at " + also + ")";
        }
        for (auto index : options.GetRange()) {
            size_t const t = static_cast<size_t>(options.GetReadType(index));
            // --profile_only: the gate its SAM was aligned at, where its header says; else nothing to say, as before.
            auto const& aligned = index < sam_gates.size() ? sam_gates[index] : std::nullopt;
            if (sams_given && !aligned) continue;
            size_t const gate = aligned ? aligned->first : options.GetSparseFirst(index);
            SparseFirstPick const pick = aligned ? aligned->second : options.GetSparseFirstPick(index);
            // At C 1 both picks look up every core: the pick then makes no difference.
            bool const same_pick = gate <= 1 || pick == database_picks[t];
            if (gate == options.ModelGate(index) && same_pick) continue;
            if (others++ == 0) {
                example = options.GetSampleId(index);
                example_gates = example + (sams_given ? " aligned gated at " : " gated at ") + level(gate, pick) +
                                ", its model at " + level(options.ModelGate(index), database_picks[t]);
            }
        }
        bool const score_differs = options.AnySparseFirst() && !sams_given &&
                                   (options.GetSparseFirstScore() != recorded_score || options.GetSparseFirstHits() != recorded_hits ||
                                    options.GetSparseFirstMinPicks() != recorded_min_picks);
        bool const gates_given = options.SparseFirstGiven() != SparseFirstGates();
        bool const picks_given = options.SparseFirstPickGiven() != SparseFirstPicks();
        std::string const given = gates_given && picks_given ? " (--sparse_first, --sparse_first_pick)" :
                                  gates_given ? " (--sparse_first)" : picks_given ? " (--sparse_first_pick)" : "";
        out << "sparse first: " << run << given
            << ", score " << options.GetSparseFirstScore() << ", " << options.GetSparseFirstHits() << " hit(s), least picks "
            << options.GetSparseFirstMinPicks() << "; the database's models were trained with " << database << ", score "
            << recorded_score << ", " << recorded_hits << " hit(s), least picks " << recorded_min_picks << std::endl;
        if (others > 0 || score_differs) {
            std::cerr << "Warning: " << (others > 0 ? std::to_string(others) + " sample(s) of this run (e.g. " + example +
                                                      (sams_given ? ") were aligned gated otherwise (their SAM "
                                                                               "headers say)"
                                                                             : ") are gated otherwise (--sparse_first, "
                                                                               "--sparse_first_pick or the map's SPARSE_FIRST "
                                                                               "or SPARSE_FIRST_PICK)")
                                                    : std::string("this run's sparse-first gates score otherwise"))
                      << " than the database's models were trained with" << (others > 0 ? " (" + example_gates + ")" : "")
                      << ": the reads a gate gives up count no failed "
                      << "candidates, so their calls are not those of training. "
                      << (sams_given ? "Align them again on this database with --sparse_first and --sparse_first_pick "
                                                  "left out (auto), or at a gate the database has a model of."
                                                : "Leave out --sparse_first and --sparse_first_pick (auto) to gate as they were "
                                                  "trained, or gate at a gate the database has a model of.") << std::endl;
        }
    }

    // A run setting's place (RunSettings, by its build_settings.tsv name): the alignment's, or the profiling's.
    static bool ProfilingSetting(std::string_view name) {
        return name == build_settings::kKeepForeignGenes || name == build_settings::kKeepSuspectCopies ||
               name == build_settings::kDepthIdentityMargin;
    }

    // Whether run setting `name` changes anything for this run's samples: the profiling's unless --no_profile; the
    // alignment's unless the run profiles SAMs it does not align (`sams_given`: --profile_only or reused SAMs, whose own
    // are checked, CheckSamSettings), each for the read types it acts on
    // (the allele scores, the site shift, the adaptive candidates and the gene end clips short reads, mate guidance
    // paired-end reads, a read type's identity floor its reads).
    static bool SettingApplies(Options const& options, std::string const& name, bool sams_given) {
        namespace b = build_settings;
        if (ProfilingSetting(name)) return !options.NoProfile();
        if (sams_given) return false;
        bool const short_reads = options.AnySample(ReadType::Paired) || options.AnySample(ReadType::Single);
        if (name == b::kAlleleScores || name == b::kWeightedSiteShift || name == b::kAdaptiveCandidates ||
            name == b::kGeneEndClips) {
            return short_reads;
        }
        if (name == b::kMateGuidance) return options.AnySample(ReadType::Paired);
        if (name.starts_with(b::kMaxScoreAni)) {
            auto const type = ReadTypeFromToken(name.substr(b::kMaxScoreAni.size()));
            return type && options.AnySample(*type);
        }
        return std::any_of(kReadTypes.begin(), kReadTypes.end(), [&options](ReadTypeInfo const& t) { return options.AnySample(t.type); });
    }

    // The unique k-mers' scheme the database's models were trained with: as build_settings.tsv records it, else cores if it
    // has unique_cores.bin (since 2026-10-14 every build writes it unless --no_unique_cores), else nothing (anchors, as
    // before).
    static std::optional<std::string> DatabaseUniqueKmers(Options const& options, build_settings::Settings const& settings) {
        if (settings.run.unique_kmers) return settings.run.unique_kmers;
        if (options.UniqueCoresDbFile().Exists()) return std::string(build_settings::kUniqueCores);
        return std::nullopt;
    }

    // The run settings (BuildSettings.h RunSettings, recorded since 2026-10-15): each not given (auto) as the database
    // records it, else protal's default (Options::FollowDatabaseRunSettings). Prints each with the database's; warns once
    // for those this run is told otherwise than the database's models were trained with (SettingApplies), and apart for a
    // run with --no_unique_cores on a database whose models counted the unique k-mers from its unique cores. A mode the
    // database records that this protal does not know (a later protal's mate guidance) is warned and left at protal's
    // default. A database built before records none: one note, and protal's defaults.
    static void ApplyRunSettings(Options& options, build_settings::Settings const& settings, std::ostream& out) {
        namespace b = build_settings;
        bool const sams_given = ProfilesExistingSams(options);
        auto const& recorded = settings.run;
        if (recorded.mate_guidance && !b::KnownMateGuidance(*recorded.mate_guidance)) {
            std::cerr << "Warning: the database records mate_guidance '" << *recorded.mate_guidance << "' (" << b::kFileName
                      << "), a mode this protal does not know (a later protal's?): paired-end reads are settled by protal's "
                      << "default (" << MateGuidanceName(kDefaultMateGuidance) << "), which the database's models may not have been trained with. "
                      << "Profile them with the protal that built the database." << std::endl;
        }
        if (recorded.gene_end_clips && !b::KnownGeneEndClips(*recorded.gene_end_clips)) {
            std::cerr << "Warning: the database records gene_end_clips '" << *recorded.gene_end_clips << "' (" << b::kFileName
                      << "), a value this protal does not know (a later protal's?): short reads are ranked with "
                      << "--gene_end_clips off, which the database's models may not have been trained with. Profile them "
                      << "with the protal that built the database." << std::endl;
        }
        if (recorded.unique_kmers && *recorded.unique_kmers != b::kUniqueCores && *recorded.unique_kmers != b::kUniqueAnchors) {
            std::cerr << "Warning: the database records unique_kmers '" << *recorded.unique_kmers << "' (" << b::kFileName
                      << "), which this protal does not know: it counts ZS, ZU and ZT "
                      << (options.UniqueCoresDbFile().Exists() && !options.NoUniqueCores() ? "from the cores each alignment covers"
                                                                                         : "from each record's anchor's seeds")
                      << ", as it would on any database." << std::endl;
        }
        options.FollowDatabaseRunSettings(recorded);
        auto used = options.RunSettingsInUse();
        used.unique_kmers = std::string(options.UniqueCoresDbFile().Exists() && !options.NoUniqueCores() ? b::kUniqueCores : b::kUniqueAnchors);
        auto database = recorded;
        database.unique_kmers = DatabaseUniqueKmers(options, settings);
        // A database that records no gene_end_clips was trained with off (every protal before ranked so): no note, and a
        // run given credit is warned.
        database.gene_end_clips = b::DatabaseGeneEndClips(recorded);
        auto const used_fields = b::RunFields(used);
        auto const database_fields = b::RunFields(database);
        auto const given = b::RunFields(options.RunSettingsGiven());

        // The block: each setting, this run's value (given or auto) and the database's.
        out << "Run settings (" << b::kFileName << "; a setting left at auto takes the database's, else protal's default):\n"
            << "  " << std::left << std::setw(26) << "setting" << std::setw(28) << "this run" << "the database's models" << '\n';
        std::vector<std::string> unrecorded;
        for (auto const& [name, value] : used_fields) {
            auto const db_value = b::FieldValue(database_fields, name);
            bool const alignment_only = sams_given && !ProfilingSetting(name);
            std::string const run_value = alignment_only ? std::string("(the SAMs')")
                                                        : value + (b::FieldValue(given, name) ? " (given)" : " (auto)");
            out << "  " << std::setw(26) << name << std::setw(28) << run_value << db_value.value_or("not recorded") << '\n';
            if (!db_value && SettingApplies(options, name, sams_given)) unrecorded.push_back(name);
        }
        out << std::right << std::flush;
        if (!unrecorded.empty()) {
            std::string list;
            for (auto const& name : unrecorded) list += (list.empty() ? "" : ", ") + name;
            out << "Note: the database records no " << list << " (built before 2026-10-15): this run uses protal's defaults "
                << "for them (or the options given), which may not be what its models were trained with"
                << (std::find(unrecorded.begin(), unrecorded.end(), std::string(b::kMateGuidance)) != unrecorded.end() &&
                    !options.RunSettingsGiven().mate_guidance
                        ? "; mate guidance: guide, the mode of every protal before (its models were trained with it)" : "") << std::endl;
        }

        // --no_unique_cores on a database whose models were trained on the unique cores.
        bool const cores_unread = !sams_given && options.NoUniqueCores() &&
                                  database.unique_kmers == std::string(b::kUniqueCores);
        if (cores_unread) {
            std::cerr << "Warning: this run counts each record's ZS, ZU and ZT from its anchor's seeds (--no_unique_cores), but "
                      << "the database's models were trained on runs that counted them from the cores each alignment covers "
                      << "(unique_kmers cores: " << Options::PROTAL_UNIQUE_CORES_FILE << "): their unique k-mer features (lu, "
                      << "lsu, lu_per_kb, lsu_per_kb, the lu and lsu gene rates, unique_seed_share, unique_seed_enrichment, "
                      << "...) differ from those of training, and so do their calls. Leave out --no_unique_cores." << std::endl;
        }
        bool const mode_unknown = recorded.mate_guidance && !b::KnownMateGuidance(*recorded.mate_guidance);  // warned above
        bool const clips_unknown = recorded.gene_end_clips && !b::KnownGeneEndClips(*recorded.gene_end_clips);  // warned above
        auto const differences = b::Differences(database_fields, used_fields, [&](std::string const& name) {
            return name != b::kUniqueKmers && !(name == b::kMateGuidance && mode_unknown) &&
                   !(name == b::kGeneEndClips && clips_unknown) && SettingApplies(options, name, sams_given);
        });
        if (!differences.empty()) {
            std::string list;
            for (auto const& d : differences) list += (list.empty() ? "" : ", ") + d;
            std::cerr << "Warning: this run is told otherwise than the database's models were trained (" << b::kFileName << "): "
                      << list << ". Its alignments and features, and so its calls, are not those of training. Leave these "
                      << "options out (auto) to run as the models were trained." << std::endl;
        }
    }

    // --profile_only: each SAM's settings line (SamSettings.h). Its gate and pick choose the sample's model, as a map's
    // SPARSE_FIRST would (Options::SetSampleGate), unless --sparse_first (of its read type) or the map's SPARSE_FIRST gives
    // the sample's own; and `sam_gates` gets them for ApplySparseFirst's warning. One warning for the SAMs aligned
    // otherwise than the database's models were trained (the traceback, the gates' score and hits, the unique k-mers'
    // scheme, the allele scores, the site shift, mate guidance, --align_top, --adaptive_candidates, --gene_end_clips, the
    // identity floor), with the first one's differences. A SAM without the line (an older protal's) is profiled as before.
    static SamGates CheckSamSettings(Options& options, build_settings::Settings const& settings, std::ostream& out) {
        namespace b = build_settings;
        namespace s = sam_settings;
        SamGates sam_gates;
        if (!ProfilesExistingSams(options)) return sam_gates;
        // What the database's models were trained with, as the SAM line names it (only what it records).
        auto expected = [&](ReadType type) {
            s::Fields fields;
            fields.emplace_back(s::kTraceback, std::string(WfaTracebackName(ResolveWfaTraceback(std::nullopt, settings.wfa_traceback))));
            if (auto const unique = DatabaseUniqueKmers(options, settings)) fields.emplace_back(s::kUniqueKmers, *unique);
            auto const run = b::RunFields(settings.run);
            auto copy = [&](std::string_view sam_key, std::string const& name) {
                if (auto value = b::FieldValue(run, name)) fields.emplace_back(std::string(sam_key), *value);
            };
            copy(s::kAlleleScores, std::string(b::kAlleleScores));
            copy(s::kWeightedSiteShift, std::string(b::kWeightedSiteShift));
            copy(s::kMateGuidance, std::string(b::kMateGuidance));
            copy(s::kAlignTop, std::string(b::kAlignTop));
            copy(s::kAdaptiveCandidates, std::string(b::kAdaptiveCandidates));
            // (only short reads' SAMs say it; a database that records none was trained with off)
            fields.emplace_back(s::kGeneEndClips, b::DatabaseGeneEndClips(settings.run));
            copy(s::kMaxScoreAni, std::string(b::kMaxScoreAni) + Info(type).token);
            return fields;
        };
        size_t lines = 0, otherwise = 0;
        std::string example;
        for (auto index : options.GetRange()) {
            auto const sam = options.SamFile(index).first;
            if (!std::filesystem::exists(sam)) continue;  // reported when it is profiled
            std::optional<s::Fields> fields;
            try {
                SamInput input(sam);
                if (input.IsOpen()) fields = s::ReadSamSettings(input.Stream());
            } catch (std::exception const&) {}  // reported when it is profiled
            if (!fields) continue;
            lines++;
            if (index >= sam_gates.size()) sam_gates.resize(index + 1);
            auto const gate = s::Get(*fields, s::kGate);
            auto const level = gate ? ParseSparseFirstLevel(*gate) : std::nullopt;
            auto const pick_text = s::Get(*fields, s::kPick);
            auto const pick = pick_text ? ParseSparseFirstPick(*pick_text) : std::nullopt;
            auto const type = options.GetReadType(index);
            if (level) {
                sam_gates[index] = std::make_pair(*level, pick.value_or(SparseFirstPick::Stride));
                options.SetSampleGate(index, *level, pick);  // unless --sparse_first or the map's SPARSE_FIRST says
            }
            // Of its other settings, those that differ from the database's (its own read type's floor; score and hits only
            // where it was gated).
            auto const sam_type = s::Get(*fields, s::kReadType);
            auto const recorded = expected(sam_type ? ReadTypeFromToken(*sam_type).value_or(type) : type);
            auto with_gate = recorded;
            if (level && *level > 0) {
                with_gate.emplace_back(s::kScore, std::to_string(settings.sparse_first_score.value_or(14)));
                with_gate.emplace_back(s::kHits, std::to_string(settings.sparse_first_hits.value_or(1)));
                if (settings.sparse_first_min_picks) with_gate.emplace_back(s::kMinPicks, std::to_string(*settings.sparse_first_min_picks));
            }
            auto const differences = b::Differences(with_gate, *fields, [](std::string const&) { return true; });
            if (differences.empty()) continue;
            if (otherwise++ == 0) {
                example = sam + ": ";
                for (size_t d = 0; d < differences.size(); d++) example += (d ? ", " : "") + differences[d];
            }
        }
        if (lines > 0) {
            out << "SAM settings: " << lines << " of " << options.GetRange().size() << " SAM file(s) say how they were "
                << "aligned (@CO protal settings); their gates choose their models" << std::endl;
        }
        if (otherwise > 0) {
            std::cerr << "Warning: " << otherwise << " SAM file(s) of this run were aligned otherwise than the database's "
                      << "models were trained (e.g. " << example << "): their features, and so their calls, are not those "
                      << "of training. Align their reads again on this database, with these options left out (auto)."
                      << std::endl;
        }
        return sam_gates;
    }

    // The database's build_settings.tsv (none: every setting unsaid); exits 8 if it cannot be read.
    build_settings::Settings DatabaseBuildSettings(Options const& options) {
        auto const file = options.BuildSettingsDbFile();
        build_settings::Settings settings;
        if (file.Exists()) {
            std::string error;
            auto const content = file.ReadAll(error);
            if (content) {
                std::istringstream is(*content);
                error = build_settings::Read(is, settings);
            }
            if (!error.empty()) {
                std::cerr << "Invalid build settings " << file.Name() << ": " << error << std::endl;
                exit(8);
            }
        }
        return settings;
    }

    // The WFA2 traceback a run aligns with (WfaTraceback.h): --wfa_traceback, else (auto) the one the database was built
    // with (build_settings.tsv; classic for a database built before 2026-10-11). Warns when the run aligns otherwise than
    // the database was built: its strain alleles, congener gaps and column weights, and the runs its models were trained
    // on, chose among equally good alignments the other way. --build records its own (WriteBuildSettings). Also how the
    // model's `fragments` feature counts (--fragment_count, else the database's), and each read type's sparse-first gate
    // (--sparse_first, else the database's). Exits 8 if build_settings.tsv cannot be read.
    static void ApplyBuildSettings(Options& options, std::ostream& out = std::cout) {
        if (options.BuildMode()) return;
        auto const settings = DatabaseBuildSettings(options);
        auto const given = options.WfaTracebackGiven();
        WfaTraceback const database = ResolveWfaTraceback(std::nullopt, settings.wfa_traceback);
        WfaTraceback const traceback = ResolveWfaTraceback(given, settings.wfa_traceback);
        options.SetWfaTraceback(traceback);
        std::string const built = std::string(WfaTracebackName(database)) +
                                  (settings.wfa_traceback ? "" : " (it has no " + Options::PROTAL_BUILD_SETTINGS_FILE + ": built before 2026-10-11)");
        out << "WFA2 traceback: " << WfaTracebackName(traceback) << (given ? " (--wfa_traceback)" : "")
            << "; the database was built with " << built << std::endl;
        if (traceback != database && !ProfilesExistingSams(options)) {
            std::cerr << "Warning: this run aligns with the " << WfaTracebackName(traceback) << " WFA2 traceback (--wfa_traceback "
                      << WfaTracebackName(traceback) << "), but the database was built with " << built << ". The two give the "
                      << "same alignment scores, but among equally good alignments they may choose differently, so the SAM "
                      << "records, and the features the database's models were trained on, can differ slightly from a run "
                      << "that aligns as the database was built. Leave out --wfa_traceback (auto) to align as it was built."
                      << std::endl;
        }
        // How the model's `fragments` feature counts (--fragment_count): given, else (auto) as the database records it
        // (mapq when it does not say). Its models were trained on runs that counted so.
        FragmentCount const recorded = settings.fragment_count.value_or(FragmentCount::Mapq);
        FragmentCount const count = options.FragmentCountGiven().value_or(recorded);
        options.SetFragmentCount(count);
        out << "fragment count: " << FragmentCountName(count) << (options.FragmentCountGiven() ? " (--fragment_count)" : "")
            << "; the database's models were trained with " << FragmentCountName(recorded) << std::endl;
        if (count != recorded) {
            std::cerr << "Warning: this run counts the model's fragments " << FragmentCountName(count) << " (--fragment_count), "
                      << "but the database's models were trained with " << FragmentCountName(recorded) << " counts: their "
                      << "calls are not those of training. Leave out --fragment_count (auto) to count as they were trained."
                      << std::endl;
        }
        ApplyRunSettings(options, settings, out);
        auto const sam_gates = CheckSamSettings(options, settings, out);
        ApplySparseFirst(options, settings, sam_gates, out);
    }

    static size_t GetTotalSystemMemory()
    {
        long pages = sysconf(_SC_PHYS_PAGES);
        long page_size = sysconf(_SC_PAGE_SIZE);
        return pages * page_size;
    }

    // The versions of the model's features' definitions it records (FeatureVersions.h; the trainer's
    // protal_feature_versions) against this protal's: one warning naming those of other definitions, or one note for a
    // model that records none. And its features' ranges in training (FeatureRanges.h; protal_feature_ranges), which the
    // profiling checks the taxa near the knob against; malformed ones are warned about and left unused.
    static void CheckModelFeatures(db::DbFile const& file, std::string const& xml, profiler::TaxonFilterObj& model) {
        if (profiler::IsPlaceholderModel(xml)) return;  // it scores every taxon 0 whatever its features
        auto const active = feature_ranges::ModelActiveFeatures(xml);
        if (auto const versions = profiler::HeaderExtension(xml, feature_versions::kExtension); !versions) {
            std::cout << "Note: the model " << file.Name() << " does not record its features' versions ("
                      << feature_versions::kExtension << ", which the trainer writes since 2026-10-15): protal cannot "
                      << "tell whether it was trained on this protal's definitions of its features (protal --feature_versions)"
                      << std::endl;
        } else {
            feature_versions::Recorded recorded;
            if (auto const error = feature_versions::Parse(*versions, recorded); !error.empty()) {
                std::cerr << "Warning: the model " << file.Name() << ": " << error << "; its features' versions are not checked"
                          << std::endl;
            } else if (auto const warning = feature_versions::Mismatch(file.Name(), recorded, active); !warning.empty()) {
                std::cerr << "Warning: " << warning << std::endl;
            }
        }
        if (auto const text = profiler::HeaderExtension(xml, feature_ranges::kExtension)) {
            feature_ranges::Ranges ranges;
            if (auto const error = feature_ranges::Parse(*text, ranges); !error.empty()) {
                std::cerr << "Warning: the model " << file.Name() << ": " << error << "; its taxa are not checked against "
                          << "its training range" << std::endl;
            } else {
                model.SetFeatureRanges(std::make_shared<feature_ranges::Ranges const>(feature_ranges::Restricted(ranges, active)));
            }
        }
    }

    // Loads the PMML presence model and its depth knobs (ParseDepthKnobCurve, or an older model's ParseDepthKnobs) and
    // checks that protal can use it (ModelContractProblemInXml); exits 2 if not. Its kind into `kind` if given
    // (profiler::ModelKindOf: gbm, forest, placeholder, or "").
    static profiler::TaxonFilterObj LoadModel(db::DbFile const& file, double knob, std::string* kind = nullptr) {
        std::string read_error;
        auto const xml = file.ReadAll(read_error);
        if (!xml) {
            std::cerr << "Cannot load the model " << file.Name() << ": " << read_error << std::endl;
            exit(2);
        }
        if (kind) *kind = profiler::ModelKindOf(*xml);
        std::optional<profiler::TaxonFilterObj> model;
        try {
            model.emplace(cpmml::Model::from_string(*xml), knob);
        } catch (std::exception const& e) {
            std::cerr << "Cannot load the model " << file.Name() << ": " << e.what() << std::endl;
            exit(2);
        }
        std::map<int, double> depth_knobs;
        profiler::DepthKnobCurve depth_knob_curve;
        profiler::FalseCallModel false_calls;
        auto problem = profiler::ModelContractProblemInXml(model.value(), *xml);
        if (problem.empty()) problem = profiler::ParseDepthKnobs(*xml, depth_knobs);
        if (problem.empty()) problem = profiler::ParseDepthKnobCurve(*xml, depth_knob_curve);
        if (problem.empty()) problem = profiler::ParseFalseCalls(*xml, false_calls);
        if (problem.empty() && !depth_knobs.empty() && !depth_knob_curve.empty()) {
            problem = "it has depth knobs twice, as bins (" + std::string(profiler::kDepthKnobsExtension) + ") and as a curve (" +
                      std::string(profiler::kDepthKnobCurveExtension) + ")";
        }
        if (!problem.empty()) {
            std::cerr << "Cannot use the model " << file.Name() << ": " << problem << std::endl;
            exit(2);
        }
        model->SetDepthKnobs(std::move(depth_knobs));
        model->SetDepthKnobCurve(std::move(depth_knob_curve));
        model->SetFalseCalls(std::move(false_calls));
        CheckModelFeatures(file, *xml, *model);
        if (profiler::IsPlaceholderModel(*xml)) {
            std::cerr << "WARNING: " << file.Name() << " is a placeholder, not a trained model: it scores every taxon 0, so "
                      << "no species is reported (--knob 0 lists every taxon with reads). Train a model for this read type "
                      << "and store it with protal --add_model FILE --read_type TYPE --db DB." << std::endl;
        }
        return *model;
    }

    std::string ModelKindOfFile(db::DbFile const& file, double knob) {
        std::string kind;
        LoadModel(file, knob, &kind);
        return kind;
    }

    // Runs protal and returns the process exit code: 0, or 1 if any sample or output failed (see
    // RunStatus). Invalid input and fatal errors still exit directly with their own codes.
    // Prints how long freeing the run's memory took: the destructors of the database (14.5M genes and their tables at GTDB
    // size), the index and the profiles run after "Run protal took". Declared before them in Run, so destroyed after them;
    // armed when the run is done.
    struct TeardownTimer {
        std::chrono::steady_clock::time_point start;
        bool armed = false;
        void Arm() { start = std::chrono::steady_clock::now(); armed = true; }
        ~TeardownTimer() {
            if (!armed) return;
            auto const ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
            std::cout << "Freeing memory took " << ms / 1000 << "s " << ms % 1000 << "ms" << std::endl;
        }
    };

    int Run(int argc, char *argv[]) {
        auto options = protal::Options::OptionsFromArguments(argc, argv);
        if (options.ShowVersion()) {
            std::cout << "protal " << protal::VersionText() << std::endl;
            exit(0);
        }
    

        PrintLogo();
        PrintProtalInformation();
        std::cout << "Total available memory is " << GetTotalSystemMemory() / (1024 * 1024 * 1024) << "GB" << std::endl;
        std::cout << std::endl;

        Benchmark bm_total("Run protal");
        bm_total.Start();



        if (options.Help()) {
            options.PrintHelp();
            exit(0);
        }

        if (options.HelpDev()) {
            options.PrintHelp(true);
            exit(0);
        }

        if (options.ShowMapHelp()) {
            options.PrintMapHelp();
            exit(0);
        }

        std::cout << "Options:\n" << options.ToString() << std::endl;

        if (options.CompressDbMode() || options.DecompressDbMode() || options.UnpackDbMode() || !options.GetAddModel().empty() ||
            !options.AddTables().empty() || !options.WriteSpeciesNeighboursFile().empty()) {
            return RunDbMaintenance(options);  // each exits 8 on failure
        }

        ApplyBuildSettings(options);

        // Directories that outputs are written into without creating them first (per-species strain
        // files, misc statistics). Created here so that a run started with -1/-2/-o, not a map, has them.
        if (!options.BuildMode()) {
            std::vector<std::string> dirs = { options.GetMiscOutputDir() };
            if (!options.NoStrains()) dirs.push_back(options.GetStrainOutputDir());
            for (auto const& dir : dirs) {
                std::error_code ec;
                if (!dir.empty()) std::filesystem::create_directories(dir, ec);
                if (ec) {
                    std::cerr << "Cannot create output directory " << dir << ": " << ec.message() << std::endl;
                    exit(2);
                }
            }
        }

        // --scratch_dir: the run's folder for intermediate files, removed when the process ends (scratch::Space), and the
        // SAMs of a run without strain MSAs in it (Options::UseScratchFolder).
        if (!options.BuildMode() && !options.ScratchDir().empty()) {
            std::string error;
            std::string const folder = scratch::Space::Get().Make(options.ScratchDir(), error);
            if (folder.empty()) {
                std::cerr << "Error: --scratch_dir: " << error << std::endl;
                exit(2);
            }
            size_t const sams = options.UseScratchFolder(folder);
            std::cout << "Scratch folder: " << folder << " (removed at the end)";
            if (sams > 0) std::cout << "; the SAM files of " << sams << " sample(s) go there and are removed once profiled (--no_strains)";
            std::cout << std::endl;
        }

        // Load protal DB into RAM: reference.fna, reference.fna.zst, or the single-file database's
        // (--build always reads the folder's files).
        auto const unique_kmers_file = options.UniqueKmersDbFile();
        int const db_threads = static_cast<int>(options.GetThreads());
        // The start-up's parts are timed like the stages (docs/claude/2026-10-04-performance-gtdb-scale: ~8 s per run
        // at r226 size were untimed): the gene tables here, the taxonomy, models and tables below, and freeing the memory
        // at the end (TeardownTimer, which prints after "Run protal took").
        TeardownTimer teardown;
        Benchmark bm_gene_tables("Loading the gene tables");
        bm_gene_tables.Start();
        // A single-file database's binary gene table (GeneTableFile.h) where it has a current one, else the text tables.
        // --build makes unique_kmers.tsv anew: one an earlier build left in the folder is not read (it would make the
        // genes it lists the only hittable ones, of which species_neighbours.tsv compares all).
        auto const unique_kmers = unique_kmers_file.Exists() && !options.BuildMode() ? std::optional(unique_kmers_file) : std::nullopt;
        auto const gene_table = options.GeneTableDbFile();
        ProtalDB db(options.SequenceDbFile(), options.SequenceMapDbFile(), unique_kmers, gene_table, db_threads);
        bm_gene_tables.Stop();
        bm_gene_tables.PrintResults();
        std::cout << db.GetGenomes().GeneTableTimes() << std::endl;

        // Skip alignment if files are present. Do not skip if either files are not there or user specified --force
        auto sam_files = options.SamFiles();
        bool all_alignments_exist = std::all_of(sam_files.begin(), sam_files.end(), [](std::string const& file){ return Utils::exists(file); });

        auto const range = options.GetRange();
        bool const piped = !options.BuildMode() && !options.ProfileOnly() &&
                           std::any_of(range.begin(), range.end(),
                                       [&options](auto index) { return ReadsArePiped(options, static_cast<int>(index)); });
        bool skip_alignment = !options.BuildMode() && all_alignments_exist && !options.Force() && !piped;

        bool run_alignment = true;
        if (!options.BuildMode() && (options.ProfileOnly() || skip_alignment) && !sam_files.empty() &&
            !options.SamFile(0).first.empty()) {
            std::cout << "All alignments are present." << std::endl;
            run_alignment = false;
        }
        bool const run_profiling = !options.BuildMode() && !options.NoProfile();
        bool const need_neighbours = !options.BuildMode() && (run_alignment || run_profiling);
        // Before any thread reads the database (RequireNcbiTaxids).
        if (run_profiling && options.NcbiTaxids()) RequireNcbiTaxids(options);

        // The rest of the database. By default the index loads on a thread of its own, beside the genome preload (all
        // threads each) and the single-threaded tables, each on a thread of its own: the taxonomy, the models, gene
        // conservation, suspect copies, species priors and the gene neighbours' table. --sequential_load loads them one
        // after another, as before. Their messages are printed in the order they always had; only the timers differ.
        bool const concurrent = !options.SequentialLoad() && !options.BuildMode();
        // The grid pick's needs, before the index load starts on a thread of its own (CheckGridIndex).
        if (run_alignment && !options.BuildMode()) CheckGridIndex(options);
        if (concurrent && run_alignment) StartIndexLoading(options, db);

        // The models the samples are profiled with (Options::SampleModelDbFile: their read types', trained at the gates
        // nearest theirs), each loaded once for all its samples; without samples, the paired-end model.
        bool const any_sample = std::any_of(kReadTypes.begin(), kReadTypes.end(), [&options](ReadTypeInfo const& t) { return options.AnySample(t.type); });
        struct ModelToLoad {
            db::DbFile file;
            ReadType type;
            size_t gate;
        };
        std::vector<ModelToLoad> model_files;
        if (run_profiling) {
            auto const add = [&model_files](db::DbFile file, ReadType type, size_t gate) {
                bool const listed = std::any_of(model_files.begin(), model_files.end(),
                                                [&](ModelToLoad const& m) { return ModelKey(m.type, m.file) == ModelKey(type, file); });
                if (!listed) model_files.push_back({std::move(file), type, gate});
            };
            for (auto index : options.GetRange()) {
                add(options.SampleModelDbFile(index), options.GetReadType(index), options.ModelGate(index));
            }
            if (!any_sample) add(options.ModelDbFile(ReadType::Paired), ReadType::Paired, 0);
        }
        std::vector<std::string> model_kinds(model_files.size());  // each one's kind (profiler::ModelKindOf), once loaded
        auto load_models = [&options, &model_files, &model_kinds]() {
            std::vector<profiler::TaxonFilterObj> loaded;
            for (size_t m = 0; m < model_files.size(); m++) loaded.push_back(LoadModel(model_files[m].file, options.GetKnob(), &model_kinds[m]));
            return loaded;
        };
        auto preload = [&options, &db]() {
            if (!options.PreloadGenomes()) return;
            std::cout << "Preload genomes" << std::endl;
            // --build's index passes take the genes from here (GenomeLoader::KeepReferenceText) rather than read reference.fna.
            if (options.BuildMode()) db.GetGenomes().KeepReferenceText(true);
            Benchmark bm_preload_genomes("Preload genomes");
            bm_preload_genomes.Start();
            db.GetGenomes().LoadAllGenomes(static_cast<int>(options.GetThreads()));
            bm_preload_genomes.Stop();
            bm_preload_genomes.PrintResults();
            std::cout << db.GetGenomes().PreloadTimes() << std::endl;
            protal::build::PrintMemory("the genome preload");
        };

        Benchmark bm_tables("Loading the taxonomy, models and tables");
        std::vector<profiler::TaxonFilterObj> loaded_models;
        std::ostringstream conservation_log, suspect_log, priors_log, ncbi_log, species_neighbours_log, gaps_log, foreign_log, weights_log,
                sites_log, neighbours_log;
        std::ostringstream alleles_log;
        // The alignment reads the strain alleles and the column weights only to score short reads' candidates with the
        // alleles (ScoreAlleles; not with --no_allele_scores), the column weights for its weighted site shift
        // (--weighted_site_shift, off by default): then they are read before the alignment, with or without profiling. Else
        // only profiling reads them, and a run without profiling leaves them unread.
        bool const alleles_for_alignment = run_alignment && !options.BuildMode() && !options.NoAlleleScores() &&
                                           (options.AnySample(ReadType::Paired) || options.AnySample(ReadType::Single));
        bool const weights_for_alignment = alleles_for_alignment && options.WeightedSiteShift() && !options.NoSiteWeights();
        bool const need_alleles = alleles_for_alignment || run_profiling;
        bool const cores_for_alignment = run_alignment && !options.BuildMode();  // the unique cores (LoadUniqueCores)
        std::ostringstream cores_log;
        std::optional<gene_neighbours::Table> neighbours;
        // In a concurrent start-up that aligns, the tables only profiling reads load beside the alignment, and the first
        // sample to be profiled waits for them (DeferredTables); the alignment waits only for what it reads: the gene
        // neighbours, the taxonomy (genera of the adaptive candidates) and the models (checked first), and the strain alleles
        // and column weights where alleles_for_alignment and weights_for_alignment. Declared after the logs, which its loads
        // write to, so that it waits for them before the logs go.
        DeferredTables deferred;
        if (concurrent) {
            bm_tables.Start();
            std::vector<std::future<void>> tables;
            auto& later = run_alignment ? deferred.futures : tables;
            std::future<std::vector<profiler::TaxonFilterObj>> models_loading;
            std::future<std::optional<gene_neighbours::Table>> neighbours_loading;
            deferred.Start();
            auto load_later = [&](auto load) { later.push_back(std::async(std::launch::async, deferred.Timed(std::move(load)))); };
            if (run_profiling) {
                tables.push_back(std::async(std::launch::async, [&]() { db.LoadTaxonomy(options.TaxonomyDbFile()); }));
                models_loading = std::async(std::launch::async, load_models);
                load_later([&]() { LoadGeneConservation(options, db.GetGenomes(), conservation_log); });
                load_later([&]() { LoadSuspectCopies(options, db.GetGenomes(), suspect_log); });
                load_later([&]() { LoadSpeciesPriors(options, db.GetGenomes(), priors_log); });
                if (options.NcbiTaxids()) load_later([&]() { LoadNcbiTaxids(options, db.GetGenomes(), ncbi_log); });
                load_later([&]() { LoadSpeciesNeighbours(options, db.GetGenomes(), species_neighbours_log); });
                load_later([&]() { LoadCongenerGaps(options, db.GetGenomes(), gaps_log); });
                load_later([&]() { LoadForeignRates(options, db.GetGenomes(), foreign_log); });
                load_later([&]() { LoadAncestrySites(options, db.GetGenomes(), sites_log); });
            }
            auto load_weights = [&]() { LoadColumnWeights(options, db.GetGenomes(), weights_log); };
            if (weights_for_alignment) tables.push_back(std::async(std::launch::async, load_weights));
            else if (run_profiling) load_later(load_weights);
            if (need_neighbours) neighbours_loading = std::async(std::launch::async, [&]() { return ReadGeneNeighbours(options, neighbours_log); });
            auto load_alleles = [&]() { LoadStrainAlleles(options, db.GetGenomes(), alleles_log); };
            if (alleles_for_alignment) tables.push_back(std::async(std::launch::async, load_alleles));
            else if (run_profiling) load_later(load_alleles);
            if (cores_for_alignment) tables.push_back(std::async(std::launch::async, [&]() { LoadUniqueCores(options, db.GetGenomes(), cores_log); }));
            if (!deferred.futures.empty()) {
                deferred.messages = [&, weights_later = !weights_for_alignment, alleles_later = !alleles_for_alignment]() {
                    return conservation_log.str() + suspect_log.str() + priors_log.str() + ncbi_log.str() + species_neighbours_log.str() +
                           gaps_log.str() + foreign_log.str() + (weights_later ? weights_log.str() : std::string()) + sites_log.str() +
                           (alleles_later ? alleles_log.str() : std::string());
                };
            }
            // The tables set their own parts of the genome loader, never its genes, which the preload fills meanwhile; those
            // loading beside the alignment set parts the alignment does not read.
            preload();
            for (auto& t : tables) t.get();
            if (models_loading.valid()) loaded_models = models_loading.get();
            if (neighbours_loading.valid()) neighbours = neighbours_loading.get();
        } else {
            preload();
            bm_tables.Start();
            if (cores_for_alignment) LoadUniqueCores(options, db.GetGenomes(), cores_log);
            if (need_alleles) LoadStrainAlleles(options, db.GetGenomes());
            if (run_profiling) {
                db.LoadTaxonomy(options.TaxonomyDbFile());
                loaded_models = load_models();
            } else if (weights_for_alignment) {
                LoadColumnWeights(options, db.GetGenomes());
            }
        }

        // Checks that would otherwise only fail after hours of alignment.
        std::vector<uint32_t> msa_taxids;
        ReadTypeModels models;
        if (run_profiling) {
            msa_taxids = ResolveMSASpecies(options, db.GetTaxonomy());
            for (size_t m = 0; m < model_files.size(); m++) {
                auto const& info = Info(model_files[m].type);
                auto const& model_file = model_files[m].file;
                // The gate it was trained at, where the database has models of more than one; its kind, and the other kinds
                // the database has models of these reads of (--model_kind).
                bool const gated = options.ModelGates(info.type).size() > 1 && options.ModelNamed(info.type).empty();
                size_t const gate = model_files[m].gate;
                std::string about = model_kinds[m];
                if (options.ModelNamed(info.type).empty() && options.ModelKind(info.type).empty()) {
                    for (auto const& kind : kModelKinds) {
                        if (kind == model_kinds[m]) continue;
                        bool has = options.DatabaseFileNamed(KindModelFile(info.model_file, kind)).Exists();
                        for (size_t g : kModelGates) has = has || options.DatabaseFileNamed(KindModelFile(GateModelFile(info.type, g), kind)).Exists();
                        if (has) about += (about.empty() ? "" : "; ") + std::string("the database also has ") + kind + " models: --model_kind " + kind;
                    }
                }
                std::cout << "Model of " << info.name << " reads"
                          << (gated ? " gated at " + (gate ? std::to_string(gate) : std::string("off")) : std::string())
                          << ": " << model_file.Name() << (about.empty() ? "" : " (" + about + ")") << std::endl;
                auto const& model = models.emplace(ModelKey(info.type, model_file), std::move(loaded_models[m])).first->second;
                if (options.FdrGiven() && options.GetFdr() > 0 && !model.HasFalseCalls()) {
                    std::cerr << "--fdr needs a model with a calibration (machine_learning_cmdline.py --fdr-calls); "
                              << model_file.Name() << " has none" << std::endl;
                    exit(2);
                }
                if (model.HasFalseCalls()) {
                    auto const& calls = model.FalseCalls();
                    bool const used = !options.KnobGiven() && options.FdrGiven() && options.GetFdr() > 0;
                    std::cout << "  calls at an expected share of false calls of " << profiler::FeatureString(calls.fdr)
                              << " (calibrated, " << calls.curve.size() << " points; training prior "
                              << profiler::FeatureString(calls.prior) << ")"
                              << (used ? "; at --fdr " + profiler::FeatureString(options.GetFdr()) + ", the depth knobs are not used"
                                       : options.KnobGiven() ? std::string("; not used, --knob is given")
                                                             : std::string("; not used (--fdr F would use them)"))
                              << std::endl;
                }
                if (!model.GetDepthKnobCurve().empty()) {
                    std::string knobs;
                    for (auto const& [x, knob] : model.GetDepthKnobCurve()) {
                        knobs += (knobs.empty() ? "" : ", ") + profiler::FeatureString(x) + ": " + profiler::FeatureString(knob);
                    }
                    std::cout << "  knobs by sample depth (log10 of the sample's fragments: knob; linear between, the "
                                 "ends' beyond): " << knobs << (options.KnobGiven() ? "; not used, --knob is given" : "")
                              << std::endl;
                } else if (!model.DepthKnobs().empty()) {
                    std::string knobs;
                    for (auto const& [bin, knob] : model.DepthKnobs()) {
                        knobs += (knobs.empty() ? "" : ", ") + std::to_string(bin) + ": " + profiler::FeatureString(knob);
                    }
                    std::cout << "  knobs by sample depth (bin b: 10^b to 10^(b+1) fragments, 2 also fewer, 6 also more): "
                              << knobs << (options.KnobGiven() ? "; not used, --knob is given"
                                                               : "; other depths --knob " + profiler::FeatureString(options.GetKnob()))
                              << std::endl;
                }
            }
            if (concurrent) {
                // Tables loading beside the alignment report when the profiling waits for them (DeferredTables::Wait).
                if (deferred.futures.empty()) {
                    std::cout << conservation_log.str() << suspect_log.str() << priors_log.str() << ncbi_log.str() << species_neighbours_log.str()
                              << gaps_log.str() << foreign_log.str() << weights_log.str() << sites_log.str();
                }
            } else {
                LoadGeneConservation(options, db.GetGenomes());
                LoadSuspectCopies(options, db.GetGenomes());
                LoadSpeciesPriors(options, db.GetGenomes());
                if (options.NcbiTaxids()) LoadNcbiTaxids(options, db.GetGenomes());
                LoadSpeciesNeighbours(options, db.GetGenomes());
                LoadCongenerGaps(options, db.GetGenomes());
                LoadForeignRates(options, db.GetGenomes());
                LoadColumnWeights(options, db.GetGenomes());
                LoadAncestrySites(options, db.GetGenomes());
            }
        }
        if (concurrent && weights_for_alignment && (!deferred.futures.empty() || !run_profiling)) std::cout << weights_log.str();
        if (concurrent && need_alleles && (alleles_for_alignment || !run_alignment)) std::cout << alleles_log.str();
        std::cout << cores_log.str();
        if (!need_alleles && run_alignment && !options.BuildMode()) {
            std::cout << "Strain alleles: not read (no profiling, and no short reads' alignments scored with them: long reads or "
                      << "--no_allele_scores)" << std::endl;
        }
        if (need_neighbours) {
            if (!concurrent) neighbours = ReadGeneNeighbours(options);
            std::cout << neighbours_log.str();
            if (neighbours) BindGeneNeighbours(options, db, std::move(*neighbours));
        }
        bm_tables.Stop();
        bm_tables.PrintResults();
        if (run_alignment && options.BenchmarkAlignment() && !options.GetRange().empty()) {
            // The benchmark takes each read's true gene from its name; without one it would stop
            // at the first read, deep inside the alignment.
            ThreadedGzIstream is{ options.GetFirstFiles(options.GetRange().front()) };
            SeqReader reader{ is };
            FastxRecord record;
            static const std::regex truth_name("^[0-9]+_[0-9]+([^0-9].*)?$");
            if (reader(record) && !std::regex_match(record.id, truth_name)) {
                std::cerr << "Error: --benchmark_alignment needs reads named <taxid>_<gene id>..., the gene each read "
                             "was simulated from; '" << record.id << "' is not (reads from simulate_metagenomes carry "
                             "no gene ids)." << std::endl;
                exit(2);
            }
        }

        // misc/cpu.tsv: each stage's wall-clock and CPU seconds from here on (CpuLog.h), the start-up's first.
        if (run_alignment || run_profiling) {
            CpuLog::Get().Open(std::filesystem::path(options.GetMiscOutputDir()) / "cpu.tsv", options.GetThreads());
            CpuLog::Get().Mark("start-up");
        }

        /*
         *  READ ALIGNMENT SECTION
         */
        // With --profile_ahead, a sample is profiled as soon as its SAM is complete, while the next sample's reads
        // are aligned (ProfilingAhead, on a quarter of the threads beside the alignment's); by default every
        // sample is profiled after the alignment of all, on all threads.
        std::optional<ProfilingContext> profiling;
        ProfilingAheadPtr ahead;
        if (run_profiling) {
            profiling.emplace(options, db, models);
            if (!deferred.futures.empty()) profiling->deferred = &deferred;
        }
        if (run_alignment && run_profiling && options.ProfileAhead() && !options.BenchmarkAlignment()) {
            ahead = MakeProfilingAhead(*profiling, std::max<size_t>(1, options.GetThreads() / 4));
        }
        auto sample_done = [&](size_t index) {
            if (!ahead) return;
            auto const& range = profiling->range;
            auto const at = std::find(range.begin(), range.end(), index);
            // Not the last sample: no alignment follows it, so the profiling stage takes it on all threads.
            if (at != range.end() && at + 1 != range.end()) SubmitProfilingAhead(*ahead, static_cast<size_t>(at - range.begin()));
        };

        if (run_alignment) RunAlignment(options, db, sample_done);

        /*
         * PROFILER
         */
        if (run_profiling) {
            Benchmark bm_profiling("Profiling");
            bm_profiling.Start();

            auto profiles = ProfileWrapper(*profiling, ahead.get());
            ahead.reset();
            // Which taxa pass, for the statistics and strains: the scores ProfileWrapper cached, each
            // from its sample's model (any model reads them), at each profile's knobs.
            auto& filter = models.begin()->second;
            bm_profiling.Stop();
            bm_profiling.PrintResults();
            CpuLog::Get().Mark("profiling");
            protal::build::PrintMemory("profiling");

            // Per taxon: its statistics in every sample it has reads in, one file per taxon (misc/<taxon>.statistics.tsv),
            // with --taxon_statistics: on a GTDB-sized database a sample has reads on thousands of taxa, and that many small
            // files took 15 s on a network file system (docs/claude/2026-10-04-performance-gtdb-scale); the profile and its
            // logs hold the same per sample.
            if (options.TaxonStatistics()) {
                Benchmark bm_statistics("Taxon statistics files");
                bm_statistics.Start();
                auto taxids = ExtractTaxa(profiles, {}, 1);
                auto& taxonomy = db.GetTaxonomy();

                for (auto taxid : taxids) {
                    TaxonStatisticsOutput stats;
                    TaxonStatisticsOutput all_stats;

                    std::string name = taxonomy.Get(taxid).scientific_name;
                    std::replace(name.begin(), name.end(), ' ', '_');

                    std::ofstream os(options.GetMiscOutputDir() + '/' + name + ".statistics.tsv");

                    stats.PrintHeader(os);
                    for (auto& profile : profiles) {
                        auto& taxa = profile.GetTaxa();
                        if (!taxa.contains(taxid)) continue;

                        auto& taxon = taxa.at(taxid);
                        bool accepted = filter.Calls(taxon, profile.Knob());
                        stats.PrintLine(os, profile.GetName(), taxon.VerticalCoverage(), taxon.TotalHits(), taxon.TotalLength(), taxon.GetMeanANI(), taxon.GetMeanMAPQ(), accepted);
                    }
                    os.close();
                    if (os.fail()) RunStatus::Get().Fail("Writing the statistics of " + name + " failed: " + options.GetMiscOutputDir() + '/' + name + ".statistics.tsv");
                }
                bm_statistics.Stop();
                bm_statistics.PrintResults();
            }

            /*
             * STRAIN PART -  RESOLVE MSAs BETWEEN SAMPLES
             */
            if (!options.NoStrains()) {
                // A sample enters a taxon's MSA with a score of its profile's MSA knob or more (--msa_knob, by
                // default the sample's knob: the samples whose profile reports the taxon; EntersMSA).
                StrainWrapper2(options, profiles, db.GetGenomes(), db.GetTaxonomy(), msa_taxids, filter);
                for (auto& profile : profiles) profile.RemoveSpillFile();  // --scratch_dir: the spill files (the folder goes at the end too)
            }
        }

        CpuLog::Get().Mark("rest of the run");  // a row only if the table was opened
        bm_total.Stop();
        bm_total.PrintResults();
        teardown.Arm();

        std::cout << std::flush;
        return RunStatus::Get().Finish();
    }
}
