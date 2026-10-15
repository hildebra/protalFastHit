// Profiling the samples, the profiling stage of protal::Run (ProfileStage.h).

#include <algorithm>
#include <chrono>
#include <condition_variable>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#include <omp.h>
#ifdef __GLIBC__
#include <malloc.h>  // malloc_trim
#endif
#include "ProfileStage.h"
#include "Benchmark.h"
#include "RunStatus.h"

namespace protal {

    // The UnreportedTaxon of each sample: <misc>/unreported_species.tsv and a warning. Without any, an
    // earlier run's file is removed.
    static void WriteUnreportedSpecies(Options const& options, std::vector<std::vector<UnreportedTaxon>> const& samples) {
        auto const dir = options.GetMiscOutputDir();
        if (dir.empty()) return;
        auto const path = dir + "/unreported_species.tsv";
        std::set<std::string> species;
        size_t rows = 0;
        for (auto const& sample : samples) {
            for (auto const& taxon : sample) species.insert(taxon.species);
            rows += sample.size();
        }
        if (rows == 0) {
            std::error_code ec;
            std::filesystem::remove(path, ec);
            return;
        }
        std::ofstream os(path, std::ios::out);
        os << "sample\tspecies\ttaxid\tscore\town_depth\thit_gene_fraction\ttop_identity\tlow_identity_share\tpasses_msa_knob\n";
        for (auto const& sample : samples) {
            for (auto const& taxon : sample) os << taxon.row << '\n';
        }
        os.close();
        if (os.fail()) RunStatus::Get().Fail("Writing " + path + " failed");
        std::string const knob = options.KnobGiven() ? "--knob " + profiler::FeatureString(options.GetKnob())
                                                     : "their sample's knob (" + profiler::FeatureString(options.GetKnob()) +
                                                       ", or the model's for the sample's depth)";
        std::string const msa_knob = options.MSAKnobGiven() ? "--msa_knob (" + profiler::FeatureString(options.GetMSAKnob()) + ")"
                                                            : "that knob";
        std::cerr << "Warning: " << species.size() << " species in " << rows << " sample profile(s) score below " << knob
                  << " and are not reported, although their own reads are strong evidence "
                  << "(depth 1x or more, reads on 90% of their genes, best reads 98% identical or more); they are "
                  << "listed in " << path << ". Their samples enter the strain MSAs if they score " << msa_knob
                  << " or more." << std::endl;
    }

    // Profiles the sample of slot idx with the model of its kind of reads, on threads_per_sample threads, into
    // ctx.profile_slots[idx] and ctx.unreported_slots[idx]; `filters` is the calling thread's own copy of the
    // models (scoring reuses a buffer; copies share the loaded model). Every taxon's score is then cached
    // (profile.ReleaseReadData scores all), so that later stages may use any model to tell which taxa pass: they
    // get the score of the sample's own model. Any thread may profile a sample of its own: a sample's profile
    // depends on its SAM and the database only, not on the thread or on the other samples.
    // A sample's profiling timers (wall clock, ProfileSample): written to its runtime table beside the alignment
    // stage's timers (misc/<sample>_runtime.tsv, made by the alignment stage of this run, else here with its header),
    // in place of the profiling rows an earlier profiling of the sample left there (a rerun once added a second set),
    // and printed in one line, `total` first and the `parts` after it. Each row has the table's columns: a wall-clock
    // timer ran on one thread, so seconds and seconds per thread are the same.
    static void WriteProfilingTimes(Options const& options, size_t index, Benchmark const& total, std::vector<Benchmark const*> const& parts) {
        auto const path = std::filesystem::path(options.GetMiscOutputDir()) / (options.GetSampleId(index) + "_runtime.tsv");
        std::string kept = "stage\tseconds\tthreads\tseconds_per_thread\n";  // the header and the other stages' rows
        {
            std::ifstream is(path);
            std::string line;
            for (bool header = true; std::getline(is, line); header = false) {
                if (!header && !line.empty() && line.rfind("Profiling: ", 0) != 0) kept += line + '\n';
            }
        }
        std::ofstream os(path, std::ios::out | std::ios::trunc);
        os << kept << std::fixed << std::setprecision(6);
        os << total.GetName() << '\t' << total.Seconds() << '\t' << total.Threads() << '\t' << total.MeanSeconds() << '\n';
        for (auto const* bm : parts) os << bm->GetName() << '\t' << bm->Seconds() << '\t' << bm->Threads() << '\t' << bm->MeanSeconds() << '\n';
        std::ostringstream line;
        line << std::fixed << std::setprecision(1) << "Profiling sample " << options.GetSampleId(index) << " took " << total.Seconds() << "s:";
        for (size_t p = 0; p < parts.size(); p++) {
            std::string name = parts[p]->GetName();
            if (name.rfind("Profiling: ", 0) == 0) name = name.substr(11);
            line << (p ? ", " : " ") << name << ' ' << parts[p]->Seconds() << 's';
        }
        #pragma omp critical(print)
        std::cout << line.str() << std::endl;
    }

    // One more profiling row of sample `index`'s runtime table (misc/<sample>_runtime.tsv), after WriteProfilingTimes wrote
    // its others: a "Profiling: " row, which WriteProfilingTimes replaces with the rest when the sample is profiled again.
    static void AppendProfilingTime(Options const& options, size_t index, Benchmark const& bm) {
        auto const path = std::filesystem::path(options.GetMiscOutputDir()) / (options.GetSampleId(index) + "_runtime.tsv");
        std::ofstream os(path, std::ios::out | std::ios::app);
        os << std::fixed << std::setprecision(6) << bm.GetName() << '\t' << bm.Seconds() << '\t' << bm.Threads() << '\t'
           << bm.MeanSeconds() << '\n';
    }

    static void ProfileSample(ProfilingContext& ctx, size_t idx, ReadTypeModels& filters, size_t threads_per_sample) {
        if (ctx.deferred) ctx.deferred->Wait();  // the profiling's tables, loaded beside the alignment
        Options& options = ctx.options;
        GenomeLoader& genomes = ctx.genomes;
        auto& taxonomy = *ctx.taxonomy;
        auto const& genera = ctx.genera;
        auto const& families = ctx.families;
        auto const& distances = ctx.distances;
        auto& profile_slots = ctx.profile_slots;
        auto& unreported_slots = ctx.unreported_slots;
        auto const& range = ctx.range;
        {
            auto i = range[idx];

            auto const read_type = options.GetReadType(i);
            auto const sample_filter = filters.find(ModelKey(read_type, options.SampleModelDbFile(i)));
            if (sample_filter == filters.end()) {  // loaded up front for every sample
                RunStatus::Get().Fail("No model loaded for the " + ReadTypeName(read_type) + " reads of sample " + options.GetSampleId(i));
                profile_slots[idx].emplace(genomes);
                return;
            }
            auto& filter = sample_filter->second;

            if (options.Verbose()) {
                auto [sam, compressed] = options.SamFile(i);
                #pragma omp critical(print)
                std::cerr << omp_get_thread_num() << " File " << i << " of " << range.size() << ":\n\t" << sam << (compressed ? " (compressed)" : "") << std::endl;
            }

            auto [sam, compressed] = options.SamFile(i);
            auto sample_name = options.GetSampleId(i);

            if (!Utils::exists(sam)) {
                RunStatus::Get().Fail("No SAM file to profile for sample " + options.GetSampleId(i) + ": " + sam);
                profile_slots[idx].emplace(genomes);
                return;
            }

            Benchmark bm_read_alignments{ "Load read alignments" };
            Benchmark bm_profile{ "Profile sample" };
            // Wall-clock timers of this sample's profiling, written to its runtime table beside the alignment's
            // stage timers (misc/<sample>_runtime.tsv) and printed in one line: the whole, and the scoring and the
            // writing here; the profiler's own steps (reading the SAM, the evidence, the SNPs) are its timers.
            Benchmark bm_sample{ "Profiling: sample" };
            Benchmark bm_score{ "Profiling: scoring" };
            Benchmark bm_write{ "Profiling: writing the profile" };
            bm_sample.Start();

            profiler::Profiler profiler(genomes);
            profiler.SetDepthIdentityMargin(options.GetDepthIdentityMargin());


            if (options.Verbose()) {
                #pragma omp critical(print)
                {
                    std::cout << "Thread " << omp_get_thread_num() << " read sam file " << sam << std::endl;
                }
                if (!Utils::exists(sam)) {
                    std::cerr << "File does not exist" << sam << std::endl;
                    exit(90);
                }
            }

            // The records whose alignment does not fit the database, for misc/<sample>.err (below).
            std::ostringstream rejected_records;
            bm_profile.Start();

            // Samples are read in parallel: a SAM is streamed, not held in memory, so only the
            // profiles themselves take memory.
            profiler::MicrobialProfile profile(genomes);
            profile.SetName(sample_name);
            profile.SetReadType(read_type);
            profile.SetGenera(genera);
            profile.SetFamilies(families);
            profile.SetCongenerDistances(distances);
            profile.SetSingletonCongener(options.GetSingletonCongener());
            profile.SetSoftFragments(options.GetFragmentCount() == FragmentCount::Soft);
            profile.SetDropForeignGenes(options.DropForeignGenes());
            // A long-read sample's strains get strain MSA rows of their own (Haplotypes.h).
            profile.SetKeepPhaseRecords(!options.NoStrains() && options.Phasing() && IsLongReadType(read_type));
            std::string sam_error = profiler.ProfileSam(sam, profile, std::optional<std::reference_wrapper<std::ostream>>{rejected_records},
                                                        options.GetSNPMinCov(), options.GetSNPMinCov(),
                                                        options.GetSNPMinAF(read_type), options.GetSNPMinMeanQual(),
                                                        options.GetSNPMinPhredSum(), options.GetSNPRequireStrand(),
                                                        threads_per_sample);
            bm_profile.Stop();
            // A SAM in the run's scratch folder (--scratch_dir with --no_strains) is read by nothing after its profile.
            if (options.IsScratchSam(i)) {
                std::error_code ec;
                std::filesystem::remove(sam, ec);
            }
            // misc/<sample>.err: the records whose alignment does not fit the database, in the run's own folder (it was
            // <sam>.err, next to a SAM that may be another run's, or in a folder protal cannot write to); a sample without
            // any gets none, and an earlier run's is removed.
            auto const rejected_file = std::filesystem::path(options.GetMiscOutputDir()) / (sample_name + ".err");
            {
                std::string const records = rejected_records.str();
                std::ostringstream().swap(rejected_records);
                std::error_code ec;
                if (records.empty()) {
                    std::filesystem::remove(rejected_file, ec);
                } else {
                    std::ofstream os(rejected_file, std::ios::out | std::ios::trunc);
                    os << records;
                    os.close();
                    if (os.fail()) RunStatus::Get().Fail("Writing the rejected records of sample " + sample_name + " failed: " + rejected_file.string());
                }
            }

            // A SAM without alignments still gets its (empty) profile files, so that every sample
            // has output; only an unreadable SAM is a failure.
            if (!sam_error.empty()) {
                RunStatus::Get().Fail("Cannot read the SAM file of sample " + options.GetSampleId(i) + " (" + sam + "): " + sam_error);
                profile_slots[idx].emplace(genomes);
                return;
            }
            // The sample's threshold, on this thread's copy of the model: --knob if given, else the model's knob for the
            // sample's depth (its fragments over all taxa) if it has depth knobs (the trainer's), else --knob's default.
            // Its taxa enter the strain MSAs at --msa_knob if given, else at the same.
            // With --fdr F above 0 and a model's calibration (machine_learning_cmdline.py --fdr-calls): the knob at which
            // the sample's expected share of false calls is at most F (context::FalseCallKnob), over its taxa the
            // singleton rule does not veto. Without --fdr the calibrated calls are not used: at GTDB r226 they called
            // 0.001-0.007 F1 below the knob curve for every read type (docs/claude/2026-10-03-r226-v5-v6-training).
            filter.SetKnob(options.GetKnob());
            bool const false_call_mode = !options.KnobGiven() && filter.HasFalseCalls() &&
                                         options.FdrGiven() && options.GetFdr() > 0;
            if (false_call_mode) {
                profile.ScoreTaxa(filter, threads_per_sample);
                std::vector<double> scores;
                for (auto const& [_, taxon] : profile.GetTaxa()) {
                    if (!taxon.Vetoed()) scores.push_back(filter.Score(taxon));
                }
                auto const& model = filter.FalseCalls();
                double const fdr = options.FdrGiven() ? options.GetFdr() : model.fdr;
                auto const calls = profiler::context::FalseCallKnob(std::move(scores), model.curve, model.prior, fdr);
                filter.SetKnob(calls.knob);
                #pragma omp critical(print)
                std::cout << "Sample " << sample_name << ": " << profile.Fragments() << " fragments, " << calls.called
                          << " taxa at an expected share of false calls of at most " << profiler::FeatureString(fdr)
                          << " (" << profiler::FeatureString(std::round(calls.expected_false * 100) / 100)
                          << " expected; prior adjusted to the sample " << profiler::FeatureString(calls.sample_prior)
                          << "), knob " << profiler::FeatureString(calls.knob) << std::endl;
            } else if (!options.KnobGiven() && filter.HasDepthKnobs()) {
                size_t const fragments = profile.Fragments();
                if (auto const depth_knob = filter.DepthKnob(fragments)) {
                    filter.SetKnob(*depth_knob);
                    #pragma omp critical(print)
                    std::cout << "Sample " << sample_name << ": " << fragments << " fragments, knob "
                              << profiler::FeatureString(*depth_knob) << " (the model's "
                              << (filter.GetDepthKnobCurve().empty() ? "for depth bin " + std::to_string(profiler::DepthKnobBin(fragments))
                                                                     : std::string("for that depth")) << ")" << std::endl;
                }
            }
            double const msa_knob = options.MSAKnobGiven() ? options.GetMSAKnob() : filter.GetKnob();
            profile.SetKnobs(filter.GetKnob(), msa_knob);

            // The outputs below score every taxon; on the sample's threads first.
            bm_score.Start();
            if (threads_per_sample > 1) profile.ScoreTaxa(filter, threads_per_sample);
            bm_score.Stop();
            // The taxa near or above the knob whose features lie outside the model's training range (FeatureRanges.h; a
            // model with the trainer's protal_feature_ranges): a line for the sample, and each taxon's in its .profile.log.
            if (auto const& ranges = filter.FeatureRanges()) {
                feature_ranges::SampleSummary summary;
                std::unordered_map<uint32_t, std::string> outside;
                for (auto const taxid : profile.SortedTaxa()) {
                    auto const& taxon = profile.GetTaxa().at(taxid);
                    if (!feature_ranges::NearKnob(filter.Score(taxon), filter.GetKnob(), taxon.Fragments())) continue;
                    auto const found = feature_ranges::Check(*ranges, profiler::TaxonFeatures(taxon));
                    summary.Add(found);
                    if (!found.hard.empty()) outside.emplace(taxid, feature_ranges::Joined(found.hard));
                }
                profile.SetOutsideTrainingRange(std::move(outside));
                if (auto const line = summary.Line(sample_name); !line.empty()) {
                    #pragma omp critical(print)
                    std::cerr << "Warning: " << line << std::endl;
                }
            }
            if (profiler.RejectedReads() > 0) {
                #pragma omp critical(print)
                std::cerr << "Warning: sample " << sample_name << ": " << profiler.RejectedReads() << " of " << profiler.Reads()
                          << " reads have an alignment that does not fit the database (a gene it lacks, a position past "
                             "a gene's end, or bases that differ from the gene); they are left out and listed in "
                          << rejected_file.string() << std::endl;
            }
            if (size_t const suspect = profile.SuspectRecords(); suspect > 0) {
                #pragma omp critical(print)
                std::cout << "Sample " << sample_name << ": " << suspect << " records on suspect gene copies left out of the "
                          << "evidence (suspect_copies.tsv; --keep_suspect_copies counts them)" << std::endl;
            }
            if (auto const [foreign, foreign_taxa] = profile.ForeignGenes(); foreign > 0) {
                #pragma omp critical(print)
                std::cout << "Sample " << sample_name << ": " << foreign << " foreign genes of " << foreign_taxa << " taxa (the genes "
                          << "next to them on their reads are mostly unlikely neighbours in the taxon's clade), "
                          << (options.DropForeignGenes() ? "left out of their depth and strain MSAs" : "kept (--keep_foreign_genes)")
                          << "; Foreign in " << options.ProfileFile(i) << ".genes.log" << std::endl;
            }


            if (options.Verbose()) {
                #pragma omp critical(print)
                {
                    std::cout << "Thread " << omp_get_thread_num() << " ";
                    bm_profile.PrintResults();
                }
            }

            // The sample's composition (Composition.h): the share of its reads its called species explain, its average
            // genome size and the unknown share that ends its profile ("?"; --no_unknown_share: the called species'
            // abundances are shares of their summed depth, as before 2026-10-08).
            {
                auto result = composition::Compute(profile.CalledSpecies(filter), profiler.Scanned(), profile.FragmentBaseShare(),
                                                   genomes.GetSpeciesPriors().WithGenomeSize() > 0);
                #pragma omp critical(print)
                std::cout << "Sample " << sample_name << ": " << composition::Summary(result) << std::endl;
                profile.SetComposition(std::move(result), !options.NoUnknownShare());
            }

            std::optional<TruthSet> truth = options.HasProfileTruths() ?
                                            std::optional<TruthSet>{ protal::GetTruth(options.ProfileTruthFile(i), taxonomy) } :
                                            std::optional<TruthSet>{};

            if (truth.has_value()) {
                std::string truth_output = options.ProfileFile(i) + ".truth_annotated";
                profile.AnnotateWithTruth(truth.value(), filter, truth_output, taxonomy);
                std::cout << "Write truth to: " << truth_output << std::endl;

                // A true species the database lacks cannot be found: it is counted apart from the misses.
                size_t tp = 0, fp = 0, fn = 0, not_in_db = 0;
                for (auto const& [taxid, _] : profile.GetTaxa()) {
                    if (!filter.Pass(profile.GetTaxa().at(taxid))) continue;
                    (truth.value().contains(taxid) ? tp : fp)++;
                }
                for (auto taxid : truth.value()) {
                    if (!genomes.GetGenomeMap().contains(taxid)) not_in_db++;
                    else if (!profile.GetTaxa().contains(taxid) || !filter.Pass(profile.GetTaxa().at(taxid))) fn++;
                }
                #pragma omp critical(print)
                std::cout << "Sample " << sample_name << ": TP " << tp << ", FP " << fp << ", FN " << fn
                          << " (and " << not_in_db << " true species not in the database)" << std::endl;

            }
            if (options.Verbose()) {
                std::cout << "Write profile to: \n" << options.ProfileFile(i) << std::endl;
            }
            {
                bool made;
#pragma omp critical(create_dir)
                made = MakeFolder(std::filesystem::path(options.ProfileFile(i)).parent_path());
                if (!made) {
                    std::cerr << "Cannot create directories for path " << options.ProfileFile(i) << std::endl;
                    exit(2);
                }
            }

            bm_write.Start();
            std::ofstream os(options.ProfileFile(i), std::ios::out);
            std::ofstream os_total(options.ProfileFile(i) + ".log", std::ios::out);
            std::ofstream os_dismissed(options.ProfileFile(i) + ".gene.log", std::ios::out);
            std::ofstream os_genes(options.ProfileFile(i) + ".genes.log", std::ios::out);
            std::ofstream os_composition(options.ProfileFile(i) + ".composition", std::ios::out);

            profile.WriteSparseProfile(taxonomy, filter, os, &os_total, &os_dismissed, threads_per_sample,
                                       options.NcbiTaxids() ? &genomes.GetNcbiTaxids() : nullptr);
            profile.WriteGeneProfile(taxonomy, filter, &os_genes, threads_per_sample);
            composition::WriteHeader(os_composition);
            composition::WriteRow(os_composition, options.GetSampleId(i), *profile.GetComposition());
            os.close();
            os_total.close();
            os_dismissed.close();
            os_genes.close();
            os_composition.close();
            bm_write.Stop();
            if (os.fail() || os_total.fail() || os_dismissed.fail() || os_genes.fail() || os_composition.fail()) {
                RunStatus::Get().Fail("Writing the profile of sample " + options.GetSampleId(i) + " failed: " + options.ProfileFile(i));
            }
            bm_sample.Stop();
            // What follows the sample's timers, by part: at r226 1.9 s of the profiling stage that no timer saw
            // (docs/claude/2026-10-12-feature-cost-vs-f1.md). A row "Profiling: after the profile" in the runtime table,
            // and a line.
            Benchmark bm_after{ "Profiling: after the profile" };
            bm_after.Start();
            std::ostringstream after_parts;
            auto part_start = std::chrono::steady_clock::now();
            auto part_done = [&after_parts, &part_start](char const* part) {
                auto const now = std::chrono::steady_clock::now();
                after_parts << (after_parts.tellp() > 0 ? ", " : "") << part << ' ' << std::fixed << std::setprecision(2)
                            << std::chrono::duration<double>(now - part_start).count() << " s";
                part_start = now;
            };
            WriteProfilingTimes(options, i, bm_sample, { &profiler.ReadTimer(), &profiler.EvidenceTimer(), &profile.m_bm_em,
                                                         &profile.m_bm_distances, &profiler.SnpTimer(), &bm_score, &bm_write });
            if (!profiler.ReadDetail().empty()) {
                #pragma omp critical(print)
                std::cout << "Sample " << options.GetSampleId(i) << " SAM read: " << profiler.ReadDetail() << std::endl;
            }
            #pragma omp critical(print)
            std::cout << "Sample " << options.GetSampleId(i) << " read EM: " << profile.EmDetail() << std::endl;
            part_done("the runtime table and lines");

            profile.SetName(options.GetSampleId(i));

            // In taxid order: the map's follows the order of the SAM's records.
            for (auto const taxid : profile.SortedTaxa()) {
                auto const& taxon = profile.GetTaxa().at(taxid);
                double const score = filter.Score(taxon);
                if (filter.Calls(taxon, filter.GetKnob()) || !profiler::StrongOwnEvidence(taxon)) continue;
                std::string const species = taxonomy.Get(taxid).scientific_name;
                std::ostringstream line;
                line << options.GetSampleId(i) << '\t' << species << '\t' << taxid << '\t' << score << '\t'
                     << taxon.VerticalCoverage() << '\t' << taxon.HitGeneFraction() << '\t' << taxon.TopIdentity() << '\t'
                     << taxon.LowIdentityShare() << '\t' << (filter.Calls(taxon, msa_knob) ? "yes" : "no");
                unreported_slots[idx].push_back({ species, line.str() });
            }
            part_done("unreported species");

            // The sample's outputs are written: only the strain stage reads its profile again, and only the taxa that
            // enter the MSAs (--msa_knob), of which it keeps what the MSAs read (Gene::KeepForMSAs); --taxon_statistics
            // reads the numbers of every taxon.
            std::optional<profiler::MSAReleaseParameters> msa;
            if (!options.NoStrains()) msa = MSAReleaseFor(options, read_type);
            profile.ReleaseReadData(filter.WithKnob(msa_knob), msa, options.TaxonStatistics(), threads_per_sample);
            part_done("releasing the read data");
            profile.ReleaseEvidence();  // the record evidence: read by nothing after the outputs; a cohort keeps the profile
            part_done("releasing the record evidence");
            // --scratch_dir: what it keeps goes to a file of its own in the run's scratch folder until the strain stage reads
            // it back per species.
            if (msa && !options.ScratchFolder().empty()) {
                auto const path = std::filesystem::path(options.ScratchFolder()) / ("strains_" + std::to_string(i) + ".bin");
                if (!profile.SpillStrainEvidence(path.string())) {
                    #pragma omp critical(print)
                    std::cerr << "Warning: sample " << sample_name << ": cannot write " << path.string()
                              << " (--scratch_dir); its strain evidence stays in memory" << std::endl;
                }
                part_done("spilling the strain evidence");
            }
            bool const last_profile = ctx.profiles_left.fetch_sub(1) == 1;  // the run's last sample to get here
#ifdef __GLIBC__
            // The sample's reads' records are freed, but glibc keeps the pages of its arenas that the blocks which stay
            // (what the profile keeps, made after the records) leave partly used; returned, a cohort's memory grows by
            // what its samples keep only (docs/claude/2026-10-07-sam-combine). malloc_trim walks every arena on this
            // thread, 0.4-0.5 s for an r226 paired-end sample (docs/claude/2026-10-12-feature-cost-vs-f1.md), so it is
            // left out where it returns nothing worth the time: after the run's last profile (no other sample's records
            // come to take the pages; the strain stage reuses the freed blocks), and while another sample's trim runs
            // (it walks every arena already; a second one would wait for their locks, and the samples' threads with it).
            if (last_profile) {
                part_done("malloc_trim skipped: the run's last profile");
            } else if (ctx.trimming.exchange(true)) {
                part_done("malloc_trim skipped: another sample's runs");
            } else {
                malloc_trim(0);
                ctx.trimming = false;
                part_done("returning freed memory (malloc_trim)");
            }
#else
            (void)last_profile;
#endif
            bm_after.Stop();
            AppendProfilingTime(options, i, bm_after);
            {
                std::ostringstream line;
                line << "Sample " << options.GetSampleId(i) << " after its profile: " << std::fixed << std::setprecision(2)
                     << bm_after.Seconds() << " s (" << after_parts.str() << ")";
                #pragma omp critical(print)
                std::cout << line.str() << std::endl;
            }

            // Each thread writes to its own pre-allocated slot — no lock needed.
            profile_slots[idx].emplace(std::move(profile));
        }
    }

    // Profiles samples while others are still being aligned (--profile_ahead). RunWrapper hands over each sample
    // whose SAM is complete (Submit); a worker thread profiles them one after another on `threads` threads of its
    // own, beside the alignment's (ProfileSample); the profiling stage then takes whatever the worker had not
    // started (ProfileWrapper: StopTaking, then Started says which). The profiles are the same as when every
    // sample is profiled after the alignment: a sample's profile depends on its SAM and the database only. This
    // can only pay where the profiling stage leaves cores idle (its serial passes, the tail of a cohort): on a
    // 6-core laptop, with the alignment on every core, the worker's CPU time came out of the alignment's and runs
    // of 4 and 8 samples took the same time either way (docs/claude/2026-10-03-performance-review), so it is off by
    // default. The last sample's profile, and any the worker is behind with, follow the alignment on all threads.
    class ProfilingAhead {
    public:
        ProfilingAhead(ProfilingContext& ctx, size_t threads) :
                m_ctx(ctx), m_threads(std::max<size_t>(threads, 1)), m_filters(ctx.models), m_started(ctx.range.size(), 0) {
            m_worker = std::thread([this] { Work(); });
        }

        ~ProfilingAhead() {
            StopTaking();
            Join();
        }

        ProfilingAhead(ProfilingAhead const&) = delete;
        ProfilingAhead& operator=(ProfilingAhead const&) = delete;

        // The sample of slot idx has its SAM: the worker profiles it when it is free.
        void Submit(size_t idx) {
            {
                std::lock_guard<std::mutex> lock(m_mutex);
                if (m_stopping) return;
                m_queue.push_back(idx);
            }
            m_wake.notify_one();
        }

        // No further sample is started; those queued but not started are the caller's (Started).
        void StopTaking() {
            {
                std::lock_guard<std::mutex> lock(m_mutex);
                m_stopping = true;
            }
            m_wake.notify_all();
        }

        // Whether the worker has started (or finished) the sample of slot idx. After StopTaking the answer is final.
        bool Started(size_t idx) {
            std::lock_guard<std::mutex> lock(m_mutex);
            return m_started[idx] != 0;
        }

        size_t Threads() const { return m_threads; }

        // Waits for the worker: after StopTaking, until the sample it is on is profiled.
        void Join() {
            if (m_worker.joinable()) m_worker.join();
        }

    private:
        void Work() {
            while (true) {
                size_t idx;
                {
                    std::unique_lock<std::mutex> lock(m_mutex);
                    m_wake.wait(lock, [this] { return m_stopping || !m_queue.empty(); });
                    if (m_stopping) return;
                    idx = m_queue.front();
                    m_queue.pop_front();
                    m_started[idx] = 1;
                }
                try {
                    ProfileSample(m_ctx, idx, m_filters, m_threads);
                } catch (std::exception const& e) {
                    RunStatus::Get().Fail("Profiling sample " + m_ctx.options.GetSampleId(m_ctx.range[idx]) + " failed: " + e.what());
                }
            }
        }

        ProfilingContext& m_ctx;
        size_t m_threads;
        ReadTypeModels m_filters;  // the worker's own copy of the models
        std::mutex m_mutex;
        std::condition_variable m_wake;
        std::deque<size_t> m_queue;
        std::vector<char> m_started;
        bool m_stopping = false;
        std::thread m_worker;
    };

    // Profiles the samples of ctx that are not yet profiled, each with the model of its kind of reads, and returns
    // every sample's profile (in the order of the samples). With `ahead`, the samples it profiled or is profiling
    // while the reads were aligned are left to it, and it is waited for.
    Profiles ProfileWrapper(ProfilingContext& ctx, ProfilingAhead* ahead) {
        Options& options = ctx.options;
        auto const& range = ctx.range;

        // Each thread scores with its own copy (firstprivate): scoring reuses a buffer. Copies share
        // the loaded model.
        ReadTypeModels filters = ctx.models;

        // The samples profiled here: all of them, or those the worker had not started when the alignment ended.
        if (ahead) ahead->StopTaking();
        std::vector<size_t> todo;
        for (size_t idx = 0; idx < range.size(); idx++) {
            if (!ahead || !ahead->Started(idx)) todo.push_back(idx);
        }
        if (ahead) {
            std::cout << "Profiled while the reads were aligned (on " << ahead->Threads() << " thread(s)): "
                      << range.size() - todo.size() << " of " << range.size() << " sample(s); " << todo.size() << " left" << std::endl;
        }

        omp_set_num_threads(options.GetThreads());

        // Samples are profiled in parallel, the largest SAM first, each next one by whichever thread is
        // free: samples differ in depth, and with the default static schedule one thread could be left
        // with several deep ones. A sample is profiled on threads in proportion to its SAM's share of
        // all the samples' bytes, at least one (ProfileSam: the same profile on any number), so that a
        // deep sample among shallow ones, or alone, is not left to one thread.
        size_t const threads = std::max<size_t>(options.GetThreads(), 1);
        std::vector<uintmax_t> sam_bytes(todo.size(), 0);
        for (size_t t = 0; t < todo.size(); t++) {
            std::error_code ec;
            auto const bytes = std::filesystem::file_size(options.SamFile(range[todo[t]]).first, ec);
            if (!ec) sam_bytes[t] = bytes;
        }
        uintmax_t const total_bytes = std::accumulate(sam_bytes.begin(), sam_bytes.end(), uintmax_t{0});
        std::vector<int> order(todo.size());
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(), [&sam_bytes](int a, int b) { return sam_bytes[a] > sam_bytes[b]; });
        auto threads_of = [&](int t) {
            if (total_bytes == 0) return std::max<size_t>(1, threads / std::max<size_t>(todo.size(), 1));
            double const share = static_cast<double>(sam_bytes[t]) / static_cast<double>(total_bytes);
            return std::clamp<size_t>(static_cast<size_t>(std::lround(share * static_cast<double>(threads))), 1, threads);
        };
        int const sample_threads = static_cast<int>(std::clamp<size_t>(todo.size(), 1, threads));
        #pragma omp parallel for schedule(dynamic, 1) num_threads(sample_threads) firstprivate(filters) shared(ctx, todo, order, threads_of)
        for (int k = 0; k < static_cast<int>(todo.size()); k++) {
            int const t = order[k];
            ProfileSample(ctx, todo[t], filters, threads_of(t));
        }
        if (ahead) ahead->Join();

        std::vector<profiler::MicrobialProfile> profiles;
        profiles.reserve(ctx.profile_slots.size());
        for (auto& slot : ctx.profile_slots) {
            if (slot.has_value()) profiles.emplace_back(std::move(slot.value()));
        }
        WriteUnreportedSpecies(options, ctx.unreported_slots);
        return profiles;
    }

    // As above, with the context made here (--profile_only, or the profiling stage without samples profiled ahead).
    Profiles ProfileWrapper(Options& options, ProtalDB& db, ReadTypeModels const& models) {
        ProfilingContext ctx(options, db, models);
        return ProfileWrapper(ctx);
    }

//     void ProfileWrapper2(Options& options, ProtalDB& db) {
//         std::cout << "ProfileWrapper2" << std::endl;
//         GenomeLoader& genomes = db.GetGenomes();
//
//         if (!db.IsTaxonomyLoaded()) db.LoadTaxonomy(options.TaxonomyDbFile());
//         auto& taxonomy = db.GetTaxonomy();
//
// //        omp_set_num_threads(6);
//
// #pragma omp parallel for default(none) shared(options, cout, taxonomy, genomes, db)
//         for (auto i : options.GetRange()) {
//         //for (auto i = 0; i < options.GetFileCount(); i++) {
//
//             Profiler::AlignmentContainer ac;
//             auto sam = options.SamFile(i);
//
//
// #pragma omp critical(read_sam)
//             ac.LoadSam(sam);
//
//             const auto& alignments = ac.AlignmentPairs();
//             Profiler::ReadFilter filter;
//             Profiler::Profiler<Profiler::ReadFilter> profiler(filter);
//             profiler.Profile(ac);//, db);
//
// //             profiler(genomes);
// //
// //            profiler.FromSam(sam);
//         }
//     }

    void ProfilingAheadDeleter::operator()(ProfilingAhead* ahead) const {
        delete ahead;
    }

    ProfilingAheadPtr MakeProfilingAhead(ProfilingContext& ctx, size_t threads) {
        return ProfilingAheadPtr(new ProfilingAhead(ctx, threads));
    }

    void SubmitProfilingAhead(ProfilingAhead& ahead, size_t idx) {
        ahead.Submit(idx);
    }

    std::shared_ptr<std::vector<uint32_t> const> TaxonomyGenera(taxonomy::IntTaxonomy const& taxonomy) {
        return profiler::GeneraOf(taxonomy);
    }
}
