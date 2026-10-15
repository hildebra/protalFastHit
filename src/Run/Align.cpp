// Aligning the samples' reads (RunWrapper), the alignment stage of protal::Run (RunStages.h).

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <vector>
#include "RunStages.h"
#include "Seedmap.h"
#include "Alignment/WFA2Wrapper2.h"
#include "AlignmentStrategy.h"
#include "ChainAnchorFinder.h"
#include "Classify.h"
#include "CpuLog.h"
#include "MemoryUse.h"
#include "RunStatus.h"
#include "SamFile.h"
#include "SeedExpectation.h"
#include "Stage1Probe.h"
#include "ThreadedGzStream.h"

namespace protal {

    // Alignments are written under a temporary name ("<sam>.partial", compressed already as the name
    // asks; SamOutput) and only get their final name once complete, so an interrupted run never leaves
    // a truncated SAM that a rerun would skip and reuse. This moves a finished file into place.
    sam_settings::Fields AlignmentSettings(Options const& options, size_t index, bool unique_cores) {
        namespace s = sam_settings;
        auto const type = options.GetReadType(index);
        bool const short_reads = !IsLongReadType(type);
        auto const run = options.RunSettingsInUse();
        s::Fields fields;
        auto add = [&fields](std::string_view key, std::string value) { fields.emplace_back(std::string(key), std::move(value)); };
        add(s::kReadType, Info(type).token);
        add(s::kGate, std::to_string(options.GetSparseFirst(index)));
        add(s::kPick, std::string(SparseFirstPickName(options.GetSparseFirstPick(index))));
        add(s::kScore, std::to_string(options.GetSparseFirstScore()));
        add(s::kHits, std::to_string(options.GetSparseFirstHits()));
        add(s::kUniqueKmers, std::string(unique_cores ? build_settings::kUniqueCores : build_settings::kUniqueAnchors));
        add(s::kTraceback, std::string(WfaTracebackName(options.GetWfaTraceback())));
        if (short_reads) add(s::kAlleleScores, build_settings::OnOff(*run.allele_scores));
        if (short_reads) add(s::kWeightedSiteShift, build_settings::OnOff(*run.weighted_site_shift));
        if (type == ReadType::Paired) add(s::kMateGuidance, *run.mate_guidance);
        add(s::kAlignTop, std::to_string(options.GetAlignTop()));
        if (short_reads) add(s::kAdaptiveCandidates, std::to_string(options.GetAdaptiveCandidates()));
        add(s::kMaxScoreAni, build_settings::Decimal(options.GetMaxScoreAni(type)));
        return fields;
    }

    static void FinishSamFile(Options& options, int index, std::string const& partial) {
        auto [sam, _] = options.SamFile(index);
        std::error_code ec;
        std::filesystem::rename(partial, sam, ec);
        if (ec) RunStatus::Get().Fail("Cannot move " + partial + " to " + sam + ": " + ec.message());
    }

    // One line per sample on the alignment stage's counts, for comparing runs (scripts/measure_performance.sh reads
    // it): the reads (or read pairs), those with an anchor (of paired reads, the mates), the candidate alignments tried (AlignAnchor calls, up to
    // --align_top per read or mate and the ties, plus the long-read rescues), those the k-mer screen refused before
    // WFA2 (AlignmentScreen.h) and those the indel bound refused in its place (ONT reads, IndelBound; 0 for others),
    // those WFA2 ran on from the anchor's exact matches or as a whole window, the alignments made (within the ANI
    // floor) and the records written.
    static void PrintAlignmentCounts(Options const& options, size_t index, ReadType read_type, Statistics const& stats,
                                     SimpleAlignmentHandler const& handler, SeedingCounts const& seeding) {
        std::cout << "Sample " << options.GetSampleId(index) << ": " << stats.reads << (read_type == ReadType::Paired ? " read pairs, " : " reads, ")
                  << stats.at_least_one_anchor << (read_type == ReadType::Paired ? " mates" : "") << " with an anchor; "
                  << handler.m_attempted_alignments << " candidate alignments tried: "
                  << handler.m_screened_alignments << " refused by the k-mer screen, " << handler.m_indel_refused << " by the indel bound, "
                  << handler.m_anchored_alignments
                  << " aligned from the anchor's exact matches and " << handler.m_whole_window_alignments << " as whole windows; "
                  << stats.total_alignments << " alignments made, " << stats.output_alignments << " records written" << std::endl;
        // And what the seeding did (SeedingCounts): the k-mer lookups, the value blocks scanned and their sizes.
        std::cout << "Sample " << options.GetSampleId(index) << " seeding: " << seeding.Text() << std::endl;
        if (seeding.sparse_reads > 0) {  // --sparse_first: the reads given up after their sampled lookups
            auto const reads = static_cast<double>(seeding.sparse_reads);
            bool const stride = options.GetSparseFirst(index) > 1 && options.GetSparseFirstPick(index) == SparseFirstPick::Stride;
            bool const grid = options.GetSparseFirst(index) > 1 && options.GetSparseFirstPick(index) == SparseFirstPick::Grid;
            std::cout << std::fixed << std::setprecision(2) << "Sample " << options.GetSampleId(index) << " sparse first (1/"
                      << options.GetSparseFirst(index) << " of the cores" << (stride ? " by stride" : grid ? " by the database's grid" : "")
                      << ", score "
                      << options.GetSparseFirstScore() << ", "
                      << options.GetSparseFirstHits() << " hit(s)): " << seeding.sparse_given_up << " of " << seeding.sparse_reads
                      << (read_type == ReadType::Paired ? " pairs" : " reads") << " given up ("
                      << 100.0 * static_cast<double>(seeding.sparse_given_up) / reads << "%), " << seeding.sparse_lookups
                      << " sampled lookups (" << static_cast<double>(seeding.sparse_lookups) / reads << " a "
                      << (read_type == ReadType::Paired ? "pair" : "read") << "); the seeding line counts the kept ones'"
                      << std::defaultfloat << std::endl;
        }
        if (handler.m_adaptive_alignments > 0) {
            std::cout << "Sample " << options.GetSampleId(index) << " adaptive candidates: " << handler.m_adaptive_alignments
                      << " of the candidate alignments tried beyond --align_top on congeners of divergent reads' best "
                      << "alignment (--adaptive_candidates)" << std::endl;
        }
        if (handler.m_allele_scored > 0) {
            std::cout << "Sample " << options.GetSampleId(index) << " strain alleles: " << handler.m_allele_scored
                      << " candidate alignments on gene copies with known alleles, " << handler.m_allele_shifted
                      << " of them with differences the species' alleles have or where they vary; " << handler.m_unsure_reads
                      << " reads with candidates of two species within " << handler.kUnsureMismatches
                      << " mismatches took those shifts, " << handler.m_settled_reads << " of them to another species "
                      << "(strain_alleles.tsv; --no_allele_scores: not)" << std::endl;
        }
    }

    // The seed expectation of each grid level (SeedExpectation.h; nullopt where the index in memory does not hold the level
    // whole, below a sparse load's): its values' unique shares by the lookups' rule (Seedmap::UniqueSharesByLevel) over
    // the whole index's (unique_kmers.tsv, summed over the genes), for ZU and ZT. A read seeded from that level alone finds
    // about those shares of unique values; the profile scales the genes' shares by them (Taxon::UniqueSeedEnrichment).
    // Empty without grid levels or without the genes' unique counts.
    static std::vector<std::optional<SeedExpectation>> SeedExpectations(Seedmap const& map, GenomeLoader& genomes,
                                                                        int threads, std::ostream& log) {
        uint64_t whole_unique = 0, whole_two = 0, whole_total = 0;
        for (auto const& [taxid, genome] : genomes.GetGenomeMap()) {
            auto const [su, lu, lsu, total] = genome.GetUniqueKmerCounts();
            whole_unique += su + lu;
            whole_two += su + lsu;
            whole_total += total;
        }
        if (!map.KeepsGridLevels() || whole_total == 0 || whole_unique == 0 || whole_two == 0) return {};
        auto const begin = std::chrono::steady_clock::now();
        auto const shares = map.UniqueSharesByLevel(threads);
        std::vector<std::optional<SeedExpectation>> expectations(shares.size());
        double const unique_share = static_cast<double>(whole_unique) / static_cast<double>(whole_total);
        double const two_share = static_cast<double>(whole_two) / static_cast<double>(whole_total);
        std::ostringstream os;
        for (size_t level = map.LoadedLevel(); level < shares.size(); level++) {
            auto const& s = shares[level];
            if (s.values == 0 || s.unique == 0 || s.unique_two == 0) continue;
            SeedExpectation e;
            e.level = static_cast<uint32_t>(level);
            e.unique = static_cast<double>(s.unique) / static_cast<double>(s.values) / unique_share;
            e.unique_two = static_cast<double>(s.unique_two) / static_cast<double>(s.values) / two_share;
            expectations[level] = e;
            os << (os.tellp() > 0 ? ", " : "") << "level " << level << " (C " << (size_t{1} << level) << ") " << std::fixed
               << std::setprecision(3) << e.unique << " / " << e.unique_two;
        }
        double const seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
        log << "Seed expectation: the unique share of a grid level's values over the whole index's (ZU / ZT): " << os.str()
            << std::fixed << std::setprecision(2) << " (" << seconds << " s)" << std::defaultfloat << std::endl;
        return expectations;
    }

    // --stage1_probe: the sample's tallies (stage1_probe::Global(), joined from the alignment threads) into
    // misc/<sample>_stage1_probe.tsv (the hash sample's rules), _stage1_probe_stride.tsv (the stride pick's),
    // _stage1_probe_taxa.tsv (both; the taxa named by their lineage, as the profile names them, if the taxonomy is loaded)
    // and, for short reads, _first_lookups_probe.tsv, with a log line of a few rules; then cleared for the next sample.
    static void WriteStage1Probe(Options const& options, int index, ProtalDB& db, bool short_reads) {
        namespace probe = stage1_probe;
        auto& tally = probe::Global();
        std::filesystem::path const dir = options.GetMiscOutputDir();
        std::error_code ec;
        std::filesystem::create_directories(dir, ec);
        std::string const prefix = (dir / options.GetSampleId(index)).string();
        {
            std::ofstream os(prefix + "_stage1_probe.tsv");
            tally.WriteRules(os);
        }
        {
            std::ofstream os(prefix + "_stage1_probe_stride.tsv");
            tally.WriteRules(os, probe::Pick::stride);
        }
        {
            std::ofstream os(prefix + "_stage1_probe_taxa.tsv");
            bool const named = db.IsTaxonomyLoaded();
            tally.WriteTaxa(os, [&](uint32_t taxid) { return named ? db.GetTaxonomy().LineageStr(static_cast<int>(taxid)) : std::string("-"); });
        }
        if (short_reads) {
            std::ofstream os(prefix + "_first_lookups_probe.tsv");
            tally.WriteFirstLookups(os);
        }
        auto score = [](uint8_t t) { return static_cast<size_t>(std::find(probe::kMinScores.begin(), probe::kMinScores.end(), t) - probe::kMinScores.begin()); };
        std::ostringstream line;
        line << "Sample " << options.GetSampleId(index) << " stage-1 probe: " << tally.Reads() << (short_reads ? " reads (pairs as one)" : " reads")
             << ", " << tally.Aligned() << " with records; " << probe::RuleText(tally, 2, score(14), 0) << "; "
             << probe::RuleText(tally, 3, score(14), 0) << "; " << probe::RuleText(tally, 4, score(14), 0) << "; "
             << probe::RuleText(tally, 3, score(15), 1);
        for (size_t level = 2; level <= 4; level++) line << "; " << probe::RuleText(tally, level, score(14), 0, probe::Pick::stride);
        if (short_reads && tally.FirstLookupMates() > 0) {
            auto const give_up = tally.GiveUp(2, 11);  // the first 4 lookups
            line << std::fixed << std::setprecision(2) << "; mates whose first 4 lookups score 11 or less: " << 100.0 * give_up[0]
                 << "% of the mates, " << 100.0 * give_up[1] << "% of those that align";
        }
        std::cout << line.str() << " (" << prefix << "_stage1_probe*.tsv)" << std::endl;
        tally.Clear();
    }

    // Aligns the samples. `sample_done` (if any) is called with a sample's index once its SAM file is complete, or
    // when the sample is skipped because its SAM exists (never one whose reads are piped: ReadsArePiped): Run
    // profiles such samples while the next ones are aligned (ProfilingAhead).
    template<typename AlignmentBenchmark=NoBenchmark>
    static void RunWrapper(Options& options, ProtalDB& db, AlignmentBenchmark benchmark=NoBenchmark{},
                           std::function<void(size_t)> const& sample_done = {}) {

        const size_t mmer_size = 15;
//        const size_t kmer_size = 27;
        const size_t kmer_size = 31;

        if (options.BuildMode()) {
            RunBuild(options, db);
        } else {

            // The index: loaded by a start-up thread beside the genome preload and the tables (Run), or here
            // (--sequential_load). Either way its messages and "Load Index took" are printed now.
            std::unique_ptr<Seedmap> index;
            if (auto loaded = db.TakeIndex()) {
                std::cout << loaded->log;
                index = std::move(loaded->map);
                loaded->timer.PrintResults();
            } else {
                Benchmark bm_load_index("Load Index");
                bm_load_index.Start();
                index = LoadIndex(options, db.GetGenomes().IndexFieldMaxima(), db.GetGenomes().MapFingerprint(), std::cout);
                bm_load_index.Stop();
                bm_load_index.PrintResults();
            }
            // The run's peak is usually here: the index, its loading threads' buffers and the preloaded genes at once.
            // Near the end of the load the threads stop one by one (zstd::LoadBudget), their buffers freed.
            if (index->ChunksOnFewerThreads() > 0) {
                std::cout << "Index: its last " << index->ChunksOnFewerThreads() << " chunks decoded on fewer threads, so that "
                          << "their buffers stay within the index's memory and what aligning takes beside it" << std::endl;
            }
            protal::build::PrintMemory("loading the index");
            Seedmap& map = *index;

            // Seeds must be sampled exactly as when the index was built.
            ClosedSyncmer minimizer{mmer_size, 7, 2, map.UsesFullSyncmerMask()};
            SimpleKmerHandler iterator{kmer_size, mmer_size, minimizer};


            KmerLookupSM kmer_lookup(map, options.GetMaxKeyUbiquity());

            GenomeLoader& genomes = db.GetGenomes();

            WFA2Wrapper2 aligner(4, 6, 2, options.GetXDrop(), options.GetWfaTraceback());
            // Long reads are aligned without X-drop: over their gene-long windows it loses the own
            // species' alignment of genes the read ends in, and leaves a relative's weaker one unique.
            WFA2Wrapper2 long_read_wfa(4, 6, 2, 0, options.GetWfaTraceback());

            // Every taxon's genus, for the adaptive candidates of short reads (--adaptive_candidates); the taxonomy is loaded
            // for profiling, and here if a run only aligns (not for long reads alone).
            std::shared_ptr<std::vector<uint32_t> const> genera;
            if (options.GetAdaptiveCandidates() > 0 && (options.AnySample(ReadType::Paired) || options.AnySample(ReadType::Single))) {
                if (!db.IsTaxonomyLoaded()) db.LoadTaxonomy(options.TaxonomyDbFile());
                genera = TaxonomyGenera(db.GetTaxonomy());
            }

            if (options.PreloadGenomes() && !genomes.AllGenomesLoaded()) {
                Benchmark bm_preload_genomes("Preload genomes");
                bm_preload_genomes.Start();
                genomes.LoadAllGenomes(static_cast<int>(options.GetThreads()));
                bm_preload_genomes.Stop();
                bm_preload_genomes.PrintResults();
            }

//            using AnchorFinder = SimpleAnchorFinder<KmerLookupSM>;
//            using AnchorFinder = HashMapAnchorFinder<KmerLookupSM>;
//            using AnchorFinder = ListAnchorFinder<KmerLookupSM>;
            using AnchorFinder = ChainAnchorFinder<KmerLookupSM>;
//            using AnchorFinder = NaiveAnchorFinder<KmerLookupSM>;


            // misc/cpu.tsv: the index (loaded beside the start-up, waited for above) and the genome preload get a row of
            // their own, so that a sample's row is its alignment only.
            CpuLog::Get().Mark("index and genomes");
            // The grid pick on the whole index: its samples gated by the grid, then seeded from every entry (GridWhole).
            bool const grid_whole = GridWhole(options);
            if (AnyGridSample(options)) {
                std::cout << "Grid pick: ";
                if (grid_whole) {
                    std::cout << "the whole index is loaded (--sparse_index off, or samples not all grid-picked), so the grid only gates "
                              << "its samples: the reads kept are seeded from every core, a gene anchored on two seeds";
                } else {
                    std::cout << "the reads kept are seeded from their grid entries alone"
                              << (options.NoSingleSeedAnchors() ? "" : ", a short read's gene anchored on one seed too")
                              << (options.SparseIndex() == "emulate" ? " (--sparse_index emulate: on the whole index, as the sparse index seeds them)" : "");
                }
                std::cout << std::endl;
            }
            // The seed expectation of the grid levels a sample is seeded from alone (a sparse index, --sparse_index emulate, or
            // --sparse_seeding with the grid pick), for its SAM's header (SeedExpectations). None with the database's unique
            // cores (UniqueCores.h): ZS, ZU and ZT count the cores each alignment covers, as the whole index's seeds would.
            std::vector<std::optional<SeedExpectation>> seed_expectations;
            if (!genomes.UniqueCores()) {
                auto const range = options.GetRange();
                bool const grid_seeded = std::any_of(range.begin(), range.end(), [&](size_t i) {
                    return options.GetSparseFirst(i) > 1 && options.GetSparseFirstPick(i) == SparseFirstPick::Grid &&
                           (!grid_whole || options.SparseSeeding());
                });
                if (grid_seeded) seed_expectations = SeedExpectations(map, genomes, static_cast<int>(options.GetThreads()), std::cout);
            }

            Benchmark bm_classify("Processing all samples");
            bm_classify.Start();
            stage1_probe::EnabledFlag() = options.Stage1Probe();  // read by Classify's loops
            bool long_genes_told = false;  // the warning of LongReadAligner::ChunksHoldEveryGene, once
            // A .sam.zst takes its records straight into the SAM, behind room for its header: for up
            // to all the database's genes, and in proportion to the read files (sam_zstd::HeaderRoom).
            size_t const database_genes = genomes.GeneCount();
            auto input_bytes = [](std::vector<std::string> const& files) {
                uint64_t total = 0;
                for (auto const& file : files) {
                    std::error_code ec;
                    if (!std::filesystem::is_regular_file(file, ec)) return uint64_t{0};  // a pipe: unknown
                    total += std::filesystem::file_size(file, ec);
                    if (ec) return uint64_t{0};
                }
                return total;
            };

            for (auto index : options.GetRange()) {
                auto const read_type = options.GetReadType(index);
                bool const single_file = read_type != ReadType::Paired;

                Benchmark bm_classify_sample("Aligning reads");
                bm_classify_sample.Start();

                // TODO implement logger in protal
                auto [sam, compressed] = options.SamFile(index);
                auto [sam_plain, _] = options.SamFile(index, true);

                if (!MakeFolder(std::filesystem::path(sam).parent_path())) {
                    std::cout << "Cannot create directories for this path " << sam << std::endl;
                    exit(32);
                }


                // std::cout << index << " Process sample " << options.GetSampleId(index) << (std::filesystem::exists(sam) ? " (sam exists)" : " (sam does not exist)") << std::endl;

                // Avoid aligning files that already exist (but reads that come through a pipe come once).
                if (!options.Force() && !ReadsArePiped(options, index) &&
                    (std::filesystem::exists(sam) || std::filesystem::exists(sam_plain))) {
                    std::cout << "Skip " << sam << " continue" << std::endl;
                    // Profile the file that is there: an earlier run may have left it uncompressed.
                    if (compressed && !std::filesystem::exists(sam)) options.UseUncompressedSamFile(index);
                    if (sample_done) sample_done(index);
                    continue;
                }
                // Its row of misc/cpu.tsv when this iteration ends, however it ends (CpuLog.h).
                CpuLog::Stage cpu_stage("aligning " + options.GetSampleId(index));
                if (options.Stage1Probe()) stage1_probe::Global().Clear();  // this sample's tallies only

                // AnchorFinder. A long read is seeded from all its k-mers: stopping at -s seeds, as for a
                // short read, would leave most of it unseeded.
                AnchorFinder anchor_finder(kmer_lookup, mmer_size, options.GetMinSuccessfulLookups(),
                                           IsLongReadType(read_type) ? SIZE_MAX : options.GetMaxSeedSize(), genomes);
                // --sparse_first: reads judged by a sample of their lookups first (ChainAnchorFinder::SetSparseFirst; long
                // reads over all their chunks, LongReadAligner), as the sample's read type is gated (or the map's
                // SPARSE_FIRST says).
                if (size_t const sparse_first = options.GetSparseFirst(index); sparse_first > 0) {
                    // The grid pick: every core looked up, the entries of the gate's grid level counting. On a sparse index
                    // (or --sparse_index emulate) a read kept is seeded from those alone, a short read's gene anchoring on a
                    // single seed too (--no_single_seed_anchors: not); on the whole index the grid only gates, and a read
                    // kept is seeded from every entry, on two seeds a gene, as other runs (GridWhole).
                    bool const grid = options.GetSparseFirstPick(index) == SparseFirstPick::Grid;
                    anchor_finder.SetSparseFirst(static_cast<uint8_t>(__builtin_ctzll(sparse_first)),
                                                 static_cast<uint32_t>(options.GetSparseFirstScore()),
                                                 options.GetSparseFirstPick(index) == SparseFirstPick::Stride, grid, grid_whole);
                    // --sparse_seeding: the reads kept seeded from those cores alone, as an index of only them would.
                    if (options.SparseSeeding()) anchor_finder.SetSeedLevel(static_cast<uint8_t>(__builtin_ctzll(sparse_first)));
                    anchor_finder.SetSingleSeedAnchors(sparse_first > 1 && grid && !grid_whole && !IsLongReadType(read_type) &&
                                                       !options.NoSingleSeedAnchors());
                }
                // AlignmentHandler approach
                double const max_score_ani = options.GetMaxScoreAni(read_type);
                SimpleAlignmentHandler alignment_handler(genomes, IsLongReadType(read_type) ? long_read_wfa : aligner, kmer_size,
                                                         options.GetAlignTop(), max_score_ani, options.FastAlign());
                // Reads from their anchors' exact matches (AnchoredAligner); long reads, whose seeds lie on
                // diagonals that their indels shift, through every link of the chain.
                alignment_handler.SetAnchoredAlignment(!options.WholeReadAlignment());
                alignment_handler.SetAnchoredIndels(IsLongReadType(read_type));
                // ZS, ZU and ZT from the gene's cores each alignment covers with exact k-mers (UniqueCores.h), found with the
                // run's k-mer handler (the build's windows), where the database has unique_cores.bin.
                if (auto const cores = genomes.UniqueCores()) alignment_handler.SetUniqueCores(cores, iterator);
                // Before WFA2, a candidate that cannot align within the budget is refused: by the k-mer screen, or for ONT
                // reads by their indel distance (IndelBound), as the screen refuses next to nothing at their identity floor
                // of 0.85 (--no_indel_bound: the screen, as before 2026-10-09). PacBio reads keep the screen: the candidates
                // it passes mostly align, and the bound on them cost more than it saved.
                bool const indel_bound = read_type == ReadType::ONT && !options.NoAlignmentScreen() && !options.NoIndelBound();
                alignment_handler.SetAlignmentScreen(!options.NoAlignmentScreen() && !indel_bound);
                alignment_handler.SetIndelBound(indel_bound);
                if (!IsLongReadType(read_type)) alignment_handler.SetAdaptiveCandidates(options.GetAdaptiveCandidates(), genera);
                if (!IsLongReadType(read_type)) alignment_handler.SetAlleleScores(!options.NoAlleleScores());
                if (!IsLongReadType(read_type)) alignment_handler.SetWeightedShift(options.WeightedSiteShift());
                // --flex_indels: the sample's lookups (the gate's too) score their flex cells with an indel per flank allowed.
                anchor_finder.SetFlexIndels(options.FlexIndels(index));



                options.SetCurrentIndex(index);
                std::string const sam_partial = sam + ".partial";  // renamed by FinishSamFile
                // The read type, and how the sample was aligned (SamSettings.h; --profile_only reads it).
                std::string const read_type_line = kSamReadTypeComment + Info(read_type).token + '\n' +
                                                   sam_settings::SamSettingsLine(AlignmentSettings(options, index, genomes.UniqueCores() != nullptr));
                // The alignment threads compress the records as the name asks (.gz, .zst). The header
                // lists the genes that the records name, so it is written once they are complete;
                // --full_sam_header lists every gene of the database, written first.
                std::optional<std::string> full_header;
                if (options.FullSamHeader()) {
                    std::ostringstream os;
                    genomes.WriteSamHeader(os);
                    os << read_type_line;
                    full_header = os.str();
                }
                // A sample's reads: one file, or several (its lanes) read one after the other; paired-end reads, as many
                // second files, paired with the first by position.
                auto const first_files = options.GetFirstFiles(index);
                auto const second_files = single_file ? std::vector<std::string>() : options.GetSecondFiles(index);
                std::vector<std::string> read_file_list = first_files;
                read_file_list.insert(read_file_list.end(), second_files.begin(), second_files.end());
                SamOutput sam_output(sam_partial, SamCompressionOf(sam), full_header,
                                     sam_zstd::HeaderRoom(database_genes, input_bytes(read_file_list), options.WriteUnmappedReads(index) ? 0 : genomes.GetGenomeMap().size()),
                                     options.ScratchFolder());
                sam_output.SetUnmappedRecords(options.WriteUnmappedReads(index));
                if (!sam_output.Ok()) {
                    RunStatus::Get().Fail("Cannot write the SAM file of sample " + options.GetSampleId(index) + ": " + sam_output.Error());
                    continue;
                }
                std::cout << "Align the " << ReadTypeName(read_type) << " reads of sample "
                          << options.GetSampleId(index) << " (-a " << max_score_ani << ")" << std::endl;

                // A gzip read file inflates on several threads by -t (ThreadedGzStream.h): a BGZF one on 1 to 4 of protal's
                // per file (a zstd file's frames, and a gzip file's members if its second begins near its start, alike),
                // another gzip file (a regular one) in pieces on 1 to 8 (a thread per 4 of -t and file, ParallelGzip.h).
                std::size_t const files = single_file ? 1 : 2;
                ThreadedGzStreambuf::SetInflateThreads(ThreadedGzStreambuf::InflateThreadsFor(options.GetThreads(), files));
                ThreadedGzStreambuf::SetGzipThreads(ThreadedGzStreambuf::GzipThreadsFor(options.GetThreads(), files));

                // Main Run Call. This is where the reads are read and alignment happens
                bool truncated = false;
                bool read_success = true;
                std::string read_problem;   // why reading failed, if it says
                size_t reads_read = 0;
                size_t bases_read = 0;  // for the header's ScannedReadsLine
                size_t placeholder_qualities = 0;  // reads (mates) with Q0 at every base, given the read type's quality
                std::vector<ThreadedGzStreambuf::FileRead> files_read;  // how each file was read, for the input line
                double reader_seconds = 0;                              // the alignment threads' time getting reads
                size_t reading_threads = 0;                             // the threads reading them beside the alignment's
                // A stream's reading threads: its files', one file at a time.
                auto stream_threads = [](ThreadedGzIstream& is) {
                    unsigned threads = 0;
                    for (auto const& file : is.rdbuf()->files_read()) threads = std::max(threads, file.threads);
                    return size_t{ threads };
                };
                std::string const read_files = options.ReadFilesText(index);
                // An input that cannot be opened reads as empty: say why instead.
                auto cannot_read = [&](ThreadedGzIstream& is) {
                    if (is.rdbuf()->is_open()) return false;
                    RunStatus::Get().Fail("Cannot read " + is.rdbuf()->current_path() + " (sample " + options.GetSampleId(index) + "): " +
                                          is.rdbuf()->open_error() + "; no SAM file was written");
                    sam_output.Discard();
                    return true;
                };
                auto read_error = [](ThreadedGzIstream& is) {
                    auto const message = is.rdbuf()->read_error_message();
                    return message.empty() ? std::string() : ": " + message;
                };
                // With --verbose, how each input is inflated as it is opened: on several threads (BGZF blocks or a zstd file's
                // frames on protal's, another gzip file in pieces), any other on one.
                ThreadedGzStreambuf::OnOpen tell_input;
                if (options.Verbose()) {
                    tell_input = [](std::string const& path, ThreadedGzStreambuf const& buf) {
                        std::string const threads = std::to_string(buf.inflate_threads()) + " thread(s)";
                        std::cout << "Input " << path << ": "
                                  << (buf.opened_as_bgzf()          ? "BGZF, inflated on " + threads
                                      : buf.opened_as_zstd_frames() ? "zstd, its frames decompressed on " + threads
                                      : buf.opened_as_gzip_members() ? "gzip of several members (not BGZF), its members inflated on " + threads
                                      : buf.parallel_gzip()         ? "gzip (not BGZF), inflated in pieces on " + threads
                                                                    : std::string("gzip, zstd or uncompressed, read on one thread"))
                                  << (buf.parallel_gzip_error().empty() ? std::string()
                                                                        : " (not in pieces: " + buf.parallel_gzip_error() + ")")
                                  << std::endl;
                    };
                }
                if (IsLongReadType(read_type)) {
                    ThreadedGzIstream is { first_files, tell_input };
                    if (cannot_read(is)) continue;
                    SeqReaderSE reader{ is, FastaQualityChar(read_type) };
                    LongReadAligner<SimpleKmerHandler<ClosedSyncmer>, AnchorFinder> long_read_aligner(
                            iterator, anchor_finder, alignment_handler, genomes, options.GetAlignTop(), max_score_ani);
                    long_read_aligner.SetBudgetMargin(options.GetLongReadBudget());
                    long_read_aligner.SetSparseFirstHits(static_cast<uint32_t>(options.GetSparseFirstHits()));
                    if (!long_read_aligner.ChunksHoldEveryGene() && !long_genes_told) {
                        long_genes_told = true;
                        std::cerr << "Warning: the database's longest gene has " << genomes.MaxGeneLength() << " bp; with its margins it "
                                  << "exceeds half of the " << kMaxLongReadChunk << " bp chunks reads longer than that are seeded in, "
                                  << "so genes that long may be missed in such reads" << std::endl;
                    }
                    ProtalLongReadOutputHandler output_handler(sam_output, options.GetMaxOut(), 1024*1024*16, genomes, 0.8);
                    auto protal_stats = protal::classify::RunLongReads(reader, options, long_read_aligner, output_handler);
                    if (options.Verbose()) {
                        protal_stats.WriteStats();
                    }
                    PrintAlignmentCounts(options, index, read_type, protal_stats, long_read_aligner.GetAlignmentHandler(),
                                         long_read_aligner.GetAnchorFinder().m_seeding);
                    reads_read = protal_stats.reads;
                    bases_read = protal_stats.bases;
                    reader_seconds = protal_stats.reader_seconds;
                    // A truncated or corrupt gzip file reads as one that ends early (ThreadedGzStream).
                    truncated = is.rdbuf()->read_failed();
                    read_problem = read_error(is);
                    read_success = reader.Success();
                    placeholder_qualities = reader.PlaceholderQualities();
                    is.close();
                    files_read = is.rdbuf()->files_read();
                    reading_threads = stream_threads(is);
                } else if (read_type == ReadType::Single) {
                    ThreadedGzIstream is { first_files, tell_input };
                    if (cannot_read(is)) continue;
                    // Reads longer than short reads are long reads given as single-end: the check
                    // before aligning sees only the first reads of the file.
                    SeqReaderSE reader{ is, FastaQualityChar(read_type), MAX_SHORT_READ_LENGTH };
                    auto align = [&](auto output_handler) {
                        return protal::classify::RunSingleEnd<
                                SimpleKmerHandler<ClosedSyncmer>,
                                AnchorFinder,
                                SimpleAlignmentHandler,
                                decltype(output_handler),
                                DEBUG_NONE,
                                AlignmentBenchmark>(
                                reader, options, anchor_finder, alignment_handler, output_handler, iterator, benchmark);
                    };
                    auto protal_stats = options.GetMAPQDebugOut() ?
                            align(ProtalSingleOutputHandler<true>(sam_output, options.GetMaxOut(), 1024*512, 1024*1024*16, genomes, 0.8)) :
                            align(ProtalSingleOutputHandler<false>(sam_output, options.GetMaxOut(), 1024*512, 1024*1024*16, genomes, 0.8));
                    if (options.Verbose()) {
                        protal_stats.WriteStats();
                    }
                    PrintAlignmentCounts(options, index, read_type, protal_stats, alignment_handler, anchor_finder.m_seeding);
                    reads_read = protal_stats.reads;
                    bases_read = protal_stats.bases;
                    reader_seconds = protal_stats.reader_seconds;
                    // A truncated or corrupt gzip file reads as one that ends early (ThreadedGzStream).
                    truncated = is.rdbuf()->read_failed();
                    read_problem = read_error(is);
                    read_success = reader.Success();
                    placeholder_qualities = reader.PlaceholderQualities();
                    if (reader.TooLong() > 0) {
                        read_problem = ": read " + reader.TooLongId() + " has " + std::to_string(reader.TooLong()) +
                                       " bp, too long for short reads; give --read_type pb or ont (or in the map's READ_TYPE "
                                       "column) for PacBio or ONT reads";
                    }
                    is.close();
                    files_read = is.rdbuf()->files_read();
                    reading_threads = stream_threads(is);
                } else {
                    // Each file inflates in a thread of its own (ThreadedGzStream.h), outside the reader lock. Of a sample's
                    // lanes, the batches end with each lane's files, so a lane whose mates differ in number fails it.
                    ThreadedGzIstream is1 { first_files, tell_input };
                    ThreadedGzIstream is2 { second_files, tell_input };
                    if (cannot_read(is1) || cannot_read(is2)) continue;
                    SeqReaderPE reader{ is1, is2, FastaQualityChar(read_type) };
                    Statistics protal_stats;

                    if (options.GetMAPQDebugOut()) {
                        using OutputHandler = ProtalPairedOutputHandler<true>;

                        OutputHandler output_handler(sam_output, options.GetMaxOut(), 1024*512, 1024*1024*16, genomes, 0.8);
                        protal_stats = protal::classify::RunPairedEnd<
                                SimpleKmerHandler<ClosedSyncmer>,
                                AnchorFinder,
                                SimpleAlignmentHandler,
                                OutputHandler,
                                DEBUG_NONE,
                                AlignmentBenchmark>(
                                reader, options, anchor_finder, alignment_handler, output_handler, iterator, genomes, benchmark);
                    } else {
                        using OutputHandler = ProtalPairedOutputHandler<false>;

                        OutputHandler output_handler(sam_output, options.GetMaxOut(), 1024*512, 1024*1024*16, genomes, 0.8);
                        protal_stats = protal::classify::RunPairedEnd<
                                SimpleKmerHandler<ClosedSyncmer>,
                                AnchorFinder,
                                SimpleAlignmentHandler,
                                OutputHandler,
                                DEBUG_NONE,
                                AlignmentBenchmark>(
                                reader, options, anchor_finder, alignment_handler, output_handler, iterator, genomes, benchmark);
                    }
                    if (options.Verbose()) {
                        protal_stats.WriteStats();
                    }
                    PrintAlignmentCounts(options, index, read_type, protal_stats, alignment_handler, anchor_finder.m_seeding);
                    reads_read = protal_stats.reads;
                    bases_read = protal_stats.bases;
                    reader_seconds = protal_stats.reader_seconds;
                    // A truncated or corrupt gzip file reads as one that ends early (ThreadedGzStream).
                    truncated = is1.rdbuf()->read_failed() || is2.rdbuf()->read_failed();
                    read_problem = read_error(is1);
                    if (read_problem.empty()) read_problem = read_error(is2);
                    read_success = reader.Success();
                    placeholder_qualities = reader.PlaceholderQualities();
                    if (!read_success && read_problem.empty() && reader.MatesOutOfStep()) {
                        read_problem = first_files.size() > 1 ? ": the two files of one of its lanes hold different numbers of reads"
                                                              : ": the two files hold different numbers of reads";
                    }
                    if (reader.NameMismatches() > 0 && read_success && !truncated) {
                        std::cerr << "Warning: sample " << options.GetSampleId(index) << ": the mates of " << reader.NameMismatches()
                                  << " read pair(s) have different names (e.g. " << reader.FirstNameMismatch()
                                  << "); are " << read_files << " in the same order?" << std::endl;
                    }
                    is1.close();
                    is2.close();
                    files_read = is1.rdbuf()->files_read();
                    files_read.insert(files_read.end(), is2.rdbuf()->files_read().begin(), is2.rdbuf()->files_read().end());
                    reading_threads = stream_threads(is1) + stream_threads(is2);
                }
                bm_classify_sample.Stop();
                bm_classify_sample.PrintResults();
                // How the sample's files were read, and how long the alignment threads spent getting reads: a share of a
                // quarter or more of the alignment's time means the input held it back (BGZF reads fastest; running.md).
                {
                    std::cout << "Sample " << options.GetSampleId(index) << " input: ";
                    for (size_t f = 0; f < files_read.size(); f++) {
                        std::cout << (f ? "; " : "") << std::filesystem::path(files_read[f].path).filename().string() << ' '
                                  << files_read[f].how;
                    }
                    double const aligning = bm_classify_sample.Seconds();
                    double const per_thread = reader_seconds / static_cast<double>(std::max<size_t>(options.GetThreads(), 1));
                    double const share = aligning > 0 ? per_thread / aligning : 0;
                    size_t const cpus = AvailableCpus();
                    std::cout << std::fixed << std::setprecision(2) << "; the " << options.GetThreads()
                              << " alignment threads spent " << per_thread << " s each getting reads, "
                              << std::setprecision(0) << 100 * share << "% of the alignment's " << std::setprecision(2)
                              << aligning << " s" << (share >= 0.25 ? ": the input held the alignment back" : "")
                              << std::defaultfloat << "; " << options.GetThreads() << " alignment and " << reading_threads
                              << " reading threads on " << cpus << " CPUs"
                              << (options.GetThreads() + reading_threads > cpus ? ", more threads than CPUs" : "") << std::endl;
                }
                if (options.Stage1Probe()) WriteStage1Probe(options, index, db, !IsLongReadType(read_type));

                if (truncated) {
                    RunStatus::Get().Fail("The FASTQ files of sample " + options.GetSampleId(index) + " are truncated or corrupt (" +
                                          read_files + read_problem + "); no SAM file was written");
                    sam_output.Discard();
                    continue;
                }
                if (!read_success) {
                    RunStatus::Get().Fail("Reading the FASTQ files of sample " + options.GetSampleId(index) + " failed (" +
                                          read_files + read_problem + "); no SAM file was written");
                    sam_output.Discard();
                    continue;
                }
                if (reads_read == 0) {
                    std::cerr << "Warning: sample " << options.GetSampleId(index) << " has no reads (" << read_files << ")" << std::endl;
                }
                // Q0 at every base means unknown qualities (IsPlaceholderQuality): taken as qualities, every base would count
                // as an error, and the reads' divergence beyond their errors (ReadExcess) would be far from the models'.
                if (placeholder_qualities > 0) {
                    std::cerr << "Note: sample " << options.GetSampleId(index) << ": " << placeholder_qualities << " of "
                              << (single_file ? reads_read : 2 * reads_read) << (single_file ? " reads" : " mates")
                              << " have Q0 at every base, which means unknown qualities (pbsim3, some BAM-to-FASTQ exports), so "
                              << "their bases were given Q" << Info(read_type).fasta_quality << ", as " << ReadTypeName(read_type)
                              << " reads without qualities (FASTA) are" << std::endl;
                }
                Benchmark bm_finish_sam("Writing the SAM header and file");
                bm_finish_sam.Start();
                auto const header_start = std::chrono::steady_clock::now();
                std::string header;
                size_t header_genes = 0;
                size_t failed_taxa = 0;
                std::string failed_line;
                if (!full_header) {
                    std::ostringstream os;
                    auto const genes = sam_output.Genes();
                    header_genes = genes.size();
                    genomes.WriteSamHeader(os, genes);
                    os << read_type_line;
                    // The reads read, which the profile's composition compares with the reads the called species explain.
                    os << ScannedReadsLine({ reads_read, single_file ? reads_read : 2 * reads_read, bases_read });
                    // How its unique seeds compare with the genes' unique values, where it was seeded from a grid level alone.
                    if (uint8_t const level = anchor_finder.GridSeedLevel(); level > 0 && level < seed_expectations.size() &&
                                                                           seed_expectations[level]) {
                        os << SeedExpectationLine(*seed_expectations[level]);
                    }
                    // The reads that aligned nowhere, counted per taxon they seeded on (--write_unmapped_reads: a record each).
                    if (!options.WriteUnmappedReads(index)) {
                        auto const failed = sam_output.FailedCandidates();
                        failed_taxa = static_cast<size_t>(std::count_if(failed.begin(), failed.end(), [](uint64_t n) { return n > 0; }));
                        failed_line = FailedCandidatesLine(failed);
                        os << failed_line;
                    }
                    header = os.str();
                }
                double const header_seconds = gene_table::SecondsSince(header_start);
                bool const written = sam_output.Finish(header);
                bm_finish_sam.Stop();
                bm_finish_sam.PrintResults();
                {
                    // Where that time went (one thread): the header's text, and how Finish placed it: in the room before
                    // the records, or with the records copied behind it (no room, or a header larger than the room).
                    auto const& placed = sam_output.Placement();
                    std::ostringstream line;
                    line << std::fixed << std::setprecision(2) << "SAM header: " << (full_header ? std::string("every gene of the database")
                                                                                                : std::to_string(header_genes) + " genes")
                         << ", " << header.size() << " bytes, written as " << placed.header_bytes << " (text built in " << header_seconds << " s); ";
                    if (failed_taxa > 0) line << "the unaligned reads' failed candidates of " << failed_taxa << " taxa in " << failed_line.size() << " bytes of it; ";
                    if (placed.copied) {
                        line << "the records (" << placed.copied_bytes << " bytes) were copied behind it in " << placed.copy_seconds << " s"
                             << (placed.room > 0 ? ", as it was larger than its room of " + std::to_string(placed.room) + " bytes" : std::string());
                    } else {
                        line << (placed.room > 0 ? "placed in its room of " + std::to_string(placed.room) + " bytes before the records"
                                                 : std::string("written before the records"));
                    }
                    std::cout << line.str() << std::endl;
                }
                if (!written) {
                    RunStatus::Get().Fail("Writing the SAM file of sample " + options.GetSampleId(index) + " failed: " + sam_output.Error());
                    continue;
                }
                FinishSamFile(options, index, sam_partial);
                if (sample_done) sample_done(index);
            }
            bm_classify.Stop();
            bm_classify.PrintResults();
            if (options.GetWfaTraceback() == WfaTraceback::Singletrack) {
                std::cout << "WFA2 Singletrack: " << WFA2Wrapper2::SingletrackFallbacks() << " alignment(s) traced back classic, "
                          << "their Singletrack backtrace having found no path through the wavefronts wf-adaptive pruned" << std::endl;
            }
            protal::build::PrintMemory("aligning");
        }
    }

    void RunAlignment(Options& options, ProtalDB& db, std::function<void(size_t)> const& sample_done) {
        using AlignmentBenchmark = CoreBenchmark;
        // Untangle Template options that need to be written out specifically.
        if (options.BenchmarkAlignment()) {
            AlignmentBenchmark alignment_benchmark{};
            if (!options.GetBenchmarkAlignmentOutputFile().empty()) {
                alignment_benchmark.SetOutput(options.GetBenchmarkAlignmentOutputFile());
            }
            RunWrapper(options, db, alignment_benchmark, sample_done);
            if (!options.GetBenchmarkAlignmentOutputFile().empty()) {
                alignment_benchmark.DestroyOutput();
            }
        } else {
            RunWrapper(options, db, NoBenchmark{}, sample_done);
        }
    }
}
