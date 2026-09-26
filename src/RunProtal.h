#pragma once

#include "Profiler.h"
#include <iostream>
#include "Benchmark.h"
#include "Options.h"
#include "Build.h"
#include "Alignment/WFA2Wrapper2.h"
#include "Classify.h"
#include "ChainAnchorFinder.h"
#include "Taxonomy.h"
#include "gzstream.h"
#include "AlignmentStrategy.h"
#include "TaxonStatisticsOutput.h"
#include "ProgressBar.h"
#include "protal_config.h"
#include "Compressor.h"
#include "RunStatus.h"

#include <atomic>
#include <iomanip>
#include <regex>
#include <ranges>
#include <unistd.h>

// #include "Profiler/ReadFilter.h"

#include <algorithm>
#include <iterator>

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


    class ProtalDB {
        GenomeLoader m_genomes;
        std::optional<taxonomy::IntTaxonomy> m_taxonomy;

    public:
        ProtalDB(std::string sequence_file, std::string map_file) :
                m_genomes(sequence_file, map_file),
                m_taxonomy() {
        }

        ProtalDB(std::string sequence_file, std::string map_file, std::string hittable_genes_file, std::string unique_kmers_file) :
                m_genomes(sequence_file, map_file),
                m_taxonomy() {
            // m_genomes.LoadHittableGenes(hittable_genes_file);
            m_genomes.LoadUniqueKmers(unique_kmers_file);
        }

        void LoadTaxonomy(std::string file) {
            m_taxonomy = taxonomy::IntTaxonomy(file);
        }

        bool IsTaxonomyLoaded() const {
            return m_taxonomy.has_value();
        }

        taxonomy::IntTaxonomy& GetTaxonomy() {
            if (!m_taxonomy.has_value()) {
                std::cerr << "Taxonomy has not been loaded. Error in code. " << std::endl;
                exit(12);
            }
            return m_taxonomy.value();
        }

        GenomeLoader& GetGenomes() {
            return m_genomes;
        }

        std::optional<taxonomy::IntTaxonomy> GetTaxonomyOptional() {
            return m_taxonomy;
        }
    };

    // Alignments are written under a temporary name ("<sam>.partial") and only get their final name
    // once complete, so an interrupted run never leaves a truncated SAM that a rerun would skip and
    // reuse. This moves a finished file into place, compressing it first when a .gz name was asked for;
    // if compression fails, the alignments are kept uncompressed (and the run reports the failure).
    static void FinishSamFile(Options& options, int index, std::string const& partial) {
        namespace fs = std::filesystem;
        auto [sam, gzipped] = options.SamFile(index);
        auto [sam_nogzip, _] = options.SamFile(index, true);
        auto const sample = options.GetSampleId(index);
        std::error_code ec;

        if (gzipped) {
            try {
                Compressor::compressInPlace(partial, options.GetThreads());
                fs::rename(partial + ".gz", sam, ec);
                if (ec) RunStatus::Get().Fail("Cannot move " + partial + ".gz to " + sam + ": " + ec.message());
                return;
            } catch (const std::exception& e) {
                fs::remove(partial + ".gz", ec);
                if (!fs::exists(partial)) {
                    RunStatus::Get().Fail("Compressing the SAM file of sample " + sample + " failed: " + e.what());
                    return;
                }
                RunStatus::Get().Fail("Compressing the SAM file of sample " + sample + " failed (" + e.what() +
                                      "); kept it uncompressed as " + sam_nogzip);
                options.SetSamFileGzip(index, false);
            }
        }
        fs::rename(partial, sam_nogzip, ec);
        if (ec) RunStatus::Get().Fail("Cannot move " + partial + " to " + sam_nogzip + ": " + ec.message());
    }

    // The index is built from --reference, but queries align reads against the database's reference.fna,
    // found through reference.map. Every --reference record must therefore be a reference.map gene of
    // the same length, with ids and length inside the index's 20-bit fields.
    static void CheckReferenceAgainstMap(Options const& options, GenomeLoader& genomes) {
        auto input = protal::build::OpenInput(options.GetSequenceFilePath());  // raw or .zst
        SeqReader reader{ input->Stream() };
        FastxRecord record;
        size_t records = 0, problems = 0;
        constexpr uint64_t max_id = (uint64_t{1} << SEEDMAP_TAXID_BITS) - 1;
        constexpr uint64_t max_gene = (uint64_t{1} << SEEDMAP_GENEID_BITS) - 1;
        constexpr uint64_t max_length = (uint64_t{1} << SEEDMAP_GENE_POS_BITS) - 1;
        while (reader(record)) {
            records++;
            std::string problem;
            size_t taxid = 0, geneid = 0;
            try {
                if (record.header.find('_') == std::string::npos) throw std::invalid_argument("no underscore");
                std::tie(taxid, geneid) = KmerUtils::ExtractHeaderInformation(record.header);
            } catch (std::exception const&) {
                problem = "header is not <taxid>_<gene id>";
            }
            if (problem.empty()) {
                if (taxid == 0 || taxid > max_id || geneid == 0 || geneid > max_gene) {
                    problem = "taxid and gene id must be between 1 and " + std::to_string(max_id);
                } else if (record.sequence.size() > max_length) {
                    problem = "gene is longer than " + std::to_string(max_length) + " bases";
                } else if (!genomes.HasGene(taxid, geneid)) {
                    problem = "not in " + options.GetSequenceMapFile();
                } else if (genomes.GeneLength(taxid, geneid) != record.sequence.size()) {
                    problem = std::to_string(record.sequence.size()) + " bases, but " +
                              std::to_string(genomes.GeneLength(taxid, geneid)) + " in reference.map";
                }
            }
            if (!problem.empty() && ++problems <= 10) {
                std::cerr << "--reference record " << record.header << ": " << problem << std::endl;
            }
        }
        if (records == 0) {
            std::cerr << "--reference " << options.GetSequenceFilePath() << " contains no sequences" << std::endl;
            exit(8);
        }
        if (problems > 0) {
            std::cerr << problems << " of " << records << " --reference records do not match the database's "
                      << "reference.map; build the index from the database's own reference.fna" << std::endl;
            exit(8);
        }
        if (records != genomes.GeneCount()) {
            std::cerr << "Warning: --reference has " << records << " sequences, reference.map lists "
                      << genomes.GeneCount() << " genes" << std::endl;
        }
    }

    template<typename AlignmentBenchmark=NoBenchmark>
    static void RunWrapper(Options& options, ProtalDB& db, AlignmentBenchmark benchmark=NoBenchmark{}) {

        const size_t mmer_size = 15;
//        const size_t kmer_size = 27;
        const size_t kmer_size = 31;

        if (options.BuildMode()) {
            CheckReferenceAgainstMap(options, db.GetGenomes());

            // New indexes compare whole s-mers (index format 2, recorded in the index header).
            ClosedSyncmer minimizer{mmer_size, 7, 2, true};
            SimpleKmerHandler iterator{kmer_size, mmer_size, minimizer};

            Benchmark bm_build("Run build");
            bm_build.Start();
            KmerPutterSM kmer_putter{};
            auto protal_stats = protal::build::Run<SimpleKmerHandler<ClosedSyncmer>, KmerPutterSM, DEBUG_NONE>(
                    options, kmer_putter, iterator);

            std::cout << "Check" << std::endl;
            protal::build::Check<SimpleKmerHandler<ClosedSyncmer>, KmerPutterSM, DEBUG_NONE>(
                    options, kmer_putter, iterator);

            // Last, as the build reads reference.fna until here.
            protal::build::CompressReference(options);

            bm_build.PrintResults();
            protal_stats.WriteStats(std::cout);

        } else {

            // Benchmark Load Time
            Benchmark bm_load_index("Load Index");
            bm_load_index.Start();
            // Load Index
            Seedmap map;
            std::string const index_file = options.ResolvedIndexFile();
            std::cout << "Load index " << index_file << std::endl;
            map.Load(index_file);
            std::cout << "Index features: " << map.FeatureDescription() << std::endl;
            if (map.HasReferenceFingerprint() &&
                !(map.GetReferenceFingerprint() == ReferenceFingerprint::Of(options.GetSequenceMapFile(), options.GetSequenceFile()))) {
                std::cerr << "index.prx was built against a different reference: " << options.GetSequenceMapFile()
                          << " or " << options.GetSequenceFile() << " changed since the index was built. Rebuild the "
                          << "index (--build) or restore the reference files it was built with." << std::endl;
                exit(8);
            }

            // Seeds must be sampled exactly as when the index was built.
            ClosedSyncmer minimizer{mmer_size, 7, 2, map.UsesFullSyncmerMask()};
            SimpleKmerHandler iterator{kmer_size, mmer_size, minimizer};
            bm_load_index.Stop();
            bm_load_index.PrintResults();


            KmerLookupSM kmer_lookup(map, options.GetMaxKeyUbiquity());

            GenomeLoader& genomes = db.GetGenomes();

            WFA2Wrapper2 aligner(4, 6, 2, options.GetXDrop());

            if (options.PreloadGenomes() && !genomes.AllGenomesLoaded()) {
                Benchmark bm_preload_genomes("Preload genomes");
                bm_preload_genomes.Start();
                genomes.LoadAllGenomes();
                bm_preload_genomes.Stop();
                bm_preload_genomes.PrintResults();
            }

//            using AnchorFinder = SimpleAnchorFinder<KmerLookupSM>;
//            using AnchorFinder = HashMapAnchorFinder<KmerLookupSM>;
//            using AnchorFinder = ListAnchorFinder<KmerLookupSM>;
            using AnchorFinder = ChainAnchorFinder<KmerLookupSM>;
//            using AnchorFinder = NaiveAnchorFinder<KmerLookupSM>;

            using OutputHandler = ProtalOutputHandler;




            Benchmark bm_classify("Processing all samples");
            bm_classify.Start();

            if (options.PairedMode()) {
                for (auto index : options.GetRange()) {
                //for (auto index = 0; index < options.GetFileCount(); index++) {

                    Benchmark bm_classify_sample("Aligning reads");
                    bm_classify_sample.Start();

                    // TODO implement logger in protal
                    auto [sam, gzipped] = options.SamFile(index);
                    auto [sam_nogzip, _] = options.SamFile(index, true);

                    auto dir = std::filesystem::path(sam).parent_path();
                    
                    if (!std::filesystem::create_directories(dir.string()) && !std::filesystem::exists(dir)) {
                        std::cout << "Cannot create directories for this path " << sam << std::endl;
                        exit(32);
                    };


                    // std::cout << index << " Process sample " << options.GetSampleId(index) << (std::filesystem::exists(sam) ? " (sam exists)" : " (sam does not exist)") << std::endl;

                    // Avoid aligning files that already exist.
                    if (!options.Force() && (std::filesystem::exists(sam) || std::filesystem::exists(sam_nogzip))) {
                        std::cout << "Skip " << sam << " continue" << std::endl;
                        // Profile the file that is there: an earlier run may have left it uncompressed.
                        if (gzipped && !std::filesystem::exists(sam)) options.SetSamFileGzip(index, false);
                        continue;
                    }

                    // AnchorFinder
                    AnchorFinder anchor_finder(kmer_lookup, mmer_size, options.GetMinSuccessfulLookups(), options.GetMaxSeedSize(), genomes);
                    // AlignmentHandler approach
                    SimpleAlignmentHandler alignment_handler(genomes, aligner, kmer_size, options.GetAlignTop(), options.GetMaxScoreAni(), options.FastAlign());



                    options.SetCurrentIndex(index);
                    std::string const sam_partial = sam_nogzip + ".partial";  // renamed by FinishSamFile
                    std::ofstream sam_output(sam_partial, std::ios::out);
                    if (!sam_output) {
                        RunStatus::Get().Fail("Cannot write the SAM file of sample " + options.GetSampleId(index) + ": " + sam_partial);
                        continue;
                    }
                    genomes.WriteSamHeader(sam_output);

                    igzstream is1 { options.GetFirstFile(index).c_str() };
                    igzstream is2 { options.GetSecondFile(index).c_str() };
                    SeqReaderPE reader{is1, is2};


                    // Main Run Call. This is where the reads are read and alignment happens
                    if (options.GetMAPQDebugOut()) {
                        using OutputHandler = ProtalPairedOutputHandler<true>;

                        OutputHandler output_handler(sam_output, options.GetMaxOut(), 1024*512, 1024*1024*16, genomes, 0.8);
                        auto protal_stats = protal::classify::RunPairedEnd<
                                SimpleKmerHandler<ClosedSyncmer>,
                                AnchorFinder,
                                SimpleAlignmentHandler,
                                OutputHandler,
                                DEBUG_NONE,
                                AlignmentBenchmark>(
                                reader, options, anchor_finder, alignment_handler, output_handler, iterator, genomes, benchmark);
                        if (options.Verbose()) {
                            protal_stats.WriteStats();
                        }
                    } else {
                        using OutputHandler = ProtalPairedOutputHandler<false>;

                        OutputHandler output_handler(sam_output, options.GetMaxOut(), 1024*512, 1024*1024*16, genomes, 0.8);
                        auto protal_stats = protal::classify::RunPairedEnd<
                                SimpleKmerHandler<ClosedSyncmer>,
                                AnchorFinder,
                                SimpleAlignmentHandler,
                                OutputHandler,
                                DEBUG_NONE,
                                AlignmentBenchmark>(
                                reader, options, anchor_finder, alignment_handler, output_handler, iterator, genomes, benchmark);
                        if (options.Verbose()) {
                            protal_stats.WriteStats();
                        }
                    }
                    bm_classify_sample.Stop();
                    bm_classify_sample.PrintResults();

                    is1.close();
                    is2.close();
                    sam_output.close();

                    if (!reader.Success()) {
                        RunStatus::Get().Fail("Reading the FASTQ files of sample " + options.GetSampleId(index) + " failed (" +
                                              options.GetFirstFile(index) + ", " + options.GetSecondFile(index) +
                                              "); no SAM file was written");
                        std::filesystem::remove(sam_partial);
                        continue;
                    }
                    if (sam_output.fail()) {
                        RunStatus::Get().Fail("Writing the SAM file of sample " + options.GetSampleId(index) + " failed: " + sam_partial);
                        std::filesystem::remove(sam_partial);
                        continue;
                    }
                    FinishSamFile(options, index, sam_partial);
                }

            } else {
                std::cout << "Single-end mode is not working currently. This will be fixed with the next version" << std::endl;
                exit(8);
                // AnchorFinder
                AnchorFinder anchor_finder(kmer_lookup, mmer_size, 4, options.GetMaxSeedSize(), genomes);
                // AlignmentHandler approach
                SimpleAlignmentHandler alignment_handler(genomes, aligner, kmer_size, options.GetAlignTop(), options.GetMaxScoreAni(), options.FastAlign());

                for (auto index : options.GetRange()) {
                    auto [sam, gzipped] = options.SamFile(index);
                    auto dir = std::filesystem::path(sam).parent_path();
                    
                    if (!std::filesystem::create_directories(dir.string()) && !std::filesystem::exists(dir)) {
                        std::cout << "Cannot create directories for this path " << sam << std::endl;
                        exit(32);
                    }

                    // Avoid aligning files that already exist.
                    if (!options.Force() && std::filesystem::exists(sam)) {
                        std::cout << "Skip " << sam << " continue" << std::endl;
                        continue;
                    }

                    std::ofstream sam_output(sam, std::ios::out);
                    genomes.WriteSamHeader(sam_output);
                    using OutputHandler = ProtalOutputHandler;
                    OutputHandler output_handler(sam_output, 1024*512, 1024*1024*16, genomes, 0.8);
                    igzstream is { options.GetFirstFile(index).c_str() };
                    SeqReader reader{ is };

                    auto protal_stats = protal::classify::Run<
                            SimpleKmerHandler<ClosedSyncmer>,
                            AnchorFinder,
                            SimpleAlignmentHandler,
                            OutputHandler,
                            DEBUG_NONE,
                            AlignmentBenchmark>(
                            reader, options, anchor_finder, alignment_handler, output_handler, iterator, benchmark);

                    if (options.Verbose()) {
                        protal_stats.WriteStats();
                    }
                    
                    // Close output streams
                    sam_output.close();
                    is.close();

                    // Handle gzip compression if needed
                    if (!options.IsSamFileGzipped(index)) {
                        try {
                            Compressor::compressInPlace(sam, options.GetThreads());
                            options.SetSamFileGzip(index, true);
                        } catch (const std::exception& e) {
                            std::cerr << "[WARNING] " << e.what() << std::endl;
                        }
                    }

                    // Clean up on read error
                    if (!reader.Success()) {
                        std::cerr << "There was an error reading the fastq file with sample " << options.GetSampleId(index) << " (" << index << ")" << std::endl;
                        std::cerr << options.GetFirstFile(index) << std::endl;
                        std::cerr << "Remove sam file: " << sam << std::endl;
                        std::filesystem::remove(sam);
                    }
                }
            }
            bm_classify.Stop();
            bm_classify.PrintResults();
        }
    }

    using Profile = profiler::MicrobialProfile;
    using Profiles = std::vector<Profile>;

//     void ProfileWrapper2(Options& options, ProtalDB& db) {
//         std::cout << "ProfileWrapper2" << std::endl;
//         GenomeLoader& genomes = db.GetGenomes();
//
//         if (!db.IsTaxonomyLoaded()) db.LoadTaxonomy(options.GetInternalTaxonomyFile());
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

    Profiles ProfileWrapper(Options& options, ProtalDB& db) {
        GenomeLoader& genomes = db.GetGenomes();

        if (!db.IsTaxonomyLoaded()) db.LoadTaxonomy(options.GetInternalTaxonomyFile());
        auto& taxonomy = db.GetTaxonomy();

        // These values do not matter anymore when a RandomForest is applied
        double min_ani = 0.95;
        double min_gene_presence = 0.50; //previously 0.5
        size_t min_total_hits = 60; //previously 70
        size_t min_mean_mapq = 10;


        using TaxonFilterObj = profiler::TaxonFilterObj;
//        TaxonFilterObj filter(min_ani, min_gene_presence, min_total_hits, min_mean_mapq);

        std::string model_path = options.GetModelPath();
        TaxonFilterObj filter(model_path, options.GetKnob());

        auto range = options.GetRange();

        // Pre-size so each thread writes to its own index slot — no emplace_back races.
        // MicrobialProfile holds a reference member so it is not assignable; use optional to allow
        // in-place construction per slot without requiring assignment.
        std::vector<std::optional<profiler::MicrobialProfile>> profile_slots(range.size());

        omp_set_num_threads(options.GetThreads());

        #pragma omp parallel for firstprivate(filter) shared(options, cout, taxonomy, profile_slots, genomes, std::cerr)//, bm_read_alignments, bm_profile)
        for (int idx = 0; idx < static_cast<int>(range.size()); idx++) {
            auto i = range[idx];

            if (options.Verbose()) {
                auto [sam, gzipped] = options.SamFile(i);
                #pragma omp critical(print)
                std::cerr << omp_get_thread_num() << " File " << i << " of " << range.size() << ":\n\t" << sam << (gzipped ? " (gzipped)" : "") << std::endl;
            }

            auto [sam, gzipped] = options.SamFile(i);
            auto sample_name = options.GetSampleId(i);

            if (!Utils::exists(sam)) {
                RunStatus::Get().Fail("No SAM file to profile for sample " + options.GetSampleId(i) + ": " + sam);
                profile_slots[idx].emplace(genomes);
                continue;
            }

            Benchmark bm_read_alignments{ "Load read alignments" };
            Benchmark bm_profile{ "Profile sample" };

            profiler::Profiler profiler(genomes);

            // New Profiler approach
            std::vector<AlignmentPair> unique_pairs;
            std::vector<std::vector<AlignmentPair>> pairs;


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

            std::string sam_error;
#pragma omp critical(load_sam)
            sam_error = profiler.FromSam(sam);

            // A SAM without alignments still gets its (empty) profile files, so that every sample
            // has output; only an unreadable SAM is a failure.
            if (!sam_error.empty()) {
                RunStatus::Get().Fail("Cannot read the SAM file of sample " + options.GetSampleId(i) + " (" + sam + "): " + sam_error);
                profile_slots[idx].emplace(genomes);
                continue;
            }


            std::ofstream erro(sam + ".err", std::ios::out);

            // profiler.PrintStats();
            bm_profile.Start();


            if (options.Verbose()) {
                #pragma omp critical(print)
                std::cout << "Thread " << omp_get_thread_num() << " run profile" << std::endl;
            }

            auto profile = profiler.Profile(sample_name, std::optional<std::reference_wrapper<std::ostream>>{erro},
                                            options.GetSNPMinCov(), options.GetSNPMinCov(),
                                            options.GetSNPMinAF(), options.GetSNPMinMeanQual(),
                                            options.GetSNPMinPhredSum(), options.GetSNPRequireStrand());
            erro.close();
            bm_profile.Stop();
            

            if (options.Verbose()) {
                #pragma omp critical(print)
                {
                    std::cout << "Thread " << omp_get_thread_num() << " ";
                    bm_profile.PrintResults();
                }
            }


            // profile.bm_add_sam.PrintResults();


            std::optional<TruthSet> truth = options.HasProfileTruths() ?
                                            std::optional<TruthSet>{ protal::GetTruth(options.ProfileTruthFile(i), taxonomy) } :
                                            std::optional<TruthSet>{};

            if (options.BenchmarkAlignment()) {
                profiler.TestSNPUtils(pairs);
                profiler.OutputErrorData(pairs);
            } else if (truth.has_value()) {
                profiler.OutputErrorData(unique_pairs, pairs, &truth.value());
            }


            if (truth.has_value()) {
                std::string truth_output = options.ProfileFile(i) + ".truth_annotated";
                profile.AnnotateWithTruth(truth.value(), filter, truth_output, taxonomy);
                std::cout << "Write truth to: " << truth_output << std::endl;

                auto filtered = profile.GetTaxa() | std::views::filter([&filter](auto a) { return filter.Pass(a.second); });
            
                // std::filter(profile.GetTaxa().begin(), profile.GetTaxa().end(), )

                auto tp = std::count_if(
                    filtered.begin(),
                    filtered.end(),
                    [&truth] (auto x) { return truth.value().contains(x.first); }
                );
                auto fp = std::ranges::distance(filtered) - tp;
                auto fn = truth.value().size() - tp;

                std::cout << "TP: " << tp << " FP: " << fp << " FN: " << fn << std::endl;

            }
            if (options.Verbose()) {
                std::cout << "Write profile to: \n" << options.ProfileFile(i) << std::endl;
            }
            {
                auto dir = std::filesystem::path(options.ProfileFile(i)).parent_path();
#pragma omp critical(create_dir)
                if (!std::filesystem::exists(dir)) {
                    if (!std::filesystem::create_directories(dir.string())) {
                        std::cerr << "Cannot create directories for path " << options.ProfileFile(i) << std::endl;
                        exit(2);
                    }
                }
            }

            std::ofstream os(options.ProfileFile(i), std::ios::out);
            std::ofstream os_total(options.ProfileFile(i) + ".log", std::ios::out);
            std::ofstream os_dismissed(options.ProfileFile(i) + ".gene.log", std::ios::out);
            std::ofstream os_genes(options.ProfileFile(i) + ".genes.log", std::ios::out);

            profile.WriteSparseProfile(taxonomy, filter, os, &os_total, &os_dismissed);
            profile.WriteGeneProfile(taxonomy, filter, &os_genes);
            os.close();
            os_total.close();
            os_dismissed.close();
            os_genes.close();
            if (os.fail() || os_total.fail() || os_dismissed.fail() || os_genes.fail()) {
                RunStatus::Get().Fail("Writing the profile of sample " + options.GetSampleId(i) + " failed: " + options.ProfileFile(i));
            }

            profile.SetName(options.GetSampleId(i));

            // Each thread writes to its own pre-allocated slot — no lock needed.
            profile_slots[idx].emplace(std::move(profile));
        }

        std::vector<profiler::MicrobialProfile> profiles;
        profiles.reserve(profile_slots.size());
        for (auto& slot : profile_slots) {
            if (slot.has_value()) profiles.emplace_back(std::move(slot.value()));
        }
        return profiles;
    }

    static std::pair<bool,bool> SharedSNP(VariantBin& a, VariantBin& b) {
        std::sort(a.begin(), a.end(), [](Variant const& va, Variant const& vb) {
            return va.Observations() > vb.Observations();
        });
        std::sort(b.begin(), b.end(), [](Variant const& va, Variant const& vb) {
            return va.Observations() > vb.Observations();
        });

        bool shared_non_ref = false;
        bool shared_ref = false;
        for (auto& var_a : a) {
            if (!var_a.GetValid()) continue;
            for (auto& var_b : b) {
                if (!var_b.GetValid()) continue;

                if (var_a.Match(var_b))  {
                    shared_ref |= var_a.IsReference() && var_b.IsReference();
                    shared_non_ref |= !var_a.IsReference() && !var_b.IsReference();
                }
            }
        }

        return { shared_ref, shared_non_ref };
    }

    static double JukesCantor(size_t mismatch, size_t length) {
        auto strain_identity = static_cast<double>(mismatch) / length;

        auto jc69 = -3.0f/4 * log(1 - (4.0f/3 * strain_identity));
        return jc69 > 1 || jc69 != jc69 ? 1 : jc69;
    }

    static bool IsSNP(VariantBin& bin, size_t min_obs, double min_frequency) {
        for (auto& var : bin) {
            size_t total_obs = std::accumulate(bin.begin(), bin.end(), 0, [](size_t acc, Variant const& var) { return acc + var.Observations(); });
            auto obs = var.Observations();
            double freq = static_cast<double>(obs) / total_obs;
            if (!var.IsReference() && freq >= min_frequency && (obs >= min_obs || var.HasFwdAndRev())) {
                return true;
            }
        }
//
//        for (auto& var : bin) {
//            std::cout << var.ToString() << '\t';
//        }
//        std::cout << endl;
//
//        Utils::Input();
        return false;
    }

    static bool HasReference(VariantBin& bin) {
        return std::any_of(bin.begin(), bin.end(), [](Variant const& var) {
            return var.IsReference();
        });
    }

    static std::pair<size_t, size_t> Distance(size_t shared_length, VariantVec& a, VariantVec& b) {
        const auto pos_getter = [](VariantBin const& bin) { return bin.front().Position(); };
        size_t shared_snps = 0;
        size_t unique_snps_a = 0;
        size_t unique_snps_b = 0;

        size_t min_observations = 2;
        double min_frequency = 0.2;

        int ai = 0, bi = 0;
//        std::cout << "a: " << a.size() << "  b: " << b.size() << std::endl;
        while (ai < a.size() && bi < b.size()) {
            auto pos_a = pos_getter(a[ai]);
            auto pos_b = pos_getter(b[bi]);

            bool has_ref_a = std::any_of(a[ai].begin(), a[ai].end(), [](Variant const& v) {return v.GetValid() && v.IsReference();});
            bool has_ref_b = std::any_of(b[bi].begin(), b[bi].end(), [](Variant const& v) {return v.GetValid() && v.IsReference();});
            bool has_snp_a = std::any_of(a[ai].begin(), a[ai].end(), [](Variant const& v) {return v.GetValid() && !v.IsReference();});
            bool has_snp_b = std::any_of(b[bi].begin(), b[bi].end(), [](Variant const& v) {return v.GetValid() && !v.IsReference();});
//            bool has_snp_b = IsSNP(b[bi], min_observations, min_frequency);

            bool same_pos = (pos_a == pos_b);
            bool advance_a = (pos_a <= pos_b);
            bool advance_b = (pos_b <= pos_a);

//            std::cout << ai << " " << pos_a << "\t->\t" << a[ai].front().ToString() << "\t" << b[bi].front().ToString() << "\t<-\t" << pos_b << " " << bi << std::endl;

            // same_pos:                 required because a shared snp needs variants detected at the same pos.
            // has_snp_a || has_snp_b:   if both have variants at the position but none are considered
            //                           genuine snps, it is assumed that both are reference alleles
            // SharedSNP(a[ai], b[bi]):  Check if they share at least one valid snp, ref or non ref.
            auto [has_shared_ref, has_shared_snp] = same_pos ? SharedSNP(a[ai], b[bi]) : std::pair<bool,bool>(false, false);
            bool is_shared_and_snp = has_shared_snp && !has_shared_ref;
            shared_snps += is_shared_and_snp;
            // has_snp_a: if all variants at this position are dodgy/fp SNPs assume its the reference allele
            // same_po
            // advance var is needed because only the smaller pos can be a snp
            bool is_unique_snp_a = (advance_a && !has_shared_ref && has_snp_a && !has_ref_a && !has_shared_snp);
            bool is_unique_snp_b = (advance_b && !has_shared_ref && has_snp_b && !has_ref_b && !has_shared_snp);
            unique_snps_a += is_unique_snp_a;
            unique_snps_b += is_unique_snp_b;

//            if (is_unique_snp_a || unique_snps_b) {
//                std::cout << std::string(60, '#') << std::endl;
//                std::cout << (is_unique_snp_a ? "SNP IN A" : "SNP IN B") << std::endl;
//
//                std::cout << "____________SNPS A noref:" << !has_ref_a << std::endl;
//                for (auto& var : a[ai]) {
//                    std::cout << var.ToString() << '\t';
//                }
//                std::cout << std::endl;
//                std::cout << "____________SNPS B noref:" << !has_ref_b << std::endl;
//                for (auto& var : b[bi]) {
//                    std::cout << var.ToString() << '\t';
//                }
//                std::cout << std::endl;
//                Utils::Input();
//            }

//            if ((has_snp_a && advance_a && !shared) || (has_snp_b && advance_b && !shared)) {
//                std::cout << "PosA: " << pos_a << "\tPosB:" << pos_b << std::endl;
//                std::cout << "SNPA: " << has_snp_a << "\tSNPB:" << has_snp_b << std::endl;
//                std::cout << "shared: " << shared << std::endl;
//                std::cout << "-------\nA is SNP? " << has_snp_a << std::endl;
//                for (auto& var : a[ai]) {
//                    std::cout << var.ToString() << '\t';
//                }
//                std::cout << endl;
//                std::cout << "B is SNP? " << has_snp_b << std::endl;
//                for (auto& var : b[bi]) {
//                    std::cout << var.ToString() << '\t';
//                }
//                std::cout << endl;
//
//                Utils::Input();
//            }

            ai += advance_a;
            bi += advance_b;
        }

        if (a.size() - ai > 0) {
//            std::cout << "A) SNPS_A: " << unique_snps_a << "\tSNPS_B: " << unique_snps_b << std::endl;
            for (; ai < a.size(); ai++) {
                bool has_var = std::any_of(a[ai].begin(), a[ai].end(), [](Variant const& v) {return v.GetValid() && !v.IsReference();});
                bool has_ref = std::any_of(a[ai].begin(), a[ai].end(), [](Variant const& v) {return v.GetValid() && v.IsReference();});
                unique_snps_a += (has_var && !has_ref);
            }
//            std::cout << "  -->  SNPS_A: " << unique_snps_a << "\tSNPS_B: " << unique_snps_b << std::endl;
        }

        if (b.size() - bi > 0) {
//            std::cout << "B) SNPS_A: " << unique_snps_a << "\tSNPS_B: " << unique_snps_b << std::endl;
            for (; bi < b.size(); bi++) {
                bool has_var = std::any_of(b[bi].begin(), b[bi].end(), [](Variant const& v) {return v.GetValid() && !v.IsReference();});
                bool has_ref = std::any_of(b[bi].begin(), b[bi].end(), [](Variant const& v) {return v.GetValid() && v.IsReference();});
                unique_snps_b += (has_var && !has_ref);
            }
        }

//        if (unique_snps_a + unique_snps_b > 0) {
//            Utils::Input();
//        }

        return std::pair{ unique_snps_a + unique_snps_b, shared_length };

//        return JukesCantor(unique_snps_a + unique_snps_b, shared_length);
    }

    static std::pair<double, size_t> Distance(profiler::Gene& g1, profiler::Gene& g2) {
        auto& strain1 = g1.GetStrainLevel();
        auto& strain2 = g2.GetStrainLevel();
        auto shared_alignment_region = SharedAlignmentRegion::GetSharedAlignmentRegion(strain1, strain2);

        auto& va = shared_alignment_region.variant_handler_a;
        auto& vb = shared_alignment_region.variant_handler_b;
        auto shared_sequence_length = shared_alignment_region.share_range.SequenceLength();


        auto [unique_snps, shared_length] = Distance(shared_sequence_length, va, vb);

        return {unique_snps, shared_sequence_length};
    }

    static std::pair<double, size_t> Distance(profiler::Taxon& t1, profiler::Taxon& t2) {

        double distance_sum = 0;
        size_t total_unique_snps = 0;
        size_t shared_region_sum = 0;

        for (auto& [gid, _] : t1.GetGenes()) {
            if (t2.GetGenes().contains(gid)) {
                //std::cout << "Stuck at " << gid << std::endl;
                auto [unique_snps, shared_region] = Distance(t1.GetGenes().at(gid), t2.GetGenes().at(gid));

                total_unique_snps += unique_snps;
                shared_region_sum += shared_region;
            }
        }

        auto distance = JukesCantor(total_unique_snps, shared_region_sum);

//        if (distance < 0.001) {
//            std::cout << "Distance: " << distance << " with " << total_unique_snps << " snps over " << shared_region_sum << " bases." << std::endl;
//        }

        return { distance, shared_region_sum };
    }


    using StrainResults = tsl::robin_map<size_t, Matrix<double>>;
    using TaxidCounts = tsl::robin_map<uint32_t, uint32_t>;
    using TaxidSet = tsl::robin_set<uint32_t>;
    using TaxidList = std::vector<uint32_t>;
    using OptionalFilter = optional<std::reference_wrapper<const profiler::TaxonFilter>>;
    using SimilarityMatrix = DoubleMatrix;

    TaxidList ExtractTaxa(Profiles const& profiles, std::optional<profiler::TaxonFilterObj> filter= {}) {
        TaxidCounts taxid_counts;
        TaxidSet taxa;
        TaxidList taxid_list;

        for (auto& profile : profiles) {
            for (auto& [id, taxon] : profile.GetTaxa()) {
                if (filter.has_value() && !filter->Pass(taxon)) continue;
                if (!taxid_counts.contains(id)) {
                    taxid_counts.insert({ id, 0 });
                }
                taxid_counts.at(id)++;
            }
        }

        for (auto& [t, c] : taxid_counts) {
            if (c > 1) {
                taxid_list.emplace_back(t);
            }
        }

        std::sort(taxid_list.begin(), taxid_list.end(), [&taxid_counts](uint32_t const& t1, uint32_t const& t2) {
            return taxid_counts.at(t1) > taxid_counts.at(t2);
        });

        return taxid_list;
    }

    static std::vector<uint32_t> ResolveMSASpecies(Options& options, taxonomy::IntTaxonomy& taxonomy) {
        const auto& msa_species = options.GetMSASpecies();
        if (msa_species.empty()) return {};

        std::vector<uint32_t> taxids;
        taxids.reserve(msa_species.size());
        std::vector<std::string> invalid;
        invalid.reserve(msa_species.size());

        for (const auto& raw_spec : msa_species) {
            std::string spec = raw_spec;
            std::replace(spec.begin(), spec.end(), ' ', '_');
            if (spec.rfind("s__", 0) != 0 || spec.find('_', 3) == std::string::npos) {
                invalid.emplace_back(raw_spec);
                continue;
            }

            uint32_t taxid = 0;
            bool found = false;
            if (taxonomy.string_to_id.contains(spec)) {
                taxid = taxonomy.Get(spec);
                found = true;
            } else {
                std::string alt = spec;
                for (size_t i = 3; i < alt.size(); i++) {
                    if (alt[i] == '_') alt[i] = ' ';
                }
                if (taxonomy.string_to_id.contains(alt)) {
                    taxid = taxonomy.Get(alt);
                    found = true;
                }
            }

            if (!found) {
                invalid.emplace_back(raw_spec);
                continue;
            }
            if (taxonomy.Get(taxid).rank != "species") {
                invalid.emplace_back(raw_spec);
                continue;
            }
            taxids.emplace_back(taxid);
        }

        if (!invalid.empty()) {
            std::cerr << "Invalid --msa_species entries (expected format s__Genus_species and present in taxonomy): ";
            std::cerr << Utils::join(invalid, ",") << std::endl;
            exit(2);
        }

        return taxids;
    }

    static std::pair<double, size_t> GetSimilarity(uint32_t taxid, Profile& profile1, Profile& profile2, Options& options, std::optional<profiler::TaxonFilterObj> filter={}) {
        size_t min_shared_region = 1000;
        if (profile1.GetTaxa().contains(taxid) && profile2.GetTaxa().contains(taxid)) {
            auto& taxon1 = profile1.GetTaxa().at(taxid);
            auto& taxon2 = profile2.GetTaxa().at(taxid);

            if (filter.has_value() && !(filter->Pass(taxon1) && filter->Pass(taxon2))) {
                return { NAN, 0 };
            }

            auto [distance, shared_region] = Distance(taxon1, taxon2);

            auto similarity = 1.0f-distance;

            if (shared_region <= min_shared_region) {
                similarity = NAN;
            }
            return { similarity, shared_region };
        }
        return { NAN, 0 };
    }

    static bool IsRowGood(std::vector<char> row, size_t min_hcov) {
        size_t count_non_n = 0;
        for (auto c : row) {
            count_non_n += (c != 'N' && c != '-');
        }
        return count_non_n >= min_hcov;
    }

    static std::vector<uint32_t> GetInformationVector(MSAVector const& msa_vector) {
        size_t info = 0;
        const auto msa_len = msa_vector.front().size();
        std::vector<uint32_t> information_vector(msa_len, 0);

        for (auto pos = 0; pos < msa_len; pos++) {
            info = 0;
            for (auto const& seq : msa_vector) info += (seq[pos] != 'N' && seq[pos] != '-');
            information_vector[pos] = info;
        }
        return information_vector;
    }

    static MSAVector ProcessMSA(MSAVector const& msa_vector, double position_coverage = 0.5) {
        auto sample_len = msa_vector.size();
        auto info_vector = GetInformationVector(msa_vector);
        MSAVector processed_msa(sample_len, std::vector<char>{});

        auto sample_size_threshold = msa_vector.size() * position_coverage;

        // std::cout << "Sample size threshold: " << sample_size_threshold << std::endl;
        
        for (auto pos = 0; pos < info_vector.size(); pos++) {
            if (info_vector[pos] > sample_size_threshold) {
                for (auto s = 0; s < sample_len; s++) {
                    processed_msa[s].emplace_back(msa_vector[s][pos]);
                }
            }
        }
        return processed_msa;
    }

    static void OutputMSA(MSAVector const& msa, std::vector<std::string> const& names, std::ostream& os=std::cout) {
        for (auto i = 0; i < msa.size(); i++) {
            auto& row = msa[i];

            if (row.empty()) continue;
            os << ">" << names[i] << std::endl;
            // os << std::string_view(&row[0], std::distance(row.begin(), row.end())) << std::endl;
            os << std::string(&row[0], std::distance(row.begin(), row.end())) << std::endl;
//            os << std::string_view(row.begin(), row.end()) << std::endl;
        }
    }

    static std::vector<size_t> GetProfilesWithTaxon(uint32_t taxid, Profiles& profiles, Options& options, std::optional<profiler::TaxonFilter>& filter) {
        std::vector<size_t> indices;

        for (auto i = 0; i < profiles.size(); i++) {
            auto& profile = profiles[i];
            if (!profile.GetTaxa().contains(taxid)) {
                continue;
            }
            auto& taxon_map = profile.GetTaxa();

            auto& taxon = taxon_map.at(taxid);
            if (filter.has_value() && !filter.value().Pass(taxon)) {
                continue;
            }
            // std::cout << taxid << " Passes " << std::endl;

            indices.emplace_back(i);
        }
        return indices;
    }

    static std::vector<uint32_t> SelectGenesForTaxon(uint32_t taxid, std::string name, std::vector<size_t>& selected_profiles, GenomeLoader& loader, Options& options, Profiles& profiles) {
        std::vector<uint32_t> selected_gene_ids;
        std::vector<uint32_t> gene_ids = loader.GetGenome(taxid).GetHittableGenes();

        auto max_gene_id = std::max_element(gene_ids.begin(), gene_ids.end());

        // std::vector<std::vector<double>> per_sample_gene_multiallelic_portion;
        // per_sample_gene_multiallelic_portion.resize(selected_profiles.size(),
        //     std::vector<double>(*max_gene_id + 1, -1) );

        std::vector<std::vector<int>> per_sample_gene_multiallelic_snps;
        per_sample_gene_multiallelic_snps.resize(selected_profiles.size(),
            std::vector<int>(*max_gene_id + 1, -2) );

        std::vector<std::vector<int>> per_sample_gene_snps;
        per_sample_gene_snps.resize(selected_profiles.size(),
            std::vector<int>(*max_gene_id + 1, -2) );

        std::vector<std::vector<int>> per_sample_gene_cov;
        per_sample_gene_cov.resize(selected_profiles.size(),
            std::vector<int>(*max_gene_id + 1, -2) );

        std::vector<std::vector<int>> per_sample_gene_noise;
        per_sample_gene_noise.resize(selected_profiles.size(),
            std::vector<int>(*max_gene_id + 1, -2) );

        for (auto si = 0; si < selected_profiles.size(); si++) {
            auto sample_index = selected_profiles[si];
            auto& taxon = profiles[sample_index].GetTaxa().at(taxid);
            auto& multiallelic_vec = per_sample_gene_multiallelic_snps[si];

            // std::cout << taxon.GetName() << " -> hittable genes" << gene_ids.size() << std::endl;
            for (auto& gene_id : gene_ids) {
                if (gene_id >= per_sample_gene_multiallelic_snps[si].size()) {
                    std::cout << gene_id << " >= " << gene_ids.size() << std::endl;
                    exit(123);
                }

                per_sample_gene_multiallelic_snps[si][gene_id] = -1;
                per_sample_gene_snps[si][gene_id] = -1;
                per_sample_gene_cov[si][gene_id] = -1;
                per_sample_gene_noise[si][gene_id] = -1;
                
                if (!taxon.GetGenes().contains(gene_id)) {
                    per_sample_gene_multiallelic_snps[si][gene_id] = 0;
                    per_sample_gene_snps[si][gene_id] = 0;
                    per_sample_gene_cov[si][gene_id] = 0;
                    per_sample_gene_noise[si][gene_id] = 0;
                    continue;
                }

                auto& gene = taxon.GetGene(gene_id);

                auto allele_counts = gene.AlleleSNPCounts(2, 60);
                auto allele_counts_nf = gene.AlleleSNPCounts(0, 0);

                // std::cout << "----Filter----\n" << allele_counts.ToString() << std::endl;
                // std::cout << "----No Filter----\n" << allele_counts_nf.ToString() << std::endl;
                
                per_sample_gene_multiallelic_snps[si][gene_id] = allele_counts.Multi();
                per_sample_gene_snps[si][gene_id] = allele_counts.AllValid();
                per_sample_gene_cov[si][gene_id] = gene.Coverage();
                per_sample_gene_noise[si][gene_id] = allele_counts.Noisy();
            }
        }

        std::ofstream os(options.GetMiscOutputDir() + '/' + name + ".snps_multiallelic.tsv", std::ios::out);
        for (auto si = 0; si < selected_profiles.size(); si++) {
            auto sample_index = selected_profiles[si];
            os << options.GetSampleId(sample_index);

            if (si >= per_sample_gene_multiallelic_snps.size()) exit(244);
            auto& vec = per_sample_gene_multiallelic_snps[si];

            double total = std::accumulate(vec.begin(), vec.end(), 0.0, [](double acc, int val) {
                return acc + (val < 0 ? 0 : val);
            });

            os << '\t' << total;
            for (auto i = 1; i < vec.size(); i++) {
                os << '\t' << vec[i];
            }
            os << std::endl;
        }
        os.close();

        std::ofstream os2(options.GetMiscOutputDir() + '/' + name + ".snps_total.tsv", std::ios::out);
        for (auto si = 0; si < selected_profiles.size(); si++) {
            auto sample_index = selected_profiles[si];
            os2 << options.GetSampleId(sample_index);

            if (si >= per_sample_gene_snps.size()) exit(244);
            auto& vec = per_sample_gene_snps[si];

            double total = std::accumulate(vec.begin(), vec.end(), 0.0, [](double acc, int val) {
                return acc + (val < 0 ? 0 : val);
            });

            os2 << '\t' << total;
            for (auto i = 1; i < vec.size(); i++) {
                os2 << '\t' << vec[i];
            }
            os2 << std::endl;
        }
        os2.close();

        std::ofstream os3(options.GetMiscOutputDir() + '/' + name + ".hcov.tsv", std::ios::out);
        for (auto si = 0; si < selected_profiles.size(); si++) {
            auto sample_index = selected_profiles[si];
            os3 << options.GetSampleId(sample_index);

            if (si >= per_sample_gene_snps.size()) exit(244);
            auto& vec = per_sample_gene_cov[si];
            
            double total = std::accumulate(vec.begin(), vec.end(), 0.0, [](double acc, int val) {
                return acc + (val < 0 ? 0 : val);
            });

            os3 << '\t' << total;
            for (auto i = 1; i < vec.size(); i++) {
                os3 << '\t' << vec[i];
            }
            os3 << std::endl;
        }
        os3.close();


        std::ofstream os4(options.GetMiscOutputDir() + '/' + name + ".snps_filtered.tsv", std::ios::out);
        for (auto si = 0; si < selected_profiles.size(); si++) {
            auto sample_index = selected_profiles[si];
            os4 << options.GetSampleId(sample_index);

            if (si >= per_sample_gene_noise.size()) exit(244);
            auto& vec = per_sample_gene_noise[si];

            double total = std::accumulate(vec.begin(), vec.end(), 0.0, [](double acc, int val) {
                return acc + (val < 0 ? 0 : val);
            });

            os4 << '\t' << total;
            for (auto i = 1; i < vec.size(); i++) {
                os4 << '\t' << vec[i];
            }
            os4 << std::endl;
        }

        os4.close();

        // Long unique k-mers are unique k-mers in index blocks with flex keys, i.e. whose 15-mer core
        // is shared with a relative. A species with relatives in the database has them in nearly
        // every gene; a gene without any is one its relatives share unchanged or lack from their
        // reference. A relative's reads of such a gene align here, and its MSA columns would show
        // them as a second strain, so it is left out. A species without relatives has hardly any
        // long unique k-mers at all (its k-mers are unique at the core already): all its genes stay.
        auto& genome = loader.GetGenome(taxid);
        auto const with_long_uniques = std::count_if(gene_ids.begin(), gene_ids.end(), [&genome](uint32_t gene_id) {
            return genome.GetGene(gene_id).HasLongUniques();
        });
        if (static_cast<size_t>(with_long_uniques) * 2 < gene_ids.size()) return gene_ids;

        std::vector<uint32_t> msa_gene_ids;
        std::copy_if(gene_ids.begin(), gene_ids.end(), std::back_inserter(msa_gene_ids), [&genome](uint32_t gene_id) {
            return genome.GetGene(gene_id).HasLongUniques();
        });
        if (msa_gene_ids.size() < gene_ids.size()) {
            std::cout << name << ": " << gene_ids.size() - msa_gene_ids.size() << " of " << gene_ids.size()
                      << " genes have no long unique k-mers and are left out of the MSA" << std::endl;
        }
        return msa_gene_ids;
    }

    // Shell-quote a path/argument for use in a std::system command line.
    static std::string ShellQuote(const std::string& s) {
        std::string out = "'";
        for (char c : s) {
            if (c == '\'') out += "'\\''";
            else out += c;
        }
        out += "'";
        return out;
    }

    static bool IsExecutableFile(std::filesystem::path const& p) {
        std::error_code ec;
        return std::filesystem::is_regular_file(p, ec) && ::access(p.c_str(), X_OK) == 0;
    }

    // Look up an executable by name in $PATH, the way a shell would.
    static std::string FindInPath(std::string const& name) {
        const char* path_env = std::getenv("PATH");
        if (!path_env || !*path_env) return "";
        std::string const path(path_env);
        for (size_t start = 0; start <= path.size(); ) {
            size_t end = path.find(':', start);
            if (end == std::string::npos) end = path.size();
            if (end > start) {
                std::filesystem::path cand = std::filesystem::path(path.substr(start, end - start)) / name;
                if (IsExecutableFile(cand)) return cand.string();
            }
            start = end + 1;
        }
        return "";
    }

    // Locate the qcmsa post-filter. qcmsa ships as an executable named 'qcmsa'
    // installed next to the protal binary (see `just install`), so that is what
    // we look for: explicit --qcmsa_script, then PROTAL_QCMSA_SCRIPT, then
    // 'qcmsa' next to the protal executable, then 'qcmsa' on $PATH. The
    // qcmsa.py script of a source checkout is only a last resort so that running
    // protal straight out of the build tree keeps working.
    // Returns "" if not found.
    static std::string FindQCMSAScript(Options& options) {
        namespace fs = std::filesystem;

        if (!options.GetQCMSAScript().empty()) {
            std::string const& given = options.GetQCMSAScript();
            // A bare command name is resolved via $PATH, a path is taken as is.
            if (given.find('/') == std::string::npos) {
                if (std::string in_path = FindInPath(given); !in_path.empty()) return in_path;
            }
            return given;
        }
        if (const char* env = std::getenv("PROTAL_QCMSA_SCRIPT"); env && *env)
            return std::string(env);

        std::error_code ec;
        fs::path exe = fs::canonical("/proc/self/exe", ec);
        if (!ec && IsExecutableFile(exe.parent_path() / "qcmsa"))
            return (exe.parent_path() / "qcmsa").string();

        if (std::string in_path = FindInPath("qcmsa"); !in_path.empty()) return in_path;

        if (!ec) {
            fs::path dir = exe.parent_path();
            for (const auto& cand : { dir / "qcmsa.py",                       // alongside the binary
                                      dir / "scripts" / "qcmsa.py",
                                      dir.parent_path() / "scripts" / "qcmsa.py" }) {
                if (fs::exists(cand)) return cand.string();
            }
        }
        return "";
    }

    // M5 step 4c: invoke the qcmsa post-filter on a species' MSA. A failure does not stop the run --
    // the raw strain outputs protal already wrote remain valid -- but it is reported and makes protal
    // exit non-zero, because the filtered <name>.msa.fna the run was asked for is missing.
    static void RunQCMSA(Options& options, const std::string& name) {
        namespace fs = std::filesystem;
        std::string script = FindQCMSAScript(options);
        if (script.empty() || !fs::exists(script)) {
            static std::atomic<bool> reported{false};
            if (!reported.exchange(true)) {
                RunStatus::Get().Fail("qcmsa not found next to the protal binary or on $PATH, so no species got "
                                      "its filtered .msa.fna (install it with 'just install', point protal at it "
                                      "with --qcmsa_script / PROTAL_QCMSA_SCRIPT, or pass --no_qcmsa)");
            }
            return;
        }
        // qcmsa is an executable with a python3 shebang. Only fall back to
        // running it through the interpreter if it is not executable (e.g. a
        // qcmsa.py picked straight out of a source checkout).
        std::string launcher = IsExecutableFile(script) ? "" : "python3 ";
        // Run qcmsa on protal's native (raw) MSA; it writes the final <name>.msa.fna.
        std::string msa = options.GetMSAOutput(name);              // .raw.msa.fna
        std::string partition = options.GetMSAPartitionOutput(name); // .raw.partition.txt
        std::string meta = options.GetSpeciesMetaOutput(name);
        if (!fs::exists(msa) || !fs::exists(partition) || !fs::exists(meta)) {
            std::cerr << "[qcmsa] WARNING: missing MSA/partition/meta for " << name
                      << "; skipping post-filter" << std::endl;
            return;
        }
        std::string prefix = options.GetStrainOutputDir() + '/' + name;
        // Remove qcmsa outputs of an earlier run first: if qcmsa writes nothing this time (e.g. no
        // sample passes its filters), a stale <name>.msa.fna must not be left behind as if it were new.
        for (auto const* ext : { ".msa.fna", ".partition.txt", ".qcmsa_summary.tsv", ".qc.png" }) {
            std::error_code ec;
            fs::remove(prefix + ext, ec);
        }
        std::ostringstream cmd;
        cmd << launcher << ShellQuote(script)
            << ' ' << ShellQuote(msa)
            << ' ' << ShellQuote(partition)
            << ' ' << ShellQuote(meta)
            << " --prefix " << ShellQuote(prefix)
            << " --reapply-hcov " << options.GetMSAMinHCOV();
        // User-supplied flags go last so they win over the defaults protal passes
        // above. Forwarded verbatim (unquoted) -- they are a flag list, not a value.
        if (!options.GetQCMSAArgs().empty()) cmd << ' ' << options.GetQCMSAArgs();
        std::cerr << "[qcmsa] " << cmd.str() << std::endl;
        int rc = std::system(cmd.str().c_str());
        if (rc != 0) {
            RunStatus::Get().Fail("qcmsa exited with code " + std::to_string(rc) + " for " + name +
                                  " (post-filter skipped; the raw MSA is still in " + msa + ")");
        } else if (!fs::exists(prefix + ".msa.fna")) {
            // Not an error (qcmsa may filter everything out), but no filtered MSA is not a success either.
            std::cerr << "[qcmsa] WARNING: " << name << ": qcmsa kept no gene or sample, so there is no "
                      << prefix << ".msa.fna (see its output above; the raw MSA is " << msa << ")" << std::endl;
        }
    }

    static void GetMSAForTaxon (uint32_t taxid, std::string taxon_name, GenomeLoader& loader, Options& options, Profiles& profiles, std::ostream* os_meta=nullptr, std::optional<profiler::TaxonFilter> filter={}) {
        // An earlier run's MSA must not survive a run that writes none (qcmsa would filter it).
        for (auto const& stale : { options.GetMSAOutput(taxon_name), options.GetMSAPartitionOutput(taxon_name) }) {
            std::error_code ec;
            std::filesystem::remove(stale, ec);
        }
        auto min_hcov = options.GetMSAMinHCOV();
        auto min_qual_sum = options.GetSNPMinPhredSum();
        auto min_cov = options.GetSNPMinCov();
        auto min_af = options.GetSNPMinAF();
        auto require_strand = options.GetSNPRequireStrand();
        auto min_mean_qual = options.GetSNPMinMeanQual();
        auto snp_max_alleles = options.GetSNPMaxAlleles();

        std::vector<size_t> profile_indices = GetProfilesWithTaxon(taxid, profiles, options, filter);

        if (profile_indices.empty()) return;

        MSAVector msa{ profile_indices.size(), std::vector<char>() };
        protal::MSARow ref_msa_row;

        // Per-sample accumulated SNP-retention statistics across all genes
        protal::MSAStats sample_stats(profile_indices.size());

        auto& genome = loader.GetGenome(taxid);
        if (!genome.IsLoaded()) genome.LoadGenomeOMP();
        std::vector<std::string> names;
        for (auto index : profile_indices) names.emplace_back(profiles[index].GetName());
        std::vector<std::string> partitions;
        size_t partition_start = 0;
        size_t previous_size = 0;

        std::vector<uint32_t> selected_genes = SelectGenesForTaxon(taxid, taxon_name, profile_indices, loader, options, profiles);

        ProgressBar prog(selected_genes.size());

        std::cout << taxid << ": " << taxon_name << " across samples " << profile_indices.size() << std::endl;
        

//        std::cout << "Process " << selected_genes.size() << std::endl;
        for (auto& geneid : selected_genes) {
//            std::cout << "GID: " << geneid << std::endl;
            // if (!loader.GetGenome(taxid).IsGeneHittable(geneid)) {
            //     continue;
            // }

            prog.UpdateAdd(1);
            MSASequenceItems items;

            // if (!genome.ValidGene(geneid)) break;
            //
            auto& gene = genome.GetGene(geneid);
            // if (!gene.IsSet()) continue;

            // Check if
            size_t samples_with_gene = 0;

            for (auto i = 0; i < profile_indices.size(); i++) {
                auto& profile = profiles[profile_indices[i]];

                auto& taxon_map = profile.GetTaxa();
                if (!taxon_map.contains(taxid)) continue;

                auto& genes = profile.GetTaxa().at(taxid).GetGenes();

                if (!genes.contains(geneid)) {
                    items.emplace_back(OptionalMSASequenceItem{});
                } else {
                    auto& gene_obs = genes.at(geneid);
                    auto& strain = gene_obs.GetStrainLevel();
                    auto& region = strain.GetSequenceRangeHandler();

                    auto ac = gene_obs.AlleleSNPCounts(min_cov, min_qual_sum);
                    region.CalculateCoverageVector();
                    auto tmp_vec = region.CalculateCoverageVector2();
                    auto counts_vcov1 = std::count_if(tmp_vec.begin(), tmp_vec.end(), [](auto val){ return(val >= 1);});
                    auto counts_vcov2 = std::count_if(tmp_vec.begin(), tmp_vec.end(), [](auto val){ return(val >= 2);});
                    // Multi-allelic positions as the MSA writes them (IUPAC codes), which qcmsa filters on.
                    size_t const multi_allelic = MultiAllelicPositions(strain.GetVariantHandler().GetVariants(), tmp_vec, min_cov,
                                                                       min_qual_sum, min_af, require_strand, min_mean_qual, snp_max_alleles);

                    double median_vcov = 0.0;
                    double mean_vcov_nonzero = 0.0;
                    double median_vcov_nonzero = 0.0;
                    double hcov = static_cast<double>(counts_vcov1) / static_cast<double>(gene_obs.m_gene_length);

                    // protal emits a raw MSA: every observed gene contributes its
                    // sequence. Coverage-based gating (horizontal coverage, depth,
                    // min samples per gene) is done by the qcmsa post-filter, which
                    // reads the hcov / mean_vcov_nonzero columns written below.
                    samples_with_gene++;
                    auto snps = SharedAlignmentRegion::GetSNPs(strain.GetVariantHandler());
                    items.emplace_back( OptionalMSASequenceItem { { std::move(snps), region } } );

                    if (os_meta) {
                        auto sorted_cov = tmp_vec;
                        sorted_cov.resize(gene_obs.m_gene_length, 0);
                        std::sort(sorted_cov.begin(), sorted_cov.end());
                        size_t n = sorted_cov.size();
                        median_vcov = n % 2 == 1 ?
                            static_cast<double>(sorted_cov[n/2]) :
                            (static_cast<double>(sorted_cov[n/2 - 1]) + static_cast<double>(sorted_cov[n/2])) / 2.0;

                        auto nonzero_begin = std::lower_bound(sorted_cov.begin(), sorted_cov.end(), 1);
                        size_t nz = std::distance(nonzero_begin, sorted_cov.end());
                        if (nz > 0) {
                            double sum = std::accumulate(nonzero_begin, sorted_cov.end(), 0.0);
                            mean_vcov_nonzero = sum / nz;
                            size_t mid = nz / 2;
                            median_vcov_nonzero = nz % 2 == 1 ?
                                static_cast<double>(*(nonzero_begin + mid)) :
                                (static_cast<double>(*(nonzero_begin + mid - 1)) + static_cast<double>(*(nonzero_begin + mid))) / 2.0;
                        }
                    }

                    if (os_meta) {
#pragma omp critical(metaout)
                        {
                            *os_meta << profile.GetName() << '\t';
                            *os_meta << geneid << '\t';
                            *os_meta << gene_obs.VerticalCoverage() << '\t';
                            *os_meta << counts_vcov1 << '\t';
                            *os_meta << counts_vcov2 << '\t';
                            *os_meta << multi_allelic << '\t';
                            *os_meta << ac.Filtered() << '\t';
                            *os_meta << (counts_vcov1 > 0 ? multi_allelic/static_cast<double>(counts_vcov1) : 0) << '\t';
                            *os_meta << (counts_vcov1 > 0 ? ac.Filtered()/static_cast<double>(counts_vcov1) : 0) << '\t';
                            *os_meta << (counts_vcov2 > 0 ? multi_allelic/static_cast<double>(counts_vcov2) : 0) << '\t';
                            *os_meta << (counts_vcov2 > 0 ? ac.Filtered()/static_cast<double>(counts_vcov2) : 0) << '\t';
                            *os_meta << median_vcov << '\t';
                            *os_meta << hcov << '\t';
                            *os_meta << gene_obs.m_gene_length << '\t';
                            *os_meta << mean_vcov_nonzero << '\t';
                            *os_meta << median_vcov_nonzero;
                            *os_meta << std::endl;
                        }
                    }
                }
            }
            if (samples_with_gene > 0) {
                previous_size = msa.front().size();

                protal::MSAStats gene_stats(items.size());
                bool result = protal::MSA(items, gene.Sequence(), msa, min_cov, min_qual_sum, min_af, require_strand, min_mean_qual, &gene_stats, &ref_msa_row, snp_max_alleles);

                if (!result) continue;

                for (size_t si = 0; si < gene_stats.size(); si++) {
                    sample_stats[si] += gene_stats[si];
                }
                // Partitions are 1-based and inclusive, as RAxML and IQ-TREE read them.
                if (msa.front().size() > partition_start) {
                    if (partitions.size() > 0) {
                        partitions.back() += std::to_string(previous_size);
                    }

                    std::string partition = "DNA, gene";
                    partition += std::to_string(geneid) + " = ";
                    partition += std::to_string(partition_start + 1) + '-';
                    partitions.emplace_back(partition);
                    partition_start = msa.front().size();
                }
            }
        }

        if (partitions.empty()) {
            std::cout << "No gene of " << taxon_name << " has enough coverage for an MSA" << std::endl;
            return;
        }

        // Prepend reference sequence as the first row so it is always present in both outputs.
        msa.insert(msa.begin(), std::move(ref_msa_row));
        names.insert(names.begin(), taxon_name + "_reference");
        sample_stats.insert(sample_stats.begin(), protal::MSASampleStats{});

        bool any_good = std::any_of(msa.begin(), msa.end(), [min_hcov](MSARow const& row){
            return IsRowGood(row, min_hcov);
        });
        if (!any_good) {
            std::cout << "No good consensus sequences found for species" << std::endl;
            return;
        }

        std::ofstream os(options.GetMSAOutput(taxon_name), std::ios::out);
        for (auto i = 0; i < msa.size(); i++) {
            auto& row = msa[i];

            if (!IsRowGood(row, min_hcov)) continue;
            os << ">" << names[i] << std::endl;

            // os << std::string_view(&row[0], std::distance(row.begin(), row.end())) << std::endl;
            os << std::string(&row[0], std::distance(row.begin(), row.end())) << std::endl;
        }
        os.close();
        if (os.fail()) RunStatus::Get().Fail("Writing the MSA of " + taxon_name + " failed: " + options.GetMSAOutput(taxon_name));
        // std::cout << " Saved MSA to " << options.GetMSAOutput(taxon_name);

        partitions.back() += std::to_string(msa.front().size());

        std::ofstream os_part(options.GetMSAPartitionOutput(taxon_name), std::ios::out);
        for (auto i = 0; i < partitions.size(); i++) {
            os_part << partitions[i] << std::endl;
        }
        os_part.close();
        if (os_part.fail()) RunStatus::Get().Fail("Writing the MSA partitions of " + taxon_name + " failed: " + options.GetMSAPartitionOutput(taxon_name));

        // std::cout << " Saved Partitions " << std::endl;

        // Strain output is the raw MSA (.raw.msa.fna) plus qcmsa's final
        // .msa.fna (+ their partition files). The genecol_filtered and processed
        // MSA variants were removed: per-gene multi-allelicity and site-level
        // filtering are handled by the qcmsa post-filter on the raw MSA.

        // Write per-sample SNP-retention statistics TSV.
        {
            std::ofstream os_stats(options.GetMSAStatsOutput(taxon_name), std::ios::out);
            os_stats << "sample"
                     // --- totals ---
                     << "\ttotal_variant_positions"
                     << "\ttotal_pass_snps"
                     << "\ttotal_filtered_snps"
                     // --- per-type pass counts + % of total_variant_positions ---
                     << "\tsnps_retained\tsnps_retained_pct"
                     << "\tinsertions_retained\tinsertions_retained_pct"
                     << "\tdeletions_retained\tdeletions_retained_pct"
                     // --- filter breakdown + % of total_variant_positions ---
                     << "\tvariants_filtered_qual_sum\tvariants_filtered_qual_sum_pct"
                     << "\tvariants_filtered_obs_cov\tvariants_filtered_obs_cov_pct"
                     << "\tvariants_filtered_af\tvariants_filtered_af_pct"
                     << "\tvariants_filtered_strand\tvariants_filtered_strand_pct"
                     // --- pass/filter summary percentages ---
                     << "\ttotal_pass_pct\ttotal_filtered_pct"
                     // --- position-level counts + % of total_positions ---
                     << "\tpositions_ref\tpositions_ref_pct"
                     << "\tpositions_below_min_cov\tpositions_below_min_cov_pct"
                     << "\tpositions_no_coverage\tpositions_no_coverage_pct"
                     // --- appended last so readers that index columns by position keep working ---
                     << "\trefs_retained\trefs_retained_pct"
                     << '\n';

            os_stats << std::fixed << std::setprecision(2);
            for (size_t si = 0; si < sample_stats.size(); si++) {
                auto const& s = sample_stats[si];
                os_stats << names[si]
                         << '\t' << s.TotalVariantPositions()
                         << '\t' << s.TotalPass()
                         << '\t' << s.TotalFiltered()
                         << '\t' << s.snps_retained              << '\t' << s.PctSnpsRetained()
                         << '\t' << s.insertions_retained        << '\t' << s.PctInsertionsRetained()
                         << '\t' << s.deletions_retained         << '\t' << s.PctDeletionsRetained()
                         << '\t' << s.variants_filtered_qual_sum << '\t' << s.PctFilteredQualSum()
                         << '\t' << s.variants_filtered_obs_cov  << '\t' << s.PctFilteredObsCov()
                         << '\t' << s.variants_filtered_af       << '\t' << s.PctFilteredAF()
                         << '\t' << s.variants_filtered_strand   << '\t' << s.PctFilteredStrand()
                         << '\t' << s.PctPass()
                         << '\t' << s.PctFiltered()
                         << '\t' << s.positions_ref              << '\t' << s.PctPositionsRef()
                         << '\t' << s.positions_below_min_cov    << '\t' << s.PctPositionsBelowMinCov()
                         << '\t' << s.positions_no_coverage      << '\t' << s.PctPositionsNoCoverage()
                         << '\t' << s.refs_retained              << '\t' << s.PctRefsRetained()
                         << '\n';
            }
            os_stats.close();
        }
    }



    static SimilarityMatrix GetSimilarityMatrixForTaxon(uint32_t taxid, Options& options, Profiles& profiles, std::optional<profiler::TaxonFilterObj> filter={}) {
        SimilarityMatrix matrix;

        size_t min_shared_length = 1000;

        size_t total_combinations = profiles.size() * profiles.size() - profiles.size() / 2;
        size_t count_combinations = 1;

        ProgressBar prog(((profiles.size() * profiles.size()) - profiles.size()) / 2);
        for (auto i = 0; i < profiles.size(); i++) {
            auto& profile1 = profiles[i];
            auto name1 = profile1.GetName();

            for (auto j = i+1; j < profiles.size(); j++) {
                prog.UpdateAdd(1);
                count_combinations++;
//                std::cout << "\rCombination " << count_combinations << " of " << total_combinations;
                auto& profile2 = profiles[j];
                auto name2 = profile2.GetName();

                auto [similarity, shared_length] = GetSimilarity(taxid, profile1, profile2, options, filter);

//                std::cout << name1 << " -- " << name2 <<  " = " << similarity << " over " << shared_length << std::endl;

                if (similarity == NAN || shared_length < min_shared_length) continue;

                if (!matrix.HasName(name1)) {
                    matrix.AddName(name1);
                }
                if (!matrix.HasName(name2)) {
                    matrix.AddName(name2);
                }

//                std::cout << taxid << " " << similarity << " " << shared_length << std::endl;

                matrix.SetValue(name1, name2, similarity, true);
                matrix.SetValue(name1, name2, static_cast<double>(shared_length), false);
            }
        }
        std::cout << std::endl;
        return matrix;
    }



    static void WriteDistanceMatrix(uint32_t id, SimilarityMatrix const& matrix, Options& options, std::string& name) {
        std::ofstream matrix_os(options.GetSimilarityMatrixOutput(name), std::ios::out);
        matrix.PrintMatrix(matrix_os, "\t", 10);
        matrix_os.close();
    }

    static void StrainWrapper2(Options& options, Profiles& profiles, GenomeLoader& loader, taxonomy::IntTaxonomy& taxonomy, std::vector<uint32_t> msa_taxids = {}, std::optional<profiler::TaxonFilterObj> filter={}) {
        Benchmark bm_strain{"Strain-level MSAs"};
        bm_strain.Start();
        std::cout << "Output " << options.GetOutputDir() << std::endl;
        auto dir = std::filesystem::path(options.GetOutputDir());
        if (!std::filesystem::create_directories(dir.string()) && !std::filesystem::exists(dir)) {
            std::cout << "Cannot create directories for this path " << dir << std::endl;
            exit(2);
        };

        auto taxids = msa_taxids.empty() ? ExtractTaxa(profiles, filter) : msa_taxids;

        auto enable_similarity_matrix = false;


        for (auto& taxid : taxids) {
            std::cout << taxonomy.Get(taxid).scientific_name << std::endl;

            std::string name = taxonomy.Get(taxid).scientific_name;
            std::replace(name.begin(), name.end(), ' ', '_');

            if (enable_similarity_matrix) {
                auto similarities = GetSimilarityMatrixForTaxon(taxid, options, profiles, filter);
                if (!similarities.AnySet()) continue;
                WriteDistanceMatrix(taxid, similarities, options, name);
            }

            std::ofstream os_meta(options.GetSpeciesMetaOutput(name));
            os_meta << "sample\tgene_id\tvertical_coverage\tcounts_vcov1\tcounts_vcov2\tmulti_allelic\tfiltered\tmulti_rate_vcov1\tfiltered_rate_vcov1\tmulti_rate_vcov2\tfiltered_rate_vcov2\tmedian_vcov\thcov\tgene_length\tmean_vcov_nonzero\tmedian_vcov_nonzero\n";
            GetMSAForTaxon(taxid, name, loader, options, profiles, &os_meta);
            os_meta.close();
            if (os_meta.fail()) RunStatus::Get().Fail("Writing the MSA metadata of " + name + " failed: " + options.GetSpeciesMetaOutput(name));

            if (options.GetRunQCMSA()) {
                RunQCMSA(options, name);
            }
//            Utils::Input();
        }
        bm_strain.Stop();
        bm_strain.PrintResults();
    }


    size_t GetTotalSystemMemory()
    {
        long pages = sysconf(_SC_PHYS_PAGES);
        long page_size = sysconf(_SC_PAGE_SIZE);
        return pages * page_size;
    }

    // Runs protal and returns the process exit code: 0, or 1 if any sample or output failed (see
    // RunStatus). Invalid input and fatal errors still exit directly with their own codes.
    static int Run(int argc, char *argv[]) {
        auto options = protal::Options::OptionsFromArguments(argc, argv);
        if (options.ShowVersion()) {
            std::cout << "protal v" << 
                protal_VERSION_MAJOR << "." <<
                protal_VERSION_MINOR << "." <<
                protal_VERSION_PATCH << std::endl;
            exit(0);
        }
    

        PrintLogo();
        PrintProtalInformation();
        std::cout << "Total available memory is " << GetTotalSystemMemory() / (1024 * 1024 * 1024) << "GB" << std::endl;
        std::cout << std::endl;

        Benchmark bm_total("Run protal");
        bm_total.Start();

        using AlignmentBenchmark = CoreBenchmark;


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

        // Load protal DB into RAM
        // reference.fna or reference.fna.zst
        ProtalDB db = options.UniqueKmersFileExists() ?
            ProtalDB(options.ResolvedSequenceFile(), options.GetSequenceMapFile(), options.GetHittableGenesMap(), options.GetUniqueKmersFile()) :
            ProtalDB(options.ResolvedSequenceFile(), options.GetSequenceMapFile());

        // Load fasta sequences of reference into RAM (advised)
        if (options.PreloadGenomes()) {
            std::cout << "Preload genomes" << std::endl;
            Benchmark bm_preload_genomes("Preload genomes");
            bm_preload_genomes.Start();
            db.GetGenomes().LoadAllGenomes();
            bm_preload_genomes.Stop();
            bm_preload_genomes.PrintResults();
        }

        // Skip alignment if files are present. Do not skip if either files are not there or user specified --force
        auto sam_files = options.SamFiles();
        bool all_alignments_exist = std::all_of(sam_files.begin(), sam_files.end(), [](std::string const& file){ return Utils::exists(file); });

        bool skip_alignment = !options.BuildMode() && all_alignments_exist && !options.Force();

        bool run_alignment = true;
        if (!options.BuildMode() && (options.ProfileOnly() || skip_alignment) && !sam_files.empty() &&
            !options.SamFile(0).first.empty()) {
            std::cout << "All alignments are present." << std::endl;
            run_alignment = false;
        }
        bool const run_profiling = !options.BuildMode() && !options.NoProfile();

        // Checks that would otherwise only fail after hours of alignment.
        std::vector<uint32_t> msa_taxids;
        if (run_profiling) {
            db.LoadTaxonomy(options.GetInternalTaxonomyFile());
            msa_taxids = ResolveMSASpecies(options, db.GetTaxonomy());
            try {
                profiler::TaxonFilterObj model_check(options.GetModelPath(), options.GetKnob());
            } catch (std::exception const& e) {
                std::cerr << "Cannot load the model " << options.GetModelPath() << ": " << e.what() << std::endl;
                exit(2);
            }
        }
        if (run_alignment && options.BenchmarkAlignment() && !options.GetRange().empty()) {
            // The benchmark takes each read's true gene from its name; without one it would stop
            // at the first read, deep inside the alignment.
            igzstream is{ options.GetFirstFile(options.GetRange().front()).c_str() };
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
        if (run_alignment && !options.BuildMode() && FindInPath("pigz").empty() &&
            std::any_of(sam_files.begin(), sam_files.end(), [](std::string const& f) { return f.ends_with(".gz"); })) {
            std::cerr << "Error: gzipped SAM output (.sam.gz) needs pigz, which is not on $PATH. "
                         "Install pigz or request uncompressed .sam files." << std::endl;
            exit(2);
        }

        /*
         *  READ ALIGNMENT SECTION
         */
        // Untangle Template options that need to be written out specifically.
        if (run_alignment) {
            if (options.BenchmarkAlignment()) {
                AlignmentBenchmark alignment_benchmark{};
                if (!options.GetBenchmarkAlignmentOutputFile().empty()) {
                    alignment_benchmark.SetOutput(options.GetBenchmarkAlignmentOutputFile());
                }
                RunWrapper(options, db, alignment_benchmark);
                if (!options.GetBenchmarkAlignmentOutputFile().empty()) {
                    alignment_benchmark.DestroyOutput();
                }
            } else {
                RunWrapper(options, db);
            }
        }

        /*
         * PROFILER
         */
        if (run_profiling) {
            Benchmark bm_profiling("Profiling");
            bm_profiling.Start();
            
            auto profiles = ProfileWrapper(options, db);
            bm_profiling.Stop();
            bm_profiling.PrintResults();


            // Remove later - keep option
            double min_ani = 0.95;
            double min_gene_presence = 0.5; //previously 0.5
            size_t min_total_hits = 60; //previously 70
            size_t min_mean_mapq = 10;

//            profiler::TaxonFilter filter(min_ani, min_gene_presence, min_total_hits, min_mean_mapq);

            std::string model = options.GetModelPath();
            // std::cout << "Model: " << model << std::endl;

            profiler::TaxonFilterObj filter(model, options.GetKnob());

            // If there is more than one profile, get per taxon output
            if (profiles.size() > 1) {
                auto taxids = ExtractTaxa(profiles);
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
                        bool accepted = filter.Pass(taxon);
                        stats.PrintLine(os, profile.GetName(), taxon.VerticalCoverage(), taxon.TotalHits(), taxon.TotalLength(), taxon.GetMeanANI(), taxon.GetMeanMAPQ(), accepted);
                    }
                    os.close();
                }
            }

            /*
             * STRAIN PART -  RESOLVE MSAs BETWEEN SAMPLES
             */
            if (!options.NoStrains()) {
                StrainWrapper2(options, profiles, db.GetGenomes(), db.GetTaxonomy(), msa_taxids, filter);
            }
        }

        bm_total.Stop();
        bm_total.PrintResults();

        std::cout << std::flush;
        return RunStatus::Get().Finish();
    }
}
