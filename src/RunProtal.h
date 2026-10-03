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
#include "ThreadedGzStream.h"
#include "AlignmentStrategy.h"
#include "TaxonStatisticsOutput.h"
#include "ProgressBar.h"
#include "protal_config.h"
#include "SamFile.h"
#include <thread>
#include <mutex>
#include <condition_variable>
#include <deque>
#include <functional>
#include <optional>
#include "RunStatus.h"
#include "BuildInfo.h"

#include <atomic>
#include <iomanip>
#include <regex>
#include <ranges>
#include <set>
#include <unistd.h>
#include <sys/wait.h>

// #include "Profiler/ReadFilter.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <numeric>

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
        // The gene tables (reference.map, unique_kmers.tsv) are read with `threads` threads.
        ProtalDB(db::DbFile sequence_file, db::DbFile map_file, int threads = 1) :
                m_genomes(std::move(sequence_file), std::move(map_file), threads),
                m_taxonomy() {
        }

        ProtalDB(db::DbFile sequence_file, db::DbFile map_file, db::DbFile const& unique_kmers_file, int threads = 1) :
                m_genomes(std::move(sequence_file), std::move(map_file), threads),
                m_taxonomy() {
            m_genomes.LoadUniqueKmers(unique_kmers_file, threads);
        }

        void LoadTaxonomy(db::DbFile const& file) {
            auto input = file.Open();
            if (!file.Exists() || !input->IsOpen()) {
                std::cerr << "Invalid taxonomy " << file.Name() << ": cannot open the file" << std::endl;
                exit(8);
            }
            m_taxonomy = taxonomy::IntTaxonomy(input->Stream(), file.Name());
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

    // Alignments are written under a temporary name ("<sam>.partial", compressed already as the name
    // asks; SamOutput) and only get their final name once complete, so an interrupted run never leaves
    // a truncated SAM that a rerun would skip and reuse. This moves a finished file into place.
    static void FinishSamFile(Options& options, int index, std::string const& partial) {
        auto [sam, _] = options.SamFile(index);
        std::error_code ec;
        std::filesystem::rename(partial, sam, ec);
        if (ec) RunStatus::Get().Fail("Cannot move " + partial + " to " + sam + ": " + ec.message());
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

    // The genes' conservation factors (Options::GeneConservationDbFile, GeneConservation.h: the database's
    // gene_conservation.tsv, or --gene_conservation's file). They give the model's conservation features
    // (conserved_fast_depth_ratio, conserved_hit_share) always, and scale the depth identity margin per gene with
    // --gene_conservation db or FILE; by default (none) every gene keeps the whole margin. Without the file the
    // features are 0 and 0.5. Exits 8 if the file cannot be read.
    static void LoadGeneConservation(Options const& options, GenomeLoader& genomes) {
        auto const file = options.GeneConservationDbFile();
        bool const scale = options.ScaleDepthMarginByConservation();
        std::string const same = "the depth identity margin is the same on every gene";
        if (!file.Exists()) {
            std::cout << "Gene conservation: the database has no " << Options::PROTAL_GENE_CONSERVATION_FILE << " (built by an "
                      << "earlier protal, or its --full_reference had no other genomes' copies of the genes): the "
                      << "conservation features are 0 and 0.5, and " << same << std::endl;
            return;
        }
        std::string error;
        auto const content = file.ReadAll(error);
        gene_conservation::Table table;
        if (content) {
            std::istringstream is(*content);
            error = table.Read(is);
        }
        if (!error.empty()) {
            std::cerr << "Invalid gene conservation factors " << file.Name() << ": " << error << std::endl;
            exit(8);
        }
        auto const [low, high] = table.Range();
        std::cout << "Gene conservation: factors " << std::setprecision(2) << low << "-" << high << std::setprecision(6) << " for "
                  << table.Genes() << " genes (" << file.Name() << "), for the conservation features; "
                  << (scale ? "they scale the depth identity margin per gene"
                            : same + " (--gene_conservation db scales it by them)") << std::endl;
        genomes.SetGeneConservation(std::move(table));
        genomes.SetScaleDepthMargin(scale);
    }

    // The database's gene neighbours (gene_neighbours.tsv, GeneNeighbours.h: how often each marker gene end faces
    // which other in a clade's genomes), for mate guidance past a gene's end, pairs of mates on neighbouring genes,
    // the genes next to a long read's genes and the profiler's adjacency features; none without the file or with
    // --no_gene_neighbours. Only these per-clade frequencies are loaded: the per-genome positions they were counted
    // from (gene_positions.tsv) stay in the database unread. Every species of the taxonomy gets the clades of its
    // lineage. Exits 8 if the file
    // cannot be read.
    static void LoadGeneNeighbours(Options const& options, ProtalDB& db) {
        auto const file = options.GeneNeighboursDbFile();
        if (!file.Exists()) {
            std::cout << "Gene neighbours: the database has none (" << Options::PROTAL_GENE_NEIGHBOURS_FILE << ", from whole "
                      << "genomes by scripts/mini_db/gene_neighbours.py before --build)" << std::endl;
            return;
        }
        if (!options.GeneNeighbours()) {
            std::cout << "Gene neighbours: not used (--no_gene_neighbours)" << std::endl;
            return;
        }
        std::string error;
        auto const content = file.ReadAll(error);
        gene_neighbours::Table table;
        if (content) {
            std::istringstream is(*content);
            error = table.Read(is);
        }
        if (!error.empty()) {
            std::cerr << "Invalid gene neighbours " << file.Name() << ": " << error << std::endl;
            exit(8);
        }
        if (!db.IsTaxonomyLoaded()) db.LoadTaxonomy(options.TaxonomyDbFile());
        auto& taxonomy = db.GetTaxonomy();
        std::vector<uint32_t> lineage;
        for (auto const& [id, node] : taxonomy.map) {
            if (node.rank != "species") continue;
            lineage.assign(1, static_cast<uint32_t>(id));
            int t = node.parent_id;
            while (t >= 0 && taxonomy.HasNode(t) && lineage.size() < 64) {
                lineage.push_back(static_cast<uint32_t>(t));
                int const parent = taxonomy.Get(t).parent_id;
                if (parent == t) break;  // the root
                t = parent;
            }
            table.SetLineage(static_cast<uint32_t>(id), lineage);
        }
        std::cout << "Gene neighbours: " << table.Rules() << " rules of " << table.Clades() << " clades from " << table.Genomes()
                  << " genomes (" << file.Name() << "), for " << table.BoundSpecies() << " species: mates looked for past "
                  << "their gene's end, pairs across neighbouring genes, long reads' neighbouring genes" << std::endl;
        db.GetGenomes().SetGeneNeighbours(std::move(table));
    }

    // Aligns the samples. `sample_done` (if any) is called with a sample's index once its SAM file is complete, or
    // when the sample is skipped because its SAM exists: Run profiles such samples while the next ones are aligned
    // (ProfilingAhead).
    template<typename AlignmentBenchmark=NoBenchmark>
    static void RunWrapper(Options& options, ProtalDB& db, AlignmentBenchmark benchmark=NoBenchmark{},
                           std::function<void(size_t)> const& sample_done = {}) {

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
                    options, kmer_putter, iterator, db.GetGenomes());

            // Last, as the build reads reference.fna until here: the single file (which compresses
            // reference.fna itself), or separate files.
            if (options.WriteBundle()) {
                protal::build::BundleDatabase(options);
            } else {
                protal::build::CompressReference(options);
                protal::build::RemoveStaleBundle(options);
            }

            bm_build.PrintResults();
            protal_stats.WriteStats(std::cout);

        } else {

            // Benchmark Load Time
            Benchmark bm_load_index("Load Index");
            bm_load_index.Start();
            // Load Index
            Seedmap map;
            auto const index_file = options.IndexDbFile();
            int const load_threads = static_cast<int>(options.GetThreads());
            bool const single_frame = !index_file.InBundle() && zstd::IsCompressed(index_file.Path()) && !zstd::IsSeekable(index_file.Path());
            std::cout << "Load index " << index_file.Name() << " (" << load_threads << " thread(s)"
                      << (single_frame ? "; a single zstd frame is read with one" : "") << ")" << std::endl;
            map.Load(index_file, load_threads);
            std::cout << "Index features: " << map.FeatureDescription() << std::endl;
            auto const map_file = options.SequenceMapDbFile(), fna_file = options.SequenceDbFile();
            if (map.HasReferenceFingerprint() && !(map.GetReferenceFingerprint() == ReferenceFingerprint::Of(map_file, fna_file))) {
                std::cerr << "index.prx was built against a different reference: " << map_file.Name() << " or "
                          << fna_file.Name() << " changed since the index was built. Rebuild the index (--build) or "
                          << "restore the reference files it was built with." << std::endl;
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
            // Long reads are aligned without X-drop: over their gene-long windows it loses the own
            // species' alignment of genes the read ends in, and leaves a relative's weaker one unique.
            WFA2Wrapper2 long_read_wfa(4, 6, 2, 0);

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


            Benchmark bm_classify("Processing all samples");
            bm_classify.Start();
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

                auto dir = std::filesystem::path(sam).parent_path();

                if (!std::filesystem::create_directories(dir.string()) && !std::filesystem::exists(dir)) {
                    std::cout << "Cannot create directories for this path " << sam << std::endl;
                    exit(32);
                };


                // std::cout << index << " Process sample " << options.GetSampleId(index) << (std::filesystem::exists(sam) ? " (sam exists)" : " (sam does not exist)") << std::endl;

                // Avoid aligning files that already exist.
                if (!options.Force() && (std::filesystem::exists(sam) || std::filesystem::exists(sam_plain))) {
                    std::cout << "Skip " << sam << " continue" << std::endl;
                    // Profile the file that is there: an earlier run may have left it uncompressed.
                    if (compressed && !std::filesystem::exists(sam)) options.UseUncompressedSamFile(index);
                    if (sample_done) sample_done(index);
                    continue;
                }

                // AnchorFinder. A long read is seeded from all its k-mers: stopping at -s seeds, as for a
                // short read, would leave most of it unseeded.
                AnchorFinder anchor_finder(kmer_lookup, mmer_size, options.GetMinSuccessfulLookups(),
                                           IsLongReadType(read_type) ? SIZE_MAX : options.GetMaxSeedSize(), genomes);
                // AlignmentHandler approach
                double const max_score_ani = options.GetMaxScoreAni(read_type);
                SimpleAlignmentHandler alignment_handler(genomes, IsLongReadType(read_type) ? long_read_wfa : aligner, kmer_size,
                                                         options.GetAlignTop(), max_score_ani, options.FastAlign());
                // Reads from their anchors' exact matches (AnchoredAligner); long reads, whose seeds lie on
                // diagonals that their indels shift, through every link of the chain.
                alignment_handler.SetAnchoredAlignment(!options.WholeReadAlignment());
                alignment_handler.SetAnchoredIndels(IsLongReadType(read_type));



                options.SetCurrentIndex(index);
                std::string const sam_partial = sam + ".partial";  // renamed by FinishSamFile
                std::string const read_type_line = kSamReadTypeComment + Info(read_type).token + '\n';
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
                std::vector<std::string> read_file_list{ options.GetFirstFile(index) };
                if (!single_file) read_file_list.push_back(options.GetSecondFile(index));
                SamOutput sam_output(sam_partial, SamCompressionOf(sam), full_header,
                                     sam_zstd::HeaderRoom(database_genes, input_bytes(read_file_list)));
                if (!sam_output.Ok()) {
                    RunStatus::Get().Fail("Cannot write the SAM file of sample " + options.GetSampleId(index) + ": " + sam_output.Error());
                    continue;
                }
                std::cout << "Align the " << ReadTypeName(read_type) << " reads of sample "
                          << options.GetSampleId(index) << " (-a " << max_score_ani << ")" << std::endl;

                // Main Run Call. This is where the reads are read and alignment happens
                bool truncated = false;
                bool read_success = true;
                std::string read_problem;   // why reading failed, if it says
                size_t reads_read = 0;
                std::string const read_files = single_file ? options.GetFirstFile(index) :
                                               options.GetFirstFile(index) + ", " + options.GetSecondFile(index);
                // An input that cannot be opened reads as empty: say why instead.
                auto cannot_read = [&](ThreadedGzIstream& is, std::string const& path) {
                    if (is.rdbuf()->is_open()) return false;
                    RunStatus::Get().Fail("Cannot read " + path + " (sample " + options.GetSampleId(index) + "): " +
                                          is.rdbuf()->open_error() + "; no SAM file was written");
                    sam_output.Discard();
                    return true;
                };
                auto read_error = [](ThreadedGzIstream& is) {
                    auto const message = is.rdbuf()->read_error_message();
                    return message.empty() ? std::string() : ": " + message;
                };
                if (IsLongReadType(read_type)) {
                    ThreadedGzIstream is { options.GetFirstFile(index).c_str() };
                    if (cannot_read(is, options.GetFirstFile(index))) continue;
                    SeqReaderSE reader{ is, FastaQualityChar(read_type) };
                    LongReadAligner<SimpleKmerHandler<ClosedSyncmer>, AnchorFinder> long_read_aligner(
                            iterator, anchor_finder, alignment_handler, genomes, options.GetAlignTop(), max_score_ani);
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
                    reads_read = protal_stats.reads;
                    // zlib-ng and libdeflate read a truncated or corrupt gzip file as one that ends early.
                    truncated = is.rdbuf()->read_failed();
                    read_problem = read_error(is);
                    read_success = reader.Success();
                    is.close();
                } else if (read_type == ReadType::Single) {
                    ThreadedGzIstream is { options.GetFirstFile(index).c_str() };
                    if (cannot_read(is, options.GetFirstFile(index))) continue;
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
                    reads_read = protal_stats.reads;
                    // zlib-ng and libdeflate read a truncated or corrupt gzip file as one that ends early.
                    truncated = is.rdbuf()->read_failed();
                    read_problem = read_error(is);
                    read_success = reader.Success();
                    if (reader.TooLong() > 0) {
                        read_problem = ": read " + reader.TooLongId() + " has " + std::to_string(reader.TooLong()) +
                                       " bp, too long for short reads; give --read_type pb or ont (or in the map's READ_TYPE "
                                       "column) for PacBio or ONT reads";
                    }
                    is.close();
                } else {
                    // Each file inflates in a thread of its own (ThreadedGzStream.h), outside the reader lock.
                    ThreadedGzIstream is1 { options.GetFirstFile(index).c_str() };
                    ThreadedGzIstream is2 { options.GetSecondFile(index).c_str() };
                    if (cannot_read(is1, options.GetFirstFile(index)) || cannot_read(is2, options.GetSecondFile(index))) continue;
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
                    reads_read = protal_stats.reads;
                    // zlib-ng and libdeflate read a truncated or corrupt gzip file as one that ends early.
                    truncated = is1.rdbuf()->read_failed() || is2.rdbuf()->read_failed();
                    read_problem = read_error(is1);
                    if (read_problem.empty()) read_problem = read_error(is2);
                    read_success = reader.Success();
                    if (reader.NameMismatches() > 0 && read_success && !truncated) {
                        std::cerr << "Warning: sample " << options.GetSampleId(index) << ": the mates of " << reader.NameMismatches()
                                  << " read pair(s) have different names (e.g. " << reader.FirstNameMismatch()
                                  << "); are " << read_files << " in the same order?" << std::endl;
                    }
                    is1.close();
                    is2.close();
                }
                bm_classify_sample.Stop();
                bm_classify_sample.PrintResults();

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
                Benchmark bm_finish_sam("Writing the SAM header and file");
                bm_finish_sam.Start();
                std::string header;
                if (!full_header) {
                    std::ostringstream os;
                    genomes.WriteSamHeader(os, sam_output.Genes());
                    os << read_type_line;
                    header = os.str();
                }
                bool const written = sam_output.Finish(header);
                bm_finish_sam.Stop();
                bm_finish_sam.PrintResults();
                if (!written) {
                    RunStatus::Get().Fail("Writing the SAM file of sample " + options.GetSampleId(index) + " failed: " + sam_output.Error());
                    continue;
                }
                FinishSamFile(options, index, sam_partial);
                if (sample_done) sample_done(index);
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

    // The model of each kind of reads (ReadType), loaded if the samples have such reads.
    using ReadTypeModels = std::array<std::optional<profiler::TaxonFilterObj>, kReadTypeCount>;

    // A taxon that a sample's profile leaves out (score below --knob) although its own reads are strong
    // evidence that it is present (profiler::StrongOwnEvidence).
    struct UnreportedTaxon {
        std::string species;
        std::string row;  // of unreported_species.tsv
    };

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

    // What profiling the samples of a run shares: the taxonomy, every taxid's genus and family (which of a read's
    // alternatives (ZA) are congeners of its taxon, and what the other taxa of a sample say of each:
    // MicrobialProfile::ApplySampleContext), the distances of the references (made as they are asked for, kept for
    // the run), the models, and a slot per sample for its profile and for the taxa it leaves out although their own
    // reads are strong evidence. Made before the alignment when samples are profiled while others are aligned
    // (ProfilingAhead), else by the profiling stage.
    struct ProfilingContext {
        Options& options;
        GenomeLoader& genomes;
        taxonomy::IntTaxonomy* taxonomy = nullptr;
        std::shared_ptr<std::vector<uint32_t> const> genera, families;
        std::shared_ptr<profiler::context::CongenerDistances> distances;
        ReadTypeModels const& models;
        std::vector<size_t> range;  // the samples; a sample's slot is its position here
        // MicrobialProfile holds a reference member so it is not assignable; a slot is filled in place.
        std::vector<std::optional<profiler::MicrobialProfile>> profile_slots;
        std::vector<std::vector<UnreportedTaxon>> unreported_slots;

        ProfilingContext(Options& options, ProtalDB& db, ReadTypeModels const& models) :
                options(options), genomes(db.GetGenomes()), models(models), range(options.GetRange()),
                profile_slots(range.size()), unreported_slots(range.size()) {
            if (!db.IsTaxonomyLoaded()) db.LoadTaxonomy(options.TaxonomyDbFile());
            taxonomy = &db.GetTaxonomy();
            genera = profiler::GeneraOf(*taxonomy);
            families = profiler::FamiliesOf(*taxonomy);
            distances = std::make_shared<profiler::context::CongenerDistances>(genomes);
        }
    };

    // Profiles the sample of slot idx with the model of its kind of reads, on threads_per_sample threads, into
    // ctx.profile_slots[idx] and ctx.unreported_slots[idx]; `filters` is the calling thread's own copy of the
    // models (scoring reuses a buffer; copies share the loaded model). Every taxon's score is then cached
    // (profile.ReleaseReadData scores all), so that later stages may use any model to tell which taxa pass: they
    // get the score of the sample's own model. Any thread may profile a sample of its own: a sample's profile
    // depends on its SAM and the database only, not on the thread or on the other samples.
    static void ProfileSample(ProfilingContext& ctx, size_t idx, ReadTypeModels& filters, size_t threads_per_sample) {
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
            auto& sample_filter = filters[static_cast<size_t>(read_type)];
            if (!sample_filter) {  // loaded up front for every kind of reads the samples have
                RunStatus::Get().Fail("No model loaded for the " + ReadTypeName(read_type) + " reads of sample " + options.GetSampleId(i));
                profile_slots[idx].emplace(genomes);
                return;
            }
            auto& filter = *sample_filter;

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

            std::ofstream erro(sam + ".err", std::ios::out);
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
            profile.SetDropForeignGenes(options.DropForeignGenes());
            // A long-read sample's strains get strain MSA rows of their own (Haplotypes.h).
            profile.SetKeepPhaseRecords(!options.NoStrains() && options.Phasing() && IsLongReadType(read_type));
            std::string sam_error = profiler.ProfileSam(sam, profile, std::optional<std::reference_wrapper<std::ostream>>{erro},
                                                        options.GetSNPMinCov(), options.GetSNPMinCov(),
                                                        options.GetSNPMinAF(read_type), options.GetSNPMinMeanQual(),
                                                        options.GetSNPMinPhredSum(), options.GetSNPRequireStrand(),
                                                        threads_per_sample);
            erro.close();
            bm_profile.Stop();

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
            // With a model's calibration (random_forest_cmdline.py --fdr-calls), unless --knob or --fdr 0: the knob at
            // which the sample's expected share of false calls is at most the model's target or --fdr's
            // (context::FalseCallKnob), over its taxa the singleton rule does not veto.
            filter.SetKnob(options.GetKnob());
            bool const false_call_mode = !options.KnobGiven() && filter.HasFalseCalls() &&
                                         (!options.FdrGiven() || options.GetFdr() > 0);
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
            if (threads_per_sample > 1) profile.ScoreTaxa(filter, threads_per_sample);
            if (profiler.RejectedReads() > 0) {
                #pragma omp critical(print)
                std::cerr << "Warning: sample " << sample_name << ": " << profiler.RejectedReads() << " of " << profiler.Reads()
                          << " reads have an alignment that does not fit the database (a gene it lacks, a position past "
                             "a gene's end, or bases that differ from the gene); they are left out and listed in "
                          << sam << ".err" << std::endl;
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

            profile.WriteSparseProfile(taxonomy, filter, os, &os_total, &os_dismissed, threads_per_sample);
            profile.WriteGeneProfile(taxonomy, filter, &os_genes, threads_per_sample);
            os.close();
            os_total.close();
            os_dismissed.close();
            os_genes.close();
            if (os.fail() || os_total.fail() || os_dismissed.fail() || os_genes.fail()) {
                RunStatus::Get().Fail("Writing the profile of sample " + options.GetSampleId(i) + " failed: " + options.ProfileFile(i));
            }

            profile.SetName(options.GetSampleId(i));

            for (auto const& [taxid, taxon] : profile.GetTaxa()) {
                double const score = filter.Score(taxon);
                if (filter.Calls(taxon, filter.GetKnob()) || !profiler::StrongOwnEvidence(taxon)) continue;
                std::string const species = taxonomy.Get(taxid).scientific_name;
                std::ostringstream line;
                line << options.GetSampleId(i) << '\t' << species << '\t' << taxid << '\t' << score << '\t'
                     << taxon.VerticalCoverage() << '\t' << taxon.HitGeneFraction() << '\t' << taxon.TopIdentity() << '\t'
                     << taxon.LowIdentityShare() << '\t' << (filter.Calls(taxon, msa_knob) ? "yes" : "no");
                unreported_slots[idx].push_back({ species, line.str() });
            }

            // The sample's outputs are written: only the strain stage reads its profile again, and only
            // the variants and read ranges of taxa that enter the MSAs (--msa_knob).
            profile.ReleaseReadData(filter.WithKnob(msa_knob), !options.NoStrains());

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
    Profiles ProfileWrapper(ProfilingContext& ctx, ProfilingAhead* ahead = nullptr) {
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
            size_t total_obs = std::accumulate(bin.begin(), bin.end(), size_t{0}, [](size_t acc, Variant const& var) { return acc + var.Observations(); });
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

    // Whether `taxon` of `profile` enters the strain MSAs: the model `filter` scores it at the profile's MSA knob or more
    // (--msa_knob, else the threshold the sample was reported at: --knob, or its model's knob for its depth or its
    // share of false calls), and the singleton rule does not veto it.
    static bool EntersMSA(profiler::TaxonFilterObj const& filter, Profile const& profile, profiler::Taxon const& taxon) {
        return filter.Calls(taxon, profile.MSAKnob());
    }

    // Taxa in at least `min_samples` profiles (entering the MSAs by `filter`, if given), the most frequent first.
    TaxidList ExtractTaxa(Profiles const& profiles, std::optional<profiler::TaxonFilterObj> filter= {}, size_t min_samples = 2) {
        TaxidCounts taxid_counts;
        TaxidSet taxa;
        TaxidList taxid_list;

        for (auto& profile : profiles) {
            for (auto& [id, taxon] : profile.GetTaxa()) {
                if (filter.has_value() && !EntersMSA(*filter, profile, taxon)) continue;
                if (!taxid_counts.contains(id)) {
                    taxid_counts.insert({ id, 0 });
                }
                taxid_counts.at(id)++;
            }
        }

        for (auto& [t, c] : taxid_counts) {
            if (c >= min_samples) {
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

            if (filter.has_value() && !(EntersMSA(*filter, profile1, taxon1) && EntersMSA(*filter, profile2, taxon2))) {
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

    static std::vector<size_t> GetProfilesWithTaxon(uint32_t taxid, Profiles& profiles, Options& options, std::optional<profiler::TaxonFilterObj> const& filter) {
        std::vector<size_t> indices;

        for (auto i = 0; i < profiles.size(); i++) {
            auto& profile = profiles[i];
            if (!profile.GetTaxa().contains(taxid)) {
                continue;
            }
            auto& taxon_map = profile.GetTaxa();

            auto& taxon = taxon_map.at(taxid);
            if (filter.has_value() && !EntersMSA(filter.value(), profile, taxon)) {
                continue;
            }
            // std::cout << taxid << " Passes " << std::endl;

            indices.emplace_back(i);
        }
        return indices;
    }

    static std::vector<uint32_t> SelectGenesForTaxon(uint32_t taxid, std::string name, std::vector<size_t>& selected_profiles, GenomeLoader& loader, Options& options, Profiles& profiles) {
        std::vector<uint32_t> selected_gene_ids;
        // The genes with reads in any of the samples. A gene without unique k-mers has reads, too (its
        // k-mers are shared, not absent); a gene without reads would add only gaps.
        std::vector<uint32_t> gene_ids;
        {
            std::set<uint32_t> observed;
            for (auto sample_index : selected_profiles) {
                auto const& taxa = profiles[sample_index].GetTaxa();
                if (!taxa.contains(taxid)) continue;
                for (auto const& [gene_id, _] : taxa.at(taxid).GetGenes()) observed.insert(gene_id);
            }
            gene_ids.assign(observed.begin(), observed.end());
        }
        if (gene_ids.empty()) return gene_ids;

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
        // them as a second strain, so it is left out. A species without relatives has fewer, as its
        // k-mers are mostly unique at the core already: in the 2026-09-29 strain audit, Dummya solo (no
        // relative in the database) had long unique k-mers in half its genes, species with congeners
        // in 98-99%. Unless 90% of a species' genes have long unique k-mers, all its genes stay.
        auto& genome = loader.GetGenome(taxid);
        auto const with_long_uniques = std::count_if(gene_ids.begin(), gene_ids.end(), [&genome](uint32_t gene_id) {
            return genome.GetGene(gene_id).HasLongUniques();
        });
        if (static_cast<size_t>(with_long_uniques) * 10 < gene_ids.size() * 9) return gene_ids;

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
        int const status = std::system(cmd.str().c_str());
        if (status != 0) {
            // std::system returns a wait status: 512 means exit code 2.
            std::string const how = status == -1 ? "could not be started"
                                  : WIFEXITED(status) ? "exited with code " + std::to_string(WEXITSTATUS(status))
                                  : WIFSIGNALED(status) ? "was killed by signal " + std::to_string(WTERMSIG(status))
                                  : "failed (status " + std::to_string(status) + ")";
            RunStatus::Get().Fail("qcmsa " + how + " for " + name +
                                  " (post-filter skipped; the raw MSA is still in " + msa + ")");
        } else if (!fs::exists(prefix + ".msa.fna")) {
            // Not an error (qcmsa may filter everything out), but no filtered MSA is not a success either.
            std::cerr << "[qcmsa] WARNING: " << name << ": qcmsa kept no gene or sample, so there is no "
                      << prefix << ".msa.fna (see its output above; the raw MSA is " << msa << ")" << std::endl;
        }
    }

    // Removes a species' strain outputs of an earlier run, so that none of them survives a run that
    // writes it no more (e.g. an MSA qcmsa now filters away, or a species with fewer samples).
    static void RemoveStrainOutputs(Options const& options, std::string const& name) {
        std::string const prefix = options.GetStrainOutputDir() + '/' + name;
        for (auto const& stale : { options.GetMSAOutput(name), options.GetMSAPartitionOutput(name),
                                   options.GetMSAStatsOutput(name), options.GetSpeciesMetaOutput(name),
                                   options.GetHaplotypesOutput(name),
                                   prefix + ".msa.fna", prefix + ".partition.txt",
                                   prefix + ".qcmsa_summary.tsv", prefix + ".qc.png" }) {
            std::error_code ec;
            std::filesystem::remove(stale, ec);
        }
    }

    // A long-read sample's strains of taxon `taxid` (Haplotypes.h): its reads' records phased at the multi-allelic sites
    // of the genes the MSA takes (`genes`, those its taxon has and does not drop), as MSAItem gives them (the taxon's
    // own reads, by --msa_identity_margin) with the SNP filters of its reads' kind and up to four alleles a site; reads
    // are cut between genes that are unlikely neighbours in the taxon's clade (the database's gene neighbours).
    static haplotypes::Phasing PhaseSample(uint32_t taxid, profiler::MicrobialProfile const& profile, std::vector<uint32_t> const& genes,
                                           Genome& genome, GenomeLoader const& loader, Options const& options) {
        auto const& taxon = profile.GetTaxa().at(taxid);
        double const min_identity = taxon.IdentityThreshold(options.GetMSAIdentityMargin());
        double const min_af = options.GetSNPMinAF(profile.GetReadType());
        auto const min_cov = options.GetSNPMinCov();
        auto const min_qual_sum = options.GetSNPMinPhredSum();
        auto const min_mean_qual = options.GetSNPMinMeanQual();
        bool const require_strand = options.GetSNPRequireStrand();
        std::vector<haplotypes::Site> sites;
        for (auto geneid : genes) {
            if (!taxon.GetGenes().contains(geneid) || taxon.DropsGene(geneid)) continue;
            auto const item = taxon.GetGenes().at(geneid).GetStrainLevel().MSAItem(min_identity, min_cov, min_af, min_mean_qual,
                                                                                 min_qual_sum, require_strand);
            auto const reference = genome.GetGene(geneid).Sequence();
            for (auto& [pos, bases] : MultiAllelicSites(item.first, item.second, min_cov, min_qual_sum, min_af, require_strand,
                                                        min_mean_qual, 4)) {
                if (pos >= reference.size()) continue;
                sites.push_back({ geneid, pos, reference[pos], std::move(bases) });
            }
        }
        return haplotypes::Phase(taxon.PhaseRecords(), sites, MaxDivergenceBin(min_identity),
                                 haplotypes::UnlikelyNeighbours(loader.GetGeneNeighbours(), taxid));
    }

    // <species>.haplotypes.tsv: per sample phased (PhaseSample), a line per block of its multi-allelic sites (the
    // sites its reads link): its genes, sites and reads, its haplotypes' reads (the most first), the haplotype of each
    // strain row (1-based) if it is phased, and how much more likely that join is than the next; and the sample's rows
    // with their shares of the reads.
    static void WriteHaplotypes(std::ostream& os, Profiles const& profiles, std::vector<size_t> const& profile_indices,
                                std::vector<haplotypes::Phasing> const& phasings, std::vector<bool> const& phased) {
        auto join = [](auto const& values, auto&& each) {
            std::ostringstream out;
            for (size_t i = 0; i < values.size(); i++) out << (i ? "," : "") << each(values[i]);
            return values.empty() ? std::string("-") : out.str();
        };
        os << "sample\trows\trow_shares\tblock\tgenes\tsites\treads\thaplotype_reads\trow_haplotypes\tlog_odds\tphased\n";
        for (size_t i = 0; i < profile_indices.size(); i++) {
            if (!phased[i]) continue;
            auto const& p = phasings[i];
            std::string const shares = join(p.shares, [](double s) { return s; });
            for (size_t b = 0; b < p.blocks.size(); b++) {
                auto const& block = p.blocks[b];
                os << profiles[profile_indices[i]].GetName() << '\t' << p.rows << '\t' << shares << '\t' << b + 1 << '\t'
                   << join(block.genes, [](uint32_t g) { return g; }) << '\t' << block.sites << '\t' << block.reads << '\t'
                   << join(block.haplotype_reads, [](size_t n) { return n; }) << '\t'
                   << join(block.row_haplotype, [](int h) { return h + 1; }) << '\t'
                   << block.log_odds << '\t' << (block.phased ? "yes" : "no") << '\n';
            }
        }
    }

    static void GetMSAForTaxon (uint32_t taxid, std::string taxon_name, GenomeLoader& loader, Options& options, Profiles& profiles, std::ostream* os_meta=nullptr, std::optional<profiler::TaxonFilterObj> const& filter={}) {
        auto min_hcov = options.GetMSAMinHCOV();
        auto min_qual_sum = options.GetSNPMinPhredSum();
        auto min_cov = options.GetSNPMinCov();
        auto require_strand = options.GetSNPRequireStrand();
        auto min_mean_qual = options.GetSNPMinMeanQual();
        auto snp_max_alleles = options.GetSNPMaxAlleles();
        uint32_t const min_depth = static_cast<uint32_t>(options.GetMSAMinDepth());

        std::vector<size_t> profile_indices = GetProfilesWithTaxon(taxid, profiles, options, filter);

        if (profile_indices.empty()) return;

        auto& genome = loader.GetGenome(taxid);
        if (!genome.IsLoaded()) genome.LoadGenomeOMP();
        std::vector<uint32_t> selected_genes = SelectGenesForTaxon(taxid, taxon_name, profile_indices, loader, options, profiles);

        // A long-read sample whose reads show two or more strains gets a row per strain, <sample>_hap1, ... (the
        // most abundant first), each called from the reads of its strain (Haplotypes.h); every other sample one, the
        // consensus of its reads.
        std::vector<haplotypes::Phasing> phasings(profile_indices.size());
        std::vector<bool> phased(profile_indices.size(), false);  // phasing tried
        if (options.Phasing()) {
            for (size_t i = 0; i < profile_indices.size(); i++) {
                auto const& profile = profiles[profile_indices[i]];
                if (!IsLongReadType(profile.GetReadType()) || profile.GetTaxa().at(taxid).PhaseRecords().empty()) continue;
                phasings[i] = PhaseSample(taxid, profile, selected_genes, genome, loader, options);
                phased[i] = true;
            }
        }
        struct Row {
            size_t sample;  // index into profile_indices
            int haplotype;  // -1: the sample's own reads
        };
        std::vector<Row> rows;
        for (size_t i = 0; i < profile_indices.size(); i++) {
            if (phasings[i].rows < 2) rows.push_back({ i, -1 });
            else for (size_t h = 0; h < phasings[i].rows; h++) rows.push_back({ i, static_cast<int>(h) });
        }
        // A strain row's reads' records, by gene.
        std::vector<std::unordered_map<uint32_t, std::vector<haplotypes::ReadRecord const*>>> row_genes(rows.size());
        for (size_t r = 0; r < rows.size(); r++) {
            if (rows[r].haplotype < 0) continue;
            auto const& records = profiles[profile_indices[rows[r].sample]].GetTaxa().at(taxid).PhaseRecords();
            for (auto index : phasings[rows[r].sample].row_records[static_cast<size_t>(rows[r].haplotype)]) {
                row_genes[r][records[index].gene].push_back(&records[index]);
            }
        }

        // Each row takes the minimum allele frequency of its sample's reads' kind.
        std::vector<double> min_afs;
        std::vector<std::string> names;
        for (auto const& row : rows) {
            auto const& profile = profiles[profile_indices[row.sample]];
            min_afs.push_back(options.GetSNPMinAF(profile.GetReadType()));
            names.push_back(row.haplotype < 0 ? profile.GetName() : profile.GetName() + "_hap" + std::to_string(row.haplotype + 1));
        }

        MSAVector msa{ rows.size(), std::vector<char>() };
        protal::MSARow ref_msa_row;

        // Per-row accumulated SNP-retention statistics across all genes
        protal::MSAStats sample_stats(rows.size());

        std::vector<std::string> partitions;
        size_t partition_start = 0;
        size_t previous_size = 0;

        ProgressBar prog(selected_genes.size());

        std::cout << taxid << ": " << taxon_name << " across samples " << profile_indices.size() << std::endl;
        for (size_t i = 0; i < profile_indices.size(); i++) {
            auto const& p = phasings[i];
            if (!phased[i]) continue;
            std::cout << "  " << profiles[profile_indices[i]].GetName() << ": " << p.reads << " long reads at " << p.blocks.size()
                      << " blocks of multi-allelic sites (" << p.cut_reads << " cut between unlikely neighbours); ";
            if (p.rows < 2) {
                std::cout << "one row" << std::endl;
            } else {
                std::cout << p.rows << " strain rows (shares";
                for (auto share : p.shares) std::cout << ' ' << share;
                std::cout << "), " << p.PhasedBlocks() << " of the blocks phased" << std::endl;
            }
        }
        if (std::find(phased.begin(), phased.end(), true) != phased.end()) {
            std::ofstream os_haplotypes(options.GetHaplotypesOutput(taxon_name));
            WriteHaplotypes(os_haplotypes, profiles, profile_indices, phasings, phased);
            os_haplotypes.close();
            if (os_haplotypes.fail()) RunStatus::Get().Fail("Writing the haplotypes of " + taxon_name + " failed: " + options.GetHaplotypesOutput(taxon_name));
        }

        // The numbers of a row's .meta.tsv line for a gene, from its item (the gene as the MSA takes it).
        struct Meta {
            double vertical_coverage = 0;
            size_t counts_vcov1 = 0, counts_vcov2 = 0, multi_allelic = 0, filtered = 0, gene_length = 0;
            double median_vcov = 0, hcov = 0, mean_vcov_nonzero = 0, median_vcov_nonzero = 0;
        };
        auto meta_of = [&](std::pair<VariantVec, CoverageVec> const& item, double vertical_coverage, size_t filtered,
                           size_t gene_length, double min_af) {
            auto const& tmp_vec = item.second;
            Meta meta;
            meta.vertical_coverage = vertical_coverage;
            meta.counts_vcov1 = std::count_if(tmp_vec.begin(), tmp_vec.end(), [](auto val){ return(val >= 1);});
            meta.counts_vcov2 = std::count_if(tmp_vec.begin(), tmp_vec.end(), [](auto val){ return(val >= 2);});
            // Multi-allelic positions as the MSA writes them (IUPAC codes), which qcmsa filters on.
            meta.multi_allelic = MultiAllelicPositions(item.first, tmp_vec, min_cov, min_qual_sum, min_af,
                                                       require_strand, min_mean_qual, snp_max_alleles);
            meta.filtered = filtered;
            meta.gene_length = gene_length;
            meta.hcov = static_cast<double>(meta.counts_vcov1) / static_cast<double>(gene_length);
            if (os_meta) {
                auto sorted_cov = tmp_vec;
                sorted_cov.resize(gene_length, 0);
                std::sort(sorted_cov.begin(), sorted_cov.end());
                size_t n = sorted_cov.size();
                meta.median_vcov = n % 2 == 1 ?
                    static_cast<double>(sorted_cov[n/2]) :
                    (static_cast<double>(sorted_cov[n/2 - 1]) + static_cast<double>(sorted_cov[n/2])) / 2.0;

                auto nonzero_begin = std::lower_bound(sorted_cov.begin(), sorted_cov.end(), 1);
                size_t nz = std::distance(nonzero_begin, sorted_cov.end());
                if (nz > 0) {
                    double sum = std::accumulate(nonzero_begin, sorted_cov.end(), 0.0);
                    meta.mean_vcov_nonzero = sum / nz;
                    size_t mid = nz / 2;
                    meta.median_vcov_nonzero = nz % 2 == 1 ?
                        static_cast<double>(*(nonzero_begin + mid)) :
                        (static_cast<double>(*(nonzero_begin + mid - 1)) + static_cast<double>(*(nonzero_begin + mid))) / 2.0;
                }
            }
            return meta;
        };

        for (auto& geneid : selected_genes) {
            prog.UpdateAdd(1);
            auto& gene = genome.GetGene(geneid);
            auto const reference = gene.Sequence();

            // Each row's gene as the MSA takes it, and the numbers of its .meta.tsv line; none for a sample without
            // the gene, or whose taxon drops it (a foreign gene), and for a strain row whose reads lack it. A sample's
            // gene from its own reads; a strain row's from its strain's (HaplotypeItem).
            MSASequenceItems items(rows.size());
            std::vector<std::optional<Meta>> metas(rows.size());
            size_t rows_with_gene = 0;

            for (size_t r = 0; r < rows.size(); r++) {
                auto& profile = profiles[profile_indices[rows[r].sample]];
                auto& taxon = profile.GetTaxa().at(taxid);
                auto& genes = taxon.GetGenes();
                if (!genes.contains(geneid) || taxon.DropsGene(geneid)) continue;
                auto& gene_obs = genes.at(geneid);
                double const min_af = options.GetSNPMinAF(profile.GetReadType());
                if (rows[r].haplotype < 0) {
                    auto& strain = gene_obs.GetStrainLevel();
                    auto ac = gene_obs.AlleleSNPCounts(min_cov, min_qual_sum);
                    // The gene from the taxon's own reads, as the MSA takes it: the reads with a base per position,
                    // what the MSA judges each position by.
                    auto item = strain.MSAItem(taxon.IdentityThreshold(options.GetMSAIdentityMargin()),
                                               min_cov, min_af,
                                               min_mean_qual, min_qual_sum, require_strand);
                    metas[r] = meta_of(item, gene_obs.VerticalCoverage(), ac.Filtered(), gene_obs.m_gene_length, min_af);
                    items[r] = std::move(item);
                } else {
                    auto const found = row_genes[r].find(geneid);
                    if (found == row_genes[r].end()) continue;
                    auto item = haplotypes::HaplotypeItem(found->second, reference, min_cov, min_af, min_mean_qual,
                                                          min_qual_sum, require_strand);
                    size_t filtered = 0;  // sites where no allele has min_cov reads and min_qual_sum (AlleleSNPCounts)
                    for (auto const& bin : item.first) {
                        filtered += std::none_of(bin.begin(), bin.end(), [&](Variant const& v) {
                            return v.Observations() >= min_cov && v.QualitySum() >= min_qual_sum;
                        });
                    }
                    double const bases = std::accumulate(item.second.begin(), item.second.end(), 0.0);
                    metas[r] = meta_of(item, bases / static_cast<double>(gene_obs.m_gene_length), filtered,
                                       gene_obs.m_gene_length, min_af);
                    items[r] = std::move(item);
                }
                // protal emits a raw MSA: every observed gene contributes its
                // sequence. Coverage-based gating (horizontal coverage, depth,
                // min samples per gene) is done by the qcmsa post-filter, which
                // reads the hcov / mean_vcov_nonzero columns written below.
                rows_with_gene++;
            }

            if (os_meta) {
                for (size_t r = 0; r < rows.size(); r++) {
                    if (!metas[r]) continue;
                    auto const& meta = *metas[r];
                    double const v1 = static_cast<double>(meta.counts_vcov1), v2 = static_cast<double>(meta.counts_vcov2);
#pragma omp critical(metaout)
                    {
                        *os_meta << names[r] << '\t';
                        *os_meta << geneid << '\t';
                        *os_meta << meta.vertical_coverage << '\t';
                        *os_meta << meta.counts_vcov1 << '\t';
                        *os_meta << meta.counts_vcov2 << '\t';
                        *os_meta << meta.multi_allelic << '\t';
                        *os_meta << meta.filtered << '\t';
                        *os_meta << (v1 > 0 ? meta.multi_allelic / v1 : 0) << '\t';
                        *os_meta << (v1 > 0 ? meta.filtered / v1 : 0) << '\t';
                        *os_meta << (v2 > 0 ? meta.multi_allelic / v2 : 0) << '\t';
                        *os_meta << (v2 > 0 ? meta.filtered / v2 : 0) << '\t';
                        *os_meta << meta.median_vcov << '\t';
                        *os_meta << meta.hcov << '\t';
                        *os_meta << meta.gene_length << '\t';
                        *os_meta << meta.mean_vcov_nonzero << '\t';
                        *os_meta << meta.median_vcov_nonzero;
                        *os_meta << std::endl;
                    }
                }
            }
            if (rows_with_gene > 0) {
                previous_size = msa.front().size();

                protal::MSAStats gene_stats(items.size());
                bool result = protal::MSA(items, reference, msa, min_cov, min_qual_sum, min_afs, require_strand, min_mean_qual,
                                          &gene_stats, &ref_msa_row, snp_max_alleles, min_depth);

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

        // The species of this run: strain outputs of other species in the directory are an earlier
        // run's (protal leaves them alone, as the directory may be shared).
        std::ostringstream species_list;
        species_list << "species\ttaxid\tsamples\traw_msa\tfiltered_msa\n";

        for (auto& taxid : taxids) {
            std::cout << taxonomy.Get(taxid).scientific_name << std::endl;

            std::string name = taxonomy.Get(taxid).scientific_name;
            std::replace(name.begin(), name.end(), ' ', '_');
            RemoveStrainOutputs(options, name);

            size_t const samples = GetProfilesWithTaxon(taxid, profiles, options, filter).size();
            if (samples == 0) {
                std::cerr << "[strains] " << name << " (--msa_species) passes in no sample, so it has no MSA" << std::endl;
            }

            if (enable_similarity_matrix) {
                auto similarities = GetSimilarityMatrixForTaxon(taxid, options, profiles, filter);
                if (!similarities.AnySet()) continue;
                WriteDistanceMatrix(taxid, similarities, options, name);
            }

            std::ofstream os_meta(options.GetSpeciesMetaOutput(name));
            os_meta << "sample\tgene_id\tvertical_coverage\tcounts_vcov1\tcounts_vcov2\tmulti_allelic\tfiltered\tmulti_rate_vcov1\tfiltered_rate_vcov1\tmulti_rate_vcov2\tfiltered_rate_vcov2\tmedian_vcov\thcov\tgene_length\tmean_vcov_nonzero\tmedian_vcov_nonzero\n";
            GetMSAForTaxon(taxid, name, loader, options, profiles, &os_meta, filter);
            os_meta.close();
            if (os_meta.fail()) RunStatus::Get().Fail("Writing the MSA metadata of " + name + " failed: " + options.GetSpeciesMetaOutput(name));

            if (options.GetRunQCMSA()) {
                RunQCMSA(options, name);
            }
            // File names: the MSAs are next to the list.
            auto written = [](std::string const& path) {
                return std::filesystem::exists(path) ? std::filesystem::path(path).filename().string() : std::string("-");
            };
            species_list << name << '\t' << taxid << '\t' << samples << '\t'
                         << written(options.GetMSAOutput(name)) << '\t'
                         << written(options.GetStrainOutputDir() + '/' + name + ".msa.fna") << '\n';
        }

        std::ofstream os_list(options.GetStrainSpeciesListOutput(), std::ios::out);
        os_list << species_list.str();
        os_list.close();
        if (os_list.fail()) RunStatus::Get().Fail("Writing the list of strain species failed: " + options.GetStrainSpeciesListOutput());
        bm_strain.Stop();
        bm_strain.PrintResults();
    }


    size_t GetTotalSystemMemory()
    {
        long pages = sysconf(_SC_PHYS_PAGES);
        long page_size = sysconf(_SC_PAGE_SIZE);
        return pages * page_size;
    }

    // Loads the PMML presence model and its depth knobs (ParseDepthKnobCurve, or an older model's ParseDepthKnobs) and
    // checks that protal can use it (ModelContractProblemInXml); exits 2 if not.
    static profiler::TaxonFilterObj LoadModel(db::DbFile const& file, double knob) {
        std::string read_error;
        auto const xml = file.ReadAll(read_error);
        if (!xml) {
            std::cerr << "Cannot load the model " << file.Name() << ": " << read_error << std::endl;
            exit(2);
        }
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
        if (profiler::IsPlaceholderModel(*xml)) {
            std::cerr << "WARNING: " << file.Name() << " is a placeholder, not a trained model: it scores every taxon 0, so "
                      << "no species is reported (--knob 0 lists every taxon with reads). Train a model for this read type "
                      << "and store it with protal --add_model FILE --read_type TYPE --db DB." << std::endl;
        }
        return *model;
    }

    // Runs protal and returns the process exit code: 0, or 1 if any sample or output failed (see
    // RunStatus). Invalid input and fatal errors still exit directly with their own codes.
    static int Run(int argc, char *argv[]) {
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

        if (options.CompressDbMode() || options.DecompressDbMode() || options.UnpackDbMode() || !options.GetAddModel().empty()) {
            // Each exits 8 on failure.
            if (options.CompressDbMode()) protal::build::CompressDatabase(options);
            else if (options.DecompressDbMode()) protal::build::DecompressDatabase(options);
            else if (options.UnpackDbMode()) protal::build::UnpackDatabase(options);
            else {
                std::vector<std::pair<std::string, std::string>> models;  // file, member
                for (auto const& [file, read_type] : options.AddModels()) {
                    LoadModel(db::DbFile::OnDisk(file), options.GetKnob());  // exits 2 if unusable, before anything is written
                    models.emplace_back(file, Info(read_type).model_file);
                }
                protal::build::AddModel(options, models);
            }
            return RunStatus::Get().Finish();
        }

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

        // Load protal DB into RAM: reference.fna, reference.fna.zst, or the single-file database's
        // (--build always reads the folder's files).
        auto const unique_kmers_file = options.UniqueKmersDbFile();
        int const db_threads = static_cast<int>(options.GetThreads());
        ProtalDB db = unique_kmers_file.Exists() ?
            ProtalDB(options.SequenceDbFile(), options.SequenceMapDbFile(), unique_kmers_file, db_threads) :
            ProtalDB(options.SequenceDbFile(), options.SequenceMapDbFile(), db_threads);

        // Load fasta sequences of reference into RAM (advised)
        if (options.PreloadGenomes()) {
            std::cout << "Preload genomes" << std::endl;
            Benchmark bm_preload_genomes("Preload genomes");
            bm_preload_genomes.Start();
            db.GetGenomes().LoadAllGenomes(static_cast<int>(options.GetThreads()));
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
        // Loaded once, for all samples: the model of each kind of reads the samples have.
        ReadTypeModels models;
        if (run_profiling) {
            db.LoadTaxonomy(options.TaxonomyDbFile());
            msa_taxids = ResolveMSASpecies(options, db.GetTaxonomy());
            bool const any_sample = std::any_of(kReadTypes.begin(), kReadTypes.end(), [&options](ReadTypeInfo const& t) { return options.AnySample(t.type); });
            for (auto const& info : kReadTypes) {
                if (options.AnySample(info.type) || (info.type == ReadType::Paired && !any_sample)) {
                    auto const model_file = options.ModelDbFile(info.type);
                    std::cout << "Model of " << info.name << " reads: " << model_file.Name() << std::endl;
                    auto const& model = models[static_cast<size_t>(info.type)].emplace(LoadModel(model_file, options.GetKnob()));
                    if (options.FdrGiven() && options.GetFdr() > 0 && !model.HasFalseCalls()) {
                        std::cerr << "--fdr needs a model with a calibration (random_forest_cmdline.py --fdr-calls); "
                                  << model_file.Name() << " has none" << std::endl;
                        exit(2);
                    }
                    if (model.HasFalseCalls()) {
                        auto const& calls = model.FalseCalls();
                        bool const used = !options.KnobGiven() && (!options.FdrGiven() || options.GetFdr() > 0);
                        std::cout << "  calls at an expected share of false calls of " << profiler::FeatureString(calls.fdr)
                                  << " (calibrated, " << calls.curve.size() << " points; training prior "
                                  << profiler::FeatureString(calls.prior) << ")"
                                  << (!used ? (options.KnobGiven() ? "; not used, --knob is given" : "; not used, --fdr 0")
                                            : options.FdrGiven() ? "; at --fdr " + profiler::FeatureString(options.GetFdr())
                                                                 : std::string("; the depth knobs are not used"))
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
            }
            LoadGeneConservation(options, db.GetGenomes());
        }
        if (!options.BuildMode() && (run_alignment || run_profiling)) LoadGeneNeighbours(options, db);
        if (run_alignment && options.BenchmarkAlignment() && !options.GetRange().empty()) {
            // The benchmark takes each read's true gene from its name; without one it would stop
            // at the first read, deep inside the alignment.
            ThreadedGzIstream is{ options.GetFirstFile(options.GetRange().front()).c_str() };
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

        /*
         *  READ ALIGNMENT SECTION
         */
        // With --profile_ahead, a sample is profiled as soon as its SAM is complete, while the next sample's reads
        // are aligned (ProfilingAhead, on a quarter of the threads beside the alignment's); by default every
        // sample is profiled after the alignment of all, on all threads.
        std::optional<ProfilingContext> profiling;
        std::unique_ptr<ProfilingAhead> ahead;
        if (run_profiling) profiling.emplace(options, db, models);
        if (run_alignment && run_profiling && options.ProfileAhead() && !options.BenchmarkAlignment()) {
            ahead = std::make_unique<ProfilingAhead>(*profiling, std::max<size_t>(1, options.GetThreads() / 4));
        }
        auto sample_done = [&](size_t index) {
            if (!ahead) return;
            auto const& range = profiling->range;
            auto const at = std::find(range.begin(), range.end(), index);
            // Not the last sample: no alignment follows it, so the profiling stage takes it on all threads.
            if (at != range.end() && at + 1 != range.end()) ahead->Submit(static_cast<size_t>(at - range.begin()));
        };

        // Untangle Template options that need to be written out specifically.
        if (run_alignment) {
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
            auto const loaded = std::find_if(models.begin(), models.end(), [](auto const& m) { return m.has_value(); });
            auto& filter = loaded->value();
            bm_profiling.Stop();
            bm_profiling.PrintResults();

            // Per taxon: its statistics in every sample it has reads in.
            {
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
            }

            /*
             * STRAIN PART -  RESOLVE MSAs BETWEEN SAMPLES
             */
            if (!options.NoStrains()) {
                // A sample enters a taxon's MSA with a score of its profile's MSA knob or more (--msa_knob, by
                // default the sample's knob: the samples whose profile reports the taxon; EntersMSA).
                StrainWrapper2(options, profiles, db.GetGenomes(), db.GetTaxonomy(), msa_taxids, filter);
            }
        }

        bm_total.Stop();
        bm_total.PrintResults();

        std::cout << std::flush;
        return RunStatus::Get().Finish();
    }
}
