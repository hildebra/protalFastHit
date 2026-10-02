//
// Created by fritsche on 14/08/22.
//

#pragma once

#include "Options.h"
#include "SequenceUtils/SeqReader.h"
#include "SequenceUtils/KmerIterator.h"
#include "Statistics.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>
#include "Constants.h"
#include "Hash/KmerLookup.h"
#include "GenomeLoader.h"
#include "AlignmentOutputHandler.h"
#include "CoreBenchmark.h"
#include "ChainAnchorFinder.h"
#include "LongReads.h"
#include "MateGuidance.h"
#include "ScoreAlignments.h"
#include "ProgressBar.h"

namespace protal::classify {

    // The seeding and alignment diagnostics of the current sample in the misc folder: the time of
    // each stage (<sample>_runtime.tsv), and histograms of seeds and anchors per read. A stage's
    // seconds are summed over the threads that ran it; per thread is that divided by threads, as
    // --verbose prints it.
    template<typename AnchorFinder>
    static void WriteAlignmentDiagnostics(protal::Options const& options, AnchorFinder& anchor_finder_global,
                                          std::vector<Benchmark*> const& stages,
                                          Utils::Histogram& seed_sizes, Utils::Histogram& anchor_sizes) {
        auto const misc_dir = std::filesystem::path(options.GetMiscOutputDir());
        auto const sample = options.GetSampleId(options.GetCurrentIndex());

        std::ofstream time_os(misc_dir / (sample + "_runtime.tsv"), std::ios::out);
        time_os << "stage\tseconds\tthreads\tseconds_per_thread\n" << std::fixed << std::setprecision(6);
        auto row = [&time_os](Benchmark const& bm) {
            time_os << bm.GetName() << '\t' << bm.Seconds() << '\t' << bm.Threads() << '\t' << bm.MeanSeconds() << '\n';
        };
        for (Benchmark* bm : { &anchor_finder_global.m_bm_seeding, &anchor_finder_global.m_bm_processing,
                               &anchor_finder_global.m_bm_pairing, &anchor_finder_global.m_bm_sorting_anchors,
                               &anchor_finder_global.m_bm_extend_anchors }) {
            row(*bm);
        }
        for (Benchmark* bm : stages) {
            row(*bm);
        }
        time_os.close();

        seed_sizes.ToTSV((misc_dir / (sample + "_seedsizes_histogram.tsv")).string());
        anchor_sizes.ToTSV((misc_dir / (sample + "_anchorsizes_histogram.tsv")).string());
    }

    // Aligns single-end reads: each read on its own, as RunPairedEnd aligns a mate, but without a
    // mate to recover anchors from or to pair alignments with.
    template<typename KmerHandler, typename AnchorFinder, typename AlignmentHandler, typename OutputHandler, DebugLevel debug, typename AlignmentBenchmark=NoBenchmark>
    static Statistics RunSingleEnd(SeqReaderSE& reader_global, protal::Options const& options, AnchorFinder& anchor_finder_global, AlignmentHandler& alignment_handler_global, OutputHandler& output_handler_global, KmerHandler& kmer_handler_global, AlignmentBenchmark benchmark_global={}) {
        constexpr bool benchmark_active = !std::is_same<AlignmentBenchmark, NoBenchmark>();

        omp_set_num_threads(options.GetThreads());

        Statistics statistics{};
        Benchmark bm_kmer_extracter_global{"Retrieve k-mers", 0};
        Benchmark bm_anchor_finder_global{"Seed- and Anchor-finding", 0};
        Benchmark bm_anchor_recovery_global{"Anchor recovery", 0};  // paired-end only; kept for the runtime table
        Benchmark bm_alignment_global{"Alignment handler", 0};
        Benchmark bm_alignment_join_sort_global{"Joining alignment pairs and sorting", 0};  // paired-end only
        Benchmark bm_output_global{"Output handler", 0};
        Benchmark bm_reader_global{"Sequence reader"};
        Benchmark bm_omp_block{"OMP Loop handler"};

        Utils::Histogram seed_sizes_global;
        Utils::Histogram anchor_sizes_global;

        std::cout << "Start parallel execution with " << options.GetThreads() << " threads" << std::endl;
        bm_omp_block.Start();
#pragma omp parallel default(none) shared(std::cout, bm_reader_global, bm_kmer_extracter_global, bm_anchor_finder_global, bm_alignment_global, bm_output_global, seed_sizes_global, anchor_sizes_global, benchmark_global, reader_global, options, kmer_handler_global, statistics, anchor_finder_global, alignment_handler_global, output_handler_global)
        {
            FastxRecord record;

            OutputHandler output_handler(output_handler_global);
            KmerHandler kmer_handler(kmer_handler_global);
            AnchorFinder anchor_finder(anchor_finder_global);
            AlignmentHandler alignment_handler(alignment_handler_global);

            SeqReaderSE reader{ reader_global };
            Statistics thread_statistics;
            thread_statistics.thread_num = omp_get_thread_num();
            AlignmentBenchmark thread_core_benchmark{ benchmark_global };

            KmerList kmers;
            SeedList seeds;
            AlignmentAnchorList anchors;
            AlignmentResultList alignment_results;

            Benchmark bm_reader{"Sequence reader", 0, Benchmark::kPerRead};
            Benchmark bm_kmer_extracter{"Retrieve k-mers", 0, Benchmark::kPerRead};
            Benchmark bm_anchor_finder{"Seed- and Anchor-finding", 0, Benchmark::kPerRead};
            Benchmark bm_alignment{"Alignment handler", 0, Benchmark::kPerRead};
            Benchmark bm_output{"Output handler", 0, Benchmark::kPerRead};

            Utils::Histogram seed_sizes;
            Utils::Histogram anchor_sizes;

            // Seeding one read ahead, as RunPairedEnd does it for pairs.
            constexpr bool seeding_ahead = requires { anchor_finder.PrepareLookups(kmers); anchor_finder.PrefetchKeys(kmers); };
            FastxRecord next;
            KmerList next_kmers;
            auto take_kmers = [&](FastxRecord const& r, KmerList& k) {
                k.clear();
                bm_kmer_extracter.Start();
                kmer_handler(std::string_view(r.sequence), k);
                bm_kmer_extracter.Stop();
                if constexpr(KmerStatisticsConcept<KmerHandler>) {
                    thread_statistics.kmers_total += kmer_handler.TotalKmers();
                    thread_statistics.kmers_accepted += kmer_handler.TotalMinimizers();
                }
                if constexpr(seeding_ahead) anchor_finder.PrefetchKeys(k);
            };

            bm_reader.Start();
            bool more = reader(record);
            bm_reader.Stop();
            if (more) take_kmers(record, kmers);
            while (more) {
                thread_statistics.reads++;
                if constexpr(seeding_ahead) {
                    bm_anchor_finder.Start(false);
                    anchor_finder.PrepareLookups(kmers);
                    bm_anchor_finder.Stop();
                }
                bm_reader.Start();
                more = reader(next);
                bm_reader.Stop();
                if (more) take_kmers(next, next_kmers);

                seeds.clear();
                anchors.clear();
                alignment_results.clear();

                bm_anchor_finder.Start();
                anchor_finder(kmers, seeds, anchors, record.sequence);
                bm_anchor_finder.Stop();

                thread_statistics.successful_lookups += anchor_finder.m_successful_lookups;
                thread_statistics.at_least_one_lookup += (anchor_finder.m_successful_lookups > 0);
                thread_statistics.total_seeds += seeds.size();
                if (!anchor_finder.Success()) {
                    thread_statistics.errors_anchor_finding++;
                }

                seed_sizes.AddObservation(seeds.size());
                anchor_sizes.AddObservation(anchors.size());
                thread_statistics.total_anchors += anchors.size();
                if (!anchors.empty()) {
                    thread_statistics.at_least_one_anchor += 1;
                    thread_statistics.best_anchor_seed_count += anchors.front().chain.size();
                }

                bm_alignment.Start();
                alignment_handler(anchors, alignment_results, record.sequence, anchor_finder.ReverseComplement(), options.GetAlignTop(), record.id);
                bm_alignment.Stop();
                thread_statistics.total_alignments += alignment_results.size();

                bm_output.Start();
                output_handler(alignment_results, record);
                bm_output.Stop();

                if constexpr (benchmark_active) {
                    thread_core_benchmark(seeds, anchors, alignment_results, record.id);
                }

                std::swap(record, next);
                std::swap(kmers, next_kmers);
            }

#pragma omp critical(statistics)
            {
                anchor_finder_global.m_bm_seeding.Join(anchor_finder.m_bm_seeding);
                anchor_finder_global.m_bm_reverse_complement.Join(anchor_finder.m_bm_reverse_complement);
                anchor_finder_global.m_bm_operator.Join(anchor_finder.m_bm_operator);
                anchor_finder_global.m_bm_processing.Join(anchor_finder.m_bm_processing);
                anchor_finder_global.m_bm_pairing.Join(anchor_finder.m_bm_pairing);
                anchor_finder_global.m_bm_sorting_anchors.Join(anchor_finder.m_bm_sorting_anchors);
                anchor_finder_global.m_bm_extend_anchors.Join(anchor_finder.m_bm_extend_anchors);

                reader_global.UpdateSuccess(reader);

                bm_reader_global.Join(bm_reader);
                bm_kmer_extracter_global.Join(bm_kmer_extracter);
                bm_anchor_finder_global.Join(bm_anchor_finder);
                bm_alignment_global.Join(bm_alignment);
                bm_output_global.Join(bm_output);
                alignment_handler_global.bm_alignment.Join(alignment_handler.bm_alignment);
                alignment_handler_global.m_bm_alignment.Join(alignment_handler.m_bm_alignment);
                alignment_handler_global.m_anchored_alignments += alignment_handler.m_anchored_alignments;
                alignment_handler_global.m_whole_window_alignments += alignment_handler.m_whole_window_alignments;

                thread_statistics.output_alignments = output_handler.alignments;
                statistics.Join(thread_statistics);

                seed_sizes_global.Join(seed_sizes);
                anchor_sizes_global.Join(anchor_sizes);

                if constexpr (benchmark_active) {
                    benchmark_global.Join(thread_core_benchmark);
                }
            }
        }
        bm_omp_block.Stop();

        if constexpr (benchmark_active) {
            if (options.Verbose()) {
                std::cout << "\n------------Alignment benchmarks------------------" << std::endl;
                benchmark_global.WriteRowStats();
                std::cout << "----------------------------------------------------\n" << std::endl;
            }
            benchmark_global.WriteRowStatsToFile(options.GetPrefix(options.GetCurrentIndex()) + "_benchmark.tsv");
        }

        if (options.Verbose()) {
            std::cout << "---------------Speed benchmarks---------------------" << std::endl;
            bm_omp_block.PrintResults();
            bm_reader_global.PrintResults();
            bm_kmer_extracter_global.PrintResults();
            bm_anchor_finder_global.PrintResults();
            std::cout << "\t";
            anchor_finder_global.m_bm_operator.PrintResults();
            std::cout << "\t\t";
            anchor_finder_global.m_bm_seeding.PrintResults();
            std::cout << "\t\t";
            anchor_finder_global.m_bm_extend_anchors.PrintResults();
            bm_alignment_global.PrintResults();
            bm_output_global.PrintResults();
            std::cout << "----------------------------------------------------\n" << std::endl;
            std::cout << "Anchors aligned from their exact matches: " << alignment_handler_global.m_anchored_alignments
                      << ", as whole reads: " << alignment_handler_global.m_whole_window_alignments << std::endl;
        }

        WriteAlignmentDiagnostics(options, anchor_finder_global,
                                  { &bm_reader_global, &bm_kmer_extracter_global, &bm_anchor_finder_global,
                                    &bm_anchor_recovery_global, &bm_alignment_global,
                                    &bm_alignment_join_sort_global, &bm_output_global },
                                  seed_sizes_global, anchor_sizes_global);
        return statistics;
    }

    // Aligns long reads (PacBio): each read into its segments, one per gene it hits (LongReads.h).
    template<typename LongReadAligner, typename OutputHandler>
    static Statistics RunLongReads(SeqReaderSE& reader_global, protal::Options const& options, LongReadAligner& aligner_global, OutputHandler& output_handler_global) {
        omp_set_num_threads(options.GetThreads());

        Statistics statistics{};
        Benchmark bm_anchor_finder_global{"Seed- and Anchor-finding", 0};  // with the alignment, per read
        Benchmark bm_anchor_recovery_global{"Anchor recovery", 0};  // paired-end only; kept for the runtime table
        Benchmark bm_alignment_global{"Alignment handler", 0};
        Benchmark bm_alignment_join_sort_global{"Joining alignment pairs and sorting", 0};  // paired-end only
        Benchmark bm_output_global{"Output handler", 0};
        Benchmark bm_reader_global{"Sequence reader"};
        Benchmark bm_omp_block{"OMP Loop handler"};
        size_t chunked_reads = 0, ambiguous_segments = 0, settled_segments = 0, inconsistent_segments = 0, neighbour_genes = 0;

        Utils::Histogram seed_sizes_global;
        Utils::Histogram anchor_sizes_global;

        std::cout << "Start parallel execution with " << options.GetThreads() << " threads" << std::endl;
        bm_omp_block.Start();
#pragma omp parallel default(none) shared(std::cout, bm_reader_global, bm_alignment_global, bm_output_global, seed_sizes_global, anchor_sizes_global, reader_global, options, statistics, aligner_global, output_handler_global, chunked_reads, ambiguous_segments, settled_segments, inconsistent_segments, neighbour_genes)
        {
            FastxRecord record;
            LongReadSegments segments;

            OutputHandler output_handler(output_handler_global);
            LongReadAligner aligner(aligner_global);
            SeqReaderSE reader{ reader_global };
            Statistics thread_statistics;
            thread_statistics.thread_num = omp_get_thread_num();

            Benchmark bm_reader{"Sequence reader", 0, Benchmark::kPerRead};
            Benchmark bm_alignment{"Alignment handler", 0, Benchmark::kPerRead};
            Benchmark bm_output{"Output handler", 0, Benchmark::kPerRead};
            Utils::Histogram seed_sizes;
            Utils::Histogram anchor_sizes;

            bm_reader.Start();
            while (reader(record)) {
                bm_reader.Stop();
                thread_statistics.reads++;

                bm_alignment.Start();
                aligner(record, segments);
                bm_alignment.Stop();

                thread_statistics.total_seeds += aligner.LastSeeds();
                thread_statistics.total_anchors += aligner.LastAnchors();
                thread_statistics.at_least_one_anchor += aligner.LastAnchors() > 0;
                seed_sizes.AddObservation(aligner.LastSeeds());
                anchor_sizes.AddObservation(aligner.LastAnchors());
                for (auto const& segment : segments) thread_statistics.total_alignments += segment.hits.size();

                bm_output.Start();
                output_handler(segments, record);
                bm_output.Stop();

                bm_reader.Start();
            }
            bm_reader.Stop();

#pragma omp critical(statistics)
            {
                auto& anchor_finder_global = aligner_global.GetAnchorFinder();
                auto& anchor_finder = aligner.GetAnchorFinder();
                anchor_finder_global.m_bm_seeding.Join(anchor_finder.m_bm_seeding);
                anchor_finder_global.m_bm_processing.Join(anchor_finder.m_bm_processing);
                anchor_finder_global.m_bm_pairing.Join(anchor_finder.m_bm_pairing);
                anchor_finder_global.m_bm_sorting_anchors.Join(anchor_finder.m_bm_sorting_anchors);
                anchor_finder_global.m_bm_extend_anchors.Join(anchor_finder.m_bm_extend_anchors);
                chunked_reads += aligner.ChunkedReads();
                ambiguous_segments += aligner.AmbiguousSegments();
                settled_segments += aligner.SettledSegments();
                inconsistent_segments += aligner.InconsistentSegments();
                neighbour_genes += aligner.NeighbourGenes();

                reader_global.UpdateSuccess(reader);
                bm_reader_global.Join(bm_reader);
                bm_alignment_global.Join(bm_alignment);
                bm_output_global.Join(bm_output);

                thread_statistics.output_alignments = output_handler.alignments;
                statistics.Join(thread_statistics);
                seed_sizes_global.Join(seed_sizes);
                anchor_sizes_global.Join(anchor_sizes);
            }
        }
        bm_omp_block.Stop();

        if (chunked_reads > 0) {
            std::cout << chunked_reads << " read(s) longer than " << kMaxLongReadChunk << " bp were seeded in chunks overlapping by "
                      << aligner_global.ChunkOverlap() << " bp" << std::endl;
        }
        std::cout << ambiguous_segments << " gene hits fit several taxa (MAPQ < " << kConfidentMapq << "); " << settled_segments
                  << " gene hits took their read's consensus taxon (another best hit or a higher MAPQ), and "
                  << inconsistent_segments << " had no hit of it or a clearly better one of another taxon (written with MAPQ 0)" << std::endl;
        if (!aligner_global.GetGenomes().GetGeneNeighbours().Empty()) {
            std::cout << "Gene neighbours: " << neighbour_genes << " genes found on reads where the database's gene neighbours put "
                      << "them next to the reads' genes" << std::endl;
        }
        if (options.Verbose()) {
            std::cout << "---------------Speed benchmarks---------------------" << std::endl;
            bm_omp_block.PrintResults();
            bm_reader_global.PrintResults();
            bm_alignment_global.PrintResults();
            bm_output_global.PrintResults();
            std::cout << "----------------------------------------------------\n" << std::endl;
        }

        WriteAlignmentDiagnostics(options, aligner_global.GetAnchorFinder(),
                                  { &bm_reader_global, &bm_anchor_finder_global, &bm_anchor_recovery_global, &bm_alignment_global,
                                    &bm_alignment_join_sort_global, &bm_output_global },
                                  seed_sizes_global, anchor_sizes_global);
        return statistics;
    }

    static void SortAlignmentPairs1(PairedAlignmentResultList &pairs) {
        std::sort(pairs.begin(), pairs.end(), [](PairedAlignment const& a, PairedAlignment const& b) {
            return ScorePairedAlignment(a) > ScorePairedAlignment(b);
        });
    }

    static void SortAlignmentPairs2(PairedAlignmentResultList &pairs) {
        std::sort(pairs.begin(), pairs.end(), [](PairedAlignment const& a, PairedAlignment const& b) {
            return PairedAlignmentComparator(a, b);
        });
    }

    static void SortAlignmentPairs3(PairedAlignmentResultList &pairs) {
        std::sort(pairs.begin(), pairs.end(), [](PairedAlignment const& a, PairedAlignment const& b) {
            return Bitscore(a) > Bitscore(b);
        });
    }

    // Every candidate of mate 1 with every candidate of mate 2 on the same gene in the right orientation, and
    // each candidate that pairs with none on its own. With `consume`, an alignment is moved into the last pair
    // it goes into and copied into the ones before (the lists are not used after this): a read on a dense
    // database has many candidates, and each copy is two strings. The pairs are the same, in the same order.
    // With across_genes (a database with gene neighbours), also mates on two genes of one taxon that are one
    // fragment across the genes' facing ends (across_genes::PairAcrossNeighbours).
    static void JoinAlignmentPairs(PairedAlignmentResultList &pairs, AlignmentResultList &read1, AlignmentResultList &read2, GenomeLoader& loader,
                                   bool consume = false, bool across_genes = false) {
        // The pairs as indices first (SIZE_MAX: no mate), then each alignment's last use.
        static thread_local std::vector<std::pair<size_t, size_t>> joined;
        static thread_local std::vector<size_t> last1, last2;
        constexpr size_t kNone = SIZE_MAX;
        joined.clear();
        std::vector<bool> selected2(read2.size(), false);
        for (size_t i = 0; i < read1.size(); i++) {
            auto const& alignment1 = read1[i];
            bool paired = false;
            for (size_t j = 0; j < read2.size(); j++) {
                auto const& alignment2 = read2[j];
                if (alignment1.Taxid() == alignment2.Taxid() && alignment1.GeneId() == alignment2.GeneId()) {
                    if (!CorrectOrientation(alignment1, alignment2)) continue;
                    joined.emplace_back(i, j);
                    selected2[j] = true;
                    paired = true;
                } else if (across_genes && across_genes::PairAcrossNeighbours(alignment1, alignment2, loader, kRescueMaxFragment)) {
                    joined.emplace_back(i, j);
                    selected2[j] = true;
                    paired = true;
                }
            }
            if (!paired) joined.emplace_back(i, kNone);
        }
        for (size_t j = 0; j < selected2.size(); j++) {
            if (!selected2[j]) joined.emplace_back(kNone, j);
        }

        last1.assign(read1.size(), kNone);
        last2.assign(read2.size(), kNone);
        for (size_t p = 0; p < joined.size(); p++) {
            if (joined[p].first != kNone) last1[joined[p].first] = p;
            if (joined[p].second != kNone) last2[joined[p].second] = p;
        }
        auto take = [&](AlignmentResultList& list, std::vector<size_t> const& last, size_t index, size_t p) {
            if (index == kNone) return AlignmentResult();
            if (consume && last[index] == p) return AlignmentResult(std::move(list[index]));
            return AlignmentResult(list[index]);
        };
        pairs.reserve(pairs.size() + joined.size());
        for (size_t p = 0; p < joined.size(); p++) {
            pairs.emplace_back(take(read1, last1, joined[p].first, p), take(read2, last2, joined[p].second, p));
        }
    }

    template<typename KmerHandler, typename AnchorFinder, typename AlignmentHandler, typename OutputHandler, DebugLevel debug, typename AlignmentBenchmark=NoBenchmark>
//    requires KmerHandlerConcept<KmerHandler> && AnchorFinderConcept<AnchorFinder>  && AlignmentHandlerConcept<AlignmentHandler>
    static Statistics RunPairedEnd(SeqReaderPE& reader_global, protal::Options const& options, AnchorFinder& anchor_finder_global, AlignmentHandler& alignment_handler_global, OutputHandler& output_handler_global, KmerHandler& kmer_handler_global, GenomeLoader& genome_loader, AlignmentBenchmark benchmark_global={}) {
        constexpr bool benchmark_active = !std::is_same<AlignmentBenchmark, NoBenchmark>();

        constexpr bool output_data = true;

//        std::ofstream adoh_os(options.GetOutputPrefix() + ".adata", std::ios::out);
//        ProtalAlignmentDataOutputHandler adoh_global(adoh_os, 1024*1024*16, 0.8);

        size_t dummy = 0;

        // Set Thread Num
        omp_set_num_threads(options.GetThreads());

        Statistics statistics{};
        Benchmark bm_kmer_extracter_global{"Retrieve k-mers", 0};
        Benchmark bm_anchor_finder_global{"Seed- and Anchor-finding", 0};
        Benchmark bm_anchor_recovery_global{"Anchor recovery", 0};
        Benchmark bm_alignment_global{"Alignment handler", 0};
        Benchmark bm_alignment_join_sort_global{"Joining alignment pairs and sorting", 0};
        Benchmark bm_output_global{"Output handler", 0};
        MateGuidanceCounts mate_guidance_global;

        Utils::Histogram seed_sizes_global;
        Utils::Histogram anchor_sizes_global;


        Benchmark bm_omp_block{"OMP Loop handler"};
        Benchmark bm_reader_global{"Sequence reader"};
        Benchmark bm_omp_before_loop_global{ "OMP before loop" };

        std::cout << "Start parallel execution with " << options.GetThreads() << " threads" << std::endl;
        bm_omp_block.Start();
#pragma omp parallel default(none) shared(std::cout, bm_reader_global, bm_kmer_extracter_global, bm_omp_before_loop_global, bm_alignment_join_sort_global, bm_anchor_recovery_global, bm_anchor_finder_global, /*adoh_global,*/ genome_loader, seed_sizes_global, anchor_sizes_global, bm_alignment_global, bm_output_global, benchmark_global, reader_global, options, dummy, kmer_handler_global, statistics, anchor_finder_global, alignment_handler_global, output_handler_global, mate_guidance_global)
        {
            bm_omp_before_loop_global.Start();
            // Private variables
            FastxRecord record1;
            FastxRecord record2;

            // TODO this is a fix to get varkit output as quickly as possible
            OutputHandler output_handler(output_handler_global);

            // Extract variables from kmi_globalF
            KmerHandler kmer_handler(kmer_handler_global);
            AnchorFinder anchor_finder1(anchor_finder_global);
            AnchorFinder anchor_finder2(anchor_finder_global);
            AlignmentHandler alignment_handler(alignment_handler_global);

            // Alignment data output handler (set with flag)
            //            ProtalAlignmentDataOutputHandler adoh(adoh_global);

            // IO
            SeqReaderPE reader{reader_global};
            Statistics thread_statistics;
            thread_statistics.thread_num = omp_get_thread_num();
            AlignmentBenchmark thread_core_benchmark{ benchmark_global };

            // Intermediate storage objects
            KmerList kmers1;
            SeedList seeds1;
            AlignmentAnchorList anchors1;
            AlignmentResultList alignment_results1;
            KmerList kmers2;
            SeedList seeds2;
            AlignmentAnchorList anchors2;
            AlignmentResultList alignment_results2;

            PairedAlignmentResultList paired_alignment_results;
            MateGuidanceCounts mate_guidance;
            bool const across = !genome_loader.GetGeneNeighbours().Empty();  // pairs across neighbouring genes

            size_t record_id = omp_get_thread_num();

            Benchmark bm_reader{"Sequence reader", 0, Benchmark::kPerRead};
            Benchmark bm_kmer_extracter{"Retrieve k-mers", 0, Benchmark::kPerRead};
            Benchmark bm_anchor_finder{"Seed- and Anchor-finding", 0, Benchmark::kPerRead};
            Benchmark bm_anchor_recovery{"Anchor recovery", 0, Benchmark::kPerRead};
            Benchmark bm_alignment{"Alignment handler", 0, Benchmark::kPerRead};
            Benchmark bm_alignment_join_sort{"Joining alignment pairs and sorting", 0, Benchmark::kPerRead};
            Benchmark bm_output{"Output handler", 0, Benchmark::kPerRead};

            Utils::Histogram seed_sizes;
            Utils::Histogram anchor_sizes;

            // Seeding one pair ahead (ChainAnchorFinder::PrefetchKeys, PrepareLookups): the next pair is read and its
            // k-mers taken while the values of this pair's lookups are fetched, and the key map blocks of its k-mers
            // are fetched while this pair is aligned. The pairs are aligned in the order they are read, as before.
            constexpr bool seeding_ahead = requires { anchor_finder1.PrepareLookups(kmers1); anchor_finder1.PrefetchKeys(kmers1); };
            FastxRecord next1;
            FastxRecord next2;
            KmerList next_kmers1;
            KmerList next_kmers2;
            auto take_kmers = [&](FastxRecord const& r1, FastxRecord const& r2, KmerList& k1, KmerList& k2) {
                k1.clear();
                k2.clear();
                bm_kmer_extracter.Start();
                kmer_handler(std::string_view(r1.sequence), k1);
                bm_kmer_extracter.Stop();
                if constexpr(KmerStatisticsConcept<KmerHandler>) {
                    thread_statistics.kmers_total += kmer_handler.TotalKmers();
                    thread_statistics.kmers_accepted += kmer_handler.TotalMinimizers();
                }
                bm_kmer_extracter.Start();
                kmer_handler(std::string_view(r2.sequence), k2);
                bm_kmer_extracter.Stop();
                if constexpr(KmerStatisticsConcept<KmerHandler>) {
                    thread_statistics.kmers_total += kmer_handler.TotalKmers();
                    thread_statistics.kmers_accepted += kmer_handler.TotalMinimizers();
                }
                if constexpr(seeding_ahead) {
                    anchor_finder1.PrefetchKeys(k1);
                    anchor_finder2.PrefetchKeys(k2);
                }
            };

            bm_omp_before_loop_global.Stop();
            bm_reader.Start();
            bool more = reader(record1, record2);
            bm_reader.Stop();
            if (more) take_kmers(record1, record2, kmers1, kmers2);
            while (more) {
                thread_statistics.reads++;
                if constexpr(seeding_ahead) {
                    bm_anchor_finder.Start(false);
                    anchor_finder1.PrepareLookups(kmers1);
                    anchor_finder2.PrepareLookups(kmers2);
                    bm_anchor_finder.Stop();
                }
                bm_reader.Start();
                more = reader(next1, next2);
                bm_reader.Stop();
                if (more) take_kmers(next1, next2, next_kmers1, next_kmers2);

                // Clear intermediate storage objects (the k-mers are taken ahead)
                seeds1.clear();
                seeds2.clear();
                anchors1.clear();
                anchors2.clear();
                alignment_results1.clear();
                alignment_results2.clear();
                paired_alignment_results.clear();

                // Calculate Anchors


                bm_anchor_finder.Start();
                anchor_finder1(kmers1, seeds1, anchors1, record1.sequence);
                anchor_finder2(kmers2, seeds2, anchors2, record2.sequence);
                bm_anchor_finder.Stop();


                thread_statistics.successful_lookups += anchor_finder1.m_successful_lookups;
                thread_statistics.at_least_one_lookup += (anchor_finder1.m_successful_lookups > 0);
                thread_statistics.successful_lookups += anchor_finder2.m_successful_lookups;
                thread_statistics.at_least_one_lookup += (anchor_finder2.m_successful_lookups > 0);
                thread_statistics.total_seeds += seeds1.size();
                thread_statistics.total_seeds += seeds2.size();

                if (!anchor_finder1.Success() || !anchor_finder2.Success()) {
                    thread_statistics.errors_anchor_finding++;
                }


                bm_anchor_recovery.Start();
                auto recover1 = anchor_finder1.RecoverAnchors(anchors1, anchor_finder2.BestAnchors(), options.GetAlignTop());
                auto recover2 = anchor_finder2.RecoverAnchors(anchors2, anchor_finder1.BestAnchors(), options.GetAlignTop());
                bm_anchor_recovery.Stop();

                seed_sizes.AddObservation(seeds1.size());
                seed_sizes.AddObservation(seeds2.size());
                anchor_sizes.AddObservation(anchors1.size());
                anchor_sizes.AddObservation(anchors2.size());

                thread_statistics.total_anchors += anchors1.size();
                thread_statistics.total_anchors += anchors2.size();

                if (!anchors1.empty()) {
                    thread_statistics.at_least_one_anchor += 1;
                    thread_statistics.best_anchor_seed_count += anchors1.front().chain.size();
                }
                if (!anchors2.empty()) {
                    thread_statistics.at_least_one_anchor += 1;
                    thread_statistics.best_anchor_seed_count += anchors2.front().chain.size();
                }

                // Do Alignment
                bm_alignment.Start();
                alignment_handler(anchors1, alignment_results1, record1.sequence, anchor_finder1.ReverseComplement(), options.GetAlignTop() + recover1, record1.id);
                alignment_handler(anchors2, alignment_results2, record2.sequence, anchor_finder2.ReverseComplement(), options.GetAlignTop() + recover2, record2.id);
                bm_alignment.Stop();

                thread_statistics.total_alignments += alignment_results1.size();
                thread_statistics.total_alignments += alignment_results2.size();

                bm_alignment_join_sort.Start();
                // This pairs SE alignments into all possible pairs.
                // The candidate lists are not used after this, but by the alignment benchmark.
                JoinAlignmentPairs(paired_alignment_results, alignment_results1,
                                   alignment_results2, genome_loader, /*consume=*/ !benchmark_active, across);

                //                SortAlignmentPairs1(paired_alignment_results);
                //                SortAlignmentPairs2(paired_alignment_results);
                SortAlignmentPairs3(paired_alignment_results);
                // A mate sure of its alignment guides the other one to its taxon and gene (MateGuidance.h).
                if (options.MateGuidance() &&
                    GuideMate(paired_alignment_results, anchors1, anchors2, record1, record2, anchor_finder1.ReverseComplement(),
                              anchor_finder2.ReverseComplement(), alignment_handler, genome_loader, mate_guidance)) {
                    SortAlignmentPairs3(paired_alignment_results);
                }
                if (across && !paired_alignment_results.empty()) {
                    auto const& [best1, best2] = paired_alignment_results.front();
                    mate_guidance.paired_across_genes += best1.IsSet() && best2.IsSet() && best1.GeneId() != best2.GeneId();
                }
                bm_alignment_join_sort.Stop();

                // Output alignments
                bm_output.Start();
                output_handler(paired_alignment_results, record1, record2, record_id);
                bm_output.Stop();

                if constexpr (benchmark_active) {
                    thread_core_benchmark(seeds1, anchors1, alignment_results1, record1.id);
                    thread_core_benchmark(seeds2, anchors2, alignment_results2, paired_alignment_results, record1.id);
                    //#pragma omp critical(write)
                    //                    thread_core_benchmark.ErrorOutput(seeds1, seeds2, anchors1, anchors2, alignment_results1, alignment_results2, paired_alignment_results, record1, record2);
                }
                // Check if this is the omp pragma that slows everything down TODO:
//                if constexpr(debug == DEBUG_VERBOSE) {
//#pragma omp critical(write)
//                    {
//                        thread_statistics.WriteStats(std::cout);
//                    }
//                }
                if constexpr(debug == DEBUG_EXTRAVERBOSE) {

                }
                record_id += options.GetThreads();

                std::swap(record1, next1);
                std::swap(record2, next2);
                std::swap(kmers1, next_kmers1);
                std::swap(kmers2, next_kmers2);
            }

#pragma omp critical(statistics)
            {
                anchor_finder_global.m_bm_seeding.Join(anchor_finder1.m_bm_seeding);
                anchor_finder_global.m_bm_seeding.Join(anchor_finder2.m_bm_seeding, false);
                anchor_finder_global.m_bm_reverse_complement.Join(anchor_finder1.m_bm_reverse_complement);
                anchor_finder_global.m_bm_reverse_complement.Join(anchor_finder2.m_bm_reverse_complement, false);
                anchor_finder_global.m_bm_operator.Join(anchor_finder1.m_bm_operator);
                anchor_finder_global.m_bm_operator.Join(anchor_finder2.m_bm_operator, false);
                anchor_finder_global.m_bm_seed_subsetting.Join(anchor_finder1.m_bm_seed_subsetting);
                anchor_finder_global.m_bm_seed_subsetting.Join(anchor_finder2.m_bm_seed_subsetting, false);
                anchor_finder_global.m_bm_processing.Join(anchor_finder1.m_bm_processing);
                anchor_finder_global.m_bm_processing.Join(anchor_finder2.m_bm_processing, false);
                anchor_finder_global.m_bm_pairing.Join(anchor_finder1.m_bm_pairing);
                anchor_finder_global.m_bm_pairing.Join(anchor_finder2.m_bm_pairing, false);
                anchor_finder_global.m_bm_sorting_anchors.Join(anchor_finder1.m_bm_sorting_anchors);
                anchor_finder_global.m_bm_sorting_anchors.Join(anchor_finder2.m_bm_sorting_anchors, false);
                anchor_finder_global.m_bm_recovering_anchors.Join(anchor_finder1.m_bm_recovering_anchors);
                anchor_finder_global.m_bm_recovering_anchors.Join(anchor_finder2.m_bm_recovering_anchors, false);
                anchor_finder_global.m_bm_extend_anchors.Join(anchor_finder1.m_bm_extend_anchors);
                anchor_finder_global.m_bm_extend_anchors.Join(anchor_finder2.m_bm_extend_anchors, false);
                anchor_finder_global.recovered_count += anchor_finder1.recovered_count;
                anchor_finder_global.recovered_count += anchor_finder2.recovered_count;
                anchor_finder_global.total_count += anchor_finder1.total_count;
                anchor_finder_global.total_count += anchor_finder2.total_count;

                reader_global.UpdateSuccess(reader);

                bm_reader_global.Join(bm_reader);

                bm_kmer_extracter_global.Join(bm_kmer_extracter);
                bm_anchor_finder_global.Join(bm_anchor_finder);
                bm_anchor_recovery_global.Join(bm_anchor_recovery);
                bm_alignment_global.Join(bm_alignment);
                bm_alignment_join_sort_global.Join(bm_alignment_join_sort);
                bm_output_global.Join(bm_output);
                mate_guidance_global.Join(mate_guidance);
                alignment_handler_global.bm_alignment.Join(alignment_handler.bm_alignment);

                alignment_handler_global.m_bm_alignment.Join(alignment_handler.m_bm_alignment);
                alignment_handler_global.m_anchored_alignments += alignment_handler.m_anchored_alignments;
                alignment_handler_global.m_whole_window_alignments += alignment_handler.m_whole_window_alignments;

                thread_statistics.output_alignments = output_handler.alignments;
                statistics.Join(thread_statistics);

                seed_sizes_global.Join(seed_sizes);
                anchor_sizes_global.Join(anchor_sizes);

                if constexpr (benchmark_active) {
                    benchmark_global.Join(thread_core_benchmark);
                }
            }
        }
        bm_omp_block.Stop();
        if (options.MateGuidance()) {
            std::cout << mate_guidance_global.guided << " fragments had one mate sure of its alignment and the other "
                      << "without a candidate of its taxon: " << mate_guidance_global.from_anchor << " aligned from its own "
                      << "anchor of that taxon, " << mate_guidance_global.rescued << " of " << mate_guidance_global.looked_for
                      << " found on the mate's gene" << std::endl;
        }
        if (!genome_loader.GetGeneNeighbours().Empty()) {
            std::cout << "Gene neighbours: " << mate_guidance_global.paired_across_genes << " fragments paired across two neighbouring "
                      << "genes; " << mate_guidance_global.rescued_on_neighbour << " guided mates found on the gene next to their "
                      << "guiding mate's" << std::endl;
        }

//        adoh_os.close();

        if constexpr (benchmark_active) {
            if (options.Verbose()) {
                std::cout << "\n------------Alignment benchmarks------------------" << std::endl;
                benchmark_global.WriteRowStats();
                std::cout << "----------------------------------------------------\n" << std::endl;
            }
            benchmark_global.WriteRowStatsToFile(options.GetPrefix(options.GetCurrentIndex()) + "_benchmark.tsv");
        }



        if (options.Verbose()) {
            std::cout << "---------------Speed benchmarks---------------------" << std::endl;
            bm_omp_block.PrintResults();
            bm_omp_before_loop_global.PrintResults();
            std::cout << "-----" << std::endl;
            bm_reader_global.PrintResults();
            bm_kmer_extracter_global.PrintResults();
            bm_anchor_finder_global.PrintResults();
            std::cout << "\t";
            anchor_finder_global.m_bm_operator.PrintResults();
            std::cout << "\t\t";
            anchor_finder_global.m_bm_reverse_complement.PrintResults();
            std::cout << "\t\t";
            anchor_finder_global.m_bm_seeding.PrintResults();
            // anchor_finder_global.m_bm_seed_subsetting.PrintResults();
            std::cout << "\t\t";
            anchor_finder_global.m_bm_processing.PrintResults();
            std::cout << "\t\t";
            anchor_finder_global.m_bm_pairing.PrintResults();
            std::cout << "\t\t";
            anchor_finder_global.m_bm_extend_anchors.PrintResults();
            std::cout << "\t\t";
            anchor_finder_global.m_bm_sorting_anchors.PrintResults();
            std::cout << "\t\t";
            anchor_finder_global.m_bm_recovering_anchors.PrintResults();
            bm_anchor_recovery_global.PrintResults();
            bm_alignment_global.PrintResults();
            bm_alignment_join_sort_global.PrintResults();
            bm_output_global.PrintResults();
            alignment_handler_global.bm_alignment.PrintResults();
            alignment_handler_global.m_bm_alignment.PrintResults();
            std::cout << "----------------------------------------------------\n" << std::endl;

            std::cout << "Reads that had anchors recovered: "
                      << static_cast<double>(anchor_finder_global.recovered_count) / anchor_finder_global.total_count
                      << std::endl;
            std::cout << "Anchors aligned from their exact matches: " << alignment_handler_global.m_anchored_alignments
                      << ", as whole reads: " << alignment_handler_global.m_whole_window_alignments << std::endl;
        }


        WriteAlignmentDiagnostics(options, anchor_finder_global,
                                  { &bm_reader_global, &bm_kmer_extracter_global, &bm_anchor_finder_global,
                                    &bm_anchor_recovery_global, &bm_alignment_global,
                                    &bm_alignment_join_sort_global, &bm_output_global },
                                  seed_sizes_global, anchor_sizes_global);
        return statistics;
    }

}