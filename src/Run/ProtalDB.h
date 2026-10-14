#pragma once

#include <future>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include "Benchmark.h"
#include "Database.h"
#include "GenomeLoader.h"
#include "Taxonomy.h"

namespace protal {
    // Hash/Seedmap.h: included only by the units that load or use the index (Load.cpp, Align.cpp, BuildStage.cpp), so
    // that the others do not compile its loading code (docs/claude/2026-10-14-compile-time).
    class Seedmap;

    // The index as a start-up thread loaded it (Run, beside the genome preload and the tables): the index, what loading it
    // printed, and its timer, for RunWrapper to print when it takes it.
    struct LoadedIndex {
        std::unique_ptr<Seedmap> map;
        std::string log;
        Benchmark timer{ "Load Index" };
    };

    class ProtalDB {
        GenomeLoader m_genomes;
        std::optional<taxonomy::IntTaxonomy> m_taxonomy;
        std::optional<std::future<LoadedIndex>> m_index;  // the index loading beside the other start-up loads (Run)

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

        // From the single-file database's binary gene table where it has a current one (gene_table, GeneTableFile.h),
        // else from reference.map and unique_kmers.tsv (if the database has one), with `threads` threads.
        ProtalDB(db::DbFile sequence_file, db::DbFile map_file, std::optional<db::DbFile> const& unique_kmers_file,
                 std::optional<db::DbFile> const& gene_table, int threads) :
                m_genomes(std::move(sequence_file), std::move(map_file), threads, gene_table,
                          unique_kmers_file ? unique_kmers_file->Size() : std::nullopt),
                m_taxonomy() {
            if (unique_kmers_file && !m_genomes.UniqueKmersFromGeneTable()) m_genomes.LoadUniqueKmers(*unique_kmers_file, threads);
        }

        // These three handle the index (a Seedmap): defined in Load.cpp, where it is a complete type.
        ~ProtalDB();

        void SetIndexLoading(std::future<LoadedIndex> index);

        // The index a start-up thread loaded, waiting for it; nullopt if none was started (RunWrapper loads it then).
        std::optional<LoadedIndex> TakeIndex();

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
}
