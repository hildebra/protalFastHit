#pragma once

// What protal::Run's stages share. Each stage is a translation unit of its own in src/Run/, including only the headers
// it uses, so that an edit recompiles the stages it touches and the stages compile in parallel
// (docs/claude/2026-10-14-compile-time):
// - Run.cpp: protal::Run, the build settings a run follows, loading the models;
// - Load.cpp: the database's tables and index (the only unit besides Align.cpp and BuildStage.cpp with Seedmap.h);
// - Align.cpp: aligning the samples (RunWrapper);
// - BuildStage.cpp: --build's index and the database maintenance modes (Build.h);
// - Profile.cpp: profiling the samples (ProfileStage.h declares it, with what Run.cpp and Strains.cpp share of it);
// - Strains.cpp: the strain-level MSAs;
// - Options.cpp (src/): parsing and checking the options.

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <future>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <optional>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include "Options.h"
#include "ProtalDB.h"
#include "SamSettings.h"

namespace protal {

    // Makes the folder `dir` and its parents; whether it is there. An empty path is the current folder: without -o the
    // outputs go there (std::filesystem::create_directories throws on an empty path, which once aborted every such run).
    inline bool MakeFolder(std::filesystem::path const& dir) {
        if (dir.empty()) return true;
        std::error_code ec;
        std::filesystem::create_directories(dir, ec);
        return std::filesystem::is_directory(dir, ec);
    }

    // Whether a sample of the run is gated with the grid pick (--sparse_first_pick grid, C 2 or more): the index's values
    // then keep their grid levels (GridLevels.h), which its lookups count by.
    inline bool AnyGridSample(Options const& options) {
        auto const range = options.GetRange();
        return std::any_of(range.begin(), range.end(), [&](size_t i) {
            return options.GetSparseFirst(i) > 1 && options.GetSparseFirstPick(i) == SparseFirstPick::Grid;
        });
    }

    // Whether a sample's reads come through a named pipe (a simulation streaming its reads into protal, as the
    // database build does for large samples): they come once, so the sample is aligned whatever SAM it has, or the
    // writer would wait for a reader forever.
    inline bool ReadsArePiped(Options const& options, int index) {
        std::error_code ec;
        for (auto const& files : {options.GetFirstFiles(index), options.GetSecondFiles(index)}) {
            for (auto const& file : files) {
                if (std::filesystem::is_fifo(file, ec)) return true;
            }
        }
        return false;
    }

    // The tables only profiling reads (gene conservation, suspect copies, species priors, NCBI taxids, species neighbours,
    // congener gaps, foreign rates, ancestry sites, and the strain alleles and column weights unless the alignment scores
    // short reads with them), loading on threads of their own beside the alignment (Run, a concurrent start-up that aligns): the first
    // sample to be profiled waits for them (Wait), which prints their messages and how long the profiling waited. At GTDB
    // r226 they held the alignment back by 12 s per run, the column weights alone (docs/claude/2026-10-11-performance-audit).
    class DeferredTables {
    public:
        std::vector<std::future<void>> futures;
        std::function<std::string()> messages;  // the tables' messages, once they are loaded

        DeferredTables() = default;
        DeferredTables(DeferredTables const&) = delete;
        DeferredTables& operator=(DeferredTables const&) = delete;
        ~DeferredTables() {
            for (auto& f : futures) {
                if (f.valid()) f.wait();
            }
        }

        void Start() { m_start = std::chrono::steady_clock::now(); }

        // A table's load as a task that notes when it finished: the last of these is when the tables were loaded.
        template<typename F>
        auto Timed(F load) {
            return [this, load = std::move(load)]() mutable {
                load();
                int64_t const done = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - m_start).count();
                int64_t last = m_last_done.load();
                while (last < done && !m_last_done.compare_exchange_weak(last, done)) {}
            };
        }

        // Waits for the tables, once (any thread; the others wait meanwhile).
        void Wait() {
            std::call_once(m_once, [this]() {
                if (futures.empty()) return;
                auto const waiting = std::chrono::steady_clock::now();
                for (auto& f : futures) f.get();
                double const waited = gene_table::SecondsSince(waiting), took = static_cast<double>(m_last_done.load()) / 1e9;
                std::ostringstream text;
                if (messages) text << messages();
                text << std::fixed << std::setprecision(2) << "Loading the profiling's tables beside the alignment took " << took
                     << " s; the profiling waited " << waited << " s for them" << std::endl;
                std::cout << text.str() << std::flush;
            });
        }

    private:
        std::once_flag m_once;
        std::chrono::steady_clock::time_point m_start = std::chrono::steady_clock::now();
        std::atomic<int64_t> m_last_done{ 0 };  // ns after m_start
    };

    // Load.cpp: the database's tables (each says there what it is for; each exits 8 if its file cannot be read) and index.
    void RequireNcbiTaxids(Options const& options);
    void LoadGeneConservation(Options const& options, GenomeLoader& genomes, std::ostream& out = std::cout);
    void LoadSuspectCopies(Options const& options, GenomeLoader& genomes, std::ostream& out = std::cout);
    void LoadSpeciesPriors(Options const& options, GenomeLoader& genomes, std::ostream& out = std::cout);
    void LoadNcbiTaxids(Options const& options, GenomeLoader& genomes, std::ostream& out = std::cout);
    void LoadSpeciesNeighbours(Options const& options, GenomeLoader& genomes, std::ostream& out = std::cout);
    void LoadCongenerGaps(Options const& options, GenomeLoader& genomes, std::ostream& out = std::cout);
    void LoadForeignRates(Options const& options, GenomeLoader& genomes, std::ostream& out = std::cout);
    void LoadStrainAlleles(Options const& options, GenomeLoader& genomes, std::ostream& out = std::cout);
    void LoadColumnWeights(Options const& options, GenomeLoader& genomes, std::ostream& out = std::cout);
    void LoadAncestrySites(Options const& options, GenomeLoader& genomes, std::ostream& out = std::cout);
    void LoadUniqueCores(Options const& options, GenomeLoader& genomes, std::ostream& out = std::cout);
    std::optional<gene_neighbours::Table> ReadGeneNeighbours(Options const& options, std::ostream& out = std::cout);
    void BindGeneNeighbours(Options const& options, ProtalDB& db, gene_neighbours::Table table, std::ostream& out = std::cout);
    bool GridWhole(Options const& options);
    void CheckGridIndex(Options const& options);
    std::unique_ptr<Seedmap> LoadIndex(Options const& options, std::tuple<uint64_t, uint64_t, uint64_t> maxima,
                                       std::optional<ReferenceFingerprint> const& fingerprint, std::ostream& log);
    // The index's load on a thread of its own, beside the genome preload and the tables (ProtalDB::TakeIndex waits for it).
    void StartIndexLoading(Options const& options, ProtalDB& db);

    // Align.cpp: aligns the samples (RunWrapper; --build: the index, RunBuild). `sample_done` (if any) is called with a
    // sample's index once its SAM file is complete, or when the sample is skipped because its SAM exists.
    void RunAlignment(Options& options, ProtalDB& db, std::function<void(size_t)> const& sample_done);
    // How sample `index` is aligned, for its SAM's header (SamSettings.h); `unique_cores`: the database's unique_cores.bin
    // is read (ZS, ZU and ZT from the cores each alignment covers).
    sam_settings::Fields AlignmentSettings(Options const& options, size_t index, bool unique_cores);

    // BuildStage.cpp: --build's index, and the database maintenance modes (--compress_db, --decompress_db, --unpack_db,
    // --write_species_neighbours, --add_tables, --add_model; returns the run's exit code, each exits 8 on failure).
    void RunBuild(Options& options, ProtalDB& db);
    int RunDbMaintenance(Options& options);

    // Run.cpp: the database's build_settings.tsv (none: every setting unsaid; exits 8 if it cannot be read), and a model
    // file's kind (profiler::ModelKindOf), loaded and checked as the profiling loads it (exits 2 if unusable).
    build_settings::Settings DatabaseBuildSettings(Options const& options);
    std::string ModelKindOfFile(db::DbFile const& file, double knob);

    // Profile.cpp: every taxid's genus (profiler::GeneraOf), for the adaptive candidates of short reads.
    std::shared_ptr<std::vector<uint32_t> const> TaxonomyGenera(taxonomy::IntTaxonomy const& taxonomy);
}
