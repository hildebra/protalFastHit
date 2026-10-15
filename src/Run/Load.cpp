// The database's tables and index, as protal::Run loads them (RunStages.h).

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <future>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <vector>
#include "RunStages.h"
#include "Seedmap.h"
#include "Stage1Probe.h"
#include "SequenceUtils/GridLevels.h"
#include "MemoryUse.h"
#include "ReferenceFingerprint.h"

namespace protal {

    ProtalDB::~ProtalDB() = default;

    void ProtalDB::SetIndexLoading(std::future<LoadedIndex> index) { m_index = std::move(index); }

    std::optional<LoadedIndex> ProtalDB::TakeIndex() {
        if (!m_index) return std::nullopt;
        auto loaded = m_index->get();
        m_index.reset();
        return loaded;
    }

    // The threads a table's reading takes (DbFile::ReadAll, the binary tables' arrays): its frames inflate in parallel. The
    // start-up's tables load at once, each on a thread of its own, so each takes a few.
    static int TableThreads(Options const& options) {
        return static_cast<int>(std::clamp<size_t>(options.GetThreads(), 1, 8));
    }

    // The genes' conservation factors (Options::GeneConservationDbFile, GeneConservation.h: the database's
    // gene_conservation.tsv, or --gene_conservation's file). They give the model's conservation features
    // (conserved_fast_depth_ratio, conserved_hit_share) always, and scale the depth identity margin per gene with
    // --gene_conservation db or FILE; by default (none) every gene keeps the whole margin. Without the file the
    // features are 0 and 0.5. Exits 8 if the file cannot be read.
    void LoadGeneConservation(Options const& options, GenomeLoader& genomes, std::ostream& out) {
        auto const file = options.GeneConservationDbFile();
        bool const scale = options.ScaleDepthMarginByConservation();
        std::string const same = "the depth identity margin is the same on every gene";
        if (!file.Exists()) {
            out << "Gene conservation: the database has no " << Options::PROTAL_GENE_CONSERVATION_FILE << " (built by an "
                      << "earlier protal, or its --full_reference had no other genomes' copies of the genes): the "
                      << "conservation features are 0 and 0.5, and " << same << std::endl;
            return;
        }
        std::string error;
        auto const content = file.ReadAll(error, TableThreads(options));
        gene_conservation::Table table;
        if (content) {
            db::TextStream is(*content);
            error = table.Read(is);
        }
        if (!error.empty()) {
            std::cerr << "Invalid gene conservation factors " << file.Name() << ": " << error << std::endl;
            exit(8);
        }
        auto const [low, high] = table.Range();
        out << "Gene conservation: factors " << std::setprecision(2) << low << "-" << high << std::setprecision(6) << " for "
                  << table.Genes() << " genes (" << file.Name() << "), for the conservation features; "
                  << (scale ? "they scale the depth identity margin per gene"
                            : same + " (--gene_conservation db scales it by them)") << std::endl;
        genomes.SetGeneConservation(std::move(table));
        genomes.SetScaleDepthMargin(scale);
    }

    // The database's suspect gene copies (suspect_copies.tsv, GeneIncongruence.h): records on them are left out of the
    // evidence; with --keep_suspect_copies they count, and the table is left unread. Exits 8 if the file cannot be read.
    void LoadSuspectCopies(Options const& options, GenomeLoader& genomes, std::ostream& out) {
        if (!options.DropSuspectCopies()) {
            out << "Suspect copies: not read (--keep_suspect_copies): every gene copy counts" << std::endl;
            return;
        }
        auto const file = options.SuspectCopiesDbFile();
        if (!file.Exists()) {
            out << "Suspect copies: the database has no " << Options::PROTAL_SUSPECT_COPIES_FILE
                      << " (built by an earlier protal, or none found): every gene copy counts" << std::endl;
            return;
        }
        std::string error;
        auto const content = file.ReadAll(error, TableThreads(options));
        gene_incongruence::Table table;
        if (content) {
            db::TextStream is(*content);
            error = table.Read(is);
        }
        if (!error.empty()) {
            std::cerr << "Invalid suspect copies " << file.Name() << ": " << error << std::endl;
            exit(8);
        }
        out << "Suspect copies: " << table.Size() << " gene copies of " << table.Species() << " species near-identical to "
                  << "another genus's copy (" << file.Name() << "), their records are left out of the evidence "
                  << "(--keep_suspect_copies counts them)" << std::endl;
        genomes.SetSuspectCopies(std::move(table));
    }

    // The species' priors (species_priors.tsv, SpeciesPriors.h: duplicated markers in the representative, CheckM
    // quality, the GTDB cluster's ANI radius and width), for the model's prior features, and their genome sizes, for the
    // profiles' composition (Composition.h). Exits 8 if the file cannot be read.
    void LoadSpeciesPriors(Options const& options, GenomeLoader& genomes, std::ostream& out) {
        auto const file = options.SpeciesPriorsDbFile();
        if (!file.Exists()) {
            out << "Species priors: the database has no " << Options::PROTAL_SPECIES_PRIORS_FILE
                      << " (converted by an earlier protal): every species' priors are unknown (-1), and the profiles have "
                      << "no unknown share (no genome sizes)" << std::endl;
            return;
        }
        std::string error;
        auto const content = file.ReadAll(error, TableThreads(options));
        species_priors::Table table;
        if (content) {
            db::TextStream is(*content);
            error = table.Read(is);
        }
        if (!error.empty()) {
            std::cerr << "Invalid species priors " << file.Name() << ": " << error << std::endl;
            exit(8);
        }
        out << "Species priors: " << table.Size() << " species (" << file.Name() << "), " << table.Informative()
                  << " with a duplicated marker, CheckM quality or a species cluster known, " << table.WithGenomeSize()
                  << " with a genome size" << (table.WithGenomeSize() > 0 ? "" : " (converted before 2026-10-08: the profiles have no unknown "
                  "share; gtdb_to_protal_db.py --priors_only and protal --add_tables add the sizes)") << std::endl;
        genomes.SetSpeciesPriors(std::move(table));
    }

    // --ncbi_taxids needs the database's ncbi_taxids.tsv: exits 2 without it. Run checks this before the start-up's loaders
    // start: an exit from a loader's thread would race the others, which may still be reading the database.
    void RequireNcbiTaxids(Options const& options) {
        if (options.NcbiTaxidsDbFile().Exists()) return;
        std::cerr << "--ncbi_taxids: the database has no " << Options::PROTAL_NCBI_TAXIDS_FILE << " (converted before "
                  << "2026-10-11): write it with scripts/mini_db/gtdb_to_protal_db.py --gtdb RELEASE --priors_only --outdir DIR "
                  << "from the GTDB release the database was converted from, and store it with protal --add_tables DIR/"
                  << Options::PROTAL_NCBI_TAXIDS_FILE << " --db DB" << std::endl;
        exit(2);
    }

    // The species' NCBI taxids (ncbi_taxids.tsv, NcbiTaxids.h) for the profiles' last columns, with --ncbi_taxids. Exits 2
    // if the database has none (before any read is aligned), 8 if the file cannot be read.
    void LoadNcbiTaxids(Options const& options, GenomeLoader& genomes, std::ostream& out) {
        RequireNcbiTaxids(options);
        auto const file = options.NcbiTaxidsDbFile();
        std::string error;
        auto const content = file.ReadAll(error, TableThreads(options));
        ncbi_taxids::Table table;
        if (content) {
            db::TextStream is(*content);
            error = table.Read(is);
        }
        if (!error.empty()) {
            std::cerr << "Invalid NCBI taxids " << file.Name() << ": " << error << std::endl;
            exit(8);
        }
        out << "NCBI taxids: " << table.Size() << " species (" << file.Name() << "), " << table.Named()
            << " of them a named NCBI species, the others an NCBI placeholder" << std::endl;
        genomes.SetNcbiTaxids(std::move(table));
    }

    // Each species' nearest congeners in the database (species_neighbours.tsv, SpeciesNeighbours.h), for the
    // database-neighbourhood features and unexpected_congener_fit_share. Exits 8 if the file cannot be read.
    void LoadSpeciesNeighbours(Options const& options, GenomeLoader& genomes, std::ostream& out) {
        auto const file = options.SpeciesNeighboursDbFile();
        if (!file.Exists()) {
            out << "Species neighbours: the database has no " << Options::PROTAL_SPECIES_NEIGHBOURS_FILE
                << " (built by an earlier protal): the database-neighbourhood features are unknown (-1)" << std::endl;
            return;
        }
        std::string error;
        auto const content = file.ReadAll(error, TableThreads(options));
        species_neighbours::Table table;
        if (content) {
            db::TextStream is(*content);
            error = table.Read(is);
        }
        if (!error.empty()) {
            std::cerr << "Invalid species neighbours " << file.Name() << ": " << error << std::endl;
            exit(8);
        }
        out << "Species neighbours: " << table.Species() << " species with " << table.Pairs() << " congeners within "
            << species_neighbours::kMaxDistance << " (" << file.Name() << ")" << std::endl;
        genomes.SetSpeciesNeighbours(std::move(table));
    }

    // A per-copy table of the database (congener_gaps.tsv, foreign_rates.tsv) read into a Table: nullopt (and a note)
    // without the file; exits 8 if it cannot be read. `what` names it in the log, `missing` says what goes unknown.
    template<typename Table>
    static std::optional<Table> ReadCopyTable(db::DbFile const& file, std::string const& what, std::string const& missing,
                                              int threads, std::ostream& out) {
        if (!file.Exists()) {
            out << what << ": the database has none (" << missing << ")" << std::endl;
            return std::nullopt;
        }
        std::string error;
        auto const content = file.ReadAll(error, threads);
        Table table;
        if (content) {
            db::TextStream is(*content);
            error = table.Read(is);
        }
        if (!error.empty()) {
            std::cerr << "Invalid " << what << " " << file.Name() << ": " << error << std::endl;
            exit(8);
        }
        out << what << ": " << table.Copies() << " gene copies of " << table.Species() << " species (" << file.Name() << ")"
            << std::endl;
        return table;
    }

    // Each gene copy's gap to its congeners' copies (congener_gaps.tsv, CongenerGaps.h), for the "gaps" features.
    void LoadCongenerGaps(Options const& options, GenomeLoader& genomes, std::ostream& out) {
        if (auto table = ReadCopyTable<congener_gaps::Table>(options.CongenerGapsDbFile(), "Congener gaps",
                                                             "built by an earlier protal: the gaps features are unknown (-1)",
                                                             TableThreads(options), out)) {
            genomes.SetCongenerGaps(std::move(*table));
        }
    }

    // Each gene copy's reads of a tiled scan and their foreign share (foreign_rates.tsv, ForeignRatesTable.h), for the
    // "foreign" features.
    void LoadForeignRates(Options const& options, GenomeLoader& genomes, std::ostream& out) {
        if (auto table = ReadCopyTable<foreign_rates::Table>(options.ForeignRatesDbFile(), "Foreign rates",
                                                             "scripts/foreign_rates.py and --add_tables make it: the foreign "
                                                             "features are unknown (-1)", TableThreads(options), out)) {
            genomes.SetForeignRates(std::move(*table));
        }
    }

    // A text table read from the single-file database's binary form of it (BinaryTable.h: column_weights.bin,
    // strain_alleles.bin) where it has a current one, else from the text; exits 8 if the text cannot be read. Returns the
    // name of what it was read from; `note` says why a binary form there was not used. `what` names the table in messages.
    template<typename Table>
    static std::string ReadTextOrBinary(db::DbFile const& text, std::optional<db::DbFile> const& binary, int threads, std::string const& what,
                                        Table& table, std::string& note) {
        if (binary) {
            std::string error;
            auto const header = binary_table::ReadHeader(*binary, Table::kMagic, error);
            if (header && text.Exists() && header->text_size != text.KnownSize()) {
                error = "made from another " + text.Name() + " (of " + std::to_string(header->text_size) + " bytes): the text is read";
            } else if (header) {
                error = table.ReadBinary(*binary, *header, threads);
                if (error.empty()) return binary->Name();
            }
            note = binary->Name() + " not used: " + error;
        }
        std::string error;
        auto const content = text.ReadAll(error, threads);
        if (content) error = table.Read(std::string_view(*content));
        if (!error.empty()) {
            std::cerr << "Invalid " << what << " " << text.Name() << ": " << error << std::endl;
            exit(8);
        }
        return text.Name();
    }

    // Each gene copy's strain alleles (strain_alleles.tsv, StrainAlleles.h; from strain_alleles.bin in a single-file database
    // that has it): for the alignment scores (--no_allele_scores) and the "alleles" features. Loaded before the alignment (the
    // scores), not only for profiling.
    void LoadStrainAlleles(Options const& options, GenomeLoader& genomes, std::ostream& out) {
        auto const file = options.StrainAllelesDbFile();
        auto const binary = options.StrainAllelesBinDbFile();
        if (!file.Exists() && !binary) {
            out << "Strain alleles: the database has none (built without --full_reference, before 2026-10-08 or with "
                << "--strain_alleles 0): alignments are scored against the references alone, the alleles features unknown (-1)" << std::endl;
            return;
        }
        strain_alleles::Table table;
        std::string note;
        auto const from = ReadTextOrBinary(file, binary, TableThreads(options), "strain alleles", table, note);
        out << "Strain alleles: " << table.Alleles() << " alleles of " << table.Copies() << " gene copies of " << table.Species()
            << " species (" << from << ")" << (options.NoAlleleScores() ? "; not in the alignment scores (--no_allele_scores)" : "")
            << (note.empty() ? "" : "; " + note) << std::endl;
        genomes.SetStrainAlleles(std::move(table));
    }

    // Each family's column weights (column_weights.tsv, ColumnWeights.h): for the "weights" features and the polarised
    // ancestry sites of small genera; left unread with --no_site_weights.
    void LoadColumnWeights(Options const& options, GenomeLoader& genomes, std::ostream& out) {
        if (options.NoSiteWeights()) {
            out << "Column weights: not read (--no_site_weights): the weights features are unknown (-1)" << std::endl;
            return;
        }
        auto const file = options.ColumnWeightsDbFile();
        auto const binary = options.ColumnWeightsBinDbFile();
        if (!file.Exists() && !binary) {
            out << "Column weights: the database has none (built before 2026-10-10 or with --column_weights 0): the weights "
                << "features are unknown (-1)" << std::endl;
            return;
        }
        column_weights::Table table;
        std::string note;
        auto const from = ReadTextOrBinary(file, binary, TableThreads(options), "column weights", table, note);
        out << "Column weights: " << table.Families() << " family-gene rows of " << table.ColumnCount() << " columns, "
            << table.Copies() << " gene copies of " << table.Species() << " species mapped (" << from << "; "
            << protal::build::HumanBytes(table.Bytes()) << " in memory)" << (note.empty() ? "" : "; " + note) << std::endl;
        genomes.SetColumnWeights(std::move(table));
    }

    // Every gene copy's ancestry sites as --build computed them (ancestry_sites.bin, AncestrySites.h Table): a run decodes a
    // copy's on its first record instead of computing them there. Without the table (a database built before 2026-10-11
    // or with --no_ancestry_sites), or with --no_site_weights (the build polarised them with the weights), a run computes
    // them as it needs them; the features are the same either way.
    void LoadAncestrySites(Options const& options, GenomeLoader& genomes, std::ostream& out) {
        auto const file = options.AncestrySitesDbFile();
        if (!file.Exists()) {
            out << "Ancestry sites: the database has none (built before 2026-10-11 or with --no_ancestry_sites): computed "
                << "for each gene copy on its first record" << std::endl;
            return;
        }
        if (options.NoSiteWeights()) {
            out << "Ancestry sites: " << file.Name() << " not read (--no_site_weights: the build polarised them with the column "
                << "weights): computed for each gene copy on its first record" << std::endl;
            return;
        }
        std::string error;
        auto content = file.ReadAll(error, TableThreads(options));
        ancestry::Table table;
        if (content) error = table.Read(std::move(*content));
        if (!error.empty()) {
            std::cerr << "Invalid ancestry sites " << file.Name() << ": " << error << std::endl;
            exit(8);
        }
        out << "Ancestry sites: " << table.Copies() << " gene copies of " << table.Species() << " species with sites (" << file.Name()
            << ", version " << table.Version() << "; " << protal::build::HumanBytes(table.Bytes()) << " in memory)" << std::endl;
        genomes.SetPrebuiltSites(std::move(table));
    }

    // Which of every gene copy's index cores are unique (unique_cores.bin, UniqueCores.h), read before the alignment: with
    // them a run counts each record's ZS, ZU and ZT from the gene's cores its alignment covers with exact k-mers, the same
    // on the whole index, a sparse index or with --sparse_seeding. Without them (a database built before 2026-10-14 or
    // with --no_unique_cores, or a run with --no_unique_cores) from the record's anchor's seeds, as before. Exits 8 if the
    // file cannot be read.
    void LoadUniqueCores(Options const& options, GenomeLoader& genomes, std::ostream& out) {
        auto const file = options.UniqueCoresDbFile();
        if (!file.Exists()) {
            out << "Unique cores: the database has none (built before 2026-10-14 or with --no_unique_cores): ZS, ZU and ZT "
                << "from each record's anchor's seeds" << std::endl;
            return;
        }
        if (options.NoUniqueCores()) {
            out << "Unique cores: " << file.Name() << " not read (--no_unique_cores): ZS, ZU and ZT from each record's anchor's "
                << "seeds" << std::endl;
            return;
        }
        std::string error;
        auto content = file.ReadAll(error, TableThreads(options));
        auto table = std::make_shared<unique_cores::Table>();
        if (content) error = table->Read(std::move(*content));
        if (!error.empty()) {
            std::cerr << "Invalid unique cores " << file.Name() << ": " << error << std::endl;
            exit(8);
        }
        out << "Unique cores: " << table->Rows() << " gene copies (" << file.Name() << "; " << protal::build::HumanBytes(table->Bytes())
            << " in memory): ZS, ZU and ZT from the cores each alignment covers" << std::endl;
        genomes.SetUniqueCores(std::move(table));
    }

    // The database's gene neighbours (gene_neighbours.tsv, GeneNeighbours.h: how often each marker gene end faces
    // which other in a clade's genomes), for mate guidance past a gene's end, pairs of mates on neighbouring genes,
    // the genes next to a long read's genes and the profiler's adjacency features; none without the file or with
    // --no_gene_neighbours. Only these per-clade frequencies are loaded: the per-genome positions they were counted
    // from (gene_positions.tsv) stay in the database unread. Every species of the taxonomy gets the clades of its
    // lineage. Exits 8 if the file
    // cannot be read.
    // ReadGeneNeighbours reads and checks the table, on any thread (messages to out); BindGeneNeighbours gives every species
    // of the taxonomy its lineage's clades and hands the table to the genomes, once the taxonomy is loaded.
    std::optional<gene_neighbours::Table> ReadGeneNeighbours(Options const& options, std::ostream& out) {
        auto const file = options.GeneNeighboursDbFile();
        if (!file.Exists()) {
            out << "Gene neighbours: the database has none (" << Options::PROTAL_GENE_NEIGHBOURS_FILE << ", from whole "
                      << "genomes by scripts/mini_db/gene_neighbours.py before --build)" << std::endl;
            return std::nullopt;
        }
        if (!options.GeneNeighbours()) {
            out << "Gene neighbours: not used (--no_gene_neighbours)" << std::endl;
            return std::nullopt;
        }
        std::string error;
        auto const content = file.ReadAll(error, TableThreads(options));
        gene_neighbours::Table table;
        if (content) {
            db::TextStream is(*content);
            error = table.Read(is);
        }
        if (!error.empty()) {
            std::cerr << "Invalid gene neighbours " << file.Name() << ": " << error << std::endl;
            exit(8);
        }
        return table;
    }

    void BindGeneNeighbours(Options const& options, ProtalDB& db, gene_neighbours::Table table, std::ostream& out) {
        auto const file = options.GeneNeighboursDbFile();
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
        out << "Gene neighbours: " << table.Rules() << " rules of " << table.Clades() << " clades from " << table.Genomes()
                  << " genomes (" << file.Name() << "), for " << table.BoundSpecies() << " species: mates looked for past "
                  << "their gene's end, pairs across neighbouring genes, long reads' neighbouring genes" << std::endl;
        db.GetGenomes().SetGeneNeighbours(std::move(table));
    }

    static void LoadGeneNeighbours(Options const& options, ProtalDB& db) {
        if (auto table = ReadGeneNeighbours(options)) BindGeneNeighbours(options, db, std::move(*table));
    }

    // Loads the database's index, held packed in the widths the reference's fields need (maxima: what
    // GenomeLoader::IndexFieldMaxima gave, before any thread changes the genes), and checks it against the reference it was
    // built from: the fingerprint the gene tables gave (GenomeLoader::MapFingerprint), else reference.map read again. Its
    // messages go to log. Exits 8 if the index was built against another reference.
    // --sparse_seeding: what an index of 1/C of the cores alone (those of level log2(C) or more, stage1_probe::Level, as
    // a gate at C samples them) would hold of the loaded one's values and keys, for each C. The key map is fixed in
    // size whatever it holds (Seedmap), so only the values would shrink in today's layout.
    static void LogSparseIndex(Seedmap const& map, std::ostream& log) {
        auto const started = std::chrono::steady_clock::now();
        auto const counts = map.CountByLevel(stage1_probe::kLevels, [](uint64_t key) { return stage1_probe::Level(key); });
        Seedmap::LevelCounts all;
        for (auto const& c : counts) all.keys += c.keys, all.slots += c.slots, all.entries += c.entries;
        double const bytes_per_slot = map.BitsPerSlot() / 8;
        auto share = [](uint64_t part, uint64_t whole) { return whole ? 100.0 * static_cast<double>(part) / static_cast<double>(whole) : 0.0; };
        log << std::fixed << std::setprecision(2) << "Sparse seeding: an index of 1/C of the cores alone would hold, of the "
            << all.entries << " entries in " << all.slots << " slots (" << static_cast<double>(all.slots) * bytes_per_slot / 1e9
            << " GB) of " << all.keys << " keys:";
        Seedmap::LevelCounts part;
        std::vector<std::string> per_c;
        for (size_t l = counts.size(); l-- > 0;) {
            part.keys += counts[l].keys, part.slots += counts[l].slots, part.entries += counts[l].entries;
            std::ostringstream os;
            os << std::fixed << std::setprecision(2) << " C " << (size_t{1} << l) << ": " << part.entries << " entries ("
               << share(part.entries, all.entries) << "%), " << static_cast<double>(part.slots) * bytes_per_slot / 1e9 << " GB of values, "
               << part.keys << " keys (" << share(part.keys, all.keys) << "%)";
            per_c.push_back(os.str());
        }
        for (auto it = per_c.rbegin(); it != per_c.rend(); ++it) log << *it << (std::next(it) == per_c.rend() ? "" : ";");
        double const seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
        log << std::defaultfloat << " (counted in " << seconds << " s); the key map stays " << std::fixed << std::setprecision(1)
            << map.KeymapBytes() / 1e9 << " GB in this layout" << std::defaultfloat << std::endl;
    }

    // The grid level of the sparse index the run loads (--sparse_index; Seedmap::SetLoadLevel): with auto or on, the
    // smallest level of its samples' gates (log2 C) if every sample is gated with the grid pick at C 2 or more, else 0, the
    // whole index; with off and emulate 0. `why`: why it is 0, for on.
    static uint8_t SparseIndexLevel(Options const& options, std::string* why = nullptr) {
        if (options.SparseIndex() == "off" || options.SparseIndex() == "emulate") return 0;
        auto const range = options.GetRange();
        uint8_t level = grid_levels::kLevels;
        for (auto i : range) {
            if (options.GetSparseFirst(i) <= 1 || options.GetSparseFirstPick(i) != SparseFirstPick::Grid) {
                if (why) *why = "sample " + options.GetSampleId(i) + " is not gated with the grid pick at C 2 or more";
                return 0;
            }
            level = std::min<uint8_t>(level, static_cast<uint8_t>(grid_levels::LevelOf(options.GetSparseFirst(i))));
        }
        if (range.empty()) {
            if (why) *why = "the run has no samples";
            return 0;
        }
        return level;
    }

    // Whether the run's grid-picked samples are seeded from every entry of the whole index (the grid only gating them, two
    // seeds a gene making an anchor, as in other runs): the whole index is loaded (--sparse_index off, or a run whose
    // samples are not all grid-picked) and not --sparse_index emulate, which seeds them as the sparse index would.
    bool GridWhole(Options const& options) {
        return SparseIndexLevel(options) == 0 && options.SparseIndex() != "emulate";
    }

    // Before the index is loaded (its load runs on a thread of its own, which must not stop the run): a run with the grid
    // pick needs an index with grid levels, and --sparse_index on a run it can load sparse. Exits 2 otherwise.
    void CheckGridIndex(Options const& options) {
        std::string why;
        bool const sparse_on = options.SparseIndex() == "on";
        if (sparse_on && SparseIndexLevel(options, &why) == 0) {
            std::cerr << "Error: --sparse_index on loads the entries of the grid pick's level alone, but " << why
                      << " (--sparse_first C, --sparse_first_pick grid)" << std::endl;
            exit(2);
        }
        if (!AnyGridSample(options)) return;
        auto const features = Seedmap::IndexFeatures(options.IndexDbFile());
        if (features && !(*features & Seedmap::kFeatureGridLevels)) {
            std::cerr << "Error: --sparse_first_pick grid counts the database's grid levels, and its index has none (built "
                      << "before 2026-10-13): rebuild the database, or gate with --sparse_first_pick stride" << std::endl;
            exit(2);
        }
    }

    std::unique_ptr<Seedmap> LoadIndex(Options const& options, std::tuple<uint64_t, uint64_t, uint64_t> maxima,
                                              std::optional<ReferenceFingerprint> const& fingerprint, std::ostream& log) {
        auto map = std::make_unique<Seedmap>();
        auto const index_file = options.IndexDbFile();
        int const load_threads = static_cast<int>(options.GetThreads());
        bool const single_frame = !index_file.InBundle() && zstd::IsCompressed(index_file.Path()) && !zstd::IsSeekable(index_file.Path());
        log << "Load index " << index_file.Name() << " (" << load_threads << " thread(s)"
            << (single_frame ? "; a single zstd frame is read with one" : "") << ")" << std::endl;
        // The values are held packed, in the widths the reference's taxids, gene ids and positions
        // need (Seedmap::PackedLayout); a position is a k-mer core's, flex_k/2 bases into its gene.
        auto const [max_taxid, max_gene, max_length] = maxima;
        Seedmap::PackedLayout layout = Seedmap::PackedLayout::For(max_taxid, max_gene, max_length + map->m_flex_k);
        // The grid pick (--sparse_first_pick grid): each entry's grid level kept, and with --sparse_index (auto: a run
        // whose every sample takes it) only the entries of the smallest gate's level loaded, about 1/C of the values.
        if (AnyGridSample(options)) layout.level_bits = GRID_LEVEL_BITS;
        map->SetLoadLevel(SparseIndexLevel(options));
        // Near the end of the load the threads' frames may hold what the alignment takes beyond the index anyway (its
        // threads' buffers: 0.95 GB at r226 on 32 threads, ~30 MB a thread), so that more threads decode to the end
        // without raising the run's peak: the r226 load took 0.67 s longer with one thread's frames allowed.
        constexpr uint64_t kAlignmentBytesPerThread = uint64_t{24} << 20;
        map->SetLoadAllowance(static_cast<uint64_t>(load_threads) * kAlignmentBytesPerThread);
        // The compact key map (the default since 2026-10-14; --direct_keys: not): the load builds it, never the direct one
        // (which --sparse_seeding's statistics read first: then it is converted after them).
        if (options.CompactKeys() && !options.SparseSeeding()) map->SetCompactLoad(true);
        map->Load(index_file, load_threads, &layout);
        log << "Index features: " << map->FeatureDescription() << std::endl;
        log << "Index in memory: " << map->MemoryDescription() << std::endl;
        if (options.SparseSeeding()) LogSparseIndex(*map, log);
        // The key map of the used cores only (CompactKeys.h), the direct one freed; the same lookups. (A raw index.prx is
        // loaded with the direct one, and converted here too.)
        if (options.CompactKeys() && !map->HasCompactKeys()) {
            auto const begin = std::chrono::steady_clock::now();
            bool const compacted = map->CompactKeymap(load_threads);
            double const seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
            if (compacted) {
                log << "Index in memory (compact key map, made in " << std::fixed << std::setprecision(2) << seconds << " s): "
                    << std::defaultfloat << map->MemoryDescription() << std::endl;
            } else {
                log << "Key map not compacted: a group of its starts would span 2^24 slots or more (the direct key map kept)" << std::endl;
            }
        }
        // The compact key map as a hash table of lines (HashedKeys.h), one line a lookup; the same lookups. The default since
        // 2026-10-15 (--compact_keys: not).
        if (options.HashedKeys() && map->HasCompactKeys()) {
            auto const begin = std::chrono::steady_clock::now();
            bool const hashed = map->HashKeymap(load_threads);
            double const seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
            if (hashed) {
                log << "Index in memory (hashed key map, made in " << std::fixed << std::setprecision(2) << seconds << " s): "
                    << std::defaultfloat << map->MemoryDescription() << std::endl;
            } else {
                log << "Key map not hashed: a core would lie too far past its home line, or its start does not fit an entry (the "
                    << "compact key map kept)" << std::endl;
            }
        }
        auto const map_file = options.SequenceMapDbFile(), fna_file = options.SequenceDbFile();
        auto const reference = fingerprint ? *fingerprint : ReferenceFingerprint::Of(map_file, fna_file);
        if (map->HasReferenceFingerprint() && !(map->GetReferenceFingerprint() == reference)) {
            std::cerr << "index.prx was built against a different reference: " << map_file.Name() << " or "
                      << fna_file.Name() << " changed since the index was built. Rebuild the index (--build) or "
                      << "restore the reference files it was built with." << std::endl;
            exit(8);
        }
        return map;
    }

    void StartIndexLoading(Options const& options, ProtalDB& db) {
        // The index's field widths from the genes, before the preload changes them.
        auto const maxima = db.GetGenomes().IndexFieldMaxima();
        auto const fingerprint = db.GetGenomes().MapFingerprint();
        db.SetIndexLoading(std::async(std::launch::async, [&options, maxima, fingerprint]() {
            LoadedIndex loaded;
            loaded.timer.Start();
            std::ostringstream log;
            loaded.map = LoadIndex(options, maxima, fingerprint, log);
            loaded.timer.Stop();
            loaded.log = log.str();
            return loaded;
        }));
    }
}
