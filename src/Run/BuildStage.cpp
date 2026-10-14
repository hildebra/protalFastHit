// --build's index and the database maintenance modes, stages of protal::Run (RunStages.h); the only unit with Build.h.

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include "RunStages.h"
#include "Build.h"
#include "Seedmap.h"
#include "RunStatus.h"

namespace protal {

    // The index is built from --reference, but queries align reads against the database's reference.fna,
    // found through reference.map. Every --reference record must therefore be a reference.map gene of
    // the same length, with ids and length inside the index's 20-bit fields. Returns whether --reference is
    // the database's own reference.fna and holds each gene of reference.map once: then the index passes
    // take the genes from the preload instead of reading it again (build::Run).
    static bool CheckReferenceAgainstMap(Options const& options, GenomeLoader& genomes) {
        auto input = protal::build::OpenInput(options.GetSequenceFilePath());  // raw or .zst
        SeqReader reader{ input->Stream() };
        FastxRecord record;
        size_t records = 0, problems = 0;
        constexpr uint64_t max_id = (uint64_t{1} << SEEDMAP_TAXID_BITS) - 1;
        constexpr uint64_t max_gene = (uint64_t{1} << SEEDMAP_GENEID_BITS) - 1;
        constexpr uint64_t max_length = (uint64_t{1} << SEEDMAP_GENE_POS_BITS) - 1;
        // Each map gene's record seen, a bit per taxid and gene id (3 MB at r226): a gene twice makes the passes read the
        // file, as they index every record. Ids too sparse for that: the file too.
        auto const [map_taxid, map_gene, map_length] = genomes.IndexFieldMaxima();
        (void)map_length;
        uint64_t const genes_per_taxid = map_gene + 1;
        bool const track = (map_taxid + 1) * genes_per_taxid <= (uint64_t{1} << 32);
        std::vector<bool> seen(track ? (map_taxid + 1) * genes_per_taxid : 0, false);
        bool repeated = !track;
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
                } else if (track) {
                    uint64_t const bit = taxid * genes_per_taxid + geneid;
                    repeated = repeated || seen[bit];
                    seen[bit] = true;
                }
            }
            if (!problem.empty() && ++problems <= 10) {
                std::cerr << "--reference record " << record.header << ": " << problem << std::endl;
            }
        }
        if (input->Stream().bad() || !reader.Success()) {  // not the end of the file: a read or decompression error
            std::cerr << "Cannot read --reference " << options.GetSequenceFilePath() << " to its end (truncated or corrupt file?)" << std::endl;
            exit(8);
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
            return false;
        }
        std::error_code ec;
        return !repeated &&
               std::filesystem::equivalent(zstd::Resolve(options.GetSequenceFilePath()), zstd::Resolve(options.GetSequenceFile()), ec);
    }

    void RunBuild(Options& options, ProtalDB& db) {
        const size_t mmer_size = 15;
        const size_t kmer_size = 31;

        bool const reference_is_mapped = CheckReferenceAgainstMap(options, db.GetGenomes());

        // New indexes compare whole s-mers (index format 2, recorded in the index header).
        ClosedSyncmer minimizer{mmer_size, 7, 2, true};
        SimpleKmerHandler iterator{kmer_size, mmer_size, minimizer};

        Benchmark bm_build("Run build");
        bm_build.Start();
        Statistics protal_stats;
        {
            KmerPutterSM kmer_putter{};
            protal_stats = protal::build::Run<SimpleKmerHandler<ClosedSyncmer>, KmerPutterSM, DEBUG_NONE>(
                    options, kmer_putter, iterator, db.GetGenomes(), reference_is_mapped);
            // Last, as the build reads reference.fna until here: the single file, which takes the index from
            // memory and frees it once written (before compressing reference.fna), or separate files.
            if (options.WriteBundle()) protal::build::BundleDatabase(options, &kmer_putter.GetMap());
        }  // the index is freed (~36 GB at r226)
        protal::build::ReleaseFreeMemory();
        if (!options.WriteBundle()) {
            protal::build::CompressReference(options);
            protal::build::RemoveStaleBundle(options);
        }
        protal::build::PrintMemory(options.WriteBundle() ? "writing " + db::kFileName : "compressing the reference");

        bm_build.PrintResults();
        protal_stats.WriteStats(std::cout);
    }

    int RunDbMaintenance(Options& options) {
        // Each exits 8 on failure.
        if (options.CompressDbMode()) protal::build::CompressDatabase(options);
        else if (options.DecompressDbMode()) protal::build::DecompressDatabase(options);
        else if (options.UnpackDbMode()) protal::build::UnpackDatabase(options);
        else if (!options.WriteSpeciesNeighboursFile().empty()) protal::build::WriteSpeciesNeighboursOnly(options);
        else if (options.GetAddModel().empty()) protal::build::AddTables(options);
        else {
            // A model of another gate than the database's own (--model_gate; build_settings.tsv) goes in as its
            // GateModelFile, which runs gated at or near that gate take; the database's gate's as the type's model. A
            // model of a kind (--model_kind) as the KindModelFile of that, which runs given --model_kind take; it must
            // be of that kind (or a placeholder). --add_tables' tables go in with them.
            auto const settings = DatabaseBuildSettings(options);
            std::vector<std::pair<std::string, std::string>> models;  // file, member
            for (auto const& model : options.AddModels()) {
                std::string kind;
                kind = ModelKindOfFile(db::DbFile::OnDisk(model.file), options.GetKnob());  // exits 2 if unusable, before anything is written
                if (model.kind && kind != *model.kind && kind != "placeholder") {
                    std::cerr << "--add_model: " << model.file << " is " << (kind.empty() ? "no gbm or forest model" : "a " + kind + " model")
                              << ", not a " << *model.kind << " model (--model_kind)" << std::endl;
                    exit(2);
                }
                size_t const own = settings.sparse_first[static_cast<size_t>(model.type)].value_or(0);
                bool const other = model.gate && *model.gate != own;
                std::string member = other ? GateModelFile(model.type, *model.gate) : Info(model.type).model_file;
                if (model.kind) member = KindModelFile(member, *model.kind);
                if (model.gate && !other) {
                    std::cout << "--model_gate " << *model.gate << " is the database's own gate of " << ReadTypeName(model.type)
                              << " reads: " << model.file << " goes in as " << member << std::endl;
                }
                if (std::any_of(models.begin(), models.end(), [&member](auto const& m) { return m.second == member; })) {
                    std::cerr << "--add_model: two models would be stored as " << member << " (the database's gate of "
                              << ReadTypeName(model.type) << " reads is " << own << "; --model_gate)" << std::endl;
                    exit(2);
                }
                models.emplace_back(model.file, member);
            }
            if (options.AddTables().empty()) protal::build::AddModel(options, models);
            else protal::build::AddTables(options, models);
        }
        return RunStatus::Get().Finish();
    }
}
