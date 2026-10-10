//
// Created by fritsche on 14/08/22.
//

#pragma once

#include "Options.h"
#include "SequenceUtils/SeqReader.h"
#include "SequenceUtils/FastaBatches.h"
#include "SequenceUtils/KmerIterator.h"
#include "Statistics.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>
#include <omp.h>
#include "Constants.h"
#include "Hash/KmerPutter.h"
#include "robin_map.h"
#include "KmerUtils.h"
#include "Benchmark.h"
#include "Zstd.h"
#include <filesystem>
#include "SequenceUtils/GenomeLoader.h"
#include "Profiling/SampleContext.h"
#include "SequenceUtils/CongenerGaps.h"
#include "SequenceUtils/StrainAllelesBuild.h"
#include "SequenceUtils/ColumnWeightsBuild.h"
#include "SequenceUtils/ColumnWeightsMsa.h"
#include "Taxonomy.h"
#ifdef __GLIBC__
#include <malloc.h>
#endif

namespace protal::build {
    // The build's resident memory now and its peak so far (VmRSS and VmHWM in /proc/self/status), after
    // a phase: which phase sets a GTDB build's peak (docs/claude/2026-10-05-build-memory). Nothing
    // where /proc is not there.
    inline void PrintMemory(std::string const& after) {
        std::ifstream status("/proc/self/status");
        std::string line;
        double rss = -1, hwm = -1;
        while (std::getline(status, line)) {
            if (line.rfind("VmRSS:", 0) == 0) rss = std::stod(line.substr(6));
            else if (line.rfind("VmHWM:", 0) == 0) hwm = std::stod(line.substr(6));
        }
        if (rss < 0 || hwm < 0) return;
        char text[96];
        std::snprintf(text, sizeof(text), "%.2f GB resident, peak %.2f GB", rss / (1024.0 * 1024.0), hwm / (1024.0 * 1024.0));
        std::cout << "Memory after " << after << ": " << text << std::endl;
    }

    // Gives the memory freed so far back to the system. glibc keeps freed small blocks in its
    // per-thread arenas: the suspect-copy scan's millions of sketches left ~12 GB resident at r226
    // scale (64 threads) that no later phase reused, under the index.
    inline void ReleaseFreeMemory() {
#ifdef __GLIBC__
        malloc_trim(0);
#endif
    }

    // The k-mer an index entry was built from, read back from its gene and encoded as the build's
    // k-mer handler encodes it: the k-mer or its reverse complement, whichever has the smaller
    // core. Entries store the position of their core, which starts flex_k/2 bases into the k-mer.
    inline bool IndexedKmer(Seedmap& map, GenomeLoader& genomes, uint64_t taxid, uint64_t geneid, uint64_t genepos, uint64_t& key) {
        size_t const k = map.m_exact_k + map.m_flex_k;
        if (genepos < map.m_flex_k_half) return false;
        auto const sequence = genomes.GetGenome(taxid).GetGeneOMP(geneid).Sequence();
        size_t const start = genepos - map.m_flex_k_half;
        if (start + k > sequence.size()) return false;
        uint64_t fwd = 0, rev = 0;
        for (size_t i = 0; i < k; i++) {
            fwd |= KmerUtils::BaseToInt(sequence[start + i], 0) << (2 * (k - 1 - i));
            rev |= KmerUtils::BaseToIntC(sequence[start + i], 0) << (2 * i);
        }
        key = map.MainKey(fwd) < map.MainKey(rev) ? fwd : rev;
        return true;
    }

    // Ambiguous bases (anything but A, C, G, T: N, the IUPAC codes) do not go into the index: a k-mer
    // whose window holds one is taken out of a record's k-mers (their second is the window's first
    // position), in the passes that count and place the values and in the uniqueness check. A sequence of
    // A, C, G, T only, nearly always, is scanned once and left alone.
    inline void DropAmbiguousKmers(std::string_view sequence, size_t k, KmerList& kmers) {
        if (kmers.empty()) return;
        unsigned ambiguous = 0;
        for (char const c : sequence) {
            unsigned char const upper = static_cast<unsigned char>(c) & 0xDF;  // lowercase is read as uppercase
            ambiguous |= !(upper == 'A' || upper == 'C' || upper == 'G' || upper == 'T');
        }
        if (!ambiguous) return;
        std::vector<uint32_t> before(sequence.size() + 1, 0);  // ambiguous bases in sequence[0, i)
        for (size_t i = 0; i < sequence.size(); i++) {
            unsigned char const upper = static_cast<unsigned char>(sequence[i]) & 0xDF;
            before[i + 1] = before[i] + !(upper == 'A' || upper == 'C' || upper == 'G' || upper == 'T');
        }
        std::erase_if(kmers, [&](KmerElement const& kmer) {
            size_t const begin = std::min(kmer.second, sequence.size());
            size_t const end = std::min(kmer.second + k, sequence.size());
            return before[end] != before[begin];
        });
    }

    // Opens a build input (--reference, --full_reference); a .zst file or sibling is fine too.
    inline std::unique_ptr<zstd::InputFile> OpenInput(std::string const& path) {
        auto input = std::make_unique<zstd::InputFile>(zstd::Resolve(path));
        if (!input->IsOpen()) {
            std::cerr << "Cannot open " << path << " (nor " << path << zstd::kExtension << ")" << std::endl;
            exit(8);
        }
        return input;
    }

    inline std::string HumanBytes(uint64_t bytes) {
        char buffer[32];
        if (bytes >= (uint64_t{1} << 30)) std::snprintf(buffer, sizeof(buffer), "%.2f GB", bytes / double(uint64_t{1} << 30));
        else if (bytes >= (uint64_t{1} << 20)) std::snprintf(buffer, sizeof(buffer), "%.2f MB", bytes / double(uint64_t{1} << 20));
        else std::snprintf(buffer, sizeof(buffer), "%llu bytes", static_cast<unsigned long long>(bytes));
        return buffer;
    }

    inline std::string FrameDescription(zstd::Params const& params) {
        return params.frame_size > 0 ? "seekable, " + HumanBytes(params.frame_size) + " frames" : "one frame";
    }

    inline void CompressionHint(protal::Options const& options) {
        if (options.CompressionParams().threads <= 1 && options.CompressionParams().level >= 16) {
            std::cout << "Note: zstd level " << options.CompressionParams().level << " with one thread compresses about "
                      << "3 MB/s; -t <threads> speeds it up almost linearly (or lower --compress_level)." << std::endl;
        }
    }

    // Writes the index as index.prx.zst (index.prx with --no_compress) via a .partial file, then
    // removes the other variant so that a stale index cannot be picked up later.
    template<typename KmerPutter>
    static void SaveIndex(protal::Options const& options, KmerPutter& putter) {
        std::string const raw = options.GetIndexFile();
        bool const compress = options.Compress();
        std::string const target = compress ? raw + zstd::kExtension : raw;
        std::string const stale = compress ? raw : raw + zstd::kExtension;
        std::string const partial = target + ".partial";
        uint64_t const size = putter.GetMap().SerializedSize();

        Benchmark bm_save("Write index");
        bm_save.Start();
        bool ok = false;
        uint64_t written = size;
        if (compress) {
            auto const params = options.CompressionParams();
            CompressionHint(options);
            std::cout << "Write " << target << " (zstd level " << params.level << ", "
                      << (params.frame_size > 0 ? "columns in " + HumanBytes(params.frame_size) + " chunks, verified" : "one frame")
                      << ", " << params.threads << " thread(s))" << std::endl;
            if (params.frame_size > 0) {
                // Column format (IndexCodec.h), chunks compressed in parallel, read back and compared.
                size_t raw_chunks = 0;
                std::string const error = putter.GetMap().SaveCompressed(partial, params, written, raw_chunks);
                ok = error.empty();
                if (!ok) std::cerr << "Error writing " << partial << ": " << error << std::endl;
                else if (raw_chunks > 0) std::cout << raw_chunks << " index chunk(s) kept as raw cells" << std::endl;
            } else {
                zstd::OStream os(partial, params, size);
                putter.Save(os);
                ok = os.Close();
                written = os.Buffer().BytesOut();
            }
        } else {
            std::cout << "Write " << target << std::endl;
            std::ofstream os(partial, std::ios::binary);
            putter.Save(os);
            os.close();
            ok = !os.fail();
        }
        std::error_code ec;
        if (ok) std::filesystem::rename(partial, target, ec);
        if (!ok || ec) {
            std::cerr << "Writing the index " << target << " failed" << (ec ? ": " + ec.message() : "") << std::endl;
            std::filesystem::remove(partial, ec);
            exit(8);
        }
        if (std::filesystem::remove(stale, ec)) std::cout << "Removed the previous " << stale << std::endl;
        bm_save.Stop();
        std::cout << "Index " << target << ": " << HumanBytes(written);
        if (compress) std::cout << " (" << HumanBytes(size) << " uncompressed, " << std::fixed << std::setprecision(1)
                                << size / double(std::max<uint64_t>(written, 1)) << "x)" << std::defaultfloat;
        std::cout << std::endl;
        bm_save.PrintResults();
    }

    // Unless --no_compress: replaces the database's reference.fna by reference.fna.zst, which is
    // decompressed and compared with reference.fna before reference.fna is removed. With
    // --no_compress, a reference.fna.zst next to reference.fna is stale and removed.
    static void CompressReference(protal::Options const& options) {
        std::string const raw = options.GetSequenceFile();
        std::string const zst = raw + zstd::kExtension;
        std::error_code ec;
        if (!std::filesystem::exists(raw)) return;  // only reference.fna.zst: nothing to do
        if (!options.Compress()) {
            if (std::filesystem::remove(zst, ec)) std::cout << "Removed the previous " << zst << std::endl;
            return;
        }
        Benchmark bm("Compress reference");
        bm.Start();
        CompressionHint(options);
        std::cout << "Write " << zst << " (zstd level " << options.CompressionParams().level << ", "
                  << FrameDescription(options.CompressionParams()) << ", verified)" << std::endl;
        std::string error;
        if (!zstd::CompressFile(raw, zst, options.CompressionParams(), true, error)) {
            std::cerr << "Compressing the reference failed: " << error << std::endl;
            exit(8);
        }
        uint64_t const before = std::filesystem::file_size(raw, ec);
        uint64_t const after = std::filesystem::file_size(zst, ec);
        std::filesystem::remove(raw, ec);
        if (ec) std::cerr << "Warning: cannot remove " << raw << ": " << ec.message() << std::endl;
        bm.Stop();
        std::cout << "Reference " << zst << ": " << HumanBytes(after) << " (" << HumanBytes(before) << " uncompressed, "
                  << std::fixed << std::setprecision(1) << before / double(std::max<uint64_t>(after, 1)) << "x); "
                  << "removed " << raw << std::defaultfloat << std::endl;
        bm.PrintResults();
    }

    // Replaces source by target (written as target.partial): renames, removes source if it differs.
    static void ReplaceFile(std::string const& partial, std::string const& target, std::string const& source) {
        std::error_code ec;
        std::filesystem::rename(partial, target, ec);
        if (ec) {
            std::cerr << "Cannot rename " << partial << " to " << target << ": " << ec.message() << std::endl;
            exit(8);
        }
        if (source != target) {
            std::filesystem::remove(source, ec);
            if (ec) std::cerr << "Warning: cannot remove " << source << ": " << ec.message() << std::endl;
        }
    }

    // --compress_db for the index: loaded in any form, written in the column format (IndexCodec.h),
    // read back and compared. Needs the index in memory.
    static void CompressIndexFile(protal::Options const& options, zstd::Params const& params) {
        std::string const raw = options.GetIndexFile();
        std::string const source = zstd::Resolve(raw);
        std::string const target = raw + zstd::kExtension;
        if (!std::filesystem::exists(source)) {
            std::cerr << "Cannot compress " << raw << ": neither it nor " << target << " exists" << std::endl;
            exit(8);
        }
        if (source == target && index_codec::IsSplitIndex(source)) {
            std::cout << target << " is already in the column format; kept" << std::endl;
            return;
        }
        Benchmark bm("Compress index.prx");
        bm.Start();
        std::cout << "Write " << target << " from " << source << " (zstd level " << params.level << ", columns in "
                  << HumanBytes(params.frame_size) << " chunks, " << params.threads << " thread(s), verified)" << std::endl;
        std::error_code ec;
        uint64_t const before = std::filesystem::file_size(source, ec);
        Seedmap map;
        map.Load(source, params.threads);
        uint64_t written = 0;
        size_t raw_chunks = 0;
        std::string const partial = target + ".partial";
        std::string const error = map.SaveCompressed(partial, params, written, raw_chunks);
        if (!error.empty()) {
            std::cerr << "Compressing " << source << " failed: " << error << std::endl;
            std::filesystem::remove(partial, ec);
            exit(8);
        }
        ReplaceFile(partial, target, source);
        bm.Stop();
        std::cout << target << ": " << HumanBytes(written) << " (was " << HumanBytes(before) << " as "
                  << std::filesystem::path(source).filename().string() << ")"
                  << (raw_chunks ? "; " + std::to_string(raw_chunks) + " chunk(s) kept as raw cells" : "") << std::endl;
        bm.PrintResults();
    }

    // The files of a database folder that go into database.protal (Database.h), in member order: the
    // index (index.prx.zst in the column format, frames copied as they are), the reference (a
    // seekable reference.fna.zst is copied the same way, reference.fna compressed), and the other
    // files queries read, compressed: reference.map, internal_taxonomy.dmp, unique_kmers.tsv,
    // gene_conservation.tsv, suspect_copies.tsv, species_neighbours.tsv and congener_gaps.tsv if the build wrote them, species_priors.tsv,
    // foreign_rates.tsv (scripts/foreign_rates.py), gene_neighbours.tsv and gene_positions.tsv if the folder has them
    // (written by scripts/mini_db/gene_neighbours.py, checked by CheckGeneNeighbours and CheckGenePositions; a run
    // reads only the first), and every presence model there is (AllModelFiles in
    // ReadType.h: model_pe.xml, model_se.xml, model_PB.xml, model_ONT.xml, and model.xml /
    // random_forest.xml of older databases).
    static std::vector<db::Source> BundleSources(protal::Options const& options) {
        namespace fs = std::filesystem;
        std::vector<db::Source> sources = {
                {Options::PROTAL_INDEX_FILE, options.ResolvedIndexFile()},
                {Options::PROTAL_SEQUENCE_FILE, options.ResolvedSequenceFile()},
                {Options::PROTAL_SEQUENCE_MAP_FILE, options.GetSequenceMapFile()},
                {Options::PROTAL_TAXONOMY_FILE, options.GetInternalTaxonomyFile()}};
        if (fs::exists(options.GetUniqueKmersFile())) sources.push_back({Options::PROTAL_UNIQUE_KMER_FILE, options.GetUniqueKmersFile()});
        if (fs::exists(options.GetGeneConservationFile())) sources.push_back({Options::PROTAL_GENE_CONSERVATION_FILE, options.GetGeneConservationFile()});
        if (fs::exists(options.GetSuspectCopiesFile())) sources.push_back({Options::PROTAL_SUSPECT_COPIES_FILE, options.GetSuspectCopiesFile()});
        if (fs::exists(options.GetSpeciesPriorsFile())) sources.push_back({Options::PROTAL_SPECIES_PRIORS_FILE, options.GetSpeciesPriorsFile()});
        if (fs::exists(options.GetSpeciesNeighboursFile())) sources.push_back({Options::PROTAL_SPECIES_NEIGHBOURS_FILE, options.GetSpeciesNeighboursFile()});
        if (fs::exists(options.GetCongenerGapsFile())) sources.push_back({Options::PROTAL_CONGENER_GAPS_FILE, options.GetCongenerGapsFile()});
        if (fs::exists(options.GetStrainAllelesFile())) sources.push_back({Options::PROTAL_STRAIN_ALLELES_FILE, options.GetStrainAllelesFile()});
        if (fs::exists(options.GetColumnWeightsFile())) sources.push_back({Options::PROTAL_COLUMN_WEIGHTS_FILE, options.GetColumnWeightsFile()});
        if (fs::exists(options.GetForeignRatesFile())) sources.push_back({Options::PROTAL_FOREIGN_RATES_FILE, options.GetForeignRatesFile()});
        if (fs::exists(options.GetGeneNeighboursFile())) sources.push_back({Options::PROTAL_GENE_NEIGHBOURS_FILE, options.GetGeneNeighboursFile()});
        if (fs::exists(options.GetGenePositionsFile())) sources.push_back({Options::PROTAL_GENE_POSITIONS_FILE, options.GetGenePositionsFile()});
        for (auto const& model : AllModelFiles()) {
            std::string const path = (fs::path(options.GetLocation().dir) / model).string();
            if (fs::exists(path)) sources.push_back({model, path});
        }
        return sources;
    }

    // The models after the other members, each group in its order: at the end of database.protal, --add_model
    // replaces them in place (db::InPlaceFrom).
    inline void ModelsLast(std::vector<db::Source>& sources) {
        auto const models = AllModelFiles();
        std::stable_partition(sources.begin(), sources.end(), [&models](db::Source const& source) {
            return std::find(models.begin(), models.end(), source.name) == models.end();
        });
    }

    // Whether names hold a model of reads of `type` (one of its ModelCandidates).
    inline bool HasModel(std::vector<std::string> const& names, ReadType type) {
        auto const candidates = ModelCandidates(type);
        return std::any_of(candidates.begin(), candidates.end(), [&](std::string const& c) {
            return std::find(names.begin(), names.end(), c) != names.end();
        });
    }

    // The read types (their tokens) whose model is among names, and those whose is not.
    inline std::pair<std::string, std::string> ModelCoverage(std::vector<std::string> const& names) {
        std::string with, without;
        for (auto const& info : kReadTypes) {
            std::string& list = HasModel(names, info.type) ? with : without;
            list += (list.empty() ? "" : ", ") + info.token;
        }
        return {with, without};
    }

    // Writes the binary gene table (GeneTableFile.h) of these tables at target: reads reference.map and unique_kmers.tsv (if
    // any) as a run does (GenomeLoader checks them) and records their sizes and reference.map's fingerprint. Exits 8 on failure.
    static void WriteGeneTableFile(db::DbFile const& fna, db::DbFile const& map, std::optional<db::DbFile> const& unique,
                                   std::string const& target, int threads) {
        Benchmark bm("Write " + gene_table_file::kFileName);
        bm.Start();
        GenomeLoader loader(fna, map, threads);
        if (unique) loader.LoadUniqueKmers(*unique, threads);
        std::string const partial = target + ".partial";
        auto const error = loader.WriteGeneTable(partial, map.Size().value_or(0), unique ? unique->Size().value_or(0) : 0,
                                                 unique.has_value());
        std::error_code ec;
        if (error.empty()) std::filesystem::rename(partial, target, ec);
        if (!error.empty() || ec) {
            std::cerr << "Writing " << target << " failed: " << (error.empty() ? ec.message() : error) << std::endl;
            std::filesystem::remove(partial, ec);
            exit(8);
        }
        bm.Stop();
        std::cout << "Wrote " << target << ": " << loader.GeneCount() << " genes, the binary form of " << map.Name()
                  << (unique ? " and " + unique->Name() : std::string()) << ", which a run from database.protal loads" << std::endl;
        bm.PrintResults();
    }

    // The index of --build as the first member of database.protal, written from memory (db::Generated): its frames
    // straight into the database file, read back from it and compared with the index, which is then freed, before the
    // other members are compressed. Before, the index was written as index.prx.zst and read back, then copied into the
    // database and both read again: ~36 GB written and ~72 GB read at r226 for an 18 GB index
    // (docs/claude/2026-10-07-database-build-audit.md, S2). The database's bytes are the same.
    inline db::Source IndexMember(Seedmap& index, zstd::Params const& params) {
        std::string error;
        auto plan = index.CompressedFrames(params, error);
        if (!plan) {
            std::cerr << "Cannot write the index: " << error << std::endl;
            exit(8);
        }
        auto const frames = std::make_shared<index_codec::FramePlan>(std::move(*plan));
        auto const timer = std::make_shared<Benchmark>("Write index");
        auto generated = std::make_shared<db::Generated>();
        generated->frames = frames->Frames();
        generated->write = [&index, frames, params, timer](zstd::FrameWriter& out) -> std::string {
            timer->Start();
            std::cout << "Write the index into " << db::kFileName << " (zstd level " << params.level << ", columns in "
                      << HumanBytes(params.frame_size) << " chunks, verified, " << params.threads << " thread(s))" << std::endl;
            std::string error;
            size_t raw_chunks = 0;
            uint64_t const before = out.Written();
            if (!index.WriteCompressedFrames(out, *frames, params, error, raw_chunks)) return error;
            uint64_t const written = out.Written() - before, size = index.SerializedSize();
            if (raw_chunks > 0) std::cout << raw_chunks << " index chunk(s) kept as raw cells" << std::endl;
            std::cout << "Index in " << db::kFileName << ": " << HumanBytes(written) << " (" << HumanBytes(size) << " uncompressed, "
                      << std::fixed << std::setprecision(1) << size / double(std::max<uint64_t>(written, 1)) << "x)"
                      << std::defaultfloat << std::endl;
            return "";
        };
        generated->verify = [&index, params](std::string const& path, zstd::SeekTable const& member) {
            return index.VerifyCompressedFrames(path, member, params.threads);
        };
        generated->done = [&index, timer]() {
            index.FreeMemory();  // ~30 GB at r226: gone before the reference is compressed
            ReleaseFreeMemory();
            timer->Stop();
            timer->PrintResults();
            PrintMemory("writing the index");
        };
        db::Source source{ Options::PROTAL_INDEX_FILE, "the index in memory" };
        source.generated = std::move(generated);
        return source;
    }

    // Packs the database folder into database.protal, checks it (db::Write), and removes the files
    // it now holds; --unpack_db writes them back. With `index` (--build), the index is taken from
    // memory (IndexMember) and freed once written; else from the folder, in the column format, as
    // --compress_db writes it.
    static void BundleDatabase(protal::Options const& options, Seedmap* index = nullptr) {
        namespace fs = std::filesystem;
        std::string const target = (fs::path(options.GetLocation().dir) / db::kFileName).string();
        auto sources = BundleSources(options);
        if (index) {
            sources.front() = IndexMember(*index, options.CompressionParams());
        } else if (!index_codec::IsSplitIndex(sources.front().path)) {
            std::cerr << "Cannot write " << target << ": the index " << sources.front().path << " is not in protal's column "
                      << "format (protal --compress_db --no_bundle converts it)" << std::endl;
            exit(8);
        }
        {
            // The binary gene table, from the text tables packed with it (GeneTableFile.h); removed with them below.
            std::string const gene_table = (fs::path(options.GetLocation().dir) / gene_table_file::kFileName).string();
            auto const unique = fs::exists(options.GetUniqueKmersFile()) ? std::optional(db::DbFile::OnDisk(options.GetUniqueKmersFile()))
                                                                          : std::nullopt;
            WriteGeneTableFile(db::DbFile::OnDisk(options.ResolvedSequenceFile()), db::DbFile::OnDisk(options.GetSequenceMapFile()),
                               unique, gene_table, static_cast<int>(std::max<size_t>(options.GetThreads(), 1)));
            sources.push_back({ gene_table_file::kFileName, gene_table });
            ModelsLast(sources);
        }
        {
            std::vector<std::string> names;
            for (auto const& source : sources) names.push_back(source.name);
            auto const [with, without] = ModelCoverage(names);
            std::cout << "Models for read types: " << (with.empty() ? "none" : with)
                      << (without.empty() ? "" : "; none for " + without + " (protal --add_model FILE --read_type TYPE --db " +
                                                 Options::ShellWord(target) + " adds one)") << std::endl;
            if (!HasModel(names, ReadType::Paired)) {
                std::cerr << "Warning: no model for paired-end reads (model_pe.xml) in " << options.GetLocation().dir
                          << "; profiling with " << target << " then needs --model" << std::endl;
            }
        }
        auto const params = options.CompressionParams();
        Benchmark bm("Write " + db::kFileName);
        bm.Start();
        std::error_code ec;
        uint64_t before = 0;
        std::string names;
        for (auto const& source : sources) {
            if (!source.generated) {
                auto const size = fs::file_size(source.path, ec);
                if (!ec) before += size;
            }
            names += (names.empty() ? "" : ", ") + (source.generated ? source.path : fs::path(source.path).filename().string());
        }
        bool const replaces = fs::exists(target, ec);
        std::cout << "Write " << target << " from " << names << " (seekable zstd files' frames as they are, the others "
                  << "compressed at zstd level " << params.level << ", " << HumanBytes(params.frame_size) << " frames; verified)" << std::endl;
        std::string error;
        auto const written = db::Write(target, sources, params, error);
        if (!written) {
            std::cerr << "Writing " << target << " failed: " << error << std::endl;
            exit(8);
        }
        for (auto const& source : sources) {
            if (source.generated) continue;
            fs::remove(source.path, ec);
            if (ec) std::cerr << "Warning: cannot remove " << source.path << ": " << ec.message() << std::endl;
        }
        // An index left beside it would take precedence over the single file.
        for (std::string const& raw : {options.GetIndexFile(), options.GetSequenceFile()}) {
            for (std::string const& path : {raw, raw + zstd::kExtension}) {
                if (fs::remove(path, ec)) std::cout << "Removed the previous " << path << std::endl;
            }
        }
        bm.Stop();
        std::cout << target << ": " << HumanBytes(*written) << " (" << HumanBytes(before) << " as separate files"
                  << (index ? ", besides the index" : "") << (replaces ? "; replaced the previous one" : "")
                  << "). Removed the separate files; " << options.ProgramName() << " --unpack_db --db "
                  << Options::ShellWord(target) << " writes them back." << std::endl;
        bm.PrintResults();
    }

    // After writing a database as separate files: a database.protal next to them is stale.
    static void RemoveStaleBundle(protal::Options const& options) {
        std::string const path = (std::filesystem::path(options.GetLocation().dir) / db::kFileName).string();
        std::error_code ec;
        if (std::filesystem::remove(path, ec)) std::cout << "Removed the previous " << path << std::endl;
    }

    // Writes what zstd::ParallelRead delivers into a file, at the same offsets.
    class FileSink : public zstd::Sink {
    public:
        explicit FileSink(int fd) : m_fd(fd) {}

        void Copy(uint64_t offset, char const* data, size_t size) override {
            while (size > 0 && !m_failed) {
                ssize_t const n = ::pwrite(m_fd, data, size, static_cast<off_t>(offset));
                if (n < 0 && errno == EINTR) continue;
                if (n <= 0) {
                    m_failed = true;
                    break;
                }
                data += n;
                size -= static_cast<size_t>(n);
                offset += static_cast<uint64_t>(n);
            }
        }

        bool Failed() const { return m_failed; }

    private:
        int m_fd;
        std::atomic<bool> m_failed{false};
    };

    // Writes the frames `frames` lists in the file at path as a seekable zstd file at target.
    inline std::string CopyFrames(std::string const& path, zstd::SeekTable const& frames, std::string const& target) {
        int const fd = ::open(path.c_str(), O_RDONLY | O_CLOEXEC);
        if (fd < 0) return "cannot open " + path + ": " + std::strerror(errno);
        zstd::FrameWriter out(target);
        std::vector<char> buffer;
        for (auto const& frame : frames.frames) {
            buffer.resize(frame.compressed_size);
            if (!zstd::PreadAll(fd, buffer.data(), frame.compressed_size, frame.compressed_offset)) {
                ::close(fd);
                return "read error in " + path;
            }
            if (!out.Add(buffer.data(), buffer.size(), frame.decompressed_size)) break;
        }
        ::close(fd);
        out.Finish();
        return out.Error();
    }

    // Writes the members of the single-file database as files into dir: the index as index.prx.zst
    // (its frames as they are) or, with raw_index, as a raw index.prx; every other member
    // uncompressed (reference.fna as --preload_genomes_off needs it). The other variant of the index
    // and the reference there (which would shadow the new file or be stale) is removed.
    static void UnpackBundle(protal::Options const& options, std::string const& dir, bool raw_index) {
        namespace fs = std::filesystem;
        auto const& bundle = *options.GetBundle();
        int const threads = static_cast<int>(std::max<size_t>(options.GetThreads(), 1));
        std::error_code ec;
        fs::create_directories(dir, ec);
        if (ec) {
            std::cerr << "Cannot create " << dir << ": " << ec.message() << std::endl;
            exit(8);
        }
        Benchmark bm("Unpack " + db::kFileName);
        bm.Start();
        for (auto const& member : bundle.Members()) {
            if (member.name == gene_table_file::kFileName) {
                std::cout << "Skip " << member.name << " (the binary form of reference.map and unique_kmers.tsv, which a folder "
                          << "does not use; --compress_db writes it again)" << std::endl;
                continue;
            }
            bool const index = member.name == Options::PROTAL_INDEX_FILE;
            bool const copy = index && !raw_index;
            std::string const target = (fs::path(dir) / (copy ? member.name + zstd::kExtension : member.name)).string();
            std::string const partial = target + ".partial";
            auto fail = [&](std::string const& what) {
                std::cerr << "Writing " << target << " from " << bundle.Path() << " failed: " << what << std::endl;
                fs::remove(partial, ec);
                exit(8);
            };
            std::cout << "Write " << target << std::endl;
            if (copy) {
                std::string error = CopyFrames(bundle.Path(), member.frames, partial);
                if (!error.empty()) fail(error);
                auto const table = zstd::ReadSeekTable(partial, error);
                if (!table || table->frames.size() != member.frames.frames.size() || table->DecompressedSize() != member.Size()) {
                    fail("it does not read back as written");
                }
            } else if (index) {
                Seedmap map;
                map.Load(db::DbFile::InBundle(bundle, member.name), threads);
                std::ofstream os(partial, std::ios::binary);
                map.Save(os);
                os.close();
                if (os.fail()) fail("write error");
            } else {
                int const fd = ::open(partial.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC, 0644);
                if (fd < 0) fail(std::strerror(errno));
                FileSink sink(fd);
                std::string error;
                db::DbFile::InBundle(bundle, member.name).ParallelRead(threads, sink, error);
                bool const closed = ::close(fd) == 0;
                if (!error.empty()) fail(error);
                if (sink.Failed() || !closed) fail("write error");
                if (fs::file_size(partial, ec) != member.Size()) fail("it has the wrong size");
            }
            fs::rename(partial, target, ec);
            if (ec) fail("cannot rename " + partial + ": " + ec.message());
            if (index || member.name == Options::PROTAL_SEQUENCE_FILE) {
                std::string const other = copy ? (fs::path(dir) / member.name).string() : target + zstd::kExtension;
                if (fs::remove(other, ec)) std::cout << "Removed the previous " << other << std::endl;
            }
        }
        bm.Stop();
        bm.PrintResults();
    }

    // --unpack_db: the single-file database's files into --unpack_dir (default: its folder); the
    // single file is kept.
    static void UnpackDatabase(protal::Options const& options) {
        std::string const dir = options.UnpackDir();
        std::cout << "Unpack " << options.GetBundle()->Path() << " into " << dir << std::endl;
        UnpackBundle(options, dir, false);
        std::cout << "Unpacked " << options.GetBundle()->Path() << " into " << dir << "; protal --db " << Options::ShellWord(dir)
                  << " uses these files (" << options.GetBundle()->Path() << " can be removed)" << std::endl;
    }

    // --add_model (the models checked already): stores each PMML file of `models` (file, member name) as member or
    // file of the database, replacing the one there. In database.protal the models are the last members, which are
    // replaced in place (db::ReplaceTail: the models and the seek table written, the rest of the file untouched); a
    // file with other members than before, or written before the models went last, is rewritten once instead, via
    // database.protal.partial, its other members' frames copied as they are and the models put last (db::Write; ~20 min
    // for a GTDB database on a network file system). Either way the result is checked. In a folder of separate files
    // each model is copied next to them.
    static void AddModel(protal::Options const& options, std::vector<std::pair<std::string, std::string>> const& models,
                         std::string const& what = "models") {
        namespace fs = std::filesystem;
        std::error_code ec;
        if (!options.IsBundle()) {
            for (auto const& [model, name] : models) {
                std::string const target = (fs::path(options.GetLocation().dir) / name).string();
                bool const replaces = fs::exists(target, ec);
                std::string const partial = target + ".partial";
                fs::copy_file(model, partial, fs::copy_options::overwrite_existing, ec);
                std::ifstream a(model, std::ios::binary), b(partial, std::ios::binary);
                bool const same = !ec && std::equal(std::istreambuf_iterator<char>(a), std::istreambuf_iterator<char>(),
                                                    std::istreambuf_iterator<char>(b), std::istreambuf_iterator<char>());
                if (same) fs::rename(partial, target, ec);
                if (!same || ec) {
                    std::cerr << "Writing " << target << " failed" << (ec ? ": " + ec.message() : "") << std::endl;
                    fs::remove(partial, ec);
                    exit(8);
                }
                std::cout << "Stored " << model << " as " << target << (replaces ? " (replaced the previous one)" : "") << std::endl;
            }
            return;
        }
        auto const& bundle = *options.GetBundle();
        std::vector<db::Source> sources;
        std::set<std::string> replaced;
        for (auto const& member : bundle.Members()) {
            auto const it = std::find_if(models.begin(), models.end(), [&member](auto const& m) { return m.second == member.name; });
            if (it != models.end()) {
                sources.push_back({member.name, it->first});
                replaced.insert(member.name);
            } else {
                sources.push_back({member.name, bundle.Path(), member.frames});
            }
        }
        for (auto const& [model, name] : models) {
            if (!replaced.contains(name)) sources.push_back({name, model});
        }
        ModelsLast(sources);
        auto const params = options.CompressionParams();
        Benchmark bm("Store the " + what + " in " + db::kFileName);
        bm.Start();
        std::string error;
        std::optional<uint64_t> written;
        if (auto const first = db::InPlaceFrom(bundle, sources)) {
            std::cout << "Replace the last " << sources.size() - *first << " of " << sources.size() << " members of " << bundle.Path()
                      << " in place (the others are not rewritten; compressed at zstd level " << params.level << "; checked, the old "
                      << "bytes kept in " << bundle.Path() << db::kJournalExtension << " until then)" << std::endl;
            written = db::ReplaceTail(bundle, sources, *first, params, error);  // its error says whether the database is unchanged
        } else {
            std::cout << "Rewrite " << bundle.Path() << " (its other members' frames copied as they are, the " << what << " compressed at "
                      << "zstd level " << params.level << ", the models last, where later --add_model runs replace them in place; "
                      << "verified)" << std::endl;
            written = db::Write(bundle.Path(), sources, params, error);
            if (!written) error += " (the database is unchanged)";
        }
        if (!written) {
            std::cerr << "Writing " << bundle.Path() << " failed: " << error << std::endl;
            exit(8);
        }
        bm.Stop();
        std::vector<std::string> names;
        for (auto const& source : sources) names.push_back(source.name);
        for (auto const& [model, name] : models) {
            std::cout << "Stored " << model << " as " << name << " in " << bundle.Path()
                      << (replaced.contains(name) ? " (replaced the previous one)" : "") << std::endl;
        }
        std::cout << bundle.Path() << ": " << HumanBytes(*written) << ". Models for read types: " << ModelCoverage(names).first
                  << std::endl;
        bm.PrintResults();
    }

    // Whether species_priors.tsv `file` is of the database's species: each of its taxids a species of the database's
    // taxonomy with the same representative genome (a table converted from another release would give other species'
    // values). Returns an error text, or "".
    static std::string PriorsMatchTaxonomy(std::string const& file, protal::Options const& options) {
        auto const taxonomy_file = options.TaxonomyDbFile();
        auto input = taxonomy_file.Open();
        if (!taxonomy_file.Exists() || !input->IsOpen()) return "the database's taxonomy cannot be read";
        taxonomy::IntTaxonomy const tax(input->Stream(), taxonomy_file.Name());
        std::ifstream is(file);
        std::string line;
        size_t number = 0;
        while (std::getline(is, line)) {
            number++;
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty() || line[0] == '#' || line.rfind("taxid", 0) == 0) continue;
            size_t const tab = line.find('\t');
            size_t const end = tab == std::string::npos ? std::string::npos : line.find('\t', tab + 1);
            int const taxid = std::stoi(line.substr(0, tab));
            std::string const rep = tab == std::string::npos ? "" : line.substr(tab + 1, end == std::string::npos ? std::string::npos : end - tab - 1);
            if (!tax.map.contains(taxid) || tax.map.at(taxid).rank != "species") {
                return "line " + std::to_string(number) + ": taxid " + std::to_string(taxid) + " is no species of the database";
            }
            auto const& known = tax.map.at(taxid).rep_genome;
            if (!known.empty() && !rep.empty() && known != rep) {
                return "line " + std::to_string(number) + ": taxid " + std::to_string(taxid) + " is " + rep + " in the table but " +
                       known + " in the database (converted from another GTDB release?)";
            }
        }
        return "";
    }

    // --add_tables: stores tables (congener_gaps.tsv, strain_alleles.tsv, foreign_rates.tsv, species_priors.tsv) in the database by their file
    // names, as AddModel stores models (a single-file database rewritten once, or replaced at its tail in place; a folder
    // gets the files copied beside the others). Each table is read and checked first; exits 8 if one is unusable.
    static void AddTables(protal::Options const& options) {
        namespace fs = std::filesystem;
        std::vector<std::pair<std::string, std::string>> tables;  // file, member
        for (auto const& file : options.AddTables()) {
            std::string const name = fs::path(file).filename().string();
            std::ifstream is(file);
            std::string error = is ? "" : "cannot be opened";
            std::string summary;
            if (error.empty() && name == Options::PROTAL_CONGENER_GAPS_FILE) {
                congener_gaps::Table table;
                error = table.Read(is);
                summary = std::to_string(table.Copies()) + " gene copies of " + std::to_string(table.Species()) + " species";
            } else if (error.empty() && name == Options::PROTAL_STRAIN_ALLELES_FILE) {
                strain_alleles::Table table;
                error = table.Read(is);
                summary = std::to_string(table.Alleles()) + " alleles of " + std::to_string(table.Copies()) + " gene copies of " +
                          std::to_string(table.Species()) + " species";
            } else if (error.empty() && name == Options::PROTAL_SPECIES_PRIORS_FILE) {
                species_priors::Table table;
                error = table.Read(is);
                if (error.empty()) error = PriorsMatchTaxonomy(file, options);
                summary = std::to_string(table.Size()) + " species, " + std::to_string(table.WithGenomeSize()) + " with a genome size";
            } else if (error.empty() && name == Options::PROTAL_COLUMN_WEIGHTS_FILE) {
                column_weights::Table table;
                error = table.Read(is);
                summary = std::to_string(table.Families()) + " family-gene rows, " + std::to_string(table.Copies()) + " gene copies of " +
                          std::to_string(table.Species()) + " species mapped";
            } else if (error.empty()) {
                foreign_rates::Table table;
                error = table.Read(is);
                summary = std::to_string(table.Copies()) + " gene copies of " + std::to_string(table.Species()) + " species";
            }
            if (!error.empty()) {
                std::cerr << "--add_tables: " << file << " is no valid " << name << ": " << error << std::endl;
                exit(8);
            }
            std::cout << file << ": " << summary << std::endl;
            tables.emplace_back(file, name);
        }
        AddModel(options, tables, "tables");
    }

    // --compress_db on a single-file database: gives it the binary gene table (GeneTableFile.h) if it has none, or one made from
    // other tables, rewriting database.protal once with its other members' frames copied as they are (as AddModel).
    static void AddGeneTable(protal::Options const& options) {
        namespace fs = std::filesystem;
        auto const& bundle = *options.GetBundle();
        int const threads = static_cast<int>(std::max<size_t>(options.GetThreads(), 1));
        auto const fna = db::DbFile::InBundle(bundle, Options::PROTAL_SEQUENCE_FILE);
        auto const map = db::DbFile::InBundle(bundle, Options::PROTAL_SEQUENCE_MAP_FILE);
        auto const unique_file = db::DbFile::InBundle(bundle, Options::PROTAL_UNIQUE_KMER_FILE);
        auto const unique = unique_file.Exists() ? std::optional(unique_file) : std::nullopt;
        if (auto const current = db::DbFile::InBundle(bundle, gene_table_file::kFileName); current.Exists()) {
            GenomeLoader check(fna, map, threads, current, unique ? unique->Size() : std::nullopt);
            if (check.FromGeneTable()) {
                std::cout << bundle.Path() << " is a single-file database with a current " << gene_table_file::kFileName << "; kept" << std::endl;
                return;
            }
        }
        std::string const table = bundle.Path() + "." + gene_table_file::kFileName;  // beside it, removed below
        WriteGeneTableFile(fna, map, unique, table, threads);
        std::vector<db::Source> sources;
        for (auto const& member : bundle.Members()) {
            if (member.name != gene_table_file::kFileName) sources.push_back({ member.name, bundle.Path(), member.frames });
        }
        sources.push_back({ gene_table_file::kFileName, table });
        ModelsLast(sources);
        auto const params = options.CompressionParams();
        std::cout << "Rewrite " << bundle.Path() << " with " << gene_table_file::kFileName << " (the other members' frames copied as "
                  << "they are, the table compressed at zstd level " << params.level << "; verified)" << std::endl;
        Benchmark bm("Rewrite " + db::kFileName);
        bm.Start();
        std::string error;
        auto const written = db::Write(bundle.Path(), sources, params, error);
        std::error_code ec;
        fs::remove(table, ec);
        if (!written) {
            std::cerr << "Writing " << bundle.Path() << " failed: " << error << " (the database is unchanged)" << std::endl;
            exit(8);
        }
        bm.Stop();
        std::cout << bundle.Path() << ": " << HumanBytes(*written) << std::endl;
        bm.PrintResults();
    }

    // --compress_db: rewrites an existing database's index and reference compressed, without
    // rebuilding it, e.g. a downloaded raw database or one compressed as a single frame: the index
    // in the column format (IndexCodec.h), reference.fna as seekable zstd, all packed into
    // database.protal unless --no_bundle. Each new file is read back and compared with the old
    // content before the old file is removed.
    static void CompressDatabase(protal::Options const& options) {
        if (options.IsBundle()) {
            AddGeneTable(options);
            return;
        }
        auto const params = options.CompressionParams();
        CompressionHint(options);
        if (params.frame_size > 0) CompressIndexFile(options, params);
        if (options.WriteBundle()) {
            BundleDatabase(options);  // compresses reference.fna itself
            return;
        }
        for (std::string const& raw : params.frame_size > 0 ? std::vector<std::string>{options.GetSequenceFile()}
                                                            : std::vector<std::string>{options.GetIndexFile(), options.GetSequenceFile()}) {
            std::string const source = zstd::Resolve(raw);
            std::string const target = raw + zstd::kExtension;
            if (!std::filesystem::exists(source)) {
                std::cerr << "Cannot compress " << raw << ": neither it nor " << target << " exists" << std::endl;
                exit(8);
            }
            if (source == target && params.frame_size > 0 && zstd::IsSeekable(source)) {
                std::cout << target << " is already seekable; kept (decompress it first to recompress)" << std::endl;
                continue;
            }
            Benchmark bm("Compress " + std::filesystem::path(raw).filename().string());
            bm.Start();
            std::cout << "Write " << target << " from " << source << " (zstd level " << params.level << ", "
                      << FrameDescription(params) << ", " << params.threads << " thread(s), verified)" << std::endl;
            std::error_code ec;
            uint64_t const before = std::filesystem::file_size(source, ec);
            std::string error;
            if (!zstd::CompressFile(source, target, params, true, error)) {
                std::cerr << "Compressing " << source << " failed: " << error << std::endl;
                exit(8);
            }
            uint64_t const after = std::filesystem::file_size(target, ec);
            if (source != target) {
                std::filesystem::remove(source, ec);
                if (ec) std::cerr << "Warning: cannot remove " << source << ": " << ec.message() << std::endl;
            }
            bm.Stop();
            std::cout << target << ": " << HumanBytes(after) << " (was " << HumanBytes(before) << " as "
                      << std::filesystem::path(source).filename().string() << ")" << std::endl;
            bm.PrintResults();
        }
    }

    // --decompress_db: writes the database's index.prx and reference.fna raw again (e.g. for
    // --preload_genomes_off or protal versions without zstd support) and removes the .zst files, or
    // database.protal after writing all its files raw next to it.
    static void DecompressDatabase(protal::Options const& options) {
        int const threads = static_cast<int>(std::max<size_t>(options.GetThreads(), 1));
        if (options.IsBundle()) {
            std::string const bundle = options.GetBundle()->Path();
            UnpackBundle(options, options.GetLocation().dir, true);
            std::error_code ec;
            std::filesystem::remove(bundle, ec);
            if (ec) std::cerr << "Warning: cannot remove " << bundle << ": " << ec.message() << std::endl;
            else std::cout << "Removed " << bundle << std::endl;
            return;
        }
        {
            std::string const raw = options.GetIndexFile(), source = zstd::Resolve(raw), partial = raw + ".partial";
            if (source == raw) {
                std::cout << raw << " is not compressed; kept" << std::endl;
            } else {
                std::cout << "Write " << raw << " from " << source << std::endl;
                Seedmap map;
                map.Load(source, threads);
                std::ofstream os(partial, std::ios::binary);
                map.Save(os);
                os.close();
                if (os.fail()) {
                    std::cerr << "Writing " << partial << " failed" << std::endl;
                    exit(8);
                }
                ReplaceFile(partial, raw, source);
            }
        }
        {
            std::string const raw = options.GetSequenceFile(), source = zstd::Resolve(raw), partial = raw + ".partial";
            if (source == raw) {
                std::cout << raw << " is not compressed; kept" << std::endl;
            } else {
                std::cout << "Write " << raw << " from " << source << std::endl;
                zstd::InputFile in(source);
                std::ofstream out(partial, std::ios::binary);
                std::vector<char> buffer(size_t{8} << 20);
                while (in.Stream()) {
                    in.Stream().read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
                    out.write(buffer.data(), in.Stream().gcount());
                }
                out.close();
                if (!in.IsOpen() || in.Stream().bad() || out.fail()) {
                    std::cerr << "Decompressing " << source << " failed" << std::endl;
                    std::error_code ec;
                    std::filesystem::remove(partial, ec);
                    exit(8);
                }
                ReplaceFile(partial, raw, source);
            }
        }
    }

    // A row per gene id of each reference taxon, for the unique k-mer table (Seedmap::CountUniqueKmers).
    inline Seedmap::GeneRows GeneRowsOf(GenomeLoader& genomes) {
        size_t max_taxid = 0;
        for (auto const& [taxid, genome] : genomes.GetGenomeMap()) max_taxid = std::max<size_t>(max_taxid, taxid);
        Seedmap::GeneRows rows;
        rows.first_row.assign(max_taxid + 2, 0);
        for (auto const& [taxid, genome] : genomes.GetGenomeMap()) rows.first_row[taxid + 1] = genome.GetGeneList().size();
        std::partial_sum(rows.first_row.begin(), rows.first_row.end(), rows.first_row.begin());
        return rows;
    }

    // The genes' conservation factors (GeneConservation.h) from the copies of the reference genes in
    // full_reference (every genome's, >taxid_geneid), read with -t threads, against the reference genes
    // in genomes. Written to target (gene_conservation.tsv) if any gene has a factor; a target left by
    // an earlier build is removed otherwise. Without other genomes' copies (no --full_reference, or one
    // genome per species) every gene keeps the whole depth identity margin.
    // gene_neighbours.tsv, which scripts/mini_db/gene_neighbours.py writes into the folder before the build and
    // BundleSources packs: read and checked against the database's genes and taxonomy, so that a query never
    // meets a bad one. Exits 8 at the first problem; says so if there is none.
    static void CheckGeneNeighbours(protal::Options const& options, GenomeLoader& genomes) {
        std::string const path = options.GetGeneNeighboursFile();
        if (!std::filesystem::exists(path)) {
            std::cout << "Gene neighbours: none (" << Options::PROTAL_GENE_NEIGHBOURS_FILE << ", which "
                      << "scripts/mini_db/gene_neighbours.py writes from whole genomes, is not in the folder)" << std::endl;
            return;
        }
        if (options.HasBuildGeneSubset()) {
            // The table was counted over every gene: a gene's neighbour there may be one the subset lacks, so
            // the gene a read really meets next (the nearest one of the subset) would be judged unlikely, and
            // long reads would look for genes that are not in the index. The folder has to be derived for
            // the subset (the frequencies counted anew from gene_positions.tsv over the subset's genes).
            std::cerr << "Cannot build with --build_gene_subset: the folder has " << Options::PROTAL_GENE_NEIGHBOURS_FILE
                      << ", whose neighbours were counted over every gene. Derive a folder of the subset's genes "
                      << "with scripts/mini_db/gtdb_to_protal_db.py --from_db <folder> --genes <list> --outdir <subset "
                      << "folder> (its gene neighbours are counted anew from " << Options::PROTAL_GENE_POSITIONS_FILE
                      << ") and build that, or remove the table (and " << Options::PROTAL_GENE_POSITIONS_FILE
                      << ") from the folder" << std::endl;
            exit(8);
        }
        auto fail = [&path](std::string const& what) {
            std::cerr << "Invalid gene neighbours " << path << ": " << what << std::endl;
            exit(8);
        };
        std::ifstream is(path);
        gene_neighbours::Table table;
        if (!is) fail("cannot open the file");
        if (auto const error = table.Read(is); !error.empty()) fail(error);

        std::vector<bool> gene_ids;
        for (auto const& [taxid, genome] : genomes.GetGenomeMap()) {
            auto const& genes = genome.GetGeneList();
            if (genes.size() >= gene_ids.size()) gene_ids.resize(genes.size() + 1, false);
            for (size_t i = 0; i < genes.size(); i++) gene_ids[i + 1] = gene_ids[i + 1] || genes[i].IsSet();
        }
        tsl::sparse_set<uint32_t> taxa;
        std::ifstream taxonomy(options.GetInternalTaxonomyFile());
        std::string line;
        while (std::getline(taxonomy, line)) {
            uint32_t id = 0;
            if (std::from_chars(line.data(), line.data() + line.size(), id).ec == std::errc()) taxa.insert(id);
        }
        std::ifstream again(path);
        size_t number = 0;
        while (std::getline(again, line)) {
            number++;
            if (line.empty() || line[0] == '#' || line.rfind("clade", 0) == 0) continue;
            uint32_t v[4] = {};
            char const* p = line.data();
            char const* const end = line.data() + line.size();
            for (auto& x : v) {
                p = std::from_chars(p, end, x).ptr;
                if (p < end) p++;
            }
            auto known = [&gene_ids](uint32_t gene) { return gene < gene_ids.size() && gene_ids[gene]; };
            if (!taxa.contains(v[0])) fail("line " + std::to_string(number) + ": clade " + std::to_string(v[0]) + " is not in the taxonomy");
            if (!known(v[1]) || (v[3] != 0 && !known(v[3]))) {
                fail("line " + std::to_string(number) + ": gene " + std::to_string(known(v[1]) ? v[3] : v[1]) + " is not in the database");
            }
        }
        std::cout << "Gene neighbours: " << table.Rules() << " rules of " << table.Clades() << " clades from "
                  << table.Genomes() << " genomes, neighbours up to " << table.MaxGap() << " bases apart ("
                  << Options::PROTAL_GENE_NEIGHBOURS_FILE << "), stored in the database" << std::endl;
    }

    // gene_positions.tsv, which scripts/mini_db/gene_neighbours.py writes beside gene_neighbours.tsv: where each gene
    // was placed in each genome the neighbours were counted from (accession, taxid, contig, contig_length, circular,
    // gene, start, end, strand, placed, kmer_share). A run does not read it; the database keeps it, checked line by
    // line against the database's genes, so that the frequencies can be traced back to the genomes and derived anew.
    // Exits 8 at the first problem; says so if there is none.
    static void CheckGenePositions(protal::Options const& options, GenomeLoader& genomes) {
        std::string const path = options.GetGenePositionsFile();
        if (!std::filesystem::exists(path)) return;
        auto fail = [&path](size_t number, std::string const& what) {
            std::cerr << "Invalid gene positions " << path << ": line " << number << ": " << what << std::endl;
            exit(8);
        };
        std::ifstream is(path);
        if (!is) fail(0, "cannot open the file");
        std::string line;
        size_t number = 0, placements = 0, traced = 0;
        tsl::sparse_set<std::string> genome_names, circular;
        tsl::sparse_set<uint32_t> species;
        std::vector<std::string_view> f;
        while (std::getline(is, line)) {
            number++;
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty() || line[0] == '#' || line.rfind("accession", 0) == 0) continue;
            f.clear();
            for (size_t at = 0;;) {
                size_t const tab = line.find('\t', at);
                f.emplace_back(std::string_view(line).substr(at, tab == std::string::npos ? std::string::npos : tab - at));
                if (tab == std::string::npos) break;
                at = tab + 1;
            }
            if (f.size() != 11) fail(number, "expected 11 fields (accession, taxid, contig, contig_length, circular, gene, start, "
                                             "end, strand, placed, kmer_share)");
            auto number_of = [&](std::string_view field, char const* what) {
                uint64_t v = 0;
                auto const r = std::from_chars(field.data(), field.data() + field.size(), v);
                if (r.ec != std::errc() || r.ptr != field.data() + field.size()) fail(number, std::string(what) + " must be a number");
                return v;
            };
            uint64_t const taxid = number_of(f[1], "taxid"), length = number_of(f[3], "contig_length");
            uint64_t const gene = number_of(f[5], "gene"), start = number_of(f[6], "start"), end = number_of(f[7], "end");
            if (f[0].empty() || f[2].empty()) fail(number, "accession and contig must not be empty");
            if (f[4] != "0" && f[4] != "1") fail(number, "circular must be 0 or 1");
            if (start < 1 || start > end || end > length) fail(number, "expected 1 <= start <= end <= contig_length");
            if (f[8] != "+" && f[8] != "-") fail(number, "strand must be + or -");
            if (f[9] != "exact" && f[9] != "trace") fail(number, "placed must be exact or trace");
            double share = -1;
            auto const r = std::from_chars(f[10].data(), f[10].data() + f[10].size(), share);
            if (r.ec != std::errc() || share < 0 || share > 1) fail(number, "kmer_share must be a number from 0 to 1");
            if (taxid > UINT32_MAX || gene > UINT32_MAX ||
                !genomes.HasGene(static_cast<uint32_t>(taxid), static_cast<uint32_t>(gene))) {
                fail(number, "species " + std::to_string(taxid) + " has no gene " + std::to_string(gene) + " in the database");
            }
            placements++;
            traced += f[9] == "trace";
            genome_names.insert(std::string(f[0]));
            if (f[4] == "1") circular.insert(std::string(f[0]));
            species.insert(static_cast<uint32_t>(taxid));
        }
        if (is.bad()) fail(number, "read error");
        std::cout << "Gene positions: " << placements << " genes (" << traced << " placed by their k-mer trace) in "
                  << genome_names.size() << " genomes of " << species.size() << " species, " << circular.size()
                  << " read as circular (" << Options::PROTAL_GENE_POSITIONS_FILE
                  << "), stored in the database, not read by queries" << std::endl;
    }

    static gene_conservation::Estimate WriteGeneConservation(protal::Options const& options, GenomeLoader& genomes,
                                                             std::string const& full_reference, std::string const& target) {
        Benchmark bm("Gene conservation");
        bm.Start();
        std::vector<uint64_t> keys;
        for (auto const& [taxid, genome] : genomes.GetGenomeMap()) {
            auto const& genes = genome.GetGeneList();
            for (size_t i = 0; i < genes.size(); i++) {
                // Genes outside --build_gene_subset get no factor: nothing can hit them.
                if (genes[i].IsSet() && options.BuildGeneAllowed(i + 1)) keys.push_back(gene_conservation::Estimator::Key(taxid, i + 1));
            }
        }
        gene_conservation::Estimator estimator(std::move(keys));
        auto add = [&genomes, &estimator](std::string const& header, std::string_view sequence) {
            auto const [taxid, geneid] = KmerUtils::ExtractHeaderInformation(header);
            if (!genomes.HasGene(taxid, geneid)) return;
            auto const slot = estimator.Take(taxid, geneid);
            if (!slot) return;
            auto const rep = genomes.GetGenome(taxid).GetGeneOMP(geneid).Sequence();
            estimator.Add(*slot, rep.View(), sequence);
        };
        int const threads = static_cast<int>(std::max<size_t>(options.GetThreads(), 1));
        std::string error;
        auto const frames = FastaFrames::Open(zstd::Resolve(full_reference), error);
        if (!error.empty()) {
            std::cerr << "Cannot read the full reference: " << error << std::endl;
            exit(8);
        }
        if (frames) {
            // A frame per thread at a time, its records in order (FastaFrames: the converter's frames, decompressed on all
            // threads instead of one). A species' copies of a gene lie in one frame (a marker file's, or a gene's), so the
            // copies it gets compared (the first kMaxCopies) are the first in the file, whatever the threads do.
            std::atomic<size_t> next{ 0 };
            std::string failure;
#pragma omp parallel num_threads(static_cast<int>(std::min<size_t>(static_cast<size_t>(threads), frames->Frames())))
            {
                FastaFrames::Reader reader(*frames);
                std::string batch, scratch, header;
                for (size_t f; (f = next.fetch_add(1)) < frames->Frames();) {
                    reader.Take(f);
                    while (reader.Next(batch, size_t{1} << 20)) {
                        ForEachFastaRecord(batch, scratch, [&](std::string_view record_header, std::string_view sequence) {
                            header.assign(record_header);
                            add(header, sequence);
                        });
                    }
                    if (!reader.Error().empty()) {
#pragma omp critical(conservation_failure)
                        if (failure.empty()) failure = reader.Error();
                        break;
                    }
                }
            }
            if (!failure.empty()) {
                std::cerr << "Cannot read the full reference: " << failure << std::endl;
                exit(8);
            }
        } else {
            auto input = OpenInput(full_reference);
            std::istream& is = input->Stream();
            bool parsed = true;
            omp_set_num_threads(threads);
#pragma omp parallel default(none) shared(is, add, parsed)
            {
                FastxRecord record;
                SeqReader reader { is };
                while (reader(record)) add(record.header, record.sequence);
#pragma omp critical(conservation_failure)
                parsed = parsed && reader.Success();
            }
            // A read or decompression error ends the stream early: the factors would rest on part of the copies.
            if (is.bad() || !parsed) {
                std::cerr << "Cannot read the full reference " << full_reference << " to its end (truncated or corrupt file?)" << std::endl;
                exit(8);
            }
        }
        auto estimate = estimator.Finish();
        std::error_code ec;
        std::filesystem::remove(target, ec);
        if (estimate.table.Empty()) {
            std::cout << "Gene conservation: no factors (" << estimate.species_with_copies << " species with other genomes' copies of "
                      << "their genes in " << full_reference << ", " << estimate.species << " of them with enough genes that differ from "
                      << "the representative's): every gene keeps the whole depth identity margin" << std::endl;
        } else {
            std::ofstream os(target);
            estimate.table.Write(os);
            os.close();
            if (!os) {
                std::cerr << "Writing " << target << " failed" << std::endl;
                exit(8);
            }
            auto const [low, high] = estimate.table.Range();
            std::cout << "Gene conservation: factors " << std::setprecision(2) << low << "-" << high << std::setprecision(6)
                      << " for " << estimate.table.Genes() << " genes, from " << estimate.species << " species (" << estimate.copies
                      << " copies compared): " << target << std::endl;
        }
        bm.Stop();
        bm.PrintResults();
        return estimate;
    }

    // The database's genera of two species or more, each its species by taxid, in the order of the genera's taxids. The
    // genus of a species is its nearest ancestor of rank genus in internal_taxonomy.dmp.
    static std::vector<std::vector<uint32_t>> Genera(protal::Options const& options, GenomeLoader& genomes) {
        std::unordered_map<uint32_t, uint32_t> parent;
        std::unordered_map<uint32_t, bool> is_genus;
        {
            std::ifstream taxonomy(options.GetInternalTaxonomyFile());
            std::string line;
            while (std::getline(taxonomy, line)) {
                std::vector<std::string> fields;
                size_t start = 0;
                for (size_t tab; (tab = line.find('\t', start)) != std::string::npos; start = tab + 1) fields.push_back(line.substr(start, tab - start));
                fields.push_back(line.substr(start));
                uint32_t id = 0, up = 0;
                if (fields.size() < 5 || std::from_chars(fields[0].data(), fields[0].data() + fields[0].size(), id).ec != std::errc() ||
                    std::from_chars(fields[1].data(), fields[1].data() + fields[1].size(), up).ec != std::errc()) continue;
                parent[id] = up;
                is_genus[id] = fields[4] == "genus";
            }
        }
        auto genus_of = [&](uint32_t taxid) -> uint32_t {
            uint32_t node = taxid;
            for (int depth = 0; depth < 32; depth++) {
                auto const it = parent.find(node);
                if (it == parent.end() || it->second == node) return 0;
                node = it->second;
                if (is_genus[node]) return node;
            }
            return 0;
        };
        std::map<uint32_t, std::vector<uint32_t>> by_genus;
        for (auto const& [taxid, genome] : genomes.GetGenomeMap()) {
            if (uint32_t const genus = genus_of(static_cast<uint32_t>(taxid))) by_genus[genus].push_back(static_cast<uint32_t>(taxid));
        }
        std::vector<std::vector<uint32_t>> genera;
        for (auto& [_, members] : by_genus) {
            std::sort(members.begin(), members.end());  // by taxid, not in the genome map's order
            if (members.size() >= 2) genera.push_back(std::move(members));
        }
        return genera;
    }

    // gene_congeners.tsv next to the database (gene_conservation::CompareCongeners): how each gene differs between the
    // representatives of congeneric species, against how it differs within species (`within`, the factors just
    // estimated). A report of the build, which queries do not read; build_gtdb_database.py keeps it in model_logs/.
    static void WriteGeneCongeners(protal::Options const& options, GenomeLoader& genomes, gene_conservation::Table const& within) {
        namespace fs = std::filesystem;
        Benchmark bm("Gene congeners");
        bm.Start();
        std::string const target = (fs::path(options.GetGeneConservationFile()).parent_path() /
                                    gene_conservation::kCongenersFileName).string();
        std::error_code ec;
        fs::remove(target, ec);
        auto const genera = Genera(options, genomes);
        if (genera.empty()) {
            std::cout << "Gene congeners: no genus with two species or more in the database; no comparison" << std::endl;
            bm.Stop();
            return;
        }
        auto genes_of = [&genomes, &options](uint32_t taxid) {
            std::vector<std::pair<uint64_t, std::string>> genes;
            auto& genome = genomes.GetGenome(taxid);
            auto const& list = genome.GetGeneList();
            for (size_t i = 0; i < list.size(); i++) {
                if (!list[i].IsSet() || !options.BuildGeneAllowed(i + 1)) continue;  // the subset's genes only
                auto const seq = genome.GetGeneOMP(i + 1).Sequence();
                genes.emplace_back(i + 1, std::string(seq.View()));
            }
            return genes;
        };
        auto const estimate = gene_conservation::CompareCongeners(genera, genes_of, within, static_cast<int>(options.GetThreads()));
        std::ofstream os(target);
        estimate.Write(os);
        os.close();
        if (!os) {
            std::cerr << "Writing " << target << " failed" << std::endl;
            exit(8);
        }
        auto pct = [](double v) {
            if (std::isnan(v)) return std::string("-");
            std::ostringstream o;
            o << std::fixed << std::setprecision(1) << 100 * v << "%";
            return o.str();
        };
        auto num = [](double v) {
            if (std::isnan(v)) return std::string("-");
            std::ostringstream o;
            o << std::fixed << std::setprecision(2) << v;
            return o.str();
        };
        std::cout << "Gene congeners: " << estimate.pairs << " pairs of species of " << estimate.genera << " genera ("
                  << estimate.species << " species); the genes' divergence between congeners correlates "
                  << num(estimate.spearman) << " with their factors (Spearman, " << estimate.correlated << " genes); genes of "
                  << "factor below 1: between-species factor " << num(estimate.conserved_between) << ", the nearest congener's "
                  << "copy identical in " << pct(estimate.conserved_identical) << " of species, nearer than "
                  << gene_conservation::kNearIdentical << " in " << pct(estimate.conserved_near) << "; the other genes: "
                  << num(estimate.fast_between) << ", " << pct(estimate.fast_identical) << ", " << pct(estimate.fast_near)
                  << ": " << target << std::endl;
        bm.Stop();
        bm.PrintResults();
    }

    // species_neighbours.tsv in the database (SpeciesNeighbours.h): every two species of a genus compared by all their
    // references' marker genes (profiler::context::SketchedTaxonDistance; a run's relative_distance compares only the
    // genes with unique k-mers, which the build does not know yet here), and each species' nearest congeners kept (at most species_neighbours::kMaxNeighbours, within
    // kMaxDistance); every species of the database gets a row, also one without a congener. The genera are taken in batches of about kNeighbourBatch
    // species, sketched and compared on all threads, so that the whole database's sketches (~35 kB a species) are never
    // held at once.
    inline constexpr size_t kNeighbourBatch = 10000;

    static void WriteSpeciesNeighbours(protal::Options const& options, GenomeLoader& genomes, std::string const& target) {
        Benchmark bm("Species neighbours");
        bm.Start();
        auto const genera = Genera(options, genomes);
        species_neighbours::Table table;
        for (auto const& [taxid, _] : genomes.GetGenomeMap()) table.Set(static_cast<uint32_t>(taxid), {});
        int const threads = std::max(1, static_cast<int>(options.GetThreads()));
        auto take = [&options](uint32_t gene) { return options.BuildGeneAllowed(gene); };
        size_t compared = 0;
        for (size_t g = 0; g < genera.size();) {
            std::vector<uint32_t> members;
            std::vector<std::pair<size_t, size_t>> ranges;  // each genus's members, [begin, end) of `members`
            while (g < genera.size() && (members.empty() || members.size() + genera[g].size() <= kNeighbourBatch)) {
                ranges.emplace_back(members.size(), members.size() + genera[g].size());
                members.insert(members.end(), genera[g].begin(), genera[g].end());
                g++;
            }
            std::vector<profiler::context::TaxonSketch> sketches(members.size());
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
            for (int64_t i = 0; i < static_cast<int64_t>(members.size()); i++) {
                sketches[static_cast<size_t>(i)] = profiler::context::ReferenceSketch(genomes, members[static_cast<size_t>(i)], take);
            }
            std::vector<std::pair<uint32_t, uint32_t>> pairs;
            for (auto const [begin, end] : ranges) {
                for (size_t i = begin; i < end; i++) {
                    for (size_t j = i + 1; j < end; j++) pairs.emplace_back(static_cast<uint32_t>(i), static_cast<uint32_t>(j));
                }
            }
            std::vector<float> distance(pairs.size());
#pragma omp parallel num_threads(threads)
            {
                std::vector<double> scratch;
#pragma omp for schedule(dynamic, 1024)
                for (int64_t p = 0; p < static_cast<int64_t>(pairs.size()); p++) {
                    auto const [i, j] = pairs[static_cast<size_t>(p)];
                    double const d = profiler::context::SketchedTaxonDistance(sketches[i], sketches[j], scratch);
                    distance[static_cast<size_t>(p)] = static_cast<float>(d);
                }
            }
            compared += pairs.size();
            std::vector<std::vector<species_neighbours::Neighbour>> near(members.size());
            for (size_t p = 0; p < pairs.size(); p++) {
                if (!(distance[p] <= species_neighbours::kMaxDistance)) continue;
                auto const [i, j] = pairs[p];
                near[i].push_back({ members[j], distance[p] });
                near[j].push_back({ members[i], distance[p] });
            }
            for (size_t i = 0; i < members.size(); i++) table.Set(members[i], std::move(near[i]));
        }
        std::ofstream os(target);
        table.Write(os);
        os.close();
        if (!os) {
            std::cerr << "Writing " << target << " failed" << std::endl;
            exit(8);
        }
        size_t within_01 = 0, within_02 = 0;
        for (auto const& [taxid, _] : genomes.GetGenomeMap()) {
            within_01 += table.Within(static_cast<uint32_t>(taxid), 0.01) > 0;
            within_02 += table.Within(static_cast<uint32_t>(taxid), 0.02) > 0;
        }
        std::cout << "Species neighbours: " << compared << " pairs of congeners compared in " << genera.size() << " genera; "
                  << table.Pairs() << " neighbours within " << species_neighbours::kMaxDistance << " listed for " << table.Species()
                  << " species, " << within_01 << " with a congener within 0.01, " << within_02 << " within 0.02: " << target << std::endl;
        bm.Stop();
        bm.PrintResults();
    }

    // --write_species_neighbours FILE: the species neighbours of a folder's references, written to FILE without a build: a
    // converted release before its databases exist (build_gtdb_database.py chooses the species to hold out by it, so that a
    // species complex is held out or kept whole), or a database of separate files. The genomes are loaded as for a build
    // (reference.fna by reference.map), the genera read from internal_taxonomy.dmp. Exits 8 when a file is missing.
    static void WriteSpeciesNeighboursOnly(protal::Options const& options) {
        namespace fs = std::filesystem;
        for (auto const& file : { options.ResolvedSequenceFile(), options.GetSequenceMapFile(), options.GetInternalTaxonomyFile() }) {
            if (!fs::is_regular_file(file)) {
                std::cerr << "--write_species_neighbours: " << file << " is missing (a converted release's folder has reference.fna, "
                          << "reference.map and internal_taxonomy.dmp)" << std::endl;
                exit(8);
            }
        }
        int const threads = static_cast<int>(std::max<size_t>(options.GetThreads(), 1));
        GenomeLoader genomes(options.SequenceDbFile(), options.SequenceMapDbFile(), threads);
        Benchmark bm("Preload genomes");
        bm.Start();
        genomes.LoadAllGenomes(threads);
        bm.Stop();
        bm.PrintResults();
        std::cout << genomes.PreloadTimes() << std::endl;
        WriteSpeciesNeighbours(options, genomes, options.WriteSpeciesNeighboursFile());
    }

    // congener_gaps.tsv in the database (CongenerGaps.h): gene by gene, every species' copy aligned against the copies of
    // its congeners (all of them up to congener_gaps::Settings::all others, else its nearest by sketch and a hashed
    // sample), its nearest congener's distance and the sample's median kept. Each gene's sequences are held while it is
    // compared (a gene's copies at r226: ~140 MB), the alignments on all threads.
    static congener_gaps::Table WriteCongenerGaps(protal::Options const& options, GenomeLoader& genomes) {
        namespace fs = std::filesystem;
        Benchmark bm("Congener gaps");
        bm.Start();
        std::string const target = options.GetCongenerGapsFile();
        std::error_code ec;
        fs::remove(target, ec);
        auto const genera = Genera(options, genomes);
        std::unordered_map<uint32_t, uint32_t> genus_of;
        for (size_t g = 0; g < genera.size(); g++) {
            for (uint32_t const taxid : genera[g]) genus_of[taxid] = static_cast<uint32_t>(g + 1);
        }
        std::vector<uint32_t> taxids;
        size_t max_gene = 0;
        for (auto const& [taxid, genome] : genomes.GetGenomeMap()) {
            if (genus_of.contains(static_cast<uint32_t>(taxid))) taxids.push_back(static_cast<uint32_t>(taxid));
            max_gene = std::max(max_gene, genome.GetGeneList().size());
        }
        std::sort(taxids.begin(), taxids.end());
        int const threads = static_cast<int>(std::max<size_t>(options.GetThreads(), 1));
        congener_gaps::Settings const settings;
        congener_gaps::Stats stats;
        std::vector<std::tuple<uint32_t, uint32_t, congener_gaps::Gap>> rows;
        std::vector<std::pair<uint32_t, congener_gaps::Gap>> found;
        std::vector<uint32_t> holders;
        std::vector<std::string> seqs;
        auto const genus = [&genus_of](uint32_t taxid) {
            auto const it = genus_of.find(taxid);
            return it == genus_of.end() ? 0u : it->second;
        };
        for (size_t g = 1; g <= max_gene; g++) {
            if (!options.BuildGeneAllowed(g)) continue;  // the subset's genes only
            holders.clear();
            for (uint32_t const taxid : taxids) {
                auto const& list = genomes.GetGenome(taxid).GetGeneList();
                if (g <= list.size() && list[g - 1].IsSet()) holders.push_back(taxid);
            }
            if (holders.size() < 2) continue;
            seqs.assign(holders.size(), std::string());
            #pragma omp parallel for schedule(dynamic, 64) num_threads(threads)
            for (int64_t i = 0; i < static_cast<int64_t>(holders.size()); i++) {
                seqs[static_cast<size_t>(i)] = std::string(genomes.GetGenome(holders[static_cast<size_t>(i)]).GetGeneOMP(g).Sequence().View());
            }
            found.clear();
            congener_gaps::ScanGene(static_cast<uint32_t>(g), holders, seqs, genus, settings, threads, found, stats);
            for (auto const& [taxid, gap] : found) rows.emplace_back(taxid, static_cast<uint32_t>(g), gap);
        }
        seqs.clear();
        seqs.shrink_to_fit();
        std::vector<double> mins;
        mins.reserve(rows.size());
        for (auto const& row : rows) mins.push_back(std::get<2>(row).Min());
        std::sort(mins.begin(), mins.end());
        auto const table = congener_gaps::Table::FromRows(std::move(rows));
        std::ofstream os(target);
        table.Write(os);
        os.close();
        if (!os) {
            std::cerr << "Writing " << target << " failed" << std::endl;
            exit(8);
        }
        auto quantile = [&mins](double q) {
            return mins.empty() ? 0.0 : mins[std::min(mins.size() - 1, static_cast<size_t>(q * static_cast<double>(mins.size())))];
        };
        std::cout << "Congener gaps: " << table.Copies() << " gene copies of " << table.Species() << " species in " << genera.size()
                  << " genera, " << stats.genes << " genes; " << stats.alignments << " copies aligned (" << stats.far
                  << " beyond " << congener_gaps::kMaxDistance << ", " << stats.unaligned << " without overlap), "
                  << stats.sampled_genera << " genus-gene groups above " << settings.all << " others sampled; the nearest "
                  << "congener's distance: quartiles " << quantile(0.25) << ", " << quantile(0.5) << ", " << quantile(0.75)
                  << ": " << target << std::endl;
        bm.Stop();
        bm.PrintResults();
        return table;
    }

    // strain_alleles.tsv in the database (StrainAllelesBuild.h): every copy of a database species' gene in the full
    // reference from a genome that may give alleles (--allele_genome_share) aligned against the representative's copy,
    // and up to --strain_alleles of each copy's alleles kept, each nearer the representative than the copy's nearest
    // congener's copy (`gaps`, congener_gaps.tsv). One pass over the full reference, as WriteGeneConservation's. No file
    // without an allele (or with --strain_alleles 0); an earlier build's is removed. Returns the table (empty without
    // one): the index passes seed the deep alleles from it (--index_alleles).
    static strain_alleles::Table WriteStrainAlleles(protal::Options const& options, GenomeLoader& genomes, congener_gaps::Table const& gaps) {
        namespace fs = std::filesystem;
        Benchmark bm("Strain alleles");
        bm.Start();
        std::string const target = options.GetStrainAllelesFile();
        std::error_code ec;
        fs::remove(target, ec);
        strain_alleles::Settings settings;
        settings.alleles = options.GetStrainAlleles();
        settings.share = options.GetAlleleGenomeShare();
        if (settings.alleles == 0) {
            std::cout << "Strain alleles: none (--strain_alleles 0)" << std::endl;
            return {};
        }
        std::string const full_reference = options.GetFullSequenceFilePath();
        strain_alleles::Collector collector(settings.candidates);
        strain_alleles::Stats stats;
        std::string failure;
        std::atomic<bool> unnamed{ false };  // a record without its genome's accession, with --allele_genome_share below 1
        // A record's allele (StrainAllelesBuild.h AddRecord); `state` is the thread's.
        struct State {
            WFA2Wrapper2 aligner{ 4, 6, 2, 0 };
            std::string rep, copy, ops, header;
            strain_alleles::Stats stats;
        };
        auto add = [&](State& state, std::string_view record_header, std::string_view sequence) {
            size_t const space = record_header.find_first_of(" \t");
            state.header.assign(record_header.substr(0, space));
            std::string_view accession;
            if (space != std::string_view::npos) {
                accession = record_header.substr(space + 1);
                accession = accession.substr(0, accession.find_first_of(" \t"));
            }
            auto const [taxid, geneid] = KmerUtils::ExtractHeaderInformation(state.header);
            if (taxid > UINT32_MAX || geneid > UINT16_MAX || !genomes.HasGene(static_cast<uint32_t>(taxid), static_cast<uint32_t>(geneid)) ||
                !options.BuildGeneAllowed(geneid)) {
                return;
            }
            if (accession.empty() && settings.share < 1) {
                unnamed = true;
                return;
            }
            auto const rep = genomes.GetGenome(taxid).GetGeneOMP(geneid).Sequence();
            state.rep.assign(rep.View());
            state.copy.assign(sequence);
            auto const gap = gaps.Find(static_cast<uint32_t>(taxid), static_cast<uint32_t>(geneid));
            strain_alleles::AddRecord(collector, state.stats, settings, static_cast<uint32_t>(taxid), static_cast<uint32_t>(geneid),
                                      accession, state.rep, state.copy, gap ? gap->Min() : -1.0, state.aligner, state.ops);
        };
        int const threads = static_cast<int>(std::max<size_t>(options.GetThreads(), 1));
        std::string error;
        auto const frames = FastaFrames::Open(zstd::Resolve(full_reference), error);
        if (!error.empty()) {
            std::cerr << "Cannot read the full reference: " << error << std::endl;
            exit(8);
        }
        if (frames) {
            std::atomic<size_t> next{ 0 };
#pragma omp parallel num_threads(static_cast<int>(std::min<size_t>(static_cast<size_t>(threads), frames->Frames())))
            {
                State state;
                FastaFrames::Reader reader(*frames);
                std::string batch, scratch;
                for (size_t f; (f = next.fetch_add(1)) < frames->Frames();) {
                    reader.Take(f);
                    while (reader.Next(batch, size_t{1} << 20)) {
                        ForEachFastaRecord(batch, scratch, [&](std::string_view record_header, std::string_view sequence) {
                            add(state, record_header, sequence);
                        });
                    }
                    if (!reader.Error().empty()) {
#pragma omp critical(alleles_failure)
                        if (failure.empty()) failure = reader.Error();
                        break;
                    }
                }
#pragma omp critical(alleles_stats)
                stats.Add(state.stats);
            }
        } else {
            auto input = OpenInput(full_reference);
            std::istream& is = input->Stream();
            bool parsed = true;
#pragma omp parallel num_threads(threads)
            {
                State state;
                FastxRecord record;
                SeqReader reader { is };
                while (reader(record)) add(state, record.header, record.sequence);
#pragma omp critical(alleles_stats)
                {
                    parsed = parsed && reader.Success();
                    stats.Add(state.stats);
                }
            }
            if (is.bad() || !parsed) failure = "truncated or corrupt file?";
        }
        if (!failure.empty()) {
            std::cerr << "Cannot read the full reference " << full_reference << " to its end: " << failure << std::endl;
            exit(8);
        }
        if (unnamed) {
            std::cerr << "--allele_genome_share " << settings.share << ": the full reference " << full_reference << " does not name "
                      << "the genome of each record (>taxid_geneid accession, as gtdb_to_protal_db.py writes it since 2026-10-08): "
                      << "convert the release again, or give --allele_genome_share 1" << std::endl;
            exit(8);
        }
        auto const table = strain_alleles::Table::FromRows(collector.Rows(settings.alleles));
        std::ostringstream counts;
        counts << "of " << stats.records << " full-reference copies of the database's species and genes: " << stats.outside_share
               << " of genomes outside --allele_genome_share " << settings.share << ", " << stats.identical
               << " identical to the representative's, " << stats.duplicates << " repeated, " << stats.unaligned
               << " covering less than " << strain_alleles::kMinCover << " of it or more than " << strain_alleles::kMaxDivergence
               << " apart, " << stats.beyond_congener << " as far as the nearest congener's copy or farther; " << stats.offered
               << " alleles offered";
        if (table.Empty()) {
            std::cout << "Strain alleles: none (" << counts.str() << ")" << std::endl;
        } else {
            std::ofstream os(target);
            table.Write(os);
            os.close();
            if (!os) {
                std::cerr << "Writing " << target << " failed" << std::endl;
                exit(8);
            }
            std::cout << "Strain alleles: " << table.Alleles() << " alleles (" << table.Edits() << " edits) of " << table.Copies()
                      << " gene copies of " << table.Species() << " species, up to " << settings.alleles << " each (" << counts.str()
                      << "): " << target << std::endl;
        }
        bm.Stop();
        bm.PrintResults();
        return table;
    }

    // column_weights.tsv in the database (ColumnWeightsBuild.h): per family and gene, how conserved each column is among
    // the family's genus references (up to --column_weights of them) and within each genus, averaged over genera, with
    // the amino acids' conservation per codon and the family's consensus base; and every copy's mapping onto the family
    // reference's columns. Gene by gene, the gene's copies held while its families are compared, the families on all
    // threads. No file with --column_weights 0; an earlier build's is removed.
    static void WriteColumnWeights(protal::Options const& options, GenomeLoader& genomes) {
        namespace fs = std::filesystem;
        Benchmark bm("Column weights");
        bm.Start();
        std::string const target = options.GetColumnWeightsFile();
        std::error_code ec;
        fs::remove(target, ec);
        column_weights::Settings settings;
        settings.genera = options.GetColumnWeights();
        if (settings.genera == 0) {
            std::cout << "Column weights: none (--column_weights 0)" << std::endl;
            return;
        }
        // Where the columns come from (--column_weights_alignment): GTDB's alignments and tree in the database folder
        // (ColumnWeightsMsa.h; auto takes them where they are), else the proteins (or the bases) aligned.
        std::string folder = fs::path(target).parent_path().string();
        if (folder.empty()) folder = ".";
        std::string const msa_folder = folder + "/" + column_weights::msa::kFolder;
        std::string const tree_path = folder + "/" + column_weights::msa::kTreeFile;
        bool const has_msa = fs::is_directory(msa_folder) && fs::exists(zstd::Resolve(tree_path));
        std::string mode = options.ColumnWeightsAlignment();
        if (mode == "auto") mode = has_msa ? "msa" : "aa";
        if (mode == "msa" && !has_msa) {
            std::cerr << "--column_weights_alignment msa: " << folder << " has no " << column_weights::msa::kFolder << "/ and "
                      << column_weights::msa::kTreeFile << " (gtdb_to_protal_db.py writes them from the release's alignments and tree)" << std::endl;
            exit(1);
        }
        bool const from_msa = mode == "msa";
        settings.codons = mode != "nt";
        column_weights::msa::Tree tree;
        if (from_msa) {
            std::string error;
            tree = column_weights::msa::Tree::Read(tree_path, error);
            if (tree.Empty()) {
                std::cerr << "Reading " << tree_path << " failed: " << error << std::endl;
                exit(1);
            }
        }
        column_weights::msa::MsaStats msa_stats;
        size_t genes_without_msa = 0;
        // Each species' genus and family from the taxonomy (as Genera does).
        std::unordered_map<uint32_t, uint32_t> parent;
        std::unordered_map<uint32_t, uint8_t> rank;  // 1 genus, 2 family
        {
            std::ifstream taxonomy(options.GetInternalTaxonomyFile());
            std::string line;
            while (std::getline(taxonomy, line)) {
                std::vector<std::string> fields;
                size_t start = 0;
                for (size_t tab; (tab = line.find('\t', start)) != std::string::npos; start = tab + 1) fields.push_back(line.substr(start, tab - start));
                fields.push_back(line.substr(start));
                uint32_t id = 0, up = 0;
                if (fields.size() < 5 || std::from_chars(fields[0].data(), fields[0].data() + fields[0].size(), id).ec != std::errc() ||
                    std::from_chars(fields[1].data(), fields[1].data() + fields[1].size(), up).ec != std::errc()) continue;
                parent[id] = up;
                rank[id] = fields[4] == "genus" ? 1 : fields[4] == "family" ? 2 : 0;
            }
        }
        auto ancestor = [&](uint32_t taxid, uint8_t wanted) -> uint32_t {
            uint32_t node = taxid;
            for (int depth = 0; depth < 32; depth++) {
                auto const it = parent.find(node);
                if (it == parent.end() || it->second == node) return 0;
                node = it->second;
                if (rank[node] == wanted) return node;
            }
            return 0;
        };
        std::vector<uint32_t> taxids;
        size_t max_gene = 0;
        for (auto const& [taxid, genome] : genomes.GetGenomeMap()) {
            taxids.push_back(static_cast<uint32_t>(taxid));
            max_gene = std::max(max_gene, genome.GetGeneList().size());
        }
        std::sort(taxids.begin(), taxids.end());
        std::vector<uint32_t> genus_all(taxids.size()), family_all(taxids.size());
        for (size_t i = 0; i < taxids.size(); i++) {
            genus_all[i] = ancestor(taxids[i], 1);
            family_all[i] = ancestor(taxids[i], 2);
        }
        int const threads = static_cast<int>(std::max<size_t>(options.GetThreads(), 1));
        column_weights::Stats stats;
        std::vector<column_weights::Family> families;
        std::vector<std::tuple<uint32_t, uint32_t, uint32_t, std::vector<column_weights::Run>>> copies;
        std::vector<uint32_t> holders, genus_of, family_of;
        std::vector<std::string> seqs;
        for (size_t g = 1; g <= max_gene; g++) {
            if (!options.BuildGeneAllowed(g)) continue;
            holders.clear();
            genus_of.clear();
            family_of.clear();
            for (size_t i = 0; i < taxids.size(); i++) {
                auto const& list = genomes.GetGenome(taxids[i]).GetGeneList();
                if (g <= list.size() && list[g - 1].IsSet() && family_all[i] != 0 && genus_all[i] != 0) {
                    holders.push_back(taxids[i]);
                    genus_of.push_back(genus_all[i]);
                    family_of.push_back(family_all[i]);
                }
            }
            if (holders.empty()) continue;
            seqs.assign(holders.size(), std::string());
            #pragma omp parallel for schedule(dynamic, 64) num_threads(threads)
            for (int64_t i = 0; i < static_cast<int64_t>(holders.size()); i++) {
                seqs[static_cast<size_t>(i)] = std::string(genomes.GetGenome(holders[static_cast<size_t>(i)]).GetGeneOMP(g).Sequence().View());
            }
            if (from_msa) {
                std::unordered_set<uint32_t> const wanted(holders.begin(), holders.end());
                size_t width = 0;
                auto const rows = column_weights::msa::ReadRows(msa_folder, static_cast<uint32_t>(g), wanted, width);
                if (rows.empty()) {
                    genes_without_msa++;
                    continue;
                }
                column_weights::msa::ScanGeneFromMsa(static_cast<uint32_t>(g), holders, seqs, genus_of, family_of, rows, width, tree,
                                                     threads, msa_stats, families, copies);
            } else {
                column_weights::ScanGene(static_cast<uint32_t>(g), holders, seqs, genus_of, family_of, settings, threads, stats, families, copies);
            }
        }
        seqs.clear();
        seqs.shrink_to_fit();
        if (families.empty()) {
            std::cout << "Column weights: none (no species with a family and a genus in the taxonomy)" << std::endl;
            return;
        }
        // How the columns' codes fall (column_weights::CodeSummary): whether the conserved columns the rate features
        // count (within, aa >= kConservedCode) and the full differences of --weighted_site_shift (among >=
        // strain_alleles::kDiscountNone) occur at all.
        auto const codes = column_weights::CodeSummary(families);
        auto const table = column_weights::Table::FromRows(std::move(families), std::move(copies));
        std::ofstream os(target);
        table.Write(os);
        os.close();
        if (!os) {
            std::cerr << "Writing " << target << " failed" << std::endl;
            exit(8);
        }
        if (from_msa) {
            auto const share = [](size_t part, size_t whole) { return whole ? static_cast<double>(part) / static_cast<double>(whole) : 0.0; };
            std::cout << "Column weights: " << table.Families() << " family-gene rows of " << table.ColumnCount() << " columns, "
                      << table.Copies() << " gene copies of " << table.Species() << " species mapped, from GTDB's alignments and tree ("
                      << msa_folder << ", " << tree_path << "): " << target << std::endl;
            std::cout << "Column weights from GTDB's alignments: " << msa_stats.copies << " copies, " << msa_stats.without_row
                      << " without a row, " << msa_stats.badly_placed << " placed below " << column_weights::msa::kMinPlaced
                      << " of their row; of the residues of those used " << share(msa_stats.placed, msa_stats.residues)
                      << " on a column (the rest insertions the alignment trimmed); " << msa_stats.genus_ancestors << " genus ancestors ("
                      << msa_stats.genus_trees << " by parsimony on the tree, the others a genus's one species), the families' ancestral bases "
                      << share(msa_stats.family_ambiguous, msa_stats.family_columns) << " ambiguous at the root, "
                      << msa_stats.polarised << " columns polarised; " << genes_without_msa << " genes without an alignment" << std::endl;
            std::cout << "Column weights codes (half nats) of the " << codes.columns << " columns: " << codes.Line(strain_alleles::kDiscountNone) << std::endl;
            bm.Stop();
            bm.PrintResults();
            return;
        }
        std::cout << "Column weights: " << table.Families() << " family-gene rows of " << table.ColumnCount() << " columns, "
                  << table.Copies() << " gene copies of " << table.Species() << " species mapped (" << stats.single_genus
                  << " rows of a family with one genus); " << stats.alignments << " alignments, " << stats.unaligned
                  << " failed; up to " << settings.genera << " genera vote per family: " << target << std::endl;
        auto const failures = [&](size_t kind) {
            double const limit = settings.MaxDivergence(kind);
            auto const& f = stats.failed[kind];
            using column_weights::MapFailure;
            std::ostringstream s;
            s << stats.made[kind] << " aligned, " << stats.made[kind] - std::accumulate(f.begin(), f.end(), size_t{ 0 }) << " kept; failed: "
              << f[static_cast<size_t>(MapFailure::kDivergence)] << " beyond " << limit << " (" << stats.far[kind] << " beyond "
              << column_weights::kFarFactor * limit << "), " << f[static_cast<size_t>(MapFailure::kGaveUp)]
              << " given up by the aligner (score beyond the bound: far or gappy), " << f[static_cast<size_t>(MapFailure::kOverlap)]
              << " overlapping less than " << column_weights::kMinOverlap << " of the shorter, " << f[static_cast<size_t>(MapFailure::kSequence)]
              << (settings.codons ? " empty, too long or too short" : " empty or too long") << "; the divergence of those compared: quartiles "
              << stats.Quantile(kind, 0.25) << ", " << stats.Quantile(kind, 0.5) << ", " << stats.Quantile(kind, 0.75) << " (bins of 0.05)";
            return s.str();
        };
        std::cout << "Column weights alignments (" << (settings.codons ? "amino acids, translated in frame" : "bases")
                  << "): genus references against their family's: " << failures(column_weights::kAmong)
                  << "; species against their genus's: " << failures(column_weights::kWithin) << "; " << stats.skipped
                  << " species' copies not aligned (their genus reference failed)" << std::endl;
        std::cout << "Column weights codes (half nats) of the " << codes.columns << " columns: " << codes.Line(strain_alleles::kDiscountNone) << std::endl;
        bm.Stop();
        bm.PrintResults();
    }

    // suspect_copies.tsv in the database and gene_incongruence.tsv beside it (gene_incongruence::Scan): every
    // species' copy of each gene sketched and compared with the other species' copies; a copy within
    // --suspect_copy_distance of another genus's copy, and farther from its congeners' or without one, is suspect
    // (contamination or a transferred gene), and a run leaves its records out (--keep_suspect_copies). The lineages
    // come from internal_taxonomy.dmp. --suspect_copy_distance 0: no scan, and no file.
    static void WriteSuspectCopies(protal::Options const& options, GenomeLoader& genomes) {
        namespace fs = std::filesystem;
        std::string const target = options.GetSuspectCopiesFile();
        std::string const report = (fs::path(target).parent_path() / gene_incongruence::kReportFileName).string();
        std::error_code ec;
        fs::remove(target, ec);
        fs::remove(report, ec);
        double const threshold = options.GetSuspectCopyDistance();
        if (threshold <= 0) {
            std::cout << "Suspect copies: not looked for (--suspect_copy_distance 0)" << std::endl;
            return;
        }
        Benchmark bm("Suspect copies");
        bm.Start();
        std::unordered_map<uint32_t, gene_incongruence::Lineage> lineages;
        if (std::string const error = gene_incongruence::LineagesFromTaxonomy(options.GetInternalTaxonomyFile(), lineages); !error.empty()) {
            std::cerr << "Suspect copies: cannot read the taxonomy " << options.GetInternalTaxonomyFile() << ": " << error << std::endl;
            exit(8);
        }
        std::vector<uint32_t> taxids;
        size_t max_gene = 0;
        for (auto const& [taxid, genome] : genomes.GetGenomeMap()) {
            taxids.push_back(static_cast<uint32_t>(taxid));
            max_gene = std::max(max_gene, genome.GetGeneList().size());
        }
        std::sort(taxids.begin(), taxids.end());
        int const threads = static_cast<int>(std::max<size_t>(options.GetThreads(), 1));
        // Gene by gene: its copies sketched (by taxid) and compared on all threads, so that only one gene's sketches are
        // held at a time (every gene's, ~13 GB at r226 with 64 threads, before) and no thread waits for the largest genes.
        gene_incongruence::Result result;
        gene_incongruence::GeneCopies copies;
        std::vector<uint32_t> holders;
        for (size_t g = 1; g <= max_gene; g++) {
            if (!options.BuildGeneAllowed(g)) continue;  // the subset's genes only
            holders.clear();
            for (uint32_t const taxid : taxids) {
                auto const& list = genomes.GetGenome(taxid).GetGeneList();
                if (g <= list.size() && list[g - 1].IsSet()) holders.push_back(taxid);
            }
            if (holders.empty()) continue;
            copies.Resize(holders.size(), gene_incongruence::kSketchSize);
            #pragma omp parallel for schedule(dynamic, 64) num_threads(threads)
            for (int64_t i = 0; i < static_cast<int64_t>(holders.size()); i++) {
                auto const seq = genomes.GetGenome(holders[static_cast<size_t>(i)]).GetGeneOMP(g).Sequence();
                copies.Set(static_cast<size_t>(i), holders[static_cast<size_t>(i)],
                           gene_conservation::BottomSketch(seq.View(), gene_incongruence::kSketchSize));
            }
            gene_incongruence::ScanGene(static_cast<uint32_t>(g), copies, lineages, threshold, threads, result);
        }
        gene_incongruence::Finish(result);
        gene_incongruence::Table table;
        for (auto const& suspect : result.suspects) table.Add(suspect);
        {
            std::ofstream os(report);
            gene_incongruence::WriteReport(os, result);
            os.close();
            if (!os) {
                std::cerr << "Writing " << report << " failed" << std::endl;
                exit(8);
            }
        }
        if (!table.Empty()) {
            std::ofstream os(target);
            table.Write(os);
            os.close();
            if (!os) {
                std::cerr << "Writing " << target << " failed" << std::endl;
                exit(8);
            }
        }
        std::cout << "Suspect copies: " << table.Size() << " of " << result.copies << " gene copies ("
                  << std::fixed << std::setprecision(2) << 100.0 * static_cast<double>(table.Size()) / std::max<double>(1, static_cast<double>(result.copies))
                  << std::defaultfloat << "%) of " << table.Species() << " species within " << threshold
                  << " of another genus's copy and farther from their congeners'"
                  << (table.Empty() ? std::string(": none") : ": " + target) << "; "
                  << result.pairs.size() << " near pairs across genera (within " << gene_incongruence::kReportDistance << ") of "
                  << result.genes << " genes, " << result.candidates << " candidate pairs, " << result.compared
                  << " of them compared in full (the others cannot be within " << gene_incongruence::kReportDistance << "): " << report
                  << std::setprecision(6) << std::endl;  // the default again: the index passes print after this
        bm.Stop();
        bm.PrintResults();
    }

    // A batch of records for the passes (PartitionedPass, UniquenessPass): the text of whole FASTA records (TextSource,
    // FramesSource), or the genes [begin, end) of the reference's order (GeneSource). Empty: nothing in this slot.
    struct BatchInput {
        std::string text;
        size_t begin = 0, end = 0;
        bool Empty() const { return text.empty() && begin == end; }
    };

    // The records of a FASTA file read in one stream (FastaBatches): one thread reads each round's batches.
    class TextSource {
    public:
        TextSource(std::istream& is, size_t batch_bytes) : m_reader(is), m_batch_bytes(batch_bytes) {}

        // Called by every thread of the pass's team: fills inputs, returns how many it filled (0: the end).
        size_t Fill(std::vector<BatchInput>& inputs) {
#pragma omp single
            {
                m_filled = 0;
                while (m_filled < inputs.size() && m_reader.Next(inputs[m_filled].text, m_batch_bytes)) m_filled++;
            }
            return m_filled;
        }

        // f(taxid, gene id, sequence) for each record of the batch.
        template<typename F>
        void ForEach(BatchInput const& input, std::string& scratch, F&& f) {
            thread_local std::string header;
            ForEachFastaRecord(input.text, scratch, [&](std::string_view record_header, std::string_view sequence) {
                header.assign(record_header);
                auto const [taxid, geneid] = KmerUtils::ExtractHeaderInformation(header);
                f(taxid, geneid, sequence);
            });
        }

        std::string const& Error() const { return m_reader.Error(); }

    private:
        FastaBatches m_reader;
        size_t m_batch_bytes;
        size_t m_filled = 0;
    };

    // Threads that decompress the full reference's frames at once in the uniqueness check (FramesSource; at most half of
    // -t): each holds a frame's window, 128 MB for the converter's, beside the index.
    inline constexpr size_t kFrameDecoders = 16;

    // The records of a FASTA file in frames that each begin a record (FastaFrames: full_reference.fna.zst as the
    // converter writes it), decompressed by `decoders` threads at once, each taking the next frame when its own ends: the
    // order of the batches is not the file's, which the uniqueness check does not need. Each decoder holds a frame's
    // window (128 MB for the converter's frames).
    class FramesSource {
    public:
        FramesSource(FastaFrames const& frames, size_t decoders, size_t batch_bytes) : m_frames(frames), m_batch_bytes(batch_bytes) {
            decoders = std::clamp<size_t>(decoders, 1, frames.Frames());
            for (size_t d = 0; d < decoders; d++) m_decoders.push_back(std::make_unique<Decoder>(frames));
        }

        size_t Decoders() const { return m_decoders.size(); }

        // Called by every thread of the pass's team: each decoder fills its share of inputs (a slot it leaves empty
        // when its frames are done). Returns inputs.size() while any decoder had records, then 0.
        size_t Fill(std::vector<BatchInput>& inputs) {
            size_t const decoders = m_decoders.size();
            size_t const share = std::max<size_t>(1, inputs.size() / decoders);
#pragma omp single
            {
                m_any.store(false, std::memory_order_relaxed);
                for (auto& input : inputs) input.text.clear();
            }
#pragma omp for schedule(static, 1)
            for (size_t d = 0; d < decoders; d++) {
                Decoder& decoder = *m_decoders[d];
                for (size_t slot = d * share; slot < std::min(inputs.size(), (d + 1) * share); slot++) {
                    while (!decoder.done) {
                        if (decoder.reader.Next(inputs[slot].text, m_batch_bytes)) break;
                        if (!decoder.reader.Error().empty()) {
                            decoder.done = true;
                            break;
                        }
                        size_t const f = m_next.fetch_add(1);
                        if (f >= m_frames.Frames()) decoder.done = true;
                        else decoder.reader.Take(f);
                    }
                    if (!inputs[slot].text.empty()) m_any.store(true, std::memory_order_relaxed);
                }
            }
            // The implicit barrier of the loop above: every decoder has filled its slots.
            return m_any.load(std::memory_order_relaxed) ? inputs.size() : 0;
        }

        template<typename F>
        void ForEach(BatchInput const& input, std::string& scratch, F&& f) {
            thread_local std::string header;
            ForEachFastaRecord(input.text, scratch, [&](std::string_view record_header, std::string_view sequence) {
                header.assign(record_header);
                auto const [taxid, geneid] = KmerUtils::ExtractHeaderInformation(header);
                f(taxid, geneid, sequence);
            });
        }

        std::string Error() const {
            for (auto const& decoder : m_decoders) {
                if (!decoder->reader.Error().empty()) return decoder->reader.Error();
            }
            return "";
        }

    private:
        struct Decoder {
            explicit Decoder(FastaFrames const& frames) : reader(frames) {}
            FastaFrames::Reader reader;
            bool done = false;
        };
        FastaFrames const& m_frames;
        size_t m_batch_bytes;
        std::vector<std::unique_ptr<Decoder>> m_decoders;
        std::atomic<size_t> m_next{ 0 };
        std::atomic<bool> m_any{ false };
    };

    // reference.fna's records from the preloaded genes (GenomeLoader::ReferenceText), in the file's order: no file read,
    // and each batch decoded by the thread that takes it.
    class GeneSource {
    public:
        GeneSource(GenomeLoader& genomes, size_t batch_bytes) : m_genomes(genomes), m_batch_bytes(batch_bytes) {}

        size_t Fill(std::vector<BatchInput>& inputs) {
#pragma omp single
            {
                auto const& order = m_genomes.FileOrder();
                m_filled = 0;
                while (m_filled < inputs.size() && m_next < order.size()) {
                    BatchInput& input = inputs[m_filled++];
                    input.text.clear();
                    input.begin = m_next;
                    size_t bytes = 0;
                    while (m_next < order.size() && (m_next == input.begin || bytes < m_batch_bytes)) {
                        bytes += m_genomes.GeneLength(order[m_next].first, order[m_next].second) + 16;  // and a header's
                        m_next++;
                    }
                    input.end = m_next;
                }
            }
            return m_filled;
        }

        template<typename F>
        void ForEach(BatchInput const& input, std::string& scratch, F&& f) {
            auto const& order = m_genomes.FileOrder();
            for (size_t g = input.begin; g < input.end; g++) {
                m_genomes.ReferenceText(g, scratch);
                f(static_cast<size_t>(order[g].first), static_cast<size_t>(order[g].second), std::string_view(scratch));
            }
        }

        std::string Error() const { return ""; }

    private:
        GenomeLoader& m_genomes;
        size_t m_batch_bytes;
        size_t m_next = 0, m_filled = 0;
    };

    // A pass over the reference in `threads` threads that updates the index as one thread would
    // (docs/claude/2026-09-29-index-build-parallel.md). The key space is cut into `ranges` ranges of
    // whole control blocks. Each round, the source gives batches of whole records, in the reference's
    // order (TextSource: one thread reads them; GeneSource: the preloaded genes); the threads parse
    // them, extract their items (extract: a record's k-mers, in order) and group each batch's items by
    // range (range_of), keeping their order within a range; then each range is applied (apply) by one
    // thread, batch by batch in reference order. A key's updates thus come in the serial order, and
    // none needs a lock: no two threads touch the same block.
    template<typename Item, typename Source, typename KmerHandler, typename Extract, typename RangeOf, typename Apply>
    void PartitionedPass(Source& source, KmerHandler const& handler_global, int threads, size_t ranges, Extract&& extract,
                         RangeOf&& range_of, Apply&& apply, Statistics& statistics) {
        struct Batch {
            std::vector<Item> items;        // grouped by range
            std::vector<uint32_t> offsets;  // items [offsets[r], offsets[r + 1]) are range r's
        };
        threads = std::max(threads, 1);
        size_t const per_round = 4 * static_cast<size_t>(threads);
        std::vector<BatchInput> inputs(per_round);
        std::vector<Batch> batches(per_round);

#pragma omp parallel num_threads(threads)
        {
            KmerHandler handler(handler_global);
            Statistics local;
            std::vector<Item> items;  // a batch's items in record order
            std::vector<uint32_t> cursor(ranges);
            std::string scratch;
            KmerList kmers;
            while (true) {
                size_t const filled = source.Fill(inputs);
                if (filled == 0) break;

#pragma omp for schedule(dynamic, 1)
                for (size_t b = 0; b < filled; b++) {
                    Batch& batch = batches[b];
                    items.clear();
                    source.ForEach(inputs[b], scratch, [&](size_t taxid, size_t geneid, std::string_view sequence) {
                        extract(handler, taxid, geneid, sequence, kmers, items, local);
                    });
                    batch.offsets.assign(ranges + 1, 0);
                    for (Item const& item : items) batch.offsets[range_of(item) + 1]++;
                    for (size_t r = 0; r < ranges; r++) batch.offsets[r + 1] += batch.offsets[r];
                    std::copy(batch.offsets.begin(), batch.offsets.end() - 1, cursor.begin());
                    batch.items.resize(items.size());
                    for (Item const& item : items) batch.items[cursor[range_of(item)]++] = item;
                }

#pragma omp for schedule(dynamic, 1)
                for (size_t r = 0; r < ranges; r++) {
                    for (size_t b = 0; b < filled; b++) {
                        Batch const& batch = batches[b];
                        for (uint32_t i = batch.offsets[r]; i < batch.offsets[r + 1]; i++) apply(batch.items[i]);
                    }
                }
            }
#pragma omp critical(statistics)
            statistics.Join(local);
        }
        if (!source.Error().empty()) {
            std::cerr << "Cannot read the reference: " << source.Error() << std::endl;
            exit(8);
        }
    }

    // What the uniqueness check did, for its log line.
    struct UniquenessCounts {
        size_t kmers = 0;     // the full reference's k-mers looked up
        size_t distinct = 0;  // distinct k-mers per key range and round, each looked up once
        size_t cells = 0;     // flex cells compared
        size_t shared = 0;    // whole k-mers found under another taxon or more than once
        size_t singles = 0;   // single entries read back from their genes
    };

    // The uniqueness check (Run) by key range, from any source of the full reference's records (TextSource, FramesSource):
    // each round's k-mers are grouped by range as in PartitionedPass, then each range's are sorted by core and flex part,
    // and each distinct whole k-mer gets the taxa that hold it (the lowest and highest taxid suffice); each core's cells
    // are then scanned once for all its k-mers of the round, instead of once per k-mer (5.3 times per value at GTDB r226,
    // ~192 cells each: 2.8e12 compared). The rule is the one a k-mer at a time applied (the serial check below, kept
    // for --serial_index_passes): an exact entry of another taxon, or several exact entries, make every exact entry
    // non-unique; a core's single entry (no flex cells) is compared with its gene. Clears commute, so the flags, and the
    // index, are the same.
    template<typename Source, typename KmerHandler>
    UniquenessCounts UniquenessPass(Source& source, KmerHandler const& handler_global, int threads, size_t ranges, int range_shift,
                                    Seedmap& map, GenomeLoader& genomes, protal::Options const& options, size_t kmer_length) {
        struct Item { uint64_t key; uint32_t taxid; };  // key: core << 32 | flex part (a whole k-mer), its record's taxid
        struct Batch {
            std::vector<Item> items;
            std::vector<uint32_t> offsets;
        };
        threads = std::max(threads, 1);
        size_t const per_round = 4 * static_cast<size_t>(threads);
        std::vector<BatchInput> inputs(per_round);
        std::vector<Batch> batches(per_round);
        UniquenessCounts counts;

#pragma omp parallel num_threads(threads)
        {
            KmerHandler handler(handler_global);
            UniquenessCounts local;
            std::vector<Item> items, group;
            std::vector<uint32_t> cursor(ranges);
            std::string scratch;
            KmerList kmers;
            struct Query { uint32_t flex; uint32_t min_taxid, max_taxid; uint32_t matches = 0, first = 0; };
            std::vector<Query> queries;
            std::vector<std::pair<uint32_t, uint32_t>> pairs;  // the further entries of queries with several: (query, entry)
            Seedmap::PackedBlock block;
            while (true) {
                size_t const filled = source.Fill(inputs);
                if (filled == 0) break;

#pragma omp for schedule(dynamic, 1)
                for (size_t b = 0; b < filled; b++) {
                    Batch& batch = batches[b];
                    items.clear();
                    if (!inputs[b].Empty()) {
                        source.ForEach(inputs[b], scratch, [&](size_t taxid, size_t geneid, std::string_view sequence) {
                            // With --build_gene_subset the database's genes are the subset: a copy of another gene is
                            // not among the sequences a read could come from, so it does not make k-mers non-unique.
                            if (options.HasBuildGeneSubset() && !options.BuildGeneAllowed(geneid)) return;
                            kmers.clear();
                            handler(sequence, kmers);
                            DropAmbiguousKmers(sequence, kmer_length, kmers);
                            for (auto const& pair : kmers) {
                                items.push_back({ map.MainKey(pair.first) << 32 | map.FlexKey(pair.first), static_cast<uint32_t>(taxid) });
                            }
                        });
                    }
                    local.kmers += items.size();
                    batch.offsets.assign(ranges + 1, 0);
                    for (Item const& item : items) batch.offsets[(item.key >> (32 + range_shift)) + 1]++;
                    for (size_t r = 0; r < ranges; r++) batch.offsets[r + 1] += batch.offsets[r];
                    std::copy(batch.offsets.begin(), batch.offsets.end() - 1, cursor.begin());
                    batch.items.resize(items.size());
                    for (Item const& item : items) batch.items[cursor[item.key >> (32 + range_shift)]++] = item;
                }

#pragma omp for schedule(dynamic, 1)
                for (size_t r = 0; r < ranges; r++) {
                    group.clear();
                    for (size_t b = 0; b < filled; b++) {
                        Batch const& batch = batches[b];
                        group.insert(group.end(), batch.items.begin() + batch.offsets[r], batch.items.begin() + batch.offsets[r + 1]);
                    }
                    if (group.empty()) continue;
                    std::sort(group.begin(), group.end(), [](Item const& a, Item const& b) { return a.key < b.key; });
                    for (size_t i = 0; i < group.size();) {
                        // The distinct k-mers of one core, each with the lowest and highest taxid holding it.
                        uint64_t const core = group[i].key >> 32;
                        queries.clear();
                        while (i < group.size() && (group[i].key >> 32) == core) {
                            uint64_t const key = group[i].key;
                            Query q{ static_cast<uint32_t>(key), group[i].taxid, group[i].taxid };
                            for (; i < group.size() && group[i].key == key; i++) {
                                q.min_taxid = std::min(q.min_taxid, group[i].taxid);
                                q.max_taxid = std::max(q.max_taxid, group[i].taxid);
                            }
                            queries.push_back(q);
                        }
                        local.distinct += queries.size();
                        uint64_t start = 0, slots = 0;
                        if (!map.Locate(core, start, slots)) continue;
                        map.BlockOfSlots(start, slots, block);
                        if (block.flex == nullptr) {
                            // No flex cells: a core of one value (or below the flex threshold, which has none to compare).
                            if (block.size != 1) continue;
                            size_t taxid = 0, geneid = 0, genepos = 0;
                            ValueEntry entry;
                            entry.value = map.EntryValue(block, 0);
                            entry.Get(taxid, geneid, genepos);
                            uint64_t indexed = 0;
                            bool read = false, have = false;
                            for (Query const& q : queries) {
                                if (q.min_taxid == taxid && q.max_taxid == taxid) continue;  // only its own taxon
                                if (!read) {
                                    read = true;
                                    local.singles++;
                                    have = IndexedKmer(map, genomes, taxid, geneid, genepos, indexed);
                                }
                                if (have && map.MainKey(indexed) == core && map.FlexKey(indexed) == q.flex) {
                                    map.ClearUniqueFlag(block, 0);
                                    local.shared++;
                                    break;
                                }
                            }
                            continue;
                        }
                        // Each cell once: the queries are sorted by their flex parts.
                        pairs.clear();
                        for (uint32_t c = 0; c < block.size; c++) {
                            uint32_t const flex = Seedmap::FlexCell(block, c);
                            auto const it = std::lower_bound(queries.begin(), queries.end(), flex,
                                                             [](Query const& q, uint32_t f) { return q.flex < f; });
                            if (it == queries.end() || it->flex != flex) continue;
                            if (it->matches++ == 0) it->first = c;
                            else pairs.emplace_back(static_cast<uint32_t>(it - queries.begin()), c);
                        }
                        local.cells += block.size;
                        for (size_t q = 0; q < queries.size(); q++) {
                            Query const& query = queries[q];
                            if (query.matches == 0 || query.matches > map.max_key_multiplicity) continue;  // as GetExact
                            if (query.matches == 1) {
                                size_t taxid = 0, geneid = 0, genepos = 0;
                                ValueEntry entry;
                                entry.value = map.EntryValue(block, query.first);
                                entry.Get(taxid, geneid, genepos);
                                if (query.min_taxid == taxid && query.max_taxid == taxid) continue;
                            }
                            local.shared++;
                            map.ClearUniqueFlag(block, query.first);
                        }
                        for (auto const& [q, c] : pairs) {
                            if (queries[q].matches <= map.max_key_multiplicity) map.ClearUniqueFlag(block, c);
                        }
                    }
                }
            }
#pragma omp critical(statistics)
            {
                counts.kmers += local.kmers;
                counts.distinct += local.distinct;
                counts.cells += local.cells;
                counts.shared += local.shared;
                counts.singles += local.singles;
            }
        }
        if (!source.Error().empty()) {
            std::cerr << "Cannot read the full reference: " << source.Error() << std::endl;
            exit(8);
        }
        return counts;
    }

    // Key ranges for PartitionedPass: at least 64 per thread for balance, and at least 64 control
    // blocks each.
    inline int IndexRangeBits(int threads, size_t main_bits) {
        int bits = 6;
        while (bits < 16 && (size_t{1} << bits) < 64 * static_cast<size_t>(std::max(threads, 1))) bits++;
        return std::min<int>(bits, static_cast<int>(main_bits) - 3 - 6);
    }

    // reference_is_mapped: --reference is the database's reference.fna and holds each gene of reference.map once
    // (CheckReferenceAgainstMap), so that the passes may take its records from the preloaded genes.
    template<typename KmerHandler, typename KmerPutter, DebugLevel debug>
    static Statistics Run(protal::Options const& options, KmerPutter& putter, KmerHandler& kmer_handler_global, GenomeLoader& genomes,
                          bool reference_is_mapped) {

        // The tables made from the genes alone first (conservation factors, congeners, suspect copies) and the checks
        // of the gene neighbours and positions: none needs the index, so they run before it is allocated. The
        // suspect-copy scan held ~20 GB at r226 scale with 64 threads, which on top of the index made the build's
        // peak (docs/claude/2026-10-05-build-memory); a bad table also stops the build before its passes.
        auto const conservation = WriteGeneConservation(options, genomes, options.GetFullSequenceFilePath(),
                                                        options.GetGeneConservationFile());
        WriteGeneCongeners(options, genomes, conservation.table);
        WriteSuspectCopies(options, genomes);
        WriteSpeciesNeighbours(options, genomes, options.GetSpeciesNeighboursFile());
        auto const gaps = WriteCongenerGaps(options, genomes);
        auto const alleles = WriteStrainAlleles(options, genomes, gaps);
        WriteColumnWeights(options, genomes);
        CheckGeneNeighbours(options, genomes);
        CheckGenePositions(options, genomes);
        ReleaseFreeMemory();
        PrintMemory("the gene tables");

        // Shared
        auto input = OpenInput(options.GetSequenceFilePath());
        std::istream& is = input->Stream();
        size_t dummy = 0;
        int read_count = 0;

        // Set Thread Num
        omp_set_num_threads(options.GetThreads());

        Statistics statistics;

        size_t main_k = putter.GetMap().m_exact_k;
        size_t main_k_bits = putter.GetMap().m_main_bits;
        size_t flex_k = putter.GetMap().m_flex_k;
        size_t flex_k_bits = putter.GetMap().m_flex_k_bits;
        size_t const kmer_length = main_k + flex_k;  // a window's length

        // The deep strain alleles' seeds (--index_alleles; StrainAllelesBuild.h AlleleSeeds): the k-mers of a copy's
        // alleles at least index_alleles from the representative that the copy lacks go into both passes under the
        // species, as non-unique entries (a strain's k-mer is not the representative's: it neither counts as a unique
        // hit nor in the copy's totals in unique_kmers.tsv, which allele_seeds keeps per gene row). Counted for the log.
        double const index_alleles = options.IndexAlleles();
        bool const seed_alleles = !alleles.Empty() && index_alleles < 1;
        auto const gene_rows = GeneRowsOf(genomes);
        std::vector<uint32_t> allele_seeds(seed_alleles ? gene_rows.Size() : 0, 0);
        std::atomic<size_t> allele_copies{ 0 }, alleles_seeded{ 0 }, allele_kmers{ 0 };
        auto note_allele_seed = [&](uint64_t value) {
            ValueEntry entry;
            entry.value = value;
            auto const [t, g, pos] = entry.Get();
            if (t + 1 < gene_rows.first_row.size() && g > 0) {
                std::atomic_ref<uint32_t>(allele_seeds[gene_rows.first_row[t] + g - 1]).fetch_add(1, std::memory_order_relaxed);
            }
            allele_kmers.fetch_add(1, std::memory_order_relaxed);
        };


        // The two passes over the reference and the value pointers run in -t threads, the passes by
        // key range (PartitionedPass), with the index one thread gives; --serial_index_passes runs
        // them on one thread as before.
        int const threads = static_cast<int>(std::max<size_t>(options.GetThreads(), 1));
        bool const serial = options.SerialIndexPasses();
        Seedmap& map = putter.GetMap();
        int const range_bits = IndexRangeBits(threads, map.m_main_bits);
        size_t const ranges = size_t{1} << range_bits;
        int const range_shift = static_cast<int>(map.m_main_bits) - range_bits;
        if (serial) omp_set_num_threads(1);
        // The passes' records: the preloaded genes in reference.fna's order (GenomeLoader::ReferenceText; the file is not
        // read again, 15 GB at r226 that the page cache often no longer held), when --reference is that file and every
        // gene is preloaded; else the file. The same records either way, so the same index.
        bool const from_genes = !serial && reference_is_mapped && !genomes.FileOrder().empty();
        std::cout << "Index passes: the reference's records " << (from_genes ? "from the preloaded genes (" + std::to_string(genomes.RawBases()) +
                     " bases of them as the reference has them, not A, C, G or T)" : "read from " + options.GetSequenceFilePath()) << std::endl;
        auto with_source = [&](auto&& pass) {
            if (from_genes) {
                GeneSource source(genomes, options.IndexBatchBytes());
                pass(source);
            } else {
                TextSource source(is, options.IndexBatchBytes());
                pass(source);
            }
        };

        // Each phase is timed in the log ("... took"): where a build at GTDB scale spends its time.
        std::cout << "Run Build" << std::endl;
        Benchmark bm_pass1("Pass 1 (count the k-mers)");
        bm_pass1.Start();
        if (!serial) {
            // Items: the k-mers' main keys, counted up by the thread that owns their range.
            std::cout << "Build: iterate records" << std::endl;
            Statistics pass;
            with_source([&](auto& source) {
                PartitionedPass<uint32_t>(source, kmer_handler_global, threads, ranges,
                    [&](KmerHandler& handler, size_t taxonomic_id, size_t gene_id, std::string_view sequence, KmerList& kmers,
                        std::vector<uint32_t>& items, Statistics& stats) {
                        if (options.HasBuildGeneSubset() && !options.BuildGeneAllowed(gene_id)) return;
                        stats.reads++;
                        kmers.clear();
                        handler(sequence, kmers);
                        DropAmbiguousKmers(sequence, kmer_length, kmers);
                        for (auto const& pair : kmers) items.push_back(static_cast<uint32_t>(map.MainKey(pair.first)));
                        if constexpr(KmerStatisticsConcept<KmerHandler>) {
                            stats.kmers_total += handler.TotalKmers();
                            stats.kmers_accepted += kmers.size();
                        }
                        if (seed_alleles) {
                            static thread_local strain_alleles::SeedScratch scratch;
                            static thread_local KmerList seeds;
                            seeds.clear();
                            strain_alleles::AlleleSeeds(alleles, index_alleles, static_cast<uint32_t>(taxonomic_id), static_cast<uint32_t>(gene_id),
                                                        sequence, kmers, kmer_length, handler, scratch, seeds);
                            for (auto const& pair : seeds) items.push_back(static_cast<uint32_t>(map.MainKey(pair.first)));
                        }
                    },
                    [range_shift](uint32_t main_key) { return static_cast<size_t>(main_key >> range_shift); },
                    [&map](uint32_t main_key) { map.CountUpKey(main_key); },
                    pass);
            });
            std::cout << "minimizers: " << pass.kmers_accepted << std::endl;
            statistics.Join(pass);
        } else {
#pragma omp parallel default(none) shared(std::cout, options, is, dummy, read_count, kmer_handler_global, statistics, putter, main_k_bits, flex_k_bits, kmer_length, seed_alleles, index_alleles, alleles)
                {
                    // Private variables
                    FastxRecord record;

                    // Extract variables from kmi_global
                    KmerHandler kmer_handler(kmer_handler_global);
                    SeqReader reader{ is };
                    Statistics thread_statistics;
                    thread_statistics.thread_num = omp_get_thread_num();

                    KmerList kmers;
                    strain_alleles::SeedScratch scratch;  // the deep alleles' seeds (--index_alleles)
                    KmerList seeds;

                    std::cout << "Build: iterate records" << std::endl;

                    while (reader(record)) {
                        thread_statistics.reads++;

                        if (options.HasBuildGeneSubset()) {
                            auto [taxonomic_id, gene_id] = KmerUtils::ExtractHeaderInformation(record.header);
                            (void)taxonomic_id;
                            if (!options.BuildGeneAllowed(gene_id)) {
                                continue;
                            }
                        }

                        // Retrieve kmers
                        kmers.clear();
                        kmer_handler(std::string_view(record.sequence), kmers);
                        DropAmbiguousKmers(record.sequence, kmer_length, kmers);
                        for (auto pair : kmers) {
                            putter.FirstPut(pair.first);
                        }
                        if (seed_alleles) {
                            auto [taxonomic_id, gene_id] = KmerUtils::ExtractHeaderInformation(record.header);
                            seeds.clear();
                            strain_alleles::AlleleSeeds(alleles, index_alleles, static_cast<uint32_t>(taxonomic_id), static_cast<uint32_t>(gene_id),
                                                        std::string_view(record.sequence), kmers, kmer_length, kmer_handler, scratch, seeds);
                            for (auto const& pair : seeds) putter.FirstPut(pair.first);
                        }

                        if constexpr(KmerStatisticsConcept<KmerHandler>) {
                            thread_statistics.kmers_total += kmer_handler.TotalKmers();
                        }

                        if constexpr(KmerStatisticsConcept<KmerHandler>) {
                            thread_statistics.kmers_accepted += kmers.size();
                        }

                        if constexpr(debug == DEBUG_VERBOSE) {
                            thread_statistics.WriteStats(std::cout);
                        }
                        if constexpr(debug == DEBUG_EXTRAVERBOSE) {

                        }
                    }

#pragma omp critical(statistics)
                    statistics.Join(thread_statistics);
                    std::cout << "minimizers: " << thread_statistics.kmers_accepted << std::endl;
                }
        }
        if (is.bad()) {  // a read or decompression error, not the end of the reference
            std::cerr << "Cannot read the reference " << options.GetSequenceFilePath() << " to its end (truncated or corrupt file?)" << std::endl;
            exit(8);
        }
        bm_pass1.Stop();
        bm_pass1.PrintResults();

        // Each key's value positions, from its count. The values are held packed from the start, in the widths
        // the reference's taxids, gene ids and positions need, as a query run holds them (Seedmap::PackedLayout;
        // 29 instead of 36 GB with the key map at r226); only the files keep the 8-byte layout. A position is a
        // k-mer core's, flex_k/2 bases into its gene.
        Benchmark bm_pointers("Value pointers");
        bm_pointers.Start();
        auto const [max_taxid, max_gene, max_length] = genomes.IndexFieldMaxima();
        Seedmap::PackedLayout const layout = Seedmap::PackedLayout::For(max_taxid, max_gene, max_length + map.m_flex_k);
        putter.InitializeForPut(serial ? 1 : threads, &layout);
        bm_pointers.Stop();
        bm_pointers.PrintResults();

        // Back to the start of the reference for the second pass
        if (!input->Rewind()) {
            std::cerr << "Cannot re-read " << input->Path() << std::endl;
            exit(8);
        }

        Benchmark bm_pass2("Pass 2 (place the values)");
        bm_pass2.Start();
        if (!serial) {
            // Items: each k-mer with its value (taxon, gene, position of its core), placed by the
            // thread that owns its range into the key's next empty slot, in reference order.
            struct Placement { uint64_t key; uint64_t value; bool allele; };
            std::cout << "After first put" << std::endl;
            with_source([&](auto& source) {
                PartitionedPass<Placement>(source, kmer_handler_global, threads, ranges,
                    [&](KmerHandler& handler, size_t taxonomic_id, size_t gene_id, std::string_view sequence, KmerList& kmers,
                        std::vector<Placement>& items, Statistics& stats) {
                        if (options.HasBuildGeneSubset() && !options.BuildGeneAllowed(gene_id)) return;
                        stats.reads++;
                        kmers.clear();
                        handler(sequence, kmers);
                        DropAmbiguousKmers(sequence, kmer_length, kmers);
                        for (auto const& pair : kmers) {
                            ValueEntry entry;
                            entry.Put(taxonomic_id, gene_id, pair.second + map.m_flex_k_half);  // as Seedmap::PutOMP
                            items.push_back({ pair.first, entry.value, false });
                        }
                        if constexpr(KmerStatisticsConcept<KmerHandler>) {
                            stats.kmers_total += handler.TotalKmers();
                            stats.kmers_accepted += handler.TotalMinimizers();
                        }
                        if (seed_alleles) {
                            static thread_local strain_alleles::SeedScratch scratch;
                            static thread_local KmerList seeds;
                            seeds.clear();
                            size_t const used = strain_alleles::AlleleSeeds(alleles, index_alleles, static_cast<uint32_t>(taxonomic_id),
                                                                            static_cast<uint32_t>(gene_id), sequence, kmers, kmer_length,
                                                                            handler, scratch, seeds);
                            allele_copies.fetch_add(used > 0, std::memory_order_relaxed);
                            alleles_seeded.fetch_add(used, std::memory_order_relaxed);
                            for (auto const& pair : seeds) {
                                ValueEntry entry;
                                entry.Put(taxonomic_id, gene_id, pair.second + map.m_flex_k_half);
                                entry.SetFlagNonUnique();
                                items.push_back({ pair.first, entry.value, true });
                            }
                        }
                    },
                    [&map, range_shift](Placement const& p) { return static_cast<size_t>(map.MainKey(p.key) >> range_shift); },
                    [&](Placement const& p) {
                        bool const placed = map.PutOwned(p.key, p.value);
                        if (placed && p.allele) note_allele_seed(p.value);
                    },
                    statistics);
            });
        } else {
#pragma omp parallel default(none) shared(std::cout, options, is, dummy, read_count, kmer_handler_global, statistics, putter, flex_k_bits, main_k_bits, kmer_length, seed_alleles, index_alleles, alleles, allele_copies, alleles_seeded, note_allele_seed)
    {
        // Private variables
        FastxRecord record;

        // Extract variables from kmi_global
        KmerHandler kmer_handler(kmer_handler_global);
        SeqReader reader { is };
        size_t kmer = 0;
        Statistics thread_statistics;
        thread_statistics.thread_num = omp_get_thread_num();

        size_t taxid = 1;
        size_t geneid = 1;
        size_t genepos = 1;

        KmerList kmers;
        strain_alleles::SeedScratch scratch;  // the deep alleles' seeds (--index_alleles)
        KmerList seeds;

        std::cout << "After first put" << std::endl;
        while (reader(record)) {
            kmer_handler.SetSequence(std::string_view(record.sequence));

            auto [taxonomic_id, gene_id] = KmerUtils::ExtractHeaderInformation(record.header);
            if (options.HasBuildGeneSubset() && !options.BuildGeneAllowed(gene_id)) {
                continue;
            }

            thread_statistics.reads++;

            // Retrieve kmers
            kmers.clear();
            kmer_handler(std::string_view(record.sequence), kmers);
            DropAmbiguousKmers(record.sequence, kmer_length, kmers);
            for (auto pair : kmers) {
                size_t pos = pair.second;
                putter.Put(pair.first, taxonomic_id, gene_id, pos);
            }
            if (seed_alleles) {
                seeds.clear();
                size_t const used = strain_alleles::AlleleSeeds(alleles, index_alleles, static_cast<uint32_t>(taxonomic_id),
                                                                static_cast<uint32_t>(gene_id), std::string_view(record.sequence), kmers,
                                                                kmer_length, kmer_handler, scratch, seeds);
                allele_copies.fetch_add(used > 0, std::memory_order_relaxed);
                alleles_seeded.fetch_add(used, std::memory_order_relaxed);
                for (auto const& pair : seeds) {
                    ValueEntry entry;
                    entry.Put(taxonomic_id, gene_id, pair.second + putter.GetMap().m_flex_k_half);
                    entry.SetFlagNonUnique();
                    if (putter.PutValue(pair.first, entry.value)) note_allele_seed(entry.value);
                }
            }

            if constexpr(KmerStatisticsConcept<KmerHandler>) {
                thread_statistics.kmers_total += kmer_handler.TotalKmers();
            }
            if constexpr(KmerStatisticsConcept<KmerHandler>) {
                thread_statistics.kmers_accepted += kmer_handler.TotalMinimizers();
            }

            if constexpr(debug == DEBUG_VERBOSE) {
                thread_statistics.WriteStats(std::cout);
            }
            if constexpr(debug == DEBUG_EXTRAVERBOSE) {

            }
        }


#pragma omp critical(statistics)
        statistics.Join(thread_statistics);
    }
        }
        if (is.bad()) {
            std::cerr << "Cannot read the reference " << options.GetSequenceFilePath() << " to its end (truncated or corrupt file?)" << std::endl;
            exit(8);
        }
        bm_pass2.Stop();
        bm_pass2.PrintResults();
        PrintMemory("pass 2");
        if (seed_alleles) {
            std::cout << "Index alleles: " << allele_kmers.load() << " seeds of " << alleles_seeded.load() << " strain alleles at least "
                      << index_alleles << " from the representative's copy, of " << allele_copies.load() << " gene copies, as non-unique "
                      << "entries under their species (--index_alleles; left out of the totals in " << Options::PROTAL_UNIQUE_KMER_FILE
                      << ")" << std::endl;
        } else if (!alleles.Empty()) {
            std::cout << "Index alleles: none (--index_alleles " << index_alleles << ")" << std::endl;
        }

        // Options falls back to --reference when no --full_reference is given, so unique_kmers.tsv is
        // always written: a database without it cannot detect anything.
        std::cout << "Check Uniqueness: " << options.GetFullSequenceFilePath() << std::endl;
        Benchmark bm_unique("Uniqueness check");
        bm_unique.Start();

        omp_set_num_threads(options.GetThreads());
        if (!serial) {
            // By key range (UniquenessPass): the full reference's frames decompressed by several threads at once when it
            // has them (FastaFrames: the converter's full_reference.fna.zst; each such thread holds a frame's 128 MB
            // window), else read in one stream.
            std::string const full_path = zstd::Resolve(options.GetFullSequenceFilePath());
            std::string error;
            auto const frames = FastaFrames::Open(full_path, error);
            if (!error.empty()) {
                std::cerr << "Cannot read the full reference: " << error << std::endl;
                exit(8);
            }
            UniquenessCounts counts;
            if (frames) {
                FramesSource source(*frames, std::clamp<size_t>(static_cast<size_t>(threads) / 2, 1, kFrameDecoders), options.IndexBatchBytes());
                std::cout << "Uniqueness check: " << frames->Frames() << " frames of " << full_path << ", decompressed on "
                          << source.Decoders() << " threads at once" << std::endl;
                counts = UniquenessPass(source, kmer_handler_global, threads, ranges, range_shift, map, genomes, options, kmer_length);
            } else {
                auto full_input = OpenInput(options.GetFullSequenceFilePath());
                TextSource source(full_input->Stream(), options.IndexBatchBytes());
                counts = UniquenessPass(source, kmer_handler_global, threads, ranges, range_shift, map, genomes, options, kmer_length);
            }
            bm_unique.Stop();
            bm_unique.PrintResults();
            std::cout << "Uniqueness check: " << counts.kmers << " k-mers, " << counts.distinct << " of them distinct in their key "
                      << "range and round and looked up once each; " << counts.cells << " flex cells compared (" << std::fixed
                      << std::setprecision(1) << counts.cells / std::max(1.0, double(counts.kmers)) << std::defaultfloat
                      << " per k-mer), " << counts.shared << " whole k-mers found under another taxon or more than once, "
                      << counts.singles << " single entries read back from their genes" << std::endl;
        } else {
        auto full_input = OpenInput(options.GetFullSequenceFilePath());
        std::istream& full_is = full_input->Stream();
        KmerLookupSM lookup_global(putter.GetMap());
        // Genes are read back to check single-entry k-mers. They are preloaded unless
        // --preload_genomes_off is given; then GetGeneOMP loads each genome on first use.
        // Counted for the log: k-mers looked up, flex parts compared, k-mers found in the index
        // under another taxon or more than once, single entries read back from their gene.
        size_t kmers_checked = 0, flex_compared = 0, kmers_shared = 0, singles_read = 0;
        bool parsed = true;

#pragma omp parallel default(none) shared(std::cout, lookup_global, options, full_is, dummy, read_count, kmer_handler_global, statistics, putter, flex_k_bits, main_k_bits, genomes, kmers_checked, flex_compared, kmers_shared, singles_read, kmer_length, parsed)
    {
        // Private variables
        FastxRecord record;

        // Extract variables from kmi_global
        KmerHandler kmer_handler(kmer_handler_global);
        SeqReader reader { full_is };
        Statistics thread_statistics;
        thread_statistics.thread_num = omp_get_thread_num();

        size_t taxid = 1;
        size_t geneid = 1;
        size_t genepos = 1;

        KmerList kmers;
        KmerLookupSM lookup(lookup_global);

        Seedmap& index = putter.GetMap();
        Seedmap::PackedBlock block;     // a core's values (the index is packed)
        std::vector<uint32_t> exact;    // the entries of those whose whole k-mer is the k-mer's
        size_t local_kmers = 0, local_compared = 0, local_shared = 0, local_singles = 0;

        while (reader(record)) {
            kmer_handler.SetSequence(std::string_view(record.sequence));

            auto [taxonomic_id, gene_id] = KmerUtils::ExtractHeaderInformation(record.header);
            // With --build_gene_subset the database's genes are the subset: a copy of another gene is
            // not among the sequences a read could come from, so it does not make k-mers non-unique
            // (nor does it cost time: most of the full reference with a small subset).
            if (options.HasBuildGeneSubset() && !options.BuildGeneAllowed(gene_id)) continue;

            thread_statistics.reads++;

            // Retrieve kmers
            kmers.clear();
            kmer_handler(std::string_view(record.sequence), kmers);
            DropAmbiguousKmers(record.sequence, kmer_length, kmers);
            local_kmers += kmers.size();

            for (auto pair : kmers) {
                // A k-mer counts where the whole of it is indexed (core and flex part): under another
                // taxon, or more than once, all its values become non-unique.
                exact.clear();
                local_compared += lookup.GetExact(pair.first, block, exact);
                if (exact.empty()) {
                    // A core that occurs once in the index has no flex keys, so GetExact cannot compare
                    // the whole k-mer; before, such entries were never checked and all stayed unique.
                    // Compare with the k-mer the entry was built from, read back from its gene, under
                    // the rule used for flex blocks: another taxon with the same k-mer makes it
                    // non-unique.
                    if (!lookup.GetSingleEntry(pair.first, block)) continue;
                    lookup.Entry(block, 0).Get(taxid, geneid, genepos);
                    if (taxonomic_id == taxid) continue;
                    local_singles++;
                    uint64_t indexed_kmer = 0;
                    if (!IndexedKmer(index, genomes, taxid, geneid, genepos, indexed_kmer) || indexed_kmer != pair.first) continue;
                    index.ClearUniqueFlag(block, 0);  // atomic
                    local_shared++;
                    continue;
                }

                lookup.Entry(block, exact.front()).Get(taxid, geneid, genepos);
                if (exact.size() == 1 && taxonomic_id == taxid) {
                    continue;
                }
                local_shared++;
                for (uint32_t const entry : exact) {
                    index.ClearUniqueFlag(block, entry);  // atomic
                }
            }
        }


    #pragma omp critical(statistics)
        {
            statistics.Join(thread_statistics);
            kmers_checked += local_kmers;
            flex_compared += local_compared;
            kmers_shared += local_shared;
            singles_read += local_singles;
            parsed = parsed && reader.Success();
        }
    }
        // A read or decompression error ends the stream early: the k-mers of the rest would stay unique.
        if (full_is.bad() || !parsed) {
            std::cerr << "Cannot read the full reference " << options.GetFullSequenceFilePath() << " to its end (truncated or corrupt "
                      << "file?)" << std::endl;
            exit(8);
        }
        bm_unique.Stop();
        bm_unique.PrintResults();
        std::cout << "Uniqueness check: " << kmers_checked << " k-mers, " << flex_compared << " flex parts compared ("
                  << std::fixed << std::setprecision(1) << flex_compared / std::max(1.0, double(kmers_checked))
                  << std::defaultfloat << " per k-mer), " << kmers_shared << " found under another taxon or more than once, "
                  << singles_read << " single entries read back from their genes" << std::endl;
        }

        // No gene is read from here on (the k-mer statistics need only how many genes each taxon has): their
        // sequences go before the index is written (~4 GB at r226).
        genomes.ReleaseGeneSequences();
        ReleaseFreeMemory();
        PrintMemory("the uniqueness check");

        std::cout << "Save unique kmer info: \n" << options.GetUniqueKmersFile() << std::endl;
        Benchmark bm_statistics("Unique k-mer statistics");
        bm_statistics.Start();
        // Via a .partial file, checked: a table cut short (a full disk) would leave the genes after the cut without
        // unique k-mers, so not hittable, in a database that is otherwise packed and verified as it is.
        std::string const unique_partial = options.GetUniqueKmersFile() + ".partial";
        std::ofstream os(unique_partial);
        auto const totals = putter.GetMap().CountUniqueKmers(os, gene_rows, static_cast<int>(options.GetThreads()),
                                                             seed_alleles ? &allele_seeds : nullptr);
        os.close();
        {
            std::error_code ec;
            if (!os) {
                std::cerr << "Writing " << unique_partial << " failed" << std::endl;
                std::filesystem::remove(unique_partial, ec);
                exit(8);
            }
            std::filesystem::rename(unique_partial, options.GetUniqueKmersFile(), ec);
            if (ec) {
                std::cerr << "Cannot rename " << unique_partial << " to " << options.GetUniqueKmersFile() << ": " << ec.message() << std::endl;
                exit(8);
            }
        }
        bm_statistics.Stop();
        bm_statistics.PrintResults();
        std::cout << "Distance-two flags: " << totals.comparisons << " flex parts compared" << std::endl;

        // Queries align against the database's reference.fna via reference.map: record which ones.
        // The fingerprint holds the uncompressed size, so it survives compressing reference.fna.
        putter.GetMap().SetReferenceFingerprint(
                ReferenceFingerprint::Of(options.GetSequenceMapFile(), options.GetSequenceFile()));
        // The single-file database (the default) takes the index straight from memory (BundleDatabase, called next);
        // separate files get index.prx.zst (or index.prx) here.
        if (!options.WriteBundle()) {
            SaveIndex(options, putter);
            PrintMemory("writing the index");
        }

        return statistics;
    }


}
