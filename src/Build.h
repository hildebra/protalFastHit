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

namespace protal::build {
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
    // gene_conservation.tsv if the build wrote one, gene_neighbours.tsv if the folder has one (written by
    // scripts/mini_db/gene_neighbours.py, checked by CheckGeneNeighbours), and every presence model there is (AllModelFiles in
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
        if (fs::exists(options.GetGeneNeighboursFile())) sources.push_back({Options::PROTAL_GENE_NEIGHBOURS_FILE, options.GetGeneNeighboursFile()});
        for (auto const& model : AllModelFiles()) {
            std::string const path = (fs::path(options.GetLocation().dir) / model).string();
            if (fs::exists(path)) sources.push_back({model, path});
        }
        return sources;
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

    // Packs the database folder into database.protal, checks it (db::Write), and removes the files
    // it now holds; --unpack_db writes them back. The index must be in the column format, as --build
    // and --compress_db write it.
    static void BundleDatabase(protal::Options const& options) {
        namespace fs = std::filesystem;
        std::string const target = (fs::path(options.GetLocation().dir) / db::kFileName).string();
        auto const sources = BundleSources(options);
        if (!index_codec::IsSplitIndex(sources.front().path)) {
            std::cerr << "Cannot write " << target << ": the index " << sources.front().path << " is not in protal's column "
                      << "format (protal --compress_db --no_bundle converts it)" << std::endl;
            exit(8);
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
            before += fs::file_size(source.path, ec);
            names += (names.empty() ? "" : ", ") + fs::path(source.path).filename().string();
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
                  << (replaces ? "; replaced the previous one" : "") << "). Removed the separate files; "
                  << options.ProgramName() << " --unpack_db --db " << Options::ShellWord(target) << " writes them back." << std::endl;
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

    // --add_model (the model checked already): stores the PMML file model as member or file `name` of
    // the database, replacing the one there. database.protal is rewritten via database.protal.partial
    // with its other members' frames copied as they are (db::Write checks it); in a folder of separate
    // files the model is copied next to them.
    static void AddModel(protal::Options const& options, std::string const& model, std::string const& name) {
        namespace fs = std::filesystem;
        std::error_code ec;
        if (!options.IsBundle()) {
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
            return;
        }
        auto const& bundle = *options.GetBundle();
        std::vector<db::Source> sources;
        bool replaces = false;
        for (auto const& member : bundle.Members()) {
            if (member.name == name) {
                sources.push_back({name, model});
                replaces = true;
            } else {
                sources.push_back({member.name, bundle.Path(), member.frames});
            }
        }
        if (!replaces) sources.push_back({name, model});
        std::string error;
        auto const written = db::Write(bundle.Path(), sources, options.CompressionParams(), error);
        if (!written) {
            std::cerr << "Writing " << bundle.Path() << " failed: " << error << " (the database is unchanged)" << std::endl;
            exit(8);
        }
        std::vector<std::string> names;
        for (auto const& source : sources) names.push_back(source.name);
        std::cout << "Stored " << model << " as " << name << " in " << bundle.Path() << (replaces ? " (replaced the previous one)" : "")
                  << "; " << HumanBytes(*written) << ". Models for read types: " << ModelCoverage(names).first << std::endl;
    }

    // --compress_db: rewrites an existing database's index and reference compressed, without
    // rebuilding it, e.g. a downloaded raw database or one compressed as a single frame: the index
    // in the column format (IndexCodec.h), reference.fna as seekable zstd, all packed into
    // database.protal unless --no_bundle. Each new file is read back and compared with the old
    // content before the old file is removed.
    static void CompressDatabase(protal::Options const& options) {
        if (options.IsBundle()) {
            std::cout << options.GetBundle()->Path() << " is already a single-file database; kept" << std::endl;
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

    static gene_conservation::Estimate WriteGeneConservation(protal::Options const& options, GenomeLoader& genomes,
                                                             std::string const& full_reference, std::string const& target) {
        Benchmark bm("Gene conservation");
        bm.Start();
        std::vector<uint64_t> keys;
        for (auto const& [taxid, genome] : genomes.GetGenomeMap()) {
            auto const& genes = genome.GetGeneList();
            for (size_t i = 0; i < genes.size(); i++) {
                if (genes[i].IsSet()) keys.push_back(gene_conservation::Estimator::Key(taxid, i + 1));
            }
        }
        gene_conservation::Estimator estimator(std::move(keys));
        auto input = OpenInput(full_reference);
        std::istream& is = input->Stream();
        omp_set_num_threads(static_cast<int>(std::max<size_t>(options.GetThreads(), 1)));
#pragma omp parallel default(none) shared(is, genomes, estimator)
        {
            FastxRecord record;
            SeqReader reader { is };
            while (reader(record)) {
                auto const [taxid, geneid] = KmerUtils::ExtractHeaderInformation(record.header);
                if (!genomes.HasGene(taxid, geneid)) continue;
                auto const slot = estimator.Take(taxid, geneid);
                if (!slot) continue;
                auto const rep = genomes.GetGenome(taxid).GetGeneOMP(geneid).Sequence();
                estimator.Add(*slot, rep.View(), record.sequence);
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

    // A pass over the reference in `threads` threads that updates the index as one thread would
    // (docs/claude/2026-09-29-index-build-parallel.md). The key space is cut into `ranges` ranges of
    // whole control blocks. Each round, one thread reads batches of whole records (FastaBatches);
    // the threads parse them, extract their items (extract: a record's k-mers, in order) and group
    // each batch's items by range (range_of), keeping their order within a range; then each range
    // is applied (apply) by one thread, batch by batch in reference order. A key's updates thus
    // come in the serial order, and none needs a lock: no two threads touch the same block.
    template<typename Item, typename KmerHandler, typename Extract, typename RangeOf, typename Apply>
    void PartitionedPass(std::istream& is, KmerHandler const& handler_global, int threads, size_t batch_bytes,
                         size_t ranges, Extract&& extract, RangeOf&& range_of, Apply&& apply, Statistics& statistics) {
        struct Batch {
            std::string text;
            std::vector<Item> items;        // grouped by range
            std::vector<uint32_t> offsets;  // items [offsets[r], offsets[r + 1]) are range r's
        };
        threads = std::max(threads, 1);
        size_t const per_round = 4 * static_cast<size_t>(threads);
        std::vector<Batch> batches(per_round);
        FastaBatches reader(is);
        size_t filled = 0;

#pragma omp parallel num_threads(threads)
        {
            KmerHandler handler(handler_global);
            Statistics local;
            std::vector<Item> items;  // a batch's items in record order
            std::vector<uint32_t> cursor(ranges);
            std::string scratch, header;
            KmerList kmers;
            while (true) {
#pragma omp single
                {
                    filled = 0;
                    while (filled < per_round && reader.Next(batches[filled].text, batch_bytes)) filled++;
                }
                if (filled == 0) break;

#pragma omp for schedule(dynamic, 1)
                for (size_t b = 0; b < filled; b++) {
                    Batch& batch = batches[b];
                    items.clear();
                    ForEachFastaRecord(batch.text, scratch, [&](std::string_view record_header, std::string_view sequence) {
                        header.assign(record_header);
                        extract(handler, header, sequence, kmers, items, local);
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
        if (!reader.Error().empty()) {
            std::cerr << "Cannot read the reference: " << reader.Error() << std::endl;
            exit(8);
        }
    }

    // Key ranges for PartitionedPass: at least 64 per thread for balance, and at least 64 control
    // blocks each.
    inline int IndexRangeBits(int threads, size_t main_bits) {
        int bits = 6;
        while (bits < 16 && (size_t{1} << bits) < 64 * static_cast<size_t>(std::max(threads, 1))) bits++;
        return std::min<int>(bits, static_cast<int>(main_bits) - 3 - 6);
    }

    template<typename KmerHandler, typename KmerPutter, DebugLevel debug>
    static Statistics Run(protal::Options const& options, KmerPutter& putter, KmerHandler& kmer_handler_global, GenomeLoader& genomes) {

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

        // Each phase is timed in the log ("... took"): where a build at GTDB scale spends its time.
        std::cout << "Run Build" << std::endl;
        Benchmark bm_pass1("Pass 1 (count the k-mers)");
        bm_pass1.Start();
        if (!serial) {
            // Items: the k-mers' main keys, counted up by the thread that owns their range.
            std::cout << "Build: iterate records" << std::endl;
            Statistics pass;
            PartitionedPass<uint32_t>(is, kmer_handler_global, threads, options.IndexBatchBytes(), ranges,
                [&](KmerHandler& handler, std::string const& header, std::string_view sequence, KmerList& kmers,
                    std::vector<uint32_t>& items, Statistics& stats) {
                    stats.reads++;
                    if (options.HasBuildGeneSubset()) {
                        auto [taxonomic_id, gene_id] = KmerUtils::ExtractHeaderInformation(header);
                        (void)taxonomic_id;
                        if (!options.BuildGeneAllowed(gene_id)) return;
                    }
                    kmers.clear();
                    handler(sequence, kmers);
                    DropAmbiguousKmers(sequence, kmer_length, kmers);
                    for (auto const& pair : kmers) items.push_back(static_cast<uint32_t>(map.MainKey(pair.first)));
                    if constexpr(KmerStatisticsConcept<KmerHandler>) {
                        stats.kmers_total += handler.TotalKmers();
                        stats.kmers_accepted += kmers.size();
                    }
                },
                [range_shift](uint32_t main_key) { return static_cast<size_t>(main_key >> range_shift); },
                [&map](uint32_t main_key) { map.CountUpKey(main_key); },
                pass);
            std::cout << "minimizers: " << pass.kmers_accepted << std::endl;
            statistics.Join(pass);
        } else {
#pragma omp parallel default(none) shared(std::cout, options, is, dummy, read_count, kmer_handler_global, statistics, putter, main_k_bits, flex_k_bits, kmer_length)
                {
                    // Private variables
                    FastxRecord record;

                    // Extract variables from kmi_global
                    KmerHandler kmer_handler(kmer_handler_global);
                    SeqReader reader{ is };
                    Statistics thread_statistics;
                    thread_statistics.thread_num = omp_get_thread_num();

                    KmerList kmers;

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
        bm_pass1.Stop();
        bm_pass1.PrintResults();

        // Each key's value positions, from its count
        Benchmark bm_pointers("Value pointers");
        bm_pointers.Start();
        putter.InitializeForPut(serial ? 1 : threads);
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
            struct Placement { uint64_t key; uint64_t value; };
            std::cout << "After first put" << std::endl;
            PartitionedPass<Placement>(is, kmer_handler_global, threads, options.IndexBatchBytes(), ranges,
                [&](KmerHandler& handler, std::string const& header, std::string_view sequence, KmerList& kmers,
                    std::vector<Placement>& items, Statistics& stats) {
                    auto [taxonomic_id, gene_id] = KmerUtils::ExtractHeaderInformation(header);
                    if (options.HasBuildGeneSubset() && !options.BuildGeneAllowed(gene_id)) return;
                    stats.reads++;
                    kmers.clear();
                    handler(sequence, kmers);
                    DropAmbiguousKmers(sequence, kmer_length, kmers);
                    for (auto const& pair : kmers) {
                        ValueEntry entry;
                        entry.Put(taxonomic_id, gene_id, pair.second + map.m_flex_k_half);  // as Seedmap::PutOMP
                        items.push_back({ pair.first, entry.value });
                    }
                    if constexpr(KmerStatisticsConcept<KmerHandler>) {
                        stats.kmers_total += handler.TotalKmers();
                        stats.kmers_accepted += handler.TotalMinimizers();
                    }
                },
                [&map, range_shift](Placement const& p) { return static_cast<size_t>(map.MainKey(p.key) >> range_shift); },
                [&map](Placement const& p) { map.PutOwned(p.key, p.value); },
                statistics);
        } else {
#pragma omp parallel default(none) shared(std::cout, options, is, dummy, read_count, kmer_handler_global, statistics, putter, flex_k_bits, main_k_bits, kmer_length)
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
        bm_pass2.Stop();
        bm_pass2.PrintResults();

        // Options falls back to --reference when no --full_reference is given, so unique_kmers.tsv is
        // always written: a database without it cannot detect anything.
        std::cout << "Check Uniqueness: " << options.GetFullSequenceFilePath() << std::endl;
        Benchmark bm_unique("Uniqueness check");
        bm_unique.Start();

        omp_set_num_threads(options.GetThreads());
        auto full_input = OpenInput(options.GetFullSequenceFilePath());
        std::istream& full_is = full_input->Stream();
        KmerLookupSM lookup_global(putter.GetMap());
        // Genes are read back to check single-entry k-mers. They are preloaded unless
        // --preload_genomes_off is given; then GetGeneOMP loads each genome on first use.
        // Counted for the log: k-mers looked up, flex parts compared, k-mers found in the index
        // under another taxon or more than once, single entries read back from their gene.
        size_t kmers_checked = 0, flex_compared = 0, kmers_shared = 0, singles_read = 0;

#pragma omp parallel default(none) shared(std::cout, lookup_global, options, full_is, dummy, read_count, kmer_handler_global, statistics, putter, flex_k_bits, main_k_bits, genomes, kmers_checked, flex_compared, kmers_shared, singles_read, kmer_length)
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

        std::vector<ValueEntry*> exact;
        size_t local_kmers = 0, local_compared = 0, local_shared = 0, local_singles = 0;

        while (reader(record)) {
            kmer_handler.SetSequence(std::string_view(record.sequence));

            auto [taxonomic_id, gene_id] = KmerUtils::ExtractHeaderInformation(record.header);

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
                local_compared += lookup.GetExact(pair.first, exact);
                if (exact.empty()) {
                    // A core that occurs once in the index has no flex keys, so GetExact cannot compare
                    // the whole k-mer; before, such entries were never checked and all stayed unique.
                    // Compare with the k-mer the entry was built from, read back from its gene, under
                    // the rule used for flex blocks: another taxon with the same k-mer makes it
                    // non-unique.
                    ValueEntry* single = lookup.GetSingleEntry(pair.first);
                    if (single == nullptr) continue;
                    single->Get(taxid, geneid, genepos);
                    if (taxonomic_id == taxid) continue;
                    local_singles++;
                    uint64_t indexed_kmer = 0;
                    if (!IndexedKmer(putter.GetMap(), genomes, taxid, geneid, genepos, indexed_kmer) || indexed_kmer != pair.first) continue;
                    single->SetFlagNonUnique();  // atomic
                    local_shared++;
                    continue;
                }

                exact.front()->Get(taxid, geneid, genepos);
                if (exact.size() == 1 && taxonomic_id == taxid) {
                    continue;
                }
                local_shared++;
                for (auto& entry : exact) {
                    entry->SetFlagNonUnique();  // atomic
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
        }
    }
        bm_unique.Stop();
        bm_unique.PrintResults();
        std::cout << "Uniqueness check: " << kmers_checked << " k-mers, " << flex_compared << " flex parts compared ("
                  << std::fixed << std::setprecision(1) << flex_compared / std::max(1.0, double(kmers_checked))
                  << std::defaultfloat << " per k-mer), " << kmers_shared << " found under another taxon or more than once, "
                  << singles_read << " single entries read back from their genes" << std::endl;

        WriteGeneConservation(options, genomes, options.GetFullSequenceFilePath(), options.GetGeneConservationFile());
        CheckGeneNeighbours(options, genomes);

        std::cout << "Save unique kmer info: \n" << options.GetUniqueKmersFile() << std::endl;
        Benchmark bm_statistics("Unique k-mer statistics");
        bm_statistics.Start();
        std::ofstream os(options.GetUniqueKmersFile());
        auto const totals = putter.GetMap().CountUniqueKmers(os, GeneRowsOf(genomes), static_cast<int>(options.GetThreads()));
        os.close();
        bm_statistics.Stop();
        bm_statistics.PrintResults();
        std::cout << "Distance-two flags: " << totals.comparisons << " flex parts compared" << std::endl;

        // Queries align against the database's reference.fna via reference.map: record which ones.
        // The fingerprint holds the uncompressed size, so it survives compressing reference.fna.
        putter.GetMap().SetReferenceFingerprint(
                ReferenceFingerprint::Of(options.GetSequenceMapFile(), options.GetSequenceFile()));
        SaveIndex(options, putter);

        return statistics;
    }


}
