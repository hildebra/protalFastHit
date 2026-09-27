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
        auto const& sequence = genomes.GetGenome(taxid).GetGeneOMP(geneid).Sequence();
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
    // files queries read, compressed: reference.map, internal_taxonomy.dmp, unique_kmers.tsv and the
    // model (model.xml, else random_forest.xml).
    static std::vector<db::Source> BundleSources(protal::Options const& options) {
        namespace fs = std::filesystem;
        std::vector<db::Source> sources = {
                {Options::PROTAL_INDEX_FILE, options.ResolvedIndexFile()},
                {Options::PROTAL_SEQUENCE_FILE, options.ResolvedSequenceFile()},
                {Options::PROTAL_SEQUENCE_MAP_FILE, options.GetSequenceMapFile()},
                {Options::PROTAL_TAXONOMY_FILE, options.GetInternalTaxonomyFile()}};
        if (fs::exists(options.GetUniqueKmersFile())) sources.push_back({Options::PROTAL_UNIQUE_KMER_FILE, options.GetUniqueKmersFile()});
        for (std::string const model : {"model.xml", "random_forest.xml"}) {
            std::string const path = (fs::path(options.GetLocation().dir) / model).string();
            if (fs::exists(path)) {
                sources.push_back({model, path});
                break;
            }
        }
        return sources;
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
        if (std::none_of(sources.begin(), sources.end(), [](db::Source const& s) { return s.name.ends_with(".xml"); })) {
            std::cerr << "Warning: no model.xml in " << options.GetLocation().dir << "; profiling with " << target
                      << " then needs --model" << std::endl;
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

    template<typename KmerHandler, typename KmerPutter, DebugLevel debug>
    static void Check(protal::Options const& options, KmerPutter& putter, KmerHandler& kmer_handler_global) {

        // Shared
        auto input = OpenInput(options.GetSequenceFilePath());
        std::istream& is = input->Stream();
        size_t dummy = 0;
        int read_count = 0;

        // Set Thread Num
        omp_set_num_threads(1);

        size_t main_k = putter.GetMap().m_exact_k;
        size_t main_k_bits = putter.GetMap().m_main_bits;
        size_t flex_k = putter.GetMap().m_flex_k;
        size_t flex_k_bits = putter.GetMap().m_flex_k_bits;

        Statistics statistics;

        KmerLookupSM lookup_global(putter.GetMap());

#pragma omp parallel default(none) shared(std::cout, lookup_global, options, is, dummy, read_count, kmer_handler_global, statistics, putter, flex_k_bits, main_k_bits)
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
            KmerLookupSM lookup(lookup_global);

            LookupList seeds;
            std::vector<ValueEntry*> max_sim_entries;
            uint32_t max_sim = 0;
            uint32_t best_possible_sim = putter.GetMap().m_flex_k;

            while (reader(record)) {
                kmer_handler.SetSequence(std::string_view(record.sequence));

                auto [taxonomic_id, gene_id] = KmerUtils::ExtractHeaderInformation(record.header);
                if (options.HasBuildGeneSubset() && !options.BuildGeneAllowed(gene_id)) {
                    continue;
                }

                thread_statistics.reads++;

                // Retrieve kmers
                kmers.clear();
                seeds.clear();
                kmer_handler(std::string_view(record.sequence), kmers);

                for (auto pair : kmers) {
                    size_t pos = pair.second;
                    max_sim = 0;
                    // std::cout << pair.first << ", " << pair.second << std::endl;
                    max_sim_entries.clear();
                    lookup.GetFlex(pair.first, max_sim_entries, max_sim);
                    if (max_sim_entries.empty()) continue;
                }
            }


#pragma omp critical(statistics)
            statistics.Join(thread_statistics);
        }
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


        // TODO Fix multithreaded database building.
        if (options.GetThreads() > 1) {
            std::cerr << "Building db with multiple threads is currently broken" << std::endl;
            omp_set_num_threads(1);
        }

        std::cout << "Run Build" << std::endl;
//        if constexpr(protal::HasFirstPut<KmerPutter>) {
        if (true) {
#pragma omp parallel default(none) shared(std::cout, options, is, dummy, read_count, kmer_handler_global, statistics, putter, main_k_bits, flex_k_bits)
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

                // E.g. if k-mers are counted before they are inserted, call Initialize for put
                // To calculate the bucket sizes
                putter.InitializeForPut();

        }

        // Back to the start of the reference for the second pass
        if (!input->Rewind()) {
            std::cerr << "Cannot re-read " << input->Path() << std::endl;
            exit(8);
        }

#pragma omp parallel default(none) shared(std::cout, options, is, dummy, read_count, kmer_handler_global, statistics, putter, flex_k_bits, main_k_bits)
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

        // Options falls back to --reference when no --full_reference is given, so unique_kmers.tsv is
        // always written: a database without it cannot detect anything.
        std::cout << "Check Uniqueness: " << options.GetFullSequenceFilePath() << std::endl;

        omp_set_num_threads(options.GetThreads());
        auto full_input = OpenInput(options.GetFullSequenceFilePath());
        std::istream& full_is = full_input->Stream();
        KmerLookupSM lookup_global(putter.GetMap());
        // Genes are read back to check single-entry k-mers. They are preloaded unless
        // --preload_genomes_off is given; then GetGeneOMP loads each genome on first use.

#pragma omp parallel default(none) shared(std::cout, lookup_global, options, full_is, dummy, read_count, kmer_handler_global, statistics, putter, flex_k_bits, main_k_bits, genomes)
    {
        // Private variables
        FastxRecord record;

        // Extract variables from kmi_global
        KmerHandler kmer_handler(kmer_handler_global);
        SeqReader reader { full_is };
        size_t kmer = 0;
        Statistics thread_statistics;
        thread_statistics.thread_num = omp_get_thread_num();

        size_t taxid = 1;
        size_t geneid = 1;
        size_t genepos = 1;

        KmerList kmers;
        KmerLookupSM lookup(lookup_global);

        LookupList seeds;
        std::vector<ValueEntry*> max_sim_entries;
        uint32_t max_sim = 0;
        uint32_t best_possible_sim = putter.GetMap().m_flex_k;

        while (reader(record)) {
            kmer_handler.SetSequence(std::string_view(record.sequence));

            auto [taxonomic_id, gene_id] = KmerUtils::ExtractHeaderInformation(record.header);

            thread_statistics.reads++;

            // Retrieve kmers
            kmers.clear();
            seeds.clear();
            kmer_handler(std::string_view(record.sequence), kmers);

            for (auto pair : kmers) {
                size_t pos = pair.second;
                max_sim = 0;
                // std::cout << pair.first << ", " << pair.second << std::endl;
                max_sim_entries.clear();
                lookup.GetFlex(pair.first, max_sim_entries, max_sim);
                if (max_sim_entries.empty()) {
                    // A core that occurs once in the index has no flex keys, so GetFlex cannot compare
                    // the whole k-mer; before, such entries were never checked and all stayed unique.
                    // Compare with the k-mer the entry was built from, read back from its gene, under
                    // the rule used for flex blocks: another taxon with the same k-mer makes it
                    // non-unique.
                    ValueEntry* single = lookup.GetSingleEntry(pair.first);
                    if (single == nullptr) continue;
                    single->Get(taxid, geneid, genepos);
                    if (taxonomic_id == taxid) continue;
                    uint64_t indexed_kmer = 0;
                    if (!IndexedKmer(putter.GetMap(), genomes, taxid, geneid, genepos, indexed_kmer) || indexed_kmer != pair.first) continue;
#pragma omp critical(SetNonUnique)
                    single->SetFlagNonUnique();
                    continue;
                }

                max_sim_entries.front()->Get(taxid, geneid, genepos);

                if (max_sim != best_possible_sim) {
                    continue;
                }
                if (max_sim_entries.size() == 1 && taxonomic_id == taxid) {
                    continue;
                }

                for (auto& entry : max_sim_entries) {
                    // std::cout << "SetNonUnique " << taxonomic_id << " != " << taxid << " entries: " << max_sim_entries.size() << " Isunique? " << entry->IsFlagUnique();
#pragma omp critical(SetNonUnique)
                    entry->SetFlagNonUnique();
                    // std::cout << " -> " << entry->IsFlagUnique() << std::endl;
                }
            }
        }


    #pragma omp critical(statistics)
        statistics.Join(thread_statistics);
    }

        std::cout << "Save unique kmer info: \n" << options.GetUniqueKmersFile() << std::endl;
        std::ofstream os(options.GetUniqueKmersFile());
        putter.GetMap().CountUniqueKmers(os, true, true);
        os.close();

        // Queries align against the database's reference.fna via reference.map: record which ones.
        // The fingerprint holds the uncompressed size, so it survives compressing reference.fna.
        putter.GetMap().SetReferenceFingerprint(
                ReferenceFingerprint::Of(options.GetSequenceMapFile(), options.GetSequenceFile()));
        SaveIndex(options, putter);

        return statistics;
    }


}
