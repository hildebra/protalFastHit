//
// Created by fritsche on 15/07/22.
//

#pragma once

#include <Constants.h>
#include "TargetClones.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <cctype>
#include <charconv>
#include <cstring>
#include <filesystem>
#include <memory>
#include <optional>
#include <sys/mman.h>
#include <unistd.h>
#include <sparse_map.h>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <fstream>
#include <err.h>
#include <KmerUtils.h>


#include <sparse_set.h>

#include "Utilities.h"
#include "Zstd.h"
#include "Database.h"
#include "ReferenceFingerprint.h"
#include "GeneTableFile.h"
#include <sysexits.h>

#include "Benchmark.h"
#include "PackedSequence.h"
#include "GeneConservation.h"
#include "GeneIncongruence.h"
#include "GeneNeighbours.h"
#include "SpeciesPriors.h"

namespace protal {
    // The database's gene tables (reference.map, unique_kmers.tsv), one line per gene (16.6M at GTDB
    // r226 size), parsed fast: the file is read in pieces of whole lines (~64 MB, reusing the
    // buffers: holding a whole table and its rows cost more in page faults than it saved), each
    // piece is cut into chunks that are parsed in parallel with std::from_chars, and the rows are
    // then added to the genomes in file order. Each chunk stops at its first problem, and a chunk's
    // problem is reported after its rows are added, so the first problem in the file is the one
    // reported, with its line number, as when the file was read line by line.
    namespace gene_table {
        inline constexpr size_t kPieceBytes = size_t{64} << 20;

        // A line's tab-separated fields, as Utils::split gives them (n tabs, n + 1 fields); the
        // first `max` of them in f.
        template<size_t max>
        struct Fields {
            std::string_view f[max];
            size_t n = 0;
            explicit Fields(std::string_view line) {
                size_t start = 0;
                for (;;) {
                    size_t const tab = line.find('\t', start);
                    std::string_view const field = line.substr(start, tab == std::string_view::npos ? std::string_view::npos : tab - start);
                    if (n < max) f[n] = field;
                    n++;
                    if (tab == std::string_view::npos) break;
                    start = tab + 1;
                }
            }
        };

        enum class NumberProblem { None, NotANumber, TooLarge };

        // A non-negative integer field: digits only (no sign or spaces).
        inline NumberProblem Number(std::string_view s, uint64_t& value) {
            if (s.empty() || !std::all_of(s.begin(), s.end(), [](char c) { return c >= '0' && c <= '9'; })) return NumberProblem::NotANumber;
            auto const [end, ec] = std::from_chars(s.data(), s.data() + s.size(), value);
            return ec == std::errc() && end == s.data() + s.size() ? NumberProblem::None : NumberProblem::TooLarge;
        }

        // Why column `column` (1-based) is not a number, as the tables' messages say it.
        inline std::string ColumnProblem(size_t column, NumberProblem problem) {
            return "column " + std::to_string(column) + (problem == NumberProblem::TooLarge ? " is too large" : " is not a non-negative integer");
        }

        // The rows of one chunk: parsed up to its first problem (line within the chunk, message).
        template<typename Row>
        struct Chunk {
            std::vector<Row> rows;       // with their line within the chunk (1-based)
            size_t lines = 0;            // lines in the chunk
            size_t problem_line = 0;     // 0: none
            std::string problem;
        };

        // Parses piece (whole lines) in chunks with `threads` threads, into chunks (resized; their
        // buffers are reused). parse(line, row, problem) gets a line without its newline and a
        // trailing '\r' (never an empty one), and returns whether it filled row; a non-empty problem
        // stops the chunk.
        template<typename Row, typename Parse>
        void ParsePiece(std::string_view piece, int threads, Parse& parse, std::vector<Chunk<Row>>& chunks) {
            size_t const parts = piece.size() < (size_t{1} << 20) ? 1 : static_cast<size_t>(std::max(threads, 1)) * 4;
            std::vector<size_t> bounds{ 0 };
            for (size_t i = 1; i < parts; i++) {
                size_t at = std::max(bounds.back(), piece.size() * i / parts);
                size_t const newline = piece.find('\n', at);
                at = newline == std::string_view::npos ? piece.size() : newline + 1;
                if (at > bounds.back() && at < piece.size()) bounds.push_back(at);
            }
            bounds.push_back(piece.size());
            chunks.resize(bounds.size() - 1);
            zstd::ParallelFor(chunks.size(), threads, [&](size_t c, size_t) -> std::string {
                auto& chunk = chunks[c];
                chunk.rows.clear();
                chunk.lines = 0;
                chunk.problem_line = 0;
                chunk.problem.clear();
                std::string_view rest = piece.substr(bounds[c], bounds[c + 1] - bounds[c]);
                Row row{};
                while (!rest.empty()) {
                    size_t const newline = rest.find('\n');
                    std::string_view line = rest.substr(0, newline);
                    rest = newline == std::string_view::npos ? std::string_view{} : rest.substr(newline + 1);
                    chunk.lines++;
                    if (!line.empty() && line.back() == '\r') line.remove_suffix(1);
                    if (line.empty()) continue;
                    if (parse(line, row, chunk.problem)) {
                        row.line = chunk.lines;
                        chunk.rows.push_back(row);
                    }
                    if (!chunk.problem.empty()) {
                        chunk.problem_line = chunk.lines;
                        break;
                    }
                }
                return {};
            });
        }

        // Where a table's load spends its wall-clock time, in seconds, for the run's log (GenomeLoader::GeneTableTimes):
        // reading and decompressing, parsing (ParsePiece), the serial grouping by genome with the genomes made
        // (AddPieceByGenome), the parallel adding, and a pass after the rows (the unique k-mers' sums per genome).
        struct Times {
            size_t rows = 0;
            double read = 0, parse = 0, group = 0, add = 0, after = 0;
        };

        inline double SecondsSince(std::chrono::steady_clock::time_point start) {
            return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        }

        // The whole content of a table, read and decompressed on `threads` threads (zstd::ParallelRead: a seekable zstd
        // member frame by frame, a raw file in chunks), for ForEachPiece. The sequential reader decompressed a 1 GB table
        // on one thread, which with the parsing and the adds parallel was most of the load at GTDB r226 size.
        class TableSink : public zstd::Sink {
        public:
            explicit TableSink(uint64_t size) : m_text(new char[size]), m_size(size) {}  // not zeroed: every byte is written
            char* Direct(uint64_t offset, size_t size) override {
                return offset + size <= m_size ? m_text.get() + offset : nullptr;
            }
            void Copy(uint64_t offset, char const* data, size_t size) override {
                if (offset + size > m_size) return;  // more than the member says it holds: ParallelRead reports the size read
                std::memcpy(m_text.get() + offset, data, size);
            }
            std::string_view Text(uint64_t read) const { return { m_text.get(), static_cast<size_t>(std::min<uint64_t>(read, m_size)) }; }
        private:
            std::unique_ptr<char[]> m_text;
            uint64_t m_size;
        };

        // Reads a gene table in pieces, parses each in parallel (ParsePiece) and hands each piece's chunks, with the
        // lines before each chunk, to add_piece(chunks, line_bases). Returns an error message if the file cannot be
        // read or decompressed, else empty. times (if given) gets the reading and parsing time and the rows.
        template<typename Row, typename Parse, typename AddPiece>
        std::string ForEachPiece(db::DbFile const& file, int threads, Parse&& parse, AddPiece&& add_piece, Times* times = nullptr,
                                XXHash64* hash = nullptr) {
            std::vector<size_t> line_bases;
            if (!file.Exists()) return "cannot open the file";
            std::vector<Chunk<Row>> chunks;
            size_t line_base = 0;    // lines before this piece
            Times ignored;
            Times& t = times ? *times : ignored;
            auto handle = [&](std::string_view piece_text) {
                if (hash) hash->add(piece_text.data(), piece_text.size());  // the pieces cover the content in order
                auto const start = std::chrono::steady_clock::now();
                ParsePiece<Row>(piece_text, threads, parse, chunks);
                t.parse += SecondsSince(start);
                line_bases.resize(chunks.size());
                for (size_t c = 0; c < chunks.size(); c++) {
                    line_bases[c] = line_base;
                    line_base += chunks[c].lines;
                    t.rows += chunks[c].rows.size();
                }
                add_piece(chunks, line_bases);
            };
            // A compressed table (a database member, a .zst file) with several threads: decompressed into memory on all
            // of them (TableSink), then parsed from there in pieces of whole lines. A raw file streams from the page cache
            // faster than it copies, and one thread has nothing to parallelise: those stream below.
            if (auto const size = file.Size(); size && threads > 1 && file.Compressed()) {
                auto const start = std::chrono::steady_clock::now();
                TableSink sink(*size);
                std::string error;
                uint64_t const read = file.ParallelRead(threads, sink, error);
                t.read += SecondsSince(start);
                if (!error.empty()) return "the file cannot be read or decompressed (" + error + ")";
                std::string_view const text = sink.Text(read);
                for (size_t pos = 0; pos < text.size();) {
                    size_t cut = text.size();
                    if (text.size() - pos > kPieceBytes) {
                        size_t const newline = text.find('\n', pos + kPieceBytes - 1);  // the first line end at or past the piece's
                        cut = newline == std::string_view::npos ? text.size() : newline + 1;
                    }
                    handle(text.substr(pos, cut - pos));
                    pos = cut;
                }
                return {};
            }
            auto input = file.Open();
            if (!input->IsOpen()) return "cannot open the file";
            std::istream* is = &input->Stream();
            std::string buffer;
            size_t kept = 0;         // bytes of an unfinished last line carried over to the next piece
            size_t piece = kPieceBytes;
            bool end = false;
            while (!end) {
                auto const start = std::chrono::steady_clock::now();
                buffer.resize(kept + piece);
                is->read(buffer.data() + kept, static_cast<std::streamsize>(piece));
                size_t const got = static_cast<size_t>(is->gcount());
                t.read += SecondsSince(start);
                if (is->bad()) return "the file cannot be read or decompressed (truncated or corrupt file?)";
                end = got < piece;
                size_t const size = kept + got;
                // Whole lines only, but for the end of the file.
                size_t cut = size;
                if (!end) {
                    size_t const newline = std::string_view(buffer.data(), size).rfind('\n');
                    if (newline == std::string_view::npos || newline < kept) {
                        kept = size;  // no line ends in this piece: read more
                        piece *= 2;
                        continue;
                    }
                    cut = newline + 1;
                }
                handle(std::string_view(buffer.data(), cut));
                std::memmove(buffer.data(), buffer.data() + cut, size - cut);
                kept = size - cut;
                piece = kPieceBytes;
            }
            return {};
        }

        // As ForEachPiece, with each chunk handed in file order to add(chunk, lines before the chunk).
        template<typename Row, typename Parse, typename Add>
        std::string ForEachChunk(db::DbFile const& file, int threads, Parse&& parse, Add&& add) {
            return ForEachPiece<Row>(file, threads, parse, [&](std::vector<Chunk<Row>> const& chunks, std::vector<size_t> const& line_bases) {
                for (size_t c = 0; c < chunks.size(); c++) add(chunks[c], line_bases[c]);
            });
        }

        // A piece's rows added by genome, each genome's rows by one thread: a genome's rows are contiguous in
        // protal's tables, so they form one run per chunk they lie in (a table that lists a genome in several
        // places gives it several runs; all of a genome's runs go to the one thread, in file order). Adding the
        // rows one by one on one thread took 6 s per run at GTDB r226 size, 24M rows per table
        // (docs/claude/2026-10-04-performance-gtdb-scale). prepare(taxid) is called on this thread for each
        // distinct genome of the piece before the rows are added, in file order (the genome map is changed there,
        // never in parallel); add(taxid, row) on any thread returns a problem with the row (empty: none). Returns
        // the earliest problem by line, a row's or a chunk's parse problem, with its line in the file (0: none),
        // as adding line by line would have met it first.
        template<typename Row, typename Prepare, typename Add>
        std::pair<size_t, std::string> AddPieceByGenome(std::vector<Chunk<Row>> const& chunks, std::vector<size_t> const& line_bases,
                                                        int threads, Prepare&& prepare, Add&& add, Times* times = nullptr) {
            auto const start = std::chrono::steady_clock::now();
            struct Run { size_t chunk, begin, end; };
            std::vector<std::vector<Run>> groups;
            std::vector<uint64_t> group_taxid;
            std::unordered_map<uint64_t, size_t> group_of;
            for (size_t c = 0; c < chunks.size(); c++) {
                auto const& rows = chunks[c].rows;
                for (size_t i = 0; i < rows.size();) {
                    size_t j = i + 1;
                    while (j < rows.size() && rows[j].taxid == rows[i].taxid) j++;
                    auto const [it, made] = group_of.try_emplace(rows[i].taxid, groups.size());
                    if (made) {
                        groups.emplace_back();
                        group_taxid.push_back(rows[i].taxid);
                        prepare(rows[i].taxid);
                    }
                    groups[it->second].push_back({ c, i, j });
                    i = j;
                }
            }
            std::vector<std::pair<size_t, std::string>> problems(groups.size());  // per genome: its first problem's line
            auto const grouped = std::chrono::steady_clock::now();
            if (times) times->group += std::chrono::duration<double>(grouped - start).count();
            struct AddTime {
                Times* times;
                std::chrono::steady_clock::time_point from;
                ~AddTime() { if (times) times->add += SecondsSince(from); }
            } add_time{ times, grouped };
            zstd::ParallelFor(groups.size(), threads, [&](size_t g, size_t) -> std::string {
                for (auto const& run : groups[g]) {
                    auto const& chunk = chunks[run.chunk];
                    for (size_t i = run.begin; i < run.end; i++) {
                        std::string problem = add(group_taxid[g], chunk.rows[i]);
                        if (!problem.empty()) {
                            problems[g] = { line_bases[run.chunk] + chunk.rows[i].line, std::move(problem) };
                            return {};
                        }
                    }
                }
                return {};
            });
            std::pair<size_t, std::string> first{ 0, {} };
            auto consider = [&first](size_t line, std::string const& problem) {
                if (line > 0 && (first.first == 0 || line < first.first)) first = { line, problem };
            };
            for (auto const& [line, problem] : problems) consider(line, problem);
            for (size_t c = 0; c < chunks.size(); c++) {
                if (!chunks[c].problem.empty()) consider(line_bases[c] + chunks[c].problem_line, chunks[c].problem);
            }
            return first;
        }
    }

    // One gene of the reference: where its sequence is (the byte in reference.fna until it is loaded,
    // then its 2-bit packed bytes, see PackedSequence.h), its length, and the k-mer counts of
    // unique_kmers.tsv. 32 bytes: a database has up to 16.6M of them.
    class Gene {
        static constexpr uint32_t kUnsetLength = 0x7fffffffu;

        uint64_t m_where = 0;            // not loaded: start byte in reference.fna; loaded: the packed bytes
        uint32_t m_id = 0;
        uint32_t m_length : 31 = kUnsetLength;
        uint32_t m_loaded : 1 = 0;
        uint32_t m_short_unique = 0;
        uint32_t m_long_unique = 0;
        uint32_t m_long_super_unique = 0;
        uint32_t m_total_kmers = 0;

    public:
        Gene(){};

        void Set(size_t id, size_t start_byte, size_t length) {
            m_id = static_cast<uint32_t>(id);
            m_where = start_byte;
            m_length = static_cast<uint32_t>(length);
            m_loaded = 0;
        }

        // The start byte in reference.fna; not kept once the gene is loaded (0).
        size_t GetStartByte() const {
            return m_loaded ? 0 : m_where;
        }

        // The gene's Bytes(GetLength()) packed bytes, in an arena or buffer that outlives the gene, to
        // be filled in place (packed::PackInto) if they are not yet.
        void SetPacked(uint8_t* data) {
            m_where = reinterpret_cast<uintptr_t>(data);
            m_loaded = 1;
        }

        uint8_t* MutablePacked() {
            return m_loaded ? reinterpret_cast<uint8_t*>(static_cast<uintptr_t>(m_where)) : nullptr;
        }

        // Forgets the packed bytes (their owner frees them, GenomeLoader::ReleaseGeneSequences): the gene is
        // not loaded and has no sequence from here.
        void ReleaseSequence() {
            if (m_loaded) {
                m_where = 0;
                m_loaded = 0;
            }
        }

        bool IsSet() const {
            return (m_id != 0 || m_length > 0) && m_length != kUnsetLength;
        }

        bool IsLoaded() const {
            return m_loaded;
        }

        size_t GetLength() const {
            return m_length;
        }

        // The gene's bases, decoded (empty while the gene is not loaded). Keep the result in a variable
        // for as long as a view of it is used.
        GeneSequence Sequence() const {
            return m_loaded ? GeneSequence(reinterpret_cast<const uint8_t*>(static_cast<uintptr_t>(m_where)), m_length)
                            : GeneSequence(std::string_view());
        }

        // As Sequence, decoding only the bases [begin, end) (cut to the gene): the view is as long as the gene
        // and indexed by gene position, but only those bases are valid (GeneSequence). For code that works
        // on a read's stretch of a gene: the stretch costs a fraction of the gene's decoding and cache misses.
        GeneSequence Window(size_t begin, size_t end) const {
            return m_loaded ? GeneSequence(reinterpret_cast<const uint8_t*>(static_cast<uintptr_t>(m_where)), m_length, begin, end)
                            : GeneSequence(std::string_view());
        }

        const size_t GetId() const {
            return m_id;
        }

        // The counts are stored in 32 bits (a gene has at most 2^20 - 1 bases).
        void SetUniqueValues(size_t short_unique, size_t long_unique, size_t long_super_unique, size_t total_kmers) {
            m_short_unique = static_cast<uint32_t>(short_unique);
            m_long_unique = static_cast<uint32_t>(long_unique);
            m_long_super_unique = static_cast<uint32_t>(long_super_unique);
            m_total_kmers = static_cast<uint32_t>(total_kmers);
        }

        [[nodiscard]] std::tuple<size_t, size_t, size_t, size_t> GetUniqueKmerCounts() const {
            return { m_short_unique, m_long_unique, m_long_super_unique, m_total_kmers };
        }

        bool HasShortUniques(const size_t threshold=0) const {
            return m_short_unique > threshold;
        }
        bool HasLongUniques(const size_t threshold=0) const {
            return m_long_unique > threshold;
        }
        bool HasLongSuperUniques(const size_t threshold=0) const {
            return m_long_super_unique > threshold;
        }

        double UniqueRate() const {
            return m_long_unique/static_cast<double>(m_total_kmers);
        }

        double SuperUniqueRate() const {
            return m_long_super_unique/static_cast<double>(m_total_kmers);
        }
    };
    static_assert(sizeof(Gene) == 32, "Gene is kept small: a database has up to 16.6M of them");


    class Genome {
    public:
        using GeneKey = size_t;
        using GenomeKey = size_t;
        using GeneList = std::vector<Gene>;

    private:
        using GeneID = uint32_t;
        GeneList m_genes;
        GenomeKey m_key;
        tsl::sparse_set<GeneID> m_hittable_genes;
        // Set once the database has said which genes are hittable (unique_kmers.tsv). Until then every
        // gene counts as hittable; after, only the listed ones, none if none are.
        bool m_hittable_known = false;
        // Whether every gene is read: set under critical(genome_loader) once they are, and read
        // without the lock by GetGeneOMP, so it is set with release and read with acquire (a thread
        // that sees it set also sees the genes). A genome is moved only while the map is built,
        // before threads use it.
        struct LoadedFlag {
            std::atomic<bool> value{ false };
            LoadedFlag() = default;
            LoadedFlag(LoadedFlag&& other) noexcept : value(other.value.load(std::memory_order_relaxed)) {}
            LoadedFlag& operator=(LoadedFlag&& other) noexcept {
                value.store(other.value.load(std::memory_order_relaxed), std::memory_order_relaxed);
                return *this;
            }
            bool Get() const { return value.load(std::memory_order_acquire); }
            void Set() { value.store(true, std::memory_order_release); }
        };
        LoadedFlag m_is_loaded;
        // Where genes are read from when they are loaded one by one: reference.fna, open (shared by all genomes,
        // used under omp critical(genome_loader)); null for a compressed reference, which only
        // GenomeLoader::LoadAllGenomes reads.
        std::ifstream* m_reader = nullptr;
        // The packed sequences of the genes loaded one by one (LoadAllGenomes keeps its genes in an arena).
        std::vector<std::unique_ptr<uint8_t[]>> m_owned;
        std::string m_scratch;
        bool m_released = false;  // the genes' sequences were freed (ReleaseSequences): none can be read again

        size_t m_short_unique = 0;
        size_t m_long_unique = 0;
        size_t m_long_super_unique = 0;
        size_t m_total_kmers = 0;



        const size_t GeneKeyToIndex(GeneKey const& key) const {
            return key - 1;
        }
    public:

//        Genome() {};
        explicit Genome(GenomeKey key) : m_key(key) {};

        GenomeKey GetKey() {
            return m_key;
        }

        void SetUniqueValues() {
            for (auto& gene : m_genes) {
                auto [su, lu, lsu, total] = gene.GetUniqueKmerCounts();
                m_short_unique += su;
                m_long_unique += lu;
                m_long_super_unique += lsu;
                m_total_kmers += total;
            }
        }

        [[nodiscard]] std::tuple<size_t, size_t, size_t, size_t> GetUniqueKmerCounts() const {
            return { m_short_unique, m_long_unique, m_long_super_unique, m_total_kmers };
        }

        void AddGene(GeneKey key, size_t id, size_t start_byte, size_t length, std::ifstream* is) {
            auto index = GeneKeyToIndex(key);
            if (index >= m_genes.size()) {
                m_genes.resize(index+1);
            }
            m_genes[index].Set(key, start_byte, length);
            m_reader = is;
        }

        void AddHittableGene(GeneID geneid) {
            m_hittable_known = true;
            m_hittable_genes.insert(geneid);
        }

        void SetHittableGenesKnown() {
            m_hittable_known = true;
        }

        bool IsGeneHittable(GeneID geneid) const {
            return !m_hittable_known || m_hittable_genes.contains(geneid);
        }

        size_t GenesWithShortUniques(size_t threshold = 0) {
            return std::count_if(m_genes.begin(), m_genes.end(), [threshold](Gene const& gene) {
                return gene.HasShortUniques(threshold);
            });
        }

        size_t GenesWithLongUniques(size_t threshold = 0) const {
            return std::count_if(m_genes.begin(), m_genes.end(), [threshold](Gene const& gene) {
                return gene.HasLongUniques(threshold);
            });
        }

        size_t GenesWithLongSuperUniques(size_t threshold = 0) const {
            return std::count_if(m_genes.begin(), m_genes.end(), [threshold](Gene const& gene) {
                return gene.HasLongSuperUniques(threshold);
            });
        }

        std::vector<uint32_t> GetHittableGenes() {
            std::vector<uint32_t> genes;
            if (!m_hittable_known) {
                genes.reserve(m_genes.size());
                for (auto const& gene : m_genes) {
                    if (gene.IsSet()) {
                        genes.emplace_back(gene.GetId());
                    }
                }
            } else {
                genes.reserve(m_hittable_genes.size());
                for (auto const& gene_id : m_hittable_genes) {
                    genes.emplace_back(gene_id);
                }
                std::sort(genes.begin(), genes.end());
            }
            return genes;
        }

        // Reads the gene from reference.fna and keeps it packed (PackedSequence.h). A failed read stops
        // protal: the stream is shared by all genes, and leaving it failed would silently turn every
        // gene read after it into NUL bytes.
        void LoadGeneData(Gene& gene) {
            if (gene.IsLoaded() || !gene.IsSet()) return;
            if (m_released) {
                errx(EX_SOFTWARE, "Gene %zu cannot be read: the genes' sequences were freed at the end of the build", gene.GetId());
            }
            size_t const length = gene.GetLength();
            if (m_reader && m_reader->is_open()) {
                m_reader->seekg(gene.GetStartByte());
                m_scratch.resize(length);
                m_reader->read(m_scratch.data(), length);
                if (!*m_reader || static_cast<size_t>(m_reader->gcount()) != length) {
                    std::cerr << "Cannot read gene " << gene.GetId() << " (bytes " << gene.GetStartByte() << "-" << gene.GetStartByte() + length
                              << ") from reference.fna: reference.map does not match the file" << std::endl;
                    exit(8);
                }
                auto packed = std::make_unique<uint8_t[]>(packed::Bytes(length));  // zeroed
                packed::Pack(m_scratch.data(), length, packed.get());
                gene.SetPacked(packed.get());
                m_owned.emplace_back(std::move(packed));
            } else if (!m_reader && length > 0) {
                errx(EX_SOFTWARE, "Gene %zu cannot be loaded on its own from a compressed reference (reference.fna.zst "
                                  "or a single-file database); the reference must be preloaded.", gene.GetId());
            }
        }

        void LoadGene(GeneKey key) {
            LoadGeneData(m_genes[GeneKeyToIndex(key)]);
        };

        // LoadGene for threads: one at a time.
        void LoadGeneOMP(GeneKey key) {
#pragma omp critical(genome_loader)
            LoadGeneData(m_genes[GeneKeyToIndex(key)]);
        };

        bool ValidGene(GeneKey key) {
            return GeneKeyToIndex(key) < m_genes.size();
        }

        // The gene list as long as `slots` genes (the largest gene id) at once, rather than as genes are added.
        void ReserveGenes(size_t slots) {
            if (m_genes.size() < slots) m_genes.resize(slots);
        }

        bool HasGene(GeneKey key) const {
            return key > 0 && GeneKeyToIndex(key) < m_genes.size() && m_genes[GeneKeyToIndex(key)].IsSet();
        }

        Gene& GetGene(GeneKey key) {
            return m_genes.at(GeneKeyToIndex(key));
        }

        const GeneList& GetGeneList() const {
            return m_genes;
        }

        GeneList& Genes() {
            return m_genes;
        }

        void MarkLoaded() {
            m_is_loaded.Set();
        }

        // Forgets the genes' sequences and frees those loaded one by one (GenomeLoader::ReleaseGeneSequences frees
        // the preload's arenas). A loaded genome stays marked loaded (GetGeneOMP returns its genes, without
        // sequences); reading a gene again stops protal.
        void ReleaseSequences() {
            for (auto& gene : m_genes) gene.ReleaseSequence();
            m_owned.clear();
            m_owned.shrink_to_fit();
            m_released = true;
        }

        size_t GeneNum() const {
            return !m_hittable_known ? std::count_if(m_genes.begin(), m_genes.end(), [](Gene const& gene) {
                return gene.IsSet();
            }) : m_hittable_genes.size();
        }

        void LoadGenome() {
            for (auto& gene : m_genes) LoadGeneData(gene);
            m_is_loaded.Set();
        };

        bool IsLoaded() const {
            return m_is_loaded.Get();
//            return std::any_of(m_genes.begin(), m_genes.end(), [](Gene const& gene){ return gene.IsLoaded(); });
        }

        void LoadGenomeOMP() {
#pragma omp critical(genome_loader)
            {
                if (!IsLoaded()) {
                    for (auto& gene : m_genes) LoadGeneData(gene);
                    m_is_loaded.Set();
                }
            }
        };

        Gene& GetGeneOMP(GeneKey key) {
            if (!m_is_loaded.Get()) {
#pragma omp critical(genome_loader)
                if (!m_is_loaded.Get()) {
                    for (auto& gene : m_genes) LoadGeneData(gene);
                    m_is_loaded.Set();
                }
            }
            return m_genes.at(GeneKeyToIndex(key));
        };
    };


    // Copies reference bytes from zstd::ParallelRead into the genes they belong to, packed (two bits per
    // base, see PackedSequence.h). The genes' packed bytes must be zero (a calloc'd arena), as
    // packed::PackInto needs. genes are sorted by start byte and do not overlap; starts[i] is genes[i]'s
    // start byte.
    class GeneSink : public zstd::Sink {
    public:
        GeneSink(std::vector<Gene*> const& genes, std::vector<uint64_t> const& starts) : m_genes(genes), m_starts(starts) {}

        void Copy(uint64_t offset, char const* data, size_t size) override {
            uint64_t const end = offset + size;
            // The last gene starting at or before offset may reach into the range.
            size_t i = static_cast<size_t>(std::upper_bound(m_starts.begin(), m_starts.end(), offset) - m_starts.begin());
            if (i > 0) i--;
            for (; i < m_genes.size() && m_starts[i] < end; i++) {
                uint64_t const begin = m_starts[i], gene_end = begin + m_genes[i]->GetLength();
                uint64_t const from = std::max(offset, begin), to = std::min(end, gene_end);
                if (from >= to) continue;
                packed::PackInto(m_genes[i]->MutablePacked(), from - begin, data + (from - offset), to - from);
            }
        }

    private:
        std::vector<Gene*> const& m_genes;
        std::vector<uint64_t> const& m_starts;
    };

    class GenomeLoader {
        using GenomeKey = Genome::GenomeKey;
        using GenomeMap = tsl::sparse_map<GenomeKey, Genome>;
        using GeneKey = Genome::GeneKey;


        db::DbFile m_reference;
        db::DbFile m_map;
        bool m_compressed = false;  // reference.fna.zst or in database.protal: genes are only read by LoadAllGenomes
        std::ifstream m_is;
        GenomeMap m_genomes;
        gene_conservation::Table m_gene_conservation;  // empty: every gene's factor is 1
        bool m_scale_depth_margin = false;  // the depth identity margin scaled by m_gene_conservation (--gene_conservation)
        gene_neighbours::Table m_gene_neighbours;      // empty: no gene's neighbours are known
        gene_incongruence::Table m_suspect_copies;     // empty: every gene copy is evidence of its species
        species_priors::Table m_species_priors;        // empty: every species' priors unknown

        int m_threads = 1;  // for reading reference.map
        gene_table::Times m_map_times, m_unique_times;  // the last loads of reference.map and unique_kmers.tsv (GeneTableTimes)
        gene_table::Times m_table_times;                // the load of gene_table.bin (GeneTableTimes)
        bool m_from_gene_table = false;         // the genes came from gene_table.bin (LoadGeneTable), not from reference.map
        bool m_unique_from_gene_table = false;  // and their unique k-mer counts too
        // reference.map's fingerprint (xxhash64 of its content, reference.fna's size), as ReferenceFingerprint::Of gives it:
        // hashed while it is parsed, or recorded in gene_table.bin; the index is checked against it.
        std::optional<ReferenceFingerprint> m_fingerprint;
        // The genomes in the order they were made (reference.map's order of first appearance), which gene_table.bin keeps (the same
        // tables give the same bytes) and its loader makes them in. No output depends on it, nor on the genome map's iteration
        // order (docs/claude/2026-10-05-order-independence).
        std::vector<GenomeKey> m_genome_order;
        struct FreeDeleter { void operator()(void* p) const { std::free(p); } };
        std::vector<std::unique_ptr<uint8_t[], FreeDeleter>> m_arenas;  // the preloaded genes' packed sequences (Gene::SetPacked)

        // Huge pages for a large arena, which the loading threads touch at random offsets (as
        // Seedmap::AdviseHugePages). Without transparent huge pages this does nothing.
        static void AdviseHugePages(void* data, size_t bytes) {
#ifdef MADV_HUGEPAGE
            constexpr size_t kMinBytes = size_t{64} << 20;
            if (bytes < kMinBytes) return;
            auto const page = static_cast<uintptr_t>(sysconf(_SC_PAGESIZE));
            auto const begin = (reinterpret_cast<uintptr_t>(data) + page - 1) & ~(page - 1);
            auto const end = (reinterpret_cast<uintptr_t>(data) + bytes) & ~(page - 1);
            if (end > begin) madvise(reinterpret_cast<void*>(begin), end - begin, MADV_HUGEPAGE);
#endif
        }

        Genome& AddOrGetGenome(GenomeKey const& key) {
            auto it = m_genomes.find(key);
            if (it == m_genomes.end()) it = m_genomes.insert({ key, Genome(key) }).first;
            return it.value();
        }

        void Open() {
            m_compressed = m_reference.Compressed();
            if (!m_compressed) m_is.open(m_reference.Path(), std::ios::in);
        }

    public:
        // reference: reference.fna, a zstd-compressed reference.fna.zst, or the member of a single-file
        // database. The byte offsets in map (reference.map) always refer to the uncompressed reference.
        // map is read with `threads` threads. gene_table: the database's gene_table.bin (GeneTableFile.h), loaded
        // instead of map where it was made from these tables (unique_size: the size of the database's
        // unique_kmers.tsv, nullopt if it has none); UniqueKmersFromGeneTable() tells whether it gave their counts too.
        GenomeLoader(db::DbFile reference, db::DbFile map, int threads = 1, std::optional<db::DbFile> const& gene_table = std::nullopt,
                     std::optional<uint64_t> unique_size = std::nullopt) :
                m_reference(std::move(reference)),
                m_map(std::move(map)),
                m_threads(threads) {
            Open();
            if (!gene_table || !LoadGeneTable(*gene_table, unique_size, m_threads)) LoadPositionMap(m_map, m_threads);
        };

        bool FromGeneTable() const { return m_from_gene_table; }
        bool UniqueKmersFromGeneTable() const { return m_unique_from_gene_table; }

        // reference.map's ReferenceFingerprint, known once the genes are loaded (nullopt for a loader that has none).
        std::optional<ReferenceFingerprint> const& MapFingerprint() const { return m_fingerprint; }

        // Loads the genes, and their unique k-mer counts if it has them, from the binary gene table (GeneTableFile.h) if it
        // was made from this database's tables: its recorded sizes of reference.map, unique_kmers.tsv (unique_size: the
        // database's, nullopt if it has none) and reference.fna are theirs. Each genome's genes are added by one thread.
        // Returns false, with nothing loaded, if the table is not this database's or not readable as one (a note says why).
        bool LoadGeneTable(db::DbFile const& table, std::optional<uint64_t> unique_size, int threads) {
            namespace gtf = gene_table_file;
            auto const start = std::chrono::steady_clock::now();
            auto const table_size = table.Size();
            auto const map_size = m_map.Size();
            auto const fna_size = m_reference.Size();
            auto skip = [&](std::string const& why) {
                std::cerr << "Note: " << table.Name() << " is not used (" << why << "); the genes are read from "
                          << m_map.Name() << std::endl;
                m_genomes.clear();
                return false;
            };
            if (!table.Exists() || !table_size || *table_size < sizeof(gtf::Header)) return skip("it cannot be read");
            if (!map_size || !fna_size) return false;  // LoadPositionMap says what is wrong
            gene_table::TableSink sink(*table_size);
            std::string error;
            uint64_t const read = table.ParallelRead(std::max(threads, 1), sink, error);
            if (!error.empty() || read != *table_size) return skip("it cannot be read: " + error);
            std::string_view const data = sink.Text(read);
            gtf::Header h;
            std::memcpy(&h, data.data(), sizeof(h));
            if (std::memcmp(h.magic, gtf::kMagic, sizeof(gtf::kMagic)) != 0) return skip("it is not a gene table");
            if (h.version != gtf::kVersion) return skip("it is of format version " + std::to_string(h.version));
            if (gtf::ExpectedSize(h) != read || h.genomes == 0) return skip("it is truncated or corrupt");
            if (h.map_size != *map_size || h.fna_size != *fna_size || (h.has_unique != 0) != unique_size.has_value() ||
                (unique_size && h.unique_size != *unique_size)) {
                return skip("it was made from other gene tables than the database's");
            }
            m_table_times = {};
            m_table_times.read = gene_table::SecondsSince(start);
            constexpr uint64_t max_id = (uint64_t{1} << SEEDMAP_TAXID_BITS) - 1;
            constexpr uint64_t max_gene = (uint64_t{1} << SEEDMAP_GENEID_BITS) - 1;
            constexpr uint64_t max_length = (uint64_t{1} << SEEDMAP_GENE_POS_BITS) - 1;
            char const* const genome_records = data.data() + sizeof(gtf::Header);
            char const* const gene_records = genome_records + h.genomes * sizeof(gtf::Genome);
            auto genome_at = [&](uint64_t i) { gtf::Genome g; std::memcpy(&g, genome_records + i * sizeof(g), sizeof(g)); return g; };
            auto gene_at = [&](uint64_t i) { gtf::Gene g; std::memcpy(&g, gene_records + i * sizeof(g), sizeof(g)); return g; };

            // The genomes, on this thread (the map is changed only here), then each one's genes on any thread.
            auto const genomes_start = std::chrono::steady_clock::now();
            m_genomes.clear();
            m_genome_order.clear();
            uint64_t next = 0;
            for (uint64_t i = 0; i < h.genomes; i++) {
                auto const g = genome_at(i);
                if (g.taxid == 0 || g.taxid > max_id || m_genomes.contains(g.taxid) || g.first != next || g.count == 0 ||
                    g.count > g.slots || g.slots > max_gene) {
                    return skip("its genome " + std::to_string(i + 1) + " is corrupt");
                }
                next += g.count;
                AddOrGetGenome(g.taxid);  // in reference.map's order, as LoadPositionMap makes them
                m_genome_order.push_back(g.taxid);
            }
            if (next != h.genes) return skip("its gene count is corrupt");
            std::vector<Genome*> genomes(h.genomes);
            for (uint64_t i = 0; i < h.genomes; i++) genomes[i] = &m_genomes.find(genome_at(i).taxid).value();  // stable now
            m_table_times.group = gene_table::SecondsSince(genomes_start);

            auto const add_start = std::chrono::steady_clock::now();
            bool const has_unique = h.has_unique != 0;
            std::atomic<bool> corrupt{ false };
            std::ifstream* const reader = m_compressed ? nullptr : &m_is;
            constexpr uint64_t kBlock = 256;
            zstd::ParallelFor((h.genomes + kBlock - 1) / kBlock, std::max(threads, 1), [&](size_t block, size_t) -> std::string {
                for (uint64_t i = block * kBlock; i < std::min<uint64_t>(h.genomes, (block + 1) * kBlock) && !corrupt; i++) {
                    auto const g = genome_at(i);
                    Genome& genome = *genomes[i];
                    genome.ReserveGenes(g.slots);
                    uint64_t last_id = 0;
                    for (uint64_t k = g.first; k < g.first + g.count; k++) {
                        auto const r = gene_at(k);
                        if (r.id == 0 || r.id > g.slots || r.id <= last_id || r.length == 0 || r.length > max_length ||
                            r.start + r.length > *fna_size) {
                            corrupt = true;
                            break;
                        }
                        last_id = r.id;
                        genome.AddGene(r.id, r.id, r.start, r.length, reader);
                        if (has_unique) {
                            genome.GetGene(r.id).SetUniqueValues(r.short_unique, r.long_unique, r.long_super_unique, r.total_kmers);
                            if (r.short_unique + r.long_unique > 0) genome.AddHittableGene(r.id);
                        }
                    }
                    if (has_unique) {
                        genome.SetHittableGenesKnown();
                        genome.SetUniqueValues();
                    }
                }
                return {};
            });
            if (corrupt) return skip("a gene record is corrupt");
            m_table_times.add = gene_table::SecondsSince(add_start);
            m_table_times.rows = h.genes;
            m_fingerprint = ReferenceFingerprint{ h.map_hash, h.fna_size };
            m_from_gene_table = true;
            m_unique_from_gene_table = has_unique;
            return true;
        }

        // Writes the genes as this loader holds them as a binary gene table (GeneTableFile.h) at path, recording the sizes
        // of the reference.map and unique_kmers.tsv they were read from (unique_size: 0 and has_unique false if none) and
        // reference.map's fingerprint (MapFingerprint, which a load from reference.map computes). Genomes in reference.map's
        // order, genes by id: the same tables give the same bytes. Returns an error message, empty on success.
        std::string WriteGeneTable(std::string const& path, uint64_t map_size, uint64_t unique_size, bool has_unique) const {
            namespace gtf = gene_table_file;
            if (!m_fingerprint) return "the loader has no reference fingerprint (its genes were not read from reference.map)";
            if (m_genome_order.size() != m_genomes.size()) return "the order the genomes were made in is not known";
            std::vector<GenomeKey> const& keys = m_genome_order;  // reference.map's order (see m_genome_order)
            std::vector<gtf::Genome> genomes;
            genomes.reserve(keys.size());
            uint64_t genes = 0;
            for (auto key : keys) {
                auto const& list = m_genomes.at(key).GetGeneList();
                uint64_t count = 0;
                for (auto const& gene : list) count += gene.IsSet();
                if (count == 0) continue;
                genomes.push_back({ key, list.size(), genes, count });
                genes += count;
            }
            gtf::Header h{};
            std::memcpy(h.magic, gtf::kMagic, sizeof(gtf::kMagic));
            h.version = gtf::kVersion;
            h.map_size = map_size;
            h.unique_size = has_unique ? unique_size : 0;
            h.has_unique = has_unique ? 1 : 0;
            h.map_hash = m_fingerprint->map_hash;
            h.fna_size = m_fingerprint->fna_size;
            h.genomes = genomes.size();
            h.genes = genes;
            std::ofstream os(path, std::ios::binary | std::ios::trunc);
            if (!os) return "cannot open " + path + " for writing";
            os.write(reinterpret_cast<char const*>(&h), sizeof(h));
            os.write(reinterpret_cast<char const*>(genomes.data()), static_cast<std::streamsize>(genomes.size() * sizeof(gtf::Genome)));
            std::vector<gtf::Gene> buffer;
            buffer.reserve(size_t{1} << 18);
            auto flush = [&]() {
                os.write(reinterpret_cast<char const*>(buffer.data()), static_cast<std::streamsize>(buffer.size() * sizeof(gtf::Gene)));
                buffer.clear();
            };
            for (auto const& g : genomes) {
                auto const& list = m_genomes.at(g.taxid).GetGeneList();
                for (auto const& gene : list) {
                    if (!gene.IsSet()) continue;
                    if (gene.IsLoaded()) return "the genes are loaded (their start bytes are not kept): write the table before LoadAllGenomes";
                    auto const [su, lu, lsu, total] = gene.GetUniqueKmerCounts();
                    buffer.push_back({ gene.GetStartByte(), static_cast<uint32_t>(gene.GetId()), static_cast<uint32_t>(gene.GetLength()),
                                       static_cast<uint32_t>(su), static_cast<uint32_t>(lu), static_cast<uint32_t>(lsu),
                                       static_cast<uint32_t>(total) });
                    if (buffer.size() == buffer.capacity()) flush();
                }
            }
            flush();
            os.close();
            if (os.fail()) return "write error on " + path;
            return {};
        }

        GenomeLoader(std::string genome_path, std::string genome_map) :
                GenomeLoader(db::DbFile::OnDisk(std::move(genome_path)), db::DbFile::OnDisk(std::move(genome_map))) {}

        GenomeLoader(const GenomeLoader& other) :
                m_reference(other.m_reference),
                m_map(other.m_map),
                m_gene_conservation(other.m_gene_conservation),
                m_scale_depth_margin(other.m_scale_depth_margin),
                m_gene_neighbours(other.m_gene_neighbours),
                m_suspect_copies(other.m_suspect_copies),
                m_species_priors(other.m_species_priors),
                m_threads(other.m_threads) {
            Open();
            LoadPositionMap(m_map, m_threads);
        }

        bool IsCompressed() const {
            return m_compressed;
        }

        // The genes' conservation factors (GeneConservation.h), which scale the depth identity margin
        // per gene; empty unless set, and then every factor is 1.
        gene_conservation::Table const& GetGeneConservation() const {
            return m_gene_conservation;
        }

        void SetGeneConservation(gene_conservation::Table table) {
            m_gene_conservation = std::move(table);
        }

        // Whether the depth identity margin is scaled per gene by the factors (--gene_conservation db or FILE); the
        // factors give the model's conservation features either way.
        bool ScaleDepthMargin() const {
            return m_scale_depth_margin;
        }

        void SetScaleDepthMargin(bool scale) {
            m_scale_depth_margin = scale;
        }

        // The gene copies that are no evidence of their species (GeneIncongruence.h: near-identical to another genus's
        // copy, contamination or transfer), whose records a run leaves out; empty unless set (a database without
        // suspect_copies.tsv, or --keep_suspect_copies).
        gene_incongruence::Table const& GetSuspectCopies() const {
            return m_suspect_copies;
        }

        void SetSuspectCopies(gene_incongruence::Table table) {
            m_suspect_copies = std::move(table);
        }

        bool IsSuspectCopy(uint32_t taxid, uint32_t geneid) const {
            return m_suspect_copies.Contains(taxid, geneid);
        }

        // What GTDB knows of each species before any read (SpeciesPriors.h: duplicated markers, CheckM quality, the
        // cluster's ANI radius and width), for the model's prior features; empty unless set (a database without
        // species_priors.tsv), and then every species' values are unknown.
        species_priors::Table const& GetSpeciesPriors() const {
            return m_species_priors;
        }

        void SetSpeciesPriors(species_priors::Table table) {
            m_species_priors = std::move(table);
        }

        // Which genes lie next to which in the species' clades (GeneNeighbours.h); empty unless set (a
        // database without gene_neighbours.tsv, or --no_gene_neighbours), and then nothing uses it.
        gene_neighbours::Table const& GetGeneNeighbours() const {
            return m_gene_neighbours;
        }

        void SetGeneNeighbours(gene_neighbours::Table table) {
            m_gene_neighbours = std::move(table);
        }

        ~GenomeLoader() {
            m_is.close();
        }

        void PrintHittableGenes() {
            for (auto& [gid, _] : m_genomes) {
                auto& genome = m_genomes.at(gid);
                auto hg = genome.GetHittableGenes();
                std::string hgstr = "";
                for (auto gene : hg) {
                    hgstr += std::to_string(gene) + ',';
                }
                std::cout << "Hittable\t" << gid << '\t' << hg.size() << '\t' << hgstr << std::endl;
            }
        }

        // unique_kmers.tsv (written by --build): taxid, gene id, then counts and rates of short,
        // long and long-super unique k-mers, and the gene's k-mer total. Every line must name a
        // gene of reference.map. Parsed with `threads` threads (gene_table).
        void LoadUniqueKmers(std::string const& file, int threads = 1) {
            LoadUniqueKmers(db::DbFile::OnDisk(file), threads);
        }

        // Where the gene tables' loads spent their time, as one line for the run's log.
        std::string GeneTableTimes() const {
            auto part = [](std::string const& name, gene_table::Times const& t, bool sums) {
                std::ostringstream os;
                os << std::fixed << std::setprecision(2) << name << " " << t.rows << " rows: reading " << t.read << " s, parsing " << t.parse
                   << " s, by genome " << t.group << " s, adding " << t.add << " s";
                if (sums) os << ", genome sums " << t.after << " s";
                return os.str();
            };
            if (m_from_gene_table) {
                std::ostringstream os;
                os << std::fixed << std::setprecision(2) << "Gene tables: " << gene_table_file::kFileName << " " << m_table_times.rows
                   << " genes" << (m_unique_from_gene_table ? " with unique k-mer counts" : "") << ": reading " << m_table_times.read
                   << " s, genomes " << m_table_times.group << " s, adding " << m_table_times.add << " s";
                return os.str();
            }
            return "Gene tables: " + part("reference.map", m_map_times, false) + "; " + part("unique_kmers.tsv", m_unique_times, true);
        }

        void LoadUniqueKmers(db::DbFile const& unique_kmers, int threads = 1) {
            std::string const& file = unique_kmers.Name();
            struct Row {
                uint64_t taxid, geneid, short_unique, long_unique, long_super_unique, total;
                size_t line;
            };
            auto parse = [](std::string_view line, Row& row, std::string& problem) {
                gene_table::Fields<9> const fields(line);
                if (fields.n < 9) {
                    problem = "expected 9 tab-separated columns, found " + std::to_string(fields.n);
                    return false;
                }
                uint64_t counts[6];
                size_t const columns[6] = { 0, 1, 2, 4, 6, 8 };  // taxid, gene, short, long, long-super, total
                for (int i = 0; i < 6; i++) {
                    auto const number = gene_table::Number(fields.f[columns[i]], counts[i]);
                    if (number != gene_table::NumberProblem::None) {
                        problem = gene_table::ColumnProblem(columns[i] + 1, number);
                        return false;
                    }
                }
                // A gene's counts are kept in 32 bits (a gene has fewer than 2^20 bases).
                for (int i = 2; i < 6; i++) {
                    if (counts[i] > UINT32_MAX) {
                        problem = gene_table::ColumnProblem(columns[i] + 1, gene_table::NumberProblem::TooLarge);
                        return false;
                    }
                }
                row = { counts[0], counts[1], counts[2], counts[3], counts[4], counts[5], 0 };
                return true;
            };
            size_t rows = 0;
            // Each genome's rows by one thread (gene_table::AddPieceByGenome); the genome map is only read.
            auto add_piece = [&](std::vector<gene_table::Chunk<Row>> const& chunks, std::vector<size_t> const& line_bases) {
                for (auto const& chunk : chunks) rows += chunk.rows.size();
                auto const [line, problem] = gene_table::AddPieceByGenome<Row>(chunks, line_bases, threads,
                    [](uint64_t) {},
                    [&](uint64_t taxid, Row const& row) -> std::string {
                        auto it = m_genomes.find(taxid);
                        if (it == m_genomes.end() || !it.value().HasGene(row.geneid)) {
                            return "gene " + std::to_string(row.taxid) + "_" + std::to_string(row.geneid) +
                                   " is not in reference.map (rebuild the database with --build)";
                        }
                        auto& taxon = it.value();
                        if (row.short_unique + row.long_unique > 0) taxon.AddHittableGene(row.geneid);
                        taxon.GetGene(row.geneid).SetUniqueValues(row.short_unique, row.long_unique, row.long_super_unique, row.total);
                        return {};
                    }, &m_unique_times);
                if (line > 0) InvalidUniqueKmers(file, line, problem);
            };
            m_unique_times = {};
            std::string const error = gene_table::ForEachPiece<Row>(unique_kmers, threads, parse, add_piece, &m_unique_times);
            if (!error.empty()) InvalidUniqueKmers(file, 0, error);
            // Without rows every taxon would fail the model silently.
            if (rows == 0) InvalidUniqueKmers(file, 0, "the file lists no genes (rebuild the database with --build)");

            // Every genome's sums over its genes, on all threads (one pass over 24M genes at GTDB r226 size).
            auto const after = std::chrono::steady_clock::now();
            std::vector<Genome*> genomes;
            genomes.reserve(m_genomes.size());
            for (auto it = m_genomes.begin(); it != m_genomes.end(); ++it) genomes.push_back(&it.value());
            zstd::ParallelFor((genomes.size() + 1023) / 1024, threads, [&](size_t block, size_t) -> std::string {
                for (size_t i = block * 1024; i < std::min(genomes.size(), (block + 1) * 1024); i++) {
                    genomes[i]->SetHittableGenesKnown();
                    genomes[i]->SetUniqueValues();
                }
                return {};
            });
            m_unique_times.after = gene_table::SecondsSince(after);

            // PrintHittableGenes();
        }

        void LoadHittableGenes(std::string const& file) {
            std::ifstream is(file, std::ios::in);

            std::vector<std::string> tokens;
            std::string line;
            while (std::getline(is, line)) {
                Utils::split(tokens, line, "\t");
                auto taxid = std::stoull(tokens[0]);
                auto gene_str = tokens[2];

                auto& taxon = m_genomes.at(taxid);

                Utils::split(tokens, gene_str, ",");
                for (auto const& g : tokens) {
                    auto gid = std::stoul(g);
                    taxon.AddHittableGene(gid);
                }
            }
            is.close();
            for (auto& [key, _] : m_genomes) m_genomes.at(key).SetHittableGenesKnown();
        }

        size_t GetLoadedGenomeCount() const {
            return std::count_if(m_genomes.begin(), m_genomes.end(), [](auto const& pair) { return pair.second.IsLoaded(); });
        }

        const Genome& GetGenome(GenomeKey const& key) const {
            return m_genomes.at(key);
        }

        PROTAL_CLONE_V3 bool HasGene(GenomeKey taxid, GeneKey gene) const {
            auto it = m_genomes.find(taxid);
            return it != m_genomes.end() && it->second.HasGene(gene);
        }

        // Length of a gene as reference.map gives it (0 if the gene is not in the map).
        PROTAL_CLONE_V3 size_t GeneLength(GenomeKey taxid, GeneKey gene) {
            return HasGene(taxid, gene) ? m_genomes.at(taxid).GetGene(gene).GetLength() : 0;
        }

        size_t GeneCount() const {
            size_t n = 0;
            for (auto const& [key, genome] : m_genomes) {
                for (auto const& gene : genome.GetGeneList()) n += gene.IsSet();
            }
            return n;
        }

        // The largest taxid (genome key), gene id and gene length of the reference: what the index's
        // fields hold, and so the widths a query run keeps them in (Seedmap::PackedLayout).
        std::tuple<uint64_t, uint64_t, uint64_t> IndexFieldMaxima() const {
            uint64_t taxid = 0, gene = 0, length = 0;
            for (auto const& [key, genome] : m_genomes) {
                taxid = std::max<uint64_t>(taxid, key);
                for (auto const& g : genome.GetGeneList()) {
                    if (!g.IsSet()) continue;
                    gene = std::max<uint64_t>(gene, g.GetId());
                    length = std::max<uint64_t>(length, g.GetLength());
                }
            }
            return { taxid, gene, length };
        }

        // Length of the longest gene of the reference.
        size_t MaxGeneLength() const {
            size_t longest = 0;
            for (auto const& [key, genome] : m_genomes) {
                for (auto const& gene : genome.GetGeneList()) {
                    if (gene.IsSet()) longest = std::max(longest, gene.GetLength());
                }
            }
            return longest;
        }

        // The genome map iterates in an order of its own (its buckets'), which nothing should depend on: code whose
        // result depends on the order of the genomes goes over SortedKeys() instead.
        GenomeMap& GetGenomeMap() {
            return m_genomes;
        }

        // The genomes' taxids in ascending order.
        std::vector<GenomeKey> SortedKeys() const {
            std::vector<GenomeKey> keys;
            keys.reserve(m_genomes.size());
            for (auto const& [key, _] : m_genomes) keys.push_back(key);
            std::sort(keys.begin(), keys.end());
            return keys;
        }

        // The header for every gene of the database (--full_sam_header), by taxid and gene id.
        void WriteSamHeader(std::ostream& os=std::cout) {
            os << "@HD\tVN:1.6\n";
            for (auto const key : SortedKeys()) {
                auto& genes = m_genomes.at(key).GetGeneList();
                for (auto i = 0; i < genes.size(); i++) {
                    if (genes[i].IsSet()) {

                        os << "@SQ\tSN:" << key << '_' << genes[i].GetId() << '\t' << "LN:" << genes[i].GetLength() << '\n';
                    }
                }
            }
        }

        // The header for the given genes only (taxid << 32 | gene id, sorted), in that order: SAM
        // needs @SQ lines only for the references that records name. Genes not in the database are
        // left out.
        void WriteSamHeader(std::ostream& os, std::vector<uint64_t> const& genes) {
            os << "@HD\tVN:1.6\n";
            std::string line;
            for (uint64_t const key : genes) {
                uint64_t const taxid = key >> 32, geneid = key & 0xffffffffu;
                if (!HasGene(taxid, geneid)) continue;
                line = "@SQ\tSN:";
                line += std::to_string(taxid);
                line += '_';
                line += std::to_string(geneid);
                line += "\tLN:";
                line += std::to_string(GeneLength(taxid, geneid));
                line += '\n';
                os << line;
            }
        }

        // Reads every gene of the reference (raw or compressed) with up to `threads` threads:
        // zstd::ParallelRead delivers the file's bytes (raw and seekable zstd files in parallel
        // chunks, other zstd files in one stream), and each piece is copied, uppercased, into the
        // genes it overlaps. Large reads instead of one seek per gene, which matters on network storage.
        void LoadAllGenomes(int threads = 1) {
            auto const keys = SortedKeys();

            std::vector<Gene*> genes;
            for (auto& key : keys) {
                auto& genome = m_genomes.at(key);
                if (genome.IsLoaded()) continue;
                for (auto& gene : genome.Genes()) {
                    if (gene.IsSet() && !gene.IsLoaded() && gene.GetLength() > 0) genes.emplace_back(&gene);
                }
            }
            std::sort(genes.begin(), genes.end(), [](Gene const* a, Gene const* b) {
                return a->GetStartByte() < b->GetStartByte();
            });

            // Give every gene its place in one arena of packed sequences first (sizing 16.6M strings one
            // by one took seconds on one thread); the threads then fill disjoint parts of it (a gene that
            // spans two chunks gets its two parts from two threads). Each gene starts on a byte. The
            // arena is zero, as PackInto needs, without being touched (calloc hands out fresh pages for
            // a block this large); the reference covers every gene byte, or protal stops below.
            std::vector<uint64_t> starts;
            starts.reserve(genes.size());
            uint64_t position = 0, packed_bytes = 0;
            for (Gene* gene : genes) {
                if (gene->GetStartByte() < position) {
                    std::cerr << "Invalid reference map " << m_map.Name() << ": gene " << gene->GetId() << " at byte "
                              << gene->GetStartByte() << " overlaps the previous gene" << std::endl;
                    exit(8);
                }
                position = gene->GetStartByte() + gene->GetLength();
                packed_bytes += packed::Bytes(gene->GetLength());
                starts.emplace_back(gene->GetStartByte());
            }
            if (!genes.empty()) {
                std::unique_ptr<uint8_t[], FreeDeleter> arena(static_cast<uint8_t*>(std::calloc(packed_bytes, 1)));
                if (!arena) {
                    std::cerr << "Cannot allocate " << packed_bytes << " bytes for the reference genes" << std::endl;
                    exit(8);
                }
                AdviseHugePages(arena.get(), packed_bytes);
                uint64_t offset = 0;
                for (Gene* gene : genes) {
                    gene->SetPacked(arena.get() + offset);
                    offset += packed::Bytes(gene->GetLength());
                }
                m_arenas.emplace_back(std::move(arena));
            }
            GeneSink sink(genes, starts);
            std::string error;
            uint64_t const size = m_reference.ParallelRead(threads, sink, error);
            if (!error.empty()) {
                std::cerr << "Cannot read the reference " << m_reference.Name() << ": " << error << std::endl;
                exit(8);
            }
            if (size < position) {
                std::cerr << "Cannot read all genes from " << m_reference.Name() << ": it holds " << size << " bytes, "
                          << m_map.Name() << " lists genes up to byte " << position << std::endl;
                exit(8);
            }
            for (auto& key : keys) {
                m_genomes.at(key).MarkLoaded();
            }
        }

        bool AllGenomesLoaded() {
            for (auto& [id, genome] : m_genomes) {
                if (!genome.IsLoaded()){
                    return false;
                }
            }
            return true;
        }

        // Frees every gene's sequence, preloaded or loaded one by one: for the end of --build, which reads no gene
        // after its uniqueness check. The genes keep their ids, lengths and k-mer counts; their sequences are empty
        // from here, and reading one again stops protal (a loaded gene no longer knows its start byte).
        void ReleaseGeneSequences() {
            for (auto it = m_genomes.begin(); it != m_genomes.end(); ++it) it.value().ReleaseSequences();
            m_arenas.clear();
            m_arenas.shrink_to_fit();
        }

        PROTAL_CLONE_V3 Genome& GetGenome(GenomeKey const& key) {
            assert(m_genomes.contains(key));
            if (!m_genomes.contains(key)) {
                std::cout << "Genomes Key: " << key << std::endl;
                exit(10);
            }
            return m_genomes.at(key);
        }


        // reference.map: taxid, gene id, start byte, end byte of the gene's sequence in reference.fna.
        // Every line is checked, so a bad map stops protal here instead of corrupting genes silently.
        // Parsed with `threads` threads (gene_table).
        void LoadPositionMap(db::DbFile const& map, int threads = 1) {
            std::string const& file_path = map.Name();
            if (!map.Exists()) InvalidMap(file_path, 0, "cannot open the file");
            // Offsets refer to the uncompressed reference, also for reference.fna.zst.
            auto const size = m_reference.Size();
            if (!size) InvalidMap(file_path, 0, "cannot read the size of " + m_reference.Name());
            uint64_t const fna_size = *size;
            struct Row {
                uint64_t taxid, geneid, start, end;
                size_t line;
            };
            std::string const& reference_name = m_reference.Name();
            auto parse = [&](std::string_view line, Row& row, std::string& problem) {
                constexpr uint64_t max_id = (uint64_t{1} << SEEDMAP_TAXID_BITS) - 1;
                constexpr uint64_t max_gene = (uint64_t{1} << SEEDMAP_GENEID_BITS) - 1;
                constexpr uint64_t max_length = (uint64_t{1} << SEEDMAP_GENE_POS_BITS) - 1;
                gene_table::Fields<4> const fields(line);
                if (fields.n != 4) {
                    problem = "expected 4 tab-separated columns, found " + std::to_string(fields.n);
                    return false;
                }
                uint64_t values[4];
                for (int i = 0; i < 4; i++) {
                    auto const number = gene_table::Number(fields.f[i], values[i]);
                    if (number != gene_table::NumberProblem::None) {
                        problem = gene_table::ColumnProblem(i + 1, number);
                        return false;
                    }
                }
                uint64_t const genome_id = values[0], gene_key = values[1], start = values[2], end = values[3];
                // taxid 0 marks empty index entries and gene ids are 1-based; both are packed into 20 bits.
                if (genome_id == 0 || genome_id > max_id) problem = "taxid must be between 1 and " + std::to_string(max_id);
                else if (gene_key == 0 || gene_key > max_gene) problem = "gene id must be between 1 and " + std::to_string(max_gene);
                else if (end <= start) problem = "end byte must be after start byte";
                else if (end > fna_size) problem = "end byte " + std::to_string(end) + " is past the end of " + reference_name + " (" + std::to_string(fna_size) + " bytes)";
                else if (end - start > max_length) problem = "gene is longer than " + std::to_string(max_length) + " bases";
                if (!problem.empty()) return false;
                row = { genome_id, gene_key, start, end, 0 };
                return true;
            };
            size_t rows = 0;
            // The genomes of a piece are made first, on this thread; then each genome's genes are added by one thread
            // (gene_table::AddPieceByGenome; the genome map is only read in parallel, and a genome's gene list is its own).
            auto add_piece = [&](std::vector<gene_table::Chunk<Row>> const& chunks, std::vector<size_t> const& line_bases) {
                for (auto const& chunk : chunks) rows += chunk.rows.size();
                auto const [line, problem] = gene_table::AddPieceByGenome<Row>(chunks, line_bases, threads,
                    [&](uint64_t taxid) {
                        if (m_genomes.contains(taxid)) return;  // listed again further on: made in an earlier piece
                        AddOrGetGenome(taxid);
                        m_genome_order.push_back(taxid);
                    },
                    [&](uint64_t taxid, Row const& row) -> std::string {
                        auto& genome = m_genomes.find(taxid).value();
                        if (genome.HasGene(row.geneid)) return "gene " + std::to_string(row.taxid) + "_" + std::to_string(row.geneid) + " is listed twice";
                        genome.AddGene(row.geneid, row.geneid, row.start, row.end - row.start, m_compressed ? nullptr : &m_is);
                        return {};
                    }, &m_map_times);
                if (line > 0) InvalidMap(file_path, line, problem);
            };
            XXHash64 hash(0);
            m_map_times = {};
            std::string const error = gene_table::ForEachPiece<Row>(map, threads, parse, add_piece, &m_map_times, &hash);
            if (!error.empty()) InvalidMap(file_path, 0, error);
            if (rows == 0) InvalidMap(file_path, 0, "the file lists no genes");
            m_fingerprint = ReferenceFingerprint{ hash.hash(), fna_size };
        }

        [[noreturn]] static void InvalidUniqueKmers(std::string const& path, size_t line_no, std::string const& reason) {
            std::cerr << "Invalid unique k-mer file " << path;
            if (line_no > 0) std::cerr << ", line " << line_no;
            std::cerr << ": " << reason << std::endl;
            exit(8);
        }

        [[noreturn]] static void InvalidMap(std::string const& path, size_t line_no, std::string const& reason) {
            std::cerr << "Invalid reference map " << path;
            if (line_no > 0) std::cerr << ", line " << line_no;
            std::cerr << ": " << reason << std::endl;
            exit(8);
        }
    };
}
