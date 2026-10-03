//
// Created by fritsche on 15/07/22.
//

#pragma once

#include <Constants.h>
#include "TargetClones.h"
#include <algorithm>
#include <atomic>
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
#include <fstream>
#include <err.h>
#include <KmerUtils.h>


#include <sparse_set.h>

#include "Utilities.h"
#include "Zstd.h"
#include "Database.h"
#include <sysexits.h>

#include "Benchmark.h"
#include "PackedSequence.h"
#include "GeneConservation.h"
#include "GeneIncongruence.h"
#include "GeneNeighbours.h"

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

        // Reads a gene table in pieces, parses each in parallel (ParsePiece) and hands its chunks, in
        // file order, to add(chunk, lines before the chunk). Returns an error message if the file
        // cannot be read or decompressed, else empty.
        template<typename Row, typename Parse, typename Add>
        std::string ForEachChunk(db::DbFile const& file, int threads, Parse&& parse, Add&& add) {
            auto input = file.Open();
            if (!file.Exists() || !input->IsOpen()) return "cannot open the file";
            std::istream& is = input->Stream();
            std::string buffer;
            std::vector<Chunk<Row>> chunks;
            size_t kept = 0;         // bytes of an unfinished last line carried over to the next piece
            size_t line_base = 0;    // lines before this piece
            size_t piece = kPieceBytes;
            bool end = false;
            while (!end) {
                buffer.resize(kept + piece);
                is.read(buffer.data() + kept, static_cast<std::streamsize>(piece));
                size_t const got = static_cast<size_t>(is.gcount());
                if (is.bad()) return "the file cannot be read or decompressed (truncated or corrupt file?)";
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
                ParsePiece<Row>(std::string_view(buffer.data(), cut), threads, parse, chunks);
                for (auto const& chunk : chunks) {
                    add(chunk, line_base);
                    line_base += chunk.lines;
                }
                std::memmove(buffer.data(), buffer.data() + cut, size - cut);
                kept = size - cut;
                piece = kPieceBytes;
            }
            return {};
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

        int m_threads = 1;  // for reading reference.map
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
        // map is read with `threads` threads.
        GenomeLoader(db::DbFile reference, db::DbFile map, int threads = 1) :
                m_reference(std::move(reference)),
                m_map(std::move(map)),
                m_threads(threads) {
            Open();
            LoadPositionMap(m_map, m_threads);
        };

        GenomeLoader(std::string genome_path, std::string genome_map) :
                GenomeLoader(db::DbFile::OnDisk(std::move(genome_path)), db::DbFile::OnDisk(std::move(genome_map))) {}

        GenomeLoader(const GenomeLoader& other) :
                m_reference(other.m_reference),
                m_map(other.m_map),
                m_gene_conservation(other.m_gene_conservation),
                m_scale_depth_margin(other.m_scale_depth_margin),
                m_gene_neighbours(other.m_gene_neighbours),
                m_suspect_copies(other.m_suspect_copies),
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
            auto add = [&](gene_table::Chunk<Row> const& chunk, size_t line_base) {
                rows += chunk.rows.size();
                for (auto const& row : chunk.rows) {
                    auto it = m_genomes.find(row.taxid);
                    if (it == m_genomes.end() || !it->second.HasGene(row.geneid)) {
                        InvalidUniqueKmers(file, line_base + row.line, "gene " + std::to_string(row.taxid) + "_" + std::to_string(row.geneid) +
                                                                       " is not in reference.map (rebuild the database with --build)");
                    }
                    auto& taxon = it.value();
                    if (row.short_unique + row.long_unique > 0) taxon.AddHittableGene(row.geneid);
                    taxon.GetGene(row.geneid).SetUniqueValues(row.short_unique, row.long_unique, row.long_super_unique, row.total);
                }
                if (!chunk.problem.empty()) InvalidUniqueKmers(file, line_base + chunk.problem_line, chunk.problem);
            };
            std::string const error = gene_table::ForEachChunk<Row>(unique_kmers, threads, parse, add);
            if (!error.empty()) InvalidUniqueKmers(file, 0, error);
            // Without rows every taxon would fail the model silently.
            if (rows == 0) InvalidUniqueKmers(file, 0, "the file lists no genes (rebuild the database with --build)");

            for (auto it = m_genomes.begin(); it != m_genomes.end(); ++it) {
                it.value().SetHittableGenesKnown();
                it.value().SetUniqueValues();
            }

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

        GenomeMap& GetGenomeMap() {
            return m_genomes;
        }

        void WriteSamHeader(std::ostream& os=std::cout) {
            os << "@HD\tVN:1.6\n";
            for (auto& [key, genome] : m_genomes) {
                auto& genes = genome.GetGeneList();
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
            std::vector<GenomeKey> keys;
            for (auto& pair : m_genomes) {
                keys.emplace_back(pair.first);
            }
            std::sort(keys.begin(), keys.end());

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
            auto add = [&](gene_table::Chunk<Row> const& chunk, size_t line_base) {
                rows += chunk.rows.size();
                for (auto const& row : chunk.rows) {
                    auto& genome = AddOrGetGenome(row.taxid);
                    if (genome.HasGene(row.geneid)) {
                        InvalidMap(file_path, line_base + row.line, "gene " + std::to_string(row.taxid) + "_" + std::to_string(row.geneid) + " is listed twice");
                    }
                    genome.AddGene(row.geneid, row.geneid, row.start, row.end - row.start, m_compressed ? nullptr : &m_is);
                }
                if (!chunk.problem.empty()) InvalidMap(file_path, line_base + chunk.problem_line, chunk.problem);
            };
            std::string const error = gene_table::ForEachChunk<Row>(map, threads, parse, add);
            if (!error.empty()) InvalidMap(file_path, 0, error);
            if (rows == 0) InvalidMap(file_path, 0, "the file lists no genes");
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
