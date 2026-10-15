// UniqueCores.h - a record's unique k-mers counted from its alignment, whatever share of the index's cores was loaded or
// looked up (docs/claude/2026-10-14-length-free-features, follow-ups).
//
// ZU and ZT counted the unique seeds of the anchor a record was aligned on, and ZS that anchor's seeds: on a grid's
// sparse index (the cores of one grid level loaded) a read finds about a quarter of them at C 4, and a short read and a
// split anchor find fewer too. A database --build made since 2026-10-14 holds unique_cores.bin: for every gene copy, in
// the order of its cores (the k-mer handler's windows over the copy's 2-bit sequence, as a run computes them), two bits:
// unique (the index entry's flag 60: its 31-mer under no other taxon and only once) and unique at distance two (flag 61
// for a core of several values, flag 60 for a core of one, as a seed counts them, KmerLookup.h). A run counts, for each
// record, the gene's cores whose whole 31-mer lies inside an exact stretch of the alignment (M columns without an
// insertion between them): those the whole index's seeds find, since a seed needs its 31-mer exact. ZS is their number
// and ZU and ZT the unique ones among them, the same for a run on the whole index, on a sparse index or with sparse
// seeding; the anchor's own seeds no longer count. A core of the 2-bit sequence the build dropped (a window over an
// ambiguous base) or whose key the index left out (too many values) has no entry and is not unique.
#pragma once

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <memory>
#include <optional>
#include <ostream>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Constants.h"

namespace protal::unique_cores {
    inline const std::string kFileName = "unique_cores.bin";
    inline const std::string kMagic = "protal unique cores 1\n";
    inline constexpr int64_t kFlank = 8;   // a core's first base lies kFlank bases into its k-mer
    inline constexpr int64_t kKmer = 31;   // the k-mer a seed matches whole
    inline constexpr uint8_t kUnique = 1;     // a core's bits: unique
    inline constexpr uint8_t kUniqueTwo = 2;  // unique at distance two
    inline constexpr size_t kCachedGenes = 4096;  // genes whose core positions a counter keeps

    struct Counts {
        uint32_t cores = 0;       // the gene's cores the alignment covers with its exact k-mers
        uint32_t unique = 0;      // of them, unique
        uint32_t unique_two = 0;  // unique at distance two
    };

    // The cores among `positions` (a gene's cores, their first bases, ascending) whose k-mer [p - kFlank, p - kFlank +
    // kKmer) lies inside an exact stretch of an alignment starting at gene position `start`: `columns` holds one operation
    // per column (AlignmentInfo::cigar: M a match, X a mismatch, I a read base, D a gene base, S/H clipped read bases),
    // and a stretch is a run of M columns with no other column between them. `flags(i)` gives core i's bits.
    template<typename Flags>
    inline Counts CountCovered(std::span<uint32_t const> positions, Flags&& flags, std::string_view columns, int64_t start) {
        Counts counts;
        int64_t r = start;  // the gene position of the next column that takes one
        int64_t run = -1;   // the first gene position of the current exact stretch; -1: none
        auto close = [&]() {
            if (run < 0) return;
            int64_t const lo = run + kFlank, hi = r - (kKmer - kFlank);  // cores p: run <= p - kFlank, p - kFlank + kKmer <= r
            run = -1;
            if (hi < lo || hi < 0) return;
            auto first = std::lower_bound(positions.begin(), positions.end(), static_cast<uint32_t>(std::max<int64_t>(lo, 0)));
            auto last = std::upper_bound(first, positions.end(), static_cast<uint32_t>(hi));
            for (auto it = first; it != last; ++it) {
                uint8_t const f = flags(static_cast<size_t>(it - positions.begin()));
                counts.cores++;
                counts.unique += (f & kUnique) != 0;
                counts.unique_two += (f & kUniqueTwo) != 0;
            }
        };
        for (char const op : columns) {
            switch (op) {
                case 'M': case '=':
                    if (run < 0) run = r;
                    r++;
                    break;
                case 'X': case 'D': case 'N':
                    close();
                    r++;
                    break;
                default:  // I, S, H: no gene base, but the stretch ends
                    close();
                    break;
            }
        }
        close();
        return counts;
    }

    // The bits of every gene copy's cores. Rows are the copies, row = first_row[taxid] + gene - 1 (Build.h GeneRowsOf);
    // a row's cores start at a byte (four cores a byte, two bits each, the first core in the lowest bits), so a row holds
    // its cores and up to three of padding.
    //
    // The file: kMagic, then little-endian 64-bit integers: the taxa (first_row's size), the rows and the flag bytes;
    // first_row (32-bit each), the rows' first bytes (rows + 1 of them, 32-bit), the flag bytes.
    class Table {
    public:
        Table() = default;

        // An empty table for these rows (first_row: GeneRows::first_row), each with room for its cores (`cores`, by row).
        Table(std::vector<uint64_t> const& first_row, std::vector<uint32_t> const& cores) {
            m_first_row.assign(first_row.begin(), first_row.end());
            m_start.assign(cores.size() + 1, 0);
            uint64_t bytes = 0;
            for (size_t row = 0; row < cores.size(); row++) {
                m_start[row] = static_cast<uint32_t>(bytes);
                bytes += (static_cast<uint64_t>(cores[row]) + 3) / 4;
            }
            m_start[cores.size()] = static_cast<uint32_t>(bytes);
            m_flags.assign(static_cast<size_t>(bytes), char{ 0 });
        }

        bool Empty() const { return m_start.size() < 2; }
        size_t Rows() const { return m_start.empty() ? 0 : m_start.size() - 1; }
        size_t Bytes() const { return m_first_row.size() * 4 + m_start.size() * 4 + m_flags.size(); }
        // Whether the rows fit 32-bit byte offsets (the file's), checked before filling a table.
        static bool Fits(std::vector<uint32_t> const& cores) {
            uint64_t bytes = 0;
            for (auto const c : cores) bytes += (static_cast<uint64_t>(c) + 3) / 4;
            return bytes <= UINT32_MAX && cores.size() < UINT32_MAX;
        }

        // The row of taxid's gene, or nullopt beyond the table.
        std::optional<size_t> Row(uint32_t taxid, uint32_t gene) const {
            if (gene == 0 || static_cast<size_t>(taxid) + 1 >= m_first_row.size()) return std::nullopt;
            uint64_t const row = static_cast<uint64_t>(m_first_row[taxid]) + gene - 1;
            if (row >= m_first_row[taxid + 1]) return std::nullopt;
            return static_cast<size_t>(row);
        }

        // Each taxon's first row (by taxid; one more than the taxa, the last the rows): RowsMismatch.
        std::vector<uint32_t> const& FirstRows() const { return m_first_row; }

        // The cores a row has room for (its cores and up to three more).
        size_t Capacity(size_t row) const { return 4 * static_cast<size_t>(m_start[row + 1] - m_start[row]); }

        uint8_t Flags(size_t row, size_t core) const {
            return static_cast<uint8_t>((static_cast<uint8_t>(m_flags[m_start[row] + core / 4]) >> (2 * (core % 4))) & 3u);
        }

        // Sets core `core` of row `row` (the build's; each row's bytes are its own, so rows may be filled on threads).
        void Set(size_t row, size_t core, uint8_t flags) {
            auto& byte = m_flags[m_start[row] + core / 4];
            auto const old = static_cast<uint8_t>(byte);
            byte = static_cast<char>((old & ~(3u << (2 * (core % 4)))) | ((flags & 3u) << (2 * (core % 4))));
        }

        // The cores with each bit set, over every row.
        std::pair<uint64_t, uint64_t> Totals() const {
            uint64_t unique = 0, two = 0;
            for (char const c : m_flags) {
                auto const byte = static_cast<uint8_t>(c);
                for (int i = 0; i < 4; i++) {
                    unique += (byte >> (2 * i)) & kUnique;
                    two += ((byte >> (2 * i)) & kUniqueTwo) != 0;
                }
            }
            return { unique, two };
        }

        void Write(std::ostream& os) const {
            std::string head = kMagic;
            auto put64 = [&](uint64_t v) { for (int i = 0; i < 8; i++) head.push_back(static_cast<char>((v >> (8 * i)) & 0xff)); };
            auto put32 = [&](uint32_t v) { for (int i = 0; i < 4; i++) head.push_back(static_cast<char>((v >> (8 * i)) & 0xff)); };
            put64(m_first_row.size());
            put64(Rows());
            put64(m_flags.size());
            for (auto const v : m_first_row) put32(v);
            for (auto const v : m_start) put32(v);
            os.write(head.data(), static_cast<std::streamsize>(head.size()));
            os.write(m_flags.data(), static_cast<std::streamsize>(m_flags.size()));
        }

        // Reads Write's bytes; "" or what is wrong with them. The flag bytes stay in the string read (its head erased), so
        // that a ~0.7 GB table at r226 is not copied.
        std::string Read(std::string content) {
            *this = Table();
            if (content.compare(0, kMagic.size(), kMagic) != 0) return "not a table of unique cores (" + kFileName + ")";
            auto get64 = [&](size_t at) { uint64_t v = 0; for (int i = 0; i < 8; i++) v |= static_cast<uint64_t>(static_cast<uint8_t>(content[at + i])) << (8 * i); return v; };
            auto get32 = [&](size_t at) { uint32_t v = 0; for (int i = 0; i < 4; i++) v |= static_cast<uint32_t>(static_cast<uint8_t>(content[at + i])) << (8 * i); return v; };
            size_t at = kMagic.size();
            if (content.size() < at + 24) return "truncated header";
            uint64_t const taxa = get64(at), rows = get64(at + 8), bytes = get64(at + 16);
            at += 24;
            size_t const left = content.size() - at;
            if (taxa > left / 4 || rows + 1 > (left - 4 * taxa) / 4 || left != 4 * taxa + 4 * (rows + 1) + bytes) {
                return "sizes do not match the file's";
            }
            Table t;
            t.m_first_row.resize(taxa);
            for (uint64_t i = 0; i < taxa; i++) t.m_first_row[i] = get32(at + 4 * i);
            at += 4 * taxa;
            t.m_start.resize(rows + 1);
            for (uint64_t i = 0; i <= rows; i++) t.m_start[i] = get32(at + 4 * i);
            at += 4 * (rows + 1);
            for (uint64_t i = 1; i < taxa; i++) {
                if (t.m_first_row[i] < t.m_first_row[i - 1]) return "rows out of order at taxon " + std::to_string(i);
            }
            if (taxa > 0 && t.m_first_row.back() != rows) return "the taxa's rows do not add up to the table's";
            for (uint64_t i = 1; i <= rows; i++) {
                if (t.m_start[i] < t.m_start[i - 1]) return "row starts out of order at row " + std::to_string(i);
            }
            if (rows > 0 && t.m_start[rows] != bytes) return "the rows' bytes do not add up to the table's";
            content.erase(0, at);
            t.m_flags = std::move(content);
            *this = std::move(t);
            return "";
        }

    private:
        std::vector<uint32_t> m_first_row;  // by taxid; one more than the taxa, the last the rows
        std::vector<uint32_t> m_start;      // each row's first flag byte; one more than the rows
        std::string m_flags;  // four cores a byte; read in place from the file's bytes (Read)
    };

    // "" if `table` has a row per gene id of each taxon of `genomes` (genes 1 to its gene count, Build.h GeneRowsOf) and
    // none besides, else the first difference: a table of another database (a copied file) would count a record's ZS,
    // ZU and ZT from another gene's cores, or, where a gene's cores do not fit its row, from the anchor's seeds (Counter),
    // mixing the two in one sample (F1 audit L6). `Genomes`: GetGenomeMap(), each genome's GetGeneList().
    template<typename Genomes>
    std::string RowsMismatch(Table const& table, Genomes& genomes) {
        size_t max_taxid = 0;
        for (auto const& [taxid, genome] : genomes.GetGenomeMap()) max_taxid = std::max<size_t>(max_taxid, taxid);
        std::vector<uint64_t> expected(max_taxid + 2, 0);
        for (auto const& [taxid, genome] : genomes.GetGenomeMap()) expected[taxid + 1] = genome.GetGeneList().size();
        for (size_t t = 1; t < expected.size(); t++) expected[t] += expected[t - 1];
        auto const& first = table.FirstRows();
        for (size_t taxid = 0; taxid + 1 < std::max(expected.size(), first.size()); taxid++) {
            auto rows = [taxid](auto const& v) -> uint64_t { return taxid + 1 < v.size() ? v[taxid + 1] - v[taxid] : 0; };
            if (rows(first) != rows(expected)) {
                return "taxon " + std::to_string(taxid) + " has " + std::to_string(rows(first)) + " rows, the database " +
                       std::to_string(rows(expected)) + " genes";
            }
        }
        return "";
    }

    // An alignment thread's counter: its own copy of the run's k-mer handler (the build's windows) finds a gene's cores,
    // whose positions it keeps for the last kCachedGenes genes; a copy starts with an empty cache. Inactive without a table.
    template<typename Handler>
    class Counter {
    public:
        Counter() = default;

        Counter(std::shared_ptr<Table const> table, Handler const& handler) : m_table(std::move(table)) {
            if (m_table) {
                m_handler = std::make_unique<Handler>(handler);
                m_mismatched = std::make_shared<std::atomic<size_t>>(0);
            }
        }

        Counter(Counter const& other) : m_table(other.m_table), m_mismatched(other.m_mismatched) {
            if (other.m_handler) m_handler = std::make_unique<Handler>(*other.m_handler);
        }

        Counter& operator=(Counter const& other) {
            if (this != &other) {
                m_table = other.m_table;
                m_mismatched = other.m_mismatched;
                m_handler = other.m_handler ? std::make_unique<Handler>(*other.m_handler) : nullptr;
                m_positions.clear();
            }
            return *this;
        }

        bool Active() const { return m_table != nullptr; }

        // Genes whose cores did not fit their row (the table is not this database's, or the handler not the build's): their
        // records keep the anchor's counts. Shared by the counter's copies (a sample's threads), each gene once per copy.
        size_t Mismatched() const { return m_mismatched ? m_mismatched->load(std::memory_order_relaxed) : 0; }

        // The counts of an alignment of `columns` from gene position `start` on taxid's gene (`gene`: its Sequence() decodes
        // it); nullopt without the table, for a gene beyond it, or one whose cores do not fit its row.
        template<typename Gene>
        std::optional<Counts> Count(uint32_t taxid, uint32_t geneid, Gene const& gene, std::string_view columns, int64_t start) {
            if (!m_table) return std::nullopt;
            auto const row = m_table->Row(taxid, geneid);
            if (!row) return std::nullopt;
            auto* positions = Positions(*row, gene);
            if (!positions) return std::nullopt;
            Table const& table = *m_table;
            size_t const r = *row;
            return CountCovered(std::span<uint32_t const>(*positions), [&](size_t i) { return table.Flags(r, i); }, columns, start);
        }

    private:
        template<typename Gene>
        std::vector<uint32_t> const* Positions(size_t row, Gene const& gene) {
            auto it = m_positions.find(row);
            if (it != m_positions.end()) return it->second.empty() && m_table->Capacity(row) > 0 ? nullptr : &it->second;
            if (m_positions.size() >= kCachedGenes) m_positions.clear();
            auto const sequence = gene.Sequence();
            (*m_handler)(std::string_view(sequence), m_kmers);
            std::vector<uint32_t> positions;
            positions.reserve(m_kmers.size());
            for (auto const& kmer : m_kmers) positions.push_back(static_cast<uint32_t>(kmer.second + kFlank));
            size_t const capacity = m_table->Capacity(row);
            bool const fits = positions.size() <= capacity && positions.size() + 4 > capacity;
            if (!fits) {
                m_mismatched->fetch_add(1, std::memory_order_relaxed);
                positions.clear();
            }
            auto& kept = m_positions[row] = std::move(positions);
            return fits ? &kept : nullptr;
        }

        std::shared_ptr<Table const> m_table;
        std::unique_ptr<Handler> m_handler;
        KmerList m_kmers;
        std::unordered_map<size_t, std::vector<uint32_t>> m_positions;
        std::shared_ptr<std::atomic<size_t>> m_mismatched;  // shared by the copies
    };
}
