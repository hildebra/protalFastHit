// ColumnWeights.h - how conserved each column of a marker gene is, from the family: the table column_weights.tsv
// (ColumnWeightsBuild.h makes it at --build), its per-copy expansion at run time, and what a record's mismatches
// weigh (the "weights" features, Profiler.h).
//
// A mismatch between a read and its reference copy is evidence of different things at different columns. At a column
// the family never changes it is a sequencing error or a read from far away; at a column the species of the genus
// change freely it is what a strain's own mutations look like. The odds that a mismatch at a column with conservation
// p is a real substitution of a close relative, against an error or a distant read's difference, scale with 1 - p, so
// the weight of a mismatch there is w = -log(1 - p): 2.3 nats at 90%, 4.6 at 99%, 6.9 at 99.9%, a ten-fold step in
// the odds per nine. Summed over a read's mismatches the weights are a log-likelihood ratio (docs/claude/2026-10-09-
// site-weighted-evidence-plan). The reads of a species the database lacks sit at the same identity as a missed
// strain's, and half of r226's false positives were reads from other genera: those put their mismatches on conserved
// columns at a share no close relative does (saturation: at 85% identity the fast columns are used up).
//
// Two conservations per column, from the same alignments (the build, ColumnWeightsBuild.h):
//  - within: the species of each genus of the family against their genus's reference copy, averaged over the family's
//    genera with equal weight per genus (no phylum-wide consensus, which would mean little), the pseudocounts over all
//    the species compared (CodeOfPooled): the rate at which a strain's own mutations hit the column;
//  - among: the genus references of the family (at most Settings::genera) against the family's reference copy: the rate
//    at long range, the homoplasy risk of a derived state and the reliability of a difference between twin species.
// And per codon column the conservation of the amino acid (aa), from the same sequences translated in the copies'
// frame (the copies are gene calls, codon 1 at position 0): a non-synonymous mismatch at a column whose amino acid the
// family keeps is the strongest single sign that a read is no close relative, since a strain's mutations are mostly
// synonymous. The family's consensus base per column also polarises the ancestry sites of a species with fewer than
// three congeners (AncestrySites.h Consensus): the species' base is derived where its congener's base is the family's.
//
// The columns are the family reference's positions; every copy of the family maps onto them through its alignment to
// its genus reference and that reference's to the family's (Mapping: a few ungapped runs per copy). The weights are
// stored in half nats, 4 bits each (Code; 0: no estimate), ~0.7 bytes per base of the reference, and expanded per
// copy once per run (Cache, as the ancestry sites are).
#pragma once

#include "BinaryTable.h"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <istream>
#include <iterator>
#include <memory>
#include <mutex>
#include <optional>
#include <ostream>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace protal::column_weights {
    inline const std::string kFileName = "column_weights.tsv";
    inline const std::string kBinaryFileName = "column_weights.bin";  // the binary form a single-file database holds beside it
    inline constexpr double kUnknown = -1;        // the features without the table
    inline constexpr uint8_t kNoCode = 0;         // a column without an estimate
    inline constexpr uint8_t kMaxCode = 15;       // 7.5 nats: 99.94% conserved, all the sample sizes support
    inline constexpr uint8_t kConservedCode = 9;  // 4.5 nats, 98.9% conserved: "a conserved column" in the rate features
    inline constexpr uint8_t kNoBase = 4;         // no consensus base
    inline constexpr size_t kMaxLength = 65535;   // columns and positions are 16-bit
    inline constexpr size_t kMaxCached = 200000;  // copies expanded and kept; the cache starts over beyond

    // A conservation share's code: round(2 w) with w = -log(1 - p), 1..kMaxCode; kNoCode without an estimate. `agree`
    // of `compared` sequences carry the consensus; the pseudocounts (k + 1)/(n + 2) cap what few sequences can claim.
    inline uint8_t CodeOf(size_t agree, size_t compared) {
        if (compared == 0) return kNoCode;
        double const p = (static_cast<double>(agree) + 1.0) / (static_cast<double>(compared) + 2.0);
        double const w = -std::log(1.0 - p);
        long const code = std::lround(2.0 * w);
        return static_cast<uint8_t>(std::clamp<long>(code, 1, kMaxCode));
    }

    // The code of a conservation share (CodeOf's pseudocounts already in it).
    inline uint8_t CodeOfShare(double share, size_t genera) {
        if (genera == 0) return kNoCode;
        double const w = -std::log(1.0 - std::clamp(share, 0.0, 0.999));
        return static_cast<uint8_t>(std::clamp<long>(std::lround(2.0 * w), 1, kMaxCode));
    }

    // The code of a share pooled over groups (the within estimate: the genera's raw shares k/n averaged with equal
    // weight per genus) with the pseudocounts over every sequence compared: (share n + 1)/(n + 2) for `compared` n
    // across the groups. The evidence accumulates across the family's genera, so a column invariant in ten genera of 24
    // species (240 compared) reaches code 11, one invariant in a genus of three code 3. Until 2026-10-10 each genus's
    // share took its own pseudocounts before the mean, which capped every column at 25/26 (code 7, 24 species per genus)
    // and left kConservedCode, and the features that count conserved columns, out of reach.
    inline uint8_t CodeOfPooled(double share, size_t compared) {
        if (compared == 0) return kNoCode;
        double const n = static_cast<double>(compared);
        return CodeOfShare((std::clamp(share, 0.0, 1.0) * n + 1.0) / (n + 2.0), 1);
    }

    // A code's weight in nats.
    inline double Nats(uint8_t code) { return 0.5 * static_cast<double>(code); }

    // A base's code (A 0, C 1, G 2, T 3), 4 for another letter.
    inline uint8_t BaseCode(char c) {
        switch (c) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            default: return 4;
        }
    }

    // The amino acid of a codon (bases as BaseCode, 16 b1 + 4 b2 + b3), as a letter; '*' a stop, 'X' with another letter.
    inline char AminoAcid(uint8_t b1, uint8_t b2, uint8_t b3) {
        static constexpr char kTable[] = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";
        if (b1 > 3 || b2 > 3 || b3 > 3) return 'X';
        return kTable[16 * b1 + 4 * b2 + b3];
    }

    // One ungapped run of a copy's alignment to its family's reference: copy positions [begin, begin + length) sit at
    // columns [column, column + length).
    struct Run {
        uint16_t begin = 0;
        uint16_t column = 0;
        uint16_t length = 0;
        bool operator==(Run const&) const = default;
    };

    // The weights of one family's copies of a gene, by column of the family reference's copy.
    struct Family {
        uint32_t family = 0;        // the family's taxid
        uint32_t gene = 0;
        uint32_t reference = 0;     // the species whose copy gives the columns
        uint16_t genera = 0;        // genera that voted in `among`
        std::vector<uint8_t> within;     // per column, a Code
        std::vector<uint8_t> among;      // per column, a Code
        std::vector<uint8_t> aa;         // per codon column (column / 3), a Code
        std::vector<uint8_t> consensus;  // per column, the family's consensus base (BaseCode), kNoBase without one
        size_t Length() const { return within.size(); }
        bool operator==(Family const&) const = default;
    };

    // A copy's columns expanded along the copy (per position): the codes of its column, kNoCode / kNoBase where the
    // copy maps to no column.
    struct Columns {
        std::vector<uint8_t> within, among, aa, consensus;
        bool Empty() const { return within.empty(); }
        size_t Known() const {
            size_t n = 0;
            for (auto const c : within) n += c != kNoCode;
            return n;
        }
    };

    inline std::string HexOf(std::vector<uint8_t> const& codes) {
        static constexpr char kDigits[] = "0123456789abcdef";
        std::string out(codes.size(), '0');
        for (size_t i = 0; i < codes.size(); i++) out[i] = kDigits[codes[i] & 15];
        return out;
    }

    inline bool HexInto(std::string_view hex, std::vector<uint8_t>& codes) {
        codes.resize(hex.size());
        for (size_t i = 0; i < hex.size(); i++) {
            char const c = hex[i];
            if (c >= '0' && c <= '9') codes[i] = static_cast<uint8_t>(c - '0');
            else if (c >= 'a' && c <= 'f') codes[i] = static_cast<uint8_t>(c - 'a' + 10);
            else return false;
        }
        return true;
    }

    inline std::string BasesOf(std::vector<uint8_t> const& bases) {
        std::string out(bases.size(), 'N');
        for (size_t i = 0; i < bases.size(); i++) out[i] = bases[i] < 4 ? "ACGT"[bases[i]] : 'N';
        return out;
    }

    // The table: per family-gene row its columns' codes, per copy (taxid, gene) its row and runs.
    //
    // The codes are packed (since 2026-10-11): within and among in one byte per column, the consensus base in 2 bits per
    // column with a bit for none, the amino acid's code in 4 bits per codon column: ~1.54 bytes a column, where a
    // Family's four vectors took 3.33 (2.0 GB at GTDB r226's 602M columns). A run reads them from column_weights.bin
    // straight into these arrays (ReadBinary) and decodes a row only for the copies its records touch (Expand, cached by
    // Cache); the text (Read) is parsed into the same arrays, without a vector per line or per copy
    // (docs/claude/2026-10-11-performance-audit).
    class Table {
    public:
        struct Copy {
            uint32_t gene = 0;
            uint32_t row = 0;     // index into the rows
            uint32_t first = 0;   // its runs: [first, first + runs) in m_runs
            uint32_t runs = 0;
        };
        // A family-gene row: whose columns, how many, and where its codes start in the packed arrays.
        struct Row {
            uint32_t family = 0;
            uint32_t gene = 0;
            uint32_t reference = 0;  // the species whose copy gives the columns
            uint16_t genera = 0;     // genera that voted in among
            uint16_t unused = 0;
            uint32_t length = 0;     // columns
            uint32_t unused2 = 0;
            uint64_t column = 0;     // its first column in m_codes, m_bases and m_no_base
            uint64_t codon = 0;      // its first codon column in m_aa
        };
        struct View {
            Row const* row = nullptr;
            size_t row_index = 0;
            Run const* runs = nullptr;
            size_t n_runs = 0;
        };

        bool Empty() const { return m_rows.empty(); }
        size_t Copies() const { return m_copies.size(); }
        size_t Families() const { return m_rows.size(); }
        size_t Species() const {
            size_t n = 0;
            for (size_t t = 0; t + 1 < m_begin.size(); t++) n += m_begin[t + 1] > m_begin[t];
            return n;
        }
        size_t ColumnCount() const { return m_columns; }
        // The bytes the table holds (its arrays), for the log.
        uint64_t Bytes() const {
            return m_rows.size() * sizeof(Row) + m_codes.size() + m_bases.size() + m_no_base.size() + m_aa.size() +
                   m_copies.size() * sizeof(Copy) + m_runs.size() * sizeof(Run) + m_begin.size() * sizeof(uint32_t);
        }

        // A row's codes at its column c (codon column k).
        uint8_t Within(Row const& r, size_t c) const { return m_codes[r.column + c] & 15; }
        uint8_t Among(Row const& r, size_t c) const { return m_codes[r.column + c] >> 4; }
        uint8_t Consensus(Row const& r, size_t c) const {
            size_t const i = r.column + c;
            return (m_no_base[i >> 3] >> (i & 7) & 1) ? kNoBase : static_cast<uint8_t>(m_bases[i >> 2] >> (2 * (i & 3)) & 3);
        }
        uint8_t Aa(Row const& r, size_t k) const {
            size_t const i = r.codon + k;
            return static_cast<uint8_t>(m_aa[i >> 1] >> (4 * (i & 1)) & 15);
        }
        static size_t Codons(size_t length) { return (length + 2) / 3; }

        // Row `index` as a Family (its codes unpacked).
        Family FamilyAt(size_t index) const {
            Row const& r = m_rows[index];
            Family f;
            f.family = r.family;
            f.gene = r.gene;
            f.reference = r.reference;
            f.genera = r.genera;
            f.within.resize(r.length);
            f.among.resize(r.length);
            f.consensus.resize(r.length);
            for (size_t c = 0; c < r.length; c++) {
                f.within[c] = Within(r, c);
                f.among[c] = Among(r, c);
                f.consensus[c] = Consensus(r, c);
            }
            f.aa.resize(Codons(r.length));
            for (size_t k = 0; k < f.aa.size(); k++) f.aa[k] = Aa(r, k);
            return f;
        }

        // Copy (taxid, gene)'s row and runs, or nullopt.
        std::optional<View> Find(uint32_t taxid, uint32_t gene) const {
            if (static_cast<size_t>(taxid) + 1 >= m_begin.size()) return std::nullopt;
            auto const first = m_copies.begin() + m_begin[taxid];
            auto const last = m_copies.begin() + m_begin[taxid + 1];
            auto const it = std::lower_bound(first, last, gene, [](Copy const& c, uint32_t g) { return c.gene < g; });
            if (it == last || it->gene != gene) return std::nullopt;
            return View{ &m_rows[it->row], it->row, m_runs.data() + it->first, it->runs };
        }

        // The copy's columns expanded along a copy of `length` positions.
        std::shared_ptr<Columns const> Expand(uint32_t taxid, uint32_t gene, size_t length) const {
            auto const view = Find(taxid, gene);
            if (!view || length == 0 || length > kMaxLength) return nullptr;
            auto out = std::make_shared<Columns>();
            out->within.assign(length, kNoCode);
            out->among.assign(length, kNoCode);
            out->aa.assign(length, kNoCode);
            out->consensus.assign(length, kNoBase);
            Row const& r = *view->row;
            size_t const codons = Codons(r.length);
            for (size_t k = 0; k < view->n_runs; k++) {
                Run const& run = view->runs[k];
                for (size_t i = 0; i < run.length; i++) {
                    size_t const p = static_cast<size_t>(run.begin) + i, c = static_cast<size_t>(run.column) + i;
                    if (p >= length || c >= r.length) break;
                    uint8_t const code = m_codes[r.column + c];
                    out->within[p] = code & 15;
                    out->among[p] = code >> 4;
                    out->aa[p] = c / 3 < codons ? Aa(r, c / 3) : kNoCode;
                    out->consensus[p] = Consensus(r, c);
                }
            }
            return out;
        }

        // The table of these families and copies (taxid, gene, family row index, runs); a later copy row of the same
        // (taxid, gene) replaces an earlier one.
        static Table FromRows(std::vector<Family> families, std::vector<std::tuple<uint32_t, uint32_t, uint32_t, std::vector<Run>>> copies) {
            std::stable_sort(copies.begin(), copies.end(), [](auto const& a, auto const& b) {
                return std::get<0>(a) != std::get<0>(b) ? std::get<0>(a) < std::get<0>(b) : std::get<1>(a) < std::get<1>(b);
            });
            Table t;
            for (auto const& f : families) {
                size_t const n = f.Length();
                t.AppendRow(f.family, f.gene, f.reference, f.genera, n);
                Row const& r = t.m_rows.back();
                for (size_t c = 0; c < n; c++) {
                    t.SetCodes(r, c, f.within[c], c < f.among.size() ? f.among[c] : kNoCode);
                    t.SetConsensus(r, c, c < f.consensus.size() ? f.consensus[c] : kNoBase);
                }
                for (size_t k = 0; k < Codons(n); k++) t.SetAa(r, k, k < f.aa.size() ? f.aa[k] : kNoCode);
            }
            uint32_t max_taxid = 0;
            for (auto const& c : copies) max_taxid = std::max(max_taxid, std::get<0>(c));
            t.m_begin.assign(copies.empty() ? 1 : static_cast<size_t>(max_taxid) + 2, 0);
            for (size_t i = 0; i < copies.size(); i++) {
                auto const& [taxid, gene, row, runs] = copies[i];
                if (row >= t.m_rows.size()) continue;
                if (!t.m_copies.empty() && i > 0 && std::get<0>(copies[i - 1]) == taxid && std::get<1>(copies[i - 1]) == gene) {
                    t.m_copies.back() = Copy{ gene, row, static_cast<uint32_t>(t.m_runs.size()), static_cast<uint32_t>(runs.size()) };
                } else {
                    t.m_copies.push_back(Copy{ gene, row, static_cast<uint32_t>(t.m_runs.size()), static_cast<uint32_t>(runs.size()) });
                }
                t.m_runs.insert(t.m_runs.end(), runs.begin(), runs.end());
                t.m_begin[taxid + 1] = static_cast<uint32_t>(t.m_copies.size());
            }
            for (size_t i = 1; i < t.m_begin.size(); i++) t.m_begin[i] = std::max(t.m_begin[i], t.m_begin[i - 1]);
            return t;
        }

        // The text: a header, then per family and gene a line
        //   F <family> <gene> <reference> <genera> <within hex> <among hex> <aa hex> <consensus bases>
        // and per copy a line
        //   C <taxid> <gene> <family> <runs: begin:column:length,...>
        // (tab-separated; the families before the copies that use them).
        void Write(std::ostream& os) const {
            static constexpr char kDigits[] = "0123456789abcdef";
            os << "#column_weights\tv1\tF: family gene reference genera within among aa consensus; C: taxid gene family runs\n";
            std::string within, among, aa, consensus;
            for (Row const& r : m_rows) {
                within.resize(r.length);
                among.resize(r.length);
                consensus.resize(r.length);
                for (size_t c = 0; c < r.length; c++) {
                    within[c] = kDigits[Within(r, c)];
                    among[c] = kDigits[Among(r, c)];
                    uint8_t const base = Consensus(r, c);
                    consensus[c] = base < 4 ? "ACGT"[base] : 'N';
                }
                aa.resize(Codons(r.length));
                for (size_t k = 0; k < aa.size(); k++) aa[k] = kDigits[Aa(r, k)];
                os << "F\t" << r.family << '\t' << r.gene << '\t' << r.reference << '\t' << r.genera << '\t' << within << '\t'
                   << among << '\t' << aa << '\t' << consensus << '\n';
            }
            for (size_t t = 0; t + 1 < m_begin.size(); t++) {
                for (uint32_t i = m_begin[t]; i < m_begin[t + 1]; i++) {
                    Copy const& c = m_copies[i];
                    os << "C\t" << t << '\t' << c.gene << '\t' << m_rows[c.row].family << '\t';
                    for (uint32_t r = 0; r < c.runs; r++) {
                        Run const& run = m_runs[c.first + r];
                        os << (r ? "," : "") << run.begin << ':' << run.column << ':' << run.length;
                    }
                    os << '\n';
                }
            }
        }

        // Reads what Write wrote; "" or what is wrong.
        std::string Read(std::istream& is) {
            std::string const text{ std::istreambuf_iterator<char>(is), std::istreambuf_iterator<char>() };
            return Read(std::string_view(text));
        }

        // Reads what Write wrote, line by line into the packed arrays; "" or what is wrong (and the table unchanged). The
        // copies come sorted from Write and go straight into place; others (a hand-made table) are sorted after, a later
        // line of the same copy replacing an earlier one.
        std::string Read(std::string_view text) {
            Table t;
            std::unordered_map<uint64_t, uint32_t> row_of;  // (family << 32 | gene) -> row
            std::vector<uint32_t> taxa;                     // each copy's taxid, in the order read
            bool in_order = true;
            size_t n = 0;
            auto number = [](std::string_view s, auto& out) {
                auto const r = std::from_chars(s.data(), s.data() + s.size(), out);
                return r.ec == std::errc() && r.ptr == s.data() + s.size();
            };
            auto hex = [](char ch) -> int {
                if (ch >= '0' && ch <= '9') return ch - '0';
                if (ch >= 'a' && ch <= 'f') return ch - 'a' + 10;
                return -1;
            };
            std::array<std::string_view, 10> f;
            for (size_t pos = 0; pos < text.size();) {
                size_t const eol = text.find('\n', pos);
                std::string_view const line = text.substr(pos, eol == std::string_view::npos ? std::string_view::npos : eol - pos);
                pos = eol == std::string_view::npos ? text.size() : eol + 1;
                n++;
                if (line.empty() || line[0] == '#') continue;
                size_t fields = 0, start = 0;
                while (true) {
                    size_t const tab = line.find('\t', start);
                    if (fields < f.size()) f[fields] = line.substr(start, tab == std::string_view::npos ? std::string_view::npos : tab - start);
                    fields++;
                    if (tab == std::string_view::npos) break;
                    start = tab + 1;
                }
                if (f[0] == "F") {
                    if (fields != 9) return "line " + std::to_string(n) + ": an F line has 9 fields";
                    uint32_t family = 0, gene = 0, reference = 0;
                    uint16_t genera = 0;
                    if (!number(f[1], family) || !number(f[2], gene) || !number(f[3], reference) || !number(f[4], genera)) {
                        return "line " + std::to_string(n) + ": family, gene, reference and genera are numbers";
                    }
                    size_t const length = f[5].size();
                    if (f[6].size() != length || f[7].size() != Codons(length) || f[8].size() != length) {
                        return "line " + std::to_string(n) + ": the columns' lengths differ";
                    }
                    if (length > kMaxLength) return "line " + std::to_string(n) + ": more than 65535 columns";
                    t.AppendRow(family, gene, reference, genera, length);
                    Row const& r = t.m_rows.back();
                    for (size_t c = 0; c < length; c++) {
                        int const w = hex(f[5][c]), a = hex(f[6][c]);
                        if (w < 0 || a < 0) return "line " + std::to_string(n) + ": the weights are hex digits";
                        t.SetCodes(r, c, static_cast<uint8_t>(w), static_cast<uint8_t>(a));
                        t.SetConsensus(r, c, BaseCode(f[8][c]));
                    }
                    for (size_t k = 0; k < f[7].size(); k++) {
                        int const code = hex(f[7][k]);
                        if (code < 0) return "line " + std::to_string(n) + ": the weights are hex digits";
                        t.SetAa(r, k, static_cast<uint8_t>(code));
                    }
                    row_of[static_cast<uint64_t>(family) << 32 | gene] = static_cast<uint32_t>(t.m_rows.size() - 1);
                } else if (f[0] == "C") {
                    if (fields != 5) return "line " + std::to_string(n) + ": a C line has 5 fields";
                    uint32_t taxid = 0, gene = 0, family = 0;
                    if (!number(f[1], taxid) || !number(f[2], gene) || !number(f[3], family)) {
                        return "line " + std::to_string(n) + ": taxid, gene and family are numbers";
                    }
                    auto const it = row_of.find(static_cast<uint64_t>(family) << 32 | gene);
                    if (it == row_of.end()) return "line " + std::to_string(n) + ": the family's F line comes first";
                    uint32_t const first = static_cast<uint32_t>(t.m_runs.size());
                    std::string_view rest = f[4];
                    while (!rest.empty()) {
                        size_t const comma = rest.find(',');
                        std::string_view const item = rest.substr(0, comma);
                        rest = comma == std::string_view::npos ? std::string_view() : rest.substr(comma + 1);
                        size_t const a = item.find(':'), b = item.rfind(':');
                        Run run;
                        if (a == std::string_view::npos || b == a || !number(item.substr(0, a), run.begin) ||
                            !number(item.substr(a + 1, b - a - 1), run.column) || !number(item.substr(b + 1), run.length)) {
                            return "line " + std::to_string(n) + ": a run is begin:column:length";
                        }
                        t.m_runs.push_back(run);
                    }
                    if (!taxa.empty() && (taxa.back() > taxid || (taxa.back() == taxid && t.m_copies.back().gene >= gene))) in_order = false;
                    t.m_copies.push_back(Copy{ gene, it->second, first, static_cast<uint32_t>(t.m_runs.size()) - first });
                    taxa.push_back(taxid);
                } else {
                    return "line " + std::to_string(n) + ": a line starts with F or C";
                }
            }
            if (!in_order) t.SortCopies(taxa);
            t.BeginFrom(taxa);
            *this = std::move(t);
            return {};
        }

        inline static constexpr char kMagic[8] = {'P', 'R', 'C', 'O', 'L', 'W', 'T', '1'};
        static constexpr uint64_t kBinaryVersion = 1;

        // column_weights.bin (BinaryTable.h): this table's arrays, made from a column_weights.tsv of `text_size` bytes.
        // False on a write error.
        bool WriteBinary(std::ostream& os, uint64_t text_size) const {
            using binary_table::Of;
            return binary_table::Write(os, kMagic, kBinaryVersion, text_size, { Of(m_rows), Of(m_codes), Of(m_bases), Of(m_no_base),
                                                                                Of(m_aa), Of(m_copies), Of(m_runs), Of(m_begin) });
        }

        // Reads column_weights.bin `file`, whose header is `header` (binary_table::ReadHeader), into this table with up to
        // `threads` threads; "" or what is wrong (and the table unchanged).
        std::string ReadBinary(db::DbFile const& file, binary_table::Header const& header, int threads) {
            if (header.version != kBinaryVersion) {
                return "version " + std::to_string(header.version) + " (this protal reads version " + std::to_string(kBinaryVersion) + ")";
            }
            if (header.sizes.size() != 8) return "it holds " + std::to_string(header.sizes.size()) + " arrays, not 8";
            size_t const units[8] = { sizeof(Row), 1, 1, 1, 1, sizeof(Copy), sizeof(Run), sizeof(uint32_t) };
            for (size_t i = 0; i < 8; i++) {
                if (header.sizes[i] % units[i] != 0) return "array " + std::to_string(i) + " is not a whole number of records";
            }
            Table t;
            t.m_rows.resize(header.sizes[0] / sizeof(Row));
            t.m_codes.resize(header.sizes[1]);
            t.m_bases.resize(header.sizes[2]);
            t.m_no_base.resize(header.sizes[3]);
            t.m_aa.resize(header.sizes[4]);
            t.m_copies.resize(header.sizes[5] / sizeof(Copy));
            t.m_runs.resize(header.sizes[6] / sizeof(Run));
            t.m_begin.resize(header.sizes[7] / sizeof(uint32_t));
            using binary_table::Into;
            auto const error = binary_table::ReadArrays(file, header, { Into(t.m_rows), Into(t.m_codes), Into(t.m_bases), Into(t.m_no_base),
                                                                         Into(t.m_aa), Into(t.m_copies), Into(t.m_runs), Into(t.m_begin) },
                                                        threads);
            if (!error.empty()) return error;
            t.m_columns = t.m_codes.size();
            for (Row const& r : t.m_rows) t.m_codons = std::max<size_t>(t.m_codons, r.codon + Codons(r.length));
            if (auto const invalid = t.Invalid(); !invalid.empty()) return invalid;
            *this = std::move(t);
            return "";
        }

    private:
        // Appends a row of `length` columns without estimates or consensus bases, which SetCodes, SetConsensus and SetAa set.
        void AppendRow(uint32_t family, uint32_t gene, uint32_t reference, uint16_t genera, size_t length) {
            Row r;
            r.family = family;
            r.gene = gene;
            r.reference = reference;
            r.genera = genera;
            r.length = static_cast<uint32_t>(length);
            r.column = m_columns;
            r.codon = m_codons;
            m_columns += length;
            m_codons += Codons(length);
            m_codes.resize(m_columns, 0);
            m_bases.resize((m_columns + 3) / 4, 0);
            m_no_base.resize((m_columns + 7) / 8, 0);
            m_aa.resize((m_codons + 1) / 2, 0);
            m_rows.push_back(r);
        }

        // Each set once per column of a row AppendRow added.
        void SetCodes(Row const& r, size_t c, uint8_t within, uint8_t among) {
            m_codes[r.column + c] = static_cast<uint8_t>((within & 15) | (among & 15) << 4);
        }
        void SetConsensus(Row const& r, size_t c, uint8_t base) {
            size_t const i = r.column + c;
            if (base > 3) m_no_base[i >> 3] |= static_cast<uint8_t>(1 << (i & 7));
            else m_bases[i >> 2] |= static_cast<uint8_t>(base << (2 * (i & 3)));
        }
        void SetAa(Row const& r, size_t k, uint8_t code) {
            size_t const i = r.codon + k;
            m_aa[i >> 1] |= static_cast<uint8_t>((code & 15) << (4 * (i & 1)));
        }

        // The copies (each of taxa[i]) sorted by taxid and gene, of copies read twice the last; taxa sorted alike.
        void SortCopies(std::vector<uint32_t>& taxa) {
            std::vector<uint32_t> order(m_copies.size());
            for (size_t i = 0; i < order.size(); i++) order[i] = static_cast<uint32_t>(i);
            std::stable_sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) {
                return taxa[a] != taxa[b] ? taxa[a] < taxa[b] : m_copies[a].gene < m_copies[b].gene;
            });
            std::vector<Copy> copies;
            std::vector<Run> runs;
            std::vector<uint32_t> sorted_taxa;
            for (size_t k = 0; k < order.size(); k++) {
                Copy const& c = m_copies[order[k]];
                uint32_t const taxid = taxa[order[k]];
                if (k + 1 < order.size() && taxa[order[k + 1]] == taxid && m_copies[order[k + 1]].gene == c.gene) continue;  // a later one
                copies.push_back(Copy{ c.gene, c.row, static_cast<uint32_t>(runs.size()), c.runs });
                runs.insert(runs.end(), m_runs.begin() + c.first, m_runs.begin() + c.first + c.runs);
                sorted_taxa.push_back(taxid);
            }
            m_copies = std::move(copies);
            m_runs = std::move(runs);
            taxa = std::move(sorted_taxa);
        }

        // m_begin from the (sorted) copies' taxids.
        void BeginFrom(std::vector<uint32_t> const& taxa) {
            uint32_t const max_taxid = taxa.empty() ? 0 : taxa.back();
            m_begin.assign(taxa.empty() ? 1 : static_cast<size_t>(max_taxid) + 2, 0);
            for (uint32_t const taxid : taxa) m_begin[taxid + 1]++;
            for (size_t i = 1; i < m_begin.size(); i++) m_begin[i] += m_begin[i - 1];
        }

        // What makes the arrays unusable (a damaged binary table), or "".
        std::string Invalid() const {
            for (Row const& r : m_rows) {
                if (r.length > kMaxLength || r.column + r.length > m_codes.size() || (r.column + r.length + 3) / 4 > m_bases.size() ||
                    (r.column + r.length + 7) / 8 > m_no_base.size() || (r.codon + Codons(r.length) + 1) / 2 > m_aa.size()) {
                    return "a row's columns lie outside the codes";
                }
            }
            if (m_begin.empty() || m_begin.front() != 0 || m_begin.back() != m_copies.size()) return "the copies' index does not cover them";
            for (size_t t = 0; t + 1 < m_begin.size(); t++) {
                if (m_begin[t + 1] < m_begin[t]) return "the copies' index is not sorted";
                for (uint32_t i = m_begin[t]; i < m_begin[t + 1]; i++) {
                    Copy const& c = m_copies[i];
                    if (c.row >= m_rows.size() || uint64_t{ c.first } + c.runs > m_runs.size()) return "a copy's row or runs lie outside the table";
                    if (i > m_begin[t] && m_copies[i - 1].gene >= c.gene) return "a species' copies are not sorted by gene";
                }
            }
            return "";
        }

        std::vector<Row> m_rows;
        std::vector<uint8_t> m_codes;    // per column: within | among << 4
        std::vector<uint8_t> m_bases;    // per column 2 bits, four to a byte: the consensus base (BaseCode 0-3)
        std::vector<uint8_t> m_no_base;  // per column 1 bit, eight to a byte: no consensus base (kNoBase)
        std::vector<uint8_t> m_aa;       // per codon column 4 bits, two to a byte: the amino acid's code
        size_t m_columns = 0, m_codons = 0;
        std::vector<Copy> m_copies;      // by taxid, then gene
        std::vector<uint32_t> m_begin;   // taxid -> first copy; size max taxid + 2
        std::vector<Run> m_runs;
    };
    static_assert(sizeof(Table::Row) == 40 && sizeof(Table::Copy) == 16 && sizeof(Run) == 6, "column_weights.bin record sizes");
    static_assert(std::is_trivially_copyable_v<Table::Row> && std::is_trivially_copyable_v<Table::Copy> &&
                  std::is_trivially_copyable_v<Run>);

    // The expanded columns of every copy a run touched, computed once (thread-safe; nullptr for a copy without a row).
    class Cache {
    public:
        template<typename Genomes>
        std::shared_ptr<Columns const> Get(uint32_t taxid, uint32_t geneid, Genomes& genomes, Table const& table) {
            if (table.Empty()) return nullptr;
            uint64_t const key = (static_cast<uint64_t>(taxid) << 32) | geneid;
            {
                std::lock_guard<std::mutex> lock(m_mutex);
                if (auto const it = m_columns.find(key); it != m_columns.end()) return it->second;
            }
            std::shared_ptr<Columns const> result;
            if (genomes.HasGene(taxid, geneid)) {
                auto const own = genomes.GetGeneOMP(taxid, geneid).Sequence();
                result = table.Expand(taxid, geneid, own.View().size());
            }
            std::lock_guard<std::mutex> lock(m_mutex);
            if (m_columns.size() >= kMaxCached) m_columns.clear();
            return m_columns.try_emplace(key, result).first->second;
        }

        void Clear() {
            std::lock_guard<std::mutex> lock(m_mutex);
            m_columns.clear();
        }

    private:
        mutable std::mutex m_mutex;
        std::unordered_map<uint64_t, std::shared_ptr<Columns const>> m_columns;
    };

    // What a record's aligned bases weigh: of its aligned columns with an estimate, their number and the sum of their
    // within codes; of its mismatches (X operations) at such columns, their number, the sums of their within and among
    // codes and how many lie at conserved columns (within >= kConservedCode); and the mismatches classified by codon
    // (the read's codon is the copy's with the read's base at the mismatch: synonymous or not), of the non-synonymous
    // ones those at codon columns whose amino acid is conserved (aa >= kConservedCode). `copy` is the reference copy
    // (for the codons); a record without SEQ ("*") counts nothing.
    struct Counts {
        uint64_t columns = 0;  // aligned bases, with or without an estimate
        uint64_t aligned = 0, within_aligned = 0;
        uint64_t mismatches = 0, within_mismatches = 0, among_mismatches = 0, conserved_mismatches = 0;
        uint64_t synonymous = 0, nonsynonymous = 0, nonsynonymous_conserved = 0;
    };

    inline Counts Count(Columns const& columns, std::string const& cigar, size_t pos, std::string const& seq, std::string_view copy) {
        Counts c;
        if (columns.Empty() || seq.empty() || seq == "*") return c;
        size_t const start = pos == 0 ? 0 : pos - 1;
        size_t ref = start, query = 0, run = 0;
        size_t const length = columns.within.size();
        for (char const ch : cigar) {
            if (ch >= '0' && ch <= '9') {
                run = run * 10 + static_cast<size_t>(ch - '0');
                continue;
            }
            if (ch == 'M' || ch == '=' || ch == 'X') {
                for (size_t i = 0; i < run && ref + i < length; i++) {
                    size_t const p = ref + i;
                    c.columns++;
                    uint8_t const w = columns.within[p];
                    if (w == kNoCode) continue;
                    c.aligned++;
                    c.within_aligned += w;
                    if (ch != 'X') continue;
                    c.mismatches++;
                    c.within_mismatches += w;
                    c.among_mismatches += columns.among[p];
                    c.conserved_mismatches += w >= kConservedCode;
                    // The codon: the copy's three bases, the read's at p.
                    size_t const codon = p - p % 3;
                    if (codon + 3 > copy.size() || query + i >= seq.size()) continue;
                    std::array<uint8_t, 3> bases{ BaseCode(copy[codon]), BaseCode(copy[codon + 1]), BaseCode(copy[codon + 2]) };
                    char const before = AminoAcid(bases[0], bases[1], bases[2]);
                    bases[p - codon] = BaseCode(seq[query + i]);
                    char const after = AminoAcid(bases[0], bases[1], bases[2]);
                    if (before == 'X' || after == 'X') continue;
                    if (before == after) {
                        c.synonymous++;
                    } else {
                        c.nonsynonymous++;
                        c.nonsynonymous_conserved += columns.aa[p] >= kConservedCode;
                    }
                }
                ref += run;
                query += run;
            } else if (ch == 'I' || ch == 'S') {
                query += run;
            } else if (ch == 'D' || ch == 'N') {
                ref += run;
            }
            run = 0;
        }
        return c;
    }
}
