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

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <istream>
#include <memory>
#include <mutex>
#include <optional>
#include <ostream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace protal::column_weights {
    inline const std::string kFileName = "column_weights.tsv";
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

    // The table: per (family, gene) its Family, per copy (taxid, gene) its family and runs.
    class Table {
    public:
        struct Copy {
            uint32_t gene = 0;
            uint32_t row = 0;     // index into the families
            uint32_t first = 0;   // its runs: [first, first + runs) in m_runs
            uint32_t runs = 0;
        };
        struct View {
            Family const* family = nullptr;
            Run const* runs = nullptr;
            size_t n_runs = 0;
        };

        bool Empty() const { return m_families.empty(); }
        size_t Copies() const { return m_copies.size(); }
        size_t Families() const { return m_families.size(); }
        size_t Species() const {
            size_t n = 0;
            for (size_t t = 0; t + 1 < m_begin.size(); t++) n += m_begin[t + 1] > m_begin[t];
            return n;
        }
        size_t ColumnCount() const {
            size_t n = 0;
            for (auto const& f : m_families) n += f.Length();
            return n;
        }
        std::vector<Family> const& AllFamilies() const { return m_families; }

        // Copy (taxid, gene)'s family and runs, or nullopt.
        std::optional<View> Find(uint32_t taxid, uint32_t gene) const {
            if (static_cast<size_t>(taxid) + 1 >= m_begin.size()) return std::nullopt;
            auto const first = m_copies.begin() + m_begin[taxid];
            auto const last = m_copies.begin() + m_begin[taxid + 1];
            auto const it = std::lower_bound(first, last, gene, [](Copy const& c, uint32_t g) { return c.gene < g; });
            if (it == last || it->gene != gene) return std::nullopt;
            return View{ &m_families[it->row], m_runs.data() + it->first, it->runs };
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
            Family const& f = *view->family;
            for (size_t r = 0; r < view->n_runs; r++) {
                Run const& run = view->runs[r];
                for (size_t i = 0; i < run.length; i++) {
                    size_t const p = static_cast<size_t>(run.begin) + i, c = static_cast<size_t>(run.column) + i;
                    if (p >= length || c >= f.Length()) break;
                    out->within[p] = f.within[c];
                    out->among[p] = f.among[c];
                    out->aa[p] = c / 3 < f.aa.size() ? f.aa[c / 3] : kNoCode;
                    out->consensus[p] = f.consensus[c];
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
            t.m_families = std::move(families);
            uint32_t max_taxid = 0;
            for (auto const& c : copies) max_taxid = std::max(max_taxid, std::get<0>(c));
            t.m_begin.assign(copies.empty() ? 1 : static_cast<size_t>(max_taxid) + 2, 0);
            for (size_t i = 0; i < copies.size(); i++) {
                auto const& [taxid, gene, row, runs] = copies[i];
                if (row >= t.m_families.size()) continue;
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
            os << "#column_weights\tv1\tF: family gene reference genera within among aa consensus; C: taxid gene family runs\n";
            for (auto const& f : m_families) {
                os << "F\t" << f.family << '\t' << f.gene << '\t' << f.reference << '\t' << f.genera << '\t' << HexOf(f.within) << '\t'
                   << HexOf(f.among) << '\t' << HexOf(f.aa) << '\t' << BasesOf(f.consensus) << '\n';
            }
            for (size_t t = 0; t + 1 < m_begin.size(); t++) {
                for (uint32_t i = m_begin[t]; i < m_begin[t + 1]; i++) {
                    Copy const& c = m_copies[i];
                    os << "C\t" << t << '\t' << c.gene << '\t' << m_families[c.row].family << '\t';
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
            std::vector<Family> families;
            std::unordered_map<uint64_t, uint32_t> row_of;  // (family << 32 | gene) -> row
            std::vector<std::tuple<uint32_t, uint32_t, uint32_t, std::vector<Run>>> copies;
            std::string line;
            size_t n = 0;
            auto number = [](std::string_view s, auto& out) {
                auto const r = std::from_chars(s.data(), s.data() + s.size(), out);
                return r.ec == std::errc() && r.ptr == s.data() + s.size();
            };
            while (std::getline(is, line)) {
                n++;
                if (line.empty() || line[0] == '#') continue;
                std::vector<std::string_view> f;
                std::string_view const v = line;
                size_t start = 0;
                for (size_t tab; (tab = v.find('\t', start)) != std::string_view::npos; start = tab + 1) f.push_back(v.substr(start, tab - start));
                f.push_back(v.substr(start));
                if (f[0] == "F") {
                    if (f.size() != 9) return "line " + std::to_string(n) + ": an F line has 9 fields";
                    Family fam;
                    if (!number(f[1], fam.family) || !number(f[2], fam.gene) || !number(f[3], fam.reference) || !number(f[4], fam.genera)) {
                        return "line " + std::to_string(n) + ": family, gene, reference and genera are numbers";
                    }
                    if (!HexInto(f[5], fam.within) || !HexInto(f[6], fam.among) || !HexInto(f[7], fam.aa)) {
                        return "line " + std::to_string(n) + ": the weights are hex digits";
                    }
                    if (fam.among.size() != fam.within.size() || fam.aa.size() != (fam.within.size() + 2) / 3 || f[8].size() != fam.within.size()) {
                        return "line " + std::to_string(n) + ": the columns' lengths differ";
                    }
                    if (fam.within.size() > kMaxLength) return "line " + std::to_string(n) + ": more than 65535 columns";
                    fam.consensus.resize(f[8].size());
                    for (size_t i = 0; i < f[8].size(); i++) fam.consensus[i] = BaseCode(f[8][i]);
                    row_of[static_cast<uint64_t>(fam.family) << 32 | fam.gene] = static_cast<uint32_t>(families.size());
                    families.push_back(std::move(fam));
                } else if (f[0] == "C") {
                    if (f.size() != 5) return "line " + std::to_string(n) + ": a C line has 5 fields";
                    uint32_t taxid = 0, gene = 0, family = 0;
                    if (!number(f[1], taxid) || !number(f[2], gene) || !number(f[3], family)) {
                        return "line " + std::to_string(n) + ": taxid, gene and family are numbers";
                    }
                    auto const it = row_of.find(static_cast<uint64_t>(family) << 32 | gene);
                    if (it == row_of.end()) return "line " + std::to_string(n) + ": the family's F line comes first";
                    std::vector<Run> runs;
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
                        runs.push_back(run);
                    }
                    copies.emplace_back(taxid, gene, it->second, std::move(runs));
                } else {
                    return "line " + std::to_string(n) + ": a line starts with F or C";
                }
            }
            *this = FromRows(std::move(families), std::move(copies));
            return {};
        }

    private:
        std::vector<Family> m_families;
        std::vector<Copy> m_copies;    // by taxid, then gene
        std::vector<uint32_t> m_begin; // taxid -> first copy; size max taxid + 2
        std::vector<Run> m_runs;
    };

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
