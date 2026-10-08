// StrainAlleles.h - strain_alleles.tsv: per species' copy of a marker gene, up to a few alleles of the species' other
// genomes, each as the part of the representative's copy it covers and its edits there (substitutions, insertions,
// deletions). StrainAllelesBuild.h makes it at --build from the full reference; a run loads it into the GenomeLoader.
//
// A run uses it twice, both times as "which of the species' known strains does this read look like":
//  - the alignment handler (AlignmentStrategy.h) gives each candidate alignment its best allele's shift: the change in
//    the read's mismatches if that allele were the reference (AlignmentInfo::allele_shift), so that a read of a known
//    strain scores on its species as it would on the strain's own copy (--no_allele_scores: not);
//  - the profiler counts, per taxon, the share of its reads' differences from the reference that an allele explains and
//    the identity it gains (Profiler.h, the "alleles" features): a strain's few differences mostly sit where the
//    species' strains differ, a novel congener's mostly do not.
// The records keep their alignment to the representative (CIGAR, identity); only the score moves.
//
// Which genomes may give alleles: --allele_genome_share of each species' genomes, chosen by a hash of the accession
// (AlleleGenome). build_gtdb_database.py simulates strains only from the other genomes, so that a simulated strain is
// never one of its species' alleles (else its reads would match an allele exactly, a gain no real sample sees), and
// the ancestry report takes its alleles from the same genomes (docs/claude/2026-10-08-r226-v18, section 5).
#pragma once

#include <algorithm>
#include <charconv>
#include <cstdint>
#include <cstdlib>
#include <istream>
#include <ostream>
#include <span>
#include <sstream>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace protal::strain_alleles {
    inline const std::string kFileName = "strain_alleles.tsv";
    inline constexpr double kUnknown = -1;          // features without the table, or without a record on a copy with alleles
    inline constexpr uint32_t kIndelTolerance = 4;  // a read's indel within this many bases of an allele's, of its kind and
                                                    // length, is the allele's (aligners place an indel in a repeat apart)

    // FNV-1a, 64 bits: the same in scripts/mini_db/gtdb_to_protal_db.py (allele_genome), so that protal, the build
    // script and the ancestry report split the genomes alike.
    inline uint64_t Fnv1a(std::string_view text) {
        uint64_t h = 1469598103934665603ull;
        for (unsigned char const c : text) {
            h ^= c;
            h *= 1099511628211ull;
        }
        return h;
    }

    // Whether the genome `accession` (the full reference's second header word, as the converter writes it: GCA_.../
    // GCF_...) may give alleles: its hash's top 53 bits as a fraction below `share`.
    inline bool AlleleGenome(std::string_view accession, double share) {
        if (share >= 1) return true;
        if (share <= 0) return false;
        return static_cast<double>(Fnv1a(accession) >> 11) * 0x1.0p-53 < share;
    }

    inline uint8_t BaseCode(char c) {
        switch (c) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            default: return 4;
        }
    }
    inline constexpr char kBases[] = "ACGT";

    enum Kind : uint8_t { kSubstitution = 0, kInsertion = 1, kDeletion = 2 };

    // An edit of the representative's copy: the base at `pos` substituted, `length` bases inserted before `pos`, or
    // `length` bases deleted from `pos` (4 bytes: kind 2 bits, base 2 bits, length 12 bits).
    struct Edit {
        uint16_t pos = 0;
        uint16_t packed = 0;

        static constexpr uint16_t kMaxLength = 4095;
        static Edit Substitution(uint16_t pos, uint8_t base) {
            return Edit{ pos, static_cast<uint16_t>(kSubstitution | (base & 3) << 2 | 1 << 4) };
        }
        static Edit Indel(Kind kind, uint16_t pos, uint16_t length) {
            return Edit{ pos, static_cast<uint16_t>(kind | std::min<uint16_t>(length, kMaxLength) << 4) };
        }
        Kind GetKind() const { return static_cast<Kind>(packed & 3); }
        uint8_t Base() const { return static_cast<uint8_t>(packed >> 2 & 3); }
        uint16_t Length() const { return static_cast<uint16_t>(packed >> 4); }
        bool operator==(Edit const&) const = default;
        bool operator<(Edit const& o) const { return pos != o.pos ? pos < o.pos : packed < o.packed; }
    };

    // An allele: the representative's bases [begin, end) it covers, and its edits there by position.
    struct Allele {
        uint16_t begin = 0;
        uint16_t end = 0;
        std::vector<Edit> edits;
        bool operator==(Allele const&) const = default;
    };

    // The edits as written in strain_alleles.tsv: "pos" + base (substitution), "pos" + "i" + length (insertion),
    // "pos" + "d" + length (deletion), comma-separated; an allele as "begin-end:edits".
    inline void WriteAllele(std::ostream& os, Allele const& allele) {
        os << allele.begin << '-' << allele.end << ':';
        for (size_t i = 0; i < allele.edits.size(); i++) {
            auto const& e = allele.edits[i];
            if (i) os << ',';
            os << e.pos;
            switch (e.GetKind()) {
                case kSubstitution: os << kBases[e.Base()]; break;
                case kInsertion: os << 'i' << e.Length(); break;
                default: os << 'd' << e.Length(); break;
            }
        }
    }

    inline std::string AlleleText(Allele const& allele) {
        std::ostringstream os;
        WriteAllele(os, allele);
        return os.str();
    }

    // One allele of WriteAllele's text; false if it is not one.
    inline bool ParseAllele(std::string_view text, Allele& allele) {
        allele = Allele{};
        auto number = [](std::string_view& s, uint32_t& value) {
            auto const [p, ec] = std::from_chars(s.data(), s.data() + s.size(), value);
            if (ec != std::errc()) return false;
            s.remove_prefix(static_cast<size_t>(p - s.data()));
            return true;
        };
        uint32_t begin = 0, end = 0;
        if (!number(text, begin) || text.empty() || text[0] != '-') return false;
        text.remove_prefix(1);
        if (!number(text, end) || text.empty() || text[0] != ':' || begin >= end || end > UINT16_MAX) return false;
        text.remove_prefix(1);
        allele.begin = static_cast<uint16_t>(begin);
        allele.end = static_cast<uint16_t>(end);
        while (!text.empty()) {
            uint32_t pos = 0;
            if (!number(text, pos) || pos > UINT16_MAX || text.empty()) return false;
            char const op = text[0];
            text.remove_prefix(1);
            if (op == 'i' || op == 'd') {
                uint32_t length = 0;
                if (!number(text, length) || length == 0 || length > Edit::kMaxLength) return false;
                allele.edits.push_back(Edit::Indel(op == 'i' ? kInsertion : kDeletion, static_cast<uint16_t>(pos),
                                                   static_cast<uint16_t>(length)));
            } else {
                uint8_t const base = BaseCode(op);
                if (base > 3) return false;
                allele.edits.push_back(Edit::Substitution(static_cast<uint16_t>(pos), base));
            }
            if (!text.empty()) {
                if (text[0] != ',' || text.size() == 1) return false;
                text.remove_prefix(1);
            }
        }
        return std::is_sorted(allele.edits.begin(), allele.edits.end());
    }

    // The table of every copy's alleles, by taxid and gene.
    class Table {
    public:
        struct Stored {
            uint32_t first_edit = 0;
            uint16_t edits = 0;
            uint16_t begin = 0;
            uint16_t end = 0;
        };
        // One allele as the table holds it: its range and a view of its edits.
        struct View {
            uint16_t begin = 0;
            uint16_t end = 0;
            std::span<Edit const> edits;
        };

        bool Empty() const { return m_entries.empty(); }
        size_t Copies() const { return m_entries.size(); }
        size_t Alleles() const { return m_alleles.size(); }
        size_t Edits() const { return m_edits.size(); }
        size_t Species() const {
            size_t n = 0;
            for (size_t t = 0; t + 1 < m_begin.size(); t++) n += m_begin[t + 1] > m_begin[t];
            return n;
        }

        // The alleles of copy (taxid, gene): an empty span without any.
        std::span<Stored const> Of(uint32_t taxid, uint32_t gene) const {
            if (static_cast<size_t>(taxid) + 1 >= m_begin.size() || gene > UINT16_MAX) return {};
            auto const first = m_entries.begin() + m_begin[taxid];
            auto const last = m_entries.begin() + m_begin[taxid + 1];
            auto const it = std::lower_bound(first, last, gene, [](Entry const& e, uint32_t g) { return e.gene < g; });
            if (it == last || it->gene != gene) return {};
            return std::span<Stored const>(m_alleles.data() + it->first_allele, it->alleles);
        }

        View Get(Stored const& s) const {
            return View{ s.begin, s.end, std::span<Edit const>(m_edits.data() + s.first_edit, s.edits) };
        }

        // The table of these rows (taxid, gene, its alleles); a later row of the same copy replaces an earlier one.
        static Table FromRows(std::vector<std::tuple<uint32_t, uint32_t, std::vector<Allele>>> rows) {
            std::stable_sort(rows.begin(), rows.end(), [](auto const& a, auto const& b) {
                return std::get<0>(a) != std::get<0>(b) ? std::get<0>(a) < std::get<0>(b) : std::get<1>(a) < std::get<1>(b);
            });
            Table t;
            uint32_t const max_taxid = rows.empty() ? 0 : std::get<0>(rows.back());
            t.m_begin.assign(rows.empty() ? 0 : static_cast<size_t>(max_taxid) + 2, 0);
            for (size_t r = 0; r < rows.size(); r++) {
                auto const& [taxid, gene, alleles] = rows[r];
                if (gene > UINT16_MAX || alleles.empty()) continue;
                if (r + 1 < rows.size() && std::get<0>(rows[r + 1]) == taxid && std::get<1>(rows[r + 1]) == gene) continue;
                t.m_entries.push_back({ static_cast<uint16_t>(gene), static_cast<uint16_t>(std::min<size_t>(alleles.size(), UINT16_MAX)),
                                        static_cast<uint32_t>(t.m_alleles.size()) });
                for (auto const& a : alleles) {
                    t.m_alleles.push_back({ static_cast<uint32_t>(t.m_edits.size()),
                                            static_cast<uint16_t>(std::min<size_t>(a.edits.size(), UINT16_MAX)), a.begin, a.end });
                    t.m_edits.insert(t.m_edits.end(), a.edits.begin(), a.edits.begin() + std::min<size_t>(a.edits.size(), UINT16_MAX));
                }
                t.m_begin[taxid + 1]++;
            }
            for (size_t i = 1; i < t.m_begin.size(); i++) t.m_begin[i] += t.m_begin[i - 1];
            return t;
        }

        // strain_alleles.tsv: a header line, then per copy "taxid<TAB>gene<TAB>allele;allele;..." (WriteAllele), by taxid and
        // gene.
        void Write(std::ostream& os) const {
            os << "taxid\tgene\talleles (begin-end:edits; an edit is pos+base, pos i length or pos d length, on the representative)\n";
            for (size_t taxid = 0; taxid + 1 < m_begin.size(); taxid++) {
                for (size_t e = m_begin[taxid]; e < m_begin[taxid + 1]; e++) {
                    auto const& entry = m_entries[e];
                    os << taxid << '\t' << entry.gene << '\t';
                    for (size_t a = 0; a < entry.alleles; a++) {
                        if (a) os << ';';
                        auto const v = Get(m_alleles[entry.first_allele + a]);
                        WriteAllele(os, Allele{ v.begin, v.end, std::vector<Edit>(v.edits.begin(), v.edits.end()) });
                    }
                    os << '\n';
                }
            }
        }

        // Reads what Write writes; lines starting with '#' or "taxid" skipped. Returns an error text, or "".
        std::string Read(std::istream& is) {
            std::vector<std::tuple<uint32_t, uint32_t, std::vector<Allele>>> rows;
            std::string line;
            size_t number = 0;
            Allele allele;
            while (std::getline(is, line)) {
                number++;
                if (!line.empty() && line.back() == '\r') line.pop_back();
                if (line.empty() || line[0] == '#' || line.rfind("taxid", 0) == 0) continue;
                std::string_view view(line);
                size_t const tab1 = view.find('\t');
                size_t const tab2 = tab1 == std::string_view::npos ? tab1 : view.find('\t', tab1 + 1);
                uint32_t taxid = 0, gene = 0;
                auto parse = [](std::string_view s, uint32_t& v) {
                    auto const [p, ec] = std::from_chars(s.data(), s.data() + s.size(), v);
                    return ec == std::errc() && p == s.data() + s.size();
                };
                if (tab2 == std::string_view::npos || !parse(view.substr(0, tab1), taxid) ||
                    !parse(view.substr(tab1 + 1, tab2 - tab1 - 1), gene) || gene > UINT16_MAX) {
                    return "line " + std::to_string(number) + ": no taxid, gene and alleles";
                }
                std::vector<Allele> alleles;
                std::string_view rest = view.substr(tab2 + 1);
                while (!rest.empty()) {
                    size_t const semi = rest.find(';');
                    std::string_view const text = rest.substr(0, semi);
                    rest = semi == std::string_view::npos ? std::string_view() : rest.substr(semi + 1);
                    if (!ParseAllele(text, allele)) {
                        return "line " + std::to_string(number) + ": '" + std::string(text) + "' is not begin-end:edits";
                    }
                    alleles.push_back(allele);
                }
                if (alleles.empty()) return "line " + std::to_string(number) + ": no allele";
                rows.emplace_back(taxid, gene, std::move(alleles));
            }
            *this = FromRows(std::move(rows));
            return "";
        }

    private:
        struct Entry {
            uint16_t gene = 0;
            uint16_t alleles = 0;
            uint32_t first_allele = 0;
        };
        std::vector<uint32_t> m_begin;  // by taxid: its entries are [m_begin[taxid], m_begin[taxid + 1])
        std::vector<Entry> m_entries;   // by taxid, then gene
        std::vector<Stored> m_alleles;
        std::vector<Edit> m_edits;
    };

    // A read's differences from the representative where it is aligned: substitutions (with the read's base), and its
    // insertions and deletions, by position on the representative; the aligned span [begin, end) and its columns.
    struct ReadDiffs {
        uint32_t begin = 0;
        uint32_t end = 0;
        uint32_t columns = 0;
        std::vector<Edit> diffs;  // a substitution's base is the read's (4: not ACGT, which no allele explains)
        std::vector<uint8_t> other_base;  // per diff: 1 if the read's base is not ACGT

        void Clear() {
            begin = end = columns = 0;
            diffs.clear();
            other_base.clear();
        }
        void Substitution(uint32_t pos, char read_base) {
            uint8_t const b = BaseCode(read_base);
            diffs.push_back(Edit::Substitution(static_cast<uint16_t>(std::min<uint32_t>(pos, UINT16_MAX)), b > 3 ? 0 : b));
            other_base.push_back(b > 3);
        }
        void Indel(Kind kind, uint32_t pos, uint32_t length) {
            diffs.push_back(Edit::Indel(kind, static_cast<uint16_t>(std::min<uint32_t>(pos, UINT16_MAX)),
                                        static_cast<uint16_t>(std::min<uint32_t>(length, Edit::kMaxLength))));
            other_base.push_back(0);
        }
    };

    // A SAM record's differences (CIGAR with =, M, X, I, D, N, S, H; 1-based `pos`; SEQ in the reference's orientation, or
    // "*": then no substitution is known and nothing is filled). A position past UINT16_MAX is no gene of the table.
    inline bool FromSamRecord(std::string_view cigar, size_t pos, std::string_view seq, ReadDiffs& out) {
        out.Clear();
        if (seq.empty() || seq == "*" || pos == 0) return false;
        uint32_t ref = static_cast<uint32_t>(pos - 1), query = 0, run = 0;
        out.begin = ref;
        for (char const op : cigar) {
            if (op >= '0' && op <= '9') {
                run = run * 10 + static_cast<uint32_t>(op - '0');
                continue;
            }
            switch (op) {
                case '=': case 'M':
                    ref += run;
                    query += run;
                    out.columns += run;
                    break;
                case 'X':
                    for (uint32_t i = 0; i < run; i++) {
                        if (query + i < seq.size()) out.Substitution(ref + i, seq[query + i]);
                    }
                    ref += run;
                    query += run;
                    out.columns += run;
                    break;
                case 'I':
                    out.Indel(kInsertion, ref, run);
                    query += run;
                    out.columns += run;
                    break;
                case 'D': case 'N':
                    out.Indel(kDeletion, ref, run);
                    ref += run;
                    out.columns += run;
                    break;
                case 'S':
                    query += run;
                    break;
                default:
                    break;
            }
            run = 0;
        }
        out.end = ref;
        return out.columns > 0;
    }

    // An alignment's differences from its per-column operations (M, X, I, D, S, H: AlignmentInfo::cigar), the read as it
    // was aligned (that strand) and the 0-based position of the first aligned reference base.
    inline void FromColumns(std::string_view ops, std::string_view read, uint32_t ref_start, ReadDiffs& out) {
        out.Clear();
        uint32_t ref = ref_start, query = 0;
        out.begin = ref;
        for (size_t i = 0; i < ops.size();) {
            char const op = ops[i];
            size_t j = i + 1;
            while (j < ops.size() && ops[j] == op) j++;
            uint32_t const run = static_cast<uint32_t>(j - i);
            switch (op) {
                case 'M': case '=':
                    ref += run;
                    query += run;
                    out.columns += run;
                    break;
                case 'X':
                    for (uint32_t k = 0; k < run; k++) {
                        if (query + k < read.size()) out.Substitution(ref + k, read[query + k]);
                    }
                    ref += run;
                    query += run;
                    out.columns += run;
                    break;
                case 'I':
                    out.Indel(kInsertion, ref, run);
                    query += run;
                    out.columns += run;
                    break;
                case 'D':
                    out.Indel(kDeletion, ref, run);
                    ref += run;
                    out.columns += run;
                    break;
                case 'S':
                    query += run;
                    break;
                default:
                    break;
            }
            i = j;
        }
        out.end = ref;
    }

    // What an allele makes of a read: the read's differences it explains, and its own differences inside the read's
    // span that the read lacks. `shift` is the change in the read's differences if the allele were the reference.
    struct Explained {
        uint32_t explained = 0;
        uint32_t contradicted = 0;
        int shift() const { return static_cast<int>(contradicted) - static_cast<int>(explained); }
    };

    inline bool Matches(Edit const& a, Edit const& r) {
        if (a.GetKind() != r.GetKind()) return false;
        if (a.GetKind() == kSubstitution) return a.pos == r.pos && a.Base() == r.Base();
        uint32_t const d = a.pos > r.pos ? a.pos - r.pos : r.pos - a.pos;
        return d <= kIndelTolerance && a.Length() == r.Length();
    }

    inline Explained Explain(Table::View const& allele, ReadDiffs const& read) {
        Explained x;
        uint32_t const lo = std::max<uint32_t>(read.begin, allele.begin);
        uint32_t const hi = std::min<uint32_t>(read.end, allele.end);
        if (lo >= hi) return x;
        // Both lists are by position; a read's diffs are few (a handful on a short read), an allele's within the span too.
        auto const first = std::lower_bound(allele.edits.begin(), allele.edits.end(), lo,
                                            [](Edit const& e, uint32_t p) { return e.pos < p; });
        size_t const n = read.diffs.size();
        // A read's diff is used once: a flag per diff (a bit up to 64 diffs, a short read's; a vector for a long read's).
        uint64_t used = 0;
        std::vector<uint8_t> used_many;
        if (n > 64) used_many.assign(n, 0);
        auto is_used = [&](size_t i) { return n > 64 ? used_many[i] != 0 : (used >> i & 1) != 0; };
        auto set_used = [&](size_t i) { if (n > 64) used_many[i] = 1; else used |= uint64_t{1} << i; };
        for (auto it = first; it != allele.edits.end() && it->pos < hi; ++it) {
            bool const substitution = it->GetKind() == kSubstitution;
            // An indel at the read's very ends (or a deletion running past them) is not inside it.
            if (!substitution && (it->pos <= read.begin || it->pos >= read.end)) continue;
            // The read's diffs are by position: those within reach of this edit (its own position for a substitution).
            uint32_t const from = substitution ? it->pos : (it->pos > kIndelTolerance ? it->pos - kIndelTolerance : 0);
            uint32_t const to = substitution ? it->pos : it->pos + kIndelTolerance;
            auto r = std::lower_bound(read.diffs.begin(), read.diffs.end(), from, [](Edit const& e, uint32_t p) { return e.pos < p; });
            bool found = false;
            for (; r != read.diffs.end() && r->pos <= to; ++r) {
                size_t const i = static_cast<size_t>(r - read.diffs.begin());
                if (is_used(i) || read.other_base[i]) continue;
                if (Matches(*it, *r)) {
                    set_used(i);
                    found = true;
                    break;
                }
            }
            if (found) x.explained++;
            else x.contradicted++;
        }
        return x;
    }

    // The best allele for a read: the least shift below 0 (of equals the first); allele -1 and shift 0 when none beats
    // the representative.
    struct Best {
        int allele = -1;
        int shift = 0;
        uint32_t explained = 0;
    };

    inline Best BestAllele(Table const& table, uint32_t taxid, uint32_t gene, ReadDiffs const& read) {
        Best best;
        auto const alleles = table.Of(taxid, gene);
        for (size_t a = 0; a < alleles.size(); a++) {
            auto const x = Explain(table.Get(alleles[a]), read);
            if (x.shift() < best.shift) {
                best.allele = static_cast<int>(a);
                best.shift = x.shift();
                best.explained = x.explained;
            }
        }
        return best;
    }
}
