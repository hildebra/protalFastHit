// ForeignRatesTable.h - foreign_rates.tsv: per species' copy of a marker gene, how many reads of a tiled scan of the full
// reference's marker genes (every genome's, a few genomes per species) landed on it (their best record, as the profiler
// counts them), and how many of those came from a genome of another species, and of another genus. scripts/foreign_rates.py
// makes it (simulate_metagenomes --tiles --tile_fasta, one protal run against the database) and protal --add_tables packs
// it into the database; a run loads it for the profiler's "foreign" features: a copy that other species' reads reach (a
// gene conserved across the genus, a strain of a congener that crosses the species boundary on it) is weaker evidence of
// its species than one that only its own reads reach (docs/claude/2026-10-07-error-read-signatures, section 7 and the
// follow-up). Every copy of the full reference has a row, reached by a read or not (reads 0): until 2026-10-08 the scan
// tiled the genomes the training samples are drawn from, and a copy was listed only when a read reached it, which told
// the models which species the simulation could draw (docs/claude/2026-10-07-congener-gaps, "The foreign features
// leak"). A copy's reads are its foreign reads plus one genome's worth of its own species' reads (not every own read:
// those scale with the species' genome count, the cluster-size rule the simulation cannot test), so ForeignShare is
// foreign / (foreign + own per genome + 1). A training database's full reference lacks its held-out species, so the
// table knows nothing of the species the models are trained to find missing.
#pragma once

#include <algorithm>
#include <charconv>
#include <cstdint>
#include <istream>
#include <ostream>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace protal::foreign_rates {
    inline const std::string kFileName = "foreign_rates.tsv";
    inline constexpr double kUnknown = -1;  // features without the table, or without a read on a scanned copy

    // A copy's reads in the scan (the foreign ones plus one genome's worth of its own; saturated at UINT16_MAX), of them
    // from other species and from other genera.
    struct Rate {
        uint16_t reads = 0;
        uint16_t foreign = 0;
        uint16_t foreign_genus = 0;
        bool operator==(Rate const&) const = default;
        // The shares, shrunk towards 0 for copies with few reads: foreign / (reads + 1).
        double ForeignShare() const { return static_cast<double>(foreign) / (static_cast<double>(reads) + 1); }
        double ForeignGenusShare() const { return static_cast<double>(foreign_genus) / (static_cast<double>(reads) + 1); }
    };

    inline uint16_t Saturated(uint64_t n) { return static_cast<uint16_t>(std::min<uint64_t>(n, UINT16_MAX)); }

    class Table {
    public:
        struct Entry {
            uint16_t gene = 0;
            Rate rate;
        };

        bool Empty() const { return m_entries.empty(); }
        size_t Copies() const { return m_entries.size(); }
        size_t Species() const {
            size_t n = 0;
            for (size_t t = 0; t + 1 < m_begin.size(); t++) n += m_begin[t + 1] > m_begin[t];
            return n;
        }

        // Copy (taxid, gene)'s rate, or nullptr (no read of the scan reached it).
        Rate const* Find(uint32_t taxid, uint32_t gene) const {
            if (static_cast<size_t>(taxid) + 1 >= m_begin.size() || gene > UINT16_MAX) return nullptr;
            auto const first = m_entries.begin() + m_begin[taxid];
            auto const last = m_entries.begin() + m_begin[taxid + 1];
            auto const it = std::lower_bound(first, last, gene, [](Entry const& e, uint32_t g) { return e.gene < g; });
            return it != last && it->gene == gene ? &it->rate : nullptr;
        }

        // The table of these rows (taxid, gene, rate); a later row of the same copy replaces an earlier one.
        static Table FromRows(std::vector<std::tuple<uint32_t, uint32_t, Rate>> rows) {
            std::stable_sort(rows.begin(), rows.end(), [](auto const& a, auto const& b) {
                return std::get<0>(a) != std::get<0>(b) ? std::get<0>(a) < std::get<0>(b) : std::get<1>(a) < std::get<1>(b);
            });
            Table t;
            uint32_t const max_taxid = rows.empty() ? 0 : std::get<0>(rows.back());
            t.m_begin.assign(rows.empty() ? 0 : static_cast<size_t>(max_taxid) + 2, 0);
            for (size_t r = 0; r < rows.size(); r++) {
                auto const& [taxid, gene, rate] = rows[r];
                if (gene > UINT16_MAX) continue;  // marker gene ids are small
                if (r + 1 < rows.size() && std::get<0>(rows[r + 1]) == taxid && std::get<1>(rows[r + 1]) == gene) continue;
                t.m_entries.push_back({ static_cast<uint16_t>(gene), rate });
                t.m_begin[taxid + 1]++;
            }
            for (size_t i = 1; i < t.m_begin.size(); i++) t.m_begin[i] += t.m_begin[i - 1];
            return t;
        }

        // foreign_rates.tsv: comment lines (#), a header line, then per species
        // "taxid<TAB>gene:reads:foreign:foreign_genus,...", genes ascending.
        void Write(std::ostream& os) const {
            os << "taxid\tgene:reads:foreign:foreign_genus\n";
            for (size_t taxid = 0; taxid + 1 < m_begin.size(); taxid++) {
                if (m_begin[taxid + 1] == m_begin[taxid]) continue;
                os << taxid << '\t';
                for (size_t e = m_begin[taxid]; e < m_begin[taxid + 1]; e++) {
                    auto const& [gene, rate] = m_entries[e];
                    if (e > m_begin[taxid]) os << ',';
                    os << gene << ':' << rate.reads << ':' << rate.foreign << ':' << rate.foreign_genus;
                }
                os << '\n';
            }
        }

        // Reads what Write writes (lines starting with '#' or "taxid" skipped; counts above UINT16_MAX saturate).
        // Returns an error text, or "".
        std::string Read(std::istream& is) {
            std::vector<std::tuple<uint32_t, uint32_t, Rate>> rows;
            std::string line;
            size_t number = 0;
            while (std::getline(is, line)) {
                number++;
                if (!line.empty() && line.back() == '\r') line.pop_back();
                if (line.empty() || line[0] == '#' || line.rfind("taxid", 0) == 0) continue;
                std::string_view view(line);
                size_t const tab = view.find('\t');
                uint32_t taxid = 0;
                std::string_view const id_text = view.substr(0, tab);
                if (auto const [p, ec] = std::from_chars(id_text.data(), id_text.data() + id_text.size(), taxid);
                    ec != std::errc() || p != id_text.data() + id_text.size() || tab == std::string_view::npos) {
                    return "line " + std::to_string(number) + ": no taxid and genes";
                }
                std::string_view rest = view.substr(tab + 1);
                while (!rest.empty()) {
                    size_t const comma = rest.find(',');
                    std::string_view const entry = rest.substr(0, comma);
                    rest = comma == std::string_view::npos ? std::string_view() : rest.substr(comma + 1);
                    uint64_t values[4] = {0, 0, 0, 0};
                    char const* p = entry.data();
                    char const* const end = entry.data() + entry.size();
                    for (size_t v = 0; v < 4; v++) {
                        auto const r = std::from_chars(p, end, values[v]);
                        if (r.ec != std::errc() || (v < 3 && (r.ptr == end || *r.ptr != ':')) || (v == 3 && r.ptr != end)) {
                            return "line " + std::to_string(number) + ": '" + std::string(entry) + "' is not gene:reads:foreign:foreign_genus";
                        }
                        p = r.ptr + (v < 3 ? 1 : 0);
                    }
                    if (values[0] > UINT16_MAX || values[2] > values[1] || values[3] > values[2]) {
                        return "line " + std::to_string(number) + ": '" + std::string(entry) +
                               "' is out of range (foreign_genus <= foreign <= reads)";
                    }
                    rows.emplace_back(taxid, static_cast<uint32_t>(values[0]),
                                      Rate{ Saturated(values[1]), Saturated(values[2]), Saturated(values[3]) });
                }
            }
            *this = FromRows(std::move(rows));
            return "";
        }

    private:
        std::vector<uint32_t> m_begin;  // by taxid: its entries are [m_begin[taxid], m_begin[taxid + 1])
        std::vector<Entry> m_entries;   // by taxid, then gene
    };
}
