// CongenerGapsTable.h - the table of congener_gaps.tsv (CongenerGaps.h makes it at --build): per species' copy of a
// marker gene, the alignment distance to its nearest congener's copy and to a typical one, and which congener's copy is
// the nearest. A run loads it into the GenomeLoader for the profiler's "gaps" features, and the ancestry sites
// (AncestrySites.h) compare a copy with that nearest congener's copy of the gene; this header holds only the table, so
// that the options and the loader need not include the aligner.
#pragma once

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <istream>
#include <optional>
#include <ostream>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace protal::congener_gaps {
    inline const std::string kFileName = "congener_gaps.tsv";
    inline constexpr uint32_t kScale = 10000;            // distances are stored in 1/kScale (0.0001)
    inline constexpr double kUnknown = -1;               // features without the table, or without an informative read
    // A copy whose nearest congener's copy is nearer than this (in 1/kScale: 0.005) leaves a read nothing to tell by: the
    // profiler counts no read on it for the gaps features.
    inline constexpr uint16_t kMinGap = 50;

    // A copy's gap: the nearest congener's distance and the median distance, in 1/kScale, how many congeners of the genus
    // have the gene, and the nearest one (its taxid; 0: not known, an older table of four fields).
    struct Gap {
        uint16_t min = 0;
        uint16_t median = 0;
        uint16_t congeners = 0;
        uint32_t nearest = 0;
        bool operator==(Gap const&) const = default;
        double Min() const { return static_cast<double>(min) / kScale; }
        double Median() const { return static_cast<double>(median) / kScale; }
    };

    inline uint16_t Scaled(double distance) {
        double const d = std::clamp(distance, 0.0, 1.0);
        return static_cast<uint16_t>(std::lround(d * kScale));
    }

    // The table of every copy's gap, by taxid and gene.
    class Table {
    public:
        // 12 bytes a copy (~175 MB at r226): the gene and the gap's fields side by side, without a Gap's padding.
        struct Entry {
            uint32_t nearest = 0;
            uint16_t gene = 0;
            uint16_t min = 0;
            uint16_t median = 0;
            uint16_t congeners = 0;
            Gap ToGap() const { return Gap{ min, median, congeners, nearest }; }
        };

        bool Empty() const { return m_entries.empty(); }
        size_t Copies() const { return m_entries.size(); }
        size_t Species() const {
            size_t n = 0;
            for (size_t t = 0; t + 1 < m_begin.size(); t++) n += m_begin[t + 1] > m_begin[t];
            return n;
        }

        // Copy (taxid, gene)'s gap, or nullopt.
        std::optional<Gap> Find(uint32_t taxid, uint32_t gene) const {
            if (static_cast<size_t>(taxid) + 1 >= m_begin.size() || gene > UINT16_MAX) return std::nullopt;
            auto const first = m_entries.begin() + m_begin[taxid];
            auto const last = m_entries.begin() + m_begin[taxid + 1];
            auto const it = std::lower_bound(first, last, gene, [](Entry const& e, uint32_t g) { return e.gene < g; });
            if (it == last || it->gene != gene) return std::nullopt;
            return it->ToGap();
        }

        // The table of these rows (taxid, gene, gap); a later row of the same copy replaces an earlier one.
        static Table FromRows(std::vector<std::tuple<uint32_t, uint32_t, Gap>> rows) {
            std::stable_sort(rows.begin(), rows.end(), [](auto const& a, auto const& b) {
                return std::get<0>(a) != std::get<0>(b) ? std::get<0>(a) < std::get<0>(b) : std::get<1>(a) < std::get<1>(b);
            });
            Table t;
            uint32_t const max_taxid = rows.empty() ? 0 : std::get<0>(rows.back());
            t.m_begin.assign(rows.empty() ? 0 : static_cast<size_t>(max_taxid) + 2, 0);
            for (size_t r = 0; r < rows.size(); r++) {
                auto const& [taxid, gene, gap] = rows[r];
                if (gene > UINT16_MAX) continue;  // marker gene ids are small
                if (r + 1 < rows.size() && std::get<0>(rows[r + 1]) == taxid && std::get<1>(rows[r + 1]) == gene) continue;
                t.m_entries.push_back({ gap.nearest, static_cast<uint16_t>(gene), gap.min, gap.median, gap.congeners });
                t.m_begin[taxid + 1]++;
            }
            for (size_t i = 1; i < t.m_begin.size(); i++) t.m_begin[i] += t.m_begin[i - 1];
            return t;
        }

        // congener_gaps.tsv: a header line, then per species "taxid<TAB>gene:min:median:congeners:nearest,..." with the
        // distances in 1/kScale, genes ascending.
        void Write(std::ostream& os) const {
            os << "taxid\tgene:min:median:congeners:nearest\n";
            for (size_t taxid = 0; taxid + 1 < m_begin.size(); taxid++) {
                if (m_begin[taxid + 1] == m_begin[taxid]) continue;
                os << taxid << '\t';
                for (size_t e = m_begin[taxid]; e < m_begin[taxid + 1]; e++) {
                    auto const& entry = m_entries[e];
                    if (e > m_begin[taxid]) os << ',';
                    os << entry.gene << ':' << entry.min << ':' << entry.median << ':' << entry.congeners << ':' << entry.nearest;
                }
                os << '\n';
            }
        }

        // Reads what Write writes, and the older tables without the nearest congener (gene:min:median:
        // congeners; nearest 0); lines starting with '#' or "taxid" skipped. Returns an error text, or "".
        std::string Read(std::istream& is) {
            std::vector<std::tuple<uint32_t, uint32_t, Gap>> rows;
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
                    uint32_t values[5] = {0, 0, 0, 0, 0};
                    char const* p = entry.data();
                    char const* const end = entry.data() + entry.size();
                    size_t fields = 0;
                    while (fields < 5) {
                        auto const r = std::from_chars(p, end, values[fields]);
                        if (r.ec != std::errc()) {
                            fields = 0;  // not a number where one belongs (also after a trailing colon)
                            break;
                        }
                        fields++;
                        p = r.ptr;
                        if (p == end) break;
                        if (*p != ':') {
                            fields = 0;
                            break;
                        }
                        p++;
                    }
                    if ((fields != 4 && fields != 5) || p != end) {
                        return "line " + std::to_string(number) + ": '" + std::string(entry) + "' is not gene:min:median:congeners[:nearest]";
                    }
                    if (values[0] > UINT16_MAX || values[1] > kScale || values[2] > kScale || values[3] > UINT16_MAX) {
                        return "line " + std::to_string(number) + ": '" + std::string(entry) + "' is out of range";
                    }
                    rows.emplace_back(taxid, values[0], Gap{ static_cast<uint16_t>(values[1]), static_cast<uint16_t>(values[2]),
                                                             static_cast<uint16_t>(values[3]), values[4] });
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
