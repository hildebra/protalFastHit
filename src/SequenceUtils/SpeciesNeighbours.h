// SpeciesNeighbours.h - each species' nearest congeners in the database, by the distance of their references' marker
// genes: the median Mash distance of the genes both references have (context::SketchedTaxonDistances, the distance
// relative_distance reads in a run). --build compares every two species of a genus and writes, per species, up to
// kMaxNeighbours of its congeners within kMaxDistance, nearest first, to species_neighbours.tsv in the database; a run
// loads it for two kinds of features (docs/claude/2026-10-06-false-positive-features):
//   - how crowded the database is around the species' reference (db_congeners_01, _02, _05: congeners within 0.01, 0.02
//     and 0.05; db_nearest_congener): a species with congeners near-identical on the marker genes is where the reads of
//     a congener the database lacks land, the soil scenarios' main false positives;
//   - the distance of two congeners a read fits, for unexpected_congener_fit_share: a read of the species fits a
//     congener at about that distance times its length more edits, so a read that fits a distant congener as well is
//     not the species' own.
// The table describes the database it was built with: a training database lacks its held-out species, so its species
// have fewer near congeners than in the finished database (as the reference's k-mer uniqueness, `ref`, has).
#pragma once

#include <algorithm>
#include <charconv>
#include <cstdint>
#include <istream>
#include <ostream>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace protal::species_neighbours {
    inline const std::string kFileName = "species_neighbours.tsv";
    inline constexpr size_t kMaxNeighbours = 16;   // congeners kept per species, the nearest
    inline constexpr double kMaxDistance = 0.15;   // and none farther
    inline constexpr double kUnknown = -1;         // the features of a database without the table

    struct Neighbour {
        uint32_t taxid = 0;
        float distance = 0;
        bool operator==(Neighbour const&) const = default;
    };

    class Table {
    public:
        bool Empty() const { return m_species == 0; }
        // The species with a row (also those without a congener within kMaxDistance).
        size_t Species() const { return m_species; }
        // The neighbours listed over all species.
        size_t Pairs() const {
            size_t n = 0;
            for (auto const& row : m_rows) n += row.size();
            return n;
        }

        // Species taxid's congeners, nearest first (empty: none within kMaxDistance, or not in the table).
        std::span<Neighbour const> Of(uint32_t taxid) const {
            if (taxid >= m_rows.size()) return {};
            return m_rows[taxid];
        }

        // Whether taxid has a row (the build compared it).
        bool Has(uint32_t taxid) const {
            return taxid < m_known.size() && m_known[taxid];
        }

        // The distance from taxid's reference to its congener's, if listed; else the least it can be: its farthest
        // listed congener's distance when the list is full, kMaxDistance otherwise.
        double DistanceAtLeast(uint32_t taxid, uint32_t congener) const {
            auto const row = Of(taxid);
            for (auto const& n : row) {
                if (n.taxid == congener) return n.distance;
            }
            return row.size() >= kMaxNeighbours ? static_cast<double>(row.back().distance) : kMaxDistance;
        }

        // Taxid's listed congeners within `distance` (at most kMaxNeighbours).
        size_t Within(uint32_t taxid, double distance) const {
            auto const row = Of(taxid);
            return static_cast<size_t>(std::count_if(row.begin(), row.end(), [distance](Neighbour const& n) { return n.distance <= distance; }));
        }

        // The nearest congener's distance; 1 without one within kMaxDistance.
        double Nearest(uint32_t taxid) const {
            auto const row = Of(taxid);
            return row.empty() ? 1.0 : static_cast<double>(row.front().distance);
        }

        // Sets taxid's row: its congeners with their distances (any order, at most one each); keeps the kMaxNeighbours
        // nearest within kMaxDistance, ties by taxid.
        void Set(uint32_t taxid, std::vector<Neighbour> neighbours) {
            neighbours.erase(std::remove_if(neighbours.begin(), neighbours.end(),
                                            [](Neighbour const& n) { return !(n.distance <= kMaxDistance); }),
                             neighbours.end());
            std::sort(neighbours.begin(), neighbours.end(), [](Neighbour const& a, Neighbour const& b) {
                return a.distance != b.distance ? a.distance < b.distance : a.taxid < b.taxid;
            });
            if (neighbours.size() > kMaxNeighbours) neighbours.resize(kMaxNeighbours);
            if (taxid >= m_rows.size()) {
                m_rows.resize(taxid + 1);
                m_known.resize(taxid + 1, false);
            }
            if (!m_known[taxid]) m_species++;
            m_known[taxid] = true;
            m_rows[taxid] = std::move(neighbours);
        }

        // species_neighbours.tsv: a header line, then per species "taxid<TAB>congener:distance,congener:distance"
        // (nearest first; the second field empty for a species without a congener within kMaxDistance).
        void Write(std::ostream& os) const {
            os << "taxid\tneighbours\n";
            for (size_t taxid = 0; taxid < m_rows.size(); taxid++) {
                if (!m_known[taxid]) continue;
                os << taxid << '\t';
                bool first = true;
                for (auto const& n : m_rows[taxid]) {
                    if (!first) os << ',';
                    first = false;
                    char buffer[32];
                    auto const [end, ec] = std::to_chars(buffer, buffer + sizeof(buffer), n.distance);
                    os << n.taxid << ':' << std::string_view(buffer, static_cast<size_t>(end - buffer));
                }
                os << '\n';
            }
        }

        // Reads what Write writes (lines starting with '#' or "taxid" skipped). Returns an error text, or "" when every
        // line was read.
        std::string Read(std::istream& is) {
            std::string line;
            size_t number = 0;
            while (std::getline(is, line)) {
                number++;
                if (!line.empty() && line.back() == '\r') line.pop_back();
                if (line.empty() || line[0] == '#' || line.rfind("taxid", 0) == 0) continue;
                std::string_view view(line);
                size_t const tab = view.find('\t');
                std::string_view const id_text = view.substr(0, tab);
                uint32_t taxid = 0;
                if (auto const [p, ec] = std::from_chars(id_text.data(), id_text.data() + id_text.size(), taxid);
                    ec != std::errc() || p != id_text.data() + id_text.size()) {
                    return "line " + std::to_string(number) + ": taxid '" + std::string(id_text) + "' is not a number";
                }
                std::vector<Neighbour> neighbours;
                std::string_view rest = tab == std::string_view::npos ? std::string_view() : view.substr(tab + 1);
                while (!rest.empty()) {
                    size_t const comma = rest.find(',');
                    std::string_view const entry = rest.substr(0, comma);
                    rest = comma == std::string_view::npos ? std::string_view() : rest.substr(comma + 1);
                    size_t const colon = entry.find(':');
                    Neighbour n;
                    if (colon == std::string_view::npos ||
                        std::from_chars(entry.data(), entry.data() + colon, n.taxid).ec != std::errc() ||
                        std::from_chars(entry.data() + colon + 1, entry.data() + entry.size(), n.distance).ec != std::errc()) {
                        return "line " + std::to_string(number) + ": '" + std::string(entry) + "' is not congener:distance";
                    }
                    neighbours.push_back(n);
                }
                Set(taxid, std::move(neighbours));
            }
            return "";
        }

    private:
        std::vector<std::vector<Neighbour>> m_rows;  // by taxid (protal's taxids are dense internal ids)
        std::vector<bool> m_known;
        size_t m_species = 0;
    };
}
