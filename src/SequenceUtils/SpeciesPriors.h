// SpeciesPriors.h - what GTDB knows of a species before any read, for the model's prior features: how many of the
// representative's single-copy marker genes the converter found twice (CheckM's contamination signature: a
// contaminating contig's genes put every present organism's reads on the species), the representative's CheckM
// completeness and contamination, and from GTDB's species clusters its ANI circumscription radius, the mean and
// minimum ANI of its genomes to the representative and how many genomes it holds (a wide or crowded cluster makes a
// cloud of reads 2-4% from the reference a strain rather than a sister species). scripts/mini_db/gtdb_to_protal_db.py
// writes species_priors.tsv (taxid, rep_genome, markers, duplicate_markers, checkm_completeness, checkm_contamination,
// ani_radius, mean_intra_ani, min_intra_ani, clustered_genomes; -1: unknown), --build packs it into database.protal, a
// run loads it and TaxonFeatures emits the per-species constants. A database without it gives every species the
// unknown values (docs/claude/2026-10-03-false-positive-fixes).
#pragma once

#include <charconv>
#include <cstdint>
#include <istream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace protal::species_priors {
    inline const std::string kFileName = "species_priors.tsv";
    inline constexpr double kUnknown = -1;

    struct Row {
        double markers = kUnknown;             // marker genes found in the representative
        double duplicate_markers = kUnknown;   // of them, found more than once
        double completeness = kUnknown;        // CheckM (v2 where the release has it), percent
        double contamination = kUnknown;
        double ani_radius = kUnknown;          // the cluster's ANI circumscription radius, percent (95 for most)
        double mean_intra_ani = kUnknown;      // mean and minimum ANI of the cluster's genomes to the representative
        double min_intra_ani = kUnknown;
        double clustered_genomes = kUnknown;   // genomes in the cluster
    };

    class Table {
    public:
        bool Empty() const { return m_rows.empty(); }
        size_t Size() const { return m_rows.size(); }

        // The species' row, or one of unknown values.
        Row const& Get(uint32_t taxid) const {
            static Row const unknown;
            auto const it = m_rows.find(taxid);
            return it == m_rows.end() ? unknown : it->second;
        }

        void Set(uint32_t taxid, Row row) {
            m_rows[taxid] = row;
        }

        // species_priors.tsv: a header line starting with "taxid", then per species the ten columns above
        // (tab-separated; lines starting with '#' are skipped; -1, "none" or an empty field: unknown). Returns an
        // error text, or "" when every line was read.
        std::string Read(std::istream& is) {
            std::string line;
            size_t number = 0;
            std::vector<std::string_view> fields;
            while (std::getline(is, line)) {
                number++;
                if (!line.empty() && line.back() == '\r') line.pop_back();
                if (line.empty() || line[0] == '#' || line.rfind("taxid", 0) == 0) continue;
                fields.clear();
                std::string_view view(line);
                size_t start = 0;
                for (size_t tab; (tab = view.find('\t', start)) != std::string_view::npos; start = tab + 1) fields.push_back(view.substr(start, tab - start));
                fields.push_back(view.substr(start));
                if (fields.size() < 10) return "line " + std::to_string(number) + " has " + std::to_string(fields.size()) + " fields, not 10";
                uint32_t taxid = 0;
                if (auto const [p, ec] = std::from_chars(fields[0].data(), fields[0].data() + fields[0].size(), taxid);
                    ec != std::errc() || p != fields[0].data() + fields[0].size()) {
                    return "line " + std::to_string(number) + ": taxid '" + std::string(fields[0]) + "' is not a number";
                }
                Row row;
                double* const values[] = { &row.markers, &row.duplicate_markers, &row.completeness, &row.contamination,
                                           &row.ani_radius, &row.mean_intra_ani, &row.min_intra_ani, &row.clustered_genomes };
                for (size_t i = 0; i < 8; i++) {
                    std::string_view const f = fields[i + 2];
                    if (f.empty() || f == "none" || f == "-1" || f == "N/A" || f == "NA") continue;
                    double v = 0;
                    if (auto const [p, ec] = std::from_chars(f.data(), f.data() + f.size(), v); ec != std::errc() || p != f.data() + f.size()) {
                        return "line " + std::to_string(number) + ": '" + std::string(f) + "' is not a number";
                    }
                    *values[i] = v;
                }
                m_rows[taxid] = row;
            }
            return "";
        }

        // The species with a known value in any column but the markers.
        size_t Informative() const {
            size_t n = 0;
            for (auto const& [_, r] : m_rows) {
                n += r.duplicate_markers > 0 || r.completeness != kUnknown || r.contamination != kUnknown || r.ani_radius != kUnknown ||
                     r.min_intra_ani != kUnknown || r.clustered_genomes != kUnknown;
            }
            return n;
        }

    private:
        std::unordered_map<uint32_t, Row> m_rows;
    };
}
