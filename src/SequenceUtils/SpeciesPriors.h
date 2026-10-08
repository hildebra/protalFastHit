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
// Since 2026-10-08 the table also gives each species' genome size (the mean of its GTDB genomes' assembly sizes, corrected
// for CheckM completeness and contamination), the representative's assembly size and the bases of its marker genes, by
// the header's names after the ten columns above: the profile's composition (Composition.h: the reads the called species
// explain, the sample's average genome size, its unknown share) uses them. A table without them gives unknown sizes.
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
        double genome_size = kUnknown;         // bp: the mean of its genomes' sizes, corrected for completeness and contamination
        double rep_genome_size = kUnknown;     // bp: the representative's assembly
        double marker_bases = kUnknown;        // bp: the representative's marker genes (a copy of each)
    };

    // The columns a table may have after the first ten, by their header names.
    inline constexpr std::string_view kGenomeSizeColumn = "genome_size";
    inline constexpr std::string_view kRepGenomeSizeColumn = "rep_genome_size";
    inline constexpr std::string_view kMarkerBasesColumn = "marker_bases";

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
        // (tab-separated; lines starting with '#' are skipped; -1, "none" or an empty field: unknown), and the genome
        // sizes' columns where the header names them (kGenomeSizeColumn, ...; at any place after the tenth). Returns an
        // error text, or "" when every line was read.
        std::string Read(std::istream& is) {
            std::string line;
            size_t number = 0;
            std::vector<std::string_view> fields;
            auto split = [&fields](std::string_view view) {
                fields.clear();
                size_t start = 0;
                for (size_t tab; (tab = view.find('\t', start)) != std::string_view::npos; start = tab + 1) fields.push_back(view.substr(start, tab - start));
                fields.push_back(view.substr(start));
            };
            // The columns of the genome sizes, by the header (SIZE_MAX: the table has none).
            size_t genome_size_column = SIZE_MAX, rep_size_column = SIZE_MAX, marker_bases_column = SIZE_MAX;
            auto value = [&number](std::string_view f, double& v) -> std::string {
                if (f.empty() || f == "none" || f == "-1" || f == "N/A" || f == "NA") return "";
                if (auto const [p, ec] = std::from_chars(f.data(), f.data() + f.size(), v); ec != std::errc() || p != f.data() + f.size()) {
                    return "line " + std::to_string(number) + ": '" + std::string(f) + "' is not a number";
                }
                return "";
            };
            while (std::getline(is, line)) {
                number++;
                if (!line.empty() && line.back() == '\r') line.pop_back();
                if (line.empty() || line[0] == '#') continue;
                split(line);
                if (line.rfind("taxid", 0) == 0) {
                    for (size_t i = 10; i < fields.size(); i++) {
                        if (fields[i] == kGenomeSizeColumn) genome_size_column = i;
                        if (fields[i] == kRepGenomeSizeColumn) rep_size_column = i;
                        if (fields[i] == kMarkerBasesColumn) marker_bases_column = i;
                    }
                    continue;
                }
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
                    if (auto error = value(fields[i + 2], *values[i]); !error.empty()) return error;
                }
                std::pair<size_t, double*> const sizes[] = { { genome_size_column, &row.genome_size }, { rep_size_column, &row.rep_genome_size },
                                                             { marker_bases_column, &row.marker_bases } };
                for (auto const& [column, target] : sizes) {
                    if (column == SIZE_MAX) continue;
                    if (column >= fields.size()) {
                        return "line " + std::to_string(number) + " has " + std::to_string(fields.size()) + " fields, but the header names " +
                               std::to_string(column + 1);
                    }
                    if (auto error = value(fields[column], *target); !error.empty()) return error;
                    if (*target != kUnknown && *target <= 0) *target = kUnknown;  // a size of 0 says nothing
                }
                m_rows[taxid] = row;
            }
            return "";
        }

        // The species with a genome size.
        size_t WithGenomeSize() const {
            size_t n = 0;
            for (auto const& [_, r] : m_rows) n += r.genome_size > 0;
            return n;
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
