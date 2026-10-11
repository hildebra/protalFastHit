// NcbiTaxids.h - each species' NCBI species taxid, for the profiles of --ncbi_taxids. scripts/mini_db/gtdb_to_protal_db.py
// writes ncbi_taxids.tsv from the GTDB release's metadata: per GTDB species the NCBI species taxid most of its genomes
// have (a tie: the representative's), NCBI's name of it, whether that is a named species (0: one of NCBI's placeholders,
// "bacterium", "uncultured Lachnospiraceae bacterium", "Streptomyces sp. IB2014 016-6", which NCBI files under the
// deepest rank it can name: 82.5% of the GTDB r226 species have only such a taxid), the genomes with it, the species'
// genomes, and how many species of the release share it (docs/claude/2026-10-11-ncbi-taxids). --build packs it into
// database.protal, --add_tables stores it in a database converted before 2026-10-11 (gtdb_to_protal_db.py --priors_only
// writes it), and a run with --ncbi_taxids loads it.
#pragma once

#include <algorithm>
#include <charconv>
#include <cstdint>
#include <istream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace protal::ncbi_taxids {
    inline const std::string kFileName = "ncbi_taxids.tsv";

    struct Row {
        uint32_t ncbi_taxid = 0;
        std::string name;    // NCBI's name of the taxid
        bool named = false;  // a named NCBI species, not a placeholder
    };

    class Table {
    public:
        bool Empty() const { return m_rows.empty(); }
        size_t Size() const { return m_rows.size(); }

        // The species' row, or nullptr if the table has none for it.
        Row const* Get(uint32_t taxid) const {
            auto const it = m_rows.find(taxid);
            return it == m_rows.end() ? nullptr : &it->second;
        }

        // The species whose taxid is a named NCBI species.
        size_t Named() const {
            size_t n = 0;
            for (auto const& [_, row] : m_rows) n += row.named;
            return n;
        }

        // ncbi_taxids.tsv: a header line starting with "taxid" that names the columns taxid, ncbi_taxid, ncbi_name and named
        // (in any order; others, such as rep_genome, species, ncbi_genomes, genomes and shared_by, are skipped), then a
        // line per species (tab-separated; lines starting with '#' are skipped). Returns an error text, or "" when every
        // line was read.
        std::string Read(std::istream& is) {
            std::string line;
            size_t number = 0;
            std::vector<std::string_view> fields;
            enum Column { kTaxid, kNcbiTaxid, kName, kNamed, kColumns };
            constexpr std::string_view names[kColumns] = { "taxid", "ncbi_taxid", "ncbi_name", "named" };
            size_t column[kColumns] = { SIZE_MAX, SIZE_MAX, SIZE_MAX, SIZE_MAX };
            size_t width = 0;  // fields a line needs: up to the last of the columns read
            auto id = [&number](std::string_view f, uint32_t& v, std::string_view what) -> std::string {
                if (auto const [p, ec] = std::from_chars(f.data(), f.data() + f.size(), v); ec != std::errc() || p != f.data() + f.size() || v == 0) {
                    return "line " + std::to_string(number) + ": " + std::string(what) + " '" + std::string(f) + "' is not a positive number";
                }
                return "";
            };
            while (std::getline(is, line)) {
                number++;
                if (!line.empty() && line.back() == '\r') line.pop_back();
                if (line.empty() || line[0] == '#') continue;
                fields.clear();
                std::string_view const view(line);
                size_t start = 0;
                for (size_t tab; (tab = view.find('\t', start)) != std::string_view::npos; start = tab + 1) fields.push_back(view.substr(start, tab - start));
                fields.push_back(view.substr(start));
                if (width == 0) {
                    if (line.rfind("taxid", 0) != 0) return "line " + std::to_string(number) + ": no header line (taxid, ncbi_taxid, ncbi_name, named, ...)";
                    for (size_t i = 0; i < fields.size(); i++) {
                        for (size_t c = 0; c < kColumns; c++) {
                            if (fields[i] == names[c]) column[c] = i;
                        }
                    }
                    for (size_t c = 0; c < kColumns; c++) {
                        if (column[c] == SIZE_MAX) return "the header names no column " + std::string(names[c]);
                        width = std::max(width, column[c] + 1);
                    }
                    continue;
                }
                if (fields.size() < width) {
                    return "line " + std::to_string(number) + " has " + std::to_string(fields.size()) + " fields, the header's columns need " +
                           std::to_string(width);
                }
                uint32_t taxid = 0;
                Row row;
                if (auto error = id(fields[column[kTaxid]], taxid, "taxid"); !error.empty()) return error;
                if (auto error = id(fields[column[kNcbiTaxid]], row.ncbi_taxid, "ncbi_taxid"); !error.empty()) return error;
                auto const named = fields[column[kNamed]];
                if (named != "0" && named != "1") return "line " + std::to_string(number) + ": named is '" + std::string(named) + "', not 0 or 1";
                row.named = named == "1";
                row.name = std::string(fields[column[kName]]);
                if (!m_rows.emplace(taxid, std::move(row)).second) return "line " + std::to_string(number) + ": taxid " + std::to_string(taxid) + " twice";
            }
            if (width == 0) return "the table is empty (no header line)";
            return "";
        }

    private:
        std::unordered_map<uint32_t, Row> m_rows;
    };
}
