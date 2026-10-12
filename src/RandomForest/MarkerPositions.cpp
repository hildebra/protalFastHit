// SPDX-License-Identifier: GPL-2.0-only
#include "MarkerPositions.h"

#include <algorithm>
#include <charconv>
#include <fstream>
#include <numeric>
#include <stdexcept>
#include <string_view>

namespace protal::sim {

namespace {

constexpr std::string_view kInsilico = "insilico_";  // insilico_strains.py's PREFIX

std::string Key(std::string_view accession, std::string_view contig) {
    std::string key;
    key.reserve(accession.size() + 1 + contig.size());
    key.append(accession);
    key += '\t';
    key.append(contig);
    return key;
}

}  // namespace

std::string MarkerPositions::AccessionOf(std::string const& genome) {
    if (!genome.starts_with(kInsilico)) return genome;
    std::string accession = genome.substr(kInsilico.size());
    auto const dot = accession.rfind('_');  // strain_name: the accession's . as _
    if (dot != std::string::npos) accession[dot] = '.';
    return accession;
}

std::shared_ptr<MarkerPositions const> MarkerPositions::Load(std::filesystem::path const& file,
                                                             std::unordered_set<std::string> const& accessions,
                                                             std::uint32_t min_overlap) {
    std::ifstream in(file);
    if (!in) throw std::runtime_error("cannot read the gene positions " + file.string());
    auto positions = std::make_shared<MarkerPositions>();
    positions->m_min_overlap = std::max<std::uint32_t>(1, min_overlap);
    std::size_t accession = 0, contig = 0, start = 0, end = 0, columns = 0;
    bool header = false;
    std::string line;
    std::vector<std::string_view> fields;
    std::size_t line_no = 0;
    while (std::getline(in, line)) {
        ++line_no;
        if (line.empty() || line[0] == '#') continue;
        fields.clear();
        std::string_view rest(line);
        while (true) {
            auto const tab = rest.find('\t');
            fields.push_back(rest.substr(0, tab));
            if (tab == std::string_view::npos) break;
            rest.remove_prefix(tab + 1);
        }
        if (!header) {
            auto column = [&](std::string_view name) {
                auto const it = std::find(fields.begin(), fields.end(), name);
                if (it == fields.end()) {
                    throw std::runtime_error(file.string() + ": no column " + std::string(name) + " in its header");
                }
                return static_cast<std::size_t>(it - fields.begin());
            };
            accession = column("accession");
            contig = column("contig");
            start = column("start");
            end = column("end");
            columns = std::max({accession, contig, start, end}) + 1;
            header = true;
            continue;
        }
        if (fields.size() < columns) {
            throw std::runtime_error(file.string() + " line " + std::to_string(line_no) + ": fewer columns than its header");
        }
        std::string const acc(fields[accession]);
        if (!accessions.empty() && !accessions.count(acc)) continue;
        std::uint64_t first = 0, last = 0;
        auto parse = [&](std::string_view text, std::uint64_t& value) {
            auto const [ptr, ec] = std::from_chars(text.data(), text.data() + text.size(), value);
            return ec == std::errc() && ptr == text.data() + text.size();
        };
        if (!parse(fields[start], first) || !parse(fields[end], last) || first < 1 || last < first) {
            throw std::runtime_error(file.string() + " line " + std::to_string(line_no) + ": start and end must be "
                                     "positions (1-based, start <= end)");
        }
        Contig& c = positions->m_contigs[Key(acc, fields[contig])];
        c.starts.push_back(first - 1);  // [first - 1, last): 0-based, half-open
        c.ends.push_back(last);
        positions->m_accessions.insert(acc);
        ++positions->m_genes;
    }
    if (!header) throw std::runtime_error(file.string() + ": no header line");
    for (auto& [_, c] : positions->m_contigs) {
        std::vector<std::size_t> order(c.starts.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
            return c.starts[a] < c.starts[b] || (c.starts[a] == c.starts[b] && c.ends[a] < c.ends[b]);
        });
        Contig sorted;
        for (std::size_t i : order) {
            sorted.starts.push_back(c.starts[i]);
            sorted.ends.push_back(c.ends[i]);
            sorted.max_ends.push_back(std::max(c.ends[i], sorted.max_ends.empty() ? 0 : sorted.max_ends.back()));
        }
        c = std::move(sorted);
    }
    return positions;
}

bool MarkerPositions::Has(std::string const& genome) const {
    return m_accessions.count(AccessionOf(genome)) > 0;
}

bool MarkerPositions::Overlaps(std::string const& genome, std::string const& contig, std::uint64_t start,
                               std::uint64_t end) const {
    if (end <= start) return false;
    std::string_view name(contig);
    if (genome.starts_with(kInsilico) && name.starts_with(genome) && name.size() > genome.size() &&
        name[genome.size()] == '_') {
        name.remove_prefix(genome.size() + 1);  // the representative's contig
    }
    auto const it = m_contigs.find(Key(AccessionOf(genome), name));
    if (it == m_contigs.end()) return false;
    Contig const& c = it->second;
    std::uint64_t const need = std::min<std::uint64_t>(m_min_overlap, end - start);
    // Genes starting at end - need or later cannot share `need` bases; of those before, only the ones whose largest
    // end so far reaches start + need can.
    auto i = static_cast<std::size_t>(std::upper_bound(c.starts.begin(), c.starts.end(), end - need) - c.starts.begin());
    while (i > 0) {
        --i;
        if (c.max_ends[i] < start + need) break;
        std::uint64_t const from = std::max(start, c.starts[i]), to = std::min(end, c.ends[i]);
        if (to > from && to - from >= need) return true;
    }
    return false;
}

}  // namespace protal::sim
