// GeneIncongruence.h - gene copies that do not belong to their species. --build compares every species' copy of each
// marker gene with the other species' copies (bottom sketches of 12-mers, as a run's congener distances) and flags a
// copy that is near-identical to a copy of another genus, family, order, class, phylum or domain while its own
// congeners' copies are farther, or it has none: a contaminating contig in a MAG, or a transferred gene. Reads on such
// a copy are no evidence of its species, and every present organism with that gene puts a perfect read on it: at GTDB
// r226 a fifth of the false positives were reads of a present species of another genus at identity 0.99, recurring
// per reference (docs/claude/2026-10-03-false-positive-anatomy). A query leaves records on suspect copies out
// (suspect_copies.tsv in the database; --keep_suspect_copies). gene_incongruence.tsv beside the database lists every
// near pair across genera, for the record and for reporting the references.
#pragma once

#include "GeneConservation.h"

#include <algorithm>
#include <array>
#include <charconv>
#include <cstdint>
#include <fstream>
#include <istream>
#include <optional>
#include <ostream>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace protal::gene_incongruence {
    inline const std::string kFileName = "suspect_copies.tsv";          // in the database: the copies a query leaves out
    inline const std::string kReportFileName = "gene_incongruence.tsv";  // beside it: every near pair across genera

    inline constexpr size_t kSketchSize = 128;              // bottom hashes per copy (gene_conservation::BottomSketch)
    inline constexpr double kDefaultSuspectDistance = 0.02;  // a copy within this of another genus's (--suspect_copy_distance)
    inline constexpr double kCongenerMargin = 0.02;          // and whose nearest congener's copy is this much farther; a copy
                                                             // without congeners: nearer to the partner than the partner's
                                                             // congeners by this much
    inline constexpr double kReportDistance = 0.05;          // pairs across genera reported in gene_incongruence.tsv
    inline constexpr size_t kProbeHashes = 8;                // a copy's smallest hashes looked up for candidates
    inline constexpr size_t kMaxBucket = 2048;               // copies sharing a hash beyond which the hash tells nothing
    inline constexpr uint32_t kGeneBits = 20;                // gene ids below 2^20 (as gene_conservation)

    // The deepest rank two species share: congeners (Genus) down to nothing but the domain, or None when their
    // domains differ or are unknown.
    enum class Rank : uint8_t { Genus = 0, Family, Order, Class, Phylum, Domain, None };

    inline constexpr std::array<std::string_view, 7> kRankNames = { "genus", "family", "order", "class", "phylum", "domain", "none" };

    inline std::string_view RankName(Rank rank) {
        return kRankNames[static_cast<size_t>(rank)];
    }

    inline std::optional<Rank> RankFromName(std::string_view name) {
        for (size_t i = 0; i < kRankNames.size(); i++) {
            if (kRankNames[i] == name) return static_cast<Rank>(i);
        }
        return std::nullopt;
    }

    // A species' ancestors by rank (taxids of internal_taxonomy.dmp; 0: unknown), genus first.
    struct Lineage {
        std::array<uint32_t, 6> ancestor{};

        uint32_t At(Rank rank) const {
            return ancestor[static_cast<size_t>(rank)];
        }
    };

    inline Rank SharedRank(Lineage const& a, Lineage const& b) {
        for (size_t r = 0; r < a.ancestor.size(); r++) {
            if (a.ancestor[r] != 0 && a.ancestor[r] == b.ancestor[r]) return static_cast<Rank>(r);
        }
        return Rank::None;
    }

    // Every node's lineage from internal_taxonomy.dmp (columns id, parent_id, external_id, name, rank, ...): the
    // nearest ancestor of each rank, the node itself included. Returns an error text, or "" when the file was read.
    inline std::string LineagesFromTaxonomy(std::istream& is, std::unordered_map<uint32_t, Lineage>& lineages) {
        std::unordered_map<uint32_t, uint32_t> parent;
        std::unordered_map<uint32_t, Rank> rank_of;
        std::string line;
        size_t nodes = 0;
        while (std::getline(is, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            std::vector<std::string_view> fields;
            std::string_view view(line);
            size_t start = 0;
            for (size_t tab; (tab = view.find('\t', start)) != std::string_view::npos; start = tab + 1) fields.push_back(view.substr(start, tab - start));
            fields.push_back(view.substr(start));
            if (fields.size() < 5) continue;
            uint32_t id = 0, up = 0;
            if (std::from_chars(fields[0].data(), fields[0].data() + fields[0].size(), id).ec != std::errc() ||
                std::from_chars(fields[1].data(), fields[1].data() + fields[1].size(), up).ec != std::errc()) continue;
            parent[id] = up;
            if (auto const rank = RankFromName(fields[4]); rank && *rank != Rank::None) rank_of[id] = *rank;
            nodes++;
        }
        if (nodes == 0) return "no taxonomy nodes";
        for (auto const& [id, _] : parent) {
            Lineage lineage;
            uint32_t node = id;
            for (int depth = 0; depth < 64; depth++) {
                if (auto const r = rank_of.find(node); r != rank_of.end() && lineage.ancestor[static_cast<size_t>(r->second)] == 0) {
                    lineage.ancestor[static_cast<size_t>(r->second)] = node;
                }
                auto const it = parent.find(node);
                if (it == parent.end() || it->second == node || it->second == 0) break;
                node = it->second;
            }
            lineages[id] = lineage;
        }
        return "";
    }

    inline std::string LineagesFromTaxonomy(std::string const& path, std::unordered_map<uint32_t, Lineage>& lineages) {
        std::ifstream is(path);
        if (!is) return "cannot read " + path;
        return LineagesFromTaxonomy(is, lineages);
    }

    // A species' copy of a gene, as its bottom sketch.
    struct Copy {
        uint32_t taxid = 0;
        std::vector<uint32_t> sketch;
    };

    // A suspect copy: near-identical to `partner`'s copy, of the shared rank `rank` (never Genus), at `distance`;
    // its nearest congener's copy at `congener_distance` (2 for none within kReportDistance).
    struct Suspect {
        uint32_t taxid = 0;
        uint32_t geneid = 0;
        uint32_t partner = 0;
        Rank rank = Rank::None;
        float distance = 2;
        float congener_distance = 2;
    };

    // A near pair of copies across genera.
    struct Pair {
        uint32_t geneid = 0;
        uint32_t taxid_a = 0;
        uint32_t taxid_b = 0;
        Rank rank = Rank::None;
        float distance = 2;
    };

    struct Result {
        std::vector<Suspect> suspects;  // by gene id, then taxid
        std::vector<Pair> pairs;        // by gene id, then taxids
        size_t copies = 0;              // copies compared
        size_t genes = 0;               // genes with two copies or more
        size_t candidates = 0;          // sketch pairs compared
    };

    // Compares the copies of every gene (copies_by_gene[geneid]): for each copy its nearest other copy among its
    // congeners and among the species of other genera (candidates: the copies sharing one of its kProbeHashes smallest
    // hashes, their distance from the sketches), on `threads` threads. A copy is suspect when a copy of another genus
    // (the partner) is within suspect_distance and: the copy's genus has other copies of the gene but its nearest one
    // is kCongenerMargin farther than the partner's (its own genus does not have this copy, the other does); or the
    // copy's genus has no other copy (a singleton genus) and the partner's genus does, with the partner's nearest
    // congener's copy kCongenerMargin farther than this copy (the copy lies inside the other genus's cluster, tighter
    // than the cluster itself: a slow gene shared by a young family does not qualify). The same for any thread count:
    // each gene is judged on its own and the results sorted.
    inline Result Scan(std::vector<std::vector<Copy>> const& copies_by_gene, std::unordered_map<uint32_t, Lineage> const& lineages,
                       double suspect_distance, int threads) {
        Result result;
        Lineage const no_lineage;
        auto lineage_of = [&](uint32_t taxid) -> Lineage const& {
            auto const it = lineages.find(taxid);
            return it == lineages.end() ? no_lineage : it->second;
        };
        size_t copies = 0, genes = 0, candidates = 0;
        #pragma omp parallel for schedule(dynamic) num_threads(std::max(threads, 1)) reduction(+:copies, genes, candidates)
        for (size_t geneid = 0; geneid < copies_by_gene.size(); geneid++) {
            auto const& copies_here = copies_by_gene[geneid];
            size_t const n = copies_here.size();
            copies += n;
            if (n < 2) continue;
            genes++;
            // Every (hash, copy) of every sketch, sorted: a hash's copies are a run.
            std::vector<std::pair<uint32_t, uint32_t>> entries;
            entries.reserve(n * kSketchSize);
            for (size_t i = 0; i < n; i++) {
                for (uint32_t const h : copies_here[i].sketch) entries.emplace_back(h, static_cast<uint32_t>(i));
            }
            std::sort(entries.begin(), entries.end());
            std::vector<Suspect> suspects;
            std::vector<Pair> pairs;
            // Copies of the gene per genus: a copy whose genus has others is judged against them.
            std::unordered_map<uint32_t, size_t> per_genus;
            std::vector<uint32_t> genus_of(n, 0);
            for (size_t i = 0; i < n; i++) {
                genus_of[i] = lineage_of(copies_here[i].taxid).At(Rank::Genus);
                if (genus_of[i] != 0) per_genus[genus_of[i]]++;
            }
            // Pass 1: each copy's nearest congener's copy and nearest other genus's copy (the partner).
            std::vector<float> congener(n, 2), foreign(n, 2);
            std::vector<size_t> partner(n, n);
            std::vector<Rank> foreign_rank(n, Rank::None);
            std::vector<uint32_t> found;
            for (size_t i = 0; i < n; i++) {
                auto const& copy = copies_here[i];
                found.clear();
                size_t const probes = std::min(kProbeHashes, copy.sketch.size());
                for (size_t p = 0; p < probes; p++) {
                    auto const [lo, hi] = std::equal_range(entries.begin(), entries.end(), std::make_pair(copy.sketch[p], uint32_t{0}),
                        [](auto const& a, auto const& b) { return a.first < b.first; });
                    if (static_cast<size_t>(hi - lo) > kMaxBucket) continue;
                    for (auto it = lo; it != hi; ++it) {
                        if (it->second != i) found.push_back(it->second);
                    }
                }
                std::sort(found.begin(), found.end());
                found.erase(std::unique(found.begin(), found.end()), found.end());
                candidates += found.size();
                auto const& lineage = lineage_of(copy.taxid);
                for (uint32_t const j : found) {
                    auto const& other = copies_here[j];
                    if (other.taxid == copy.taxid) continue;
                    auto const d = static_cast<float>(gene_conservation::SketchDistance(copy.sketch, other.sketch));
                    if (d > kReportDistance) continue;
                    Rank const rank = SharedRank(lineage, lineage_of(other.taxid));
                    if (rank == Rank::Genus) {
                        congener[i] = std::min(congener[i], d);
                        continue;
                    }
                    if (d < foreign[i] || (d == foreign[i] && other.taxid < copies_here[partner[i]].taxid)) {
                        foreign[i] = d;
                        partner[i] = j;
                        foreign_rank[i] = rank;
                    }
                    if (copy.taxid < other.taxid) pairs.push_back({ static_cast<uint32_t>(geneid), copy.taxid, other.taxid, rank, d });
                    else pairs.push_back({ static_cast<uint32_t>(geneid), other.taxid, copy.taxid, rank, d });
                }
            }
            // Pass 2: the verdicts.
            for (size_t i = 0; i < n; i++) {
                if (foreign[i] > suspect_distance || partner[i] >= n) continue;
                bool const has_congeners = genus_of[i] != 0 && per_genus[genus_of[i]] > 1;
                bool suspect;
                if (has_congeners) {
                    suspect = congener[i] > foreign[i] + kCongenerMargin;  // 2 when none within kReportDistance
                } else {
                    size_t const p = partner[i];
                    bool const partner_has_congeners = genus_of[p] != 0 && per_genus[genus_of[p]] > 1;
                    suspect = partner_has_congeners && congener[p] < 2 && foreign[i] + kCongenerMargin <= congener[p];
                }
                if (suspect) {
                    suspects.push_back({ copies_here[i].taxid, static_cast<uint32_t>(geneid), copies_here[partner[i]].taxid,
                                         foreign_rank[i], foreign[i], congener[i] });
                }
            }
            #pragma omp critical(gene_incongruence_merge)
            {
                result.suspects.insert(result.suspects.end(), suspects.begin(), suspects.end());
                result.pairs.insert(result.pairs.end(), pairs.begin(), pairs.end());
            }
        }
        result.copies = copies;
        result.genes = genes;
        result.candidates = candidates;
        std::sort(result.suspects.begin(), result.suspects.end(), [](Suspect const& a, Suspect const& b) {
            return a.geneid != b.geneid ? a.geneid < b.geneid : a.taxid < b.taxid;
        });
        auto pair_key = [](Pair const& p) { return std::tie(p.geneid, p.taxid_a, p.taxid_b); };
        std::sort(result.pairs.begin(), result.pairs.end(), [&](Pair const& a, Pair const& b) { return pair_key(a) < pair_key(b); });
        result.pairs.erase(std::unique(result.pairs.begin(), result.pairs.end(), [&](Pair const& a, Pair const& b) { return pair_key(a) == pair_key(b); }),
                           result.pairs.end());
        return result;
    }

    // The suspect copies of a database (suspect_copies.tsv): a run leaves records on them out.
    class Table {
    public:
        static uint64_t Key(uint32_t taxid, uint32_t geneid) {
            return (static_cast<uint64_t>(taxid) << kGeneBits) | geneid;
        }

        void Add(Suspect const& suspect) {
            if (m_keys.insert(Key(suspect.taxid, suspect.geneid)).second) m_rows.push_back(suspect);
        }

        bool Contains(uint32_t taxid, uint32_t geneid) const {
            return !m_keys.empty() && m_keys.count(Key(taxid, geneid)) > 0;
        }

        bool Empty() const { return m_rows.empty(); }
        size_t Size() const { return m_rows.size(); }
        std::vector<Suspect> const& Rows() const { return m_rows; }

        // The species with a suspect copy.
        size_t Species() const {
            std::unordered_set<uint32_t> species;
            for (auto const& row : m_rows) species.insert(row.taxid);
            return species.size();
        }

        // suspect_copies.tsv: a header line starting with "taxid", then per copy taxid, geneid, partner_taxid,
        // shared_rank, distance, congener_distance (tab-separated; lines starting with '#' are skipped). Returns an
        // error text, or "" when every line was read.
        std::string Read(std::istream& is) {
            std::string line;
            size_t number = 0;
            while (std::getline(is, line)) {
                number++;
                if (!line.empty() && line.back() == '\r') line.pop_back();
                if (line.empty() || line[0] == '#' || line.rfind("taxid", 0) == 0) continue;
                std::vector<std::string_view> fields;
                std::string_view view(line);
                size_t start = 0;
                for (size_t tab; (tab = view.find('\t', start)) != std::string_view::npos; start = tab + 1) fields.push_back(view.substr(start, tab - start));
                fields.push_back(view.substr(start));
                if (fields.size() < 6) return "line " + std::to_string(number) + " has " + std::to_string(fields.size()) + " fields, not 6";
                Suspect row;
                auto number_of = [&](std::string_view f, auto& value) {
                    auto const [ptr, ec] = std::from_chars(f.data(), f.data() + f.size(), value);
                    return ec == std::errc() && ptr == f.data() + f.size();
                };
                auto const rank = RankFromName(fields[3]);
                bool const no_congener = fields[5] == "none";
                if (no_congener) row.congener_distance = 2;
                if (!number_of(fields[0], row.taxid) || !number_of(fields[1], row.geneid) || !number_of(fields[2], row.partner) ||
                    !rank || !number_of(fields[4], row.distance) || (!no_congener && !number_of(fields[5], row.congener_distance))) {
                    return "line " + std::to_string(number) + " is not taxid, geneid, partner_taxid, shared_rank, distance, congener_distance: " + line;
                }
                if (row.geneid >= (uint32_t{1} << kGeneBits)) return "line " + std::to_string(number) + ": gene id " + std::to_string(row.geneid) + " too large";
                row.rank = *rank;
                Add(row);
            }
            return "";
        }

        void Write(std::ostream& os) const {
            os << "taxid\tgeneid\tpartner_taxid\tshared_rank\tdistance\tcongener_distance\n";
            for (auto const& row : m_rows) {
                os << row.taxid << '\t' << row.geneid << '\t' << row.partner << '\t' << RankName(row.rank) << '\t' << row.distance << '\t';
                if (row.congener_distance >= 2) os << "none";
                else os << row.congener_distance;
                os << '\n';
            }
        }

    private:
        std::unordered_set<uint64_t> m_keys;
        std::vector<Suspect> m_rows;
    };

    // gene_incongruence.tsv: every near pair across genera (distance at most kReportDistance), with which of the two
    // copies are suspect.
    inline void WriteReport(std::ostream& os, Result const& result) {
        std::unordered_set<uint64_t> suspect;
        for (auto const& s : result.suspects) suspect.insert(Table::Key(s.taxid, s.geneid));
        os << "geneid\ttaxid_a\ttaxid_b\tshared_rank\tdistance\tsuspect_a\tsuspect_b\n";
        for (auto const& p : result.pairs) {
            os << p.geneid << '\t' << p.taxid_a << '\t' << p.taxid_b << '\t' << RankName(p.rank) << '\t' << p.distance << '\t'
               << (suspect.count(Table::Key(p.taxid_a, p.geneid)) ? 1 : 0) << '\t' << (suspect.count(Table::Key(p.taxid_b, p.geneid)) ? 1 : 0) << '\n';
        }
    }
}
