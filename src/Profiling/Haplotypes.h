// Haplotypes.h - a long-read sample's strains of a species as rows of their own in the species' strain MSA.
//
// A sample's strain MSA row is the consensus of its reads: where two strains of the species are in the sample, an
// IUPAC code at each site where both strains' alleles pass, which qcmsa reads as a mixture. A long read covers
// several marker genes (about 15 of bac120's markers lie in the S10-spc-alpha cluster, ~14 kb), and at each of the
// sample's multi-allelic sites on them it shows the allele of the one strain it comes from. Phase() clusters the reads
// by their alleles into haplotypes, block by block (a block: the sites that reads link, directly or through other
// sites), and joins the blocks' haplotypes into rows by their shares of the reads, which are the strains' shares of
// the sample in every block: the most abundant strain is the first row in each. Each row is then called from its own
// reads, as a sample's row from the sample's (HaplotypeItem): its coverage is its strain's, and its bases are its
// strain's also where the sample's SNP filters did not see the strain's allele. A block whose haplotypes' shares do not
// tell which row each belongs to (strains of about equal abundance, or few reads) gives no row its reads.
//
// A read is cut where two of its genes next to each other are unlikely neighbours in the species' clade (the
// database's gene neighbours: a chimeric read, or a gene from elsewhere), so that its parts link nothing across.
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <map>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <vector>
#include "GeneNeighbours.h"
#include "Strain.h"
#include "VariantHandler.h"

namespace protal::haplotypes {
    // A long read's record on one gene of a taxon, as the taxon keeps it for phasing (Taxon::AddSam).
    struct ReadRecord {
        uint32_t link = 0;        // the read (all its records share it)
        uint32_t gene = 0;
        uint32_t read_start = 0;  // where on the read the record lies, in the read's own orientation
        uint32_t read_end = 0;
        bool forward = true;      // on the gene's strand
        uint8_t divergence = 0;   // DivergenceBin of the alignment
        ReadAlleles alleles;
    };

    // A multi-allelic site of a sample's gene (MultiAllelicSites in Strain.h): the bases that pass there, the most
    // observed first.
    struct Site {
        uint32_t gene = 0;
        uint32_t pos = 0;
        char reference = 'N';
        std::vector<char> alleles;
    };

    struct Settings {
        uint32_t min_reads = 3;       // reads a haplotype needs
        double min_share = 0.05;      // and this share of its block's reads
        size_t max_haplotypes = 4;    // per block, and rows
        double min_log_odds = 3.0;    // ln 20: a block's best join of its haplotypes to the rows against the next best
        double min_consensus = 0.75;  // a haplotype's base at a site: this share of its reads there, else none
        // Genes the phased blocks need between them: strains differ all along their genomes, while a second allele on
        // one or two genes only is a gene from elsewhere (another genome's copy, a foreign gene the run kept).
        size_t min_genes = 3;
    };

    // A block of a sample: its sites, its reads (or parts of reads), its haplotypes and how they join the rows.
    struct Block {
        std::vector<uint32_t> genes;
        size_t sites = 0;
        size_t reads = 0;
        std::vector<size_t> haplotype_reads;  // reads of each haplotype, the most first
        std::vector<int> row_haplotype;       // each row's haplotype; empty unless phased
        double log_odds = 0;                  // the best join against the next best
        bool phased = false;
    };

    // The phasing of one sample's species.
    struct Phasing {
        size_t rows = 0;             // haplotype rows; 0: not phased (one row, the sample's calls)
        std::vector<double> shares;  // each row's share of the reads, the largest first
        std::vector<Block> blocks;
        size_t reads = 0;            // reads (or parts) with an allele at a site
        size_t cut_reads = 0;        // reads cut between unlikely neighbours
        // Each row's reads: the records (indices into the records phased) of the reads of its haplotype in each
        // phased block, from which its MSA row is called (HaplotypeItem).
        std::vector<std::vector<uint32_t>> row_records;
        // gene -> position -> each row's haplotype's base there (0: none), at the phased sites: what the rows' reads
        // were told apart by.
        std::unordered_map<uint32_t, std::unordered_map<uint32_t, std::vector<char>>> calls;

        // The base of row `row`'s haplotype at position `pos` of gene `gene`, 0 if none (no phased site there).
        char Call(uint32_t gene, uint32_t pos, size_t row) const {
            auto const g = calls.find(gene);
            if (g == calls.end()) return 0;
            auto const p = g->second.find(pos);
            return p == g->second.end() || row >= p->second.size() ? 0 : p->second[row];
        }

        // The phased sites of gene `gene` where row `row`'s haplotype has a base.
        size_t Calls(uint32_t gene, size_t row) const {
            auto const g = calls.find(gene);
            if (g == calls.end()) return 0;
            size_t n = 0;
            for (auto const& [_, bases] : g->second) n += row < bases.size() && bases[row] != 0;
            return n;
        }

        size_t PhasedBlocks() const {
            return static_cast<size_t>(std::count_if(blocks.begin(), blocks.end(), [](Block const& b) { return b.phased; }));
        }
    };

    // Whether a read is cut between its records a and b (b after a on the read, a's taxon's): two genes next to each
    // other on it, at most the table's max_gap apart, whose facing ends are unlikely neighbours in the taxon's clade.
    using Breaks = std::function<bool(ReadRecord const&, ReadRecord const&)>;

    inline Breaks UnlikelyNeighbours(gene_neighbours::Table const& table, uint32_t taxid) {
        if (table.Empty()) return {};
        uint32_t const max_gap = table.MaxGap() > 0 ? static_cast<uint32_t>(table.MaxGap()) : 3000;
        return [&table, taxid, max_gap](ReadRecord const& a, ReadRecord const& b) {
            if (a.gene == b.gene || b.read_start > a.read_end + max_gap) return false;
            auto const verdict = table.Assess(taxid, a.gene, gene_neighbours::EndAhead(a.forward), b.gene,
                                              gene_neighbours::EndAhead(!b.forward)).verdict;
            return verdict == gene_neighbours::Verdict::Unlikely;
        };
    }

    namespace detail {
        // A read (or part of one) by its alleles at the sites: (site, allele index), ascending sites.
        using Unit = std::vector<std::pair<uint32_t, uint8_t>>;

        struct UnionFind {
            std::vector<uint32_t> parent;
            explicit UnionFind(size_t n) : parent(n) { std::iota(parent.begin(), parent.end(), 0u); }
            uint32_t Find(uint32_t x) {
                while (parent[x] != x) x = parent[x] = parent[parent[x]];
                return x;
            }
            void Unite(uint32_t a, uint32_t b) {
                a = Find(a);
                b = Find(b);
                if (a != b) parent[std::max(a, b)] = std::min(a, b);
            }
        };

        // The reads of one block clustered into haplotypes (Cluster).
        class Clustering {
        public:
            Clustering(std::vector<Unit const*> units, size_t sites, Settings const& settings) :
                    m_units(std::move(units)), m_sites(sites), m_settings(settings), m_assign(m_units.size(), -1) {}

            // Seeds a haplotype with the read of the most sites (the one most like the block's majority alleles of
            // those), and lets it grow by the reads that share sites with it and agree there (at most kJoinConflict
            // of them otherwise), better than with any other haplotype, pass after pass, so that a haplotype extends
            // along the block through the reads that overlap. When no read joins, the read of the most sites that
            // conflicts with every haplotype it shares sites with (at two sites or more, one in a block of one site,
            // and too many to join it) seeds the next, up to max_haplotypes. Reads that conflict a little then
            // join the haplotype they fit best; a haplotype of fewer than min_reads reads, or less than min_share of
            // the block's, gives its reads to the others, and two that differ at fewer than two sites are one. Reads
            // that fit no haplotype better than another (they show only alleles two of them share) are left out.
            // A haplotype seeded from a read of a strain thus grows along that strain only: clusters made by
            // splitting reads at a site and moving them would mix two strains' halves of a block longer than a read.
            void Cluster() {
                if (m_units.empty()) return;
                m_counts.clear();
                m_k = 0;
                std::vector<std::array<uint32_t, 4>> global(m_sites, std::array<uint32_t, 4>{});
                for (auto const* u : m_units) {
                    for (auto const& [s, a] : *u) global[s][a]++;
                }
                auto majority_agreement = [&](Unit const& u) {
                    size_t n = 0;
                    for (auto const& [s, a] : u) n += global[s][a] == *std::max_element(global[s].begin(), global[s].end());
                    return n;
                };
                // The best seed among the units that `ok`: the most sites, then the most like the majority.
                auto seed = [&](auto&& ok) {
                    int best = -1;
                    std::pair<size_t, size_t> best_key{ 0, 0 };
                    for (size_t u = 0; u < m_units.size(); u++) {
                        if (m_assign[u] >= 0 || !ok(u)) continue;
                        std::pair<size_t, size_t> const key{ m_units[u]->size(), majority_agreement(*m_units[u]) };
                        if (best < 0 || key > best_key) {
                            best = static_cast<int>(u);
                            best_key = key;
                        }
                    }
                    return best;
                };
                uint32_t const needed = m_sites == 1 ? 1 : 2;
                int first = seed([](size_t) { return true; });
                if (first < 0) return;
                AddHaplotype(static_cast<size_t>(first));
                while (true) {
                    Grow();
                    if (m_k >= m_settings.max_haplotypes) break;
                    int const next = seed([&](size_t u) {
                        bool shares = false;
                        for (size_t h = 0; h < m_k; h++) {
                            auto const m = MatchOf(*m_units[u], static_cast<int>(h));
                            if (m.shared == 0) continue;
                            shares = true;
                            uint32_t const conflicts = m.shared - m.agree;
                            if (conflicts < needed || conflicts <= kJoinConflict * m.shared) return false;
                        }
                        return shares;
                    });
                    if (next < 0) break;
                    AddHaplotype(static_cast<size_t>(next));
                }
                Absorb();
                // Too small, or too like another: its reads go to the others.
                for (bool changed = true; changed && m_k > 1;) {
                    changed = false;
                    size_t assigned = 0;
                    std::vector<size_t> sizes(m_k, 0);
                    for (auto a : m_assign) {
                        if (a >= 0) {
                            sizes[static_cast<size_t>(a)]++;
                            assigned++;
                        }
                    }
                    int drop = -1;
                    for (size_t h = 0; h < m_k; h++) {
                        if (sizes[h] < m_settings.min_reads || sizes[h] < m_settings.min_share * assigned) {
                            if (drop < 0 || sizes[h] < sizes[static_cast<size_t>(drop)]) drop = static_cast<int>(h);
                        }
                    }
                    for (size_t a = 0; drop < 0 && a < m_k; a++) {
                        for (size_t b = a + 1; drop < 0 && b < m_k; b++) {
                            if (Differences(static_cast<int>(a), static_cast<int>(b)) < needed) {
                                drop = static_cast<int>(sizes[a] < sizes[b] ? a : b);
                            }
                        }
                    }
                    if (drop < 0) break;
                    for (size_t u = 0; u < m_units.size(); u++) {
                        if (m_assign[u] == drop) Move(u, -1);
                    }
                    Compact();
                    Absorb();
                    changed = true;
                }
                // A few rounds of every read to the haplotype it fits strictly best.
                for (int round = 0; round < 3; round++) {
                    bool moved = false;
                    for (size_t u = 0; u < m_units.size(); u++) {
                        if (m_assign[u] < 0) continue;
                        int const best = Best(*m_units[u], m_assign[u]);
                        if (best >= 0 && best != m_assign[u]) {
                            Move(u, best);
                            moved = true;
                        }
                    }
                    if (!moved) break;
                }
                Compact();
            }

            size_t Haplotypes() const { return m_k; }

            // Unit u's haplotype (as Size and Consensus number them), -1 for none.
            int Of(size_t u) const { return m_assign[u]; }

            // The haplotypes, the most reads first.
            std::vector<size_t> Order() const {
                std::vector<size_t> order(m_k);
                std::iota(order.begin(), order.end(), 0);
                std::stable_sort(order.begin(), order.end(), [this](size_t a, size_t b) { return Size(a) > Size(b); });
                return order;
            }

            size_t Size(size_t c) const {
                return static_cast<size_t>(std::count(m_assign.begin(), m_assign.end(), static_cast<int>(c)));
            }

            // Haplotype c's allele at site s: the one min_consensus of its reads there show, else -1.
            int Consensus(size_t c, uint32_t s) const {
                auto const& counts = m_counts[Index(static_cast<int>(c), s)];
                uint32_t const total = counts[0] + counts[1] + counts[2] + counts[3];
                if (total == 0) return -1;
                auto const top = std::max_element(counts.begin(), counts.end());
                if (*top < m_settings.min_consensus * total) return -1;
                return static_cast<int>(top - counts.begin());
            }

        private:
            // A read with this much of its shared sites against a haplotype's alleles does not join it while it grows.
            static constexpr double kJoinConflict = 0.2;

            struct Match {
                uint32_t shared = 0;  // sites of the read where the haplotype has a majority allele
                uint32_t agree = 0;   // of those, where the read shows it
            };

            size_t Index(int c, uint32_t s) const { return static_cast<size_t>(c) * m_sites + s; }

            // The plain majority of haplotype c at site s (-1: none or a tie), which a read is matched against.
            int Majority(int c, uint32_t s) const {
                auto const& counts = m_counts[Index(c, s)];
                auto const top = std::max_element(counts.begin(), counts.end());
                if (*top == 0 || std::count(counts.begin(), counts.end(), *top) > 1) return -1;
                return static_cast<int>(top - counts.begin());
            }

            Match MatchOf(Unit const& unit, int c) const {
                Match m;
                for (auto const& [s, a] : unit) {
                    int const major = Majority(c, s);
                    if (major < 0) continue;
                    m.shared++;
                    m.agree += major == a;
                }
                return m;
            }

            // The haplotype a read fits strictly best (agreeing sites less conflicting ones, with a shared site), or
            // `keep` if that one ties for best; -1 if none.
            int Best(Unit const& unit, int keep = -1) const {
                int best = -1, best_score = 0, keep_score = INT32_MIN;
                bool tie = false;
                for (int c = 0; c < static_cast<int>(m_k); c++) {
                    auto const m = MatchOf(unit, c);
                    if (m.shared == 0) continue;
                    int const score = 2 * static_cast<int>(m.agree) - static_cast<int>(m.shared);
                    if (c == keep) keep_score = score;
                    if (best < 0 || score > best_score) {
                        best = c;
                        best_score = score;
                        tie = false;
                    } else if (score == best_score) {
                        tie = true;
                    }
                }
                if (!tie) return best;
                return keep >= 0 && keep_score == best_score ? keep : -1;
            }

            void Move(size_t u, int to) {
                if (m_assign[u] >= 0) {
                    for (auto const& [s, a] : *m_units[u]) m_counts[Index(m_assign[u], s)][a]--;
                }
                m_assign[u] = to;
                if (to >= 0) {
                    for (auto const& [s, a] : *m_units[u]) m_counts[Index(to, s)][a]++;
                }
            }

            void AddHaplotype(size_t seed) {
                m_counts.resize((m_k + 1) * m_sites, std::array<uint32_t, 4>{});
                Move(seed, static_cast<int>(m_k));
                m_k++;
            }

            // Passes over the reads left: each joins the haplotype it fits strictly best if it agrees with it at all
            // but kJoinConflict of the sites they share; until a pass adds none.
            void Grow() {
                for (bool grew = true; grew;) {
                    grew = false;
                    for (size_t u = 0; u < m_units.size(); u++) {
                        if (m_assign[u] >= 0) continue;
                        int const best = Best(*m_units[u]);
                        if (best < 0) continue;
                        auto const m = MatchOf(*m_units[u], best);
                        if (m.shared - m.agree > kJoinConflict * m.shared) continue;
                        Move(u, best);
                        grew = true;
                    }
                }
            }

            // The reads left join the haplotype they fit strictly best, conflicts or not.
            void Absorb() {
                for (size_t u = 0; u < m_units.size(); u++) {
                    if (m_assign[u] >= 0) continue;
                    int const best = Best(*m_units[u]);
                    if (best >= 0) Move(u, best);
                }
            }

            // Renumbers the haplotypes without the empty ones.
            void Compact() {
                std::vector<int> map(m_k, -1);
                int next = 0;
                for (size_t c = 0; c < m_k; c++) {
                    if (std::find(m_assign.begin(), m_assign.end(), static_cast<int>(c)) != m_assign.end()) map[c] = next++;
                }
                std::vector<std::array<uint32_t, 4>> counts(static_cast<size_t>(next) * m_sites, std::array<uint32_t, 4>{});
                for (size_t c = 0; c < m_k; c++) {
                    if (map[c] < 0) continue;
                    std::copy_n(m_counts.begin() + static_cast<std::ptrdiff_t>(c * m_sites), m_sites,
                                counts.begin() + static_cast<std::ptrdiff_t>(static_cast<size_t>(map[c]) * m_sites));
                }
                for (auto& a : m_assign) {
                    if (a >= 0) a = map[static_cast<size_t>(a)];
                }
                m_counts = std::move(counts);
                m_k = static_cast<size_t>(next);
            }

            // Sites where haplotypes a and b both have an allele (Consensus) and they differ.
            uint32_t Differences(int a, int b) const {
                uint32_t n = 0;
                for (uint32_t s = 0; s < m_sites; s++) {
                    int const x = Consensus(static_cast<size_t>(a), s), y = Consensus(static_cast<size_t>(b), s);
                    n += x >= 0 && y >= 0 && x != y;
                }
                return n;
            }

            std::vector<Unit const*> m_units;
            size_t m_sites;
            Settings m_settings;
            std::vector<int> m_assign;  // -1: in no haplotype
            size_t m_k = 0;
            std::vector<std::array<uint32_t, 4>> m_counts;  // (haplotype, site) -> reads of each allele
        };

        // Every map of `rows` rows onto `haplotypes` haplotypes that gives each haplotype a row (rows >= haplotypes).
        inline void Joins(size_t rows, size_t haplotypes, std::vector<std::vector<int>>& out) {
            out.clear();
            std::vector<int> join(rows, 0);
            while (true) {
                std::vector<bool> used(haplotypes, false);
                for (int h : join) used[static_cast<size_t>(h)] = true;
                if (std::all_of(used.begin(), used.end(), [](bool u) { return u; })) out.push_back(join);
                size_t i = 0;
                while (i < rows && ++join[i] == static_cast<int>(haplotypes)) join[i++] = 0;
                if (i == rows) break;
            }
        }
    }

    // The haplotype rows of one sample's species from its long reads' records (`records`, the taxon's own; those of
    // more than max_divergence are left out, as from the strain MSA) at the sample's multi-allelic sites (`sites`).
    // `breaks` (may be empty) says where a read is cut (UnlikelyNeighbours).
    //
    // Per block, the reads are clustered (detail::Clustering). The rows are as many as the blocks most often have
    // haplotypes (each block weighted by its sites), their shares the blocks' with that many haplotypes, pooled by
    // rank. The anchor, the block of that many haplotypes with the most sites, gives the rows its haplotypes by rank:
    // within a block the reads tell the strains apart, whatever their shares. Any other block's haplotypes join the
    // rows by the multinomial likelihood of its haplotypes' reads: each row gets one haplotype (and a block of fewer
    // haplotypes than rows gives one to several rows: strains alike there), its most likely join if that is
    // min_log_odds more likely than the next; else the block is not phased. A row's reads are those of its haplotype
    // in each phased block (row_records: all their records, on genes with sites or not). Not phased at all (rows 0)
    // with fewer than two rows, or with phased blocks on fewer than min_genes genes.
    inline Phasing Phase(std::vector<ReadRecord> const& records, std::vector<Site> const& sites, uint8_t max_divergence,
                         Breaks const& breaks = {}, Settings const& settings = {}) {
        Phasing result;
        if (sites.empty() || records.empty()) return result;

        // The sites of each gene, by position.
        std::unordered_map<uint32_t, std::vector<std::pair<uint32_t, uint32_t>>> of_gene;  // gene -> (pos, site)
        for (uint32_t s = 0; s < sites.size(); s++) of_gene[sites[s].gene].emplace_back(sites[s].pos, s);
        for (auto& [_, list] : of_gene) std::sort(list.begin(), list.end());

        // Each read's records in read order, cut where `breaks` says; each part a unit of its alleles at the sites, with
        // its records (indices into `records`).
        std::vector<uint32_t> order;
        order.reserve(records.size());
        for (uint32_t i = 0; i < records.size(); i++) {
            if (records[i].divergence <= max_divergence) order.push_back(i);
        }
        std::stable_sort(order.begin(), order.end(), [&records](uint32_t a, uint32_t b) {
            return records[a].link != records[b].link ? records[a].link < records[b].link : records[a].read_start < records[b].read_start;
        });
        std::vector<detail::Unit> units;
        std::vector<std::vector<uint32_t>> unit_records;
        detail::Unit unit;
        std::vector<uint32_t> in_unit;
        auto finish = [&]() {
            if (!unit.empty()) {
                std::sort(unit.begin(), unit.end());
                unit.erase(std::unique(unit.begin(), unit.end(), [](auto const& a, auto const& b) { return a.first == b.first; }), unit.end());
                units.push_back(std::move(unit));
                unit_records.push_back(std::move(in_unit));
            }
            unit.clear();
            in_unit.clear();
        };
        for (size_t i = 0; i < order.size(); i++) {
            auto const& r = records[order[i]];
            if (i > 0 && records[order[i - 1]].link != r.link) {
                finish();
            } else if (i > 0 && breaks && breaks(records[order[i - 1]], r)) {
                finish();
                result.cut_reads++;
            }
            in_unit.push_back(order[i]);
            auto const g = of_gene.find(r.gene);
            if (g == of_gene.end()) continue;
            auto it = std::lower_bound(g->second.begin(), g->second.end(), std::pair<uint32_t, uint32_t>{ r.alleles.start, 0 });
            for (; it != g->second.end() && it->first < r.alleles.end; ++it) {
                auto const& site = sites[it->second];
                char const base = r.alleles.BaseAt(it->first, site.reference);
                auto const a = std::find(site.alleles.begin(), site.alleles.end(), base);
                if (base == 0 || a == site.alleles.end() || a - site.alleles.begin() >= 4) continue;
                unit.emplace_back(it->second, static_cast<uint8_t>(a - site.alleles.begin()));
            }
        }
        finish();
        result.reads = units.size();

        // Blocks: the sites linked by the units.
        detail::UnionFind links(sites.size());
        for (auto const& u : units) {
            for (size_t i = 1; i < u.size(); i++) links.Unite(u[0].first, u[i].first);
        }
        std::map<uint32_t, std::vector<uint32_t>> block_sites;  // root -> sites, by the block's first site
        for (uint32_t s = 0; s < sites.size(); s++) block_sites[links.Find(s)].push_back(s);
        std::map<uint32_t, std::vector<size_t>> block_units;
        for (size_t u = 0; u < units.size(); u++) block_units[links.Find(units[u][0].first)].push_back(u);

        // Each block's haplotypes: their reads, the most first, their alleles at the block's sites, and each unit's.
        struct Clustered {
            std::vector<uint32_t> sites;
            std::vector<size_t> reads;               // per haplotype, the most first
            std::vector<std::vector<int>> alleles;   // per haplotype, per site of the block: allele index or -1
            std::vector<size_t> units;               // the block's units
            std::vector<int> haplotype;              // each unit's haplotype (as `reads`), -1 for none
        };
        std::vector<Clustered> clustered;
        for (auto const& [root, members] : block_units) {
            auto const& in_block = block_sites.at(root);
            std::unordered_map<uint32_t, uint32_t> local;
            for (uint32_t i = 0; i < in_block.size(); i++) local[in_block[i]] = i;
            std::vector<detail::Unit> renumbered;
            renumbered.reserve(members.size());
            for (auto const u : members) {
                detail::Unit r;
                for (auto const& [s, a] : units[u]) r.emplace_back(local.at(s), a);
                renumbered.push_back(std::move(r));
            }
            std::vector<detail::Unit const*> pointers;
            for (auto const& r : renumbered) pointers.push_back(&r);
            detail::Clustering clustering(std::move(pointers), in_block.size(), settings);
            clustering.Cluster();
            Clustered c;
            c.sites = in_block;
            c.units = members;
            auto const ranked = clustering.Order();
            std::vector<int> rank_of(ranked.size(), -1);
            for (size_t k = 0; k < ranked.size(); k++) {
                rank_of[ranked[k]] = static_cast<int>(k);
                c.reads.push_back(clustering.Size(ranked[k]));
                std::vector<int> alleles(in_block.size());
                for (uint32_t s = 0; s < in_block.size(); s++) alleles[s] = clustering.Consensus(ranked[k], s);
                c.alleles.push_back(std::move(alleles));
            }
            for (size_t j = 0; j < members.size(); j++) {
                int const h = clustering.Of(j);
                c.haplotype.push_back(h < 0 ? -1 : rank_of[static_cast<size_t>(h)]);
            }
            clustered.push_back(std::move(c));
        }

        // The rows: as many as the blocks most often have haplotypes, weighted by their sites.
        std::map<size_t, size_t> weight;
        for (auto const& c : clustered) weight[c.reads.size()] += c.sites.size();
        size_t rows = 0, best_weight = 0;
        for (auto const& [k, w] : weight) {
            if (w > best_weight) {
                rows = k;
                best_weight = w;
            }
        }
        rows = std::min(rows, settings.max_haplotypes);
        std::vector<double> pooled(rows, 0);
        for (auto const& c : clustered) {
            if (c.reads.size() != rows) continue;
            for (size_t r = 0; r < rows; r++) pooled[r] += static_cast<double>(c.reads[r]);
        }
        double const total = std::accumulate(pooled.begin(), pooled.end(), 0.0);

        // The anchor: the block with as many haplotypes as rows and the most sites (then reads). Its haplotypes are the
        // rows, by rank; the others join them.
        size_t anchor = clustered.size();
        for (size_t b = 0; b < clustered.size(); b++) {
            auto const& c = clustered[b];
            if (rows < 2 || c.reads.size() != rows) continue;
            auto key = [&clustered](size_t i) {
                auto const& x = clustered[i];
                return std::pair(x.sites.size(), std::accumulate(x.reads.begin(), x.reads.end(), size_t{0}));
            };
            if (anchor == clustered.size() || key(b) > key(anchor)) anchor = b;
        }

        std::vector<std::vector<int>> joins;
        result.row_records.assign(rows >= 2 ? rows : 0, {});
        for (size_t b = 0; b < clustered.size(); b++) {
            auto const& c = clustered[b];
            Block block;
            block.sites = c.sites.size();
            block.reads = std::accumulate(c.reads.begin(), c.reads.end(), size_t{0});
            block.haplotype_reads = c.reads;
            for (auto s : c.sites) {
                if (std::find(block.genes.begin(), block.genes.end(), sites[s].gene) == block.genes.end()) block.genes.push_back(sites[s].gene);
            }
            size_t const haplotypes = std::min(c.reads.size(), rows);
            if (rows >= 2 && haplotypes >= 2 && total > 0) {
                detail::Joins(rows, haplotypes, joins);
                double best = -INFINITY, next = -INFINITY;
                std::vector<int> const* best_join = nullptr;
                for (auto const& join : joins) {
                    std::vector<double> q(haplotypes, 0);
                    for (size_t r = 0; r < rows; r++) q[static_cast<size_t>(join[r])] += pooled[r] / total;
                    double ll = 0;
                    for (size_t h = 0; h < haplotypes; h++) ll += static_cast<double>(c.reads[h]) * std::log(std::max(q[h], 1e-12));
                    if (ll > best) {
                        next = best;
                        best = ll;
                        best_join = &join;
                    } else if (ll > next) {
                        next = ll;
                    }
                }
                block.log_odds = std::isfinite(next) ? best - next : INFINITY;
                std::vector<int> by_rank(rows);
                std::iota(by_rank.begin(), by_rank.end(), 0);
                if (b == anchor || (best_join && block.log_odds >= settings.min_log_odds)) {
                    block.phased = true;
                    block.row_haplotype = b == anchor ? by_rank : *best_join;
                    for (uint32_t i = 0; i < c.sites.size(); i++) {
                        auto const& site = sites[c.sites[i]];
                        std::vector<char> bases(rows, 0);
                        for (size_t r = 0; r < rows; r++) {
                            int const allele = c.alleles[static_cast<size_t>(block.row_haplotype[r])][i];
                            if (allele >= 0) bases[r] = site.alleles[static_cast<size_t>(allele)];
                        }
                        if (std::any_of(bases.begin(), bases.end(), [](char b) { return b != 0; })) {
                            result.calls[site.gene][site.pos] = std::move(bases);
                        }
                    }
                    for (size_t r = 0; r < rows; r++) {
                        for (size_t j = 0; j < c.units.size(); j++) {
                            if (c.haplotype[j] != block.row_haplotype[r]) continue;
                            auto const& of_unit = unit_records[c.units[j]];
                            result.row_records[r].insert(result.row_records[r].end(), of_unit.begin(), of_unit.end());
                        }
                    }
                }
            }
            result.blocks.push_back(std::move(block));
        }
        std::vector<uint32_t> phased_genes;
        for (auto const& block : result.blocks) {
            if (block.phased) phased_genes.insert(phased_genes.end(), block.genes.begin(), block.genes.end());
        }
        std::sort(phased_genes.begin(), phased_genes.end());
        phased_genes.erase(std::unique(phased_genes.begin(), phased_genes.end()), phased_genes.end());
        if (rows >= 2 && phased_genes.size() < settings.min_genes) {
            for (auto& block : result.blocks) {
                block.phased = false;
                block.row_haplotype.clear();
            }
        }
        if (rows >= 2 && result.PhasedBlocks() > 0) {
            result.rows = rows;
            for (auto p : pooled) result.shares.push_back(total > 0 ? p / total : 0);
            for (auto& r : result.row_records) std::sort(r.begin(), r.end());
        } else {
            result.calls.clear();
            result.row_records.clear();
        }
        return result;
    }

    // A strain row's gene as the strain MSA takes a sample's (StrainLevelContainer::MSAItem): from the records of the
    // row's reads on the gene (Phasing::row_records), the reads with a base at each position per strand, and their SNPs
    // by their base qualities, each site completed and filtered as a sample's (VariantHandler::PostProcessSNPBin, with
    // the SNP filters of the sample's reads). Insertions are not kept: the row has none, and a deletion is a gap.
    inline std::pair<VariantVec, CoverageVec> HaplotypeItem(std::vector<ReadRecord const*> const& gene_records,
                                                            std::string_view reference, size_t min_observations,
                                                            double min_frequency, size_t min_avg_quality,
                                                            size_t min_phred_sum, bool require_strand) {
        size_t const n = reference.size();
        CoverageVec forward(n, 0), reverse(n, 0);
        std::map<uint32_t, VariantBin> bins;
        for (auto const* r : gene_records) {
            auto& coverage = r->forward ? forward : reverse;
            auto const& a = r->alleles;
            size_t k = 0;
            for (uint32_t pos = a.start; pos < a.end && pos < n; pos++) {
                while (k < a.no_base.size() && a.no_base[k] < pos) k++;
                if (k < a.no_base.size() && a.no_base[k] == pos) continue;
                coverage[pos]++;
            }
            for (size_t s = 0; s < a.snps.size(); s++) {
                auto const [pos, base] = a.snps[s];
                if (pos >= n) continue;
                Qual const quality = s < a.snp_quals.size() ? a.snp_quals[s] : 0;
                VariantHandler::GetVariant(bins[pos], pos, base, reference[pos]).AddObservation(quality, r->forward, r->divergence);
            }
        }
        VariantVec out;
        out.reserve(bins.size());
        for (auto& [pos, bin] : bins) {
            VariantHandler::PostProcessSNPBin(bin, forward[pos], reverse[pos], min_observations, min_observations, min_frequency,
                                              min_avg_quality, min_phred_sum, require_strand);
            out.push_back(std::move(bin));
        }
        CoverageVec informative(n, 0);
        for (size_t pos = 0; pos < n; pos++) informative[pos] = forward[pos] + reverse[pos];
        return { std::move(out), std::move(informative) };
    }
}
