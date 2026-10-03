// GeneNeighbours.h - how often each marker gene end faces which other in the genomes of a clade, and how far
// apart: gene_neighbours.tsv, a gene-by-gene frequency table per clade (family, order, class, phylum,
// domain), counted by scripts/mini_db/gene_neighbours.py from every genome at hand when the database is built
// (the representatives and other strains downloaded to simulate training data: genes placed by their exact
// sequence or k-mer trace; complete replicons, and representatives in one sequence, read as circular, across
// the origin), and packed into the database by --build. The database also holds where each gene lies in each
// of those genomes (kPositionsFileName), which a run does not load: the frequencies are derived from it, and a
// copy without some species derives them anew. A database without them works as before.
//
// Gene order is conserved across bacteria (the ribosomal protein operons, rpoBC, ...), so this helps to
// follow a fragment or a long read from one gene into the next; it does not tell congeners apart. A run
// uses it to look for a mate past the end of its guiding mate's gene, on the gene there (MateGuidance.h),
// to pair mates on two neighbouring genes (classify::JoinAlignmentPairs), to look on a long read for the
// genes its genes are next to (LongReads.h), and for the profiler's adjacent_expected_share,
// adjacent_unlikely_share and adjacent_support (Profiler.h). --no_gene_neighbours turns all of it off.
//
// Each gene has two ends in its coding orientation: 5' (position 0) and 3' (its last base). A line
//   clade gene end partner partner_end species informative gap_median gap_min gap_max
// says that in `species` of the clade's `informative` species (those whose genomes go on for the build's
// max_gap past that end), the end faces partner's partner_end gap bases away (negative: they overlap): the
// observed frequency species / informative (Rule::Share); partner 0 (partner_end 0): no marker within
// max_gap. For a species, a pairing's share is that of its nearest clade with data on the gene's end, each
// clade's share smoothed towards its parent's: a clade with few species leans on the clades above it, the
// fewer the more (Table::Assess, Table::Partners).
#pragma once

#include <algorithm>
#include <atomic>
#include <charconv>
#include <cstdint>
#include <istream>
#include <map>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace protal::gene_neighbours {
    inline const std::string kFileName = "gene_neighbours.tsv";
    // Where each gene lies in each genome the frequencies were counted from (gene_neighbours.py), kept in the
    // database (checked and packed by --build) but not read by a run.
    inline const std::string kPositionsFileName = "gene_positions.tsv";

    // A gene's end in its coding orientation.
    enum class End : uint8_t { Five = 5, Three = 3 };

    // The end of a gene a read in orientation `forward` (to the gene's coding strand) points to, i.e. runs
    // towards from its start.
    inline End EndAhead(bool forward) { return forward ? End::Three : End::Five; }

    // Whether two genes whose ends a and b face each other lie on the same strand: the 3' end of one faces
    // the 5' end of the other.
    inline bool SameStrand(End a, End b) { return a != b; }

    // The orientation on gene b of a read in orientation forward_on_a on gene a, the two facing each other
    // with ends a and b.
    inline bool OrientationOnPartner(bool forward_on_a, End a, End b) {
        return SameStrand(a, b) ? forward_on_a : !forward_on_a;
    }

    // The bases of a read from where it starts on a gene to the gene's end `end` it runs towards, on an
    // alignment at gene position start, length bases long, on a gene of gene_length: the stretch of a
    // fragment that lies on this gene.
    inline int64_t ReachToEnd(int64_t start, int64_t length, int64_t gene_length, End end) {
        return end == End::Three ? gene_length - start : start + length;
    }

    // The gene positions [first, last) of partner gene b (partner_length bases) that lie from `from` to `to`
    // bases past the end of gene a facing it (0: the first base past it), with the genes gap_min to gap_max
    // bases apart; first >= last if none.
    inline std::pair<int64_t, int64_t> PartnerStretch(int64_t from, int64_t to, int64_t gap_min, int64_t gap_max,
                                                      End b, int64_t partner_length) {
        int64_t const near = std::max<int64_t>(0, from - gap_max);   // from b's facing end
        int64_t const far = std::min(partner_length, to - gap_min);
        if (near >= far) return { 0, 0 };
        return b == End::Five ? std::pair<int64_t, int64_t>{ near, far }
                              : std::pair<int64_t, int64_t>{ partner_length - far, partner_length - near };
    }

    // What one clade's genomes show at one end of one gene: one partner and in how many of its species.
    struct Rule {
        uint32_t partner = 0;      // gene id; 0: no marker within max_gap
        End partner_end = End::Five;
        uint32_t species = 0;      // species with this partner at this end
        uint32_t informative = 0;  // species in which this end was informative
        int32_t gap_median = 0;
        int32_t gap_min = 0;
        int32_t gap_max = 0;

        bool HasPartner() const { return partner != 0; }
        // The share of the clade's species informative at this end that have this partner there.
        double Share() const { return informative == 0 ? 0 : static_cast<double>(species) / informative; }
    };

    // A pairing's share in a species' clades, from the top down: the top clade with data on the gene's end has
    // its own share; each one below has (species + kPriorSpecies * its parent's share) / (informative +
    // kPriorSpecies), so that a clade with the end informative in kPriorSpecies species counts as much as the
    // clades above it, one with fewer less, one with many more. The nearest clade's share says: expected from
    // kExpectedShare, unlikely up to kUnlikelyShare (0: never seen), rare in between. If even the top clade has
    // the end informative in fewer than kMinInformative species, too few to call a pairing they lack unlikely,
    // a pairing seen in one of them is expected and another unknown.
    inline constexpr double kPriorSpecies = 3;
    inline constexpr uint32_t kMinInformative = 5;
    inline constexpr double kExpectedShare = 0.2;
    inline constexpr double kUnlikelyShare = 0.05;

    enum class Verdict : uint8_t {
        Unknown,   // no clade of the species has data on this end, or too little
        Expected,  // the two ends face each other in kExpectedShare of the species' clades or more (smoothed)
        Unlikely,  // in kUnlikelyShare or fewer (none: never seen)
        Rare,      // in between
    };

    struct Assessment {
        Verdict verdict = Verdict::Unknown;
        Rule const* rule = nullptr;  // the nearest clade's rule for the pairing among those that saw it (with Expected)
        uint32_t clade = 0;          // the nearest clade with data on the end
        double share = 0;            // the pairing's smoothed share there
    };

    // A partner of one end of a gene in a species (Table::Partners).
    struct Partner {
        Rule const* rule = nullptr;  // that of the nearest of the species' clades that saw it
        double share = 0;            // smoothed, as Assessment::share
        Verdict verdict = Verdict::Unknown;
    };
    // Bases of leeway on a clade's gap range, for a species whose genes are a little further or closer.
    inline constexpr int64_t kGapSlack = 100;

    class Table {
    public:
        bool Empty() const { return m_rules.empty(); }
        size_t Rules() const { return m_rules.size(); }
        // Differs between tables and changes with every Read and SetLineage: what a cache of a table's answers
        // (across_genes::Neighbours) checks before it trusts them.
        uint64_t Generation() const { return m_generation; }
        size_t Clades() const { return m_clades.size(); }
        size_t Genomes() const { return m_genomes; }   // from the file's comment line, 0 if it has none
        size_t MaxGap() const { return m_max_gap; }
        size_t BoundSpecies() const { return m_chain_of.size(); }

        // gene_neighbours.tsv: a header line starting with "clade", lines starting with '#' (one may say
        // "genomes=N max_gap=N") and ten tab-separated numbers per line. Returns the first problem with its
        // line ("line 3: ..."), or an empty string.
        std::string Read(std::istream& is) {
            std::string line;
            size_t number = 0;
            std::vector<std::pair<uint64_t, Rule>> rows;
            while (std::getline(is, line)) {
                number++;
                if (!line.empty() && line.back() == '\r') line.pop_back();
                if (line.empty()) continue;
                if (line[0] == '#') {
                    ReadComment(line);
                    continue;
                }
                if (line.rfind("clade", 0) == 0) continue;
                auto problem = [&](std::string const& what) { return "line " + std::to_string(number) + ": " + what; };
                int64_t v[10];
                std::string_view rest(line);
                for (size_t i = 0; i < 10; i++) {
                    size_t const tab = rest.find('\t');
                    auto const field = rest.substr(0, tab);
                    auto const r = std::from_chars(field.data(), field.data() + field.size(), v[i]);
                    if (r.ec != std::errc() || r.ptr != field.data() + field.size()) {
                        return problem("expected ten numbers separated by tabs (clade, gene, end, partner, partner_end, "
                                       "species, informative, gap_median, gap_min, gap_max)");
                    }
                    if ((tab == std::string_view::npos) != (i == 9)) {
                        return problem("expected ten fields separated by tabs");
                    }
                    if (tab != std::string_view::npos) rest.remove_prefix(tab + 1);
                }
                auto const [clade, gene, end, partner, partner_end, species, informative, median, low, high] = v;
                if (clade <= 0 || clade > UINT32_MAX || gene <= 0 || gene >= kGeneLimit || partner < 0 || partner >= kGeneLimit) {
                    return problem("clade and gene ids must be positive numbers, partner a gene id or 0");
                }
                if (end != 3 && end != 5) return problem("end must be 3 or 5");
                if ((partner == 0) != (partner_end == 0) || (partner != 0 && partner_end != 3 && partner_end != 5)) {
                    return problem("partner_end must be 3 or 5, or 0 with partner 0");
                }
                if (species <= 0 || informative < species || informative > UINT32_MAX) {
                    return problem("species must be at least 1 and at most informative");
                }
                if (low > median || median > high || low < INT32_MIN || high > INT32_MAX) {
                    return problem("expected gap_min <= gap_median <= gap_max");
                }
                Rule rule;
                rule.partner = static_cast<uint32_t>(partner);
                rule.partner_end = partner_end == 3 ? End::Three : End::Five;
                rule.species = static_cast<uint32_t>(species);
                rule.informative = static_cast<uint32_t>(informative);
                rule.gap_median = static_cast<int32_t>(median);
                rule.gap_min = static_cast<int32_t>(low);
                rule.gap_max = static_cast<int32_t>(high);
                rows.emplace_back(Key(static_cast<uint32_t>(clade), static_cast<uint32_t>(gene), end == 3 ? End::Three : End::Five), rule);
            }
            if (is.bad()) return "read error";
            std::stable_sort(rows.begin(), rows.end(), [](auto const& a, auto const& b) {
                if (a.first != b.first) return a.first < b.first;
                return std::pair(a.second.partner, a.second.partner_end) < std::pair(b.second.partner, b.second.partner_end);
            });
            for (size_t i = 1; i < rows.size(); i++) {
                auto const& [ka, a] = rows[i - 1];
                auto const& [kb, b] = rows[i];
                if (ka == kb && a.partner == b.partner && a.partner_end == b.partner_end) {
                    return "clade " + std::to_string(kb >> 32) + ", gene " + std::to_string((kb >> 1) & (kGeneLimit - 1)) +
                           ": a partner is listed twice";
                }
                if (ka == kb && a.informative != b.informative) {
                    return "clade " + std::to_string(kb >> 32) + ", gene " + std::to_string((kb >> 1) & (kGeneLimit - 1)) +
                           ": the lines of one gene end differ in informative";
                }
            }
            m_generation = NextGeneration();
            m_keys.clear();
            m_rules.clear();
            m_clades.clear();
            for (auto& [key, rule] : rows) {
                m_keys.push_back(key);
                m_rules.push_back(rule);
                uint32_t const clade = static_cast<uint32_t>(key >> 32);
                if (m_clades.empty() || m_clades.back() != clade) m_clades.push_back(clade);
            }
            m_chain_of.clear();
            m_chains.clear();
            m_chain_index.clear();
            return {};
        }

        // Which clades of the rules a species belongs to: its ancestors, nearest first (any ranks; those without
        // rules are left out). A species not set has no rules.
        void SetLineage(uint32_t species, std::vector<uint32_t> const& ancestors) {
            std::vector<uint32_t> chain;
            for (auto const clade : ancestors) {
                if (std::binary_search(m_clades.begin(), m_clades.end(), clade)) chain.push_back(clade);
            }
            if (chain.empty()) return;
            auto [it, added] = m_chain_index.try_emplace(chain, static_cast<uint32_t>(m_chains.size()));
            if (added) m_chains.push_back(std::move(chain));
            m_chain_of[species] = it->second;
            m_generation = NextGeneration();
        }

        // The partners that species taxid's clades saw at end `end` of `gene`, none (partner 0) among them, each
        // with the rule of the nearest clade that saw it and its smoothed share and verdict (Assess); empty if no
        // clade has data on that end. The nearest clade's partners first.
        void Partners(uint32_t taxid, uint32_t gene, End end, std::vector<Partner>& out) const {
            out.clear();
            auto const* chain = ChainOf(taxid);
            if (!chain) return;
            for (auto const clade : *chain) {
                for (auto const& rule : Range(clade, gene, end)) {
                    bool const known = std::any_of(out.begin(), out.end(), [&rule](Partner const& p) {
                        return p.rule->partner == rule.partner && p.rule->partner_end == rule.partner_end;
                    });
                    if (!known) out.push_back({ &rule });
                }
            }
            for (auto& p : out) {
                auto const assessment = Assess(taxid, gene, end, p.rule->partner, p.rule->partner_end);
                p.share = assessment.share;
                p.verdict = assessment.verdict;
            }
        }

        // Whether, in species taxid, end `end` of `gene` faces end partner_end of `partner`: by the pairing's
        // share smoothed over the species' clades (kPriorSpecies), at the nearest clade with data on that end
        // (Expected, Rare, Unlikely: kExpectedShare, kUnlikelyShare); if the top clade with data has the end
        // informative in fewer than kMinInformative species, Expected if a clade saw it, else Unknown. Unknown
        // if no clade has data.
        Assessment Assess(uint32_t taxid, uint32_t gene, End end, uint32_t partner, End partner_end) const {
            auto const* chain = ChainOf(taxid);
            if (!chain) return {};
            Assessment a;
            bool populated = false;
            for (auto clade = chain->rbegin(); clade != chain->rend(); ++clade) {  // from the top
                auto const range = Range(*clade, gene, end);
                if (range.empty()) continue;
                Rule const* seen = nullptr;
                for (auto const& rule : range) {
                    if (rule.partner == partner && rule.partner_end == partner_end) seen = &rule;
                }
                double const species = seen ? seen->species : 0.0;
                double const informative = range.front().informative;
                if (a.clade == 0) {
                    a.share = species / informative;
                    populated = range.front().informative >= kMinInformative;
                } else {
                    a.share = (species + kPriorSpecies * a.share) / (informative + kPriorSpecies);
                }
                a.clade = *clade;
                if (seen) a.rule = seen;
            }
            if (a.clade == 0) return {};
            if (!populated) {
                if (!a.rule) return {};
                a.verdict = Verdict::Expected;
                return a;
            }
            a.verdict = a.share >= kExpectedShare ? Verdict::Expected
                      : a.share <= kUnlikelyShare ? Verdict::Unlikely : Verdict::Rare;
            return a;
        }

        // As Assess, for two genes of a read whose ends may face either way: Expected if any of their four end
        // pairings is, Unlikely if one is and none is expected, else Unknown.
        Verdict AssessGenes(uint32_t taxid, uint32_t gene, uint32_t partner) const {
            bool unlikely = false;
            for (End a : { End::Five, End::Three }) {
                for (End b : { End::Five, End::Three }) {
                    auto const verdict = Assess(taxid, gene, a, partner, b).verdict;
                    if (verdict == Verdict::Expected) return verdict;
                    unlikely |= verdict == Verdict::Unlikely;
                }
            }
            return unlikely ? Verdict::Unlikely : Verdict::Unknown;
        }

    private:
        static constexpr int64_t kGeneLimit = int64_t{1} << 30;

        static uint64_t Key(uint32_t clade, uint32_t gene, End end) {
            return (static_cast<uint64_t>(clade) << 32) | (static_cast<uint64_t>(gene) << 1) | (end == End::Three ? 1u : 0u);
        }

        std::span<Rule const> Range(uint32_t clade, uint32_t gene, End end) const {
            auto const key = Key(clade, gene, end);
            auto const [first, last] = std::equal_range(m_keys.begin(), m_keys.end(), key);
            return { m_rules.data() + (first - m_keys.begin()), static_cast<size_t>(last - first) };
        }

        std::vector<uint32_t> const* ChainOf(uint32_t taxid) const {
            auto const it = m_chain_of.find(taxid);
            return it == m_chain_of.end() ? nullptr : &m_chains[it->second];
        }

        void ReadComment(std::string const& line) {
            auto value = [&line](std::string_view name) -> size_t {
                auto const at = line.find(name);
                if (at == std::string::npos) return 0;
                size_t v = 0;
                char const* p = line.data() + at + name.size();
                std::from_chars(p, line.data() + line.size(), v);
                return v;
            };
            if (auto const g = value("genomes=")) m_genomes = g;
            if (auto const m = value("max_gap=")) m_max_gap = m;
        }

        static uint64_t NextGeneration() {
            static std::atomic<uint64_t> next{ 1 };
            return next++;
        }

        uint64_t m_generation = NextGeneration();
        std::vector<uint64_t> m_keys;   // (clade, gene, end) of each rule, sorted
        std::vector<Rule> m_rules;      // by key, then partner
        std::vector<uint32_t> m_clades; // sorted
        std::unordered_map<uint32_t, uint32_t> m_chain_of;  // species -> m_chains
        std::vector<std::vector<uint32_t>> m_chains;        // clades with rules, nearest first
        std::map<std::vector<uint32_t>, uint32_t> m_chain_index;
        size_t m_genomes = 0;
        size_t m_max_gap = 0;
    };
}
