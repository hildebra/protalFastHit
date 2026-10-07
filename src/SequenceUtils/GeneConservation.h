// GeneConservation.h - how fast each gene diverges within species compared with the species' other
// genes: one factor per gene id (a marker, the same in every species), 1 for a gene of typical
// conservation, below 1 for conserved genes (ribosomal proteins), above 1 for fast ones. With
// --gene_conservation db it scales the identity margin of a species' depth per gene (GeneMargin;
// Taxon::OwnIdentityThreshold in Profiler.h): a strain's reads on a fast gene sit further below the
// species' best reads than on a conserved gene. By default the margin is the same on every gene: summed over
// three simulated worlds that did better, as the scaled margin also admits more of a missing relative's reads
// (docs/claude/2026-10-01-gene-scaled-margin). Every database built with other genomes' copies in
// --full_reference holds the factors, so the scaled margin, or other uses of the genes' rates, need no rebuild.
//
// --build estimates the factors from --full_reference (every genome's copy of each gene, >taxid_geneid)
// against the representatives' genes (reference.fna) and writes gene_conservation.tsv:
//   geneid <tab> factor <tab> species      (species: how many the factor was estimated from)
// For each species and gene, the mean k-mer (Mash, k = 12) distance of other genomes' copies to the
// representative's (at most kMaxCopies, none further than kMaxCopyDistance), over that of the species'
// median gene, so that how far a species' strains are from its representative cancels out; the factor
// of a gene is the median of this ratio over the species with kMinGenes genes or more and a median gene
// at least kMinTypical from the representative, scaled so that the median gene has factor 1, shrunk
// towards 1 by kPrior species and kept within kMinFactor..kMaxFactor. On two simulated worlds whose genes
// evolve at different rates, the estimates correlated 0.985 and 0.999 with the true rates.
#pragma once

#include <algorithm>
#include <atomic>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <istream>
#include <memory>
#include <optional>
#include <ostream>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

namespace protal::gene_conservation {
    inline const std::string kFileName = "gene_conservation.tsv";

    inline constexpr size_t kK = 12;
    inline constexpr double kMaxCopyDistance = 0.25;  // further: not a strain's copy of the gene
    inline constexpr uint32_t kMaxCopies = 16;        // copies compared per species and gene
    inline constexpr size_t kMinGenes = 10;
    inline constexpr double kMinTypical = 0.002;
    inline constexpr size_t kMinSpecies = 2;
    inline constexpr double kPrior = 3;
    inline constexpr double kMinFactor = 0.25;
    inline constexpr double kMaxFactor = 4;
    inline constexpr uint64_t kGeneIdLimit = uint64_t{1} << 20;  // the index's gene id field

    // The part of the depth identity margin that is the same on every gene: read errors and the spread
    // of read identities (about two standard deviations of a 100-base read's identity at 2-3%
    // differences). The rest scales with the gene's factor.
    inline constexpr double kFixedMargin = 0.03;

    // The depth identity margin on a gene with conservation `factor`: kFixedMargin (or the whole margin
    // if it is smaller) plus the rest times the factor. `margin` on a gene of factor 1; a margin of 1 or
    // more lets every read count on every gene.
    inline double GeneMargin(double margin, double factor) {
        if (margin >= 1) return margin;
        double const fixed = std::min(margin, kFixedMargin);
        return fixed + (margin - fixed) * factor;
    }

    // The distinct k-mers (kK bases, 2 bits each) of the A/C/G/T stretches of seq, sorted.
    inline std::vector<uint32_t> Kmers(std::string_view seq) {
        std::vector<uint32_t> kmers;
        if (seq.size() >= kK) kmers.reserve(seq.size() - kK + 1);
        uint32_t constexpr mask = (uint32_t{1} << (2 * kK)) - 1;
        uint32_t code = 0;
        size_t run = 0;
        for (char const c : seq) {
            uint32_t base = 0;
            switch (c) {
                case 'A': case 'a': base = 0; break;
                case 'C': case 'c': base = 1; break;
                case 'G': case 'g': base = 2; break;
                case 'T': case 't': base = 3; break;
                default: run = 0; continue;
            }
            code = ((code << 2) | base) & mask;
            if (++run >= kK) kmers.push_back(code);
        }
        std::sort(kmers.begin(), kmers.end());
        kmers.erase(std::unique(kmers.begin(), kmers.end()), kmers.end());
        return kmers;
    }

    // A bijection of 32-bit k-mer codes (lowbias32): k-mers of similar codes get unrelated hashes, so that the
    // smallest hashes are a random sample of a gene's k-mers.
    inline uint32_t MixKmer(uint32_t x) {
        x ^= x >> 16;
        x *= 0x7feb352dU;
        x ^= x >> 15;
        x *= 0x846ca68bU;
        x ^= x >> 16;
        return x;
    }

    // The `size` smallest hashes (MixKmer) of the distinct k-mers of seq (kK bases, 2 bits each, of its A/C/G/T
    // stretches, as Kmers takes them), sorted; all of them for a gene with fewer distinct k-mers: Mash's bottom
    // sketch, for a run's congener distances (SampleContext.h) and the build's suspect copies (GeneIncongruence.h).
    //
    // The k-mers are hashed as they come, repeats included (MixKmer is a bijection: distinct k-mers are distinct
    // hashes and a repeat repeats its hash). The 2 x size smallest hashes are selected in linear time and sorted;
    // whenever they hold `size` distinct values, the smallest distinct ones are among them, as every hash outside is
    // at least as large as any inside. Only a gene of many repeats sorts all its hashes. The first version collected
    // and sorted every distinct k-mer first (test SketchesEqualTheReferencesToTheLastHash).
    inline std::vector<uint32_t> BottomSketch(std::string_view seq, size_t size) {
        thread_local std::vector<uint32_t> hashes;
        hashes.clear();
        uint32_t constexpr mask = (uint32_t{1} << (2 * kK)) - 1;
        uint32_t code = 0;
        size_t run = 0;
        for (char const c : seq) {
            uint32_t base = 0;
            switch (c) {
                case 'A': case 'a': base = 0; break;
                case 'C': case 'c': base = 1; break;
                case 'G': case 'g': base = 2; break;
                case 'T': case 't': base = 3; break;
                default: run = 0; continue;
            }
            code = ((code << 2) | base) & mask;
            if (++run >= kK) hashes.push_back(MixKmer(code));
        }
        size_t const n = hashes.size();
        size_t const probe = std::min(n, 2 * size);
        auto const probe_end = hashes.begin() + static_cast<std::ptrdiff_t>(probe);
        if (probe < n) std::nth_element(hashes.begin(), probe_end, hashes.end());
        std::sort(hashes.begin(), probe_end);
        size_t distinct = 0;
        for (size_t i = 0; i < probe; i++) distinct += i == 0 || hashes[i] != hashes[i - 1];
        auto end = probe_end;
        if (distinct < size && probe < n) {  // too many repeats among the smallest: all of them, sorted
            std::sort(hashes.begin(), hashes.end());
            end = hashes.end();
        }
        end = std::unique(hashes.begin(), end);
        size_t const keep = std::min(size, static_cast<size_t>(end - hashes.begin()));
        return std::vector<uint32_t>(hashes.begin(), hashes.begin() + static_cast<std::ptrdiff_t>(keep));
    }

    // SketchDistance of two sketches whose s smallest union hashes hold `shared` in both: the one place the distance
    // is computed, so that a bound on `shared` (GeneIncongruence.h) says exactly which distances a pair can reach.
    inline double SketchDistanceOf(size_t shared, size_t s) {
        if (s == 0 || shared == 0) return 1;
        double const jaccard = static_cast<double>(shared) / static_cast<double>(s);
        // max with 0: identical sketches give -log(1) = -0.0, which the reports would print as "-0".
        return std::min(1.0, std::max(0.0, -std::log(2 * jaccard / (1 + jaccard)) / static_cast<double>(kK)));
    }

    // The Mash distance of two genes from their bottom sketches (BottomSketch): of the s smallest hashes of their
    // union (s the smaller sketch's size), the share in both estimates their k-mers' Jaccard index j;
    // -ln(2j / (1 + j)) / k, 1 if they share none, as MashDistance on all k-mers.
    inline double SketchDistance(std::vector<uint32_t> const& a, std::vector<uint32_t> const& b) {
        size_t const s = std::min(a.size(), b.size());
        if (s == 0) return 1;
        size_t i = 0, j = 0, seen = 0, shared = 0;
        while (seen < s && i < a.size() && j < b.size()) {
            if (a[i] < b[j]) {
                i++;
            } else if (b[j] < a[i]) {
                j++;
            } else {
                shared++;
                i++;
                j++;
            }
            seen++;
        }
        return SketchDistanceOf(shared, s);
    }

    // The Mash distance of two sorted k-mer sets, an estimate of the share of bases that differ, from
    // their Jaccard index j: -ln(2j / (1 + j)) / k; 1 if they share none.
    inline double MashDistance(std::vector<uint32_t> const& a, std::vector<uint32_t> const& b) {
        size_t shared = 0;
        auto i = a.begin();
        auto j = b.begin();
        while (i != a.end() && j != b.end()) {
            if (*i < *j) {
                ++i;
            } else if (*j < *i) {
                ++j;
            } else {
                ++shared;
                ++i;
                ++j;
            }
        }
        if (shared == 0) return 1;
        double const jaccard = static_cast<double>(shared) / static_cast<double>(a.size() + b.size() - shared);
        return std::min(1.0, -std::log(2 * jaccard / (1 + jaccard)) / static_cast<double>(kK));
    }

    // The median of values (the mean of the two middle ones for an even count); 0 if empty.
    inline double Median(std::vector<double> values) {
        if (values.empty()) return 0;
        size_t const mid = values.size() / 2;
        std::nth_element(values.begin(), values.begin() + mid, values.end());
        double const upper = values[mid];
        if (values.size() % 2 == 1) return upper;
        return (*std::max_element(values.begin(), values.begin() + mid) + upper) / 2;
    }

    // The factors by gene id; a gene without one has factor 1.
    class Table {
    public:
        double Factor(uint64_t geneid) const {
            return geneid < m_factors.size() && m_factors[geneid] > 0 ? m_factors[geneid] : 1.0;
        }

        bool Empty() const { return m_genes == 0; }
        size_t Genes() const { return m_genes; }

        // Whether gene `geneid` has a factor of its own.
        bool Has(uint64_t geneid) const { return geneid < m_factors.size() && m_factors[geneid] > 0; }

        void Set(uint64_t geneid, double factor, uint32_t species = 0) {
            if (geneid >= m_factors.size()) {
                m_factors.resize(geneid + 1, 0);
                m_species.resize(geneid + 1, 0);
            }
            if (m_factors[geneid] <= 0) m_genes++;
            m_factors[geneid] = static_cast<float>(factor);
            m_species[geneid] = species;
        }

        // The lowest and highest factor (1, 1 without any).
        std::pair<double, double> Range() const {
            double low = 0, high = 0;
            for (float const f : m_factors) {
                if (f <= 0) continue;
                low = low == 0 ? f : std::min<double>(low, f);
                high = std::max<double>(high, f);
            }
            return m_genes == 0 ? std::pair<double, double>{1, 1} : std::pair<double, double>{low, high};
        }

        // gene_conservation.tsv: lines "geneid <tab> factor [<tab> species]", a header line starting
        // with "geneid" and lines starting with '#' skipped. Returns the first problem with its line
        // ("line 3: ..."), or an empty string.
        std::string Read(std::istream& is) {
            std::string line;
            size_t number = 0;
            while (std::getline(is, line)) {
                number++;
                if (!line.empty() && line.back() == '\r') line.pop_back();
                if (line.empty() || line[0] == '#' || line.rfind("geneid", 0) == 0) continue;
                auto problem = [&](std::string const& what) { return "line " + std::to_string(number) + ": " + what; };
                std::vector<std::string_view> fields;
                std::string_view rest(line);
                while (true) {
                    size_t const tab = rest.find('\t');
                    fields.push_back(rest.substr(0, tab));
                    if (tab == std::string_view::npos) break;
                    rest.remove_prefix(tab + 1);
                }
                if (fields.size() < 2) return problem("expected a gene id and a factor, separated by a tab");
                uint64_t geneid = 0;
                auto const id = std::from_chars(fields[0].data(), fields[0].data() + fields[0].size(), geneid);
                if (id.ec != std::errc() || id.ptr != fields[0].data() + fields[0].size()) return problem("the gene id is not a number");
                if (geneid >= kGeneIdLimit) return problem("gene id " + std::to_string(geneid) + " is too large");
                double factor = 0;
                auto const f = std::from_chars(fields[1].data(), fields[1].data() + fields[1].size(), factor);
                if (f.ec != std::errc() || f.ptr != fields[1].data() + fields[1].size() || !std::isfinite(factor)) {
                    return problem("the factor is not a number");
                }
                if (factor <= 0 || factor > 100) return problem("the factor must be above 0 and at most 100");
                uint32_t species = 0;
                if (fields.size() >= 3) std::from_chars(fields[2].data(), fields[2].data() + fields[2].size(), species);
                if (geneid < m_factors.size() && m_factors[geneid] > 0) return problem("gene " + std::to_string(geneid) + " is listed twice");
                Set(geneid, factor, species);
            }
            if (is.bad()) return "read error";
            return {};
        }

        void Write(std::ostream& os) const {
            os << "geneid\tfactor\tspecies\n";
            char buffer[32];
            for (size_t id = 0; id < m_factors.size(); id++) {
                if (m_factors[id] <= 0) continue;
                std::snprintf(buffer, sizeof(buffer), "%.4f", m_factors[id]);
                os << id << '\t' << buffer << '\t' << m_species[id] << '\n';
            }
        }

    private:
        std::vector<float> m_factors;     // 0: none
        std::vector<uint32_t> m_species;  // species each factor was estimated from
        size_t m_genes = 0;
    };

    // How each gene differs between congeneric species, beside how it differs within species (the factors): --build
    // compares the representatives' copies (reference.fna) of species of one genus and writes gene_congeners.tsv next
    // to the database, a report that queries do not read (docs/claude/2026-10-01-conservation-pattern):
    //   geneid  within_factor  between_factor  pairs  species  identical_share  near_identical_share
    // between_factor: the median over pairs of species of the gene's k-mer distance over that of the pair's median
    // gene, scaled so that the median gene has 1, as the within-species factors are; pairs: the pairs it was estimated
    // from (kMinGenes shared genes, the median gene kMinTypical to kMaxPairTypical apart); species: the species whose
    // copy was compared with a congener's; identical_share and near_identical_share: of those, the share whose nearest
    // congener's copy is identical (the same k-mers) or nearer than kNearIdentical, so that most of its reads fit both.
    // The factors assume that a gene conserved within species is conserved between them; a relative the database lacks
    // aligns best where it differs least from its congeners, and a read of a gene the same in two congeners fits both.
    inline const std::string kCongenersFileName = "gene_congeners.tsv";
    inline constexpr size_t kMaxSpeciesPerGenus = 64;  // of a larger genus, a fixed sample (by a hash of the taxid)
    inline constexpr size_t kCongenerSteps = 4;        // each species against the next 4 of its genus's sample, cyclically
    inline constexpr double kMaxPairTypical = 0.2;     // a pair whose median gene is further apart: saturated, no ratios
    inline constexpr double kNearIdentical = 0.01;     // about one difference in 100 bases

    struct CongenerRow {
        uint64_t geneid = 0;
        double within = 0;   // its factor (0: none)
        double between = std::nan("");  // its between-species factor (NaN: fewer than kMinSpecies pairs)
        size_t pairs = 0, species = 0, identical = 0, near = 0;
    };

    struct CongenerEstimate {
        std::vector<CongenerRow> genes;  // the genes compared between congeners, by gene id
        size_t genera = 0, species = 0, pairs = 0;
        double spearman = std::nan("");  // between the within- and between-species factors, over the genes with both
        size_t correlated = 0;           // the genes it is over
        // The genes with a factor below 1 (conserved) and the others with one: the median between-species factor,
        // and the share of compared species whose nearest congener's copy is identical, or nearer than kNearIdentical.
        double conserved_between = std::nan(""), fast_between = std::nan("");
        double conserved_identical = std::nan(""), fast_identical = std::nan("");
        double conserved_near = std::nan(""), fast_near = std::nan("");

        void Write(std::ostream& os) const {
            os << "geneid\twithin_factor\tbetween_factor\tpairs\tspecies\tidentical_share\tnear_identical_share\n";
            for (auto const& g : genes) {
                os << g.geneid << '\t';
                if (g.within > 0) os << g.within; else os << "NA";
                os << '\t';
                if (!std::isnan(g.between)) os << g.between; else os << "NA";
                os << '\t' << g.pairs << '\t' << g.species << '\t'
                   << (g.species ? static_cast<double>(g.identical) / static_cast<double>(g.species) : 0.0) << '\t'
                   << (g.species ? static_cast<double>(g.near) / static_cast<double>(g.species) : 0.0) << '\n';
            }
        }
    };

    // Spearman's rank correlation of x and y (ties ranked by their mean); NaN for fewer than 3 values.
    inline double Spearman(std::vector<double> const& x, std::vector<double> const& y) {
        size_t const n = x.size();
        if (n < 3 || y.size() != n) return std::nan("");
        auto ranks = [n](std::vector<double> const& v) {
            std::vector<size_t> order(n);
            for (size_t i = 0; i < n; i++) order[i] = i;
            std::sort(order.begin(), order.end(), [&v](size_t a, size_t b) { return v[a] < v[b]; });
            std::vector<double> r(n);
            for (size_t i = 0; i < n;) {
                size_t j = i;
                while (j + 1 < n && v[order[j + 1]] == v[order[i]]) j++;
                for (size_t k = i; k <= j; k++) r[order[k]] = (static_cast<double>(i) + static_cast<double>(j)) / 2;
                i = j + 1;
            }
            return r;
        };
        auto const rx = ranks(x), ry = ranks(y);
        double const mean = (static_cast<double>(n) - 1) / 2;
        double num = 0, dx = 0, dy = 0;
        for (size_t i = 0; i < n; i++) {
            num += (rx[i] - mean) * (ry[i] - mean);
            dx += (rx[i] - mean) * (rx[i] - mean);
            dy += (ry[i] - mean) * (ry[i] - mean);
        }
        return dx > 0 && dy > 0 ? num / std::sqrt(dx * dy) : std::nan("");
    }

    // Compares the genes of congeneric species. genera: the species (taxids) of each genus; genes_of(taxid): a
    // species' representative genes as (geneid, sequence) pairs, called once per species compared, from any thread;
    // within: the within-species factors.
    template<typename GenesOf>
    CongenerEstimate CompareCongeners(std::vector<std::vector<uint32_t>> const& genera, GenesOf&& genes_of,
                                      Table const& within, int threads) {
        struct GeneSums {
            std::vector<double> ratios;
            size_t pairs = 0, species = 0, identical = 0, near = 0;
        };
        std::vector<GeneSums> sums;
        CongenerEstimate estimate;
#pragma omp parallel for schedule(dynamic, 1) num_threads(std::max(threads, 1)) default(none) shared(genera, genes_of, sums, estimate)
        for (size_t g = 0; g < genera.size(); g++) {
            std::vector<uint32_t> members = genera[g];
            if (members.size() < 2) continue;
            if (members.size() > kMaxSpeciesPerGenus) {
                auto const hash = [](uint32_t t) { return (static_cast<uint64_t>(t) * 0x9E3779B97F4A7C15ull) >> 32; };
                std::sort(members.begin(), members.end(), [&hash](uint32_t a, uint32_t b) { return hash(a) < hash(b); });
                members.resize(kMaxSpeciesPerGenus);
            }
            std::sort(members.begin(), members.end());
            size_t const n = members.size();
            std::vector<std::vector<std::pair<uint64_t, std::vector<uint32_t>>>> kmers(n);  // by geneid
            for (size_t i = 0; i < n; i++) {
                for (auto const& [geneid, seq] : genes_of(members[i])) kmers[i].emplace_back(geneid, Kmers(seq));
                std::sort(kmers[i].begin(), kmers[i].end(), [](auto const& a, auto const& b) { return a.first < b.first; });
            }
            std::vector<std::vector<std::pair<uint64_t, double>>> nearest(n);  // by geneid: the nearest congener's distance
            for (size_t i = 0; i < n; i++) {
                for (auto const& [geneid, _] : kmers[i]) nearest[i].emplace_back(geneid, 2.0);
            }
            std::vector<std::pair<size_t, size_t>> pair_list;
            for (size_t i = 0; i < n; i++) {
                for (size_t step = 1; step <= kCongenerSteps && step < n; step++) {
                    size_t const j = (i + step) % n;
                    pair_list.emplace_back(std::min(i, j), std::max(i, j));
                }
            }
            std::sort(pair_list.begin(), pair_list.end());
            pair_list.erase(std::unique(pair_list.begin(), pair_list.end()), pair_list.end());
            std::vector<std::pair<uint64_t, double>> local_ratios;  // geneid, ratio
            size_t used_pairs = 0;
            std::vector<double> distances;
            std::vector<std::tuple<uint64_t, size_t, size_t, double>> shared;  // geneid, index in i, index in j, distance
            for (auto const [i, j] : pair_list) {
                shared.clear();
                for (size_t a = 0, b = 0; a < kmers[i].size() && b < kmers[j].size();) {
                    if (kmers[i][a].first < kmers[j][b].first) a++;
                    else if (kmers[i][a].first > kmers[j][b].first) b++;
                    else {
                        shared.emplace_back(kmers[i][a].first, a, b, MashDistance(kmers[i][a].second, kmers[j][b].second));
                        a++;
                        b++;
                    }
                }
                if (shared.size() < kMinGenes) continue;
                for (auto const& [geneid, a, b, d] : shared) {
                    nearest[i][a].second = std::min(nearest[i][a].second, d);
                    nearest[j][b].second = std::min(nearest[j][b].second, d);
                }
                distances.clear();
                for (auto const& t : shared) distances.push_back(std::get<3>(t));
                double const typical = Median(distances);
                if (typical < kMinTypical || typical > kMaxPairTypical) continue;
                used_pairs++;
                for (auto const& [geneid, a, b, d] : shared) local_ratios.emplace_back(geneid, d / typical);
            }
#pragma omp critical(congener_sums)
            {
                estimate.genera++;
                estimate.species += n;
                estimate.pairs += used_pairs;
                for (auto const& [geneid, ratio] : local_ratios) {
                    if (geneid >= sums.size()) sums.resize(geneid + 1);
                    sums[geneid].ratios.push_back(ratio);
                    sums[geneid].pairs++;
                }
                for (size_t i = 0; i < n; i++) {
                    for (auto const& [geneid, d] : nearest[i]) {
                        if (d > 1.5) continue;  // no congener's copy compared
                        if (geneid >= sums.size()) sums.resize(geneid + 1);
                        sums[geneid].species++;
                        sums[geneid].identical += d == 0;
                        sums[geneid].near += d < kNearIdentical;
                    }
                }
            }
        }
        std::vector<double> medians;
        for (auto const& s : sums) {
            if (s.ratios.size() >= kMinSpecies) medians.push_back(Median(s.ratios));
        }
        double const typical = Median(medians);
        std::vector<double> x, y, conserved_between, fast_between;
        size_t c_species = 0, c_identical = 0, c_near = 0, f_species = 0, f_identical = 0, f_near = 0;
        for (uint64_t id = 0; id < sums.size(); id++) {
            auto const& s = sums[id];
            if (s.species == 0 && s.ratios.empty()) continue;
            CongenerRow row;
            row.geneid = id;
            row.within = within.Has(id) ? within.Factor(id) : 0;
            row.between = s.ratios.size() >= kMinSpecies && typical > 0 ? Median(s.ratios) / typical : std::nan("");
            row.pairs = s.pairs;
            row.species = s.species;
            row.identical = s.identical;
            row.near = s.near;
            estimate.genes.push_back(row);
            if (row.within <= 0) continue;
            bool const conserved = row.within < 1;
            (conserved ? c_species : f_species) += s.species;
            (conserved ? c_identical : f_identical) += s.identical;
            (conserved ? c_near : f_near) += s.near;
            if (std::isnan(row.between)) continue;
            x.push_back(row.within);
            y.push_back(row.between);
            (conserved ? conserved_between : fast_between).push_back(row.between);
        }
        estimate.spearman = Spearman(x, y);
        estimate.correlated = x.size();
        if (!conserved_between.empty()) estimate.conserved_between = Median(conserved_between);
        if (!fast_between.empty()) estimate.fast_between = Median(fast_between);
        auto share = [](size_t k, size_t n) { return n ? static_cast<double>(k) / static_cast<double>(n) : std::nan(""); };
        estimate.conserved_identical = share(c_identical, c_species);
        estimate.conserved_near = share(c_near, c_species);
        estimate.fast_identical = share(f_identical, f_species);
        estimate.fast_near = share(f_near, f_species);
        return estimate;
    }

    struct Estimate {
        Table table;
        size_t species = 0;              // species the factors were estimated from
        size_t species_with_copies = 0;  // species with another genome's copy of a gene
        size_t copies = 0;               // copies compared
    };

    // Collects the distances of the genes' copies to the representatives' genes, from any number of
    // threads, and estimates the factors from them (Finish).
    class Estimator {
    public:
        static uint64_t Key(uint64_t taxid, uint64_t geneid) {
            return (taxid << 20) | geneid;
        }

        // keys: the reference genes, as Key(taxid, geneid).
        explicit Estimator(std::vector<uint64_t> keys) : m_keys(std::move(keys)) {
            std::sort(m_keys.begin(), m_keys.end());
            m_keys.erase(std::unique(m_keys.begin(), m_keys.end()), m_keys.end());
            m_slots = std::make_unique<Slot[]>(m_keys.size());
        }

        // The slot of gene (taxid, geneid) for one more copy: none if it is not a reference gene or
        // kMaxCopies of its copies were taken already. Thread-safe.
        std::optional<size_t> Take(uint64_t taxid, uint64_t geneid) {
            if (geneid >= kGeneIdLimit) return std::nullopt;
            uint64_t const key = Key(taxid, geneid);
            auto const it = std::lower_bound(m_keys.begin(), m_keys.end(), key);
            if (it == m_keys.end() || *it != key) return std::nullopt;
            size_t const slot = static_cast<size_t>(it - m_keys.begin());
            if (m_slots[slot].taken.fetch_add(1, std::memory_order_relaxed) >= kMaxCopies) return std::nullopt;
            return slot;
        }

        // Records the copy of the gene in `slot` against the representative's, rep. A copy identical to
        // it (read as the 2-bit gene store does: other bases as A) is the representative's own, or a
        // genome identical to it there. Thread-safe.
        void Add(size_t slot, std::string_view rep, std::string_view copy) {
            Slot& s = m_slots[slot];
            if (SameAsStored(rep, copy)) {
                s.zeros.fetch_add(1, std::memory_order_relaxed);
                return;
            }
            double const distance = MashDistance(Kmers(rep), Kmers(copy));
            if (distance > kMaxCopyDistance) return;
            s.others.fetch_add(1, std::memory_order_relaxed);
            s.sum.fetch_add(static_cast<uint32_t>(std::lround(distance * 1e6)), std::memory_order_relaxed);
        }

        Estimate Finish() const {
            Estimate estimate;
            std::vector<std::vector<double>> ratios;
            std::vector<std::pair<uint64_t, double>> genes;  // geneid, distance
            std::vector<double> distances;
            for (size_t begin = 0; begin < m_keys.size();) {
                uint64_t const taxid = m_keys[begin] >> 20;
                size_t end = begin;
                genes.clear();
                for (; end < m_keys.size() && (m_keys[end] >> 20) == taxid; end++) {
                    Slot const& s = m_slots[end];
                    uint32_t const zeros = s.zeros.load(), others = s.others.load();
                    // One identical copy is the representative's own.
                    uint32_t const copies = others + zeros - (zeros > 0 ? 1 : 0);
                    estimate.copies += others + zeros;
                    if (copies == 0) continue;
                    genes.emplace_back(m_keys[end] & (kGeneIdLimit - 1), static_cast<double>(s.sum.load()) / 1e6 / copies);
                }
                begin = end;
                if (genes.empty()) continue;
                estimate.species_with_copies++;
                if (genes.size() < kMinGenes) continue;
                distances.clear();
                for (auto const& [id, d] : genes) distances.push_back(d);
                double const typical = Median(distances);
                if (typical < kMinTypical) continue;
                estimate.species++;
                for (auto const& [id, d] : genes) {
                    if (id >= ratios.size()) ratios.resize(id + 1);
                    ratios[id].push_back(d / typical);
                }
            }
            std::vector<double> medians(ratios.size(), 0), all;
            for (size_t id = 0; id < ratios.size(); id++) {
                if (ratios[id].size() < kMinSpecies) continue;
                medians[id] = Median(ratios[id]);
                all.push_back(medians[id]);
            }
            double const typical = Median(all);
            if (typical <= 0) return estimate;
            for (size_t id = 0; id < ratios.size(); id++) {
                if (ratios[id].size() < kMinSpecies) continue;
                double const n = static_cast<double>(ratios[id].size());
                double const factor = (n * medians[id] / typical + kPrior) / (n + kPrior);
                estimate.table.Set(id, std::clamp(factor, kMinFactor, kMaxFactor), static_cast<uint32_t>(ratios[id].size()));
            }
            return estimate;
        }

    private:
        struct Slot {
            std::atomic<uint32_t> taken{0}, zeros{0}, others{0};
            std::atomic<uint32_t> sum{0};  // the others' distances in millionths
        };

        static bool SameAsStored(std::string_view rep, std::string_view copy) {
            if (rep.size() != copy.size()) return false;
            for (size_t i = 0; i < copy.size(); i++) {
                char c = static_cast<char>(copy[i] & ~0x20);  // upper case
                if (c != 'C' && c != 'G' && c != 'T') c = 'A';
                if (c != rep[i]) return false;
            }
            return true;
        }

        std::vector<uint64_t> m_keys;
        std::unique_ptr<Slot[]> m_slots;
    };
}
