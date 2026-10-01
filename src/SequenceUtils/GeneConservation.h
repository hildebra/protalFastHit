// GeneConservation.h - how fast each gene diverges within species compared with the species' other
// genes: one factor per gene id (a marker, the same in every species), 1 for a gene of typical
// conservation, below 1 for conserved genes (ribosomal proteins), above 1 for fast ones. It scales the
// identity margin of a species' depth per gene (GeneMargin; Taxon::OwnIdentityThreshold in Profiler.h):
// a strain's reads on a fast gene sit further below the species' best reads than on a conserved gene.
//
// --build estimates the factors from --full_reference (every genome's copy of each gene, >taxid_geneid)
// against the representatives' genes (reference.fna) and writes gene_conservation.tsv:
//   geneid <tab> factor <tab> species      (species: how many the factor was estimated from)
// For each species and gene, the mean k-mer (Mash, k = 12) distance of other genomes' copies to the
// representative's (at most kMaxCopies, none further than kMaxCopyDistance), over that of the species'
// median gene, so that how far a species' strains are from its representative cancels out; the factor
// of a gene is the median of this ratio over the species with kMinGenes genes or more and a median gene
// at least kMinTypical from the representative, scaled so that the median gene has factor 1, shrunk
// towards 1 by kPrior species and kept within kMinFactor..kMaxFactor. On a simulated world whose genes
// evolve at different rates, the estimates correlated 0.99 with the true rates
// (docs/claude/2026-09-30-depth-margin-stress).
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
