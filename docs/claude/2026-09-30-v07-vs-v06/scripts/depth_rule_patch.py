#!/usr/bin/env python3
"""depth_rule_patch.py PROFILER_H - an experimental build: the identity threshold of a taxon's depth (Profiler.h
OwnIdentityThreshold, and per gene in VerticalCoverage) from $PROTAL_DEPTH_RULE instead of the fixed "98th
percentile minus --depth_identity_margin", to compare rules on the benchmarks' SAMs (--profile_only). Not for
the repository: applied to a copy of the source by dynamic_margin.sh and the stress test's run.sh.

Rules (identities as AlignmentIdentity computes them; the gene median = the median over the taxon's genes
with 3 or more reads of each gene's median read identity):
  top:M              98th percentile - M (protal's rule; M = --depth_identity_margin by default)
  genemed:M          gene median - M
  genemedz:Z         gene median - Z * sqrt(p (1 - p) / L), p = 1 - gene median (at least 0.002), L = the
                     median aligned length of the taxon's reads: Z binomial standard deviations of a read's
                     identity at the divergence the gene median shows
  scaled:M0:K        98th percentile - (M0 + K * (1 - gene median))
  geneq:Q:M          the median over genes of each gene's Q-quantile read identity - M (Q = 0.5: genemed)
  floor:M:F          min(gene median - M, 98th percentile - F): the gene median's threshold, never closer to
                     the top than F (a minor, more distant strain beside a dominant close one)
Per gene, with each gene's conservation factor r (mean 1; $PROTAL_GENE_RATES: lines "geneid factor", 1 for
genes it lacks):
  gtop:M0:K          98th percentile - (M0 + r * K): the margin wide on fast genes, narrow on conserved ones
  gmedc:M:K          1 - r * (D + K) - M, D = the median over genes of (1 - the gene's median identity) / r:
                     the strain's divergence from the reference, measured on every gene and scaled back to it
"""
import sys

path = sys.argv[1]
src = open(path).read()
anchors = [
    """            double OwnIdentityThreshold() const {
                if (m_depth_identity_margin >= 1 || PresentGenes() == 0) return 0;
                return TopIdentity() - m_depth_identity_margin;
            }""",
    """            double OwnIdentityThreshold() const {
                return IdentityThreshold(m_depth_identity_margin);
            }""",
]
anchor = next((a for a in anchors if src.count(a) == 1), None)
assert anchor, "OwnIdentityThreshold is not as expected"
replacement = """            // [experiment] the median over genes (3+ reads) of each gene's `quantile` read identity.
            double GeneMedianIdentity(double quantile = 0.5) const {
                std::vector<double> medians;
                for (auto const& [id, gene] : m_genes) {
                    if (gene.m_read_identities.size() < 3) continue;
                    std::vector<float> ids;
                    for (auto const& read : gene.m_read_identities) ids.push_back(read.first);
                    size_t const at = std::min(ids.size() - 1, static_cast<size_t>(quantile * static_cast<double>(ids.size())));
                    std::nth_element(ids.begin(), ids.begin() + at, ids.end());
                    medians.push_back(ids.at(at));
                }
                if (medians.empty()) return TopIdentity();
                std::nth_element(medians.begin(), medians.begin() + medians.size() / 2, medians.end());
                return medians.at(medians.size() / 2);
            }

            double MedianReadLength() const {
                std::vector<uint32_t> lengths;
                for (auto const& [id, gene] : m_genes) {
                    for (auto const& read : gene.m_read_identities) lengths.push_back(read.second);
                }
                if (lengths.empty()) return 150;
                std::nth_element(lengths.begin(), lengths.begin() + lengths.size() / 2, lengths.end());
                return lengths.at(lengths.size() / 2);
            }

            static std::string DepthRule() {
                char const* env = std::getenv("PROTAL_DEPTH_RULE");
                return env ? env : "";
            }

            static double RuleField(std::string const& rule, size_t n) {  // n-th ':'-separated field
                size_t start = 0;
                for (size_t i = 0; i < n; i++) start = rule.find(':', start) + 1;
                return std::stod(rule.substr(start, rule.find(':', start) - start));
            }

            // [experiment] a gene's conservation factor ($PROTAL_GENE_RATES), 1 if it has none.
            static double GeneRate(size_t geneid) {
                static std::unordered_map<size_t, double> const rates = [] {
                    std::unordered_map<size_t, double> r;
                    if (char const* p = std::getenv("PROTAL_GENE_RATES")) {
                        std::ifstream is(p);
                        size_t g = 0;
                        double f = 0;
                        while (is >> g >> f) r[g] = f;
                    }
                    return r;
                }();
                auto const it = rates.find(geneid);
                return it == rates.end() ? 1.0 : it->second;
            }

            // [experiment] the strain's divergence from the reference: the median over genes (3+ reads) of
            // (1 - the gene's median read identity) / its conservation factor.
            double ScaledDivergence() const {
                if (m_scaled_divergence) return *m_scaled_divergence;
                std::vector<double> d;
                for (auto const& [id, gene] : m_genes) {
                    if (gene.m_read_identities.size() < 3) continue;
                    std::vector<float> ids;
                    for (auto const& read : gene.m_read_identities) ids.push_back(read.first);
                    std::nth_element(ids.begin(), ids.begin() + ids.size() / 2, ids.end());
                    d.push_back((1.0 - ids.at(ids.size() / 2)) / GeneRate(id));
                }
                double v = 0;
                if (!d.empty()) {
                    std::nth_element(d.begin(), d.begin() + d.size() / 2, d.end());
                    v = d.at(d.size() / 2);
                }
                m_scaled_divergence = v;
                return v;
            }

            // [experiment] the threshold of one gene's reads: per gene for gtop and gmedc, else the taxon's.
            double GeneThreshold(size_t geneid, double taxon_threshold) const {
                if (m_depth_identity_margin >= 1) return taxon_threshold;
                std::string const rule = DepthRule();
                if (rule.rfind("gtop:", 0) == 0) return TopIdentity() - (RuleField(rule, 1) + GeneRate(geneid) * RuleField(rule, 2));
                if (rule.rfind("gmedc:", 0) == 0) return 1 - GeneRate(geneid) * (ScaledDivergence() + RuleField(rule, 2)) - RuleField(rule, 1);
                return taxon_threshold;
            }

            double OwnIdentityThreshold() const {
                if (m_depth_identity_margin >= 1 || PresentGenes() == 0) return 0;
                std::string const rule = DepthRule();
                auto field = [&rule](size_t n) { return RuleField(rule, n); };
                if (rule.rfind("genemedz:", 0) == 0) {
                    double const med = GeneMedianIdentity(), p = std::max(0.002, 1 - med);
                    return med - field(1) * std::sqrt(p * (1 - p) / MedianReadLength());
                }
                if (rule.rfind("genemed:", 0) == 0) return GeneMedianIdentity() - field(1);
                if (rule.rfind("geneq:", 0) == 0) return GeneMedianIdentity(field(1)) - field(2);
                if (rule.rfind("floor:", 0) == 0) return std::min(GeneMedianIdentity() - field(1), TopIdentity() - field(2));
                if (rule.rfind("scaled:", 0) == 0) return TopIdentity() - (field(1) + field(2) * (1 - GeneMedianIdentity()));
                if (rule.rfind("top:", 0) == 0) return TopIdentity() - field(1);
                return TopIdentity() - m_depth_identity_margin;
            }"""
src = src.replace(anchor, replacement)
member = "            mutable std::optional<double> m_top_identity;  // cached by TopIdentity\n"
assert src.count(member) == 1, "m_top_identity is not as expected"
src = src.replace(member, member + "            mutable std::optional<double> m_scaled_divergence;  // [experiment] cached by ScaledDivergence\n")
reset = "                m_top_identity.reset();\n"
assert src.count(reset) == 1, "Changed() is not as expected"
src = src.replace(reset, reset + "                m_scaled_divergence.reset();\n")
loop = "                        size_t const bases = gene.MappedLength(min_identity);\n"
assert src.count(loop) == 1, "VerticalCoverage is not as expected"
src = src.replace(loop, "                        size_t const bases = gene.MappedLength(GeneThreshold(geneid, min_identity));\n")
src = src.replace("#pragma once\n", "#pragma once\n#include <cmath>\n#include <cstdlib>\n#include <fstream>\n#include <unordered_map>\n", 1)
open(path, "w").write(src)
print(f"patched {path}")
