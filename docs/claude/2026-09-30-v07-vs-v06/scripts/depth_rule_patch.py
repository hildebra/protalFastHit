#!/usr/bin/env python3
"""depth_rule_patch.py PROFILER_H - an experimental build: the identity threshold of a taxon's depth (Profiler.h
OwnIdentityThreshold) from $PROTAL_DEPTH_RULE instead of the fixed "98th percentile minus
--depth_identity_margin", to compare dynamic rules on the benchmark's SAMs (--profile_only). Not for the
repository: applied to a copy of the source by dynamic_margin.sh.

Rules (identities as AlignmentIdentity computes them; the gene median = the median over the taxon's genes
with 3 or more reads of each gene's median read identity):
  top:M              98th percentile - M (protal's rule; M = --depth_identity_margin by default)
  genemed:M          gene median - M
  genemedz:Z         gene median - Z * sqrt(p (1 - p) / L), p = 1 - gene median (at least 0.002), L = the
                     median aligned length of the taxon's reads: Z binomial standard deviations of a read's
                     identity at the divergence the gene median shows
  scaled:M0:K        98th percentile - (M0 + K * (1 - gene median))
"""
import sys

path = sys.argv[1]
src = open(path).read()
anchor = """            double OwnIdentityThreshold() const {
                if (m_depth_identity_margin >= 1 || PresentGenes() == 0) return 0;
                return TopIdentity() - m_depth_identity_margin;
            }"""
assert src.count(anchor) == 1, "OwnIdentityThreshold is not as expected"
replacement = """            // [experiment] the median over genes (3+ reads) of each gene's median read identity.
            double GeneMedianIdentity() const {
                std::vector<double> medians;
                for (auto const& [id, gene] : m_genes) {
                    if (gene.m_read_identities.size() < 3) continue;
                    std::vector<float> ids;
                    for (auto const& read : gene.m_read_identities) ids.push_back(read.first);
                    std::nth_element(ids.begin(), ids.begin() + ids.size() / 2, ids.end());
                    medians.push_back(ids.at(ids.size() / 2));
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

            double OwnIdentityThreshold() const {
                if (m_depth_identity_margin >= 1 || PresentGenes() == 0) return 0;
                char const* env = std::getenv("PROTAL_DEPTH_RULE");
                std::string const rule = env ? env : "";
                auto field = [&rule](size_t n) {  // n-th ':'-separated field
                    size_t start = 0;
                    for (size_t i = 0; i < n; i++) start = rule.find(':', start) + 1;
                    return std::stod(rule.substr(start, rule.find(':', start) - start));
                };
                if (rule.rfind("genemedz:", 0) == 0) {
                    double const med = GeneMedianIdentity(), p = std::max(0.002, 1 - med);
                    return med - field(1) * std::sqrt(p * (1 - p) / MedianReadLength());
                }
                if (rule.rfind("genemed:", 0) == 0) return GeneMedianIdentity() - field(1);
                if (rule.rfind("scaled:", 0) == 0) return TopIdentity() - (field(1) + field(2) * (1 - GeneMedianIdentity()));
                if (rule.rfind("top:", 0) == 0) return TopIdentity() - field(1);
                return TopIdentity() - m_depth_identity_margin;
            }"""
src = src.replace(anchor, replacement)
if "#include <cstdlib>" not in src:
    src = src.replace("#pragma once\n", "#pragma once\n#include <cstdlib>\n#include <cmath>\n", 1)
open(path, "w").write(src)
print(f"patched {path}")
