#!/usr/bin/env python3
"""rule_patch.py PROFILER_H - an experimental build: the identity threshold of each gene's reads for a taxon's
depth (Profiler.h, Taxon::OwnIdentityThreshold(geneid)) from $PROTAL_DEPTH_RULE instead of protal's rule, to
compare rules on the same SAMs (--profile_only). Not for the repository: applied to a copy of the source by
run_rules.sh.

r is the gene's conservation factor (the database's gene_conservation.tsv, or --gene_conservation; 1 without);
the gene median = the median over the taxon's genes with 3 or more reads of each gene's median read identity.
  (unset)       protal's rule: 98th percentile - (0.03 + (M - 0.03) * r), M = --depth_identity_margin
  gtop:F:M      98th percentile - (F + (M - F) * r): the fixed part F, the margin M on a gene of factor 1
  genemed:M     gene median - M, on every gene
  gmedc:M:K     1 - r * (D + K) - M, D = the median over genes of (1 - the gene's median identity) / r: the
                strain's divergence from the reference measured on every gene, scaled back to each gene
Two rules for relatives' reads, the user's question of 2026-10-01:
  gsplit:T:F:M  each gene's threshold as gtop:F:M; then the depth from one side of the genes when the median
                depths of the conserved genes (r < 1) and of the fast ones (r >= 1; 5 genes with reads or more
                on each side) differ: T >= 1: the fast genes alone if the conserved ones are more than T times
                as deep; T < 1: the conserved genes alone if they are less than T times as deep as the fast
                ones. The first was the premise (a relative the database lacks aligns mostly on the conserved
                genes); split_diag.py found the reverse, its reads on the fast genes, where they are unique.
  gcong:N:R     protal's rule; then each taxon with a detected congener (another species of its genus with a
                model score of 0.5 or more, under protal's rule) whose depth is at least R times its own gets
                the margin N instead of --depth_identity_margin (scaled per gene as protal's rule is). In
                RunProtal.h, after a sample's SAM is profiled.
"""
import os
import sys

path = sys.argv[1]
src = open(path).read()
anchor = """            double OwnIdentityThreshold(uint64_t geneid) const {
                if (m_depth_identity_margin >= 1 || PresentGenes() == 0) return 0;
                return TopIdentity() - gene_conservation::GeneMargin(m_depth_identity_margin, GeneFactor(geneid));
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
                    d.push_back((1.0 - ids.at(ids.size() / 2)) / GeneFactor(id));
                }
                double v = 0;
                if (!d.empty()) {
                    std::nth_element(d.begin(), d.begin() + d.size() / 2, d.end());
                    v = d.at(d.size() / 2);
                }
                m_scaled_divergence = v;
                return v;
            }

            static double RuleField(std::string const& rule, size_t n) {  // n-th ':'-separated field
                size_t start = 0;
                for (size_t i = 0; i < n; i++) start = rule.find(':', start) + 1;
                return std::stod(rule.substr(start, rule.find(':', start) - start));
            }

            double OwnIdentityThreshold(uint64_t geneid) const {
                if (m_depth_identity_margin >= 1 || PresentGenes() == 0) return 0;
                static std::string const rule = [] { char const* e = std::getenv("PROTAL_DEPTH_RULE"); return std::string(e ? e : ""); }();
                double const r = GeneFactor(geneid);
                if (rule.rfind("gtop:", 0) == 0) {
                    double const fixed = RuleField(rule, 1), margin = RuleField(rule, 2);
                    return TopIdentity() - (fixed + (margin - fixed) * r);
                }
                if (rule.rfind("gsplit:", 0) == 0) {
                    double const fixed = RuleField(rule, 2), margin = RuleField(rule, 3);
                    return TopIdentity() - (fixed + (margin - fixed) * r);
                }
                if (rule.rfind("genemed:", 0) == 0) return GeneMedianIdentity() - RuleField(rule, 1);
                if (rule.rfind("gmedc:", 0) == 0) return 1 - r * (ScaledDivergence() + RuleField(rule, 2)) - RuleField(rule, 1);
                return TopIdentity() - gene_conservation::GeneMargin(m_depth_identity_margin, r);
            }"""
src = src.replace(anchor, replacement)
member = "            mutable std::optional<double> m_top_identity;  // cached by TopIdentity\n"
assert src.count(member) == 1, "m_top_identity is not as expected"
src = src.replace(member, member + "            mutable std::optional<double> m_scaled_divergence;  // [experiment] cached by ScaledDivergence\n")
reset = "                m_top_identity.reset();\n"
assert src.count(reset) == 1, "Changed() is not as expected"
src = src.replace(reset, reset + "                m_scaled_divergence.reset();\n")

# gsplit: the depth from the fast genes when the conserved ones are much deeper.
blend = "                    m_vcov = BlendedDepth(Median(vcovs), own_bases, expected_length, vcovs.size(), expected_genes);\n"
assert src.count(blend) == 1, "VerticalCoverage is not as expected"
src = src.replace(blend, blend + "                    if (auto split = SplitDepth()) m_vcov = *split;  // [experiment]\n")
split = """            // [experiment] gsplit:T:F:M - the depth of one side of the genes alone, see rule_patch.py.
            std::optional<double> SplitDepth() const {
                static std::string const rule = [] { char const* e = std::getenv("PROTAL_DEPTH_RULE"); return std::string(e ? e : ""); }();
                if (rule.rfind("gsplit:", 0) != 0 || m_depth_identity_margin >= 1) return std::nullopt;
                std::vector<double> slow, fast;
                size_t slow_bases = 0, fast_bases = 0;
                for (auto const& [geneid, gene] : m_genes) {
                    size_t const bases = gene.MappedLength(OwnIdentityThreshold(geneid));
                    if (bases == 0 || gene.m_gene_length == 0) continue;
                    double const depth = static_cast<double>(bases) / static_cast<double>(gene.m_gene_length);
                    if (GeneFactor(geneid) < 1) {
                        slow.push_back(depth);
                        slow_bases += bases;
                    } else {
                        fast.push_back(depth);
                        fast_bases += bases;
                    }
                }
                if (slow.size() < 5 || fast.size() < 5) return std::nullopt;
                std::sort(slow.begin(), slow.end());
                std::sort(fast.begin(), fast.end());
                double const t = RuleField(rule, 1);
                bool use_fast;
                if (t >= 1 && Median(slow) > t * Median(fast)) use_fast = true;
                else if (t < 1 && Median(slow) < t * Median(fast)) use_fast = false;
                else return std::nullopt;
                size_t expected_length = 0, expected_genes = 0;
                for (auto gene_id : m_genome->GetHittableGenes()) {
                    if (!m_genome->HasGene(gene_id) || (GeneFactor(gene_id) >= 1) != use_fast) continue;
                    expected_length += m_genome->GetGene(gene_id).GetLength();
                    expected_genes++;
                }
                auto const& side = use_fast ? fast : slow;
                return BlendedDepth(Median(side), use_fast ? fast_bases : slow_bases, expected_length, side.size(), expected_genes);
            }

            // Depth of the taxon from its own reads (see OwnIdentityThreshold), estimated by BlendedDepth.
"""
anchor_vcov = "            // Depth of the taxon from its own reads (see OwnIdentityThreshold), estimated by BlendedDepth.\n"
assert src.count(anchor_vcov) == 1, "the VerticalCoverage comment is not as expected"
src = src.replace(anchor_vcov, split)
src = src.replace("#pragma once\n", "#pragma once\n#include <cstdlib>\n", 1)
open(path, "w").write(src)
print(f"patched {path}")

# gcong: after a sample's SAM is profiled, a narrower margin for taxa with a detected, deeper congener.
run_path = os.path.join(os.path.dirname(os.path.dirname(path)), "RunProtal.h")
run = open(run_path).read()
truth = """            std::optional<TruthSet> truth = options.HasProfileTruths() ?"""
assert run.count(truth) == 1, "RunProtal.h is not as expected"
gcong = """            // [experiment] gcong:N:R, see rule_patch.py.
            if (char const* env = std::getenv("PROTAL_DEPTH_RULE"); env && std::string(env).rfind("gcong:", 0) == 0) {
                std::string const rule(env);
                size_t const a = rule.find(':') + 1, b = rule.find(':', a);
                double const narrow = std::stod(rule.substr(a, b - a)), ratio = std::stod(rule.substr(b + 1));
                auto& taxa = profile.GetTaxa();
                std::vector<std::tuple<size_t, int, double>> detected;  // taxid, genus, depth
                for (auto it = taxa.begin(); it != taxa.end(); ++it) {
                    if (filter.Score(it.value()) >= 0.5) {
                        detected.emplace_back(it->first, taxonomy.GetParent(static_cast<int>(it->first)).id, it.value().VerticalCoverage());
                    }
                }
                std::vector<size_t> narrowed;
                for (auto it = taxa.begin(); it != taxa.end(); ++it) {
                    int const genus = taxonomy.GetParent(static_cast<int>(it->first)).id;
                    double const depth = it.value().VerticalCoverage();
                    for (auto const& [d, g, dd] : detected) {
                        if (d != it->first && g == genus && dd >= ratio * depth) {
                            narrowed.push_back(it->first);
                            break;
                        }
                    }
                }
                for (auto t : narrowed) taxa.at(t).SetDepthIdentityMargin(narrow);
                #pragma omp critical(print)
                std::cout << "[experiment] gcong: " << detected.size() << " taxa detected, " << narrowed.size()
                          << " with a detected congener get margin " << narrow << std::endl;
            }

"""
run = run.replace(truth, gcong + truth)
if "#include <cstdlib>" not in run:
    run = run.replace("#include <atomic>\n", "#include <atomic>\n#include <cstdlib>\n#include <tuple>\n", 1)
open(run_path, "w").write(run)
print(f"patched {run_path}")
