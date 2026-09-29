//
// Created by fritsche on 07/10/22.
//

#pragma once

#include "Strain.h"
#include "AlignmentUtils.h"
#include <vector>
#include <sparse_map.h>
#include <numeric>
#include <unordered_set>
#include <ProfilerDefinitions.h>

#include "LineSplitter.h"
#include "Taxonomy.h"
#include "InternalReadAlignment.h"
#include "Constants.h"
#include "SNPUtils.h"
#include "ScoreAlignments.h"
#include "gzstream/gzstream.h"
#include "cPMML.h"
#include "sparse_map.h"
#include "Benchmark.h"
#include "RunStatus.h"
#include <algorithm>
#include <charconv>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <set>
#include <string_view>
#include <string>
#include <ranges>
#include <cmath>

namespace protal {
    inline bool IsDigit(std::string &test) {
        return false;
    }

    using TruthSet = tsl::robin_set<uint32_t>;

    // Reads the true species of a sample. Each line names one species either by a GTDB-style lineage
    // in any tab-separated field ("d__...;...;s__Genus species", as simulate_metagenomes writes it) or
    // by an internal taxid in the first field. A header line, blank lines and '#' comments are skipped;
    // lineages without a species (e.g. unclassified genomes) are reported and skipped.
    inline TruthSet GetTruth(std::string const& file_path, taxonomy::IntTaxonomy& taxonomy) {
        std::ifstream is(file_path, std::ios::in);
        if (!is) {
            std::cerr << "Cannot read truth file " << file_path << std::endl;
            exit(8);
        }
        TruthSet truths;
        std::vector<std::string> tokens;
        std::vector<std::string> ranks;
        std::string delim = "\t";
        auto is_lineage = [](std::string const& s) {
            return s.rfind("d__", 0) == 0 || s.rfind("k__", 0) == 0 || s.rfind("s__", 0) == 0 || s.find(";s__") != std::string::npos;
        };
        auto is_integer = [](std::string const& s) {
            return !s.empty() && std::all_of(s.begin(), s.end(), [](char c) { return std::isdigit(static_cast<unsigned char>(c)); });
        };

        std::string line;
        size_t line_no = 0;
        while (std::getline(is, line)) {
            line_no++;
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty() || line[0] == '#') continue;
            LineSplitter::Split(line, delim, tokens);

            auto lineage = std::find_if(tokens.begin(), tokens.end(), is_lineage);
            if (lineage == tokens.end()) {
                if (is_integer(tokens.front())) {
                    auto taxid = std::stoul(tokens.front());
                    if (!taxonomy.map.contains(taxid)) {
                        std::cerr << "Truth file " << file_path << ", line " << line_no << ": taxid " << taxid << " is not in the taxonomy" << std::endl;
                        continue;
                    }
                    truths.insert(taxid);
                } else if (line_no > 1) {
                    std::cerr << "Truth file " << file_path << ", line " << line_no << ": no lineage or taxid, ignored" << std::endl;
                }
                continue;
            }

            LineSplitter::Split(*lineage, ";", ranks);
            std::string species;
            for (auto& e : ranks) {
                if (e.rfind("s__", 0) == 0 && e.size() > 3) species = e;
            }
            if (species.empty()) {
                std::cerr << "No species classification in truth file " << file_path << ", line " << line_no << std::endl;
                continue;
            }
            if (!taxonomy.string_to_id.contains(species)) {
                std::cerr << "Species unknown to taxonomy: " << species << std::endl;
                continue;
            }

            truths.insert(taxonomy.Get(species));
        }
        return truths;
    }


    namespace profiler {
        // Forward declaration
        class TaxonFilterForest;
        class TaxonFilter;

        struct MicrobialProfile;
        struct Gene;
        struct Taxon;

        using GeneMap = tsl::sparse_map<uint32_t, Gene>;
        using TaxonMap = tsl::sparse_map<uint32_t, Taxon>;
        using GeneRef = protal::Gene;

        using TaxonFilterObj = TaxonFilterForest;

        // Data structures

        struct AlleleCounts {
            uint32_t filtered = 0;
            uint32_t mono = 0;
            uint32_t bi = 0;
            uint32_t tri = 0;
            uint32_t tetra = 0;

            uint32_t Multi() {
                return bi + tri + tetra;
            }

            uint32_t Filtered() {
                return filtered;
            }

            uint32_t Noisy() {
                return Multi() + filtered;
            }

            uint32_t AllValid() {
                return Multi() + mono;
            }

            uint32_t All() {
                return Multi() + mono;
            }

            std::string ToString() {
                return "Filtered\t" + std::to_string(filtered) + '\n' +
                       "Mono\t" + std::to_string(mono) + '\n' +
                       "Bi\t" + std::to_string(bi) + '\n' +
                       "Tri\t" + std::to_string(tri) + '\n' +
                       "Tetra\t" + std::to_string(tetra);
            }
        };

        class Gene {
        public:
//            using SNPs = std::vector<SNP>;

            size_t m_mapped_reads = 0;
            size_t m_mapped_length = 0;
            size_t m_total_kmers = 0;
            size_t m_mapq_sum = 0;
            size_t m_unique_mers = 0;
            size_t m_unique_two_mers = 0;
            size_t m_unique_mer_reads = 0;
            size_t m_unique_two_mer_reads = 0;
            double m_ani_sum = 0;

            size_t m_gene_length = 0;
            StrainLevelContainer m_strain_level;
            // Identity and aligned reference length of every read, for depth from a taxon's own reads.
            std::vector<std::pair<float, uint32_t>> m_read_identities;
            double m_identity_bases = 0;  // aligned reference bases times their read's identity
            size_t m_fragments = 0;       // reads, a pair counting once
            size_t m_last_read = SIZE_MAX;

            // Aligned reference bases of the reads with at least `min_identity`.
            size_t MappedLength(double min_identity) const {
                size_t bases = 0;
                for (auto const& [identity, length] : m_read_identities) {
                    if (identity >= min_identity) bases += length;
                }
                return bases;
            }


            GeneRef* m_gene_ref = nullptr;

            Gene(protal::Gene& gene_ref) :
                    m_gene_ref(&gene_ref),
                    m_strain_level(gene_ref) {
            };

            size_t Coverage(size_t above=0) {
                auto cov_vec = GetStrainLevel().GetSequenceRangeHandler().CalculateCoverageVector2();
                auto cov = std::count_if(cov_vec.begin(), cov_vec.end(), [above](const uint32_t e){ return e > above; });
                return cov;
            }

            AlleleCounts AlleleSNPCounts(uint32_t min_cov, size_t min_qual_sum) {
                std::vector<size_t> alleles(5, 0);
                GetAlleles(alleles, min_cov, min_qual_sum);
                AlleleCounts counts;


                counts.filtered = alleles[0];
                counts.mono = alleles[1];
                counts.bi = alleles[2];
                counts.tri = alleles[3];
                counts.tetra = alleles[4];

                return counts;
            }

            void AddRead(GenePos pos) {
                m_mapped_reads++;
            }

            void GetAlleles(std::vector<size_t>& allele_counts, uint32_t min_cov=0, uint32_t min_qual_sum=0) const {
                auto& vars = m_strain_level.GetVariantHandler().GetVariants();

                for (const auto& [pos, _] : vars) {
                    auto const& var = vars.at(pos);
                    size_t total_cov = std::accumulate(var.begin(), var.end(), size_t{0}, [](size_t acc, Variant const& var) {
                        return acc + var.Observations();
                    });

                    size_t valid_alleles = 0;
                    for (auto& v : var) {
                        valid_alleles += (v.Observations() >= min_cov & v.QualitySum() >= min_qual_sum);
                    }
                    if (valid_alleles >= allele_counts.size()) {
                        allele_counts.resize(valid_alleles+1, 0);
                    }
                    allele_counts[valid_alleles]++;
                }
            }

            const StrainLevelContainer& GetStrainLevel() const {
                return m_strain_level;
            }

            StrainLevelContainer& GetStrainLevel() {
                return m_strain_level;
            }

            std::string GetStatisticsString(std::string delimiter="\t") {
                std::string s = "";

                auto alleles = AlleleSNPCounts(0, 0);


                auto alleles_filtered = AlleleSNPCounts(2, 60);

                auto cov_vec = m_strain_level.GetSequenceRangeHandler().CalculateCoverageVector2();
                auto cov_sum = std::accumulate(cov_vec.begin(), cov_vec.end(), size_t{0});
                auto cov = m_strain_level.GetSequenceRangeHandler().CoveredPortion();
                s += std::to_string(m_mapped_reads) + '\t';
                s += std::to_string(cov_sum) + '\t';
                s += std::to_string(m_gene_length) + '\t';
                s += std::to_string(m_ani_sum) + '\t';
                s += std::to_string(m_mapq_sum) + '\t';
                s += std::to_string(m_ani_sum/static_cast<double>(m_mapped_reads)) + '\t';
                s += std::to_string(m_mapq_sum/static_cast<double>(m_mapped_reads)) + '\t';
                s += std::to_string(m_strain_level.GetSequenceRangeHandler().CoveredPortion()) + '\t';
                s += std::to_string(m_strain_level.GetSequenceRangeHandler().CoveredPortion()/static_cast<double>(m_gene_length)) + '\t';
                s += std::to_string(alleles.mono) + '\t';
                s += std::to_string(alleles.bi) + '\t';
                s += std::to_string(alleles.tri) + '\t';
                s += std::to_string(alleles.tetra) + '\t';
                s += std::to_string(alleles_filtered.filtered) + '\t';
                s += std::to_string(alleles_filtered.mono) + '\t';
                s += std::to_string(alleles_filtered.bi) + '\t';
                s += std::to_string(alleles_filtered.tri) + '\t';
                s += std::to_string(alleles_filtered.tetra);
                return s;
            }



            bool AddSam(SamEntry const& sam, size_t read_id, double ani=0.0, bool no_strain=true) {
                if (!no_strain) {
                    auto successful = m_strain_level.AddSam(sam, read_id, true);
                    if (!successful) {
                        return false;;
                    }
                }


                // if (extract_snps) {
                //     ExtractSNPs(sam, m_gene_ref->Sequence(), m_snps, 0, 0, read_id);
                // }

                auto const [identity, length] = AlignmentIdentity(sam.m_cigar);
                m_read_identities.emplace_back(static_cast<float>(identity), static_cast<uint32_t>(length));
                m_identity_bases += identity * static_cast<double>(length);
                if (read_id != m_last_read) {
                    m_fragments++;
                    m_last_read = read_id;
                }

                m_mapped_reads++;
                m_mapped_length += length;
                m_total_kmers += (length - 30) * 0.2;
                m_mapq_sum += sam.m_mapq;
                m_ani_sum += ani;

                m_unique_mers += sam.m_uniques;
                m_unique_two_mers += sam.m_uniques_two;

                m_unique_mer_reads += sam.m_uniques > 0;
                m_unique_two_mer_reads += sam.m_uniques_two > 0;

                return true;
            }

            void SetLength(size_t length) {
                m_gene_length = length;
            }

            size_t LongUniques() const {
                return m_unique_mers;
            }

            size_t ReadsWithLongUniques() const {
                return m_unique_mer_reads;
            }

            size_t LongSuperUniques() const {
                return m_unique_two_mers;
            }

            size_t ReadsWithLongSuperUniques() const {
                return m_unique_two_mer_reads;
            }

            size_t TotalKmers() const {
                return m_total_kmers;
            }

            double UniqueRate() const {
                return m_unique_mers / static_cast<double>(m_total_kmers);
            }

            double SuperUniqueRate() const {
                return m_unique_two_mers / static_cast<double>(m_total_kmers);
            }


            double VerticalCoverage() const {
                return static_cast<double>(m_mapped_length)/static_cast<double>(m_gene_length);
            }

//            SNPs& GetSNPs() {
//                return m_snps;
//            }
//
//            const SNPs& GetSNPs() const {
//                return m_snps;
//            }

            std::string ToString() const {
                std::string str = "{";
                str += std::to_string(m_mapped_reads) + ",";
                str += std::to_string(static_cast<double>(m_mapq_sum)/m_mapped_reads) + ",";
                str += std::to_string(m_ani_sum/m_mapped_reads) + "}";


                return str;
            }

            std::string ToString2() const {
                std::string str;
                str += std::to_string(m_mapped_reads) + "\t";
                str += std::to_string(static_cast<double>(m_mapq_sum)/m_mapped_reads) + "\t";
                str += std::to_string(m_ani_sum/m_mapped_reads) + "\t";
//                str += std::to_string(m_snps.size());
                return str;
            }
        };

        class Taxon {
        private:
            size_t m_id;
            std::string m_name;
            size_t m_total_hits = 0;
            size_t m_total_kmers = 0;
            size_t m_unique_mers = 0;
            size_t m_unique_mer_reads = 0;
            size_t m_unique_hits = 0;
            double m_ani_sum = 0;
            size_t m_mapq_sum = 0;
            mutable double m_vcov = -1;  // cached by VerticalCoverage
            mutable double m_low_identity_share = 0;
            double m_depth_identity_margin = 1;  // every read counts towards depth unless set
            mutable std::optional<double> m_model_score;  // cached by TaxonFilterForest::Score
            mutable std::optional<double> m_top_identity;  // cached by TopIdentity
            size_t m_fragments = 0;  // reads with an accepted alignment, a pair counting once
            size_t m_last_read = SIZE_MAX;
            GeneMap m_genes;

            Genome* m_genome;  // the database's genome, shared by all samples
            size_t m_genome_gene_count = 0;

        public:

            Taxon(Genome& genome) : m_genome(&genome), m_genome_gene_count(genome.GeneNum()) {}

            // Drops what is computed from the reads, when a read is added.
            void Changed() {
                m_model_score.reset();
                m_top_identity.reset();
                m_vcov = -1;
            }

            void AddHit(GeneId geneid, GenePos genepos, double ani, bool unique) {
                Changed();
                if (!m_genes.contains(geneid)) {
                    auto& g = m_genome->GetGene(geneid);
                    g.LoadOMP();
                    m_genes.insert( { geneid, profiler::Gene(g) } );
                    m_genes.at(geneid).SetLength(m_genome->GetGene(geneid).Sequence().length());
                }

                m_genes.at(geneid).AddRead(genepos);
                m_total_hits++;
                m_ani_sum += ani;
                m_unique_hits += unique;
            }

            std::vector<size_t> GetAlleles(size_t min_cov=0, size_t min_qual_sum=0) const {
                std::vector<size_t> alleles(5, 0);
                for (auto& [id, _] : GetGenes()) {
                    auto& gene = GetGenes().at(id);
                    gene.GetAlleles(alleles, min_cov, min_qual_sum);
                }
                alleles.resize(5,0);
                return alleles;
            }

            const GeneMap& GetGenes() const {
                return m_genes;
            }

            GeneMap& GetGenes() {
                return m_genes;
            }

            std::vector<uint32_t> SortedGeneIds() const {
                std::vector<uint32_t> ids;
                ids.reserve(m_genes.size());
                for (auto const& [id, _] : m_genes) ids.push_back(id);
                std::sort(ids.begin(), ids.end());
                return ids;
            }

            bool AddSam(GeneId geneid, SamEntry const& sam, double score, bool unique, size_t read_id, bool no_strain=true) {
                Changed();
                bool const new_gene = !m_genes.contains(geneid);
                if (new_gene) {
                    auto& g = m_genome->GetGene(geneid);
                    g.LoadOMP();
                    m_genes.insert( { geneid, profiler::Gene(g) } );
                    m_genes.at(geneid).SetLength(m_genome->GetGene(geneid).Sequence().length());
                }

                bool success = m_genes.at(geneid).AddSam(sam, read_id, score, no_strain);
                if (!success) {
                    // A gene is present only with at least one read.
                    if (new_gene) m_genes.erase(geneid);
                    return false;
                }

                m_unique_mers += sam.m_uniques;
                m_unique_mer_reads += sam.m_uniques > 0;
                m_total_hits++;
                if (read_id != m_last_read) {
                    m_fragments++;
                    m_last_read = read_id;
                }
                m_total_kmers += (sam.m_seq.length() - 30) * 0.2; //TODO: store that info in sam
                m_ani_sum += score;
                m_mapq_sum += sam.m_mapq;
                m_unique_hits += unique;
                return true;
            }


            size_t GetGenomeGeneNumber() const {
                return m_genome->GeneNum();
            }

            const Genome& GetGenome() const {
                return *m_genome;
            }

            Genome& GetGenome() {
                return *m_genome;
            }

            size_t LongUniques() const {
                return std::accumulate(m_genes.begin(), m_genes.end(), size_t{0}, [](auto acc, auto const& pair) { return acc + pair.second.LongUniques(); });
            }

            size_t LongSuperUniques() const {
                return std::accumulate(m_genes.begin(), m_genes.end(), size_t{0}, [](auto acc, auto const& pair) { return acc + pair.second.LongSuperUniques(); });
            }

            size_t GenesWithLongUniques(const size_t threshold=0) const {
                return std::count_if(m_genes.begin(), m_genes.end(), [threshold](auto const& pair) { return pair.second.LongUniques() > threshold; });
            }

            size_t GenesWithLongSuperUniques(const size_t threshold=0) const {
                return std::count_if(m_genes.begin(), m_genes.end(), [threshold](auto const& pair) { return pair.second.LongSuperUniques() > threshold; });
            }

            double GetLongUniqueGeneRate(const size_t threshold=0) const {
                auto lu_genes = GenesWithLongUniques(threshold);
                auto lu_genes_ref = m_genome->GenesWithLongUniques(threshold);
                return lu_genes == 0 || lu_genes_ref == 0 ? 0 : lu_genes/static_cast<double>(lu_genes_ref);
            }

            double GetLongSuperUniqueGeneRate(const size_t threshold=0) const {
                auto lsu_genes = GenesWithLongSuperUniques(threshold);
                auto lsu_genes_ref = m_genome->GenesWithLongSuperUniques(threshold);
                return lsu_genes == 0 || lsu_genes_ref == 0 ? 0 : lsu_genes/static_cast<double>(lsu_genes_ref);
            }

            size_t PresentGenes() const {
                return m_genes.size();
            }

            double Uniqueness() const {
                return m_total_hits == 0 ? 0 : static_cast<double>(m_unique_hits)/static_cast<double>(m_total_hits);
            }

            double GetMeanANI() const {
                return m_total_hits == 0 ? 0 : m_ani_sum/m_total_hits;
            }

            double GetMeanMAPQ() const {
                // Integer division, as the model was trained with.
                return m_total_hits == 0 ? 0 : m_mapq_sum/m_total_hits;
            }

            size_t TotalHits() const {
                return m_total_hits;
            }

            size_t TotalLength() const {
                return std::accumulate(m_genes.begin(), m_genes.end(), size_t{0}, [](size_t acc, std::pair<uint32_t, Gene> const& pair){
                    return acc + pair.second.m_mapped_length;
                });
            }

            double GetGeneVariance(size_t min_gene_hits = 1) const {
                std::vector<double> vec = VerticalCoverageVector();
                std::vector<double> vec_filtered;
                std::copy_if(vec.begin(), vec.end(), std::back_inserter(vec_filtered),
                 [min_gene_hits](double value) { return value >= min_gene_hits; });

                auto var = Variance(vec_filtered);
                return var;
            }

            size_t UniqueHits() const {
                return m_unique_hits;
            }

            void SetId(size_t id) {
                m_id = id;
            }

            void SetName(std::string name) {
                m_name = name;
            }

            std::vector<double> VerticalCoverageVector() const {
                std::vector<double> vcovs;
                for (auto& [geneid, gene] : m_genes) {
                    vcovs.emplace_back(gene.VerticalCoverage());
                }
                return vcovs;
            }



            double Variance(std::vector<double> v) const {
                if (v.size() < 2) return 500.0f;
                double sum = std::accumulate(v.begin(), v.end(), 0.0);
                double mean = sum / v.size();

                std::vector<double> diff(v.size());
                std::transform(v.begin(), v.end(), diff.begin(),
                               [mean](double x){ return x - mean; });
                double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
                double stdev = std::sqrt(sq_sum / v.size());
                return stdev;
            }

            double StandardDeviation(std::vector<double> v) const {
                return std::sqrt(Variance(v));
            }

            double VCovStdDev() const {
                auto v = VerticalCoverageVector();
                double sum = std::accumulate(v.begin(), v.end(), 0.0);
                double mean = sum / v.size();

                std::vector<double> diff(v.size());
                std::transform(v.begin(), v.end(), diff.begin(),
                               [mean](double x){ return x - mean; });
                double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
                double stdev = std::sqrt(sq_sum / v.size());
                return stdev;
            }

            static double Median(std::vector<double> const& v) {
                if (v.empty()) return 0;
                size_t mid_index = v.size()/2;
                return v.size() % 2 == 1 ?
                    v.at(mid_index) :
                    (v.at(mid_index - 1) + v.at(mid_index)) / 2;
            }

            static double SmoothStep(double x) {
                x = std::clamp(x, 0.0, 1.0);
                return x * x * (3 - 2 * x);
            }

            // Depth estimate from two estimators that agree at high coverage. The median depth over
            // genes with reads is robust to outlier genes but biased upwards at low coverage, because
            // genes that drew no read are left out (by about 1/(1-e^-n) for n reads per gene). The
            // aligned bases over the length of all expected (hittable) genes, zeros included, is
            // unbiased at any coverage. The weight of the median rises smoothly with the fraction of
            // expected genes hit (0.80 to 0.95) or with the median depth itself (0.5x to 1x),
            // whichever is higher, so that the estimate changes continuously, without a step.
            static double BlendedDepth(double median_depth, size_t mapped_bases, size_t expected_length,
                                       size_t hit_genes, size_t expected_genes) {
                if (expected_length == 0 || expected_genes == 0) return median_depth;
                double const depth_all_genes = static_cast<double>(mapped_bases) / static_cast<double>(expected_length);
                double const hit_fraction = std::min(1.0, static_cast<double>(hit_genes) / static_cast<double>(expected_genes));
                // A high median depth is trusted only once enough genes are hit: two reads on one
                // short gene give a high median depth, too.
                double const weight = std::max(SmoothStep((hit_fraction - 0.80) / 0.15),
                                               SmoothStep((median_depth - 0.5) / 0.5) * SmoothStep((hit_fraction - 0.25) / 0.25));
                return (1 - weight) * depth_all_genes + weight * median_depth;
            }

            void SetDepthIdentityMargin(double margin) {
                m_depth_identity_margin = margin;
                m_vcov = -1;
                m_model_score.reset();  // the depth is a feature
            }

            // The model score is computed once per taxon (it is asked for by every writer) and dropped
            // whenever the taxon changes.
            std::optional<double> ModelScore() const { return m_model_score; }
            void SetModelScore(double score) const { m_model_score = score; }
            void InvalidateModelScore() { m_model_score.reset(); }

            // Frees the per-read data once the taxon's outputs are computed: every read's identity
            // (after the depth is cached) and, unless `keep_strain_data`, the genes' variants and read
            // ranges. Depth, counters and a cached model score stay valid; features do not, so the
            // taxon must be scored first (see MicrobialProfile::ReleaseReadData).
            void ReleaseReadData(bool keep_strain_data) {
                VerticalCoverage();
                for (auto it = m_genes.begin(); it != m_genes.end(); ++it) {
                    auto& gene = it.value();
                    std::vector<std::pair<float, uint32_t>>{}.swap(gene.m_read_identities);
                    if (!keep_strain_data) gene.GetStrainLevel().Clear();
                }
            }

            // The identity of the taxon's best-matching reads: the 98th percentile, by aligned bases,
            // of its reads' identities (indels count as differences). 0 without reads.
            double TopIdentity() const {
                if (m_top_identity) return *m_top_identity;
                std::vector<std::pair<float, uint32_t>> reads;
                size_t total = 0;
                for (auto const& [id, gene] : m_genes) {
                    for (auto const& read : gene.m_read_identities) {
                        reads.push_back(read);
                        total += read.second;
                    }
                }
                double top = 0;
                if (!reads.empty()) {
                    std::sort(reads.begin(), reads.end());
                    size_t cumulative = 0;
                    top = reads.back().first;
                    for (auto const& [identity, length] : reads) {
                        cumulative += length;
                        if (cumulative >= 0.98 * static_cast<double>(total)) {
                            top = identity;
                            break;
                        }
                    }
                }
                m_top_identity = top;
                return top;
            }

            // The lowest identity of a read that counts towards the taxon's depth: `margin` below
            // TopIdentity. A present species' own reads form this top cluster; reads of relatives
            // (absent from the database, or much more abundant) align at lower identity and would
            // inflate its depth. They still count for detection: the model's features use every read.
            double OwnIdentityThreshold() const {
                if (m_depth_identity_margin >= 1 || PresentGenes() == 0) return 0;
                return TopIdentity() - m_depth_identity_margin;
            }

            // Fragments with an accepted alignment: a read pair counts once, as it is one draw from the
            // genome (TotalHits counts mates).
            size_t Fragments() const {
                return m_fragments;
            }

            // Aligned bases weighted by their read's identity, over all aligned bases: identity with
            // indels as differences (GetMeanANI counts them as matches).
            double BaseIdentity() const {
                double weighted = 0;
                size_t bases = 0;
                for (auto const& [id, gene] : m_genes) {
                    weighted += gene.m_identity_bases;
                    bases += gene.m_mapped_length;
                }
                return bases == 0 ? 0 : weighted / static_cast<double>(bases);
            }

            // `count` per 1000 aligned bases.
            double PerAlignedKb(double count) const {
                size_t const bases = TotalLength();
                return bases == 0 ? 0 : count * 1000.0 / static_cast<double>(bases);
            }

            // Reference bases covered by at least one read.
            size_t CoveredBases() const {
                size_t covered = 0;
                for (auto const& [id, gene] : m_genes) covered += gene.GetStrainLevel().GetSequenceRangeHandler().CoveredPortion();
                return covered;
            }

            // Genes with reads, as a fraction of the hittable genes.
            double HitGeneFraction() const {
                size_t const hittable = m_genome->GeneNum();
                return hittable == 0 ? 0 : static_cast<double>(PresentGenes()) / static_cast<double>(hittable);
            }

            // Genes with reads over the number expected if the fragments fell on the hittable genes in
            // proportion to their length: about 1 for a species that is present, less when the reads
            // gather on a few genes (those shared with a relative, say). It does not depend on depth
            // or the size of the genome.
            double GenePresenceRatio() const {
                auto [lengths, total] = HittableGeneLengths();
                if (total == 0 || m_fragments == 0) return 0;
                double expected = 0;
                for (auto const& [id, length] : lengths) {
                    expected += 1 - std::pow(1 - static_cast<double>(length) / total, static_cast<double>(m_fragments));
                }
                return expected > 0 ? static_cast<double>(PresentGenes()) / expected : 0;
            }

            // Pearson chi-square per degree of freedom of the fragments per hittable gene against their
            // share of its length: about 1 when the reads spread over the genome at random, whatever the
            // depth; large when they gather on some genes.
            double GeneDispersion() const {
                auto [lengths, total] = HittableGeneLengths();
                if (lengths.size() < 2 || total == 0 || m_fragments == 0) return 0;
                double chi_square = 0;
                for (auto const& [id, length] : lengths) {
                    double const expected = static_cast<double>(m_fragments) * static_cast<double>(length) / total;
                    double const observed = m_genes.contains(id) ? static_cast<double>(m_genes.at(id).m_fragments) : 0;
                    chi_square += (observed - expected) * (observed - expected) / expected;
                }
                return chi_square / static_cast<double>(lengths.size() - 1);
            }

            // Coefficient of variation of the depth of the genes with reads (VCovStdDev over their mean).
            double DepthCV() const {
                auto v = VerticalCoverageVector();
                if (v.empty()) return 0;
                double const mean = std::accumulate(v.begin(), v.end(), 0.0) / static_cast<double>(v.size());
                return mean > 0 ? VCovStdDev() / mean : 0;
            }

            // The hittable genes of the genome with their lengths, and the lengths' sum.
            std::pair<std::vector<std::pair<uint32_t, size_t>>, double> HittableGeneLengths() const {
                std::vector<std::pair<uint32_t, size_t>> lengths;
                double total = 0;
                for (auto id : m_genome->GetHittableGenes()) {
                    if (!m_genome->HasGene(id)) continue;
                    size_t const length = m_genome->GetGene(id).GetLength();
                    if (length == 0) continue;
                    lengths.emplace_back(id, length);
                    total += static_cast<double>(length);
                }
                return { lengths, total };
            }

            // Depth of the taxon from its own reads (see OwnIdentityThreshold), estimated by BlendedDepth.
            double VerticalCoverage(bool force=false) const {
                if (m_vcov == -1 || force) {
                    double const min_identity = OwnIdentityThreshold();
                    std::vector<double> vcovs;
                    size_t own_bases = 0, all_bases = 0;
                    for (auto& [geneid, gene] : m_genes) {
                        size_t const bases = gene.MappedLength(min_identity);
                        all_bases += gene.m_mapped_length;
                        if (bases == 0 || gene.m_gene_length == 0) continue;
                        vcovs.emplace_back(static_cast<double>(bases) / static_cast<double>(gene.m_gene_length));
                        own_bases += bases;
                    }
                    m_low_identity_share = all_bases == 0 ? 0 : 1 - static_cast<double>(own_bases) / static_cast<double>(all_bases);
                    std::sort(vcovs.begin(), vcovs.end());

                    if (vcovs.empty()) {
                        m_vcov = 0;
                        return m_vcov;
                    }

                    size_t expected_length = 0;
                    size_t expected_genes = 0;
                    for (auto gene_id : m_genome->GetHittableGenes()) {
                        if (!m_genome->HasGene(gene_id)) continue;
                        expected_length += m_genome->GetGene(gene_id).GetLength();
                        expected_genes++;
                    }
                    m_vcov = BlendedDepth(Median(vcovs), own_bases, expected_length, vcovs.size(), expected_genes);
                }
                return m_vcov;
            }

            // Share of the taxon's aligned bases below OwnIdentityThreshold: reads of relatives, e.g. of
            // a species the database lacks.
            double LowIdentityShare() const {
                VerticalCoverage();
                return m_low_identity_share;
            }

            Gene& GetGene(int geneid) {
                if (!m_genes.contains(geneid)) {
                    std::cerr << geneid << " not in genes " << std::endl;
                    exit(9);
                }
                return m_genes.at(geneid);
            }

            double GetAbundance(double total_vertical_coverage) {
                return VerticalCoverage()/total_vertical_coverage;
            }


            std::string ToString() const {
                std::string str;

                str += std::to_string(m_id);
                str += "\t{ ";
                str += "Genes: " + std::to_string(m_genes.size()) + ", ";
                str += "Total Hits: " + std::to_string(m_total_hits) + ", ";
                str += "MeanMAPQ: " + std::to_string(GetMeanMAPQ()) + ", ";
                str += "VCOV: " + std::to_string(VerticalCoverage()) + ", ";
                str += " }";

                return str;
            }

            std::string ToString() {
                std::string str;

                str += std::to_string(m_id);
                str += "\t{ ";
                str += "Genes: " + std::to_string(m_genes.size()) + ", ";
                str += "Total Hits: " + std::to_string(m_total_hits) + ", ";
                str += "MeanMAPQ: " + std::to_string(GetMeanMAPQ()) + ", ";
                str += "VCOV: " + std::to_string(VerticalCoverage()) + ", ";
                str += " }";

                return str;
            }

            std::string& GetName() {
                return m_name;
            }

            const std::string& GetName() const {
                return m_name;
            }

            std::string ToString(taxonomy::IntTaxonomy& taxonomy) const {
                std::string str;
                auto& genes = m_genome->GetGeneList();
                auto total_gene_length = std::accumulate(genes.begin(), genes.end(), size_t{0}, [](size_t acc, protal::Gene const& gene){
                    return acc + (gene.IsSet() ? gene.GetLength() : 0);
                });

                str += taxonomy.Get(m_id).scientific_name;
                str += "\t" + std::to_string(m_id) + "\t";
                str += "{ ";
                str += "Genes: " + std::to_string(m_genes.size()) + ", ";
                str += "HittableGenes: " + std::to_string(m_genome->GeneNum()) + ", ";
                str += "Total Hits: " + std::to_string(m_total_hits) + ", ";
                str += "Total MGenome: " + std::to_string(total_gene_length) + ", ";
                str += "MeanANI: " + std::to_string(m_ani_sum/(double)m_total_hits) + ", ";
                str += "MeanMAPQ: " + std::to_string(GetMeanMAPQ()) + ", ";
                str += "VCOV: " + std::to_string(VerticalCoverage()) + ", ";
                str += " }";

//                for (auto& [id, gene] : m_genes) {
//                    auto& variant_handler = gene.GetStrainLevel().GetVariantHandler();
//                    str += "\tGene\t" + std::to_string(id) + "\t" + gene.ToString2() + '\n';
//                    for (auto& [key, bin] : variant_handler.GetVariants()) {
//                        str += std::to_string(key) + '\t' + variant_handler.VariantBinToString(bin) + '\n';
////                        str += "\t\t" + snp.ToString() + '\n';
//                    }
//                }

                return str;
            }
        };



        class TaxonFilter {
            double m_min_mean_ani = 0.94;
            double m_min_gene_presence = 0.7;
            size_t m_min_reads = 50;
            size_t m_min_mean_mapq = 5;
            double m_min_uniqueness = 0.1;
        public:
            TaxonFilter() {}
            TaxonFilter(double min_mean_ani, double min_gene_presence, double min_reads, size_t min_mean_mapq) :
                    m_min_mean_ani(min_mean_ani),
                    m_min_gene_presence(min_gene_presence),
                    m_min_reads(min_reads),
                    m_min_mean_mapq(min_mean_mapq) {
            }

            static double ExpectedGenes(size_t marker_genes_count, size_t mapped_reads) {
                return (1 - pow((double) (marker_genes_count - 1)/marker_genes_count, mapped_reads)) * marker_genes_count;
            }

            static double ExpectedGenesUniqueWeighted(Taxon const& taxon) {
                auto const& genome = taxon.GetGenome();
                auto const& genes = genome.GetGeneList();

                std::vector<double> weights;
                weights.reserve(genes.size());

                for (auto const& gene : genes) {
                    if (!gene.IsSet()) continue;
                    if (!genome.IsGeneHittable(gene.GetId())) continue;

                    auto [short_unique, long_unique, long_super_unique, total_kmers] = gene.GetUniqueKmerCounts();
                    (void) short_unique;
                    (void) long_super_unique;
                    (void) total_kmers;

                    // Use only reference long-unique k-mer counts as weight.
                    auto weight = static_cast<double>(long_unique);
                    if (weight > 0.0) {
                        weights.emplace_back(weight);
                    }
                }

                if (weights.empty()) {
                    return ExpectedGenes(taxon.GetGenomeGeneNumber(), taxon.LongUniques());
                }

                auto mapped_unique_hits = taxon.LongUniques();
                auto total_weight = std::accumulate(weights.begin(), weights.end(), 0.0);
                if (total_weight <= 0.0) {
                    return ExpectedGenes(taxon.GetGenomeGeneNumber(), mapped_unique_hits);
                }

                double expected = 0.0;
                for (auto weight : weights) {
                    auto p = weight / total_weight;
                    expected += 1.0 - std::pow(1.0 - p, static_cast<double>(mapped_unique_hits));
                }
                return expected;
            }

            static double ExpectedGenePresence(Taxon const& taxon) {
                return ExpectedGenes(taxon.GetGenomeGeneNumber(), taxon.TotalHits());
            }

            static double ExpectedGenePresenceRatio(Taxon const& taxon) {
                return static_cast<double>(taxon.PresentGenes()) / ExpectedGenes(taxon.GetGenomeGeneNumber(), taxon.TotalHits());
            }

            static double ExpectedGenePresenceUniqueWeighted(Taxon const& taxon) {
                return ExpectedGenesUniqueWeighted(taxon);
            }

            static size_t PresentGenesUniqueWeighted(Taxon const& taxon) {
                auto const& genes = taxon.GetGenes();
                return std::count_if(genes.begin(), genes.end(), [](auto const& pair) {
                    return pair.second.LongUniques() > 0;
                });
            }

            static double ExpectedGenePresenceRatioUniqueWeighted(Taxon const& taxon) {
                auto expected = ExpectedGenePresenceUniqueWeighted(taxon);
                if (expected <= 0.0) return 0.0;
                return static_cast<double>(PresentGenesUniqueWeighted(taxon)) / expected;
            }

            bool Formula1(Taxon const& taxon) const {
                bool pass = // Currently mix of expected gene presence, min ani and min gene presence
                        taxon.GetMeanANI() >= m_min_mean_ani &&
                        taxon.TotalHits() >= m_min_reads &&
                        ExpectedGenePresenceRatio(taxon) >= m_min_gene_presence &&
                        taxon.GetMeanMAPQ() >= m_min_mean_mapq;
                return pass;
            }

            void PrintViolatingValues(Taxon const& taxon) const {
                if (taxon.GetMeanANI() < m_min_mean_ani) {
                    std::cout << "MeanANI:      " << taxon.GetMeanANI() << " < " << m_min_mean_ani << " (min)" << std::endl;
                }
                if (taxon.TotalHits() < m_min_reads) {
                    std::cout << "TotalHits:    " << taxon.TotalHits() << " < " << m_min_reads << " (min)" << std::endl;
                }
                if (ExpectedGenePresenceRatio(taxon) < m_min_gene_presence) {
                    std::cout << "GenePresence: " << ExpectedGenePresenceRatio(taxon) << " < " << m_min_gene_presence << " (min)" << std::endl;
                }
                if (taxon.GetMeanMAPQ() < m_min_mean_mapq) {
                    std::cout << "MeanMAPQ: " << taxon.GetMeanMAPQ() << " < " << m_min_mean_mapq << " (min)" << std::endl;
                }
            }

            bool Pass(Taxon const& taxon) const {
                return Formula1(taxon);
            }
        };

        class HittableGeneMask {
        private:
            tsl::sparse_map<Profiler::TaxonID, tsl::sparse_set<Profiler::GeneID>> m_hittable_map;

        public:
            void Load(std::string const& file) {
                std::ifstream is(file, std::ios::in);

                std::vector<std::string> tokens;
                std::string line;
                while (std::getline(is, line)) {
                    Utils::split(tokens, line, "\t");
                    auto taxid = std::stoull(tokens[0]);
                    auto gene_str = tokens[2];

                    Utils::split(tokens, gene_str, ",");
                    for (auto const& g : tokens) {
                        auto gid = std::stoul(g);
                        m_hittable_map[taxid].insert(gid);
                    }
                }
            }

            size_t HittableGenes(Profiler::TaxonID taxid) {
                return m_hittable_map.at(taxid).size();
            }

            bool IsHittable(Profiler::TaxonID taxid, Profiler::GeneID geneid) {
                return m_hittable_map[taxid].contains(geneid);
            }
        };

        // The features of a taxon, in the column order of the training dump (<profile>.truth_annotated).
        // The random forest is given them from here, and the dump is written from here, so that a
        // model is trained on exactly the quantities it is later scored with.
        using TaxonFeatureList = std::vector<std::pair<std::string, double>>;

        inline TaxonFeatureList TaxonFeatures(Taxon const& taxon) {
            auto a = taxon.GetAlleles();
            auto af = taxon.GetAlleles(2, 60);
            double const a_sum = std::accumulate(a.begin(), a.end(), 0.0);
            double const af_sum = std::accumulate(af.begin(), af.end(), 0.0);
            auto rate = [](double part, double whole) { return part == 0 || whole == 0 ? 0.0 : part / whole; };
            auto [su, lu, lsu, all] = taxon.GetGenome().GetUniqueKmerCounts();

            TaxonFeatureList f;
            f.reserve(80);
            f.emplace_back("present_genes", taxon.PresentGenes());
            f.emplace_back("total_hits", taxon.TotalHits());
            f.emplace_back("unique_hits", taxon.UniqueHits());
            f.emplace_back("mean_ani", taxon.GetMeanANI());
            f.emplace_back("expected_gene_presence", TaxonFilter::ExpectedGenePresence(taxon));
            f.emplace_back("expected_gene_presence_ratio", TaxonFilter::ExpectedGenePresenceRatio(taxon));
            f.emplace_back("expected_gene_presence_unique_weighted", TaxonFilter::ExpectedGenePresenceUniqueWeighted(taxon));
            f.emplace_back("expected_gene_presence_ratio_unique_weighted", TaxonFilter::ExpectedGenePresenceRatioUniqueWeighted(taxon));
            f.emplace_back("uniqueness", taxon.Uniqueness());
            f.emplace_back("mean_mapq", taxon.GetMeanMAPQ());
            f.emplace_back("variance1", taxon.GetGeneVariance(1));
            f.emplace_back("variance2", taxon.GetGeneVariance(5));
            // Variant positions by number of alleles (all, and with >= 2 observations and a quality
            // sum >= 60), and as fractions of all variant positions.
            for (size_t i = 0; i < 5; i++) f.emplace_back("A" + std::to_string(i), a[i]);
            for (size_t i = 0; i < 5; i++) f.emplace_back("AF" + std::to_string(i), af[i]);
            for (size_t i = 0; i < 5; i++) f.emplace_back("RAF" + std::to_string(i), rate(af[i], af_sum));
            for (size_t i = 0; i < 5; i++) f.emplace_back("RA" + std::to_string(i), rate(a[i], a_sum));
            f.emplace_back("stddev", taxon.VCovStdDev());
            f.emplace_back("hittable", taxon.GetGenomeGeneNumber());
            f.emplace_back("lu", taxon.LongUniques());
            f.emplace_back("lu_genes", taxon.GenesWithLongUniques());
            f.emplace_back("lsu", taxon.LongSuperUniques());
            f.emplace_back("lsu_genes", taxon.GenesWithLongSuperUniques());
            // Unique k-mers of the species' reference, from the index.
            f.emplace_back("su_genome", su);
            f.emplace_back("lu_genome", lu);
            f.emplace_back("lsu_genome", lsu);
            f.emplace_back("total_genome", all);
            f.emplace_back("su_rate_ref", rate(su, all));
            f.emplace_back("lu_rate_ref", rate(lu, all));
            f.emplace_back("lsu_rate_ref", rate(lsu, all));
            f.emplace_back("lu_gene_rate", taxon.GetLongUniqueGeneRate());
            f.emplace_back("lsu_gene_rate", taxon.GetLongSuperUniqueGeneRate());
            f.emplace_back("lu_gene_rate2", taxon.GetLongUniqueGeneRate(1));
            f.emplace_back("lsu_gene_rate2", taxon.GetLongSuperUniqueGeneRate(1));
            f.emplace_back("lu_gene_rate3", taxon.GetLongUniqueGeneRate(5));
            f.emplace_back("lsu_gene_rate3", taxon.GetLongSuperUniqueGeneRate(5));
            f.emplace_back("lsu_per_read", rate(taxon.LongSuperUniques(), taxon.TotalHits()));
            f.emplace_back("lu_per_read", rate(taxon.LongUniques(), taxon.TotalHits()));

            // Features for a model that holds across databases, depths and libraries: fractions and
            // ratios instead of counts of genes and k-mers (archaea have fewer marker genes and
            // k-mers than bacteria), rates per aligned kb instead of per read (reads differ in length),
            // fragments instead of mates, and identity with indels as differences. The model shipped
            // with protal does not use them; they are in the training dump for the next one.
            size_t const covered = taxon.CoveredBases();
            auto per_covered_kb = [covered](double count) { return covered == 0 ? 0.0 : count * 1000.0 / static_cast<double>(covered); };
            f.emplace_back("fragments", taxon.Fragments());
            f.emplace_back("depth", taxon.VerticalCoverage());
            f.emplace_back("hit_gene_fraction", taxon.HitGeneFraction());
            f.emplace_back("gene_presence_ratio", taxon.GenePresenceRatio());
            f.emplace_back("gene_dispersion", taxon.GeneDispersion());
            f.emplace_back("depth_cv", taxon.DepthCV());
            f.emplace_back("identity", taxon.BaseIdentity());
            f.emplace_back("top_identity", taxon.TopIdentity());
            f.emplace_back("low_identity_share", taxon.LowIdentityShare());
            f.emplace_back("lu_per_kb", taxon.PerAlignedKb(taxon.LongUniques()));
            f.emplace_back("lsu_per_kb", taxon.PerAlignedKb(taxon.LongSuperUniques()));
            f.emplace_back("variant_sites_per_kb", per_covered_kb(af[1] + af[2] + af[3] + af[4]));
            f.emplace_back("multiallelic_sites_per_kb", per_covered_kb(af[2] + af[3] + af[4]));
            return f;
        }

        // A feature value as the model and the dump get it: the shortest text that reads back as
        // the same double (fixed six decimals turned small rates into 0).
        inline std::string FeatureString(double value) {
            char buffer[64];
            auto [end, ec] = std::to_chars(buffer, buffer + sizeof(buffer), value);
            return std::string(buffer, end);
        }

        class TaxonFilterForest {
            cpmml::Model m_model{};
            double m_knob = 0.5;

            mutable std::unordered_map<std::string, std::string> m_sample;


        public:

            TaxonFilterForest(const std::string& path, double knob = 0.5) : m_model(path), m_knob(knob) {}

            // A model already parsed, e.g. with cpmml::Model::from_string from a single-file database.
            TaxonFilterForest(cpmml::Model model, double knob) : m_model(std::move(model)), m_knob(knob) {}

            TaxonFilterForest(const TaxonFilterForest& other) :
                    m_model(other.m_model), m_knob(other.m_knob), m_sample() {
            }

            TaxonFilterForest(const TaxonFilterForest&& other) :
                    m_model(other.m_model),
                    m_knob(other.m_knob),
                    m_sample(other.m_sample) {
            }

            // Returns the model's probability for TRUE (0–1).
            double Score(Taxon const& taxon) const {
                if (auto cached = taxon.ModelScore()) return *cached;
                double const score = ScoreFeatures(taxon);
                taxon.SetModelScore(score);
                return score;
            }

            double ScoreFeatures(Taxon const& taxon) const {
                m_sample.clear();
                for (auto const& [name, value] : TaxonFeatures(taxon)) m_sample[name] = FeatureString(value);

                auto dist = m_model.score(m_sample).distribution();
                auto it = dist.find("TRUE");
                if (it != dist.end() && std::isfinite(it->second)) {
                    return it->second;
                }
                // Fallback: model lacks probability support (no ScoreDistributions or
                // all record counts are zero). Treat predict() output as hard 0/1.
                return m_model.predict(m_sample) == "TRUE" ? 1.0 : 0.0;
            }

            bool Pass(Taxon const& taxon) const {
                return Score(taxon) >= m_knob;
            }

            double GetKnob() const { return m_knob; }
        };

        // A placeholder model (scripts/placeholder_models.py) fills a read type's slot in the database
        // until a trained model replaces it: it scores every taxon 0, and protal warns when it loads one.
        inline constexpr std::string_view kPlaceholderModelMarker = "<Annotation>protal:placeholder</Annotation>";

        inline bool IsPlaceholderModel(std::string const& xml) {
            auto const header_end = xml.find("</Header>");
            return header_end != std::string::npos && xml.rfind(kPlaceholderModelMarker, header_end) != std::string::npos;
        }

        // Why protal cannot use the PMML model `model` parsed from `xml`, or an empty string. The
        // model must take its inputs from TaxonFeatures (a missing one would stop the run after the
        // alignment), predict the label TRUE (its probability is the taxon's score; another label
        // scores every taxon 0), and score a taxon.
        inline std::string ModelContractProblemInXml(TaxonFilterForest const& model, std::string const& xml) {
            auto attribute = [](std::string_view tag, std::string const& name) -> std::string {
                auto const key = " " + name + "=\"";
                auto start = tag.find(key);
                if (start == std::string_view::npos) return {};
                start += key.size();
                auto end = tag.find('"', start);
                return end == std::string_view::npos ? std::string{} : std::string(tag.substr(start, end - start));
            };
            auto tags = [&xml](std::string const& element, size_t from = 0, size_t to = std::string::npos) {
                std::vector<std::string_view> found;
                std::string const open = "<" + element + " ";
                for (size_t pos = xml.find(open, from); pos != std::string::npos && pos < to; pos = xml.find(open, pos + 1)) {
                    auto end = xml.find('>', pos);
                    if (end == std::string::npos) break;
                    found.emplace_back(std::string_view(xml).substr(pos, end - pos + 1));
                }
                return found;
            };

            Genome no_genome(0);
            std::set<std::string> provided;
            for (auto const& [name, _] : TaxonFeatures(Taxon(no_genome))) provided.insert(name);

            std::string target;
            std::vector<std::string> missing;
            for (auto tag : tags("MiningField")) {
                auto const name = attribute(tag, "name");
                auto usage = attribute(tag, "usageType");
                if (usage.empty()) usage = "active";
                if (usage == "predicted" || usage == "target") target = name;
                else if (usage == "active" && !provided.contains(name)) missing.push_back(name);
            }
            if (!missing.empty()) {
                std::string list;
                for (auto const& name : missing) list += (list.empty() ? "" : ", ") + name;
                return "it needs " + std::to_string(missing.size()) + " input(s) protal does not compute: " + list;
            }
            if (target.empty()) return "it names no predicted field";

            auto const data_field = xml.find("<DataField name=\"" + target + "\"");
            auto const data_end = data_field == std::string::npos ? std::string::npos : xml.find("</DataField>", data_field);
            bool has_true = false;
            if (data_end != std::string::npos) {
                for (auto tag : tags("Value", data_field, data_end)) has_true |= attribute(tag, "value") == "TRUE";
            }
            if (!has_true) return "its predicted field '" + target + "' has no value TRUE (protal reports the probability of TRUE)";

            try {
                double const score = model.ScoreFeatures(Taxon(no_genome));
                if (!std::isfinite(score)) return "it gives a taxon a score that is not a number";
            } catch (std::exception const& e) {
                return std::string("scoring a taxon fails: ") + e.what();
            }
            return {};
        }

        // ModelContractProblemInXml for the model file at `path`.
        inline std::string ModelContractProblem(TaxonFilterForest const& model, std::string const& path) {
            std::ifstream is(path);
            if (!is) return "cannot read " + path;
            std::string const xml((std::istreambuf_iterator<char>(is)), std::istreambuf_iterator<char>());
            return ModelContractProblemInXml(model, xml);
        }

        class MicrobialProfile {
        public:
            Benchmark bm_add_sam{"Add Sam profile"};
            MicrobialProfile(GenomeLoader& genome_loader) : m_genome_loader(genome_loader) {}

            // See Taxon::OwnIdentityThreshold; 1 or more lets every read count towards depth.
            void SetDepthIdentityMargin(double margin) {
                m_depth_identity_margin = margin;
                for (auto& [id, _] : m_taxa) m_taxa.at(id).SetDepthIdentityMargin(margin);
            }

            // Frees the per-read data of every taxon once the profile's outputs are written (see
            // Taxon::ReleaseReadData), after scoring it. The strain stage reads the variants and read
            // ranges of taxa that pass `filter` only, so only theirs are kept, and only if
            // `keep_strain_data`.
            void ReleaseReadData(TaxonFilterObj const& filter, bool keep_strain_data) {
                for (auto it = m_taxa.begin(); it != m_taxa.end(); ++it) {
                    auto& taxon = it.value();
                    bool const pass = filter.Pass(taxon);  // caches the score
                    taxon.ReleaseReadData(keep_strain_data && pass);
                }
            }

            void AddRead(InternalReadAlignment const& ira, bool unique=true) {
                if (!m_taxa.contains(ira.taxid)) {
                    auto& genome = m_genome_loader.GetGenome(ira.taxid);
                    m_taxa.insert( { ira.taxid, Taxon(genome) } );
                    m_taxa.at(ira.taxid).SetName(std::to_string(ira.taxid));
                }
                auto& taxon = m_taxa.at(ira.taxid);
                taxon.AddHit(ira.geneid, ira.genepos, ira.alignment_ani, unique);
            }

            bool AddSam(int taxid, int geneid, SamEntry const& sam, double score, bool unique=true, int read_id=0, bool no_strain=true) {
                // A record on a gene this database does not have, or reaching past the gene's end (a
                // SAM aligned against another database), is rejected rather than read out of bounds.
                if (!m_genome_loader.HasGene(taxid, geneid) ||
                    sam.m_pos - 1 + AlignmentLengthRef(sam.m_cigar) > m_genome_loader.GeneLength(taxid, geneid)) {
                    return false;
                }
                // Reads on genes without unique k-mers are ignored: they are not evidence of the taxon.
                if (!m_genome_loader.GetGenome(taxid).IsGeneHittable(geneid)) {
                    return true;
                }

                if (!m_taxa.contains(taxid)) {
                    auto &genome = m_genome_loader.GetGenome(taxid);

                    m_taxa.insert( { taxid, Taxon(genome) } );
                    m_taxa.at(taxid).SetId(taxid);
                    m_taxa.at(taxid).SetDepthIdentityMargin(m_depth_identity_margin);
                    m_taxa.at(taxid).SetName(std::to_string(taxid));
                }

                // Debug
                if (!m_taxa.contains(taxid)) {
                    std::cout << "Error: " << taxid << std::endl;
                    std::cout << sam.ToString() << std::endl;
                }
                auto& taxon = m_taxa.at(taxid);
                unique = sam.m_mapq > 20;
                bm_add_sam.Start();
                bool success = taxon.AddSam(geneid, sam, score, unique, read_id, no_strain);
                bm_add_sam.Stop();
                // A taxon exists only with at least one read (its means divide by the read count).
                if (!success && taxon.TotalHits() == 0) m_taxa.erase(taxid);
                return success;
            }

            void SetName(std::string name) {
                m_name = name;
            }

            const std::string& GetName() const {
                return m_name;
            }

            std::string& GetName() {
                return m_name;
            }

            void PostProcessSNPs(size_t min_observations=2, size_t min_observations_fwdrev=2, double min_frequency=0.0, size_t min_avg_quality=15, size_t min_phred_sum=0, bool require_strand=false) {
                for (auto& [tid, _] : m_taxa) {
                    auto& taxon = m_taxa.at(tid);
                    for (auto& [gid, __] : taxon.GetGenes()) {
                        auto& gene = taxon.GetGenes().at(gid);
                        auto& strain_handler = gene.GetStrainLevel();
                        strain_handler.PostProcess(min_observations, min_observations_fwdrev, min_frequency, min_avg_quality, min_phred_sum, require_strand);
                    }
                    taxon.InvalidateModelScore();
                }
            }

            std::unordered_set<size_t> GetKeySet(std::optional<TaxonFilter> filter={}) {
                std::unordered_set<size_t> keys;
                for (auto& [id, taxon] : m_taxa) {
                    if (filter.has_value() && !filter->Pass(taxon)) continue;
                    keys.insert(id);
                }
                return keys;
            }

            std::string ToString(taxonomy::IntTaxonomy& taxonomy, std::optional<TaxonFilterObj> filter={}) {
                std::string str;

                str += m_name + '\t';
                str += std::to_string(m_taxa.size());
                if (filter.has_value()) {
                    str += " (Filtered)";
                }
                str += '\n';
                for (auto const& [id, _] : m_taxa) {
                    auto& taxon = m_taxa.at(id);
                    if (filter.has_value() && !filter->Pass(taxon)) continue;
                    taxon.VerticalCoverage();
                    str += taxon.ToString(taxonomy) + '\n';
                }

                return str;
            }



            std::string ToString() {
                std::string str;

                str += m_name + '\t';
                str += std::to_string(m_taxa.size());
                for (auto const& [id, _] : m_taxa) {
                    str += m_taxa.at(id).ToString() + '\n';
                }

                return str;
            }

            TaxonMap& GetTaxa() {
                return m_taxa;
            }

            const TaxonMap& GetTaxa() const {
                return m_taxa;
            }

            // Writes the training data for the random forest: per taxon the truth, the model's call,
            // and the features exactly as the model is given them (TaxonFeatures).
            void AnnotateWithTruth(TruthSet const& set, TaxonFilterObj& filter, std::string& output, taxonomy::IntTaxonomy& taxonomy) {
                std::ofstream os(output, std::ios::out);

                Genome no_genome(0);
                os << "truth\tprediction\tprobability\ttaxon\ttaxon_name";
                for (auto const& [name, _] : TaxonFeatures(Taxon(no_genome))) os << '\t' << name;
                os << '\n';

                for (auto key : SortedTaxa()) {
                    auto const& taxon = m_taxa.at(key);
                    double probability = filter.Score(taxon);
                    os << set.contains(key) << '\t'
                       << (probability >= filter.GetKnob()) << '\t'
                       << FeatureString(probability) << '\t'
                       << key << '\t'
                       << taxonomy.Get(key).scientific_name;
                    for (auto const& [_, value] : TaxonFeatures(taxon)) os << '\t' << FeatureString(value);
                    os << '\n';
                }

                os.close();
                if (os.fail()) RunStatus::Get().Fail("Writing the truth annotation failed: " + output);
            }

            // Taxon ids in ascending order: outputs list taxa (and genes) independently of the order of
            // the reads in the SAM.
            std::vector<uint32_t> SortedTaxa() const {
                std::vector<uint32_t> ids;
                ids.reserve(m_taxa.size());
                for (auto const& [id, _] : m_taxa) ids.push_back(id);
                std::sort(ids.begin(), ids.end());
                return ids;
            }

            // Summed depth of the taxa that pass the model: abundances are fractions of it.
            double PassingDepth(TaxonFilterObj const& filter) {
                double total = 0;
                for (auto const& [id, _] : m_taxa) {
                    auto& taxon = m_taxa.at(id);
                    if (filter.Pass(taxon)) total += taxon.VerticalCoverage();
                }
                return total;
            }

            // .genes.log: one line per gene hit, with its taxon's call. TaxAbundance is 0 for taxa the model
            // rejects; MAPQ and ANI are means over the gene's reads.
            void WriteGeneProfile(taxonomy::IntTaxonomy& taxonomy, TaxonFilterObj const& filter, std::ostream* os) {
                *os << "Predicted\tProbability\tTaxID\tLineage\tTaxVCOV\tTaxAbundance\tGeneID\tGeneRefLength\t"
                    << "TotalReads\tTotalMappedLength\tMAPQ\tUniqueMers\tUniqueTwoMers\tUniqueMerReads\tUniqueTwoMerReads\tANI\t"
                    << "VCov\tVCovExp\tHCovExp\tHCovObs\tHCovObsRel\tConsistency\n";

                double const total_vcov = PassingDepth(filter);

                for (auto tax_id : SortedTaxa()) {
                    auto& taxon = m_taxa.at(tax_id);
                    double probability = filter.Score(taxon);
                    bool prediction = probability >= filter.GetKnob();

                    auto predicted_vcov = taxon.VerticalCoverage();
                    double const predicted_abundance = prediction && total_vcov > 0 ? predicted_vcov / total_vcov : 0;

                    for (auto gene_id : taxon.SortedGeneIds()) {
                        auto& gene = taxon.GetGenes().at(gene_id);
                        auto vcov = gene.VerticalCoverage();
                        auto expected_vcov = static_cast<double>(gene.m_mapped_length) / static_cast<double>(gene.m_gene_length);
                        auto expected_hcov = 1 - exp(-expected_vcov);
                        auto hcov_obs = gene.GetStrainLevel().GetSequenceRangeHandler().CoveredPortion();
                        auto hcov_obs_rel = static_cast<double>(hcov_obs) / static_cast<double>(gene.m_gene_length);
                        auto coverage_consistency_ratio = expected_hcov / hcov_obs_rel;
                        double const reads = static_cast<double>(gene.m_mapped_reads);

                        *os
                            << prediction << '\t'
                            << probability << '\t'
                            << tax_id << '\t'
                            << taxonomy.LineageStr(tax_id) << '\t'
                            << predicted_vcov << '\t'
                            << predicted_abundance << '\t'
                            << gene_id << '\t'
                            << gene.m_gene_length << '\t'
                            << gene.m_mapped_reads << '\t'
                            << gene.m_mapped_length << '\t'
                            << (reads > 0 ? gene.m_mapq_sum / reads : 0) << '\t'
                            << gene.m_unique_mers << '\t'
                            << gene.m_unique_two_mers << '\t'
                            << gene.m_unique_mer_reads << '\t'
                            << gene.m_unique_two_mer_reads << '\t'
                            << (reads > 0 ? gene.m_ani_sum / reads : 0) << '\t'
                            << vcov << '\t'
                            << expected_vcov << '\t'
                            << expected_hcov << '\t'
                            << hcov_obs << '\t'
                            << hcov_obs_rel << '\t'
                            << coverage_consistency_ratio
                            << '\n';
                    }
                }
            }

            // The profile (taxa that pass the model), .profile.log (every taxon with its call, features
            // and per-gene coverage; one column per gene id of the database, so every file has the same
            // columns) and .gene.log (statistics per gene hit).
            void WriteSparseProfile(taxonomy::IntTaxonomy& taxonomy, TaxonFilterObj const& filter, std::ostream &os_filtered=std::cout, std::ostream* os_total=nullptr, std::ostream* os_genes=nullptr) {
                bool one_pass = false;

                std::string gene_covs_str = "";
                std::string gene_cov_ratios_str = "";

                size_t max_gene_id = 0;
                for (auto& [key, genome] : m_genome_loader.GetGenomeMap()) {
                    max_gene_id = std::max(max_gene_id, genome.GetGeneList().size());
                }
                std::vector<size_t> gene_covs(max_gene_id + 1, 0);
                std::vector<double> gene_cov_ratios(max_gene_id + 1, 0.0);

                if (os_total) {
                    *os_total << "Predicted\tProbability\tRepGenome\tLineage\tAbundance\tVCovStdDev\tGeneVariance\tGeneVariance5\t"
                              << "Name\tTaxID\tSummary\tVCov\tLowIdentityShare\tMeanGeneCov\tMeanGeneCovRatio";
                    for (size_t id = 0; id <= max_gene_id; id++) *os_total << "\tGeneCov" << id;
                    for (size_t id = 0; id <= max_gene_id; id++) *os_total << "\tGeneCovRatio" << id;
                    *os_total << '\n';
                }
                if (os_genes) {
                    *os_genes << "Sample\tTaxID\tLineage\tName\tGeneID\tGene\tReads\tCoverageSum\tGeneLength\tANISum\tMAPQSum\t"
                              << "MeanANI\tMeanMAPQ\tCoveredBases\tCoveredFraction\tMono\tBi\tTri\tTetra\t"
                              << "FilteredNoAllele\tFilteredMono\tFilteredBi\tFilteredTri\tFilteredTetra\n";
                }

                double const total_vcov = PassingDepth(filter);

                for (auto key : SortedTaxa()) {
                    // this is necessary as taxon cannot be constant
                    auto& taxon = m_taxa.at(key);
                    double probability = filter.Score(taxon);
                    bool prediction = probability >= filter.GetKnob();
                    one_pass |= prediction;
                    double const vcov = taxon.VerticalCoverage();
                    double const abundance = prediction && total_vcov > 0 ? vcov / total_vcov : 0;

                    gene_covs_str.clear();
                    gene_cov_ratios_str.clear();
                    std::fill(gene_covs.begin(), gene_covs.end(), 0);
                    std::fill(gene_cov_ratios.begin(), gene_cov_ratios.end(), 0.0);

                    for (auto id : taxon.SortedGeneIds()) {
                        auto& gene = taxon.GetGenes().at(id);
                        size_t cov = gene.GetStrainLevel().GetSequenceRangeHandler().CoveredPortion();
                        double ratio = static_cast<double>(cov) / gene.m_gene_length;
                        gene_covs[id] = cov;
                        gene_cov_ratios[id] = ratio;

                        if (os_genes) {
                            *os_genes << m_name << '\t';
                            *os_genes << key << '\t';
                            *os_genes << taxonomy.LineageStr(key) << '\t';
                            *os_genes << taxon.GetName() << '\t';
                            *os_genes << id << '\t';
                            *os_genes << "Gene" + std::to_string(id) << '\t';
                            *os_genes << gene.GetStatisticsString() << '\n';
                        }
                    }

                    gene_covs_str = std::accumulate(gene_covs.begin(), gene_covs.end(), std::string{}, [](std::string acc, size_t val) {
                        return(acc + "\t" + std::to_string(val));
                    });
                    gene_cov_ratios_str = std::accumulate(gene_cov_ratios.begin(), gene_cov_ratios.end(), std::string{}, [](std::string acc, double val) {
                        return(acc + "\t" + std::to_string(val));
                    });

                    double mean_gene_covs = std::accumulate(gene_covs.begin(), gene_covs.end(), size_t{0}) /
                                            static_cast<double>(taxon.GetGenes().size());
                    double mean_gene_cov_ratios = std::accumulate(gene_cov_ratios.begin(), gene_cov_ratios.end(), 0.0) /
                                            static_cast<double>(taxon.GetGenes().size());

                    auto node = taxonomy.Get(key);

                    if (prediction) {
                        os_filtered << node.rep_genome << '\t' << taxonomy.LineageStr(key) << '\t' << abundance << std::endl;
                    }

                    if (os_total) {
                        *os_total << (prediction ? "1" : "0") << "\t" << probability << "\t" << node.rep_genome << '\t' << taxonomy.LineageStr(key) << '\t' << abundance;
                        *os_total << '\t' << taxon.VCovStdDev();
                        *os_total << '\t' << taxon.GetGeneVariance();
                        *os_total << '\t' << taxon.GetGeneVariance(5);
                        *os_total << '\t' << taxon.ToString(taxonomy);
                        *os_total << '\t' << vcov << '\t' << taxon.LowIdentityShare();
                        *os_total << '\t' << mean_gene_covs << '\t' << mean_gene_cov_ratios << gene_covs_str;
                        *os_total << gene_cov_ratios_str << std::endl;
                    }
                }


                if (!one_pass) {
                    std::cout << "No taxon passes the model in sample " << m_name << std::endl;
                }
            }


        private:
            std::string m_name;
            mutable TaxonMap m_taxa;
            GenomeLoader &m_genome_loader;
            double m_depth_identity_margin = 1;
        };


        using SamPairList = std::vector<std::vector<AlignmentPair>>;
        using SamPairs = std::vector<AlignmentPair>;
        using SamPairList_ptr = std::vector<std::vector<AlignmentPair>*>;
        using SamPairs_ptr = std::vector<AlignmentPair*>;

        template<typename ScoringSystem=ScoreAlignments>
        class Profiler {

        public:
            Profiler(GenomeLoader& genome_loader) :
                    m_genome_loader(genome_loader) {}

        private:
            GenomeLoader& m_genome_loader;


            SamPairs_ptr m_pairs_unique_ptr;
            SamPairList_ptr m_pairs_nonunique_ptr;

            ScoringSystem m_score;

            CigarInfo m_info1;
            CigarInfo m_info2;

            size_t m_min_alignment_length = 50;
            size_t m_min_mapq = 4;
            double m_depth_identity_margin = 1;
            size_t m_reads = 0;
            size_t m_rejected_reads = 0;

            // Variants are recorded with or without --no_strains: the model's allele features come
            // from them, so a profile must not depend on whether strain MSAs are written.
            static constexpr bool kNoStrain = false;
            std::string m_id = {};

        public:

            SamPairs m_pairs_unique;
            SamPairList m_pairs_nonunique;
            SamPairs m_pairs_nonunique_best;
            Benchmark m_post_process_bm{"Post-processing"};

            // See Taxon::OwnIdentityThreshold.
            void SetDepthIdentityMargin(double margin) {
                m_depth_identity_margin = margin;
            }

            static std::pair<double, int> ScorePairedAlignment(OptIRA const& a, OptIRA const& b, double divide_penalty=2.2, double alone_penalty=2) {
                int score = 0;
                double ani = 0;
                if (a.has_value() && b.has_value()) {
                    score = (a.value().alignment_score + b.value().alignment_score) / divide_penalty;
                    ani = (a.value().alignment_ani + b.value().alignment_ani) / 2;
                } else if (a.has_value()) {
                    score = a.value().alignment_score - alone_penalty;
                    ani = a.value().alignment_ani;
                } else if (b.has_value()) {
                    score = b.value().alignment_score - alone_penalty;
                    ani = b.value().alignment_ani;
                }
                return { ani, score };
            }

            static std::pair<double, int> ScorePairedAlignment(OptIRAPair const& pair, double divide_penalty=2.2, double alone_penalty=2) {
                return ScorePairedAlignment(pair.first, pair.second);
            }

            static bool SameReadId(OptIRA const &ira1, OptIRA const &ira2, OptIRAPair const& pair) {
                assert(pair.first.has_value() || pair.second.has_value());
                assert(ira1.has_value() || ira2.has_value());

                if (ira1.has_value() && ira2.has_value() && ira1.value().readid != ira2.value().readid) {
                    exit(12);
                }

                auto& rid = ira1.has_value() ? ira1.value().readid : ira2.value().readid;
                auto& rid2 = pair.first.has_value() ? pair.first.value().readid : pair.second.value().readid;

                return rid == rid2;
            }

            // Both mates carry the same QNAME (PairQName), so the whole name identifies the read; the
            // uniqueness test in FromSam compares whole names too.
            static bool SameRead(AlignmentPair& pair, AlignmentPair& other) {
                return pair.Any().m_qname == other.Any().m_qname;
            }

//            void Process(std::optional<InternalReadAlignment> const &ira,
//                         std::optional<InternalReadAlignment> const &ira2,
//                         bool unique) {
//                assert(ira.has_value() || ira2.has_value());
//                if (unique) {
//                    ProcessUnique(ira, ira2);
//                } else {
////                    ProcessNonUnique(ira, ira2);
//                }
//            }

//            void Process(size_t read_id, SamEntry &sam1, SamEntry &sam2, bool has_sam1, bool has_sam2, bool unique) {
//                CigarInfo info;
//                size_t min_alignment_length = 70;
//                if (has_sam1) {
//                    CompressedCigarInfo(sam1.m_cigar, info);
//                    if (info.clipped_alignment_length < min_alignment_length) has_sam1 = false;
//                }
//                if (has_sam2) {
//                    CompressedCigarInfo(sam2.m_cigar, info);
//                    if (info.clipped_alignment_length < min_alignment_length) has_sam2 = false;
//                }
//                if (!has_sam1 && !has_sam2) return;
//                OptIRA oira1 = has_sam1 ? OptIRA(InternalReadAlignment(read_id, sam1)) : OptIRA();
//                OptIRA oira2 = has_sam2 ? OptIRA(InternalReadAlignment(read_id, sam2)) : OptIRA();
//                Process(oira1, oira2, unique);
//            }

            void PrintStats() {
                std::cout << "Unique Sams: " << m_pairs_unique.size() << std::endl;
                std::cout << "NonUnique Sams: " << m_pairs_nonunique.size() << std::endl;
                std::cout << "NonUnique Best Sams: " << m_pairs_nonunique_best.size() << std::endl;
            }

            bool HasReads() const {
                return !(m_pairs_unique.empty() && m_pairs_nonunique.empty() && m_pairs_nonunique_best.empty());
            }

            // Throws SamFormatError if a SAM header line names a reference sequence (@SQ) that is not a
            // gene of this database with the same length: the SAM was aligned against another database,
            // and its records would be compared with the wrong genes.
            void CheckReference(std::string const& line) {
                if (line.rfind("@SQ\t", 0) != 0) return;
                std::string_view name, length;
                for (size_t start = 4; start <= line.size();) {
                    size_t end = line.find('\t', start);
                    if (end == std::string::npos) end = line.size();
                    std::string_view field(line.data() + start, end - start);
                    if (field.rfind("SN:", 0) == 0) name = field.substr(3);
                    if (field.rfind("LN:", 0) == 0) length = field.substr(3);
                    start = end + 1;
                }
                auto number = [](std::string_view s, size_t& value) {
                    auto [end, error] = std::from_chars(s.data(), s.data() + s.size(), value);
                    return error == std::errc{} && end == s.data() + s.size() && !s.empty();
                };
                size_t const underscore = name.find('_');
                size_t const suffix = name.find('_', underscore == std::string_view::npos ? name.size() : underscore + 1);
                size_t taxid = 0, geneid = 0, sq_length = 0;
                if (underscore == std::string_view::npos || !number(name.substr(0, underscore), taxid) ||
                    !number(name.substr(underscore + 1, suffix - underscore - 1), geneid)) {
                    return;  // not a protal gene: its records are skipped as such
                }
                std::string const gene(name);
                if (!m_genome_loader.HasGene(taxid, geneid)) {
                    throw SamFormatError("gene " + gene + " of the SAM header (@SQ) is not in the database: the SAM "
                                         "was aligned against another database");
                }
                size_t const db_length = m_genome_loader.GeneLength(taxid, geneid);
                if (number(length, sq_length) && sq_length != db_length) {
                    throw SamFormatError("gene " + gene + " is " + std::string(length) + " bp in the SAM header (@SQ) but " +
                                         std::to_string(db_length) + " bp in the database: the SAM was aligned against "
                                         "another database");
                }
            }

            // Reads a SAM file (plain or gzipped) and hands each read's group of candidate alignments to
            // `on_group`: adjacent records with one QNAME are one read's candidates. Returns an error
            // message if the file cannot be read, is truncated or was aligned against another database
            // (CheckReference); a SAM without alignments is not an error.
            template<typename OnGroup>
            std::string ReadSamGroups(std::string const& file_path, OnGroup&& on_group) {
                if (std::filesystem::exists(file_path) && std::filesystem::file_size(file_path) == 0) {
                    return "the file is empty (not even a SAM header)";
                }
                igzstream file{ file_path.c_str() };
                if (!file.good()) return "cannot open the file";
                SamReader reader(file, [this](std::string const& line) { CheckReference(line); });
                // zlib reads a truncated or corrupt gzip file as one that ends early.
                auto truncated = [&file]() {
                    return "the file is truncated or corrupt (" + file.rdbuf()->read_error_message() + ")";
                };

                SamEntry sam1;
                SamEntry sam2;
                bool has_sam1 = false, has_sam2 = false;
                std::vector<AlignmentPair> group;

                try {
                    while (reader.Next(sam1, sam2, has_sam1, has_sam2)) {
                        AlignmentPair pair(
                                has_sam1 ? std::optional<SamEntry>{ sam1 } : std::optional<SamEntry>{},
                                has_sam2 ? std::optional<SamEntry>{ sam2 } : std::optional<SamEntry>{});
                        if (!group.empty() && !SameRead(pair, group.front())) {
                            on_group(group);
                            group.clear();
                        }
                        group.emplace_back(std::move(pair));
                    }
                } catch (SamFormatError const& e) {
                    // A truncated file's last line is cut short, too: the truncation is the cause.
                    if (file.rdbuf()->read_failed()) return truncated();
                    return e.what();
                }
                if (file.rdbuf()->read_failed()) return truncated();
                if (!group.empty()) on_group(group);

                for (auto const& [reason, count] : reader.Skipped()) {
                    std::cerr << file_path << ": skipped " << count << " record(s): " << reason << std::endl;
                }
                if (reader.Records() == 0) {
                    std::cerr << file_path << " contains no usable alignments" << std::endl;
                } else if (reader.RecordsWithoutTags() > 0) {
                    std::cerr << "Warning: " << reader.RecordsWithoutTags() << " of " << reader.Records() << " records in "
                              << file_path << " have no ZU tag (protal's unique k-mer count). A SAM file not "
                              << "written by protal lacks it, and the model then rejects most taxa." << std::endl;
                }
                return {};
            }

            // The best alignment of a multi-mapped read: its primary one (no 0x100; protal writes it first).
            static AlignmentPair& BestOfGroup(std::vector<AlignmentPair>& group) {
                auto best = std::find_if(group.begin(), group.end(), [](AlignmentPair& pair) {
                    return !Flag::IsNotPrimaryAlignment(pair.Any().m_flag);
                });
                return best == group.end() ? group.front() : *best;
            }

            // Loads the alignments of a SAM file into m_pairs_unique (reads with one candidate),
            // m_pairs_nonunique and m_pairs_nonunique_best (see BestOfGroup). ProfileSam profiles a file
            // without holding it; this keeps it, for inspection. Returns an error message as ReadSamGroups.
            std::string FromSam(std::string file_path, bool truth_in_header=false) {
                m_pairs_unique.clear();
                m_pairs_nonunique.clear();
                m_pairs_nonunique_best.clear();
                auto error = ReadSamGroups(file_path, [this](std::vector<AlignmentPair>& group) {
                    if (group.size() == 1) {
                        m_pairs_unique.emplace_back(std::move(group.front()));
                    } else {
                        m_pairs_nonunique_best.emplace_back(BestOfGroup(group));
                        m_pairs_nonunique.emplace_back(std::move(group));
                    }
                });
                if (!error.empty()) return error;

                if (truth_in_header) {
                    OutputErrorData(m_pairs_unique, m_pairs_nonunique);
                }
                return {};
            }


            std::string ErrorLine(int id, int length1, int score1, bool forward1, int length2, int score2, bool forward2, bool possibly_true, bool same_gene, int distance, bool paired, int ref_length) {
                std::string line;
                line += std::to_string(id) + '\t';
                line += std::to_string(length1) + '\t';
                line += std::to_string(score1) + '\t';
                line += std::to_string(forward1) + '\t';
                line += std::to_string(length2) + '\t';
                line += std::to_string(score2) + '\t';
                line += std::to_string(forward2) + '\t';
                line += std::to_string(possibly_true) + '\t';
                line += std::to_string(same_gene) + '\t';
                line += std::to_string(distance) + '\t';
                line += std::to_string(paired) + '\t';
                line += std::to_string(ref_length);
                return line;
            }

            void TestSNPUtils(std::vector<AlignmentPair> &pairs) {
//                std::cout << "TestSNPUtils" << std::endl;
                SNPList snps;
                size_t sum = 0;
                Benchmark bm("SNPs");
                for (auto &pair: pairs) {
                    auto &any = pair.Any();
                    auto &[tid, gid] = ExtractTaxidGeneid(any.m_rname);
                    auto &gene = m_genome_loader.GetGenome(tid).GetGeneOMP(gid);

                    gene.LoadOMP();
                    auto &ref = gene.Sequence();
//                        std::cout << "Extract: " << any.m_qname << " ---> " <<  any.m_rname << std::endl;

                    PrintAlignment(any, ref);
                    ExtractSNPs(any, ref, snps, tid, gid);
//                    for (auto &snp: snps) {
//                        std::cout << snp.ToString() << std::endl;
//                    }
//                    sum += snps.size();
//                    std::cout << ".........>" << std::endl;
//                    Utils::Input();
                }
                bm.PrintResults();
                std::cout << "Sum: " << sum << std::endl;

            }

            void TestSNPUtils(std::vector<std::vector<AlignmentPair>> &pair_lists) {
//                std::cout << "TestSNPUtils" << std::endl;
                SNPList snps;
                size_t sum = 0;
                Benchmark bm("SNPs");
                for (auto& pairlist : pair_lists) {
                    for (auto &pair: pairlist) {
                        auto& any = pair.Any();
                        auto &[ tid, gid ] = ExtractTaxidGeneid(any.m_rname);
                        auto& gene = m_genome_loader.GetGenome(tid).GetGene(gid);
                        gene.LoadOMP();
                        auto& ref = gene.Sequence();
//                        std::cout << "Extract: " << any.m_qname << " ---> " <<  any.m_rname << std::endl;

                        PrintAlignment(any, ref);
                        ExtractSNPs(any, ref, snps, tid, gid);
                        for (auto& snp : snps) {
                            std::cout << snp.ToString() << std::endl;
                        }
                        sum += snps.size();
                    }
                }
                bm.PrintResults();
                std::cout << "Sum: " << sum << std::endl;

            }

            void OutputErrorData(std::vector<std::vector<AlignmentPair>> &pair_lists) {
                int pairlist_id = 0;
                CigarInfo info1;
                CigarInfo info2;

                size_t count_truth = 0;
                size_t count_false = 0;
                for (auto& pairlist : pair_lists) {
                    for (auto& pair : pairlist) {
                        auto &[ tid, gid ] = ExtractTaxidGeneid(pair.Any().m_rname);
                        auto &[ atid, agid ] = ExtractTaxidGeneid(pair.Any().m_qname);
//                        std::cerr << pair.Any().m_qname << " " << atid << " " << agid << std::endl;
//                        std::cerr << "CAPTURE\t";
                        if (pair.IsPair()) {
                            auto &[ tid1, gid1 ] = ExtractTaxidGeneid(pair.First().m_rname);
                            auto &[ tid2, gid2 ] = ExtractTaxidGeneid(pair.Second().m_rname);
                            auto &[ ttid1, tgid1 ] = ExtractTaxidGeneid(pair.First().m_qname);
                            auto &[ ttid2, tgid2 ] = ExtractTaxidGeneid(pair.Second().m_qname);
                            CompressedCigarInfo(pair.First().m_cigar, info1);
                            CompressedCigarInfo(pair.Second().m_cigar, info2);
                            bool same_gene = gid1 == gid2;
                            int insert = std::max(pair.First().m_pos, pair.Second().m_pos) - std::min(pair.First().m_pos, pair.Second().m_pos);
                            count_truth += tid == ttid1;
                            count_false += tid != ttid1;
//                            std::cerr << ErrorLine(pairlist_id, info1.clipped_alignment_length, info1.Score(), Flag::IsRead1ReverseComplement(pair.First().m_flag), info2.clipped_alignment_length, info2.Score(), Flag::IsRead2ReverseComplement(pair.Second().m_flag),  tid == ttid1, same_gene, insert, true, m_genome_loader.GetGenome(tid1).GetGene(gid1).GetLength()) << std::endl;
                        } else if (pair.HasFirst()) {
                            CompressedCigarInfo(pair.First().m_cigar, info1);
                            auto &[ ttid, tgid ] = ExtractTaxidGeneid(pair.First().m_qname);
                            count_truth += tid == ttid;
                            count_false += tid != ttid;
//                            std::cerr << ErrorLine(pairlist_id, info1.clipped_alignment_length, info1.Score(), Flag::IsRead1ReverseComplement(pair.First().m_flag), -1, -1, false, tid == ttid, false, -1, false, m_genome_loader.GetGenome(tid).GetGene(gid).GetLength()) << std::endl;
                        } else if (pair.HasSecond()) {
                            CompressedCigarInfo(pair.Second().m_cigar, info2);
                            auto &[ ttid, tgid ] = ExtractTaxidGeneid(pair.Second().m_qname);
                            count_truth += tid == ttid;
                            count_false += tid != ttid;
//                            std::cerr << ErrorLine(pairlist_id, -1, -1, false, info2.clipped_alignment_length, info2.Score(), Flag::IsRead2ReverseComplement(pair.Second().m_flag), tid == ttid, false, -1, false, m_genome_loader.GetGenome(tid).GetGene(gid).GetLength()) << std::endl;
                        }
                    }

                    pairlist_id++;
                }
                std::cout << "Truth: " << count_truth << std::endl;
                std::cout << "False: " << count_false << std::endl;
            }

//            void OutputErrorData(TruthSet const& set, std::vector<std::vector<AlignmentPair>> &pair_lists) {
//                int pairlist_id = 0;
//                CigarInfo info1;
//                CigarInfo info2;
//                for (auto& pairlist : pair_lists) {
//                    for (auto& pair : pairlist) {
//                        auto &[ tid, gid ] = ExtractTaxidGeneid(pair.Any().m_rname);
//                        bool maybe_true = set.contains(tid);
//                        std::cerr << "CAPTURE\t";
//                        if (pair.IsPair()) {
//                            auto &[ tid1, gid1 ] = ExtractTaxidGeneid(pair.First().m_rname);
//                            auto &[ tid2, gid2 ] = ExtractTaxidGeneid(pair.Second().m_rname);
//                            CompressedCigarInfo(pair.First().m_cigar, info1);
//                            CompressedCigarInfo(pair.Second().m_cigar, info2);
//                            bool same_gene = gid1 == gid2;
//                            int insert = std::max(pair.First().m_pos, pair.Second().m_pos) - std::min(pair.First().m_pos, pair.Second().m_pos);
//
//                            std::cerr << ErrorLine(pairlist_id, info1.clipped_alignment_length, info1.Score(), Flag::IsRead1ReverseComplement(pair.First().m_flag), info2.clipped_alignment_length, info2.Score(), Flag::IsRead2ReverseComplement(pair.Second().m_flag),  maybe_true, same_gene, insert, true, m_genome_loader.GetGenome(tid1).GetGene(gid1).GetLength()) << std::endl;
//                        } else if (pair.HasFirst()) {
//                            CompressedCigarInfo(pair.First().m_cigar, info1);
//                            std::cerr << ErrorLine(pairlist_id, info1.clipped_alignment_length, info1.Score(), Flag::IsRead1ReverseComplement(pair.First().m_flag), -1, -1, false, maybe_true, false, -1, false, m_genome_loader.GetGenome(tid).GetGene(gid).GetLength()) << std::endl;
//                        } else if (pair.HasSecond()) {
//                            CompressedCigarInfo(pair.Second().m_cigar, info2);
//                            std::cerr << ErrorLine(pairlist_id, -1, -1, false, info2.clipped_alignment_length, info2.Score(), Flag::IsRead2ReverseComplement(pair.Second().m_flag), maybe_true, false, -1, false, m_genome_loader.GetGenome(tid).GetGene(gid).GetLength()) << std::endl;
//                        }
//                    }
//
//                    pairlist_id++;
//                }
//            }


            std::pair<bool,bool> IsSamCorrect(SamEntry &entry) {
                auto &[ rtid, rgid ] = ExtractTaxidGeneid(entry.m_rname);
                auto &[ qtid, qgid ] = ExtractTaxidGeneid(entry.m_rname);
                return { rtid == qtid, rtid == qtid && rgid == qgid };
            }

            std::pair<bool,bool> IsCorrect(SamEntry &entry, TruthSet* set=nullptr) {
                if (set != nullptr) {
                    auto &[ tid, gid ] = ExtractTaxidGeneid(entry.m_rname);
                    return { set->contains(tid), false };
                } else {
                    return IsSamCorrect(entry);
                }
            }

            void OutputErrorData(SamPairs &pairs, SamPairList &pair_lists, TruthSet* set=nullptr) {
                int pairlist_id = 0;
                CigarInfo info1;
                CigarInfo info2;

                size_t count_truth = 0;
                size_t count_false = 0;
                for (auto& pairlist : pair_lists) {
                    for (auto& pair : pairlist) {
//                        std::cerr << "CAPTURE\t";
                        if (pair.IsPair()) {
                            auto &[ tid1, gid1 ] = ExtractTaxidGeneid(pair.First().m_rname);
                            auto &[ tid2, gid2 ] = ExtractTaxidGeneid(pair.Second().m_rname);
                            auto [true_taxid1, true_geneid1] = IsCorrect(pair.First(), set);
                            auto [true_taxid2, true_geneid2] = IsCorrect(pair.Second(), set);
                            auto &[read_truth_tid, read_truth_gid] = ExtractTaxidGeneid(pair.First().m_qname);
                            bool same_gene = gid1 == gid2;
                            CompressedCigarInfo(pair.First().m_cigar, info1);
                            CompressedCigarInfo(pair.Second().m_cigar, info2);

                            count_truth += read_truth_tid == tid1;
                            count_false += read_truth_tid != tid1;

                            int insert = std::max(pair.First().m_pos, pair.Second().m_pos) - std::min(pair.First().m_pos, pair.Second().m_pos);
//                            std::cerr << ErrorLine(pairlist_id, info1.clipped_alignment_length, info1.Score(), Flag::IsRead1ReverseComplement(pair.First().m_flag), info2.clipped_alignment_length, info2.Score(), Flag::IsRead2ReverseComplement(pair.Second().m_flag), true_taxid1, same_gene, insert, true, m_genome_loader.GetGenome(tid1).GetGene(gid1).GetLength()) << std::endl;
                        } else if (pair.HasFirst()) {
                            auto &[ tid, gid ] = ExtractTaxidGeneid(pair.First().m_rname);
                            auto [true_taxid, true_geneid] = IsCorrect(pair.First(), set);
                            auto &[read_truth_tid, read_truth_gid] = ExtractTaxidGeneid(pair.First().m_qname);
                            CompressedCigarInfo(pair.First().m_cigar, info1);
                            count_truth += read_truth_tid == tid;
                            count_false += read_truth_tid != tid;
//                            std::cerr << ErrorLine(pairlist_id, info1.clipped_alignment_length, info1.Score(), Flag::IsRead1ReverseComplement(pair.First().m_flag), -1, -1, false, true_taxid, false, -1, false, m_genome_loader.GetGenome(tid).GetGene(gid).GetLength()) << std::endl;
                        } else if (pair.HasSecond()) {
                            auto &[ tid, gid ] = ExtractTaxidGeneid(pair.Second().m_rname);
                            auto [true_taxid, true_geneid] = IsCorrect(pair.Second(), set);
                            auto &[read_truth_tid, read_truth_gid] = ExtractTaxidGeneid(pair.Second().m_qname);
                            CompressedCigarInfo(pair.Second().m_cigar, info2);
                            count_truth += read_truth_tid == tid;
                            count_false += read_truth_tid != tid;

//                            std::cerr << ErrorLine(pairlist_id, -1, -1, false, info2.clipped_alignment_length, info2.Score(), Flag::IsRead2ReverseComplement(pair.Second().m_flag), true_taxid, false, -1, false, m_genome_loader.GetGenome(tid).GetGene(gid).GetLength()) << std::endl;
                        }
                        break;
                    }

                    pairlist_id++;
                }
                for (auto &pair: pairs) {
//                    std::cerr << "CAPTURE\t";
                    if (pair.IsPair()) {
                        auto &[tid1, gid1] = ExtractTaxidGeneid(pair.First().m_rname);
                        auto &[tid2, gid2] = ExtractTaxidGeneid(pair.Second().m_rname);
                        auto &[read_truth_tid, read_truth_gid] = ExtractTaxidGeneid(pair.First().m_qname);
                        auto [true_taxid1, true_geneid1] = IsCorrect(pair.First(), set);
                        auto [true_taxid2, true_geneid2] = IsCorrect(pair.Second(), set);
                        bool same_gene = gid1 == gid2;

                        std::cout << tid1 << " " << read_truth_tid << std::endl;

                        CompressedCigarInfo(pair.First().m_cigar, info1);
                        CompressedCigarInfo(pair.Second().m_cigar, info2);
                        count_truth += read_truth_tid == tid1;
                        count_false += read_truth_tid != tid1;
                        int insert = std::max(pair.First().m_pos, pair.Second().m_pos) -
                                     std::min(pair.First().m_pos, pair.Second().m_pos);
//                        std::cerr << ErrorLine(pairlist_id, info1.clipped_alignment_length, info1.Score(),
//                                               Flag::IsRead1ReverseComplement(pair.First().m_flag),
//                                               info2.clipped_alignment_length, info2.Score(),
//                                               Flag::IsRead2ReverseComplement(pair.Second().m_flag), true_taxid1,
//                                               same_gene, insert, true,
//                                               m_genome_loader.GetGenome(tid1).GetGene(gid1).GetLength()) << std::endl;
                    } else if (pair.HasFirst()) {
                        auto [true_taxid, true_geneid] = IsCorrect(pair.First(), set);
                        auto &[ tid, gid ] = ExtractTaxidGeneid(pair.First().m_rname);
                        auto &[read_truth_tid, read_truth_gid] = ExtractTaxidGeneid(pair.First().m_qname);
                        CompressedCigarInfo(pair.First().m_cigar, info1);
                        count_truth += read_truth_tid == tid;
                        count_false += read_truth_tid != tid;
//                        std::cerr << ErrorLine(pairlist_id, info1.clipped_alignment_length, info1.Score(),
//                                               Flag::IsRead1ReverseComplement(pair.First().m_flag), -1, -1, false,
//                                               true_taxid, false, -1, false,
//                                               m_genome_loader.GetGenome(tid).GetGene(gid).GetLength()) << std::endl;
                    } else if (pair.HasSecond()) {
                        auto &[ tid, gid ] = ExtractTaxidGeneid(pair.Second().m_rname);
                        auto [true_taxid, true_geneid] = IsCorrect(pair.Second(), set);
                        auto &[read_truth_tid, read_truth_gid] = ExtractTaxidGeneid(pair.Second().m_qname);
                        CompressedCigarInfo(pair.Second().m_cigar, info2);
                        count_truth += read_truth_tid == tid;
                        count_false += read_truth_tid != tid;
//                        std::cerr
//                                << ErrorLine(pairlist_id, -1, -1, false, info2.clipped_alignment_length, info2.Score(),
//                                             Flag::IsRead2ReverseComplement(pair.Second().m_flag), true_taxid, false,
//                                             -1, false, m_genome_loader.GetGenome(tid).GetGene(gid).GetLength())
//                                << std::endl;
                    }
                }


                std::cout << "Truth: " << count_truth << std::endl;
                std::cout << "False: " << count_false << std::endl;

                pairlist_id++;
            }

            Benchmark m_add_sam{"Add sam.."};
            Benchmark m_cigar_info{"Compressed cigar info"};


            bool ProcessMAPQ(MicrobialProfile& profile, AlignmentPair& ap, int read_id=0) {
                bool valid_sam = true;

                bool take_first = false;
                bool take_second = false;
                m_cigar_info.Start();
                if (ap.HasFirst()) CompressedCigarInfo(ap.First().m_cigar, m_info1);
                if (ap.HasSecond()) CompressedCigarInfo(ap.Second().m_cigar, m_info2);
                m_cigar_info.Stop();

                // MAPQ is judged per mate: the mates of a pair aligned together share one MAPQ, those
                // of a fragment split over two genes each have their own.
                if (ap.HasFirst() && ap.First().m_mapq >= m_min_mapq && m_info1.clipped_alignment_length > m_min_alignment_length) {
                    take_first = true;
                }
                if (ap.HasSecond() && ap.Second().m_mapq >= m_min_mapq && m_info2.clipped_alignment_length > m_min_alignment_length) {
                    take_second = true;
                }
                if (!take_first && !take_second) return true;


                if (take_first && take_second) {
                    auto [tid1, geneid1] = ExtractTaxidGeneid(ap.First().m_rname);
                    auto [tid2, geneid2] = ExtractTaxidGeneid(ap.Second().m_rname);

                    m_add_sam.Start();
                    valid_sam &= profile.AddSam(tid1, geneid1, ap.First(), m_info1.Ani(), true, read_id, kNoStrain);
                    valid_sam &= profile.AddSam(tid2, geneid2, ap.Second(), m_info2.Ani(), true, read_id, kNoStrain);
                    m_add_sam.Stop();

                } else if (take_first) {
                    auto [tid, geneid] = ExtractTaxidGeneid(ap.First().m_rname);

                    m_add_sam.Start();
                    valid_sam &=profile.AddSam(tid, geneid, ap.First(), m_info1.Ani(), true, read_id, kNoStrain);
                    m_add_sam.Stop();



                } else if (take_second) {
                    auto [tid, geneid] = ExtractTaxidGeneid(ap.Second().m_rname);

                    m_add_sam.Start();
                    valid_sam &=profile.AddSam(tid, geneid, ap.Second(), m_info2.Ani(), true, read_id, kNoStrain);
                    m_add_sam.Stop();
                }
                return valid_sam;
            }

            void ProcessUnique(MicrobialProfile& profile, AlignmentPair& ap) {
                bool take_first = false;
                bool take_second = false;
                if (ap.HasFirst()) CompressedCigarInfo(ap.First().m_cigar, m_info1);
                if (ap.HasSecond()) CompressedCigarInfo(ap.Second().m_cigar, m_info2);

                if (ap.HasFirst() && ap.HasSecond() && (m_info1.clipped_alignment_length + m_info2.clipped_alignment_length) > m_min_alignment_length) {
                    take_first = true;
                    take_second = true;
                } else if (ap.HasFirst() && m_info1.clipped_alignment_length > m_min_alignment_length) {
                    take_first = true;
                } else if (ap.HasSecond() && m_info2.clipped_alignment_length > m_min_alignment_length) {
                    take_second = true;
                }

                if (!take_first && !take_second) return;

                if (take_first && take_second) {
                    auto [tid1, geneid1] = ExtractTaxidGeneid(ap.First().m_rname);
                    auto [tid2, geneid2] = ExtractTaxidGeneid(ap.Second().m_rname);

                    double score = m_score.Score(ap.First(), ap.Second());
                    double score1 = m_score.Score(ap.First());
                    double score2 = m_score.Score(ap.Second());
                    double score_diff = std::abs(score1 - score2);

                    profile.AddSam(tid1, geneid1, ap.First(), score, true);
                    profile.AddSam(tid2, geneid2, ap.Second(), score, true);

                } else if (take_first) {
                    auto [tid, geneid] = ExtractTaxidGeneid(ap.First().m_rname);
                    double score = m_score.Score(ap.First());

                    profile.AddSam(tid, geneid, ap.First(), score, true);
                } else if (take_second) {
                    auto [tid, geneid] = ExtractTaxidGeneid(ap.Second().m_rname);
                    double score = m_score.Score(ap.Second());

                    profile.AddSam(tid, geneid, ap.Second(), score, true);
                }
            }

            using OptionalRefOstream = std::optional<std::reference_wrapper<std::ostream>>;
            // Profiles a SAM file one read at a time, without holding the file in memory: a read with
            // one candidate alignment, or the best of a multi-mapped read's (BestOfGroup), is added to
            // `profile`. Records the profile rejects (genes outside the database, alignments that do
            // not match their gene) go to `erroneous_sam_out`. Returns an error message as ReadSamGroups.
            std::string ProfileSam(std::string const& file_path, MicrobialProfile& profile, OptionalRefOstream erroneous_sam_out={}, size_t snp_min_cov=2, size_t snp_min_obs_fwdrev=2, double snp_min_af=0.0, size_t snp_min_mean_qual=15, size_t snp_min_phred_sum=0, bool snp_require_strand=false) {
                profile.SetDepthIdentityMargin(m_depth_identity_margin);
                size_t read_id = 0;
                m_rejected_reads = 0;
                auto error = ReadSamGroups(file_path, [&](std::vector<AlignmentPair>& group) {
                    auto& pair = group.size() == 1 ? group.front() : BestOfGroup(group);
                    bool const valid = ProcessMAPQ(profile, pair, read_id++);
                    m_rejected_reads += !valid;
                    if (!valid && erroneous_sam_out.has_value()) {
                        auto& os = erroneous_sam_out.value().get();
                        if (pair.first.has_value()) os << pair.first.value().ToString() << '\n';
                        if (pair.second.has_value()) os << pair.second.value().ToString() << '\n';
                    }
                });
                if (!error.empty()) return error;

                m_reads = read_id;
                m_post_process_bm.Start();
                profile.PostProcessSNPs(snp_min_cov, snp_min_obs_fwdrev, snp_min_af, snp_min_mean_qual, snp_min_phred_sum, snp_require_strand);
                m_post_process_bm.Stop();
                return {};
            }

            // Reads of the last ProfileSam, and those with an alignment the database rejected (a gene it
            // lacks, a position past a gene's end, or bases that do not match the gene).
            size_t Reads() const { return m_reads; }
            size_t RejectedReads() const { return m_rejected_reads; }
        };
    }
}
