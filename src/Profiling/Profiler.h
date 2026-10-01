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
#include "SamFile.h"
#include "cPMML.h"
#include "sparse_map.h"
#include "Benchmark.h"
#include "RunStatus.h"
#include "ReadType.h"
#include "SamChunks.h"
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <charconv>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <set>
#include <sstream>
#include <unordered_map>
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
                auto const [identity, length] = AlignmentIdentity(sam.m_cigar);
                if (!no_strain) {
                    // The read's identity lets the strain MSA keep the taxon's own reads only.
                    if (!m_strain_level.AddSam(sam, read_id, true, identity)) return false;
                }

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

        inline constexpr MAPQ_t kLowMapq = 10;  // a record below fits another candidate nearly as well
        inline constexpr int kAlternativeFitEdits = 1;  // an alternative this many edits worse fits as well
        inline constexpr double kHighExcess = 0.02;  // a record's divergence beyond its base qualities above this is high

        // A record's differences per aligned base (X + I + D over M + X + I + D, as AlignmentIdentity counts them) less
        // the mean error probability of its bases by their qualities (10^(-Q/10)): how far its genome differs from the
        // reference beyond what its sequencing errors explain. A relative's reads exceed their errors by several
        // percent, a present species' own reads by about one (docs/claude/2026-10-01-f1-opportunities). None without
        // base qualities ("*") or aligned bases.
        inline std::optional<float> ReadExcess(SamEntry const& sam) {
            if (sam.m_qual.empty() || sam.m_qual == "*") return std::nullopt;
            static std::array<double, 256> const error = [] {
                std::array<double, 256> e{};
                for (int c = 0; c < 256; c++) e[c] = c < 33 ? 1.0 : std::pow(10.0, -(c - 33) / 10.0);
                return e;
            }();
            size_t matches = 0, differences = 0, run = 0;
            for (char const c : sam.m_cigar) {
                if (c >= '0' && c <= '9') {
                    run = run * 10 + static_cast<size_t>(c - '0');
                    continue;
                }
                if (c == 'M' || c == '=') matches += run;
                else if (c == 'X' || c == 'I' || c == 'D') differences += run;
                run = 0;
            }
            if (matches + differences == 0) return std::nullopt;
            double expected = 0;
            for (unsigned char const c : sam.m_qual) expected += error[c];
            expected /= static_cast<double>(sam.m_qual.size());
            return static_cast<float>(static_cast<double>(differences) / static_cast<double>(matches + differences) - expected);
        }

        // The median of values (the mean of the two middle ones for an even count); 0 if empty.
        inline double MedianOf(std::vector<double> values) {
            if (values.empty()) return 0;
            size_t const mid = values.size() / 2;
            std::nth_element(values.begin(), values.begin() + mid, values.end());
            double const upper = values[mid];
            if (values.size() % 2 == 1) return upper;
            return (*std::max_element(values.begin(), values.begin() + mid) + upper) / 2;
        }

        // What a read's record tells beyond its own alignment (MicrobialProfile::AddSam).
        struct ReadEvidence {
            size_t link = SIZE_MAX;  // the read across its records (both mates, a long read's genes); SIZE_MAX: none
        };

        // A taxon's best records of all reads, before the profiler's MAPQ and length filters drop any (the reads that
        // fit another taxon as well have MAPQ near 0 and would never be counted): MicrobialProfile::NoteRecord.
        struct RecordEvidence {
            size_t records = 0;
            size_t low_mapq = 0;  // MAPQ below kLowMapq
            size_t congener_fit = 0;  // another species of the genus fits the read within kAlternativeFitEdits (ZA)
            size_t other_genus_fit = 0;  // a species of another genus does
            size_t adjacent = 0;  // genes next to each other on its reads (MicrobialProfile::NoteLinkedRecord)
            size_t adjacent_expected = 0;  // of these, whose ends face each other in its clade (gene_neighbours::Verdict::Expected)
            size_t adjacent_unlikely = 0;  // that never do in a clade with data on them (Verdict::Unlikely)
            std::vector<float> excess;  // each record's ReadExcess (records with base qualities)

            RecordEvidence& operator+=(RecordEvidence const& other) {
                records += other.records;
                low_mapq += other.low_mapq;
                congener_fit += other.congener_fit;
                other_genus_fit += other.other_genus_fit;
                adjacent += other.adjacent;
                adjacent_expected += other.adjacent_expected;
                adjacent_unlikely += other.adjacent_unlikely;
                excess.insert(excess.end(), other.excess.begin(), other.excess.end());
                return *this;
            }
        };

        // Calls on_alternative(taxid, edits more) for each entry of a ZA tag ("12:0,40:3"; "*" or empty: none).
        template<typename F>
        inline void ForEachAlternative(std::string const& tag, F&& on_alternative) {
            if (tag.empty() || tag == "*") return;
            char const* p = tag.data();
            char const* const end = p + tag.size();
            while (p < end) {
                uint32_t taxid = 0;
                int more = 0;
                auto r1 = std::from_chars(p, end, taxid);
                if (r1.ec != std::errc() || r1.ptr == end || *r1.ptr != ':') return;
                auto r2 = std::from_chars(r1.ptr + 1, end, more);
                if (r2.ec != std::errc()) return;
                on_alternative(taxid, more);
                p = r2.ptr;
                if (p < end && *p == ',') p++;
            }
        }

        // The bases clipped (H, S) at the start of a CIGAR, or at its end (trailing): where a long read's record
        // starts on the read is the clip at its start, or at its end if it is reverse.
        inline uint32_t Clip(std::string_view cigar, bool trailing) {
            uint32_t total = 0;
            if (!trailing) {
                size_t i = 0;
                while (i < cigar.size()) {
                    uint32_t n = 0;
                    size_t j = i;
                    while (j < cigar.size() && cigar[j] >= '0' && cigar[j] <= '9') n = n * 10 + static_cast<uint32_t>(cigar[j++] - '0');
                    if (j == i || j >= cigar.size() || (cigar[j] != 'H' && cigar[j] != 'S')) break;
                    total += n;
                    i = j + 1;
                }
                return total;
            }
            size_t end = cigar.size();
            while (end > 0 && (cigar[end - 1] == 'H' || cigar[end - 1] == 'S')) {
                size_t j = end - 1;
                uint32_t n = 0, scale = 1;
                while (j > 0 && cigar[j - 1] >= '0' && cigar[j - 1] <= '9') {
                    n += static_cast<uint32_t>(cigar[j - 1] - '0') * scale;
                    scale *= 10;
                    j--;
                }
                if (j == end - 1) break;
                total += n;
                end = j;
            }
            return total;
        }

        // The read's bases a CIGAR aligns (M, I, =, X).
        inline uint32_t QueryBases(std::string_view cigar) {
            uint32_t total = 0, n = 0;
            for (char const c : cigar) {
                if (c >= '0' && c <= '9') {
                    n = n * 10 + static_cast<uint32_t>(c - '0');
                    continue;
                }
                if (c == 'M' || c == 'I' || c == '=' || c == 'X') total += n;
                n = 0;
            }
            return total;
        }

        class Taxon {
        private:
            size_t m_id;
            std::string m_name;
            size_t m_total_hits = 0;
            RecordEvidence m_records;  // all of its best records, before the filters (MicrobialProfile::ApplyRecordEvidence)
            size_t m_links = 0;  // reads with a record here (a pair once, a long read once with all its genes)
            size_t m_linked = 0;  // of those, with two or more records here: both mates, or two of a long read's genes
            size_t m_last_link = SIZE_MAX;
            size_t m_link_records = 0;  // records of m_last_link here
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
            gene_conservation::Table const* m_conservation = nullptr;  // none: every gene's factor is 1
            bool m_scale_margin = false;  // the depth identity margin scaled by the factors (--gene_conservation db)
            double m_excess_median = 0;  // see ExcessMedian
            double m_excess_high_share = 0;

        public:

            // conservation: the genes' conservation factors (GenomeLoader::GetGeneConservation), for the features and,
            // with scale_margin, to scale the depth identity margin per gene.
            Taxon(Genome& genome, gene_conservation::Table const* conservation = nullptr, bool scale_margin = false) :
                    m_genome(&genome), m_genome_gene_count(genome.GeneNum()), m_conservation(conservation),
                    m_scale_margin(scale_margin) {}

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
                    m_genome->LoadGeneOMP(geneid);
                    m_genes.insert( { geneid, profiler::Gene(g) } );
                    m_genes.at(geneid).SetLength(m_genome->GetGene(geneid).GetLength());
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

            bool AddSam(GeneId geneid, SamEntry const& sam, double score, bool unique, size_t read_id, bool no_strain=true,
                        ReadEvidence const& evidence = {}) {
                Changed();
                bool const new_gene = !m_genes.contains(geneid);
                if (new_gene) {
                    auto& g = m_genome->GetGene(geneid);
                    m_genome->LoadGeneOMP(geneid);
                    m_genes.insert( { geneid, profiler::Gene(g) } );
                    m_genes.at(geneid).SetLength(m_genome->GetGene(geneid).GetLength());
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
                if (evidence.link != SIZE_MAX) {
                    if (evidence.link != m_last_link) {
                        m_links++;
                        m_last_link = evidence.link;
                        m_link_records = 0;
                    }
                    m_linked += ++m_link_records == 2;
                }
                return true;
            }

            // Takes the counts of the taxon's best records (MicrobialProfile::ApplyRecordEvidence), and of their excesses
            // only the median and the high share.
            void SetRecordEvidence(RecordEvidence const& records) {
                m_records = records;
                m_records.excess.clear();
                m_records.excess.shrink_to_fit();
                m_excess_median = 0;
                m_excess_high_share = 0;
                if (!records.excess.empty()) {
                    std::vector<double> values(records.excess.begin(), records.excess.end());
                    m_excess_high_share = static_cast<double>(std::count_if(values.begin(), values.end(),
                        [](double v) { return v > kHighExcess; })) / static_cast<double>(values.size());
                    m_excess_median = MedianOf(std::move(values));
                }
            }

            // The median ReadExcess of the taxon's best records (all reads, before the filters), and the share above
            // kHighExcess; both 0 without base qualities.
            double ExcessMedian() const { return m_excess_median; }
            double ExcessHighShare() const { return m_excess_high_share; }

            // The conservation pattern of the genes its reads hit, by their factors (gene_conservation.tsv): log2 of the
            // median depth of its hit genes with factor below 1 (conserved) over that of the others, each + 0.001 (0 if
            // either has none); and the conserved share of its hit genes (0.5 without genes or factors). A relative the
            // database lacks makes a species' fast genes deeper (docs/claude/2026-10-01-gene-scaled-margin).
            std::pair<double, double> ConservationPattern() const {
                if (!m_conservation || m_conservation->Empty() || m_genes.empty()) return { 0.0, 0.5 };
                std::vector<double> conserved, fast;
                for (auto const& [id, gene] : m_genes) {
                    (m_conservation->Factor(id) < 1 ? conserved : fast).push_back(gene.VerticalCoverage());
                }
                double const share = static_cast<double>(conserved.size()) / static_cast<double>(m_genes.size());
                if (conserved.empty() || fast.empty()) return { 0.0, share };
                return { std::log2((MedianOf(std::move(conserved)) + 1e-3) / (MedianOf(std::move(fast)) + 1e-3)), share };
            }

            // Shares of the taxon's best records of all reads (also those the filters left out): with MAPQ below
            // kLowMapq, and whose read another species of the genus, or one of another genus, fits within
            // kAlternativeFitEdits edits (the ZA tag).
            double LowMapqShare() const { return m_records.records == 0 ? 0 : m_records.low_mapq / static_cast<double>(m_records.records); }
            double CongenerFitShare() const { return m_records.records == 0 ? 0 : m_records.congener_fit / static_cast<double>(m_records.records); }
            double OtherGenusFitShare() const { return m_records.records == 0 ? 0 : m_records.other_genus_fit / static_cast<double>(m_records.records); }
            // Share of the taxon's reads with two or more records on it: both mates of a pair (on one gene or
            // two), or two or more genes of a long read. Single-end reads have one record each: 0.
            double LinkedShare() const { return m_links == 0 ? 0 : m_linked / static_cast<double>(m_links); }
            // Of the genes next to each other on the taxon's reads (a pair's mates on two genes, a long read's consecutive
            // genes; MicrobialProfile::NoteLinkedRecord), the share whose ends face each other in the taxon's clade, and
            // the share whose never do there; both 0 without the database's gene neighbours.
            double AdjacentExpectedShare() const { return m_records.adjacent == 0 ? 0 : m_records.adjacent_expected / static_cast<double>(m_records.adjacent); }
            double AdjacentUnlikelyShare() const { return m_records.adjacent == 0 ? 0 : m_records.adjacent_unlikely / static_cast<double>(m_records.adjacent); }


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

            // The lowest identity of a read on gene `geneid` that counts towards the taxon's depth: the
            // depth identity margin below TopIdentity, the same on every gene unless conservation factors
            // are given (--gene_conservation db): then gene_conservation::GeneMargin, 0.03 + 0.05 x the
            // factor at 0.08, as a strain's reads on a fast gene sit further below the best ones. A
            // present species' own reads form this top cluster; reads of relatives (absent from the
            // database, or much more abundant) align at lower identity and would inflate its depth. They
            // still count for detection: the model's features use every read.
            double OwnIdentityThreshold(uint64_t geneid) const {
                if (m_depth_identity_margin >= 1 || PresentGenes() == 0) return 0;
                return TopIdentity() - gene_conservation::GeneMargin(m_depth_identity_margin, GeneFactor(geneid));
            }

            // The conservation factor of gene `geneid` the depth identity margin is scaled by: 1 unless the margin is
            // scaled (--gene_conservation db or FILE) and factors are given.
            double GeneFactor(uint64_t geneid) const {
                return m_scale_margin && m_conservation ? m_conservation->Factor(geneid) : 1.0;
            }

            // The lowest identity of a read within `margin` of TopIdentity; 0 (every read) for a margin of 1
            // or more. The strain MSA takes its reads with a stricter margin than the depth
            // (--msa_identity_margin): a relative's reads that a wider margin admits add false alleles.
            double IdentityThreshold(double margin) const {
                if (margin >= 1 || PresentGenes() == 0) return 0;
                return TopIdentity() - margin;
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
                    std::vector<double> vcovs;
                    size_t own_bases = 0, all_bases = 0;
                    for (auto& [geneid, gene] : m_genes) {
                        size_t const bases = gene.MappedLength(OwnIdentityThreshold(geneid));
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

            // Share of the taxon's aligned bases below their gene's OwnIdentityThreshold: reads of
            // relatives, e.g. of a species the database lacks.
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
            // Evidence from the reads' other candidates and from their other records: the share of records
            // with low MAPQ, of records whose read a congener or a species of another genus fits as well (ZA,
            // within kAlternativeFitEdits), and of reads with two or more records on the taxon (both mates, or
            // two genes of a long read).
            f.emplace_back("low_mapq_share", taxon.LowMapqShare());
            f.emplace_back("congener_fit_share", taxon.CongenerFitShare());
            f.emplace_back("other_genus_fit_share", taxon.OtherGenusFitShare());
            f.emplace_back("linked_share", taxon.LinkedShare());
            // Whether the genes next to each other on its reads are neighbours in the taxon's clade (the database's gene
            // neighbours; 0 without them): reads of the taxon itself, or of a relative of the same gene order, mostly
            // are; reads of genes that crossed from elsewhere are not.
            f.emplace_back("adjacent_expected_share", taxon.AdjacentExpectedShare());
            f.emplace_back("adjacent_unlikely_share", taxon.AdjacentUnlikelyShare());
            // How far its reads differ from the reference beyond their base qualities' errors (ReadExcess: the median,
            // and the share above kHighExcess), and the depth of its conserved hit genes against its fast ones
            // (ConservationPattern): a relative's reads exceed their errors and land on the fast genes
            // (docs/claude/2026-10-01-f1-opportunities).
            f.emplace_back("excess_median", taxon.ExcessMedian());
            f.emplace_back("excess_high_share", taxon.ExcessHighShare());
            auto const [conserved_ratio, conserved_share] = taxon.ConservationPattern();
            f.emplace_back("conserved_fast_depth_ratio", conserved_ratio);
            f.emplace_back("conserved_hit_share", conserved_share);
            return f;
        }

        // A feature value as the model and the dump get it: the shortest text that reads back as
        // the same double (fixed six decimals turned small rates into 0).
        inline std::string FeatureString(double value) {
            char buffer[64];
            auto [end, ec] = std::to_chars(buffer, buffer + sizeof(buffer), value);
            return std::string(buffer, end);
        }

        // A sample's depth bin for a model's depth knobs: floor(log10) of the sample's fragments over all its taxa, 2 to
        // 6 (below 100 fragments 2, a million or more 6), as scripts/random_forest_cmdline.py --depth-knobs bins the
        // training samples.
        inline int DepthKnobBin(size_t fragments) {
            int bin = 0;
            for (; fragments >= 10; fragments /= 10) bin++;
            return std::clamp(bin, 2, 6);
        }

        class TaxonFilterForest {
            cpmml::Model m_model{};
            double m_knob = 0.5;
            std::map<int, double> m_depth_knobs;  // see DepthKnob

            mutable std::unordered_map<std::string, std::string> m_sample;


        public:

            TaxonFilterForest(const std::string& path, double knob = 0.5) : m_model(path), m_knob(knob) {}

            // A model already parsed, e.g. with cpmml::Model::from_string from a single-file database.
            TaxonFilterForest(cpmml::Model model, double knob) : m_model(std::move(model)), m_knob(knob) {}

            TaxonFilterForest(const TaxonFilterForest& other) :
                    m_model(other.m_model), m_knob(other.m_knob), m_depth_knobs(other.m_depth_knobs), m_sample() {
            }

            TaxonFilterForest(const TaxonFilterForest&& other) :
                    m_model(other.m_model),
                    m_knob(other.m_knob),
                    m_depth_knobs(other.m_depth_knobs),
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

            // The threshold Pass compares with: a sample's (ProfileWrapper sets each sample's on its thread's copy).
            void SetKnob(double knob) { m_knob = knob; }

            // The model's knobs by sample depth (ParseDepthKnobs), bin (DepthKnobBin) -> knob; empty for most models.
            void SetDepthKnobs(std::map<int, double> knobs) { m_depth_knobs = std::move(knobs); }
            std::map<int, double> const& DepthKnobs() const { return m_depth_knobs; }

            // The model's knob for a sample of `fragments` fragments over all its taxa: the F1-optimal threshold the
            // trainer found for training samples of that depth bin, if the model has one for it.
            std::optional<double> DepthKnob(size_t fragments) const {
                auto const it = m_depth_knobs.find(DepthKnobBin(fragments));
                if (it == m_depth_knobs.end()) return std::nullopt;
                return it->second;
            }

            // The same model with another threshold (e.g. --msa_knob). Scores are cached per taxon,
            // so both give a taxon the same score.
            TaxonFilterForest WithKnob(double knob) const {
                TaxonFilterForest copy(*this);
                copy.m_knob = knob;
                return copy;
            }
        };

        // Whether a taxon's own reads are strong evidence that it is present, whatever its model
        // score: 1x depth or more from its own reads, reads on 90% of its genes, best reads at least
        // 98% identical to the reference, and at most half of its aligned bases from reads of lower
        // identity (a relative's). In the model's test sets on the toy database (2026-09-30), this
        // held for 9% of the present taxa that scored below 0.5, and for none of 14,390 absent ones.
        inline bool StrongOwnEvidence(Taxon const& taxon) {
            return taxon.VerticalCoverage() >= 1.0 && taxon.HitGeneFraction() >= 0.9 &&
                   taxon.TopIdentity() >= 0.98 && taxon.LowIdentityShare() <= 0.5;
        }

        // A placeholder model (scripts/placeholder_models.py) fills a read type's slot in the database
        // until a trained model replaces it: it scores every taxon 0, and protal warns when it loads one.
        inline constexpr std::string_view kPlaceholderModelMarker = "<Annotation>protal:placeholder</Annotation>";

        inline bool IsPlaceholderModel(std::string const& xml) {
            auto const header_end = xml.find("</Header>");
            return header_end != std::string::npos && xml.rfind(kPlaceholderModelMarker, header_end) != std::string::npos;
        }

        // A model's knobs by sample depth, which scripts/random_forest_cmdline.py --depth-knobs writes into the header:
        // <Extension name="protal_depth_knobs" value="2:0.31,3:0.42"/>, depth bin (DepthKnobBin) : knob.
        inline constexpr std::string_view kDepthKnobsExtension = "protal_depth_knobs";

        // Reads the depth knobs from the header of the PMML `xml` into `knobs` (empty without them); an error message if
        // they are malformed: a bin outside 2-6, given twice, or a knob outside 0-1.
        inline std::string ParseDepthKnobs(std::string const& xml, std::map<int, double>& knobs) {
            knobs.clear();
            auto const header_end = xml.find("</Header>");
            if (header_end == std::string::npos) return {};
            std::string const name = "name=\"" + std::string(kDepthKnobsExtension) + "\"";
            for (size_t pos = xml.find("<Extension ", 0); pos < header_end; pos = xml.find("<Extension ", pos + 1)) {
                auto const end = xml.find('>', pos);
                if (end == std::string::npos || end > header_end) break;
                std::string_view const tag(xml.data() + pos, end - pos + 1);
                if (tag.find(" " + name) == std::string_view::npos) continue;
                auto const key = tag.find(" value=\"");
                auto const close = key == std::string_view::npos ? key : tag.find('"', key + 8);
                if (close == std::string_view::npos) return "its depth knobs have no value";
                std::string const value(tag.substr(key + 8, close - key - 8));
                if (value.empty()) return "its depth knobs are malformed (no bins: bin 2-6 : knob 0-1, each bin once)";
                std::istringstream items(value);
                std::string item;
                while (std::getline(items, item, ',')) {
                    auto const colon = item.find(':');
                    int bin = 0;
                    double knob = -1;
                    bool parsed = colon != std::string::npos;
                    if (parsed) {
                        auto const* b = item.data();
                        auto const [bin_end, bin_ec] = std::from_chars(b, b + colon, bin);
                        auto const [knob_end, knob_ec] = std::from_chars(b + colon + 1, b + item.size(), knob);
                        parsed = bin_ec == std::errc() && bin_end == b + colon && knob_ec == std::errc() && knob_end == b + item.size();
                    }
                    if (!parsed || bin < 2 || bin > 6 || !(knob >= 0 && knob <= 1) || knobs.contains(bin)) {
                        knobs.clear();
                        return "its depth knobs are malformed ('" + item + "' in \"" + value + "\": bin 2-6 : knob 0-1, each bin once)";
                    }
                    knobs[bin] = knob;
                }
                return {};
            }
            return {};
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

        // The counts of the taxa's best records that MicrobialProfile::NoteRecord and NoteLinkedRecord collect before the
        // filters, kept apart from the taxa until ApplyRecordEvidence hands them over: one per profile, and one per chunk
        // of a SAM profiled on several threads (Profiler::ProfileSam), whose counts are then added up. They are sums, so
        // the order the chunks are added in does not matter.
        class RecordEvidenceCollector {
        public:
            explicit RecordEvidenceCollector(GenomeLoader& genome_loader) : m_genome_loader(&genome_loader) {}

            void SetGenera(std::shared_ptr<std::vector<uint32_t> const> genera) {
                m_genera = std::move(genera);
            }

            uint32_t GenusOf(uint32_t taxid) const {
                return m_genera && taxid < m_genera->size() ? (*m_genera)[taxid] : 0;
            }

            // A collector without counts that counts as this one does (the same genera).
            RecordEvidenceCollector Empty() const {
                RecordEvidenceCollector empty(*m_genome_loader);
                empty.m_genera = m_genera;
                return empty;
            }

            // See MicrobialProfile::NoteRecord.
            void NoteRecord(uint32_t taxid, SamEntry const& sam) {
                auto& e = m_counts[taxid];
                e.records++;
                e.low_mapq += sam.m_mapq < kLowMapq;
                if (auto const excess = ReadExcess(sam)) e.excess.push_back(*excess);
                if (!m_genera) return;
                auto const genus = GenusOf(taxid);
                bool congener = false, other = false;
                ForEachAlternative(sam.m_alternatives, [&](uint32_t alternative, int more) {
                    if (more > kAlternativeFitEdits) return;
                    if (genus != 0 && GenusOf(alternative) == genus) congener = true;
                    else other = true;
                });
                e.congener_fit += congener;
                e.other_genus_fit += other;
            }

            // See MicrobialProfile::NoteLinkedRecord.
            void NoteLinkedRecord(uint32_t taxid, uint32_t geneid, SamEntry const& sam, size_t link) {
                if (link == SIZE_MAX || m_genome_loader->GetGeneNeighbours().Empty()) return;
                if (link != m_link) {
                    FinishLink();
                    m_link = link;
                }
                bool const reverse = Flag::IsReverseComplement(sam.m_flag);
                // Where on the read it starts: the clips before it in the read's orientation (those at the CIGAR's end if
                // it is reverse); the reader moved the hard clips out of the CIGAR.
                uint32_t const start = (reverse ? sam.m_hard_clip_end : sam.m_hard_clip_start) + Clip(sam.m_cigar, reverse);
                m_link_records.push_back({ taxid, geneid, !reverse, start, start + QueryBases(sam.m_cigar) });
                m_link_paired = Flag::IsPaired(sam.m_flag);
            }

            // See MicrobialProfile::FinishLink.
            void FinishLink() {
                auto const& table = m_genome_loader->GetGeneNeighbours();
                if (m_link_records.size() >= 2) {
                    using gene_neighbours::EndAhead;
                    if (!m_link_paired) {
                        std::stable_sort(m_link_records.begin(), m_link_records.end(),
                                         [](LinkRecord const& a, LinkRecord const& b) { return a.start < b.start; });
                    }
                    uint32_t const max_gap = table.MaxGap() > 0 ? static_cast<uint32_t>(table.MaxGap()) : 3000;
                    for (size_t i = 1; i < m_link_records.size(); i++) {
                        auto const& a = m_link_records[i - 1];
                        auto const& b = m_link_records[i];
                        if (a.gene == b.gene) continue;
                        if (!m_link_paired && b.start > a.end + max_gap) continue;
                        auto const end_a = EndAhead(a.forward);
                        auto const end_b = EndAhead(m_link_paired ? b.forward : !b.forward);
                        auto credit = [&](uint32_t taxid) {
                            auto const verdict = table.Assess(taxid, a.gene, end_a, b.gene, end_b).verdict;
                            auto& e = m_counts[taxid];
                            e.adjacent++;
                            e.adjacent_expected += verdict == gene_neighbours::Verdict::Expected;
                            e.adjacent_unlikely += verdict == gene_neighbours::Verdict::Unlikely;
                        };
                        credit(a.taxid);
                        if (b.taxid != a.taxid) credit(b.taxid);
                    }
                }
                m_link_records.clear();
            }

            // Adds the counts of `other`, whose last link is finished (FinishLink).
            void Add(RecordEvidenceCollector const& other) {
                for (auto const& [taxid, counts] : other.m_counts) m_counts[taxid] += counts;
            }

            // The counts of taxon taxid, or nullptr if it has none.
            RecordEvidence const* Find(uint32_t taxid) const {
                auto const found = m_counts.find(taxid);
                return found == m_counts.end() ? nullptr : &found->second;
            }

        private:
            // The best records of the current link, of every taxon (NoteLinkedRecord): taxon, gene, orientation, where on
            // the read they start and end.
            struct LinkRecord {
                uint32_t taxid;
                uint32_t gene;
                bool forward;
                uint32_t start;
                uint32_t end;
            };
            GenomeLoader* m_genome_loader;
            std::shared_ptr<std::vector<uint32_t> const> m_genera;  // taxid -> genus (0: none)
            std::unordered_map<uint32_t, RecordEvidence> m_counts;
            std::vector<LinkRecord> m_link_records;
            size_t m_link = SIZE_MAX;
            bool m_link_paired = false;  // the current link is a pair's mates, not a long read's genes
        };

        class MicrobialProfile {
        public:
            MicrobialProfile(GenomeLoader& genome_loader) : m_genome_loader(genome_loader), m_evidence(genome_loader) {}

            // Every taxid's genus (GeneraOf), which tells a read's alternatives (ZA) of the same genus from those
            // of another. Without it congener_fit_share and other_genus_fit_share are 0.
            void SetGenera(std::shared_ptr<std::vector<uint32_t> const> genera) {
                m_genera = std::move(genera);
                m_evidence.SetGenera(m_genera);
            }

            uint32_t GenusOf(uint32_t taxid) const {
                return m_genera && taxid < m_genera->size() ? (*m_genera)[taxid] : 0;
            }

            // Counts a read's best record (a mate, a single read, a long read's gene) for its taxon, whatever the
            // filters later make of it: its MAPQ, and whether its alternatives (ZA) hold a congener or a species of
            // another genus within kAlternativeFitEdits. ApplyRecordEvidence hands the counts to the taxa.
            void NoteRecord(uint32_t taxid, SamEntry const& sam) {
                m_evidence.NoteRecord(taxid, sam);
            }

            // Gives every taxon the counts NoteRecord collected for it (the features low_mapq_share,
            // congener_fit_share, other_genus_fit_share); taxa without a counted record keep none.
            void ApplyRecordEvidence() {
                FinishLink();
                for (auto it = m_taxa.begin(); it != m_taxa.end(); ++it) {
                    auto const* found = m_evidence.Find(static_cast<uint32_t>(it->first));
                    it.value().SetRecordEvidence(found ? *found : RecordEvidence{});
                }
            }

            // Notes a read's best record (as NoteRecord: before the filters) of taxon taxid on gene geneid for the adjacency
            // counts, if the database has gene neighbours; link: the read across its records (both mates, a long read's
            // genes). When the next link begins (or ApplyRecordEvidence), FinishLink judges the last one's records.
            void NoteLinkedRecord(uint32_t taxid, uint32_t geneid, SamEntry const& sam, size_t link) {
                m_evidence.NoteLinkedRecord(taxid, geneid, sam, link);
            }

            // The genes of the last link's records next to each other: a pair's mates, which run towards each other, on
            // two genes, or each two consecutive genes of a long read in read order, at most the table's max_gap apart on
            // the read (genes further apart are no neighbours, whatever lies between); the ends that face each other
            // follow from the records' orientations. Whatever taxa the two records are of (a read's gene that its taxon
            // lacks is written on another), each of them is credited with whether its clade has those ends facing each
            // other (gene_neighbours::Table::Assess).
            void FinishLink() {
                m_evidence.FinishLink();
            }

            // The counts NoteRecord and NoteLinkedRecord collect.
            RecordEvidenceCollector& Evidence() {
                return m_evidence;
            }

            RecordEvidenceCollector const& Evidence() const {
                return m_evidence;
            }

            // Drops what the reads added (the taxa and the record evidence), keeping the settings (name, read type,
            // genera, depth identity margin): the profile is as before its first read.
            void ClearReads() {
                m_taxa = TaxonMap();
                m_evidence = m_evidence.Empty();
            }

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

            // Scores every taxon with `filter` on `threads` threads, each with a copy of it (scoring reuses a buffer;
            // copies share the loaded model). The taxa cache their scores and the depth and top identity the scores
            // are computed from, so the outputs find them as if they had scored each taxon themselves.
            void ScoreTaxa(TaxonFilterObj const& filter, size_t threads) {
                std::vector<Taxon const*> taxa;
                taxa.reserve(m_taxa.size());
                for (auto const& [id, taxon] : m_taxa) taxa.push_back(&taxon);
                size_t const blocks = std::min(taxa.size(), 8 * std::max<size_t>(threads, 1));
                sam_chunks::ParallelFor(blocks, threads, [&](size_t b) {
                    TaxonFilterObj model(filter);
                    for (size_t i = b; i < taxa.size(); i += blocks) model.Score(*taxa[i]);
                });
            }

            void AddRead(InternalReadAlignment const& ira, bool unique=true) {
                if (!m_taxa.contains(ira.taxid)) {
                    auto& genome = m_genome_loader.GetGenome(ira.taxid);
                    m_taxa.insert( { ira.taxid, Taxon(genome, &m_genome_loader.GetGeneConservation(), m_genome_loader.ScaleDepthMargin()) } );
                    m_taxa.at(ira.taxid).SetName(std::to_string(ira.taxid));
                }
                auto& taxon = m_taxa.at(ira.taxid);
                taxon.AddHit(ira.geneid, ira.genepos, ira.alignment_ani, unique);
            }

            // What AddSam does with a record (CheckSam).
            enum class SamCheck { kReject, kSkip, kAdd };

            // Whether AddSam takes a record: kReject for one on a gene this database does not have, or reaching past the
            // gene's end (a SAM aligned against another database), which is rejected rather than read out of bounds;
            // kSkip for one on a gene without unique k-mers, which is no evidence of the taxon. Only reads the database.
            SamCheck CheckSam(int taxid, int geneid, SamEntry const& sam) const {
                if (!m_genome_loader.HasGene(taxid, geneid) ||
                    sam.m_pos - 1 + AlignmentLengthRef(sam.m_cigar) > m_genome_loader.GeneLength(taxid, geneid)) {
                    return SamCheck::kReject;
                }
                if (!m_genome_loader.GetGenome(taxid).IsGeneHittable(geneid)) return SamCheck::kSkip;
                return SamCheck::kAdd;
            }

            // The taxon `taxid`, made if the profile has none yet.
            Taxon& TaxonOf(int taxid) {
                if (!m_taxa.contains(taxid)) {
                    auto &genome = m_genome_loader.GetGenome(taxid);

                    m_taxa.insert( { taxid, Taxon(genome, &m_genome_loader.GetGeneConservation(), m_genome_loader.ScaleDepthMargin()) } );
                    m_taxa.at(taxid).SetId(taxid);
                    m_taxa.at(taxid).SetDepthIdentityMargin(m_depth_identity_margin);
                    m_taxa.at(taxid).SetName(std::to_string(taxid));
                }
                return m_taxa.at(taxid);
            }

            // AddSam's work on its taxon: false if the record does not fit its gene. The taxon is not removed then;
            // AddCheckedSam does that.
            static bool AddToTaxon(Taxon& taxon, int geneid, SamEntry const& sam, double score, int read_id, bool no_strain, size_t link) {
                bool const unique = sam.m_mapq > 20;
                ReadEvidence evidence;
                evidence.link = link;
                return taxon.AddSam(geneid, sam, score, unique, read_id, no_strain, evidence);
            }

            // AddSam for a record CheckSam takes (kAdd).
            bool AddCheckedSam(int taxid, int geneid, SamEntry const& sam, double score, int read_id, bool no_strain, size_t link) {
                auto& taxon = TaxonOf(taxid);
                bool success = AddToTaxon(taxon, geneid, sam, score, read_id, no_strain, link);
                // A taxon exists only with at least one read (its means divide by the read count).
                if (!success && taxon.TotalHits() == 0) m_taxa.erase(taxid);
                return success;
            }

            // link: the read across its records (both mates, a long read's genes), for linked_share; SIZE_MAX: none.
            // `unique` is unused: a record with MAPQ above 20 is unique.
            bool AddSam(int taxid, int geneid, SamEntry const& sam, double score, bool unique=true, int read_id=0, bool no_strain=true,
                        size_t link=SIZE_MAX) {
                switch (CheckSam(taxid, geneid, sam)) {
                    case SamCheck::kReject: return false;
                    case SamCheck::kSkip: return true;
                    case SamCheck::kAdd: break;
                }
                return AddCheckedSam(taxid, geneid, sam, score, read_id, no_strain, link);
            }

            void SetName(std::string name) {
                m_name = name;
            }

            // The kind of the sample's reads, whose SNP filters its strain MSA rows take.
            void SetReadType(ReadType type) {
                m_read_type = type;
            }

            ReadType GetReadType() const {
                return m_read_type;
            }

            // The threshold its taxa were reported at (--knob, or its model's knob for the sample's depth), and the one
            // they enter the strain MSAs at (--msa_knob, else the same): ProfileWrapper sets both.
            void SetKnobs(double knob, double msa_knob) {
                m_knob = knob;
                m_msa_knob = msa_knob;
            }

            double Knob() const {
                return m_knob;
            }

            double MSAKnob() const {
                return m_msa_knob;
            }

            // The fragments of all its taxa: the sample's depth, by which a model's depth knobs apply (DepthKnobBin).
            size_t Fragments() const {
                size_t fragments = 0;
                for (auto const& [_, taxon] : m_taxa) fragments += taxon.Fragments();
                return fragments;
            }

            const std::string& GetName() const {
                return m_name;
            }

            std::string& GetName() {
                return m_name;
            }

            // The SNP filters on every gene, on `threads` threads (a taxon each).
            void PostProcessSNPs(size_t min_observations=2, size_t min_observations_fwdrev=2, double min_frequency=0.0, size_t min_avg_quality=15, size_t min_phred_sum=0, bool require_strand=false,
                                 size_t threads=1) {
                std::vector<Taxon*> taxa;
                taxa.reserve(m_taxa.size());
                for (auto it = m_taxa.begin(); it != m_taxa.end(); ++it) taxa.push_back(&it.value());
                sam_chunks::ParallelFor(taxa.size(), threads, [&](size_t i) {
                    auto& taxon = *taxa[i];
                    for (auto& [gid, __] : taxon.GetGenes()) {
                        auto& gene = taxon.GetGenes().at(gid);
                        auto& strain_handler = gene.GetStrainLevel();
                        strain_handler.PostProcess(min_observations, min_observations_fwdrev, min_frequency, min_avg_quality, min_phred_sum, require_strand);
                    }
                    taxon.InvalidateModelScore();
                });
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
            // On `threads` threads, each taxon's lines apart, then written in the order of the taxa: the same file.
            void WriteGeneProfile(taxonomy::IntTaxonomy& taxonomy, TaxonFilterObj const& filter, std::ostream* os, size_t threads=1) {
                *os << "Predicted\tProbability\tTaxID\tLineage\tTaxVCOV\tTaxAbundance\tGeneID\tGeneRefLength\t"
                    << "TotalReads\tTotalMappedLength\tMAPQ\tUniqueMers\tUniqueTwoMers\tUniqueMerReads\tUniqueTwoMerReads\tANI\t"
                    << "VCov\tVCovExp\tHCovExp\tHCovObs\tHCovObsRel\tConsistency\n";

                double const total_vcov = PassingDepth(filter);

                auto const keys = SortedTaxa();
                std::vector<std::string> lines(keys.size());
                sam_chunks::ParallelFor(keys.size(), threads, [&](size_t k) {
                    auto const tax_id = keys[k];
                    TaxonFilterObj model(filter);  // scoring reuses a buffer
                    std::ostringstream out;
                    auto& taxon = m_taxa.at(tax_id);
                    double probability = model.Score(taxon);
                    bool prediction = probability >= model.GetKnob();

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

                        out
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
                    lines[k] = out.str();
                });
                for (auto const& text : lines) *os << text;
            }

            // The profile (taxa that pass the model), .profile.log (every taxon with its call, features
            // and per-gene coverage; one column per gene id of the database, so every file has the same
            // columns) and .gene.log (statistics per gene hit). On `threads` threads, each taxon's lines
            // apart, then written in the order of the taxa: the same files.
            void WriteSparseProfile(taxonomy::IntTaxonomy& taxonomy, TaxonFilterObj const& filter, std::ostream &os_filtered=std::cout, std::ostream* os_total=nullptr, std::ostream* os_genes=nullptr,
                                    size_t threads=1) {
                size_t max_gene_id = 0;
                for (auto& [key, genome] : m_genome_loader.GetGenomeMap()) {
                    max_gene_id = std::max(max_gene_id, genome.GetGeneList().size());
                }

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

                struct Lines {
                    std::string filtered, total, genes;
                    bool prediction = false;
                };
                auto const keys = SortedTaxa();
                std::vector<Lines> lines(keys.size());
                sam_chunks::ParallelFor(keys.size(), threads, [&](size_t k) {
                    auto const key = keys[k];
                    TaxonFilterObj model(filter);  // scoring reuses a buffer
                    std::ostringstream filtered, total, genes;
                    std::vector<size_t> gene_covs(max_gene_id + 1, 0);
                    std::vector<double> gene_cov_ratios(max_gene_id + 1, 0.0);

                    // this is necessary as taxon cannot be constant
                    auto& taxon = m_taxa.at(key);
                    double probability = model.Score(taxon);
                    bool prediction = probability >= model.GetKnob();
                    lines[k].prediction = prediction;
                    double const vcov = taxon.VerticalCoverage();
                    double const abundance = prediction && total_vcov > 0 ? vcov / total_vcov : 0;

                    for (auto id : taxon.SortedGeneIds()) {
                        auto& gene = taxon.GetGenes().at(id);
                        size_t cov = gene.GetStrainLevel().GetSequenceRangeHandler().CoveredPortion();
                        double ratio = static_cast<double>(cov) / gene.m_gene_length;
                        gene_covs[id] = cov;
                        gene_cov_ratios[id] = ratio;

                        if (os_genes) {
                            genes << m_name << '\t';
                            genes << key << '\t';
                            genes << taxonomy.LineageStr(key) << '\t';
                            genes << taxon.GetName() << '\t';
                            genes << id << '\t';
                            genes << "Gene" + std::to_string(id) << '\t';
                            genes << gene.GetStatisticsString() << '\n';
                        }
                    }

                    std::string gene_covs_str, gene_cov_ratios_str;
                    for (auto val : gene_covs) gene_covs_str += "\t" + std::to_string(val);
                    for (auto val : gene_cov_ratios) gene_cov_ratios_str += "\t" + std::to_string(val);

                    double mean_gene_covs = std::accumulate(gene_covs.begin(), gene_covs.end(), size_t{0}) /
                                            static_cast<double>(taxon.GetGenes().size());
                    double mean_gene_cov_ratios = std::accumulate(gene_cov_ratios.begin(), gene_cov_ratios.end(), 0.0) /
                                            static_cast<double>(taxon.GetGenes().size());

                    auto const& node = taxonomy.Get(key);

                    if (prediction) {
                        filtered << node.rep_genome << '\t' << taxonomy.LineageStr(key) << '\t' << abundance << '\n';
                    }

                    if (os_total) {
                        total << (prediction ? "1" : "0") << "\t" << probability << "\t" << node.rep_genome << '\t' << taxonomy.LineageStr(key) << '\t' << abundance;
                        total << '\t' << taxon.VCovStdDev();
                        total << '\t' << taxon.GetGeneVariance();
                        total << '\t' << taxon.GetGeneVariance(5);
                        total << '\t' << taxon.ToString(taxonomy);
                        total << '\t' << vcov << '\t' << taxon.LowIdentityShare();
                        total << '\t' << mean_gene_covs << '\t' << mean_gene_cov_ratios << gene_covs_str;
                        total << gene_cov_ratios_str << '\n';
                    }
                    lines[k].filtered = filtered.str();
                    lines[k].total = total.str();
                    lines[k].genes = genes.str();
                });

                bool one_pass = false;
                for (auto const& taxon_lines : lines) {
                    os_filtered << taxon_lines.filtered;
                    if (os_total) *os_total << taxon_lines.total;
                    if (os_genes) *os_genes << taxon_lines.genes;
                    one_pass |= taxon_lines.prediction;
                }
                os_filtered.flush();

                if (!one_pass) {
                    std::cout << "No taxon passes the model in sample " << m_name << std::endl;
                }
            }


        private:
            std::string m_name;
            ReadType m_read_type = ReadType::Paired;
            double m_knob = 0.5;      // see Knob
            double m_msa_knob = 0.5;  // see MSAKnob
            mutable TaxonMap m_taxa;
            GenomeLoader &m_genome_loader;
            double m_depth_identity_margin = 1;
            std::shared_ptr<std::vector<uint32_t> const> m_genera;  // taxid -> genus (0: none); see SetGenera
            RecordEvidenceCollector m_evidence;  // NoteRecord, NoteLinkedRecord
        };

        // Every taxon's genus (its taxid; 0 for a taxon without one), for MicrobialProfile::SetGenera.
        inline std::shared_ptr<std::vector<uint32_t> const> GeneraOf(taxonomy::IntTaxonomy const& taxonomy) {
            size_t max_id = 0;
            for (auto const& [id, node] : taxonomy.map) max_id = std::max<size_t>(max_id, static_cast<size_t>(id));
            auto genera = std::make_shared<std::vector<uint32_t>>(max_id + 1, 0);
            for (auto const& [id, node] : taxonomy.map) {
                int t = static_cast<int>(id);
                for (int steps = 0; steps < 64 && taxonomy.map.contains(t); steps++) {
                    auto const& n = taxonomy.map.at(t);
                    if (n.rank == "genus") {
                        (*genera)[id] = static_cast<uint32_t>(n.id);
                        break;
                    }
                    if (n.parent_id < 0 || n.parent_id == t) break;
                    t = n.parent_id;
                }
            }
            return genera;
        }


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
            void CheckReference(std::string const& line) const {
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

            // SamReader's counts of a SAM's records, which ReportSamRecords reports; added up over the chunks of a SAM
            // profiled on several threads.
            struct SamRecordCounts {
                size_t records = 0;
                size_t without_tags = 0;
                size_t primary = 0;
                size_t primary_without_alternatives = 0;
                std::map<std::string, size_t> skipped;

                void Add(SamReader const& reader) {
                    records += reader.Records();
                    without_tags += reader.RecordsWithoutTags();
                    primary += reader.PrimaryRecords();
                    primary_without_alternatives += reader.PrimaryRecordsWithoutAlternatives();
                    for (auto const& [reason, count] : reader.Skipped()) skipped[reason] += count;
                }

                void Add(SamRecordCounts const& other) {
                    records += other.records;
                    without_tags += other.without_tags;
                    primary += other.primary;
                    primary_without_alternatives += other.primary_without_alternatives;
                    for (auto const& [reason, count] : other.skipped) skipped[reason] += count;
                }
            };

            // The records skipped, and a SAM without usable records or without protal's tags, once it is read.
            static void ReportSamRecords(std::string const& file_path, SamRecordCounts const& counts) {
                for (auto const& [reason, count] : counts.skipped) {
                    std::cerr << file_path << ": skipped " << count << " record(s): " << reason << std::endl;
                }
                if (counts.records == 0) {
                    std::cerr << file_path << " contains no usable alignments" << std::endl;
                } else if (counts.without_tags > 0) {
                    std::cerr << "Warning: " << counts.without_tags << " of " << counts.records << " records in "
                              << file_path << " have no ZU tag (protal's unique k-mer count). A SAM file not "
                              << "written by protal lacks it, and the model then rejects most taxa." << std::endl;
                } else if (counts.primary > 0 && counts.primary_without_alternatives == counts.primary) {
                    std::cerr << "Warning: no record in " << file_path << " has a ZA tag (the read's alternative alignments "
                              << "to other taxa), which older protal versions do not write: every taxon's congener_fit_share "
                              << "and other_genus_fit_share is 0. Align the reads again (--force) for a model that uses them."
                              << std::endl;
                }
            }

            // Hands the reader's records to `on_group` as groups of candidate alignments: adjacent records with one
            // QNAME are one read's candidates, and a supplementary record (0x800) starts a group of its own, one part
            // of a long read. When the records end, the last group is left in `group`.
            template<typename OnGroup>
            static void CollectGroups(SamReader& reader, std::vector<AlignmentPair>& group, OnGroup&& on_group) {
                SamEntry sam1;
                SamEntry sam2;
                bool has_sam1 = false, has_sam2 = false;
                while (reader.Next(sam1, sam2, has_sam1, has_sam2)) {
                    AlignmentPair pair(
                            has_sam1 ? std::optional<SamEntry>{ std::move(sam1) } : std::optional<SamEntry>{},
                            has_sam2 ? std::optional<SamEntry>{ std::move(sam2) } : std::optional<SamEntry>{});
                    if (!group.empty() && (!SameRead(pair, group.front()) || Flag::IsSupplementaryAlignment(pair.Any().m_flag))) {
                        on_group(group);
                        group.clear();
                    }
                    group.emplace_back(std::move(pair));
                }
            }

            // Reads a SAM file (plain, gzip or zstd; SamInput) and hands each read's group of candidate
            // alignments to `on_group` (CollectGroups). Returns an error message if the file cannot be read,
            // is truncated or was aligned against another database (CheckReference); a SAM without
            // alignments is not an error.
            template<typename OnGroup>
            std::string ReadSamGroups(std::string const& file_path, OnGroup&& on_group) {
                if (std::filesystem::exists(file_path) && std::filesystem::file_size(file_path) == 0) {
                    return "the file is empty (not even a SAM header)";
                }
                SamInput input(file_path);
                if (!input.IsOpen()) return "cannot open the file";
                // A file cut at a block or frame boundary lacks its format's end marker.
                if (!input.Problem().empty()) return "the file is truncated or corrupt (" + input.Problem() + ")";
                std::istream& file = input.Stream();
                SamReader reader(file, [this](std::string const& line) { CheckReference(line); });
                // The gzip readers and zstd read a truncated or corrupt file as one that ends early.
                auto truncated = [&input]() {
                    return "the file is truncated or corrupt (" + input.ReadError() + ")";
                };

                std::vector<AlignmentPair> group;
                try {
                    CollectGroups(reader, group, on_group);
                } catch (SamFormatError const& e) {
                    // A truncated file's last line is cut short, too: the truncation is the cause.
                    if (input.ReadFailed()) return truncated();
                    return e.what();
                }
                if (input.ReadFailed()) return truncated();
                if (!group.empty()) on_group(group);

                SamRecordCounts counts;
                counts.Add(reader);
                ReportSamRecords(file_path, counts);
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

                    auto const ref = gene.Sequence();
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
                        auto& gene = m_genome_loader.GetGenome(tid).GetGeneOMP(gid);
                        auto const ref = gene.Sequence();
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



            // A record that ProcessMAPQ adds to a taxon (MicrobialProfile::AddCheckedSam), with what ProfileSam's chunks
            // need: the read's number and link in the file, and whether the record fit its gene.
            struct SamAddition {
                int taxid = 0;
                int geneid = 0;
                SamEntry const* sam = nullptr;
                double ani = 0;
                int read_id = 0;
                size_t link = SIZE_MAX;
                bool ok = true;
            };

        private:
            std::vector<SamAddition> m_additions;  // ProcessMAPQ's

        public:
            // ProcessMAPQ before it adds anything to a taxon: the read's evidence (NoteRecord, NoteLinkedRecord) into
            // `evidence`, and its records that pass the MAPQ and length filters into `additions` if CheckSam takes them.
            // Returns false if CheckSam rejects one (the read does not fit the database). Reads `profile` only, so that
            // the chunks of a SAM can be prepared on several threads.
            bool PrepareMAPQ(MicrobialProfile const& profile, AlignmentPair& ap, size_t link, RecordEvidenceCollector& evidence,
                             CigarInfo& info1, CigarInfo& info2, std::vector<SamAddition>& additions) const {
                bool valid_sam = true;

                bool take_first = false;
                bool take_second = false;
                if (ap.HasFirst()) CompressedCigarInfo(ap.First().m_cigar, info1);
                if (ap.HasSecond()) CompressedCigarInfo(ap.Second().m_cigar, info2);

                // Every best record counts for its taxon's MAPQ and alternative evidence, also one the MAPQ filter
                // below leaves out: a read that fits another taxon as well has MAPQ near 0.
                if (ap.HasFirst()) evidence.NoteRecord(ExtractTaxidGeneid(ap.First().m_rname).first, ap.First());
                if (ap.HasSecond()) evidence.NoteRecord(ExtractTaxidGeneid(ap.Second().m_rname).first, ap.Second());
                for (auto* sam : { ap.HasFirst() ? &ap.First() : nullptr, ap.HasSecond() ? &ap.Second() : nullptr }) {
                    if (!sam) continue;
                    auto const [taxid, geneid] = ExtractTaxidGeneid(sam->m_rname);
                    evidence.NoteLinkedRecord(static_cast<uint32_t>(taxid), static_cast<uint32_t>(geneid), *sam, link);
                }

                // MAPQ is judged per mate: the mates of a pair aligned together share one MAPQ, those
                // of a fragment split over two genes each have their own.
                if (ap.HasFirst() && ap.First().m_mapq >= m_min_mapq && info1.clipped_alignment_length > m_min_alignment_length) {
                    take_first = true;
                }
                if (ap.HasSecond() && ap.Second().m_mapq >= m_min_mapq && info2.clipped_alignment_length > m_min_alignment_length) {
                    take_second = true;
                }
                if (!take_first && !take_second) return true;

                auto take = [&](SamEntry& sam, CigarInfo const& info) {
                    auto [tid, geneid] = ExtractTaxidGeneid(sam.m_rname);
                    int const taxid = static_cast<int>(tid), gene = static_cast<int>(geneid);
                    switch (profile.CheckSam(taxid, gene, sam)) {
                        case MicrobialProfile::SamCheck::kReject: valid_sam = false; break;
                        case MicrobialProfile::SamCheck::kSkip: break;
                        case MicrobialProfile::SamCheck::kAdd: {
                            SamAddition addition;
                            addition.taxid = taxid;
                            addition.geneid = gene;
                            addition.sam = &sam;
                            addition.ani = info.Ani();
                            additions.push_back(addition);
                            break;
                        }
                    }
                };
                if (take_first) take(ap.First(), info1);
                if (take_second) take(ap.Second(), info2);
                return valid_sam;
            }

            // Adds a read's best records to the profile: false if one does not fit the database (a gene it lacks, a
            // position past a gene's end, or bases that differ from the gene). link: the read across its groups (a long
            // read's genes; see ProfileSam).
            bool ProcessMAPQ(MicrobialProfile& profile, AlignmentPair& ap, int read_id=0, size_t link=SIZE_MAX) {
                m_additions.clear();
                bool valid_sam = PrepareMAPQ(profile, ap, link, profile.Evidence(), m_info1, m_info2, m_additions);
                for (auto const& addition : m_additions) {
                    valid_sam &= profile.AddCheckedSam(addition.taxid, addition.geneid, *addition.sam, addition.ani, read_id, kNoStrain, link);
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

            // The bytes of SAM text per chunk when a SAM is profiled on several threads (ProfileSam); 0: by the
            // number of threads. Small chunks let the tests cut a small SAM into many.
            void SetChunkBytes(size_t bytes) {
                m_chunk_bytes = bytes;
            }

        private:
            size_t m_chunk_bytes = 0;

            // A read of a chunk (ProfileSamParallel): its group of candidates, which holds its records, the one of
            // them that counts (BestOfGroup), its link within the chunk, its records to add (SamAdditions
            // [first_addition, first_addition + additions) of the chunk) and whether PrepareMAPQ found it valid.
            struct ChunkRead {
                std::vector<AlignmentPair> group;
                size_t pair = 0;
                size_t link = 0;
                size_t first_addition = 0;
                size_t additions = 0;
                bool valid = true;
            };

            // What ParseChunk makes of a chunk. The reads are in a deque: the additions point to their records.
            struct ParsedChunk {
                explicit ParsedChunk(RecordEvidenceCollector evidence) : evidence(std::move(evidence)) {}
                std::deque<ChunkRead> reads;
                std::vector<SamAddition> additions;
                RecordEvidenceCollector evidence;
                SamRecordCounts counts;
                size_t links = 0;  // reads with their own QNAME
                size_t lines = 0;  // the line the chunk ends with
                bool started = false;  // a record was read (a read begun)
                std::string error;  // a SamFormatError; the reads before it count
            };

            // ProfileSamParallel's work on a chunk, on any thread: its reads (CollectGroups) with their links within the
            // chunk and PrepareMAPQ. A chunk starts with a read of its own, so a read is numbered and linked in the
            // chunk as ProfileSam numbers and links it in the file, less the reads and links of the chunks before.
            // `last_read`: the chunk's last read counts once the chunk ends; not if the file is truncated there,
            // as ReadSamGroups leaves it out then.
            void ParseChunk(MicrobialProfile const& profile, sam_chunks::Chunk& chunk, ParsedChunk& out, bool last_read) const {
                sam_chunks::TextStreambuf buffer(chunk.text);
                std::istream is(&buffer);
                SamReader reader(is, [this](std::string const& line) { CheckReference(line); }, chunk.first_line);
                CigarInfo info1, info2;
                size_t link = 0;
                std::string last_qname;
                auto on_group = [&](std::vector<AlignmentPair>& group) {
                    auto& read = out.reads.emplace_back();
                    read.group = std::move(group);
                    auto& pair = read.group.size() == 1 ? read.group.front() : BestOfGroup(read.group);
                    read.pair = static_cast<size_t>(&pair - read.group.data());
                    auto const& qname = pair.Any().m_qname;
                    if (out.reads.size() > 1 && qname != last_qname) link++;
                    last_qname = qname;
                    read.link = link;
                    read.first_addition = out.additions.size();
                    read.valid = PrepareMAPQ(profile, pair, link, out.evidence, info1, info2, out.additions);
                    read.additions = out.additions.size() - read.first_addition;
                };
                std::vector<AlignmentPair> group;
                try {
                    CollectGroups(reader, group, on_group);
                    if (last_read && !group.empty()) on_group(group);
                } catch (SamFormatError const& e) {
                    out.error = e.what();
                }
                out.started = !out.reads.empty() || !group.empty();
                out.evidence.FinishLink();
                out.links = out.reads.empty() ? 0 : link + 1;
                out.lines = reader.Lines();
                out.counts.Add(reader);
                std::string().swap(chunk.text);  // the records hold their own copies
            }

            // The records of one taxon in a wave, in file order (PlanChunks), which one thread adds to it (RunJob).
            struct TaxonJob {
                Taxon* taxon = nullptr;
                std::vector<SamAddition*> additions;
            };

            // What has to be done in file order when the parsed chunks [0, used) of a wave are added to the profile: each
            // read's number and link in the file (read_offset, link_offset: those of the chunks before), and every taxon
            // made as ProfileSam makes them. Returns each taxon's records, the taxa with most first. No taxon is made or
            // removed until the jobs are done, so the references hold.
            std::vector<TaxonJob> PlanChunks(MicrobialProfile& profile, std::vector<ParsedChunk>& parsed, size_t used,
                                             size_t& read_offset, size_t& link_offset) const {
                std::unordered_map<uint32_t, size_t> job_of;
                std::vector<std::pair<uint32_t, std::vector<SamAddition*>>> records;
                for (size_t c = 0; c < used; c++) {
                    auto& chunk = parsed[c];
                    for (size_t r = 0; r < chunk.reads.size(); r++) {
                        auto const& read = chunk.reads[r];
                        for (size_t k = read.first_addition; k < read.first_addition + read.additions; k++) {
                            auto& addition = chunk.additions[k];
                            addition.read_id = static_cast<int>(read_offset + r);
                            addition.link = link_offset + read.link;
                            profile.TaxonOf(addition.taxid);
                            auto const [it, made] = job_of.try_emplace(static_cast<uint32_t>(addition.taxid), records.size());
                            if (made) records.emplace_back(static_cast<uint32_t>(addition.taxid), std::vector<SamAddition*>{});
                            records[it->second].second.push_back(&addition);
                        }
                    }
                    read_offset += chunk.reads.size();
                    link_offset += chunk.links;
                }
                std::sort(records.begin(), records.end(), [](auto const& a, auto const& b) { return a.second.size() > b.second.size(); });
                std::vector<TaxonJob> jobs;
                jobs.reserve(records.size());
                for (auto& [taxid, additions] : records) jobs.push_back({ &profile.GetTaxa().at(taxid), std::move(additions) });
                return jobs;
            }

            // Adds a taxon's records of a wave to it. False if a record that does not fit its gene leaves the taxon
            // without records: ProfileSam would remove the taxon (and maybe make it anew) at that point.
            static bool RunJob(TaxonJob& job) {
                for (SamAddition* addition : job.additions) {
                    addition->ok = MicrobialProfile::AddToTaxon(*job.taxon, addition->geneid, *addition->sam, addition->ani,
                                                                 addition->read_id, kNoStrain, addition->link);
                    if (!addition->ok && job.taxon->TotalHits() == 0) return false;
                }
                return true;
            }

            // ProfileSam on `threads` threads. The SAM is read in a thread of its own and cut into chunks of whole
            // reads (SamChunks.h), which are parsed and prepared (ParseChunk) in waves. Of a wave's records, what
            // depends on their order is done in file order (PlanChunks), and then the taxa add their records on all
            // threads (RunJob) while the next wave is parsed: every taxon gets the same records in the same order as
            // from ProfileSam on one thread, so the profile is the same. The rejected reads' records go to `rejected`,
            // in file order. Sets `serial` and stops if a record that does not fit its gene would leave its taxon
            // without records (RunJob); ProfileSam then profiles the file on one thread.
            std::string ProfileSamParallel(std::string const& file_path, MicrobialProfile& profile, std::string& rejected,
                                           size_t threads, bool& serial) {
                if (std::filesystem::exists(file_path) && std::filesystem::file_size(file_path) == 0) {
                    return "the file is empty (not even a SAM header)";
                }
                SamInput input(file_path);
                if (!input.IsOpen()) return "cannot open the file";
                // A file cut at a block or frame boundary lacks its format's end marker.
                if (!input.Problem().empty()) return "the file is truncated or corrupt (" + input.Problem() + ")";
                // The gzip readers and zstd read a truncated or corrupt file as one that ends early.
                auto truncated = [&input]() {
                    return "the file is truncated or corrupt (" + input.ReadError() + ")";
                };

                // Chunks of 1 MB, four per thread in a wave (up to 128 MB: a wave's records are held in memory), so that
                // threads of different speed share a wave's work.
                size_t const bytes = m_chunk_bytes > 0 ? m_chunk_bytes : size_t{1} << 20;
                size_t const per_wave = std::max(threads, std::min<size_t>(4 * threads, 128));
                sam_chunks::ChunkReader reader(input.Stream(), bytes, 2 * per_wave);
                SamRecordCounts counts;
                size_t read_offset = 0, link_offset = 0, lines = 0;
                std::string error;
                bool error_at_end = false;  // in the last chunk
                bool failed = false;        // the stream could not be read to its end
                bool end = false;           // the last chunk is taken
                // Where the last read's rejected records begin in `rejected`, if it was rejected: when reading stops at
                // an error or at a truncation, ReadSamGroups leaves out the read it has begun, and that is the last read
                // of an earlier chunk if the chunk with the error began none.
                size_t last_mark = 0;
                bool last_rejected = false;

                std::vector<sam_chunks::Chunk> wave, next;     // this wave's chunks, the next one's
                std::vector<ParsedChunk> parsed, upcoming, done;  // and their parsed reads; the last wave's, to free
                auto fetch = [&]() {
                    next.clear();
                    upcoming.clear();
                    sam_chunks::Chunk chunk;
                    while (next.size() < per_wave && !end && reader.Next(chunk)) {
                        end = chunk.last;
                        next.push_back(std::move(chunk));
                    }
                    // The reader has read all of the stream once it hands over the last chunk.
                    if (!next.empty() && next.back().last) failed = reader.Bad() || input.ReadFailed();
                    upcoming.reserve(next.size());
                    for (size_t i = 0; i < next.size(); i++) upcoming.emplace_back(profile.Evidence().Empty());
                };
                auto parse = [&](size_t i) { ParseChunk(profile, next[i], upcoming[i], !(next[i].last && failed)); };
                auto release = [](ParsedChunk& chunk) {
                    std::deque<ChunkRead>().swap(chunk.reads);
                    std::vector<SamAddition>().swap(chunk.additions);
                };

                fetch();
                sam_chunks::ParallelFor(next.size(), threads, parse);
                while (!next.empty()) {
                    std::swap(wave, next);
                    std::swap(parsed, upcoming);
                    // The chunks up to the first with an error count, as the reads before an error do in ProfileSam.
                    size_t used = wave.size();
                    for (size_t i = 0; i < wave.size(); i++) {
                        if (parsed[i].error.empty()) continue;
                        used = i + 1;
                        error = parsed[i].error;
                        error_at_end = wave[i].last;
                        break;
                    }
                    for (size_t i = 0; i < used; i++) {
                        counts.Add(parsed[i].counts);
                        profile.Evidence().Add(parsed[i].evidence);
                        lines = parsed[i].lines;
                    }
                    auto jobs = PlanChunks(profile, parsed, used, read_offset, link_offset);

                    // The taxa add this wave's records while the next wave is parsed and the last one is freed.
                    if (error.empty()) fetch();
                    else {
                        next.clear();
                        upcoming.clear();
                    }
                    std::atomic<bool> emptied{ false };
                    size_t const n_jobs = jobs.size(), n_parse = next.size();
                    sam_chunks::ParallelFor(n_jobs + n_parse + done.size(), threads, [&](size_t t) {
                        if (t < n_jobs) {
                            if (!emptied && !RunJob(jobs[t])) emptied = true;
                        } else if (t < n_jobs + n_parse) {
                            parse(t - n_jobs);
                        } else {
                            release(done[t - n_jobs - n_parse]);
                        }
                    });
                    done.clear();
                    if (emptied) {
                        serial = true;
                        break;
                    }

                    for (size_t i = 0; i < used; i++) {
                        auto& chunk = parsed[i];
                        for (auto& read : chunk.reads) {
                            bool valid = read.valid;
                            for (size_t k = read.first_addition; k < read.first_addition + read.additions; k++) valid &= chunk.additions[k].ok;
                            last_mark = rejected.size();
                            last_rejected = !valid;
                            if (valid) continue;
                            m_rejected_reads++;
                            auto& pair = read.group[read.pair];
                            if (pair.first.has_value()) rejected += pair.first.value().ToString() + '\n';
                            if (pair.second.has_value()) rejected += pair.second.value().ToString() + '\n';
                        }
                    }
                    size_t const stop = !error.empty() ? used - 1 : wave.back().last && failed ? wave.size() - 1 : SIZE_MAX;
                    if (stop != SIZE_MAX && !parsed[stop].started && last_rejected) {
                        rejected.resize(last_mark);
                        m_rejected_reads--;
                    }
                    std::swap(done, parsed);
                    if (!error.empty()) break;
                }
                // What is left of the waves, freed on all threads.
                for (auto* waves : { &done, &parsed, &upcoming }) {
                    sam_chunks::ParallelFor(waves->size(), threads, [&](size_t i) { release((*waves)[i]); });
                }
                reader.Stop();
                if (serial) return {};
                if (!error.empty()) {
                    // A truncated file's last line is cut short, too: the truncation is the cause.
                    if (error_at_end && input.ReadFailed()) return truncated();
                    return error;
                }
                if (failed) return input.ReadFailed() ? truncated() : "read error after line " + std::to_string(lines);
                m_reads = read_offset;
                ReportSamRecords(file_path, counts);
                return {};
            }

        public:
            using OptionalRefOstream = std::optional<std::reference_wrapper<std::ostream>>;
            // Profiles a SAM file one read at a time, without holding the file in memory: a read with
            // one candidate alignment, or the best of a multi-mapped read's (BestOfGroup), is added to
            // `profile`. Records the profile rejects (genes outside the database, alignments that do
            // not match their gene) go to `erroneous_sam_out`. Returns an error message as ReadSamGroups;
            // the profile is then incomplete. On `threads` threads (ProfileSamParallel) the profile is the
            // same, and so are an error and the rejected records.
            std::string ProfileSam(std::string const& file_path, MicrobialProfile& profile, OptionalRefOstream erroneous_sam_out={}, size_t snp_min_cov=2, size_t snp_min_obs_fwdrev=2, double snp_min_af=0.0, size_t snp_min_mean_qual=15, size_t snp_min_phred_sum=0, bool snp_require_strand=false,
                                   size_t threads=1) {
                profile.SetDepthIdentityMargin(m_depth_identity_margin);
                m_rejected_reads = 0;
                if (threads > 1) {
                    std::string rejected;
                    bool serial = false;
                    auto const error = ProfileSamParallel(file_path, profile, rejected, threads, serial);
                    if (!serial) {
                        if (erroneous_sam_out.has_value()) erroneous_sam_out.value().get() << rejected;
                        if (!error.empty()) return error;
                        profile.ApplyRecordEvidence();
                        m_post_process_bm.Start();
                        profile.PostProcessSNPs(snp_min_cov, snp_min_obs_fwdrev, snp_min_af, snp_min_mean_qual, snp_min_phred_sum, snp_require_strand, threads);
                        m_post_process_bm.Stop();
                        return {};
                    }
                    // Only the serial profile removes a taxon whose records do not fit, at the point it does: again on one thread.
                    profile.ClearReads();
                    m_rejected_reads = 0;
                }
                size_t read_id = 0;
                // One link per read: a pair's group, or all groups of a long read (its genes, one supplementary
                // record each, which follow each other with the read's name).
                size_t link = 0;
                std::string last_qname;
                auto error = ReadSamGroups(file_path, [&](std::vector<AlignmentPair>& group) {
                    auto& pair = group.size() == 1 ? group.front() : BestOfGroup(group);
                    auto const& qname = pair.Any().m_qname;
                    if (read_id > 0 && qname != last_qname) link++;
                    last_qname = qname;
                    bool const valid = ProcessMAPQ(profile, pair, read_id++, link);
                    m_rejected_reads += !valid;
                    if (!valid && erroneous_sam_out.has_value()) {
                        auto& os = erroneous_sam_out.value().get();
                        if (pair.first.has_value()) os << pair.first.value().ToString() << '\n';
                        if (pair.second.has_value()) os << pair.second.value().ToString() << '\n';
                    }
                });
                if (!error.empty()) return error;

                m_reads = read_id;
                profile.ApplyRecordEvidence();
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
