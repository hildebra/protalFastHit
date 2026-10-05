//
// Created by fritsche on 21/02/23.
//

#pragma once

#include <omp.h>

#include "Matrix.h"
#include <optional>
#include <unordered_set>
#include "VariantHandler.h"
#include "SequenceRangeHandler.h"
#include <string_view>

namespace protal {
    // Alleles ranked as calls: by observations, then by quality sum, then by the allele itself (Variant::AlleleBefore),
    // not by the order the reads first showed them in.
    static bool RanksBefore(Variant const& a, Variant const& b) {
        if (a.Observations() != b.Observations()) return a.Observations() > b.Observations();
        if (a.QualitySum() != b.QualitySum()) return a.QualitySum() > b.QualitySum();
        return a.AlleleBefore(b);
    }

    // The call at a site from its base alleles (the reference, SNPs, a deletion; not insertions,
    // which come before the base): the first that `passes`, ranked by RanksBefore. With none
    // passing, the top allele is returned and `pass` is false: the site is written as N.
    template<typename Passes>
    static Variant const* BaseCall(VariantBin const& bin, Passes&& passes, bool& pass) {
        std::vector<Variant const*> alleles;
        for (auto const& v : bin) if (!v.IsINS()) alleles.push_back(&v);
        pass = false;
        if (alleles.empty()) return nullptr;
        std::stable_sort(alleles.begin(), alleles.end(), [](Variant const* a, Variant const* b) { return RanksBefore(*a, *b); });
        for (auto const* a : alleles) {
            if (passes(*a)) {
                pass = true;
                return a;
            }
        }
        return alleles.front();
    }

    // The insertion written before a site's base: the top passing one, or none.
    template<typename Passes>
    static Variant const* InsertionCall(VariantBin const& bin, Passes&& passes) {
        Variant const* best = nullptr;
        for (auto const& v : bin) {
            if (v.IsINS() && passes(v) && (!best || RanksBefore(v, *best))) best = &v;
        }
        return best;
    }


    struct SharedAlignmentRegion;
    using VariantVec = std::vector<VariantBin>;

    // This should be in alignment utils but there are weird circular dependencies
    static size_t AlignmentLengthRef(const std::string& cigar) {
        if (cigar.empty()) return 0;

        int score = 0;
        int count = 0;
        char c;
        int digit_start = -1;

        size_t ins = 0;
        size_t del = 0;
        size_t soft = 0;
        size_t mismatch = 0;
        size_t match = 0;

        for (auto i = 0; i < cigar.length(); i++) {
            if (!std::isdigit(cigar[i])) {
                count = CigarCount(cigar, digit_start, static_cast<size_t>(i));
                c = cigar[i];
                digit_start = -1;
                if (c == 'X') mismatch += count;
                if (c == 'M') match += count;
                if (c == 'I') ins += count;
                if (c == 'D') del += count;
            } else if (digit_start == -1) {
                digit_start = i;
            }
        }
        return del + mismatch + match;
    }

    // Identity of an alignment: matches over aligned columns (M, X, I and D; protal writes M only
    // for exact matches), and the reference length it covers (M, X and D).
    static std::pair<double, size_t> AlignmentIdentity(const std::string& cigar) {
        size_t matches = 0, differences = 0, ref_length = 0;
        size_t count = 0;
        for (char c : cigar) {
            if (std::isdigit(static_cast<unsigned char>(c))) {
                count = count * 10 + (c - '0');
                continue;
            }
            if (c == 'M') matches += count;
            if (c == 'X' || c == 'I' || c == 'D') differences += count;
            if (c == 'M' || c == 'X' || c == 'D') ref_length += count;
            count = 0;
        }
        double const columns = static_cast<double>(matches + differences);
        return { columns > 0 ? matches / columns : 0.0, ref_length };
    }

    class StrainLevelContainer {
        friend SharedAlignmentRegion;
        // Needs two things, variant handler which stores and handles SNPs and INDELS
        // And sequence range handler which stores which regions are covered by the alignment
        VariantHandler m_variant_handler;
        SequenceRangeHandler m_sequence_range_handler;

        std::vector<SNP> m_snp_tmp;
        const Gene& m_reference;
        // The fragment of the last alignment added and the reference interval its reads cover here:
        // its other mate skips that interval, so that a fragment counts once.
        size_t m_last_read_id = SIZE_MAX;
        std::pair<size_t, size_t> m_last_interval{ 0, 0 };

        void AddReadRange(size_t start, size_t end, size_t read_id, bool forward, uint8_t divergence) {
            if (end <= start) return;
            SequenceRange range(start, end);
            ReadInfo rinfo;
            rinfo.read_id = read_id;
            rinfo.length = static_cast<uint32_t>(end - start);
            rinfo.start = static_cast<uint32_t>(start);
            rinfo.forward = forward;
            rinfo.divergence = divergence;
            range.AddReadInfo(rinfo);
            m_sequence_range_handler.Merge(std::move(range));
        }

    public:
        StrainLevelContainer(Gene const& reference) :
                m_reference(reference), m_variant_handler(reference) {
        };

        const VariantHandler& GetVariantHandler() const {
            return m_variant_handler;
        }

        // Drops the gene's variants and read ranges and frees their memory.
        void Clear() {
            m_variant_handler.Clear();
            m_sequence_range_handler.Clear();
            std::vector<SNP>{}.swap(m_snp_tmp);
        }

        const SequenceRangeHandler& GetSequenceRangeHandler() const {
            return m_sequence_range_handler;
        }

        SequenceRangeHandler& GetSequenceRangeHandler() {
            return m_sequence_range_handler;
        }

        void PostProcess(size_t min_observations=2, size_t min_observations_fwdrev=2, double min_frequency=0.2, size_t min_avg_quality=15, size_t min_phred_sum=0, bool require_strand=false) {
            auto forward = m_sequence_range_handler.CalculateCoverageVector2(SequenceRange::kForward);
            auto reverse = m_sequence_range_handler.CalculateCoverageVector2(SequenceRange::kReverse);
            m_variant_handler.PostProcessSNPs(forward, reverse, min_observations, min_observations_fwdrev, min_frequency, min_avg_quality, min_phred_sum, require_strand);
        }

        // Reads per position that have a base there: the coverage less the reads with an N or inside
        // a deletion. This is the depth SNP calls and the MSA are judged by.
        CoverageVec InformativeCoverage() const {
            auto cov = m_sequence_range_handler.CalculateCoverageVector2();
            for (auto const& [pos, reads] : m_variant_handler.GetNoBase()) {
                if (pos >= cov.size()) continue;
                uint32_t const none = static_cast<uint32_t>(reads.size());
                cov[pos] = cov[pos] > none ? cov[pos] - none : 0;
            }
            return cov;
        }

        // The gene as the strain MSA takes it: its variant bins (sorted by position) and informative
        // coverage from the reads with at least `min_identity` only, the taxon's own reads (see
        // Taxon::OwnIdentityThreshold): reads of relatives would add false SNPs and mixtures. The
        // reference allele of each bin is inferred again from these reads. The profile's variants
        // (and the model's features) keep every read.
        std::pair<VariantVec, CoverageVec> MSAItem(double min_identity, size_t min_observations, double min_frequency,
                                                   size_t min_avg_quality, size_t min_phred_sum, bool require_strand) const {
            uint8_t const max_divergence = MaxDivergenceBin(min_identity);
            auto forward = m_sequence_range_handler.CalculateCoverageVector2(SequenceRange::kForward, max_divergence);
            auto reverse = m_sequence_range_handler.CalculateCoverageVector2(SequenceRange::kReverse, max_divergence);
            for (size_t pos = 0; pos < forward.size(); pos++) {
                uint32_t const f = m_variant_handler.NoBase(pos, SequenceRange::kForward, max_divergence);
                uint32_t const r = m_variant_handler.NoBase(pos, SequenceRange::kReverse, max_divergence);
                forward[pos] = forward[pos] > f ? forward[pos] - f : 0;
                if (pos < reverse.size()) reverse[pos] = reverse[pos] > r ? reverse[pos] - r : 0;
            }

            VariantVec bins;
            for (auto const& [pos, bin] : m_variant_handler.GetVariants()) {
                VariantBin own;
                for (auto const& v : bin) {
                    if (v.IsReference()) continue;
                    auto filtered = v.WithMaxDivergence(max_divergence);
                    if (filtered.Observations() > 0) own.push_back(std::move(filtered));
                }
                if (own.empty()) continue;
                size_t const f = pos < forward.size() ? forward[pos] : 0;
                size_t const r = pos < reverse.size() ? reverse[pos] : 0;
                VariantHandler::PostProcessSNPBin(own, f, r, min_observations, min_observations, min_frequency,
                                                  min_avg_quality, min_phred_sum, require_strand);
                bins.push_back(std::move(own));
            }
            std::sort(bins.begin(), bins.end(), [](VariantBin const& a, VariantBin const& b) {
                return a.front().Position() < b.front().Position();
            });

            CoverageVec informative(forward.size(), 0);
            for (size_t pos = 0; pos < informative.size(); pos++) {
                informative[pos] = forward[pos] + (pos < reverse.size() ? reverse[pos] : 0);
            }
            return { std::move(bins), std::move(informative) };
        }

        void AddToSequenceRange(SamEntry const& sam, size_t read_id) {
            size_t start = sam.m_pos-1;
            size_t length = AlignmentLengthRef(sam.m_cigar);
            size_t end = start + length;

            SequenceRange query_range(start, end);
            ReadInfo rinfo;
            rinfo.read_id = read_id;
            rinfo.length = length;
            rinfo.start = start;
            rinfo.forward = !Flag::IsReverseComplement(sam.m_flag);
            query_range.AddReadInfo(rinfo);
            m_sequence_range_handler.Merge(std::move(query_range));
        }

        // Adds an alignment of a read (read_id: its fragment, the same for both mates). With
        // read_variants, its variants are recorded and its range is the part VariantHandler::AddAlignment
        // trusts, less what the fragment's other mate already covered; an alignment that does not fit
        // the gene adds nothing and returns false. Without, only its whole range is added. With read_variants
        // and `alleles`, what the read shows of the gene is also written there (ReadAlleles, for phasing).
        bool AddSam(SamEntry const& sam, size_t read_id, bool read_variants=false, double identity=1.0,
                    ReadAlleles* alleles = nullptr) {
            uint8_t const divergence = DivergenceBin(identity);
            if (!read_variants) {
                AddToSequenceRange(sam, read_id);
                return true;
            }
            size_t skip_begin = 0, skip_end = 0;
            if (read_id == m_last_read_id) std::tie(skip_begin, skip_end) = m_last_interval;
            auto const interval = m_variant_handler.AddAlignment(sam, skip_begin, skip_end, divergence, alleles);
            if (!interval) {
                return false;
            }
            auto const [start, end] = *interval;
            bool const forward = !Flag::IsReverseComplement(sam.m_flag);
            AddReadRange(start, std::min(end, skip_begin), read_id, forward, divergence);
            AddReadRange(std::max(start, skip_end), end, read_id, forward, divergence);

            if (read_id == m_last_read_id) {
                m_last_interval = { std::min(start, m_last_interval.first), std::max(end, m_last_interval.second) };
            } else {
                m_last_read_id = read_id;
                m_last_interval = { start, end };
            }
            return true;
        }
    };

    class SharedAlignmentRegion {
    public:
        const SequenceRangeHandler share_range;
        VariantVec variant_handler_a;
        VariantVec variant_handler_b;
//        const VariantVec variant_handler_a;
//        const VariantVec variant_handler_b;

        SharedAlignmentRegion(SequenceRangeHandler const&& shared_range,
                              VariantVec const&& variant_handler_a,
                              VariantVec const&& variant_handler_b) :
                share_range(shared_range),
                variant_handler_a(variant_handler_a),
                variant_handler_b(variant_handler_b) {};

        SharedAlignmentRegion(SequenceRangeHandler const& shared_range,
                              VariantVec const& variant_handler_a,
                              VariantVec const& variant_handler_b) :
                share_range(shared_range),
                variant_handler_a(variant_handler_a),
                variant_handler_b(variant_handler_b) {};

        static VariantVec GetSNPs(VariantHandler const& a) {
            VariantVec var;

            for (auto& [pos, varbin] : a.GetVariants()) {
                var.emplace_back(varbin);
            }

            std::sort(var.begin(), var.end(), [](VariantBin const& var1, VariantBin const& var2) {
                return var1.front().Position() < var2.front().Position();
            });

            return var;
        }

        static VariantVec GetSNPsInRegion(CoverageVec const& coverage, VariantHandler const& a) {
            VariantVec var;

//            std::cout << "GREP GetSNPsInRegion: " << a.GetVariants().size() << std::endl;
            for (auto& [pos, varbin] : a.GetVariants()) {
//                std::cout << "Pos: " << pos << " cov: " << coverage[pos] << std::endl;
                if (coverage[pos] != 0) {
                    var.emplace_back(varbin);
                    varbin.front().Position();
                }
            }
            std::sort(var.begin(), var.end(), [](VariantBin const& var1, VariantBin const& var2) {
                return var1.front().Position() < var2.front().Position();
            });

            return var;
        }

        static SharedAlignmentRegion GetSharedAlignmentRegion(StrainLevelContainer& a, StrainLevelContainer& b) {
            size_t min_cov = 3;
            size_t min_len = 30;
            SequenceRangeHandler ranges_a, ranges_b;
//            a.GetSequenceRangeHandler().CalculateCoverageVector();
//            b.GetSequenceRangeHandler().CalculateCoverageVector();

//            std::cout << "GetSharedAlignmentRegion: " << a.GetVariantHandler().GetVariants().size() << " " << b.GetVariantHandler().GetVariants().size() << std::endl;
            auto intersection = a.GetSequenceRangeHandler().Intersect(b.GetSequenceRangeHandler(), min_cov, 10);

            auto variant_vec_a = GetSNPsInRegion(intersection.GetCoverageVector(), a.GetVariantHandler());
            auto variant_vec_b = GetSNPsInRegion(intersection.GetCoverageVector(), b.GetVariantHandler());

            return SharedAlignmentRegion(intersection, variant_vec_a, variant_vec_b);
        };
    };



    using DoubleMatrix = Matrix<double>;


    static bool VariantPass(Variant const& call, VariantBin const& bin,
                            size_t min_var_qual_sum=50, size_t min_var_cov=3,
                            double min_frequency=0.0, uint32_t coverage=0,
                            bool require_strand=false, size_t min_mean_qual=0) {
        if (call.Observations() < min_var_cov) return false;
        // OR logic: passes if either quality gate holds
        bool quality_ok = call.QualitySum() >= min_var_qual_sum ||
                          (min_mean_qual > 0 && call.MeanQuality() >= min_mean_qual);
        if (!quality_ok) return false;
        if (min_frequency > 0.0 && coverage > 0) {
            double freq = static_cast<double>(call.Observations()) / coverage;
            if (freq < min_frequency) return false;
        }
        if (require_strand && !PassesStrandFilter(call, bin)) return false;
        return true;
    }

    // The bases of a site's SNP alleles (the reference included) that pass, at most max_alleles, the
    // most observed first (RanksBefore): more than one is written as an IUPAC code.
    template<typename Passes>
    static std::vector<char> PassingBases(VariantBin const& bin, Passes&& passes, size_t max_alleles) {
        std::vector<Variant const*> alleles;
        for (auto const& v : bin) if (v.IsSNP() && passes(v)) alleles.push_back(&v);
        std::stable_sort(alleles.begin(), alleles.end(), [](Variant const* a, Variant const* b) { return RanksBefore(*a, *b); });
        std::vector<char> bases;
        for (auto const* a : alleles) {
            if (bases.size() >= max_alleles) break;
            if (std::find(bases.begin(), bases.end(), a->GetVariant()) == bases.end()) bases.push_back(a->GetVariant());
        }
        return bases;
    }

    // Positions that MSA() writes as an IUPAC code, with the same rule and parameters, each with its passing bases (at
    // most max_alleles, the most observed first): the site's call (BaseCall) is a single-base allele that passes, and
    // at least one other base passes too; sites inside a deletion the MSA writes as gaps do not count. `bins` and
    // `coverage` are the gene as the MSA takes it (StrainLevelContainer::MSAItem: bins sorted by position, informative
    // coverage). A long-read sample's strains are phased at these sites (Haplotypes.h).
    static std::vector<std::pair<uint32_t, std::vector<char>>> MultiAllelicSites(
            std::vector<VariantBin> const& bins, CoverageVec const& coverage, uint32_t min_cov, uint32_t min_qual_sum,
            double min_frequency, bool require_strand, size_t min_mean_qual, size_t max_alleles) {
        std::vector<std::pair<uint32_t, std::vector<char>>> sites;
        if (max_alleles < 2) return sites;
        size_t deleted_until = 0;
        for (auto const& bin : bins) {
            if (bin.empty()) continue;
            size_t const pos = bin.front().Position();
            if (pos < deleted_until) continue;
            uint32_t const cov = pos < coverage.size() ? coverage[pos] : 0;
            if (cov < min_cov) continue;
            auto passes = [&](Variant const& v) {
                return VariantPass(v, bin, min_qual_sum, min_cov, min_frequency, cov, require_strand, min_mean_qual);
            };
            bool pass = false;
            auto const* call = BaseCall(bin, passes, pass);
            if (!call || !pass) continue;
            if (call->IsDEL()) {
                deleted_until = pos + call->GetStructuralSize();
                continue;
            }
            auto bases = PassingBases(bin, passes, max_alleles);
            if (bases.size() > 1) sites.emplace_back(static_cast<uint32_t>(pos), std::move(bases));
        }
        return sites;
    }

    // How many positions MSA() writes as an IUPAC code (MultiAllelicSites with snp_max_alleles). The .meta.tsv
    // reports this count, so that qcmsa filters on what the MSA holds.
    static size_t MultiAllelicPositions(std::vector<VariantBin> const& bins, CoverageVec const& coverage, uint32_t min_cov,
                                        uint32_t min_qual_sum, double min_frequency, bool require_strand,
                                        size_t min_mean_qual, size_t snp_max_alleles) {
        return MultiAllelicSites(bins, coverage, min_cov, min_qual_sum, min_frequency, require_strand, min_mean_qual,
                                 snp_max_alleles).size();
    }

    using MSAVector = std::vector<std::vector<char>>;
    using MSARow = std::vector<char>;
    static void AddInsertionGap(MSARow& msa_row, size_t ins_count) {
        for (auto j = ins_count; j > 0; j--) msa_row.emplace_back('-');
    }

    // Returns IUPAC ambiguity code for a sorted, deduplicated set of allele bases.
    static char IUPACCode(std::vector<char>& alleles) {
        static const std::unordered_map<std::string, char> iupac = {
            {"AG", 'R'}, {"CT", 'Y'}, {"AT", 'W'}, {"CG", 'S'}, {"AC", 'M'}, {"GT", 'K'},
            {"CGT", 'B'}, {"ACT", 'H'}, {"AGT", 'D'}, {"ACG", 'V'}, {"ACGT", 'N'}
        };
        std::string key(alleles.begin(), alleles.end());
        auto it = iupac.find(key);
        return (it != iupac.end()) ? it->second : (alleles.empty() ? 'N' : alleles.front());
    }

    using VariantVecRef = std::reference_wrapper<VariantVec>;
    // A sample's gene for the MSA: its variant bins (sorted by position) and its informative coverage
    // (StrainLevelContainer::InformativeCoverage).
    using OptionalMSASequenceItem = std::optional<std::pair<VariantVec, CoverageVec>>;
    using MSASequenceItems = std::vector<OptionalMSASequenceItem>;
    using CoverageVecs = std::vector<CoverageVec>;
    using OptionalVariant = std::optional<Variant>;
    using OptionalVariantBin = std::optional<VariantBin>;

    // Per-sample statistics on SNP/variant retention during MSA construction.
    // One entry per sample (aligned with MSASequenceItems indices).
    struct MSASampleStats {
        // Variant outcomes at positions where a variant was called
        size_t snps_retained = 0;              ///< Variant passed all filters and is a SNP
        size_t refs_retained = 0;              ///< Consensus is the reference allele and passed all filters
        size_t insertions_retained = 0;        ///< Variant passed all filters and is an insertion
        size_t deletions_retained = 0;         ///< Variant passed all filters and is a deletion
        size_t variants_filtered_qual_sum = 0; ///< Variant failed quality gate (phred_sum < min AND mean_qual < min)
        size_t variants_filtered_obs_cov = 0;  ///< Variant failed minimum observation count filter (--snp_min_cov)
        size_t variants_filtered_af = 0;       ///< Variant failed allele frequency filter (--snp_min_af)
        size_t variants_filtered_strand = 0;   ///< Variant failed strand-bias filter (--snp_no_strand not set)

        // Position-level outcomes (one count per reference position, excluding deletion continuations)
        size_t positions_ref = 0;              ///< Covered >=min_cov, no variant called → reference base used
        size_t positions_below_min_cov = 0;    ///< Coverage present but below min_depth (--msa_min_depth) → gap
        size_t positions_no_coverage = 0;      ///< No coverage or sample absent → gap

        // Filled externally after ProcessMSA (vertical coverage filter)

        size_t TotalPass() const {
            return snps_retained + refs_retained + insertions_retained + deletions_retained;
        }

        size_t TotalFiltered() const {
            return variants_filtered_qual_sum + variants_filtered_obs_cov
                 + variants_filtered_af + variants_filtered_strand;
        }

        size_t TotalVariantPositions() const {
            return TotalPass() + TotalFiltered();
        }

        // All reference positions visited for this sample (covered + uncovered).
        size_t TotalPositions() const {
            return positions_ref + TotalVariantPositions()
                   + positions_below_min_cov + positions_no_coverage;
        }

        // Percentage helpers — return 0 when denominator is zero.
        double PctPass() const {
            auto d = TotalVariantPositions();
            return d > 0 ? 100.0 * TotalPass() / d : 0.0;
        }

        double PctFiltered() const {
            auto d = TotalVariantPositions();
            return d > 0 ? 100.0 * TotalFiltered() / d : 0.0;
        }

        double PctFilteredQualSum() const {
            auto d = TotalVariantPositions();
            return d > 0 ? 100.0 * variants_filtered_qual_sum / d : 0.0;
        }

        double PctFilteredObsCov() const {
            auto d = TotalVariantPositions();
            return d > 0 ? 100.0 * variants_filtered_obs_cov / d : 0.0;
        }

        double PctFilteredAF() const {
            auto d = TotalVariantPositions();
            return d > 0 ? 100.0 * variants_filtered_af / d : 0.0;
        }

        double PctFilteredStrand() const {
            auto d = TotalVariantPositions();
            return d > 0 ? 100.0 * variants_filtered_strand / d : 0.0;
        }

        double PctSnpsRetained() const {
            auto d = TotalVariantPositions();
            return d > 0 ? 100.0 * snps_retained / d : 0.0;
        }

        double PctRefsRetained() const {
            auto d = TotalVariantPositions();
            return d > 0 ? 100.0 * refs_retained / d : 0.0;
        }

        double PctInsertionsRetained() const {
            auto d = TotalVariantPositions();
            return d > 0 ? 100.0 * insertions_retained / d : 0.0;
        }

        double PctDeletionsRetained() const {
            auto d = TotalVariantPositions();
            return d > 0 ? 100.0 * deletions_retained / d : 0.0;
        }

        double PctPositionsRef() const {
            auto d = TotalPositions();
            return d > 0 ? 100.0 * positions_ref / d : 0.0;
        }

        double PctPositionsBelowMinCov() const {
            auto d = TotalPositions();
            return d > 0 ? 100.0 * positions_below_min_cov / d : 0.0;
        }

        double PctPositionsNoCoverage() const {
            auto d = TotalPositions();
            return d > 0 ? 100.0 * positions_no_coverage / d : 0.0;
        }

        MSASampleStats& operator+=(MSASampleStats const& o) {
            snps_retained              += o.snps_retained;
            refs_retained              += o.refs_retained;
            insertions_retained        += o.insertions_retained;
            deletions_retained         += o.deletions_retained;
            variants_filtered_qual_sum += o.variants_filtered_qual_sum;
            variants_filtered_obs_cov  += o.variants_filtered_obs_cov;
            variants_filtered_af       += o.variants_filtered_af;
            variants_filtered_strand   += o.variants_filtered_strand;
            positions_ref              += o.positions_ref;
            positions_below_min_cov    += o.positions_below_min_cov;
            positions_no_coverage      += o.positions_no_coverage;
            return *this;
        }
    };
    using MSAStats = std::vector<MSASampleStats>;

    static auto FindUnequalIndex(std::vector<OptionalVariant> const& column, std::vector<bool> const& column_pass) {
        auto count = std::count_if(column.begin(), column.end(), [](OptionalVariant const& v) { return v.has_value(); });
        if (count <= 1) return -1;

        char baseline = 'X';
        for (auto i = 0; i < column.size(); i++) {
            if (column_pass[i] && column[i].has_value() && column[i].value().GetVariant() != 'N') {
                if (baseline == 'X') {
                    baseline = column[i].value().GetVariant();
                } else {
                    if (column[i].value().GetVariant() != baseline)
                        return i;
                }
            }
        }

//        auto baseline = std::find_if(column.begin(), column.end(), [](OptionalVariant const& v) { return v.has_value() && v.value().GetVariant() != 'N'; });
//        for (auto& v : column) {
//            if (v.has_value() && v.value().GetVariant() != 'N' && !v.value().Match(baseline->value())) {
//                return true;
//            }
//        }
        return -1;
    }

    static auto WrongPos(std::vector<OptionalVariant> const& column) {
        auto count = std::count_if(column.begin(), column.end(), [](OptionalVariant const& v) { return v.has_value(); });
        if (count <= 1) return false;
        auto baseline = std::find_if(column.begin(), column.end(), [](OptionalVariant const& v) { return v.has_value(); });
        auto pos = baseline->value().Position();
        auto correct = std::all_of(column.begin(), column.end(), [pos](OptionalVariant const& var) {
            return !var.has_value() || var.value().Position() == pos;
        });
        std::cout << pos << ", " << count << " -> " << correct << std::endl;
        return !correct;
    }

    static auto IsVariantBinBad(VariantBin const& bin) {
        auto baseline = bin.front().Position();
        for (auto& v : bin) {
            if (baseline != v.Position()) return true;
        }
        return false;
    }

    static auto PrintColumn(std::vector<OptionalVariant> const& column) {
        std::cout << "Bad column" << std::endl;
        for (auto& var : column) {
            if (var.has_value())
                std::cout << var->ToString() << std::endl;
            else std::cout << " no " << std::endl;
        }
    }



    // The items' informative coverage vectors (empty for a sample without the gene).
    static CoverageVecs LoadCoverageVectors(MSASequenceItems const& items, std::string_view const reference) {
        CoverageVecs covs;
        for (auto i = 0; i < items.size(); i++) {
            covs.emplace_back(items[i].has_value() ? items[i]->second : CoverageVec());
            if (covs.back().size() > reference.size()) {
                std::cerr << "Coverage of MSA item " << i << " is longer than its gene: " << covs.back().size()
                          << " > " << reference.size() << std::endl;
            }
        }
        return covs;
    }

    static size_t GetValidBases(CoverageVecs const& covs, uint32_t min_cov) {
        return std::count_if(covs.begin(), covs.end(), [min_cov](auto const& cv) {
            return std::count_if(cv.begin(), cv.end(), [min_cov](auto c){ return c >= min_cov; });
        });
    }

    static std::vector<bool> HasVariantVector(MSASequenceItems const& items, std::string_view const reference) {
        std::vector<bool> has_variant(reference.length(), false);
        for (auto& item : items) {
            if (!item.has_value()) continue;
            for (const auto& var : item.value().first) {
                if (var.front().Position() < has_variant.size()) has_variant[var.front().Position()] = true;
            }
        }
        return has_variant;
    }

    // The MSA of a gene over samples (items), each with its own minimum allele frequency
    // (min_frequencies, one per item: the SNP filter of its reads' kind). Per sample and position:
    //   - fewer than min_depth reads with a base (informative coverage): '-';
    //   - no variant: the reference base;
    //   - else the site's call (BaseCall): its base, an IUPAC code where more than one SNP allele
    //     passes (at most snp_max_alleles), '-' over a passing deletion, or N if no allele passes;
    //   - a passing insertion before the base fills insertion columns, which all other rows (and the
    //     reference row) get as '-'.
    // An allele passes with min_cov reads, or with fewer when it has every read of the site: reads
    // that agree need no second one, while a mixture needs min_cov reads per allele.
    static bool MSA(MSASequenceItems const& items, std::string_view const reference, MSAVector& msa, uint32_t min_cov, uint32_t min_qual_sum, std::vector<double> const& min_frequencies, bool require_strand=false, size_t min_mean_qual=0, MSAStats* stats = nullptr, MSARow* ref_row = nullptr, size_t snp_max_alleles = 1, uint32_t min_depth = 1) {
        if (min_frequencies.size() != items.size()) {
            std::cerr << "MSA: " << min_frequencies.size() << " minimum allele frequencies for " << items.size() << " samples" << std::endl;
            return false;
        }
        if (msa.size() != items.size()) {
            std::cerr << msa.size() << " != " << items.size() << " <- items" << std::endl;
            std::cerr << "Msa object must be of the same length as items" << std::endl;
            return false;
        }

        if (items.empty()) return false;

        min_depth = std::max<uint32_t>(min_depth, 1);
        CoverageVecs covs = LoadCoverageVectors(items, reference);
        if (GetValidBases(covs, min_depth) == 0) return false;

        // Where this gene's columns begin, to take them back should the rows get out of step.
        std::vector<size_t> row_starts;
        for (auto const& row : msa) row_starts.push_back(row.size());
        size_t const ref_row_start = ref_row ? ref_row->size() : 0;

        std::vector<size_t> indices(items.size(), 0);
        std::vector<uint16_t> pause_timer(items.size(), 0);
        std::vector<VariantBin const*> bins(items.size(), nullptr);
        std::vector<Variant const*> base_calls(items.size(), nullptr);
        std::vector<Variant const*> ins_calls(items.size(), nullptr);
        std::vector<bool> base_pass(items.size(), false);

        const char LACKING_COVERAGE = '-';
        const char NO_PASS = 'N';

        for (size_t rpos = 0; rpos < reference.length(); rpos++) {
            char const ref = reference[rpos];
            size_t max_ins = 0;

            // The calls of every sample at this position.
            for (size_t i = 0; i < items.size(); i++) {
                bins[i] = nullptr;
                base_calls[i] = nullptr;
                ins_calls[i] = nullptr;
                base_pass[i] = false;
                auto const& cov = covs[i];
                uint32_t const depth = rpos < cov.size() ? cov[rpos] : 0;
                if (!items[i].has_value() || depth == 0) continue;

                auto const& variants = items[i].value().first;
                while (indices[i] < variants.size() && variants[indices[i]].front().Position() < rpos) indices[i]++;
                if (indices[i] == variants.size() || variants[indices[i]].front().Position() != rpos) continue;

                auto const& bin = variants[indices[i]];
                auto passes = [&](Variant const& v) {
                    uint32_t const reads = v.Observations() >= depth ? 1 : min_cov;
                    return VariantPass(v, bin, min_qual_sum, reads, min_frequencies[i], depth, require_strand, min_mean_qual);
                };
                bins[i] = &bin;
                bool pass = false;
                base_calls[i] = BaseCall(bin, passes, pass);
                base_pass[i] = pass;
                ins_calls[i] = InsertionCall(bin, passes);
                if (ins_calls[i] && pause_timer[i] == 0) max_ins = std::max<size_t>(max_ins, ins_calls[i]->GetStructuralSize());
            }

            // Reference row: gaps for any insertion slots, then the reference base itself.
            if (ref_row) {
                AddInsertionGap(*ref_row, max_ins);
                ref_row->emplace_back(ref);
            }

            bool bad = false;
            for (size_t i = 0; i < items.size(); i++) {
                auto& msa_row = msa[i];
                size_t const before = msa_row.size();
                auto const& cov = covs[i];
                uint32_t const depth = rpos < cov.size() ? cov[rpos] : 0;

                if (pause_timer[i] > 0) {
                    // Inside a deletion.
                    AddInsertionGap(msa_row, max_ins);
                    msa_row.emplace_back('-');
                    pause_timer[i]--;
                } else if (!items[i].has_value() || depth < min_depth) {
                    if (stats) {
                        auto& s = (*stats)[i];
                        if (!items[i].has_value() || depth == 0) s.positions_no_coverage++;
                        else s.positions_below_min_cov++;
                    }
                    AddInsertionGap(msa_row, max_ins);
                    msa_row.emplace_back(LACKING_COVERAGE);
                } else {
                    // Insertion columns: the sample's passing insertion, left-aligned, or gaps.
                    if (auto const* ins = ins_calls[i]) {
                        if (stats) (*stats)[i].insertions_retained++;
                        for (auto c : ins->GetStructural()) msa_row.emplace_back(c);
                        AddInsertionGap(msa_row, max_ins - ins->GetStructuralSize());
                    } else {
                        AddInsertionGap(msa_row, max_ins);
                    }

                    auto const* call = base_calls[i];
                    if (!bins[i] || !call) {
                        // No variant: every read with a base shows the reference.
                        if (stats) (*stats)[i].positions_ref++;
                        msa_row.emplace_back(ref);
                    } else if (!base_pass[i]) {
                        // No allele passes: an ambiguous base. Each filter the top allele fails is counted.
                        if (stats) {
                            auto& s = (*stats)[i];
                            if (call->Observations() < min_cov) s.variants_filtered_obs_cov++;
                            bool const qual_fails = call->QualitySum() < min_qual_sum &&
                                                    !(min_mean_qual > 0 && call->MeanQuality() >= min_mean_qual);
                            if (qual_fails) s.variants_filtered_qual_sum++;
                            if (min_frequencies[i] > 0.0 && static_cast<double>(call->Observations()) / depth < min_frequencies[i])
                                s.variants_filtered_af++;
                            if (require_strand && !PassesStrandFilter(*call, *bins[i])) s.variants_filtered_strand++;
                        }
                        msa_row.emplace_back(NO_PASS);
                    } else if (call->IsDEL()) {
                        if (stats) (*stats)[i].deletions_retained++;
                        pause_timer[i] = call->GetStructuralSize() - 1;
                        msa_row.emplace_back('-');
                    } else {
                        if (stats) {
                            if (call->IsReference()) (*stats)[i].refs_retained++;
                            else (*stats)[i].snps_retained++;
                        }
                        char base = call->GetVariant();
                        if (snp_max_alleles > 1) {
                            auto const& bin = *bins[i];
                            auto passes = [&](Variant const& v) {
                                return VariantPass(v, bin, min_qual_sum, min_cov, min_frequencies[i], depth, require_strand, min_mean_qual);
                            };
                            auto alleles = PassingBases(bin, passes, snp_max_alleles);
                            std::sort(alleles.begin(), alleles.end());
                            if (alleles.size() > 1) base = IUPACCode(alleles);
                        }
                        msa_row.emplace_back(base);
                    }
                }

                if (msa_row.size() != before + max_ins + 1 || (i > 0 && msa[i].size() != msa[i - 1].size())) bad = true;
            }

            if (bad) {
                std::cerr << "Error: the MSA rows differ in length at reference position " << rpos
                          << "; the gene is left out of the MSA" << std::endl;
                for (size_t i = 0; i < msa.size(); i++) msa[i].resize(row_starts[i]);
                if (ref_row) ref_row->resize(ref_row_start);
                return false;
            }
        }

        return true;
    }

    // The MSA with one minimum allele frequency for all samples.
    static bool MSA(MSASequenceItems const& items, std::string_view const reference, MSAVector& msa, uint32_t min_cov, uint32_t min_qual_sum, double min_frequency=0.0, bool require_strand=false, size_t min_mean_qual=0, MSAStats* stats = nullptr, MSARow* ref_row = nullptr, size_t snp_max_alleles = 1, uint32_t min_depth = 1) {
        return MSA(items, reference, msa, min_cov, min_qual_sum, std::vector<double>(items.size(), min_frequency), require_strand,
                   min_mean_qual, stats, ref_row, snp_max_alleles, min_depth);
    }
}
