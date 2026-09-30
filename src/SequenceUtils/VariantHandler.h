//
// Created by fritsche on 21/02/23.
//

#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <optional>
#include <utility>
#include "SNPUtils.h"
#include "Variant.h"
#include "robin_map.h"
#include "GenomeLoader.h"
#include "SequenceRange.h"

namespace protal {
    // The reads at a site, per strand: those with a base there (the non-insertion alleles of its bin).
    inline std::pair<uint32_t, uint32_t> SiteStrandCounts(VariantBin const& bin) {
        uint32_t forward = 0, reverse = 0;
        for (auto const& v : bin) {
            if (v.IsINS()) continue;
            forward += v.ObservationsForward();
            reverse += v.ObservationsReverse();
        }
        return { forward, reverse };
    }

    // Strand-bias filter. An allele seen on one strand only fails when that is unlikely among the
    // site's reads: the chance that its k reads, drawn from the site's n reads, all come from the s
    // reads of that strand is prod_{i<k} (s-i)/(n-i). Below kStrandBiasP the allele fails. Where all
    // reads are on one strand (common at low depth), an allele on that strand passes. The reference
    // allele is tested like any other.
    inline constexpr double kStrandBiasP = 0.05;
    inline bool PassesStrandFilter(Variant const& allele, VariantBin const& bin) {
        if (allele.HasFwdAndRev()) return true;
        size_t const k = allele.Observations();
        if (k == 0) return true;
        auto [site_forward, site_reverse] = SiteStrandCounts(bin);
        size_t s = allele.ObservationsForward() > 0 ? site_forward : site_reverse;
        size_t n = site_forward + site_reverse;
        s = std::max(s, k);  // an insertion's reads are counted by their base, which may lie elsewhere
        n = std::max(n, s);
        double p = 1.0;
        for (size_t i = 0; i < k && p >= kStrandBiasP; i++) {
            p *= static_cast<double>(s - i) / static_cast<double>(n - i);
        }
        return p >= kStrandBiasP;
    }

    class VariantHandler {
        Variants m_variants;
        // Reads that cover a position without a base there: an N in the read or the reference, or
        // inside a deletion after its first position. They support no allele, so they are not taken
        // for reference support. One byte per read: its strand (bit 7) and divergence (bits 0-6).
        tsl::robin_map<VariantPos, std::vector<uint8_t>> m_no_base;
        std::string_view m_reference;  // the gene's sequence, which outlives the handler

    public:
        Benchmark bm_next_compressed_cigar{"Next compressed cigar"};

        explicit VariantHandler(std::string_view reference) : m_reference(reference) {};

        // Drops all variants and frees their memory.
        void Clear() {
            Variants{}.swap(m_variants);
            tsl::robin_map<VariantPos, std::vector<uint8_t>>{}.swap(m_no_base);
        }

        // Phred score of a Sanger/Illumina 1.8+ (Phred+33) quality character. Characters below the
        // offset clamp to 0 instead of wrapping around in the unsigned Qual type.
        static Qual PhredScore(char quality_char) {
            constexpr int phred_offset = 33;
            return quality_char > phred_offset ? static_cast<Qual>(quality_char - phred_offset) : 0;
        }

        bool HasVariantBin(VariantPos position) {
            return m_variants.contains(position);
        }

        VariantBin& GetVariantBin(VariantPos position) {
            return m_variants.at(position);
        }

        bool HasReference(VariantBin& variant_bin, VariantPos pos) {
            return std::any_of(variant_bin.begin(), variant_bin.end(), [](Variant const& var) {
                return var.IsReference();
            });
        }

        static Variant& GetVariant(VariantBin& variant_bin, VariantPos pos, Base snp, Base ref) {
            for (auto& variant :  variant_bin) {
                if (variant.Match(pos, snp, ref)) {
                    return variant;
                }
            }
            variant_bin.emplace_back( Variant(pos, snp, ref) );
            return variant_bin.back();
        }


        std::string_view GetReference() const {
            return m_reference;
        }

        static Variant& GetVariant(VariantBin& variant_bin, VariantType type, VariantPos pos, Base ref, std::string& structural) {
            for (auto& variant :  variant_bin) {
                if (variant.Match(type, pos, structural)) {
                    return variant;
                }
            }
            variant_bin.emplace_back( Variant(type, pos, ref, structural) );
            return variant_bin.back();
        }

        const Variants& GetVariants() const {
            return m_variants;
        }

        // Reads without a base at `pos` (see m_no_base): one strand, or both; with max_divergence,
        // only reads of at most that divergence.
        uint32_t NoBase(VariantPos pos, int strand = SequenceRange::kBothStrands, uint8_t max_divergence = 127) const {
            auto it = m_no_base.find(pos);
            if (it == m_no_base.end()) return 0;
            uint32_t n = 0;
            for (auto read : it->second) {
                bool const forward = read & 0x80;
                if ((read & 0x7f) > max_divergence) continue;
                if (strand == SequenceRange::kBothStrands || forward == (strand == SequenceRange::kForward)) n++;
            }
            return n;
        }

        const tsl::robin_map<VariantPos, std::vector<uint8_t>>& GetNoBase() const {
            return m_no_base;
        }

        void AddSNP(VariantPos position, Base snp, Base ref, bool on_forward, Qual quality, uint8_t divergence = 0) {
            if (!HasVariantBin(position)) m_variants.insert({position, VariantBin() });
            auto& variant_bin = GetVariantBin(position);
            auto& variant = GetVariant(variant_bin, position, snp, ref);
            variant.AddObservation(quality, on_forward, divergence);
        }

        void AddINDEL(VariantType type, VariantPos position, Base ref, std::string&& structural, bool on_forward, Qual quality, uint8_t divergence = 0) {
            auto& variant_bin = HasVariantBin(position) ? GetVariantBin(position) : m_variants[position];
            auto& variant = GetVariant(variant_bin, type, position, ref, structural);
            variant.AddObservation(quality, on_forward, divergence);
        }

        void AddNoBase(VariantPos position, bool on_forward, uint8_t divergence = 0) {
            m_no_base[position].push_back(static_cast<uint8_t>((on_forward ? 0x80 : 0) | (divergence & 0x7f)));
        }

        // A terminal segment of an alignment, before its first (after its last) run of kAnchor
        // matching bases, is not trusted when it holds this many mismatched bases or more: reads
        // that start or end at an insertion are written with such clusters of mismatches.
        static constexpr int kAnchor = 5;
        static constexpr int kUntrustedMismatches = 2;

        // Records the variants one alignment supports, and returns the reference interval
        // [start, end) its trusted part covers, or nothing if the alignment does not fit the gene
        // (a position past its end, or M at a mismatch: protal writes M only for exact matches).
        // Nothing is recorded for an alignment that does not fit.
        //
        // SNPs come from X ops. A base that is N in the read or the reference supports no allele,
        // nor do the positions of a deletion after its first one (the deletion itself is the allele
        // at its first position). Insertions and deletions are called only between trusted aligned
        // bases: at either end of an alignment they are artefacts. A deletion has no base of its own
        // and gets the lower quality of its two flanking bases. Positions in [skip_begin, skip_end),
        // the part the fragment's other mate already covered, are left out: a fragment counts once.
        // Every record carries the read's divergence from the gene (DivergenceBin).
        std::optional<std::pair<size_t, size_t>> AddAlignment(SamEntry const& sam, size_t skip_begin = 0, size_t skip_end = 0, uint8_t divergence = 0) {
            std::vector<std::pair<int, char>> ops;
            {
                int cpos = 0, count = 0;
                char op = ' ';
                while (NextCompressedCigar(cpos, sam.m_cigar, count, op)) ops.emplace_back(count, op);
            }
            auto const n_ops = static_cast<std::ptrdiff_t>(ops.size());
            auto is_aligned = [](std::pair<int, char> const& o) { return o.second == 'M' || o.second == 'X'; };
            std::ptrdiff_t first = std::find_if(ops.begin(), ops.end(), is_aligned) - ops.begin();
            std::ptrdiff_t last = n_ops - 1 - (std::find_if(ops.rbegin(), ops.rend(), is_aligned) - ops.rbegin());
            if (first >= n_ops) return std::nullopt;  // no aligned base

            // Trim untrusted ends: a terminal segment with kUntrustedMismatches or more mismatches
            // before the first (after the last) anchor.
            {
                std::ptrdiff_t anchor = first;
                while (anchor <= last && !(ops[anchor].second == 'M' && ops[anchor].first >= kAnchor)) anchor++;
                if (anchor <= last) {
                    int mismatches = 0;
                    for (auto k = first; k < anchor; k++) mismatches += ops[k].second == 'X' ? ops[k].first : 0;
                    if (mismatches >= kUntrustedMismatches) first = anchor;
                }
                anchor = last;
                while (anchor >= first && !(ops[anchor].second == 'M' && ops[anchor].first >= kAnchor)) anchor--;
                if (anchor >= first) {
                    int mismatches = 0;
                    for (auto k = anchor + 1; k <= last; k++) mismatches += ops[k].second == 'X' ? ops[k].first : 0;
                    if (mismatches >= kUntrustedMismatches) last = anchor;
                }
            }

            struct Snp { VariantPos pos; Base base; Base ref; Qual qual; };
            struct Indel { VariantType type; VariantPos pos; Base ref; std::string structural; Qual qual; };
            std::vector<Snp> snps;
            std::vector<Indel> indels;
            std::vector<VariantPos> no_base;
            auto skipped = [&](size_t pos) { return pos >= skip_begin && pos < skip_end; };

            bool const is_fwd = !sam.IsReversed();
            size_t qpos = 0;
            size_t rpos = sam.m_pos - 1;
            size_t start = 0, end = 0;

            bm_next_compressed_cigar.Start();
            for (std::ptrdiff_t k = 0; k < n_ops; k++) {
                auto const [count, op] = ops[k];
                bool const consumes_ref = !(op == 'I' || op == 'S');
                bool const consumes_query = op != 'D';
                if ((consumes_ref && rpos + count > m_reference.size()) || (consumes_query && qpos + count > sam.m_seq.size())) {
                    bm_next_compressed_cigar.Stop();
                    return std::nullopt;
                }
                bool const trusted = k >= first && k <= last;
                bool const between_trusted_bases = k > first && k < last;
                if (k == first) start = rpos;

                if (op == 'M') {
                    for (auto i = 0; i < count; i++) {
                        char const base = sam.m_seq[qpos + i];
                        char const ref = m_reference[rpos + i];
                        if (base == 'N' || ref == 'N') {
                            if (trusted && !skipped(rpos + i)) no_base.push_back(rpos + i);
                        } else if (base != ref) {
                            // The record goes to the sample's .err file; the first few are shown.
                            static std::atomic<size_t> shown{0};
                            constexpr size_t kShown = 3;
                            size_t const n = shown++;
                            if (n < kShown) {
                                #pragma omp critical(print)
                                {
                                    std::cerr << "Alignment does not match its gene (M at a mismatch; a SAM aligned against "
                                                 "another database?):\n" << sam.ToString() << std::endl;
                                    PrintAlignment(sam, m_reference, std::cerr);
                                    if (n + 1 == kShown) std::cerr << "Further such alignments are not shown." << std::endl;
                                }
                            }
                            bm_next_compressed_cigar.Stop();
                            return std::nullopt;
                        }
                    }
                } else if (op == 'X' && trusted) {
                    for (auto i = 0; i < count; i++) {
                        if (skipped(rpos + i)) continue;
                        char const base = sam.m_seq[qpos + i];
                        char const ref = m_reference[rpos + i];
                        if (base == 'N' || ref == 'N') {
                            no_base.push_back(rpos + i);
                        } else {
                            snps.push_back({ static_cast<VariantPos>(rpos + i), base, ref, PhredScore(sam.m_qual[qpos + i]) });
                        }
                    }
                } else if (op == 'I' && between_trusted_bases && !(rpos > skip_begin && rpos < skip_end)) {
                    size_t qual_sum = 0;
                    for (auto i = 0; i < count; i++) {
                        qual_sum += PhredScore(sam.m_qual[qpos + i]);
                    }
                    indels.push_back({ VariantType::INS, static_cast<VariantPos>(rpos), m_reference[rpos],
                                       sam.m_seq.substr(qpos, count), static_cast<Qual>(qual_sum / count) });
                } else if (op == 'D' && between_trusted_bases && !(rpos > skip_begin && rpos + count < skip_end)) {
                    Qual const flank = std::min(PhredScore(sam.m_qual[qpos - 1]), PhredScore(sam.m_qual[qpos]));
                    indels.push_back({ VariantType::DEL, static_cast<VariantPos>(rpos), m_reference[rpos],
                                       std::string(m_reference.substr(rpos, count)), flank });
                    for (auto i = 1; i < count; i++) {
                        if (!skipped(rpos + i)) no_base.push_back(rpos + i);
                    }
                }

                qpos += consumes_query * count;
                rpos += consumes_ref * count;
                if (k == last) end = rpos;
            }
            bm_next_compressed_cigar.Stop();

            for (auto const& s : snps) AddSNP(s.pos, s.base, s.ref, is_fwd, s.qual, divergence);
            for (auto& d : indels) AddINDEL(d.type, d.pos, d.ref, std::move(d.structural), is_fwd, d.qual, divergence);
            for (auto pos : no_base) AddNoBase(pos, is_fwd, divergence);
            return std::make_pair(start, end);
        }

        // AddAlignment without a mate: whether the alignment fits the gene.
        bool AddVariantsFromSam(SamEntry const& sam) {
            return AddAlignment(sam).has_value();
        }

        static std::string VariantBinToMinimalString(const VariantBin& variant_bin) {
            std::string str;
            for (auto& variant : variant_bin) {
                str += variant.ToMinimalString() + '\t';
            }
            return str;
        }

        static std::string VariantBinToString(const VariantBin& variant_bin) {
            std::string str;
            for (auto& variant : variant_bin) {
                str += variant.ToString() + '\t';
            }
            return str;
        }

        static bool FilterSNPs(Variant const& var, VariantBin const& bin, size_t coverage, size_t min_observations, size_t min_observations_fwdrev, double min_frequency, size_t min_avg_quality, size_t min_phred_sum=0, bool require_strand=false) {
            auto observations = var.Observations();
            auto frequency = coverage == 0 ? 0.0 : static_cast<double>(observations) / coverage;
            auto mean_qual = var.MeanQuality();

            bool count_ok = (observations >= min_observations ||
                             (var.HasFwdAndRev() && observations >= min_observations_fwdrev));
            bool freq_ok = frequency >= min_frequency;
            // OR logic: passes if either quality gate holds
            bool quality_ok = mean_qual >= min_avg_quality || var.QualitySum() >= min_phred_sum;
            bool strand_ok = !require_strand || PassesStrandFilter(var, bin);

            return count_ok && freq_ok && quality_ok && strand_ok;
        }

        static void FilterSNPs(VariantBin& bin, size_t coverage, size_t min_observations, size_t min_observations_fwdrev, double min_frequency, size_t min_avg_quality, size_t min_phred_sum=0, bool require_strand=false) {
            for (auto& var : bin) {
                var.SetValid(FilterSNPs(var, bin, coverage, min_observations, min_observations_fwdrev, min_frequency, min_avg_quality, min_phred_sum, require_strand));
            }
        }

        // Completes a site's bin: adds the reference allele, with the reads of each strand that have a
        // base here (`forward`, `reverse`) less those carrying another base allele (a SNP or the
        // deletion), and marks each allele valid or not. Insertions sit before the base and are not
        // subtracted: their reads have a base here, too.
        static void PostProcessSNPBin(VariantBin& bin, size_t forward, size_t reverse, size_t min_observations=5, size_t min_observations_fwdrev=3, double min_frequency=0.2, size_t min_avg_quality=15, size_t min_phred_sum=0, bool require_strand=false) {
            auto var_pos = bin.front().Position();
            auto [allele_forward, allele_reverse] = SiteStrandCounts(bin);
            size_t const ref_forward = forward > allele_forward ? forward - allele_forward : 0;
            size_t const ref_reverse = reverse > allele_reverse ? reverse - allele_reverse : 0;

            // Reference allele (also stored in insertions, deleteions), unless every read carries a variant
            if (ref_forward + ref_reverse > 0) {
                char ref = bin.front().Reference();
                auto& variant = GetVariant(bin, var_pos, ref, ref);
                variant.SetObservations(ref_forward, ref_reverse);
            }

            std::stable_sort(bin.begin(), bin.end(), [](Variant const& a, Variant const& b) {
                return a.Observations() > b.Observations();
            });
            FilterSNPs(bin, forward + reverse, min_observations, min_observations_fwdrev, min_frequency, min_avg_quality, min_phred_sum, require_strand);
        }

        // Completes every bin (PostProcessSNPBin) from the reads per position and strand.
        void PostProcessSNPs(CoverageVec const& forward, CoverageVec const& reverse, size_t min_observations=2, size_t min_observations_fwdrev=2, double min_frequency=0.2, size_t min_avg_quality=15, size_t min_phred_sum=0, bool require_strand=false) {
            for (auto& [variant_pos, variant_bin] : m_variants) {
                auto& bin = m_variants.at(variant_pos);
                auto informative = [&](CoverageVec const& cov, int strand) -> size_t {
                    size_t const c = variant_pos < cov.size() ? cov[variant_pos] : 0;
                    size_t const none = NoBase(variant_pos, strand);
                    return c > none ? c - none : 0;
                };
                PostProcessSNPBin(bin, informative(forward, SequenceRange::kForward), informative(reverse, SequenceRange::kReverse),
                                  min_observations, min_observations_fwdrev, min_frequency, min_avg_quality, min_phred_sum, require_strand);
            }
        }
    };
}
