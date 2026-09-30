//
// Created by fritsche on 21/02/23.
//

#pragma once

#include <algorithm>
#include <atomic>
#include "SNPUtils.h"
#include "Variant.h"
#include "robin_map.h"
#include "GenomeLoader.h"
#include "SequenceRange.h"

namespace protal {
    class VariantHandler {
        Variants m_variants;
        tsl::robin_map<VariantPos, uint32_t> m_uncalled;  // reads with an N at a position
        std::string_view m_reference;  // the gene's sequence, which outlives the handler

    public:
        Benchmark bm_next_compressed_cigar{"Next compressed cigar"};

        explicit VariantHandler(std::string_view reference) : m_reference(reference) {};

        // Drops all variants and frees their memory.
        void Clear() {
            Variants{}.swap(m_variants);
            tsl::robin_map<VariantPos, uint32_t>{}.swap(m_uncalled);
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

        Variant& GetVariant(VariantBin& variant_bin, VariantPos pos, Base snp, Base ref) {
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

        Variant& GetVariant(VariantBin& variant_bin, VariantType type, VariantPos pos, Base ref, std::string& structural) {
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

        void AddSNP(VariantPos position, Base snp, Base ref, bool on_forward, Qual quality) {
            if (!HasVariantBin(position)) m_variants.insert({position, VariantBin() });
            auto& variant_bin = GetVariantBin(position);
            auto& variant = GetVariant(variant_bin, position, snp, ref);
            variant.AddObservation(quality, on_forward);
        }

        void AddINDEL(VariantType type, VariantPos position, Base ref, std::string&& structural, bool on_forward, Qual quality) {
            auto& variant_bin = HasVariantBin(position) ? GetVariantBin(position) : m_variants[position];
            auto& variant = GetVariant(variant_bin, type, position, ref, structural);
            variant.AddObservation(quality, on_forward);
        }

        // Records the variants one alignment supports. SNPs come from X ops (protal writes M only for
        // exact matches). A base that is N in the read or the reference supports no allele: it is
        // counted as uncalled, so that it is not taken for reference support either. Insertions and
        // deletions are called only between aligned bases: at either end of an alignment they are
        // artefacts, and a trailing insertion would lie past the read's coverage. A deletion has no
        // base of its own and gets the lower quality of its two flanking bases.
        bool ExtractVariants(SamEntry const& sam, size_t read_id = 0) {
            std::vector<std::pair<int, char>> ops;
            {
                int cpos = 0, count = 0;
                char op = ' ';
                while (NextCompressedCigar(cpos, sam.m_cigar, count, op)) ops.emplace_back(count, op);
            }
            auto is_aligned = [](std::pair<int, char> const& o) { return o.second == 'M' || o.second == 'X'; };
            std::ptrdiff_t const first_aligned = std::find_if(ops.begin(), ops.end(), is_aligned) - ops.begin();
            std::ptrdiff_t const last_aligned = static_cast<std::ptrdiff_t>(ops.size()) - 1 -
                                                (std::find_if(ops.rbegin(), ops.rend(), is_aligned) - ops.rbegin());

            bool const is_fwd = !sam.IsReversed();
            size_t qpos = 0;
            size_t rpos = sam.m_pos - 1;

            bm_next_compressed_cigar.Start();
            for (std::ptrdiff_t k = 0; k < static_cast<std::ptrdiff_t>(ops.size()); k++) {
                auto const [count, op] = ops[k];
                bool const consumes_ref = !(op == 'I' || op == 'S');
                bool const consumes_query = op != 'D';
                if ((consumes_ref && rpos + count > m_reference.size()) || (consumes_query && qpos + count > sam.m_seq.size())) {
                    bm_next_compressed_cigar.Stop();
                    return false;
                }
                bool const between_aligned_bases = k > first_aligned && k < last_aligned;

                if (op == 'M') {
                    for (auto i = 0; i < count; i++) {
                        char const base = sam.m_seq[qpos + i];
                        char const ref = m_reference[rpos + i];
                        if (base == 'N' || ref == 'N') {
                            m_uncalled[rpos + i]++;
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
                            return false;
                        }
                    }
                } else if (op == 'X') {
                    for (auto i = 0; i < count; i++) {
                        char const base = sam.m_seq[qpos + i];
                        char const ref = m_reference[rpos + i];
                        if (base == 'N' || ref == 'N') {
                            m_uncalled[rpos + i]++;
                        } else {
                            AddSNP(rpos + i, base, ref, is_fwd, PhredScore(sam.m_qual[qpos + i]));
                        }
                    }
                } else if (op == 'I' && between_aligned_bases) {
                    size_t qual_sum = 0;
                    for (auto i = 0; i < count; i++) {
                        qual_sum += PhredScore(sam.m_qual[qpos + i]);
                    }
                    AddINDEL(VariantType::INS, rpos, m_reference[rpos], sam.m_seq.substr(qpos, count), is_fwd, qual_sum / count);
                } else if (op == 'D' && between_aligned_bases) {
                    Qual const flank = std::min(PhredScore(sam.m_qual[qpos - 1]), PhredScore(sam.m_qual[qpos]));
                    AddINDEL(VariantType::DEL, rpos, m_reference[rpos], std::string(m_reference.substr(rpos, count)), is_fwd, flank);
                }

                qpos += consumes_query * count;
                rpos += consumes_ref * count;
            }
            bm_next_compressed_cigar.Stop();
            (void)read_id;
            return true;
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

        bool FilterSNPs(Variant& var, size_t coverage, size_t min_observations, size_t min_observations_fwdrev, double min_frequency, size_t min_avg_quality, size_t min_phred_sum=0, bool require_strand=false) {
            auto observations = var.Observations();
            auto frequency = coverage == 0 ? 0.0 : static_cast<double>(observations) / coverage;
            auto mean_qual = var.MeanQuality();

            bool count_ok = (observations >= min_observations ||
                             (var.HasFwdAndRev() && observations >= min_observations_fwdrev));
            bool freq_ok = frequency >= min_frequency;
            // OR logic: passes if either quality gate holds
            bool quality_ok = mean_qual >= min_avg_quality || var.QualitySum() >= min_phred_sum;
            bool strand_ok = !require_strand || var.PassesStrandFilter();

            return count_ok && freq_ok && quality_ok && strand_ok;
        }

        void FilterSNPs(VariantBin& bin, size_t coverage, size_t min_observations, size_t min_observations_fwdrev, double min_frequency, size_t min_avg_quality, size_t min_phred_sum=0, bool require_strand=false) {
            for (auto& var : bin) {
                auto is_valid = FilterSNPs(var, coverage, min_observations, min_observations_fwdrev, min_frequency, min_avg_quality, min_phred_sum, require_strand);
                var.SetValid(is_valid);
            }
        }

        void PostProcessSNPBin(VariantBin& bin, size_t coverage, size_t min_observations=5, size_t min_observations_fwdrev=3, double min_frequency=0.2, size_t min_avg_quality=15, size_t min_phred_sum=0, bool require_strand=false) {
            auto var_pos = bin.front().Position();

            // Total variant var_observations
            auto var_observations = std::accumulate(bin.begin(), bin.end(), size_t{0}, [](size_t acc, Variant const & a) {
                return acc + a.Observations();
            });

            // Sort SNPs based on
            std::sort(bin.begin(), bin.end(), [](Variant const& a, Variant const& b) {
                return a.Observations() > b.Observations();
            });

            // Average quality for bin.
            auto qual = std::accumulate(bin.begin(), bin.end(), size_t{0}, [](size_t acc, Variant const & a) {
                return acc + a.MeanQuality();
            });
            qual /= bin.size();

//            if (coverage < var_observations) {
//                std::cout << "Cov: " << coverage << std::endl;
//                std::cout << "Obs: " << var_observations << std::endl;
//                for (auto& var : bin) {
//                    std::cout << var.ToString() << std::endl;
//                }
//                std::cout << coverage - var_observations << std::endl;
//                Utils::Input();
//            }

            // Reference allele (also stored in insertions, deleteions), unless every read carries a variant
            if (coverage > var_observations) {
                char ref = bin.front().Reference();

                auto& variant = GetVariant(bin, var_pos, ref, ref);
                variant.SetObservations(coverage - var_observations);
            }

            FilterSNPs(bin, coverage, min_observations, min_observations_fwdrev, min_frequency, min_avg_quality, min_phred_sum, require_strand);
        }

        void PostProcessSNPs(CoverageVec const& coverage, size_t min_observations=2, size_t min_observations_fwdrev=2, double min_frequency=0.2, size_t min_avg_quality=15, size_t min_phred_sum=0, bool require_strand=false) {
            // Iterate all variant positions.
            for (auto& [variant_pos, variant_bin] : m_variants) {
                auto& bin = m_variants.at(variant_pos);
                size_t cov = variant_pos < coverage.size() ? coverage[variant_pos] : 0;
                // Reads with an N here cover the position but support no allele.
                auto uncalled = m_uncalled.find(variant_pos);
                if (uncalled != m_uncalled.end()) cov = cov > uncalled->second ? cov - uncalled->second : 0;

                PostProcessSNPBin(bin, cov, min_observations, min_observations_fwdrev, min_frequency, min_avg_quality, min_phred_sum, require_strand);
            }
        }

        bool AddVariantsFromSam(SamEntry const& sam, size_t read_id = 0) {
            return ExtractVariants(sam, read_id);
        }

    };
}
