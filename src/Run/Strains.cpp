// The strain-level MSAs over the samples' profiles, the strain stage of protal::Run (ProfileStage.h).

#include <algorithm>
#include <atomic>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>
#include <sys/wait.h>
#include "ProfileStage.h"
#include "Benchmark.h"
#include "ProgressBar.h"
#include "RunStatus.h"

namespace protal {

    static std::pair<bool,bool> SharedSNP(VariantBin& a, VariantBin& b) {
        std::sort(a.begin(), a.end(), [](Variant const& va, Variant const& vb) {
            return va.Observations() > vb.Observations();
        });
        std::sort(b.begin(), b.end(), [](Variant const& va, Variant const& vb) {
            return va.Observations() > vb.Observations();
        });

        bool shared_non_ref = false;
        bool shared_ref = false;
        for (auto& var_a : a) {
            if (!var_a.GetValid()) continue;
            for (auto& var_b : b) {
                if (!var_b.GetValid()) continue;

                if (var_a.Match(var_b))  {
                    shared_ref |= var_a.IsReference() && var_b.IsReference();
                    shared_non_ref |= !var_a.IsReference() && !var_b.IsReference();
                }
            }
        }

        return { shared_ref, shared_non_ref };
    }

    static double JukesCantor(size_t mismatch, size_t length) {
        auto strain_identity = static_cast<double>(mismatch) / length;

        auto jc69 = -3.0f/4 * log(1 - (4.0f/3 * strain_identity));
        return jc69 > 1 || jc69 != jc69 ? 1 : jc69;
    }

    static bool IsSNP(VariantBin& bin, size_t min_obs, double min_frequency) {
        for (auto& var : bin) {
            size_t total_obs = std::accumulate(bin.begin(), bin.end(), size_t{0}, [](size_t acc, Variant const& var) { return acc + var.Observations(); });
            auto obs = var.Observations();
            double freq = static_cast<double>(obs) / total_obs;
            if (!var.IsReference() && freq >= min_frequency && (obs >= min_obs || var.HasFwdAndRev())) {
                return true;
            }
        }
//
//        for (auto& var : bin) {
//            std::cout << var.ToString() << '\t';
//        }
//        std::cout << endl;
//
//        Utils::Input();
        return false;
    }

    static bool HasReference(VariantBin& bin) {
        return std::any_of(bin.begin(), bin.end(), [](Variant const& var) {
            return var.IsReference();
        });
    }

    static std::pair<size_t, size_t> Distance(size_t shared_length, VariantVec& a, VariantVec& b) {
        const auto pos_getter = [](VariantBin const& bin) { return bin.front().Position(); };
        size_t shared_snps = 0;
        size_t unique_snps_a = 0;
        size_t unique_snps_b = 0;

        size_t min_observations = 2;
        double min_frequency = 0.2;

        int ai = 0, bi = 0;
//        std::cout << "a: " << a.size() << "  b: " << b.size() << std::endl;
        while (ai < a.size() && bi < b.size()) {
            auto pos_a = pos_getter(a[ai]);
            auto pos_b = pos_getter(b[bi]);

            bool has_ref_a = std::any_of(a[ai].begin(), a[ai].end(), [](Variant const& v) {return v.GetValid() && v.IsReference();});
            bool has_ref_b = std::any_of(b[bi].begin(), b[bi].end(), [](Variant const& v) {return v.GetValid() && v.IsReference();});
            bool has_snp_a = std::any_of(a[ai].begin(), a[ai].end(), [](Variant const& v) {return v.GetValid() && !v.IsReference();});
            bool has_snp_b = std::any_of(b[bi].begin(), b[bi].end(), [](Variant const& v) {return v.GetValid() && !v.IsReference();});
//            bool has_snp_b = IsSNP(b[bi], min_observations, min_frequency);

            bool same_pos = (pos_a == pos_b);
            bool advance_a = (pos_a <= pos_b);
            bool advance_b = (pos_b <= pos_a);

//            std::cout << ai << " " << pos_a << "\t->\t" << a[ai].front().ToString() << "\t" << b[bi].front().ToString() << "\t<-\t" << pos_b << " " << bi << std::endl;

            // same_pos:                 required because a shared snp needs variants detected at the same pos.
            // has_snp_a || has_snp_b:   if both have variants at the position but none are considered
            //                           genuine snps, it is assumed that both are reference alleles
            // SharedSNP(a[ai], b[bi]):  Check if they share at least one valid snp, ref or non ref.
            auto [has_shared_ref, has_shared_snp] = same_pos ? SharedSNP(a[ai], b[bi]) : std::pair<bool,bool>(false, false);
            bool is_shared_and_snp = has_shared_snp && !has_shared_ref;
            shared_snps += is_shared_and_snp;
            // has_snp_a: if all variants at this position are dodgy/fp SNPs assume its the reference allele
            // same_po
            // advance var is needed because only the smaller pos can be a snp
            bool is_unique_snp_a = (advance_a && !has_shared_ref && has_snp_a && !has_ref_a && !has_shared_snp);
            bool is_unique_snp_b = (advance_b && !has_shared_ref && has_snp_b && !has_ref_b && !has_shared_snp);
            unique_snps_a += is_unique_snp_a;
            unique_snps_b += is_unique_snp_b;

//            if (is_unique_snp_a || unique_snps_b) {
//                std::cout << std::string(60, '#') << std::endl;
//                std::cout << (is_unique_snp_a ? "SNP IN A" : "SNP IN B") << std::endl;
//
//                std::cout << "____________SNPS A noref:" << !has_ref_a << std::endl;
//                for (auto& var : a[ai]) {
//                    std::cout << var.ToString() << '\t';
//                }
//                std::cout << std::endl;
//                std::cout << "____________SNPS B noref:" << !has_ref_b << std::endl;
//                for (auto& var : b[bi]) {
//                    std::cout << var.ToString() << '\t';
//                }
//                std::cout << std::endl;
//                Utils::Input();
//            }

//            if ((has_snp_a && advance_a && !shared) || (has_snp_b && advance_b && !shared)) {
//                std::cout << "PosA: " << pos_a << "\tPosB:" << pos_b << std::endl;
//                std::cout << "SNPA: " << has_snp_a << "\tSNPB:" << has_snp_b << std::endl;
//                std::cout << "shared: " << shared << std::endl;
//                std::cout << "-------\nA is SNP? " << has_snp_a << std::endl;
//                for (auto& var : a[ai]) {
//                    std::cout << var.ToString() << '\t';
//                }
//                std::cout << endl;
//                std::cout << "B is SNP? " << has_snp_b << std::endl;
//                for (auto& var : b[bi]) {
//                    std::cout << var.ToString() << '\t';
//                }
//                std::cout << endl;
//
//                Utils::Input();
//            }

            ai += advance_a;
            bi += advance_b;
        }

        if (a.size() - ai > 0) {
//            std::cout << "A) SNPS_A: " << unique_snps_a << "\tSNPS_B: " << unique_snps_b << std::endl;
            for (; ai < a.size(); ai++) {
                bool has_var = std::any_of(a[ai].begin(), a[ai].end(), [](Variant const& v) {return v.GetValid() && !v.IsReference();});
                bool has_ref = std::any_of(a[ai].begin(), a[ai].end(), [](Variant const& v) {return v.GetValid() && v.IsReference();});
                unique_snps_a += (has_var && !has_ref);
            }
//            std::cout << "  -->  SNPS_A: " << unique_snps_a << "\tSNPS_B: " << unique_snps_b << std::endl;
        }

        if (b.size() - bi > 0) {
//            std::cout << "B) SNPS_A: " << unique_snps_a << "\tSNPS_B: " << unique_snps_b << std::endl;
            for (; bi < b.size(); bi++) {
                bool has_var = std::any_of(b[bi].begin(), b[bi].end(), [](Variant const& v) {return v.GetValid() && !v.IsReference();});
                bool has_ref = std::any_of(b[bi].begin(), b[bi].end(), [](Variant const& v) {return v.GetValid() && v.IsReference();});
                unique_snps_b += (has_var && !has_ref);
            }
        }

//        if (unique_snps_a + unique_snps_b > 0) {
//            Utils::Input();
//        }

        return std::pair{ unique_snps_a + unique_snps_b, shared_length };

//        return JukesCantor(unique_snps_a + unique_snps_b, shared_length);
    }

    static std::pair<double, size_t> Distance(profiler::Gene& g1, profiler::Gene& g2) {
        auto& strain1 = g1.GetStrainLevel();
        auto& strain2 = g2.GetStrainLevel();
        auto shared_alignment_region = SharedAlignmentRegion::GetSharedAlignmentRegion(strain1, strain2);

        auto& va = shared_alignment_region.variant_handler_a;
        auto& vb = shared_alignment_region.variant_handler_b;
        auto shared_sequence_length = shared_alignment_region.share_range.SequenceLength();


        auto [unique_snps, shared_length] = Distance(shared_sequence_length, va, vb);

        return {unique_snps, shared_sequence_length};
    }

    static std::pair<double, size_t> Distance(profiler::Taxon& t1, profiler::Taxon& t2) {

        double distance_sum = 0;
        size_t total_unique_snps = 0;
        size_t shared_region_sum = 0;

        for (auto& [gid, _] : t1.GetGenes()) {
            if (t2.GetGenes().contains(gid)) {
                //std::cout << "Stuck at " << gid << std::endl;
                auto [unique_snps, shared_region] = Distance(t1.GetGenes().at(gid), t2.GetGenes().at(gid));

                total_unique_snps += unique_snps;
                shared_region_sum += shared_region;
            }
        }

        auto distance = JukesCantor(total_unique_snps, shared_region_sum);

//        if (distance < 0.001) {
//            std::cout << "Distance: " << distance << " with " << total_unique_snps << " snps over " << shared_region_sum << " bases." << std::endl;
//        }

        return { distance, shared_region_sum };
    }


    using StrainResults = tsl::robin_map<size_t, Matrix<double>>;
    using TaxidCounts = tsl::robin_map<uint32_t, uint32_t>;
    using TaxidSet = tsl::robin_set<uint32_t>;
    using OptionalFilter = optional<std::reference_wrapper<const profiler::TaxonFilter>>;
    using SimilarityMatrix = DoubleMatrix;

    // Whether `taxon` of `profile` enters the strain MSAs: the model `filter` scores it at the profile's MSA knob or more
    // (--msa_knob, else the threshold the sample was reported at: --knob, or its model's knob for its depth or its
    // share of false calls), and the singleton rule does not veto it.
    static bool EntersMSA(profiler::TaxonFilterObj const& filter, Profile const& profile, profiler::Taxon const& taxon) {
        return filter.Calls(taxon, profile.MSAKnob());
    }

    // Taxa in at least `min_samples` profiles (entering the MSAs by `filter`, if given), the most frequent first.
    TaxidList ExtractTaxa(Profiles const& profiles, std::optional<profiler::TaxonFilterObj> filter, size_t min_samples) {
        TaxidCounts taxid_counts;
        TaxidSet taxa;
        TaxidList taxid_list;

        for (auto& profile : profiles) {
            for (auto& [id, taxon] : profile.GetTaxa()) {
                if (filter.has_value() && !EntersMSA(*filter, profile, taxon)) continue;
                if (!taxid_counts.contains(id)) {
                    taxid_counts.insert({ id, 0 });
                }
                taxid_counts.at(id)++;
            }
        }

        for (auto& [t, c] : taxid_counts) {
            if (c >= min_samples) {
                taxid_list.emplace_back(t);
            }
        }

        // Ties by taxid: the counts' map iterates in the order the taxa were met.
        std::sort(taxid_list.begin(), taxid_list.end(), [&taxid_counts](uint32_t const& t1, uint32_t const& t2) {
            auto const c1 = taxid_counts.at(t1), c2 = taxid_counts.at(t2);
            return c1 != c2 ? c1 > c2 : t1 < t2;
        });

        return taxid_list;
    }

    std::vector<uint32_t> ResolveMSASpecies(Options& options, taxonomy::IntTaxonomy& taxonomy) {
        const auto& msa_species = options.GetMSASpecies();
        if (msa_species.empty()) return {};

        std::vector<uint32_t> taxids;
        taxids.reserve(msa_species.size());
        std::vector<std::string> invalid;
        invalid.reserve(msa_species.size());

        for (const auto& raw_spec : msa_species) {
            std::string spec = raw_spec;
            std::replace(spec.begin(), spec.end(), ' ', '_');
            if (spec.rfind("s__", 0) != 0 || spec.find('_', 3) == std::string::npos) {
                invalid.emplace_back(raw_spec);
                continue;
            }

            uint32_t taxid = 0;
            bool found = false;
            if (taxonomy.string_to_id.contains(spec)) {
                taxid = taxonomy.Get(spec);
                found = true;
            } else {
                std::string alt = spec;
                for (size_t i = 3; i < alt.size(); i++) {
                    if (alt[i] == '_') alt[i] = ' ';
                }
                if (taxonomy.string_to_id.contains(alt)) {
                    taxid = taxonomy.Get(alt);
                    found = true;
                }
            }

            if (!found) {
                invalid.emplace_back(raw_spec);
                continue;
            }
            if (taxonomy.Get(taxid).rank != "species") {
                invalid.emplace_back(raw_spec);
                continue;
            }
            taxids.emplace_back(taxid);
        }

        if (!invalid.empty()) {
            std::cerr << "Invalid --msa_species entries (expected format s__Genus_species and present in taxonomy): ";
            std::cerr << Utils::join(invalid, ",") << std::endl;
            exit(2);
        }

        return taxids;
    }

    static std::pair<double, size_t> GetSimilarity(uint32_t taxid, Profile& profile1, Profile& profile2, Options& options, std::optional<profiler::TaxonFilterObj> filter={}) {
        size_t min_shared_region = 1000;
        if (profile1.GetTaxa().contains(taxid) && profile2.GetTaxa().contains(taxid)) {
            auto& taxon1 = profile1.GetTaxa().at(taxid);
            auto& taxon2 = profile2.GetTaxa().at(taxid);

            if (filter.has_value() && !(EntersMSA(*filter, profile1, taxon1) && EntersMSA(*filter, profile2, taxon2))) {
                return { NAN, 0 };
            }

            auto [distance, shared_region] = Distance(taxon1, taxon2);

            auto similarity = 1.0f-distance;

            if (shared_region <= min_shared_region) {
                similarity = NAN;
            }
            return { similarity, shared_region };
        }
        return { NAN, 0 };
    }

    static bool IsRowGood(std::vector<char> row, size_t min_hcov) {
        size_t count_non_n = 0;
        for (auto c : row) {
            count_non_n += (c != 'N' && c != '-');
        }
        return count_non_n >= min_hcov;
    }

    static std::vector<uint32_t> GetInformationVector(MSAVector const& msa_vector) {
        size_t info = 0;
        const auto msa_len = msa_vector.front().size();
        std::vector<uint32_t> information_vector(msa_len, 0);

        for (auto pos = 0; pos < msa_len; pos++) {
            info = 0;
            for (auto const& seq : msa_vector) info += (seq[pos] != 'N' && seq[pos] != '-');
            information_vector[pos] = info;
        }
        return information_vector;
    }

    static MSAVector ProcessMSA(MSAVector const& msa_vector, double position_coverage = 0.5) {
        auto sample_len = msa_vector.size();
        auto info_vector = GetInformationVector(msa_vector);
        MSAVector processed_msa(sample_len, std::vector<char>{});

        auto sample_size_threshold = msa_vector.size() * position_coverage;

        // std::cout << "Sample size threshold: " << sample_size_threshold << std::endl;
        
        for (auto pos = 0; pos < info_vector.size(); pos++) {
            if (info_vector[pos] > sample_size_threshold) {
                for (auto s = 0; s < sample_len; s++) {
                    processed_msa[s].emplace_back(msa_vector[s][pos]);
                }
            }
        }
        return processed_msa;
    }

    static void OutputMSA(MSAVector const& msa, std::vector<std::string> const& names, std::ostream& os=std::cout) {
        for (auto i = 0; i < msa.size(); i++) {
            auto& row = msa[i];

            if (row.empty()) continue;
            os << ">" << names[i] << std::endl;
            // os << std::string_view(&row[0], std::distance(row.begin(), row.end())) << std::endl;
            os << std::string(&row[0], std::distance(row.begin(), row.end())) << std::endl;
//            os << std::string_view(row.begin(), row.end()) << std::endl;
        }
    }

    static std::vector<size_t> GetProfilesWithTaxon(uint32_t taxid, Profiles& profiles, Options& options, std::optional<profiler::TaxonFilterObj> const& filter) {
        std::vector<size_t> indices;

        for (auto i = 0; i < profiles.size(); i++) {
            auto& profile = profiles[i];
            if (!profile.GetTaxa().contains(taxid)) {
                continue;
            }
            auto& taxon_map = profile.GetTaxa();

            auto& taxon = taxon_map.at(taxid);
            if (filter.has_value() && !EntersMSA(filter.value(), profile, taxon)) {
                continue;
            }
            // std::cout << taxid << " Passes " << std::endl;

            indices.emplace_back(i);
        }
        return indices;
    }

    static std::vector<uint32_t> SelectGenesForTaxon(uint32_t taxid, std::string name, std::vector<size_t>& selected_profiles, GenomeLoader& loader, Options& options, Profiles& profiles, std::ostream& out = std::cout) {
        std::vector<uint32_t> selected_gene_ids;
        // The genes with reads in any of the samples. A gene without unique k-mers has reads, too (its
        // k-mers are shared, not absent); a gene without reads would add only gaps.
        std::vector<uint32_t> gene_ids;
        {
            std::set<uint32_t> observed;
            for (auto sample_index : selected_profiles) {
                auto const& taxa = profiles[sample_index].GetTaxa();
                if (!taxa.contains(taxid)) continue;
                for (auto const& [gene_id, _] : taxa.at(taxid).GetGenes()) observed.insert(gene_id);
            }
            gene_ids.assign(observed.begin(), observed.end());
        }
        if (gene_ids.empty()) return gene_ids;

        auto max_gene_id = std::max_element(gene_ids.begin(), gene_ids.end());

        // std::vector<std::vector<double>> per_sample_gene_multiallelic_portion;
        // per_sample_gene_multiallelic_portion.resize(selected_profiles.size(),
        //     std::vector<double>(*max_gene_id + 1, -1) );

        std::vector<std::vector<int>> per_sample_gene_multiallelic_snps;
        per_sample_gene_multiallelic_snps.resize(selected_profiles.size(),
            std::vector<int>(*max_gene_id + 1, -2) );

        std::vector<std::vector<int>> per_sample_gene_snps;
        per_sample_gene_snps.resize(selected_profiles.size(),
            std::vector<int>(*max_gene_id + 1, -2) );

        std::vector<std::vector<int>> per_sample_gene_cov;
        per_sample_gene_cov.resize(selected_profiles.size(),
            std::vector<int>(*max_gene_id + 1, -2) );

        std::vector<std::vector<int>> per_sample_gene_noise;
        per_sample_gene_noise.resize(selected_profiles.size(),
            std::vector<int>(*max_gene_id + 1, -2) );

        for (auto si = 0; si < selected_profiles.size(); si++) {
            auto sample_index = selected_profiles[si];
            auto& taxon = profiles[sample_index].GetTaxa().at(taxid);
            auto& multiallelic_vec = per_sample_gene_multiallelic_snps[si];

            // std::cout << taxon.GetName() << " -> hittable genes" << gene_ids.size() << std::endl;
            for (auto& gene_id : gene_ids) {
                if (gene_id >= per_sample_gene_multiallelic_snps[si].size()) {
                    out << gene_id << " >= " << gene_ids.size() << std::endl;
                    exit(123);
                }

                per_sample_gene_multiallelic_snps[si][gene_id] = -1;
                per_sample_gene_snps[si][gene_id] = -1;
                per_sample_gene_cov[si][gene_id] = -1;
                per_sample_gene_noise[si][gene_id] = -1;
                
                if (!taxon.GetGenes().contains(gene_id)) {
                    per_sample_gene_multiallelic_snps[si][gene_id] = 0;
                    per_sample_gene_snps[si][gene_id] = 0;
                    per_sample_gene_cov[si][gene_id] = 0;
                    per_sample_gene_noise[si][gene_id] = 0;
                    continue;
                }

                auto& gene = taxon.GetGene(gene_id);

                auto allele_counts = gene.AlleleSNPCounts(2, 60);
                auto allele_counts_nf = gene.AlleleSNPCounts(0, 0);

                // std::cout << "----Filter----\n" << allele_counts.ToString() << std::endl;
                // std::cout << "----No Filter----\n" << allele_counts_nf.ToString() << std::endl;
                
                per_sample_gene_multiallelic_snps[si][gene_id] = allele_counts.Multi();
                per_sample_gene_snps[si][gene_id] = allele_counts.AllValid();
                per_sample_gene_cov[si][gene_id] = gene.Coverage();
                per_sample_gene_noise[si][gene_id] = allele_counts.Noisy();
            }
        }

        std::ofstream os(options.GetMiscOutputDir() + '/' + name + ".snps_multiallelic.tsv", std::ios::out);
        for (auto si = 0; si < selected_profiles.size(); si++) {
            auto sample_index = selected_profiles[si];
            os << options.GetSampleId(sample_index);

            if (si >= per_sample_gene_multiallelic_snps.size()) exit(244);
            auto& vec = per_sample_gene_multiallelic_snps[si];

            double total = std::accumulate(vec.begin(), vec.end(), 0.0, [](double acc, int val) {
                return acc + (val < 0 ? 0 : val);
            });

            os << '\t' << total;
            for (auto i = 1; i < vec.size(); i++) {
                os << '\t' << vec[i];
            }
            os << std::endl;
        }
        os.close();

        std::ofstream os2(options.GetMiscOutputDir() + '/' + name + ".snps_total.tsv", std::ios::out);
        for (auto si = 0; si < selected_profiles.size(); si++) {
            auto sample_index = selected_profiles[si];
            os2 << options.GetSampleId(sample_index);

            if (si >= per_sample_gene_snps.size()) exit(244);
            auto& vec = per_sample_gene_snps[si];

            double total = std::accumulate(vec.begin(), vec.end(), 0.0, [](double acc, int val) {
                return acc + (val < 0 ? 0 : val);
            });

            os2 << '\t' << total;
            for (auto i = 1; i < vec.size(); i++) {
                os2 << '\t' << vec[i];
            }
            os2 << std::endl;
        }
        os2.close();

        std::ofstream os3(options.GetMiscOutputDir() + '/' + name + ".hcov.tsv", std::ios::out);
        for (auto si = 0; si < selected_profiles.size(); si++) {
            auto sample_index = selected_profiles[si];
            os3 << options.GetSampleId(sample_index);

            if (si >= per_sample_gene_snps.size()) exit(244);
            auto& vec = per_sample_gene_cov[si];
            
            double total = std::accumulate(vec.begin(), vec.end(), 0.0, [](double acc, int val) {
                return acc + (val < 0 ? 0 : val);
            });

            os3 << '\t' << total;
            for (auto i = 1; i < vec.size(); i++) {
                os3 << '\t' << vec[i];
            }
            os3 << std::endl;
        }
        os3.close();


        std::ofstream os4(options.GetMiscOutputDir() + '/' + name + ".snps_filtered.tsv", std::ios::out);
        for (auto si = 0; si < selected_profiles.size(); si++) {
            auto sample_index = selected_profiles[si];
            os4 << options.GetSampleId(sample_index);

            if (si >= per_sample_gene_noise.size()) exit(244);
            auto& vec = per_sample_gene_noise[si];

            double total = std::accumulate(vec.begin(), vec.end(), 0.0, [](double acc, int val) {
                return acc + (val < 0 ? 0 : val);
            });

            os4 << '\t' << total;
            for (auto i = 1; i < vec.size(); i++) {
                os4 << '\t' << vec[i];
            }
            os4 << std::endl;
        }

        os4.close();

        // Long unique k-mers are unique k-mers in index blocks with flex keys, i.e. whose 15-mer core
        // is shared with a relative. A species with relatives in the database has them in nearly
        // every gene; a gene without any is one its relatives share unchanged or lack from their
        // reference. A relative's reads of such a gene align here, and its MSA columns would show
        // them as a second strain, so it is left out. A species without relatives has fewer, as its
        // k-mers are mostly unique at the core already: in the 2026-09-29 strain audit, Dummya solo (no
        // relative in the database) had long unique k-mers in half its genes, species with congeners
        // in 98-99%. Unless 90% of a species' genes have long unique k-mers, all its genes stay.
        auto& genome = loader.GetGenome(taxid);
        auto const with_long_uniques = std::count_if(gene_ids.begin(), gene_ids.end(), [&genome](uint32_t gene_id) {
            return genome.GetGene(gene_id).HasLongUniques();
        });
        if (static_cast<size_t>(with_long_uniques) * 10 < gene_ids.size() * 9) return gene_ids;

        std::vector<uint32_t> msa_gene_ids;
        std::copy_if(gene_ids.begin(), gene_ids.end(), std::back_inserter(msa_gene_ids), [&genome](uint32_t gene_id) {
            return genome.GetGene(gene_id).HasLongUniques();
        });
        if (msa_gene_ids.size() < gene_ids.size()) {
            out << name << ": " << gene_ids.size() - msa_gene_ids.size() << " of " << gene_ids.size()
                      << " genes have no long unique k-mers and are left out of the MSA" << std::endl;
        }
        return msa_gene_ids;
    }

    // Shell-quote a path/argument for use in a std::system command line.
    static std::string ShellQuote(const std::string& s) {
        std::string out = "'";
        for (char c : s) {
            if (c == '\'') out += "'\\''";
            else out += c;
        }
        out += "'";
        return out;
    }

    static bool IsExecutableFile(std::filesystem::path const& p) {
        std::error_code ec;
        return std::filesystem::is_regular_file(p, ec) && ::access(p.c_str(), X_OK) == 0;
    }

    // Look up an executable by name in $PATH, the way a shell would.
    static std::string FindInPath(std::string const& name) {
        const char* path_env = std::getenv("PATH");
        if (!path_env || !*path_env) return "";
        std::string const path(path_env);
        for (size_t start = 0; start <= path.size(); ) {
            size_t end = path.find(':', start);
            if (end == std::string::npos) end = path.size();
            if (end > start) {
                std::filesystem::path cand = std::filesystem::path(path.substr(start, end - start)) / name;
                if (IsExecutableFile(cand)) return cand.string();
            }
            start = end + 1;
        }
        return "";
    }

    // Locate the qcmsa post-filter. qcmsa ships as an executable named 'qcmsa'
    // installed next to the protal binary (see `just install`), so that is what
    // we look for: explicit --qcmsa_script, then PROTAL_QCMSA_SCRIPT, then
    // 'qcmsa' next to the protal executable, then 'qcmsa' on $PATH. The
    // qcmsa.py script of a source checkout is only a last resort so that running
    // protal straight out of the build tree keeps working.
    // Returns "" if not found.
    static std::string FindQCMSAScript(Options& options) {
        namespace fs = std::filesystem;

        if (!options.GetQCMSAScript().empty()) {
            std::string const& given = options.GetQCMSAScript();
            // A bare command name is resolved via $PATH, a path is taken as is.
            if (given.find('/') == std::string::npos) {
                if (std::string in_path = FindInPath(given); !in_path.empty()) return in_path;
            }
            return given;
        }
        if (const char* env = std::getenv("PROTAL_QCMSA_SCRIPT"); env && *env)
            return std::string(env);

        std::error_code ec;
        fs::path exe = fs::canonical("/proc/self/exe", ec);
        if (!ec && IsExecutableFile(exe.parent_path() / "qcmsa"))
            return (exe.parent_path() / "qcmsa").string();

        if (std::string in_path = FindInPath("qcmsa"); !in_path.empty()) return in_path;

        if (!ec) {
            fs::path dir = exe.parent_path();
            for (const auto& cand : { dir / "qcmsa.py",                       // alongside the binary
                                      dir / "scripts" / "qcmsa.py",
                                      dir.parent_path() / "scripts" / "qcmsa.py" }) {
                if (fs::exists(cand)) return cand.string();
            }
        }
        return "";
    }

    // M5 step 4c: invoke the qcmsa post-filter on a species' MSA. A failure does not stop the run --
    // the raw strain outputs protal already wrote remain valid -- but it is reported and makes protal
    // exit non-zero, because the filtered <name>.msa.fna the run was asked for is missing.
    // log: where its messages go. qcmsa's own output (progress bars, its report per species) is not shown: it goes to a file
    // beside its outputs, which is printed (its last lines, with the command) only if qcmsa fails, and removed. parallel
    // (the strain stage runs species at once, StrainWrapper2): qcmsa's numerical libraries run on one thread each.
    static void RunQCMSA(Options& options, const std::string& name, std::ostream& log = std::cerr, bool parallel = false) {
        namespace fs = std::filesystem;
        std::string script = FindQCMSAScript(options);
        if (script.empty() || !fs::exists(script)) {
            static std::atomic<bool> reported{false};
            if (!reported.exchange(true)) {
                RunStatus::Get().Fail("qcmsa not found next to the protal binary or on $PATH, so no species got "
                                      "its filtered .msa.fna (install it with 'just install', point protal at it "
                                      "with --qcmsa_script / PROTAL_QCMSA_SCRIPT, or pass --no_qcmsa)");
            }
            return;
        }
        // qcmsa is an executable with a python3 shebang. Only fall back to
        // running it through the interpreter if it is not executable (e.g. a
        // qcmsa.py picked straight out of a source checkout).
        std::string launcher = IsExecutableFile(script) ? "" : "python3 ";
        // Run qcmsa on protal's native (raw) MSA; it writes the final <name>.msa.fna.
        std::string msa = options.GetMSAOutput(name);              // .raw.msa.fna
        std::string partition = options.GetMSAPartitionOutput(name); // .raw.partition.txt
        std::string meta = options.GetSpeciesMetaOutput(name);
        if (!fs::exists(msa) || !fs::exists(partition) || !fs::exists(meta)) {
            log << "[qcmsa] WARNING: missing MSA/partition/meta for " << name
                << "; skipping post-filter" << std::endl;
            return;
        }
        std::string prefix = options.GetStrainOutputDir() + '/' + name;
        // Remove qcmsa outputs of an earlier run first: if qcmsa writes nothing this time (e.g. no
        // sample passes its filters), a stale <name>.msa.fna must not be left behind as if it were new.
        for (auto const* ext : { ".msa.fna", ".partition.txt", ".qcmsa_summary.tsv", ".qc.png" }) {
            std::error_code ec;
            fs::remove(prefix + ext, ec);
        }
        std::ostringstream cmd;
        if (parallel) cmd << "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 ";
        // --scratch_dir: the temporary files of qcmsa's libraries in the run's scratch folder.
        if (!options.ScratchFolder().empty()) cmd << "TMPDIR=" << ShellQuote(options.ScratchFolder()) << ' ';
        cmd << launcher << ShellQuote(script)
            << ' ' << ShellQuote(msa)
            << ' ' << ShellQuote(partition)
            << ' ' << ShellQuote(meta)
            << " --prefix " << ShellQuote(prefix)
            << " --reapply-hcov " << options.GetMSAMinHCOV();
        // User-supplied flags go last so they win over the defaults protal passes
        // above. Forwarded verbatim (unquoted) -- they are a flag list, not a value.
        if (!options.GetQCMSAArgs().empty()) cmd << ' ' << options.GetQCMSAArgs();
        std::string const command = cmd.str();
        // What qcmsa prints, read and removed below: in the run's scratch folder with --scratch_dir.
        std::string const output = options.ScratchFolder().empty() ? prefix + ".qcmsa_output.partial"
                                 : (fs::path(options.ScratchFolder()) / (name + ".qcmsa_output.partial")).string();
        cmd << " > " << ShellQuote(output) << " 2>&1";
        int const status = std::system(cmd.str().c_str());
        std::string text;
        {
            std::ifstream is(output, std::ios::binary);
            text.assign(std::istreambuf_iterator<char>(is), std::istreambuf_iterator<char>());
            std::error_code ec;
            fs::remove(output, ec);
        }
        if (status != 0) {
            // std::system returns a wait status: 512 means exit code 2.
            std::string const how = status == -1 ? "could not be started"
                                  : WIFEXITED(status) ? "exited with code " + std::to_string(WEXITSTATUS(status))
                                  : WIFSIGNALED(status) ? "was killed by signal " + std::to_string(WTERMSIG(status))
                                  : "failed (status " + std::to_string(status) + ")";
            // Its last lines (progress bars redraw a line with carriage returns: each redraw is a line here).
            std::replace(text.begin(), text.end(), '\r', '\n');
            std::vector<std::string> lines;
            std::istringstream is(text);
            for (std::string line; std::getline(is, line);) {
                if (!line.empty()) lines.push_back(line);
            }
            constexpr size_t kShown = 20;
            log << "[qcmsa] " << name << ": " << command << std::endl;
            for (size_t i = lines.size() > kShown ? lines.size() - kShown : 0; i < lines.size(); i++) log << "[qcmsa]   " << lines[i] << '\n';
            log << std::flush;
            RunStatus::Get().Fail("qcmsa " + how + " for " + name +
                                  " (post-filter skipped; the raw MSA is still in " + msa + "; its output above)");
        } else if (!fs::exists(prefix + ".msa.fna")) {
            // Not an error (qcmsa may filter everything out), but no filtered MSA is not a success either.
            log << "[qcmsa] WARNING: " << name << ": qcmsa kept no gene or sample, so there is no " << prefix
                << ".msa.fna (the raw MSA is " << msa << "; " << command << " shows why)" << std::endl;
        }
    }

    // Removes a species' strain outputs of an earlier run, so that none of them survives a run that
    // writes it no more (e.g. an MSA qcmsa now filters away, or a species with fewer samples).
    static void RemoveStrainOutputs(Options const& options, std::string const& name) {
        std::string const prefix = options.GetStrainOutputDir() + '/' + name;
        for (auto const& stale : { options.GetMSAOutput(name), options.GetMSAPartitionOutput(name),
                                   options.GetMSAStatsOutput(name), options.GetSpeciesMetaOutput(name),
                                   options.GetHaplotypesOutput(name),
                                   prefix + ".msa.fna", prefix + ".partition.txt",
                                   prefix + ".qcmsa_summary.tsv", prefix + ".qc.png", prefix + ".qcmsa_output.partial" }) {
            std::error_code ec;
            std::filesystem::remove(stale, ec);
        }
    }

    // A long-read sample's strains of taxon `taxid` (Haplotypes.h): its reads' records phased at the multi-allelic sites
    // of the genes the MSA takes (`genes`, those its taxon has and does not drop), as MSAItem gives them (the taxon's
    // own reads, by --msa_identity_margin) with the SNP filters of its reads' kind and up to four alleles a site; reads
    // are cut between genes that are unlikely neighbours in the taxon's clade (the database's gene neighbours).
    static haplotypes::Phasing PhaseSample(uint32_t taxid, profiler::MicrobialProfile const& profile, std::vector<uint32_t> const& genes,
                                           Genome& genome, GenomeLoader const& loader, Options const& options) {
        auto const& taxon = profile.GetTaxa().at(taxid);
        double const min_identity = taxon.IdentityThreshold(options.GetMSAIdentityMargin());
        double const min_af = options.GetSNPMinAF(profile.GetReadType());
        auto const min_cov = options.GetSNPMinCov();
        auto const min_qual_sum = options.GetSNPMinPhredSum();
        auto const min_mean_qual = options.GetSNPMinMeanQual();
        bool const require_strand = options.GetSNPRequireStrand();
        std::vector<haplotypes::Site> sites;
        for (auto geneid : genes) {
            if (!taxon.GetGenes().contains(geneid) || taxon.DropsGene(geneid)) continue;
            auto const item = taxon.GetGenes().at(geneid).MSAItem(MSAItemFor(options, profile.GetReadType(), min_identity));
            auto const reference = genome.GetGene(geneid).Sequence();
            for (auto& [pos, bases] : MultiAllelicSites(item.first, item.second, min_cov, min_qual_sum, min_af, require_strand,
                                                        min_mean_qual, 4)) {
                if (pos >= reference.size()) continue;
                sites.push_back({ geneid, pos, reference[pos], std::move(bases) });
            }
        }
        return haplotypes::Phase(taxon.PhaseRecords(), sites, MaxDivergenceBin(min_identity),
                                 haplotypes::UnlikelyNeighbours(loader.GetGeneNeighbours(), taxid));
    }

    // <species>.haplotypes.tsv: per sample phased (PhaseSample), a line per block of its multi-allelic sites (the
    // sites its reads link): its genes, sites and reads, its haplotypes' reads (the most first), the haplotype of each
    // strain row (1-based) if it is phased, and how much more likely that join is than the next; and the sample's rows
    // with their shares of the reads.
    static void WriteHaplotypes(std::ostream& os, Profiles const& profiles, std::vector<size_t> const& profile_indices,
                                std::vector<haplotypes::Phasing> const& phasings, std::vector<bool> const& phased) {
        auto join = [](auto const& values, auto&& each) {
            std::ostringstream out;
            for (size_t i = 0; i < values.size(); i++) out << (i ? "," : "") << each(values[i]);
            return values.empty() ? std::string("-") : out.str();
        };
        os << "sample\trows\trow_shares\tblock\tgenes\tsites\treads\thaplotype_reads\trow_haplotypes\tlog_odds\tphased\n";
        for (size_t i = 0; i < profile_indices.size(); i++) {
            if (!phased[i]) continue;
            auto const& p = phasings[i];
            std::string const shares = join(p.shares, [](double s) { return s; });
            for (size_t b = 0; b < p.blocks.size(); b++) {
                auto const& block = p.blocks[b];
                os << profiles[profile_indices[i]].GetName() << '\t' << p.rows << '\t' << shares << '\t' << b + 1 << '\t'
                   << join(block.genes, [](uint32_t g) { return g; }) << '\t' << block.sites << '\t' << block.reads << '\t'
                   << join(block.haplotype_reads, [](size_t n) { return n; }) << '\t'
                   << join(block.row_haplotype, [](int h) { return h + 1; }) << '\t'
                   << block.log_odds << '\t' << (block.phased ? "yes" : "no") << '\n';
            }
        }
    }

    static void GetMSAForTaxon (uint32_t taxid, std::string taxon_name, GenomeLoader& loader, Options& options, Profiles& profiles, std::ostream* os_meta=nullptr, std::optional<profiler::TaxonFilterObj> const& filter={}, std::ostream& out = std::cout) {
        auto min_hcov = options.GetMSAMinHCOV();
        auto min_qual_sum = options.GetSNPMinPhredSum();
        auto min_cov = options.GetSNPMinCov();
        auto require_strand = options.GetSNPRequireStrand();
        auto min_mean_qual = options.GetSNPMinMeanQual();
        auto snp_max_alleles = options.GetSNPMaxAlleles();
        uint32_t const min_depth = static_cast<uint32_t>(options.GetMSAMinDepth());

        std::vector<size_t> profile_indices = GetProfilesWithTaxon(taxid, profiles, options, filter);

        if (profile_indices.empty()) return;

        auto& genome = loader.GetGenome(taxid);
        if (!genome.IsLoaded()) genome.LoadGenomeOMP();
        std::vector<uint32_t> selected_genes = SelectGenesForTaxon(taxid, taxon_name, profile_indices, loader, options, profiles, out);

        // A long-read sample whose reads show two or more strains gets a row per strain, <sample>_hap1, ... (the
        // most abundant first), each called from the reads of its strain (Haplotypes.h); every other sample one, the
        // consensus of its reads.
        std::vector<haplotypes::Phasing> phasings(profile_indices.size());
        std::vector<bool> phased(profile_indices.size(), false);  // phasing tried
        if (options.Phasing()) {
            for (size_t i = 0; i < profile_indices.size(); i++) {
                auto const& profile = profiles[profile_indices[i]];
                if (!IsLongReadType(profile.GetReadType()) || profile.GetTaxa().at(taxid).PhaseRecords().empty()) continue;
                phasings[i] = PhaseSample(taxid, profile, selected_genes, genome, loader, options);
                phased[i] = true;
            }
        }
        struct Row {
            size_t sample;  // index into profile_indices
            int haplotype;  // -1: the sample's own reads
        };
        std::vector<Row> rows;
        for (size_t i = 0; i < profile_indices.size(); i++) {
            if (phasings[i].rows < 2) rows.push_back({ i, -1 });
            else for (size_t h = 0; h < phasings[i].rows; h++) rows.push_back({ i, static_cast<int>(h) });
        }
        // A strain row's reads' records, by gene.
        std::vector<std::unordered_map<uint32_t, std::vector<haplotypes::ReadRecord const*>>> row_genes(rows.size());
        for (size_t r = 0; r < rows.size(); r++) {
            if (rows[r].haplotype < 0) continue;
            auto const& records = profiles[profile_indices[rows[r].sample]].GetTaxa().at(taxid).PhaseRecords();
            for (auto index : phasings[rows[r].sample].row_records[static_cast<size_t>(rows[r].haplotype)]) {
                row_genes[r][records[index].gene].push_back(&records[index]);
            }
        }

        // Each row takes the minimum allele frequency of its sample's reads' kind.
        std::vector<double> min_afs;
        std::vector<std::string> names;
        for (auto const& row : rows) {
            auto const& profile = profiles[profile_indices[row.sample]];
            min_afs.push_back(options.GetSNPMinAF(profile.GetReadType()));
            names.push_back(row.haplotype < 0 ? profile.GetName() : profile.GetName() + "_hap" + std::to_string(row.haplotype + 1));
        }

        MSAVector msa{ rows.size(), std::vector<char>() };
        protal::MSARow ref_msa_row;

        // Per-row accumulated SNP-retention statistics across all genes
        protal::MSAStats sample_stats(rows.size());

        std::vector<std::string> partitions;
        size_t partition_start = 0;
        size_t previous_size = 0;

        out << taxid << ": " << taxon_name << " across samples " << profile_indices.size() << std::endl;
        for (size_t i = 0; i < profile_indices.size(); i++) {
            auto const& p = phasings[i];
            if (!phased[i]) continue;
            out << "  " << profiles[profile_indices[i]].GetName() << ": " << p.reads << " long reads at " << p.blocks.size()
                      << " blocks of multi-allelic sites (" << p.cut_reads << " cut between unlikely neighbours); ";
            if (p.rows < 2) {
                out << "one row" << std::endl;
            } else {
                out << p.rows << " strain rows (shares";
                for (auto share : p.shares) out << ' ' << share;
                out << "), " << p.PhasedBlocks() << " of the blocks phased" << std::endl;
            }
        }
        if (std::find(phased.begin(), phased.end(), true) != phased.end()) {
            std::ofstream os_haplotypes(options.GetHaplotypesOutput(taxon_name));
            WriteHaplotypes(os_haplotypes, profiles, profile_indices, phasings, phased);
            os_haplotypes.close();
            if (os_haplotypes.fail()) RunStatus::Get().Fail("Writing the haplotypes of " + taxon_name + " failed: " + options.GetHaplotypesOutput(taxon_name));
        }

        // The numbers of a row's .meta.tsv line for a gene, from its item (the gene as the MSA takes it).
        struct Meta {
            double vertical_coverage = 0;
            size_t counts_vcov1 = 0, counts_vcov2 = 0, multi_allelic = 0, filtered = 0, gene_length = 0;
            double median_vcov = 0, hcov = 0, mean_vcov_nonzero = 0, median_vcov_nonzero = 0;
        };
        auto meta_of = [&](std::pair<VariantVec, CoverageVec> const& item, double vertical_coverage, size_t filtered,
                           size_t gene_length, double min_af) {
            auto const& tmp_vec = item.second;
            Meta meta;
            meta.vertical_coverage = vertical_coverage;
            meta.counts_vcov1 = std::count_if(tmp_vec.begin(), tmp_vec.end(), [](auto val){ return(val >= 1);});
            meta.counts_vcov2 = std::count_if(tmp_vec.begin(), tmp_vec.end(), [](auto val){ return(val >= 2);});
            // Multi-allelic positions as the MSA writes them (IUPAC codes), which qcmsa filters on.
            meta.multi_allelic = MultiAllelicPositions(item.first, tmp_vec, min_cov, min_qual_sum, min_af,
                                                       require_strand, min_mean_qual, snp_max_alleles);
            meta.filtered = filtered;
            meta.gene_length = gene_length;
            meta.hcov = static_cast<double>(meta.counts_vcov1) / static_cast<double>(gene_length);
            if (os_meta) {
                auto sorted_cov = tmp_vec;
                sorted_cov.resize(gene_length, 0);
                std::sort(sorted_cov.begin(), sorted_cov.end());
                size_t n = sorted_cov.size();
                meta.median_vcov = n % 2 == 1 ?
                    static_cast<double>(sorted_cov[n/2]) :
                    (static_cast<double>(sorted_cov[n/2 - 1]) + static_cast<double>(sorted_cov[n/2])) / 2.0;

                auto nonzero_begin = std::lower_bound(sorted_cov.begin(), sorted_cov.end(), 1);
                size_t nz = std::distance(nonzero_begin, sorted_cov.end());
                if (nz > 0) {
                    double sum = std::accumulate(nonzero_begin, sorted_cov.end(), 0.0);
                    meta.mean_vcov_nonzero = sum / nz;
                    size_t mid = nz / 2;
                    meta.median_vcov_nonzero = nz % 2 == 1 ?
                        static_cast<double>(*(nonzero_begin + mid)) :
                        (static_cast<double>(*(nonzero_begin + mid - 1)) + static_cast<double>(*(nonzero_begin + mid))) / 2.0;
                }
            }
            return meta;
        };

        for (auto& geneid : selected_genes) {
            auto& gene = genome.GetGene(geneid);
            auto const reference = gene.Sequence();

            // Each row's gene as the MSA takes it, and the numbers of its .meta.tsv line; none for a sample without
            // the gene, or whose taxon drops it (a foreign gene), and for a strain row whose reads lack it. A sample's
            // gene from its own reads; a strain row's from its strain's (HaplotypeItem).
            MSASequenceItems items(rows.size());
            std::vector<std::optional<Meta>> metas(rows.size());
            size_t rows_with_gene = 0;

            for (size_t r = 0; r < rows.size(); r++) {
                auto& profile = profiles[profile_indices[rows[r].sample]];
                auto& taxon = profile.GetTaxa().at(taxid);
                auto& genes = taxon.GetGenes();
                if (!genes.contains(geneid) || taxon.DropsGene(geneid)) continue;
                auto& gene_obs = genes.at(geneid);
                double const min_af = options.GetSNPMinAF(profile.GetReadType());
                if (rows[r].haplotype < 0) {
                    auto ac = gene_obs.AlleleSNPCounts(min_cov, min_qual_sum);
                    // The gene from the taxon's own reads, as the MSA takes it: the reads with a base per position,
                    // what the MSA judges each position by (kept so when the sample was profiled, Gene::KeepForMSAs).
                    auto item = gene_obs.MSAItem(MSAItemFor(options, profile.GetReadType(),
                                                            taxon.IdentityThreshold(options.GetMSAIdentityMargin())));
                    metas[r] = meta_of(item, gene_obs.VerticalCoverage(), ac.Filtered(), gene_obs.m_gene_length, min_af);
                    items[r] = std::move(item);
                } else {
                    auto const found = row_genes[r].find(geneid);
                    if (found == row_genes[r].end()) continue;
                    auto item = haplotypes::HaplotypeItem(found->second, reference, min_cov, min_af, min_mean_qual,
                                                          min_qual_sum, require_strand);
                    size_t filtered = 0;  // sites where no allele has min_cov reads and min_qual_sum (AlleleSNPCounts)
                    for (auto const& bin : item.first) {
                        filtered += std::none_of(bin.begin(), bin.end(), [&](Variant const& v) {
                            return v.Observations() >= min_cov && v.QualitySum() >= min_qual_sum;
                        });
                    }
                    double const bases = std::accumulate(item.second.begin(), item.second.end(), 0.0);
                    metas[r] = meta_of(item, bases / static_cast<double>(gene_obs.m_gene_length), filtered,
                                       gene_obs.m_gene_length, min_af);
                    items[r] = std::move(item);
                }
                // protal emits a raw MSA: every observed gene contributes its
                // sequence. Coverage-based gating (horizontal coverage, depth,
                // min samples per gene) is done by the qcmsa post-filter, which
                // reads the hcov / mean_vcov_nonzero columns written below.
                rows_with_gene++;
            }

            if (os_meta) {
                for (size_t r = 0; r < rows.size(); r++) {
                    if (!metas[r]) continue;
                    auto const& meta = *metas[r];
                    double const v1 = static_cast<double>(meta.counts_vcov1), v2 = static_cast<double>(meta.counts_vcov2);
#pragma omp critical(metaout)
                    {
                        *os_meta << names[r] << '\t';
                        *os_meta << geneid << '\t';
                        *os_meta << meta.vertical_coverage << '\t';
                        *os_meta << meta.counts_vcov1 << '\t';
                        *os_meta << meta.counts_vcov2 << '\t';
                        *os_meta << meta.multi_allelic << '\t';
                        *os_meta << meta.filtered << '\t';
                        *os_meta << (v1 > 0 ? meta.multi_allelic / v1 : 0) << '\t';
                        *os_meta << (v1 > 0 ? meta.filtered / v1 : 0) << '\t';
                        *os_meta << (v2 > 0 ? meta.multi_allelic / v2 : 0) << '\t';
                        *os_meta << (v2 > 0 ? meta.filtered / v2 : 0) << '\t';
                        *os_meta << meta.median_vcov << '\t';
                        *os_meta << meta.hcov << '\t';
                        *os_meta << meta.gene_length << '\t';
                        *os_meta << meta.mean_vcov_nonzero << '\t';
                        *os_meta << meta.median_vcov_nonzero;
                        *os_meta << std::endl;
                    }
                }
            }
            if (rows_with_gene > 0) {
                previous_size = msa.front().size();

                protal::MSAStats gene_stats(items.size());
                bool result = protal::MSA(items, reference, msa, min_cov, min_qual_sum, min_afs, require_strand, min_mean_qual,
                                          &gene_stats, &ref_msa_row, snp_max_alleles, min_depth);

                if (!result) continue;

                for (size_t si = 0; si < gene_stats.size(); si++) {
                    sample_stats[si] += gene_stats[si];
                }
                // Partitions are 1-based and inclusive, as RAxML and IQ-TREE read them.
                if (msa.front().size() > partition_start) {
                    if (partitions.size() > 0) {
                        partitions.back() += std::to_string(previous_size);
                    }

                    std::string partition = "DNA, gene";
                    partition += std::to_string(geneid) + " = ";
                    partition += std::to_string(partition_start + 1) + '-';
                    partitions.emplace_back(partition);
                    partition_start = msa.front().size();
                }
            }
        }

        if (partitions.empty()) {
            out << "No gene of " << taxon_name << " has enough coverage for an MSA" << std::endl;
            return;
        }

        // Prepend reference sequence as the first row so it is always present in both outputs.
        msa.insert(msa.begin(), std::move(ref_msa_row));
        names.insert(names.begin(), taxon_name + "_reference");
        sample_stats.insert(sample_stats.begin(), protal::MSASampleStats{});

        bool any_good = std::any_of(msa.begin(), msa.end(), [min_hcov](MSARow const& row){
            return IsRowGood(row, min_hcov);
        });
        if (!any_good) {
            out << "No good consensus sequences found for species" << std::endl;
            return;
        }

        std::ofstream os(options.GetMSAOutput(taxon_name), std::ios::out);
        for (auto i = 0; i < msa.size(); i++) {
            auto& row = msa[i];

            if (!IsRowGood(row, min_hcov)) continue;
            os << ">" << names[i] << std::endl;

            // os << std::string_view(&row[0], std::distance(row.begin(), row.end())) << std::endl;
            os << std::string(&row[0], std::distance(row.begin(), row.end())) << std::endl;
        }
        os.close();
        if (os.fail()) RunStatus::Get().Fail("Writing the MSA of " + taxon_name + " failed: " + options.GetMSAOutput(taxon_name));
        // std::cout << " Saved MSA to " << options.GetMSAOutput(taxon_name);

        partitions.back() += std::to_string(msa.front().size());

        std::ofstream os_part(options.GetMSAPartitionOutput(taxon_name), std::ios::out);
        for (auto i = 0; i < partitions.size(); i++) {
            os_part << partitions[i] << std::endl;
        }
        os_part.close();
        if (os_part.fail()) RunStatus::Get().Fail("Writing the MSA partitions of " + taxon_name + " failed: " + options.GetMSAPartitionOutput(taxon_name));

        // std::cout << " Saved Partitions " << std::endl;

        // Strain output is the raw MSA (.raw.msa.fna) plus qcmsa's final
        // .msa.fna (+ their partition files). The genecol_filtered and processed
        // MSA variants were removed: per-gene multi-allelicity and site-level
        // filtering are handled by the qcmsa post-filter on the raw MSA.

        // Write per-sample SNP-retention statistics TSV.
        {
            std::ofstream os_stats(options.GetMSAStatsOutput(taxon_name), std::ios::out);
            os_stats << "sample"
                     // --- totals ---
                     << "\ttotal_variant_positions"
                     << "\ttotal_pass_snps"
                     << "\ttotal_filtered_snps"
                     // --- per-type pass counts + % of total_variant_positions ---
                     << "\tsnps_retained\tsnps_retained_pct"
                     << "\tinsertions_retained\tinsertions_retained_pct"
                     << "\tdeletions_retained\tdeletions_retained_pct"
                     // --- filter breakdown + % of total_variant_positions ---
                     << "\tvariants_filtered_qual_sum\tvariants_filtered_qual_sum_pct"
                     << "\tvariants_filtered_obs_cov\tvariants_filtered_obs_cov_pct"
                     << "\tvariants_filtered_af\tvariants_filtered_af_pct"
                     << "\tvariants_filtered_strand\tvariants_filtered_strand_pct"
                     // --- pass/filter summary percentages ---
                     << "\ttotal_pass_pct\ttotal_filtered_pct"
                     // --- position-level counts + % of total_positions ---
                     << "\tpositions_ref\tpositions_ref_pct"
                     << "\tpositions_below_min_cov\tpositions_below_min_cov_pct"
                     << "\tpositions_no_coverage\tpositions_no_coverage_pct"
                     // --- appended last so readers that index columns by position keep working ---
                     << "\trefs_retained\trefs_retained_pct"
                     << '\n';

            os_stats << std::fixed << std::setprecision(2);
            for (size_t si = 0; si < sample_stats.size(); si++) {
                auto const& s = sample_stats[si];
                os_stats << names[si]
                         << '\t' << s.TotalVariantPositions()
                         << '\t' << s.TotalPass()
                         << '\t' << s.TotalFiltered()
                         << '\t' << s.snps_retained              << '\t' << s.PctSnpsRetained()
                         << '\t' << s.insertions_retained        << '\t' << s.PctInsertionsRetained()
                         << '\t' << s.deletions_retained         << '\t' << s.PctDeletionsRetained()
                         << '\t' << s.variants_filtered_qual_sum << '\t' << s.PctFilteredQualSum()
                         << '\t' << s.variants_filtered_obs_cov  << '\t' << s.PctFilteredObsCov()
                         << '\t' << s.variants_filtered_af       << '\t' << s.PctFilteredAF()
                         << '\t' << s.variants_filtered_strand   << '\t' << s.PctFilteredStrand()
                         << '\t' << s.PctPass()
                         << '\t' << s.PctFiltered()
                         << '\t' << s.positions_ref              << '\t' << s.PctPositionsRef()
                         << '\t' << s.positions_below_min_cov    << '\t' << s.PctPositionsBelowMinCov()
                         << '\t' << s.positions_no_coverage      << '\t' << s.PctPositionsNoCoverage()
                         << '\t' << s.refs_retained              << '\t' << s.PctRefsRetained()
                         << '\n';
            }
            os_stats.close();
        }
    }



    static SimilarityMatrix GetSimilarityMatrixForTaxon(uint32_t taxid, Options& options, Profiles& profiles, std::optional<profiler::TaxonFilterObj> filter={}) {
        SimilarityMatrix matrix;

        size_t min_shared_length = 1000;

        size_t total_combinations = profiles.size() * profiles.size() - profiles.size() / 2;
        size_t count_combinations = 1;

        ProgressBar prog(((profiles.size() * profiles.size()) - profiles.size()) / 2);
        for (auto i = 0; i < profiles.size(); i++) {
            auto& profile1 = profiles[i];
            auto name1 = profile1.GetName();

            for (auto j = i+1; j < profiles.size(); j++) {
                prog.UpdateAdd(1);
                count_combinations++;
//                std::cout << "\rCombination " << count_combinations << " of " << total_combinations;
                auto& profile2 = profiles[j];
                auto name2 = profile2.GetName();

                auto [similarity, shared_length] = GetSimilarity(taxid, profile1, profile2, options, filter);

//                std::cout << name1 << " -- " << name2 <<  " = " << similarity << " over " << shared_length << std::endl;

                if (similarity == NAN || shared_length < min_shared_length) continue;

                if (!matrix.HasName(name1)) {
                    matrix.AddName(name1);
                }
                if (!matrix.HasName(name2)) {
                    matrix.AddName(name2);
                }

//                std::cout << taxid << " " << similarity << " " << shared_length << std::endl;

                matrix.SetValue(name1, name2, similarity, true);
                matrix.SetValue(name1, name2, static_cast<double>(shared_length), false);
            }
        }
        std::cout << std::endl;
        return matrix;
    }



    static void WriteDistanceMatrix(uint32_t id, SimilarityMatrix const& matrix, Options& options, std::string& name) {
        std::ofstream matrix_os(options.GetSimilarityMatrixOutput(name), std::ios::out);
        matrix.PrintMatrix(matrix_os, "\t", 10);
        matrix_os.close();
    }

    void StrainWrapper2(Options& options, Profiles& profiles, GenomeLoader& loader, taxonomy::IntTaxonomy& taxonomy, std::vector<uint32_t> msa_taxids, std::optional<profiler::TaxonFilterObj> filter) {
        Benchmark bm_strain{"Strain-level MSAs"};
        bm_strain.Start();
        // Its two parts, over all species (scripts/measure_performance.sh reads them): building the MSAs and the qcMSA filter.
        Benchmark bm_msa{ "Building the strain MSAs" };
        Benchmark bm_qcmsa{ "qcMSA" };
        size_t species = 0, raw_msas = 0, filtered_msas = 0;
        // The folder the MSAs go to (made at the start of the run as well; -o may be the current folder, "").
        std::cout << "Strain outputs: " << options.GetStrainOutputDir() << std::endl;
        if (!MakeFolder(options.GetStrainOutputDir())) {
            std::cout << "Cannot create directories for this path " << options.GetStrainOutputDir() << std::endl;
            exit(2);
        }

        auto taxids = msa_taxids.empty() ? ExtractTaxa(profiles, filter) : msa_taxids;

        // The species of this run: strain outputs of other species in the directory are an earlier
        // run's (protal leaves them alone, as the directory may be shared).
        std::ostringstream species_list;
        species_list << "species\ttaxid\tsamples\traw_msa\tfiltered_msa\n";

        // The species are independent: their MSAs and qcMSA run over the threads, each species on one (at GTDB r226, 114
        // species took 29 s one after another, 23 s of it qcMSA; docs/claude/2026-10-04-performance-gtdb-scale). A
        // species' messages and qcMSA's output are held until it is done and printed in species order; the list of species
        // is written in that order. On one thread everything is printed as it comes, as before.
        size_t const n = taxids.size();
        size_t const workers = std::max<size_t>(1, std::min<size_t>(options.GetThreads(), n));
        bool const parallel = workers > 1;
        // First, on this thread, the samples each species enters: that scores its taxa with the model (the score is cached
        // in each taxon), so that the workers only read the profiles.
        std::vector<std::string> names(n);
        std::vector<size_t> samples(n);
        for (size_t i = 0; i < n; i++) {
            names[i] = taxonomy.Get(taxids[i]).scientific_name;
            std::replace(names[i].begin(), names[i].end(), ' ', '_');
            samples[i] = GetProfilesWithTaxon(taxids[i], profiles, options, filter).size();
        }
        // File names: the MSAs are next to the list.
        auto written = [](std::string const& path) {
            return std::filesystem::exists(path) ? std::filesystem::path(path).filename().string() : std::string("-");
        };
        std::vector<std::string> rows(n), outs(n), errs(n);
        std::vector<char> done(n, 0);
        std::mutex mutex;
        size_t next_print = 0;
        if (options.GetRunQCMSA() && n > 0) {
            std::string const script = FindQCMSAScript(options);
            std::cout << "qcMSA: " << (script.empty() ? std::string("not found") : script) << " on each species' raw MSA (its "
                      << "output is shown only if it fails)" << std::endl;
        }
        zstd::ParallelFor(n, static_cast<int>(workers), [&](size_t i, size_t) -> std::string {
            uint32_t const taxid = taxids[i];
            std::string const& name = names[i];
            std::ostringstream out_buffer, err_buffer;
            std::ostream& out = parallel ? static_cast<std::ostream&>(out_buffer) : std::cout;
            std::ostream& err = parallel ? static_cast<std::ostream&>(err_buffer) : std::cerr;
            out << taxonomy.Get(taxid).scientific_name << std::endl;
            RemoveStrainOutputs(options, name);
            if (samples[i] == 0) err << "[strains] " << name << " (--msa_species) passes in no sample, so it has no MSA" << std::endl;

            Benchmark msa{ "Building the strain MSAs" }, qc{ "qcMSA" };
            msa.Start();
            std::ofstream os_meta(options.GetSpeciesMetaOutput(name));
            os_meta << "sample\tgene_id\tvertical_coverage\tcounts_vcov1\tcounts_vcov2\tmulti_allelic\tfiltered\tmulti_rate_vcov1\tfiltered_rate_vcov1\tmulti_rate_vcov2\tfiltered_rate_vcov2\tmedian_vcov\thcov\tgene_length\tmean_vcov_nonzero\tmedian_vcov_nonzero\n";
            // With --scratch_dir the samples' evidence of this species is read back from their spill files for its MSA
            // only (the other species' workers read theirs: each species is its own taxa in the profiles).
            auto const in_msa = GetProfilesWithTaxon(taxid, profiles, options, filter);
            bool loaded = true;
            for (auto s : in_msa) {
                if (profiles[s].LoadStrainEvidence(taxid)) continue;
                RunStatus::Get().Fail("Cannot read the strain evidence of " + name + " in sample " + profiles[s].GetName() +
                                      " back from " + profiles[s].SpillPath() + " (--scratch_dir)");
                loaded = false;
            }
            if (loaded) GetMSAForTaxon(taxid, name, loader, options, profiles, &os_meta, filter, out);
            for (auto s : in_msa) profiles[s].DropStrainEvidence(taxid);
            os_meta.close();
            if (os_meta.fail()) RunStatus::Get().Fail("Writing the MSA metadata of " + name + " failed: " + options.GetSpeciesMetaOutput(name));
            msa.Stop();
            if (options.GetRunQCMSA()) {
                qc.Start();
                RunQCMSA(options, name, err, parallel);
                qc.Stop();
            }
            std::string const raw_msa = options.GetMSAOutput(name);
            std::string const filtered_msa = options.GetStrainOutputDir() + '/' + name + ".msa.fna";
            std::ostringstream row;
            row << name << '\t' << taxid << '\t' << samples[i] << '\t' << written(raw_msa) << '\t' << written(filtered_msa) << '\n';
            bool const has_raw = std::filesystem::exists(raw_msa), has_filtered = std::filesystem::exists(filtered_msa);

            std::lock_guard<std::mutex> lock(mutex);
            bm_msa.Join(msa, false);  // the sums over the species (not their wall time: they overlap)
            bm_msa.AddObservation();
            if (options.GetRunQCMSA()) {
                bm_qcmsa.Join(qc, false);
                bm_qcmsa.AddObservation();
            }
            species++;
            raw_msas += has_raw;
            filtered_msas += has_filtered;
            rows[i] = row.str();
            outs[i] = out_buffer.str();
            errs[i] = err_buffer.str();
            done[i] = 1;
            for (; next_print < n && done[next_print]; next_print++) {
                std::cout << outs[next_print] << std::flush;
                std::cerr << errs[next_print] << std::flush;
                std::string().swap(outs[next_print]);
                std::string().swap(errs[next_print]);
            }
            return {};
        });
        for (auto const& row : rows) species_list << row;

        std::ofstream os_list(options.GetStrainSpeciesListOutput(), std::ios::out);
        os_list << species_list.str();
        os_list.close();
        if (os_list.fail()) RunStatus::Get().Fail("Writing the list of strain species failed: " + options.GetStrainSpeciesListOutput());
        bm_strain.Stop();
        bm_strain.PrintResults();
        std::cout << "Strain-level MSAs of " << species << " species: " << raw_msas << " raw MSAs, " << filtered_msas << " filtered by qcMSA"
                  << (options.GetRunQCMSA() ? "" : " (--no_qcmsa)") << "; " << workers << " species at a time, the times below summed over "
                  << "the species" << std::endl;
        bm_msa.PrintResults();
        if (options.GetRunQCMSA()) bm_qcmsa.PrintResults();
    }
}
