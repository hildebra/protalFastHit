// FeatureVersions.h - the version of each presence-model feature's definition in this protal (since 2026-10-15). A model
// reaches its features by name (profiler::TaxonFeatures), so a model trained before a feature's definition changed would
// be given the new values under the old name without notice. Each feature here has a version (1: as first defined, one
// more for each change of what it computes from the same reads and database), the date of its last change and a note.
// `protal --feature_versions` prints the table as a TSV (WriteTable); the trainer stores it in each model it writes:
//   <Extension name="protal_feature_versions" value="name:version,name:version,..."/>
// and a run that loads a model checks the model's features against this table (Mismatch, RunProtal's CheckModelFeatures):
// one warning per model naming the features whose definition changed since its training; a model without the extension
// gets one note. A test checks that the table holds every feature TaxonFeatures writes and nothing else
// (test_RunSettings.cpp). When a feature's definition changes, raise its version here and give the date and the change.
#pragma once

#include <charconv>
#include <map>
#include <ostream>
#include <string>
#include <string_view>
#include <vector>

namespace protal::feature_versions {
    inline constexpr std::string_view kExtension = "protal_feature_versions";

    struct FeatureVersion {
        std::string_view name;
        int version;               // 1: as first defined
        std::string_view changed;  // the date of its definition (YYYY-MM-DD): its introduction, or its last change
        std::string_view note;
    };

    // Every feature of profiler::TaxonFeatures, in its order. The versions start on 2026-10-15: 1 for a feature computed as
    // when it was introduced, 2 for one whose computation from the same reads and database changed since (however often;
    // the note names the commits). Not counted: changes that act only through the database's tables or index (a model is
    // trained on its own database), opt-in options, and run settings build_settings.tsv records since 2026-10-15 (the
    // depth margin, foreign genes, the adaptive candidates). From here on, each change raises the version by one.
    inline constexpr FeatureVersion kFeatureVersions[] = {
        { "present_genes", 2, "2024-02-01", "genes with records; records on genes the database does not mark hittable left out since 1375534" },
        { "total_hits", 2, "2024-02-01", "records on the taxon; non-hittable genes' left out since 1375534" },
        { "unique_hits", 2, "2024-02-01", "records with MAPQ above 20 (not ZU); non-hittable genes' left out since 1375534" },
        { "mean_ani", 2, "2024-02-01", "mean identity of the records; non-hittable genes' left out since 1375534" },
        { "expected_gene_presence", 2, "2024-02-01", "genes expected hit by its records; hittable genes only since 1375534" },
        { "expected_gene_presence_ratio", 1, "2024-05-27", "present genes over expected_gene_presence" },
        { "expected_gene_presence_unique_weighted", 2, "2026-10-14", "genes expected hit by its ZU draws; ZU with short uniques (11cdf6b), from unique_cores.bin (3ef6ad6)" },
        { "expected_gene_presence_ratio_unique_weighted", 2, "2026-10-14", "genes with ZU over their expectation; ZU with short uniques (11cdf6b), from unique_cores.bin (3ef6ad6)" },
        { "uniqueness", 1, "2024-05-27", "share of records with MAPQ above 20" },
        { "mean_mapq", 2, "2024-02-01", "mean MAPQ of the records; non-hittable genes' left out since 1375534" },
        { "variance1", 2, "2026-10-02", "variance of the genes' depth; fragment bases, a pair's overlap once (4a0a8e8)" },
        { "variance2", 2, "2026-10-02", "variance of the depth of genes with 5+ reads; fragment bases (4a0a8e8)" },
        { "A0", 2, "2026-09-30", "variant positions with 0 alleles; allele rules ca2df11, --no_strains 050e6cd, once per fragment 692db41" },
        { "A1", 2, "2026-09-30", "variant positions with 1 allele; ca2df11, 050e6cd, 692db41" },
        { "A2", 2, "2026-09-30", "variant positions with 2 alleles; ca2df11, 050e6cd, 692db41" },
        { "A3", 2, "2026-09-30", "variant positions with 3 alleles; ca2df11, 050e6cd, 692db41" },
        { "A4", 2, "2026-09-30", "variant positions with 4 alleles; ca2df11, 050e6cd, 692db41" },
        { "AF0", 2, "2026-09-30", "A0 with 2+ observations and quality sum 60+; Phred+33 2da2748, ca2df11, 050e6cd, 692db41" },
        { "AF1", 2, "2026-09-30", "A1 filtered as AF0; 2da2748, ca2df11, 050e6cd, 692db41" },
        { "AF2", 2, "2026-09-30", "A2 filtered as AF0; 2da2748, ca2df11, 050e6cd, 692db41" },
        { "AF3", 2, "2026-09-30", "A3 filtered as AF0; 2da2748, ca2df11, 050e6cd, 692db41" },
        { "AF4", 2, "2026-09-30", "A4 filtered as AF0; 2da2748, ca2df11, 050e6cd, 692db41" },
        { "RAF0", 2, "2026-09-30", "AF0 over all AF; 2da2748, ca2df11, the dump's class-1 bug 050e6cd, 692db41" },
        { "RAF1", 2, "2026-09-30", "AF1 over all AF; 2da2748, ca2df11, 050e6cd, 692db41" },
        { "RAF2", 2, "2026-09-30", "AF2 over all AF; 2da2748, ca2df11, 050e6cd, 692db41" },
        { "RAF3", 2, "2026-09-30", "AF3 over all AF; 2da2748, ca2df11, 050e6cd, 692db41" },
        { "RAF4", 2, "2026-09-30", "AF4 over all AF; 2da2748, ca2df11, 050e6cd, 692db41" },
        { "RA0", 2, "2026-09-30", "A0 over all A; ca2df11, the dump's class-1 bug 050e6cd, 692db41" },
        { "RA1", 2, "2026-09-30", "A1 over all A; zero test eabf090, ca2df11, 050e6cd, 692db41" },
        { "RA2", 2, "2026-09-30", "A2 over all A; eabf090, ca2df11, 050e6cd, 692db41" },
        { "RA3", 2, "2026-09-30", "A3 over all A; eabf090, ca2df11, 050e6cd, 692db41" },
        { "RA4", 2, "2026-09-30", "A4 over all A; eabf090, ca2df11, 050e6cd, 692db41" },
        { "stddev", 2, "2026-10-02", "SD of the vertical coverage; fragment bases (4a0a8e8)" },
        { "hittable", 2, "2026-09-27", "the reference's hittable genes; a genome without unique k-mers has none since 0d04d2c" },
        { "lu", 2, "2026-10-14", "sum of ZU; short uniques (11cdf6b), from the cores the alignment covers in unique_cores.bin (3ef6ad6)" },
        { "lu_genes", 2, "2026-10-14", "genes with ZU above 0; 11cdf6b, 3ef6ad6" },
        { "lsu", 2, "2026-10-14", "sum of ZT; short uniques (11cdf6b), from unique_cores.bin (3ef6ad6)" },
        { "lsu_genes", 2, "2026-10-14", "genes with ZT above 0; 11cdf6b, 3ef6ad6" },
        { "su_genome", 1, "2024-05-27", "the reference's unique k-mers of a single-value core (unique_kmers.tsv)" },
        { "lu_genome", 2, "2026-10-06", "the reference's unique k-mers; plus the short ones since 11cdf6b" },
        { "lsu_genome", 2, "2026-10-06", "the reference's unique k-mers at distance two; plus the short ones since 11cdf6b" },
        { "total_genome", 1, "2024-05-27", "the reference's k-mers" },
        { "su_rate_ref", 1, "2024-05-27", "su_genome over total_genome" },
        { "lu_rate_ref", 2, "2026-10-06", "lu_genome over total_genome; with the short uniques (11cdf6b)" },
        { "lsu_rate_ref", 2, "2026-10-06", "lsu_genome over total_genome; with the short uniques (11cdf6b)" },
        { "lu_gene_rate", 2, "2026-10-14", "genes with ZU over genes with unique k-mers; e5335fd, 11cdf6b, unique_cores.bin 3ef6ad6" },
        { "lsu_gene_rate", 2, "2026-10-14", "genes with ZT over genes with such k-mers; e5335fd, 11cdf6b, 3ef6ad6" },
        { "lu_gene_rate2", 2, "2026-10-14", "lu_gene_rate of genes with 1+ such k-mers; e5335fd, 11cdf6b, 3ef6ad6" },
        { "lsu_gene_rate2", 2, "2026-10-14", "lsu_gene_rate of genes with 1+ such k-mers; e5335fd, 11cdf6b, 3ef6ad6" },
        { "lu_gene_rate3", 2, "2026-10-14", "lu_gene_rate of genes with 5+ such k-mers; e5335fd, 11cdf6b, 3ef6ad6" },
        { "lsu_gene_rate3", 2, "2026-10-14", "lsu_gene_rate of genes with 5+ such k-mers; e5335fd, 11cdf6b, 3ef6ad6" },
        { "lsu_per_read", 2, "2026-10-14", "lsu over total_hits; 11cdf6b, 3ef6ad6" },
        { "lu_per_read", 2, "2026-10-14", "lu over total_hits; 11cdf6b, 3ef6ad6" },
        { "fragments", 1, "2026-09-27", "fragments the MAPQ and length filters keep (--fragment_count soft: plus the dropped by EM share)" },
        { "depth", 2, "2026-10-02", "vertical coverage within the depth margin; fragment bases once (4a0a8e8), foreign genes left out (c5475ff)" },
        { "hit_gene_fraction", 1, "2026-09-27", "share of the reference's genes hit" },
        { "gene_presence_ratio", 1, "2026-09-27", "genes hit over those expected at the depth" },
        { "gene_dispersion", 1, "2026-09-27", "chi-square per degree of freedom of the fragments per hittable gene against its length share" },
        { "depth_cv", 2, "2026-10-02", "CV of the genes' depth; fragment bases (4a0a8e8)" },
        { "identity", 1, "2026-09-27", "base identity of the records, indels as differences" },
        { "top_identity", 2, "2026-10-02", "98th percentile identity; weighted by fragment bases (4a0a8e8)" },
        { "low_identity_share", 2, "2026-10-02", "share of aligned bases below the depth margin; fragment bases (4a0a8e8)" },
        { "lu_per_kb", 2, "2026-10-14", "lu per aligned kb; short uniques (11cdf6b), unique_cores.bin (3ef6ad6)" },
        { "lsu_per_kb", 2, "2026-10-14", "lsu per aligned kb; short uniques (11cdf6b), unique_cores.bin (3ef6ad6)" },
        { "variant_sites_per_kb", 2, "2026-09-30", "variant sites per covered kb; variants once per fragment (692db41)" },
        { "multiallelic_sites_per_kb", 2, "2026-09-30", "sites of 2+ alleles per covered kb; 692db41" },
        { "low_mapq_share", 1, "2026-10-01", "share of records of low MAPQ" },
        { "congener_fit_share", 1, "2026-10-01", "share of records a congener fits as well (ZA)" },
        { "other_genus_fit_share", 1, "2026-10-01", "share of records another genus fits as well (ZA)" },
        { "linked_share", 1, "2026-10-01", "share of reads with 2+ records on the taxon" },
        { "adjacent_expected_share", 2, "2026-10-02", "share of neighbouring genes on reads that the clade expects; smoothing and thresholds d11381f" },
        { "adjacent_unlikely_share", 2, "2026-10-02", "share of neighbouring genes on reads unlikely in the clade; d11381f" },
        { "adjacent_support", 1, "2026-10-02", "mean frequency of the reads' gene pairings in the clade" },
        { "excess_median", 1, "2026-10-01", "median divergence beyond the base qualities' errors" },
        { "excess_high_share", 1, "2026-10-01", "share of reads whose excess is above 0.02" },
        { "conserved_fast_depth_ratio", 2, "2026-10-02", "depth of conserved over fast genes; fragment bases (4a0a8e8)" },
        { "conserved_hit_share", 1, "2026-10-01", "share of hit genes that are conserved" },
        { "conserved_fast_record_ratio", 1, "2026-10-01", "conserved over fast genes' depth of every best record before the filters" },
        { "conserved_fast_kept_ratio", 2, "2026-10-02", "the same of the filtered records; fragment bases (4a0a8e8)" },
        { "excess_scaled_median", 1, "2026-10-03", "median excess over the genes' conservation factors" },
        { "excess_conserved_fast_ratio", 1, "2026-10-03", "excess on conserved over fast genes" },
        { "third_position_share", 1, "2026-10-03", "share of mismatches at codon third positions" },
        { "mate_lost_share", 1, "2026-10-03", "share of paired fragments whose expected mate has no record on it" },
        { "genus_skew", 1, "2026-10-03", "log ratio of its fragments to the genus's most abundant species'" },
        { "family_skew", 1, "2026-10-03", "the same against the family's" },
        { "genus_share", 1, "2026-10-03", "its share of the genus's fragments" },
        { "genus_spill", 1, "2026-10-03", "fragments the genus's top species would spill onto it, over its own" },
        { "relative_skew", 1, "2026-10-03", "log ratio to its nearest abundant relative by reference distance" },
        { "relative_distance", 1, "2026-10-03", "distance to that relative's reference" },
        { "relative_spill", 1, "2026-10-03", "fragments that relative would spill onto it, over its own" },
        { "relative_close_share", 1, "2026-10-03", "share of its fragments from close relatives" },
        { "em_own_share", 2, "2026-10-05", "share of its best records the read EM leaves it; stop rule 0f4e664, 100 sweeps f2bafd9" },
        { "em_kept_own_share", 2, "2026-10-05", "the same of the filtered records; 0f4e664, f2bafd9" },
        { "sample_log_fragments", 1, "2026-10-03", "log10 of the sample's fragments over all taxa" },
        { "sample_log_taxa", 1, "2026-10-07", "log10 of the sample's taxa with fragments" },
        { "sample_low_identity", 1, "2026-10-07", "the sample's fragments' share on low-identity bases" },
        { "sample_identity", 1, "2026-10-07", "the sample's taxa's fragment-weighted median identity" },
        { "sample_thin_share", 1, "2026-10-16", "the sample's taxa of 1-2 fragments over its taxa with fragments" },
        { "sample_thin_close_share", 1, "2026-10-16", "of those, the share within 0.01 of the sample's identity" },
        { "fragments_all", 1, "2026-10-03", "fragments before the filters" },
        { "em_fragments", 2, "2026-10-05", "em_own_share times fragments_all; EM 100 sweeps f2bafd9" },
        { "failed_candidate_rate", 2, "2026-10-14", "reads seeded on it that did not align, over those tried; single-seed anchors out (ce2940e, 3ef6ad6)" },
        { "rep_duplicate_share", 1, "2026-10-03", "share of the representative's single-copy markers found twice (species_priors.tsv)" },
        { "rep_completeness", 1, "2026-10-03", "the representative's CheckM completeness" },
        { "rep_contamination", 1, "2026-10-03", "the representative's CheckM contamination" },
        { "cluster_ani_radius", 1, "2026-10-03", "the species cluster's ANI radius" },
        { "cluster_mean_ani", 1, "2026-10-03", "the cluster's mean intra-species ANI" },
        { "cluster_min_ani", 1, "2026-10-03", "the cluster's minimum intra-species ANI" },
        { "cluster_genomes_log10", 1, "2026-10-03", "log10 of the cluster's genomes" },
        { "genus_top_fragments", 1, "2026-10-03", "fragments of the genus's most abundant other species" },
        { "read_consensus_share", 1, "2026-10-06", "share of long reads' genes their consensus taxon gave (ZR)" },
        { "read_inconsistent_share", 1, "2026-10-06", "share of genes clearly another taxon's (ZR)" },
        { "seed_crowding", 1, "2026-10-06", "mean log2 of the taxa its reads seeded on nearly as well as on the best (ZN)" },
        { "unexpected_congener_fit_share", 1, "2026-10-06", "share of reads a congener fits better than the references' distance allows" },
        { "split_fragment_share", 1, "2026-10-06", "share of its fragments with a best record on a congener too" },
        { "congener_gene_overlap", 1, "2026-10-06", "how its genes' reads and its congener's complement each other" },
        { "gene_divergence_dispersion", 1, "2026-10-06", "dispersion of its genes' divergence" },
        { "breadth_ratio", 1, "2026-10-06", "bases its hit genes' reads cover over the Lander-Waterman expectation at their depth" },
        { "failed_gene_share", 2, "2026-10-14", "share of genes with more failed reads than aligned; single-seed anchors out (ce2940e, 3ef6ad6)" },
        { "fixed_difference_rate", 1, "2026-10-06", "sites of a fixed non-reference allele per covered site, by conservation factor" },
        { "polymorphic_site_rate", 1, "2026-10-06", "sites with a second allele per covered site" },
        { "db_congeners_01", 1, "2026-10-06", "database congeners within 0.01 (species_neighbours.tsv)" },
        { "db_congeners_02", 1, "2026-10-06", "database congeners within 0.02" },
        { "db_congeners_05", 1, "2026-10-06", "database congeners within 0.05" },
        { "db_nearest_congener", 1, "2026-10-06", "distance to the nearest database congener" },
        { "ancestry_sites_per_record", 2, "2026-10-10", "ancestry sites per record; consensus sites 45bc47b 551b9ec 4b13640 f566cd0, allele-shared 371ac84" },
        { "ancestry_agreement", 2, "2026-10-10", "share of sites where the read has the reference's base; 45bc47b 551b9ec 4b13640 f566cd0 371ac84" },
        { "ancestry_congener_share", 2, "2026-10-10", "share of sites with the congeners' base; 45bc47b 551b9ec 4b13640 f566cd0 371ac84" },
        { "ancestry_indel_sites_per_record", 2, "2026-10-10", "indel ancestry sites per record; f566cd0, 371ac84" },
        { "ancestry_indel_congener_share", 2, "2026-10-10", "share of indel sites taken the congeners' way; f566cd0, 371ac84" },
        { "gap_informative_share", 1, "2026-10-07", "share of reads in the gaps to its congeners' copies (congener_gaps.tsv)" },
        { "gap_within_min_share", 1, "2026-10-07", "share of reads within the smallest gap" },
        { "gap_within_median_share", 1, "2026-10-07", "share of reads within the median gap" },
        { "gap_position", 1, "2026-10-07", "where its reads lie between it and its congeners" },
        { "foreign_scanned_share", 2, "2026-10-08", "share of scan reads from other species on its copies (foreign_rates.tsv); scan and share 551b9ec" },
        { "foreign_copy_share", 2, "2026-10-08", "share of its copies other species reach; 551b9ec" },
        { "foreign_genus_copy_share", 2, "2026-10-08", "share of its copies other genera reach; 551b9ec" },
        { "untried_candidate_rate", 1, "2026-10-07", "reads whose seeds fit it but never tried it (ZC), over those plus its reads with a record" },
        { "allele_copy_share", 1, "2026-10-08", "share of its copies with strain alleles (strain_alleles.tsv)" },
        { "allele_explained_share", 1, "2026-10-08", "share of the reads' differences its alleles explain" },
        { "allele_identity_gain", 1, "2026-10-08", "identity gained against the best allele" },
        { "polymorphic_known_share", 1, "2026-10-09", "share of polymorphic sites with a known allele's base" },
        { "polymorphic_novel_share", 1, "2026-10-09", "share of polymorphic sites with a base no allele has" },
        { "ancestry_fixed_gain", 2, "2026-10-10", "gain at the ancestry sites fixed within the species; consensus sites f566cd0" },
        { "ancestry_fixed_agreement", 1, "2026-10-10", "agreement at the fixed ancestry sites" },
        { "allele_sites_per_kb", 1, "2026-10-09", "allele sites per covered kb" },
        { "ancestry_fixed_share", 2, "2026-10-10", "share of ancestry sites fixed within the species; f566cd0" },
        { "conserved_mismatch_ratio", 1, "2026-10-10", "mismatches at conserved columns against the expectation (column_weights.tsv)" },
        { "conserved_mismatch_rate", 1, "2026-10-10", "mismatches at conserved columns per kb (0 with weights built before 0d961e1)" },
        { "ancestry_agreement_weighted", 2, "2026-10-10", "ancestry agreement weighted by the columns' conservation; allele-shared sites 371ac84" },
        { "nonsynonymous_share", 1, "2026-10-10", "share of mismatches that change the amino acid" },
        { "nonsynonymous_conserved_rate", 1, "2026-10-10", "nonsynonymous mismatches at conserved amino acids per kb (0 with weights before 0d961e1)" },
        { "column_weight_coverage", 1, "2026-10-10", "share of its family's columns with weights" },
        { "allele_identity", 1, "2026-10-13", "identity against each read's best strain allele" },
        { "allele_top_identity", 1, "2026-10-13", "top_identity against the best strain allele" },
        { "unique_seed_share", 2, "2026-10-14", "ZU over ZS; from the cores the alignment covers in unique_cores.bin (3ef6ad6)" },
        { "unique2_seed_share", 2, "2026-10-14", "ZT over ZS; unique_cores.bin (3ef6ad6)" },
        { "unique_seed_enrichment", 2, "2026-10-14", "unique_seed_share over the genes' unique share; unique_cores.bin, no seed expectation (3ef6ad6)" },
        { "unique2_seed_enrichment", 2, "2026-10-14", "the same at distance two; unique_cores.bin (3ef6ad6)" },
        { "ancestry_sites_per_kb", 1, "2026-10-14", "ancestry sites per aligned kb" },
        { "ancestry_indel_sites_per_kb", 1, "2026-10-14", "indel ancestry sites per aligned kb" },
        { "top_identity_bound", 1, "2026-10-14", "the best reads' identity bound given their aligned bases" },
        { "allele_top_identity_bound", 1, "2026-10-14", "the same against the best allele" },
        { "sample_read_length", 1, "2026-10-14", "the sample's mean read length" },
        { "sample_read_length_sd", 1, "2026-10-14", "the SD of the sample's read lengths" },
        { "aligned_length_ratio", 1, "2026-10-14", "its reads' aligned length over the sample's read length" },
    };

    inline FeatureVersion const* Find(std::string_view name) {
        for (auto const& feature : kFeatureVersions) {
            if (feature.name == name) return &feature;
        }
        return nullptr;
    }

    // The table as protal --feature_versions prints it: feature, version, changed, note.
    inline void WriteTable(std::ostream& os) {
        os << "feature\tversion\tchanged\tnote\n";
        for (auto const& f : kFeatureVersions) os << f.name << '\t' << f.version << '\t' << f.changed << '\t' << f.note << '\n';
    }

    // The extension's value as the trainer writes it from the table: "name:version,name:version,...".
    inline std::string ExtensionValue() {
        std::string value;
        for (auto const& f : kFeatureVersions) value += (value.empty() ? "" : ",") + std::string(f.name) + ":" + std::to_string(f.version);
        return value;
    }

    using Recorded = std::map<std::string, int, std::less<>>;

    // The versions of a model's extension into `recorded` (empty first); an error text if it is malformed: an item that
    // is not name:version with a version of 1 or more, or a feature twice.
    inline std::string Parse(std::string_view text, Recorded& recorded) {
        recorded.clear();
        while (!text.empty()) {
            size_t const comma = text.find(',');
            std::string_view const item = text.substr(0, comma);
            text = comma == std::string_view::npos ? std::string_view() : text.substr(comma + 1);
            size_t const colon = item.find(':');
            int version = 0;
            bool ok = colon != std::string_view::npos && colon > 0 && !(comma != std::string_view::npos && text.empty());
            if (ok) {
                std::string_view const number = item.substr(colon + 1);
                auto const [end, error] = std::from_chars(number.data(), number.data() + number.size(), version);
                ok = !number.empty() && error == std::errc() && end == number.data() + number.size() && version >= 1 &&
                     !recorded.contains(item.substr(0, colon));
            }
            if (!ok) {
                recorded.clear();
                return "its features' versions (" + std::string(kExtension) + ") are malformed ('" + std::string(item) +
                       "': name:version, the version 1 or more, each feature once)";
            }
            recorded.emplace(std::string(item.substr(0, colon)), version);
        }
        return {};
    }

    // Of the model's input features `features`, those whose version `recorded` (the model's) differs from this protal's
    // (or that this protal does not know), and those it records none of.
    struct Comparison {
        std::vector<std::string> differ;      // "name" of each, in the order of `features`
        std::vector<int> model_versions;      // and its version in the model
        std::vector<int> protal_versions;     // and in this protal (0: unknown to it)
        std::vector<std::string> unrecorded;  // features the model records no version of
    };

    inline Comparison Compare(Recorded const& recorded, std::vector<std::string> const& features) {
        Comparison c;
        for (auto const& name : features) {
            auto const it = recorded.find(name);
            if (it == recorded.end()) {
                c.unrecorded.push_back(name);
                continue;
            }
            auto const* own = Find(name);
            int const version = own ? own->version : 0;
            if (version == it->second) continue;
            c.differ.push_back(name);
            c.model_versions.push_back(it->second);
            c.protal_versions.push_back(version);
        }
        return c;
    }

    // The warning of model `model` (its name) for Compare's result, or "" if its features are this protal's: "the model X
    // was trained on other definitions of features a, b (versions in the model 1, 1; in this protal 2, 2): retrain it, or
    // expect calls that differ from its training", and the features it records no version of.
    inline std::string Mismatch(std::string const& model, Recorded const& recorded, std::vector<std::string> const& features) {
        auto const c = Compare(recorded, features);
        if (c.differ.empty() && c.unrecorded.empty()) return {};
        auto list = [](std::vector<std::string> const& items) {
            std::string text;
            for (auto const& item : items) text += (text.empty() ? "" : ", ") + item;
            return text;
        };
        auto versions = [](std::vector<int> const& items) {
            std::string text;
            for (int item : items) text += (text.empty() ? "" : ", ") + (item ? std::to_string(item) : std::string("none"));
            return text;
        };
        std::string text = "the model " + model;
        if (!c.differ.empty()) {
            text += " was trained on other definitions of feature(s) " + list(c.differ) + " (versions in the model " +
                    versions(c.model_versions) + "; in this protal " + versions(c.protal_versions) + ", protal --feature_versions)";
            if (!c.unrecorded.empty()) text += ", and";
        }
        if (!c.unrecorded.empty()) text += " records no version of feature(s) " + list(c.unrecorded);
        return text + ": retrain it, or expect calls that differ from its training";
    }
}
