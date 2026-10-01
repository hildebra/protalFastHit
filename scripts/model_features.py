"""Which columns of protal's training dump (<profile>.truth_annotated) a model may use.

The dump holds every feature protal computes for a taxon (TaxonFeatures in
src/Profiling/Profiler.h). The model shipped with protal uses absolute counts of genes, k-mers and
reads, which differ between databases (archaea have 52 marker genes, bacteria 119), depths and
read lengths. NORMALIZED_FEATURES are those that do not: fractions, ratios and rates per aligned
kb, counting fragments rather than mates. `fragments` and `depth` say how much evidence there is.
The dump also has gene_dispersion, which is left out here: on simulated data it separated present
from absent taxa hardly at all (AUC 0.43) and grew with depth.

The last four are about the reads' other candidates (docs/claude/2026-10-01-alignment-features): the reads' MAPQ
(mean, and the share below 10) and the shares of reads that another species of the genus or a species of another
genus fits as well (within an edit, from protal's ZA tag). A dump of an older protal lacks them. The dump's
linked_share (reads with two records on the taxon: both mates, or two genes of a long read) is left out: it lowered
the paired-end F1 on the tuning world.

Four more (docs/claude/2026-10-01-f1-opportunities): how far the reads differ from the reference beyond
what their base qualities explain (excess_median, excess_high_share; 0 without qualities), and the depth of the
taxon's conserved hit genes against its fast ones by the database's gene_conservation.tsv (conserved_fast_depth_ratio,
conserved_hit_share; 0 and 0.5 without it). A relative's reads exceed their errors and land on the fast genes.
"""

# Columns of the dump that describe the row, not the taxon's evidence; columns that
# collect_training_data.py adds start with "meta_".
NON_FEATURE_COLUMNS = {"truth", "truth_raw", "prediction", "probability", "taxon", "taxon_name", "dataset", "total_hits"}

NORMALIZED_FEATURES = [
    "fragments",
    "depth",
    "hit_gene_fraction",
    "gene_presence_ratio",
    "depth_cv",
    "identity",
    "top_identity",
    "low_identity_share",
    "uniqueness",
    "lu_per_kb",
    "lsu_per_kb",
    "variant_sites_per_kb",
    "multiallelic_sites_per_kb",
    "RAF0", "RAF1", "RAF2", "RAF3", "RAF4",
    "lu_gene_rate", "lsu_gene_rate",
    "lu_gene_rate2", "lsu_gene_rate2",
    "lu_gene_rate3", "lsu_gene_rate3",
    "mean_mapq", "low_mapq_share",
    "congener_fit_share", "other_genus_fit_share",
    "excess_median", "excess_high_share",
    "conserved_fast_depth_ratio", "conserved_hit_share",
]


def feature_columns(columns, feature_set="all"):
    """The feature columns among `columns`: all of them, or the normalized ones."""
    available = [c for c in columns if c not in NON_FEATURE_COLUMNS and not c.startswith("meta_")]
    if feature_set == "all":
        return available
    if feature_set == "normalized":
        missing = [c for c in NORMALIZED_FEATURES if c not in available]
        if missing:
            raise RuntimeError("the training data lacks normalized features (written by an older protal?): " + ", ".join(missing))
        return list(NORMALIZED_FEATURES)
    raise ValueError(f"unknown feature set {feature_set!r}")
