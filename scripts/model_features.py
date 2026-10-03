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
conserved_hit_share; 0 and 0.5 without it). A relative's reads exceed their errors, and align best on the conserved
genes, where it differs least from its congeners (docs/claude/2026-10-01-conservation-pattern). The last one,
conserved_fast_record_ratio, is the same depth ratio over every best record before the MAPQ filter, which drops most
of a relative's reads on conserved genes (they fit several congeners equally): the drop from conserved to fast genes
shows there, and blurs after the filter; conserved_fast_kept_ratio is the ratio of the records the filters keep, so the
two together say how much of the conserved genes' reads were ambiguous.

ADJACENCY_FEATURES, by the database's gene neighbours (docs/claude/2026-10-02-gene-neighbour-frequencies): of the genes
next to each other on a taxon's reads (a pair's mates on two genes, a long read's consecutive genes), the shares that
are expected and unlikely neighbours in the taxon's clade, and the mean frequency of their pairings there.
"normalized+adjacency" is the normalized features and these; "normalized" leaves them out, to test them on real data:
on a synthetic world they changed test F1 within noise. Without gene neighbours in the
database (or with protal --no_gene_neighbours) they are 0, 0 and 0.5, as for a taxon with no reads across genes.

RELATIVE_FEATURES, from the other taxa of the taxon's sample (protal's MicrobialProfile::ApplySampleContext,
docs/claude/2026-10-02-amplicon-denoising): how its fragments compare with its genus's and family's most abundant other
species (genus_skew, family_skew, genus_share), against the fragments they would spill onto it by rank (genus_spill) and
by the distance between their references (relative_skew, relative_distance, relative_spill), how its reads lie on the
genes where its reference and the likeliest source's are most alike (relative_close_share), and the share of its reads an
abundance-weighted assignment over the reads' alternatives leaves to it (em_own_share, em_kept_own_share). They let a
model call a thin taxon with no relative in the sample and reject one beside an abundant congener whose reads it holds;
trained on samples whose species are drawn uniformly, where congeners hardly ever share a sample, the model learns to
reject every taxon beside an abundant congener, so train them on samples with congener groups
(collect_training_data.py --congeners SHARE:MIN-MAX, build_gtdb_database.py's default). They are the set
"normalized+adjacency+relatives", opt-in: on the benchmark world they raised cross-validated F1, AP and log loss, but cost
0.001-0.004 of the test sets' F1 at the knobs protal calls with, and missed more minor congeners
(docs/claude/2026-10-03-denoising-implementation); at GTDB r226 they lowered the paired-end log loss by 27% and the false
positives by 29% with species held out, yet on the test set at the knob curve they did no better than the four distance
features (docs/claude/2026-10-03-r226-v5-v6-training).

SAMPLE_FEATURES ("depth"): the sample's depth, log10 of its fragments over all its taxa, the number protal reads a
model's knob curve at. Without it no feature says how deep the sample is, and what a taxon of one perfect read is worth
depends on nothing else: at GTDB r226 half of the paired-end false positives were such reads, equal to true single-read
species in every feature and apart from them only by the sample's depth. As a feature it halved the false positives at
knob 0.5 and raised the test F1 by 0.004 (pe) to 0.009 (se) over the distance set, above its knob curve; a knob curve
fitted on top of it corrects twice and loses, so the trainer fits none when the feature is in the set
(docs/claude/2026-10-03-false-positive-anatomy).

DIVERGENCE_FEATURES ("divergence", the same report): excess_scaled_median, the reads' divergence beyond their base
qualities with each read's excess divided by its gene's conservation factor, the genome's divergence from the reference
as the species definition (95% ANI) measures it rather than the marker genes', which compress it;
excess_conserved_fast_ratio, log2 of that divergence on the conserved genes over the fast ones (a species' own reads
follow the factors, a relative's reads that align only where the gene is conserved do not); third_position_share, the
share of the mismatches at third codon positions (errors fall on all three alike, a strain's differences mostly on the
third); and mate_lost_share, of the paired fragments whose mate was expected on the taxon (both mates with a record, or
room for the fragment inside the gene), the share whose mate has no record on it (a relative's read that fits a
conserved stretch has a mate that fits nowhere). A dump of a protal before them lacks the five columns.

UNFILTERED_FEATURES ("unfiltered", docs/claude/2026-10-03-false-positive-fixes): fragments_all, the taxon's reads
with a best record before the MAPQ and length filters (a pair or a long read once); em_fragments, of them what the
abundance-weighted assignment leaves to the taxon (em_own_share x fragments_all: the fragments a divergent strain
would have had, had its reads not tied with a congener's reference and fallen to MAPQ 0, which is where 159 of the
328 r226 misses lost their evidence); and failed_candidate_rate, of the reads that seeded on the taxon strongly
enough to be aligned against it, the share that did not align to it (protal's ZF tag, also on an unmapped record
when a read aligned nowhere: a relative the database lacks seeds on its nearest species and fails there, a present
species' reads align).

PRIORS_FEATURES ("priors", the same report): per-species constants from GTDB, written by the converter
(species_priors.tsv, protal's SpeciesPriors.h; -1 unknown): the share of the representative's single-copy markers
found twice (rep_duplicate_share, CheckM's contamination signature: a contaminating contig's genes put every
present organism's reads on the species), its CheckM completeness and contamination, and from GTDB's species
clusters the ANI circumscription radius, mean and minimum intra-species ANI and the cluster's size (log10): a wide
or crowded cluster makes a cloud of reads a few percent from the reference a strain rather than a sister species.
On a synthetic world they are all unknown and do nothing.

A set's name is its groups joined by "+": normalized, adjacency, relatives or distance, depth, divergence,
unfiltered, priors; "all" is every feature column of the dump. The trainer's default set (DEFAULT_FEATURE_SET) is
"normalized+adjacency+distance+depth+divergence+unfiltered+priors". At r226 the distance features added 0.004 F1
(paired-end) on the test set at the knob curve and 0.007 at knob 0.5 over normalized+adjacency, single-end and long
reads within noise (docs/claude/2026-10-03-r226-v5-v6-training); a table of a protal before them (be35d15) needs
--features normalized+adjacency, one before the depth and divergence features --features
normalized+adjacency+distance, one before the unfiltered and priors features
--features normalized+adjacency+distance+depth+divergence. genus_top_fragments, the fragments of the taxon's most
abundant congener, is an input of the singleton rule (random_forest_cmdline.py --singleton-congener), not a feature
of these sets: it counts reads.
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
    "conserved_fast_record_ratio", "conserved_fast_kept_ratio",
]


ADJACENCY_FEATURES = ["adjacent_expected_share", "adjacent_unlikely_share", "adjacent_support"]

RELATIVE_FEATURES = ["genus_skew", "family_skew", "genus_share", "genus_spill", "relative_skew", "relative_distance",
                     "relative_spill", "relative_close_share", "em_own_share", "em_kept_own_share"]
# Of them, those by the distance of the references: on the benchmark world the best paired-end set at every knob
# (+0.001 to +0.006 test F1 over normalized+adjacency), single-end within noise, but more minor congeners missed
# (docs/claude/2026-10-03-denoising-implementation); at r226 +0.004 paired-end at the knob curve, the rest within noise
# (docs/claude/2026-10-03-r226-v5-v6-training). In the default set.
DISTANCE_FEATURES = ["relative_skew", "relative_distance", "relative_spill", "relative_close_share"]

# The sample's depth (see above): with it in a set, the trainer fits no knob curve.
SAMPLE_FEATURES = ["sample_log_fragments"]

# The reads' divergence by gene conservation and codon position, and the mates (see above).
DIVERGENCE_FEATURES = ["excess_scaled_median", "excess_conserved_fast_ratio", "third_position_share", "mate_lost_share"]

# The reads before the filters and the reads that failed on the taxon (see above).
UNFILTERED_FEATURES = ["fragments_all", "em_fragments", "failed_candidate_rate"]

# What GTDB knows of the species before any read (see above).
PRIORS_FEATURES = ["rep_duplicate_share", "rep_completeness", "rep_contamination", "cluster_ani_radius", "cluster_mean_ani",
                   "cluster_min_ani", "cluster_genomes_log10"]

# The groups a set's name may join with "+", in the order they are listed.
FEATURE_GROUPS = {"normalized": NORMALIZED_FEATURES, "adjacency": ADJACENCY_FEATURES, "relatives": RELATIVE_FEATURES,
                  "distance": DISTANCE_FEATURES, "depth": SAMPLE_FEATURES, "divergence": DIVERGENCE_FEATURES,
                  "unfiltered": UNFILTERED_FEATURES, "priors": PRIORS_FEATURES}

# The sets worth naming (--features takes any groups joined by "+", and "all").
FEATURE_SETS = ("normalized", "normalized+adjacency", "normalized+adjacency+relatives", "normalized+adjacency+distance",
                "normalized+adjacency+distance+depth", "normalized+adjacency+distance+depth+divergence",
                "normalized+adjacency+distance+divergence+unfiltered+priors",
                "normalized+adjacency+distance+depth+divergence+unfiltered+priors",
                "normalized+adjacency+relatives+depth+divergence+unfiltered+priors", "all")
DEFAULT_FEATURE_SET = "normalized+adjacency+distance+depth+divergence+unfiltered+priors"


def feature_set_columns(feature_set):
    """The columns a set's name stands for, in group order; ValueError for a name that is no groups joined by "+"."""
    groups = feature_set.split("+")
    unknown = [g for g in groups if g not in FEATURE_GROUPS]
    if unknown or len(set(groups)) != len(groups) or "normalized" not in groups:
        raise ValueError(f"unknown feature set {feature_set!r}: groups joined by '+' from " + ", ".join(FEATURE_GROUPS) +
                         " (normalized among them), or all")
    chosen = []
    for group in FEATURE_GROUPS:
        if group in groups:
            chosen += [c for c in FEATURE_GROUPS[group] if c not in chosen]
    return chosen


def feature_set_name(feature_set):
    """argparse type for --features: the name checked (ValueError, which argparse reports, for an unknown one)."""
    if feature_set != "all":
        feature_set_columns(feature_set)
    return feature_set


def has_sample_depth(columns):
    """Whether the chosen feature columns hold the sample's depth (SAMPLE_FEATURES): then no knob curve is fitted."""
    return any(c in columns for c in SAMPLE_FEATURES)


def feature_columns(columns, feature_set="all"):
    """The feature columns among `columns`: all of them, or those of the set's groups (feature_set_columns), which the
    table must have."""
    available = [c for c in columns if c not in NON_FEATURE_COLUMNS and not c.startswith("meta_")]
    if feature_set == "all":
        return available
    chosen = feature_set_columns(feature_set)
    missing = [c for c in chosen if c not in available]
    if missing:
        raise RuntimeError(f"the training data lacks {feature_set} features (written by an older protal?): " +
                           ", ".join(missing))
    return list(chosen)
