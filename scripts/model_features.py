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
features (docs/claude/2026-10-03-r226-v5-v6-training); beside the depth, divergence and unfiltered features
and without the priors (the r226 v10 tables) they changed the test F1 by +0.0003 (pe), +0.0007 (se), +0.0013 (pb) and 0
(ont), within noise (docs/claude/2026-10-04-r226-v10-evaluation).

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
enough to be aligned against it, the share that did not align to it (protal's ZF tag, or the SAM header's counts
when a read aligned nowhere: a relative the database lacks seeds on its nearest species and fails there, a present
species' reads align).

PRIORS_FEATURES ("priors", the same report): per-species constants from GTDB, written by the converter
(species_priors.tsv, protal's SpeciesPriors.h; -1 unknown): the share of the representative's single-copy markers
found twice (rep_duplicate_share, CheckM's contamination signature: a contaminating contig's genes put every
present organism's reads on the species; 0 for every species of a GTDB release, whose marker files hold one copy
per genome, so there rep_contamination carries the signal), its CheckM completeness and contamination, and from GTDB's species
clusters the ANI circumscription radius, mean and minimum intra-species ANI and the cluster's size (log10): a wide
or crowded cluster makes a cloud of reads a few percent from the reference a strain rather than a sister species.
On a synthetic world they are all unknown and do nothing.

REF_FEATURES ("ref", 2026-10-06; in the dump since 2024): the shares of the species' reference k-mers that are unique in
the database's index (su_rate_ref, lu_rate_ref, lsu_rate_ref: super-unique, locally unique, locally super-unique), a
property of the reference and the database, before any read: a reference with few unique k-mers has close relatives
in the database. Unlike the priors they do not depend on how many genomes GTDB has of a species. On the r226 v13
tables they added 0.002-0.003 of F1 in soil to gradient boosting (pe +0.003 soil and shallow soil, pb +0.003, ont
+0.002; the design's test set -0.001 to +0.003), and nothing to a forest (docs/claude/2026-10-06-r226-v13-soil). They
are computed against the database that profiles the sample: the training database lacks its held-out species, so a
reference whose close relatives were held out is more unique there than in the finished database, which the
simulations cannot show (real samples can).

COMPLEXITY_FEATURES ("complexity", 2026-10-06): the sample's complexity, the same for each of its taxa, from all of
them before any is called (protal's context::SampleComplexityOf): sample_log_taxa, log10 of its taxa with fragments;
sample_low_identity, its fragments' share on low-identity bases (the taxa's low_identity_share weighted by their
fragments: how much of the sample is relatives the database lacks); sample_identity, its taxa's fragment-weighted median
identity (of those with 10 fragments or more, if any). A soil sample has thousands of taxa and much of its reads on
relatives the database lacks, a gut sample hundreds; without these a model learns the kind of community from the
sample's depth (sample_log_fragments), which each simulated scenario draws from a narrow band around its preset. At
r226 (v14) pe's shallowest shallow-soil sample, held out whole and below every other soil sample's depth, scored F1
0.758 against 0.920 with species held out; with these three the six shallow samples held out scored 0.935 against
0.908, the hold-out samples +0.001 to +0.003, the design's test set as before. For the other read types, with whole
samples held out: se soil +0.002, PacBio soil and shallow soil +0.002, gut +0.005, Nanopore +0.001 to +0.005, each
read type's design test set -0.0017 (Nanopore) to +0.001 (PacBio) (docs/claude/2026-10-06-r226-v14). In the default
set since 2026-10-06. A dump of a protal before them lacks the three columns. In the r226 v15 build: pe shallow soil
with samples held out 0.908 -> 0.936, gut +0.004 to +0.006 with species held out, no cost on the design's test set at
a fixed knob (docs/claude/2026-10-07-r226-v15).

Three groups against the false positives of complex communities (soil: 59-76% beside a novel species near-identical
to the taxon on the marker genes; docs/claude/2026-10-06-false-positive-features), added 2026-10-06 and not yet
measured at r226 (the next build's tables decide; train without them with --features
normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity):
CONSISTENCY_FEATURES ("consistency"), whether the taxon's reads are its own: long reads' segment records whose hit the
read's consensus taxon gave them (read_consensus_share, ZR:i:1) or that are clearly another taxon's gene
(read_inconsistent_share, ZR:i:2); the mean log2 of the taxa a read's seeds could not tell apart (seed_crowding, ZN,
also those beyond ZA's four alternatives); the share of records whose read a congener fits better than the two
references' distance lets a read of the taxon fit it (unexpected_congener_fit_share, by the database's
species_neighbours.tsv); fragments with a best record on another species of the genus too (split_fragment_share);
and how the taxon's hit genes and those of the congener its reads fit most often complement each other
(congener_gene_overlap, below 1 for one novel species split over two references).
SHAPE_FEATURES ("shape"), how its reads lie on its genes: the dispersion of the genes' divergence around one
factor-scaled genome divergence (gene_divergence_dispersion), covered bases over those expected at the genes' depth
(breadth_ratio), the share of its genes with more failed reads than records (failed_gene_share, ZF genes), and its
consensus' fixed differences (factor-scaled) and polymorphic sites per covered site (fixed_difference_rate,
polymorphic_site_rate).
NEIGHBOURHOOD_FEATURES ("neighbourhood"), the database around the reference (species_neighbours.tsv, written by --build):
congeners within 0.01, 0.02 and 0.05 marker distance and the nearest one's distance (-1 without the table). Like ref,
they describe the database in use: a training database lacks its held-out species.
ANCESTRY_FEATURES ("ancestry", 0.7.9), which side its reads take where the reference differs from its congeners' consensus
(AncestrySites.h; since 2026-10-08 nine in ten of three or more congeners compared, the nearest congener alone with
fewer, as 0.7.9 had everywhere): the sites per record, and the shares where the read has the reference's and the
congeners' base; and the indel sites (the reference's own gain or loss of three bases or more against its congeners,
2026-10-08): per record, and the share where the read deletes or inserts what the congeners do.
GAP_FEATURES ("gaps"), where its reads lie in the gaps to its congeners' copies (congener_gaps.tsv, written by --build: per
copy the alignment distance to the nearest and the median congener's copy): the share of its kept records on copies with a
congener at least 0.005 away (gap_informative_share), of those the shares closer to the reference than the nearest and
than the median congener is (gap_within_min_share, gap_within_median_share), and their median divergence over the nearest
congener's distance (gap_position); -1 without the table or without such a record. FOREIGN_FEATURES ("foreign", in no
named set): how far other species' reads reach its copies in a tiled scan of the genomes at hand (foreign_rates.tsv,
scripts/foreign_rates.py and protal --add_tables): the share of its kept records on scanned copies and their mean shares
of reads from other species and other genera; -1 without the table. The scan reads the simulation's own genomes, so at
r226 v17 these features told the species the samples were drawn from (scanned share ~1) from the rest (~0), the label
in disguise (docs/claude/2026-10-07-congener-gaps, "The foreign features leak"): train on them only to measure that.
UNTRIED_FEATURES ("untried"): the reads whose seeds fit it as well as the taxa they were aligned against but never tried
it, beyond --align_top (ZC), over those plus its reads (untried_candidate_rate).
ALLELE_FEATURES ("alleles", since 2026-10-08), its reads against its species' known strain alleles (strain_alleles.tsv,
written by --build from the full reference): of its kept records' differences from the reference the share their copies'
best alleles explain, and the identity those alleles gain (allele_explained_share, allele_identity_gain); -1 without the
table, 0 for a species without alleles as for reads they do not explain. A strain of the species has its few
differences where its known strains differ, a novel congener mostly elsewhere. The share of its records on copies with
alleles (allele_copy_share, a column of the dump) is in no group: whether a species has alleles is whether GTDB has
other genomes of it, the cluster size the priors carry, which the simulation cannot test. The training database's alleles come from genomes
the simulations never draw strains from (build_gtdb_database.py --allele-genome-share), so that no simulated strain is
its species' own allele.
POLYMORPHIC_FEATURES ("polymorphic", since 2026-10-09), its reads at its species' polymorphic sites (where one of
the species' alleles in strain_alleles.tsv differs from the reference) and at the ancestry sites the species does not
vary at (docs/claude/2026-10-09-r226-v19, section 7): of the polymorphic sites its kept records cover, the share where
the read has a base (or indel) one of the alleles has and the share where it has another (polymorphic_known_share,
polymorphic_novel_share); the share of the species' base at the ancestry sites fixed within the species less that
at all its ancestry sites, weighted by the alleles covering each (ancestry_fixed_gain); and that share at the fixed
sites itself (ancestry_fixed_agreement, since 2026-10-09 evening: a strain near 1 however deep it sits in the species, a
novel congener at the fraction of the species' stem it shares; the models could only rebuild it from
ancestry_agreement and the gain, docs/claude/2026-10-09-ancestry-true-positive-test). A missed real strain sides with
the congener at the congener sites more often than a novel congener does, but less at the fixed ones (v19: AUC 0.59-0.67
against 0.37-0.46). -1 without the table, 0 without a site (a species without alleles too). The polymorphic sites per
kb and the share of the ancestry sites that are fixed (allele_sites_per_kb, ancestry_fixed_share) are columns of the
dump in no group, as allele_copy_share. The shape group's polymorphic_site_rate is the sample's own read variants.
WEIGHT_FEATURES ("weights", since 2026-10-10), its mismatches against how conserved each column of the gene is in the
species' family (column_weights.tsv, written by --build: per column the weight -log(1 - p) of the conservation within
the family's genera, averaged over genera, and among the genus references, and of the amino acid per codon; the odds
that a mismatch at a column is a close relative's real substitution scale with 1 - p, so 99% counts ten times 90%):
the mean weight of its mismatches over the mean weight of its aligned columns (conserved_mismatch_ratio: 1 for
mismatches spread as errors are, below 1 for a relative whose mutations avoid the conserved columns, above 1 for a
read from afar whose fast columns are saturated), its mismatches at columns conserved at ~99% per kb
(conserved_mismatch_rate), the ancestry agreement with each site weighted by its among weight
(ancestry_agreement_weighted; -1 without a site), and of the mismatches classified by codon the share that change the
amino acid (nonsynonymous_share) and those at codon columns whose amino acid is conserved per kb
(nonsynonymous_conserved_rate). -1 without the table (--no_site_weights). Half of r226's false positives were reads
from other genera (docs/claude/2026-10-03-false-positive-anatomy); these see them. column_weight_coverage (the share
of the aligned bases at columns with an estimate) is a column of the dump in no group.

A set's name is its groups joined by "+": normalized, adjacency, relatives or distance, depth, divergence,
unfiltered, ref, complexity, consistency, shape, neighbourhood, ancestry, gaps, foreign, untried, alleles, polymorphic, weights, priors; "all" is every
feature column of the dump. The trainer's default set (DEFAULT_FEATURE_SET) is "normalized+adjacency+distance+depth+
divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood+ancestry+gaps+untried+alleles+polymorphic+weights"
(without weights before 2026-10-10, without polymorphic before 2026-10-09, without alleles before 2026-10-08; with foreign from the congener-gaps merge to 2026-10-08, when the r226 v17 build showed its leak; without gaps,
foreign and untried in 0.7.9 before the congener-gaps merge; without ancestry before 0.7.9; without consistency, shape
and neighbourhood before 2026-10-07, when the false-positive groups were merged; without ref and complexity before
2026-10-06; the set without them is the one every model of the r226 v12 and v13
builds chose with --features auto); the priors are opt-in ("+priors": at r226 +0.007 to
+0.009 of test F1, all of it the cluster size, a bet that a divergent read cloud on a one-genome species is a relative
the database lacks, which the simulation cannot test; for a database of a densely sampled environment,
docs/claude/2026-10-03-r226-v9-evaluation). At r226 the distance features added 0.004 F1
(paired-end) on the test set at the knob curve and 0.007 at knob 0.5 over normalized+adjacency, single-end and long
reads within noise (docs/claude/2026-10-03-r226-v5-v6-training); a table of a protal before them (be35d15) needs
--features normalized+adjacency, one before the depth and divergence features --features
normalized+adjacency+distance, one before the unfiltered and priors features
--features normalized+adjacency+distance+depth+divergence, one before the sample's complexity (2026-10-06, the r226
v14 tables and older) --features normalized+adjacency+distance+depth+divergence+unfiltered+ref, one before the
false-positive groups (the r226 v15 tables and older) --features
normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity. genus_top_fragments, the fragments of the taxon's most
abundant congener, is an input of the singleton rule (machine_learning_cmdline.py --singleton-congener), not a feature
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

# The reference's k-mer uniqueness in the database's index (see above).
REF_FEATURES = ["su_rate_ref", "lu_rate_ref", "lsu_rate_ref"]

# The sample's complexity: its taxa, its share of low-identity bases and its median identity (see above).
COMPLEXITY_FEATURES = ["sample_log_taxa", "sample_low_identity", "sample_identity"]

# What GTDB knows of the species before any read (see above).
PRIORS_FEATURES = ["rep_duplicate_share", "rep_completeness", "rep_contamination", "cluster_ani_radius", "cluster_mean_ani",
                   "cluster_min_ani", "cluster_genomes_log10"]

# Against the false positives of complex communities (see above).
CONSISTENCY_FEATURES = ["read_consensus_share", "read_inconsistent_share", "seed_crowding", "unexpected_congener_fit_share",
                        "split_fragment_share", "congener_gene_overlap"]
SHAPE_FEATURES = ["gene_divergence_dispersion", "breadth_ratio", "failed_gene_share", "fixed_difference_rate",
                  "polymorphic_site_rate"]
NEIGHBOURHOOD_FEATURES = ["db_congeners_01", "db_congeners_02", "db_congeners_05", "db_nearest_congener"]
# Which side the reads take where the reference differs from its nearest congener's (AncestrySites.h, 0.7.9): a
# strain of the species carries the species' base at those sites, a congener that branched off below some of them the
# congener's; the one place in the reads where a missed strain and a novel congener at the same identity differ
# (docs/claude/2026-10-07-error-read-signatures). -1 without a site (no species neighbours, or no congener with the gene).
ANCESTRY_FEATURES = ["ancestry_sites_per_record", "ancestry_agreement", "ancestry_congener_share",
                     "ancestry_indel_sites_per_record", "ancestry_indel_congener_share"]

# Where a taxon's reads lie in the gaps to its congeners' copies (congener_gaps.tsv), how far other species' reads reach its
# copies in a tiled scan of the genomes (foreign_rates.tsv), and the reads whose seeds fit it as well as the taxa they were
# aligned against but never tried it (ZC): see above. In the default set untested at r226.
GAP_FEATURES = ["gap_informative_share", "gap_within_min_share", "gap_within_median_share", "gap_position"]
FOREIGN_FEATURES = ["foreign_scanned_share", "foreign_copy_share", "foreign_genus_copy_share"]
UNTRIED_FEATURES = ["untried_candidate_rate"]
# Its reads against its species' known strain alleles (strain_alleles.tsv): see above. In the default set since 2026-10-08,
# untested at r226.
ALLELE_FEATURES = ["allele_explained_share", "allele_identity_gain"]
# Its reads at its species' polymorphic sites and at the ancestry sites fixed within it: see above. In the default set
# since 2026-10-09, untested at r226.
POLYMORPHIC_FEATURES = ["polymorphic_known_share", "polymorphic_novel_share", "ancestry_fixed_gain",
                        "ancestry_fixed_agreement"]
# Its mismatches against the conservation of the columns in its species' family (column_weights.tsv, written by
# --build since 2026-10-10): see above. In the default set, untested at r226. column_weight_coverage is in no group.
WEIGHT_FEATURES = ["conserved_mismatch_ratio", "conserved_mismatch_rate", "ancestry_agreement_weighted",
                   "nonsynonymous_share", "nonsynonymous_conserved_rate"]

# The groups a set's name may join with "+", in the order they are listed.
FEATURE_GROUPS = {"normalized": NORMALIZED_FEATURES, "adjacency": ADJACENCY_FEATURES, "relatives": RELATIVE_FEATURES,
                  "distance": DISTANCE_FEATURES, "depth": SAMPLE_FEATURES, "divergence": DIVERGENCE_FEATURES,
                  "unfiltered": UNFILTERED_FEATURES, "ref": REF_FEATURES, "complexity": COMPLEXITY_FEATURES,
                  "consistency": CONSISTENCY_FEATURES, "shape": SHAPE_FEATURES, "neighbourhood": NEIGHBOURHOOD_FEATURES,
                  "ancestry": ANCESTRY_FEATURES, "gaps": GAP_FEATURES, "foreign": FOREIGN_FEATURES, "untried": UNTRIED_FEATURES,
                  "alleles": ALLELE_FEATURES, "polymorphic": POLYMORPHIC_FEATURES, "weights": WEIGHT_FEATURES,
                  "priors": PRIORS_FEATURES}

# The sets worth naming (--features takes any groups joined by "+", and "all").
FEATURE_SETS = ("normalized", "normalized+adjacency", "normalized+adjacency+relatives", "normalized+adjacency+distance",
                "normalized+adjacency+distance+depth", "normalized+adjacency+distance+depth+divergence",
                "normalized+adjacency+distance+depth+divergence+unfiltered",
                "normalized+adjacency+distance+depth+divergence+unfiltered+ref",
                "normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity",
                "normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency",
                "normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape",
                "normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood",
                "normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood+ancestry",
                "normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood+ancestry+"
                "gaps+untried",
                "normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood+ancestry+"
                "gaps+untried+alleles",
                "normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood+ancestry+"
                "gaps+untried+alleles+polymorphic",
                "normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood+ancestry+"
                "gaps+untried+alleles+polymorphic+weights",
                "normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood+ancestry+"
                "gaps+untried+alleles+polymorphic+weights+foreign",
                "normalized+adjacency+distance+divergence+unfiltered+priors",
                "normalized+adjacency+distance+depth+divergence+unfiltered+priors",
                "normalized+adjacency+relatives+depth+divergence+unfiltered+priors", "all")
# The priors are not in the default set: their gain at r226 (+0.007 to +0.009 of test F1) is the GTDB cluster size
# alone, a rule that a divergent read cloud on a one-genome species is a relative the database lacks, which the
# simulation cannot test (a one-genome species has no strain to simulate from) and which rejects the strains of
# single-MAG species that dominate environments GTDB has sampled sparsely (docs/claude/2026-10-03-r226-v9-evaluation).
# consistency, shape, neighbourhood, ancestry, gaps and untried are in it untested at r226 (the last two:
# docs/claude/2026-10-07-congener-gaps): the next build's training says what they are worth. foreign is not in it: the
# r226 v17 models learnt from its first scan (of the genomes the samples are drawn from) which species the simulation
# could draw (same report, "The foreign features leak"); the scan reads every species' marker genes alike since
# 2026-10-08, and the named set with foreign lets --features auto say whether that scan earns anything.
# alleles (since 2026-10-08) and polymorphic (since 2026-10-09) are in it untested at r226 too: the next build says what
# the strain alleles are worth.
DEFAULT_FEATURE_SET = ("normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+"
                       "neighbourhood+ancestry+gaps+untried+alleles+polymorphic+weights")
# --features auto: the trainer scores each of these sets with species held out and keeps the best, the default unless
# another beats it by AUTO_MIN_GAIN of F1 at the knob (machine_learning_cmdline.choose_feature_set). Without the priors
# (above) and "all"; auto+priors adds the sets with the priors.
AUTO_FEATURE_SETS = ("auto", "auto+priors")
AUTO_CANDIDATES = tuple(s for s in FEATURE_SETS if s != "all" and "priors" not in s.split("+")) + \
    ("normalized+adjacency+relatives+depth+divergence+unfiltered",)
AUTO_PRIORS_CANDIDATES = AUTO_CANDIDATES + tuple(s for s in FEATURE_SETS if "priors" in s.split("+"))
AUTO_MIN_GAIN = 0.002


def auto_candidates(feature_set):
    """The sets --features auto (or auto+priors) chooses among."""
    return AUTO_PRIORS_CANDIDATES if feature_set == "auto+priors" else AUTO_CANDIDATES


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
    if feature_set != "all" and feature_set not in AUTO_FEATURE_SETS:
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
