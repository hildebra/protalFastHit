# The features of protal's presence model

protal decides whether a species is present with a model of decision trees (gradient-boosted trees
since 2026-10-06, a random forest before; [databases.md](databases.md#training)). For every species with reads it
computes the quantities below from the sample's alignments, the database and the sample's other
taxa, hands them to the model and reports the species whose probability reaches `--knob`. This
page lists every feature: since which version protal computes it, what it measures, how much the
GTDB r226 models lean on it for each read type, and in which situations it matters. How the models
are trained is in [databases.md](databases.md#the-presence-model); the reports the evidence comes from are
in [`docs/claude/`](claude/README.md).

**Where they come from.** `TaxonFeatures` in `src/Profiling/Profiler.h` computes them, in the column
order of the training dump (`<profile>.truth_annotated`), so a model is scored on exactly the
quantities it was trained on. Since 0.7.6 the sums behind them (the identity-weighted bases, the
records' ANIs) are exact, so a feature no longer varies in its last digit with the order in which
multi-threaded alignment wrote the SAM. `scripts/model_features.py` groups them into named sets; a model's
training picks a set with `--features`. The groups, in the order a set's name joins them with `+`:

| group | since | features | in the default set |
|---|---|---|---|
| `normalized` | 0.7.0 (grown in 0.7.1 and 0.7.2) | 34 fractions, ratios and rates that do not depend on the database, the depth or the read length | yes |
| `adjacency` | 0.7.1 (`adjacent_support` 0.7.3) | gene neighbours on the taxon's reads | yes |
| `distance` | 0.7.4 | the taxon against its most abundant congeners, by the references' distance (4 of the 10 `relatives`) | yes |
| `relatives` | 0.7.4 | all 10 features from the sample's other taxa | opt-in, in place of `distance` |
| `depth` | 0.7.5 | the sample's depth | yes |
| `divergence` | 0.7.5 | the reads' divergence by gene conservation and codon position, and the mates | yes |
| `unfiltered` | 0.7.5 | the reads before the MAPQ filter and the reads that failed on the taxon | yes |
| `ref` | 2026-10-06 (in the dump since 2024) | the reference's k-mer uniqueness in the database | yes, since 2026-10-06 |
| `complexity` | 2026-10-06 | the sample's complexity: its taxa, its share of low-identity bases, its median identity | yes, since 2026-10-06 |
| `consistency` | 2026-10-06 | whether the taxon's reads are its own: long reads' consensus tags, the seeds' crowding, congeners that fit better than their distance allows, fragments split with a congener, genes complementing a congener's | yes, untested at r226 |
| `shape` | 2026-10-06 | how the reads lie on the genes: divergence dispersion, breadth against depth, genes where reads fail, fixed and polymorphic sites | yes, untested at r226 |
| `neighbourhood` | 2026-10-06 | the database's congeners near the reference (`species_neighbours.tsv`) | yes, untested at r226 |
| `ancestry` | 0.7.9 | which side the reads take where the reference differs from its nearest congener's | yes, untested at r226 |
| `gaps` | 2026-10-07 | where the reads lie in the gaps to the congeners' copies of their genes (`congener_gaps.tsv`) | yes, untested at r226 |
| `foreign` | 2026-10-07 | how far other species' reads reach the taxon's gene copies in a tiled scan of the genomes (`foreign_rates.tsv`) | no: its scan reads the simulation's genomes, which it tells the models ([below](#the-foreign-features-leak)); in the default set from the merge to 2026-10-08 |
| `untried` | 2026-10-07 | the reads whose seeds fit the taxon as well as the taxa they were aligned against, but never tried it (`ZC`) | yes, untested at r226 |
| `priors` | 0.7.5 | what GTDB knows of the species before any read | opt-in (`+priors`) since 0.7.6; in 0.7.5's default |

The default set is
`normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood+ancestry+gaps+untried`
(78 features; with `foreign` 81, from the congener-gaps merge to 2026-10-08; without `gaps` and `untried` 73, in
0.7.9 before the congener-gaps merge: train such a table with
`--features normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood+ancestry`;
without `ancestry` 70 and without the three groups against false positives 55, both before 2026-10-07; without `ref`
and `complexity` too 49, before 2026-10-06). The r226 v17 build was the first to train the newest groups; its
evaluation is in progress, and its `foreign` features leaked ([below](#the-foreign-features-leak)). A training table of the r226 v15 build or older lacks their columns: train it with
`--features normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity`, one of v14 or older with
`--features normalized+adjacency+distance+depth+divergence+unfiltered+ref`. 0.7.6
adds no feature; it changes how the models are trained (below: the priors opt-in, in-silico strains,
more training depths). 0.6.0a
shipped one model on absolute counts (genes, k-mers and mates); those columns are still in the dump
([below](#in-the-dump-but-in-no-set)) but in no 0.7 set. A training table of an older protal
lacks the columns of later groups: train it with the groups it has (before 0.7.5
`normalized+adjacency+distance`, before 0.7.4 `normalized+adjacency`). Since 2026-10-05 the trainer
chooses the set itself by default (`--features auto`): among the named sets without the priors, with
species held out, keeping the default set unless another is 0.002 of F1 better, and saying why
([databases.md](databases.md#training)). The build did so from 2026-10-05 to 2026-10-06, when every
model of the r226 v12 and v13 builds kept the default set; since then it trains on the default set
(with `ref`), as choosing doubled a boosted model's training.

**How to read the importance columns.** The numbers are the forests' Gini importances (scikit-learn's
`feature_importances_`, which sum to 1 over a model; a boosted model's `varimp.tsv` holds its splits'
gains instead, also summing to 1, not comparable cell for cell) of the two latest forest trainings at GTDB r226, read
from their `trained_model*.varimp.tsv` in `local/v10` and `local/v9`: **v10**
([report](claude/2026-10-04-r226-v10-evaluation/README.md)) is the build of `281a4ba`
(2026-10-04, after 0.7.5) with the default set 0.7.6 trains, without priors; **v9** (`26b065c`,
2026-10-03, 0.7.5's default) is the same design and samples with `+priors` (56 features). Neither
has 0.7.6's in-silico strains or extra training depths; v11 will. Each cell gives
pe / se / pb / ont for v10, with v9 in brackets where it differs by 0.01 or more. An importance is
how much a feature's splits reduce impurity, not how much F1 would fall without it: correlated
features share the credit (the nine unique-k-mer rates, the three excess features), which is why
single rows move between trainings while their groups do not ([group sums](#importance-by-group)).
The feature-set ablations below say what a group is worth in F1.

**Where a feature matters.** The last column names the situations from protal's error analyses
([false positives](claude/2026-10-03-false-positive-anatomy/README.md),
[the fixes and the F1 list](claude/2026-10-03-false-positive-fixes/README.md),
[v9](claude/2026-10-03-r226-v9-evaluation/README.md), and v10's own training log):

- *FP at depth*: a single perfect read on the wrong species in a deep sample (half of r226's paired-end false positives);
- *missing relative*: reads of a species the database lacks landing on its nearest relative (most other false positives; in v10, 320 of 430 held-out false positives are congeners of a held-out species);
- *minor congener*: a present species beside a present congener ten times as abundant whose reads it could hold;
- *thin strain*: a present strain of 1-10 fragments whose genome is 2-4% from the reference (v10: strains miss at 4.6% against 0.9% for representatives, 9.3% at 1-10 fragments; most of r226's misses);
- *shallow sample*: samples of a few thousand reads, where most present taxa have one or two fragments;
- *abundance*: the forest decides presence only, but a few features are the quantities the abundance rests on, so an error in them shows in both;
- *long reads*: features whose signal differs with the error profile of PacBio HiFi or Nanopore reads;
- *contamination*: gene copies that belong to another genus (a contaminated or transferred contig).

## Normalized features (`normalized`)

### Evidence and coverage

| feature | since | what it measures | importance pe / se / pb / ont | matters for |
|---|---|---|---|---|
| `fragments` | 0.7.0 | fragments with an accepted alignment after the MAPQ and length filters: a pair once, a long read once | 0.001 / 0.001 / 0.001 / 0.001 | the amount of evidence; low on its own because the rates below carry it. Shallow samples, thin strains |
| `depth` | 0.7.0 | vertical coverage from the taxon's own reads (those within `--depth_identity_margin` of its best), fragment bases counted once, foreign genes left out | 0.002 / 0.003 / 0.004 / 0.003 | **abundance**: this is the number the profile reports; thin strains |
| `hit_gene_fraction` | 0.7.0 | genes with reads as a fraction of the species' hittable genes | 0.009 / 0.010 / 0.003 / 0.003 (v9 0.020 / 0.009 / 0.007 / 0.005) | missing relative (its reads pile on the conserved genes), shallow samples |
| `gene_presence_ratio` | 0.7.0 | genes hit over the number expected if the fragments fell on the genes by length (about 1 when present, less when the reads gather on a few genes); independent of depth and genome size | 0.016 / 0.013 / 0.012 / 0.011 | missing relative, contamination (reads on one or two genes) |
| `depth_cv` | 0.7.0 | coefficient of variation of the hit genes' depths | 0.003 / 0.003 / 0.001 / 0.001 | missing relative, contamination |
| `variant_sites_per_kb` | 0.7.0 | variant positions (2 or more observations, quality sum 60 or more) per covered kb | 0.002 / 0.001 / 0.005 / 0.001 | a relative's reads add variants: missing relative; strain mixtures |
| `multiallelic_sites_per_kb` | 0.7.0 | positions with 3 or more alleles per covered kb | 0.000 for all | strain mixtures; nearly unused |
| `RAF0` to `RAF4` | 0.6.0a | the share of the quality-filtered variant positions with 0, 1, 2, 3 or 4 alleles (`RAF1` is positions with one alternative allele) | `RAF0`, `RAF1` 0.001-0.002, the rest 0.000 | strain mixtures; nearly unused since the rates per kb exist |

### Read identity

| feature | since | what it measures | importance pe / se / pb / ont | matters for |
|---|---|---|---|---|
| `identity` | 0.7.0 | aligned bases weighted by their read's identity over all aligned bases, indels counted as differences | **0.088 / 0.098 / 0.084 / 0.033** (v9 0.106 / 0.085 / 0.079 / 0.024) | missing relative (reads 3-10% off), thin strain (missed strains have median identity 0.967 against 0.987 for found ones), **abundance** (the margin rule starts from the reads' identities) |
| `top_identity` | 0.7.0 | the 98th percentile of the reads' identities, weighted by aligned bases: what the taxon's best reads look like | 0.039 / 0.055 / **0.136** / 0.036 (v9 0.053 / 0.052 / 0.190 / 0.070) | among the top three for **PacBio HiFi**, where a present species' best reads reach 1.0 and a relative's do not; **abundance** (the depth threshold is `top_identity` minus the margin) |
| `low_identity_share` | 0.7.0 | the share of aligned bases below their gene's own-identity threshold (`top_identity` minus the margin): reads of relatives | 0.001 for all | minor congener and missing relative in principle; the model prefers the direct measures |

### Unique k-mers

The index flags a k-mer that occurs in one species only of the database ("long unique", `lu`: the
whole k-mer), and among those the ones with no k-mer of another species within one difference
("long super unique", `lsu`). A read of the species itself carries many of them; a relative's read,
however well it aligns, carries few. `unique_kmers.tsv` holds the reference's counts per gene, so
the rates below are relative to what the species could show.

Since 2026-10-06 the counts include the k-mers whose 15-base core has a single value in the index
("short unique", `su`). Such a core has no flex cells, so the lookup cannot compare the rest of the
read's k-mer with the entry's; before, these seeds were never counted as unique, though the build
had checked their whole k-mer against the other taxa. The anchor now compares the whole k-mer with
the gene and counts the seed in `ZU` (and in `ZT`: no other value has its core) when it matches; the
reference side (`lu_rate_ref`, `lsu_rate_ref`, the gene rates' denominators) counts them too. At
GTDB r226 they are 5.1M of the index's 2.13 billion unique values (0.24%), so the r226 features
barely move; a small database has many more, and its models must be retrained.

| feature | since | what it measures | importance pe / se / pb / ont | matters for |
|---|---|---|---|---|
| `lu_per_kb` | 0.7.0 | long unique k-mer hits per aligned kb | **0.163 / 0.137 / 0.138 / 0.151** (v9 0.121 / 0.110 / 0.094 / 0.117), the top feature of every v10 model but PacBio's, where it ties for first | the strongest single evidence of a species' own reads; missing relative, minor congener (a congener's reads hit few of the taxon's unique k-mers) |
| `lsu_per_kb` | 0.7.0 | long super unique k-mer hits per aligned kb | 0.038 / 0.052 / 0.031 / 0.053 (v9 0.021 / 0.041 / 0.045 / 0.048) | as above, stricter: long reads with errors still hit them |
| `lu_gene_rate`, `lsu_gene_rate` | 0.6.0a | of the reference's genes that have (super) unique k-mers, the share with at least one hit in the reads | 0.014 / 0.025 / 0.010 / 0.011 and 0.043 / 0.034 / 0.016 / 0.025 | whether the unique evidence spreads over the genome or sits on one gene: contamination, missing relative |
| `lu_gene_rate2`, `lsu_gene_rate2` | 0.6.0a | the same with more than 1 hit per gene | 0.022 / 0.018 / 0.009 / 0.022 and 0.028 / 0.025 / 0.014 / 0.027 (v9 `lu_gene_rate2` 0.045 / 0.042 / 0.024 / 0.025) | as above, at a little depth |
| `lu_gene_rate3`, `lsu_gene_rate3` | 0.6.0a | the same with more than 5 hits per gene | 0.052 / 0.030 / 0.023 / 0.047 and 0.060 / 0.033 / 0.037 / 0.016 (v9 0.062 / 0.020 / 0.032 / 0.058 and 0.054 / 0.030 / 0.025 / 0.033) | the depth-robust form; together the nine unique-k-mer features are a quarter to a third of every model's importance |
| `uniqueness` | 0.6.0a | mates with a unique k-mer hit over all mates | 0.001-0.002 | superseded by the rates per kb |

### The reads' other candidates (0.7.1)

From every read's best record, also those the MAPQ filter leaves out; the alternatives come from the
`ZA` tag protal writes (the other taxa among a read's aligned candidates within 5 edits). A SAM of
an older protal lacks the tag (protal warns, and the fit shares are 0): align it again (`--force`)
for a model that uses them.

| feature | since | what it measures | importance pe / se / pb / ont | matters for |
|---|---|---|---|---|
| `mean_mapq` | 0.6.0a | mean MAPQ of the taxon's records (integer division, as the first model was trained) | 0.002 / 0.003 / 0.003 / 0.002 | minor congener, missing relative (a read that fits two species gets MAPQ 0); low because the EM and distance features say it better |
| `low_mapq_share` | 0.7.1 | the share of records with MAPQ below 10 | 0.001 / 0.002 / 0.006 / 0.007 (v9 pb 0.011) | minor congener; more for long reads, whose alignments span several genes |
| `congener_fit_share` | 0.7.1 | the share of records whose read another species of the genus fits within one edit | 0.003 / 0.005 / 0.003 / 0.004 | minor congener, missing relative (a relative's reads fit several congeners alike) |
| `other_genus_fit_share` | 0.7.1 | the same for a species of another genus | 0.000 for all | contamination, transferred genes; rarely non-zero |

### Divergence beyond the base qualities, and gene conservation (0.7.2)

| feature | since | what it measures | importance pe / se / pb / ont | matters for |
|---|---|---|---|---|
| `excess_median` | 0.7.2 | median over the taxon's best records of the record's differences per aligned base minus the mean error probability of its bases (10^(−Q/10)); 0 without qualities | 0.049 / 0.070 / **0.115** / **0.111** (v9 0.062 / 0.062 / 0.100 / 0.093) | **missing relative**: present taxa exceed their predicted errors by about 0.00 (r226 medians −0.001), absent congeners of a missing species by 0.05-0.06; **long reads**, whose errors differ most from read to read, so that identity alone cannot tell error from divergence |
| `excess_high_share` | 0.7.2 | the share of records with an excess above 0.02 | 0.068 / 0.095 / 0.111 / **0.138** (v9 0.050 / 0.071 / 0.101 / 0.122) | second feature for **Nanopore**; missing relative |
| `conserved_fast_depth_ratio` | 0.7.2 | log2 of the median depth of the hit genes with conservation factor below 1 over the other hit genes (after the filters); 0 without `gene_conservation.tsv` | 0.001 / 0.001 / 0.000 / 0.000 | the pattern it was built for (a relative's reads pile on the conserved genes) is erased by the MAPQ filter: every class has median 0 at r226 |
| `conserved_hit_share` | 0.7.2 | the share of the hit genes that are conserved (factor below 1); 0.5 without the table | 0.001 / 0.001 / 0.002 / 0.001 | missing relative, in principle |
| `conserved_fast_record_ratio` | 0.7.2 | the same depth ratio over every best record before the filters | 0.002 / 0.002 / 0.003 / 0.002 | missing relative: median +0.32 (pe) for absent congeners of a missing species against −0.07 for present taxa, with wide quartiles, so the forest rarely splits on it |
| `conserved_fast_kept_ratio` | 0.7.2 | the same over the records the filters keep | 0.003 / 0.003 / 0.002 / 0.002 | with the record ratio, how many of the taxon's reads on conserved genes were ambiguous |

The 0.7.2 features gained paired-end F1 +0.003 and Nanopore +0.007 per sample on the benchmark
world, mostly through fewer false positives (paired-end 2.12 → 1.50 per sample), the excess pair
carrying the gain ([report](claude/2026-10-01-features-depth-knobs/README.md)); at r226 the excess
features led the importances from the first realistic run on
([v3](claude/2026-10-02-r226-v3-training/README.md)).

## Gene neighbours (`adjacency`, 0.7.1 and 0.7.3)

From the pairs of genes next to each other on the taxon's reads (a pair's mates on two genes, a long
read's consecutive genes), judged by the database's gene neighbours per clade. Without gene
neighbours, or with `--no_gene_neighbours`, they are 0, 0 and 0.5.

| feature | since | what it measures | importance pe / se / pb / ont | matters for |
|---|---|---|---|---|
| `adjacent_expected_share` | 0.7.1 | the share of the pairings whose gene ends face each other in 20% or more of the clade's species | 0.000 / 0.000 / 0.001 / 0.001 | genes that crossed from elsewhere: contamination, transferred genes |
| `adjacent_unlikely_share` | 0.7.1 | the share of pairings seen in 5% of the clade's species or less | 0.000 / 0.000 / 0.001 / 0.000 | as above |
| `adjacent_support` | 0.7.3 | the mean frequency of the pairings in the clade, pulled towards 0.5 by one pairing | 0.000 / 0.000 / 0.002 / 0.001 | as above; long reads carry more pairings |

They separate present from absent taxa hardly at all (AUC 0.39-0.72 on the benchmark world,
[report](claude/2026-10-01-gene-neighbours-run/README.md)): a missing relative's reads pair across
the same genes as the species' own. At r226 the group changes F1 within ±0.001 (the feature-set
table below). They stay in the default set for real data, where gene order differs more between
clades than in the simulated worlds; single-end reads cannot have them (0 importance). A database
whose `gene_neighbours.tsv` has lines of a species' own gene order judges that order expected, so
train on the table the model will run with.

## The sample's other taxa (`distance` and `relatives`, 0.7.4)

Computed once a sample's reads are in (`MicrobialProfile::ApplySampleContext`), after the way
amplicon denoisers judge a candidate against the more abundant sequences that could have produced
it ([report](claude/2026-10-02-amplicon-denoising/README.md)). The four `relative_*` features are
in the default set; the other six are opt-in with `+relatives`. They need training samples in
which congeners share a sample (`build_gtdb_database.py --congeners`, the default): trained on
uniform draws, a model learns that an abundant congener means absence and misses 69-73% of the
minor congeners.

| feature | set | what it measures | importance pe / se / pb / ont | matters for |
|---|---|---|---|---|
| `relative_skew` | distance | log10 of the taxon's fragments (+1) over those of the congener most likely to spill reads onto it (+1); the likeliest source is the congener among the four most abundant with the largest fragments × 0.01 × 10^(−d/0.05), d the median Mash distance of the two references over their shared genes | 0.022 / 0.021 / 0.011 / 0.011 (v9 0.028 / 0.033 / 0.028 / 0.013) | **minor congener** (a thin taxon beside a close abundant congener), **missing relative** (its reads land on several congeners, which then look like each other's spill) |
| `relative_distance` | distance | that congener's distance d (1 without one) | 0.015 / 0.020 / 0.013 / 0.007 (v9 0.024 / 0.034 / 0.022 / 0.027) | tells a close congener's spill from a distant one's: minor congener |
| `relative_spill` | distance | log10 of the taxon's fragments over the spill its abundant congeners would produce | 0.011 / 0.011 / 0.004 / 0.005 | minor congener |
| `relative_close_share` | distance | of the taxon's fragments on the genes it shares with that congener, the share on the half where the two references are most alike, over that half's length share (about 1 for the species' own reads, up to 2 for spilled ones) | 0.003 / 0.003 / 0.002 / 0.001 | minor congener; also the singleton rule's read test |
| `genus_skew`, `family_skew` | relatives | log10 ratio of fragments against the genus's, or the family's other genera's, most abundant species | not in the shipped models | minor congener by rank alone |
| `genus_share`, `genus_spill` | relatives | the taxon's share of its genus's fragments; its fragments over 0.001 × the congeners' plus 0.0001 × the family's | not in the shipped models | as above |
| `em_own_share`, `em_kept_own_share` | relatives | the share of the taxon's best records (all, and those the filters keep) that an abundance-weighted EM over the reads' alternatives leaves to it: a read fitting a taxon m edits worse goes to it in proportion to that taxon's records × (θ/(1−θ))^m, θ that taxon's own divergence in the sample | not in the shipped models (among the top five when trained) | **minor congener**, **thin strain** (its reads tie with a congener's reference and fall to MAPQ 0; the EM gives them back); also the singleton rule's input |

What they are worth: on the benchmark world the four distance features were the best paired-end
set at every knob (+0.001 to +0.006 test F1), single-end within noise; the full relatives set
raised the cross-validated F1 and the log loss but cost up to 0.005 of test F1 at protal's knobs and
missed 1.6-2.3 times as many minor congeners
([report](claude/2026-10-03-denoising-implementation/README.md)). At r226 the distance features add
+0.004 paired-end F1 at the knob curve (+0.007 at knob 0.5) over `normalized+adjacency`; all ten
cut the paired-end log loss by 27% with species held out but did no better at the knob
([v5/v6](claude/2026-10-03-r226-v5-v6-training/README.md)). In v10's held-out table the relatives
set still beats the distance set when the priors are in (pe F1 0.9706 against 0.9691, se 0.9683
against 0.9667). Without the priors, beside the depth, divergence and unfiltered groups, it does
not: retrained on v10's tables, relatives in place of distance changes the test F1 by +0.0003 (pe),
+0.0007 (se), +0.0013 (pb) and 0 (ont), within noise
([v10 follow-up](claude/2026-10-04-r226-v10-evaluation/README.md)). So the default stays `distance`.

## The sample's depth (`depth`, 0.7.5)

| feature | since | what it measures | importance pe / se / pb / ont | matters for |
|---|---|---|---|---|
| `sample_log_fragments` | 0.7.5 | log10 of the sample's fragments over all its taxa, the number a knob curve was read at | 0.017 / 0.022 / 0.007 / 0.008 (v9 0.022 / 0.020 / 0.007 / 0.011) | **FP at depth**: a taxon of one perfect read is a present species in a sample of 5,000 pairs and spill-over in one of 5M; without this feature nothing told them apart. **Shallow samples** on the other side of the same split |

At r226 it halved the false positives at knob 0.5 (pe 218 → 122, se 231 → 118, single-fragment
ones 108 → 39) for a third more misses among short reads, net test F1 +0.004 (pe) to +0.009 (se)
and the log loss −11 to −19% ([anatomy](claude/2026-10-03-false-positive-anatomy/README.md)); v8
against v7 confirmed it on identical samples at GTDB scale
([v7/v8](claude/2026-10-03-r226-v7-v8-evaluation/README.md)). On the 900-species benchmark world,
with few false positives to remove, it costs 0.005-0.007 through misses in deep samples, because
the forest learns the training design's prior for how many rare species a deep sample holds. Two
consequences: the trainer fits no knob curve when the feature is in the set (fitted on top, it
corrects twice and loses 0.010), and a model cannot extrapolate past the deepest training sample,
so the build's design should reach the depths you profile (r226 v9 and v10 train to 30M pairs).
Since 0.7.6 the design also has points at 2,000, 50,000 and 200,000 pairs: v10's test set had 42%
of its paired-end false positives at 50,000 pairs and called 15% of the absent taxa at 2,000, depths
between the points it was trained at.

## The reads' divergence by gene and codon, and the mates (`divergence`, 0.7.5)

| feature | since | what it measures | importance pe / se / pb / ont | matters for |
|---|---|---|---|---|
| `excess_scaled_median` | 0.7.5 | the median excess (as `excess_median`) with each record's excess divided by its gene's conservation factor: the genome's divergence from the reference as the 95% ANI species boundary measures it, rather than the marker genes', which compress a 6-10% divergence into 2-4% | **0.083 / 0.079 / 0.138 / 0.119** (v9 0.054 / 0.071 / 0.073 / 0.108); the top feature of v10's PacBio model | **missing relative** against **thin strain**: places the read cloud on the genome's ANI scale, where a sister species (below 95%) and a strain (above) separate; long reads |
| `excess_conserved_fast_ratio` | 0.7.5 | log2 of that divergence on the conserved genes over the fast ones (0 with under 200 aligned bases on either kind) | 0.001 for all | a species' own reads follow the factors, a relative's reads that align only where the gene is conserved do not; rarely split on |
| `third_position_share` | 0.7.5 | the share of mismatches at third codon positions of the reference (1/3 without mismatches) | 0.005 / 0.004 / 0.033 / **0.110** (v9 0.012 / 0.010 / 0.026 / 0.081) | **long reads**, Nanopore above all: sequencing errors fall on the three positions alike, a strain's differences mostly on the synonymous third, so the share tells error from divergence where the qualities do not |
| `mate_lost_share` | 0.7.5 | of the paired fragments whose mate was expected on the taxon (both mates with a record, or room for the fragment inside the gene), the share whose mate has no record on it; 0 for single-end and long reads | 0.009 / 0 / 0 / 0 | **missing relative** with paired reads: a read that fits a conserved stretch has a mate that fits nowhere on the reference |

On the benchmark world the group gains +0.001 (pe) at knob 0.5 with the log loss down 6-8% and the
misses down (119 → 112) at unchanged false positives; `excess_scaled_median` was the top feature
there too ([report](claude/2026-10-03-false-positive-fixes/README.md)). At r226 (v10, species held
out) it takes the paired-end false positives per sample from 2.34 to 2.19 and F1 from 0.9618 to
0.9641.

## The reads before the filters, and the reads that failed (`unfiltered`, 0.7.5)

| feature | since | what it measures | importance pe / se / pb / ont | matters for |
|---|---|---|---|---|
| `fragments_all` | 0.7.5 | the taxon's reads with a best record before the MAPQ and length filters, a pair or a long read once | 0.001 / 0.002 / 0.001 / 0.001 | thin strains whose reads fell to MAPQ 0; at r226 it is `fragments` plus ~15% and adds nothing |
| `em_fragments` | 0.7.5 | `em_own_share` × `fragments_all`: the fragments a divergent strain would have had, had its reads not tied with a congener's reference | 0.002 / 0.002 / 0.001 / 0.001 | thin strain; adds nothing at r226 (the misses blamed on the MAPQ filter were ties, which `em_own_share` already covered) |
| `failed_candidate_rate` | 0.7.5 | of the reads that seeded on the taxon strongly enough to be aligned against it, the share that did not align to it (protal's `ZF` tag, and the SAM header's counts of the reads that aligned nowhere) | **0.119 / 0.110** / 0.011 / 0.021 (v9 0.125 / 0.162 / 0.012 / 0.039) | **missing relative**: it seeds on its nearest species and fails there (median 0.94 for absent taxa, 0.13 for present ones, 0 for present taxa of one or two fragments). The second feature of the short-read models |

At r226 the group ranks rather than calls: in v10's held-out table it takes the paired-end log loss
from 0.0319 to 0.0309 and the false positives per sample from 2.19 to 2.14 at the same F1
(0.9641 → 0.9643), because what `failed_candidate_rate` separates the depth and distance features
already mostly did at the knob ([v9](claude/2026-10-03-r226-v9-evaluation/README.md)). Low for
long reads, which seed on many genes and rarely fail outright.

## The reference's uniqueness (`ref`, 2026-10-06)

How unique the species' reference is in the database that profiles the sample, before any read: the
shares of its marker genes' k-mers that are unique in the index. protal has written them since 2024
(0.6.0a's model used them with the absolute counts below); no 0.7 set took them until 2026-10-06.

| feature | what it measures | matters for |
|---|---|---|
| `su_rate_ref` | the share of the reference's k-mers that are super-unique (in no other taxon) | **missing relative**: a reference with few unique k-mers has close relatives in the database, and a species the database lacks lands on such references |
| `lu_rate_ref`, `lsu_rate_ref` | the shares that are locally unique and locally super-unique | the same, by the index's local uniqueness |

Unlike the priors they do not depend on how many genomes GTDB has of a species. On the r226 v13
tables they added 0.002-0.003 of soil's F1 to gradient boosting (paired-end +0.003 in soil and
shallow soil, PacBio +0.003, Nanopore +0.002; the design's test set -0.001 to +0.003) and nothing to
a forest ([report](claude/2026-10-06-r226-v13-soil/README.md)): among the soil samples' most
confidently called false positives, the reference's uniqueness was lower (`lu_rate_ref` 0.64 against
0.77 of true positives at the same score and divergence). They are computed against the database in
use: the training database lacks its held-out species, so a reference whose close relatives were
held out is more unique there than in the finished database. The simulations cannot show what that
shift does; real samples can.

## The sample's complexity (`complexity`, 2026-10-06)

What all the taxa of the sample say together, before any is called; the same for every taxon of the
sample (`context::SampleComplexityOf` in `src/Profiling/SampleContext.h`).

| feature | what it measures | matters for |
|---|---|---|
| `sample_log_taxa` | log10 of the sample's taxa with fragments (the profile's rows) | **the kind of community at any depth**: a soil sample has thousands of taxa, a gut one hundreds, a host-dominated one tens |
| `sample_low_identity` | the sample's fragments' share on low-identity bases: the taxa's `low_identity_share` weighted by their fragments | how much of the sample is relatives the database lacks (60% of a soil sample's species in the r226 scenarios) |
| `sample_identity` | the fragment-weighted median identity of the sample's taxa with 10 fragments or more (of every taxon with fragments if none has 10) | the reads' quality and how close the community is to the database |

Without them a model can tell the kind of community only by the sample's depth
(`sample_log_fragments`), and the simulated scenarios draw each kind's samples from a narrow band of
depths around its preset (0.5-2× since 2026-10-06). At r226 (v14) the paired-end model scored its
shallowest shallow-soil sample, held out whole and below every other soil sample's depth, at F1
0.758 against 0.920 with species held out: it took the sample for a design sample of that depth. With
these three features the six shallow samples held out scored 0.935 against 0.908, the hold-out
samples +0.001 to +0.003, the design's test set as before
([report](claude/2026-10-06-r226-v14/README.md)). The other read types, refitted the same way:
with whole samples held out se soil +0.002, PacBio soil and shallow soil +0.002 and gut +0.005,
Nanopore +0.001 to +0.005; on the design's test set -0.0006 (se, at 0.5) to +0.001 (PacBio), Nanopore
-0.0017. A real soil sample at a depth the simulations did not draw is the case they are for. In the
default set since 2026-10-06; the next GTDB build is the first to train with them. A training table
of a protal before them lacks the three columns.

## Against false positives in complex communities (`consistency`, `shape`, `neighbourhood`, 2026-10-06)

In the r226 v13 soil scenarios 95-99% of the false positives were database species beside a congener
the database lacks, 59-76% of them beside one near-identical on the marker genes, and the model told
them apart no better on its own training rows: what was missing was information
([soil report](claude/2026-10-06-r226-v13-soil/README.md)). These 15 features add what the alignment
knew and threw away, and what the database knows of a reference's neighbourhood
([report](claude/2026-10-06-false-positive-features.md)). They are in the default set but untested at
GTDB scale; no importances yet. Three need SAM tags of protal since 2026-10-06 (`ZN`, the gene in
`ZF`; `ZR` is older but was never read); a SAM of an older protal gives them 0.

| feature | group | what it measures | matters for |
|---|---|---|---|
| `read_consensus_share` | consistency | of the taxon's best records, the share of long reads' gene records whose hit or MAPQ the read's consensus taxon gave them (`ZR:i:1`); 0 for short reads | **missing relative** with long reads: a novel species near-identical to the taxon on some genes and to a congener on others splits its reads; a present taxon wins its genes on their own |
| `read_inconsistent_share` | consistency | the share of its records on a gene that is clearly another taxon's than the read's consensus (`ZR:i:2`) | as above, from the side of the taxon that holds the stray genes |
| `seed_crowding` | consistency | the mean log2 of the taxa whose anchors were at least 0.8 as long as the read's longest (`ZN`), also those beyond ZA's four alternatives and align_top that were never aligned | **missing relative** in crowded genera: a read that ten species seed on alike is evidence of none |
| `unexpected_congener_fit_share` | consistency | the share of its records whose read a congener (ZA) fits with so few edits more that a read of the taxon's reference would do so with probability below 0.01 (Poisson, the two references' distance in `species_neighbours.tsv` times the aligned bases); 0 without the table | **missing relative**: ambiguity that the congener's distance does not explain, where `congener_fit_share` also counts the ambiguity of crowded genera |
| `split_fragment_share` | consistency | of its fragments (a pair, a long read), the share with a best record on another species of its genus too | **missing relative**: one novel species split over two references, mate by mate or gene by gene |
| `congener_gene_overlap` | consistency | its hit genes shared with the congener its reads' alternatives name most often, over the overlap expected by their numbers of hit genes ((both + 0.5) / (expected + 0.5)); 1 without such a congener | **missing relative**: a novel species' genes nearer one reference here and the other there give two complementary taxa (below 1) |
| `gene_divergence_dispersion` | shape | Pearson's chi-square over its genes of their records' differences against one genome divergence scaled by each gene's conservation factor plus the expected errors, per degree of freedom; 1 with fewer than 3 genes of 100 aligned bases | **missing relative** against **thin strain**: a strain diverges on every gene by its factor; a sister species near-identical on some markers and diverged on others does not |
| `breadth_ratio` | shape | the bases its genes' kept reads cover over those expected at their depth, L (1 - e^-c) per gene (as inStrain's breadth over expected breadth) | **missing relative** with divergent sources: its reads align where the gene is conserved |
| `failed_gene_share` | shape | of its genes with records or failed reads (the gene in `ZF`, from reads that aligned elsewhere), the share with more failed reads than records | **missing relative**: it fails on the genes where it differs most from the reference |
| `fixed_difference_rate` | shape | sites where a non-reference allele has 80% of at least 4 reads, per covered site weighted by the genes' factors: the consensus' divergence in genome units, free of sequencing errors | **long reads**, whose errors the excess features must subtract; thin strain against relative |
| `polymorphic_site_rate` | shape | sites where a second allele has 20% of the reads and 2 reads or more, per covered site | **minor congener** (the taxon's strain plus a relative's spill-over), strain mixtures |
| `db_congeners_01`, `db_congeners_02`, `db_congeners_05` | neighbourhood | the database's congeners within 0.01, 0.02 and 0.05 of the reference (at most 16); -1 without `species_neighbours.tsv` | **missing relative**: the soil's false positives sit where the database has congeners near-identical on the markers; unlike `ref` it says how near |
| `db_nearest_congener` | neighbourhood | the nearest congener's distance; 1 without one within 0.15, -1 without the table | as above |

Like `ref`, the neighbourhood describes the database in use: a training database lacks its held-out
species, so its species have fewer near congeners than in the finished database. The table's
distance is the median Mash distance of all the marker genes two references have; a run's
`relative_distance` takes it over the genes with unique k-mers only (genes identical in two congeners
have none), so for near-identical congeners the table's distance is the smaller one.

## Which side the reads take (`ancestry`, 0.7.9)

The reads of a false positive are mostly a species the database lacks, landing on its nearest congener at
the same identity as a missed strain's reads on its own species, and nothing in the alignments told them
apart ([report](claude/2026-10-07-error-read-signatures/README.md)). What does differ is *where* the
mismatches fall. Where a species' gene copy differs from its nearest congener's copy, the species has its
derived states (and the congener its own): a strain of the species carries the species' base at those
sites; a species that branched off the lineage below some of them carries the congener's base there. At
run time protal compares each hit gene's copy with the nearest congener's copy in the database: that
gene's nearest by alignment (`congener_gaps.tsv`, since 2026-10-08, the same congener the `gaps` features
measure), else the species' nearest (`species_neighbours.tsv`, nearest first), the first of three whose copy
pairs along their shared 12-mers; keeps the differing positions with the congener's base, and reads each best
record's base at the sites it covers (`AncestrySites.h`; a few hundred bytes per copy, once per run).

| feature | since | what it measures | importance pe / se / pb / ont | matters for |
|---|---|---|---|---|
| `ancestry_sites_per_record` | 0.7.9 | the sites the best records cover, per record | untested | how much the next two can say |
| `ancestry_agreement` | 0.7.9 | of the covered sites, the share where the read has the reference's base; -1 without a site | untested | a strain near 1, a novel congener at the fraction of the branch it shares |
| `ancestry_congener_share` | 0.7.9 | the share where the read has the congener's base; -1 without a site | untested | a novel congener's reads, or a congener's spilling over |

A species without a congener in the database (none for the gene in `congener_gaps.tsv`, none within 0.15
in `species_neighbours.tsv`, or a database without both tables) has no sites: the shares are -1 and the
model falls back on the other features. Species with a second genome could refine the sites further (a site
where the species' own strains vary is no evidence either way); the build does not store that yet.

## The gene copies' gaps and foreign reads, and the untried candidates (`gaps`, `foreign`, `untried`, 2026-10-07)

The reads behind the errors of the r226 v15 build ([report](claude/2026-10-07-error-read-signatures/README.md))
showed false positives fed by species the database lacks, landing on their nearest congener at ~0.97, and misses that
are strains of wide species, their reads at the same identity; and a strain's own species often never aligned against,
ranked below its congeners by the seeds. Three groups add what the reads alone cannot say:

| feature | group | what | against |
|---|---|---|---|
| `gap_informative_share` | gaps | of the taxon's kept records, the share on gene copies whose nearest congener's copy is at least 0.005 away (`congener_gaps.tsv`); 0 without such records, -1 without the table | how much of the evidence can tell the species from its congeners at all |
| `gap_within_min_share`, `gap_within_median_share` | gaps | of those records, the shares whose divergence is below the copy's distance to its nearest congener's copy, and below the median congener's; -1 without such records | **missing relative** (a congener the database lacks lies about as far from the reference as its congeners do) against a **strain** (within the species' gap) |
| `gap_position` | gaps | their median divergence over the nearest congener's distance (0: identical to the reference, 1: as far as the nearest congener) | as above |
| `foreign_scanned_share` | foreign | the share of the kept records on copies the tiled scan reached (`foreign_rates.tsv`); -1 without the table | none: whether the species' genome was scanned ([leak](#the-foreign-features-leak)) |
| `foreign_copy_share`, `foreign_genus_copy_share` | foreign | those copies' mean shares of the scan's reads from other species and from other genera, foreign / (reads + 1); -1 without such records | **contamination, transferred genes, conserved genes**: a copy other species' reads reach is weak evidence of its species |
| `untried_candidate_rate` | untried | the reads whose seeds fit the taxon as well as the taxa they were aligned against (`ZN`'s crowd) but never were aligned against it (beyond `--align_top`, their `ZC` tag), over those plus its reads | **missed strain**: its reads went to a congener without trying it |

`congener_gaps.tsv` is written by `--build`: every species' copy of each marker gene aligned (WFA2, protal's
scores, the ends partly free) against the copies of its genus: all of them up to 24 others, else its 4 nearest by
their k-mer sketches and 16 others drawn by a hash of the pair, whose median gives the median. A read covers a part of
a gene whose divergence varies along it, so the test of one read against the whole gene's gap is noisy: the features
count over a taxon's reads. `foreign_rates.tsv` is made by `scripts/foreign_rates.py` (reads of 150 bases every 500
of every genome at hand, `simulate_metagenomes --tiles`, aligned once; a training database's held-out species left out)
and stored with `protal --add_tables`; `build_gtdb_database.py --foreign-rates` does both ([databases.md](databases.md)).
Both tables describe the database in use, as `ref` does.

### The foreign features leak

The `foreign` group is in no named feature set since 2026-10-08, so neither the default nor `--features auto` trains on
it. The genomes at hand that the scan reads are the genomes the training samples are drawn from. A copy is in the
table only if a scan read landed on it. A species whose genome was scanned has all its copies there
(`foreign_scanned_share` about 1); any other species has only the copies other species' reads reached (about 0, and
`foreign_copy_share` near 1 or -1). The features therefore tell the species the simulation can draw from the rest. In
the r226 v17 build they ranked second and third in every model, and separated present taxa from absent congeners at the
same identity with an AUC of 0.89-0.91 (docs/claude/2026-10-08-r226-v17). In use, the shipped database's scan covers
the same ~50,000 genomes, so the 83% of GTDB species without one would be pushed towards absent. A scan without the
leak would take every species' copies alike (the full reference's marker genes, the held-out species left out) and
count only other species' reads per copy position
([report](claude/2026-10-07-congener-gaps/README.md#the-foreign-features-leak-2026-10-08)).

## The species' priors (`priors`, 0.7.5, opt-in)

Per-species constants from GTDB, written by the converter into `species_priors.tsv` (−1 unknown).
On a synthetic world they are all unknown and do nothing. Not in v10's models; the importances
are v9's.

| feature | since | what it measures | importance pe / se / pb / ont (v9) | matters for |
|---|---|---|---|---|
| `rep_duplicate_share` | 0.7.5 | the share of the representative's single-copy markers found twice | 0.000 for all | **contamination**; 0 for every species of a GTDB release (its marker files hold one copy per genome), so the signal is `rep_contamination`'s |
| `rep_completeness`, `rep_contamination` | 0.7.5 | the representative's CheckM completeness and contamination | 0.001-0.002 | contamination; nothing measurable at r226 |
| `cluster_ani_radius` | 0.7.5 | the species cluster's ANI circumscription radius | 0.000-0.001 | 95 for nearly every species, so nothing |
| `cluster_mean_ani`, `cluster_min_ani` | 0.7.5 | mean and minimum intra-species ANI of the cluster (unknown for one-genome clusters) | 0.004 / 0.005 / 0.007 / 0.006 each | with the size: a wide or crowded cluster makes a read cloud a few percent from the reference a strain, not a sister species |
| `cluster_genomes_log10` | 0.7.5 | log10 of the cluster's genomes | 0.004 / 0.004 / 0.006 / 0.004 | **missing relative** against **thin strain** on one-genome species |

The group is worth +0.005 of paired-end and single-end F1 with species held out (v10's table:
0.9643 → 0.9691, 0.9616 → 0.9667) and +0.007 to +0.009 on v9's test set, all of it the
cluster size: the forest learns that a divergent read cloud on a one-genome species is a relative
the database lacks, never a strain of it, because the simulation has no second genome to make such
a strain from. That is a fair bet for environments GTDB has sampled densely and rejects the strains
of single-MAG species elsewhere, so the priors are opt-in since 0.7.6, until checked on real samples
([v9](claude/2026-10-03-r226-v9-evaluation/README.md)). Without them v10 calls absent one-genome
species 2.9 times as often and misses more strains of multi-genome species
([v10](claude/2026-10-04-r226-v10-evaluation/README.md)). 0.7.6's training gives every one-genome
species an in-silico strain: a codon-aware mutated copy of its representative, as far from it as
real strains are (`scripts/insilico_strains.py`). In such training data the cluster size can no
longer stand for "has no strains", so v11 shows what the priors are worth beyond that rule.

## Importance by group

Sum of the v10 importances per group and read type, with the normalized group split by kind
(v9 in brackets where the priors moved it):

| group | pe | se | pb | ont |
|---|---|---|---|---|
| unique k-mers (9) | 0.421 (0.364) | 0.357 (0.309) | 0.279 (0.264) | 0.353 (0.326) |
| identity (3) | 0.127 (0.159) | 0.154 (0.138) | 0.221 (0.270) | 0.069 (0.094) |
| excess beyond base qualities, 0.7.2 (2) | 0.116 | 0.166 (0.133) | 0.227 (0.200) | 0.249 (0.215) |
| coverage and alleles (12) | 0.034 | 0.033 | 0.028 | 0.022 |
| other candidates, 0.7.1 (4) | 0.007 | 0.010 | 0.012 | 0.013 |
| conservation pattern, 0.7.2 (4) | 0.006 | 0.007 | 0.008 | 0.006 |
| `adjacency` (3) | 0.001 | 0.000 | 0.004 | 0.002 |
| `distance` (4) | 0.051 (0.068) | 0.055 (0.083) | 0.029 (0.055) | 0.024 (0.050) |
| `depth` (1) | 0.017 (0.022) | 0.022 | 0.007 | 0.008 (0.011) |
| `divergence` (4) | 0.097 (0.073) | 0.084 | 0.173 (0.099) | 0.231 (0.190) |
| `unfiltered` (3) | 0.123 | 0.114 (0.166) | 0.013 | 0.023 (0.041) |
| `priors` (7) | — (0.015) | — (0.017) | — (0.025) | — (0.020) |

The short-read models lean on the unique k-mers, the identity and the failed candidates; the
long-read models on the divergence features and `top_identity`, where the reads' errors must be
told from a genome's divergence base by base. With the priors in (v9), the forest shifts credit
from the unique k-mers and the divergence features to the distance features and, for single-end
reads, to `failed_candidate_rate`.

## What the groups are worth: the r226 feature-set ablation

The trainer refits the forest on every named set and scores it with the species held out (v10,
`local/v10/classifier_training*.log`, "Feature sets", rows held out by species; knob 0.5; v9
gives the same numbers for the sets both have, so the shared rows are one table). F1 and false
positives per sample:

| set | pe F1 / FP | se F1 / FP | pb F1 / FP | ont F1 / FP |
|---|---|---|---|---|
| `normalized` | 0.9486 / 3.27 | 0.9405 / 3.75 | 0.9654 / 1.11 | 0.9626 / 1.34 |
| `+adjacency` | 0.9485 / 3.28 | 0.9414 / 3.71 | 0.9664 / 1.05 | 0.9622 / 1.33 |
| `+adjacency+distance` | 0.9546 / 3.00 | 0.9494 / 3.35 | 0.9665 / 1.02 | 0.9659 / 1.26 |
| `+adjacency+relatives` | 0.9590 / 2.57 | 0.9550 / 2.80 | 0.9662 / 1.02 | 0.9658 / 1.27 |
| `+adjacency+distance+depth` | 0.9618 / 2.34 | 0.9592 / 2.51 | 0.9679 / 0.97 | 0.9663 / 1.25 |
| `+adjacency+distance+depth+divergence` | 0.9641 / 2.19 | 0.9595 / 2.41 | 0.9676 / 0.96 | 0.9681 / 1.15 |
| **the default** (`+unfiltered`; v10's models, 0.7.6's default) | **0.9643 / 2.14** | **0.9616 / 2.29** | **0.9692 / 0.91** | **0.9685 / 1.11** |
| the default + `priors` (v9's models, 0.7.5's default) | 0.9691 / 1.79 | 0.9667 / 2.04 | 0.9735 / 0.79 | 0.9722 / 0.98 |
| `+relatives` in place of `distance`, + `priors` | 0.9706 / 1.69 | 0.9683 / 1.94 | 0.9729 / 0.78 | 0.9720 / 0.99 |
| `all` (every dump column, 104) | 0.9702 / 1.72 | 0.9667 / 1.99 | 0.9728 / 0.83 | 0.9713 / 1.04 |

On the independent test set (78 paired-end samples of 500 to 5M pairs, sigma 2.0) v10's models
score F1 0.9630 (pe), 0.9596 (se), 0.9675 (pb) and 0.9689 (ont) at knob 0.5, false positives
2.15 / 2.26 / 0.90 / 1.10 per sample; v9's, with the priors, 0.9699 / 0.9668 / 0.9727 / 0.9741.

Each group added in 0.7.4 and 0.7.5 lowers the false positives and raises F1 for short reads; the
adjacency group changes nothing, and for long reads the gains are smaller and mostly in the log
loss. `all`, which includes 0.6.0a's absolute counts, does no better than the named sets, so
nothing is left in the dump that the sets miss. Held out by rows instead of by species the numbers
are the same within 0.002, so no set has learnt the training species.

## In the dump, but in no set

The training dump also holds columns that no 0.7 set uses:

- **0.6.0a's absolute counts**, the inputs of the model 0.6.0a shipped (`present_genes`,
  `total_hits`, `unique_hits`, `mean_ani`, `expected_gene_presence` and its ratio, `variance1`,
  `variance2`, `stddev`, `hittable`, `lu`, `lsu`, their gene counts, the reference's unique k-mer
  counts `su_genome`, `lu_genome`, `lsu_genome`, `total_genome` (their rates are `ref`'s), `lu_per_read`,
  `lsu_per_read`, the allele counts `A0`-`A4`, `AF0`-`AF4`, `RA0`-`RA4`). They depend on the
  database (archaea have 52 marker genes, bacteria 119), the depth and the read length, which is
  why 0.6.0a found a third of the archaea present with probabilities pinned near 0.5 and why 0.7.0
  retrained on the normalized features. protal still scores a database's old `model.xml` with them;
  on v10's test set that model finds 1.3% of the archaea and 75% of the bacteria (F1 0.788
  against 0.963).
- `gene_dispersion` (0.7.0): a chi-square of the reads over the genes; it separated present from
  absent taxa hardly at all (AUC 0.43) and grew with depth.
- `linked_share` (0.7.1): the share of reads with two records on the taxon (both mates, or two
  genes of a long read); it lowered the paired-end F1 on the tuning world and gained 0.001 at r226.
- `genus_top_fragments` (0.7.4): the fragments of the genus's most abundant other species, an input
  of the singleton rule (`--singleton_congener`), not a feature.

## Situations and the features that address them

| situation | features that carry it | what the evidence says |
|---|---|---|
| **False positives at depth** (a perfect single read on the wrong species) | `sample_log_fragments`; with it `identity`, `top_identity`, `lu_per_kb` | r226: 52% of paired-end false positives; the depth feature halved them. The remaining ones are the species boundary on marker genes: 3+ reads of a strain 95-97% from the reference |
| **A relative the database lacks** (reads on the nearest species) | `failed_candidate_rate`, `excess_median`, `excess_high_share`, `excess_scaled_median`, `lu_per_kb`, `identity`, `relative_skew`, `mate_lost_share`, `gene_presence_ratio` | the main source of false positives in every error analysis (v10: 320 of 430 held-out false positives are congeners of a held-out species); the training database leaves out whole clades so that the model sees them |
| **A minor congener beside an abundant one** | `relative_skew`, `relative_distance`, `relative_spill`, `relative_close_share`, `em_own_share` (opt-in), `congener_fit_share`, `low_mapq_share` | the distance features gain here only when trained with congener groups; the full relatives set over-rejects (1.6-2.3× the minor-congener misses) |
| **Thin strains** (1-10 fragments, 2-4% from the reference: the false negatives) | `identity`, `excess_scaled_median`, `em_own_share`/`em_fragments`, `cluster_genomes_log10` (opt-in), `fragments` | v10: strains miss at 4.6% against 0.9% for representatives, 9.3% at 1-10 fragments; missed strains' identity median 0.967. Half the misses had most records below MAPQ 10. Since 0.7.6 one-genome species are trained with in-silico strains too. The remaining lever is strain alleles in the index, not a feature |
| **Shallow samples** (a few thousand reads) | `sample_log_fragments`, `fragments`, `hit_gene_fraction` | at 500-1,000 pairs the depth feature lets a single read count; the knob curve did the same more coarsely |
| **A complex community at a depth the training lacked** (soil at 2M pairs) | `sample_log_taxa`, `sample_low_identity`, `sample_identity` | r226 v14: with the depth alone, a shallow-soil sample below the simulated soil depths was scored like a design sample (F1 0.758); with the sample's complexity 0.935 over the six held out |
| **Abundance estimates** | `depth`, `top_identity`, `identity` (the margin rule), `conserved_*`, and the foreign-gene exclusion behind `depth` | not a forest decision: abundance is the depth from the reads within `--depth_identity_margin` (0.08) of `top_identity`, foreign genes left out, fragment bases once ([version_changes.md](version_changes.md)) |
| **Long reads** (errors differ read by read) | `excess_scaled_median`, `excess_median`, `excess_high_share`, `top_identity` (PacBio), `third_position_share` (Nanopore), `low_mapq_share` | PacBio HiFi models were trained on reads with real qualities since 0.7.3; the excess features need base qualities and are 0 without them |
| **Contamination and transferred gene copies** | the suspect-copy scan at build (`suspect_copies.tsv`, not a feature), then `gene_presence_ratio`, `depth_cv`, `other_genus_fit_share`, `adjacent_*`, `rep_contamination` | 47% of r226's cross-genus false positives were reads on copies near-identical to another genus's; the scan removes them from the evidence before the features see them (v10: 5,589 of 14.5M copies) |
| **Archaea and databases of another size** | the whole normalized group (fractions, ratios, rates per kb) in place of 0.6.0a's counts | archaea recall 0.24 (0.6.0a) → 0.74 (0.7.0+) on the benchmark world; 0.013 → 0.978 on v10's r226 test set |

## Sources

Importances and ablations: `local/v10/model_logs/trained_model*.varimp.tsv`,
`local/v10/classifier_training*.log` (the build of 2026-10-04, `281a4ba`;
[report](claude/2026-10-04-r226-v10-evaluation/README.md)) and the same files in `local/v9`
(git-ignored). Definitions: `src/Profiling/Profiler.h` (`TaxonFeatures`, the `Taxon`
accessors), `src/Profiling/SampleContext.h`, `scripts/model_features.py`. Evidence: the reports
linked above, all in [`docs/claude/`](claude/README.md).
