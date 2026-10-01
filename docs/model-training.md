# Training protal's presence model

protal decides whether a species is present with a random forest (PMML, `model_pe.xml` in the
database for paired-end reads; `model.xml` in databases of earlier versions). It scores every species with reads, and reports those whose probability of `TRUE` is at
least `--knob` (0 to 1, default 0.5).

The model's features come from the reads and the database, so a model belongs to the kind of
database it was trained on. [building-a-database.md](building-a-database.md#build-and-train-in-one-command)
runs the whole pipeline below in one command; this page describes its parts.

Features are counted per read, and a pair of reads counts twice, so a model also belongs to the
kind of reads it was trained on. Single-end samples are profiled with `model_se.xml` of the
database (or `--model_se`; [running.md](running.md#single-end-reads)), a model trained the same way
on single-end samples, and PacBio and Nanopore samples with `model_PB.xml` and `model_ONT.xml`
(`collect_training_data.py --read_types`).

## What protal expects of a model

At start, before any alignment, protal checks that the model
- takes only inputs protal computes: the features of `TaxonFeatures` in
  `src/Profiling/Profiler.h`, the columns of the training dump below;
- predicts a field with the value `TRUE`, whose probability is the species' score;
- scores a species.

A model that fails one of these stops the run with exit code 2, naming the problem.

## Training data

With a truth file (`--profile_truth`, one per sample, comma-separated, or a `PROFILE_TRUTH` column
in the map), protal writes `<profile>.truth_annotated`. This file has one row per species with
reads: whether the species was present (`truth`), the current model's call and probability, and
every feature. The features are written exactly as the model is given them.

A truth file names one species present per line: by GTDB lineage in any tab-separated field
(`d__...;s__Genus species`) or by internal taxid in the first field. `simulate_metagenomes
--protal_metafile` writes such files and a map that uses them ([simulation.md](simulation.md)).

`collect_training_data.py` produces such data at scale. It simulates metagenomes with
`simulate_metagenomes` over a grid of read setups (length, ART profile, fragment size) and depths,
profiles them against a database, and joins all dumps into `training_data.tsv`:

    python3 scripts/collect_training_data.py --db DB --genome_table genomes.tsv -o training \
        --archaea 2 --species_per_sample 10-40 -t 8

With `--read_types`, it collects other reads of the same communities too: `se`, the paired-end
samples' first reads alone, profiled as single-end reads (no new simulation); `pb` and `ont`,
long reads simulated with pbsim3 from each paired-end sample's manifest, every genome given its
share of `--long_read_bases` by abundance times length (point i of the long reads replays the
communities of paired-end point i; contigs under 100 bases are left out of a genome before pbsim3
reads it, since pbsim3 stops at them). All samples are profiled in one protal run, each with its read
type's model and settings (the map's `READ_TYPE` column), and each read type gets its table:
`training_data.tsv` (pe), `training_data_se.tsv`, `training_data_pb.tsv`, `training_data_ont.tsv`.

| Option | Default | |
|---|---|---|
| `--db`, `--genome_table`, `-o` | required | database, simulator genome table (accession, GTDB taxonomy, FASTA path), output directory |
| `--protal`, `--simulator` | on `$PATH` | the binaries |
| `--samples` | 4 | samples per design point |
| `--read_pairs` | `1000,5000,20000,100000,500000` | depths, one design point each |
| `--read_setups` | `100:HS20:300:40,150:HSXt:350:50,250:MSv3:550:50` | length : ART profile : fragment mean : fragment SD, one design point each; the profile `file=R1.txt+R2.txt` (or `file=P.txt`) uses quality profiles of `art_profiler_illumina` instead of a built-in one |
| `--species_per_sample` | `5-30` | N or MIN-MAX |
| `--strains_per_species` | one strain | probabilities of a second, third, ... strain of a species, e.g. `0.3,0.1` |
| `--abundance` | the simulator's (Poisson-lognormal, sigma 1.3) | `lognormal:SIGMA`, `powerlaw:ALPHA` or `negbin:R:P` |
| `--read_types` | `pe` | `pe`, `se`, `pb`, `ont`, comma-separated |
| `--long_read_bases` | `300000,1500000,6000000,30000000,150000000` | bases per long-read sample, one design point each |
| `--pb_setup`, `--ont_setup` | `errhmm:ERRHMM-SEQUEL:15000:3000:0.999`, `qshmm:QSHMM-ONT-HQ:8000:6000:0.97:39/24/36` | pbsim3 method : model : length mean : length SD : accuracy mean (: error mix, for qshmm) |
| `--pbsim`, `--pbsim_models` | `pbsim`, its data folder | pbsim3 and the folder of its `.model` files |
| `--archaea` | 0 | archaeal species per sample |
| `--congeners` | 0 | species of one genus in every sample of a design point (the genus drawn per point): relatives share real samples, but hardly ever uniform draws from many genera |
| `--novel_species` | | species the database lacks (e.g. those a training database leaves out), optionally with the rank they were held out at and the clade (`heldout_species.txt`), for `meta_novel_*` |
| `--novel_clades` | 0 | species of held-out clades (ranks above species in `--novel_species`) in every sample, per rank: each design point takes one clade of each rank, in turn |
| `--taxonomy` | | the database's `internal_taxonomy.dmp`, for `meta_rep_genome`, `meta_relative_rank`, `meta_novel_level` and `meta_neighbour_rank` |
| `-t`, `--seed` | 4, 1 | threads of the protal run; seed |
| `--jobs` | `-t` | design points (and long-read genomes) simulated at a time (ART and pbsim3 simulate one genome at a time) |

The design points are simulated in parallel, then all their samples are profiled in one protal run,
which loads the database once. A rerun resumes: a design point already simulated from the same
inputs (simulator, genome table, seed and design) and profiled against the same database with the
same protal is skipped; one made from other inputs is simulated or profiled again (its folder holds
the key of what made it, `simulated.json` and `profiled.json`). Columns it adds start with
`meta_` and say where each row comes from: design point, sample, read type (`meta_read_type`), read
length (the mean for long reads), depth (`meta_read_pairs`: read pairs, reads for se, bases for pb
and ont), the taxon's
domain (`meta_domain`, from the genome table's lineages; `unknown` for species the table lacks),
how many species of `--novel_species` the sample holds (`meta_novel_species`), whether the taxon
shares a genus with one of them (`meta_novel_congener`: the taxa their reads land on), and whether a
present species was simulated from its representative, the database's reference, or another
genome (`meta_rep_genome` 1 or 0). With the ranks in `--novel_species` and `--taxonomy` also: the
sample's novel species by the rank they were held out at (`meta_novel_levels`, e.g.
`species:2,family:1`), the deepest rank a taxon shares with a species simulated in the sample
(`meta_relative_rank`; `species` for those species themselves), and, for an absent taxon whose
closest simulated species is a novel one, the rank that one was held out at (`meta_novel_level`):
its reads are the likely source of the taxon's. For a present taxon, `meta_neighbour_rank` is the
deepest rank it shares with another species of its sample: a congener's reads fit it nearly as
well, so it may be missed. The training report breaks its errors down by these.

The genome table should include species that the database lacks. Their reads land on relatives
the database has, and those species are the false positives the model must learn to reject: build
the training database with some species left out (`gtdb_to_protal_db.py --exclude_species`, which
keeps their taxids, so the model applies to the full database) and simulate from all of them.
`build_gtdb_database.py --holdout` does this, and `--holdout-clades` leaves whole families,
classes and phyla out too (below). Include archaea,
too (`--archaea`), since they have fewer marker genes than bacteria. Include genomes other than the
database's references (other strains of its species): real strains differ from the reference by up
to a few percent, and a model that has only seen reads of the reference itself may call them
absent. `scripts/mini_db/simulate_gtdb_release.py --strain_divergence 0.002-0.015
--species_divergence 0.015-0.04` makes a small world with such strains and close relatives.

## Features

The shipped model was trained on absolute counts: genes, k-mers and mates. These depend on the
database, on the domain (archaea have 52 marker genes, bacteria 119), on depth and on read length.
On simulated data it finds about a third of the archaea present, with probabilities pinned near 0.5.

`model_features.py` lists `NORMALIZED_FEATURES`, which do not depend on these:
- fractions of the hittable genes (`hit_gene_fraction`, `gene_presence_ratio`: genes hit against
  the number expected from the number of fragments);
- rates per aligned or covered kb (`lu_per_kb`, `lsu_per_kb`, `variant_sites_per_kb`,
  `multiallelic_sites_per_kb`);
- read identity with indels counted as differences (`identity`, `top_identity`);
- the share of low-identity reads, i.e. reads of relatives (`low_identity_share`);
- the depth's coefficient of variation across genes (`depth_cv`), uniqueness and allele frequency
  classes (`uniqueness`, `RAF0`-`RAF4`, the `*_gene_rate*` columns);
- `fragments` and `depth`, which say how much evidence there is;
- what the reads' other candidates and other records say: `mean_mapq` and `low_mapq_share` (the
  share of the taxon's reads with MAPQ below 10, i.e. nearly as good a second candidate);
  `congener_fit_share` and `other_genus_fit_share` (the share of its reads that another species of
  its genus, or a species of another genus, fits within one edit: the reads of a relative the
  database lacks fit several of its congeners about as well); `linked_share` (the share of its
  reads with two records on it: both mates of a pair, on one gene or two, or two genes of a long
  read). These count every read's best record, also those the profiler's MAPQ filter leaves out.

The dump also has `adjacent_expected_share` and `adjacent_unlikely_share`: of the genes next to
each other on a taxon's reads (a pair's mates on two genes, a long read's consecutive genes within
3 kb of each other on the read; every read's best records, also those the filters leave out), the
shares whose ends face each other in the taxon's clade and that never do there, by the database's
gene neighbours ([running.md](running.md#options-the-website-does-not-list); both 0 without them).
They are not among the normalised features the models train on: reads of a congener the database
lacks pair across the same genes as the species' own, and `linked_share`, which they refine,
lowered the paired-end F1 when it was tried ([report](claude/2026-10-01-gene-neighbours-run/README.md)).

protal writes the alternatives as the `ZA` tag of a read's best record (`ZA:Z:<taxid>:<edits
more>,...`, the other taxa among the read's aligned candidates with at most 5 edits more, or `*`);
a SAM file of an older protal lacks it, protal warns, and the two fit shares are then 0, so such
files are aligned again (`--force`) for a model that uses them. A dump of an older protal lacks the
five columns, and the trainer stops with `--features normalized`.

## Training

    python3 scripts/random_forest_cmdline.py --truth-file training/training_data.tsv \
        --output-prefix training/model

The trainer needs Python 3 with numpy, pandas, joblib and scikit-learn; no Java. `model_pmml.py`
writes the forest as PMML itself, so that protal's probabilities equal scikit-learn's bit for bit:
scikit-learn compares a feature as float32 with its thresholds and cPMML as a double, so each
threshold is written as the largest double that rounds to a float32 at or below it; leaf counts are
written so that cPMML's count over total is scikit-learn's leaf probability; and the trees are
averaged in file order, as protal does. The trainer checks this on every training row and fails if
one differs.

| Option | Default | |
|---|---|---|
| `--truth-file`, `--output-prefix` | required | the training table; the prefix of the outputs |
| `--features` | `normalized` | `normalized`: only `NORMALIZED_FEATURES`; `all`: every feature column of the dump |
| `--reference-pmml` | | train on the input fields of an existing model instead |
| `--ntree`, `--maxnodes`, `--min-samples-leaf`, `--max-features` | 64, 128, 1, `sqrt` | the forest (`--maxnodes 0`: no limit on leaves) |
| `--knob` | 0.5 | the threshold protal will use; calls and their errors are counted at it |
| `--folds` | 5 | folds of the held-out evaluations |
| `--evaluation` | `full` | `basic`: the held-out evaluations only; `none`: fit and export only |
| `--taxonomy` | | the database's `internal_taxonomy.dmp`: domains the table's `meta_domain` lacks, and the lineages for holding out whole clades |
| `--test-file` | | an independent test table (the collector with another design and seed), scored by the fitted forest: the report's section "Independent test set" (metrics, FN and FP rates by depth and by rank, the threshold with the highest F1 there), `<prefix>.test_predictions.tsv.gz`, and a warning when it scores clearly worse than cross-validation |
| `--seed`, `--threads` | 1, 4 | |

Rows of one sample share its reads, and rows of one species share its reference, so a random split
of rows scores a model on samples and species it was trained on. The trainer scores each row with
forests that saw neither its sample ("by sample") nor its species ("by species"), and out of bag
(scikit-learn draws each tree's rows by their class weight, so when absent taxa far outnumber
present ones, some present rows are drawn for every tree and have no out-of-bag score; they are left
out of that estimate and counted in the report). By species is the estimate that matters for a large database: of GTDB's ~130,000 species a
training set holds a few thousand, so most species protal meets in real samples were never in
training. With `--taxonomy` it also holds out whole genera, families, orders, classes and phyla (a forest
that saw no taxon of the row's clade): how far the model carries to parts of the tree the training
data barely cover. The report also compares with the model that profiled the training samples (the dump's
`probability`, e.g. the shipped model), and `--evaluation full` adds studies of whether the data
and settings suffice: the other feature set; the procedure this trainer used before (a grid search
over `max_features`, then a 512-tree forest on only the top features, judged on random rows); the
number of leaves and of trees; and a learning curve with fewer training samples.

| Output | |
|---|---|
| `<prefix>.xml` | the model |
| `<prefix>.report.txt` | the evaluation (also printed): data summary with warnings, held-out results by domain, depth and evidence, the hardest taxa, the threshold, the studies |
| `<prefix>.metrics.json` | the report's numbers |
| `<prefix>.predictions.tsv.gz` | every taxon's probability out of bag, by sample, by species, by rows and by clade held out, with its main features |
| `<prefix>.thresholds.tsv` | precision, sensitivity and F1 by threshold, species held out |
| `<prefix>.varimp.tsv` | feature importances |
| `<prefix>.joblib` | the fitted scikit-learn forest |

On simulated worlds, 64 trees scored as well as 256 or 512 (the forest is 8 times smaller and
loads faster in protal), the leaf limit did not bind, and the grid search, which took most of the
old trainer's time, chose a few top features and did no better on species held out.

### Species and clades the database lacks

Real samples hold organisms the database lacks at every depth: a species of a genus it has, or of a
genus, family, order, class or phylum it has none of. Their reads land on the closest relatives the
database has, or nowhere, and those relatives are the false positives to avoid. `build_gtdb_database.py`
makes such samples: its training database lacks many single species (`--holdout`, 20%) and whole
clades of every rank (`--holdout-clades`, from phylum down to genus), every sample has species of
held-out clades (`--novel-clades-per-sample`), and the absent taxa their reads reach are negatives
the forest is trained on, like any other row.

The report's section "False positives and false negatives by taxonomic rank" gives rates in % at
the knob, for the model scored with species held out and for the collection model, in rows from
species to phylum:

- **False positives from what the training database lacks**: absent taxa whose closest species in
  their sample is one the database lacks (`meta_novel_level`), by the rank it was held out at: the
  novel species simulated, those taxa, how many are called (FP), the FP rate and FP per 100 novel
  species; and the other absent taxa for comparison. A second table splits the same taxa by the rank
  they share with that species (`meta_relative_rank`).
- **False negatives** of present taxa by the deepest rank they share with another species of their
  sample (`meta_neighbour_rank`): whether congeners, or relatives further up, cost sensitivity.
- **With the taxon's clade held out of training** (the cross-validation above, species to phylum):
  FN rate, FP rate, FP per sample and F1 when the forest has seen nothing of that species, genus,
  family, order, class or phylum.

The summary repeats the FP and FN rates by rank.

`check_model_parity.py` re-profiles saved training samples (`--profile_only` on the SAMs a
`collect_training_data.py` folder keeps) with a model and checks that protal's probabilities are
the model file's, and that protal computes the features as it did when the training data was
collected (another protal version may not):

    python3 scripts/check_model_parity.py --db DB --model training/model.xml --training training

`--read_type se` (or `pb`, `ont`) checks a model of other reads, on the collection's samples of that
read type (protal is given it with `--model_se`, `--model_pb` or `--model_ont`).

`gradient_boosted_cmdline.py` and `hist_gradient_boosted_cmdline.py` train gradient-boosted trees
with the same inputs; they still export with sklearn2pmml, which needs Java. The histogram variant
is faster, but its PMML is not guaranteed to load in cPMML: check its scores before you use it. The
R scripts (`random_forest_cmdline.R`, `random_forest.Rmd`) are the older caret pipeline the Python
trainer was modelled on.

## Using a new model

Try it first without changing the database:

    protal --db DB --model training/model.xml --map test.map -t 16

or, cheaper, on existing alignments with `--profile_only`. To make it the database's default,
add it with `protal --add_model training/model.xml --read_type pe --db DB` (`se`, `pb` or `ont` for a
model of those reads), which checks it and replaces the database's model of that read type.

Train the production model on simulations from the database's own genomes (for example GTDB), with
held-out species as negatives. Choose `--knob` on held-out samples like the ones it will profile;
`<prefix>.thresholds.tsv` is a start.
