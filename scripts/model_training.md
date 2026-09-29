# Training protal's presence model

protal decides whether a species is present with a random forest (PMML) from the database. It
scores every species with reads, and reports those whose probability of `TRUE` is at least
`--knob` (0 to 1, default 0.5).

A database holds one model per read type, and `--read_type` picks it: `model_pe.xml` (paired-end,
the default; older databases: `model.xml`), `model_se.xml` (single-end reads < 500 bp),
`model_PB.xml` (PacBio), `model_ONT.xml` (Oxford Nanopore). Train each on simulated reads of its
type, and store it in the database with

    protal --add_model training/model.xml --read_type se --db DB -t 8

which checks the model as below and replaces the one stored for that read type. Before `--build`,
the same names in the database folder go into `database.protal` too.

## What protal expects of a model

At start, before any alignment, protal checks that the model
- takes only inputs protal computes: the features of `TaxonFeatures` in
  `src/Profiling/Profiler.h`, the columns of the training dump below;
- predicts a field with the value `TRUE`, whose probability is the species' score;
- scores a species.

A model that fails one of these stops the run with exit code 2, naming the problem.

## Training data

With a truth file (`--profile_truth`, or a `PROFILE_TRUTH` column in the map), protal writes
`<profile>.truth_annotated`. This file has one row per species with reads: whether the species was
present (`truth`), the current model's call and probability, and every feature. The features are
written exactly as the model is given them.

`collect_training_data.py` produces such data at scale. It simulates metagenomes with
`simulate_metagenomes` over a grid of read setups (length, ART profile, fragment size) and depths,
profiles them against a database, and joins all dumps into `training_data.tsv`:

    python3 scripts/collect_training_data.py --db DB --genome_table genomes.tsv -o training \
        --archaea 2 --species_per_sample 10-40 -t 8

The genome table should include species that the database lacks. Their reads land on relatives
the database has, and those species are the false positives the model must learn to reject.
Include archaea, too (`--archaea`), since they have fewer marker genes than bacteria.

## Features

The shipped model was trained on absolute counts: genes, k-mers and mates. These depend on the
database, on the domain (archaea have 52 marker genes, bacteria 119), on depth and on read length.
On simulated data it finds about a third of the archaea present, with probabilities pinned near 0.5.

`model_features.py` lists `NORMALIZED_FEATURES`, which do not depend on these:
- fractions of the hittable genes;
- the gene presence expected from the number of fragments;
- rates per aligned or covered kb;
- read identity with indels counted as differences;
- the share of low-identity reads (reads of relatives).

Train on them with `--features normalized`:

    python3 scripts/random_forest_cmdline.py --truth-file training/training_data.tsv \
        --output-prefix training/model --features normalized

Train the production model on simulations from the database's own genomes (for example GTDB), with
held-out species as negatives. Choose `--knob` on held-out samples like the ones it will profile.
