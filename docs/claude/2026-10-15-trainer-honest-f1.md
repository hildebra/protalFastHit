# The trainer's F1 made honest: marker-truth headline, folds by sample and species, no test-set veto

**Date:** 2026-10-15.

**Code:** branch `f1-trainer`, on `3305b2c`; the changes are in the commit that adds this report:
`scripts/machine_learning_cmdline.py`, `scripts/model_pmml.py`, their tests in `scripts/test_model_pmml.py` and
`scripts/mini_db/test_gtdb_pipeline.py`, [databases.md](../databases.md#the-headline-the-decisions-folds-and-the-features-ranges-since-2026-10-15).

**Why:** the [F1 audit](2026-10-15-f1-audit-since-0.6.0a.md)'s findings E2, E3, E4, C1 and C4 (its section 5, items 3,
5 and 7).

**Data:** the paired-end tables of the `c4test` build (`~/c4test/db/work/{training,test}/training_data.tsv` and their
`.unseen.tsv`, copied to `~/f1trainer/c4`): 8,215 training rows of 200 samples and 524 species, no scenarios; a test
table of 1,304 rows. **Machine:** WSL on the shared laptop, CPUs 0-3 (`taskset -c 0-3 nice -n 5`), Python 3.14.7,
scikit-learn 1.9.1. Times are from a loaded machine and are rough.

## What changed

- **The headline (E3).** The first section of every model's report, and the first key of its `metrics.json`
  (`headline`), is the F1 over every truth species the marker genes can see, as protal calls by default: a present
  species with reads on a marker gene counts also when it has no row (`<table>.unseen.tsv`). It is given held out
  (`held_out`: the decisions' scheme; `held_out_species`: species alone) and on the test set (`test`), each with
  `marker_truth_f1`, `_tp`, `_fp`, `_fn`, `unseen` and the row-based `rows_f1`, `rows_tp`, `rows_fp`, `rows_fn`,
  `rows_scored`. A missing unseen table or marker column is named in `headline.notes`. The summary's other F1s say
  "rows only".
- **Folds that hold out sample and species (E4), `--cv-scheme samples+species`, the default.** Samples and species are
  dealt into `--folds` (k) groups each, balanced by rows, in an order drawn by `--seed`. Fold (i, j) scores the rows of
  sample group i and species group j with a model fitted on the rows of neither. Rotation r pairs i with
  j = (i + r) mod k. `--cv-rotations 2` gives 2k fits that score 2/k of the rows, each once, every fit on
  (1 − 1/k)² = 64% of the rows. This estimate now takes the decisions: the global knob, the knob curve, the calibration
  of `--fdr-calls`, and `--features auto`. Species alone is still computed, reported, and its knob shown beside the new
  one; `--cv-scheme species` restores the old decisions. Unseen species count in the cells the folds score, through
  the table's taxon of their name, or a seeded hash for a species without any row.
- **`--features auto` by the weighted F1.** The candidates are compared by the F1 with each row weighted as the fits
  weigh it (the scenarios' rows 0.25), with "F1 unweighted" beside it.
- **No veto by the test set (E2).** `choose_global_knob` no longer drops a knob that loses on the test set; the test
  set's F1 at the chosen knob, at 0.5 and at the median threshold is reported (`global_knob.test`).
- **Feature versions and ranges in every model (C1, C4).** `model_pmml._write_head` writes, for forests and boosted
  models alike:
  - `<Extension name="protal_feature_ranges" value="name:min:p01:p99:max,..."/>`, every model feature in the model's
    order, over all rows the final model is fitted on (unweighted), numbers as C's `%.9g`;
  - with `--feature-versions FILE` (the table of `protal --feature_versions`: columns feature, version, changed, note,
    a header line), `<Extension name="protal_feature_versions" value="name:version,..."/>` for the model's features
    in the file; the rest are left out and named in the report.

  Feature names must match `[A-Za-z0-9_]+`, versions `[A-Za-z0-9_.+-]+`; the trainer stops otherwise. Read back by
  `read_feature_ranges` and `read_feature_versions`, and checked after every export.

## Measured on c4test

Commands (the build's own trainer options; `$C=~/f1trainer/c4`):

```bash
COMMON="--truth-file $C/training/training_data.tsv --seed 1 --threads 4 --scenario-weight 0.25 \
  --taxonomy ~/c4test/db/work/internal_taxonomy.dmp --evaluation basic --depth-knobs \
  --test-file $C/test/training_data.tsv"
FOREST="--model forest --ntree 64 --maxnodes 512"; GBM="--model gbm --maxnodes 63"
SET=normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood+ancestry+gaps+untried+alleles+polymorphic+weights+strain_identity+length_free+read_length
python <HEAD's machine_learning_cmdline.py> $COMMON $FOREST --features $SET --output-prefix $C/out/old
python scripts/machine_learning_cmdline.py $COMMON $FOREST --features $SET --cv-scheme species --output-prefix $C/out/new_species
python scripts/machine_learning_cmdline.py $COMMON $FOREST --features $SET --output-prefix $C/out/new_crossed
# the same three with --features auto; and $GBM with --cv-scheme species and without
```

**Cost.** The crossed folds: 5 groups each, 2 rotations, 10 fits, 3,276 of 8,215 rows scored (39.9%). Species held
out: 5 fits, every row. The forest's evaluation took 59 s against 45 s before. `--features auto` took 605 s against
310 s (22 candidates × 10 fits instead of × 5).

**HEAD and `--cv-scheme species` agree exactly:** the same evaluation, knob 0.47 and test F1, so the fallback is the
old procedure. The old test veto did not act on this table: the knob gained on the test set.

**Forest (the build's default model):**

| | species held out (old) | samples and species held out (new) |
|---|---|---|
| rows scored | 8,215 | 3,276 |
| F1 at 0.5, rows only | 0.9730 | 0.9731 |
| on the same 3,276 rows: F1 at 0.5, AP, log loss | 0.9734, 0.9962, 0.0847 | 0.9731, 0.9940, 0.1049 |
| best threshold on those rows | 0.47 | 0.58 |
| global knob (median, support) | 0.47 (0.47, 95.0%) | none: 0.5 (0.57, 63.5%) |
| test F1 at the model's knob, rows only | 0.9603 | 0.9588 |
| test marker-truth F1 (4 unseen) | 0.9685 | 0.9670 |
| held-out marker-truth F1 | 0.9775 (at 0.47) | 0.9775 (0.9768 by species at 0.5) |

`--features auto` chose the default set under both schemes, by an F1 of 0.9730 (species) and 0.9731 (crossed)
against a mean of 0.9715 for the other 21 sets.

**Gradient boosting:**

| | species held out (old) | samples and species held out (new) |
|---|---|---|
| on the same 3,276 rows: F1 at 0.5, AP | 0.9752, 0.9970 | 0.9723, 0.9965 |
| global knob (median, support) | none: 0.5 (0.64, 87%) | 0.63 (0.63, 97%, +0.0016) |
| test F1 at the model's knob, rows only (FP, FN) | 0.9639 (15, 36) | 0.9659 (10, 38) |
| test marker-truth F1 | 0.9722 | 0.9742 |

- **What it shows:** on the same rows, holding the sample out costs the boosted model 0.0029 of F1 and the forest
  0.0003. That is the per-sample learning E4 describes, and the v13 report saw it too: the forest exploits it less.
  Part of the drop is the smaller training share (64% against 80%).
- **The knob moved under both models,** once to a value that won on the test set (gbm 0.63: +0.002), once back to
  0.5 (forest: −0.0015 against 0.47). Both differences are within the ±0.002 that earlier reports called noise.
- **Without scenarios,** this table does not show the shallow-soil effect that v13 measured (0.881 against 0.931).
  The scenario hold-outs of an r226 build are where the new scheme should matter most.

**cPMML ignores the extensions:** `check_model_parity.py` with this branch's protal on the training database of
`c4test` (80 samples of `rl150_p10000` and `rl150_p50000`, 2,630 taxa):
- the forest and the boosted model with `protal_feature_ranges`, and the forest with both extensions: probabilities
  equal the model file's for every taxon;
- the features equal the training data's;
- the trees of the model with versions are byte for byte those of the model without.

## Checks

`scripts/test_model_pmml.py`, 58 tests, all pass, among them these new ones:
- `FeatureExtensionsTest`: the ranges' `%.9g` text, bad names, versions and numbers, the versions table read, and
  placeholders without extensions.
- `ForestExportTest.test_feature_versions_and_ranges`; `BoostedExportTest.test_header_extensions` extended.
- `CrossedFoldsTest`: test and training rows share no sample or species; training is every row of neither group;
  rotations score r/k of the rows, each once; the groups are balanced and drawn by the seed; unseen species fall in
  the cells as their rows would.
- `TrainerHeadlineTest`:
  - the headline is first in `report.txt` and `metrics.json`;
  - marker truth = the rows' calls plus the unseen species in the scored cells as FN;
  - the note for a missing test unseen table;
  - `--cv-scheme species`, `--feature-versions`, and a bad versions file stopping the run.
- `TrainerDepthKnobsTest.test_one_knob_for_every_sample`: the knob stands against a test set it loses on.

`scripts/mini_db/test_gtdb_pipeline.py` checks the headline and the ranges of a GTDB build's models. With this branch's
protal and `simulate_metagenomes`: all 8 of its tests pass (1,896 s at load 13-15), including the build's parity checks
of the new models, and so do the 39 of `test_gtdb_build.py`.

## Left open

- **A crossed estimate at r226,** on tables with scenarios. This is the measurement that tells whether the new
  scheme changes builds.
- **The build's summaries** (`build_gtdb_database.py`, edited by another session) still print the species-held-out
  row F1. They should print `headline.held_out.marker_truth_f1` and `headline.test.marker_truth_f1`.
- **protal itself** (C++) does not yet read `protal_feature_versions` or `protal_feature_ranges`.
- **The website** is not affected (trainer only).
