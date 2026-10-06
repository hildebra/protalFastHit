# Training time with gradient boosting: where it goes, and the build's new defaults

**Question.** With gradient-boosted models as the default (`1750475`), how long does training take, does
`--features auto` cost much, and how can the build train faster?

**Data and commands.** The r226 v13 tables (`local/v13`: the paired-end training table, 224,011 rows of
321 design and 12 scenario samples; the test table; `internal_taxonomy.dmp`), the trainer of `1750475`
(`scripts/machine_learning_cmdline.py`) as the build ran it then, on 6 threads in WSL
([`time_trainer.sh`](time_trainer.sh)):

    READ_TYPE=pe FEATURES=auto EVAL=full bash time_trainer.sh
    # = machine_learning_cmdline.py --truth-file v13/training/training_data.tsv --test-file v13/test/training_data.tsv
    #     --taxonomy v13/internal_taxonomy.dmp --features auto --evaluation full --depth-knobs --threads 6

The machine (6 cores) was shared with other sessions' builds, tests and trainings throughout (load average
7-14), so wall-clock times are inflated, by different amounts in different stages; CPU time is the steadier
measure. The report: [`pe-auto-full.report.txt`](pe-auto-full.report.txt), its timing:
[`pe-auto-full.time.txt`](pe-auto-full.time.txt).

## Where the time goes

The trainer's work is model fits; each fold of a cross-validation is one on four fifths of the rows (5 folds).
By stage, with `--features auto` (9 candidate sets) and `--evaluation full`:

| stage | fits | wall time here |
|---|---|---|
| feature choice (`--features auto`): 9 sets x (5 folds + 1 fit on all rows for the test sets) | 54 | 50 min |
| final model | 1 | 25 s |
| evaluation (`basic`): rows, samples (5 each), species (reused from the choice), 5 clade ranks (25) | 35 | 63 min |
| study "Feature sets" (`full`): 12 sets by rows, 4 more by species | 80 | 44 min |
| study "Model size" (`full`): 5 sizes x 5 folds, rounds from 5 fits of 1,000 | 35 (~40 fit-equivalents) | 22 min |
| study "More training samples?" (`full`): a quarter and half of the samples (5 folds each), all | 15 | 27 min |
| knob, test set, scenarios, export | 0 | 0.5 min |
| **total** | **~220** | **3 h 28 min** (11.7 CPU-hours; 347% of the 6 cores) |

One fit on all 224,011 rows took 25 s (52 features, 500 rounds of 63 leaves); the stages' per-fit times here
(35-110 s) carry the machine's other load.

- **Does `--features auto` cost much?** With `--evaluation full` hardly: the "Feature sets" study scores every
  set with species held out anyway, and reuses the choice's scores, so the choice adds its 9 fits on all rows
  (~4%). With `--evaluation basic` it doubles the training: 54 fits beside the ~41 the rest needs.
- **What does the build read of it?** `summary.txt`, `build_metadata.tsv` and the console take the scores with
  species held out, the test set, the knob and the scenarios' rows: the evaluation with species (and, for the
  scenarios' table, samples) held out. The studies and the clade ranks are for audits.
- **What auto chose:** exactly the default set with the reference's uniqueness
  (`normalized+adjacency+distance+depth+divergence+unfiltered+ref`, F1 0.9554 with species held out; the relatives
  variant 0.9559, within 0.002), as every model of the v12 and v13 builds chose the default set.

## What the studies said about the model

- **Size** (species held out): 500 rounds of 63 leaves has the lowest log loss (0.0459); 250 rounds 0.0487
  (F1 -0.003), 100 rounds 0.0567; 750 and 1,000 rounds and 127 leaves leave the log loss where it is (0.0455-0.0462)
  and only shift the calls at 0.5 from sensitivity to precision. Halving the rounds would halve a fit's time at a
  cost; the settings stay.
- **More samples of the same design** do not help: log loss 0.0502 with a quarter of the samples, 0.0466 with
  half, 0.0459 with all.

## The build's defaults since `14632a8`

`build_gtdb_database.py` trains on the default set (`--features`, the set auto chose) with `--evaluation basic`:
1 fit + 5 (species) + 35 = ~41 fits for the paired-end model against ~220, about a fifth of the fits (measured
below: 2.6 times less CPU time on the shared machine). `--features auto` and `--evaluation full` remain for audits. The trainer
also takes auto's fit of the chosen set as its final model, and the learning curve's all-samples point from the
evaluation's species folds (the same models; 6 fits fewer).

## The new defaults, measured

The trainer of `14632a8` with the build's new defaults on the same table and machine
(`FEATURES=normalized+adjacency+distance+depth+divergence+unfiltered+ref EVAL=basic bash time_trainer.sh`;
[`pe-default-basic.report.txt`](pe-default-basic.report.txt), [`pe-default-basic.time.txt`](pe-default-basic.time.txt);
the load during it in [`load_during_new_defaults.txt`](load_during_new_defaults.txt), 7-13):

| | old defaults (auto, full) | new defaults (default set, basic) |
|---|---|---|
| fits | ~220 | 41 |
| wall time (6 shared cores) | 3 h 28 min | 1 h 17 min |
| CPU time | 11.7 h | 4.4 h |

2.6 times less CPU time, not the 5 times the fit counts promise: the machine's other load inflated this run's
fits more (its final fit took 55 s, the earlier run's 25 s; under oversubscription OpenMP threads that wait burn
CPU time too). On a node of its own the ratio should approach the fit counts'.

## One fit's threads, and cheaper settings

[`gbm_speed.py`](../2026-10-06-r226-v13-soil/gbm_speed.py) ([`speed_pe.tsv`](speed_pe.tsv)): one fit on four fifths of
the paired-end rows (the default set, the sample's depth rounded to 0.25 as varied scenario depths will make it,
scenario rows at 0.25), and settings scored with species held out (5 folds) and on the test table:

| threads | 1 | 2 | 3 | 6 |
|---|---|---|---|---|
| one fit (s) | 88 | 52 | 52 | 236 |

| setting | F1 species held out | design test | soil hold-out | shallow soil hold-out |
|---|---|---|---|---|
| 500 rounds at 0.05 (default) | 0.9542 | 0.9691 | 0.9554 | 0.9265 |
| 250 rounds at 0.1 | 0.9537 | 0.9693 | 0.9551 | 0.9257 |
| 150 rounds at 0.15 | 0.9529 | 0.9673 | 0.9549 | 0.9242 |
| 500 rounds, 63 bins | 0.9532 | 0.9674 | 0.9545 | 0.9260 |

- **Threads.** Two threads were 1.7 times as fast as one; three no faster than two, and six four times slower
  than three, the machine's other load (load ~9) leaving boosting's OpenMP threads waiting on each other. How far
  one fit scales on a node of its own is not measured here; the shape says that a fit uses a few threads well and
  many poorly.
- **Fewer rounds at a higher rate.** 250 rounds at a rate of 0.1 score as the default within 0.001 everywhere
  (one seed: -0.0005 with species held out, +0.0002 on the design's test set, -0.0003 soil, -0.0008 shallow soil)
  at half the trees, so half the time of every fit and of the model's scoring in protal. (At the default rate, the
  "Model size" study's 250 rounds were worse: the rate makes up for the rounds.) 150 rounds at 0.15 and 63 bins
  lose 0.001-0.002.

The times of the settings' rows (96-301 s per fold) follow the machine's load, not the settings.

## Further cuts, not made

In order of what they are worth:
1. **250 rounds at a rate of 0.1** (`--rounds 250 --learning-rate 0.1` as the defaults): half of every fit, F1
   within 0.001 (one seed here; worth a second seed or a build to confirm).
2. **The folds side by side**: the 5 folds of a cross-validation fitted at once in processes of a few threads each
   (5 x 3 of a trainer's 16) instead of one after another on all 16, as one fit uses few threads well. Likely 2-4
   times faster on the cluster (an estimate from the shape above); the scores stay the same.
3. **The clade ranks into `--evaluation full`**: 25 of the 41 fits are the five clade ranks, which the summary does
   not read; at r226 v13 they scored within 0.0004 of each other (F1 0.9474-0.9478).

Together, with the new defaults: 16 fits of 250 rounds, five at a time, a few minutes per model instead of the
hours of `--features auto` and `--evaluation full` (an estimate).
