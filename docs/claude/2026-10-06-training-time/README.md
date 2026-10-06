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
1 fit + 5 (species) + 35 = ~41 fits for the paired-end model against ~220, about a fifth (an estimate from the
fit counts; the measurement below). `--features auto` and `--evaluation full` remain for audits. The trainer
also takes auto's fit of the chosen set as its final model, and the learning curve's all-samples point from the
evaluation's species folds (the same models; 6 fits fewer).

## Still to measure (running on 2026-10-06)

- The trainer with the build's new defaults on the same table and machine (`FEATURES=<default set> EVAL=basic`).
- One fit's time at 1, 2, 3 and 6 threads, and cheaper settings' F1 (250 rounds at 0.1, 150 at 0.15, 63 bins)
  with species held out and on the test sets ([`gbm_speed.py`](../2026-10-06-r226-v13-soil/gbm_speed.py)).

Further cuts, not made: the clade ranks (25 of the ~41 fits, audit only) into `--evaluation full`; the folds
fitted side by side in processes rather than one after another with all threads (worth it if one fit scales
poorly with threads; the measurement above says).
