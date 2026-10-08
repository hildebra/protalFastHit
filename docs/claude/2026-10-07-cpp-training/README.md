# Training the presence model in C++: possible, and how much faster?

**Question.** The build trains the presence models (gradient boosting by default, or a random forest) in Python
(`scripts/machine_learning_cmdline.py`, scikit-learn). Could protal's C++ do the training as well, and would it be
faster?

**Short answer.** protal's C++ has no training code to extend: it only scores models (cPMML). Training in C++
would mean linking a tree library such as LightGBM or XGBoost, or writing one. Those libraries fit protal's boosted
model **1.6-1.8 times faster** than scikit-learn on the r226 v14 paired-end table (five folds with species held out,
4 threads), with identical scores. That would save about 30% of the trainer's time. Fitting is two thirds of the
trainer, and scikit-learn already does it in compiled code (Cython and OpenMP); Python only drives it. In a build
this is a minute or two: the four models train side by side in a few minutes of a build that takes hours. The same
gain is available without moving training into C++: the trainer could call LightGBM from Python. Porting the
trainer's 2,100 lines of evaluation into protal is not worth that gain.

## Data and commands

- **Code:** `9ecf7ee` (`audit-fixes`; the trainer as of `652ed53`, defaults of `15006c2`: 250 rounds at a learning rate
  of 0.1, at most 63 leaves, at least 20 rows per leaf, L2 1, classes balanced, scenario rows at weight 0.25).
- **Data:** the r226 v14 build's paired-end tables (`local/v14/training/training_data.tsv`: 369,803 rows from 333
  samples; `local/v14/test/training_data.tsv`), copied to the WSL file system. The features are the v14 build's set
  (`normalized+adjacency+distance+depth+divergence+unfiltered+ref`, 52 features).
- **Machine:** WSL on the Core Ultra 7 258V, every run pinned to cores 0-3 (`taskset -c 0-3 nice -n 5`). The machine
  is shared with other sessions; their load is noted where it changed a result. Software: Python 3.14,
  scikit-learn 1.9.1, numpy 2.5.3, and in a separate venv (the build environment is untouched) LightGBM 4.7.0 and
  XGBoost 3.4.1 (`xgboost-cpu`) from PyPI.
- **Scripts:** [`prof_trainer.sh`](prof_trainer.sh) (the trainer under cProfile, then without it),
  [`profile_breakdown.py`](profile_breakdown.py) (the profile by callee), [`bench_libraries.py`](bench_libraries.py)
  (the same table, features, row weights and settings fitted in scikit-learn, LightGBM and XGBoost),
  [`bench.sh`](bench.sh) and [`bench2.sh`](bench2.sh) (the benchmark runs).

## What the C++ side has

- protal scores the model with **cPMML** (`lib/cPMML`), a library that only scores PMML ("High-Performance PMML
  Scoring"; it has no training code).
- `src/RandomForest/` is the training-data **simulator** (`MetagenomeSimulator`, `CommunityProfileDesigner`,
  `ArtIlluminaWrapper`), despite its name. No forest or boosting is built there.
- So "the C++ implementation trains" means new code in protal: a tree library linked in (LightGBM is MIT-licensed,
  XGBoost Apache-2.0; both use OpenMP, as protal already does), or a histogram-boosting fit written from scratch.

## Where the trainer's time goes

The trainer with the build's defaults (`--evaluation basic --depth-knobs`, the test table), on 4 threads (folds one
after another, which is what `--fold-jobs` defaults to at 4 threads), under cProfile
([`trainer_profiled.log`](trainer_profiled.log), [`profile_breakdown.txt`](profile_breakdown.txt)):
290.5 s, 905 CPU-seconds.

| part | seconds | share | runs in |
|---|---|---|---|
| 16 boosting fits (1 final, 15 folds: rows, samples, species) | 196 | 67% | scikit-learn's Cython/OpenMP code, driven node by node by its Python grower |
| scikit-learn's predictions (folds, test set, all rows) | 14.5 | 5% | compiled |
| reading the two tables (`pandas.read_csv`, `round_trip`) | 23 | 8% | pandas' C parser |
| writing `predictions`/`calls`/`test_predictions` `.tsv.gz` (`to_csv`, gzip) | 26 | 9% | pandas formatting, zlib |
| the PMML parity check (`model_pmml` scores all 369,803 rows in numpy) | 20 | 7% | numpy |
| metrics, report, scenarios and the rest | 11 | 4% | Python, pandas |

(Run again without the profiler, while another session's `rsync` shared the cores: 400 s,
[`trainer_unprofiled.log`](trainer_unprofiled.log). Both runs: PMML and scikit-learn identical on every row.)

The fits dominate, and scikit-learn already runs them in compiled code. A C++ trainer can only beat it with a faster
fitting algorithm, or by dropping the Python loop that scikit-learn's grower runs per node. The benchmark below
measures both.

## The fit in scikit-learn, LightGBM and XGBoost

Same rows, features, weights and settings (LightGBM takes each by name; XGBoost has no minimum of rows per leaf and
keeps its `min_child_weight` of 1). [`fit_times.tsv`](fit_times.tsv): one fit on all 369,803 rows, each twice with
the libraries' order reversed:

| threads | scikit-learn | LightGBM | XGBoost |
|---|---|---|---|
| 1 | 29.0-30.2 s | 13.4-18.8 s | 17.2-20.4 s |
| 2 | 15.5-16.8 s | 9.2-9.8 s | 10.2-11.1 s |
| 4 | 11.7-13.8 s | 8.7-8.9 s | 10.4 s (19.5 s once) |
| training log loss | 0.02874 | 0.02871 | 0.02874 |

What the trainer actually does is five folds with species held out ([`folds1.tsv`](folds1.tsv),
[`folds2.tsv`](folds2.tsv); 4 threads, fit and predict, two runs, quiet machine):

| | scikit-learn | LightGBM | XGBoost |
|---|---|---|---|
| 5 folds | 48.8 s, 54.4 s | 31.7 s, 30.5 s | 29.4 s, 31.9 s |
| F1 at 0.5 | 0.9520 | 0.9521 | 0.9520 |
| AP / AUC | 0.9917 / 0.9977 | 0.9917 / 0.9977 | 0.9917 / 0.9977 |
| log loss | 0.0472 | 0.0471 | 0.0472 |

scikit-learn's row is exactly the trainer's own "species" row (F1 0.9520, AP 0.9917, AUC 0.9977, log loss 0.0472),
so the benchmark reproduces the trainer.

**Where the difference comes from.** On 4,000 rows drawn from the table, the 250 trees cost almost nothing per row,
so what remains is each library's fixed cost per tree and node ([`small.tsv`](small.tsv), quiet machine; a first
run during another session's tests, [`small_contended.tsv`](small_contended.tsv), had scikit-learn at 9-12 s on 4
threads):

| threads | scikit-learn | LightGBM | XGBoost |
|---|---|---|---|
| 1 | 1.8-1.9 s | 0.31 s | 0.35-0.38 s |
| 2 | 2.2-2.3 s | 0.25 s | 0.27-0.31 s |
| 4 | 3.8-4.4 s | 0.23 s | 0.24 s |

scikit-learn grows each tree from Python, node by node, with an OpenMP region per histogram. That costs ~1.5 s per
fit on one thread, and the cost grows with the threads (~4 s on 4). This is why one scikit-learn fit uses a few
threads well and many poorly (docs/claude/2026-10-06-training-time) and why the trainer fits folds side by side.
The C++ libraries have almost no fixed cost. Their histograms are also faster: on one thread the full fit is 1.6-2.2
times faster, beyond the fixed part. **Caveat:** LightGBM's OpenMP threads spin while they wait, so it degrades badly
when other work shares its cores. Its first run of the five folds fell on another session's tests on the same cores
and took 342 s against 69 s for scikit-learn ([`quality_first_run.tsv`](quality_first_run.tsv)). Rerun on a quiet
machine, it took 31 s. On a shared node it needs cores of its own.

**The forest** (`--model forest`: 64 trees of at most 256 leaves, √52 features per split; [`forest.tsv`](forest.tsv)):
five folds with species held out took 88.5 s in scikit-learn, 15.0 s in LightGBM and 15.5 s in XGBoost (5.7-5.9 times
faster), with the same F1 (0.9431 / 0.9436 / 0.9438) and AP. scikit-learn's forest searches exact splits, whereas the
C++ libraries use 255-bin histograms. Their probabilities differ, though: log loss 0.168 against 0.063, because a
histogram forest's leaf is a Newton step, not a share of the classes. The knob and the calibration would have to be
chosen afresh. Boosting is the default since 2026-10-06, so this gain is for a model the build no longer trains.

## What training in protal (C++) would take, and what it would gain

- **The gain.** With fits at LightGBM's speed (×0.6), the trainer above goes from ~290 s to ~210 s. C++ would also
  read and write the tables and score the parity check several times faster. Those parts (~70 s here) can be cut in
  Python too (below). A C++ trainer would end near ~150-170 s on this table, about what Python with LightGBM and
  those fixes would reach. Its fits would be the same library's.
- **In a build.** At r226 v14 (cluster, 16 threads per trainer, the older defaults: 500 rounds, clade ranks in
  `basic`), the four trainers ran side by side and took 7-15 min (pe 413 s, se 892 s, pb 664 s, ont 716 s; their
  `classifier_training*.log`). Since `15006c2` a trainer fits 16 times instead of 41, each with half the trees, so the
  next build's training should take a few minutes (an estimate; the next build's logs will tell). Fits 1.6-1.8 times
  faster save one or two of those minutes, in a build of several hours (simulation, the training database, protal on
  ~1,000 samples).
- **The cost of the port.** Linking LightGBM into protal is the small part. The trainer is ~2,100 lines (66
  functions) on pandas and scikit-learn: grouped folds by rows, samples, species and clades; AP, AUC, log loss and
  thresholds; the knob and its depth curve; isotonic calibration (`--fdr-calls`); scenario hold-outs; the feature-set
  choice; the studies; the reports, `metrics.json`, `calls.tsv.gz`; and the PMML writer that the parity check and
  41 tests in `test_model_pmml.py` guard. All of it would be rewritten in C++ and tested again. protal's binary
  would carry a training library that only the build uses.

## Recommendations

1. **Do not move training into protal's C++.** The gain is ~30% of a trainer that now takes minutes, and the port
   would duplicate the trainer's evaluation, which keeps changing (scenarios, features, knobs).
2. **If training time ever matters, swap the estimator, not the language:** a `--model lgbm` in the trainer, with
   LightGBM's C++ fit called from Python. It scores the same here and fits 1.6-1.8 times faster. It needs:
   - a PMML writer for LightGBM's trees in `model_pmml.py`, ~60 lines. `Booster.dump_model()` gives each split's
     feature, `<=` threshold, the default direction and the leaf values, and LightGBM folds its starting score into
     the first tree, as `write_boosted` does: the same `modelChain` of regression trees into a logit (checked on a
     toy model: the dumped leaves, summed per row, give LightGBM's raw score to within 9e-16);
   - the parity check on every row, unchanged;
   - `lightgbm` in `envs/protal-db-build.yaml` and in CI;
   - threads matched to the cores a trainer really has. Fold jobs already do this; LightGBM needs it more than
     scikit-learn does.
3. **Cheaper cuts in the Python that is there** (~70 s of 290 s here, no new dependency):
   - writing the gz tables (26 s): a lower gzip level, or zstd through `scripts/compressed.py`;
   - reading the two tables (23 s): a faster parser, if it reads back the same doubles as `round_trip`;
   - the PMML check of every row (20 s): `model_pmml`'s numpy traversal could be faster. The check itself is worth
     keeping: it is what proves that protal scores what was trained.

None of this changes protal's behaviour or options, so neither `docs/` nor the website needs an update.
