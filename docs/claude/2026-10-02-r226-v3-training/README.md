# The second and third GTDB r226 builds: is the training sound, and did it improve?

2026-10-02. Data: the training logs of two SLURM runs of `build_gtdb_database.py` on GTDB r226 with the
tuned design of `460c9a1` (scripts and binary between `460c9a1` and `9fd81c8`; their `build_metadata.tsv`
was not written): v2 (job 23884920, 16 threads) and v3 (job 23887614, 64 threads), copied from
`OUTDIR/{convert,training_db,training_db_index,index_and_package}.log`, `classifier_training*.log` and
`model_logs/holdout.txt` into `local/protal_r226_v{2,3}_logs.tgz` (git-ignored; extracted to
`local/protal0.7.3_r226_v{2,3}/`). Compared with the first r226 run, v1 (job 23865669, scripts `ff57266`,
[2026-10-02-r226-build-evaluation](../2026-10-02-r226-build-evaluation/README.md), `local/protal0.7_r226_v1/`).
Numbers up to the tuning suggestions are from the reports in those logs (`sed`/`grep` of their sections). The
sections after them use a second upload of v3's OUTDIR, `local/protal_r226_v3_training.tgz` (predictions,
thresholds, parity, simulation logs and the training and test tables; extracted to the same folder), and the
scripts in this folder: `retrain.sh` (the four models trained again on v3's tables with suggestions 1 and 4),
`knob_noise.py` and `strains_missed.py` (with their outputs as `*_output.txt`).

## Verdict

The training is sound and the models are usable. Both runs stopped at step 7/8, the parity check, but
not because of the models: v2 failed on the pb model and v3 on the pe model, with byte-identical
training reports (the four `classifier_training*.log` differ only in timings and paths). The check
failed on `adjacent_support` differing in its last digits (about 1e-15) between a sample profiled
among others and alone on more threads, which `9fd81c8` fixed (an order-independent sum, and a
tolerance of 1e-12 relative in `check_model_parity.py`).

## Performance against v1

F1 at knob 0.5 with species held out of training (cross-validation), on the independent test set, and on
the test set with the model's knob curve over sample depth, as protal calls by default (v1 had no knobs):

| read type | v1 held out | v1 test | v3 held out | v3 test | v3 test, knob curve |
|---|---|---|---|---|---|
| pe | 0.952 | 0.948 | 0.958 | 0.953 | **0.962** (FP 305 → 159, FN 162 → 212) |
| se | 0.947 | 0.947 | 0.953 | 0.949 | **0.957** (FP 326 → 146, FN 176 → 259) |
| pb | 0.917 | 0.886 | 0.972 | 0.975 | 0.974 (FP 23 → 22, FN 36 → 39) |
| ont | 0.925 | 0.904 | 0.973 | 0.970 | 0.972 (FP 40 → 34, FN 31 → 32) |

- **Short reads: a real, modest gain.** v3's test set adds deep samples (1M pairs 12, 5M pairs 6), which
  are easier, so compare the depths both test sets have (500 to 1M pairs, 72 samples each): pe F1 0.948 →
  0.953 at knob 0.5 (FN 191 → 161, FP 283 → 261), se 0.947 → 0.949. FN rates fell at every depth (pe at
  500 pairs 22.7% → 4.3%, 200k 4.0% → 3.4%, 1M 1.42% → 1.03%), FP rates stayed similar. The knob curve adds
  +0.009 (pe, se) out of sample, trading 50 (pe) and 83 (se) more misses for 146 and 180 fewer false
  positives.
- **False positives no longer grow with depth.** v1 extrapolated about 130 per sample at 10M pairs; v3's
  training samples at 10M pairs have 6.7 at knob 0.5 (40 of 12,236 absent taxa in 6 samples, FP rate
  0.33%), at 2M 9.9, at 500k 6.3; with the curve (knobs 0.78-0.94 there) fewer. The deep design points
  replaced v1's extrapolation with data.
- **Long reads: not comparable.** v1's long-read data came from the old pbsim per-contig collector
  (reads of ~100 bp at the shallow points), so its 0.886 / 0.904 measured broken samples. v3's 0.975 / 0.970
  are the first numbers on realistic long reads (HiFi made by the collector, ONT by `pbsim --strategy templ`).
- **What the model uses:** the 0.7.2 divergence features lead the importances (pe: `excess_high_share`
  0.151, `excess_median` 0.149, `identity` 0.116; pb, ont: `excess_median` first). The previous procedure
  (grid search, top features, 512 trees) is again clearly worse: pe 0.933, se 0.923 against 0.958, 0.953.
- **Reproducible:** v2 (16 threads) and v3 (64 threads) gave the same trees and numbers.

## Where the errors are (v3, pe, species held out, knob 0.5)

- **Misses:** 292 of 355 are strains (another genome than the representative; FN rate 4.69% against
  1.16% of representatives), 243 of them with 1-10 fragments. Missed strains' identity median 0.968 (found:
  0.986), as in v1. Long reads miss a few strains even with more than 100 fragments (pb 10 of 419, ont 8 of
  408; representatives 0): their reads differ by ~3% beyond the base qualities, which is what a missing
  congener's reads look like.
- **False positives:** 381 of 642 are taxa beside a species held out of the training database (FP rate
  1.78%), 32 beside held-out clades, 229 beside species it has (1.23%; v1 2.25%). The highest-scored ones have `top_identity` 1.0 and a few
  fragments: species whose marker genes are nearly identical to a present one's.
- **Clades held out** (genus to phylum): F1 0.951-0.952, so the model holds on parts of the tree the
  training lacks.

## Tuning suggestions

1. **Leaves per read type.** 512 leaves beat 256 for short reads (pe log loss 0.0573 → 0.0538, FP 642 → 582,
   F1 0.9578 → 0.9592; se 0.0708 → 0.0679, F1 0.9527 → 0.9540). For long reads 256 does not bind (trees grow
   to ~210-256 leaves), and 128 leaves calibrate better (pb log loss 0.0986 → 0.0872, ont 0.0875 → 0.0855) at
   the same F1. So 512 for pe and se, 128 for pb and ont (a per-read-type `--maxnodes`).
2. **More long-read samples.** The ONT learning curve is still falling (log loss 0.116, 0.098, 0.088 at a
   quarter, half and all of 126 samples), PacBio's nearly so (0.121, 0.105, 0.099). `--long-read-samples 36`
   (the most communities a depth index has: 12 samples of 3 setups) gives 186 per type, 1.5 times the
   samples at the five shallower long-read points (the two deep ones keep their 4 and 2).
3. **A larger long-read test set.** 22 samples (2,108 and 2,507 taxa) cannot tell the knob curve from 0.5
   (pb −0.001, ont +0.002): 8 samples per long-read test point would.
4. **Steadier knob ends.** The deepest point of the pe, se and ont curves rests on 6, 6 and 2 window samples
   and falls below its neighbour (pe 0.94 → 0.78, se 0.92 → 0.75, ont 0.86 → 0.55); protal keeps it for any
   deeper sample. A point should need a minimum number of samples in its window (e.g. 10), its bin merged
   with the next otherwise; or the deepest training point gets more samples (`10000000:4`).
5. **Not worth it:** more short-read samples (pe and se learning curves are flat), all 78 dump columns
   (+0.001 F1, with counts that do not carry across databases), the previous-procedure study (up to 45 s per
   model; worse every time, `--previous-procedure` can stay off).

## Follow-up: suggestions 1-4 implemented

On `6d2e241`: `build_gtdb_database.py --maxnodes` takes `N` and `TYPE:N` items (default `512,pb:128,ont:128`,
recorded per model as `classifier_max_leaves`), `--long-read-samples` defaults to 36, the new
`--test-long-read-samples` to 8, and a knob point needs 6 samples in its window (`DEPTH_KNOB_MIN_SAMPLES` in
`random_forest_cmdline.py`; `depth_knob_windows`: a bin with fewer joins the next deeper one, bins left at the
deep end join the point before; the report's knob table gives each point's samples' depth range).

The minimum is 6, not the 10 suggested above, because of how noisy the knobs are (next section): the pe and se
10M-pair point rests on 6 samples, but its knob is well determined (bootstrap 0.68-0.87, where the 2M point
is 0.93-0.95) and folding it into the 2M point cost 0.002-0.007 F1 on those samples. The long-read deep points
of 4 and 2 samples are the noisy ones (bootstrap 0.33-0.91), and with 6 they become one point. On v3's data
this changes nothing for pe and se (se's two deepest bins, of 2 and 4 samples, share a window of 6), and merges
the deepest two points of pb and ont.

### Trained again on v3's tables

`retrain.sh` trained the four models on v3's training and test tables with the new leaves and knob minimum,
as the build runs the trainer (64 trees, seed 1); 2:15, 2:00, 0:33 and 0:36 on 6 threads. More long-read
samples (suggestions 2 and 3) need a new run. F1, species held out and on the test set:

| read type | leaves | held out, v3 → new | log loss | FP per sample | test at 0.5 | test at the knob curve |
|---|---|---|---|---|---|---|
| pe | 256 → 512 | 0.9578 → 0.9592 | 0.0573 → 0.0538 | 3.24 → 2.94 | 0.9533 → 0.9542 | 0.9621 → 0.9618 |
| se | 256 → 512 | 0.9527 → 0.9540 | 0.0708 → 0.0679 | 3.42 → 3.20 | 0.9486 → 0.9471 | 0.9574 → 0.9582 |
| pb | 256 → 128 | 0.9724 → 0.9717 | 0.0986 → 0.0872 | 1.06 → 1.00 | 0.9753 → 0.9736 | 0.9744 → 0.9736 |
| ont | 256 → 128 | 0.9731 → 0.9730 | 0.0875 → 0.0855 | 1.20 → 1.17 | 0.9696 → 0.9701 | 0.9717 → 0.9689 |

- **Short reads:** fewer false positives (pe 642 → 582 held out, test 3.91 → 3.60 per sample) and better
  calibration for 12-22 more misses; F1 +0.001 held out, the same on the test set within its noise (se at 0.5
  −0.0015, at the curve +0.0008). A knob minimum of 10 gave pe 0.9624 on the test set.
- **Long reads:** better calibrated (log loss −12% pb, −2% ont) at the same F1. ONT at the knob curve lost
  0.003 (34 + 32 → 27 + 45 errors); all of it at the two 3 Gb test samples, which the merged deep point (6
  samples, knob 0.93) now calls with a stricter knob than v3's (0.72, between its 4- and 2-sample points).
  The errors of those two samples are flat in the knob: 21, 19, 20, 18, 18, 17 and 23 at 0.5, 0.6, 0.7, 0.8,
  0.85, 0.9 and 0.93 (pb: 15, 14, 15, 13, 15, 20, 22), so the test set cannot tell the curves apart; the
  next run's 8 test samples per point and 1.5 times the training samples will.

The changes are worth keeping: they lower the false positives and the log loss everywhere, cost nothing
measurable in F1, and they are what the next run trains with.

## The knobs: how well the training data pin them

`knob_noise.py` (output in `knob_noise_output.txt`): per half decade of sample depth, the threshold of the
highest F1 on the scores with species held out, the thresholds within 0.002 F1 of it, and the 10th-90th
percentile of the best threshold over 200 bootstrap resamples of the bin's samples.

- **Shallow samples (below ~300 fragments):** flat; any knob from 0.05 to ~0.3 is within 0.002 F1. The curve's
  low knobs there do no harm and little good.
- **Middle and deep short-read samples:** well pinned. pe 10^4-10^5 fragments: best 0.94, bootstrap 0.88-0.95
  (+0.019 to +0.031 F1 over 0.5); 10M pairs: 0.78, bootstrap 0.68-0.87, so the drop from 0.94 is real, not
  noise (+0.020 F1 over 0.5). se the same (0.91-0.93, then 0.67-0.85).
- **Deep long reads:** not pinned. pb 4.5-5.0 (4 samples) bootstrap 0.74-0.77 but 5.0-5.5 (2 samples)
  0.36-0.83; ont 0.65-0.91 and 0.33-0.90. These are the points the minimum of 6 merges.

## Long-read strains missed with many fragments

`strains_missed.py` (output in `strains_missed_output.txt`): the present taxa with more than 100 fragments
that the model misses with species held out (p below 0.5): pb 10 of 761, ont 8 of 742, all strains (another
genome than the representative), and mostly the same species in both (*Blattabacterium cuenoti_B* with 10,224
and 8,818 fragments, *Vulcanococcus* sp019264575, *Streptomyces* sp041435355, *Mycobacterium phocaicum*,
*Methylocystis rosea*, *Alloprevotella* sp004556745, UMGS1470 sp900552105). Their medians against found
strains (pb): identity 0.969 against 0.998, `excess_high_share` 0.69 against 0.009, `congener_fit_share`
0.15 against 0.0005, `low_mapq_share` 0.47 against 0.005: the simulated genome is ~3% from the representative,
at the species boundary, and its reads look like those of a congener the database lacks. Not a tuning
problem: no threshold separates these from the false positives they resemble, and they are 1-2% of the deep
present taxa.

## Parity

`parity.log` confirms the failure: 24 pe samples (`rl100_p1000`, `rl250_p500000`) re-profiled, 4,587 taxa,
largest relative difference `adjacent_support` 8.07e-16, below the 1e-12 tolerance `9fd81c8` added.

## Run time (v3, 64 threads)

1:54 to the parity check (v2 on 16 threads: 3:38). The training data's simulation is now the critical
path: it ran from 0:13 to 1:08 in the background, the collection waited 27 minutes for it after the
training database was built at 0:41; collection 34 min, test set 9 min, training 2:17 for all four. Peak
scratch 120 GB (estimate ~100 GB), 210 GB free on the node.

The simulation log shows where the 55 minutes go: the 21 paired-end design points ran together and took
15:33; only then did the 14 long-read points start (252 samples, 64 at a time), and they took 39:27, of which
the last 17 minutes are the two ONT 6 Gb samples (pbsim3, one process each) after everything else had
finished (the test set: its two ONT 3 Gb samples 8 of its 16 minutes). Starting the long-read points with the
paired-end ones, deepest ONT first, or splitting a deep ONT sample into chunks simulated in parallel would cut
the simulation by 15-30 minutes and with it the collection's wait. With 36 long-read samples the shallow
long-read points grow, but the deep ones (4 and 2 samples) and so this path stay as they are.
