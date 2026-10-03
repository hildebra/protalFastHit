# The r226 v4 and v2 reruns: did suggestions 1-4 help?

2026-10-03. Data: two SLURM runs of `build_gtdb_database.py` on GTDB r226, uploaded as
`local/protal_r226_v4_training.tgz` and `local/protal_r226_v2_training.tgz` (git-ignored; extracted to
`local/r226_v4/` and `local/r226_v2_new/`): the four models' predictions, test predictions, thresholds,
`metrics.json` and `varimp.tsv`, the parity logs, the collection and simulation logs, and the training and test
tables. No `classifier_training*.log` reports; the numbers below come from the `metrics.json` files and from the
predictions. Compared with v3 (job 23887614, [2026-10-02-r226-v3-training](../2026-10-02-r226-v3-training/README.md),
`local/protal0.7.3_r226_v3/`).

- **v4** (OUTDIR `protal0.7.3_r226_v4`, started 18:27 on 2026-10-02): the build of `6d2e241` with its defaults
  (256 leaves for every read type, 24 long-read samples, 4 per long-read test point: 126 long-read training and 22
  test samples per type). Its trainer ran at 20:13, after the checkout had moved to `6ca86ee`, so its knob points
  already need 6 samples (the pb and ont curves end in one deep point); this is inferred from the curves.
- **v2-new** (OUTDIR `protal0.7.3_r226_v2` again, started 19:36): `6ca86ee`, with suggestions 1-4 of the v3
  report (512 leaves for pe and se, 128 for pb and ont; 36 long-read samples, i.e. 186 per type; 8 per long-read
  test point, i.e. 42; knob points of 6 samples), on 52 rather than 64 threads.

`run_all.sh` in this folder makes every output here from the extracted runs (`*_output.txt`): `compare_runs.py`
(the headline numbers of each `metrics.json`), `table_diff.py` (which columns of two runs' tables differ),
`common_rows.py` (two runs' species-held-out scores on the rows both have), long-read models trained locally
on v4's and v2-new's tables and all scored on v2-new's test set (`local/r226_cross/{A,B,C}`), and
`curve_shrink.py` (test-set F1 at knob 0.5, at the knob curve, and at the curve with weak points set back to
0.5).

## Verdict

Both runs are sound and finished their parity check: protal's probabilities and features equal the training
data's for all four models in each run (`parity*.log`). The local retraining of v2-new's long-read models
reproduces the HPC run's test predictions byte for byte.

Suggestions 1-4 did what they were meant to do except raise F1, which stayed within noise:

- **512 leaves for short reads:** fewer false positives and better calibration on the same tables. With species
  held out, pe false positives 659 → 575 (−13%), log loss 0.0574 → 0.0536, F1 0.9577 → 0.9595; se 671 → 628,
  0.0706 → 0.0676, 0.9529 → 0.9539. The gain is in deep samples (pe F1 at 500k pairs +0.004, 2M +0.004, 10M
  +0.007); the shallowest lose a little calibration (pe log loss at 1k pairs 0.262 → 0.287). On the test set F1
  stays put: pe 0.9528 → 0.9533 at knob 0.5 and 0.9611 → 0.9599 at the curve, se 0.9488 → 0.9485 and 0.9574 →
  0.9571; v3 and v4 (same settings, features differing as below) differ by as much. Test false positives per
  sample at 0.5: pe 4.01 → 3.65, se 4.12 → 3.95.
- **128 leaves and 1.5 times the long-read samples:** the same F1 on one test set. Scored on v2-new's 42 test
  samples (1,936 pb and 1,928 ont present taxa), at knob 0.5: pb 0.9744 (v4's tables, 256 leaves) → 0.9746 (128
  leaves) → 0.9749 (v2-new's tables); ont 0.9691 → 0.9681 → 0.9684. On the 126 samples both runs train on
  (species held out), pb improves (F1 0.9725 → 0.9736, log loss 0.0944 → 0.0878, false positives 130 → 112),
  ont does not (0.9726 → 0.9726, 0.0875 → 0.0865). The ONT learning curve, the reason for 36 samples, is now
  nearly flat (log loss 0.120, 0.098, 0.095 at a quarter, half and all of 186 samples).
- **The larger long-read test set** showed what 22 samples could not: the long-read knob curves do not carry
  over to new samples. On v2-new's test set its curves cost pb 0.0014 and ont 0.0037 against knob 0.5 (ont 61 +
  61 errors at 0.5, 64 + 72 at the curve; most at 250 Mb, 54 → 65). Over the five long-read model sets of this
  report the curve's mean effect on the test set is 0.000 (pb) and −0.001 (ont), at most ±0.004, where it adds
  +0.007 to +0.009 for short reads in every run.
- **Knob points of 6 samples** merged the long-read 1.5 and 6 Gb points into one, as intended; on the 3 Gb test
  samples, which that point calls, the curve makes as many errors as knob 0.5 (ont 20 and 20, pb 17 and 15).

## The models against v3

F1 with species held out / on the test set at knob 0.5 / on the test set at the knob curve (as protal calls by
default). The long-read test sets differ: v3 and v4 22 samples, v2-new 42 (v4's 22 among them).

| read type | v3 | v4 | v2-new |
|---|---|---|---|
| pe | 0.9578 / 0.9533 / 0.9621 | 0.9577 / 0.9528 / 0.9611 | 0.9595 / 0.9533 / 0.9599 |
| se | 0.9527 / 0.9486 / 0.9574 | 0.9529 / 0.9488 / 0.9574 | 0.9539 / 0.9485 / 0.9571 |
| pb | 0.9724 / 0.9753 / 0.9744 | 0.9725 / 0.9737 / 0.9732 | 0.9725 / 0.9749 / 0.9735 |
| ont | 0.9731 / 0.9696 / 0.9717 | 0.9726 / 0.9714 / 0.9698 | 0.9702 / 0.9684 / 0.9647 |

v2-new's long-read held-out F1 is over 186 samples, 60 of them new at the shallow points; on the 126 both runs
have it is pb 0.9736 and ont 0.9726 (above).

**v3 → v4 is the binary's doing, and neutral.** The runs simulated the same samples (55,694 pe rows, the same
keys), but v4's binary has `4a0a8e8` (a fragment's bases on a gene counted once): `conserved_fast_record_ratio`
and `conserved_fast_kept_ratio` differ in 52,928 of 55,694 pe rows (median 0.6-0.8% relative),
`conserved_fast_depth_ratio` in 4,979 and `conserved_hit_share` in 3,892; and long reads are aligned through
their chain (`29369eb`). The models are as good: pe 0.9578 → 0.9577 held out, se 0.9527 → 0.9529, pb 0.9724 →
0.9725, ont 0.9731 → 0.9726. v4's and v2-new's short-read tables are identical, and v4's long-read training
and test tables are rows of v2-new's, unchanged, so v4 → v2-new isolates the settings.

## The knob curve: points that gain little

Each knob point's F1 gain over 0.5 on its training window (species held out) is large at the shallow and deep
ends of the short-read curves (+0.009 to +0.032) and small in between and for long reads (ont +0.000 to
+0.010, 3.2-3.7 log10 fragments +0.000 to +0.003 for every read type). Setting the points that gain less than
0.002 back to 0.5 (`curve_shrink.py`), test-set F1 against the full curve:

| read type | runs | full curve | points gaining ≥ 0.002 | change |
|---|---|---|---|---|
| pe | v3, v4, v2-new | 0.9621, 0.9611, 0.9599 | 0.9617, 0.9606, 0.9610 | −0.0004, −0.0005, +0.0011 |
| se | v3, v4, v2-new | 0.9574, 0.9574, 0.9571 | 0.9588, 0.9587, 0.9588 | +0.0014, +0.0013, +0.0017 |
| pb | v3, v4, A, B, v2-new | 0.9744, 0.9732, 0.9751, 0.9756, 0.9735 | 0.9767, 0.9732, 0.9757, 0.9757, 0.9765 | +0.0023, 0, +0.0006, +0.0001, +0.0030 |
| ont | v3, v4, A, B, v2-new | 0.9717, 0.9698, 0.9673, 0.9684, 0.9647 | 0.9712, 0.9698, 0.9678, 0.9684, 0.9652 | −0.0005, 0, +0.0005, 0, +0.0005 |

(A, B: v4's long-read tables at 256 and 128 leaves, scored on v2-new's test set.) The long-read rows score
overlapping test samples (v4's are among v2-new's), so they show how the curve varies between fits, not between
samples. For ONT no threshold helps much: with 0.01 its curves are 0.5 nearly everywhere and score 0.9681-0.9714,
equal to or 0.001 above knob 0.5.

## Run time

v2-new (52 threads) against v4 (64): training-data simulation 16:35 + 44:28 against 15:49 + 39:25 (the 120 more
long-read samples, on fewer threads), collection 35:10 for 768 samples against 29:49 for 648, test set 13:46
against 8:49. The critical path is unchanged: the long-read points start after the paired-end ones and the two
ONT 6 Gb samples finish last, 18 minutes after the rest.

## Suggestions

1. **A knob point needs a gain:** keep 0.5 where the point's F1 at its knob beats 0.5 by less than 0.002 on its
   window. In the trainer, next to `DEPTH_KNOB_MIN_SAMPLES`. On these test sets: se +0.0015 in all three runs,
   pb +0.001 (up to +0.003), pe and ont ±0.001.
2. **No knob curve for ONT** (`--depth-knob-read-types pe,se,pb`), until a run shows it helps: on the 42-sample
   test set it costs 0.004, and its points gain at most 0.010 in training.
3. **Keep** 512/128 leaves (fewer false positives, better calibration, no F1 cost), 8 samples per long-read test
   point (it made these comparisons possible), and 36 long-read training samples (pb's calibration improved; it
   costs ~6 minutes of simulation and ~5 of collection).
4. **The simulation's critical path** (from the v3 report) stands: start the deepest ONT samples first, with the
   paired-end points, or split them into chunks simulated in parallel; 15-30 minutes.

## Files

None more are needed to judge the models. To check that the runs completed and how long they took, the SLURM
stdout and stderr, `build_metadata.tsv` and `model_logs/` of each OUTDIR would do.
