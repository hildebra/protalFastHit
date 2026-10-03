# The r226 v5 and v6 runs: do the distance features, the read EM and the calibrated calls help?

2026-10-03. Data: `local/protal_r226_v5_{logs,training}.tgz` and `local/protal_r226_v6_{logs,training}.tgz`
(git-ignored; extracted to `local/v5/`, `local/v6/`). Both are builds of `27423c6` on GTDB r226: 768 samples each
(198 pe, 198 se, 186 pb, 186 ont; 78 pe/se and 42 pb/ont independent test samples), congener groups in the samples
(`--congeners 0.25:2-5`, new since `be35d15`), 512 leaves for short reads and 128 for long reads, the same held-out species
and clades, the same simulated samples (same row keys).

- **v5**: `--features normalized+adjacency+distance --call-mode fdr` (the run suggested by the
  [denoising report](../2026-10-03-denoising-implementation/README.md)).
- **v6**: the defaults (`normalized+adjacency`, knob curve), plus `--previous-procedure`, which also trains the old
  procedure (grid search, top features, 512 trees) for comparison. Its `model_logs` lack `build_metadata.tsv` and
  `summary.txt`; its trainer commands are in `classifier_training*.log`.

`compare_models.py` and `run_all.sh` make `compare_models_output.txt`. `run_all.sh` also retrains the models on v5's
tables with `+relatives` and with the default set (`local/v5_relatives`, `local/v5_adjacency`, v5's trainer arguments),
which puts every feature set on one table and one test set; the HPC runs did not train the relatives set as a model, only
cross-validated it in the report's feature-set table.

## Verdict

**Training: yes.** The distance features lower the log loss by 16-17% for short reads and cut false positives; the full
relatives set does more (-27%). **Predicted performance of protal: a small gain for paired-end only, and only against
what protal calls with by default.** `+distance` adds 0.004 F1 for pe at the knob curve (single-end 0.0005, long reads
none). The calibrated calls (`--call-mode fdr`) are worse than the knob curve for all four models, and the relatives set
is not better at the curve either. The EM acts through those features; there is no other test of it here.

## Training (species held out, knob 0.5; each row on the same tables)

| read type | adjacency (v6) / distance (v5) / relatives | F1 | AP | log loss | FP/sample |
|---|---|---|---|---|---|
| pe | | 0.9563 / 0.9625 / 0.9689 | 0.9899 / 0.9925 / 0.9942 | 0.0569 / 0.0477 / 0.0415 | 3.07 / 2.76 / 2.19 |
| se | | 0.9508 / 0.9579 / 0.9650 | 0.9884 / 0.9915 / 0.9926 | 0.0705 / 0.0586 / 0.0515 | 3.35 / 3.01 / 2.40 |
| pb | | 0.9719 / 0.9733 / 0.9735 | 0.9954 / 0.9959 / 0.9960 | 0.1002 / 0.0963 / 0.0932 | 1.00 / 0.96 / 0.91 |
| ont | | 0.9690 / 0.9716 / 0.9716 | 0.9946 / 0.9954 / 0.9955 | 0.0971 / 0.0902 / 0.0884 | 1.33 / 1.21 / 1.17 |

- The gain is in species beside a congener. pe, v6 → v5, false negatives by the rank shared with the closest other species:
  family 4.0 → 2.8%, order 3.2 → 2.0%, class 3.0 → 1.9%, genus unchanged (4.6%); strains 5.3 → 4.5%; false positives of
  absent taxa that are not congeners of a missing species 217 → 163.
- With whole clades held out v5 pe stays at F1 0.9604-0.9625 from species to phylum; v6 falls to 0.9507-0.9527 at
  genus to phylum.
- More training samples would still help: log loss falls with 25%, 50%, 100% of the samples in every model (pe v5
  0.0564, 0.0534, 0.0477).
- **The old procedure** (v6, species held out): pe F1 0.931 against 0.956 (1,206 against 608 false positives); se 0.922 against
  0.951; pb 0.971 against 0.972; ont 0.969 against 0.969 (a tie). As in the first r226 run: the grid search kept 7 of 37
  features for pe and doubled the nodes. The new trainer is better for short reads and equal for long reads, at scale.

## Predicted performance (independent test set)

F1 (false positives + false negatives) for what protal would call with: knob 0.5, the model's knob curve (the default of
`--call-mode curve` models), the target share of false calls (the default of `--call-mode fdr` models).

| read type | model | knob 0.5 | knob curve | target share (fdr) |
|---|---|---|---|---|
| pe | v6 adjacency | 0.9542 (243+180) | 0.9605 (146+214) | - |
| | v5 distance | 0.9615 (218+139) | **0.9642** (160+168) | 0.9573 (206+187) |
| | relatives (local) | 0.9625 (195+151) | 0.9616 (151+200) | 0.9516 (153+285) |
| se | v6 adjacency | 0.9513 (252+188) | 0.9588 (123+242) | - |
| | v5 distance | 0.9570 (231+158) | 0.9593 (137+225) | 0.9556 (186+211) |
| | relatives (local) | 0.9608 (195+158) | 0.9601 (138+217) | 0.9569 (180+206) |
| pb | v6 adjacency | 0.9714 (39+59) | 0.9773 (32+46) | - |
| | v5 distance | 0.9714 (38+60) | 0.9761 (33+49) | 0.9711 (40+59) |
| | relatives (local) | 0.9711 (38+61) | 0.9760 (29+53) | 0.9729 (41+52) |
| ont | v6 adjacency | 0.9698 (53+53) | 0.9708 (39+63) | - |
| | v5 distance | 0.9704 (50+54) | 0.9696 (42+64) | 0.9685 (42+68) |
| | relatives (local) | 0.9686 (48+62) | 0.9673 (41+73) | 0.9688 (43+66) |

- **Noise.** Retraining v6's feature set on v5's tables (`v5-adj-local` in the output) moves the test F1 by 0.0002-0.0014 at
  0.5 and by 0.0003-0.0020 at the curve (pe 0.9605 → 0.9591, ont 0.9708 → 0.9722). Differences below about 0.002 are
  not findings; that rules out se and long reads at the curve.
- **Distance features:** pe +0.0037 at the curve and +0.0073 at 0.5 (fewer false positives and fewer false negatives at 0.5;
  test log loss 0.0560 → 0.0484); se +0.0005 at the curve, +0.0057 at 0.5; long reads within noise (pb and ont -0.001
  at the curve). At 0.5 they mostly supply the calibration that v6's knob curve already supplies.
- **Relatives:** best at knob 0.5 for short reads (pe 0.9625 with 195 false positives against 243; se 0.9608) and best AP and
  log loss, but the curve brings them nothing (pe -0.0009, se -0.0007) where it adds +0.003 to +0.008 to the other two sets.
  They do not beat `+distance` at the curve anywhere, nor at oracle thresholds (below): their gain in ranking with species
  held out does not carry over to the test samples (as on the benchmark world).
- **The calibrated calls are worse than the curve for every v5 model:** pe -0.0069, se -0.0037, pb -0.0050, ont -0.0011.
  By depth (pe test): at 500 pairs they miss 30 of 53 present species (curve 3), at 2,000 pairs 18 of 177 (curve 3); from
  50,000 pairs they make more false positives (76 against 50; 57 against 39 at 1M); only at 10,000-200,000 pairs do they
  miss fewer present species (41 against 59 at 50,000). One target (0.03 per sample) cannot suit samples with 50 and
  with 2,000 present taxa; the F1 over targets on species held out is flat (0.9582 at 0.02, 0.9587 at 0.03, 0.9540 at 0.05).
- **Against the earlier runs.** v6 against `v2-new` (`6ca86ee`, the same settings, but samples without congener groups: a
  different, easier world): held-out F1 -0.003 (pe) and -0.003 (se); test F1 at the curve +0.0006, +0.0017, +0.0038 (pb),
  +0.0061 (ont). The new binary did not harm the old feature set. v5 against v2-new: pe test 0.9533 → 0.9615 at 0.5 and
  0.9599 → 0.9642 at the curve, false positives per sample 3.65 → 2.79 (two changes at once, world and features).

## Where the remaining errors are

`error_budget.py` (output `error_budget_output.txt`) sorts each model's test-set errors at its knob curve and gives the
F1 if a kind of error were avoided; it also gives the F1 at oracle thresholds, the best one per test depth, chosen on the
test set itself (an upper bound for any threshold rule). For v5 (distance):

| read type | F1 at the curve | FN that are strains (another genome than the representative) | FP whose closest simulated species the database lacks | FP of 1 fragment | F1 at oracle thresholds by depth |
|---|---|---|---|---|---|
| pe | 0.9642 (160 FP, 168 FN) | 140 (83%), +0.016 if avoided | 123 (77%), +0.013 | 90 (56%) | 0.9670 (+0.0028) |
| se | 0.9593 (137, 225) | 181 (80%), +0.021 | 102 (74%), +0.011 | 82 (60%) | 0.9650 (+0.0057) |
| pb | 0.9761 (33, 49) | 47 (96%), +0.014 | 29 (88%), +0.008 | 18 (55%) | 0.9798 (+0.0037) |
| ont | 0.9696 (42, 64) | 58 (91%), +0.017 | 36 (86%), +0.010 | 25 (60%) | 0.9754 (+0.0058) |

- **Misses are strains the database lacks**: 80-96% of the false negatives are species simulated from another genome than
  the GTDB representative the database holds, most with 1-10 fragments (the missed strains' reads 0.94-0.97 identical).
  Of the representatives, 28-44 short-read and 2-7 long-read species are missed.
- **False positives are relatives of species the database lacks**: 74-88% of them are taxa whose closest simulated
  species is one held out of the database (55-67% a congener); over half have a single fragment.
- **Thresholds have little left**: oracle thresholds by depth would add 0.003-0.006, and they are chosen on the test set.
- **The relatives set has no hidden gain**: even at oracle thresholds it stays below the distance set (pe 0.9645 against
  0.9670, se 0.9637 against 0.9650, pb 0.9788 against 0.9798, ont 0.9737 against 0.9754). Its better log loss and AP with
  species held out do not carry over to the test samples' design.

## A defect found on the way: gene-conservation factors differ between builds

v5 and v6 simulated the same samples, yet four pe columns differ in 4,800-50,800 of 56,761 rows: `conserved_fast_depth_ratio`,
`conserved_hit_share`, `conserved_fast_record_ratio`, `conserved_fast_kept_ratio` (median relative difference 1-4%; checked
with the `table_diff.py` of the [v4/v2 report](../2026-10-03-r226-v4-v2-rerun/README.md)). The two builds' conservation
estimates compared 31,615,841 and 31,615,837 copies (training database; 64 and 52 threads), and the final databases 29,509 and
29,510 species. `GeneConservation::Estimator::Take` (`src/SequenceUtils/GeneConservation.h`) gives each gene's first
`kMaxCopies` (16) copies to whichever thread asks first, so for genes with more copies the chosen copies, and the factors,
probably depend on thread timing (the count differed with the thread number here; I did not run two builds with the same
thread number). The effect is small (the noise above includes it) and the parity check cannot see it (it re-profiles
with the build's own database), but a rebuild of the same inputs is not reproducible. A fix: take the 16 copies with the
lowest genome ids (or the lowest hash of the genome) per gene.

## Suggestions

1. **Next r226 build: `--features normalized+adjacency+distance` with `--call-mode curve` (the default).** pe gains 0.004 at the
   curve and 0.007 at 0.5; long reads are unchanged. Do not use `--call-mode fdr`; try it again only with targets by depth.
2. **Relatives:** keep opt-in. Their gain with species held out (log loss -27%, false positives -29% for pe) does not reach
   the test set, even at oracle thresholds (above).
3. **More training samples** (log loss still falling in all four models).
4. **A deterministic gene-conservation estimate** (above).
5. The earlier suggestions (a knob point needs a gain of 0.002; no ONT curve) are still not implemented; ONT's curve again
   costs 0.001-0.003 against knob 0.5 for the distance and relatives models.

Not covered: real samples, and the singleton rule and the EM alone (the runs trained no model on the EM features
without the rest).

The website is not affected (no protal option changes).
