# The r226 v12 build: what the four scenarios did to the models

**Data.** The r226 v12 build (SLURM 23974400, 2026-10-05, `3910d98` scripts with the collector's later fixes;
`local/protal_r226_v12_logs.tgz` and `_training.tgz`, extracted to `local/v12`): 24,959 species to simulate from
(6,000 with strains, 18,978 with their representative only, each given an in-silico strain), 9,524 held out of the
training database, the design of 0.7.6 and the four default scenarios (`gut`, `soil`, `soil_shallow`, `host`; 3 hold-in
and 2 hold-out samples each), `--features auto`. **Script.** [`scenario_ablation.py`](scenario_ablation.py), run
locally on the build's training and test tables:

    python3 docs/claude/2026-10-05-r226-v12-scenarios/scenario_ablation.py --training local/v12/training \
        --test local/v12/test --out docs/claude/2026-10-05-r226-v12-scenarios/results -t 6 --seeds 1,2,3,4,5
    # the weights: --out .../results_weights --variants "with scenarios,scenarios x0.5,scenarios x0.25,scenarios x0.1,design only"

Each read type's forest is fitted as the v12 trainer fitted it (the feature set all four trainers chose,
`normalized+adjacency+distance+depth+divergence+unfiltered`; 64 trees, 512 leaves for pe and se, 128 for pb and ont,
balanced classes, `max_features` sqrt) on variants of the training table, over 5 seeds, and scored as protal calls (knob
0.5; no knob curve, the sample's depth is a feature) on the build's test tables: the design's independent test set and
each scenario's hold-out samples. Outputs: `results/` and `results_weights/` (`runs.tsv` per seed, `summary.tsv`,
`by_depth.tsv`).

## The build

| read type | species held out F1 | test set F1 (v10) | gut | soil | soil_shallow | host |
|---|---|---|---|---|---|---|
| pe | 0.960 | 0.970 (0.963) | 0.995 | 0.959 | 0.937 | 1.000 |
| se (Ultima in the scenarios) | 0.971 | 0.972 (0.960) | | 0.967 | | 0.988 |
| pb | 0.966 | 0.981 (0.968) | 0.980 | 0.961 | 0.960 | 0.976 |
| ont | 0.960 | 0.974 (0.969) | 0.983 | 0.953 | 0.950 | 0.953 |

Scenario columns: hold-out samples (2 each), at knob 0.5. What the numbers rest on:
- **The scenarios' rows are most of the training data**: pe 146,852 of 218,383 rows (67%), se 52%, pb 80%, ont 80%,
  nearly all of them soil's. The "species held out" figures and FP per sample of the report (pe 6.84) are therefore
  mostly soil's (~200 false positives in a sample of ~10,000 species) and are not comparable with earlier builds; the
  design's independent test set is.
- **The hold-out samples score as the hold-in ones with species held out** (pe soil 0.959 vs 0.964, soil_shallow 0.937
  vs 0.941; pb and ont within 0.003): cross-validation by species predicts unseen samples of a scenario. The in-sample
  scores are 0.01-0.02 higher in soil.
- **The soil tables are crowded**: 12,880 species (7,728 held out) for samples of 9,000-11,000, so a sample takes ~78%
  of its table and the hold-out samples share most species with the hold-in ones (the held-out side is the same in
  both). That the hold-out matches the species-held-out estimate says it did not inflate the score much.
- **host** has 2 hold-out samples of ~21 present species each: one call is 0.01-0.02 of F1 (the ONT warning,
  0.953 vs 0.984, is that).
- **v10 to v12 on the test set** (+0.005 to +0.013) is not the scenarios' doing (below): v12 also has in-silico
  strains, the 2k/50k/200k depths and three times the species; and its test set is drawn from that larger pool.

## What the scenarios did

F1 change from removing the scenarios' rows (design only - with scenarios), mean of 5 seeds (SD 0.0002-0.002 on the
design test set and gut, up to 0.016 on soil and host):

| | pe | se | pb | ont |
|---|---|---|---|---|
| design test set | +0.0044 | 0.0000 | +0.0021 | +0.0025 |
| gut hold-out | +0.0010 | | +0.0079 | +0.0079 |
| host hold-out | 0.0000 | 0.0000 | +0.0138 | +0.0253 |
| soil hold-out | **-0.0097** | **-0.0095** | **-0.0268** | **-0.0607** |
| soil_shallow hold-out | **-0.0998** | | **-0.0681** | **-0.1453** |

- **For soil they are what makes the models work**: trained on the design alone the models miss 6-32% of a soil
  sample's species against 3-6% (pe shallow soil: 1,900 of 7,137 present missed against 363; ont 1,637 of 5,203
  against 288) and score 0.80-0.96 against 0.94-0.97. Ten thousand species, most with a few reads and 60% of them
  absent from the database, are not in the design (20-200 species per sample).
- **Everywhere else they cost a little**: 0.002-0.004 on the design's test set, up to 0.008 for long reads on gut.
- **The cost is a shift of the threshold, not lost ranking.** On the pe test set the scenarios' model calls more:
  false positives 153 against 85, false negatives 87 against 119, mostly at 200k-1M read pairs (FP 35 and 47 against
  13 and 15). The F1 at each model's best threshold is about the same (0.973 against 0.975); that threshold is 0.65
  with the scenarios and 0.43 without. A taxon of a few fragments in a deep sample is often present in soil and seldom
  in the design's samples, and no feature tells the two kinds of sample apart: `sample_log_taxa` (log10 of the taxa
  with reads in the sample, which protal could compute) as an extra feature recovered only 0.001 (pe) and nothing
  elsewhere.

## Down-weighting the scenarios' rows

F1 change against the v12 training (sample weight of the scenarios' rows, times the balanced class weights):

| | weight 0.5 | 0.25 | 0.1 | 0 (design only) |
|---|---|---|---|---|
| pe design test | +0.0017 | +0.0030 | +0.0040 | +0.0044 |
| pe soil / shallow | -0.0005 / +0.0003 | -0.0009 / -0.0003 | -0.0014 / -0.0015 | -0.0097 / -0.0998 |
| se design / soil | +0.0003 / -0.0009 | -0.0002 / -0.0019 | 0.0000 / -0.0026 | 0.0000 / -0.0095 |
| pb design test | +0.0007 | +0.0010 | +0.0004 | +0.0021 |
| pb gut / soil / shallow | -0.0005 / -0.0008 / -0.0005 | +0.0014 / -0.0006 / -0.0022 | +0.0011 / -0.0017 / -0.0052 | +0.0079 / -0.0268 / -0.0681 |
| ont design test | +0.0013 | +0.0016 | +0.0013 | +0.0025 |
| ont gut / soil / shallow | +0.0014 / +0.0004 / -0.0005 | +0.0017 / 0.0000 / -0.0016 | +0.0026 / -0.0005 / -0.0028 | +0.0079 / -0.0607 / -0.1453 |

A weight of 0.25 keeps the scenarios' gain (soil within 0.002) and recovers half to two thirds of their cost on the
design's test set (pe +0.003 of 0.0044, ont +0.0016 of 0.0025, pb +0.001 of 0.002) and part of gut's for long reads;
0.1 recovers nearly all of pe's at -0.0015 in soil.

## Other observations

- **The pe model's knob**: its best threshold is 0.698 with species held out (F1 +0.0046 over knob 0.5) and 0.699 on
  the test set (+0.0029), so a knob chosen on the training data alone would have gained on unseen samples. For se
  (0.635 vs 0.437), pb (0.496 vs 0.586) and ont (0.556 vs 0.471) the two disagree or sit near 0.5. With the sample's
  depth as a feature the trainer fits no knob at all; one global knob chosen with species held out, kept only when it
  gains 0.002 (as the depth knobs were), would have set pe's to ~0.70.
- **Strains are mostly in-silico now**: of the pe training rows' present strains 20,001 are in-silico and 6,784 real
  (the download's 6,000 species with strains against 18,978 one-genome species). The real ones are missed 2.4 times as
  often with species held out (6.1% against 2.5%; 1-10 fragments 13.0% against 6.2%), as in the test set (4.3% against
  2.0%). The in-silico strains (substitutions only, every gene by its factor) are easier than real ones, and they
  dominate what the models learn of strains. The build's warning "nearly all simulated species will be the database's
  own reference" (14.2%) is computed before the in-silico strains (52.2% after) and is misleading.
- **Long reads beat Illumina on shallow soil at the same bases** (1.5 Gb: pb 0.960, ont 0.950, pe 0.937), Ultima beat
  Illumina on soil (0.967 against 0.959 at 6 Gb).

## Recommendations

1. Keep the scenarios: they are the only training data like soil, and without them soil's F1 falls by 0.01-0.15.
2. Down-weight their rows to ~0.25 in the trainer (a `--scenario-weight` option; the collector's `meta_scenario` says
   which rows), for +0.001-0.003 on the design's test set at no cost in soil.
3. Let the trainer choose one global knob with species held out when it gains 0.002 or more (pe: ~0.70, +0.003 on the
   test set).
4. Watch the real strains: more species with real strains in the download (`--species`), or fewer in-silico ones, so
   that the models learn strains from real genomes too; fix the genome-table warning to count the in-silico strains.
5. More host samples (they are cheap: 2-50 species) to make the host scores more than a few calls.

## Follow-up: recommendations 2 and 3 implemented

`random_forest_cmdline.py --scenario-weight` (default 0.25; `build_gtdb_database.py --scenario-weight`) weighs the
scenarios' rows in every forest the trainer fits (the final one, the held-out estimates, the feature sets' choice, the
studies). With `--depth-knobs` and the sample's depth a feature the trainer now chooses one knob for every sample
(`choose_global_knob`: the threshold of the highest F1 with species held out, the rows weighted as in the forests, kept
if it gains 0.002 over 0.5), written as a knob curve of one point, which protal already reads as that knob at every
depth. The four v12 models retrained with the trainer itself on the v12 tables (`--features
normalized+adjacency+distance+depth+divergence+unfiltered --depth-knobs --evaluation basic`, weight 0.25 against 1;
reports and metrics in [`retrained/`](retrained)):

| | knob chosen (0.25 / 1) | design test set F1, v12 as built | weight 1 + knob | weight 0.25 + knob (the new default) |
|---|---|---|---|---|
| pe | 0.69 / 0.70 | 0.970 | 0.973 | 0.973 |
| se | - / - (best 0.60 and 0.64 gain 0.0004 and 0.0009) | 0.972 | 0.972 | 0.972 |
| pb | - / - | 0.981 | 0.981 | 0.982 |
| ont | - / - | 0.974 | 0.974 | 0.975 |

Scenario hold-outs, pe at its knob, weight 0.25 (weight 1): gut 0.991 (0.993), soil 0.962 (0.962), soil_shallow
0.938 (0.940), host 1.000 (1.000); v12 as built 0.995, 0.959, 0.937, 1.000. Long reads at weight 0.25: gut pb 0.984 (v12
0.980), ont 0.984 (0.983), soil and shallow soil within 0.002 of v12. For pe the weight and the knob do the same thing,
calling less freely, and do not add up: either gives the test set 0.973. With both the knob chosen with species held out
(0.69) is above the test set's best (0.58), as soil's rows, a third of the weight, prefer ~0.7, and pe's gut hold-out
loses 0.004 against v12 as built (9 of its 756 present species missed against 2; 2 samples).
