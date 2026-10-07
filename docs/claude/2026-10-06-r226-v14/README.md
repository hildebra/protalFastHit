# The r226 v14 build: gradient boosting and scenario samples at varied depths

**Data.** The r226 v14 build (SLURM 24005369, q512n17, 2026-10-06 14:40-19:49, 64 threads; `local/protal_r226_v14_logs.tgz`
and `_training.tgz`, extracted to `local/v14`): the same 24,980 species, 9,030 in-silico strains and 10,600 held-out
species as v13 (same seed). It is the first build with gradient-boosted presence models (500 rounds at a learning rate
of 0.05, 63 leaves), with the `ref` features in the default set (`normalized+adjacency+distance+depth+divergence+unfiltered+ref`,
no `--features auto`) and `--evaluation basic`. Its scenario samples have depths of their own: 6 hold-in and 3 hold-out
samples per scenario (v13: 3 and 2), each at 0.5-2x its preset's depth (training 0.52-1.96x, test 0.51-1.91x).
Scenario rows are weighted 0.25, with one global knob per read type, as in v13. Compared with the v13 build
(`local/v13`, [its report](../2026-10-06-r226-v13-soil/README.md)).

**Which commit.** `build_metadata.tsv` says `46e89f6` for the scripts, protal and the simulator, but the run used older
scripts. Its trainer fitted 500 rounds at 0.05 and the metadata says "500 rounds"; those were the defaults before
`15006c2` (250 at 0.1, committed 18:30). The build started at 14:40, after `14632a8` (12:54, build defaults) and
`cdce3c0` (14:09), and at 19:28 the trainer still read the old file. `provenance()` in `build_gtdb_database.py` runs
`git rev-parse HEAD` and `protal --version` when it writes the metadata, at the end of the run. A checkout pulled (and
rebuilt) during a five-hour run is therefore what the metadata records, and `git status` is clean after a pull, so no
"(scripts changed since)" either. In effect v14 is `14632a8`/`cdce3c0`.

Scripts here, run in WSL (6 cores) with the scikit-learn venv of the v13 report and the trainer of the working tree
(only its table loader and feature groups, which differ from v13's commit by `ref` alone):

    python3 per_sample.py --builds v13=local/v13,v14=local/v14 --out . > per_sample.txt
    python3 design_by_depth.py --builds v13=local/v13,v14=local/v14 > design_by_depth.txt
    python3 samples_held_out.py --build local/v14 > samples_held_out_v14.txt      # and local/v13
    python3 knob_by_taxa.py --build local/v14 --read-types pe,se,pb,ont > knob_by_taxa_v14.txt
    python3 ../2026-10-06-r226-v13-soil/soil_errors.py --build local/v14 --save v14_rows.tsv.gz > soil_errors_v14.txt
    python3 ../2026-10-06-r226-v13-soil/soil_partition.py --rows v14_rows.tsv.gz > soil_partition_v14.txt
    python3 ../2026-10-06-r226-v13-soil/fp_consistency.py --build local/v14 > fp_consistency_v14.txt
    python3 cross_fit.py --builds v13=local/v13,v14=local/v14 --out cross --models forest,gbm+ref,gbm,gbm250 -t 6
    python3 summarize_cross.py cross/sets.tsv > cross/summary.txt
    python3 samples_cv.py --build local/v14 --other local/v13 --read-types pe --out samples_cv -t 3

The joined rows (`v14_rows.tsv.gz`, 80 MB) stay in WSL `~/soil13`, beside v13's; everything else is here.
`per_sample.py` reproduces the build's summary table, `cross_fit.py` reproduces v13's forest (design test 0.9652, soil
hold-out 0.9535, shallow 0.9253, knob 0.70) and v14's boosted model (design test 0.9649, knob 0.82), and
`samples_cv.py` reproduces v14's "samples held out" estimate (soil 0.9610, shallow 0.9076).

## Summary

- **The design's test set is the same in both builds**, row for row: 22,305 pe taxa with the same labels and the same
  feature values (likewise se, pb, ont). On it, any difference comes from the model alone. pe is unchanged (0.9652 →
  0.9649 at its knob), se +0.002, pb +0.002, ont +0.004.
- **The soil scenarios gained +0.003 to +0.016** (pe soil 0.954 → 0.965, shallow soil 0.925 → 0.941; se soil 0.951 →
  0.965; pb 0.947 → 0.958 and 0.947 → 0.955; ont 0.941 → 0.944 and 0.934 → 0.945). Gut rose too (pb 0.986 → 0.992,
  ont 0.979 → 0.988). The hold-out samples are not the same in the two builds: they are re-simulated at other depths
  and are partly other communities. Refits on both builds' tables, scored on both builds' samples, put the gain in
  boosting (+0.003 to +0.010), which works only with v14's varied depths (fitted on v13's samples it fails at other
  depths: ont soil 0.928, pb soil 0.942), plus +0.003 to +0.005 for easier samples (pe, se, pb soil); for a forest
  v14's training samples change little: [attribution](#where-the-gain-comes-from).
- **The false positives fell most.** Per pe deep-soil sample they halved (152 → 71), the congeners of a novel species
  that is near-identical on the markers among them (108 → 51). The misses of strains far from their reference barely
  moved (108 → 98 per sample); they are now the largest class in deep soil (ceiling +0.013).
- **The sample's depth still stands for its scenario.** With whole samples held out, pe's shallowest shallow-soil
  sample (4.86 log10 fragments, below every other soil sample's depth) scores 0.758, against 0.920 with species held
  out. A real soil sample shallower than the simulations would be scored like it. Adding the sample's complexity
  (taxa with reads, low-identity share, median identity) as features fixes that (0.908 → 0.935 over the six
  shallow samples held out) and gains +0.001 to +0.003 on the hold-out samples:
  [below](#the-samples-depth-still-stands-for-its-scenario).
- The knob: pe 0.82 and se 0.73 cost the design's test set 0.003 and 0.0024 against 0.5, and gain soil 0.004-0.013.
  A knob per bin of the sample's taxa count gains nothing: with species held out the design's training samples prefer
  the high knob too.

## The design's test set: the same rows, a different model

| read type | v13 (knob): F1, FP, FN | v14 (knob): F1, FP, FN | v14 at 0.5 | v14 best (at) | v13 best (at) |
|---|---|---|---|---|---|
| pe | 0.70: 0.9652, 93, 220 | 0.82: 0.9649, 63, 251 | 0.9678 | 0.9687 (0.58) | 0.9656 (0.68) |
| se | 0.50: 0.9592, 186, 179 | 0.73: 0.9612, 88, 253 | 0.9636 | 0.9636 (0.50) | 0.9607 (0.56) |
| pb | 0.50: 0.9767, 55, 43 | 0.50: 0.9791, 54, 34 | 0.9791 | 0.9814 (0.63) | 0.9769 (0.54) |
| ont | 0.50: 0.9720, 66, 55 | 0.50: 0.9759, 56, 48 | 0.9759 | 0.9762 (0.44) | 0.9728 (0.56) |

Boosting is better on the design's samples for every read type (best F1 +0.003 to +0.005). For pe and se the global
knob gives that back. It is chosen with species held out over the training rows, weighted as in the fit, where the soil
scenarios carry about half the weight (pe: 271,078 soil and shallow-soil rows at 0.25 against 76,870 design rows at 1).
The design's test set prefers 0.5-0.6.
By depth ([`design_by_depth.txt`](design_by_depth.txt)) the knob's extra misses are taxa at 50k-1M read pairs (pe 220
→ 251 misses at the knob, 128 → 133 at 0.5); at 5M pairs, where the design's samples overlap shallow soil's depths,
pe has 13 false positives at the knob in both builds.

## The scenarios, as the builds report them

The hold-out samples (scored by the final model; per sample, as v14 has 3 and v13 2):

| | v13 F1 | v14 F1 | FP per sample | FN per sample |
|---|---|---|---|---|
| pe soil | 0.9535 | 0.9646 | 152 → 71 | 216 → 199 |
| pe shallow soil | 0.9253 | 0.9408 | 166 → 172 | 364 → 262 |
| se soil | 0.9511 | 0.9645 | 217 → 90 | 176 → 182 |
| pb soil | 0.9469 | 0.9577 | 161 → 109 | 216 → 193 |
| pb shallow soil | 0.9474 | 0.9553 | 89 → 103 | 168 → 134 |
| ont soil | 0.9412 | 0.9439 | 216 → 150 | 213 → 254 |
| ont shallow soil | 0.9339 | 0.9453 | 142 → 141 | 202 → 162 |
| gut pe / pb / ont | 0.9950 / 0.9856 / 0.9788 | 0.9958 / 0.9924 / 0.9878 | | |

The samples differ ([`per_sample.txt`](per_sample.txt)). Hold-out sample 1 of each soil scenario is the same community
in both builds (3,652 of 3,654 present pe taxa shared) re-simulated at another depth; sample 2 shares 80% of its taxa;
v14's sample 3 is new. Between communities F1 varies more than between builds: pe soil sample 2 is the hardest in both
(v13 0.946, v14 0.952), sample 1 the easiest (0.962 at 5.78 log10 fragments → 0.975 at 5.90 in v14). ont's soil score
barely rose because its sample 2 has 410 misses (287 of them with 1-2 fragments) against 170-180 in the others.

## Where the gain comes from

[`cross_fit.py`](cross_fit.py) refits each model on each build's training table and scores it on the shared design
test set and on both builds' hold-out samples, the knob chosen as the trainer chooses it (species held out, rows
weighted as in the fit). F1 at that knob, on v13's / v14's hold-out samples ([`cross/summary.txt`](cross/summary.txt);
every set, F1 at 0.5 and the best threshold in [`cross/sets.tsv`](cross/sets.tsv), each hold-out sample in
[`cross/samples.tsv`](cross/samples.tsv)):

| | forest, v13's table (= v13) | forest, v14's table | boosting, v13's table | boosting, v14's table (= v14) | 250 rounds at 0.1, v14's table | boosting without ref, v14's table |
|---|---|---|---|---|---|---|
| pe soil | 0.9535 / 0.9567 | 0.9541 / 0.9567 | 0.9583 / 0.9638 | 0.9589 / 0.9646 | 0.9582 / 0.9650 | 0.9567 / 0.9619 |
| pe shallow soil | 0.9253 / 0.9303 | 0.9284 / 0.9321 | **0.8862** / 0.9190 | 0.9348 / 0.9408 | 0.9343 / 0.9407 | 0.9311 / 0.9367 |
| pe design test (knob) | 0.9652 (0.70) | 0.9656 (0.71) | 0.9684 (0.66) | 0.9649 (0.82) | 0.9655 (0.79) | 0.9661 (0.76) |
| se soil | 0.9511 / 0.9539 | 0.9526 / 0.9552 | 0.9581 / 0.9637 | 0.9607 / 0.9645 | 0.9565 / 0.9611 | 0.9539 / 0.9584 |
| se design test (knob) | 0.9592 (0.5) | 0.9594 (0.5) | 0.9630 (0.5) | 0.9612 (0.73) | 0.9635 (0.5) | 0.9608 (0.5) |
| pb soil | 0.9469 / 0.9515 | 0.9474 / 0.9508 | 0.9541 / **0.9416** | 0.9552 / 0.9577 | 0.9555 / 0.9572 | 0.9521 / 0.9548 |
| pb shallow soil | 0.9474 / 0.9470 | 0.9477 / 0.9493 | 0.9537 / **0.9432** | 0.9546 / 0.9553 | 0.9533 / 0.9555 | 0.9508 / 0.9533 |
| pb design test | 0.9767 | 0.9760 | 0.9800 | 0.9791 | 0.9789 | 0.9786 |
| ont soil | 0.9412 / 0.9411 | 0.9414 / 0.9407 | 0.9479 / **0.9280** | 0.9494 / 0.9439 | 0.9486 / 0.9435 | 0.9480 / 0.9425 |
| ont shallow soil | 0.9339 / 0.9351 | 0.9355 / 0.9354 | 0.9440 / **0.9230** | 0.9437 / 0.9453 | 0.9429 / 0.9444 | 0.9395 / 0.9423 |
| ont design test | 0.9720 | 0.9718 | 0.9759 | 0.9759 | 0.9743 | 0.9738 |

- **The gain is boosting.** On v14's table, boosting against the forest: soil +0.005 to +0.009 (ont on v14's samples
  +0.003), shallow soil +0.006 to +0.010, design test +0.003 to +0.008 at 0.5 (pe and se give part of it back at their
  knobs, above), gut +0.007 to +0.012 for long reads.
- **Boosting needs v14's samples.** Fitted on v13's samples, all at the preset's depth, it scores samples at that depth
  well (v13's: pb soil 0.954, ont soil 0.948) and fails at other depths (v14's: pb soil 0.942, ont soil 0.928, ont
  shallow 0.923, all below the forest); pe fails even on v13's own shallow samples (0.886). It learns the depth as the
  scenario's label. Fitted on v14's six samples at 0.5-2x it holds on both builds' samples. For a forest v14's table
  changes little (-0.001 to +0.003).
- **The samples.** For the same model v14's hold-out samples are easier than v13's by +0.003 to +0.005 (pe, se, pb
  soil; pe shallow), not at all for ont and pb shallow. Split up, pe soil's +0.011 from v13 to v14 is +0.003 samples,
  +0.000 training table (with the forest) and +0.008 model; pe shallow soil's +0.016 is +0.005, +0.002, +0.009; se
  soil's +0.013 is +0.003, +0.001, +0.009; pb soil's +0.011 is +0.005, -0.001, +0.007; ont shallow's +0.011 is +0.001,
  +0.000, +0.010.
- **The `ref` features** add +0.001 to +0.004 to boosting in soil and shallow soil for every read type (design test
  -0.001 to +0.002; se compared at 0.5, as its knobs differ: +0.002 to +0.003). They are computed against the training
  database, which lacks the held-out species, so a reference whose close relatives were held out looks more unique there
  than in the shipped database (see `model_features.py`): real samples have to confirm them.
- **250 rounds at 0.1** (the trainer's default since `15006c2`) score within 0.0013 of 500 at 0.05 everywhere but ont's
  design test (-0.0016). For se its knob stays at 0.5 (500 rounds: 0.73), which costs soil 0.004 and gains the design
  test 0.002: the knob is chosen by a 0.002 margin and can flip between near-identical models.

## What the remaining soil errors are

[`soil_partition_v14.txt`](soil_partition_v14.txt), the v13 report's classes (first match wins for the misses), per
hold-out sample (v13 → v14):

| per hold-out sample | pe soil | pe shallow | se soil | pb soil | pb shallow | ont soil | ont shallow |
|---|---|---|---|---|---|---|---|
| FN, 1-2 fragments | 59 → 41 | 224 → 132 | 44 → 35 | 89 → 99 | 111 → 78 | 103 → 146 | 139 → 107 |
| FN, divergent strain (excess_scaled_median > 0.015) | 108 → 98 | 98 → 91 | 93 → 91 | 116 → 78 | 55 → 51 | 93 → 83 | 55 → 47 |
| FN beside a called congener with >= 5x its fragments | 34 → 34 | 22 → 18 | 30 → 30 | 6 → 7 | 1 → 1 | 8 → 11 | 3 → 3 |
| FP, congener of a novel species near-identical on the markers | 108 → 51 | 127 → 140 | 127 → 64 | 114 → 68 | 62 → 58 | 165 → 107 | 100 → 99 |
| FP, congener of a more divergent novel species | 43 → 19 | 34 → 28 | 87 → 26 | 42 → 38 | 26 → 43 | 45 → 38 | 42 → 38 |

F1 of v14's hold-out samples without a class (ceilings): pe soil 0.9646 → 0.9777 without the divergent strains'
misses, 0.9700 without the 1-2-fragment misses, 0.9711 without the near-identical relatives' false positives; pe
shallow 0.9408 → 0.9538, 0.9595, 0.9591; pb soil 0.9577 → 0.9690, 0.9720, 0.9669.

- In deep soil the near-identical relatives' false positives halved (pe, se, pb, ont -35 to -53%), and the divergent
  strains are now the largest class (+0.011 to +0.013 if gone). These are real strains of wide GTDB clusters, 2-5%
  from the one reference the index holds; as the v13 report found, the present features do not separate them from a
  sister species.
- In shallow soil boosting recovered taxa with 1-2 fragments (pe -41%, pb -30%), but the near-identical relatives'
  false positives did not fall (pe 127 → 140 per sample). There the three classes are of about equal weight.
- The false positives still recur ([`fp_consistency_v14.txt`](fp_consistency_v14.txt)): 67% of pe's soil false positives
  are in taxa that are false positives in another soil sample (v13 60%). The bound of a per-taxon prior from the hold-in
  samples is smaller now: pe soil 0.9646 → 0.9680 (+0.003; v13 +0.008), se +0.004, pb +0.002, ont +0.003; shallow
  +0.003 to +0.005. A prior learned from the design's samples instead gives +0.001 (pe, se) and nothing for long reads.

## The sample's depth still stands for its scenario

The v13 report found that a scenario sample's fixed depth identified the sample, so boosting learned per-sample
offsets. v14 varies the depths, and the hold-out samples no longer collapse (pe shallow 0.886 when boosting is fitted
on v13's table, 0.941 now). But six samples at 0.5-2x one preset still cover a narrow band of depths, and soil-like
communities exist in the training data only there. The trainer's "samples held out" estimate (5 folds grouped by sample) shows it
([`samples_held_out_v14.txt`](samples_held_out_v14.txt)):

| pe shallow soil, hold-in sample | log10 fragments | F1, samples held out | F1, species held out | FN, samples / species held out |
|---|---|---|---|---|
| 6 | 4.86 | 0.758 | 0.920 | 1,194 / 308 |
| 1 | 5.02 | 0.904 | 0.941 | 681 / 295 |
| 5 | 5.08 | 0.946 | 0.943 | 275 / 279 |
| 4 | 5.27 | 0.901 | 0.928 | 549 / 306 |
| 2 | 5.38 | 0.925 | 0.937 | 443 / 287 |
| 3 | 5.45 | 0.958 | 0.962 | 144 / 201 |

Held out, sample 6 sits below every other soil sample's depth, and the model scores it like the design's samples of
that depth: its present taxa's median score falls from 0.998 to 0.959, and 1,194 of its 3,090 present taxa (39%) drop
below the knob. Over the six shallow samples the estimate is 0.908, against 0.940 with species held out and 0.941 on
the hold-out samples, whose depths (0.78-1.43x) lie inside the training band. Deep soil shows it mildly (0.961 against
0.962). se, pb and ont do not: pooled over a scenario's samples, samples held out are within 0.002 of species held out
(se soil 0.963 against 0.965, pb shallow 0.952 against 0.953, ont shallow 0.944 against 0.945); per sample within 0.013.
A real soil sample of 2M read pairs, or one with fewer fragments for other reasons, would be scored like sample 6.

[`samples_cv.py`](samples_cv.py) refits pe's boosted model on v14's table with the same sample folds, without the depth
feature and with the sample's complexity, three features protal could compute from the profile before calling
(`soil_experiments.add_context` of the v13 report: `sample_log_taxa`, log10 of the taxa with reads; `sample_low_identity`,
the sample's fragment share on low-identity bases; `sample_identity`, its taxa's fragment-weighted median identity).
All at pe's knob 0.82 ([`samples_cv/samples_cv.tsv`](samples_cv/samples_cv.tsv)):

| pe | shallow soil, samples held out | soil, samples held out | gut, samples held out | design, samples held out | soil / shallow hold-out (v14) | soil / shallow hold-out (v13's samples) | design test (best) |
|---|---|---|---|---|---|---|---|
| v14 as built | 0.9076 | 0.9610 | 0.9929 | 0.9699 | 0.9646 / 0.9408 | 0.9589 / 0.9348 | 0.9649 (0.9687) |
| no depth feature | 0.9368 | 0.9638 | 0.9909 | 0.9696 | 0.9652 / 0.9398 | 0.9613 / 0.9330 | 0.9655 (0.9671) |
| + sample complexity | 0.9348 | 0.9636 | 0.9940 | 0.9693 | 0.9659 / 0.9419 | 0.9618 / 0.9358 | 0.9647 (0.9696) |
| + sample complexity, no depth | 0.9365 | 0.9643 | 0.9940 | 0.9690 | 0.9673 / 0.9420 | 0.9620 / 0.9363 | 0.9622 (0.9673) |

- Without the depth feature the failure is gone (shallow 0.908 → 0.937) at no cost on the hold-out samples, but gut
  and host held out lose a little, and the design's test set needs a higher knob (0.9638 at 0.5, against 0.9678).
- With the sample's complexity added and the depth kept, the model can tell a soil-like sample by its complexity
  rather than by its depth: shallow held out 0.935, gut 0.994 (best of all four), hold-out samples +0.001 (v14) to
  +0.003 (v13's samples), design test as built (best F1 0.9696, the highest).
- Both without the depth: best in soil, but the design's test set loses 0.003 at the knob.

The v12 report tried `sample_log_taxa` in a forest and found nothing; it was scored on samples whose depth already
told the scenario, where complexity adds no information. Its value is where the depth does not: samples held out,
and real samples at depths the simulations did not draw. One seed, pe only; the knob was not re-chosen per variant.

## The knob

[`knob_by_taxa_v14.txt`](knob_by_taxa_v14.txt): per bin of log10(taxa with reads), the knob of the highest species-held-out
F1 on the training rows of that bin, kept if it gains 0.002 over the global one. For pe every bin keeps 0.82, the
design's (< 1,000 taxa) included: with species held out the design's training samples prefer the high knob too, and
only the independent test set (another design) would rather have 0.58. se moves one bin (3.0-3.5) to 0.52, ont one
to 0.90; on the test sets this changes nothing (se, pe) or loses (ont design 0.9759 → 0.9741). The cost on the design's
test set stays: pe 0.9678 → 0.9649, se 0.9636 → 0.9612. At 0.5 soil would lose more (pe 0.9646 → 0.9587, shallow
0.9408 → 0.9276, se 0.9645 → 0.9606).

## The build

5 h 09 min in all: genome table and conversion 0:11, in-silico strains 0:03, `training_db` 0:34 with `protal_db` 0:27
in the background, training data collected in 3:56 (the training simulations 4:23 and the test set's 2:19 in the
background; the test set's profiling, 1:54, takes turns with it), training 0:16 (pe 7:43, se 15:41, pb 11:58, ont
12:48 on 16 threads each; pe's 413 s are 337 s of held-out estimates and 9 s for the final fit), parity 0:27,
relatives 1:01, `--add_model` 2:49. Scratch peak 252.7 GB of 263.1 GB free: the doubled scenario samples leave 10 GB to spare on
this node. With `15006c2` (250 rounds at 0.1, folds side by side) training would take about half; `fd85db9` would make
`--add_model` seconds.

## What to do next

1. **The next build at the current commit**, so that its metadata is true and it gets 250 rounds at 0.1, the folds side
   by side, `--add_model` in place, and the reads behind every model error (`--error-reads`, `652ed53`). Record the
   commit and versions when the build starts (and warn if they change before it ends), not when it writes the
   metadata.
2. **The sample's complexity as features** (a "sample" group: taxa with reads, low-identity share, median identity),
   computed by protal from the profile before calling and written to the dump. It removes the depth's role as a
   scenario label (pe shallow soil, samples held out, +0.027) at no cost on the samples tested. To check on all read
   types with species and samples held out before making it a default. A cheaper partial fix: wider scenario depths
   (e.g. 0.25-4x), which widen the band but keep the depth as a label inside it.
3. **Deep soil's misses are now strains 2-5% from their reference** (+0.011 to +0.013 if gone): strain alleles in the
   index, or the cluster's other genomes' markers, are the lever, not the model.
4. **Shallow soil's false positives from near-identical relatives did not fall** (pe 140 per sample): the
   complex-community features of branch `fp-features` (untested at r226) are aimed at them.
5. The global knob costs pe and se 0.003 on the design's test set; no per-taxa-count knob recovers it.

## Follow-up: the sample's complexity for every read type, and in protal

The same refit with whole samples held out for se, PacBio and Nanopore (`samples_cv.py --read-types <type> --knob <the
build's> --variants v14,+sample`, [`run_samples_cv_other.sh`](run_samples_cv_other.sh); tables
`samples_cv/samples_cv_<type>.tsv`), v14 as built → with the sample's complexity, F1 at each build's knob:

| | design, samples held out | gut, samples held out | soil, samples held out | shallow soil, samples held out | design test (best) | soil hold-out, v14's / v13's samples | shallow hold-out, v14's / v13's | gut hold-out, v14's / v13's |
|---|---|---|---|---|---|---|---|---|
| pe (0.82) | 0.9699 → 0.9693 | 0.9929 → 0.9940 | 0.9610 → 0.9636 | 0.9076 → 0.9348 | 0.9649 → 0.9647 (0.9687 → 0.9696) | +0.0013 / +0.0029 | +0.0011 / +0.0010 | -0.0005 / -0.0025 |
| se (0.73) | 0.9674 → 0.9671 | | 0.9629 → 0.9648 | | 0.9612 → 0.9590 (0.9636 → 0.9631) | +0.0012 / +0.0021 | | |
| pb (0.5) | 0.9756 → 0.9762 | 0.9896 → 0.9949 | 0.9564 → 0.9583 | 0.9516 → 0.9540 | 0.9791 → 0.9802 (0.9814 → 0.9809) | +0.0014 / +0.0015 | +0.0011 / -0.0038 | +0.0030 / +0.0025 |
| ont (0.5) | 0.9728 → 0.9743 | 0.9874 → 0.9927 | 0.9500 → 0.9521 | 0.9435 → 0.9475 | 0.9759 → 0.9742 (0.9762 → 0.9753) | +0.0057 / +0.0015 | -0.0010 / -0.0018 | +0.0046 / +0.0050 |

With whole samples held out every read type gains in every scenario (pe shallow soil +0.027, the rest +0.002 to
+0.005, gut for long reads +0.005); on the hold-out samples mostly +0.001 to +0.006, with a few small losses; on the
design's test set -0.0017 (ont) to +0.001 (pb), se -0.002 at its knob of 0.73 (-0.0006 at 0.5; a model trained with
them chooses its own knob). One seed each.

Implemented in `096567f`: protal computes the three features for every taxon of a sample
(`context::SampleComplexityOf` in `src/Profiling/SampleContext.h`, called in `MicrobialProfile::ApplySampleContext`;
`sample_identity` is the fragment-weighted median identity of the taxa with 10 fragments or more, as
`soil_experiments.add_context` computed it, of every taxon with fragments when none has 10) and writes them to the
training dump after `sample_log_fragments`; `model_features.py` has them as the group `complexity`, in the default set.
Profiles of the mini database's e2e samples are identical to fd85db9's on 1 and 3 threads, the dumps identical but
for the three columns. A training table of v14 or older needs `--features normalized+adjacency+distance+depth+divergence+unfiltered+ref`.
Since `d3d269a` the build records the versions it starts with (`build_metadata.tsv`, with the end's where they
differ, and a warning on the console).
