# The sample's tail as features: thin species in deep, species-rich, uneven samples

Date: 2026-10-16/17. Branch `audit-fixes` on `da69b40`; the code is in the commit that adds this report.

**Why.** At r226 v31 a deep (1.38×), species-rich (10,521 species), uneven (lognormal sigma 2.5) test soil sample lost
85% of its species with 1-2 marker reads. Their own evidence (identity, allele identity, fragments) looked like that of
the thin species the training samples call, but their median score was 0.11 against 0.89. The models had learnt from
`sample_log_fragments` and `sample_log_taxa` that a thin row in a deep, rich sample is a congener's spill: every deep,
rich training soil had even abundances ([v31 report](../2026-10-16-r226-v31/README.md), sections 2 and 8). Dropping
`sample_log_taxa` recovered half of se's lost thin species but little of pe's. `234675b` adds deep, rich, uneven
scenario samples to the training design; this report adds a feature that describes the sample's tail directly, the v31
report's other option.

**Data.** r226 v31's share archive (`local/protal0.8.5_r226_v31_share.tar.gz`, its `work/` tables unpacked to
`local/v31/`, not in git), and r226 v27's tables (`local/v27/`, which hold the same training scenario simulations).
**Machine:** WSL, 4 cores (`taskset -c 0-3 nice -n 5`), Python 3.14.7 and scikit-learn 1.9.1 as on the cluster.

## Summary

- **Two features** of the whole sample, the same for each of its taxa (`context::SampleComplexityOf`, beside
  `sample_log_taxa`), in every dump and table since this commit:
  - `sample_thin_share`: the sample's taxa with 1-2 fragments over its taxa with fragments;
  - `sample_thin_close_share`: of those thin taxa, the share whose identity is within 0.01 of `sample_identity`, the
    sample's own species' identity. A rare species of the sample reads as its abundant species do; a congener's spill
    onto a species the sample lacks reads some percent lower.
- **They are a group of their own, `tail`, not in the default set.** `--features <default>+tail` (a named set, so
  `--features auto` can also consider it) trains with them.
- **On v31's tables they recover half of the lost thin species and cost nothing measurable elsewhere** (pe and se
  retrained as built and with `+tail`, the tail columns computed from the tables exactly as protal computes them;
  section 3):
  - pe's problem soil sample: thin species missed 675 → 342 (median score 0.15 → 0.86); shallow-soil sample 2:
    454 → 178. se's problem soil sample: 546 → 247.
  - The soil hold-out F1 rises by 0.019 (pe) and 0.018 (se), shallow soil by 0.016 (pe). The price is more false
    positives in those samples (pe +126 and +145), far fewer than the misses saved.
  - With species held out and on the design's test set, F1 is level or better (pe +0.0001..+0.0006, se +0.0020 and
    +0.0030 held out). At the same fixed knob, every set is within -0.0002 of the model without them or better,
    except pe's gut hold-out at knob 0.5 (-0.0021, 7 more false positives in 4 samples).
  - `sample_thin_close_share` is the 5th (pe) and 4th (se) most important of the 100 features. `sample_thin_share`
    ranks 40th-43rd: alone it tells little, spill being thin too.
- **The margin of 0.01 is better than 0.03** in every comparison above (pe gut -0.0006 against -0.0024; se's problem
  sample 247 against 387 thin misses).
- **Recommendation:** train the next r226 build with `--features <default>+tail` and compare it with v31 on the
  scenario samples and with the default set on the real mocks, then decide on the default set. The new deep, rich,
  uneven training samples (`234675b`) give the models what these features describe from another side, and that build
  has both. Real samples carry contaminant DNA that looks like a tail to these features and is present too (caveats):
  check `tail` on the real mocks with the contaminants counted as present.

## 1. Where the tail shows: v27's training scenarios

`tail_from_tables.py summary` on r226 v27's training tables (the same scenario simulations as v31's). The table's
`fragments` and `identity` are the values protal summarises (`Fragments()`, `BaseIdentity()`), and its rows are the
sample's taxa: the taxa with fragments equal `10^sample_log_taxa` in all 359 pe and 339 se samples. So both features
come out of a table exactly as protal would compute them.

| pe soil sample | sigma | sample_log_fragments | present taxa of 1-2 fragments | thin share | close share at 0.03 | at 0.01 | at 0.005 |
|---|---:|---:|---:|---:|---:|---:|---:|
| 9 (even, deep, rich) | 1.0 | 5.91 | 14 | 0.530 | 0.136 | 0.053 | 0.039 |
| 5 (even, deep, rich) | 1.0 | 5.70 | 24 | 0.551 | 0.147 | 0.061 | 0.048 |
| 4 | 2.5 | 5.89 | 193 | 0.580 | 0.145 | 0.062 | 0.044 |
| 8 | 2.5 | 5.91 | 387 | 0.575 | 0.224 | 0.116 | 0.101 |
| 10 | 1.5 | 5.25 | 747 | 0.585 | 0.230 | 0.121 | 0.092 |

- The thin share is 0.47-0.71 in every scenario sample: most thin rows are spill, whatever the community.
- The close share follows the thin present species. At 0.01 the deep, even soils have 0.05-0.06, those with hundreds
  of thin species 0.12 (2.2×); at 0.03 1.6×, where reads of a congener within 3% pass for the sample's own. A tighter
  margin separates more, but one or two reads of a noisier read type then fall outside it more often: 0.01 was chosen
  before the retraining, and the retraining agrees (section 3).
- v31's problem sample (soil 4 of the test set) has a close share of 0.118 at 0.01 (`tail_test_samples.txt`), as v27's
  training soils with hundreds of thin species, against 0.054 for the deep test soil 2 (sigma 1.5).

## 2. Implementation

- `src/Profiling/SampleContext.h`: `kThinFragments` (2), `kThinIdentityMargin` (0.01); `SampleComplexity::thin_share`
  and `thin_close_share`, computed by `SampleComplexityOf` from every taxon's fragments and identity after the sample's
  median identity (order-free).
- `src/Profiling/Profiler.h`: `SampleEvidence::sample_thin_share`, `sample_thin_close_share`, written after
  `sample_identity` (`TaxonFeatures`).
- `src/Profiling/FeatureVersions.h`: both at version 1, 2026-10-16.
- `scripts/model_features.py`: `TAIL_FEATURES`, group `tail` (not in `DEFAULT_FEATURE_SET`), the named set
  `<default>+tail` in `FEATURE_SETS`, and its description.
- Tests: `SampleContext.TheSamplesTailIsItsThinTaxaAndHowManyOfThemReadAsItsOwnSpecies` (thin and close counts, the
  margin, order), `test_FalsePositiveFeatures.cpp` (both features of a two-taxon sample), e2e (the same in every row of
  a sample, the thin share = the rows of 1-2 fragments), `test_model_pmml.py` (the group, not in the default set, the
  named set). WSL, 4 cores: ctest 635 passed, e2e 192, script suites 133 (`test_model_pmml`, `test_gtdb_build`,
  `test_collector`), GTDB pipeline 8/8.
- Docs: [features.md](../../features.md#the-samples-tail-tail-2026-10-16-opt-in) (group table, section, situations),
  `version_changes.md` (Since 0.8.5).

## 3. Retrained on v31's tables

**How.** [retrain_with_tail.sh](retrain_with_tail.sh) adds the two columns to v31's pe and se training and test tables
with [tail_from_tables.py](tail_from_tables.py) `add` (margin 0.01, protal's; a first run used 0.03), then trains a GBM
on them twice with the build's trainer options (`trained_model*.report.txt` line 3) but `--threads 4` and without the
C 1 / C 8 extra test tables: the default set (as built) and the default set + `tail`. The as-built models do not read
the new columns, so the second run reused them. [compare_tail.py](compare_tail.py) → [tail.txt](tail.txt) compares
the two on the same rows (each at its own knob, as protal calls; 95% interval over 1,000 resamples of the samples);
[fixed_knobs.py](fixed_knobs.py) → [fixed_knobs.txt](fixed_knobs.txt) scores both at the same knobs. The margin-0.03
run: `tail_margin003.txt`, `fixed_knobs_margin003.txt`, `tail_test_samples_margin003.txt`. Runtime: 6.5 min per pe
model, 3.8 min per se model.

As in the v31 report's section 8, the local as-built model is a refit of the build's (0.34% of its calls differ, a
thread count's rounding), so both sides here are local.

**pe** (knobs: 0.69 as built, 0.71 with `tail`):

| set | F1 as built | F1 with tail | difference | 95% interval |
|---|---:|---:|---:|---|
| species held out, design | 0.9706 | 0.9709 | +0.0003 | [-0.0005, +0.0011] |
| species held out, scenarios (hold-in) | 0.9506 | 0.9507 | +0.0001 | [-0.0005, +0.0008] |
| test, design | 0.9674 | 0.9681 | +0.0006 | [-0.0009, +0.0021] |
| test, scenarios (hold-out) | 0.9337 | 0.9459 | +0.0122 | [-0.0002, +0.0254] |
| test, gut hold-out | 0.9931 | 0.9925 | -0.0006 | [-0.0029, +0.0007] |
| test, moderate hold-out | 0.9683 | 0.9700 | +0.0017 | [-0.0019, +0.0085] |
| test, soil hold-out | 0.9161 | 0.9346 | +0.0185 | [-0.0033, +0.0478] |
| test, soil_shallow hold-out | 0.9157 | 0.9321 | +0.0164 | [+0.0007, +0.0322] |

**se** (knobs: 0.5 as built, 0.64 with `tail`):

| set | F1 as built | F1 with tail | difference | 95% interval |
|---|---:|---:|---:|---|
| species held out, design | 0.9649 | 0.9669 | +0.0020 | [+0.0010, +0.0030] |
| species held out, scenarios (hold-in) | 0.9519 | 0.9549 | +0.0030 | [+0.0018, +0.0043] |
| test, design | 0.9617 | 0.9608 | -0.0010 | [-0.0030, +0.0013] |
| test, scenarios (hold-out) | 0.9350 | 0.9470 | +0.0119 | [+0.0025, +0.0258] |
| test, moderate hold-out | 0.9660 | 0.9697 | +0.0037 | [+0.0013, +0.0064] |
| test, soil hold-out | 0.9148 | 0.9323 | +0.0175 | [+0.0018, +0.0359] |

At the same knob (`fixed_knobs.txt`) se's design test set is +0.0005 (0.5) and +0.0019 (0.69): its -0.0010 above is
the higher knob the trainer chose with `tail`. pe's gut hold-out is -0.0021 at 0.5 and 0.0000 at 0.69.

**The hold-out samples one by one** (thin: present species with 1-2 marker reads):

| sample | close share | thin | thin missed, as built → with tail | median score of thin | FP |
|---|---:|---:|---|---|---|
| pe soil 4 (sigma 2.5, deep, 10,521 species) | 0.118 | 823 | 675 → 342 | 0.15 → 0.86 | 47 → 173 |
| pe soil_shallow 2 | 0.110 | 598 | 454 → 178 | 0.18 → 0.92 | 48 → 193 |
| pe soil_shallow 4 | 0.156 | 297 | 90 → 68 | 0.87 → 0.97 | 37 → 63 |
| pe moderate 3 | 0.080 | 59 | 35 → 19 | 0.57 → 0.82 | 5 → 13 |
| pe soil 2 (sigma 1.5, deep, even-ish) | 0.054 | 29 | 11 → 20 | 0.85 → 0.28 | 48 → 25 |
| pe moderate 2 | 0.060 | 36 | 12 → 24 | 0.83 → 0.29 | 34 → 19 |
| se soil 4 | 0.182 | 691 | 546 → 247 | 0.10 → 0.86 | 81 → 214 |
| se moderate 4 | 0.536 | 709 | 122 → 78 | 0.87 → 0.97 | 41 → 69 |
| se soil 2 | 0.132 | 30 | 7 → 10 | 0.85 → 0.83 | 124 → 73 |

- The features act both ways: a sample with a tail of its own species (high close share) gets its thin species
  back, one whose thin rows are spill (low close share) loses a few and calls fewer false positives. The gain is far
  larger: pe's soil and shallow-soil hold-out samples miss 804 fewer species for 291 more false positives.
- The false positives that come with the recovered species in pe's soil 4 and shallow-soil 2 (+126, +145) are congener
  spill those samples' tails now let through. These two samples have no counterpart in v31's training data; the next
  build's corner samples do (`234675b`).

**Margin 0.03** (`tail_margin003.txt`): pe soil 4 675 → 350, shallow soil 2 454 → 180, se soil 4 546 → 387; soil hold-out
+0.0161 (pe) and +0.0145 (se); pe gut hold-out -0.0024 [-0.0048, +0.0000]. 0.01 is better on every one.

## Caveats

- Four hold-out samples per scenario: the intervals of the scenario sets are wide, and the gain rests on the two
  samples the features were made for. The design's test set and the held-out estimates over 386-436 samples say only
  that the features cost nothing measurable.
- The features were computed from the tables, not by a protal run. The computation is the same (the counts check
  against `sample_log_taxa` in every sample), but the tables round to 9 significant digits, so a taxon exactly at the
  margin could fall otherwise.
- The models were trained on v31's training samples, which lack deep, rich, uneven samples. The next build trains on
  such samples (`234675b`), which may teach the models the same from the other features; whether `tail` still adds
  something then is what that build should say.
- **Contaminant DNA reads like a tail, and is present.** Another session's report written meanwhile
  ([mock regression](../2026-10-16-mock-regression/README.md)) finds most of 0.8.5's extra calls on real mocks to be thin,
  high-identity rows of DNA the libraries carry from species the mock's list lacks; the simulations had none. Such a row
  is a close thin taxon here, as a rare member of the community is, so a model trained with `tail` may call more of
  them. That is not wrong: their DNA is in the sample, so they are present species, true positives (the user,
  2026-10-17). Only a mock's list counts them as false positives, so the mocks must be scored with them counted as
  present before `tail` is judged there. The simulations gain contaminants on 2026-10-17
  ([contaminants](../2026-10-17-contaminants/README.md)); for the models to learn that a thin contaminant is present,
  they belong in the truth (as written there they are left out of it).

## Next

- Build r226 with `--features normalized+...+read_length+tail` (the named set `<default>+tail`), or train both sets
  and compare them on the scenario hold-out samples and the deep, rich, uneven training samples (held out whole).
- Add `tail` to the default set if that build agrees (the user's decision); the website then needs the two features.
