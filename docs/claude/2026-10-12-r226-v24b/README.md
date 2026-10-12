# The r226 v24b build (0.8.2, the larger download): lower F1 from a harder, more realistic sample mix

**Data.**
- **v24b:** SLURM 24171430, q512n11, 2026-10-11 23:05 to 2026-10-12 00:58, 84 threads. protal, the simulator and the
  scripts at `4a42f3b` (0.8.2; `build_metadata.tsv`), seed 1, no `--index-alleles`. The output directory is the
  earlier v24's (`/hpc-home/hildebra/DB/protal/protal0.8.1_r226_v24`), so the share archive has the same name; it
  replaced the earlier archive in `local/` and is unpacked to `local/v24b`. The earlier v24's logs stay in `local/v24`.
  The run uses the larger download (`refGenos`: 85,548 genomes of 45,482 species), Singletrack
  (`build_settings.tsv`; protal's runs report it), the weighted site shift in the training runs, and no indexed
  alleles (`Index alleles: none`).
- **Compared with** v24 ([2026-10-11-r226-v24](../2026-10-11-r226-v24/README.md); `7d5df4b`, 47,403 genomes, shift
  off, classic traceback, `--index-alleles 0.01`) and v23 ([2026-10-11-r226-v23](../2026-10-11-r226-v23/README.md);
  `8d8c89b`). The three builds share the seed, design and scenarios, but have other genomes and hold-outs.
- **Commands** (WSL, cores 0-3; the outputs are in this folder):
  - [builds_by_class.py](builds_by_class.py) (v24's script, any build directory name): `local <read type> v24 v24b`
    and `v23 v24b` → [builds_by_class.txt](builds_by_class.txt).
  - [builds_best_knob.py](builds_best_knob.py): finer classes (18: present taxa by genome and fragments, absent taxa by
    held-out congener distance), every model's rates at every build's mix at the best knob, and a bootstrap over
    samples → [builds_best_knob.txt](builds_best_knob.txt).
  - [refit.sh](refit.sh) with [thin_table.py](thin_table.py): v24b's pe model refitted on its own table (`work/`) with
    the build's trainer (`git archive 4a42f3b scripts`, species folds only, as in the v23 refits). The table is used
    whole, with the real strains thinned to v23's share, or with as many in-silico strains removed. Evaluated by
    [refit_classes.py](refit_classes.py) → [refits.txt](refits.txt).
  - [depth_by_class.py](depth_by_class.py) → [depth_by_class.txt](depth_by_class.txt);
    [strain_identity.py](strain_identity.py) → [strain_identity.txt](strain_identity.txt).
  - [shift_depth_pairs.py](shift_depth_pairs.py) on the small real world's builds with and without the shift
    (`~/rw_ws_{off,on}_s{1,2,3}`, [weighted_shift_runs.sh](../2026-10-11-r226-v23/weighted_shift_runs.sh)) →
    [shift_depth_pairs.txt](shift_depth_pairs.txt).
  - The aligning CPU per sample, from `logs/protal_cpu_{training,test}.tsv` of both builds (awk, section 5).

**Question (the user's, 2026-10-12).** The latest database build's results: how do they compare, and what changed?

## Summary

- **Every headline F1 is lower.** Species held out: pe 0.9608 → 0.9513, se 0.9661 → 0.9591, pb 0.9794 → 0.9714,
  ont 0.9742 → 0.9649. Test sets: −0.0036 to −0.0059. The scenarios' complex communities lost the most (pe
  soil_shallow hold-out 0.945 → 0.934).
- **About 80% of the drop is the harder sample mix.** The larger download doubled the real strains, and
  `--insilico-strains auto` now strains every one-genome species. Representatives (the database's own genome, the
  easiest class) fell from 62% to 29% of the present taxa; real strains rose from 18% to 36%. At v24b's mix and each
  model's best knob, v24's error rates would give pe 0.9578, against v24b's 0.9562 (v24 at its own mix: 0.9650).
- **At one mix, v24b's models are within noise of v24's and v23's.** At v24b's mix and the best knob, the differences
  in favour of v24 are pe +0.0017 [−0.0042, +0.0071], se +0.0001, pb +0.0012, ont +0.0024 [−0.0009, +0.0056] (95%
  bootstrap intervals over samples). ONT is the only one near the edge (v24 better in 92% of resamples, v23 in 96%).
  ONT is also where Singletrack changes the most alignments (47-49% of records in ties).
- **Real strains in training make the model call strains more readily.** At knob 0.5, false positives beside a
  held-out congener rose by about half in every read type (pe 1.29% → 1.92%). A refit with the real strains thinned
  to v23's share brings them back to 1.43%; removing as many in-silico strains does not (1.80%). But at each model's
  best knob the thinned refit gains only +0.0004 F1 (0.9567 against 0.9563). The real strains shift the probability
  scale, and the knob the build chooses absorbs that. Keep the training mix.
- **The weighted site shift lowers short-read depth estimates.** pe species depth fell 1-2 points within every class
  (representatives −1.5% → −3.2%), while PacBio's did not move. In the small real world, where both arms profiled the
  same samples, the shift lowers the depth of 19-22% of the called species and raises 1% (mean −1.1 to −1.4%). Its
  F1 effect at r226 is not separable here (the pe rate difference is within noise).
- **The build is faster and lighter.** 1:53 against 2:06 with twice the genomes. The collector's peak memory fell
  from 118.4 GB to 47.9 GB (0.8.2 releases each sample's evidence). PacBio aligning CPU fell 31% per collection.
  `database.protal` is 26.1 GB against 26.7 GB (no indexed alleles).

## 1. The headline numbers

F1 from each build's `console.log`: species held out at knob 0.5, the test set at the knob protal calls at (pe 0.74
in v23/v24, 0.72 in v24b; se 0.7 in v24b; ont 0.53 in v24; the others 0.5).

| | pe held out | pe test | se held out | se test | pb held out | pb test | ont held out | ont test |
|---|---|---|---|---|---|---|---|---|
| v23 | 0.9595 | 0.9709 | 0.9650 | 0.9662 | 0.9791 | 0.9828 | 0.9741 | 0.9820 |
| v24 | 0.9608 | 0.9697 | 0.9661 | 0.9677 | 0.9794 | 0.9878 | 0.9742 | 0.9791 |
| v24b | 0.9513 | 0.9652 | 0.9591 | 0.9630 | 0.9714 | 0.9842 | 0.9649 | 0.9732 |

pe with species held out: FP 3,780 → 4,795, FN 2,281 → 2,747. Both errors rose in every read type: FP pb 832 →
1,206, ont 1,114 → 1,590; FN pb 1,255 → 1,708, ont 1,557 → 2,058.

Scenario hold-out F1 (pe): gut 0.979 → 0.977, moderate 0.976 → 0.970, soil 0.965 → 0.955, soil_shallow 0.945 →
0.934, host 1.000 in both.

What the simulations drew from:

| | v23 | v24 | v24b |
|---|---|---|---|
| genomes / species in the genome table | 47,403 / 34,956 | 47,403 / 34,956 | 85,548 / 45,482 |
| species held out of the training database | 12,142 | 12,142 | 16,120 |
| simulated species another genome than the representative | 53.3% | 29.0% | 56.3% |
| of them real strains / in-silico strains | 14.5% / 38.8% | 14.5% / 14.5% | 30.4% / 25.9% |
| present taxa, pe held out: representative / real strain / in-silico | 31% / 15% / 54% | 62% / 18% / 19% | 29% / 36% / 35% |

## 2. Mix against rates

Each model's error rates by class (species held out) are applied to each build's class mix. There are 18 classes:
present taxa by genome simulated × 1-10, 11-100, >100 fragments, and absent taxa by the distance of the held-out
congener beside them, or by their neighbourhood. Each model is taken at the knob that is best for the mix
([builds_best_knob.txt](builds_best_knob.txt)). The coarser 7-class version at knob 0.5 is in
[builds_by_class.txt](builds_by_class.txt).

| F1, best knob | pe | se | pb | ont |
|---|---|---|---|---|
| v24 at its own mix | 0.9650 | 0.9680 | 0.9795 | 0.9744 |
| v24's rates at v24b's mix | 0.9578 | 0.9613 | 0.9726 | 0.9675 |
| v23's rates at v24b's mix | 0.9582 | 0.9619 | 0.9730 | 0.9682 |
| v24b at its own mix | 0.9562 | 0.9610 | 0.9716 | 0.9649 |
| the mix (v24 → v24b) | −0.0072 | −0.0067 | −0.0069 | −0.0069 |
| the rates (v24 → v24b) | −0.0016 | −0.0003 | −0.0010 | −0.0026 |
| v24 − v24b, 300 bootstraps over samples | +0.0017 [−0.0042, +0.0071] | +0.0001 [−0.0062, +0.0076] | +0.0012 [−0.0013, +0.0037] | +0.0024 [−0.0009, +0.0056] |
| v23 − v24b | +0.0020 [−0.0038, +0.0078] | +0.0011 [−0.0053, +0.0069] | +0.0016 [−0.0013, +0.0048] | +0.0032 [−0.0003, +0.0070] |

The bootstrap resamples each build's training samples independently: the builds' samples are not paired. The
intervals are wide (pe ±0.006) because samples range from 1,000 to 30 M read pairs. The rate parts are therefore
not distinguishable from noise for pe, se and pb. ONT's is the closest to a real loss.

Error rates at knob 0.5 by the coarse classes, pe (the other read types move the same way,
[builds_by_class.txt](builds_by_class.txt)):

| class | v23 | v24 | v24b |
|---|---|---|---|
| FN real strain, alleles | 6.51% | 6.51% | 4.95% |
| FN real strain, no alleles | 4.88% | 6.01% | 5.46% |
| FN in-silico strain | 2.75% | 3.35% | 3.14% |
| FN representative | 1.84% | 1.90% | 2.22% |
| FP beside a held-out congener | 1.41% | 1.29% | 1.92% |
| FP beside present congeners only | 0.61% | 0.59% | 0.65% |
| FP no congener in the sample | 1.22% | 1.24% | 1.32% |

The congener false positives rose within every distance bin (pe, the trainer's report): 0.01-0.02 23.9% →
28.8%, 0.02-0.05 5.0% → 7.3%, 0.05-0.15 0.70% → 0.94%. The FP that v24b calls have lower read identity: median
0.978 against 0.982, and 54% below 0.98 against 45% ([strain_identity.txt](strain_identity.txt)). The model accepts
more divergent evidence. The real strains themselves are hardly more divergent than v24's: 12.1% of those with more
than 10 fragments have identity below 0.98, against 10.3%.

## 3. What the real strains in training do (refits, pe)

v24b's pe model was refitted on its own training table, once whole and twice with present-strain rows removed at
random ([refits.txt](refits.txt), [builds_best_knob.txt](builds_best_knob.txt)). Rates are on the rows each refit
kept; F1 is at v24b's full mix.

| refit | rows removed | FP beside held-out congener (0.5) | FN real strain, alleles / none (0.5) | F1 at 0.5 | F1 at the best knob |
|---|---|---|---|---|---|
| full (reproduces the build: 0.9514 vs 0.9513) | none | 1.93% | 4.91% / 5.47% | 0.9514 | 0.9563 (0.76) |
| real strains to v23's share (31.7% kept) | 18,780 | 1.43% | 6.71% / 5.79% | 0.9549 | 0.9567 (0.67) |
| as many in-silico strains (30.1% kept) | 18,770 | 1.80% | 4.90% / 5.76% | 0.9518 | 0.9556 (0.68) |

- The real strains, not strains in general, push the model toward calling. With them thinned, the congener false
  positives fall back to v23's level and the real strains are missed as often as in v24.
- The knob undoes nearly all of it. At its own best knob for v24b's mix, the thinned refit's errors resemble the
  full refit's: congener FP 0.94% against 0.96%, real strains with alleles 8.3% against 7.6%. F1 differs by
  +0.0004. Fewer positive training rows lower the probabilities, as a higher knob would.
- So the training mix does not need changing. The build's knob choice (one knob from bootstrap resamples, pe 0.72
  in v24b) is what makes the difference. It was 0.0030 above 0.5 here, as in v24 (0.0031).

## 4. Composition: the shift lowers short-read depth

The species' depth / true depth − 1 of the present, called species (design samples, median), by class
([depth_by_class.txt](depth_by_class.txt)):

| pe training samples | v24 | v24b | pb training samples | v24 | v24b |
|---|---|---|---|---|---|
| representative | −1.5% | −3.2% | representative | +22.9% | +23.1% |
| in-silico strain | −1.6% | −2.7% | in-silico strain | +23.6% | +23.7% |
| real strain | −6.0% | −6.4% | real strain | +18.8% | +19.2% |
| all | −2.3% | −4.0% | all | +22.3% | +22.2% |

- pe's depth fell 1-2 points within every class, and the larger share of real strains (−6%) adds to it. PacBio did
  not move, and the shift applies only to short reads (`RunProtal.h`: `SetWeightedShift` unless a long read type).
- The small real world isolates the shift, with the same samples in both arms
  ([shift_depth_pairs.txt](shift_depth_pairs.txt)). With the shift on, 18.8-22.4% of the called pe species have a
  lower depth and 1.0-1.2% a higher one. The mean change is −1.1% to −1.4%, and the median error moves from −1.3..−1.8%
  to −1.7..−2.3%. PacBio is identical. The shift moves or ambiguates some of a species' own reads, so they leave its
  depth.
- The console's summary agrees: pe species' depth −4.8% (v24 +1.2%), explained share −0.011 (+0.018). PacBio's
  explained share is +0.19 in both builds: the known PacBio depth bias, unchanged.

## 5. Speed and resources

| | v24 | v24b |
|---|---|---|
| whole build | 2:06:28 | 1:52:56 |
| training database | 0:17:33, 42.1 GB | 0:21:58, 39.5 GB |
| training data collected | 1:08:11, peak 118.4 GB | 1:00:54, peak 47.9 GB |
| test data collected | 1:10:18, peak 66.7 GB | 1:03:25, peak 45.5 GB |
| protal_db in the background | 1:28:17, 46.3 GB | 1:15:27, 43.6 GB |
| `database.protal` | 26.7 GB | 26.1 GB |
| genome store on scratch | 46.3 GB (estimate) | 83.0 GB |
| protal CPU, training / test runs | 102,932 / 41,274 s | 87,711 / 35,347 s |
| start-up of a cohort run | 62 s (training run1) | 22 s |

Aligning CPU per collection, v24b / v24 (the same samples by name, but other genomes in them):

| | pe | se | pb | ont |
|---|---|---|---|---|
| training, total | 1.03 | 0.91 | 0.69 | 0.95 |
| training, per-sample median | 1.00 | 0.89 | 0.76 | 0.97 |
| test, total | 1.02 | 1.03 | 0.69 | 0.95 |

PacBio gains the most, as the WFA2 speed report expected (Singletrack and the four WFA2 patches). The samples hold
more divergent strains now, which costs alignment work, so pe and se staying flat is not a measure of the patches. The
collector's peak memory fell most likely because 0.8.2 releases each sample's evidence once it is profiled
(`0a1beaa`), and the runs no longer load indexed alleles.

## 6. What to do

1. **v24b is usable as the database.** Its models are within noise of v24's at the same mix, and its training mix
   (71% of the present taxa strains, half of those real) is closer to real samples, where most taxa are not
   the GTDB representative. Its lower F1 is the more honest estimate.
2. **Compare later builds with v24b's mix,** not v24's: the same download and seed keep the mix fixed.
   [builds_best_knob.py](builds_best_knob.py) compares two builds at one mix and their best knobs.
3. **ONT and Singletrack:** ONT's rates are the only ones near a real loss (−0.0024 to −0.0032; 92-96% of
   resamples). A local test (small real world with ONT reads, classic against Singletrack, three seeds) would tell
   whether Singletrack costs ONT accuracy before a cluster build with `--wfa-traceback classic` is worth running.
4. **The weighted site shift's depth cost:** about −1% of short-read depth, with no measured F1 gain at r226. Before
   keeping it as the default, find which reads it takes from their species (lower MAPQ after the shift, or reads
   moved to a congener). Then either leave those reads in the depth and use the shift only for the presence
   features, or turn the shift off. **Done 2026-10-12 (the user's decision):** off by default in protal and in
   `build_gtdb_database.py`; `--weighted_site_shift` turns it on, and v24b should be profiled with it.
5. **Name the next output directory after the version.** This run wrote 0.8.2 into `protal0.8.1_r226_v24`, the
   earlier v24's directory, so that database was replaced on the cluster unless it was moved first.
