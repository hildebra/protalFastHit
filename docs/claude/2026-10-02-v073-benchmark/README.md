# protal 0.7.3 against 0.6.0a, 0.7.0, 0.7.1 and 0.7.2: detection, abundance, speed and strains

- **Date**: 2026-10-02.
- **Code**: 0.7.3 = `b146951` (Version 0.7.3, branch `audit-fixes`), from `git archive`, Release, WSL; the two commits
  after it (`29369eb` long reads through their chain, `87b857e` prefetched k-mer lookups) are not in it. The other
  versions are the v0.7.2 benchmark's builds (`../2026-10-02-v072-benchmark`): 0.6.0a = `014f4a9`, 0.7.0 = `4b21427`,
  0.7.1 = `1c11a00`, 0.7.2 = `fd1e719`.
- **World and databases**: the v0.7.1 benchmark's (`../2026-10-01-v071-benchmark`, WSL `~/bench071`): 900 simulated
  species (135 in no database), markers in operon-like clusters, genes of different conservation. 0.7.3 with the
  databases and four models its own `build_gtdb_database.py` trained on the release (`V073`): the design and seed of
  the other versions' pipelines and 0.7.1's simulator (so the same training communities), `--long-read-samples 4` (the
  number of long-read samples per point the earlier pipelines made; 0.7.3's default is 24), and otherwise 0.7.3's
  defaults: knob curves over the sample's depth for all four read types, PacBio HiFi training reads (`hifi_reads.py`),
  256 leaves per tree. "missing" = its training database, 290 species left out.
- **Samples**: the v0.7.1 benchmark's paired-end (26, 2 of them 5M pairs) and single-end (24) samples; long reads made
  again by 0.7.3's collector (`samples_lr073`, below) and the v0.7.2 benchmark's `samples_lr072`; the v0.7.2
  benchmark's strain samples.
- **Runs**: `scripts/profile.sh`: all five versions again, one sample after the other with every version in turn, and
  0.7.3 also at `--knob 0.5` (its knob curve off), `-t 6`, no qcmsa, timed: 552 runs. `scripts/long_reads.sh` (256),
  `scripts/strain_runs.sh` (0.7.3, and 0.7.3 `--no_phasing` for the long reads).
- **Scores**: the v0.7.1 benchmark's `score.py` through `scripts/score.py` (sets `short`, `lr073`, `lr072`); strains by
  the v0.7.2 benchmark's `strain_score.py` through `scripts/strain_score.py`; `scripts/speed.py` (time by depth).
- **Machine**: WSL2, Intel Core Ultra 7 258V, 6 threads, 23 GB, idle at the start (load 0.05), then shared with this
  session's build and tests of the parity fix at nice 19 (15:06-15:55, during the 500,000-pair and 5M-pair samples
  and half of the long reads; load before the runs 3-7). Wall times carry that noise.
- **Reproducibility**: 0.6.0a, 0.7.0, 0.7.1 and 0.7.2 made the same calls on every sample as in the v0.7.2 benchmark
  (same F1, Bray-Curtis and false positives to the fourth digit).

## Summary

- **Paired-end**: 0.7.3 has the best F1 with the full database, 0.913 against 0.911 (0.7.2), with the fewest false
  positives (1.35 per sample against 1.88), but with species missing it is below 0.7.2 (0.902 against 0.906; neither
  difference significant per sample). Its model alone (at `--knob 0.5`) equals 0.7.2's (+0.0001, +0.0002): the
  difference is the knob curve, +0.002 with the full database and −0.004 with species missing (not significant).
- **Single-end**: 0.7.3 is the best version, F1 0.893 against 0.888 (full) and 0.886 against 0.884 (missing); here
  the knob curve helps (+0.0064 per sample over `--knob 0.5` with the full database, significant; 0 with species
  missing).
- **Long reads**: on PacBio HiFi samples (0.7.3's collector, reads with qualities by length) 0.7.3 is the best version,
  F1 0.826 / 0.800 (full / missing) against 0.817 / 0.793 for 0.7.2 (+0.009 per sample, significant, and +0.007), with
  2 fewer false positives per sample and closer abundance. Nanopore: as 0.7.2 (0.872 / 0.858). 0.7.3's knob curves
  cost long reads little (PacBio −0.0002 and −0.0019, Nanopore −0.0016 and −0.0024 against `--knob 0.5`, none
  significant) but add 0.1-0.75 false positives per sample; 0.7.2's depth knobs cost up to 0.015. On the old pbsim3
  PacBio samples (every base quality 0) 0.7.3's PacBio model is worse than 0.7.2's (−0.006 and −0.010): it was trained
  on reads with real qualities.
- **Strains**: paired-end and Nanopore MSAs of 0.7.3 equal 0.7.2's (SNP F1 0.967 and 0.939, true tree for 11 and 8 of
  12 species). PacBio: 0.7.3's phasing splits 8 of the 96 single-strain samples into two haplotype rows (none for
  Nanopore); without them SNP F1 is 0.903 against 0.908, and the trees are right for 7 of 12 species instead of 6,
  partly because they have fewer leaves. With `--no_phasing` 0.7.3 equals 0.7.2.
- **Speed**: 0.7.3 is as fast as 0.7.2 or faster at every depth: per paired-end sample −1.3 s (−3.9, +0.1) with the
  full database, single-end −0.1 s; long reads alike. Memory 3.3-3.5 GB for every 0.7 version.
- **0.7.3's pipeline stopped at its parity check** (`adjacent_support` 1.4e-15 apart): a sum whose order depended on
  the threads. Fixed in `9fd81c8` (below); this benchmark's copy of the scripts accepted a relative difference up to
  1e-12.

## Paired-end and single-end reads

Means over all samples of a read type, at each version's default (`results/short/overall.md`; archaea over all present
archaea, every depth):

| reads | database | version | F1 | precision | recall | FP per sample | Bray-Curtis | median abs log2 | archaea recall |
|---|---|---|---|---|---|---|---|---|---|
| pe | full | 0.6.0a | 0.634 | 0.993 | 0.556 | 1.00 | 0.345 | 0.758 | 0.242 |
| pe | full | 0.7.0 | 0.903 | 0.976 | 0.869 | 2.96 | 0.131 | 0.399 | 0.736 |
| pe | full | 0.7.1 | 0.908 | 0.980 | 0.873 | 2.12 | 0.087 | 0.261 | 0.747 |
| pe | full | 0.7.2 | 0.911 | 0.982 | 0.875 | 1.88 | 0.085 | 0.262 | 0.749 |
| pe | full | 0.7.3 | **0.913** | 0.982 | 0.876 | 1.35 | **0.084** | 0.263 | 0.742 |
| pe | full | 0.7.3, `--knob 0.5` | 0.911 | 0.983 | 0.875 | 1.85 | 0.085 | 0.261 | 0.749 |
| pe | missing | 0.6.0a | 0.629 | 0.982 | 0.560 | 1.50 | 0.357 | 0.789 | 0.225 |
| pe | missing | 0.7.0 | 0.895 | 0.963 | 0.865 | 2.00 | 0.131 | 0.408 | 0.697 |
| pe | missing | 0.7.1 | 0.900 | 0.967 | 0.869 | 1.54 | 0.093 | 0.295 | 0.704 |
| pe | missing | 0.7.2 | 0.906 | 0.977 | 0.869 | 1.08 | 0.088 | 0.293 | 0.700 |
| pe | missing | 0.7.3 | 0.902 | 0.963 | 0.872 | 1.46 | 0.089 | 0.292 | 0.708 |
| pe | missing | 0.7.3, `--knob 0.5` | **0.906** | 0.975 | 0.871 | 1.27 | 0.089 | 0.291 | 0.708 |
| se | full | 0.7.2 | 0.888 | 0.982 | 0.843 | 1.75 | 0.107 | 0.333 | 0.700 |
| se | full | 0.7.3 | **0.893** | 0.984 | 0.847 | 1.08 | 0.106 | 0.336 | 0.708 |
| se | full | 0.7.3, `--knob 0.5` | 0.886 | 0.976 | 0.845 | 2.25 | 0.107 | 0.334 | 0.708 |
| se | missing | 0.7.2 | 0.884 | 0.977 | 0.839 | 1.33 | 0.109 | 0.367 | 0.648 |
| se | missing | 0.7.3 | **0.886** | 0.978 | 0.838 | 1.21 | 0.111 | 0.378 | 0.652 |
| se | missing | 0.7.3, `--knob 0.5` | **0.886** | 0.976 | 0.842 | 1.46 | 0.109 | 0.370 | 0.652 |

(0.7.0 and 0.7.1 single-end as in the v0.7.2 benchmark.) Paired differences per sample, mean (95% interval;
`results/short/paired.md`):

| reads | database | comparison | F1 | FP per sample | Bray-Curtis | wall (s) |
|---|---|---|---|---|---|---|
| pe | full | 0.7.3 less 0.7.2 | +0.0020 (−0.0046, +0.0076) | −0.54 (−1.12, 0.00) | −0.0007 (−0.0020, +0.0002) | −1.3 (−3.9, +0.1) |
| pe | full | 0.7.3 `--knob 0.5` less 0.7.2 | +0.0001 (−0.0009, +0.0014) | −0.04 | −0.0002 | −1.0 |
| pe | missing | 0.7.3 less 0.7.2 | −0.0038 (−0.0152, +0.0059) | +0.38 (−0.12, +0.92) | +0.0007 (−0.0007, +0.0023) | −0.1 |
| pe | missing | 0.7.3 `--knob 0.5` less 0.7.2 | +0.0002 (−0.0019, +0.0024) | +0.19 | +0.0006 | +0.5 |
| se | full | 0.7.3 less 0.7.2 | +0.0047 (−0.0004, +0.0095) | −0.67 (−1.50, +0.04) | −0.0006 | −0.1 (−0.2, −0.0) |
| se | full | 0.7.3 less 0.7.3 `--knob 0.5` | **+0.0064 (+0.0032, +0.0100)** | −1.17 (−2.12, −0.38) | −0.0008 | +0.1 |
| se | missing | 0.7.3 less 0.7.2 | +0.0019 (−0.0043, +0.0085) | −0.12 | +0.0015 (−0.0005, +0.0037) | −0.2 |

By depth (`results/short/summary.md`): archaea at 10,000 pairs or more 91% (0.7.3) and 93% (0.7.2) with the full
database, 88% and 88% with species missing; species at 0.01-0.1% at 10,000 pairs found 81-85% by every 0.7 version;
0.7.3's abundance error equals 0.7.2's. 0.7.3's knob curves (`results/pipeline_0.7.3/depth_knob_curves.txt`, log10 of
the sample's fragments: knob): paired-end 2.45: 0.26, 2.52: 0.26, 3.81: 0.66, 4.81: 0.71; single-end 2.39: 0.40,
3.73: 0.52, 4.74: 0.82; PacBio 1.80: 0.18, 2.04: 0.18, 3.20: 0.27, 4.23: 0.43; Nanopore 1.86: 0.13, 3.19: 0.35,
4.20: 0.63 (each from the pipeline's 3 design points of 4 samples).

## Speed and memory

Mean wall / CPU seconds per run (load average before), both databases, `-t 6` (`results/speed.md`):

| reads | depth | 0.6.0a | 0.7.0 | 0.7.1 | 0.7.2 | 0.7.3 |
|---|---|---|---|---|---|---|
| pe | 1,000 pairs | 4.0 / 4 | 1.1 / 4 | 1.1 / 4 | 1.1 / 4 | 0.9 / 3 |
| pe | 10,000 pairs | 3.6 / 4 | 1.1 / 4 | 1.1 / 4 | 1.0 / 4 | 0.9 / 3 |
| pe | 500,000 pairs | 11.9 / 36 | 7.6 / 32 | 6.5 / 25 | 5.1 / 25 | 5.2 / 25 |
| pe | 5M pairs | 93.1 / 398 | 76.0 / 339 | 81.0 / 311 | 58.4 / 307 | 50.4 / 268 |
| se | 500,000 reads | - | 4.0 / 16 | 3.7 / 14 | 2.9 / 12 | 2.7 / 11 |
| pb | 90 Mb | - | 10.1 / 54 | 11.1 / 60 | 10.6 / 59 | 10.2 / 57 |
| ont | 90 Mb | - | 10.4 / 54 | 10.4 / 54 | 10.0 / 55 | 10.0 / 55 |

The 5M-pair samples ran during this session's test build (load 5-7): 0.6.0a to 0.7.2 took 3-62% longer on them than
in the v0.7.2 benchmark, at 34-82% more CPU time (not traced). The four deep runs took 35-99 s with 0.7.2 and 0.7.3
(each the faster on two), 45-159 s with the earlier versions, and one version up to 3 times as long on one run as on
another (0.7.1 53 and 159 s), so these means say little. Peak memory 3.3-3.5 GB for every 0.7 version. 0.7.3's pipeline took 10 min
up to the parity check on the idle machine, and 1:42 min resumed (0.7.2's: 21.5 min under load 9).

## Strains

The v0.7.2 benchmark's strain samples (12 species of 8 strains, one strain per sample, 2-40x;
`../2026-10-02-v072-benchmark/scripts/make_strains.py`), 0.7.3 against its finished database
(`results/strains/summary.md`):

| reads | version | MSA columns (raw / qcmsa) | rows after qcmsa | SNP precision | SNP recall | SNP F1 | miscalled SNPs per 100 kb | true topology |
|---|---|---|---|---|---|---|---|---|
| pe | 0.7.1, 0.7.2, 0.7.3 | 114,695 / 114,695 | 8.0 | 0.998 | 0.939 | 0.967 | 1.9 | 11 of 12 |
| pb | 0.7.0-0.7.2, 0.7.3 `--no_phasing` | 115,220 / 114,112 | 8.0 | 0.953 | 0.868 | 0.908 | 45.8 | 6 of 12 |
| pb | 0.7.3 | 115,185 / 113,922 | 7.3 | 0.951 | 0.861 | 0.903 | 48.1 | 7 of 12 |
| ont | 0.7.0-0.7.3, with or without phasing | 114,681 / 114,110 | 8.0 | 0.994 | 0.889 | 0.939 | 5.7 | 8 of 12 |

0.7.3 phases long reads into strain rows (`Haplotypes.h`, `c5475ff`). Every sample here holds one strain, so a split
is wrong; with PacBio 0.7.3 split 8 of the 96 sample rows in two (`str_sN_hap1`, `_hap2`) in 5 species, 4 of them
with a log-odds of 0.6 or less (`results/strains/pb_phased_samples.txt`). The scorer finds no sample of that name,
so those rows are left out: the 0.903 is over the other rows, and one species (Damviia corcomia) has a right tree
with its 6 remaining leaves where all 8 give a wrong one. Nanopore split none. The PacBio strain reads are pbsim3's errhmm reads
given Q30, not the HiFi reads of `hifi_reads.py`, which may matter for the phasing.

## Long reads

`scripts/long_reads.sh` makes the long-read samples again with 0.7.3's collector (`samples_lr073`: the design and seed
of the v0.7.2 benchmark's `samples_lr072`). Its Nanopore reads are the same as `samples_lr072`'s (every version gives
the same results on both); its PacBio reads are HiFi reads of `hifi_reads.py` (qualities by read length) instead of
pbsim3's quality-0 reads. 0.7.0 to 0.7.2 were trained on quality-0 PacBio reads, 0.7.3 on HiFi reads; real HiFi reads
carry qualities, so `samples_lr073` is the one like real data. F1 (false positives per sample)
(`results/long_reads_0.7.3_samples/overall.md`):

| reads | database | 0.7.0 | 0.7.1 | 0.7.2 | 0.7.2, `--knob 0.5` | 0.7.3 | 0.7.3, `--knob 0.5` |
|---|---|---|---|---|---|---|---|
| pb (HiFi) | full | 0.818 (3.00) | 0.824 (1.38) | 0.817 (2.75) | 0.824 (1.75) | **0.826** (0.62) | **0.826** (0.00) |
| pb (HiFi) | missing | 0.792 (2.50) | 0.798 (1.75) | 0.793 (2.62) | 0.799 (1.62) | 0.800 (0.88) | **0.802** (0.12) |
| ont | full | 0.865 (1.25) | 0.870 (0.88) | 0.872 (0.38) | 0.873 (0.38) | 0.872 (0.62) | **0.874** (0.25) |
| ont | missing | 0.843 (2.38) | 0.854 (1.75) | 0.845 (2.50) | 0.860 (1.25) | 0.858 (1.12) | **0.861** (1.00) |

Per sample (`results/long_reads_0.7.3_samples/paired.md`): 0.7.3 less 0.7.2 PacBio **+0.0089 (+0.0020, +0.0212)**
and +0.0070 (−0.0003, +0.0170), false positives −2.1 and −1.75, Bray-Curtis −0.0024 and −0.0105 (all significant);
Nanopore +0.0000 and +0.0128 (−0.0008, +0.0284). Both at `--knob 0.5`: PacBio +0.0028 and +0.0032, Nanopore +0.0005
and +0.0005. 0.7.3's knob curves less `--knob 0.5`: PacBio −0.0002 (−0.0010, +0.0006) and −0.0019 (−0.0080, +0.0034),
Nanopore −0.0016 (−0.0094, +0.0042) and −0.0024 (−0.0065, +0.0013), with +0.1 to +0.75 false positives per sample
(PacBio significant); 0.7.2's depth knobs there: PacBio −0.0062 and −0.0057, Nanopore −0.0011 and −0.0148.

On the v0.7.2 benchmark's samples (`results/long_reads_0.7.2_samples/`), with quality-0 PacBio reads, 0.7.3's PacBio
model is worse than 0.7.2's: −0.0064 (−0.0107, −0.0025) and −0.0103 (−0.0192, −0.0024) per sample (F1 0.783 and
0.751 against 0.790 and 0.761). Its features (`excess_median`, which compares a read's differences with what its
qualities expect) meet reads unlike its training reads.

## 0.7.3's pipeline, and the parity check

The pipeline (`results/pipeline_0.7.3/`) first stopped at its step 7/8, the check that protal scores the models as the
trainer does and computes the features as during the collection:

    PROBLEM: protal computes features differently than when the training data was collected (largest relative
    difference): adjacent_support 1.4e-15

The same on every rerun. `adjacent_support` summed the smoothed shares of a taxon's gene links as doubles, read by
read on one thread and chunk by chunk on several; the collector profiles a design point's samples together (each on
part of the threads), the check its 8 samples together on all 6, so the sums were added in another order. Profiling
0.7.3's training samples together on 6 threads, alone on 1 and alone on 6, only `adjacent_support` differed (up to
1.2e-14, pe, pb and ont). `9fd81c8` sums the shares as integers (units of 2^-32): after it every numeric column of the
dumps is equal in all three ways, and `check_model_parity.py` notes a difference of 1e-12 or less as rounding instead
of failing. Here the pipeline was resumed with that tolerance in its copy of the scripts (`scripts/profile.sh`); with
it all four parity checks passed, with the probabilities equal. The models are 0.7.3's, unchanged.

Its own evaluation (`results/pipeline_0.7.3/summary.txt`, at knob 0.5; 0.7.2's in brackets): species held out pe
0.9827 (0.9835), se 0.9814 (0.9796), pb 0.9903 (0.9721), ont 0.9866 (0.9884); independent test set pe 0.9726
(0.9739), se 0.9687 (0.9690), pb 0.9828 (0.9809), ont 0.9827 (0.9827). The PacBio sets are HiFi reads now.

## Not covered

Real samples, a GTDB-scale database, strain mixtures, long-read strain samples of HiFi reads, the commits after
0.7.3 (`29369eb` changes long-read alignment). Deep-sample timings on a busy machine.

## Reproduce

```
bash docs/claude/2026-10-02-v073-benchmark/scripts/profile.sh     # 0.7.3 built, its pipeline, the pe and se runs
bash docs/claude/2026-10-02-v073-benchmark/scripts/long_reads.sh  # samples_lr073 by 0.7.3's collector, long-read runs
bash docs/claude/2026-10-02-v073-benchmark/scripts/strain_runs.sh
python3 docs/claude/2026-10-02-v073-benchmark/scripts/score.py ~/bench071   # short, lr073, lr072
python3 docs/claude/2026-10-02-v073-benchmark/scripts/strain_score.py ~/bench071
python3 docs/claude/2026-10-02-v073-benchmark/scripts/speed.py ~/bench071 results_v073 results_lr073
```
