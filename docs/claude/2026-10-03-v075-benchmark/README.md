# protal 0.7.5 against 0.7.3 (and 0.6.0a to 0.7.2): detection, abundance, speed and strains

- **Date**: 2026-10-03.
- **Code**: 0.7.5 = `b892133` (Version 0.7.5, branch `audit-fixes`), from `git archive`, Release, WSL. 0.7.3 =
  `b146951` and 0.6.0a = `014f4a9`, the v0.7.3 benchmark's builds (`../2026-10-02-v073-benchmark`). 0.7.4 was set in
  `8b2e376` the same day and never had a database of its own; this report's 0.7.5 holds 0.7.4's changes too.
- **World and databases**: the v0.7.1 benchmark's world (`../2026-10-01-v071-benchmark`, WSL `~/bench071`): 900
  simulated species (135 in no database), markers in operon-like clusters, genes of different conservation. 0.7.5 with
  the databases and four models its own `build_gtdb_database.py` trained on the release (`V075`, 14.6 min, 4.3 GB):
  the design and seed of the earlier pipelines (`--samples 4 --read-pairs 1000,20000,200000 --long-read-samples 4`,
  the test points of the v0.7.1 pipelines) and otherwise 0.7.5's defaults: the depth and divergence features in the
  default set (no knob curve), suspect gene copies, species priors, the mixed training design (lognormal sigma 1.3
  and 2.0, congeners), HiFi PacBio reads, 0.7.5's own `simulate_metagenomes` (0.7.1's, which the earlier pipelines
  used so that their training communities were the same, fails `--test` on what 0.7.5's collector asks of it), and
  its holdout design: **349 species left out instead of 290** (2 phylum, 4 class, 6 order, 8 family and 12 genus
  clades of 171 species, and 178 species alone), so the training communities and the "missing" database are not
  those of 0.7.3's pipeline. 0.7.3 with the v0.7.3 benchmark's `V073`; 0.6.0a with its own databases.
- **What is not here**: the 0.7.0, 0.7.1 and 0.7.2 pipelines' databases are gone from `~/bench071` (as are all
  earlier run folders), so those versions did not run again; their rows below are the v0.7.3 benchmark's, which found
  them reproducible to the fourth digit across its reruns (0.6.0a and 0.7.3 reproduce here to the fourth digit as
  well: 0.6338 / 0.9126 against 0.634 / 0.913).
- **Samples**: the v0.7.1 benchmark's paired-end (26, 2 of them 5M pairs) and single-end (24) samples; the v0.7.3
  benchmark's long-read samples (`samples_lr073`: PacBio HiFi reads with qualities by length, Nanopore by pbsim3; 3 and
  90 Mb, 4 samples each); the v0.7.2 benchmark's strain samples (12 species, 8 strains, one per sample).
- **Runs**: `scripts/profile.sh` (0.6.0a, 0.7.3 and 0.7.5, one sample after the other with every version in turn,
  `-t 6`, no qcmsa, timed; 252 runs), `scripts/long_reads.sh` (32), `scripts/strain_runs.sh` (0.7.3 and 0.7.5, the
  long reads also `--no_phasing`; 10). The missing-database runs of 0.7.5 are scored against its own held-out list
  (`scripts/score.py`), the others against 0.7.3's. `scripts/run_all.sh` runs everything; `scripts/collect_results.sh`
  scores and copies the tables into `results/`.
- **Machine**: WSL2, Intel Core Ultra 7 258V, 6 threads, 23 GB, otherwise idle (load 3-6 during the runs from the
  runs themselves). Wall times are this laptop's; repeated runs vary by 10-20%.

## Summary

- **Paired-end, full database**: 0.7.5 has the best F1 of every version, 0.9145 against 0.9126 (0.7.3), with higher
  precision (0.988 against 0.982) and the fewest false positives (1.27 per sample); per sample +0.0019 (−0.0024,
  +0.0075), not significant. Abundance as 0.7.3. **Single-end, full**: equal (0.8930 against 0.8928).
- **Species missing**: 0.7.5 equals 0.7.3 on F1 (pe 0.9021 against 0.9017, se 0.8811 against 0.8858, neither
  significant) but makes more false positives (pe 1.77 against 1.46, se 1.71 against 1.21 per sample) and its
  Bray-Curtis is worse on paired-end (+0.0088, significant). Part of that is the design: 0.7.5's training database
  lacks 349 species against 0.7.3's 290, so its samples hold more species with no reference whose reads land on
  relatives; what part is the model cannot be told from these runs (section "Species missing").
- **Long reads**: PacBio HiFi with species missing is 0.7.5's clearest gain, F1 0.8114 against 0.8001 (+0.0113,
  (+0.0004, +0.0246), significant), false positives 0.12 against 0.88 per sample and closer abundance (Bray-Curtis
  −0.011, significant); with the full database equal F1 (0.826) at 0.12 against 0.62 false positives. Nanopore: equal
  F1 (0.873 / 0.860), fewer false positives (0.38 / 0.75 against 0.62 / 1.12).
- **Strains**: paired-end MSAs and SNPs identical to 0.7.3's (SNP F1 0.967, true tree for 11 of 12 species).
  PacBio and Nanopore within 0.002 of 0.7.3 (phasing splits a few single-strain samples as in 0.7.3; with
  `--no_phasing` 0.907 against 0.908); PacBio trees right for 8 of 12 species against 7.
- **Speed and memory**: short reads as 0.7.3 (5M pairs 29.7 against 30.5 s, 500k pairs 5.7 against 6.1 s; single-end
  500k 3.3 against 3.0 s, +0.3 s per run, significant); long reads 20% faster (90 Mb: Nanopore 7.2 against 9.1 s,
  PacBio 4.7 against 6.0 s, −0.6 to −1.1 s per run, significant). Peak memory 3.42-3.45 GB against 3.45-3.48 GB:
  the packed index saves 70 MB of this world's 213 MB of values; the fixed 3.2 GB key map is the rest (at GTDB r226
  the saving is 7.7 GB, `../2026-10-03-memory-audit`).

## Paired-end and single-end reads

Means over all samples of a read type (`results/short/overall.md`; 0.7.0-0.7.2 and 0.7.3 `--knob 0.5` from the
v0.7.3 benchmark's table; archaea over all present archaea, every depth):

| reads | database | version | F1 | precision | recall | FP per sample | Bray-Curtis | median abs log2 | archaea recall |
|---|---|---|---|---|---|---|---|---|---|
| pe | full | 0.6.0a | 0.634 | 0.993 | 0.556 | 1.00 | 0.345 | 0.758 | 0.242 |
| pe | full | 0.7.0 | 0.903 | 0.976 | 0.869 | 2.96 | 0.131 | 0.399 | 0.736 |
| pe | full | 0.7.1 | 0.908 | 0.980 | 0.873 | 2.12 | 0.087 | 0.261 | 0.747 |
| pe | full | 0.7.2 | 0.911 | 0.982 | 0.875 | 1.88 | 0.085 | 0.262 | 0.749 |
| pe | full | 0.7.3 | 0.913 | 0.982 | 0.876 | 1.35 | 0.084 | 0.263 | 0.742 |
| pe | full | **0.7.5** | **0.915** | **0.988** | 0.876 | **1.27** | 0.085 | 0.261 | 0.745 |
| pe | missing | 0.6.0a | 0.629 | 0.982 | 0.560 | 1.50 | 0.357 | 0.789 | 0.225 |
| pe | missing | 0.7.0 | 0.895 | 0.963 | 0.865 | 2.00 | 0.131 | 0.408 | 0.697 |
| pe | missing | 0.7.1 | 0.900 | 0.967 | 0.869 | 1.54 | 0.093 | 0.295 | 0.704 |
| pe | missing | 0.7.2 | 0.906 | 0.977 | 0.869 | 1.08 | 0.088 | 0.293 | 0.700 |
| pe | missing | 0.7.3 | 0.902 | 0.963 | 0.872 | 1.46 | 0.089 | 0.292 | 0.708 |
| pe | missing | 0.7.5 (349 species out) | 0.902 | 0.966 | 0.872 | 1.77 | 0.098 | 0.289 | 0.706 |
| se | full | 0.7.2 | 0.888 | 0.982 | 0.843 | 1.75 | 0.107 | 0.333 | 0.700 |
| se | full | 0.7.3 | 0.893 | 0.984 | 0.847 | 1.08 | 0.106 | 0.336 | 0.708 |
| se | full | **0.7.5** | **0.893** | **0.987** | 0.847 | 1.08 | 0.106 | 0.330 | 0.702 |
| se | missing | 0.7.2 | 0.884 | 0.977 | 0.839 | 1.33 | 0.109 | 0.367 | 0.648 |
| se | missing | 0.7.3 | 0.886 | 0.978 | 0.838 | 1.21 | 0.111 | 0.378 | 0.652 |
| se | missing | 0.7.5 (349 species out) | 0.881 | 0.962 | 0.846 | 1.71 | 0.113 | 0.356 | 0.647 |

Paired differences per sample, mean (95% bootstrap interval; `results/short/paired.md`):

| reads | database | comparison | F1 | FP per sample | Bray-Curtis | wall (s) |
|---|---|---|---|---|---|---|
| pe | full | 0.7.5 less 0.7.3 | +0.0019 (−0.0024, +0.0075) | −0.08 (−0.54, +0.38) | +0.0005 (−0.0001, +0.0012) | −0.0 (−0.4, +0.3) |
| se | full | 0.7.5 less 0.7.3 | +0.0002 (−0.0037, +0.0047) | +0.00 (−0.58, +0.67) | +0.0004 (−0.0015, +0.0024) | +0.3 (+0.1, +0.5) |
| pe | missing | 0.7.5 less 0.7.3 | +0.0004 (−0.0096, +0.0120) | +0.31 (−0.65, +1.50) | **+0.0088 (+0.0016, +0.0162)** | −0.1 (−0.4, +0.1) |
| se | missing | 0.7.5 less 0.7.3 | −0.0047 (−0.0149, +0.0051) | +0.50 (−0.21, +1.21) | +0.0023 (−0.0054, +0.0095) | +0.1 (−0.1, +0.3) |
| pe | full | 0.7.3 less 0.6.0a | +0.2789 (+0.1854, +0.3704) | +0.35 (−0.54, +1.12) | −0.2605 (−0.3371, −0.1823) | −6.4 (−9.5, −3.7) |

By depth and abundance: `results/short/summary.md`.

### Species missing

0.7.5's missing-database runs are not like-for-like with 0.7.3's: its pipeline holds out 349 species (its holdout
design: whole clades at every rank plus species alone) where 0.7.3's held out 290, and the scorer counts each
version's own held-out species as unfindable. A sample then has more species with no reference, whose reads land on
relatives and make false positives, and whose share is renormalised away from the true abundances, which is where
the extra 0.3-0.5 false positives per sample and the worse Bray-Curtis can come from. The v0.7.3 benchmark saw
0.7.3's knob curve itself add +0.38 false positives per sample with species missing, and 0.7.5 has no knob curve,
so the model's part may well be favourable; telling the two apart needs 0.7.5's pipeline run with 0.7.3's held-out
list (`build_gtdb_database.py --holdout-species V073/heldout_species.txt`) or 0.7.3's with 0.7.5's. Not done here.

## Long reads

The v0.7.3 benchmark's long-read samples (`samples_lr073`), both versions against both databases
(`results/long_reads/overall.md`, `paired.md`; 0.7.0-0.7.2 from the v0.7.3 benchmark):

| reads | database | 0.7.0 | 0.7.1 | 0.7.2 | 0.7.3 | 0.7.5 | 0.7.5 less 0.7.3 per sample |
|---|---|---|---|---|---|---|---|
| pb (HiFi) | full | 0.818 (3.00) | 0.824 (1.38) | 0.817 (2.75) | 0.826 (0.62) | 0.826 (**0.12**) | −0.0001 (−0.0010, +0.0009); FP −0.50; Bray-Curtis **−0.0007 (−0.0015, −0.0001)** |
| pb (HiFi) | missing | 0.792 (2.50) | 0.798 (1.75) | 0.793 (2.62) | 0.800 (0.88) | **0.811** (**0.12**) | **+0.0113 (+0.0004, +0.0246)**; FP **−0.75 (−1.38, −0.25)**; Bray-Curtis **−0.0112 (−0.0187, −0.0048)** |
| ont | full | 0.865 (1.25) | 0.870 (0.88) | 0.872 (0.38) | 0.872 (0.62) | 0.873 (0.38) | +0.0007 (−0.0060, +0.0093); FP −0.25 |
| ont | missing | 0.843 (2.38) | 0.854 (1.75) | 0.845 (2.50) | 0.858 (1.12) | 0.860 (0.75) | +0.0015 (−0.0103, +0.0147); FP −0.38 |

F1 (false positives per sample). 0.7.5's precision is 0.998-0.999 on PacBio and 0.980-0.997 on Nanopore; its recall
is 0.7.3's within 0.003. Here 0.7.5's missing database (349 species out) does not cost it.

## Strains

The v0.7.2 benchmark's strain samples, each version against its finished database (`results/strains/summary.md`):

| reads | version | MSA columns (raw / qcmsa) | rows after qcmsa | SNP precision | SNP recall | SNP F1 | miscalled SNPs per 100 kb | true topology |
|---|---|---|---|---|---|---|---|---|
| pe | 0.7.3, 0.7.5 | 114,695 / 114,695 | 8.0 | 0.998 | 0.939 | 0.967 | 1.9 | 11 of 12 |
| pb | 0.7.3 | 115,188 / 113,924 | 7.4 | 0.951 | 0.862 | 0.904 | 47.5 | 7 of 12 |
| pb | 0.7.5 | 115,171 / 114,064 | 7.6 | 0.948 | 0.861 | 0.903 | 50.9 | 8 of 12 |
| pb | 0.7.3 `--no_phasing` | 115,220 / 114,112 | 8.0 | 0.953 | 0.868 | 0.908 | 45.8 | 6 of 12 |
| pb | 0.7.5 `--no_phasing` | 115,184 / 114,077 | 8.0 | 0.951 | 0.867 | 0.907 | 47.9 | 7 of 12 |
| ont | 0.7.3 (with or without phasing) | 114,681 / 114,110 | 8.0 | 0.994 | 0.889 | 0.939 | 5.7 | 8 of 12 |
| ont | 0.7.5 (with or without phasing) | 114,679 / 114,108 | 8.0 | 0.994 | 0.889 | 0.938 | 5.9 | 8 of 12 |

Paired-end is byte-for-byte 0.7.3's MSA. The long-read MSAs differ by a few columns (the long-read alignment through
every chain link, `29369eb`, came after 0.7.3) with the same SNP F1 within 0.002. Every sample holds one strain, so
the phasing's split rows (PacBio, as in 0.7.3) are wrong and are left out by the scorer. Strain runs took 21 / 37 / 50
s (pe / pb / ont) against 19 / 47 / 65 s.

## Speed and memory

Mean wall / CPU seconds per run, both databases, `-t 6` (`results/speed.md`; 0.7.0-0.7.2 from the v0.7.3
benchmark, measured under load then):

| reads | depth | 0.6.0a | 0.7.0 | 0.7.1 | 0.7.2 | 0.7.3 | 0.7.5 |
|---|---|---|---|---|---|---|---|
| pe | 1,000 pairs | 3.1 / 3 | 1.1 / 4 | 1.1 / 4 | 1.1 / 4 | 0.9 / 3 | 1.0 / 3 |
| pe | 10,000 pairs | 3.9 / 4 | 1.1 / 4 | 1.1 / 4 | 1.0 / 4 | 1.3 / 4 | 1.5 / 4 |
| pe | 500,000 pairs | 13.2 / 36 | 7.6 / 32 | 6.5 / 25 | 5.1 / 25 | 6.1 / 28 | 5.7 / 25 |
| pe | 5M pairs | 59.8 / 276 | 76.0 / 339 | 81.0 / 311 | 58.4 / 307 | 30.5 / 167 | 29.7 / 158 |
| se | 500,000 reads | - | 4.0 / 16 | 3.7 / 14 | 2.9 / 12 | 3.0 / 13 | 3.3 / 13 |
| pb | 90 Mb | - | 10.1 / 54 | 11.1 / 60 | 10.6 / 59 | 6.0 / 32 | 4.7 / 25 |
| ont | 90 Mb | - | 10.4 / 54 | 10.4 / 54 | 10.0 / 55 | 9.1 / 49 | 7.2 / 39 |

The 5M-pair and long-read rows of 0.7.3 are much faster than in the v0.7.3 benchmark (30.5 against 50.4 s; 6.0
against 10.2 s): that benchmark's deep runs shared the machine with a build; the 0.7.3 and 0.7.5 columns here ran in
turn on the same otherwise idle machine and are the comparison that holds. Per sample (`paired.md`): paired-end
−0.0 s (−0.4, +0.3) and −0.1 s, single-end +0.3 s (+0.1, +0.5) with the full database (the new features' profiling
on the smallest runs), long reads −0.6 to −1.1 s on every set (significant). Peak RSS 3.42-3.45 GB (0.7.5) against
3.45-3.48 GB (0.7.3) and 3.21-3.32 GB (0.6.0a).

## 0.7.5's pipeline

`results/pipeline_0.7.5/`: 14:37 wall and 4.3 GB peak for the conversion, two databases, 96 training and 24 test
samples, four models and their parity checks (all passed: "the same probabilities and features"). Its own evaluation
at knob 0.5 (0.7.3's in brackets, on 0.7.3's design, so not the same samples): species held out pe 0.9819 (0.9827),
se 0.9755 (0.9814), pb 0.9852 (0.9903), ont 0.9790 (0.9866); independent test set pe 0.9669 (0.9726), se 0.9671
(0.9687), pb 0.9779 (0.9828), ont 0.9730 (0.9827); false positives per sample 0.11-0.92 (0.08-1.03). The mixed design's
sigma-2.0 communities are the harder half of these sets (the r226 reports found the same), which the benchmark
samples above do not show.

## Not covered

Real samples, a GTDB-scale database (the r226 runs of the memory audit and the v7/v8 evaluation cover 0.7.3-era code
there), 0.7.5's missing database with 0.7.3's held-out list, strain mixtures, long-read strains of HiFi reads.

## Reproduce

```
bash docs/claude/2026-10-03-v075-benchmark/scripts/run_all.sh          # build, pipeline, all runs, scores (~45 min here)
bash docs/claude/2026-10-03-v075-benchmark/scripts/collect_results.sh  # scores again and copies the tables into results/
```

`~/bench071` must hold the v0.7.1 benchmark's world and samples, the v0.7.3 benchmark's `V073`, `samples_lr073` and
binaries, and 0.6.0a's `db060_*`; `profile.sh` builds 0.7.5 and its pipeline (`V075`) if they are not there.
