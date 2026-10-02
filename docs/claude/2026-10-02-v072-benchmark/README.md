# protal 0.6.0a, 0.7.0, 0.7.1 and 0.7.2: detection, abundance, speed and strains

- **Date**: 2026-10-02.
- **Code**: 0.6.0a = tag `0.6.0a` (`014f4a9`), 0.7.0 = `4b21427`, 0.7.1 = `1c11a00`, 0.7.2 = `fd1e719` (Version 0.7.2,
  branch `audit-fixes`; the uncommitted work of other sessions is not in it). Each from `git archive`, Release, WSL;
  0.6.0a, 0.7.0 and 0.7.1 are the v0.7.1 benchmark's builds.
- **World and databases**: the v0.7.1 benchmark's (`../2026-10-01-v071-benchmark`, WSL `~/bench071`): 900 simulated
  species (135 in no database), markers in operon-like clusters, genes of different conservation. Each 0.7 version
  with the finished database and four models its own `build_gtdb_database.py` trained on the release (0.7.2's new:
  `V072`, the same design and seed as 0.7.0's and 0.7.1's; its long-read models with depth knobs, its default), and the
  training database without 290 species ("missing"); 0.6.0a with its own builds and shipped model.
- **Samples**: the v0.7.1 benchmark's (26 paired-end incl. 2 of 5M pairs, their first reads as 24 single-end, 8
  PacBio, 8 Nanopore). Long reads also made again by 0.7.2's collector (`scripts/long_reads.sh`, see below). For the
  strains, new samples (`scripts/make_strains.py`, below).
- **Runs**: `scripts/profile.sh`: all four versions run again, one sample after the other with every version in
  turn, `-t 6`, no qcmsa, timed (`/usr/bin/time -v`); paired-end for all four, single-end, PacBio and Nanopore for the
  0.7 versions; 480 runs. `scripts/strain_runs.sh`, `scripts/long_reads.sh`.
- **Scores**: the v0.7.1 benchmark's `score.py` (wrapped by `scripts/score.py`, `scripts/long_reads_score.py`): F1,
  precision, recall, Bray-Curtis, median |log2| of abundance, archaea, low-abundance taxa, false positives, time and
  memory; paired differences per sample with a 95% bootstrap interval. Strains: `scripts/strain_score.py`.
- **Machine**: WSL2, Intel Core Ultra 7 258V, 6 threads, 23 GB, shared with other sessions' builds and tests: the load
  average before the runs had a median of 9.0 and a maximum of 24.3, so wall times carry noise; versions ran each
  sample back to back, so the paired differences are fairer than the means.
- **Reproducibility**: the 316 runs of 0.6.0a, 0.7.0 and 0.7.1 made the same calls (TP, FP, FN) as in the v0.7.1
  benchmark a day earlier, run for run.

## Summary

- **Paired-end**: 0.7.2 is the best version, by a little over 0.7.1: F1 0.911 against 0.908 (+0.002 per sample,
  interval +0.0002 to +0.005) with every species in the database and 0.906 against 0.900 with species missing
  (+0.005, not significant), with fewer false positives (1.9 against 2.1 and 1.1 against 1.5 per sample) and a
  slightly closer abundance (Bray-Curtis −0.0015 and −0.0048 per sample, both significant). 0.6.0a stays far behind
  (F1 0.63, Bray-Curtis 0.34). Single-end: 0.7.2 as 0.7.1 (F1 0.888), Bray-Curtis −0.004 with species missing.
- **Speed**: 0.7.2 is as fast as 0.7.1 or a little faster (paired-end runs 6.7 s against 7.1 s with the full
  database, 5.9 s against 7.1 s with the training database; CPU 22.6 s against 23.5 s), twice as fast as 0.6.0a
  (14.7 s). Memory 3.4-3.5 GB for every 0.7 version. Under the machine's load these are a few seconds per run.
- **Strains** (paired-end, 12 species of 8 strains each): the 0.7 versions call SNPs with F1 0.967 (precision 0.998,
  recall 0.94), 0.6.0a with 0.869 (precision 1.000, recall 0.77, 0.35-0.45 at 2-3x coverage); 0.7.1 and 0.7.2
  recover the true strain tree for 11 of 12 species, 0.7.0 for 10, 0.6.0a for 4. 0.7.1 and 0.7.2 give the same MSAs.
- **Long reads**: the benchmark's long-read samples, made by 0.7.1's collector, are half reads under 1 kb (pbsim3's
  sampling per genome), as 0.7.0's and 0.7.1's training reads were; 0.7.2 trains on reads of the intended length. On
  samples made again by 0.7.2's collector, 0.7.2's models at `--knob 0.5` are the best of the 0.7 versions (PacBio
  +0.004 F1, Nanopore +0.003 to +0.006 over 0.7.1), but **its default depth knobs cost it** (PacBio −0.004 and
  −0.007, Nanopore −0.001 and −0.015 with species missing), as on the v0.7.1 benchmark
  (`../2026-10-01-features-depth-knobs`).

## Paired-end and single-end reads

Means over all samples of a read type, at the model's knob 0.5 (archaea recall over all present archaea, every
depth):

| reads | database | version | F1 | precision | recall | FP per sample | Bray-Curtis | median abs log2 | archaea recall |
|---|---|---|---|---|---|---|---|---|---|
| pe | full | 0.6.0a | 0.634 | 0.993 | 0.556 | 1.00 | 0.345 | 0.758 | 0.242 |
| pe | full | 0.7.0 | 0.903 | 0.976 | 0.869 | 2.96 | 0.131 | 0.399 | 0.736 |
| pe | full | 0.7.1 | 0.908 | 0.980 | 0.873 | 2.12 | 0.087 | 0.261 | 0.747 |
| pe | full | 0.7.2 | **0.911** | 0.982 | 0.875 | 1.88 | **0.085** | 0.262 | 0.749 |
| pe | missing | 0.6.0a | 0.629 | 0.982 | 0.560 | 1.50 | 0.357 | 0.789 | 0.225 |
| pe | missing | 0.7.0 | 0.895 | 0.963 | 0.865 | 2.00 | 0.131 | 0.408 | 0.697 |
| pe | missing | 0.7.1 | 0.900 | 0.967 | 0.869 | 1.54 | 0.093 | 0.295 | 0.704 |
| pe | missing | 0.7.2 | **0.906** | 0.977 | 0.869 | 1.08 | **0.088** | 0.293 | 0.700 |
| se | full | 0.7.0 | 0.874 | 0.975 | 0.830 | 2.88 | 0.145 | 0.467 | 0.682 |
| se | full | 0.7.1 | 0.888 | 0.981 | 0.844 | 1.92 | 0.107 | 0.331 | 0.696 |
| se | full | 0.7.2 | 0.888 | 0.982 | 0.843 | 1.75 | 0.107 | 0.333 | 0.700 |
| se | missing | 0.7.0 | 0.866 | 0.962 | 0.826 | 2.04 | 0.146 | 0.467 | 0.648 |
| se | missing | 0.7.1 | 0.881 | 0.973 | 0.836 | 1.42 | 0.113 | 0.361 | 0.648 |
| se | missing | 0.7.2 | 0.884 | 0.977 | 0.839 | 1.33 | 0.109 | 0.367 | 0.648 |

Paired differences per sample, mean (95% interval):

| reads | database | comparison | F1 | FP per sample | Bray-Curtis |
|---|---|---|---|---|---|
| pe | full | 0.7.0 less 0.6.0a | +0.269 (+0.175, +0.360) | +1.96 | −0.214 |
| pe | full | 0.7.1 less 0.7.0 | +0.0055 (+0.0023, +0.0087) | −0.85 | −0.044 |
| pe | full | 0.7.2 less 0.7.1 | **+0.0022 (+0.0002, +0.0046)** | −0.23 (−0.69, +0.19) | −0.0015 (−0.0029, −0.0002) |
| pe | missing | 0.7.2 less 0.7.1 | +0.0054 (−0.0010, +0.0138) | −0.46 (−0.96, −0.04) | −0.0048 (−0.0089, −0.0016) |
| se | full | 0.7.2 less 0.7.1 | +0.0001 (−0.0037, +0.0041) | −0.17 | −0.0005 |
| se | missing | 0.7.2 less 0.7.1 | +0.0034 (−0.0034, +0.0105) | −0.08 | −0.0042 (−0.0080, −0.0010) |

By depth and abundance (`results/short_and_old_long/summary.md`): at 1,000 read pairs every 0.7 version finds 37-48%
of the species at 0.1-1% (0.6.0a 0-3%), at 10,000 pairs 81-85% of those at 0.01-0.1% (0.6.0a 1-8%); 0.7.1 and 0.7.2
estimate their abundance alike (median |log2| 0.52 at 10,000 pairs, 2x150), and 2-4 times closer than 0.7.0 above
that (0.07 against 0.23 at 500,000 pairs). Archaea: 0.7.x finds 82-100% at 10,000 pairs or more, 0.6.0a 3-85%. False
positives per sample by depth range from 0 to 7.5 (0.7.2: 0.8 to 3.2 with the full database), most of them congeners
of a species in the sample.

## Speed and memory

Means per run (full database; the training database alike), `-t 6`, wall / CPU / peak RSS:

| reads | 0.6.0a | 0.7.0 | 0.7.1 | 0.7.2 |
|---|---|---|---|---|
| pe, all 26 samples | 14.7 s / 38 s / 3.32 GB | 7.6 s / 25 s / 3.48 GB | 7.1 s / 24 s / 3.44 GB | 6.7 s / 23 s / 3.44 GB |
| pe, 2 x 5M pairs (2x150) | 86.4 s / 298 s / 3.39 GB | 54.3 s / 191 s / 3.54 GB | 49.3 s / 176 s / 3.48 GB | 52.6 s / 179 s / 3.49 GB |
| se, 24 samples | - | 2.2 s / 6 s | 2.2 s / 6 s | 2.2 s / 7 s |
| pb, 8 samples | - | 6.0 s / 21 s | 6.9 s / 23 s | 6.9 s / 23 s |
| ont, 8 samples | - | 5.7 s / 28 s | 6.0 s / 29 s | 5.3 s / 27 s |

0.7.2 less 0.7.1 per paired-end sample: −0.4 s (−1.5, +1.0) with the full database, −1.2 s (−2.4, −0.3) with the
training database. The 5M-pair samples took 3 s longer with 0.7.2 at about the same CPU time, presumably the load.
0.7.2's pipeline took 21.5 min against
24.3 min for 0.7.1's, the same design (peak 4.9 GB).

## Strains

`scripts/make_strains.py`: 12 species of the databases (10 bacterial, 2 archaeal), each grown into 8 strains from its
representative (the databases' reference) along a random coalescent tree, root-to-tip 0.4-1.2% substitutions per site
at a gene of rate 1, each marker gene at its rate in the world (ribosomal proteins slowest), substitutions only, so
that every strain SNP against the reference is known. Sample j holds strain j of every species (MSAs of 8 rows, whose
true tree is the strain tree) at the species' coverage (2 to 40x, log-uniform: 4 species at 2-6x) times 0.4-1.6 per
sample, and 15 background species at 0.5-5x. Paired-end 2x150 (ART HSXt); PacBio HiFi-like (pbsim3 errhmm, 15 kb,
99.9%) and Nanopore (qshmm, 8 kb, 97%) from the same genomes. Each version profiles the 8 samples in one run, with its
qcmsa, against its finished database.

`scripts/strain_score.py`: SNP calls from the raw MSA (qcmsa drops columns, so its columns no longer map to the genes'
positions): at each gene position, a sample's base is a called SNP if it is A/C/G/T and differs from the reference;
TP if it is the strain's base, FP (miscalled) otherwise, FN a strain SNP not called so. Phylogeny: p-distances of the
sample rows, neighbour joining, its splits against the true tree's (Robinson-Foulds, normalised), and the Spearman
correlation of the p-distances with the true patristic distances.

| reads | version | species with an MSA | MSA columns (raw / qcmsa) | SNP precision | SNP recall | SNP F1 | miscalled SNPs per 100 kb | true topology | RF | Spearman |
|---|---|---|---|---|---|---|---|---|---|---|
| pe | 0.6.0a | 12 of 12 | 114,125 / 108,841 | 1.000 | 0.768 | 0.869 | 0.0 | 4 of 12 | 0.154 | 0.826 |
| pe | 0.7.0 | 12 of 12 | 114,640 / 114,416 | 0.998 | 0.937 | 0.967 | 2.0 | 10 of 12 | 0.033 | 0.922 |
| pe | 0.7.1 | 12 of 12 | 114,695 / 114,695 | 0.998 | 0.939 | **0.967** | 1.9 | **11 of 12** | 0.017 | 0.922 |
| pe | 0.7.2 | 12 of 12 | 114,695 / 114,695 | 0.998 | 0.939 | **0.967** | 1.9 | **11 of 12** | 0.017 | 0.922 |
| pb | 0.7.x | 12 of 12 | 115,220 / 114,112 | 0.953 | 0.868 | 0.908 | 45.8 | 6 of 12 | 0.117 | 0.903 |
| ont | 0.7.x | 12 of 12 | 114,681 / 114,110 | 0.994 | 0.889 | 0.939 | 5.7 | 8 of 12 | 0.083 | 0.926 |

The MSA lengths are the 12 species' mean (52 to 119 genes each). 0.6.0a calls no SNP wrongly but misses a quarter of
them, most at the low-coverage species (recall 0.35-0.45 at 2-3x, where the 0.7 versions have 0.79-0.84; from 9x on
0.6.0a has 0.88-0.998, the 0.7 versions 0.98-1.00). Its trees are wrong for 8 of the 12 species, 3 of them at 37-39x
where it calls 99% of the SNPs; why was not traced. One species (Damviia corcomia, 39x), whose tree has an internal
branch of a few SNPs, no version gets right; Vitultulia yaldamvenum (39x) only 0.7.1 and 0.7.2
(`results/strains/species.tsv`). The three 0.7 versions give the same long-read results (0.7.0 for PacBio differs in
the fourth digit). Nanopore calls SNPs nearly as well as paired-end reads; PacBio has 24 times as many miscalled SNPs
as paired-end reads (not traced).

Two things had to be fixed for this to measure anything:
- **0.6.0a writes no MSA unless its strains folder exists** (it does not create it; `[qcmsa] WARNING: missing
  MSA/partition/meta`): `strain_runs.sh` creates the folder before every run. 0.6.0a's partition files are 0-based,
  0.7's 1-based.
- **pbsim3's errhmm model gives every PacBio base quality 0** (`!`), which every SNP filter rejects: no PacBio SNP
  was called. The PacBio strain reads get Q30 (their accuracy) instead (`make_strains.py --pb-qualities`). The
  benchmark's and the pipelines' PacBio samples have Q0 as well, which is why `excess_median` is about −0.97 there.

## Long reads, and the samples they are judged on

0.7.1's collector had pbsim3 sample reads once per genome, which cut each contig's last read to the contig's share of
the bases: the benchmark's 3 Mb samples are half reads under 1 kb (PacBio mean 2.3 kb, Nanopore 1.7 kb, against 15
and 8 kb meant), and so were 0.7.0's and 0.7.1's training samples. 0.7.2's collector draws each read's genome, length
and start itself (`c40cec3`), so its models are trained on reads of the intended length (PacBio mean 14 kb). On the
old samples 0.7.2's PacBio model is the worst of the three (F1 0.953 at `--knob 0.5` against 0.969 for 0.7.1, full
database; `results/short_and_old_long/overall.md`): its features meet reads it was not trained on.

`scripts/long_reads.sh` makes the samples again with 0.7.2's collector (the same design and seed; 3 Mb now holds
about 200 PacBio or 375 Nanopore reads, hence the lower F1 of every version) and runs the 0.7 versions on them, 0.7.2
also at
`--knob 0.5` (`results/long_reads_0.7.2_samples/`):

| reads | database | 0.7.0 | 0.7.1 | 0.7.2 (depth knobs) | 0.7.2, `--knob 0.5` |
|---|---|---|---|---|---|
| pb | full | 0.789 | 0.790 | 0.790 | **0.794** |
| pb | missing | 0.758 | 0.763 | 0.761 | **0.768** |
| ont | full | 0.865 | 0.870 | 0.872 | **0.873** |
| ont | missing | 0.843 | 0.854 | 0.845 | **0.860** |

Per sample: 0.7.2 at 0.5 less 0.7.1 PacBio +0.0039 (+0.0011, +0.0076) and +0.0051 (−0.0073, +0.0160), Nanopore +0.0029
(−0.0003, +0.0072) and +0.0060 (+0.0006, +0.0124); the depth knobs less 0.5 PacBio −0.0038 (−0.0105, +0.0006) and
−0.0073 (−0.0125, −0.0024), Nanopore −0.0011 (−0.0023, 0.0000) and −0.0148 (−0.0293, −0.0013), with more false
positives (0 to +1.25 per sample). 0.7.2's knobs (`results/pipeline_0.7.2/depth_knobs.txt`): PacBio 2: 0.27,
3: 0.21, 4: 0.60; Nanopore 2: 0.13, 3: 0.46, 4: 0.76, each bin on 2-6 training samples; on the pipeline's own test
set they left PacBio as it was (F1 0.981) and cost Nanopore 0.007 (0.983 to 0.976).

## 0.7.2's pipeline

`results/pipeline_0.7.2/summary.txt` against 0.7.1's (`../2026-10-01-v071-benchmark`): F1 of species held out pe
0.9835 (0.7.1 0.9791), se 0.9796 (0.9747), pb 0.9721 (0.9729), ont 0.9884 (0.9750); on its independent test set pe
0.9739 (0.9682), se 0.9690 (0.9642), pb 0.9809 (0.9658), ont 0.9827 (0.9503); at knob 0.5, the long-read test sets
made by each version's own collector.

## Not covered

Real samples, a GTDB-scale database, strain mixtures (one strain per species and sample here), and indels in the
strains. Timing on a loaded, shared machine.

## Reproduce

```
bash docs/claude/2026-10-02-v072-benchmark/scripts/profile.sh       # 0.7.2 built, its pipeline, all runs
python3 docs/claude/2026-10-02-v072-benchmark/scripts/make_strains.py  # strain samples (env python: numpy)
bash docs/claude/2026-10-02-v072-benchmark/scripts/strain_runs.sh
python3 docs/claude/2026-10-02-v072-benchmark/scripts/strain_score.py
bash docs/claude/2026-10-02-v072-benchmark/scripts/knob_runs.sh     # 0.7.2's long reads at --knob 0.5
bash docs/claude/2026-10-02-v072-benchmark/scripts/long_reads.sh    # long-read samples by 0.7.2's collector
python3 docs/claude/2026-10-02-v072-benchmark/scripts/score.py; python3 docs/claude/2026-10-02-v072-benchmark/scripts/long_reads_score.py
```
