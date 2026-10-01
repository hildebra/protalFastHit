# protal 0.6.0a, 0.7.0 and 0.7.1 on an operon world: detection, abundance, archaea, rare taxa, false positives

- **Date**: 2026-10-01.
- **Code**: 0.6.0a = the tag `0.6.0a` (`014f4a9`); 0.7.0 = `4b21427`; 0.7.1 = `1c11a00` (branch `audit-fixes`,
  Version 0.7.1). Each exported with `git archive` and built in Release in WSL (`scripts/build_versions.sh`).
- **Machine**: WSL2 Ubuntu 24.04, Intel Core Ultra 7 258V, 6 threads, 23 GB; shared with other sessions (the
  load average before each run is in `results/runs.tsv`).
- **World** (`scripts/world.sh`): 900 GTDB-like species (108 archaeal), 3 genomes each, strains 0.4-4% from
  their representative, congeneric species 3-12% apart at a gene of factor 1, genes of different conservation
  (`simulate_gtdb_release.py --gene_rates categories`) and, new today, the markers laid out in operon-like
  clusters (`--operons`); 135 species are in no database (`release_p`, the tuning study's `make_release_p.py`),
  so the release has 765.
- **Databases and models**: each 0.7 version's own `build_gtdb_database.py` on the release
  (`scripts/pipelines.sh`, the V2 design of 2026-09-30 at reduced size, seed 1): the finished database with its
  four trained models (pe, se, pb, ont; 0.7.1's with the four alignment features of today in its model, and its
  database with `gene_neighbours.tsv` and `gene_conservation.tsv`), and the training database, which leaves out
  290 of the 765 species (whole clades of every rank and 119 single species; the same in both pipelines). 0.6.0a:
  its own builds of the same converted release, full and without those species (`scripts/db060.sh`), with the
  model it ships.
- **Samples** (`scripts/samples.sh`, 0.7.1's collector, in the design of the pipelines' independent test sets):
  from all 900 species, 10-300 species per sample, Poisson-lognormal abundances of sigma 2, a second strain in
  half of the species and a third in a fifth, 6 archaeal species per sample; paired-end 2x150 (HiSeq X profile)
  and 2x100 (HiSeq 2500) at 1,000, 10,000 and 500,000 read pairs, 4 samples each, and 2 samples of 5M pairs
  2x150; single-end: the first reads; PacBio and Nanopore (pbsim3) of the same communities, 3 and 90 Mb, 4 each.
- **Runs** (`scripts/profile.sh`): every sample by every version against both databases, one timed run each, 6
  threads, no qcmsa; 0.6.0a profiles paired-end reads only. `scripts/margin_check.sh` profiles 0.7.0's SAMs again
  with the margin 0.08. `scripts/score.py` scores them (`results/runs.tsv`, `species.tsv.zst`, `summary.md`).
- **Scores**: at each model's knob (0.5). Species no database has, and against the missing-species database the
  species it leaves out, cannot be found: they are not counted as present, the others' true abundances are
  renormalised, and their reads that land on relatives make false positives.

## Summary

0.7.1 is the best of the three on every read type. With paired-end reads it detects species as 0.7.0 does
(F1 0.908 against 0.903) and 0.6.0a does not (0.634). It halves 0.7.0's abundance error and beats 0.6.0a's:
Bray-Curtis 0.087 against 0.131 and 0.345. It has fewer false positives than 0.7.0 (2.1 per sample against
3.0). The abundance gain is the depth margin of 0.08: 0.7.0's own alignments profiled with 0.08 give 0.087
too. The detection and false-positive gains are 0.7.1's binary and model together: 0.7.0 at 0.08 keeps
F1 0.902 and 3.5 false positives per sample, and 0.7.1's binary with 0.7.0's model makes 6.5.

Archaea: 0.7.x finds 92–93% of the archaeal species at 10,000 read pairs or more, 0.6.0a 35% (its shipped
model). Low-abundance taxa: at 10,000 pairs 0.7.x finds 80–85% of the species with 1–10 read pairs, 0.6.0a
1–8%. 0.7.1 reconstructs their abundance better (median |log2| 0.52 against 0.55 for 0.7.0). Above
that, 0.7.1's abundances are 2–4 times as close (0.07 against 0.23 at 500,000 pairs). Single-end, PacBio
and Nanopore reads: 0.7.1 improves F1 by 0.002–0.014 and Bray-Curtis by 0.002–0.038 over 0.7.0.

Means over all samples of a read type (paired-end: 26, the 4 depths of both read lengths and the deep ones; the
others 24, 8, 8), at the model's knob 0.5:

| reads | database | version | F1 | Bray-Curtis | FP per sample | precision | FP share of abundance | archaea found (10,000 pairs or more) |
|---|---|---|---|---|---|---|---|---|
| pe | full | 0.6.0a | 0.634 | 0.345 | 1.00 | 0.993 | 0.48% | 130/376 = 0.35 |
| pe | full | 0.7.0 | 0.903 | 0.131 | 2.96 | 0.976 | 0.17% | 345/376 = 0.92 |
| pe | full | 0.7.0, margin 0.08 | 0.902 | 0.087 | 3.50 | 0.973 | 0.37% | 345/376 = 0.92 |
| pe | full | 0.7.1 | **0.908** | **0.087** | 2.12 | 0.980 | 0.31% | 348/376 = 0.93 |
| pe | full | 0.7.1, 0.7.0's model | 0.888 | 0.087 | 6.54 | 0.946 | 0.35% | 346/376 = 0.92 |
| pe | missing | 0.6.0a | 0.628 | 0.357 | 1.50 | 0.981 | 1.42% | 58/182 = 0.32 |
| pe | missing | 0.7.0 | 0.895 | 0.131 | 2.00 | 0.963 | 0.55% | 162/182 = 0.89 |
| pe | missing | 0.7.0, margin 0.08 | 0.895 | 0.095 | 2.31 | 0.960 | 0.97% | 163/182 = 0.90 |
| pe | missing | 0.7.1 | **0.900** | **0.093** | 1.54 | 0.967 | 0.77% | 163/182 = 0.90 |
| pe | missing | 0.7.1, 0.7.0's model | 0.882 | 0.095 | 3.50 | 0.936 | 1.09% | 164/182 = 0.90 |
| se | full | 0.7.0 | 0.874 | 0.145 | 2.88 | 0.975 | 0.20% | 290/328 = 0.88 |
| se | full | 0.7.1 | **0.888** | **0.107** | 1.92 | 0.981 | 0.19% | 292/328 = 0.89 |
| se | missing | 0.7.0 | 0.866 | 0.146 | 2.04 | 0.962 | 0.65% | 136/159 = 0.86 |
| se | missing | 0.7.1 | **0.880** | **0.113** | 1.42 | 0.973 | 0.88% | 134/159 = 0.84 |
| pb | full | 0.7.0 | 0.967 | 0.150 | 1.62 | 0.989 | 0.15% | 185/196 = 0.94 |
| pb | full | 0.7.1 | **0.969** | **0.124** | 1.25 | 0.991 | 0.07% | 186/196 = 0.95 |
| pb | missing | 0.7.0 | 0.959 | 0.142 | 2.50 | 0.976 | 0.98% | 90/96 = 0.94 |
| pb | missing | 0.7.1 | **0.962** | **0.120** | 1.50 | 0.983 | 0.81% | 91/96 = 0.95 |
| ont | full | 0.7.0 | 0.977 | 0.133 | 2.38 | 0.987 | 0.29% | 182/196 = 0.93 |
| ont | full | 0.7.1 | **0.982** | **0.112** | 1.88 | 0.989 | 0.08% | 184/196 = 0.94 |
| ont | missing | 0.7.0 | 0.960 | 0.119 | 1.75 | 0.976 | 0.89% | 82/96 = 0.85 |
| ont | missing | 0.7.1 | **0.964** | **0.117** | 1.62 | 0.978 | 0.88% | 83/96 = 0.86 |

"Missing": the training database, without 290 of the 765 species (whole clades and single species), on top of
the 135 no database has. "0.7.0, margin 0.08": 0.7.0's SAMs profiled again with `--depth_identity_margin 0.08`
(`scripts/margin_check.sh`). All tables per scenario are in `results/summary.md`.

## Paired-end reads by depth

F1 / Bray-Curtis, 2x150, mean of the samples (4; 2 at 5M):

| read pairs | database | 0.6.0a | 0.7.0 | 0.7.1 | 0.7.1, 0.7.0's model |
|---|---|---|---|---|---|
| 1,000 | full | 0.242 / 0.588 | 0.749 / 0.192 | **0.755 / 0.155** | 0.749 / 0.157 |
| 10,000 | full | 0.602 / 0.444 | 0.949 / 0.094 | **0.954 / 0.063** | 0.932 / 0.063 |
| 500,000 | full | 0.964 / 0.040 | **0.984** / 0.076 | **0.984 / 0.028** | 0.941 / 0.028 |
| 5,000,000 | full | 0.986 / 0.043 | 0.980 / 0.090 | **0.992 / 0.030** | 0.976 / 0.030 |
| 1,000 | missing | 0.236 / 0.601 | 0.737 / 0.199 | **0.741 / 0.164** | 0.726 / 0.173 |
| 10,000 | missing | 0.595 / 0.420 | 0.932 / 0.086 | **0.942 / 0.068** | 0.906 / 0.070 |
| 500,000 | missing | 0.940 / 0.078 | 0.969 / 0.073 | **0.982 / 0.029** | 0.943 / 0.029 |
| 5,000,000 | missing | 0.974 / 0.059 | 0.991 / 0.090 | **0.993 / 0.032** | 0.989 / 0.032 |

2x100 shows the same (`results/summary.md`); at 500,000 pairs F1 0.956 / 0.980 / 0.985 and Bray-Curtis
0.108 / 0.105 / 0.041 for 0.6.0a / 0.7.0 / 0.7.1.

- **Abundance**: 0.7.0's error at depth (0.076–0.090) is the margin of 0.04 it ships with, which cuts strains
  ([2026-09-30 benchmark](../2026-09-30-v07-vs-v06/README.md)). With 0.08, 0.7.0's own SAMs give 0.7.1's
  abundances (0.028 at 500,000 pairs, 0.030 at 5M), closer than 0.6.0a's (0.040, 0.043).
- **Detection at low depth** is the trained models': 0.6.0a's shipped model calls few species below 500,000
  pairs (F1 0.24 at 1,000, 0.60 at 10,000). At 5M pairs it is close (0.986 against 0.992).
- **0.7.1's binary with 0.7.0's model** has the same abundances as 0.7.1 but detects worse than either. Its
  features (depth from the wider margin, the reads' new candidate counts) are not those the model was trained
  on. A model has to be retrained with its binary, as the pipeline does.

## Archaea

Archaeal species found, of those present that a database has (paired-end, 2x150, full database):

| read pairs | 0.6.0a | 0.7.0 | 0.7.1 |
|---|---|---|---|
| 1,000 | 0/81 = 0.00 | 29/81 = 0.36 | 31/81 = 0.38 |
| 10,000 | 4/115 = 0.03 | 94/115 = 0.82 | 95/115 = 0.83 |
| 500,000 | 29/46 = 0.63 | 45/46 = 0.98 | 46/46 = 1.00 |
| 5,000,000 | 41/48 = 0.85 | 48/48 = 1.00 | 48/48 = 1.00 |

0.6.0a's model misses archaea, as known since round 4; the models trained for this database find them about as
well as bacteria. Archaeal false positives are 0–5 per point (of 4 samples) for every version. Long reads: 0.7.1
finds 115/115 archaea at 90 Mb with PacBio and Nanopore, 71/81 and 69/81 at 3 Mb (0.7.0: 115, 114, 70, 68).

## Low-abundance taxa

Recall / median |log2(reported / true)| of the species found, by true relative abundance (paired-end 2x150,
full database). At 10,000 pairs, 0.01–0.1% is 1–10 read pairs; the simulator's truth lists no species below
0.01%.

| read pairs | abundance | n | 0.6.0a | 0.7.0 | 0.7.1 |
|---|---|---|---|---|---|
| 1,000 | 0.1–1% | 445 | 0.03 / 2.96 | 0.47 / 0.80 | 0.48 / 0.74 |
| 1,000 | ≥ 1% | 71 | 0.59 / 1.50 | 0.99 / 0.37 | 0.99 / 0.24 |
| 10,000 | 0.01–0.1% | 370 | 0.08 / 3.77 | 0.84 / 0.55 | 0.85 / 0.52 |
| 10,000 | 0.1–1% | 321 | 0.63 / 1.79 | 0.98 / 0.38 | 0.98 / 0.28 |
| 10,000 | ≥ 1% | 69 | 0.80 / 0.64 | 1.00 / 0.20 | 1.00 / 0.10 |
| 500,000 | 0.01–0.1% | 74 | 0.85 / 0.60 | 0.99 / 0.23 | 0.99 / 0.10 |
| 500,000 | 0.1–1% | 203 | 0.97 / 0.12 | 1.00 / 0.23 | 1.00 / 0.07 |
| 500,000 | ≥ 1% | 63 | 1.00 / 0.08 | 1.00 / 0.29 | 1.00 / 0.07 |
| 5,000,000 | 0.01–0.1% | 149 | 1.00 / 0.07 | 1.00 / 0.21 | 1.00 / 0.06 |

- With few reads, 0.7.x finds 84–85% of the species with 1–10 read pairs, 0.6.0a 8%.
- 0.7.0's abundance error does not shrink with depth (0.21–0.29 at 500,000 and 5M pairs): it is the margin's
  bias, not noise. 0.7.1's falls to 0.06–0.10. 0.6.0a is nearly as close at depth (0.07–0.12), for the species it reports.
- With species missing the picture is the same; recall at 10,000 pairs and 0.01–0.1% is 0.81 for 0.7.1, 0.80
  for 0.7.0, 0.03 for 0.6.0a.

## False positives

Per sample at the knob, mean over the point's samples, and how many of them are congeners of a species in
the sample (paired-end 2x150; `results/summary.md` has all points):

| read pairs | database | 0.6.0a | 0.7.0 | 0.7.1 | 0.7.1, 0.7.0's model |
|---|---|---|---|---|---|
| 10,000 | full | 0.2 (1 of 1) | 3.2 (11 of 13) | 2.2 (9 of 9) | 6.5 (24 of 26) |
| 500,000 | full | 2.2 (9 of 9) | 2.2 (7 of 9) | 2.2 (5 of 9) | 10.2 (33 of 41) |
| 5,000,000 | full | 1.5 (3 of 3) | 7.5 (8 of 15) | 3.0 (3 of 6) | 9.0 (13 of 18) |
| 10,000 | missing | 0.2 (1 of 1) | 4.0 (12 of 16) | 3.0 (9 of 12) | 6.2 (20 of 25) |
| 500,000 | missing | 3.5 (14 of 14) | 2.2 (7 of 9) | 1.5 (6 of 6) | 5.5 (20 of 22) |
| 5,000,000 | missing | 2.5 (5 of 5) | 0.5 (1 of 1) | 0.5 (1 of 1) | 2.0 (3 of 4) |

- Most false positives of every version are congeners of a species in the sample, often one no database has:
  its reads land on the relative.
- They carry little abundance with the 0.7 models (0.1–0.9% of the reported total, means per read type). With
  0.6.0a they carry up to 4.1% against the missing-species database at 500,000 pairs: it calls the relatives of
  missing species with their borrowed reads.
- 0.6.0a has the fewest false positives at low depth because it calls almost nothing there.

## Gene neighbours and operons

The world's markers lie in operon-like clusters, and 0.7.1's database has the gene neighbours its pipeline
found in the genomes (301,504 rules of 260 clades). 0.7.1 paired about 3% of fragments across two
neighbouring genes: 15,888 of 500,000 at 2x150 and 152,340 of 5M. It found 877 and 8,512 guided mates on
the gene next door. On Nanopore reads at 90 Mb it placed 117 genes where the neighbour rules put them. Their
share in the gains above is not separated here (`--no_gene_neighbours` was not run); the gene-neighbours run
of today found calls unchanged on an operon world ([report](../2026-10-01-gene-neighbours-run/README.md)).

## Speed and memory

Mean wall time / CPU / peak RSS per run, full database (6 threads):

| scenario | 0.6.0a | 0.7.0 | 0.7.1 |
|---|---|---|---|
| pe 2x150, 500,000 pairs | 9.5 s / 33 s / 3.38 GB | 6.1 s / 25 s / 3.54 GB | 5.0 s / 20 s / 3.49 GB |
| pe 2x150, 5M pairs | 64.2 s / 285 s / 3.39 GB | 42.8 s / 191 s / 3.54 GB | 40.8 s / 180 s / 3.49 GB |
| pe 2x100, 500,000 pairs | 8.0 s / 23 s / 3.37 GB | 4.9 s / 18 s / 3.56 GB | 4.4 s / 16 s / 3.50 GB |
| se 150 bp, 500,000 reads | – | 2.8 s / 11 s / 3.48 GB | 2.6 s / 10 s / 3.44 GB |
| PacBio, 90 Mb | – | 8.3 s / 43 s / 3.48 GB | 8.6 s / 44 s / 3.54 GB |
| Nanopore, 90 Mb | – | 12.7 s / 60 s / 3.52 GB | 14.0 s / 69 s / 3.61 GB |

- 0.7.1 is 5–18% faster than 0.7.0 on short reads (today's sampled timers, reverse complement, 2-bit syncmers
  and gene windows) and 1.6–1.9× faster than 0.6.0a. On long reads it is 4% (PacBio) and 10% (Nanopore)
  slower than 0.7.0; not analysed here.
- The machine was shared: the load average before the runs had a median of 3.7 (another session), so wall
  times vary more than usual; CPU times are steadier. The laptop ran on battery from about 16:55 to 17:10,
  during the end of 0.7.0's training pipeline and 0.6.0a's builds, before any run timed here.
- The pipelines (both training databases, the finished ones, 96 training and 48 test samples, four models):
  0.7.0 21 min, 0.7.1 24 min, with the sample simulation running alongside the latter.

## Not covered

Real samples, a GTDB-scale database, and depth beyond 5M pairs. The models are trained on this world and
scored on it with another design and seed; that flatters 0.7.x against 0.6.0a's shipped model, which was
trained on other databases (as in the [2026-09-30 benchmark](../2026-09-30-v07-vs-v06/README.md)).
