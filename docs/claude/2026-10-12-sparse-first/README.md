# Sparse-first lookups: an early exit from seeding (--sparse_first), benchmarked locally and at r226

- **Date**: 2026-10-12.
- **Question (the user's)**: implement the sparse-first lookup order the stage-1 probe suggested, and benchmark it
  locally on synthetic data; it is an early-fail mechanism for the seeding phase. Then: use it for PacBio at least
  (C = 8), if the sparse lookup can be made cheaper.
- **Code**: the local benchmarks ran on the working tree on `ccdbf52` (with the other sessions' uncommitted hunks of the
  time), committed as `e751c7f`; the r226 run on `5222386`. `--sparse_first C`, `--sparse_first_score T` (14),
  `--sparse_first_hits M` (1); a development option, off by default.
- **Machine**: WSL2 on the Core Ultra 7 258V, cores 0-3, Release build, 4 threads unless said; the r226 run on a
  cluster nodes ([At r226](#at-r226-cluster-runs)).

## What it does

A read (a pair as one; a long read over all its chunks) is looked up by a sample of its k-mers first: those whose
15-base core's hash falls in 1/C of the range (`stage1_probe::Level`, the probe's sampling). A sampled lookup is a hit
if its block's best flex cell shares T or more of the 16 flank bases with the read (a key of one value, which has no
flex cell, is scored from the gene's packed bases as its flex cell would be). With fewer than M hits the read is given
up before the rest of its seeding, its anchors, the k-mer screen and the alignment: it is written as a read without
candidates (no record, no failed candidates). A read kept is looked up and aligned as without the option.

- **Short reads** (`Classify.h`): seeded one read (pair) ahead, as the normal path is. The next pair's sampled cores'
  key map blocks are fetched while the pair before it aligns (`PrefetchSparseKeys`); at the top of a pair's turn its
  sampled lookups are made and their blocks fetched (`SparsePrepare`); the next pair is read and its k-mers taken; then
  the sampled lookups are scored (`SparseScore`, flex cells only: no entries, no seeds). A pair kept is then looked up in
  full (`FindSeeds`).
- **Long reads** (`LongReads.h`): the read's chunks' sampled lookups in read order (`SparseHits`), each core's key map
  block fetched 8 k-mers ahead, each block fetched as soon as the lookup is found and scored 4 lookups later, stopping
  at the M-th hit.
- Counts: `Sample X sparse first (1/C of the cores, score T, M hit(s)): N of R pairs given up (p%), L sampled lookups
  (l a pair)`; the seeding line then counts the kept reads' seeding only.
- Not with `--stage1_probe` (which measures what this would give up).

**Exact where it keeps a read**: with every core sampled (C = 1) and any score a hit (T = 0), only pairs without any
k-mer in the index are given up, and those have no seeds anyway: on the archaea world's 100,000 pairs the SAM records
(as a sorted set) and the profile were identical to a run without the option ([sf_exact.sh](scripts/sf_exact.sh)).

## Making the sampled lookup cheap

The first version was slower than no gate (archaea world, C = 1: aligning +65%), with fewer instructions: callgrind
counted 3.64 G instructions against plain's 4.15 G for 100,000 pairs on one thread. Two causes:

- **Waits**: it made a pair's sampled lookups and scored them at once. The normal path fetches a pair's blocks while
  the next pair is read; the gate's waited. Fixed by the split above (SparsePrepare before the next read, SparseScore
  after), and for long reads by the streaming order.
- **Keys of one value**: 45% of the archaea world's blocks (r226: 1.2%). The first version decoded the gene's window
  and compared strings, both orientations; now the 31-mer's 2-bit code is assembled from the packed gene, the strand
  of the smaller core taken (as the index keys a k-mer), and its flank bits scored against the read's as a flex cell is
  (`KmerLookupSM::FlexSimilarity`): ~67 instructions a call.

After both ([sf_callgrind.sh](scripts/sf_callgrind.sh), the paired-end loop on one thread, 100,000 pairs):

| arm | instructions | vs plain | the largest parts |
|---|---|---|---|
| plain | 4,152 M | | k-mers 943 M (23%), lookups (`LookupKmers`) 383 M, `GetFromLookup` 288 M, reverse complements 246 M |
| C = 1, T = 14 | 3,346 M | −19% | k-mers 943 M, `SingleValueScore` 579 M (17%: this world's single values), `SparsePrepare` 408 M |
| C = 2, T = 14 | 2,860 M | −31% | k-mers 943 M (33%), `SparsePrepare` 292 M, `SingleValueScore` 290 M |

Taking the k-mers (the closed syncmers) is the same in every arm: a read must be scanned to find its sampled cores.

## Results

### Paired-end, archaea world: speed

[sf_speed.sh](scripts/sf_speed.sh): every r226 archaeal species' marker genes (the
[seeding profile](../2026-10-12-seeding-profile/README.md)'s world), 1M pairs of whole-genome reads of 120 species with
strains and congeners (2.9% of the pairs have records, as in a real metagenome); aligning only, 4 threads, medians of 3
alternated runs.

| arm | aligning s | wall s | CPU s | pairs given up | sampled lookups | k-mers looked up in full | records |
|---|---|---|---|---|---|---|---|
| plain | 3.43 | 4.62 | 15.9 | | | 44.7 M | 35,692 |
| C = 1, T = 14 | 3.94 | 5.71 | 19.6 | 96.2% | 19.2 M | 1.69 M | 35,692 (0%) |
| C = 2, T = 14 | 2.93 (−15%) | 4.55 | 14.9 | 96.3% | 9.6 M | 1.64 M | 35,555 (−0.38%) |
| C = 4, T = 14 | 2.15 (−37%) | 3.75 | 11.2 | 96.4% | 4.8 M | 1.59 M | 35,013 (−1.90%) |

Locally a lookup is cheap: this world's index is 0.4 GB of values (43% of the k-mers' cores in it, 4 cells a block),
mostly in cache. At r226 a lookup takes ~160 ns, ~100 of it waiting for memory, and seeding is two-thirds of aligning
([seeding profile](../2026-10-12-seeding-profile/README.md)), so the same counts should pay much more there. C = 1
samples every lookup and is slower here in wall time though it runs 19% fewer instructions, probably from the packed
gene bytes each key of one value reads (a further cache miss for 45% of this world's lookups, 1.2% of r226's).

### Paired-end, small real world: accuracy

[sf_accuracy.sh](scripts/sf_accuracy.sh): the 45 paired-end test samples of `~/realworld_sp_s1` (a build of 6 GTDB
families, 530 real genomes; 15 samples each of 10,000, 100,000 and 1,000,000 pairs), simulated again with the recorded
commands and this simulator (new reads, new truth), profiled with that build's database and models, which were trained
without `--sparse_first`; 4 threads, a whole run of all 45 per arm.

| arm | TP | FP | FN | precision | recall | F1 | records | run wall / CPU s |
|---|---|---|---|---|---|---|---|---|
| plain | 1487 | 13 | 236 | 0.9913 | 0.8630 | 0.9227 | 848,074 | 34.0 / 122.0 |
| C = 1, T = 14 | 1490 | 18 | 233 | 0.9881 | 0.8648 | 0.9223 | 847,977 | 39.2 / 143.2 |
| C = 2, T = 14 | 1490 | 18 | 233 | 0.9881 | 0.8648 | 0.9223 | 844,698 | 36.1 / 131.6 |
| C = 4, T = 14 | 1491 | 18 | 232 | 0.9881 | 0.8654 | 0.9226 | 832,363 | 33.8 / 122.3 |

- **The calls hardly change**: F1 within 0.0004, with 3-4 more true and 5 more false species.
- **`failed_candidate_rate` changes a lot**: the reads given up no longer count as failed candidates. Its mean over
  the true taxa falls from 0.573 to 0.206 (C = 1), over the false ones from 0.923 to 0.757; the model, trained on the
  former, calls a little more readily. It is the short-read models' second feature, so a database to be run with
  `--sparse_first` should train its models with it (the build cannot pass protal options to its training runs yet).
- This world's database is 27 MB: no arm can be faster here, and the run includes the profiling.

### HiFi, small real world: accuracy

The 24 HiFi test samples of that build made again as its collector makes them ([sf_hifi_sim.sh](scripts/sf_hifi_sim.sh):
each replays a regenerated paired-end community, 12 of 3 Mb and 12 of 30 Mb, HiFi model 15,000 ± 3,000 bp):

| arm | TP | FP | FN | F1 | records | CPU s |
|---|---|---|---|---|---|---|
| plain | 497 | 1 | 364 | 0.7314 | 20,044 | 14.1 |
| C = 4, T = 14 | 497 | 1 | 364 | 0.7314 | 20,036 (−0.04%) | 11.3 |
| C = 8, T = 14 | 497 | 1 | 364 | 0.7314 | 20,020 (−0.12%) | 11.2 |
| C = 16, T = 14 | 496 | 1 | 365 | 0.7305 | 19,981 (−0.31%) | 10.9 |

The calls are the same up to C = 8 (C = 16 loses one species). These genomes are small and dense in marker genes:
58-62% of the reads were given up (r226's gut sample: 82% at C = 8). `failed_candidate_rate` falls here too (true taxa
0.731 → 0.539 at C = 8), without changing a call.

### HiFi, archaea world: speed

300 Mb of HiFi reads (20,027 reads of 15,000 ± 3,000 bp) of the archaea world's 1M-pair community, made by
[sf_hifi_sim.sh](scripts/sf_hifi_sim.sh); `--read_type pb`, aligning only, 4 threads, medians of 3 alternated runs:

| arm | aligning s | wall s | CPU s | reads given up | sampled lookups (a read) | k-mers looked up in full | records |
|---|---|---|---|---|---|---|---|
| plain | 4.97 | 6.12 | 22.1 | | | 55.7 M | 11,093 |
| C = 4, T = 14 | 2.77 (−44%) | 3.76 | 13.5 | 71.4% | 4.77 M (238) | 16.4 M | 11,090 (−0.03%) |
| C = 8, T = 14 | 2.54 (−49%) | 3.58 | 12.5 | 72.4% | 2.41 M (120) | 15.8 M | 11,083 (−0.09%) |
| C = 16, T = 14 | 2.21 (−56%) | 3.26 | 11.8 | 72.9% | 1.22 M (61) | 15.5 M | 11,047 (−0.41%) |

At C = 8 the candidate alignments tried fell from 750,891 to 221,016 (the k-mer screen's refusals from 705,046 to
195,848), the alignments made from 15,041 to 15,028. A long read makes hundreds of lookups, so even 1/8 of them find
its gene, and the reads given up are most of the run's seeding and screening. This holds with lookups as cheap as
here; at r226 (C = 8 gave up 82% of the HiFi reads with 99.2% of the reads with records kept) the lookups saved are
dearer.

Tables of every run: [results/](results/).

## At r226 (cluster runs)

The probe at r226 ([report](../2026-10-12-stage1-probe-r226.md)) predicted the share of reads kept and the recall of
each rule: paired-end C = 2, T = 14: 9.4% of the pairs kept, 99.63% of those with records; HiFi C = 8: 17.9% kept,
99.17%. The user ran `scripts/measure_performance.sh --test_sparse` twice on the same real samples, with protal
`5222386` (`e751c7f`, this option, and the next commit's profiling changes), the r226 v24 database
(`protal0.8.1_r226_v24`, models trained without the option; on NFS, as is the run folder) and 32 threads:

- **Run 1**: SLURM 24200987, node t512n10: 2 × Xeon E5-2683 v4 (Broadwell, 2.1 GHz, AVX2, 2 NUMA nodes), all 32 cores.
  Tables: [results/r226_runs.tsv](results/r226_runs.tsv) (timed runs), [results/r226_sparse.tsv](results/r226_sparse.tsv)
  (levels; the species columns recomputed with the fixed comparison, below),
  [results/r226_species_changes.tsv](results/r226_species_changes.tsv) (each species lost or gained).
- **Run 2**: SLURM 24208987, node q512n10: AMD EPYC 9634 (Zen 4, up to 3.7 GHz, AVX-512, 1 NUMA node), 32 of its 84
  cores; with an ONT sample. Only its printed summary is here.
- **Run 3**: SLURM 24209738, the same node and samples again; its times repeat run 2's (aligning within 2%; HiFi
  +13..+29% again). Tables: [results/r226_epyc_runs.tsv](results/r226_epyc_runs.tsv),
  [results/r226_epyc_sparse.tsv](results/r226_epyc_sparse.tsv),
  [results/r226_epyc_species_changes.tsv](results/r226_epyc_species_changes.tsv), and the damage table below
  ([results/r226_damage.tsv](results/r226_damage.tsv)). Its paired-end C = 2 run took 63.9 s of wall time with a
  normal alignment: writing the profile took 5.5 s (0.6 s otherwise) and the profiling stage 29.5 s, the NFS run
  folder's stall, not the gate.
- **Samples**: a real gut sample, 50.6 M Illumina pairs; a real HiFi run, 497,656 reads; a real soil ONT run (SIP),
  3.1 M reads; all `.fastq.gz`. Three timed runs each without the option (medians below), then one run at each level.

Both runs gave up the same reads and wrote the same records at each level: the gate is deterministic, and the machines
differ in time only. The probe was right: C = 2 gave up 90.6% of the pairs and lost 0.37% of the records (the probe:
9.4% kept, 99.63% recall); HiFi C = 8 gave up 82.1% of the reads (17.9% kept).

### Run 1, Broadwell

Paired-end (wall, aligning and CPU against the timed runs' medians; CPU = the whole run's user seconds; cycles from
perf; Bray-Curtis as written, and over the species alone, scaled to sum 1):

| arm | wall s | aligning s | user CPU s | cycles | pairs given up | records | species lost / gained | Bray-Curtis (species alone) |
|---|---|---|---|---|---|---|---|---|
| plain | 83.6 | 61.5 | 2,075 | 4.79 T | | 2,223,279 | | |
| C = 1 | 59.0 (−29%) | 36.2 (−41%) | 1,379 (−34%) | −31% | 90.1% | −0.01% | 3 / 2 | 0.0001 (0.0001) |
| C = 2 | 63.3 (−24%) | 41.2 (−33%) | 1,315 (−37%) | −45% | 90.6% | −0.37% | 3 / 2 | 0.0005 (0.0007) |
| C = 4 | 61.2 (−27%) | 39.6 (−36%) | 1,177 (−43%) | −51% | 91.3% | −2.13% | 3 / 5 | 0.0030 (0.0029) |
| C = 8 | 60.3 (−28%) | 39.8 (−35%) | 1,079 (−48%) | −56% | 91.9% | −7.10% | 9 / 8 | 0.0189 (0.0079) |

HiFi:

| arm | wall s | aligning s | user CPU s | cycles | reads given up | records | species lost / gained | Bray-Curtis (species alone) |
|---|---|---|---|---|---|---|---|---|
| plain | 25.7 | 15.6 | 672 | 1.67 T | | 103,403 | | |
| C = 4 | 21.8 (−15%) | 10.8 (−31%) | 444 (−34%) | −43% | 78.8% | −0.19% | 0 / 2 | 0.0002 (0.0004) |
| C = 8 | 21.2 (−17%) | 10.8 (−31%) | 415 (−38%) | −47% | 82.1% | −0.54% | 0 / 2 | 0.0006 (0.0006) |
| C = 16 | 21.2 (−18%) | 10.7 (−32%) | 398 (−41%) | −50% | 83.9% | −1.34% | 0 / 2 | 0.0021 (0.0026) |
| C = 32 | 21.3 (−17%) | 10.7 (−32%) | 386 (−43%) | −51% | 85.5% | −3.59% | 1 / 3 | 0.0048 (0.0043) |

- **CPU falls with C, wall time does not**: from C = 1 to C = 8 the paired-end run's cycles fall from −31% to −56% (IPC
  1.46 → 1.74: the lookups left out were the memory-bound part), but the aligning time stays at 36-41 s; HiFi alike
  (cycles −43% to −51%, aligning 10.7-10.8 s at every level). At 32 threads the alignment has hit a floor the gate
  does not touch. The likeliest one: the input. protal inflates a gzip file that is not BGZF on one thread
  (`ThreadedGzStreambuf`; BGZF on several), and all 32 threads take their reads from those streams under one lock;
  whether these files are BGZF the logs do not say (protal tells it with `--verbose` only). In the timed runs a thread
  already spent 8.1 s of the 61.5 s getting reads (paired-end; HiFi 1.0 s of 15.6 s). C = 1's aligning time (36.2 s,
  the shortest though it does the most work) is one run each, so it may be noise around that floor.
- **Species**: every species lost or gained is near the detection limit (relative abundance at most 4.2e-5; the plain
  profile's median is 1.2e-4). At C = 1, which loses 0.01% of the records, three species are lost and two gained all
  the same: the models, trained without the option, see other `failed_candidate_rate` values (the reads given up
  count no failed candidates), and taxa at the edge flip. The unknown share hardly moves (paired-end 0.4635 → 0.4639 at
  C = 2, 0.4830 at C = 8, where half the Bray-Curtis is that shift; HiFi 0.2577 → 0.2578 at C = 8).

### Run 2, EPYC, with ONT

From the run's printed summary (changes against the timed runs' medians; Bray-Curtis as written):

| sample | arm | wall s | aligning s | given up | records | Bray-Curtis |
|---|---|---|---|---|---|---|
| paired-end | plain | 35.6 | 23.6 | | 2,223,279 | |
| | C = 1 | 34.7 (−2%) | 23.1 (−2%) | 90.1% | −0.01% | 0.0001 |
| | C = 2 | 34.4 (−3%) | 23.4 (−1%) | 90.6% | −0.37% | 0.0005 |
| | C = 4 | 35.1 (−1%) | 23.8 (+1%) | 91.3% | −2.13% | 0.0030 |
| | C = 8 | 35.2 (−1%) | 24.5 (+4%) | 91.9% | −7.10% | 0.0189 |
| HiFi | plain | 12.6 | 7.75 | | 103,403 | |
| | C = 4 | 13.9 (+10%) | 8.8 (+14%) | 78.8% | −0.19% | 0.0002 |
| | C = 8 | 14.3 (+14%) | 9.5 (+22%) | 82.1% | −0.54% | 0.0006 |
| | C = 16 | 14.9 (+18%) | 9.9 (+27%) | 83.9% | −1.34% | 0.0021 |
| | C = 32 | 14.8 (+18%) | 10.0 (+29%) | 85.5% | −3.59% | 0.0048 |
| ONT | plain | 50.9 | 44.0 | | 110,185 | |
| | C = 4 | 42.7 (−16%) | 36.1 (−18%) | 82.3% | −0.23% | 0.1447 |
| | C = 8 | 45.0 (−12%) | 38.1 (−14%) | 87.8% | −0.78% | 0.1447 |
| | C = 16 | 47.1 (−7%) | 39.4 (−11%) | 91.3% | −2.30% | 0.1514 |
| | C = 32 | 48.9 (−4%) | 41.1 (−7%) | 93.8% | −6.80% | 0.1454 |

- **Paired-end gains nothing here**: the plain run aligns in 23.6 s, and no level is faster. Run 1's levels stalled at
  36-41 s; scaled by the clock (2.1 GHz → up to 3.7 GHz) that floor would be ~21-23 s, so on this node the plain run
  is already at it. That fits a floor set by one thread (the input's inflating, the reader's lock), which the gate
  cannot move, but it is not measured.
- **HiFi is slower with the gate** (+14% to +29% aligning), the more so the higher C; on run 1 it was 31% faster.
  Unexplained. The gate's own cost should fall with C (fewer sampled lookups) while the reads given up rise; that the
  time rises with C instead points at something per read given up that costs more than the plain path's work on that
  read on this machine, or at the reader (threads done with their reads sooner queue on the input's lock more). The
  gate's timer and the reader's time per thread (below) tell the two apart.
- **ONT gains 18% at C = 4**, less the higher C (−7% at C = 32), the same rising pattern. Its plain run is not at a
  floor (the alignment handler 41.4 s and the indel bound 22.6 s a thread of 44.0 s).
- **ONT's profile moves**: Bray-Curtis 0.145 at every level, already at C = 4, where 0.23% of the records are lost: one
  species (below).

### How much the gate changes the profiles

No truth for real samples, so the yardstick is the plain run's profile (protal is deterministic: both nodes wrote the
same plain profiles). Over the species alone (each profile's species scaled to sum 1, so a change of the unknown
share does not count): Bray-Curtis; the mass of the species lost (in the plain profile) and gained (in the level's);
Bray-Curtis over the species both call (scaled over them); the largest relative change of a shared species of 0.1%
or more; the unknown share (run 3, [results/r226_damage.tsv](results/r226_damage.tsv)):

| sample | C | records | Bray-Curtis | mass lost / gained | Bray-Curtis, shared | largest change, shared ≥ 0.1% | unknown share |
|---|---|---|---|---|---|---|---|
| paired-end | 1 | −0.01% | 0.0001 | 0 / 0 | 0.0001 | 0.2% | 0.464 → 0.464 |
| | 2 | −0.37% | 0.0007 | 0 / 0 | 0.0007 | 2.4% | 0.464 → 0.464 |
| | 4 | −2.13% | 0.0029 | 0 / 0 | 0.0029 | 3.8% | 0.464 → 0.466 |
| | 8 | −7.10% | 0.0079 | 0.0001 / 0.0001 | 0.0078 | 10.5% | 0.464 → 0.483 |
| HiFi | 4 | −0.19% | 0.0004 | 0 / 0.0001 | 0.0004 | 4.0% | 0.258 → 0.258 |
| | 8 | −0.54% | 0.0006 | 0 / 0.0001 | 0.0006 | 5.0% | 0.258 → 0.258 |
| | 16 | −1.34% | 0.0026 | 0 / 0.0001 | 0.0026 | 12.2% | 0.258 → 0.259 |
| | 32 | −3.59% | 0.0043 | 0 / 0.0001 | 0.0043 | 12.8% | 0.258 → 0.263 |
| ONT | 4 | −0.23% | 0.2056 | 0.0006 / 0.2056 | 0.0002 | 8.1% | 0.926 → 0.901 |
| | 8 | −0.78% | 0.2059 | 0.0001 / 0.2059 | 0.0011 | 8.8% | 0.926 → 0.901 |
| | 16 | −2.30% | 0.2173 | 0.0001 / 0.2173 | 0.0031 | 13.7% | 0.926 → 0.899 |
| | 32 | −6.80% | 0.2065 | 0.0019 / 0.2065 | 0.0201 | 25.1% | 0.926 → 0.903 |

- **What losing reads does** is small at low C and grows with it, as the records lost do: paired-end C ≤ 2 and HiFi
  C ≤ 8 move the species both call by Bray-Curtis ≤ 0.0007, no species of 0.1% or more by more than 5%, and lose or
  gain only species of 4.1e-5 or less (relative abundance as written; the plain paired-end profile's median is 1.2e-4); paired-end C = 8 moves a shared species by up to 10%
  and the unknown share by 0.02, HiFi C ≥ 16 by up to 13%.
- **ONT: one call flips.** JADGHX01 sp019683695 (Steroidobacteraceae; a genus without a name, as most of this soil's)
  is not in the plain profile and is called at every level, with 2.0% of the sample (20% of what is classified): the
  whole of ONT's Bray-Curtis. Everything else scales down by the same factor, and the shared species move little
  (Bray-Curtis 0.0002 at C = 4). With records nearly the same (−0.23%), the flip is the model's: the reads given up no
  longer count as failed candidates, the failed-candidate features fall, and a taxon the plain model refused passes.
  The plain run's log warns that 4 species with strong evidence (depth 1x or more, reads on 90% of their genes, best
  reads 98% identical or more) scored below the knob, and the C = 4 run's that 3 did; JADGHX01 is probably one of
  the 4, which `misc/unreported_species.tsv` of a kept plain run (`KEEP=1`) would show. Whether the plain profile or
  the gated one is right, these data cannot say; either way the gate changed a call worth 2% of the sample through
  the models, not through the reads.
- **For comparison**, the small real world's truth ([above](#paired-end-small-real-world-accuracy)): F1 0.9227 →
  0.9223-0.9226 for paired-end C = 1-4, HiFi unchanged to C = 8, with models trained without the gate.

### What to use, and the next run

**What to use**: no default, and not for HiFi until run 2's slowdown is understood. Where the alignment is
compute-bound (run 1, ONT on run 2) the gate cuts 18-41% of the aligning time; where a run is at the input floor
(paired-end on run 2) it saves CPU only, which pays when several samples share a node. Records: paired-end C ≤ 2,
HiFi C ≤ 8 and ONT C ≤ 8 lose under 1%, and at paired-end C ≤ 2 and HiFi C ≤ 8 the profiles barely move (above).
Before a default: models trained with the option; ONT's flipped call shows what models trained without it do.

**The script's bug**: its first version printed "lost 0, gained 0" for every level. An awk array element read before
its membership test is created by the read, so every species counted as in both profiles. Fixed; run 1's numbers above
are recomputed from its profiles.

**Next run**: `measure_performance.sh` now writes the alignment's CPU seconds, the cores it kept busy and a thread's
seconds getting reads (`aligning_cpu_s`, `aligning_cores`, `reader_s_per_thread`) to `runs.tsv` and `sparse.tsv`, the
gate's own seconds a thread (`gate_s_per_thread`, protal's new "Sparse first" timer), each level's stage timers
(`sparse_stages.tsv`), Bray-Curtis over the species alone and both unknown shares, and in `environment.txt` each read
file's format (BGZF, gzip, zstd, plain). protal's `misc/cpu.tsv` now has a row `index and genomes` before the first
sample's: before, that sample's `aligning` row held the index load and the genome preload. The timer and the row need
a protal with this report's second commit. Few cores busy and much reader time would confirm the input floor; if the
files are plain gzip, the same reads as BGZF (`zcat R1.fastq.gz | bgzip -@ 16 > R1.bgzf.fq.gz`) as further samples
would show what the gate gives without it:

```bash
sbatch --exclusive -c 32 --mem 100G --wrap "bash scripts/measure_performance.sh --test_sparse runs/sparse DB.protal \
    pe:R1.fq.gz:R2.fq.gz pe:R1.bgzf.fq.gz:R2.bgzf.fq.gz pb:hifi.fq.gz pb:hifi.bgzf.fq.gz ont:ont.fq.gz"
```

## Writing out the reads kept

The gate decides per read before any seeding, so the reads it keeps (9.4% of a r226 gut sample's pairs at C = 2, 18%
of the HiFi reads at C = 8) could be written to a file as they are kept: BGZF blocks (ISA-L level 1, as `.sam.gz`
records are compressed) from each thread's buffer, the two mates' files in step. About 1 GB for the r226 gut sample.
What it would serve: running protal again on the same database (other options, another model, strain MSAs) on a tenth
of the reads; the input of a two-stage design. What it would not: another database, whose index keeps other reads (a
read near-exact to a genome only the new release has). Not implemented.
