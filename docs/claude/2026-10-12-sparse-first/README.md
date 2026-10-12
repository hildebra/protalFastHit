# Sparse-first lookups: an early exit from seeding (--sparse_first), benchmarked on synthetic data

- **Date**: 2026-10-12.
- **Question (the user's)**: implement the sparse-first lookup order the stage-1 probe suggested, and benchmark it
  locally on synthetic data; it is an early-fail mechanism for the seeding phase. Then: use it for PacBio at least
  (C = 8), if the sparse lookup can be made cheaper.
- **Code**: the working tree on `ccdbf52` (with the other sessions' uncommitted hunks of the time), not committed.
  `--sparse_first C`, `--sparse_first_score T` (14), `--sparse_first_hits M` (1); a development option, off by default.
- **Machine**: WSL2 on the Core Ultra 7 258V, cores 0-3, Release build, 4 threads unless said.

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

## For r226

The probe at r226 ([report](../2026-10-12-stage1-probe-r226.md)) gives the share of reads kept and the recall of each
rule: paired-end C = 2, T = 14: 9.4% of the pairs kept, 99.63% of those with records; HiFi C = 8: 17.9% kept, 99.17%.
What the sampled lookups cost there: a lookup without its entries and seeds is the key map (11% of a lookup) and the
flex cells (60%), and C = 2 makes half of a pair's lookups, so the seeding of the 90.6% of pairs given up falls to
~35% of theirs, and their anchors, k-mer screen and alignments (the screen refused 90% of r226's 176 M candidates) go.
Only a cluster run can say what that is in seconds: `scripts/measure_performance.sh --test_sparse` runs each sample at
four levels after its timed runs (C = 1, 2, 4, 8 for short reads, 4, 8, 16, 32 for long reads) and compares each with
them (`sparse.tsv`: times, reads given up, records, species lost and gained, Bray-Curtis):

```bash
sbatch --exclusive -c 32 --mem 100G --wrap "bash scripts/measure_performance.sh --threads 8 --test_sparse runs/sparse \
    DB.protal pe:R1.fq.gz:R2.fq.gz pb:hifi.fq.gz"
```

## Writing out the reads kept

The gate decides per read before any seeding, so the reads it keeps (9.4% of a r226 gut sample's pairs at C = 2, 18%
of the HiFi reads at C = 8) could be written to a file as they are kept: BGZF blocks (ISA-L level 1, as `.sam.gz`
records are compressed) from each thread's buffer, the two mates' files in step. About 1 GB for the r226 gut sample.
What it would serve: running protal again on the same database (other options, another model, strain MSAs) on a tenth
of the reads; the input of a two-stage design. What it would not: another database, whose index keeps other reads (a
read near-exact to a genome only the new release has). Not implemented.
