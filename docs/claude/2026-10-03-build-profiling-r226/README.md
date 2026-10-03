# Where the r226 database build spends its time now, and the easiest wins

2026-10-03, branch `audit-fixes` at `4cae2f4` (the scripts and simulator unchanged since `be35d15`). The question:
where does the whole `build_gtdb_database.py` workflow lose time, including the random-forest training? Should
`simulate_metagenomes` be parallelised, or run several simulations in one command? What else would make it faster?

Sources:

- **The r226 v5 run** (job 23897575, `27423c6`, `-t 64`, `--scratch` on the node's SSD, 84-CPU node). Its logs are in
  `local/v5` and `local/protal_r226_v5_logs.tgz`, which are git-ignored. The step times come from the logs' own
  clocks (`*_simulation.log`, `training_data.log`, `test_data.log`, `training_db_index.log`, the trainer's `time:`
  line) and from the files' modification times in the tarball. The console log (stdout) was not uploaded, so the
  steps before the conversion's end (start-up, genome table, gene neighbours) are not timed here.
- **Unit costs measured today** on the WSL laptop, on the padded real-sized genomes of the 2026-10-01 report
  (`~/bprof/b2/db/genomes.tsv`: 2,295 genomes, median 3.4 Mb, 78 contigs). The binaries were built at v0.7.3
  (`~/protal-fp/build`), with ART 2.5.8, pbsim3 3.0.5 and Python 3.12. The machine was shared with other sessions
  (load 1.6-18 during the runs), so CPU seconds are more reliable than wall seconds. Work data: WSL `~/bprof26`.

## Summary

**The build waits on the simulations, and within them on scheduling, not on the simulator's speed.** After the
held-out species are chosen, v5 took 1:55. The training database was ready after 20 minutes. The collector then
waited another ~40-44 minutes for the training data's simulations, which ran 1:00 in all (paired-end 15:30, then
long reads 44:23). Then came protal on the training samples (31:23), protal on the test set (9:29, run after the
training samples' run), the four models (78 s), parity, and `--add_model`.

The simulations are slow for three reasons, none of which is the per-command cost of `simulate_metagenomes`:

1. **The deep long-read samples are single tasks, queued last.** The long-read samples go to the pool in design
   order, shallow first. All 360 samples of the ten shallow points were done 16:49 into the long reads. The 12 deep
   samples (1.5 Gb ×4 per read type, 6 Gb ×2 per read type) then ran alone for another 27.5 minutes, on a node that
   was otherwise idle because both database builds had finished. One 6 Gb Nanopore sample is a single pbsim3 run
   at ~0.4 s per Mb (measured: ~39 min on the laptop).
2. **The long reads start only after every paired-end point is simulated** (15:30). They need only those points'
   manifests, and a point's manifest exists when that point is done.
3. **The deepest paired-end samples run on one thread each.** The three 10M-pair points (2 samples, 2 threads
   each) ended at 10:52, 11:42 and 15:30. The other 18 points were all done by 10:00. The collector also hands out
   only 42 of the 64 threads (33 for the test set).

**Easiest wins, ranked by v5 wall time saved** (all five implemented the same day, [below](#follow-up-2026-10-03-fixes-1-5-implemented)):

| # | Change | Where | Effort | Saves (v5, est.) |
|---|---|---|---|---|
| 1 | Long-read samples **largest first**, and deep samples **split into chunks** (e.g. ≤ 250 Mb of templates each), simulated in parallel and the gzip members concatenated | `collect_training_data.simulate_long`, `long_read_sample` | sort: 1 line; chunks: ~50 lines + tests | long reads 44 → ~10 min |
| 2 | Long reads **not behind the paired-end points**: start a long-read point as soon as its community points are simulated (one pool for both) | `collect_training_data.main` | small | ~15 min more; the simulations then end before the training database |
| 3 | The paired-end points get **all `--threads`**, and `simulate_metagenomes` simulates **a sample's genomes on threads** (ART per genome side by side, appended in order) | collector's `threads_of`; `MetagenomeSimulator::write_all_reads` | small; medium | the 15:30 tail → ~5 min (margin once 1 and 2 are in) |
| 4 | The training database **on `--scratch`** | `build_gtdb_database.py` | small | of its 15:36 build, 6:15 is writing `database.protal` to `/hpc-home`; ~4-5 min once 1 and 2 put the build on the critical path; faster loads too |
| 5 | **One protal run** for both collections (or the test set profiled as soon as it is simulated) | `build_gtdb_database.collect` | small | one index load (~3 min?), and the test set's 9:29 if it overlaps idle time |

With 1 and 2, the simulations of both collections (~15 min if 3 is in too) finish before the training database
(20 min). The post-holdout part of v5 would then be about 20 (database) + ~38 (one protal run) + ~8 (training,
parity, `--add_model`), so **~1:05-1:10 instead of 1:55**. The next critical step is protal on 1,008 samples. Its
per-sample times are in the collector's `profile_all/protal.log` on scratch, which was not uploaded ("Aligning reads
took" per sample).

**The random-forest training is not worth optimising.** The four models train in parallel in 78 s at most (pe:
fit 1.1 s, evaluation 19 s, feature-set study 24 s, capacity 19 s). The previous-procedure grid, 82% of the trainer
before, has been off by default since `d11381f`.

**On the question asked: `simulate_metagenomes` already runs the samples of one command in parallel** (`-t`, one
sample per thread, `write_all_reads`), and the collector already runs all the design points at once. Putting
several simulations in one command would save almost nothing: a command's own cost (reading a 17k-row genome table
that has lengths, designing the communities) is milliseconds. The time is spent per genome of a sample (85 ms) and
per read (~96 s per million 150 bp pairs on one thread). What is missing is parallelism *within* a sample, and the
same for the long-read samples (#1, #3). One command for all points would let one thread pool balance everything,
but the collector can do that with the commands it has. A decompressed-genome cache shared across samples would
save at most the 29 ms per genome that is not ART, and genomes recur only ~2.3 times per collection, so it is not
worth it.

## The v5 timeline (after the conversion)

| Clock | Step | Took | Note |
|---|---|---:|---|
| 09:46 | conversion done | (4:20 of steps in `convert.log`) | earlier steps not in the uploaded logs |
| 09:48 | held-out species chosen; both collections' simulations and the finished database's build start in the background | | simulations at `nice` 10, `-t 64` each |
| 09:53 | training database's files | ~5 min | `training_db.log` |
| 10:08 | training database built | 15:36 | preload 17 s, passes 49 s, uniqueness 3:15, gene conservation 1:37, write index 1:44, **write `database.protal` 6:15** |
| 10:09 | finished database built (background) | 15:57 | |
| ~10:03 | training data: paired-end simulated | 15:30 | 21 points, 42 of 64 threads handed out; 18 points by 10:00, then the 10M-pair points |
| ~10:16 | test set simulated | 10:34 + 17:41 | |
| ~10:21 | training data: the 360 shallow long-read samples done | 16:49 into the long reads | |
| ~10:48 | training data: the last long-read sample (6 Gb ONT) done | 44:23 | 12 deep samples alone for 27.5 min |
| 10:52 → 11:24 | training data: protal, 768 samples, one run | 31:23 | 4:00 in: sample 60; 12-18 min: the three 10M-pair samples; 25-29 min: the 6 Gb ONT sample; then ~2.4 min profiling |
| 11:24 → 11:34 | test set: protal, 240 samples | 9:29 | after the training run, although its samples were ready at ~10:16 |
| 11:36 → 11:37 | four models, in parallel | 78 s max | 13 threads each |
| 11:38 | parity (4 runs) | ~1 min | |
| 11:43 | gene congeners, `--add_model` (one rewrite), metadata | ~5 min | |

From the holdout to the end: **1:55**. The step that waits on others is the collection's wait for its simulations
(10:08 → ~10:50, ~42 min), followed by the two protal runs in series (41 min).

## Unit costs (laptop)

`scripts/sim_cost.sh`: one `simulate_metagenomes` run of 2 samples per depth, 1 thread, the build's design (20-200
species, strains 0.3,0.1, 2 archaea), 150 bp HSXt, every ART call timed:

| Read pairs per sample | Wall | CPU | Genomes (2 samples) | ART: calls, wall, CPU | Not ART (wall) |
|---:|---:|---:|---:|---|---:|
| 1,000 | 38.5 s | 31.8 s | 453 | 453, 25.2 s, 19.7 s | 13.3 s |
| 500,000 | 130.5 s | 98.3 s | 402 | 402, 76.8 s, 64.9 s | 53.7 s |

- **Per genome of a sample: ~85 ms**, of which 56 ms is ART (process start, profile, reading the genome) and 29 ms
  is the simulator decompressing the genome for ART (31 ms alone for a 3 MB genome) and appending its reads. At
  1,000 pairs the median genome gets 2-3 pairs (ART `-f` median 0.00022), so a shallow sample is ~19 s of fixed
  cost. On the HPC it was ~40 s per sample (12 samples of a 1,000-pair point in ~8 min on one thread).
- **Per read: ~96 s per million 150 bp pairs**, of which ~54 s is ART and ~42 s is the simulator reading ART's
  FASTQ back and compressing it as BGZF at level 6 (`Bgzf.h`, `kLevel = 6`). A 10M-pair sample of 250 bp is about
  20-25 min on one thread.

`scripts/long_cost.py`: the collector's `long_read_sample` with each step timed, for one community (209 genomes),
30 Mb:

| Type | Step | Wall | CPU |
|---|---|---:|---:|
| ONT | templates (Python, 4,639 reads) | 10.3 s | 10.3 s |
| ONT | pbsim3 `--strategy templ` | 11.8 s | 23.6 s (children) |
| ONT | renaming + gzip (Python) | 1.2 s | 1.2 s |
| PacBio | templates (2,570 reads) | 8.6 s | 8.6 s |
| PacBio | `hifi_reads.simulate` | 4.3 s | 4.2 s |

pbsim3 runs at ~0.4 s wall (0.8 CPU-s) per Mb and HiFi at ~0.14 s per Mb, so a 6 Gb ONT sample is ~40 min in one
pbsim3 run and a 6 Gb PacBio sample ~15 min. The templates are mostly a fixed cost per sample:
`scripts/gil_profile.py` (cProfile) shows that 62% of their time is zlib decompression and 27% is bytes
`split`/`upper`/`join`. **Every genome is decompressed twice** (424 decompressions for 209 genomes), once per
drawing round in `long_read_templates`.

The collector runs the long-read samples on a thread pool, which raises the question of the GIL.
`scripts/gil_check.py`, 4 PacBio samples of 30 Mb: one after the other 61.3 s, on 4 threads 20.9 s, in 4
processes 16.3 s. The threads mostly do run side by side (zlib and most numpy release the GIL), but processes are
22% faster at 4 workers already, and the share that holds the GIL limits 64 threads more. (A counter thread beside
one sample, `scripts/gil_share.py`, gave 0.03-0.74 depending on the machine's load; too noisy to quote.)

Estimated CPU of v5's training long reads (24.8 Gb each of PacBio and ONT, 372 samples): pbsim3 ~5.4 h, renaming
~0.3 h, HiFi ~1 h, templates ~1 h, so ~8 CPU h. That is ~7.5 min on 64 cores, against the 44 minutes it took.

## The changes in more detail

1. **Long reads: largest first, deep samples in chunks.** In `simulate_long`, submit the tasks sorted by
   `bases` (descending). For a sample above a chunk size, draw its templates once, write them into k files, run one
   pbsim3 (or `hifi_reads`) per file with seeds derived from the sample's, then concatenate the renamed gzip
   outputs (concatenated gzip members are valid gzip; protal reads them). Read names stay `g<genome>x_<n>` with n
   global. The samples change (other random streams), so the long-read keys change and those points are simulated
   again. Sorting alone is the 1-line part. It lets the 6 Gb samples start at once, so they end after ~40 min on
   the laptop's speed (on the HPC the v5 ONT one took ≤ 44 min). Only chunks remove the tail.
   Also worth doing while there: cache the contigs read in the first drawing round (the second round re-reads
   them; ~half the template time), and a `ProcessPoolExecutor` instead of threads (tasks are plain dicts and
   `long_read_sample` is a module-level function).
2. **Long reads beside the paired-end points.** A long-read unit's communities are the manifests of its paired-end
   points (`unit_communities`). Submit each long-read unit's samples when its community points finish, into one
   pool of `--jobs` slots that both kinds share. The paired-end and long-read simulations then overlap instead of
   adding up.
3. **Paired-end: all threads, and threads within a sample.** `threads_of` gives each point
   `round(threads × cost share)` but at least 1. The shallow points get 1 and the total stays well under
   `--threads` (42 of 64). Handing out the remainder to the costliest points uses the node. Within
   `simulate_metagenomes`, `write_all_reads` gives a sample one thread. Running the ART calls of a sample's genomes
   on its threads (each into its own temporary FASTQ, appended in assignment order) keeps the reads byte-identical
   and divides a deep sample's ~55% ART share. The BGZF compression (most of the other ~45%) can run on a block
   pipeline too, or at level 1, since these files are read once by protal and deleted (the bytes change, the reads
   do not; not measured here).
4. **The training database on scratch.** It is read only by the collections and the parity check, and it is
   written at zstd level 3 to `/hpc-home` (6:15 of its 15:36 build was writing `database.protal`). On the node's SSD
   the write and every later load are faster. `--scratch` already holds the samples. It needs at least the
   database's ~22 GB more space there.
5. **One protal run for both collections.** Once the simulations finish before the database, the test set's samples
   can go into the training run's map (`profile_all`), with the tables split afterwards. That saves one index load
   and the run's start-up. While the simulations still end late, profiling the test set as soon as it is simulated
   uses the idle node instead.

Not worth it at r226 now: the trainer (78 s); `--add_model` (one rewrite since `460c9a1`); the per-genome fixed cost
of shallow paired-end samples (an in-process Illumina simulator for genomes of a few pairs would cut the 85 ms, but
those points run in parallel and end well before the deep ones); merging simulator commands (above).

Not measured, and worth a look in the next run's logs: the steps before the holdout (v1: start-up 4:32, genome
table 2:43, gene neighbours 15:39 before the 9.7× fix; the console output has the step lines), and protal's
per-sample alignment times at r226 (`profile_all/protal.log` on scratch).

## Follow-up (2026-10-03): fixes 1-5 implemented

On `4cae2f4` (uncommitted when measured): `src/RandomForest/MetagenomeSimulator.{cpp,h}`,
`src/simulate_metagenomes_main.cpp` (help text), `scripts/collect_training_data.py`, `scripts/build_gtdb_database.py`,
`scripts/mini_db/test_mini_db.py`, `docs/model-training.md`, `docs/simulation.md`, `docs/building-a-database.md`.

1. **Long reads, largest first and in chunks.** All the collector's simulations are now jobs in one queue
   (`Scheduler`), which runs the ready job of highest priority on the cores it needs (`--jobs`, default `-t`) and
   fills the cores left with smaller jobs. A long-read sample's priority is its estimated time (`long_read_seconds`:
   10 s plus 0.45 s per Mb for pbsim3, 0.15 for HiFi). A pbsim3 job takes 2 cores, because pbsim3 keeps about two
   busy. A sample above `--long_read_chunk` (250 Mb) is split into k = ceil(bases / 250 Mb) chunks
   (`long_read_chunks`). Each chunk draws a k-th of the bases with a seed of its own and names its reads every k-th
   from c + 1, so the names stay `g<genome>x_<n>`, unique, but not contiguous. The chunks' gzip files are then
   concatenated (`join_chunks`) into the sample's file, a gzip of several members. Samples of 250 Mb or less make
   the same reads as before. A chunked unit's key gets `chunk`, so only those points are simulated again.
2. **Long reads beside the paired-end points.** For each paired-end point a long-read point replays, a design
   run comes first (`design`: the point's own `simulate_metagenomes` command with `--test` into
   `points/<name>/design`, 0.03 s). The long-read unit's samples are made from that manifest
   (`long_unit_jobs(..., source)`), with truth paths into the point's `sim` folder. Once the point's reads are
   simulated, `same_design` checks that its communities are the design's, and the design folder is removed at the
   end. The paired-end points take threads in proportion to their estimated work (`pe_threads`), out of a share of
   the cores equal to their share of all the estimated work. The long reads get the rest, plus the cores the
   paired-end points free; the paired-end points are queued before any long-read sample.
3. **Threads within a paired-end sample.** `write_reads` runs a sample's ART calls on up to `threads` threads,
   each genome in a folder of its own, at most 2 × threads genomes ahead. One thread appends the genomes in the
   sample's order, so the files are the same for any thread count. `write_all_reads` gives a sample the threads
   beyond the samples (the first samples take the remainder). `pe_threads` can give a point up to 16 threads per
   sample.
4. **The training database on `--scratch`** (`SCRATCH/training_db`). The build's `--scratch` help and
   building-a-database.md now say 175 GB.
5. **One protal run for both collections.** Once the training data's simulations (and the test set's) are done,
   `build_gtdb_database.py` runs the test collector with `--prepare_profiling`. That writes
   `test/profile_all/samples.map`, with its units' folders in `samples.map.units`, and stops. The training
   collector then runs with `--also_profile` on that map: one protal run of both collections' samples, after which
   it marks the test units profiled. The test collector then finds everything profiled and writes its tables.

**Same outputs** (`scripts/sim_identity.sh`; one or two samples of 20,000 pairs from the padded genomes):

- The new simulator's reads are byte-identical to the old binary's, at `-t 1`, at `-t 4` (2 samples, 2 genome
  threads each) and at `-t 6` (one sample, 6 genome threads).
- The design run's manifest and `protal.meta` equal the real run's. One sample on 6 threads took 2.4 s instead of
  15.0 s.

**A collection before and after** (`scripts/collect_ab.sh`). `--simulate_only`, 6 cores, the padded real-sized
genomes, pe, pb and ont. Paired-end: 1,000, 20,000, 500,000 (2 samples) and 2,000,000 (1 sample) pairs of 150 bp,
3 samples per point. Long reads: 0.3, 6, 30 and 600 Mb (1 sample). Old: the collector of `4cae2f4` and its
simulator. Load 1-5 from other sessions.

| | Wall | CPU (user) | Paired-end done | Long reads done |
|---|---:|---:|---:|---:|
| before | 6:05 | 821 s | 1:25 | 4:40 after the paired-end points (the 600 Mb ONT sample last, alone) |
| after | **3:11** | 880 s | 1:45 (sharing the cores) | 3:11, beside them |

- All 18 paired-end read files and all 16 long-read samples of 30 Mb or less have the same reads before and
  after. They are compared decompressed, because Python's gzip writes the time into the header.
- The two 600 Mb samples are 3 chunks each now. They have other reads of the same size (PacBio 600.0 Mb in 47,529
  reads, before 47,312; ONT 596.9 Mb in 85,325, before 85,546) and no read name twice.
- The 6% more CPU is the chunks' templates: each chunk reads its genomes again.

At r226's v5 design the same mechanisms apply at 64 cores. The 12 deep long-read samples become ~140 chunks
started first; the long reads no longer wait 15:30 for the paired-end points; the 10M-pair samples get up to 32
threads. The estimate above (simulations ending before the 20-minute training database, ~1:05-1:10 after the
holdout instead of 1:55) is what the next r226 run's `*_simulation.log`, `training_data.log` ("and N samples of
another collection in one protal run") and console will show.

**Tests.** `scripts/mini_db/test_mini_db.py`, 50 tests OK on `4cae2f4` plus the changes, built from
`git archive` + the changed files in WSL `~/bt26`:

- `test_long_read_replay` now also checks chunked samples: names unique and every 5th per chunk, the bases, reads
  of their genomes, the same reads on 1 and 3 slots, the chunk files gone.
- New `test_scheduler`: priority, backfill, dependencies, new jobs, a failure stopping the queue, `pe_threads`.
- `GtdbBuildTest` checks the training database in `SCRATCH/training_db` (not in OUTDIR), the test set's 2 samples
  in the training data's protal run (10 samples), and `test_data.log` writing tables without a protal run.

The website is not affected: none of this changes how protal itself is run.

## Reproducing

From this folder, in WSL (padded genomes of the 2026-10-01 report in `~/bprof/b2/db/genomes.tsv`;
`SIM=~/protal-fp/build/simulate_metagenomes`; pbsim3 from `~/micromamba/envs/protal-db-build`):

```
bash scripts/run_all.sh ~/bprof26          # sim_cost.sh (1,000 and 500,000 pairs), long_cost.py (30 Mb)
python3 scripts/gil_check.py ~/bprof26/sim/p1000/manifest.tsv ~/bprof26/gil --samples 4
python3 scripts/gil_profile.py ~/bprof26/sim/p1000/manifest.tsv ~/bprof26/gilp --top 14
python3 scripts/gil_share.py ~/bprof26/sim/p1000/manifest.tsv ~/bprof26/gil2      # noisy under load
# the follow-up (OLD: the simulator before; NEW: built from git archive + the changes in ~/bt26)
bash scripts/sim_identity.sh ~/bprof26/identity
bash scripts/collect_ab.sh ~/bprof26/ab 6     # OLD_SCRIPTS: the collector of 4cae2f4 (git show HEAD:scripts/...)
```
