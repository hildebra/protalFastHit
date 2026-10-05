# Why the collector uses 3 of 64 cores since the scenarios: one GIL and every genome read per chunk

2026-10-05, branch `audit-fixes` at `22c445d` (`collect_training_data.py` and `scenarios.py` as committed in
`3910d98`). The question: the r226 build with the scenarios (v12, still running) is noticeably slower;
`collect_training_data.py` runs with 64 threads but uses ~300% CPU. Is it an I/O bottleneck? The user also pointed
to the v10/v11 logs for the build's timings.

Sources:

- **The r226 v10 and v11 runs** (`local/v10`, `local/protal_r226_v11_logs.tgz` and `..._training.tgz`, unpacked into
  `local/v11`, all git-ignored). v11: job 23954399, `735617c`, `-t 64`, `--scratch` on the node's SSD, genomes and
  OUTDIR on `/hpc-home`. Both are **before the scenarios** (`3910d98`), so they show the build without them. Times
  come from the logs' own clocks and the files' modification times.
- **v12** (the first build with the default scenarios) is still running; its logs are not here. What this report
  says about it is measured locally and projected, and the projection should be checked against its logs.
- **Local measurements** in WSL (Core Ultra 7 258V, 6 vCPUs, Python 3.12.3, numpy 1.26.4, pbsim3 3.0.5 from
  `~/micromamba/envs/protal-db-build`), on 400 synthetic genomes of real size (`scripts/make_genomes.py`: random
  sequence, median 3.4 Mb, ~78 contigs, 80 bases per line, gzip level 6, ~1.1 MB each). Work dir `~/colprof`. The
  machine was shared with other sessions (load 1-9), so CPU seconds are firmer than wall seconds.

## Summary

**It is the GIL first, and the volume of genome reads behind it.** The collector runs every simulation job on a
thread of one Python process (`Scheduler`, a `ThreadPoolExecutor`). That was fine while the jobs were mostly
`simulate_metagenomes` and pbsim3 (other processes). The scenarios' reads are mostly Python: a long-read or Ultima
sample's templates are drawn by `long_read_templates`, and a 6 Gb sample is 24 chunks of 250 Mb, **each of which
reads, decompresses and parses every genome of the community again**, once per drawing round. A soil community has
~3,700 genomes (2,383-2,912 species at the r226 pool, x1.4 for strains). Reading one genome costs 21 ms, of which
13.8 ms is zlib (releases the GIL) and 7.2 ms is the file read and the Python parse (holds it).

- Locally, the same jobs on 6 slots use **2.1-2.5 cores** and burn 30-45% more CPU than on one slot (GIL
  contention). In worker processes they use 4.8-5.3. A process whose threads hold the GIL most of the time cannot
  use more than ~1 core for that work, plus what zlib, numpy and gzip do beside it: **the 300%**.
- The default scenarios at r226 make **2.0 M genome reads** in the two collections (~2.3 TB of gzipped FASTA
  read from `/hpc-home`), against ~0.1 M for all of v11's long reads. The GIL-held part alone is **3.0 core-hours
  in the training collector and 2.0 in the test collector** (laptop core), each serial in its process, plus the
  contention. **The scenarios' drawn reads therefore take at least ~4 h** (training) and ~2.7 h (test), side by
  side. All of v11's simulations took 46 min.
- **I/O is the second limit, not the first.** At the current pace a collector reads at most ~140 genomes a second
  (one per 7.2 ms of GIL), ~150-300 MB/s for the two from `/hpc-home`. If the GIL were removed alone, 64 workers
  would ask ~3,000 genome reads/s (~3.5 GB/s) of the network file system, which would then be the limit. The
  thread-state check [below](#checking-the-running-v12) tells which one binds right now.

**The fix, ranked** (details [below](#fixes-ranked)):

| # | Change | Effect at r226 (both collections) | Effort |
|---|---|---|---|
| 1 | **One pass over a sample's genomes for all its chunks** (each genome read once per drawing round per sample, not per chunk): a draw job per sample, then the chunks' reads as jobs. Same reads, byte for byte (prototype checked) | genome reads 2.04 M → 0.23 M; GIL-held CPU 5.0 → 1.4 h; all CPU 31.7 → 21.2 h; reads from `/hpc-home` ~2.3 TB → ~260 GB | ~80 lines + tests |
| 2 | **Python jobs in worker processes** (`ProcessPoolExecutor`; the `Scheduler` stays the planner) | removes the ~3-core ceiling (locally 2.1-2.5 → 4.8-5.3 of 6 cores) | ~30 lines |
| 3 | The simulated genomes **on `--scratch`** (copied once), if the check shows D states | the remaining ~270 GB of genome reads from the node's SSD | small, in `build_gtdb_database.py` |
| 4 | Long-read priorities that count the community's genomes | a soil Ultima chunk is ~150 s, estimated 52 s | small |

With 1 and 2, the scenarios' drawn reads are ~21 laptop core-hours in both collections: **~20-40 min on the node**
instead of ≥4 h, with the longest single piece a soil sample's draw (~3.5 min) followed by its chunks. Fix 2
without fix 1 would turn the GIL limit into a network-FS limit. Fix 1 without fix 2 leaves 1.4 h of GIL-held work
(0.8 h in the training collector), most of it the Nanopore renaming loop.

Not the problem: `hifi_reads.py` (numpy releases the GIL: 3.0-3.8x on 6 threads, as fast as processes), the host's
paired-end fragments (0.02 h), the scenarios' paired-end reads (`simulate_metagenomes`, C++ threads), the trainer.

## The v10 and v11 timelines (before the scenarios)

v11 (`735617c`), clock times from the logs and file times:

| Clock | Step | Took | Note |
|---|---|---:|---|
| 08:25:30 | conversion done | | `convert.log` |
| 08:27:36 | held-out species chosen | | 8,668 species; 2,621 of the genome table's 7,998 to simulate |
| ~08:31:30 | both collections' simulations start (nice 10), the training database's files written | | |
| 08:56 | test set simulated | 24:19 | |
| 08:59:23 | finished database built (background) | | `index_and_package.log` |
| 09:00:25 | training database built | | `training_db_index.log` |
| 09:17 | training data simulated | 45:43 | **17 min after the database**: the long reads' tail |
| 09:17 → 09:48 | protal, 1,230 + 240 samples, one run | 31:22 | tail 13-24 min: the deep samples |
| 09:49 → 09:52 | four models, in parallel | ~2.5 min | |
| 09:55:04 | parity, gene congeners, metadata | | |

From the holdout to the end: **1:27:28** (v10: 1:35, simulations 41:39, protal 40:31).

In the training simulations of v11:

- The 33 paired-end points ran while both database builds (64 threads each, nice 0) held the node. The shallow ones
  (12 samples, 1,000-5,000 pairs, 1 thread) took 9-14 min instead of the ~2.5 min an idle core needs (~150 genomes
  x 12 samples x 85 ms): they ran on what the builds left (nice 10), as intended. The last, `rl250_p30000000` (5
  threads), ended at 28:47.
- The long reads were the tail: PacBio 1.5 Gb at 18:30, 6 Gb at 30:59, the 180 shallow PacBio samples at
  31:43-31:49; Nanopore 1.5 Gb at 36:39, 6 Gb at 42:42, the 144 shallow ones last, at 45:30-45:43. After the
  database (29 min in), only these jobs ran, on threads of the collector, and their Python part (templates,
  renaming) is the same GIL-bound work as the scenarios', at a smaller scale (~0.1 M genome reads).

## Measurements

### Unit costs, one core (`scripts/unit_costs.py`, `results/unit_costs_py312.txt`)

| Step | Cost | GIL |
|---|---|---|
| a genome read (3.36 MB of FASTA): file | 0.4 ms | held |
| ... `gzip.decompress` | 13.8 ms | released |
| ... parse (`split`, `join`, `upper` per record) | 6.8 ms | held |
| ... the same parse as one `translate` per record (same contigs) | 5.0 ms | held |
| template drawing per read, besides the genome reads: Ultima (300 bp) / HiFi (15 kb) / ONT (8 kb) | 3.7 / 26.9 / 16.6 µs | held |
| `hifi_reads.simulate`: Ultima flow model / HiFi | 223 / 112 ms per Mb | mostly released (numpy) |
| pbsim3 `--strategy templ` (ONT) | 591 ms CPU per Mb, in its own process (~2 cores) | none |
| pbsim3's reads renamed and re-gzipped (Python loop) | 31 ms per Mb (222 µs per read) | held, but for zlib |
| a host fragment (`host_pe_chunk`'s loop) | 2.5 µs | held |

A genome is read 1.3-1.7 times per sample or chunk (the second and third drawing rounds re-read the genomes that get
reads again). Within one chunk's templates, the genome reads are 95-98% of the time.

### The collector's jobs on 6 slots (`scripts/gil_scaling.py`, `results/scaling.jsonl`)

The collector's own `Scheduler`, `long_read_chunks`, `long_read_sample` and `join_chunks` (imported from `22c445d`),
2 Ultima samples of 60 Mb in chunks of 10 Mb, communities of 300 genomes (lognormal sigma 1.5): 12 chunks. Modes:
*threads* is the collector as it is; *processes* runs each job in a worker process; *onepass* draws a sample's
templates for all its chunks in one pass over its genomes, then the chunks' reads; *onepass+processes* both. Two
alternated repetitions, load 3-9 from other sessions:

| Mode | Slots | Wall | CPU (all processes) | Cores used |
|---|---:|---:|---:|---:|
| threads | 1 | 91.2 s | 91.0 s | 1.00 |
| threads | 6 | 46.7 / 62.2 s | 118.1 / 133.4 s | 2.53 / 2.14 |
| processes | 6 | 28.9 / 29.6 s | 138.0 / 155.7 s | 4.78 / 5.27 |
| onepass | 6 | 15.5 / 26.2 s | 40.5 / 55.1 s | 2.62 / 2.10 |
| onepass+processes | 6 | 14.4 / 20.1 s | 45.9 / 62.5 s | 3.19 / 3.10 |

One ONT sample of 60 Mb (6 chunks, pbsim3 on 2 cores each): threads 40.5 s (49.6 s CPU in the collector, 48.9 s in
pbsim3), processes 25.0 s, onepass 20.7 s, onepass+processes 26.2 s (one sample: its draw is one job, then six
pbsim3 runs; noisy).

- On threads, 6 slots give 1.5-2x, not 6x, and the CPU rises 30-45% over one slot: threads waiting for the GIL
  and handing it over. On the node, 64 slots cannot do better than the same ~1 core of GIL-held work plus what runs
  beside it, which is the 300% the build shows.
- In processes the same jobs use 4.8-5.3 of the 6 cores (the rest was taken by other sessions).
- One pass needs a third to a half of the CPU (the genomes are read 12 times less here). With two samples it is
  two draw jobs first, so it cannot fill 6 slots at this size; at r226 there are 5 samples per scenario and read
  type, 50 draw jobs in all.

**Same reads.** `gil_scaling.py --check`: the onepass prototype's reads are byte-identical (decompressed) to the
collector's, for Ultima and ONT (one sample of 12 Mb in 4 chunks, 120 genomes), and the processes and
onepass+processes paths give identical PacBio reads. Each chunk keeps its own random stream, consumed in the same
order: its round's plan, then its draws genome by genome in sorted order; one pass only visits each genome once for
all chunks that need it in that round.

### `hifi_reads.py` on threads (`scripts/hifi_scaling.py`, `results/hifi_scaling.txt`)

6 template files of 20 Mb, one after the other, on 6 threads, in 6 processes: Ultima 3.76x and 2.37x on threads,
2.76x and 2.44x in processes; HiFi 3.14x and 3.00x on threads, 2.06x and 2.35x in processes. numpy releases the GIL
for most of `mutate`; threads are not what limits it.

## The scenarios at r226 (`scripts/scenario_load.py`, `results/scenario_load.txt`)

The default scenarios (`gut`, `soil`, `soil_shallow`, `host`), 3 hold-in samples each in the training collection
and 2 hold-out in the test collection, with the v11 pool (5,377 species the training database has, 2,621 it lacks):
soil and soil_shallow scale to 2,383-2,912 species (`fit_species`), gut keeps 350-450, host 2-50. The drawing is
modelled after `long_read_templates` (vectorised; the genome reads per chunk and round counted), the costs are the
unit costs above (laptop core):

| Unit (per collection: training / test) | Chunks | Genomes | Genome reads now | One pass | GIL-held CPU h now | One pass | All CPU h now | One pass |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| soil Ultima 20M reads (6 Gb) | 72 / 48 | 3,706 | 320,064 / 208,481 | 19,840 / 12,949 | 0.70 / 0.46 | 0.10 / 0.07 | 3.04 / 2.00 | 1.29 / 0.86 |
| soil PacBio 6 Gb | 72 / 48 | 3,706 | 272,963 / 178,808 | 29,901 / 20,078 | 0.55 / 0.36 | 0.07 / 0.05 | 2.16 / 1.42 | 0.74 / 0.50 |
| soil Nanopore 6 Gb | 72 / 48 | 3,706 | 324,897 / 206,089 | 30,933 / 19,699 | 0.82 / 0.52 | 0.23 / 0.15 | 5.02 / 3.28 | 3.30 / 2.20 |
| soil_shallow PacBio / Nanopore 1.5 Gb | 18 + 18 / 12 + 12 | 3,706 | 144,132 / 94,107 | 43,003 / 28,299 | 0.33 / 0.22 | 0.13 / 0.09 | 1.76 / 1.16 | 1.17 / 0.78 |
| gut PacBio / Nanopore 6 Gb | 72 + 72 / 48 + 48 | 560 | 169,436 / 105,069 | 12,907 / 8,201 | 0.51 / 0.32 | 0.20 / 0.13 | 4.68 / 3.07 | 3.77 / 2.50 |
| host Ultima, PacBio, Nanopore 3 Gb (90% host) | 108 / 72 | 36 | 9,811 / 6,285 | 1,146 / 726 | 0.13 / 0.09 | 0.10 / 0.06 | 2.47 / 1.65 | 2.42 / 1.62 |
| **total** | **504 / 336** | | **1,241,303 / 798,839** | **137,730 / 89,952** | **3.04 / 1.97** | **0.83 / 0.55** | **19.13 / 12.59** | **12.70 / 8.45** |

Of the one-pass CPU, 11.3 / 7.6 h are the reads themselves (pbsim3 for Nanopore, numpy for the others); the
remaining GIL-held time is mostly the Nanopore renaming loop (49.5 Gb x 31 ms per Mb = 0.43 h in training, an
upper bound since zlib runs beside it). The host's 27 M + 18 M paired-end host fragments add 0.02 h of GIL-held
Python and ~0.7 h of ART.

The scenarios' paired-end reads are `simulate_metagenomes` runs (C++, threads): 165 M pairs in training and 110 M
in test, ~4.4 and ~2.9 core-hours at 96 µs per pair, plus ~85 ms per genome and sample (5 min per soil sample).
They are not GIL-bound and fit in the same window once the drawn reads are.

**What v12 should show**, if this is right: in `training_data_simulation.log`, the `sc_*` long-read and Ultima
units ending hours after the design's points, the training collection ≥ ~4 h and the test collection ≥ ~2.7 h after
the start (more if a node core is slower than the laptop's for one thread).

## Checking the running v12

On the node, for each collector (`R` running, `S` sleeping, here mostly on the GIL's lock, `D` waiting on I/O):

```bash
for p in $(pgrep -f collect_training_data.py); do echo "== $p"; grep -h '^State' /proc/$p/task/*/status | sort | uniq -c; done
```

- **GIL-bound**: 1-3 threads `R`, the rest `S`, about 300% in `top`.
- **I/O-bound**: many threads `D`. Then `cat /proc/<pid>/io` twice a minute apart: `rchar` grows by GB per minute
  (genome reads; on a network file system `read_bytes` may stay 0).
- Optional, if `py-spy` can be installed (`pip install --user py-spy`): `py-spy dump --pid <pid>` shows every
  thread's stack (expected: most in `read_contigs` or `long_read_templates`), `py-spy top --gil --pid <pid>` the
  functions holding the GIL.

## Fixes, ranked

1. **One pass over a sample's genomes for all its chunks.** `templates_once` in `scripts/gil_scaling.py` is the
   prototype: per round, every chunk that still needs bases plans its reads with its own random stream, then the
   genomes any chunk drew are visited once in sorted order and each chunk draws its reads from them. The same
   reads as now (checked), so no long-read key changes. In the collector, `long_unit_jobs` would add one draw job per
   sample (Python, ~1 core, priority by genomes x 21 ms + reads x per-read cost), then one job per chunk for its
   reads (`hifi_reads.simulate`, or pbsim3 and the renaming), then the join as now. `long_read_sample` splits into
   its two halves. Memory: the draw job holds the names of all chunks; for Ultima only their count is needed
   (`hifi_reads` names reads after the templates), for Nanopore 750 k names per 6 Gb sample.
2. **Python jobs in worker processes.** Keep the `Scheduler` (threads that only wait), and have the jobs that run
   Python (`long_read_sample` or its two halves, `host_pe_chunk`) submit to a `ProcessPoolExecutor` of `--jobs`
   workers and wait for the result. The tasks are plain dicts (picklable); `scenarios.Host.of` opens its memory map
   per process. Use the `spawn` or `forkserver` start method, since the collector has threads when the pool starts.
3. **The genomes on scratch.** `build_gtdb_database.py` could copy the genome files of the simulated table
   (`genomes_simulated.tsv`, ~19 k genomes, ~22 GB gzipped) to `SCRATCH/genomes` once and give the collectors a
   table pointing there. After fix 1 this is ~260 GB of reads from the node's SSD instead of `/hpc-home`. Worth it
   if the check above shows `D` states, or once fix 2 lets 64 workers read at once.
4. **Priorities that know the community.** `long_read_seconds` estimates 10 s + bases x rate: a soil Ultima chunk
   (4,400 genome reads, ~150 s) comes out at 52 s, a gut Nanopore chunk too high next to it. With fix 1 the draw
   job's estimate is genomes x 21 ms x ~1.8 + reads x per-read cost.

Not worth doing: the `translate` parse alone (6.8 → 5.0 ms per genome, still under the GIL), `hifi_reads.py`
(scales on threads), the host fragments (0.02 h), a decompressed-genome cache (fix 1 leaves 1.8 reads per genome and
sample).

**After that, the next step on the critical path is protal on the scenarios' samples.** They add ~290 Gb of reads
to profile (per hold-in/hold-out sample set: gut 18 Gb, soil 24, soil_shallow 4.5, host 12; 5 sets) to v11's ~190
Gb, including 25 single samples of 6 Gb (gut and soil PacBio and Nanopore, soil Ultima); deep samples were the
tail of v11's protal run (13-24 min). v12's `profile_all/protal.log` on scratch ("Aligning reads took" per sample)
will show how much.

The website is not affected: nothing here changes how protal is run.

## Follow-up (2026-10-05): v12 ran out of scratch; fixes 1 and 2 and block profiling implemented

v12 (job 23962682, `--scratch` with 263.1 GB free, a genome pool of 24,959 species and 53,178 genomes with the
in-silico strains, so soil at its full 9,000-11,000 species) stopped 2:42 into the test set's simulations with
`No space left on device` (soil PacBio chunks, a gut PacBio join). The user asked whether the scenarios are deeper
than real studies, then for fixes 1 and 2, then to profile blocks of samples as they are done and remove their reads.

### Disk per read type (`scripts/disk_costs.py`, `scripts/zstd_levels.sh`, `results/disk_costs.txt`, `results/zstd_levels.txt`)

Bytes on disk as the collector writes them (synthetic genomes: random sequence, so the sequence half compresses a
little worse than real genomes):

| Reads | As written | Per sample of the scenarios |
|---|---|---|
| paired-end 100 / 150 / 250 bp (`simulate_metagenomes`, BGZF level 6) | 192 / 174 / 379 bytes per pair | 20M pairs of 150 bp: 3.5 GB |
| Ultima 300 bp (`hifi_reads.py`, gzip level 1) | 1.11 GB per Gb | 6 Gb: 6.7 GB |
| PacBio HiFi (`hifi_reads.py`, gzip level 1) | 1.05 GB per Gb | 6 Gb: 6.3 GB |
| Nanopore (pbsim3 renamed, gzip level 1) | 0.89 GB per Gb, and while a chunk runs pbsim3's files 1.4 and the plain templates 1.0 GB per Gb | 6 Gb: 5.3 GB |

One sample of each preset is ~52 GB (gut 15.1, soil 21.8, soil_shallow 3.8, host 11.2), 80% of it long and Ultima
reads; the default 3 hold-in + 2 hold-out samples ~260 GB. With the design (training ~99 GB, v12's test set ~27
GB), the training database (23.7 GB) and two host copies (6.2 GB), v12 needed ~420-450 GB at once, against 263 GB.

**Are the scenarios deeper than real studies? Mostly not.** 20M read pairs (6 Gb) is at the deep end of typical
human gut shotgun studies (often 2-10 Gb per sample) but common for current NovaSeq runs; for soil, 6 Gb of
Illumina or Ultima reads is shallow (soil studies commonly sequence tens of Gb per sample); 6 Gb of HiFi or Nanopore
reads is typical to low for long-read metagenomes (often 10-30 Gb per sample, multiplexed). These ranges are general
knowledge, not measured here. The disk problem is the number of deep samples on the disk at once, not their depth:
making them shallower would move the scenarios' scores away from the studies they stand for, while not keeping
their reads (below) costs nothing. If disk is still short: `--scenario-samples 2` saves ~52 GB.

**zstd** (one thread; ratio to the plain FASTQ; speeds in MB/s of plain FASTQ):

| Reads | As written | zstd -1 | zstd -3 | zstd -6 | zstd -12 | zstd -19 |
|---|---|---|---|---|---|---|
| paired-end 150 bp R1 | 0.250 (BGZF 6) | 0.259, 384/1089 | 0.247, 340/1116 | 0.246, 96/1240 | 0.220, 26/1338 | 0.185, 1.2/1110 |
| Ultima | 0.542 (gzip 1, 70/167) | 0.480, 305/874 | 0.473, 120/899 | 0.487 | 0.487 | 0.448, 0.7/222 |
| HiFi | 0.525 (gzip 1, 45/105) | 0.468, 206/696 | 0.462, 100/396 | 0.472 | 0.478 | 0.449, 0.6/207 |
| Nanopore | 0.446 (gzip 1, 46/63) | 0.451, 100/574 | 0.430, 67/370 | 0.410, 24/442 | 0.401, 7/364 | 0.367, 0.8/489 |

At the fast levels zstd saves ~0% (paired-end), 11-13% (Ultima, HiFi) and 4% (Nanopore): ~7% of a scenario set.
It compresses 2-7x and decompresses 4-9x faster than the gzip levels used now, which would save the simulators
~1.5-2 of the ~21 core-hours of drawn reads and some of protal's input decoding. It needs protal to read `.fq.zst`
(its reader takes plain, gzip and BGZF), `simulate_metagenomes` to write zstd, and the Python writers
`compression.zstd` (Python 3.14) or the `zstd` command. Worth doing for speed and for users' `.fq.zst` files, not as
the disk fix; not done here.

### What was implemented (written on `22c445d`, committed with this report on top of `186c8ed`)

`scripts/collect_training_data.py`, `scripts/build_gtdb_database.py`, `scripts/hifi_reads.py` (a docstring),
`scripts/mini_db/test_mini_db.py`, `docs/databases.md`:

1. **One pass over a sample's genomes for all its chunks** (`draw_templates`, `draw_chunks`, `make_reads`). A
   chunked sample is a draw job (one core, group `draw`, at most a quarter of `--jobs` at once), which writes every
   chunk's templates gzipped, then a job per chunk for its reads (priority above any new sample, so that drawn
   templates do not pile up), then the join, which removes each chunk's file once appended. pbsim3's reads are
   renamed from the templates' headers, so no job holds a sample's read names. Samples of one chunk are one job,
   as before. The reads are the same: **byte-identical** to `22c445d`'s for 21 samples of Ultima, HiFi and
   Nanopore, in one piece and in 4 chunks, one with a host share, on threads and in processes
   (`scripts/ab_collector.py`; wall on 6 slots: 102.9 s before, 80.7 s one pass on threads, 54.9 s one pass in
   processes, on 150-genome communities: the gain grows with the community and the number of chunks).
2. **The Python work in worker processes** (`Workers`: a `ProcessPoolExecutor` with the `forkserver` start method,
   started by the collector's main; `Workers.call` runs inline where none is started, as in the tests). The
   drawing, the chunks' reads, single long-read samples and the host's paired-end fragments run there; the
   Scheduler's threads wait for them.
3. **Profiling as the simulations go on** (`--follow`, the build's `--profile-blocks`, default 20 GB). The
   `--simulate_only` run writes `OUT/simulating.json` while it runs, marked `prepared` once its scenario tables and
   host genome are written (a follower waits for that, so that the two never write them at once). Once the training
   database is built, the build starts a follower per collection (`collect_training_data.py --follow`); it profiles the design points whose
   simulations have ended (their `simulated.json`) in a protal run once their reads reach `--profile_block` GB (all
   that are ready once the simulations have ended), the two followers' runs taking turns on a lock file. After each
   run the reads of every point whose read types are all profiled are removed (`sim/reads_removed.txt` lists them);
   SAMs, profiles and dumps stay (the parity check and `trace_relatives.py` read the SAMs). Then the follower writes
   its tables as before. Points it cannot profile once the simulations have ended (reads removed by an earlier run
   that profiled against another database, say) it simulates itself; a run without `--follow` does the same before
   profiling. protal's start-up at r226 is ~14 s warm (2026-10-04 report), so a run per block costs little.
4. **Room on the disk** (`--min_free`, the build's `--keep-free`, default 30 GB). A simulation job that writes
   reads starts only while the disk has room for it (an estimate: `PE_BYTES` etc.) and for the jobs running, and
   one that starts new work (a point, a sample's drawing) only if `--keep-free` GB are left besides; finishing work
   (a chunk's reads, the host's chunks) may use them. While nothing can start, the simulations wait for the
   followers to remove reads.

**Tests** (WSL, Python 3.12; build tests with the `0e800b2` binaries of `~/protal-scen` and scikit-learn from
`protal-db-build`): new `test_one_pass_and_workers` (chunk templates drawn together equal each drawn alone, a host
genome too; fewer genome reads; the same reads inline and in workers), `test_room_on_the_disk`, `test_follow` (a
run while the simulations go on, the reads removed, re-simulation after another protal), `test_scheduler` (group
limits), `GtdbBuildTest.test_g_profiled_as_simulated` (test_f's scenario build followed in blocks of a few kB:
**the same tables** as test_f's one protal run, every read removed, SAMs kept, a rerun profiles nothing);
`test_a` and `test_f` run with `--profile-blocks 0` (they check the one-run flow and read the reads afterwards).
Fixes 1 and 2 alone passed `GtdbBuildTest` (6 tests) and the 59 tests of `test_mini_db.py`.

## Reproducing

In WSL (`~/colprof`; pbsim3 and its models from `~/micromamba/envs/protal-db-build`):

```
bash scripts/run_all.sh ~/colprof /mnt/c/Users/hildebra/Documents/locDev/protal 22c445d
```

It makes the 400 genomes (`make_genomes.py`, ~440 MB), the unit costs (`unit_costs.py`), the identity check and the
scaling runs (`gil_scaling.py`), `hifi_scaling.py`, and the projection (`scenario_load.py`, ~4 min; its unit costs
are constants at the top of the script, from `unit_costs.py`). The follow-up's measurements:

```
python3 scripts/disk_costs.py --scripts ~/colprof/scripts --genomes ~/colprof/genomes/genomes.tsv \
    --simulator SIMULATE_METAGENOMES --out ~/colprof/disk --pbsim_models ~/micromamba/envs/protal-db-build/data
bash scripts/zstd_levels.sh ~/colprof/disk        # after copying a pe150_R1.fq.gz there (see the script)
python3 scripts/ab_collector.py --scripts OLD_SCRIPTS --genomes G --out ab/old --pbsim_models M
python3 scripts/ab_collector.py --scripts NEW_SCRIPTS --genomes G --out ab/new --pbsim_models M --processes
python3 scripts/ab_collector.py --compare ab/old ab/new
```
