# The GTDB build's order of steps after faster builds, simulators and streaming

Date: 2026-10-07. Branch `audit-fixes` at `54b1aa1` (`scripts/build_gtdb_database.py`, `collect_training_data.py`,
`machine_learning_cmdline.py`; protal's profiling stage in `src/RunProtal.h`). Read-only apart from a one-second
planning model ([`scripts/plan_model.py`](scripts/plan_model.py)): no build, protal or simulator was run (the local
cores are reserved).

**Since implemented** (same day, [section 8](#8-implemented-changes-1-3-and-the-scratch-estimate)): changes 1-3 and
the scratch estimate, which chooses `--stream-above` (`auto`, the new default). Change 4 is open.

**Question.** Several changes landed today: faster `protal --build` (S1-S4 of the
[database build audit](../2026-10-07-database-build-audit.md), `54b1aa1`); C++ Illumina and long-read simulators
([Illumina](../2026-10-07-illumina-model/README.md), [long reads](../2026-10-07-long-read-simulator/README.md)); the
~50 GB genome store on scratch ([simulator audit](../2026-10-07-simulate-metagenomes-audit/README.md)); samples above
`--stream-above` (2 GB) streamed into protal through named pipes. Scratch is also ~500 GB larger. Is the order in which
`build_gtdb_database.py` runs its programs still right? Where do programs compete for cores? This matters most for
gradient boosting, which slows down badly when its OpenMP threads share cores.

**Sources.**
- **r226 v15** (`local/v15`, SLURM job 24027266, `-t 84`, scripts of `ce85bd7`): `training_data*.log`,
  `test_data*.log`, the trainers' reports and the file times in `model_logs/`. Its `console.log`, index logs and
  per-run `protal.log` files are not in `local/`. Times below are derived from the collectors' clocks.
- **r226 v14** (`local/v14`, job 24005369, `-t 64`): the same logs plus `training_db_index.log`.
- **The default design at `54b1aa1`**, from the collector's own planning functions (`units_of`,
  `streamed_simulations`, `pe_bytes`, `sample_bytes`, `pe_point_seconds`, `long_read_seconds`):
  `plan_model.py`, and `plan_model.py --v15` for v15's design (4 scenarios of 6 + 3 samples, depth factors 1/2-2).

## Summary

**The order is no longer right.** It was built when the simulations were the critical path. In v15 they were:
the training samples took 4:04 to simulate, one soil point on one thread. After today's changes they should take
minutes to tens of minutes (~14 core-hours by the collector's own laptop constants, which are high for paired-end reads
now that genomes come from the store). The critical path is now:

> training database build → **protal profiling** (~670 GB of reads, ~2.6 h of protal runs at v15's rate) → ~2 min of
> training → parity → `--add_model`

Four parts of the current order lose wall time on that path. The first three are the largest:

| # | What happens now | Cost on the critical path (estimate) | Change |
|---|---|---|---|
| 1 | At the default `--stream-above 2`, 758 of 825 GB (the model's sizes) are streamed in **44 protal runs of their own** (25 training, 19 test). Each run pays protal's start-up and ends in its own profiling stage. In v14/v15, runs with deep long-read samples spent 2-5 min after the last sample began aligning | ~45-110 min more than batched runs (profiling stages ~30-90, start-ups ~15-20) | Choose streaming from the measured free scratch: with ~1 TB nothing needs streaming. If some must be, merge streamed points into shared runs and try `--profile_ahead` |
| 2 | The **training database** (the gate of profiling) is built beside the finished database's build (both `-t T`, nice 0) and both collections' simulations. v14: 33:28 against the finished database's 26:21. v15: profiling began ~61 min after both builds started | ~5-15 min with the faster builds (more in v15) | Build the training database alone. Build the finished one after it at idle CPU priority, during profiling |
| 3 | After the models, the database waits for `trace_relatives` and `error_reads` (v15: 10:55:36 → 11:10:38) before `--add_model` | ~15 min to "database ready" | `--add_model` right after parity. The reports run after it, in parallel |
| 4 | Setup is serial: genome table, conversion, gene neighbours, **in-silico strains**, then the derive (v15: 4.5 min). The in-silico strains feed only the simulations | the in-silico step's minutes plus ~3 min (P4/P5 of the build audit) | In-silico strains beside the derive and training build; conversion beside the genome table |

**Gradient boosting.** In v15 each trainer took 68-120 s. The four ran side by side on 84 cores, 21 threads each (5
fold jobs × 4 threads), with nothing else on the node. Nothing competes with them today. The point is to keep it so:
with change 2 the finished database's build may still be running when training starts, and it must then pause.
Pipelined steps (parity, error reads) may use only a finished trainer's threads.

**Scratch.** The script prints the free space at the start and the peak at the end, but its choices ignore the free
space. It should measure the free space at the start, estimate what the plan needs and pick the streaming from that.
At the end it should report the peak against that estimate (section 5). With the default streaming, the written reads
are only ~54 GB (model 66 GB × 0.81, the calibration in section 2). The extra 500 GB is unused unless streaming is
reduced.

Together (estimates, to be checked on the next build): from the training database's start to "database ready", ~4-5.5 h
with today's order and ~3-3.3 h with changes 1-3. Most of the difference, and most of the uncertainty, is change 1.

## 1. The order today

`main()` in `build_gtdb_database.py` at `54b1aa1`, with the defaults (`--profile-blocks 20`, `--stream-above 2`, no
`--one-build-at-a-time`):

| Step | Needs | Threads, priority | v15 time |
|---|---|---|---|
| genome table (`make_genome_table`) | the release | `-t` | ~3 min (v11, build audit P4) |
| conversion + `gene_neighbours.py` (into `OUTDIR/protal_db`) | the release; neighbours also the table | `-t` | 3.7 min + ~1-2 min (v14 `convert.log`) |
| in-silico strains → `genomes_simulated.tsv` | `gene_positions.tsv`, taxonomy | `-t` | not in `local/` |
| held-out species; derive `SCRATCH/training_db` | taxonomy, table | `-t` | 06:36:51 → 06:41:19 (4.5 min) |
| **both collections' simulations** (`--simulate_only`) | `genomes_simulated.tsv`, held-out species | 2 × `-t` slots, nice 10 | training 4:04:22, test 2:52:47 |
| **finished database build** (background) | `protal_db` | `-t`, nice 0 | not in `local/` (v14 26:21) |
| **training database build** | its files | `-t`, nice 0 | ~61 min (to the followers' first runs, ~07:42) |
| followers: protal runs in blocks of ≥ 20 GB, streamed points one run each, one run at a time (lock) | training database, simulated reads | `-t` each, nice 0; streamed runs' simulator `-t`, nice 0 | training 5 runs, 93:53 of protal; test 4 runs, 45:46 (v15 streamed nothing: streaming came after it) |
| 4 trainers side by side | both collections' tables | `-t`/4 each | 10:53:11 → ~10:55 (68-120 s) |
| parity, 4 protal runs one after another | models, training database | `-t` | done by 10:55:36 |
| `trace_relatives` (pe), `error_reads` (per read type, one after another) | SAMs | `-t` | → 10:58:08; → ~11:10:38 |
| wait for the finished build; `--add_model` | finished database | `-t` | seconds |

v15's timeline: the derive started at 06:36:51. The simulations and both builds started at 06:41:19 and profiling
at ~07:42. The test simulations ended at 09:34 and the training ones at 10:45:41. The last protal run ended ~10:52,
the trainers started at 10:53:11 and the summary was written at 11:10:40. protal ran 140 min of the ~190 min the
followers were active. The rest was waiting for the simulations: from 2:12:48 to 3:02:50 on the training follower's
clock, it logged "0 of them simulated" every 10 min.

## 2. What the changes do to the critical path

| Stage | v15 | next build, expected | Basis |
|---|---|---|---|
| simulations, both collections | 4:04 (training), 2:53 (test); ART, Python long reads | ~14 core-hours by the collector's constants (pe 4.8 + 1.8 with pre-store genome costs, drawn 5.4 + 2.2): tens of minutes, mostly inside streamed runs | `plan_model.py`; per genome < 0.4 ms from the store, per pair ~8.8-18 µs (simulator audit §7.3) |
| reads protal profiles | 405.5 + 186.0 GB (logged) | 589 + 236 GB (model) × 0.81 ≈ **477 + 191 GB** | the model gives 503 + 231 GB for v15's design against the logged 405.5 + 186.0 (ratio 0.81 for both). The next build has 10 + 4 scenario samples and a fifth scenario, depths 1/8-2× |
| protal time (reads ÷ v15's rate, 4.3 and 4.1 GB/min) | 94 + 46 min | **~110 + 47 min**, plus the cost of more runs (§3.1) | v15's protal runs at 84 threads |
| training database build | ~61 min (shared with the finished build and the simulations, which filled ~350 GB of scratch in their first 40 min) | S1-S4: the build audit estimated ~10-12 min for a build; shared, more | `54b1aa1`; v14 phase table in the build audit |
| trainers | 68-120 s | similar (+ moderate scenario rows) | v15 reports' `time:` lines |
| post-training steps | ~15 min (error reads ~12) | similar | `model_logs/` file times |

**Profiling is now the long pole.** Its length is set by the reads and by protal's speed, which the build script
cannot change. What the script controls is when profiling starts (the training database) and how many runs it is cut
into.

## 3. Findings and changes, in order of gain

### 3.1 Streaming at 2 GB turns most of the reads into 44 protal runs

`plan_model.py` at the defaults:

| `--stream-above` | streamed runs (training + test) | of them long-read or Ultima | streamed (model GB) | written to scratch (model GB; × 0.81 real) |
|---|---|---|---|---|
| 2 (default) | 25 + 19 = **44** | 15 + 15 | 543 + 215 | 46 + 20 (~54) |
| 4 | 19 + 15 = 34 | 13 + 11 | 483 + 193 | 106 + 43 (~121) |
| 8 | 7 + 6 = **13** | 5 + 5 | 266 + 109 | 322 + 127 (~364) |
| 16 or none | **0** | 0 | 0 | 589 + 236 (~668) |

Every streamed simulation is a protal run of its own (`stream_run`): its samples, one at a time, then their
profiling stage. Batched runs share one profiling stage. Two costs grow with the number of runs:

- **The profiling stage at the end of a run.** protal profiles every sample after aligning all of them
  (`ProfileWrapper`). Samples run in parallel, each with threads in proportion to its SAM, so a few deep samples set
  how long it takes. From the last "aligning sample" mark to the run's end:

  | Run | Samples | Last aligning mark → end |
  |---|---|---|
  | v14 training 8 | 6 pb + 6 ont (2 scenario points), 68.5 GB | 21:00 → 26:10 |
  | v14 test 3 | 3 pb + 3 ont, 36.9 GB | 8:00 → 13:14 |
  | v15 test 2 | 3 pe + 3 pb + 3 ont, 47.8 GB | 12:00 → 17:19 |
  | v15 training 3 | 8 pb + 120 ont, 95.1 GB | 24:00 → 28:01 |
  | v15 training 4 | 6 each of pe, se, pb, ont, 137.3 GB | 22:00 → 24:23 |
  | v14 training 6 | 6 ont, 34.8 GB | 20:00 → 22:37 |
  | v14 training 3 / 2 | 6 pe / 5 pe + 5 se | 3:00 → 3:06 / 3:00 → 3:59 |

  The last sample's alignment (~1-2 min for a 5-7 GB sample) is inside these spans. The profiling stage itself is
  therefore ~1-4 min in runs with deep long-read samples and under a minute for paired-end runs. Streaming gives ~30
  long-read and Ultima runs. If each pays 1-3 min instead of a handful of batched runs sharing it, that is **~30-90 min
  more**.
- **Start-up.** protal's start-up at r226 is ~15-30 s (index load 3.2 s, preload 3.2 s, the rest per the
  [GTDB-scale report](../2026-10-04-performance-gtdb-scale/README.md)). The simulator also loads the genome table and
  its design. That is ~44 × 0.3-0.5 min ≈ **15-20 min**.

Streaming is no longer needed to save time: protal, not the simulator, bounds a streamed run. It saves only space.
That space is now there: a scratch of ~1 TB holds all ~668 GB at once, besides the store, the training database and
the SAMs. Changes:

1. **Choose the streaming from the measured free space** (section 5): stream the largest simulations only until the
   rest fits. With ~1 TB this streams nothing. With v15's 500 GB it is about `--stream-above 8` (13 runs, ~364 GB
   written over the run, removed as profiled). The batched runs then pick up whatever is ready. The simulations take
   ~15 times less CPU than protal (~14 against ~220 core-hours), so the first block is ready within minutes and later blocks grow by themselves.
2. If many points must still be streamed (small scratch): **one protal run for several streamed simulations**, their
   samples in map order and a simulator companion per point (`stream_run` takes a list). The profiling stages
   then overlap, as in batched runs.
3. **`--profile_ahead` for runs with few deep samples.** Each sample is profiled while the next one aligns, on a
   quarter of the threads. It was written for a node where the profiling stage leaves cores idle, and left opt-in
   because a 6-vCPU laptop showed no gain ([performance review, item 4](../2026-10-03-performance-review/README.md)).
   It adds 25% threads beside the alignment, so A/B one streamed scenario point on the cluster before making it the
   collector's default.
4. Optionally, a **cap on a block** (e.g. 150 GB) besides its 20 GB minimum. With fast simulations the first run
   would otherwise take everything ready (hundreds of GB). A failure there re-profiles all of it, and its reads stay
   on the disk until it ends.

### 3.2 The training database's build shares the node with the finished database's

`build_gtdb_database.py:2221-2240` starts the finished build in the background (`-t T`, nice 0) and then builds the
training database (`-t T`, nice 0) in the foreground. Both collections' simulations run at nice 10 with `T` slots each.
In v14 the two builds' CPU-bound phases ran at the same time and took the same time: suspect copies 9:24 / 9:45,
uniqueness 3:06 / 3:29. Each had half the node. The training database's write phases were 2-4× slower than the
finished one's (12:22 / 3:31 writing the index): it writes to the scratch that the simulations fill. In v15,
profiling began ~61 min after the builds started. That was the simulations' heaviest hour: ART, and the scratch went from
499.9 to 149.7 GB free by the first "waits for room" at 0:40.

The finished database is read only by `--add_model`, ~3-4 hours later. So:

1. **The training database is built alone at normal priority.** The finished database's build starts after it, under
   `SCHED_IDLE` (`os.sched_setscheduler(0, os.SCHED_IDLE, ...)` in `Job`'s `preexec_fn`; nice 19 as a fallback). It
   writes to `OUTDIR` and adds no scratch I/O. It then takes only the cycles protal leaves idle: start-ups, the
   profiling stages, the waits between runs. At r226 it needs ~5-10 core-hours after S1-S4 (estimate). Profiling leaves
   tens of core-hours idle at 84 cores over ~2.5 h (estimate). Memory: ~35 GB beside protal's ~35 GB.
2. Its inputs must be read before it packs `protal_db`: the derive reads it, and so does the in-silico step if that
   moves (3.4). So it starts once both have read it, as now.
3. If it has not ended when training starts, **pause it** (`os.killpg(pid, SIGSTOP)`, `SIGCONT` after the trainers;
   ~2 min). Under `SCHED_IDLE` the trainers' threads would preempt it anyway. But memory bandwidth and cache are not
   scheduled, and boosting's OpenMP threads wait on the slowest at every barrier.
4. **Simulations during the training build: lowest I/O priority too** (`ionice -c 3` in front of the collector's
   command). Node-local NVMe often runs the `none` I/O scheduler, which ignores I/O priorities
   (`cat /sys/block/<dev>/queue/scheduler`). If so, start the simulations after the training database is built: the
   first 20 GB block takes about a minute at full speed. The genome store's first fill (~54k FASTAs from the network
   file system) adds a few minutes. Compare `training_db_index.log`'s write phases with v14's (12:22 index write)
   before choosing.

`--one-build-at-a-time` stays as it is (both builds in sequence). It costs the finished build's whole time at the end.

### 3.3 The database waits for the error-read reports

After parity (done by 10:55:36 in v15), `trace_relatives` (→ 10:58:08) and `error_reads` (→ ~11:10:38) ran before the
script waited for the finished build and ran `--add_model`. Neither report feeds the database. Instead:

- **parity → `--add_model` → "Ready protal database"**, then the reports. The database is ready ~15 min sooner.
- `error_reads` runs per read type, one after another, each on `-t` threads (`build_gtdb_database.py:740-775`). The
  four read types are independent. Run them side by side on `-t`/4 threads each, with `trace_relatives` beside them,
  and measure the gain: they read SAMs on scratch, and v15's pe extraction was killed for memory with 84 workers (it has a memory
  budget since `2543408`).

### 3.4 Setup: the in-silico strains are on the critical path but feed only the simulations

`insilico_strains.py` writes `genomes_simulated.tsv`, which only the collections read. The held-out species, the
derive and the training database use `genomes.tsv`, and the scenario notes are text. Running it beside the derive and
the training build, at nice 19, takes its minutes off the path. The simulations then start once it ends. Two smaller
items from the build audit apply as before:
- P4: one `os.walk` instead of six recursive globs in the genome-table step, ~1-2 min. The converter does not need the
  genome table, so it can run beside that step.
- P5: `gene_neighbours.py` while the converter writes the full reference, ~1:13.

The audit's P1 (start the simulations before the derive) gains nothing now that they are off the critical path.

## 4. Cores: where programs compete

| Phase | Running together today | Threads on `T` cores | Effect | After the changes |
|---|---|---|---|---|
| builds | training build (`T`, nice 0) + finished build (`T`, nice 0) + 2 simulation collectors (`T` slots each, nice 10) | ~4T | the critical-path build gets about half the node; its scratch writes compete with the simulations' | training build alone; the others idle-class |
| profiling | protal (`T`) + a streamed run's simulator (`T`, nice 0, a few % of protal's CPU) + background simulations (nice 10) + the finished build if still running (`T`, nice 0) | 2-3T | protal shares with an equal-priority build if the build is slow | finished build `SCHED_IDLE`. Fewer streamed runs, so fewer simulator companions |
| training | 4 trainers × (5 fold jobs × 4 threads) = 80 of 84 | ≤ T | none observed (68-120 s each) | keep alone: pause the finished build; nothing else at normal priority |
| reports | `trace_relatives`, `error_reads` one after another, `-t` each | T | none | side by side, threads split, after `--add_model` |

Rules for the build:
- Programs on the critical path (training build, protal, trainers) run at normal priority with `T` threads. Nothing
  else runs at normal priority beside them. Background work is idle-class (CPU and, where the device honours it, I/O).
- Never put OpenMP programs side by side at normal priority with more threads in total than cores. The trainers'
  split (`T`/4 each, fold jobs × fold threads ≤ `T`/4) already respects that. The
  [training-time report](../2026-10-06-training-time/README.md) saw one boosted fit 4× slower on 6 threads than on 3
  under load. A step pipelined beside trainers (parity of a finished read type, say) takes only that trainer's
  threads.
- **Check `-t` against the allocation**: `len(os.sched_getaffinity(0))` is what SLURM gave the job. A `-t` above it
  oversubscribes every stage at once. A one-line warning at the start (the trainer's report prints
  `os.cpu_count()`, the node's CPUs).

## 5. Scratch: what the script should measure

Today: `Scratch` (`build_gtdb_database.py:249-264`) records the file system's use at the start and its peak while the
script waits, and prints the free space at the start and the peak at the end. Nothing is decided from it. The
collector's `--min_free` (`--keep-free 30`) only makes a simulation wait.

Proposed:
1. **At the start, estimate the plan's need** and print it beside the free space:
   - the genome store: 0.25 bytes per base of the genomes in `genomes_simulated.tsv` not yet stored (the table has
     the lengths), ~50 GB at r226;
   - the training database: ~22 GB, plus what its build holds before packing (L7 of the build audit; S2 made it
     smaller). Measure it once;
   - the reads written: the collector's per-simulation sizes (as `plan_model.py` computes them) for the streaming
     chosen. Its constants were 1/0.81 of v15's real sizes, so they are safe as they are;
   - what stays: SAMs, dumps and the host genome. at the end of v15's profiling ~62 GB more were in use than when its collectors
     started (499.9 → 437.5 GB free), the training database's growth included;
   - `--keep-free`.
2. **`--stream-above auto` (the new default)**: no streaming if everything fits at once; else stream the largest
   simulations, one by one, until the written rest fits. A number keeps today's behaviour. Print the choice: "scratch
   N GB free; the plan needs ~M GB; K simulations streamed (above X GB)".
3. **Fail early** when even full streaming does not fit (store + training database + SAMs + keep-free), instead of
   hours later (v12's ENOSPC).
4. **At the end, the peak against the estimate, and its parts**: store, training database, reads now, SAMs. The
   next build's estimate can then be calibrated.
5. Copy each protal run's `protal.log` (small) from `SCRATCH/*/profile_all/run*/` to `OUTDIR`. Node-local scratch
   goes with the job, and these logs hold the stage timers that would confirm §3.1.

The peak includes other jobs' files if several share the node's disk. A `du` of the run's own folders at the end
tells them apart.

## 6. The proposed order

1. Genome table ∥ conversion → gene neighbours (P5 beside the full reference's write).
2. Held-out species → derive → **training database build** (alone, `T` threads).
   ∥ in-silico strains (nice 19) → simulations of both collections (idle CPU and I/O priority; or after step 2's build
   if the scratch ignores I/O priority). Streaming chosen from the measured scratch.
3. After the training build: **profiling** (followers as now; fewer streamed runs; `--profile_ahead` if the A/B
   gains). ∥ the **finished database build**, `SCHED_IDLE`.
4. **Trainers** alone (the finished build paused if still running).
5. Parity → (finished build done) → **`--add_model`**: database ready.
6. Reports: `trace_relatives` ∥ `error_reads` per read type ∥ the share archive.

## 7. What the next r226 build should show

- `console.log`: each step's start and end (written since `d5ffc20`).
- `*_simulation.log`: the last "in all" time. Expect well under an hour, even with today's order.
- `training_data.log`, `test_data.log`: the number of "protal run N" lines and how many say "streamed". Time from the
  last "aligning sample" mark to "protal profiled". The sum of the runs' times against ~157 min.
- `training_db_index.log`: "Run build took" and the write phases, against v14 (33:28; 12:22 writing the index).
- If kept: `profile_all/run*/protal.log` "Processing all samples took" and "Profiling took" per run.
- `classifier_training*.log` or the reports' `time:` lines: still ~1-2 min per trainer.
- The scratch line at the end.

The website is not affected: nothing here changes how protal is run. If `--stream-above auto` or the new order is
implemented, `docs/databases.md` (build stages, scratch size) needs updating.

## 8. Implemented: changes 1-3 and the scratch estimate

The user asked to implement changes 1-3 and the scratch estimate, with change 1 choosing `--stream-above` from that
estimate. Commit `73607bb` (`build_gtdb_database.py`, `collect_training_data.py`, `trace_relatives.py`, the tests,
`docs/databases.md`). Change 4 (setup in parallel) is not part of it.

**The scratch estimate and the streaming (change 1).**
- `disk_needs()` counts what the run puts on the samples' disk besides the reads it writes. Each part is counted only if
  it is on that file system:
  - the genome store's growth: 0.25 bytes a base of `genomes_simulated.tsv`'s genomes, less what the store holds;
  - each database still to be built there: 1.5 times its `reference.fna`;
  - the host genome of each collection that has not prepared it;
  - the SAMs and profiles: a tenth of the collector's read estimate;
  - `--keep-free`.
- `simulation_sizes()` uses the collector's `simulation_bytes()` to give each simulation's largest sample (what
  `--stream_above` compares) and all its reads.
- **Before the builds**, the run says the free space, these needs and the reads' total. With `--stream-above auto` it
  stops if the needs alone exceed the free space (nothing fits even when every point is streamed).
- **Once the training database is built**, it measures again. `stream_threshold()` picks the value that streams the
  fewest simulations, largest samples first, until the others' reads fit at once. The value lies halfway between two
  simulations' largest samples, so equal ones go together, and is passed to both collectors exactly (`repr`).
- The choice is printed (`--stream-above auto: N of the M simulations streamed ...`) and recorded in
  `build_metadata.tsv` (`samples_streamed`). The end of the run gives the peak beside the estimate and what the
  store, the training database and the samples hold.
- At the defaults (`plan_model.py --free`, with a 50 GB store and `--keep-free 30`; 162 GB besides the reads):

  | free once the training database is built | `--stream-above` chosen | streamed | written |
  |---|---|---|---|
  | 1,000 GB or more | none | 0 of 116 | 824 GB (estimate) |
  | 900 GB | 12.1 GB | 2 | 734 GB |
  | 700 GB | 9.75 GB | 9 | 533 GB |
  | 500 GB | 5.7 GB | 21 | 306 GB |
  | 300 GB | 3.64 GB | 36 | 125 GB |

- **Streamed runs merged.** `stream_run()` takes several simulations: every streamed simulation whose communities are
  there joins one protal run, a simulator companion per point, in map order. Companions waiting for later pipes block
  at their first open and hold little; the pipeline opens a sample's outputs before it loads genomes.
- **`--profile-block-max`** (200 GB; collector `--profile_block_max`, default 0) caps a run, at least one point. It
  caps written blocks too, by whole simulations, first fit (`first_batch()`).
- **`--profile-ahead`** (collector `--profile_ahead`) passes protal's `--profile_ahead`. It is off by default until an
  A/B on a cluster node.
- Each protal run's log is copied to `OUTDIR/<collection>/protal_runs/` with `--scratch`.

**The order (change 2).**
- The training database is built alone: no finished build and no simulation beside it.
- After it, the finished database's build starts in the background under `SCHED_IDLE`
  (`Job(idle=True)`; nice 19 where that class is missing). Then both collections' simulations start at nice 10
  (before `--rank-genes`, as before).
- Before the trainers start, the finished build is paused (`Job.pause()`, SIGSTOP to its group) and resumed after
  them. `Job.kill()` sends SIGCONT after its SIGTERM, so a paused build still stops at once.
- The cost: the first protal block waits for the first simulations and the genome store's first fill. At r226 that is
  a few minutes (estimate), against v15's ~61 min training build beside the simulations.

**The reports after the database (change 3).**
- Parity, then `--add_model` and "Ready protal database". Then `reports()`: `trace_relatives.py` and one
  `error_reads.py` per read type, side by side.
- Each run gets an even share of `-t`, and each error-read run the matching share of `error_reads.py`'s memory budget
  (`--memory`).
- Logs: `trace_relatives.log` and `error_reads_<read type>.log`. The single `error_reads.log` is gone.
- The genomes' contig names are read once into `genome_contigs.tsv.gz` before the parallel runs. Every pe point's
  manifest of both collections is read, scenarios included.
- `trace_relatives.genome_contigs()` now locks the cache (`CACHE.lock`): shared to read, exclusive to append one whole
  gzip member. Runs side by side then neither mix their writes nor read a half-written member.

**Tests.**
- `test_collector`: `simulation_bytes` agrees with `streamed_simulations` at four values; `first_batch`.
- `test_gtdb_build`: `stream_threshold` and `stream_spec`; `disk_needs` part by part; a `Job` under `SCHED_IDLE`
  paused, resumed and killed while paused.
- `test_gtdb_pipeline`:
  - test_a: the order (training database, then the finished build at idle, then the simulations; "Ready" before the
    reports), the room line, the streaming line, the protal run log copied;
  - test_h: both streamed points of the training data in one protal run;
  - the rest unchanged.
- Results, from a snapshot of the working tree (`git stash create`, `git archive`), 4 pinned WSL cores,
  `PROTAL_TESTS_REQUIRED=1`. The binaries together equal `54b1aa1`: protal from `~/dbaudit/new`, the simulator from
  `~/simaudit/mine`.
  - `test_gtdb_pipeline`, `test_collector` and `test_gtdb_build`: 40 OK in 190 s.
  - `test_error_reads` and `test_trace_relatives`: 16 OK.
- The first pipeline run failed test_b. Its check that no training table was written stood for "stops at once", and
  it no longer holds: the finished build now starts after the training database, and at mini scale the collection can
  end before a background failure is noticed. test_b's build now fails at once, and the test checks that no model is
  trained. The script also looks at its background jobs every second instead of every 5 s.

**What the next r226 build should show** (besides section 7):
- the `room on ...` and `--stream-above auto: ...` lines, and their estimate against the peak at the end;
- `training_db_index.log` against v14's 33:28 (S1-S4 and no neighbours on the node);
- the streamed runs' "protal run N: K simulations streamed" lines and how long their profiling stage takes
  (`OUTDIR/<collection>/protal_runs/`);
- the line "... build in the background paused meanwhile" at training;
- "Ready protal database" before the reports.

## Reproducing

```bash
python3 docs/claude/2026-10-07-build-ordering/scripts/plan_model.py          # the default design at 54b1aa1
python3 docs/claude/2026-10-07-build-ordering/scripts/plan_model.py --v15    # v15's design, to calibrate the sizes
python3 docs/claude/2026-10-07-build-ordering/scripts/plan_model.py --free 300,500,700,900,1000   # --stream-above auto (section 8)
grep -h "protal run [0-9]*:\|protal profiled\|aligning sample" local/v15/training_data.log local/v15/test_data.log
grep -h "^time:" local/v15/model_logs/trained_model*.report.txt
```

`plan_model.py` imports `scripts/collect_training_data.py` and `scripts/scenarios.py`. It reads and writes no files
(`--v15` writes a temporary scenario file) and runs in about a second.
