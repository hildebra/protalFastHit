# Wall time on 32 threads: where the 33 s run went, and what each feature costs against its F1

**Question (the user's, 2026-10-12).** On 32 threads the r226 paired-end run takes 41.6-48.2 s against the
best of ~33 s; should features be cut, weighing what each adds to F1 against what it costs in wall time?

**Data.**
- The new run: `local/Perf0.8.2e` (the user's cluster job 24198970, `scripts/measure_performance.sh`, 32 threads,
  node q512n19, protal 0.8.2 `61f33f1`, database `protal0.8.1_r226_v24`, samples pe1 = PEARL_AGE B10 paired-end
  and pb2 = a HiFi barcode, 3 runs each, then a cohort run of both).
- The 33 s run: `local/Perf0.7.5_v13` (thirteenth cluster run, protal 0.7.7 `230bf64`, r226 v11, node q512n10,
  same samples and script).
- For the trend: `local/Perf0.8.1` (0.8.1 `7d5df4b`, 32 threads) and `local/Perf0.8.2d` (0.8.2 `ccdbf52`,
  8 threads).
- F1 of the feature groups: refits of the r226 v23 build with groups left out
  ([2026-10-11-r226-v23](2026-10-11-r226-v23/README.md), sections 4-5: species held out, knob 0.5, paired
  bootstrap), the strain study ([2026-10-11-strain-false-negatives](2026-10-11-strain-false-negatives/README.md))
  and [2026-10-10-r226-v22](2026-10-10-r226-v22/README.md).
- Code read at `61f33f1`. Nothing was rerun for this report: the numbers come from the runs' `runs.tsv`,
  `stages.tsv` and `logs/*.log`.

## Summary

- **Aligning is not where the time went.** pe1's best run (41.65 s) aligned in 23.6 s, 1.1 s *faster* than the
  33.09 s run. The 8.6 s difference is profiling (+6.0 s) and finishing the SAM on the network file system (+4.0 s).
- **The SAM parse in profiling stopped scaling with threads**: 12 thread-seconds at 0.7.7, 119 at 0.8.2 on 32
  threads, 32 on 8 threads. Its parallel part takes ~4-5.6 s of wall time whatever the thread count. The cause is
  very likely the run-wide caches of ancestry sites and column weights (added after 0.7.7). Every 1 MB chunk of the
  SAM starts with empty caches of its own and takes up to four global mutexes per gene copy it meets. System time
  rose from 24 to 95 s, which fits contended mutexes sleeping in the kernel. This costs ~3.4 s of wall time and is
  a scaling bug, not what the features cost.
- **1.8 s of profiling falls outside every timer**: the end of `ProfileSample` after its sample timer stops.
  New there since 0.7.7: `malloc_trim(0)` (`899d82f`, 2026-10-07), the release of the evidence and the read data,
  and the scan for unreported species.
- **SAM finish**: 0.74 s at 0.7.7 on q512n10, but 1.6-6.7 s on q512n19. One run's alignment also stalled 7 s
  writing the SAM. Both depend on the file server, not on protal's version.
- **What the features themselves cost is small at 32 threads.** An upper bound for everything in profiling that is
  new since 0.7.7: ~0.7 s of parse once it scales, +0.4 s evidence and +0.3 s profile writing, so ~1.4 s of 41.6.
  In the alignment, the allele scores cost nothing measurable and the adaptive candidates are 0.12% of the
  candidates tried.
- **What they are worth**: everything since the v9 core set is worth +0.0055 pe F1 (se +0.0051, pb +0.0042, ONT
  +0.0058). The three groups built on the ancestry sites are worth +0.0022 pe (ONT +0.0031). The strain alleles
  turn ~7% more strains into calls in the strain study.
- **Recommendation: cut no feature for wall time.** Fix the parse contention (−3 s expected) and time and trim the
  profiling tail (up to −1.8 s). For benchmarks, write the outputs to node-local disk. Together these should bring
  pe1 back to ~33-37 s with every feature kept (37 s if the SAM close on this file server stays as it was).
  Feature cuts would buy ~1 s at best. The candidates are column weights (1.47 GB resident, ~0 F1 alone) and
  polymorphic (~0 except ONT), and the reason to cut them would be memory and simplicity, not speed.

## 1. The 41.65 s run against the 33.09 s run

Best runs (run 3 of each; the 0.8.2 median, 48.18 s, is run 2, whose alignment stalled on the SAM output: "Output
handler" 10.97 s per thread against 0.68 in run 3).

| | 0.7.7, v11 DB, q512n10 | 0.8.2, v24 DB, q512n19 | difference |
|---|---|---|---|
| wall | 33.09 | 41.65 | **+8.56** |
| start-up before aligning (index 2.67 / 2.75) | 3.83 | 3.49 | −0.34 |
| aligning | 24.77 | 23.64 | −1.13 |
| writing the SAM header and file (`sam_finish_s`) | 0.74 | 4.76 | +4.01 |
| profiling | 3.37 | 9.32 | +5.96 |
| of it: reading the SAM | 1.85 | 5.30 | +3.45 |
| of it: record evidence and sample context | 1.00 | 1.37 | +0.37 |
| of it: writing the profile | 0.32 | 0.63 | +0.31 |
| of it: read EM, congener distances, SNPs, scoring | 0.96 | 0.96 | 0.00 |
| of it: outside the sample's timers | 0.10 | 1.91 | +1.81 |
| process exit, freeing | 0.38 | 0.44 | +0.06 |

The start-up row is "Run protal took" minus the stages below it. The profiling parts overlap a little, so they add
up to more than the "Profiling: sample" timer (3.27 and 7.42 s); the differences are what counts. The aligning
sub-timers moved little: seeding 12.9 → 11.6 s per thread, alignment handler 2.8 → 3.1, joining pairs and sorting
0.45 → 0.64.

PacBio (pb2, run 3 of each) lost 1.9 s (12.25 → 14.11 s): aligning 7.93 → 7.46, writing the SAM 0.30 → 1.61,
profiling 0.35 → 1.05 (reading the SAM 0.2 → 0.3 s, writing the profile 0.0 → 0.2 s, the untimed tail 0.25 s).

## 2. The SAM parse no longer scales

Each sample's log line `Sample X SAM read:` gives the parse time summed over the threads, and the wall time of the
threads' part:

| run | threads | pe1 parsing (thread-s) | pe1 parallel part (wall s) | pb2 parsing (thread-s) |
|---|---|---|---|---|
| 0.7.7 (v13) | 32 | 12.08 | 0.90 | 0.78 |
| 0.7.7 (v13), cohort | 28 / 4 | 9.32 | 0.85 | 0.77 |
| 0.8.1 | 32 | 164.75 | 5.69 | |
| 0.8.2 (`ccdbf52`) | 8 | 32.10 | 5.61 | 2.11 |
| 0.8.2 (`61f33f1`) | 32 | 118.71 | 4.18 | 5.29 |
| 0.8.2 (`61f33f1`), cohort | 28 / 4 | 108.35 | 4.42 | 4.11 |

At 0.7.7 the parse scaled perfectly (pb2 0.77 vs 0.78 thread-s on 4 and 32 threads). Since 0.8.1 each thread
spends ~3.7-4.0 s whatever their number, so four times the threads do four times the work in about the same wall
time. In the cohort, pb2's 4 threads parse slower (4.11 thread-s) than pb2 alone on 8 (2.11), because pe1's 28
threads are parsing at the same time. Something process-wide is serialising the parse. Over the same versions
system time went from 23.6-24.6 s to 92-97 s at 32 threads (16-18 s at 8 threads), and user time stayed the same
(954 s).

**The likely cause** (from the code; a `perf` profile of a `--profile_only` run would confirm it):
- `RecordEvidenceCollector::NoteRecord` (`Profiler.h:3156`) looks up each record's gene copy in its collector's
  own maps (`m_ancestry_sites`, `m_species_sites`, `m_column_weights`). On a miss it asks the `GenomeLoader`'s
  run-wide caches.
- There is one collector per SAM chunk (`ProfileSam`, `Profiler.h:5580`: `profile.Evidence().Empty()` per
  chunk), and the pe1 SAM has 1,617 chunks of ~1,375 records. So almost every copy a chunk meets is a miss, up to
  ~2 M misses for the sample.
- Each miss takes up to four acquisitions of three global mutexes:
  - `AncestrySitesOf` (`ancestry::Cache::Lookup`, `AncestrySites.h:644`);
  - `SpeciesAncestrySitesOf`, which calls `AncestrySitesOf` again and then takes `m_species_sites_mutex`
    (`GenomeLoader.h:1261`);
  - `ColumnWeightsOf` (`column_weights::Cache::Get`, `ColumnWeights.h:611`).
- With 32 threads on the same mutexes, each acquisition becomes a futex sleep and wake-up, which is system time.
  None of these lookups existed at 0.7.7: the column weights and the per-copy ancestry sites came in 0.8.0-0.8.1.
- The comment at `Profiler.h:3187` assumed the opposite ("a chunk's records hit the same taxa and genes over and
  over"). At r226 a 1 MB chunk spreads over hundreds to thousands of copies.

**Fix (outputs unchanged: a cached value depends on its key only).**
- Shard the three caches by key hash (e.g. 64 shards, each with its own mutex and 1/64 of `kMaxCached`).
- Let `SpeciesAncestrySitesOf` take the sites the caller already has instead of looking them up again.

Expected: the parse back to about 1 s of wall time on 32 threads. 0.7.7 took 0.9 s. On top of that comes the new
features' own work: at most (32 − 12)/32 ≈ 0.6 s, because the 8-thread run still has some contention. So reading
the SAM goes from 5.3 to ~2-2.5 s, which is −3 s of pe1's wall time. PacBio gains ~0.1 s (its parse is 0.2 s of
wall time even with the contention).

## 3. The rest of the gap

- **Profiling outside the sample's timers, 1.9 s** (0.1 s at 0.7.7; 1.95 s in the cohort too). After
  `bm_sample.Stop()` `ProfileSample` does the following (`RunProtal.h` at `61f33f1`):
  - it scores every taxon of the profile again for `unreported_species.tsv` (`filter.Score` before the `continue`);
  - it calls `ReleaseReadData` and `ReleaseEvidence` (2026-10-11);
  - it calls `malloc_trim(0)` (`RunProtal.h:1735`, since `899d82f`, for cohorts' memory).

  `malloc_trim` walks every glibc arena and returns its free pages one `madvise` at a time. After a 30 GB run that
  can take seconds of system time. Put a timer on this tail first. Then skip `malloc_trim` when no other sample
  of the process remains to be profiled. If the unreported-species scan is part of it, check `StrongOwnEvidence`
  before scoring.
- **Writing the SAM header and file**: 0.74-0.81 s on q512n10 (0.7.7, also 0.17-0.95 s for 0.8.2 at 8 threads),
  1.63-6.73 s on q512n19. The SAMs are the same size (1.66 vs 1.70 GB of text, header placed in its room both
  times), so the time is the network file system's flush at `close`. Run 2's alignment stalled the same way. For
  benchmarks, give `measure_performance.sh` an `OUT_DIR` on node-local disk and copy it back afterwards. For
  users, the close could run on a thread of its own while the profiling reads the SAM from the page cache, with
  its result checked before the run reports success. That saves 1-6 s whenever the output folder is on NFS.
- **Cold start** (not in the comparison): the first run of each sample loaded the index in 935 s (369 s and 227 s
  in earlier runs). That is the 30 GB database read from NFS at ~32 MB/s, and on a cluster it is the largest wall
  time a user sees on a fresh node. Copying the database to node-local disk, or reading it once ahead of the job,
  avoids it.

## 4. What each feature costs against what it adds

F1: change when the group is left out of the model (v23 refits, species held out, pe / se / pb / ont; negative =
the group helps). Wall time: pe1 on 32 threads once the parse scales, estimated from the v23 report's instruction
shares (section 5 there, adjusted for the prebuilt `ancestry_sites.bin`, which removed 47% of the profiling
instructions) and the 8-thread parse above. Memory as the run logs it.

| feature (group) | where it runs | wall time, pe1, 32 threads | memory | F1 left out (pe / se / pb / ont) |
|---|---|---|---|---|
| ancestry sites: the ancestry group and the fixed-site columns (8 columns) | profiling: decode a copy's sites once, walk each record's sites | ~0.1-0.3 s | `ancestry_sites.bin`, 11.9 M copies | −0.0015 / −0.0009 / −0.0001 / −0.0019 |
| ancestry + polymorphic + weights together | profiling | ~0.3-0.5 s | as above + column weights | −0.0022 / −0.0016 / −0.0013 / −0.0031 |
| column weights (weights, 5 columns) | profiling: expand a copy's columns once, count each record | ~0.1-0.2 s | **1.47 GB** | alone −0.0001 / −0.0001 / −0.0003 / +0.0000 |
| polymorphic sites (part of the alleles' table) | profiling | < 0.05 s | none of its own | alone −0.0002 / +0.0000 / +0.0000 / −0.0006 |
| strain alleles (alleles; also the alignment's allele scores) | alignment (0 measurable) and profiling (~1% of its instructions) | < 0.1 s | 9.2 M alleles | alone ~0; with the site groups −0.0026 / −0.0021 / −0.0020 / −0.0032; strain calls 0.830 → 0.903 |
| congener distances (distance) | profiling, per sample | 0.23 s | sketch cache | −0.0005 / −0.0006 / +0.0000 / −0.0001 |
| congener gaps (gaps), untried candidates | profiling | < 0.02 s | `congener_gaps.tsv`, 13 M copies | +0.0003 / −0.0001 / +0.0001 / +0.0000 |
| gene neighbours (adjacency; also mate guidance) | alignment and profiling | ~0 (0.06% of profiling) | 1.8 M rules | +0.0000 / +0.0000 / −0.0002 / +0.0002 |
| adaptive candidates | alignment: 206,699 of 176 M candidates (0.12%) | ~0.03 s | none | not measured alone |
| everything since the v9 core set | | ≤ ~1.4 s in profiling, ~0 in alignment | 1.47 GB of column weights plus the other tables (not measured apart) | **−0.0055 / −0.0051 / −0.0042 / −0.0058** |

How to read it:
- The site groups stand in for each other. Each alone is worth ~0, but the three together are worth about as
  much as v22 → v23 gained. Cutting one of them saves at most a few tenths of a second and moves F1 within noise.
  Cutting all of them saves ~0.5 s and costs 0.002-0.003 F1.
- PacBio and ONT spend 1-4% of their CPU in profiling (pb2: 0.9 s of 13.5), so no feature cut changes their wall
  time noticeably.
- Every cut needs retrained models (the r226 v24 pe model uses features of every table), so a cut belongs with the
  next database build, not with a protal release on its own.

## 5. Recommendations, in order

1. **Shard the run-wide caches** (section 2). Expected −3 s on pe1 at 32 threads (−0.1 s pb2), and the parse
   scales with threads again. Profiles stay byte-identical. Check that the "SAM read" line's parsing goes back to
   ~15-25 thread-seconds and system time drops toward 25 s.
2. **Time the profiling tail and skip `malloc_trim` when no sample is left** (section 3). Up to −1.8 s per run.
3. **Benchmarks: write `OUT_DIR` to node-local disk.** This removes 1-6 s of noise in the SAM finish and in the
   alignment's output. A background close of the SAM would give users the same on NFS output folders.
4. **Cut no feature for wall time.** If memory or build time matters, the candidates are:
   - column weights (1.47 GB resident, about 3 minutes of the training database's build, ~0 F1 alone);
   - congener gaps (~0 F1).

   Decide with a refit without them at the next build. Expect less than 0.3 s of wall time back.

After 1-3 the 32-thread pe1 run should be at ~33-34 s with every feature kept (about 37 s if the SAM close on
q512n19 stays as it was), against 41.65 s now.

## 6. Follow-up (the same day): what was implemented

The user asked for recommendations 1 and 2 (the timer, not yet the `malloc_trim` change), a `--scratch_dir` for
protal's intermediate files (replacing `--strain_spill`), and a folder for the benchmark's runs, sylph's in the same
folder. Committed on `e751c7f` in one commit; the benchmark script's `RUN_DIR` went into `e751c7f` itself, which
another session committed with the whole script.

- **Sharded caches.** `src/Utilities/ShardedCache.h`: 64 shards with a lock each, keyed by the mixed key's top 6 bits;
  each shard starts over at 1/64 of the old cap. It holds `ancestry::Cache` (`AncestrySites.h`),
  `column_weights::Cache` (`ColumnWeights.h`) and `GenomeLoader`'s species sites. `SpeciesAncestrySitesOf` takes the
  sites the caller already has (`NoteRecord`), one lookup fewer per new copy. A cached value depends on its key only,
  so the profiles cannot change; the unit and e2e suites pass on the snapshot (HEAD + these changes). The contention
  itself shows only on many cores: the next cluster run should show the "SAM read" line's parsing at ~15-25
  thread-seconds and the system time near 25 s.
- **The profiling's tail.** A runtime-table row "Profiling: after the profile" and a line per sample,
  `Sample <name> after its profile: T s (the runtime table and lines, unreported species, releasing the read data,
  releasing the record evidence, [spilling the strain evidence,] returning freed memory (malloc_trim))`. Whether to skip
  `malloc_trim` for the last sample waits for those numbers.
- **`--scratch_dir DIR`** (`src/Utilities/ScratchSpace.h`, `Options::UseScratchFolder`): a folder
  `protal_<host>_<pid>_<ms>` made in DIR at the start and removed when the process ends (normal end or `exit()`; a
  killed run leaves it). In it:
  - each sample's strain evidence, read back per species (what `--strain_spill DIR` did; that option is hidden and now
    means `--scratch_dir`);
  - with `--no_strains`, the SAMs protal names: written there and removed once their sample is profiled. Not a map's
    `SAM` column, not `--no_profile`, not `--profile_only`, and not a SAM an earlier run left in the output folder,
    which is profiled as before;
  - the records of a `--sam_format gz`/`sam` SAM before its header (a zstd SAM whose header outgrows its room still
    renames them beside the SAM);
  - qcmsa's captured output, with `TMPDIR` pointing there.
- **`scripts/measure_performance.sh`**: `RUN_DIR` (environment, default OUT_DIR) holds protal's runs in
  `RUN_DIR/protal` (each run's `-o` folder and the cohort's SAMs) and sylph's profiles in `RUN_DIR/sylph` (copied to
  `OUT_DIR/sylph` at the end); `environment.txt` names its file system. The `--test_sparse` runs go to `RUN_DIR/protal`
  as well.
