# The test suites made faster, their checks kept

**Date:** 2026-10-12.

**Code:** `audit-fixes` at `4a42f3b` (protal 0.8.2; `714649b` after it changes only reports). The changes are in the
commit that adds this report.

**Earlier reviews:** [2026-10-06 audit](../2026-10-06-test-suite-audit/README.md) and
[2026-10-07 fixes](../2026-10-07-test-suite-fixes/README.md). On 2026-10-07 every suite ran in about 4 minutes on
4 cores. Since then 120 unit tests, 30 end-to-end tests, the ancestry true-positive test, scenarios, streaming, strain
alleles and NCBI taxids were added, about 12,000 lines of tests.

**Machine:** WSL Ubuntu 24.04 on the shared Windows laptop (Core Ultra 7 258V). Every build and suite ran on CPUs
0-3 (`taskset -c 0-3 nice -n 5`). Python 3.14.7 with scikit-learn 1.9.1, numpy 2.5.3 and pandas 3.0.6 (the
micromamba environment `protal-db-build`). CI has Ubuntu's Python 3.12 with scikit-learn 1.4.1, not installed here.
Other sessions compiled, ran callgrind and trained models on the same machine during most of the runs, so wall times
vary up to 2-2.5×. The comparison below alternates the old and the new tree to even that out.

**Scripts:** [`scripts/`](scripts/). `run_all.sh` runs every suite as CI does, timed. `timed_unittest.py` records
each test's and each class's wall time (`setUpClass` included). `ab_all.sh` alternates the old and the new tree. The
rest are the single measurements quoted below. They use the scratch paths of the session that wrote them.

## Answer in short

The tests took about 15 minutes on 4 cores at `4a42f3b`. Two causes made up most of it:
- scikit-learn 1.9.1's slow binning of weighted fits, in one trainer test (5 minutes);
- the GTDB build end to end (4 minutes).

The changes keep every check, and add the two script tests that CI never ran. What they remove is repeated work:
- builds made twice: the mini database in CI, the release of two end-to-end classes;
- databases rebuilt where a copy serves: the GTDB build's streamed and followed collections;
- level-19 compression of throwaway test databases;
- the trainer test's slow binning path;
- serial waits: the GTDB runs on 2 of 4 threads, a determinism check generated after its first world.

| | before (quiet) | after |
|---|---|---|
| Trainer and other script tests | 344 s | 78-118 s (loaded rounds; 1062 s for the old tree in round 1) |
| `TrainerScenariosTest` | 312 s | 69 s (load 6.8) |
| Mini database + example in CI | 20 + 20 s, two builds | 16-21 s, one build (91 s for the old tree in round 1) |
| `GtdbBuildTest` `test_g` + `test_h` | 31 + 31 s | 11 + 11 s (load ~5) |
| `GtdbBuildTest` in all | 252 s | 226 s (load ~5) |
| Unit tests | 7-11 s | unchanged |

Taken together, the new tree ran every suite in 15.5 minutes at load averages of 3-5 (round 2). The old tree took
36 minutes at 5-7 (round 1), and about 15 minutes when the machine was quiet. A quiet run of the new tree was not
made: the laptop went on battery. On a quiet machine, the savings above add up to about 5-6 minutes, from 15 to about
9-10. The trainer test is the largest share, but it is local only: CI's scikit-learn 1.4.1 never had that cost.

Found on the way: scikit-learn 1.8 and later makes every weighted fit of the trainer slow in real GTDB builds too:
268 s against 35 s for one fit of 200,000 rows (section 1). That fix belongs in the trainer and is left as a task of
its own.

## Where the time went

The suites at `4a42f3b`, run as CI runs them, on a quiet machine (load average below 1.5):

| Suite | Tests | Time | The slowest classes |
|---|---|---|---|
| C++ unit tests (`ctest -j4`) | 516 | 11.5 s (6.9 s in an earlier quiet run; 15.5 s serially in one process) | `IlluminaSimulation` 2.7 s, `ThreadedGzStream` 1.8 s, `ProfileSam` 1.4 s (serial) |
| Mini database (`build_mini_db.sh`) | | 20 s | `protal --build -t 1` at zstd level 19: 13 s |
| End-to-end | 155 | 80 s | `CompressedDatabaseTest` 20 s, `BuildIndexTest` 9 s, `GeneNeighboursTest` 9 s, `WfaTracebackTest` 8 s |
| `examples/mini_db/run.sh` | | 20 s | a second build of the same mini database |
| `scripts/mini_db/test_*.py` | 114 | 435 s | `GtdbBuildTest` 252 s, ancestry `PipelineTest` 74 s, `MiniDbTest` 30 s, `WorldTest` 23 s, `GeneSubsetTest` 22 s |
| Other script and trainer tests | 98 | 344 s | `TrainerScenariosTest` 312 s, of it `test_features_auto` 245 s |

That is about 15 minutes. Three causes account for most of it:
- scikit-learn's binning of weighted fits (section 1);
- the GTDB build end to end, which built 13 databases and ran each pipeline on 2 of the 4 cores (section 2);
- databases built at zstd level 19 or on one thread, and builds made twice (sections 3-5).

## 1. The trainer's tests: scikit-learn's binning of weighted fits

`TrainerScenariosTest` trains on tables with scenario rows, which the trainer weighs 0.25 (`--scenario-weight`).
cProfile of the trainer's `--features auto` run on the test's tables (`prof_auto.py`): 402 s, of which 395 s are in
`sklearn/ensemble/_hist_gradient_boosting/binning.py::_find_binning_thresholds`, in 1.86 million calls of
`sklearn.utils.stats._weighted_percentile`. With a sample weight and a column of more than 255 distinct values,
scikit-learn 1.8 and later computes the column's 254 bin thresholds one percentile at a time, each call sorting the
column again, although `_weighted_percentile` takes an array of percentiles. One fit of the test's 720 rows took
3.2 s.

**Change:** the test's synthetic features are rounded to steps of 0.1 (about 70 distinct values per column). A column
of 255 values or fewer is binned at the midpoints between them, weighted or not. The rows keep their weights. Every
assertion is about the procedure (which set is chosen and why, the scenarios' counts, that weight 1 gives another
model, parity with the PMML), not about values.

| | before | after |
|---|---|---|
| `test_features_auto` | 245 s | 22 s (under load 6.8) |
| `TrainerScenariosTest` | 312 s | 69 s (under load 6.8) |

CI's scikit-learn 1.4.1 bins without the weights, so CI did not pay this; a machine with a current scikit-learn
(the build environment here) did.

**The same cost in real trainings.** One fit on one thread, 90 features, 30 rounds (`weighted_binning.py`):

| rows | unweighted | weighted (a fifth at 0.25) |
|---|---|---|
| 2,000 | 1.2 s | 25 s |
| 20,000 | 5.8 s | 53 s |
| 200,000 | 35 s | 268 s |

A GTDB build's trainers weigh the scenario rows (the default since `f12cac1`) and fit many times (folds, the feature
sets of `--features auto`, the final model). With scikit-learn 1.8 or later, which `protal-db-build` has, every such
fit pays this. The fix belongs in the trainer: compute the same thresholds with one call per column. It is left for a
task of its own, as it changes the training code. It must give bit-identical models, and CI's 1.4.1 has no weighted
binning to patch. (Done the same day: [2026-10-12-trainer-weighted-binning](../2026-10-12-trainer-weighted-binning.md).
A weighted fit now takes what an unweighted one does, with the same thresholds and probabilities.)

## 2. The GTDB build end to end (`GtdbBuildTest`)

Its runs used 1.6 of the 4 cores on average (`keep_gtdb.py` under `/usr/bin/time`: 392 CPU-s in 245 s). Each
pipeline run converts the release, builds a training and a finished database, simulates, profiles and trains, so
every test that runs one costs about 30 s.

**Changes:**
- `test_g` (profiled as simulated) and `test_h` (streamed) run in copies of the full build without its collections
  (`copy_of_full`: the training data, the test set and the genome store removed). The build script's stages are keyed
  by content, not by OUTDIR (`test_c` already relied on this), so the run takes the conversion, both databases and
  their foreign scans as they are, and simulates and profiles again. The in-silico strains are made again: their
  stage key names the scratch folder. What the tests check is unchanged: the same tables as the full build's, the
  reads removed, the streams' folders, the genome store. Both now also assert that nothing was built. The full
  build's state is kept as `full_snapshot`, which `test_c` copies too (it used to run in the only copy, which it
  changes).
- The runs use 4 threads instead of 2. Alternated twice (`test_f` and `test_i`, setup included): 88 s against 79 s.
- `test_b` and `test_i` do not scan the full references for foreign rates; they never looked at the scans.

| | before | after |
|---|---|---|
| `test_g_profiled_as_simulated` | 31 s | 12.5 s |
| `test_h_streamed` | 31 s | 15.7 s |

## 3. The ancestry true-positive test (`test_ancestry_world.py`)

`PipelineTest` runs `ancestry_truth_test.py run --quick`. Its step times (`.done` markers): the two `protal --build`
runs took 35 and 21 s, at `--build`'s level 19; the two worlds 16 s each, one after the other on one core each.
`WorldTest` generated its world twice in a row, to compare them.

**Changes:**
- `--quick` (the tests' mode) builds at zstd level 1: the same database content. The test world's build spent 24 s
  of its 35 s writing `database.protal` at level 19 (`build.log`).
- `WorldTest` starts the second generation in `setUpClass`, beside the first; the determinism test waits for it.

**Tried and taken back:** `--quick` making the two worlds side by side. It needs 6.5 GB of protal memory at once
(`memwatch.sh`), and in both loaded rounds the new `PipelineTest` was slower than the old (165 s against 136 s in
round 1). Without a quiet comparison it stays one world after the other.

## 4. The end-to-end tests

**Changes:**
- `GeneNeighboursTest` and `WfaTracebackTest` each simulated a release of 200 kb genomes, converted and built it, and
  simulated reads (8.7 and 8.0 s). `WfaTracebackTest` needs only a database built with the default traceback and
  reads that align. Both now share one module-level `OperonRelease`, made on first use. Five builds remain in the module.
- `CompressedDatabaseTest` (20 s) writes a raw index of 3 GB three times and hashes it three times. Timed step by step
  (`cdb_steps.sh`): `--decompress_db` 3.4 s, sha256 2.2-2.8 s, each `--compress_db` 1.2-1.8 s. sha256 is the fastest
  hash Python has here (`hashspeed.sh`: crc32 2.0 s, sha1 3.2 s, md5 7.0 s, blake2b 7.8 s). Now the raw index is
  hashed on a thread while the copies are made (about 2 s). The copies compress at level 1 instead of 3 (about 1 s):
  the 1 MB frames are what the test needs.

  **Tried and taken back:** the two round trips side by side, each on 2 threads. Peak protal memory went from 3 to
  6.2 GB and `/tmp` to 6.5 GB, and under load the class was not faster (48 and 69 s against 42 s in the loaded
  rounds).

## 5. The mini database and CI

`build_mini_db.sh` built on one thread at `--build`'s level 19: 13 s (`minibuild.sh`, alternated twice). On 4 threads
at level 3: 4 s, the file 1.6 instead of 1.1 MB (level 1 is no faster and larger). CI built it twice: once for the
end-to-end tests and once in `examples/mini_db/run.sh`, from the same seed.

**Changes:**
- `build_mini_db.sh` builds with `-t 4 --compress_level 3`; `PROTAL_BUILD_ARGS` still overrides both.
- CI runs `examples/mini_db/run.sh data/mini_db` first: it builds the mini database where the end-to-end tests look
  for it and checks the accuracy, in one step instead of two builds. In the first round of the comparison: 21 s
  instead of 57 + 34 s.
- `GeneSubsetTest`'s three builds use 4 threads.

## 6. Tests CI did not run

`scripts/test_ancestry_sites.py` (8 tests) and `scripts/test_foreign_rates.py` (2 tests) ran in no CI step, and the
first in no `just` recipe. They take under a second and need numpy only. Both are now in CI's script step and in
`just mini-db-test`.

## What was left as it was

- **The C++ unit tests** (516, 7-11 s with `ctest -j4`, 14-16 s in one process). The slowest test takes 1 s; the
  heaviest are threads-against-one comparisons and the Illumina simulator's pipes. The 51 death tests cost about
  40 ms each (fork and exit). No test does work twice, and trimming read counts would save under 2 s.
- **`MiniDbTest`** (30 s): many conversions with different options, each 1-1.5 s, some of which could run side by side
  (about 5 s). `test_genome_length_range` (4.6 s, genomes up to 400 kb) asserts on three seeded draws; a smaller
  range would make it depend on the draws.
- **Running the Python modules side by side in CI.** `GtdbBuildTest` builds its finished database in the background
  at the idle scheduling class. Any other load starves that build: with another session's training running, it took
  2 minutes instead of about 15 s.
- **`TrainerDepthKnobsTest`** (19 s): four trainer runs of different options, 3-8 s each, a second or two of each
  importing scikit-learn.

## How the numbers were taken

```bash
bash scripts/run_all.sh SRC BUILD RESULTS unit mini e2e example minidb scripts    # the old tree, as CI ran it
bash scripts/run_all.sh SRC BUILD RESULTS unit ciexample e2e minidb scripts       # the new tree, as CI runs it now
bash scripts/ab_all.sh                                                            # both, alternated: old new new old
```

The baseline at `4a42f3b` ran from 07:22 to 07:38 at load averages of 1-3 (the build before it, then quiet). The
comparison ran from 08:30, while other sessions compiled, ran callgrind and their own end-to-end tests:
- old round 1 at load averages of 5-7;
- new rounds 1 and 2 at 3-6;
- old round 2, stopped in its end-to-end tests at 10-12.

Every suite of both new rounds passed. Per class, the loaded rounds are noise: `DownloadTest`, whose code did not
change, took 8 s quiet and 12-32 s in them. Only the large differences above stand out of that.

**The website** is not affected: nothing here changes how protal runs.
