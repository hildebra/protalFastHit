# The test audit's recommendations implemented: what changed, what it costs now, what is left

**Date:** 2026-10-07.

**Code:** `audit-fixes` from `ce85bd7`. The changes are four commits:
- `049f7a1`: C++ unit tests;
- `ac32e2d`: end-to-end tests;
- `7086c09`: script tests;
- `51d111d`: CI and docs.

The six unverified WIP branches of 2026-10-06 (the audit's follow-up section) were the starting point. They were merged, built, run, fixed and finished here.

**Audit:** [2026-10-06-test-suite-audit](../2026-10-06-test-suite-audit/README.md). Its section 7 lists recommendations 1-8; the user asked for all of them, and for long tests to be cut.

**Machine:** WSL Ubuntu 24.04 on the shared Windows machine. Every build and test ran on CPUs 0-3 (`taskset -c 0-3`, nice 5), as the user asked: at most 4 cores. Python 3.12 with scikit-learn 1.4.1, pandas, joblib and numpy, the versions Ubuntu 24.04's packages have, which CI installs.

**Verification:** [`scripts/verify.sh`](scripts/verify.sh) runs every suite as CI does on the content of `51d111d` (`git archive`), with `PROTAL_TESTS_REQUIRED=1`. Its output is [`results/verify.out`](results/verify.out); the sanitizer run, by [`scripts/asan.sh`](scripts/asan.sh), is [`results/asan.out`](results/asan.out). Both scripts use paths in the scratch directory and in `~/testaudit` and `~/ta` of the session that wrote this.

## Answer in short

All eight recommendations are implemented, and every suite now runs in CI. On the content of `51d111d`, run as CI runs them:
- the unit tests: 396 in 5 s;
- the end-to-end tests: 125 in 49 s;
- the accuracy example: 15 s;
- the mini database and GTDB build tests: 71 in 136 s;
- the other script and trainer tests: 67 in 20 s.

Nothing was skipped and all passed, about 4 minutes of tests in all on 4 cores. With the audit's numbers, the same suites took about 30 minutes, 20 of them `test_model_pmml.py` with scikit-learn 1.9.1, and most of them ran nowhere automatically. A few things are left; they are listed at the end.

| Suite | Before (audit, or same machine where noted) | Now | In CI |
|---|---|---|---|
| C++ unit tests (`ctest -j4`) | 372 tests, 54 s; one test peaked at 13.2 GB (two raw 3 GB indexes) | 396 tests, 5-9 s; that test builds 3 MB key maps, 0.3 s | yes |
| ASan + UBSan (Debug, `-j4`) | not timed in the audit | 95 s, now with `-DPROTAL_NO_CLONES` | yes |
| End-to-end | 134 tests, 154-161 s (same machine and binaries) | 125 tests, 49 s | **new** |
| `examples/mini_db/run.sh` (accuracy) | not run since 2026-09-26 | 15 s, abundances within 0.002 | **new** |
| `scripts/mini_db/test_*.py` (GTDB build end to end included) | 70 tests: 65 s, plus `GtdbBuildTest` 246-354 s (here, under load) | 71 tests, 136 s in all (`GtdbBuildTest` 149-196 s under the same load) | **new** for the GTDB build |
| Other script tests and the trainer's | 52 tests; `test_model_pmml.py` 55 s here (20.5 min with scikit-learn 1.9.1 in the audit) | 67 tests, 20 s; `test_model_pmml.py` 28 s | **new** for the trainer |

Without a prerequisite a test is skipped. With `PROTAL_TESTS_REQUIRED=1`, as in CI, it fails and names what is missing. Before, the end-to-end suite reported "OK" having run nothing, and two Python suites errored.

## 1. Run what exists, and make missing prerequisites fail loudly (recommendation 1)

- **CI's Release job** now also runs:
  - the mini database build;
  - the end-to-end tests;
  - `examples/mini_db/run.sh`;
  - every script test, including the GTDB build end to end;
  - the trainer's tests.

  It installs `python3-pandas`, `python3-sklearn`, `python3-joblib` and `art-nextgen-simulation-tools`, and sets `PROTAL_TESTS_REQUIRED=1`.
- **`PROTAL_TESTS_REQUIRED=1`** (the audit proposed `PROTAL_E2E_REQUIRED`; one variable serves every suite) turns a missing prerequisite into an error. It is implemented in `tests/e2e/test_protal_e2e.py` and in `scripts/prerequisites.py`.
- **Three modes were checked:**
  1. everything present with REQUIRED: all pass;
  2. system `python3` without scikit-learn or `$PROTAL`: clean skips;
  3. REQUIRED with something missing: failures that name it.

  The end-to-end tests that need no database (the version, the simulator, qcmsa's contract, small builds, gene neighbours) run without one.
- **The gates the audit found wrong are fixed:**
  - `GtdbBuildTest` errored without `$PROTAL`.
  - `TrainerDepthKnobsTest` errored without scikit-learn.
  - The trace tests were gated on tools they do not use.
- **`test_model_pmml.py`** fits on one OpenMP/BLAS thread and trains each table once.

## 2. Assertions that could not fail (recommendation 2)

Every row of the audit's section 2 is fixed or removed. The confirmed ones:
- **`test_SampleContext.cpp`:** the veto checks now run with the singleton rule on. They cover the own read, the identity branch with 6 against 4 mismatches, and rules 0 and 201. The references' distance (0.04562) and `relative_spill` are pinned by number.
- **`test_CommunityDesign.cpp`:** bound at 30 (twice the expected 14.7), and grouped pairs above twice the uniform ones.
- **`test_Zstd.cpp`:** what a truncated or corrupt file gives is a strict prefix of the content.
- **Pipe tests:** the stream is closed before `join()`, and SIGPIPE is ignored, so a reader that stops early fails its test.
- **End to end:**
  - the truth line asserts `TP 3, FP 0, FN 0` (and a thin sample's `TP 2, FP 0, FN 1`);
  - the MSA-knob tests are independent of each other and need species;
  - the knob tests use a sample whose profile the knob changes;
  - no check passes on an empty glob.

The constants behind default features (`kSpillAtZero`, `kSpillDecade`, the edit ratio's bounds) are now pinned by their numbers, not by their symbols.

## 3. Accuracy and the default model (recommendation 3)

- **Accuracy:**
  - The end-to-end suite checks truth counts and abundances (within 0.01 of the simulated pairs).
  - Reads of nothing in the database give an empty profile.
  - `run.sh` runs in CI and exits 1 when accuracy is off; it was checked with `--knob 1`.
- **A gradient-boosted `modelChain`** (`tests/e2e/data/model_gbm_small.xml`) is scored by protal end to end:
  - every leaf is reached;
  - each probability equals one computed by hand (1e-12);
  - the trainer's scorer agrees;
  - the fixture is what the exporter writes.
- **The feature contract:** every group of `FEATURE_GROUPS`, the default set and the model's fields are checked against protal's feature dump.
- **Golden vectors** (`tests/data/golden_model_rules.tsv`): the prior adjusted to a sample, the calls at a share of false calls, and the knob by depth. `test_GoldenModelRules.cpp` checks protal on them, `test_model_pmml.py` the trainer. Before, each side had hand-copied numbers of its own.

## 4. Cost (recommendation 4, and the user's "remove lengthy tests")

- **C++:**
  - The packed-index tests build `Seedmap`s of 10-base cores through a test-only constructor (`Seedmap(size_t exact_k)`, a 3 MB key map instead of 3 GB). `ABuildIntoThePackedLayoutWritesTheSameIndex` went from 53 s and 13.2 GB to 0.3 s.
  - The flex-scan bench moved to [`flex_scan_bench.cpp`](../2026-10-06-performance-profiling/scripts/flex_scan_bench.cpp).
  - The gene-table bench was removed; it stays in `git show ce85bd7:tests/test_InputValidation.cpp`.
  - `IndexCodec.ColumnsCompressBetterThanRawCells` uses level 3 on a smaller index. Its ratio still holds.
- **End to end:**
  - one paired-end alignment of the module's samples, which tests that change only profiling profile again (`--profile_only`);
  - builds at `--compress_level 1` on 4 threads;
  - two builds where five were (ambiguous bases, uniqueness, parallel batches);
  - the parallel build compared by its compressed index, not two 3 GB raw ones hashed;
  - `--add_model` with the small model;
  - no sleeps.
- **Python:** in `GtdbBuildTest` each build serves every check of what it does. Old `test_d` and `test_e` were folded into `test_a` and `test_b`, the seed-2 run works on a copy of the full build, and every build but the scenarios passes `--error-reads none`. **Dropped on purpose:**
  - the reduced database's rerun;
  - the `--genes` pipeline run (the converter's `--genes` is tested on its own);
  - the scenarios in follow mode.

## 5. Dead code and duplicates (recommendation 5)

- **Removed with its tests:** `Compressor.h`, `bgzf::CompressFile`, `zstd::MemoryReader`, `AlignmentEdits(std::string)`, `HoldsPairedReads`, `CalculateCoverageVector()`, `ModelContractProblem(model, path)` and `scripts/index_reference.py`.
- **Moved:** `SamStreamSink`, which only tests used, into `tests/TestSamSink.h`.
- **Rewritten:** `BgzfWriter.PiecesOfAnySizeGiveTheSameFile` replaces the test whose oracle was the dead `CompressFile`. It compares piece sizes with each other, walks the BGZF blocks and reads them back with zlib-ng.
- **Duplicates:** every row of the audit's duplicate table is merged or removed, including the two `IndexHeader` tests, now one in `test_Index.cpp`. The end-to-end duplicates are merged into parameterised tests.

## 6. Reorganisation (recommendation 6)

- **Shared helpers:**
  - `tests/TestUtil.h` holds the scratch directory, whole files, random and mutated sequences, and `TestData`.
  - `tests/TestReference.h` holds taxa of genes written as `reference.fna`/`reference.map` and loaded.
  - They replace 13 copies of the scratch directory and 10 of the reference writer. Each object now has its own directory; two `TinyReference`s used to share one.
- **Suites:**
  - `GeneConservation` and `ModelFeatures` have files of their own.
  - `test_Operations.cpp` is now `test_RunStatusAndBgzf.cpp`.
  - Tests in the wrong file moved (`SeqReader`, `ExtractVariants`, a `SamReader` test).
  - Suites at different layers that shared a name were renamed (`SamPairs`, `NeighbourFeatures`, `FlexNeighbours`).
- **Python:**
  - `test_mini_db.py` is split into six files by subject, with `mini_db_fixtures.py`.
  - The console regexes with step counters are now checks on files.

## 7. Gaps closed (recommendation 7)

- **The sanitizer job builds with `-DPROTAL_NO_CLONES`.** It runs the plain x86-64 copies of the dual-compiled functions under ASan and UBSan, at no extra CI time; the Release job runs the x86-64-v3 copies.
- **Thread identity:** 10 genes per taxon, congeners and a conservation table, so the distance, conservation and failed-candidate features differ from their defaults and are compared across threads and chunks.
- **The alignment screen through `operator()`,** with the shared `ReadKmers` strands and about half the reads reverse.
- **Error paths:**
  - `Bundle::Open`'s directory checks;
  - `/dev/full`;
  - `LoadFromMap`;
  - an edited `reference.map` of the same size (pinned as accepted);
  - `SortKey` at its maxima.
- **Code that had no unit test:**
  - `LongReadAligner`;
  - the `Haplotypes` settings;
  - `recurrent_calls.py`, `prevalence_calls.py` and `protal_profile_utils`;
  - `strain_test/db_gene_counts.py` and `strain_report.py`.

## 8. Documentation (recommendation 8)

`docs/development.md` now describes the unit tests as they are, the end-to-end tests (written for the mini database, with what they need), the script tests by file, the redesigned GTDB build test, and CI. The analysis script `test_errors.py` is renamed `error_breakdown.py`, so that `pytest` does not collect it.

## A bug the tests found

**`build_gtdb_database.py --rank-genes`** ranked a folder unpacked from `database.protal`. That folder lacks the build's report `gene_congeners.tsv`, which sits beside the bundle, not in it, so `gene_ranking.tsv` had NA for `between_factor` and `identical_share`. The report is now copied in, and `test_a` compares that ranking's whole text with `--n-genes`'.

## What is left, or narrower than before

- **Not tested anywhere:**
  - an interrupted run resumed;
  - `CoverageVector` with a read outside its range: Debug builds assert, a Release build would still write out of bounds;
  - the duplicate-column and empty-key branches of `LoadFromMap`.
- **Helpers still copied:**
  - `MakeSam` (4 files, different fields);
  - `Random` (3);
  - `RandomBases`;
  - `ExactRuns` (2);
  - `MateGuidance`'s two tests stay in `test_SamRoundTrip.cpp`, beside the helpers they use.
- **Full 15-base `Seedmap`s** (3 GB `calloc`, untouched) remain in a few tests that do not need 2^30 keys: `test_GeneWindows`, `test_Index`, `test_InputValidation` and `test_Zstd`. They cost no memory, but under ASan they are the slow tests' common factor.
- **End to end:**
  - the parallel-index check runs on the tiny build reference instead of the mini database;
  - the "an ambiguous base leaves the other genes alone" comparison is gone (the `AmbiguousKmers.*` unit tests check the windows);
  - `Failure.test_truncated_index` is dropped as a near-duplicate.
- **`test_model_pmml.py` with scikit-learn 1.9.1** was not timed here: the machine has only 1.4.1, which CI installs too.

**The website** is not affected: nothing changes how protal runs.
