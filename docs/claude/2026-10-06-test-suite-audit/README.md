# The test suites: what they check, what is redundant, what runs nowhere

**Date:** 2026-10-06.

**Code:**
- The test files of the working tree on `audit-fixes`: HEAD `3bcf535` plus other sessions' uncommitted edits.
- `tests/test_SharedSeeds.cpp` is untracked. `tests/CMakeLists.txt` already lists it but is uncommitted too, so the two must be committed together.
- The measurements are from `git archive 3bcf535`, which has 372 C++ tests instead of the working tree's 376.

**Machine:** WSL Ubuntu 24.04, Python 3.12 with numpy, no scikit-learn, art_illumina installed. The machine was shared with other sessions (load average 7-15), so times are indicative only.

**Scripts:**
- [`scripts/run_all_tests.sh`](scripts/run_all_tests.sh) builds HEAD and runs every suite with timings.
- [`scripts/mem_tests.sh`](scripts/mem_tests.sh) measures the peak memory of the slowest unit tests.
- Both use paths in the scratch directory and in `~/testaudit` of the session that wrote this.

**Results:** the timings and logs are in [`results/`](results/).

**Method:** six reviewers each read one group of test files in full, with the code under test. Their findings were checked against the code; the ones marked ✓ below were confirmed by hand, the rest come from a reviewer's reading with file and line. One of the six also ran `test_model_pmml.py` with scikit-learn 1.9.1 and `GtdbBuildTest` without `$PROTAL`.

## Answer in short

**The C++ unit tests mostly make sense.** Most of them compare the code with an independent oracle:
- libzstd and zlib-ng's reader;
- WFA2 alignments;
- brute-force definitions of syncmers, flex neighbours and the alignment screen;
- the text loaders against the binary gene table;
- serial against parallel profiling.

Their seeds are fixed. The AVX2 paths are switched at run time and compared with the scalar ones.

**Redundancy is moderate, not the main problem.** About 20 of 376 tests duplicate another test outright. On top of that there are:
- eight suites split over two or three files;
- the same helpers (`ScratchDir`, `TinyReference`, `Reference`) copied into up to 13 files;
- in the end-to-end suite, the same sample aligned about 12 times.

**The larger problems are these:**
1. **What runs where.**
   - Nothing runs automatically: the end-to-end suite (133 tests), the model tests (37), the GTDB build test (7) or the accuracy example.
   - CI builds `protal` and `simulate_metagenomes` and then never runs them.
   - The end-to-end suite reports "OK" having run nothing when `PROTAL_TEST_DB` is unset.
   - Two Python suites fail with errors instead of skipping when a prerequisite is missing.
2. **About 15 assertions are vacuous or wrong.** They pass whatever the code does, or pass only because of their seeds (section 2).
3. **No automated test checks that a profile is right, or loads the default model.**
   - The mini database ships the 2024 random forest (`scripts/random_forest.xml`), so the gradient-boosted `modelChain` default of `1750475` is loaded only by `GtdbBuildTest`, which runs nowhere.
   - Truth TP/FP/FN numbers and abundances are checked only by `examples/mini_db/run.sh`, last run on record 2026-09-26 with 0.6.0.
4. **Cost hot spots.**
   - One unit test peaks at 13.2 GB of memory.
   - The `PackedIndex` suite takes 68% of the unit-test time.
   - `test_model_pmml.py` takes 20 minutes. Its gradient-boosting fits use every OpenMP thread; with one thread it takes 6.

## 1. Inventory and measurements

| Suite | Tests | Time here | In CI | Needs; when it is missing |
|---|---|---|---|---|
| C++ `protal_tests` (`tests/*.cpp`, 34 files) | 372 at HEAD, 376 in the working tree; 3 always skipped (2 benchmarks, 1 for AddressSanitizer only) | ctest `-j4`: 54 s wall (177 test-seconds, `PackedIndex` 120 of them). One process: 534 s, `PackedIndex` 416 s (memory pressure on the shared machine) | yes: Release, and Debug with ASan + UBSan | gtest |
| `tests/e2e/test_protal_e2e.py` | 133 | 143 s, plus 20 s to build the mini database | **no** (`just e2e`) | `PROTAL_TEST_DB`, both binaries, zstd CLI, numpy, Linux. Without the database: `Ran 0 tests … OK (skipped=1)`, exit 0 |
| `scripts/mini_db/test_mini_db.py` | 63 | 65 s | yes, except `GtdbBuildTest` and the 2 tests that need `$PROTAL` | numpy, art_illumina |
| ├ `GtdbBuildTest` | 7 | minutes (3 of the tests: 114-153 s on 2026-09-30) | **no** | skipped without scikit-learn, but **7 errors (`KeyError: 'PROTAL'`) with scikit-learn and no `$PROTAL`** ✓ (the gate at :2855 checks only the Python modules) |
| `scripts/test_insilico_strains.py`, `scripts/test_trace_relatives.py` | 11 | 0.1 s | yes | numpy |
| `scripts/test_model_pmml.py` | 37 | 20.5 min with scikit-learn; 6.0 min with `OMP_NUM_THREADS=1` (`BoostedExportTest.test_same_probabilities` 189 s → 1.1 s) | **no** (`just model-test`) | without scikit-learn: **1 error** ✓, because `TrainerDepthKnobsTest` has no skip gate (`ModuleNotFoundError: joblib`, [results](results/model_pmml_without_sklearn.log)) |
| `examples/mini_db/run.sh` | accuracy check | not run | **no** | protal |

The slowest tests ([ctest](results/ctest_j4.tsv), [one process](results/gtest_one_process.tsv), [peak memory](results/peak_memory.tsv), [e2e](results/e2e_durations.txt)):

| Test | Time (ctest) | Peak memory |
|---|---|---|
| `PackedIndex.ABuildIntoThePackedLayoutWritesTheSameIndex` | 53 s | **13.2 GB** |
| `PackedIndex.TaxidAndGeneAsOneNumberComeBackForEveryGeneCount` | 22 s | 3.1 GB |
| `PackedIndex.TheColumnFormatPacksChunkByChunkToTheSameValues` | 15 s | 6.3 GB |
| `PackedIndex.AValueOutsideTheLayoutStopsTheColumnLoad` | 11 s | 3.1 GB |
| `IndexCodec.ColumnsCompressBetterThanRawCells` (zstd level 19, a ratio) | 6 s | 0.1 GB |
| e2e `BuildAmbiguousBasesTest` (5 builds with `--no_compress`) | 22 s | |
| e2e `ParallelIndexBuildTest` (three raw indexes hashed) | 17 s | |

**Why the 13.2 GB test is so big.** It saves two raw indexes of about 3 GB each to `/tmp` (`test_PackedIndex.cpp:278-283`). It then holds both as strings to compare them, while both `Seedmap`s are still alive. A GitHub `ubuntu-24.04` runner has 16 GB of memory. The CI sanitizer job runs the same test as a Debug build with AddressSanitizer.

## 2. Assertions that are vacuous or wrong (fix first)

### C++ unit tests

| Where | Problem |
|---|---|
| `test_SampleContext.cpp:524-533` ✓ | The `own` and `close` profiles never call `SetSingletonCongener`, and the default is 0 (`SampleContext.h:188`). So `EXPECT_FALSE(Vetoed())` holds whatever the rule does. The exemption it claims to test, a read that fits its own reference, is untested. So is the identity branch (`BaseIdentity < 0.95`, `Profiler.h:2801`). |
| `test_CommunityDesign.cpp:98` ✓ | `EXPECT_LE(uniform_pairs, 15)` sits at the mean. 20 species drawn from 200 genera of 5 give 14.67 genera with two or more species over 20 samples, and 52% of samples have one. The test passes because of its seeds; one more random draw in the designer can break it. The comment "hardly ever has two of a genus" is wrong. Bound it at about 30, or assert grouped > 2 × uniform. |
| `test_Zstd.cpp:185` ✓ | `EXPECT_NE(read, data + "garbage")` cannot fail, because the reader never produces the garbage. Assert that `read` is a prefix of `data`. |
| `test_ReaderAndTimers.cpp:412-428`, `:304-317` ✓ | If the reader stops early, the writer thread blocks on the full FIFO. The stream is destroyed only after `writer.join()`, so CI hangs instead of reporting a failure. Close the stream before `join()`. |
| `test_ReaderAndTimers.cpp:545` | In the `SeqReader` branch `batch_global.Success()` is vacuous, and the block readers' `Success()` is never checked. |
| `test_StrainOutput.cpp:408-409` | Asserts only "some IUPAC code" and "some plain base", so an `N` would pass. Pin `R` and `G`. |
| `test_StrainOutput.cpp:143-144` | Tests `CalculateCoverageVector()`, which has no production caller and returns `CalculateCoverageVector2`; it cannot fail. |
| `test_SNPCalling.cpp:270` | Asserts that two uint16 counters overflow alike: a pin of an overflow, not of behaviour. |
| `test_Parsing.cpp:680-733` | `NormalizeCigarReference` is a verbatim copy of `NormalizeCigar`'s slow path (`SamHandler.h:316-352`). For every CIGAR off the fast path the test compares the code with itself. |
| `test_Parsing.cpp:789-796` | Pins `std::stoul`/`stoi` exception types that no caller catches. |
| Tautological lines | `test_SharedSeeds.cpp:62` (its argument is `const&`), `test_AnchoredAlignment.cpp:526-527` (const strings compared with themselves), `test_AlignmentScreen.cpp:373`, `test_WFA2Wrapper.cpp:88/:90` (a value read back just after it was set), `test_PackedSequence.cpp:84-92` (checks the test's own helper), `test_SampleContext.cpp:196-198` (a `std::map` passed twice unchanged). |
| `test_GeneTableFile.cpp:252` ✓ | Pins an unordered map's iteration order and `bucket_count`, which, by its own comment, no output depends on. |
| `test_WFA2Wrapper.cpp:70-71` | Pins a WFA2-lib bug: the test fails once a library upgrade fixes it. Skip when the bug does not reproduce. |

### End-to-end tests

| Where | Problem |
|---|---|
| `:331` ✓ | The only truth comparison in the suite. Its regex `TP \d+, FP \d+, FN \d+` never checks the numbers. |
| `:479` | `MSAKnob.test_msa_samples_mirror_the_profiles` passes with no species reported. Both `MSAKnob` tests write to one `out` folder, so the second depends on the first. |
| `:880-909`, `:965` | `DepthKnobs` and `FalseCalls` compare profiles at knob 0 with the default. Every candidate in `sa` is a true species (`:836`), so the profiles are the same either way; only log lines show that the knob was applied. |
| `:241`, `:285`, `:550` | Pass on an empty glob, or only check that files exist. |
| `:2120` | Skips the whole format-identity test, gzip part included, without the zstd CLI, although `sam_text()` could read the file. |

### Python tests

| Where | Problem |
|---|---|
| `test_mini_db.py:3280` ✓ | `test_g_profiled_as_simulated` reads `test_f`'s output: it fails when run alone. |
| `test_mini_db.py:2855` ✓ | The `GtdbBuildTest` gate (section 1). |
| `test_model_pmml.py:378` ✓ | `TrainerDepthKnobsTest` has no skip gate (section 1). |
| `test_mini_db.py:2419` | `TraceRelativesTest` is gated on `$PROTAL`, `$SIMULATE` and art_illumina, none of which it uses, so it is skipped in CI. |

## 3. Tests of code nothing in production uses

| Test | Code | Status |
|---|---|---|
| `Compressor.*` ×3 (`test_Operations.cpp:72, :90, :106`) ✓ | `Utilities/Compressor.h` | Included nowhere in `src/`. `bgzf::CompressFile` is called only from it. |
| `Compressor.TheStreamingWriterGivesTheSameBytes` (`:118`) | `bgzf::Writer`, production code | Its oracle is the dead `CompressFile`. Compare piece sizes with each other and round-trip through zlib-ng instead. The `$(touch INJECTED)` and mtime checks are left over from the pigz era. |
| `test_SamRoundTrip.cpp:476-478` | `AlignmentEdits(std::string)` | No production caller. |
| `test_SamRoundTrip.cpp:316` | `HoldsPairedReads` | Only tests call it; production uses `ReadsOfSam`. |
| `test_SamFile.cpp:447` | `SamStreamSink` | No production user. |
| `ModelFeatures.TheModelMustFitProtal` | `ModelContractProblem(model, path)` | Production calls `ModelContractProblemInXml`. |
| `test_ReaderAndTimers.cpp:378-380` | `igzstream`'s `read_failed` | Tests the third-party library itself. |
| `zstd::MemoryReader` | | Exists only for a test helper. |
| `GenomeLoaderPacked.AGeneTakesThirtyTwoBytes` | | Repeats a `static_assert` (`GenomeLoader.h:492`). |

`scripts/index_reference.py` is referenced nowhere ✓. A reviewer also called `gradient_boosted_cmdline.py` and `hist_gradient_boosted_cmdline.py` dead, but `docs/databases.md` documents them.

## 4. Redundancy

### True duplicates in the C++ unit tests (merge or remove)

| Test | Duplicates | Action |
|---|---|---|
| `ZstdSeekable.GenomeLoaderReadsFramesInParallel` (`test_Zstd.cpp:507`) | `GenomeLoaderPacked.GenesCutByFrames…` (`test_PackedSequence.cpp:426`, a superset) | add its lowercase base to `:426`, remove |
| `test_Zstd.cpp:561` (writes seekable files with the default frame size) | `:507` | keep only its death test |
| `IndexHeader.FeaturesRoundTrip…` (`test_Index.cpp:92`) | `IndexHeader.ReferenceFingerprintRoundTrips` (`test_InputValidation.cpp:242`) | one test in `test_Index.cpp` |
| `Syncmers.ScanReusesItsBuffersAcrossLengths` (`:98`) | `:59` reuses one handler across lengths | remove |
| `ReverseComplement.TwiceGivesTheReadBack` (`:38`) | follows from `:19` and `:27` | remove |
| `AlignmentScreen.PackedWindowsGiveTheDecodedWindowsAnswer` (`:144`) | `:217` (uncommitted) | move its k = 9-10 and IUPAC cases into `:217` |
| `ProfileSam` `test_ProfileThreads.cpp:373` | `:389` | one test looping over `with_neighbours` |
| `GeneWindows` `:85` | `:153` (`:85` only catches uninitialised memory under ASan) | run `:85`'s cases inside the fill loop |
| `ThreadedGzStream` `:77`; `:383`; `:89` | `:253`; `:394`; `:560` and `:604` | remove |
| `UniqueKmers` `test_InputValidation.cpp:288-289` | `:531-534` | merge |
| `ReadExcess` cases in `test_Parsing.cpp:590-597` and `test_FalsePositiveFeatures.cpp:169-172` | `DivergenceBeyondTheBaseQualities` | keep them in one place |
| `ReadEvidence` `:142` (`FailedCandidatesLine`) | `SingleOutputHandler.CountsUnalignedReads…` (`test_SamRoundTrip.cpp:434`) | keep the handler half in one, the reader half in the other |
| `ReadEvidence` `:282` (`FeatureString`) | `ModelFeatures.NamesAreUnique…` (`test_StrainOutput.cpp:826`) | merge |
| `PhredScore` `:27` and `:34` | each other | merge |
| `IndexCodec` `:190` | re-implements the `RoundTrip` helper | fold into the helper |
| `PackedIndex` `:128-129` (exact layout numbers) | `:379-381` and the minimality check at `:398-403` | remove the pins |

### Organisation

- **Suites split over files ✓:** `IndexHeader`, `UniqueKmers` (test_Index / test_InputValidation), `SeqReader` (test_LongReads / test_ReaderAndTimers), `SamReader` (3 files), `ExtractVariants`, `ModelFeatures` (test_StrainOutput / test_InputValidation), `ProfileSam`, `MicrobialProfile`.
- **Tests in a file named for something else:**
  - 7 `GeneConservation` tests in `test_Database.cpp`.
  - `MateGuidance` in `test_SamRoundTrip.cpp`; `AcrossGenes.MateGuidanceLooksPastTheGenesEnd` is in `test_GeneNeighbours.cpp`.
  - A `SamReader` test in `test_GeneNeighbours.cpp:357`.
  - Two simulator tests in `test_ReadEvidence.cpp:257, :269`.
  - Feature tests in `test_Parsing.cpp`.
  - `ModelFeatures` and `Abundance` (about 350 lines) in `test_StrainOutput.cpp`.
  - `test_Operations.cpp` is named after a unit that no longer exists.
- **Copied helpers:**
  - `ScratchDir`/`TempDir` in 13 files.
  - `TinyReference` 4 times.
  - `Reference`/`Record`/`Profile`/`RandomSequence`/`Mutated` in 3-6 files.
  - `RandomReference`/`ExactRuns`/`Derive` 3 times; `TestData` twice; `SmallIndex` twice.

  One `tests/TestUtil.h` would remove about 150-300 lines. It would also fix the per-process directory of `test_StrainOutput.cpp`'s `TinyReference`, which two live objects share (`:208, :230, :266, :292`).

### End-to-end tests

**Within the file:**
- The same paired-end sample with default options is aligned about 12 times (GeneConservation, DepthKnobs, FalseCalls, BadInput ×2, SamInput, ReadTypeModel, ModelContract). Tests that change only profiling could use `--profile_only` on one shared SAM.
- Three FASTA-quality tests (`:2330, :2602, :2701`), three mixed-map tests (`:2314, :2617, :2715`) and three `--profile_only` read-type tests (`:2294, :2546, :2689`) could each be one parameterised test.
- `Qcmsa.test_no_filtered_msa_is_reported` (`:2201`) is contained in `:2219`.
- `Failure.test_truncated_index` (`:1357`) ≈ `CompressedDatabase.test_corrupt_seekable_index` (`:1665`).
- `BuildAmbiguousBasesTest` (22 s) runs five builds where one dirty build with all four codes would do. The unit test `AmbiguousKmers.*` already checks the window arithmetic.

**Against the unit tests:**
- `BadInput`, `SamInput`'s broken SAMs, `FailFast`, the corrupt databases and the malformed knob files repeat unit tests at the command line, often with the same messages.
- That is the right layer for exit codes and for isolating samples from each other. Keep them, but they need only check the exit code and the isolation, leaving the message text to the unit tests.

### Python tests

- **The same rules tested twice.** `test_model_pmml.py` and the C++ tests (`test_SampleContext.cpp:375-440`, `test_StrainOutput.cpp:767-819`) test the prior adjustment, the false-call calibration and the knob curve, each on hand-copied vectors. Nothing checks that Python and C++ give the same numbers. One golden-vector file read by both would.
- **Duplicated tests.** `TraceRelativesTest` (`test_mini_db.py`) and `scripts/test_trace_relatives.py` test the same script: merge them into the latter. Two classes are named `GeneNeighboursTest` (mini_db and e2e), at different layers: rename one.
- **Inside `test_mini_db.py` (3,431 lines):**
  - The `reference.map` offset check is copied 3 times, and the fake mirror and `datasets` stub 4 times.
  - `GtdbBuildTest` makes about 9 full pipeline runs, has four "a rerun does nothing" checks, and its `test_e` repeats `test_a` and `test_d`.
  - About 23 of `MiniDbTest`'s 42 tests do not use its fixture.
  - About 30 regular expressions on console text include step counters (`5/9`), which have broken before.
- **`examples/mini_db/run.sh` is not redundant.** It is the only check of accuracy against the truth.

## 5. What no automated test checks

1. **The end-to-end suite and the accuracy example.** Both binaries are built in CI and never run, as the round-5 audit noted (B4, [2026-09-30](2026-09-30-round5-audit/README.md)).
2. **Accuracy.** No automated test checks TP/FP/FN counts, abundances against the truth, or that reads of nothing in the database give an empty profile.
3. **The default model.**
   - No CI test loads a gradient-boosted `modelChain` PMML.
   - The mini database copies `scripts/random_forest.xml` as `model_pe.xml` ✓ (`gtdb_to_protal_db.py:638`), and its se/pb/ont models are that file again.
4. **The feature-name contract between protal and the trainer.**
   - e2e `:819` ✓ checks only the normalized and adjacency groups. The distance, depth, divergence, unfiltered and ref names, all in the default set, are not checked until a training fails.
   - The constants behind default features (`kPriorSpecies`, `kSpillAtZero`, `kSpillDecade`) are tested only through the same symbols, so a change passes every test and invalidates the shipped models.
5. **Thread identity of every default feature.**
   - `test_ProfileThreads.cpp`'s taxa have 4 genes, fewer than `kMinSharedGenes` (10). They also have no conservation table, so the distance and conservation features are trivial there.
   - The thread loops in `test_FalsePositiveFeatures.cpp` and `test_ReadEvidence.cpp` profile a few-kB SAM in one 1 MB chunk.
6. **The default (non-AVX2) copies of the `target_clones` functions.** CI runners have AVX2, so `AnchoredAligner::Align`, `AlignAnchor` and the others run only their x86-64-v3 copy. A build with `-DPROTAL_NO_CLONES` ✓ (`TargetClones.h`) would test the other.
7. **The alignment screen as production runs it.** That is `operator()` with the shared `ReadKmers` strands; the handler test calls `AlignAnchor` without it. Also untested: reverse-strand alignments through the handler.
8. **Error paths:**
   - `Bundle::Open`'s directory checks;
   - a full disk when writing a SAM (`/dev/full`);
   - `LoadFromMap`'s errors;
   - an edited `reference.map` of the same size, which is accepted (`GenomeLoader.h:920`);
   - `CoverageVector` with a read outside its range, which writes out of bounds;
   - `SortKey` at its maxima;
   - an interrupted run resumed.
9. **Code without a unit test:** `LongReadAligner` (end to end only), and the `Haplotypes` settings (`max_haplotypes`, `min_share`, `min_reads`).
10. **Scripts without a test:** `recurrent_calls.py`, `prevalence_calls.py`, `protal_profile_utils` (which the conda recipe ships) and `scripts/strain_test/*`.

## 6. Documentation that is out of date

- **`docs/development.md:17-20`** lists 13 areas for 34 test files: list all of them, or say "for example".
- **`docs/development.md:44-54`** says the end-to-end tests run against "any database". In fact they hard-code the mini database's taxa, taxids and genes: against the v0.7.3 database 13 of 102 failed ([2026-10-03](2026-10-03-performance-review/README.md)). They also need numpy and Linux, and cover PacBio, ONT and phasing besides paired-end and single-end reads.
- **`docs/development.md:56-58`**, and the justfile comment on `mini-db-test`, promise that `GtdbBuildTest` skips when it cannot run: it errors without `$PROTAL`.
- **`docs/claude/2026-10-01-f1-opportunities/scripts/test_errors.py`** is an analysis script, not a test. `pytest` would collect it by its name; rename it.

## 7. Recommendations, ranked

1. **Run what exists, and make missing prerequisites fail loudly.**
   - Add a CI job: `build_mini_db.sh` (20 s here), the end-to-end suite (143 s here), and `examples/mini_db/run.sh`.
   - Add a switch (`PROTAL_E2E_REQUIRED=1`) that turns the module-level skips into errors in CI. Move `Version` and `Simulator` out of the database gate.
   - Fix the gates of `GtdbBuildTest`, `TrainerDepthKnobsTest` and `TraceRelativesTest`.
   - Set `OMP_NUM_THREADS=1` in `test_model_pmml.py` and drop its second `--features auto` training (`:840`). Then add `python3-sklearn python3-pandas python3-joblib` to CI and run it.
2. **Fix the assertions of section 2,** starting with the `Vetoed` checks, the community-design bound, the hanging pipe tests and the e2e truth line (assert `TP 3, FP 0, FN 0`).
3. **Test accuracy and the default model.**
   - Give the mini database a small committed GBM PMML, with sklearn's probabilities for it, scored by protal end to end.
   - Add a sample of foreign reads that must give an empty profile.
   - Check every default feature group's names at e2e `:819`.
   - Share golden vectors between `test_model_pmml.py` and the C++ tests.
4. **Cut the cost.**
   - In `PackedIndex.ABuildIntoThePackedLayoutWritesTheSameIndex`, compare the raw files in chunks or by hash, and free `wide` once its files are written. That removes the two ~3 GB strings, and the second index for part of the test (not measured).
   - Gate `IndexCodec.ColumnsCompressBetterThanRawCells` like `FlexScan.Bench`. Move both benchmarks (`FlexScan.Bench`, `GeneTables.BenchLoadOfLargeTables`) out of the unit binary.
   - In the end-to-end suite, run one shared baseline and use `--profile_only`, and run one ambiguous-base build.
5. **Remove the tests of dead code, and the dead code** (`Compressor.h`, `bgzf::CompressFile`, `AlignmentEdits(std::string)`, `HoldsPairedReads`, `CalculateCoverageVector()`, `index_reference.py`). Merge or remove the duplicates of section 4.
6. **Reorganise.** Add a shared `tests/TestUtil.h` and regroup the split suites (`GeneConservation` gets its own file). Rename `test_Operations.cpp`. Split `test_mini_db.py` into download, collector and scenarios, gene neighbours, and GTDB build. Replace console regexes with checks on files.
7. **Close the gaps of section 5,** starting with a `-DPROTAL_NO_CLONES` CI build, a thread test with ≥10 genes per taxon and a conservation table, and the reverse-strand and `operator()` screen tests.
8. **Fix the documentation of section 6.**

**When the website would change:** none of this changes how protal runs, so the website is not affected.

## Follow-up: implementation started, stopped unfinished (2026-10-06)

Six agents started on recommendations 1-8, one area each, on branches from `a20bd2e`. The work was
stopped before any of it was built or run, because of the session's token budget. Each agent's partial
edits are saved as one unverified WIP commit on its own branch. None of it is on `audit-fixes`.

| Branch | Commit | Area |
|---|---|---|
| `worktree-agent-a9f8eac5e10212542` | `e27ce02` | C++: I/O, compression, database, packed index (`Compressor.h` removed) |
| `worktree-agent-ab680925d62414681` | `b9f3f6a` | C++: parsing, input validation, strain output, SNPs, haplotypes, exact sums |
| `worktree-agent-a1ed36a11b648c8f2` | `a230218` | C++: packing, index, syncmers, seeds, screen, anchored alignment, WFA2, long reads |
| `worktree-agent-af4d1a274c207ebe1` | `dba2066` | C++: gene neighbours, windows, profile threads, sample context, FP features, read evidence |
| `worktree-agent-a662f5866a84af2fd` | `c7eacd0` | end-to-end suite (with a `tests/e2e/data/` fixture) |
| `worktree-agent-abfb50c7edf922063` | `55bf6a2` | Python tests: `test_mini_db.py` split into five files, gates, `prerequisites.py` |

To continue:
1. Build and run each branch.
2. Finish its items from the lists above.
3. Merge the branches. They touch disjoint files, apart from `tests/CMakeLists.txt`.
4. Do the reorganisation pass (shared test helpers, regrouped suites).
5. Change CI and the docs.

The Python branch splits `test_mini_db.py`. Since `a20bd2e`, `652ed53` has changed that file, so the
split has to be redone on top of `652ed53`.

**Done on 2026-10-07:** the WIP branches were merged, built, run and finished, and recommendations 1-8 are on
`audit-fixes` (`049f7a1`, `ac32e2d`, `7086c09`, `51d111d`). See
[2026-10-07-test-suite-fixes](../2026-10-07-test-suite-fixes/README.md) for what changed, what the suites cost
now, and what is left.
