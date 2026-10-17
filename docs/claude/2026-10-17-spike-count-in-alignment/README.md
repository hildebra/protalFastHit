# The spike check from the records as they are written, not from a second read of the SAM

Date: 2026-10-17. Branch `audit-fixes`, on `1e60c1a` (the working tree; not committed).

- **Why:** the speed report of 0.8.5 ([2026-10-17-perf-0.8.5.md](../2026-10-17-perf-0.8.5.md), section 4) found that
  the spike-gene check (`--spike_genes auto`, new in 0.8.5) reads every sample's SAM a second time before profiling
  it.
  - At GTDB r226 that cost 0.72 s of a paired-end sample's 5.1 s of profiling, 3.7% of a gated pe run, on every
    sample, spiky or not.
  - For a HiFi sample it cost 0.12 s, for an ONT sample 0.15 s.
- **The task:** count the records per gene as the alignment writes them, and keep the SAM pass only for a SAM the run
  did not write.
- **Checked on:**
  - the unit and e2e suites (WSL, cores 0-3), with a new e2e test of a spiky sample;
  - two real mock samples on the reduced r226 database (`~/Mocks`), against protal at `1e60c1a` (before the change).

## What changed

| file | change |
|---|---|
| `src/IO/SamFile.h` | `SamOutput::SetObserver`: a function that sees each block of records `Write` takes, on the writing thread, outside the lock. |
| `src/Profiling/SpikeGenes.h` | `Tally` (a block's records counted by `CountRecords` on its thread, then added under a lock); `Found` (a check and its spike genes); `AlignedChecks` (the checks of the SAMs the run wrote, by path, each taken once). |
| `src/Profiling/SpikeCheck.h` (new) | `CountedRecords`, `FindSpikeGenes` (moved out of `MicrobialProfile`, which used only the gene tables) and `Check`, so that the alignment can make the check without `Profiler.h`. |
| `src/Run/Align.cpp` | In a run that profiles with the check on, the sample's `SamOutput` gets a `Tally` as its observer. Once the SAM is complete and renamed, its check goes to `AlignedChecks`, before `sample_done` (`--profile_ahead` may profile the sample at once). |
| `src/Run/Profile.cpp`, `Profiler.h` | `ProfileSample` takes the SAM's check if this run wrote it (`Profiler::UseSpikeCheck`); `ProfileSam` uses it, else counts the SAM as before (`CountGeneRecords`, then `Check`). |

- **Why the result is the same:**
  - The tally applies the same parser (`spike_genes::CountRecords`) to the same text that goes into the file, block by
    block. The blocks are whole reads' records.
  - Nothing changes the records after `Write`; `Finish` adds only the header, whose lines the parser skips.
  - So the counts are the counts of the finished SAM, and the same function (`Check`) turns them into the check and
    the spike genes.
- **When the SAM pass still runs:**
  - `--profile_only`;
  - a SAM that already existed (the alignment skips the sample);
  - any SAM the run did not write. There is no entry for its path, and protal counts it from the file as before.
- **Memory:** only the check is kept between the alignment and the profiling: the sample check's numbers and its spike
  genes (usually none). The counts themselves (about 77,000 genes for an r226 pe sample) are freed when the SAM is
  done.
- **Outputs:** the log, profile and spike files are unchanged. The runtime row "Profiling: spike genes" is now ~0 s for
  a sample whose SAM the run wrote. The check itself runs at the end of the sample's alignment, outside the stage
  timers, and takes milliseconds.

## Tests

- **Unit** (641 pass):
  - `SamFile.AnObserverSeesTheRecordsTheFileHolds`: SAMs written from 4 threads, with records of flags 0, 4, 256,
    2048 and 99:
    - in every format (plain, BGZF, zstd, zstd with header room, the header first);
    - the observer saw exactly the bytes after the header;
    - the tally equals `CountRecords` over the file as `SamInput` reads it back;
    - an empty `Write` (genes only) shows the observer nothing.
  - `SpikeGenes.TallyOfBlocksFromThreadsCountsAsTheWholeText`: 4 threads adding blocks give `CountRecords` of all the
    text, and `Take` empties the tally.
  - `SpikeGenes.AlignedChecksAreTakenOnce`.
- **e2e** (195 pass): `SpikeCheckTest` puts a pile of 150 pairs on one gene of each of two species of the mini
  database, plus one pair on each of 12 of their other genes and an evenly covered third species. The sample is in
  spike mode with one spike gene per piled species. Three runs of it give the same spike line (but for the file's
  path), `.profile` and `.profile.spike_genes.tsv`:
  - aligned and profiled in one run: spike timer 0.000001 s, the check handed over;
  - `--profile_only` of that SAM: spike timer 0.013 s, the SAM pass;
  - `--profile_ahead` with a second sample after it, so that the spiky sample is profiled while the next one aligns.

## Before and after on real mock samples

`scripts/time_spikes.sh`:
- the reduced r226 database (`~/Mocks/reduced_db_r226`);
- two D6331 paired-end samples: Dohlman et al. SRR38158090, the spiky one, and HKU SRR23926922;
- 4 threads on cores 0-3;
- 3 rounds, the two binaries' order alternating.

Seconds:

| sample | build | spike genes (runtime row) | profiling of the sample | aligning |
|---|---|---|---|---|
| Dohlman | `1e60c1a` | 0.056 / 0.082 / 0.065 | 0.147 / 0.200 / 0.167 | 1.02 / 1.40 / 1.17 |
| | working tree | 0.0000 (all three) | 0.102 / 0.108 / 0.105 | 1.09 / 1.18 / 1.16 |
| HKU | `1e60c1a` | 0.055 / 0.053 / 0.053 | 0.311 / 0.289 / 0.283 | 1.77 / 1.59 / 1.61 |
| | working tree | 0.0000 (all three) | 0.262 / 0.303 / 0.236 | 1.84 / 1.68 / 1.61 |

- **Identical outputs:**
  - the SAM records (sorted; 114,393 and 91,044);
  - `.profile`, `.profile.composition` and `.profile.spike_genes.tsv`;
  - the spike line. Dohlman is in spike mode: 8 of its 10 most abundant taxa spiked, 78.9% of its reads, 20 spike
    genes, 80,219 records left out.
- **`--profile_only` of the new SAMs** (the SAM pass) gives the same spike line and files, in 0.076 s (Dohlman) and
  0.052 s (HKU).
- **The profiling of the sample is shorter by the spike pass**, 0.04-0.09 s on these small samples. Aligning varies
  more between rounds than the tally could cost (below).

**What the tally costs the alignment threads** (`scripts/tally_bench.cpp`: the samples' SAM text in 1 MB blocks through
`Tally::Add` on one thread, `-O3 -march=x86-64-v3`, best of 5):
- HKU: 127 ms per GB of SAM text, 51 ns a record;
- Dohlman: 84 ms per GB, 35 ns a record.

At r226 pe1 (1.73 GB of SAM text, 2.22 M records), that is ~0.1-0.2 s of CPU spread over 32 alignment threads, ~5 ms of
wall time, against the 0.72 s pass it replaces. The lock is held only to add a block's few hundred genes.

**Expected on the cluster** (`measure_performance.sh`, as in Perf0.8.5): "Profiling: spike genes" from 0.72 to ~0 s for
pe, from 0.12 and 0.15 s for HiFi and ONT. Gated pe at C 4: from ~19.3 to ~18.6 s.

## Files

- [scripts/](scripts/) (WSL, cores 0-3):
  - `build_test.sh`: the working tree's build, unit tests, mini database and e2e suite;
  - `build_ref.sh`: protal at `1e60c1a` from `git archive`;
  - `time_spikes.sh`: the before/after runs and the comparison;
  - `check_bench.sh` and `tally_bench.cpp`: the sorted SAM comparison and the tally's cost;
  - `e2e_spike.sh`: `SpikeCheckTest` alone, verbose, with its outputs.
- Run outputs in WSL `~/spikecount`, removed after the report.
