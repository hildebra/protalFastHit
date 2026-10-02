# Prefetching the k-mer lookups

- **Date**: 2026-10-02.
- **Code**: `audit-fixes` at `b7322c1`, plus this change (`scripts/prefetch.patch`), committed as `87b857e`. The output checks and timings ran on the change before its last step, which only puts the new lookups back under the seeding timers (below).
- **Machine**: WSL2 Ubuntu 24.04 on an Intel Core Ultra 7 258V (4 fast and 4 low-power cores, 6 vCPUs), GCC 13.3.
  - No other sessions were running.
  - The laptop still changed speed from run to run, by up to 2×, and later in the day ran about 1.8× slower throughout (throttling, or the low-power cores).
- **Data**: the v0.7.1 benchmark world's database (`~/bench071/V071`). `rl150_p500000_s_1` gives 500k pairs; `rl150_p5000000_s_1` gives 5M.
- **Question**: round 3 (`../2026-10-02-performance-round3/README.md`) found seeding at 19% of short-read time for ~6% of the instructions, so it waits on memory. Can prefetching the lookups help, across reads in particular?

## Summary

1. **Why seeding waits.** A lookup reads its k-mer's 24-byte control block in the key map (3.2 GB), then the values that block points to. Both miss the caches. `Seedmap::Get` branches on what it loaded (an empty block or key returns), and a mispredicted branch on a cache miss throws away the following lookups' loads. So a read's lookups waited one after another, not together.
2. **What pays, measured in one process** (`scripts/pfbench.cpp`; the real index and the k-mers of 400,000 real reads; variants alternated per batch of 1,000 reads; thread CPU time):

   | variant | seeding time |
   |---|---|
   | as protal | 1.000 |
   | the read's own control blocks prefetched first (round 2's prototype) | 0.865–0.870 |
   | the control blocks prefetched two reads ahead | 0.847–0.881 |
   | pipelined within the read: the block 8 k-mers ahead, the values 16 lookups ahead | **0.758** |
   | two stages: blocks ahead, then the next read's lookups made and their values prefetched | **0.688–0.696** |

   The values a block points to are the larger half. Prefetching them needs the block first, hence the two stages.
3. **In protal** (both of the above):
   - `ChainAnchorFinder::FindSeeds` pipelines within the read, which covers every read type, long reads included.
   - The paired-end and single-end loops seed one read ahead. They prefetch the next pair's control blocks while a pair is aligned (`PrefetchKeys`). At the start of the next pair they make its lookups and prefetch their values (`PrepareLookups`), while the pair after it is read and its k-mers taken.
   - The lookups, their order and the reads' order stay as they were. SAM text is identical at `-t 1` (paired-end, single-end, Nanopore), and the sorted SAM and the profile at `-t 6`.
4. **Whole runs: about 4–5% less alignment time, within this laptop's noise.**
   - 5M pairs at six threads: fastest of 6 alternated runs 35.06 → 33.76 s (−3.7%), median 35.98 → 34.32 s (−4.6%), CPU time median −4.6%. In the first, quiet round of `pf_check.sh`, though, the change was 2.7% slower (19.16 → 19.68 s), one run each.
   - 500k pairs at one thread: −5% in that quiet round (9.04 → 8.58 s), −4.5% fastest of 8 later, the medians equal (17.08 and 17.18 s).
   - Single-end reads at one thread: −4.5% in the quiet round (3.65 → 3.48 s).
   - That is about what the benchmark implies: seeding is 15–19% of the time and got ~30% cheaper. A quiet machine (or `perf stat` cycles on the HPC) would pin it down better than these runs.
5. **At GTDB size the gain should be larger.** The values array grows to ~35 GB (round 2's memory report), so more of the value reads miss, and their prefetch is the part that pays most here.

## Measurements

**pfbench** (`results/pfbench.txt`): 22.9 k-mers per read, about 21.5 seeds. Two runs on different days' states of the machine gave the same ratios within 0.02. Each variant makes the same lookups and seeds for a read. The seeds per read in the table differ slightly only because each variant gets other batches of reads.

**Whole runs**:

- `results/pf_check.txt`: the output checks, and five alternated rounds of 500k pairs (`-t 1`), 500k single-end reads (`-t 1`) and 5M pairs (`-t 6`). After round 1 the run times rose 2× in both builds alike, so only round 1 compares.
- `results/fastest_of.txt`, `runs_fastest_of.tsv`: 8 alternated runs of 500k pairs and 6 of 5M pairs, later, with the machine slower throughout.
- Stage timers. In the change as timed, the "Seeding" timer dropped by two thirds: `PrepareLookups` ran outside it, so the lookups were no longer counted. In the commit, `PrepareLookups` counts as part of the read's "Seeding", "Seed-finding operator" and "Seed- and Anchor-finding" (`Benchmark::Start(false)`: timed, not counted as a call of its own), so `misc/<prefix>_runtime.tsv` compares with earlier runs.
  - `results/stage_timers.txt`, 500k pairs at `-t 1`, three alternated runs each: Seeding 1.53–1.76 s with the change, 2.23–3.17 s without. That is about what pfbench found.
- `results/suite.txt`: the committed change passes the unit tests (295) and the end-to-end tests (128).

## The change

- `Seedmap::PrefetchKey`: the 32 bytes `Get` reads of a key's control block and the next one (one or two cache lines).
- `KmerLookupSM::PrefetchKey` and `PrefetchValues`: the first line of a lookup's values and of its flex cells.
- `ChainAnchorFinder`:
  - `LookupKmers` prefetches the block of the k-mer `kPrefetchKmers` (8) ahead.
  - `FindSeeds` prefetches the values of the lookup `kPrefetchLookups` (16) ahead, and takes the lookups `PrepareLookups` made for this read instead of making them.
  - `PrefetchKeys` and `PrepareLookups` are for the loops. A prepared read must be the next one the finder seeds; the loops call the finder on every read and mate, and `operator()` always seeds.
- `classify::RunPairedEnd` and `RunSingleEnd` read one pair (read) ahead into a second set of records and k-mer lists and swap them in after each pair. The pairs are aligned in the order read, so the output order does not change.
- Long reads take only the in-read pipelining. A long read's thousands of k-mers give the prefetches lead time by themselves.

## How it was run

```bash
bash scripts/pfbench_run.sh     # pfbench against the head07c build's libraries, pinned (ROUNDS, READS, CPU)
bash scripts/pf_check.sh        # both builds; outputs compared; five alternated rounds of whole runs
bash scripts/pf_time2.sh        # fastest of 8 and 6 alternated runs
bash scripts/tcheck.sh          # the stage timers of the committed change against HEAD's build
```

The scripts use the session's scratch paths (`build_wt.sh` and `check.sh` built a clean `git archive` plus the patch in `~/mt-work/<name>`); adjust them to rerun.
