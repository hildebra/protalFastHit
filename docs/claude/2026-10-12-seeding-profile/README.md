# Profiling the seed and anchor finding of short reads at GTDB r226 scale

- **Date**: 2026-10-12.
- **Code**: `714649b` (0.8.2 and the performance audit's report), unchanged. Built from `git archive` in WSL, Release
  with `-g` (the same code, with debug info for callgrind's source lines).
- **Question**: at r226 on 8 threads a paired-end run is now 80% aligning, and two thirds of that is "Seed- and
  Anchor-finding" (the 16th cluster run, [performance audit section 11](../2026-10-11-performance-audit/README.md)).
  Where does that time go, and what would cut it?
- **Data**:
  - r226: the 16th cluster run's logs and stage timers (0.8.2 on the r226 v24 database, 8 threads, the real gut sample
    B10 of 50.6M pairs; `../2026-10-11-performance-audit/results/cluster_v16_082/`).
  - Local: every GTDB r226 archaeal species (6,968 species, the marker genes of 17,245 genomes), converted by HEAD's
    `build_gtdb_database.py` and built with the pipeline's finished-database arguments ([sa_setup.sh](scripts/sa_setup.sh));
    1M simulated pairs of 150 bp from 120 species with strains and congeners, and their first 100,000.
- **Machine**: WSL2 on the Core Ultra 7 258V (AVX2, no AVX-512, no `perf`), jobs pinned to cores 0-3. Other sessions
  ran builds and tests meanwhile (load 4-7); the instruction counts do not mind, the local wall times do.
- **Status**: the local part and the cluster job are done ([results](#results-at-r226-the-cluster-job-2026-10-12));
  whole blocks are prefetched since ([implemented](#implemented-whole-blocks-prefetched-2026-10-12)); the give-up rule
  waits for the stage-1 probe's counts.

## Summary

1. **What the 62 s per thread are** (r226, pe, 8 threads): the k-mer lookups ("Seeding") 44.2 s (71%), extending the
   anchors 8.3 s (13%), sorting the seeds 4.3 s (7%), pairing 2.3 s, recovering anchors from the mate 2.1 s, sorting the
   anchors 0.6 s. Taking the k-mers ("Retrieve k-mers", 6.0 s) is timed apart.
2. **A lookup costs ~175 ns of a thread** (the benchmark at r226; the run's own timer gives 159 ns): 2.22G lookups,
   nearly every k-mer's core found in the index (99.2%), every block scanned (the stop at 128 seeds never comes), 70.4
   flex cells and 2.0 seeds per lookup. Split ([results](#results-at-r226-the-cluster-job-2026-10-12)): **the key map
   19.5 ns** (11%; prefetched a pair ahead, it is nearly free), **the flex cells 105 ns** (60%), the entries, the seeds
   and the sort 50 ns (29%). With warm caches the whole lookup takes 75 ns: **~100 ns of each lookup waits on
   memory**, mostly for the block's flex cells. On 32 threads a lookup costs 3% more: the limit is each thread's
   misses in flight, not the node's memory.
3. **Prefetching more pays a little, exactly**: every line of a block prefetched instead of the first (P2) −10%, the
   lookups made a pair earlier as well −11%, four pairs interleaved −5%; the same seeds. The misses are not late, they
   are many: ~10 lines a lookup, 22 lookups a mate.
4. **95% of the lookups are chance matches.** For 94.9% of the lookups that scan flex cells (85% of the cells scanned)
   no entry is within 5 bases of the read's 16 flank bases (best score ≤ 11); only 1-3% of the lookups into blocks
   under 256 cells match exactly (15% of those into larger blocks). Nearly every 15-base core is in the r226 index, and
   most of a metagenome's reads come from outside the marker genes: the sample's 101M mates gave 2.2M SAM records.
   Their seeds still pair: 79% of the mates get an anchor, and the k-mer screen refuses 90% of the 177M candidates.
   **Most of the seed and anchor finding (and of the screen) is spent on reads that will never align.**
5. **So**:
   - exact, now: prefetch whole blocks (capped), **−10% of the seeding**, ~4.4 s per thread of a pe run on 8 threads
     (~4% of it);
   - rejected: large blocks grouped by their first flex bases (to skip groups when the best is 15-16): 75% of the
     lookups there score ≤ 11;
   - **the large lever, not exact: give up on a mate whose lookups score low** (e.g. no best score above 11 among its
     first few). It would cut most of the seeding, the seed sort, the pairing and the screen's refusals. It first needs
     a count of how the lookups of mates that end up aligned score (none of them should be given up), then the outputs
     on the full sample and the r226 F1 benchmark.
6. **Exact, locally** (callgrind on all r226 archaea, 100,000 pairs, one thread; that world is not r226's: 43% of its
   cores in the index, 4 cells per block): **82 instructions to make a lookup (per k-mer), 226 to score a block and emit
   its seeds** (3.2 cells), 47 of them the AVX2 scan. The seed and anchor finding is 30% of the alignment loop's
   instructions, taking the k-mers 24%, and ~17 heap allocations per pair cost up to ~10% (`operator new` with its
   `malloc`, and `free`: a `resize` per read in the k-mer handler, the timers, the column-weights cache).

## Where the 62 s go at r226

From the 16th cluster run's `stages.tsv` (run 2, pe, seconds per thread of 8):

| stage | s per thread | what it does |
|---|---|---|
| **Seed- and Anchor-finding** | **62.2** | `ChainAnchorFinder::operator()` and the mate's anchor recovery |
| – Seeding | **44.2** | `PrepareLookups` (the lookups made, their values' first lines prefetched) and `FindSeeds` (the blocks scored, smallest first, until 128 seeds) |
| – Extending Anchors | 8.3 | each anchor's gene record and bases around it (prefetched 8 anchors ahead), exact matches extended |
| – Sorting Seeds | 4.3 | `SharedSeeds::Sort`: seeds counted by taxon and gene, those sharing one sorted as 128-bit keys |
| – Pairing | 2.3 | `FindPairs`: anchors from runs of seeds on one gene's diagonal |
| – Anchor recovery | 2.1 | the mate's best genes looked for among this mate's lookups |
| – Sorting Anchors | 0.6 | |
| Retrieve k-mers (apart) | 6.0 | closed syncmers of each mate (AVX2) |

The run's seeding line: 2,242,166,093 k-mers looked up, 99.2% in the index, 2,221,730,360 blocks scanned with 156.4G
flex cells (70.4 per block; blocks of 16-255 cells: 63% of blocks and 56% of cells; 256-4,095: 5.5% and 41%),
4.37G seeds (12.6% of them sharing their taxon and gene with another seed of the mate, the only ones that can make an
anchor), 211M anchors. Per thread: 278M lookups in 44.2 s, **159 ns each**. On 32 threads (the 13th run, 0.7.8) it
was 186 ns: the lookup's cost per thread hardly depends on how many threads share the node.

## How a lookup goes, and what it touches

Every closed syncmer of a mate is a 31-mer: a 15-base core and 8 flank bases on each side (`Seedmap`). Per lookup:

1. **The key map** (`Seedmap::Locate`): the core's 24-byte control block in the 3.2 GB key map gives the slots of its
   values. One cache line, sometimes two. `PrefetchKeys` fetches them while the previous pair aligns (~6 µs ahead).
2. **The flex cells** (`flex_scan::ScanAvx512`): one 32-bit cell per entry, the entry's 16 flank bases; each is scored
   by the bases it shares with the read's. 70 cells are 280 bytes, ~5 lines. `PrepareLookups` prefetches the first
   line just before the pair is seeded (~1.5 µs ahead: the next pair's reading and k-mers); the rest come by the
   hardware prefetcher, or by waiting.
3. **The entries of the best cells** (`Seedmap::EntryValue`): taxon, gene and position of each best cell, ~42 bits
   each, after the flex cells. Two per lookup on average, anywhere in the block's ~370 bytes of entries: ~1-2 lines.
   The first line of the entries is prefetched with the flex cells; the best ones are known only after the scan.

So a lookup reads ~8 lines from ~30 GB of index: ~1.1 TB over the run's seeding, ~26 GB/s on 8 threads. That is far
from the node's bandwidth but near what one thread gets from its outstanding misses (~10 lines in flight, ~100 ns
each: ~6 GB/s a thread, ~3.2 GB/s here). The arithmetic is ~400 instructions (below), ~60 ns at the run's IPC of 1.9.

## Exact, locally: instructions of the seed and anchor finding

Callgrind on 100,000 pairs, one thread, `--no_profile` ([sa_cg.sh](scripts/sa_cg.sh),
[callgrind_local_100k.txt](results/callgrind_local_100k.txt)). The alignment loop (`RunPairedEnd`) is 4.17 G
instructions (41.7k a pair):

| part | instructions | share of the loop | per unit |
|---|---|---|---|
| seed and anchor finding (`ChainAnchorFinder::operator()`) | 1.24 G | 30% | |
| – `FindSeeds` | 0.59 G | 14% | 335 per lookup (in the index) |
| – – `GetFromLookup` | 0.40 G | 10% | 226 per lookup (3.2 cells); the AVX2 scan 47 |
| – `SharedSeeds::Sort` | 0.14 G | 3% | |
| – `FindPairs`, `ExtendAnchor`, `RecoverAnchors` | 0.15 G | 4% | |
| making the lookups (`PrepareLookups`, `LookupKmers`) | 0.39 G | 9% | **82 per k-mer** |
| taking the k-mers (`SimpleKmerHandler`, AVX2) | 1.02 G | 24% | 5,100 per pair |
| the alignment handler | 0.55 G | 13% | |
| reading the pairs (`SeqReaderPE`) | 0.51 G | 12% | |
| reverse complements | 0.24 G | 6% | |
| heap: `operator new` (with its `malloc`) and `free` | ~0.45 G | up to ~10% | 1.7M allocations, ~17 a pair |

At r226 a lookup scans 70 cells with AVX-512 (~9 instructions per 16 cells) and emits 2 seeds: ~400 instructions with
the key map's 82, as estimated above. The cluster's warm-cache variant measures it.

## Why the local world is not r226 here

| | local (all r226 archaea) | r226 (16th run) |
|---|---|---|
| k-mers whose core is in the index | 43% | 99.2% |
| flex cells per block | 4.0 | 70.4 |
| blocks of 256-4,095 cells | 0.0% of blocks, 2.3% of cells | 5.5% of blocks, 41% of cells |
| the index in memory | 3.2 GB key map, 0.4 GB values | 3.2 GB key map, 26.6 GB values |

r226's 2.1G index entries cover nearly every 15-base core, so nearly every lookup reads values, and many cores are
shared by thousands of species. The local world has neither. Its seedbench runs ([seedbench_local_run1.tsv](results/seedbench_local_run1.tsv),
[run 2](results/seedbench_local_run2.tsv), the second at a load of 7) only show that the benchmark works.

## The benchmark (`seedbench`)

[seedbench.cpp](scripts/seedbench.cpp), built as an extra target in a copy of the source tree with protal's flags
([seedbench_build.sh](scripts/seedbench_build.sh)); not part of protal. It loads the database's gene tables and index as
protal does (the packed layout from the reference's maxima), takes the k-mers of the first N pairs as protal takes them,
and runs protal's paired-end seeding on T threads: per pair, `PrepareLookups` for both mates; a spin of 1.5 µs (the
next pair's reading and k-mers); the next pair's keys prefetched and both mates' `FindSeeds`; a spin of 6 µs (the rest
of the pair). The spins stand for the work between, without touching memory. Variants alternate per batch of 500 pairs;
each is timed inside the seeding only, as protal's timer is.

| variant | what it measures |
|---|---|
| P0 as protal | the reference |
| P1 values a pair ahead | the lookups made one pair earlier (keys two pairs ahead): do the values arrive too late? |
| P2 whole blocks prefetched | every line of a block prefetched, not only the first: do the later flex lines and the entries wait? |
| P3 | P1 and P2 |
| P4 four pairs together | four pairs' lookups made and prefetched at once: more misses in flight |
| B0 the key map only | the lookups made, no values read |
| B1 key map and scans | the flex cells scored, no entries read, no seeds |
| B2 warm caches | each pair seeded twice, the second timed: the arithmetic |

P0-P4 emit the same seeds (a hash over all seeds of the first 20,000 pairs, per variant). B1 also counts each block's
best score by the block's size: the share of lookups (and of cells) whose best is 16, 15, 14, ...

## The cluster job

On one node of the kind the runs use, with the r226 v24 database and the B10 sample
([seedbench_cluster.sh](scripts/seedbench_cluster.sh)): it builds protal at the checkout's HEAD with `-g` and seedbench
beside it, then

1. seedbench on the first 1M pairs, on 8 and on 32 threads (`seedbench_8.tsv`, `seedbench_32.tsv`);
2. protal on the first 2M pairs, aligning only, 8 threads, under `perf stat` (cycles, instructions, cache and TLB
   misses), and the index's huge pages once it is loaded (`AnonHugePages` of the process);
3. the same run under `perf record`: the hottest symbols, and the hottest instructions of the hottest three with
   their source lines (`perf.data` is kept to annotate again).

```bash
cd ~/dev/protalFastHit && git pull
D=docs/claude/2026-10-12-seeding-profile/scripts
sbatch --exclusive -c 32 --mem 100G --wrap "bash $D/seedbench_cluster.sh $PWD \
    /hpc-home/hildebra/DB/protal/protal0.8.1_r226_v24/protal_db/database.protal \
    /hpc-home/hildebra/data/raw/Illumina/PEARL_AGE/00_fastq/B10_R1_001.fastq.gz \
    /hpc-home/hildebra/data/raw/Illumina/PEARL_AGE/00_fastq/B10_R2_001.fastq.gz runs/seedbench"
```

**How to read it**, with P0 = 1:

- B2 (warm) is the arithmetic. If it is near 1, the waits are gone already and only fewer instructions help.
- B0 is the key map's share, B1 − B0 the flex cells', P0 − B1 the entries' and the seeds'.
- P1 or P3 well below 1: the values come too late; making the lookups a pair earlier is a small change to the
  paired-end and single-end loops (`Classify.h`), exact.
- P2 below 1: the later lines of a block wait; prefetch more of each block (or the entries of the best cells as soon
  as they are known), exact.
- P4 below P1: more misses in flight per thread pay; interleave several pairs' lookups, exact.
- The best-score counts of the 256-4,095 blocks: if most lookups there score 16 (or 15), the large blocks could be
  grouped by their first flex bases at load time, so that a lookup with a high best score skips the groups that cannot
  reach it. That cuts 41% of the cells streamed. Exact, if the ties are emitted in entry order. If the best scores
  there are low (reads of species the database lacks), it would not pay.

## Results at r226: the cluster job (2026-10-12)

SLURM 24188138 on `q512n10` (the 16th run's node), `3931e48`, the r226 v24 database: seedbench on the first 1M pairs of
B10, `perf` on the first 2M ([cluster_run1.txt](results/cluster_run1.txt)).

| variant | 8 threads: ns per lookup | relative | 32 threads: ns per lookup | relative |
|---|---|---|---|---|
| P0 as protal | 174.6 | 1.000 | 180.4 | 1.000 |
| P1 the lookups a pair earlier | 164.5 | 0.942 | 171.8 | 0.952 |
| **P2 whole blocks prefetched** | **157.4** | **0.902** | 164.1 | 0.910 |
| P3 P1 and P2 | 156.1 | 0.894 | 162.9 | 0.903 |
| P4 four pairs interleaved | 166.5 | 0.954 | 174.2 | 0.965 |
| B0 the key map only | 19.5 | 0.112 | 20.4 | 0.113 |
| B1 the key map and the scans | 124.4 | 0.713 | 130.8 | 0.725 |
| B2 warm caches (the arithmetic) | 74.8 | 0.428 | 69.1 | 0.383 |

P0-P4 emitted the same seeds. The benchmark's 174.6 ns are close to the run's own 159 ns (its timer, its gaps).

- **The key map: 19.5 ns** (11%). The pair-ahead prefetch of `87b857e` hides it.
- **The flex cells: 105 ns** (B1 − B0, 60%). The AVX-512 arithmetic of ~70 cells is ~10-15 ns of it; the rest waits for
  ~5 lines a lookup.
- **The entries and seeds: 50 ns** (P0 − B1, 29%): the best cells' entries (lines of the block found only after the
  scan), their decoding, the seeds emitted, the lookups sorted by block size.
- **The arithmetic: 75 ns** (B2, 43%), with the data in cache. **~100 ns of a lookup (57%) wait on memory.**
- **Prefetching more**: P2 (every line of each of the next 16 lookups' blocks, up to 32 lines, flex cells and entries,
  instead of their first lines) −10%; P1 −6%; both −11%; four pairs interleaved −5%. The misses are not late, they are
  many: 44 lookups a pair, ~10 lines each, and a core keeps only so many in flight. On 32 threads a lookup costs 3%
  more and the variants pay the same: the limit is per thread, not the node's memory.

**The blocks' best scores** (8 threads; the share of the lookups whose best flex cell shares 16, 15, ... of the read's
16 flank bases):

| block size | lookups | 16 | 15 | 14 | 13 | 12 | ≤ 11 | cells in lookups of ≤ 11 |
|---|---|---|---|---|---|---|---|---|
| 2-15 | 4.97M | 1.0% | 0.2% | 0.1% | 0.1% | 0.2% | **98.5%** | 98.2% |
| 16-255 | 10.39M | 3.1% | 0.5% | 0.2% | 0.3% | 0.8% | **95.0%** | 93.4% |
| 256-4,095 | 0.93M | 15.0% | 3.0% | 1.9% | 1.2% | 4.3% | **74.6%** | 73.3% |
| all | 16.29M | | | | | | **94.9%** | 85.2% of all cells |

Sixteen random bases share four with the read on average, and the best of 70 such cells ~9. So 95% of the lookups are
matches of the core by chance: at r226 nearly every 15-base core is in the index (99.2% of the k-mers), and a
metagenome's reads come mostly from outside the marker genes (101M mates, 2.2M SAM records written). Those lookups
cost as much as any other. Their seeds pair by chance often enough that 79% of the mates get an anchor (79.5M of
101.2M), and the k-mer screen then refuses 90% of the 176M candidates (159M).

**`perf stat`** (protal on 2M pairs, the index load included, which dominates): IPC 2.42; 183M L2 TLB misses (page
walks), ~2 per lookup if all were the seeding's. Whether the index lies on huge pages was not measured: the script read
`smaps_rollup` when protal's (buffered) log reached the index line, which was at the run's end, after the index was
freed. THP is `always` on the node and both arrays are `madvise`d; the script now polls `smaps_rollup` during the run.

**`perf record`** (the same run): two thirds of the samples are the index load. Of the seeding's symbols,
`flex_scan::ScanAvx512` is the largest (6.9% of all samples: its loads of the flex cells wait), then `PrepareLookups`
1.6%, `GetFromLookup` 1.3%, the finder's `operator()` 1.2%, `FindSeeds` 1.0%, `SharedSeeds::Sort` 0.9%, the k-mers
(`ScanWindowsAvx512`, `FillAvx512`) 1.4%.

**What it means for a pe run on 8 threads** (seeding 44.2 s per thread, split as the benchmark splits it): the key map
~5 s, the flex cells ~26.5 s, the entries and seeds ~12.6 s; about 25 s of the 44 are waits.

- **Whole blocks prefetched (P2)**: −10% of the seeding, ~−4.4 s per thread, exact. A change to
  `KmerLookupSM::PrefetchValues` (it needs the entry width, which the map has).
- **Giving up on chance matches**: if the 95% of lookups that match by chance come from mates that never align, and such
  a mate could be given up after ~4 of its 22 lookups, its seeding would fall by ~80%: the seeding by up to ~70%
  (~30 s per thread), with the seed sort, the pairing, the extension and the screen of their anchors. An upper bound,
  not exact: what decides it is how the lookups of the mates that do end up aligned score (counted per mate: the best
  score of each of its first lookups, and whether it was written). A rule that gives up none of those would change no
  alignment; the full sample's outputs and the r226 F1 benchmark then confirm it.

## Implemented: whole blocks prefetched (2026-10-12)

P2 in protal: `KmerLookupSM::PrefetchValues` fetches every cache line of a lookup's block, its flex cells and its entries,
up to 32 lines (`kPrefetchLines`), where it fetched the first line of each. Its callers are unchanged: `PrepareLookups`
for every lookup of the next pair, `FindSeeds` for the lookup 16 ahead, for short and long reads alike. To know a
block's end, `Seedmap::PackedBlock` carries its byte span (`span`, set by `BlockOfSlots` from the entry width; it fits
the struct's padding, so `LookupPointer` stays 48 bytes). seedbench's P0 keeps the old prefetch, so that it still
measures against it.

Checked (WSL, `a315c72` with the change against `a315c72`, [p2_check.sh](scripts/p2_check.sh)): unit tests 516 run, 513
passed, 3 skipped; e2e 155 OK on a fresh mini database; the same outputs as before (SAM records sorted, profiles) for
the archaea world's 1M pairs and the small real world's pe, se, PacBio and ONT samples. The local world cannot show the
gain (4 cells a block, mostly a single line).

At r226 the benchmark predicts −10% of the seeding: ~4.4 s per thread of a pe run on 8 threads (~4% of it), ~1.2 s
on 32. The next cluster run shows it: `measure_performance.sh` as for the 16th run, on this commit.

**The per-mate count for giving up** is part of another session's `--stage1_probe` (in progress when this was written,
uncommitted): per mate, the best score among its first 1, 2, 4, 8 and all lookups, by whether the mate aligns, only its
fragment does, or neither. Its output on the r226 sample decides the give-up rule; this report takes it up then.

## Options, with what the earlier rounds did

From the GTDB-scale report's list (`../2026-10-04-performance-gtdb-scale`, "The seeding: how a lookup works"):

| option | status |
|---|---|
| the flex scan vectorised | done: AVX2 `cbdb069`, the ties by masks, AVX-512 `f53db5e` |
| prefetch: the key map a pair ahead, the values' first lines before seeding | done `87b857e` |
| huge pages for the key map and the values (`madvise`) | done; whether r226's 30 GB get them is still to see (the job read them at the wrong moment) |
| more threads (the node has 84 cores) | a user's choice; a lookup costs 3% more on 32 threads than on 8 (the benchmark) |
| more misses in flight (several pairs interleaved; the lookups a pair earlier) | measured: −5% (P4), −6% (P1); not worth more than P2 |
| **more of each block prefetched** | **measured −10% (P2), exact; implemented (whole blocks, up to 32 lines)** |
| large blocks grouped by flex bases, the best found without a full scan | rejected: 75% of the lookups into large blocks score ≤ 11 and would scan them all |
| small blocks next to their key (one miss instead of two) | not worth it: the key map is 11% of a lookup |
| fewer lookups: **giving up on a mate whose lookups match by chance**; sparser syncmers; a lower seed target | not exact; the first is the large lever (up to ~70% of the seeding), to count first (the mates that align: how their lookups score) |

## How it was run

In WSL under `~/seed` (kept until the cluster results are in): [sa_setup.sh](scripts/sa_setup.sh) (the build, the
archaea world's database, the samples), [sa_cg.sh](scripts/sa_cg.sh) (the 1M-pair run's seeding line, callgrind on
100,000 pairs, seedbench), [sa_bench_compile.sh](scripts/sa_bench_compile.sh) and [sa_bench_run.sh](scripts/sa_bench_run.sh)
(seedbench again with the monotonic clock). Results in [results/](results/).
