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
- **Status**: the local part is done. The r226 part needs a cluster job ([below](#the-cluster-job)), because the
  local world cannot stand in for r226 here.

## Summary

1. **What the 62 s per thread are** (r226, pe, 8 threads): the k-mer lookups ("Seeding") 44.2 s (71%), extending the
   anchors 8.3 s (13%), sorting the seeds 4.3 s (7%), pairing 2.3 s, recovering anchors from the mate 2.1 s, sorting the
   anchors 0.6 s. Taking the k-mers ("Retrieve k-mers", 6.0 s) is timed apart.
2. **A lookup costs ~159 ns of a thread**: 2.22G lookups, nearly every k-mer's core found in the index (99.2%), every
   block scanned (the stop at 128 seeds never comes), 70.4 flex cells and 2.0 seeds per lookup. Its arithmetic is at
   most ~60 ns (~400 instructions: the key map's control block, the AVX-512 scan of ~70 cells, two entries decoded and
   emitted), so **about two thirds is waiting on memory**: a line of the 3.2 GB key map, ~5 lines of flex cells and
   ~2 lines of entries in the ~18 GB of values, per lookup. Earlier rounds prefetch the key map a pair ahead and the
   first line of each block's values just before seeding (`87b857e`), and vectorised the scan (`f53db5e`).
3. **The local world cannot answer which of those waits dominate.** Its k-mer cores are in the index for 43% of the
   k-mers (r226: 99.2%), and its blocks hold 4 cells (r226: 70). Its 3.2 GB key map and 0.4 GB of values are not r226's
   3.2 + 18 GB either. What it gives exactly (callgrind, 100,000 pairs, one thread): **82 instructions to make a lookup
   (per k-mer), 226 to score a block and emit its seeds** (with 3.2 cells), 47 of them the AVX2 scan; the seed and
   anchor finding is 30% of the alignment loop's instructions, taking the k-mers 24%.
4. **A benchmark of the seeding in parts and variants** ([seedbench.cpp](scripts/seedbench.cpp)) runs protal's own
   lookups on a real database and real reads, in protal's pipeline. It times: the key map alone; adding the scans;
   the whole seeding; the whole seeding with warm caches (the arithmetic); and four exact variants that put more misses
   in flight. Every exact variant gives the same seeds (checked by a hash over all seeds). It also counts how often a
   block's best score is 16 (the read's 16 flank bases occur exactly), which decides whether sorting the large blocks
   by flex cell would pay. Locally it works, but here its numbers say nothing about r226.
5. **Seen in passing** (exact, locally): ~17 heap allocations per pair (`operator new` with its `malloc`, and `free`:
   ~0.45 G instructions, up to ~10% of the loop; a `resize` per read in the k-mer handler, the timers, the column-weights
   cache), and taking the k-mers is 24% of the
   loop's instructions (6.0 s per thread at r226, compute-bound). Both are smaller than the lookups' waits.

**Next**: the cluster job below (one node, ~15 minutes). Its four numbers (the key map alone, plus the scans, the
whole, warm) split the 159 ns, and the variants and the best-score counts say which of the options at the end pays.

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

So a lookup reads ~8 lines from ~21 GB of index: ~1.1 TB over the run's seeding, ~26 GB/s on 8 threads. That is far
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
| the index in memory | 3.2 GB key map, 0.4 GB values | 3.2 GB key map, ~18 GB values |

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

## Options, with what the earlier rounds did

From the GTDB-scale report's list (`../2026-10-04-performance-gtdb-scale`, "The seeding: how a lookup works"):

| option | status |
|---|---|
| the flex scan vectorised | done: AVX2 `cbdb069`, the ties by masks, AVX-512 `f53db5e` |
| prefetch: the key map a pair ahead, the values' first lines before seeding | done `87b857e` |
| huge pages for the key map and the values (`madvise`) | done; whether r226's 21 GB get them is step 2 of the job |
| more threads (the node has 84 cores) | a user's choice; the lookup's cost per thread barely rises with threads (159 ns on 8, 186 ns on 32) |
| more misses in flight (several pairs interleaved) | open: P1, P4 |
| more of each block prefetched | open: P2 |
| large blocks grouped by flex bases, the best found without a full scan | open: the best-score counts |
| small blocks next to their key (one miss instead of two) | open; 31% of r226's blocks have 2-15 cells, but a lookup's flex lines matter more there than its key line |
| fewer lookups (sparser syncmers, a lower seed target) | not exact: changes sensitivity; after the exact ones |

## How it was run

In WSL under `~/seed` (kept until the cluster results are in): [sa_setup.sh](scripts/sa_setup.sh) (the build, the
archaea world's database, the samples), [sa_cg.sh](scripts/sa_cg.sh) (the 1M-pair run's seeding line, callgrind on
100,000 pairs, seedbench), [sa_bench_compile.sh](scripts/sa_bench_compile.sh) and [sa_bench_run.sh](scripts/sa_bench_run.sh)
(seedbench again with the monotonic clock). Results in [results/](results/).
