# Speed of protal 0.8.4 at GTDB r226 (cluster run Perf0.8.4a), and where to improve it

Date: 2026-10-15. Branch `audit-fixes`.

- **Cluster run:** protal 0.8.4 at `5369c99`, on the user's cluster as `local/Perf0.8.4a` (`environment.txt`):
  - SLURM 24243431 on q512n20 (AMD EPYC 9634, 84 cores, 502 GB), 32 threads, transparent huge pages `always`.
  - Database: r226 v24 (`protal0.8.1_r226_v24`, 30 GB, on NFS). Run folder also on NFS.
  - Samples: pe1 (PEARL_AGE B10, 50.6 M pairs, gzip of 290,592 members per file), pb2 (one HiFi barcode, 497,656
    reads), ont3 (soil ONT, 3.10 M reads). The HiFi and ONT files are single-member gzip.
- **What the run did:** `scripts/measure_performance.sh` with three repeats per sample, the cohort (all three from
  their SAMs, strain MSAs, qcMSA), the comparison profiler on the same reads, and `--test_sparse` (one run per gate:
  pe 1 2 4 8, long reads 4 8 16 32).
- **Earlier runs compared:** all on 32 threads and r226 v24 unless noted.
  - `local/Perf0.8.2g`: 0.8.2 `5222386`, q512n10.
  - `local/Perf0.8.3a`: 0.8.3 `9c27756`, q512n11.
  - `local/Perf0.7.5_v13`: 0.7.7 `230bf64`, on the r226 v11 database.
- **Local trials:** at `4ec3f65` (HEAD; the same alignment and input code as `5369c99`), WSL on cores 0-3, with data
  from `~/Mocks`. The scripts are in [scripts/](scripts/) and the raw outputs in [results/](results/). They are
  described in section 6.

Below, "warm" means repeats 2 and 3. Repeat 1 of each sample read its input, and for pe1 also the database, from NFS
for the first time (section 4).

## Summary

- **0.8.4 is the fastest 0.8.x so far on paired-end reads.**
  - pe1 warm: 35.3 and 35.8 s. 0.8.2 took 38.0-38.2 s, 0.8.3 took 67-68 s (the gzip-members regression, fixed in
    `0ab85f4`), and 0.7.7 took 33.1-33.6 s on a smaller database.
  - HiFi takes 12.2 s, as in 0.8.2/0.8.3 (11.6-13.8 s).
  - ONT takes 51.8-54.3 s (0.8.3: 50.5-50.9 s).
- **Peak memory dropped.** It is 35.1-36.2 GB, against 37.5-38.2 GB in 0.8.2 and 39.4-45.5 GB in 0.8.3. The key map is
  0.5 GB instead of 3.2 GB, and the compressed input is no longer mapped whole.
- **One regression: every k-mer lookup costs more.**
  - Per lookup, the seeding takes +18.5% for pe (166 → 197 ns), +13% for PacBio and +10% for ONT. The gate's lookups
    take +21% (119 → 144 ns).
  - The work is identical: the same k-mers, blocks, flex cells, seeds and anchors.
  - The node is not slower: the comparison profiler ran at the same cycle counts there as on the other nodes.
  - Each lookup now misses the cache 1.3 (pe) to 2.0 (long reads) times more often. That is the compact key map's
    second, dependent cache line, now the default (`ab6bf25`).
  - Cost: about 2.1 s per thread of pe seeding (≈ 1.5-2 s of pe wall time, ~5%), ~0.3 s for HiFi, ~1 s for ONT.
  - To confirm: one A/B run with `--direct_keys` (section 7).
- **The database gates are where 0.8.4 is fast.** A database built since `4dfc549` gates pe at C 4, pb at C 16 and ONT
  at C 8 by default. At those gates, on this database:
  - pe: 21.2 s (−40%), records −0.20%.
  - HiFi: 9.2 s (−24%), records −0.36%.
  - ONT: 32.8 s (−38%), records −0.34%.
  - The stride pick keeps far more pe records than 0.8.3's hash pick did: −0.20% against −2.1% at C 4.
- **Gated long-read runs are limited by the gzip input, not the CPU.**
  - ONT aligning stays at 26.3-28.1 s from C 4 to C 32. Meanwhile, the time threads spend waiting for reads grows from
    2.0 to 7.3 s per thread.
  - HiFi aligning stays at ~4 s from C 8 on, with 43% of it waiting for reads at C 16.
  - Locally, ParallelGzip reads long-read FASTQ at 0.9 GB/s on one thread and 1.4-1.6 GB/s on four. Retuning its knobs
    gained nothing (section 6.1).
  - A model fitted to the local measurements predicts ~+40% input speed with 16 inflate threads instead of 8, which
    `PROTAL_INFLATE_THREADS` can test on the cluster without a code change.
- **Against the comparison profiler on the same node:**
  - Ungated, protal is 1.9× faster on pe, as fast on HiFi, and 1.25× slower on ONT.
  - At the database gates it is faster on all three: 3.2×, 1.3× and 1.3×.
  - It uses ~35 GB against 4-5 GB.

## 1. Wall time and memory

Warm repeats (2 / 3). Times in seconds. Peak memory in GB (max RSS).

| sample | build | wall | aligning | profiling | index load | peak | user CPU |
|---|---|---|---|---|---|---|---|
| pe1 | 0.7.7 `230bf64` (v11 db) | 33.6 / 33.1 | 25.5 / 24.8 | 3.2 / 3.4 | 2.7 | 34.2 | 954 / 955 |
| | 0.8.2 `5222386` | 38.0 / 38.2 | 23.6 / 24.5 | 7.9 / 6.8 | 2.7 | 37.7 | 944 / 947 |
| | 0.8.3 `9c27756` | 67.2 / 68.0 | 58.2 / 58.3 | 4.6 / 4.7 | 2.8 | 45.0 | 1,078 / 1,079 |
| | **0.8.4 `5369c99`** | **35.3 / 35.8** | 25.5 / 25.7 | 4.75 / 4.86 | 3.07 / 3.09 | **35.1** | 1,031 / 1,039 |
| pb2 | 0.7.7 (v11 db) | 12.3 / 12.3 | 7.9 / 7.9 | 0.38 / 0.35 | 2.7 | 34.5 | 355 / 356 |
| | 0.8.2 | 13.8 / 12.7 | 7.7 / 7.6 | 0.83 / 0.80 | 2.7 | 37.5 | 325 / 324 |
| | 0.8.3 | 11.8 / 11.6 | 6.74 / 6.73 | 0.59 / 0.58 | 2.8 | 39.4 | 327 / 326 |
| | **0.8.4** | **12.2 / 15.0** | 6.71 / 6.90 | 0.84 / 1.97 | 3.08 | **35.5** | 346 / 345 |
| ont3 | 0.8.2 | 52.0 / 51.4 | 44.0 / 44.2 | 2.97 / 2.44 | 2.7 | 38.1 | 1,520 / 1,521 |
| | 0.8.3 | 50.5 / 50.9 | 44.2 / 44.6 | 1.74 / 1.74 | 2.8 | 45.4 | 1,548 / 1,550 |
| | **0.8.4** | **54.3 / 51.8** | 45.4 / 45.0 | 4.15 / 1.72 | 3.07 | **36.2** | 1,604 / 1,594 |

- **Noise from writing to NFS.** pb2 repeat 3 (15.0 s) and ont3 repeat 2 (54.3 s) lost 1.2-2.4 s writing their
  outputs to the NFS run folder:
  - pb2 rep 3: "writing the profile" 1.50 s against 0.31 s, and the SAM finish 1.62 s against 0.53 s.
  - ont3 rep 2: "writing the profile" 2.73 s against 0.31 s.
  
  Without those outliers, HiFi takes 12.2 s and ONT ~52 s.
- **Index load +0.3 s** (2.74-2.78 → 3.07-3.09 s): the compact key map is built while the index loads.
- **Memory −2.6 GB against 0.8.2, −9.9 GB against 0.8.3 (pe).**
  - The key map is 0.5 GB compact, against 3.2 GB direct ("key map 0.5 GB (compact: 98629120 cores …)" in the log).
  - The genome preload peaks 2.8 GB lower (28.6 against 31.4 GB).
  - 0.8.3's extra 7.3 GB of mapped gzip input is gone.
- **pe user CPU is +87 s against 0.8.2,** even though the wall time is 2.6 s lower:
  - Seeding added 69 thread-s (section 3).
  - Profiling saved 2-3 s, mostly from 0.8.3's changes: "writing the profile" went from 1.33 to 0.36 s and "after the
    profile" from 1.93 to 0.28 s.
  - Untimed time fell from 3.8 to 2.0 s.

## 2. Where the time goes

Seconds per alignment thread, mean of the warm repeats.

| stage | pe1 0.8.2 | pe1 0.8.3 | pe1 0.8.4 | pb2 0.8.2 | pb2 0.8.4 | ont3 0.8.2 | ont3 0.8.4 |
|---|---|---|---|---|---|---|---|
| Seeding | 11.55 | 11.58 | **13.72** | 2.30 | **2.54** | 10.09 | **10.90** |
| Sorting seeds | 1.09 | 1.11 | 1.08 | 0.35 | 0.34 | 1.26 | 1.22 |
| Pairing | 0.55 | 0.56 | 0.55 | 0.13 | 0.12 | 3.13 | 1.07 |
| Extending anchors | 2.09 | 2.15 | 2.13 | 0.55 | 0.55 | 2.03 | 2.06 |
| Retrieve k-mers | 1.49 | 1.50 | 1.50 | – | – | – | – |
| (Seed- and anchor-finding, total) | 16.05 | 16.16 | 18.23 | – | – | – | – |
| Alignment handler | 3.09 | 3.16 | 3.10 | 6.48 | 6.57 | 42.48 | 42.68 |
| of it: k-mer screen / indel bound | 1.37 | 1.40 | 1.37 | 0.46 | 0.44 | 22.58 | 22.77 |
| Sequence reader (waiting for reads) | 1.27 | **39.25** | 1.16 | 1.49 | 0.30 | 0.49 | 0.38 |
| Output handler | 0.62 | 0.82 | 0.67 | 0.11 | 0.10 | 0.20 | 0.18 |

- **pe:** seed and anchor finding is 71% of the 25.5 s of aligning, and the seeding alone is 54%.
- **pe profiling (4.75 s):**
  - Reading the SAM: 2.6 s. Of that, parsing is 29.1 thread-s on 32 threads, ~13 µs a record. "Adding the chunks'
    counts" takes 0.79 s on one thread. The reading thread waits 1.26 s for room.
  - Record evidence and sample context: 1.38 s.
  - Read EM: 0.62 s.
  - Congener distances: 0.24 s.
  - Writing the profile: 0.3-0.4 s.
- **HiFi:** the alignment handler (WFA2) is 6.6 of 6.7 s per thread.
- **ONT:** the indel bound alone is half of the alignment handler (it refuses 32.7 M of 35.1 M candidates).
- **Fixed costs per run:** start-up (index load 3.1 s, the tables beside it) ~3.8 s, the SAM finish ~0.6 s, and freeing
  memory ~0.3 s. A map of many samples pays the start-up once.

## 3. The lookups are slower: identical work, more cache misses

| | pe1 0.8.2 | pe1 0.8.3 | pe1 0.8.4 | pb2 0.8.2 | pb2 0.8.4 | ont3 0.8.2 | ont3 0.8.4 |
|---|---|---|---|---|---|---|---|
| k-mers found in the index (M) | 2,224.5 | 2,224.5 | 2,224.5 | 470.0 | 470.0 | 1,581.4 | 1,581.4 |
| flex cells scanned (G) | 156.36 | 156.36 | 156.36 | 34.89 | 34.89 | 253.26 | 253.26 |
| seeding ns per lookup (thread-s / lookups) | 166 | 167 | **197** | 157 | **173** | 204 | **221** |
| cycles (G) | 3,437 | 3,890 | 3,742 | 1,142 | 1,228 | 5,563 | 5,817 |
| cache misses (G) | 30.33 | 30.87 | **33.14** | 8.62 | **9.48** | 32.22 | **35.31** |
| extra misses per lookup vs 0.8.2 | | | **+1.26** | | **+1.83** | | **+1.95** |

- **The work is the same.** Seeds (4,373,336,840 pe), anchors and every other count are identical in all three builds.
- **The node is as fast as the others.** The comparison profiler ran there at the same cycle counts:
  - pe: 941-951 G on q512n20, 949-980 G on q512n11, 942-947 G on q512n19.
  - HiFi: 66-67 G against 65-66 G.
  - ONT: 247 G against 246-247 G.
- **The other stages did not change:** sorting, extension, the screen and the alignment handler are all within ±2%.
- **The gate's lookups at C 1** (the same 2,224 M sampled lookups in 0.8.3 and 0.8.4): 8.24 → 10.04 s per thread, 119 →
  144 ns a lookup.
- **What changed in the lookup between `9c27756` and `5369c99`:**
  - the compact key map, built in `47e3cc2`, fixed in `27f4c74`, the default since `ab6bf25`;
  - the direct map's one-line read (D2);
  - grid-level checks per lookup (`LookupPointer` gained `min_level`/`level_entries`; the loaded level is compared on
    every block);
  - the lookup loop compiled for x86-64-v3 too.
  
  Of these, only the compact map adds a cache line per lookup:
  - a lookup reads the presence/rank line (153 MB array);
  - then, depending on it, a 3-byte start in a 296 MB array and a group base in a 12 MB one;
  - the direct map reads one line of its 3.2 GB.
  
  +1.3 to +2.0 misses per lookup matches that second line, plus part of the bases array.
- **Why the local tests saw no cost.** The tests of 2026-10-14 found compact = direct on real reads. They used the
  reduced r226 database, whose starts fit in the cache; at full r226 they do not. The synthetic r226-density test
  then put the cost at +12-13 ns a lookup, and r226 measures +31 ns (pe).
  - The synthetic test timed the key map step alone. In a real run, each lookup also reads ~10 lines of its block, so
    an extra outstanding miss costs more than it does in isolation.

What it costs, if the cause is confirmed (section 7):
- Ungated: ~2.1 s per thread of pe seeding (≈ 1.5-2 s of wall time, 4-6%), ~0.3 s per thread HiFi, ~0.8 s ONT.
- Gated pe at C 4: ~25 ns × (594 M gate lookups + 194 M lookups of kept reads) / 32 threads ≈ 0.6 s per thread, ~6%
  of its 10 s of aligning.
- The saving bought with it: 2.7 GB of 35 GB (7.7%).

## 4. Cold runs, the cohort, and the comparison profiler

- **Cold runs.** pe1 repeat 1 loaded the 30 GB database from NFS in 129 s (0.23 GB/s; wall 165.6 s).
  - Earlier first runs took 369 s (0.7.7), 935 s (0.8.2e) and 453 s (0.8.3): NFS throughput varies 4×.
  - The first read of pb2's and ont3's input from NFS cost +5.8 s and +25 s. The alignment threads waited 2.6 and 31 s
    per thread (reader share 0.20 and 0.43).
  - This is the network file system, not protal. It is a reason to keep databases and inputs on node-local disks
    (`--scratch_dir` covers protal's own intermediates only).
- **Cohort** (profile-only from the three SAMs, strain MSAs of 130 species, qcMSA): 17.1 s wall, against 16.0 s in
  0.8.3 and 18.4 s in 0.8.2. Peak 13.1 GB.
  - qcMSA, summed over the species, took 119.6 s against 64.6 s in 0.8.3. `scripts/qcmsa.py` did not change between
    the builds (`git diff 9c27756 5369c99 -- scripts/qcmsa.py` is empty), so this is the node's Python environment.
  - The strain stage took 6.2 s against 4.6 s of wall time.
- **Comparison profiler,** same node, reads and threads:

  | | its wall (warm) | its peak | protal ungated | protal at the database gate |
  |---|---|---|---|---|
  | pe | 66.4-70.5 s | 4.5 GB | 35.5 s (1.9× faster) | C 4: 21.2 s (3.2×) |
  | HiFi | 12.0-12.2 s | 3.8 GB | 12.2 s (equal) | C 16: 9.2 s (1.3×) |
  | ONT | 42.1-42.2 s | 5.0 GB | ~52 s (1.25× slower) | C 8: 32.8 s (1.3×) |

## 5. The database gates (`--test_sparse`, one run each)

The v24 database has no gates and only ungated models, so these runs gate by `--sparse_first` and score with the
ungated models (protal warns). A database built since `4dfc549` gates by default and has models for its gates.

| sample | gate | wall s | aligning s | records vs ungated | given up | species lost / gained | waiting for reads, s/thread |
|---|---|---|---|---|---|---|---|
| pe1 | none | 35.3 | 25.5 | – | – | – | 1.16 |
| | C 1 | 27.0 | 17.1 | −0.014% | 90.1% | 3 / 2 | 1.08 |
| | C 2 | 23.2 | 12.4 | −0.03% | 90.5% | 3 / 3 | 1.15 |
| | **C 4** | **21.2** | 10.0 | **−0.20%** | 91.0% | 3 / 3 | 1.21 |
| | C 8 | 20.3 | 8.6 | −1.19% | 92.4% | 3 / 2 | 1.19 |
| pb2 | none | 12.2 | 6.7 | – | – | – | 0.30 |
| | C 4 | 9.6 | 4.32 | −0.03% | 77.9% | 0 / 1 | 1.36 |
| | C 8 | 9.2 | 4.05 | −0.12% | 81.2% | 0 / 1 | 1.54 |
| | **C 16** | **9.2** | 3.95 | **−0.36%** | 83.4% | 0 / 1 | **1.71 (43%)** |
| | C 32 | 9.0 | 3.98 | −0.88% | 84.7% | 0 / 1 | 1.73 |
| ont3 | none | ~52 | 45.0 | – | – | – | 0.38 |
| | C 4 | 34.2 | 28.1 | −0.10% | 81.4% | 3 / 17 | 1.95 |
| | **C 8** | **32.8** | 26.9 | **−0.34%** | 87.1% | 4 / 18 | 3.77 |
| | C 16 | 32.5 | 26.5 | −1.15% | 90.8% | 3 / 23 | 5.75 |
| | C 32 | 32.3 | 26.3 | −3.68% | 93.3% | 5 / 28 | **7.32 (28%)** |

- **Even C 1 pays on short reads.** It makes the same 2.22 G lookups, but drops the 90% of pairs without a hit before
  the seed sort, pairing, screen and alignment: −24% wall time for 305 fewer records of 2.22 M.
- **Stride pick against hash pick.** pe C 4 now keeps 99.80% of the records; 0.8.3's hash pick at C 4 kept 97.9%.
  - Same at C 8: 98.8% against 92.9%.
- **ONT's calls differ from the ungated run** at every gate (Bray-Curtis 0.14, 17-28 species gained), as they did in
  0.8.2 and 0.8.3 (0.1447 at C 4 in both). That is a scoring question for gated ONT models, not a speed one.
- **The long reads stop gaining at C 4 (ONT) or C 8 (HiFi).**
  - The gate's CPU savings continue: ONT alignment CPU 572 → 320 s from C 4 to C 32.
  - The alignment threads instead wait for the input: ONT busy cores 20.3 → 12.1, protal's log saying "the input
    held the alignment back".
  - Estimated from the k-mers looked up (6.8 bp a lookup, as for pe1), ont3 is ~11 Gbp, ~22 GB of FASTQ, read in
    ~26 s: ~0.85 GB/s.
  - pb2 is ~3.2 Gbp, ~6.5 GB, read in ~4 s: ~1.6 GB/s.
  - Both are single-member gzip, inflated in pieces on 8 threads.

## 6. What to improve, and the local trials

Ranked by what they would save on the cluster with gated databases (the default from now on).

### 6.1 Single-member gzip input (gated long reads): the largest lever

**Potential:**
- ONT at C 8 spends 26.9 s aligning, but its alignment CPU (443 s over 32 threads) needs ~14 s. With an input twice as
  fast, the run would take ~20 s instead of 32.8 s.
- HiFi at C 16 would take ~7.9 s instead of 9.2 s.
- Gated pe waits 12-14% of its aligning for its (multi-member) input too.

**Local trials.** The data is the local ONT mock (ENA ERR15459768) as 2.0 GB of FASTQ: 186,729 reads, 1.0 Gbp, mean
5,355 bp. It was compressed one member at a time by gzip -6, gzip -1 and pigz -6, and ENA's own file was used too. The
HiFi mock SRR13128012 (2.0 GB of FASTQ) was the second file. Everything was read through `ThreadedGzIstream` with
`bench_pgz` ([scripts/bench_pgz.cpp](scripts/bench_pgz.cpp), from the gzip-members report).

1. **Compressors** ([results/pgz_compressors.txt](results/pgz_compressors.txt)).
   - All of them decode in pieces (none falls back to in-order decoding).
   - Long-read FASTQ compresses only 2.4× (the qualities), so one ISA-L thread gives 0.6-0.9 GB/s of it. That is about
     half the short-read figure of the ParallelGzip report.
   - 4 threads give 1.5-1.7× one thread's rate: 0.9-1.2 GB/s in the machine's slow state, 1.4-1.6 GB/s in its fast
     one. Compare within a row of the result files.
2. **Piece size and gap** ([results/pgz_piece_gap_sweep.txt](results/pgz_piece_gap_sweep.txt)).
   - Gaps of 1-3 pieces (workers closer to the reader) are slower than the default T + 2. For example ONT, 4
     threads, 1 MB pieces: gap 1 0.69 / 1.00 GB/s against the default's 1.11 / 1.34.
   - 512 KB and 2 MB pieces are within the noise of 1 MB. No setting to change.
3. **The lead's waits** ([scripts/parallel_gzip_lead_experiment.diff](scripts/parallel_gzip_lead_experiment.diff),
   [results/pgz_lead_nowait.txt](results/pgz_lead_nowait.txt)). Timers were added, plus a policy in which the reader's
   thread decodes a piece itself instead of waiting for the worker decoding it.
   - The policy gives the same speed (ONT 4 threads: 1.54-1.57 against 1.56-1.58 GB/s).
   - The lead hardly waits (0.03 s per file).
   - It spends 0.94 of 1.29 s decoding in order the ~41% of the pieces that no worker reached. Each speculative piece
     costs a worker two passes (labels A and B), so 3 workers keep up with only ~59% of the file.
4. **A model** that fits these numbers: throughput = 1 / (f / r1 + o), where
   - r1 = 0.87 GB/s is the lead's in-order rate;
   - o = 0.17 s/GB is the lead's copy, CRC and window tails;
   - f = 1 / (1 + W / 2.2) is the share the lead decodes, with W workers.
   
   It gives 1.56 GB/s for W = 3 (measured 1.56), ~2.2 for W = 7 (8 threads, protal's choice on 32 threads), ~3.1 for
   W = 15, and a ceiling of ~5.9 GB/s.
5. **zstd as the alternative input** ([scripts/zstd_pgz.sh](scripts/zstd_pgz.sh),
   [scripts/pzstd_pgz.sh](scripts/pzstd_pgz.sh)).
   - One zstd frame (`zstd -T4`) reads at 1.1 GB/s on any number of threads.
   - Multi-frame files from `split | zstd` and from `pzstd` read on one thread too (0.8-1.2 GB/s at 4 threads). Their
     frames declare no content size (`zstd -lv`: no decompressed size; frame descriptor 0x04), and protal's
     parallel-frames path takes only frames that declare one (`ThreadedGzStream.h`:
     `ZSTD_getFrameContentSize(...) <= kMaxZstdFrame`).

**Suggestions:**
- **(a) Cluster test, no code change:** `PROTAL_INFLATE_THREADS=16` for the long-read gated runs. The model predicts
  +40% input speed. The cluster has 84 CPUs and the runs use 32 + 8.
- **(b) Report it in the log** (cheap): ParallelGzip's counts (pieces ahead, in order, streamed) and the input's GB/s in
  each sample's "input" line. The cluster's ~0.85 GB/s for ont3 is below even the local 4-thread rate. Only the counts
  can tell whether that file decodes in order (a compressor whose blocks the finder misses) or is just slow.
- **(c)** If (a) pays: `GzipThreadsFor` up to 16 when the CPUs allow, rather than 8.
- **(d) Larger work:** speculative pieces in one pass instead of two. This needs a decoder that writes window
  references as 16-bit symbols, so ISA-L's stock API cannot do it; the ceiling then moves to the lead's copy and CRC
  (o).
- **(e) Small change:** take zstd frames without a declared size into the parallel path. Their ends can be found by
  walking the block headers without decompressing, and the output buffer can grow as needed. Today such a `.fq.zst`
  reads no faster than gzip.

### 6.2 The compact key map's second line

- **First, confirm** with an A/B run (section 7).
- **If confirmed, the choice is between:**
  - the direct map again for whole-index runs: +2.7 GB, −2.1 s per thread of pe seeding;
  - a layout with one line per lookup.
- **What the earlier key-map report found:**
  - It tried a 128-core bucket line, which matched the direct map on synthetic cores. Real cores cluster (17.5% of the
    reduced r226's in overflowing lines), so it lost on real reads.
  - Its "starts a read ahead" variant (Ca) saved ~2 ns against protal's (Cf), so a different prefetch order will not
    close a 31 ns gap.
- **A design not yet tried:** buckets chosen by a hash of the core, which breaks the clustering.
  - Each 64-byte line holds 6 entries of core, start and slot count; an overflow goes to the next line.
  - That is ~1.3-1.5 GB at a load of 0.7-0.8: between compact (0.45 GB) and direct (3.2 GB), with one random line
    per lookup.
  - Worth a `keybench` layout only once the A/B shows the cost is real.

### 6.3 Short-read profiling

- **Size:** 4.4-4.75 s per pe sample at r226, 21% of a gated pe run.
- **Local measurement:** a profile-only run of a 1 M-pair D6331 HKU sample on the reduced r226 database parses at
  6.9 µs a record on 4 threads (r226 on the cluster: ~13 µs). See [scripts/trial_b.sh](scripts/trial_b.sh) and
  [results/profile_only_sam_read.txt](results/profile_only_sam_read.txt).
- **callgrind** (one thread, [scripts/cg_parse.sh](scripts/cg_parse.sh),
  [results/callgrind_profile_only_inclusive.txt](results/callgrind_profile_only_inclusive.txt)): the SAM read is 2.86 G
  instructions, ~57 k a record.
  - PrepareMAPQ: 1.43 G. Of it, `RecordEvidenceCollector::NoteRecord` is 1.20 G (13.6 k a call), within which
    `column_weights::Count` is 0.19 G and the ancestry sites 0.09 G.
  - `SamReader::Advance`: 0.82 G, of it `SamFromTokens` 0.24 G.
  - `Taxon::AddSam`: 0.37 G.
  - `VariantHandler::AddAlignment`: 0.27 G.
  - `memset`: 0.28 G in 393,624 calls (~2.8 KB zeroed per call, ~8 a record).
  - libgcc's software popcount: 0.17 G, 538,776 calls from code that is not cloned for x86-64-v3, as the key map had
    on 2026-10-14.
- **At r226, two serial pieces** are worth looking at first:
  - "adding the chunks' counts" (0.79 s on one thread; 0.01 s locally): the evidence of 1,646 chunks over 143,614
    taxa;
  - the reading thread's 1.26 s waiting for room.
  
  Together with a cheaper parse, ~1-1.5 s per pe sample.

### 6.4 Smaller items

- **Start-up, ~3.8 s per run:** index load 3.1 s on 32 threads. It is paid once per run, so it matters for one-sample
  runs like these timings, not for a map of many samples.
- **NFS:** outputs on NFS added 1.2-2.4 s to two of the nine warm runs, and the first read of a database or an input
  costs 5-130 s. This is cluster practice rather than protal code.

## 7. Next cluster run

The same script and samples, ideally on the same node:

```bash
PROTAL_ARGS=--direct_keys bash scripts/measure_performance.sh --threads 32 ... OUT_direct DB pe:R1:R2 pb:R ont:R
PROTAL_INFLATE_THREADS=16 bash scripts/measure_performance.sh --threads 32 --test_sparse ... OUT_inflate16 DB pe:R1:R2 pb:R ont:R
```

- **The first run settles section 3:** seeding per lookup back to ~166 ns, misses −1.3 per lookup, if the compact key
  map is the cause.
- **The second run settles 6.1(a):** ONT aligning at C 8 below 26.9 s, and HiFi at C 16 below 3.95 s, if more
  inflate threads help.
- A run on a database built with gates (r226 v27 or later) gives the true default speeds and its own models.

## Files

- **[scripts/](scripts/):** the WSL scripts of the trials (they read from the session's scratch folder, `SP=`), the
  bench source, and the ParallelGzip experiment as a diff on `4ec3f65`.
  - `setup_pgz.sh`: the bench and the test files.
  - `run_pgz.sh`: compressors.
  - `sweep_pgz.sh`: piece size and gap.
  - `nowait_pgz.sh`: the lead's waits.
  - `zstd_pgz.sh`, `pzstd_pgz.sh`: zstd input.
  - `build_protal.sh`, `trial_b.sh`, `cg_parse.sh`: profiling.
- **[results/](results/):** their outputs.
- **Removed after the trials:** the WSL tree `~/perf085` (the build, ~7 GB of test FASTQ and gzip/zstd copies, the
  callgrind output). The cluster data stays in `local/Perf0.8.4a` (not in git).
