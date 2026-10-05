# Performance at GTDB scale: protal 0.7.5 on the r226 database with real samples

- **Date**: 2026-10-04.
- **Code**: protal v0.7.5, commit `281a4ba` (branch `audit-fixes`), the binary the cluster job ran.
- **Measurement**: `scripts/measure_performance.sh` (SLURM job 23933082, node `q512n10`): three runs per sample, one
  sample at a time, 32 threads, `--no_qcmsa`, protal's outputs removed after each run. The script's stdout, `runs.tsv`,
  `stages.tsv`, `environment.txt` and the six protal logs are copied unchanged into [`results/`](results/).
- **Machine**: AMD EPYC 9634 (84 cores, one socket, one NUMA node, 384 MiB L3, 3.7 GHz max), 502 GB RAM, no swap,
  transparent huge pages `always`, `perf stat` allowed. The job had 32 threads of the 84 CPUs; whether the node was
  otherwise idle is not recorded (`free` showed 11 GB used at the start).
- **Database**: `protal0.7.3_r226_v10/protal_db/database.protal`, 27 GB on NFS. 143,614 species, 168 genes; the index
  2,903,557,366 entries of 42 bits with 32-bit flex cells (27.3 GB in memory) and a 3.2 GB key map; gene neighbours
  (731,664 rules), species priors, suspect copies and gene conservation tables loaded as the log says.
- **Samples**: `pe1`, a real paired-end gut sample (PEARL_AGE B10, two gzipped FASTQ files); `pb2`, a real PacBio HiFi
  sample (one demultiplexed barcode of a Sequel IIe cell). **The script does not record how many reads a sample has**,
  and protal prints the count only with `--verbose`, so no per-read or per-base cost is given below; see the gaps.
- **Earlier rounds**, all on a 765-species world in WSL at 1-6 threads: [the review of 2026-10-03](../2026-10-03-performance-review/README.md)
  (whose item 7 asked for this measurement), [pass 2](../2026-10-03-performance-pass2/README.md),
  [round 3](../2026-10-02-performance-round3/README.md), [the prefetch](../2026-10-02-prefetch/README.md),
  [the memory audit](../2026-10-03-memory-audit/README.md). This report says what the GTDB-sized database and real
  reads change about their picture.

## Summary

1. **The runs.** A warm paired-end run takes 131 s wall (median; the first, cold run 284 s), of which aligning 62 s,
   profiling 36 s, index load 3.2 s, genome preload 3.2 s, and **22-26 s that no stage timer covers** (start-up tables,
   model, output writing, teardown). The HiFi run takes 104 s: aligning 84 s, profiling 3.4 s, 9 s untimed. Both runs
   peak at 38.0 GB RSS. Run-to-run noise on the warm runs is 1-2% of wall time; instruction counts agree to 1e-5.
2. **Memory is as the memory audit predicted.** 38.0 GB for a full r226 run with the packed index (the audit expected
   ~35 GB for the index values and genes after packing, from 43 GB): the index is 27.3 + 3.2 GB as loaded, the
   genes, tables and per-run buffers the remaining 7.5 GB. No major faults on warm runs.
3. **Short reads: the k-mer lookups grew from 14% to 34% of the alignment stage.** Per thread of the 61.9 s: seeding
   20.8 s (34%), the alignment handler 27.2 s (44%), sorting seeds 3.4, extending anchors 3.6, k-mers 2.0, output 1.5,
   reader 1.0; anchor finding as a whole 29.2 s (47%). On the small world seeding was 14% and the alignment handler 34%.
   The lookups into a 27 GB index are what scales with the database (the prefetch report found seeding waiting on
   memory already on a 3.7 GB index); the whole run still reaches an IPC of 3.17 because the WFA2 kernels are dense.
4. **Short reads: at least 44.9 million fragments were aligned against and failed everywhere.** The SAM holds
   44,884,669 unmapped records; protal writes one only for a fragment whose anchors were tried (up to `align_top` 3 per
   mate, with ties) and none reached 90% ANI. Every such fragment paid its seeding and its WFA2 attempts in full, and
   a failing alignment costs the whole score budget. On the simulated world almost every fragment came from a
   database species and aligned; on real reads at GTDB scale the failing candidates are a large part of the 27.2 s
   of the alignment handler. How large needs the attempted/aligned counts (gap 1). The records themselves are evidence
   (the profiler counts the reads that fail per taxon), so the cost to attack is the alignment, not the output.
5. **Short reads: the profiling stage takes 36 s on about 3-6 of 32 cores.** Of the run's 2,277 CPU-seconds, the
   alignment accounts for ~1,920-1,980 (its per-thread stage sums equal its wall time: all 32 threads busy) and the
   index load for up to 100; what is left for profiling, start-up and output is ~200-360 CPU-s over ~62 s of wall.
   The stage is 27% of the warm wall time, as on the deep simulated sample, and more threads will not shorten it.
   The script does not split it; the review's split (parse, EM, congener distances, scoring, writing) has to be
   remeasured on this SAM (gap 2). What is known to be serial: the one thread that reads and cuts the SAM (0.9 GB/s;
   with 32 threads five zstd threads decompress), the order-dependent planning of each wave, the EM, the writers.
6. **Long reads: 91% of the alignment stage is WFA2.** The alignment handler is 83.1 of the 83.7 s per thread; the
   seeding, sorting, pairing and extending inside it sum to 7.0 s, so the whole-window and chain alignments are ~76 s
   per thread. The run's IPC is 3.85 and its cache-miss rate 1.9 per 1,000 instructions: compute-bound in the
   wavefront kernels. 430,087 reads seeded on taxa and aligned to none. The small-world picture (88% alignment in
   round 3) holds; the levers are the alignment itself: the free ends from the last link (round 3's open item 7,
   which leaves fewer reads to the whole window), WFA2's AVX2 kernels retested on long windows (they were rejected on
   short-read flanks, where the kernels are not the cost), and a screen against candidates that cannot reach 90% ANI.
7. **The cold start is the NFS read of the database: 131 s for 27 GB, about 200 MB/s**, one sequential stream; the
   reads' first pass from NFS added 3 s (paired-end) and 5 s (HiFi) to the reader. A warm load takes 3.2 s. A job that
   runs one sample on a cold node spends half its time loading; a map of samples pays it once.
8. **Threads.** The alignment stage keeps all 32 threads busy for both read types (pb2 averages 27 cores over the whole
   run), so a 64- or 84-thread run would roughly halve it; the profiling stage, the start-up and the cold load would
   not shrink. At 38 GB a node of 502 GB holds a dozen concurrent runs.

## What I would do, in order

| # | change | where | expected | outputs | effort |
|---|---|---|---|---|---|
| 1 | **Record the counts**: protal prints the reads processed, the anchors tried and the alignments made (the `Statistics` it has, today only with `--verbose`) in every run, and the script copies them into `runs.tsv`; the HiFi run also the whole-window/anchored split | `RunProtal.h`, `Classify.h`, `scripts/measure_performance.sh` | per-read costs become comparable across machines and databases; items 2 and 4 get their denominators | none | small |
| 2 | **A screen before WFA2 for candidates that cannot reach the ANI floor**: the q-gram bound (a read and a gene window within e edits share at least n − k + 1 − e·k k-mers; e from the score budget at the smallest penalty, n the part the free ends leave to be aligned) counted on the window with a small table, so a candidate below the bound is a failed candidate without an alignment; applies to short-read anchors and long-read windows alike | `AlignmentStrategy.h` (`SimpleAlignmentHandler::operator()`), `LongReads.h` (`Align`) | removes the WFA2 cost of the 44.9M fragments and 430k reads that fail everywhere, and of the failed candidates of mapped reads; the size needs item 1's counts, I would guess a quarter to a half of the alignment handler on real short reads | identical if the bound is derived conservatively (a candidate it rejects would have failed); to be shown by a test against WFA2 as the flank shortcuts were | medium |
| 3 | **Time the profiling stage on this SAM** (`KEEP=1`, then `--profile_only` with the review's timing build, or stage timers in `Profiler.h` where only `Post-processing` and `SNPs` have them) and parallelise the largest serial part; on maps, try `--profile_ahead`, which was written for a node where the stage leaves cores idle | `Profiler.h` (`ProfileSamParallel`, `ApplySampleContext`) | up to −36 s of 131 s per deep short-read sample | identical | small to measure, medium to fix |
| 4 | **Long reads**: (a) the window's free ends from the last link's diagonal (round 3's item 7), (b) WFA2 with its AVX2 kernels measured on `pb2`, (c) item 2's screen, (d) a per-segment budget from the best candidate's score | `AnchoredAlignment.h`, `WFA2Wrapper2.h`, `LongReads.h` | (a)+(b) unknown here, to measure; (d) saves the far candidates' WFA2 work | (a), (d) change some alignments; (b), (c) identical | medium |
| 5 | **Seeding at GTDB scale**: profile the lookups on this node (`perf record` on `LookupKmers`/`GetFromLookup`); then deeper prefetch than one read ahead, or lookups batched across reads by key block, and check whether `max key ubiquity` 256 walks long runs for little | `ChainAnchorFinder.h`, `Seedmap.h` | seeding is 34% of the stage; the prefetch report's gain of 4-5% was on a 3.7 GB index | identical | medium |
| 6 | **Timers for the untimed 22-26 s**: start-up (taxonomy, species priors, gene neighbours, suspect copies, conservation, model) and the output after profiling, printed like the other stages, so the next run of the script shows where they go; the paired-end run has 13 s more of it than the HiFi run, so part depends on the sample (output size or teardown) | `RunProtal.h` | knowing what 17% (pe) and 9% (pb) of a run is; the start-up is per run, not per sample | none | small |
| 7 | **Cold load**: read the database file with several threads over disjoint ranges (NFS gives more with parallel streams than the 200 MB/s of one), or stage it on node-local disk at the start of a job as the build does with `--scratch` | `Database.h`, `IndexCodec.h`, docs | 131 s → the node's parallel NFS or local read rate; matters for single-sample jobs only | none | small to medium |

Not worth it now: the index decode (3.2 s warm at 32 threads), the genome preload (3.2 s), the SAM header (0.7 s),
the strain stage (one sample, 4 ms), the output handler (1.5 s per thread, and sampled one call in 61).

## Measurements

### Whole runs (`results/runs.tsv`)

| run | wall | user | sys | max RSS | major faults | load index | aligning | profiling | instructions | IPC | cache misses |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pe1 run 1 (cold) | 283.5 s | 2,234 s | 62.0 s | 38.1 GB | 57 | 134.2 s | 63.5 s | 35.8 s | 26.03 T | 3.17 | 26.0 G |
| pe1 run 2 | 128.6 s | 2,242 s | 40.9 s | 38.0 GB | 0 | 3.2 s | 61.9 s | 36.4 s | 26.03 T | 3.16 | 26.1 G |
| pe1 run 3 | 131.5 s | 2,236 s | 40.5 s | 38.0 GB | 17 | 3.2 s | 61.8 s | 36.0 s | 26.03 T | 3.17 | 26.0 G |
| pb2 run 1 | 105.5 s | 2,778 s | 30.0 s | 38.0 GB | 0 | 3.2 s | 84.0 s | 3.6 s | 39.31 T | 3.84 | 74.7 G |
| pb2 run 2 | 103.7 s | 2,773 s | 28.3 s | 38.0 GB | 0 | 3.2 s | 83.6 s | 3.4 s | 39.31 T | 3.85 | 74.8 G |
| pb2 run 3 | 103.5 s | 2,775 s | 28.5 s | 38.0 GB | 0 | 3.2 s | 83.7 s | 3.4 s | 39.31 T | 3.85 | 74.5 G |

`perf stat` counts the whole process (all threads). Instructions per CPU-second: 11.4 G (pe1), 14.0 G (pb2), at 3.7
GHz an IPC of 3.1-3.8 as reported. Cache misses per 1,000 instructions: 1.0 (pe1), 1.9 (pb2). The paired-end run was
on average about 17 cores busy of 32, the HiFi run 27.

### A run's wall time by stage (protal's own timers in the logs; run 3 of each)

| | pe1 | pb2 |
|---|---|---|
| Preload genomes | 3.15 s | 3.17 s |
| Load index (32 threads, warm) | 3.20 s | 3.21 s |
| Aligning reads | 61.85 s | 83.70 s |
| SAM header and file | 0.66 s | 0.33 s |
| Profiling | 36.03 s | 3.40 s |
| Strain-level MSAs | 0.004 s | 0.004 s |
| **untimed** (Run protal minus the above) | **25.9 s** (run 2: 22.5 s) | **9.2 s** (run 2: 9.4 s) |
| Run protal | 130.8 s | 103.0 s |
| outside protal's timer (`/usr/bin/time` minus Run protal) | 0.6 s | 0.5 s |

The untimed time lies before "Preload genomes" (taxonomy, the species priors of 143,614 species, 731,664 gene-neighbour
rules, suspect copies, conservation factors, the model) and after "Profiling" (writing the profile and its logs,
freeing 38 GB). The start-up is the same database for both runs, so it is at most the HiFi run's 9 s; the paired-end
run's extra 13-17 s depends on the sample.

### Stage timers inside the alignment stage (`results/stages.tsv`, seconds per thread, median of 3)

| stage | pe1 | share of 61.9 s | small world, 1 thread (review) | pb2 | share of 83.7 s |
|---|---|---|---|---|---|
| Seed- and Anchor-finding (encloses the next six) | 29.15 | 47% | 3.4 s = 35% | – | – |
| Seeding (k-mer lookups) | 20.75 | 34% | 1.37 s = 14% | 4.51 | 5% |
| Sorting Seeds | 3.37 | 5% | 0.70 s = 7% | 1.32 | 2% |
| Pairing | 0.66 | 1% | 0.38 s = 4% | 0.15 | 0.2% |
| Sorting Anchors | 0.17 | 0.3% | 0.09 s = 1% | 0.04 | – |
| Extending Anchors | 3.57 | 6% | 0.65 s = 7% | 1.00 | 1% |
| Retrieve k-mers | 1.98 | 3% | 0.58 s = 6% | – | – |
| Anchor recovery | 0.51 | 1% | 0.19 s = 2% | – | – |
| Alignment handler | 27.19 | 44% | 3.28 s = 34% | 83.09 | 99% (encloses seeding to extending for long reads) |
| Joining alignment pairs and sorting | 0.49 | 1% | 0.99 s = 10% | – | – |
| Output handler | 1.53 | 2% | 0.82 s = 8% | 0.06 | 0.1% |
| Sequence reader (wait for input) | 0.97 (cold run: 3.70) | 2% | 0.37 s = 4% | 0.49 (cold run: 5.44) | 0.6% |

The per-thread stage sums match the stage's wall time (pe1 59.9 of 61.9 s; pb2 83.7 of 83.7 s), so no thread waited
on the reader or the output. For long reads the aligner seeds inside the alignment handler (`LongReadAligner::operator()`),
so its seeding, sorting, pairing and extending are part of the 83.1 s; the rest, 76 s per thread, is the alignment of
candidate windows (`Align` → `SimpleAlignmentHandler::AlignAnchor`: through the chain where it allows, else the whole
window, both WFA2) and the neighbour-gene rescue. The per-read timers time one call in 61; the output handler's
estimate is the least reliable (the review's note 6).

### What the logs say about the samples

- pe1: 44,884,669 unmapped records skipped by the profiler (fragments with anchors tried and no alignment at 90% ANI
  anywhere; `UnmappedRecord` writes nothing for a fragment without failed candidates); 233,255 fragments with one sure
  mate, of which 29,756 found on the mate's gene; 66,966 fragments paired across neighbouring genes; 21,631 records on
  suspect gene copies; 47 foreign genes; 5 species under the knob with strong reads.
- pb2: 430,087 unmapped records; 33,621 gene hits of several taxa (MAPQ < 4), 31,258 settled by their read's
  consensus, 1,552 written with MAPQ 0; 472 neighbour genes found; 356 suspect-copy records; 58 foreign genes.
- Both: "Index in memory: 2,903,557,366 entries of 42 bits ... 27.3 GB; key map 3.2 GB" (the packed layout of `d25fb66`).

### Cold against warm

The first paired-end run loaded the index in 134.2 s against 3.2 s afterwards: 27 GB from NFS at about 200 MB/s, one
sequential stream decoded in parallel; its reader stage was 3.7 s per thread against 1.0 (the FASTQ files from NFS),
and it took 57 major faults. The HiFi run's first pass, with the database already cached, only paid the reads: reader
5.4 s against 0.5 s. Everything else was within noise of the warm runs, so nothing but the file reads is cold.

## Gaps

1. **Read counts.** Not recorded; protal prints "Processed reads", "Anchors" and "Total alignments" only under
   `--verbose`, and the script does not count the FASTQ records. Without them the costs here cannot be put per
   fragment or per base, nor compared to the simulated world's 171k instructions per pair. Item 1 fixes it.
2. **The profiling stage's inner split** (parse, EM, congener distances, scoring, writing) was not measured; the
   script removes the SAM (`KEEP=1` keeps it), and the timing build of the review is not in the repository.
3. **The untimed 9-26 s** per run (table above) have no timers.
4. **Attempted against aligned candidates**: the number of WFA2 calls that fail is the size of item 2; it needs item 1
   or a counting build on the cluster.
5. One sample of each type, one node, one database; no multi-sample map, so `--profile_ahead` and the strain stage
   were not measured; whether the node was shared is not known.

## Follow-up the same day: items 1-3 implemented

Built from `c94bd0e` + this work in WSL (`scripts/build_work.sh`; the reference `c94bd0e` from `git archive`), checked
on the v0.7.3 benchmark world's database with the 500k-pair sample at one thread, the 5M-pair sample at six and
90 Mb of simulated PacBio reads at six (`scripts/compare.sh`, `scripts/deep.sh`; `results/followup_wsl.txt`). The
code is uncommitted in the working tree. Wall times on this laptop vary by 5-25% between identical runs (the second
deep round ran 25% slower than the first for both binaries), so only the counts and the identity checks are firm here.

### 1. The counts, in every run

Protal prints one line per sample after the alignment stage, from the `Statistics` it kept and the handler's
counters (`PrintAlignmentCounts`, `RunProtal.h`; the long-read aligner's handler copies are joined over the threads
like the others, `SimpleAlignmentHandler::JoinCounts`):

```
Sample s: 4995812 read pairs, 4738945 with an anchor; 12782429 candidate alignments tried: 720480 refused by the
k-mer screen, 11792557 aligned from the anchor's exact matches and 269392 as whole windows; 10456457 alignments made,
2351020 records written
```

`scripts/measure_performance.sh` records them as eight columns of `runs.tsv` (`reads`, `anchored_reads`, `tried`,
`screened`, `from_anchors`, `whole_windows`, `made`, `written`) and prints them per sample. The HPC run of this report
would have shown how many of its candidates fail; the next one will.

### 2. The k-mer screen before WFA2 (`src/Alignment/AlignmentScreen.h`)

- In `SimpleAlignmentHandler::AlignAnchor`, after the window and its budget are known and before either alignment
  method: the read's bases that every alignment must cover (the read without the free ends the window allows) are
  scanned for k-mers that occur in the gene window. If fewer occur than any alignment of score under the budget
  leaves intact, the candidate fails without WFA2; its taxon stays a failed candidate (the `ZF` tag and the
  profiler's failed-candidate counts are unchanged). The bound is the q-gram lemma with WFA2's penalties: an
  operation of score s touches at most s / c k-mers, c the smallest score per touched k-mer over mismatches
  (4 per k k-mers), insertion runs (6 + 2g per k − 1 + g) and deletion runs (6 + 2g per k − 1); with k ≤ 8 that is
  k/4 per score unit. A read k-mer with a base other than ACGT counts as present (WFA2 matches equal characters).
- k is 6 for windows of up to 512 bases and 7 above (`AlignmentScreen::KFor`). Measured on 19,309 random candidates of
  2-25% divergence with the budgets and free ends AlignAnchor gives them (test `RefusesOnlyWhatWFA2Fails`): WFA2
  failed 10,247; k = 6 refused 47% of the short reads' failures and 36% of the long windows', k = 7 36% and 39%,
  k = 8 (the first choice) 19% and 24%, k = 5 40% and 15%. Nothing refused was aligned by WFA2, for any k. The q-gram
  bound is weak near the ANI floor because random edits leave more intact k-mers than the lemma's worst case
  (evenly spread edits): it refuses candidates below roughly 83% identity, not the 88-90% ones.
- The counting is a generation-stamped table of 4^k 16-bit stamps per handler (8 KB at k = 6, 32 KB at k = 7),
  one pass over the window and one over the read, with early exits either way; `--no_alignment_screen` turns it
  off (a dev option, for measuring).
- Tests (`tests/test_AlignmentScreen.cpp`, 4 tests): the bound's values; a read of its own gene passes and random reads
  are refused; the WFA2 check above; and through the handler, 4,460 random anchors give the same outcome, score,
  CIGAR and start with the screen on and off, from the anchor's exact matches and as whole windows (278 refused).
  All 338 unit tests pass (one skipped, as before).
- Whole runs, outputs: SAM text and every profile output identical to `c94bd0e` on the 500k-pair sample (456,293
  records) and the 5M-pair sample (4,663,232 records, compared as sets: at six threads the record order differs
  between any two runs); the PacBio SAM identical as a set of records (36,359) with the screen and without it, the
  profile outputs identical.
- Whole runs, counts and time: on the simulated samples nearly every candidate is from a database species and
  aligns, so the screen refuses 5% of the short-read candidates (71,332 of 1,423,103; 720,480 of 12,782,429 on the
  deep sample) and the alignment stage's time does not change within this machine's noise. On the PacBio sample it
  refuses 20% (27,043 of 134,051), the ones WFA2 would have run through the whole budget of a gene-long window:
  aligning 2.85 s → 1.81 s against the same binary with the screen off (−36%), 3.27 s on `c94bd0e`. What it saves on
  the real samples of this report (44.9M fragments and 430k long reads failing everywhere) is for the next cluster
  run to show; the counts line reports it directly.

### 3. Timers in the profiling stage

Wall-clock timers around the profiler's steps (`Profiler::m_bm_read`, `m_bm_evidence`, the renamed SNP timer;
`MicrobialProfile::m_bm_em`, `m_bm_distances`; scoring and writing in `ProfileSample`) are appended to the sample's
`misc/<sample>_runtime.tsv` after the alignment stage's rows (`WriteProfilingTimes`), so `measure_performance.sh`
collects them into `stages.tsv` without a change, and printed in one line:

```
Profiling sample s took 6.1s: reading the SAM 4.6s, record evidence and sample context 1.0s, read EM 0.7s,
congener distances 0.2s, SNP post-processing 0.2s, scoring 0.2s, writing the profile 0.1s
```

On the 5M-pair sample at six threads (4.66M records), reading and parsing the SAM, with the taxa adding their
records, is 75% of the stage (4.6 of 6.1 s); the read EM, which the review of 2026-10-03 found at 7.5 s, is 0.7 s
after its rewrite, the distances 0.2 s, the rest under 0.3 s each. So the part to parallelise further is the SAM's
reading (`ProfileSamParallel`: one thread reads and cuts the stream, below eight threads it also decompresses; the
chunks are parsed in waves on all threads, and what depends on record order runs in file order). Whether the same
holds for the 36 s at 32 threads on the r226 run, where the SAM holds 45M unmapped records and the taxa are 143,614,
the next `measure_performance.sh` run on the cluster will show in its `stages.tsv`; that is the measurement item 3
asked for, and the parallelisation follows it. `--profile_ahead` for maps is unchanged and untested here (one sample).

Also timed: the per-taxon statistics files written after the profiling stage ("Taxon statistics files took"), one
file per taxon with reads, a candidate for the untimed seconds on a network file system.

### 4. Long reads: (a) the window's right end from the last link, (b) WFA2's AVX2 kernels, (c) the screen, (d) a budget from the best candidate

Measured as above (`scripts/chain2.sh`, `scripts/avx.sh`; `results/followup_longreads.txt`, `results/followup_avx2.txt`):
90 Mb of simulated PacBio (6,569 reads) and Nanopore (12,011 reads) at six threads, two rounds each where times are
compared, and the 500k-pair sample at one thread to show short reads untouched.

- **(a) On by default.** For chains with indels (long reads, `AnchoredAligner::IndelsAllowed`), `AlignAnchor` places the
  read's end on the gene by the last link's diagonal (`AlignmentOrientation::Update` with two diagonals; one diagonal
  gives the window it gave before), so the bases past the gene's end are the free ones by the end's own diagonal.
  Before, the first link's diagonal placed the end; where the indels had shifted it by more than the 9-base dovetail,
  the chain's right flank did not fit and the read went to the whole-window alignment, which forced the overhang
  into the alignment. Test `LongReadsThroughTheirChain`: 56 of 60 reads through their chain (51 before; the 4 left are
  chains the anchored aligner does not handle for other reasons). PacBio: whole-window alignments 100 → 38 of 134,051
  candidates, 14 of 36,359 records changed, the profile
  identical. Nanopore: 1,605 of 35,747 records changed (4.5%; its 4,001 whole-window alignments remain, chains the
  anchored aligner does not handle for other reasons), the profile's values shift in 244 lines. Short reads: SAM and
  outputs identical. The change is the one round 3 described; it makes the long-read models' inputs differ slightly
  for Nanopore, so the next retrain absorbs it.
- **(b) WFA2's AVX2 extend kernels: no measurable gain on long reads either; left off.** The run-time dispatch of
  `1fced40` (protal's own AVX2 kernels, the same alignments) built on the working tree: PacBio aligning 2.05 s against
  2.14 s mean over three alternated rounds (−4%, within this machine's ±10%), user CPU 13.7 against 14.2 s; short reads
  within noise; SAMs identical. The report of 2026-10-02 found the same on short-read flanks. The extend kernels are
  a third of WFA2's instructions and the wavefronts of a 10% budget are not long enough for SIMD to pay.
- **(c)** is item 2 above: on PacBio it refuses 20% of the candidates and the alignment stage takes a third less; on
  Nanopore it refuses none (the reads' error rate puts every candidate near the floor, where the bound is weak).
- **(d) `--long_read_budget E`, opt-in, off by default.** Once a candidate of a read's gene has aligned, the gene's other
  candidates get a WFA2 budget of the best so far plus E edits (as mismatches, 4 each); one that would cost more fails
  there and counts as a failed candidate (`LongReadAligner::operator()`, the cap through `Align` into `AlignAnchor`).
  PacBio with E = 20: aligning 2.34 s against 2.63-3.18 s (−10-25%), 59,086 alignments made instead of 86,889; but
  28,759 of 36,359 records change (79%): the segments' MAPQ (from the best against the second), the `ZA` alternatives
  and the `ZF` failed candidates, and 2,600 fewer gene hits take their read's consensus taxon, so the profile changes
  (240 lines). E = 50 changes 22,060 records; its time here is within noise of none. Nanopore: E = 50 −15%, 17,564
  records changed. The saving is real but modest at this budget's size, and a cut candidate is not the far relative
  one assumes: a longer window with more penalties can outrank a shorter one by bitscore, so the cap can change which
  hit is best. It stays off until a benchmark with retrained models says the gain is worth the different MAPQ
  distribution; the option is there for that run (`build_gtdb_database.py` would pass it through `PROTAL_ARGS`-style
  settings, not done here).

## The second cluster run: `50ecb4b` (items 1-4) on the same node, database and samples

SLURM job 23937247 on `q512n10`, 2026-10-04 20:13, 32 threads, three runs per sample as before; the database and the
reads were in the page cache from the first run (no cold load). Files in [`results_v2/`](results_v2/).

| | pe1 before | pe1 now | pb2 before | pb2 now |
|---|---|---|---|---|
| wall (median) | 131.5 s | **106.7 s** (−19%) | 103.7 s | **31.9 s** (−69%) |
| user CPU | 2,236 s | 1,485 s (−34%) | 2,775 s | 489 s (−82%) |
| instructions | 26.0 T | 12.6 T (−51%) | 39.3 T | 4.80 T (−88%) |
| IPC / cache misses | 3.17 / 26.0 G | 2.31 / 25.9 G | 3.85 / 74.7 G | 2.67 / 8.5 G |
| aligning | 61.9 s | 38.7 s (−37%) | 83.7 s | 11.8 s (−86%) |
| – alignment handler per thread | 27.2 s | 3.9 s (−86%) | 83.1 s | 11.4 s (−86%) |
| – seeding per thread (share of aligning) | 20.8 s (34%) | 20.6 s (53%) | 4.5 s (5%) | 4.5 s (38%) |
| profiling | 36.0 s | 36.3 s | 3.4 s | 3.4 s |
| – reading the SAM | – | 33.8 s (94%) | – | 3.2 s |
| taxon statistics files (newly timed) | – | 14.9 s | – | 1.6 s |
| untimed (Run protal minus the stages) | 22-26 s | 8.6 s | 9 s | 7.5 s |
| max RSS | 38.0 GB | 38.1 GB | 38.0 GB | 38.0 GB |

The counts line (`runs.tsv`):

| | pe1 | pb2 |
|---|---|---|
| reads | 50,607,792 pairs | 497,656 |
| with an anchor | 79,531,387 mates (of 101.2M) | 497,655 reads |
| candidate alignments tried | 175,963,654 | 9,679,514 |
| refused by the k-mer screen | 159,137,128 (**90.4%**) | 9,254,056 (**95.6%**) |
| WFA2 ran, from the anchor / whole window | 16,764,300 / 62,226 | 424,449 / 1,009 |
| alignments made (within the ANI floor) | 10,114,232 | 268,742 |
| records written | 2,221,613 | 103,394 |
| unmapped records (fragments that failed everywhere) | 44,884,669 (as before) | 430,086 (430,087 before) |

What this says:

1. **The screen was the right lever, and it is exact at scale.** Nine in ten short-read candidates and 96 of 100
   long-read candidates on these real samples cannot reach 90% identity, and the screen refuses them for a pass over
   the read and the window. The alignment handler's time fell by 86% for both read types; per read pair the
   alignment stage now costs 24 µs of CPU (39 before), per HiFi read 0.76 ms (5.4 ms). The outputs are the same by
   construction, and the run bears it out: the paired-end run writes the same 44,884,669 unmapped records, the same
   233,255 one-sure-mate fragments and the same 66,966 cross-gene pairs as before. The HiFi run differs in one read
   (430,086 unmapped records against 430,087) and two settled gene hits (31,260 against 31,258): item 4a's placing of
   the window's right end, which gave one read an alignment it had lost to a whole-window alignment. Instructions per
   run fell by half and by 88%; the IPC fell with them, since what is left is the memory-bound seeding.
2. **The paired-end run is now bound by two things the screen does not touch: the SAM reading in profiling (33.8 s,
   94% of the stage, 32% of the run) and the seeding (20.6 of 38.7 s per thread of aligning, 53%).** The HiFi run
   has become short enough that its fixed costs are half of it: index load 3.2 s, genome preload 3.2 s, start-up ~7.5 s,
   taxon statistics files 1.6 s, profiling 3.4 s, against 11.8 s of aligning.
3. **The profiling stage's split answers item 3's question: it is the SAM reading, 33.8 of 36.3 s**, with the read EM
   at 1.3 s and everything else under 0.3 s. The SAM holds 47.1M records, 95% of them the unmapped records of
   fragments that failed everywhere (each with its `ZF` tag of failed taxa, the evidence the profiler counts per
   taxon). `ProfileSamParallel` reads and cuts the stream on one thread, parses the chunks on all 32 and adds the
   records per taxon in parallel; 190 MB/s of SAM text through a pipeline with that much parallelism points at the
   single reading thread or at what is serial per wave. The experiment below locates it.
4. **The per-taxon statistics files cost 14.9 s on the paired-end run (14% of it)**, one `misc/<taxon>.statistics.tsv`
   per taxon with reads, written after profiling: thousands of small files on NFS. They were most of the untimed
   time of the first run. One table instead of a file per taxon (an output-format change, so a website update), or
   writing them only for the taxa the profile reports, removes almost all of it. The untimed rest, 7.5-8.6 s, is the
   start-up (taxonomy, priors, neighbours, suspect copies, conservation, model) and the teardown, now worth its own
   timers (item 6).
5. **The next levers, in order**: the SAM reading of profiling (−30 s of 107, once located), the taxon statistics
   files (−14 s), the seeding at GTDB scale (item 5: 53% of aligning; the k-mer lookups into the 27 GB index), the
   start-up (item 6), and for HiFi runs the fixed costs. `--long_read_budget` is not needed for this sample: after the
   screen the long-read alignment is 11.8 s of a 32 s run.

### The SAM reading, located and fixed (`scripts/unmapped.sh`, `unmapped2.sh`, `verify_fold.sh`; `results/followup_profiling_unmapped.txt`)

The cluster's paired-end SAM has 47.1M records, 95% of them unmapped records with `ZF` tags naming the taxa the
fragment seeded on and failed against. In WSL, the deep sample's SAM (4.66M records, 393k unmapped) with 45M synthetic
unmapped records appended, `--profile_only` at six threads:

| SAM | records | failed taxa named per record | reading the SAM, 6 threads | 1 thread |
|---|---|---|---|---|
| deep sample as is | 4.66M | few (760 taxa in the database) | 3.6 s | – |
| + 45M unmapped, `ZF` over 760 taxids | 49.7M | 1-4 | 6.7 s | 18.1 s |
| + 45M unmapped, `ZF` over 140,000 taxids (as r226) | 49.7M | 1-8 | **50.5 s** | **21.2 s** |

The record count alone costs little (3 s for 45M records at six threads). What the GTDB-sized database adds is the
number of distinct failed taxa: a 1 MB chunk of unmapped records names tens of thousands of them, and for each the
chunk's `RecordEvidenceCollector` made a full `RecordEvidence` (a `robin_map` and four vectors, ~300 bytes) in an
`unordered_map`, which the main thread merged chunk by chunk into the profile's collector between waves: 9,000 chunks
times ~30,000 taxa, serial. Six threads were slower than one. A first repair, a light hash map of counts per chunk
folded into the records at the end, cut the per-chunk records but kept 200M serial hash-map increments and merges
(41 s at six threads, 38 s at one, under load). The fix in the commit: the failed candidates are counted in vectors
indexed by taxid (`FailedCandidateCounts` in `SamHandler.h`, grown as needed; protal's taxids are dense internal
ids), merged by vector addition and folded into the taxa's records once, before the counts are read
(`RecordEvidenceCollector::FoldFailedCandidates`, called at the start of `ApplyRecordEvidence`). Per chunk that is
574 KB at r226 size, 147 MB per wave of 128 chunks at 32 threads.

| | before (`50ecb4b`) | after |
|---|---|---|
| 140k-taxid SAM, 6 threads: profiling / reading the SAM | 56.8 s whole run | 11.7 s; 10.5 / 8.5 s |
| 140k-taxid SAM, 1 thread | – | 23.2 s; 22.2 / 18.6 s |
| deep sample, 6 threads | 8.1 s whole run | 6.8 s; 5.7 / 4.1 s |
| deep sample, 1 thread | – | 13.6 s; 12.6 / 9.4 s |
| profile outputs | | identical to `50ecb4b` on both SAMs, and the same at 1 and 6 threads |

(The whole-run times of this round are 1.3-2× those of the earlier rounds: other work on the machine; compare within
the round.) The reading now scales from one to six threads as the deep sample's does (2.2×). The 28 profiling,
evidence and screen tests pass. On the cluster this should take most of the 33.8 s out of the paired-end run's
profiling stage; the next `measure_performance.sh` run shows how much.

### The taxon statistics files behind a flag, and the start-up timed (`scripts/check_startup.sh`)

- The per-taxon `misc/<taxon>.statistics.tsv` files (one per taxon with reads: coverage, reads, ANI, MAPQ and whether
  it is reported, per sample) are written only with `--taxon_statistics`. They cost 14.9 s of the 107 s paired-end
  run on NFS; the per-sample profile files hold the same. The sample's own `.profile.log`, `.profile.gene.log` and
  `.profile.genes.log` stay: they are three files per sample and take 0.2 s together ("writing the profile" above),
  and `.profile.truth_annotated` is written only with a truth file. On the 500k-pair sample, 306 files with the flag
  and none without; every other output identical (the 6-thread SAMs identical as sets of records). The website's
  output page lists the statistics files and needs the flag added.
- Timers for what was untimed: "Loading the gene tables took" (the `ProtalDB` construction: `reference.map`, 14.5M
  lines at r226 size, parsed on all threads), "Loading the taxonomy, models and tables took" (taxonomy, the models'
  XML, gene conservation, suspect copies, species priors, gene neighbours), and "Freeing memory took", printed after
  "Run protal took" by an object destroyed after the database (`TeardownTimer`). Locally these are 0.1, 0.13 and
  0.005 s; the next cluster run says how the ~8 s of start-up at r226 size divide.

## The third cluster run: `7693f28` (the SAM-read fix, the statistics flag, the start-up timers)

SLURM job 23941225 on `q512n10`, 2026-10-04 21:42, 32 threads, three runs per sample, page cache warm. Files in
[`results_v3/`](results_v3/).

| | pe1 first run | pe1 second | pe1 third | pb2 first | pb2 second | pb2 third |
|---|---|---|---|---|---|---|
| wall (median) | 131.5 s | 106.7 s | **60.9 s** | 103.7 s | 31.9 s | **27.5 s** |
| user CPU | 2,236 s | 1,485 s | 1,419 s | 2,775 s | 489 s | 486 s |
| aligning | 61.9 s | 38.7 s | 38.3 s | 83.7 s | 11.8 s | 11.8 s |
| profiling | 36.0 s | 36.3 s | **6.4 s** | 3.4 s | 3.4 s | **0.5 s** |
| – reading the SAM | – | 33.8 s | 4.2 s | – | 3.2 s | 0.3 s |
| taxon statistics files | – | 14.9 s | 0 (off) | – | 1.6 s | 0 (off) |
| loading the gene tables | – | – | 6.0 s | – | – | 6.1 s |
| taxonomy, models and tables | – | – | 1.4 s | – | – | 1.3 s |
| preload genomes + index load | 6.4 s | 6.4 s | 6.4 s | 6.4 s | 6.4 s | 6.4 s |
| freeing memory | – | – | 0.16 s | – | – | 0.14 s |
| untimed | 22-26 s | 8.6 s | **0.7 s** | 9 s | 7.5 s | **0.2 s** |

The counts are identical to the second run's, as they must be (nothing in this commit touches the alignment), and the
SAMs hold the same records.

1. **The profiling stage went from 36.3 to 6.4 s** on the paired-end sample, the SAM reading from 33.8 to 4.2 s: the
   dense counting of failed candidates did at scale what it did on the synthetic SAM. What is left of the stage is the
   reading at 4.2 s (47M records on 32 threads, now in line with the local scaling), the record evidence 1.6 s (the EM
   1.2 s) and the writing 0.5 s. The HiFi stage is 0.5 s.
2. **Every second of a run is now timed.** The start-up divides into loading the gene tables 6.0 s, the genome
   preload 3.2 s, the index 3.2 s, and the taxonomy, models and tables 1.4 s: 13.8 s of fixed cost per run, which is
   half of the 27.5 s HiFi run and a quarter of the paired-end one. Freeing the memory takes 0.15 s, so the teardown
   is not it, and nothing is left untimed (0.2-0.7 s).
3. **The gene tables are the one surprise: 6 s.** `GenomeLoader::LoadPositionMap` and `LoadUniqueKmers` parse
   `reference.map` and `unique_kmers.tsv` (14.5M lines each at r226 size) in parallel chunks, but add the rows to the
   genomes in file order on one thread: per row a genome lookup, a duplicate check and a gene insertion, twice over.
   Adding per genome in parallel (the rows of a genome are contiguous in the file), or a binary table of the genes in
   the database file, would bring this to about a second. That is the next start-up item; the other 7.8 s are the
   loads of the sequences and the index, already parallel.
4. **Where the paired-end run's 61 s go now**: aligning 38.3 s (63%: seeding 20.6 s per thread, extending anchors
   3.6, sorting seeds 3.2, the alignment handler 3.9, output 1.8, reader 1.0), the fixed costs 13.8 s (23%), profiling
   6.4 s (11%), the SAM header 1.2 s. The seeding is the lever that is left for throughput; the gene tables for the
   fixed costs. The HiFi run at 27.5 s is half fixed costs, and the rest is its 11.8 s of aligning.
5. **Over the day**: the paired-end run from 131 to 61 s (−54%), the HiFi run from 104 to 27 s (−74%), with the same
   alignments and profiles but for item 4a's one HiFi read, on 32 of the node's 84 cores. The alignment stage kept all
   its threads busy in every run, so 84 threads would shorten it further; the fixed costs would not move.

### The gene tables loaded in parallel (`scripts/gene_table_check.sh`, `build_head_tests.sh`; `results/followup_gene_tables.txt`)

`reference.map` and `unique_kmers.tsv` (14.5M rows each at r226 size) were parsed in parallel chunks but their rows
were added to the genomes on one thread, and a database member was decompressed on one thread: 6.0 s per run on the
cluster. Two changes, outputs unchanged (the same genes, lengths, offsets and unique counts; the 10 gene-table tests
pass, with a new one whose genomes are spread over the whole file):

- `gene_table::AddPieceByGenome`: a piece's rows grouped by genome (a genome's rows are contiguous in protal's tables,
  so a group is a run, or one per chunk where a chunk cuts it; a table that lists a genome in several places still
  gives one group), the genomes made on the main thread, then each genome's rows added by one thread; the earliest
  problem by line is the one reported, as before (test `TheFirstProblemInTheFileIsReportedWithItsLine`).
- A compressed table (a database member or a `.zst` file) with several threads is decompressed into memory on all of
  them (`zstd::ParallelRead`, the reader the genome preload uses) and parsed from there; a raw file keeps streaming,
  which from the page cache is faster than a copy.

The bench (`GeneTables.BenchLoadOfLargeTables`, `PROTAL_GENE_TABLE_TAXA=30000`: 5.04M rows per table, a fifth of r226,
on this laptop's six vCPUs):

| both tables, 6 threads | `c76e838` | now |
|---|---|---|
| members of a single-file database (64 MB frames), the path a run takes | 709 ms | **417 ms** (−41%) |
| raw files on disk | 496 ms | 421 ms |
| one thread (either way) | 1.15-1.18 s | 1.24-1.27 s |

What is left at six threads is the parallel parse and adds at about 40 ns per row, and the serial grouping of the
rows (a few ns per row). On the cluster's 32 threads the 6.0 s should become 1-2 s; the next run tells. **A binary
table is not needed for that**: the parse it would save is parallel now, the members are 150 MB compressed, and a
new member kind would touch the build, the packer and `--unpack_db` for perhaps half a second more.

## The fourth cluster run (`8c37dd6`), the binary gene table and the concurrent start-up

The fourth run (SLURM job 23943078, [`results_v4/`](results_v4/)): paired-end 60.9 → 58.5 s, HiFi 27.5 → 24.6 s; the gene
tables 6.0 → 3.4 s, not the 1-2 s expected. The bench at r226's species count (`PROTAL_GENE_TABLE_TAXA=143614`, 168 genes each: 24.1M rows per table, 1.7× r226's 14.5M; with the split the loads now print; [`results/followup_gene_table_bin.txt`](results/followup_gene_table_bin.txt))
showed why, on six threads from a single-file database:

| part (both tables) | 1 thread | 6 threads |
|---|---|---|
| reading and decompressing | 1.08 s | 0.39 s |
| parsing the text | 2.44 s | 0.73 s |
| grouping by genome, making the genomes (serial) | 0.18 s | 0.19 s |
| adding the genes | 0.81 s | 0.33 s |
| whole | 4.63 s | 1.76 s |

Parsing is the largest part and scales about 3× on six threads; with the serial grouping that fits 3.4 s on the cluster's
32 slower cores. Two more findings: the index load read and decompressed the whole `reference.map` a second time on one
thread for the reference fingerprint (`ReferenceFingerprint::Of`, inside "Load Index took"), and the start-up's
single-threaded parts ran one after another.

**The binary gene table** ([`GeneTableFile.h`](../../../src/SequenceUtils/GeneTableFile.h)): `gene_table.bin`, a member of
`database.protal` beside the two text tables, holding per genome its genes (start byte, id, length, the four unique k-mer
counts) and the fingerprint of `reference.map`. The packer writes it from the text tables it packs, read and checked as a
run reads them (`--build`; `--compress_db` on a folder, and on a single-file database without a current one, rewriting it
once with the other members' frames copied); `--unpack_db` leaves it out. A run from `database.protal` loads it on all
threads with no parsing and no grouping, and takes the fingerprint from it; a table whose recorded sizes of
`reference.map`, `unique_kmers.tsv` and `reference.fna` are not the database's, or that is corrupt, is not used (a note
says why) and the text tables are read. At r226 size (14.5M genes) it is ~465 MB (an estimated ~150 MB compressed, under 1% of the file).

| the bench (24.1M genes), six threads, from `database.protal` | text tables | `gene_table.bin` |
|---|---|---|
| loading the genes and their unique k-mer counts | 1.76 s | **0.43 s** (reading 0.18, genomes 0.02, adding 0.22) |
| the fingerprint for the index check | `reference.map` read again on one thread | from the table |

The text path also hashes `reference.map` while it parses it, so a folder database does not read it twice either.

One trap found on the way: the genome map's iteration order follows the order genomes are made in, and some outputs
break ties by it (a run whose map was built in taxid order changed one value of `.profile.gene.log` and 13 PacBio
alignments). The table therefore keeps `reference.map`'s genome order and its loader makes the genomes in that order, with
no `reserve`; a test checks that the text and binary loaders iterate the map identically (`TheGenomeMapIsBuiltInTheSameOrder`).
Where outputs depend on that order is worth finding; it is outside this change.

Follow-up 2026-10-05 ([order independence](../2026-10-05-order-independence/README.md)): the genome map iterates in
taxid order whatever the insertion order (identity hash, dense taxids), and no output depended on it; the varying
`.profile.gene.log` value was a double sum in SAM record order, now exact.

**Concurrent start-up** (the user's request; `--sequential_load` keeps the old order for comparison): after the gene
tables, the index loads on a thread of its own beside the genome preload, and the single-threaded tables each on one
(taxonomy, the models, gene conservation, suspect copies, species priors, the gene neighbours' table; the neighbours bind
their lineages once the taxonomy is in). The index's field widths are taken from the genes before the preload changes
them; the tables set their own members of the genome loader, never its genes. Messages are buffered and printed in the
order they always had; only the timers move.

**Outputs.** Every check identical to the reference build `c94bd0e` but for what earlier commits already changed:
paired-end SAM and profile from the single file with `gene_table.bin`, sequential and concurrent, and from the unpacked
folder; PacBio identical between the text and binary tables and between the two loading orders (its difference from
`c94bd0e` is item 4a's). One value of `.profile.gene.log` varies in its last digit between repeated runs of the same build
in either loading order (against the first sequential run: one of three other sequential runs, three of four concurrent ones), while profiling one SAM four times gives the
same output: it depends on the SAM's record order, which multi-threaded alignment varies. That predates this work.
All 345 unit tests pass (two skipped: the bench and the one skipped before), including 7 new `GeneTableFile` tests.

**Time, locally**: on this database (84k genes, a 3.7 GB index) the start-up is under a second and the two loading orders
are within the machine's noise (whole runs 3.8-5.8 s either way); the concurrent load can only show at r226 size, where
the parts it overlaps are the preload (3.2 s), the index (3.2 s) and the tables (1.4 s), at most ~4.6 s together.

**What the next cluster run needs**: add the table to the r226 database once (it rewrites the 27 GB file through
`database.protal.partial`, so it needs that much free space beside it), then the script twice:

```bash
protal --compress_db --db /hpc-home/hildebra/DB/protal/protal0.7.3_r226_v10/protal_db/database.protal -t 32 --compress_level 9
bash scripts/measure_performance.sh Perf0.7.5_v5 $DB pe:...:... pb:...
PROTAL_ARGS="--sequential_load" bash scripts/measure_performance.sh Perf0.7.5_v5_sequential $DB pe:...:... pb:...
```

## The fifth cluster run (2026-10-05): `gene_table.bin` and the concurrent start-up at r226 size

SLURM job 23960603 on `q512n10`, protal at `735617c` (which holds `48e8cd0`), 32 threads, on a new database,
`protal0.7.5_r226_v11`, built by the 0.7.5 pipeline and so packed with `gene_table.bin`. It is the same reference as v10:
the same 2,903,557,366 index entries, 5,589 suspect copies, 143,614 species priors and 731,664 gene-neighbour rules.
Files in [`results_v5/`](results_v5/).

| | pe1 fourth run (v10) | pe1 fifth run (v11) | pb2 fourth | pb2 fifth |
|---|---|---|---|---|
| wall (median of 3) | 58.5 s | **51.4 s** | 24.6 s | **17.7 s** |
| loading the gene tables | 3.37 s | **0.46 s** (`gene_table.bin`: reading 0.08, genomes 0.07, adding 0.31) | 3.39 s | 0.46 s |
| genome preload | 3.21 s | 3.80 s, beside the index and the tables | 3.22 s | 3.73 s |
| taxonomy, models and tables | 1.39 s | done within the preload (their timer, 3.84 s, spans it) | 1.30 s | idem |
| index load | 3.23 s, after the others | 2.54 s, beside them (no second read of `reference.map`) | 3.22 s | 2.54 s |
| start-up before aligning | 11.2 s | **4.3 s** | 11.1 s | 4.2 s |
| aligning / profiling | 38.4 / 6.2 s | 37.9 / 6.4 s | 11.9 / 0.53 s | 11.8 / 0.49 s |
| max RSS | 38.2 GB | 36.2 GB | 38.1 GB | 36.2 GB |

The counts are the same as in every run since `50ecb4b`. The 1.9 GB less peak memory is the text tables' in-memory copies,
which the binary table does not need. The preload took 0.6 s longer beside the index decode, which shares the cores'
memory bandwidth; the overlap still saves 4 s.

**The gene count.** The table holds 14,527,738 genes: 101 per species of the 168 markers, as a species carries its own
domain's markers. Earlier reports extrapolated 16.6M from the benchmark world, and the sections above assumed 24M (all 168
per species). So the bench's 24.1M rows per table overstated the per-table times by about 1.7×; the comments that said
24M are corrected. The count is the database's, not a bug.

**Where the time goes now.** Paired-end, 51.4 s: aligning 37.9 s (74%; seeding 20.4 s per thread, 54% of it), profiling
6.4 s (13%), start-up 4.3 s (8%), the SAM header and file 1.6 s, the rest under a second. HiFi, 17.7 s: aligning 11.8 s
(67%), start-up 4.2 s (24%). Over the two days the paired-end run went from 131 s to 51 s (−61%) and the HiFi run from
104 s to 18 s (−83%), on the same 32 threads, with the same alignments and profiles but for item 4a's long-read windows and the order-independence fix of `3c2d56c`.

**The cold first run**, 222.6 s, read the 27 GB file from NFS for the first time: the index (174 s) and the preload
(47.6 s) streamed their members at the same time. It is not comparable with the first day's cold run, whose page cache
held part of the old database.

**What is left, by size**: the seeding at GTDB scale (20.4 s per thread of the paired-end run, memory latency in the
k-mer lookups; item 5); the preload, now the start-up's critical path (3.8 s; the genes stored 2-bit packed in the
database would skip the packing, a format change); the SAM header and file (1.6 s on the paired-end run, up from 0.5 s,
to be looked at); and the cold start (item 7: parallel reads or a copy to node-local disk).

## The seeding: how a lookup works, and the options (2026-10-05)

**How it works now** (`ChainAnchorFinder::FindSeeds`, `KmerLookupSM::GetFromLookup`, `Seedmap::GetPacked`). Every
closed syncmer of a read is a k-mer whose 15-base core is looked up in the key map (3.2 GB): its control block, then the
core's slot range, one or two dependent cache misses. The core's values are a block in the 27 GB value array: first one
32-bit flex cell per entry (the k-mer's other 16 bases), then the 42-bit entries. The lookup scores every flex cell of
the block against the read's (`Similarity`: matching bases, a popcount), keeps the entries of the best score, and drops
the lookup if more than `--max_key_ubiquity` (256) entries tie. A read's lookups are taken smallest block first; a short
read stops once it has 128 seeds and enough lookups, a long read takes all. The prefetch of `87b857e` fetches key-map
blocks 8 k-mers ahead, value blocks 16 lookups ahead, and the next read's key-map blocks while this read aligns.

At r226 the seeding is 20.4 s per thread of the paired-end run's 37.9 s of aligning (6.4 µs of one thread per mate), its
IPC low. What it waits on is the open question: the key map's misses (latency; many small blocks), or the scans of large
blocks (the cores of conserved genes, shared by thousands of species; a read hits a block in proportion to its size).

**Options, the measurement first:**

| # | option | exact | effort | pays when |
|---|---|---|---|---|
| 1 | **Count first** (done, below): k-mers looked up, cores in the index, blocks scanned, flex cells, blocks and cells by block size, seeds | yes | – | decides 3-8 |
| 2 | **More threads**: the stage is latency-bound and kept all 32 threads busy; the node has 84 cores. `THREADS=84` in the script | yes | none | always; likely the largest single gain |
| 3 | **The flex scan with AVX2**: one pass instead of two (score, then collect the ties), 8 cells per instruction (equal-base mask, popcount by nibble table); the scalar loop pushes every score into a vector | yes | small | cells per block in the hundreds |
| 4 | **Blocks sorted by flex cell at load** (a load-time permutation, no database change): a read whose 16 bases occur exactly (the common case for a species in the database) finds its ties by binary search, a contiguous run; only reads without an exact match scan the block. Emitted in the original entry order, so the same seeds | yes | medium | most cells in blocks ≥ 256 |
| 5 | **Huge pages for the index's arrays**: with 4 KB pages every lookup also misses the TLB (a page walk per miss); THP is `always` on the node, but whether the 30 GB are backed by 2 MB pages is not known (`AnonHugePages` in `/proc/<pid>/smaps` during a run); `madvise(MADV_HUGEPAGE)` on them as on the gene arena | yes | small | latency-bound |
| 6 | **More misses in flight per thread**: interleave the lookups of 2-4 reads (the prefetch report's in-process bench gave 12-15% on a 3.7 GB index; at 27 GB the latency is longer) | yes | medium | latency-bound |
| 7 | **Small blocks next to their key** (a load-time layout: blocks of 1-15 entries stored in or beside the control block, flex cells and entries interleaved): one miss instead of two or three | yes | medium-large, more memory | most blocks small |
| 8 | **Fewer lookups**: sparser syncmers, a lower seed target for short reads, or a cap for long reads | no: sensitivity | small, a benchmark | after 2-7 |

The local check of option 1, on the 765-species world (where the picture is not GTDB's): paired-end 22.8M k-mers, 8.5M
cores in the index, 9.3 flex cells per block, 64% of the cells in blocks of 16-255 entries; HiFi 14.6 per block. The r226
run will show the block sizes that matter.

## The SAM header (2026-10-05)

How it is built, and what is parallel:

- During the alignment each thread compresses its own records (zstd, `SamOutput::Write`); under one mutex the genes they
  name go into a `std::unordered_set` and the compressed bytes into the file. A `.sam.zst` is written with room for the
  header at its start: up to all of the database's genes at 8 bytes each, at most 1 MB (`sam_zstd::HeaderRoom`).
- After the alignment, on one thread ("Writing the SAM header and file took"): the set is copied and sorted
  (`SamOutput::Genes`); `GenomeLoader::WriteSamHeader` writes one `@SQ` line per gene (three hash-map lookups each); the
  header is compressed in one call; it is written into the room (`PlaceHead`). A header larger than the room makes
  `MoveRecordsBehind` copy all the records behind it, on one thread.
- So nothing of it is parallel. Locally it is negligible: 36,522 genes give a 127 KB compressed header (5M pairs), well
  inside the room. At r226 the room is capped at 1 MB, about 290,000 genes at that compression; a gut sample against
  14.5M genes may name more, and then the 1.6 s are mostly the copy of a SAM with 47M records. The run now prints which
  ("SAM header: N genes, B bytes, written as C (text built in T s); the records (R bytes) were copied behind it in X s, as
  it was larger than its room of M bytes", or "placed in its room"), and `measure_performance.sh` records the genes and
  the copy's time.
- If it is the copy: raise the room's cap (it is a hole in the file where the file system has them; at worst its zeros
  are written once), sized from the genes the run actually names. If it is the text: the per-gene lookups can be done in
  parallel chunks. Both exact; for the next run's numbers to decide.

## MSA and qcMSA in the measurement (2026-10-05)

The strain stage needs a species in two samples, so the script's single-sample runs never had MSAs to time. Now, with two
samples or more, it ends with a cohort run: all samples in one protal run from a map whose SAM column names the
per-sample runs' last SAMs, so protal profiles them again without aligning ("All alignments are present"), then builds
the strain MSAs and runs qcMSA (on by default; `QCMSA=0` without, `COHORT=0` no cohort, `COHORT_REPEATS` runs).
`cohort.tsv` records the profiling, the strain stage, building the MSAs and qcMSA separately (new timers in
`StrainWrapper2`: "Building the strain MSAs took", "qcMSA took", and "Strain-level MSAs of S species: R raw MSAs, F
filtered by qcMSA"). Each sample's profiling timers from the cohort go to `stages.tsv` as `cohort:<sample>`.

Locally (pe 500k pairs and HiFi 90 Mb from the 765-species world, six threads): 14 species, the MSAs built in 0.43 s, qcMSA
4.35 s (one python process per species, one after another): qcMSA is most of the strain stage, and the first thing to
parallelise there (its species are independent) once the r226 numbers are in.

## The sixth cluster run (2026-10-05): the seeding's counts, the SAM header, the strain stage

SLURM job 23962556, protal v0.7.6 at `959d46c` (with `f29595f`'s counters), the r226 v11 database, 32 threads, three runs
per sample and one cohort run. Files in [`results_v6/`](results_v6/). Wall times as in the fifth run: paired-end 50.8 s,
HiFi 18.0 s.

**The seeding** (the "seeding" lines):

| | paired-end (101.2M mates) | HiFi (497,656 reads) |
|---|---|---|
| k-mers looked up | 2,242,166,093 (22.2 per mate) | 473,852,030 (952 per read) |
| cores in the index | 99.2% | 99.2% |
| blocks scanned | 2,221,730,360: the stop at 128 seeds almost never comes | all |
| flex cells scanned | 156.4 G, **70.4 per block** | 34.9 G, 74.2 per block |
| cells in blocks of 16-255 / 256-4095 entries | 55.5% / 40.9% (none larger) | 54.1% / 42.6% |
| seeds | 4.37 G | 0.96 G |
| seeding time per thread | 20.5 s: ~295 ns per lookup, ~4 ns per flex cell | 4.6 s |

So the seeding is not only memory latency: every lookup scores ~70 flex cells in a scalar loop that also pushes each
score into a vector and then passes over it again for the ties, and 96% of the cells are in blocks of 16 to 4,095
entries. That favours option 3 (one AVX2 pass over the cells) for every lookup, and option 4 (blocks sorted by flex cell,
the ties found by binary search when the read's 16 bases occur exactly) for the 41% of cells in blocks of 256 or more;
how often a lookup has an exact flex match is the one count option 4 still needs. Option 2 (`THREADS=84`) stays the
cheapest test.

**The SAM header** is not the cause of the ~1 s: paired-end 80,168 genes, 1.9 MB of text built in 0.02 s, 371 KB
compressed, placed in its 1 MB room (no copy); HiFi 13,870 genes, 66 KB. "Writing the SAM header and file" still takes
1.10 s and 1.07 s. What is left in it is the seek table and closing the file; on NFS `close()` waits for the client to
write back and commit all the SAM's data written during the alignment, so the second is most likely that. Writing the
SAM to node-local disk, or closing it on a thread of its own while the next sample aligns, would show and remove it.

**The strain stage** (the cohort run, both samples from their SAMs, 40.1 s): profiling 6.2 s, strain-level MSAs 29.0 s
for 114 species (114 raw MSAs, 102 kept by qcMSA): building the MSAs 5.4 s and **qcMSA 23.3 s**, one species after
another (qcMSA a python process each, 0.2 s per species). The species are independent: run the MSAs and the qcMSA
processes over the threads and the stage should take a few seconds. The environment line's "qcmsa not found" was the
script's own check, which missed protal's last fallback (`scripts/qcmsa.py` beside a source checkout's `build/`); fixed.

## The AVX2 flex scan and the parallel strain stage (2026-10-05)

Checked locally (WSL, six vCPUs, the 765-species world; [`results/followup_flex_strain.txt`](results/followup_flex_strain.txt)).

**The flex scan** ([`FlexScan.h`](../../../src/Hash/FlexScan.h)): `KmerLookupSM::GetFromLookup` scored each cell of a
block in a scalar loop that pushed every score into a vector, and passed over it again for the ties. Now
`flex_scan::Score` scores 8 cells per step with AVX2 (the 8 cells as words shifted by the block's bit offset, the equal
bases by XOR, popcount by nibble table), writes the scores to a byte array and returns the best score and its count;
the ties are then taken as before, in entry order. It is chosen at run time where the CPU has AVX2;
`PROTAL_FLEX_SCAN=scalar` keeps the scalar loop. Exact:

- test `FlexScan.Avx2ScoresAsTheScalarScan`: 2,880 random blocks of 1-5,000 cells at every bit shift, keys matching
  cells exactly, partly or not at all; the same scores, best and count;
- runs: SAMs and profiles identical with the AVX2 and the scalar scan (paired-end 500k pairs, PacBio 90 Mb).

The bench (`PROTAL_FLEX_BENCH=1`): 3.60 → 1.07 ns per cell for blocks of 70 cells (r226's mean; 3.4×), 3.32 → 1.00 for
1,000, 3.42 → 1.94 for 16. If the scalar scan was most of the r226 lookup's ~295 ns (70 cells × ~3.6 ns), the seeding's
20.5 s per thread could fall by up to half; locally the blocks average 9 cells and the seeding is within noise either
way. The next cluster run shows it.

**The strain stage** (`StrainWrapper2`): species ran one after another. Now they run over the threads, each species on
one: its MSA, then qcMSA. The samples each species enters are found first on one thread (which caches the model's
scores in the taxa, so the workers only read the profiles); each worker writes its species' files; a species' messages
are held until it and every species before it are done, then printed in species order, and the species list is written
in its order. The qcmsa processes run with their numerical libraries on one thread each. "Building the strain MSAs took"
and "qcMSA took" are now sums over the species (they overlap); "Strain-level MSAs took" is the stage's wall time.

| cohort of pe 500k + PacBio 90 Mb, 14 species | 1 thread | 6 threads |
|---|---|---|
| strain-level MSAs (wall) | 7.1 s | **0.74 s** |
| building the MSAs (summed) | 0.36 s | 0.23 s |
| qcMSA (summed) | 6.7 s | 3.4 s (one library thread per process is cheaper per species) |
| strain outputs | | identical (113 files) |

**The log**: qcMSA's own output (its progress bars and report, per species) is no longer shown: it goes to a file that is
printed, its last 20 lines with the command, only if qcmsa fails; the stage says once which qcmsa it runs. The MSA
building's per-gene progress bar is gone too (it wrote to `std::cout` from every worker). The cohort's log shrank from
125 KB to 7.7 KB, and protal's own lines start lines again (qcmsa's bars ended without a newline, which hid
"Strain-level MSAs took" from a line-start search; the script's parser now finds a timer anywhere on a line).

At r226 (114 species, 29 s one after another, 32 threads) the stage should take a few seconds; the next cohort run shows.

## How it was run

On the cluster (the user's job; the paths are the cluster's):

```bash
sbatch --exclusive -c 32 --mem 100G --wrap "bash scripts/measure_performance.sh \
    /hpc-home/hildebra/dev/protalFastHit/runs/Perf0.7.5 \
    /hpc-home/hildebra/DB/protal/protal0.7.3_r226_v10/protal_db/database.protal \
    pe:.../PEARL_AGE/00_fastq/B10_R1_001.fastq.gz:.../B10_R2_001.fastq.gz \
    pb:.../m64319e_230622_130841.bc1002--bc1002.hifi_reads.fastq.gz"
```

The output folder was copied to `local/Perf0.7.5/` and from there, unchanged, into `results/`. The derived numbers
above are arithmetic on `runs.tsv`, `stages.tsv` and the "took" lines of the logs.
