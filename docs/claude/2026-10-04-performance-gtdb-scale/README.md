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
