# Multithreading: locks, hand-offs and what does not run in parallel

- **Date**: 2026-10-01.
- **Code**: branch `audit-fixes` at `1c11a00` (Version 0.7.1; `src/` is unchanged at `066a7e7`, the HEAD when
  this was written). Measured with two builds of it: the plain Release build of the 0.7.1 benchmark of the same
  day (`protal-base`), and the same commit with `scripts/mtaudit.patch` (`protal-instr`): wall-clock counters
  around every lock and hand-off on the alignment and profiling paths, printed as `MTAUDIT` lines when protal
  exits. The patch is for measuring only and is not meant for the tree.
- **Machine**: WSL2 Ubuntu 24.04 on an Intel Core Ultra 7 258V (4 fast and 4 low-power cores; WSL gets 6
  vCPUs), 23 GB, gcc 13. It is shared with other sessions. The scaling runs (18:08–18:18) had only this
  session's load (load average 1.6–4.4 before each run). From 18:19 other sessions ran a Python experiment and
  training-data simulations (load 6–13). The input benchmarks made after that are marked, and the ones at load
  10–13 are left out.
- **Data** (read only, from the 0.7.1 benchmark, `docs/claude/2026-10-01-v071-benchmark`): the 765-species
  database `~/bench071/V071/protal_db` (single file `database.protal`, with gene neighbours); the sample
  `rl150_p5000000_s_1` (5M read pairs 2x150, BGZF as the simulator writes it, 407 MB per file).
- **Runs**: whole protal runs (`-1 R1 -2 R2 --no_qcmsa`) at `-t` 1, 2, 3, 4 and 6 with `protal-instr`, and at
  6 twice with each build, niced (`scripts/scaling.sh`); `readbench`, the input path alone (`scripts/readbench.cpp`):
  protal's `ThreadedGzIstream` and `SeqReaderPE` on T threads that do nothing with the reads. Commands and builds
  are in `scripts/`, the numbers in `results/`.

## Summary

The locks are not what limits protal's threads. During alignment, threads wait on the reader lock for 1.3–1.5%
of their time at 6 threads and on the SAM writer for none. The inflating threads never starved the aligners. The
statistics merges and the critical sections on error paths do not matter.

What limits a run on many threads is the work that runs on **one** thread:

1. **Profiling a sample runs on one thread.** It takes 13–20 s of a 5M-pair run, 28–31% of the wall time at
   6 threads. At 32 threads it would be about 60% of the run. Parallelising it is the one change that would make
   a large difference: about 2x end-to-end for a single-sample job on 16 or more threads.
2. **Each FASTQ file is inflated by one thread.** The input path delivers at most 1.15M read pairs/s from
   standard gzip, and 1.4–1.7M from BGZF. Alignment here runs at 25–36k pairs/s per thread, so input becomes
   the limit at roughly 30–60 threads, on this database. A denser database (GTDB) costs more per pair and pushes
   that further out.
3. Smaller serial pieces are the index load (0.7–1 s at 6 threads), copying the SAM records behind the header
   (1.4–1.7 s at 6 threads), and the profiling of several samples, whose static schedule balances unevenly sized
   samples poorly.

## The locks and hand-offs

| Where | What it guards | How often | Measured cost |
|---|---|---|---|
| `critical(reader)`, `SeqReaderPE`/`SE` (`SeqReader.h`) | taking the next 32 records of each file from the inflated buffer (`TakeLines`: `memchr` and one copy) | once per 32 pairs | 1.3–1.5% of thread time waiting at 6 threads, 0.3% at 4; held 6% of the time at 6 threads (below) |
| `ThreadedGzStreambuf` mutex and condition variable (`ThreadedGzStream.h`) | 1 MB blocks between each file's inflating thread and the lock holder | once per MB | consumers waited ≤ 4 ms per run for inflated data; the inflating threads idle 90% of the time |
| `SamOutput::m_mutex` (`SamFile.h`) | `write()` of a compressed block and the gene set | once per 16 MB of SAM per thread (≈ 205 per run) | waiting ≤ 0.07 s per run; held 2.6–8.5 ms per flush; compression (10 s summed at 6 threads) is outside the lock |
| `critical(statistics)` (`Classify.h`, `Build.h`) | merging each thread's counters and timers | once per thread per sample | — |
| `critical(genome_loader)` (`GenomeLoader.h`) | loading a gene from `reference.fna`; with preloaded genomes (the default) only entered by `LoadGeneOMP` | profiling: once per gene per taxon per sample | negligible; see correctness note |
| `critical(err_out / errout / debug_out / invalid_align / print)` | printing inconsistent alignments and debug output | only on errors | — |
| `critical(print / create_dir)`, `RunStatus` mutex (`RunProtal.h`) | profiling samples in parallel: messages, output folders, failures | per sample | — |
| `critical(put)`, `omp single` in `PartitionedPass`, atomics (`Seedmap.h`, `Build.h`) | the index build | — | not re-measured (see "The index build") |

Nothing else synchronises threads on the alignment path. Each thread has its own k-mer handler, anchor finder,
lookup copy, WFA aligner and output buffer, and the index and genes are read only.

### The reader lock in numbers

156,121 entries per run, each taking 32 pairs (both files):

| `-t` | aligning (s) | speed-up | held (s) | held per entry | waiting (s) | waiting per entry | lock busy | waiting, share of thread time |
|---|---|---|---|---|---|---|---|---|
| 1 | 137.5 | 1.00 | 0.87 | 5.6 µs | 0.03 | 0.2 µs | 0.6% | 0.02% |
| 2 | 70.8 | 1.94 | 0.90 | 5.8 µs | 0.05 | 0.3 µs | 1.3% | 0.03% |
| 3 | 52.9 | 2.60 | 1.05 | 6.7 µs | 0.12 | 0.7 µs | 2.0% | 0.07% |
| 4 | 45.2 | 3.04 | 1.42 | 9.1 µs | 0.57 | 3.6 µs | 3.1% | 0.31% |
| 6 | 33.1 | 4.15 | 2.01 | 12.9 µs | 2.89 | 18.5 µs | 6.1% | 1.45% |
| 6 (repeat) | 31.2 | 4.41 | 1.82 | 11.7 µs | 2.41 | 15.5 µs | 5.8% | 1.29% |

The plain build aligned in 29.9 and 31.1 s at 6 threads, within the noise of the instrumented one, so the
counters cost little.

Two things stand out. Holding the lock takes twice as long per entry at 6 threads as at 1: the inflated block is
in another core's cache, and 6 aligners plus 2 inflating threads share 6 vCPUs, so a holder can be descheduled.
Also, at 6 threads the mean wait (15–19 µs) is longer than the hold, although the lock is busy only 6% of the
time. The hand-off itself (waking a sleeping waiter, a descheduled holder) costs more than the copy it guards.
This does not matter at 6 threads. It tells what to expect on many cores: the cost grows faster than the number
of threads.

Alignment loses 27–31% efficiency at 6 threads. The waits measured here account for at most 1.5 points of that.
The rest comes from the hardware: low-power cores beyond the fourth, shared cache and memory bandwidth, lower
clocks with all cores busy. A lock change cannot recover it.

### The input ceiling

`readbench` reads the 5M pairs through protal's input path on T threads that do nothing else, so the input path
alone sets the rate (`results/readbench*.tsv`):

| input | 1 thread | 2 | 3 | 4 | 6 |
|---|---|---|---|---|---|
| BGZF (simulator), pass 1, load 3–5 | 1.70M pairs/s | 1.48M | 1.25M | 1.29M | 0.85M |
| uncompressed, pass 1 | 1.78M | 2.97M | 2.86M | 3.26M | 2.30M |
| BGZF, load ≈ 6 | 1.47M | 1.46M | | | |
| standard gzip (`pigz -6`, one member, the zlib-ng path), load ≈ 6 | 1.15M | 1.17M | | | |

With gzip input, more readers bring nothing. Readers that do nothing else contend so hard that the rate drops.
The limit is the single inflating thread per file: about 450 MB/s per file for standard gzip at load 6, and
0.5–0.9 GB/s for BGZF in the alignment runs (`inflate_busy`). Uncompressed input goes up to 3–4.5M pairs/s.
With 6 readers it collapses again: 8 busy threads on 6 vCPUs.

Real FASTQ from sequencers is standard gzip, so its ceiling is about 1.15M pairs/s on this machine. Alignment
on this database runs at 36k pairs/s on one thread and 25–27k per thread at 6 threads. Input would therefore cap
a run at about 45 threads, and the reader lock's hand-offs may cost noticeably before that. On GTDB a pair costs
more to align, and the cap moves to more threads. This needs measuring on a many-core node with the real
database.

### The SAM writer

Each thread fills a 16 MB buffer, compresses it with zstd outside the lock (6.5–10.4 s summed per run), and
holds the lock only to append it and add its genes to the set. Waiting for that lock totalled 0.0002–0.07 s per
run. On a slow network file system, `write()` under the lock would show up first. A dedicated writer thread
would take it off the aligners, but nothing measured here calls for one.

## What runs on one thread

A 5M-pair run at 6 threads (`i6b`, 47.5 s):

| Stage | Time | Threads |
|---|---|---|
| index load | 0.84 s | all (15.2 s on one thread) |
| alignment | 31.2 s | all |
| copying the records behind the SAM header (`Writing the SAM header and file`) | 1.37 s | one |
| profiling | 13.3 s | one |

The serial part is 15.6 s, 33% of the run.

### Profiling a sample

`ProfileWrapper` profiles samples in parallel, one sample per thread. A single sample, as in a per-sample cluster
job, runs on one thread while the others idle. Its parts, from the counters in `ProfileSam`:

| part | `-t 1` run | `-t 6` run (`i6b`) |
|---|---|---|
| reading and parsing the SAM (zstd, `SamReader`, the `AlignmentPair` copies) | 7.6 s | 4.9 s |
| the reads' work (`BestOfGroup`, CIGAR, `NoteRecord`, gene neighbours, `AddSam`) | 8.3 s | 5.3 s |
| of it `MicrobialProfile::AddSam` (genes, variants, read ranges) | 5.4 s | 3.5 s |
| `ApplyRecordEvidence` and `PostProcessSNPs` | 1.0 s | 0.8 s |
| after `ProfileSam`: model scores, outputs | 2.6 s | 2.3 s |
| **profiling, total** | **19.6 s** | **13.3 s** |

The same single-threaded work took 13–20 s depending on which core and clock it got. The machine's timings vary
up to 2x.

Projected to 32 threads (the 6-thread efficiency of about 70% assumed): alignment ≈ 6 s, the serial parts ≈ 15 s,
so profiling is about 60% of the run. Going from 32 to 64 threads would save less than 2 s, because input caps
the alignment and profiling does not change. With profiling at about 3 s, the 32-thread run would drop from about
21 s to about 11 s.

A design that keeps the outputs byte-identical:

1. **Parse in parallel.** A `.sam.zst` is a seekable list of independent 1 MB frames (`SamFile.h`). Cut it at
   frame boundaries into ranges, one per thread. Each range starts at the first read whose name differs from the
   line before, and the previous range reads on to the end of its read. Each range yields its records' work
   items (taxon, gene, the record, its read and link numbers within the range). It also yields the
   order-independent integer counts: `NoteRecord`, the gene-neighbour credits, which stay within a read.
   BGZF `.sam.gz` splits the same way; plain SAM splits at byte offsets.
2. **Apply in order, by taxon.** Split the taxa over the threads. Each thread applies its taxa's items range
   after range, in file order. Every taxon then sees its reads in the same order as now, so the floating-point
   sums, read identities and variants come out the same. Read numbers become the range's offset plus the local
   number, and the offset is known once the earlier ranges are parsed.
3. Post-process, score and write each taxon in parallel.

The largest taxon holds 13% of the records here, which bounds step 2 at about 7.5x. Step 2 alone, with a serial
parser, would stop at the 5–7.6 s of parsing. Parsing could also get cheaper serially: `ReadSamGroups` copies
each record into a `std::optional<SamEntry>`, and `SamReader` splits each line into a vector of strings. This was
not measured in detail.

### Several samples

The per-sample loop in `ProfileWrapper` (`RunProtal.h`) is `#pragma omp parallel for` with the default static
schedule. With more samples than threads and samples of different depth, a thread can be left with several deep
ones. `schedule(dynamic, 1)` fixes that and leaves the outputs unchanged.

### The SAM's records behind its header

`SamOutput::Finish` copies the records file behind the header once the alignment is done (1.4–1.7 s at 6
threads, 2.1 s at 1), because the header lists only the genes the records name. Profiling could read the records
file while the final SAM is written, or the copy could overlap the next sample's alignment. That saves about 3%
at 6 threads and 7% at 32.

## Recommendations, by what they save

| # | Change | Saves | Effort | Outputs |
|---|---|---|---|---|
| 1 | Profile one sample on all threads (the design above) | 10–16 s per 5M-pair sample; about 2x end-to-end at 16 or more threads for single-sample jobs | medium–large: the range splitter, the ordered taxon shards, a merge of the integer counts, tests that the profiles are byte-identical | identical by construction |
| 1a | First step: post-process, score and write taxa in parallel, and shard `AddSam` by taxon behind the serial parser | about half of profiling | medium | identical |
| 2 | `schedule(dynamic, 1)` for the per-sample profiling loop | idle threads in multi-sample runs | one line | identical |
| 3 | Input for 40 or more threads: larger reader batches (32 → 128–256 pairs) to amortise the hand-offs; several inflating threads for BGZF input (its blocks are independent); for standard gzip, which cannot be split, a faster inflater (ISA-L, declined earlier for portability) or zstd FASTQ input | lifts the ceiling of 1.15–1.7M pairs/s; nothing below about 30 threads | the batch size is one constant; parallel BGZF is small; zstd input is medium | identical |
| 4 | Overlap the SAM records copy with profiling or the next sample | 1.4–2 s per sample | small–medium | identical |
| — | SAM writer, statistics merges, error-path sections | nothing measurable | — | — |

Items 1 and 3 should be checked on a many-core node with the GTDB database before more work goes into item 3.
There, a pair costs more to align and the profile has many more taxa, so the balance between alignment, input
and profiling will differ from this laptop's.

## Correctness notes (not speed)

- `Genome::GetGeneOMP` (`GenomeLoader.h`) checks a plain `bool m_is_loaded` outside `critical(genome_loader)`
  and sets it inside, which is double-checked locking without an atomic. With `--preload_genomes_off`, a thread
  may see the flag set before the genes it guards, and it is a data race ThreadSanitizer would report. It should
  be a `std::atomic<bool>` (release when set, acquire when read). With preloaded genomes (the default) nothing
  writes it during alignment.
- With `--preload_genomes_off`, one global `critical(genome_loader)` serialises every gene and genome load, disk
  reads included, across all genomes. A lock per genome would let different genomes load in parallel. This only
  matters for that option.
- `Taxon::AddSam` and `AddHit` call `LoadGeneOMP` for each new gene of a taxon, which enters the global critical
  section even when the gene is loaded. That is once per gene per taxon per sample, too rare to cost anything.

## Asides

- `-t` counts the aligners only. A paired run also has one inflating thread per file, and profiling has a zstd
  read-ahead thread. They are busy about 10% of the time during alignment, so they cost little, but on an
  `N`-core allocation they share the cores with the `N` aligners.
- The output handler takes 16–18% of thread time at 6 threads (30–34 s summed). zstd is a third of that; the
  rest is SAM formatting and `ExtractSNPs`. The handlers use `ExtractSNPs` only as a consistency check and drop
  its SNP list. This is CPU, not locking, and was not split further here.
- **The index build** was not re-measured. What stays serial there is described in
  [Parallelising the index build](../2026-09-29-index-build-parallel.md) and
  [Faster index builds](../2026-09-30-index-build-gains/README.md). In the partitioned passes, one thread reads
  each round's batches in an `omp single` while the others wait at its barrier. Reading the next round while the
  current one is processed would hide that.

## How it was run

```bash
bash scripts/build.sh          # ~/mt-audit/bin/protal-instr, protal-base, readbench/readbench
bash scripts/scaling.sh        # 8 runs; then scripts/summary.sh > results/scaling_summary.tsv
bash scripts/readbench.sh      # BGZF and uncompressed, 1-6 readers, two passes
bash scripts/readbench_gzip.sh # standard gzip (pigz -6) against BGZF, 1-2 readers
```

The work folder is `~/mt-audit` in WSL; set `MT_AUDIT` for another. It holds the outputs and the uncompressed and
pigz copies of the sample (7 GB). `results/run_<name>.txt` has each run's stage table (`_runtime.tsv`), its
`MTAUDIT` counters, its stage times and `/usr/bin/time` lines. In the counters, `reader_hold` counts two entries
per batch for paired reads, because the patch adds the hold twice in `SeqReaderPE::LoadBlockOMP`, once at ~0 s;
the sum is right. The second and third passes of `readbench_gzip.sh` ran at load 10–13 and are not used.

## Follow-up: the quick fixes and parallel profiling, implemented

Later on 2026-10-01, on branch `audit-fixes`, the same machine and data as above. Three commits:

- **`f8b1b6b`**: recommendation 2, `schedule(dynamic, 1)` for the per-sample profiling loop. Also the
  correctness note: `Genome`'s loaded flag is now an atomic, set with release and read with acquire.
- **`2f637f9`**: recommendation 1, one sample profiled on all threads with the same profile.
  - A thread of its own reads and decompresses the SAM. It cuts the text into chunks of whole reads
    (`src/Profiling/SamChunks.h`): a chunk ends before a record whose QNAME differs from the one before.
  - Waves of chunks (1 MB each, four per thread) are parsed and prepared on all threads
    (`Profiler::ParseChunk`, `PrepareMAPQ`). Preparing means the record evidence into a collector per chunk,
    the MAPQ and length filters, and the database checks.
  - In file order come the reads' numbers and links and the taxa made (`PlanChunks`). Then each taxon adds its
    records on one thread, in file order (`RunJob`), while the next wave is parsed and the last one freed.
  - Every taxon therefore gets the same records in the same order, and the same maps are built in the same
    order, so every floating-point sum and output is the same.
  - In one case the parallel path cannot copy the serial one: a record that does not fit its gene would leave
    its taxon without records, and the serial profiler removes that taxon on the spot. There the file is
    profiled again on one thread; that happens only with a SAM aligned against another database.
  - `PostProcessSNPs`, the model scores (`ScoreTaxa`) and the lines of `.profile.log`, `.gene.log` and
    `.genes.log` are computed per taxon on the sample's threads and written in taxon order.
  - Three serial costs went too: `AlignmentPair` copied every record twice; `IntTaxonomy::LineageStr` copied
    each node on the way to the root, with its children, for every gene line; the gene coverage columns were
    rebuilt from scratch for each value.
- **`f5ee645`**: samples start largest SAM first, and each gets threads in proportion to its SAM's share of
  all the samples' bytes, at least one. Measuring `2f637f9` showed why this is needed: with more samples than
  threads, each sample got one thread, so a deep sample among shallow ones gained nothing.

### Checks

- **Byte-identical outputs** (`scripts/followup/compare_profiles.sh`, `results/followup/identity_*.txt`).
  Against 0.7.1 (`1c11a00`), every output file of `--profile_only` is the same: profiles, `.profile.log`,
  `.gene.log`, `.genes.log`, `misc/` statistics, `.err`, and the strain MSAs of the multi-sample runs.
  - Each sample alone at 1, 3 and 6 threads: eight SAMs of the full database (paired-end 5M, 500k ×2, 10k
    and 1k pairs, single-end, PacBio, Nanopore) and three of the training database.
  - The eight together (3,543 files) and the three together (647 files) at 6 and 16 threads; with
    `f5ee645` also at 3.
  - 43 runs in all. The logs of the 37 runs with `2f637f9` are the same apart from the output path.
- **Tests on a clean tree** (`git archive` plus the patch; `scripts/followup/check.sh`): unit tests
  258/258 and e2e 119/119 for each commit.
  - The 7 new unit tests (`tests/test_ProfileThreads.cpp`) cover the cut rule and the chunk reader. They also
    check that a profile is the same, map order and floating-point values included, on 1 to 8 threads and
    in chunks of 1 byte to 4 MB.
  - Their SAMs have pairs, single mates, secondary candidates, long reads in parts, low MAPQ, alternatives,
    records of genes the database lacks or past a gene's end, records that do not fit their gene, unusable
    records, blank lines and CRLF.
  - They also check the same errors, line numbers and rejected records for malformed and foreign SAMs,
    plain, BGZF, gzip and zstd SAMs cut short, and the fallback to one thread.
- **ThreadSanitizer** finds nothing on the 16 profiling tests (`scripts/followup/tsan.sh`).

### Speed

`results/followup/speed.tsv` (`scripts/followup/speed.sh`, `speed_b2.sh`). The 0.7.1 binary is alternated with
this change's, niced, at load 0.7–2.3 before each run (C ran at 1.8–4.5). Medians of three:

**Profiling one 5M-pair SAM (`--profile_only`), seconds of the profiling stage:**

| `-t` | 0.7.1 | `2f637f9` | speed-up | wall 0.7.1 | wall `2f637f9` |
|---|---|---|---|---|---|
| 1 | 11.87 | 11.14 | 1.07 | 12.65 | 11.80 |
| 2 | 11.39 | 6.31 | 1.80 | 12.04 | 7.00 |
| 3 | 11.62 | 4.93 | 2.36 | 12.27 | 5.65 |
| 4 | 11.95 | 4.46 | 2.68 | 12.61 | 5.24 |
| 6 | 11.82 | 3.90 | 3.03 | 12.47 | 4.70 |

Peak memory goes from 1.18 GB to 1.20–1.27 GB, for the waves held in memory.

**Eight samples together at `-t 6`** (the full database's eight SAMs, among them the 5M-pair one): 0.7.1
took 12.8 s (12.5–14.0). `2f637f9` took 12.4 and 13.3 s, no gain, because the deep sample still had one
thread. `f5ee645` takes **5.4 s** (5.36, 5.41, 6.78), with peak memory 2.22 → 2.29–2.36 GB.

**Whole run (alignment and profiling) of the 5M-pair sample at `-t 6`, two runs each:**

| | 0.7.1 | `2f637f9` |
|---|---|---|
| wall | 42.2 / 42.6 s | **33.6 / 33.7 s** (−21%) |
| aligning | 27.1 / 26.1 s | 26.8 / 26.6 s |
| profiling | 13.6 / 15.1 s | 5.0 / 5.5 s |
| peak memory | 3.66 GB | 3.66 GB |

Profiling takes longer inside a whole run (5.0–5.5 s) than alone (3.9 s). In the whole runs the index sits in
memory beside it, and the machine was more loaded.

### What is left

`results/followup/breakdown.txt` (`scripts/followup/proftimers.patch` on `2f637f9`). Profiling the 5M-pair SAM
at 6 threads takes 3.6–3.8 s:

| Part | Time |
|---|---|
| combined phase (parse, add, free), parallel | 2.6–2.7 s |
| in file order (`PlanChunks`) | 0.11 s |
| rejected reads | 0.06 s |
| SNP post-processing | 0.23 s |
| scoring | 0.18 s |
| writing | 0.32 s |

The reader thread decompresses for 1.6 s and waits for room for 0.7–0.9 s; the parsers never wait for it.

On many threads the combined phase shrinks until the reader thread bounds it: about 1.6 s for this SAM. So
profiling a sample would take about 2 s at 32 threads, against 13–20 s before. That puts the 32-thread
single-sample run projected in the summary at about 10 s instead of about 21 s.

What would come next:
- **Decompress the frames of a `.sam.zst` in parallel.** Its seek table lists them, so the reader thread
  stops being the limit.
- **Parse SAM lines more cheaply.** `SamReader` splits each line into a vector of strings, and parsing is
  half of the CPU time.
- **The input path for aligning on 40 or more threads**, from the summary above: larger reader batches,
  several inflate threads for BGZF.
- **The serial copy of the SAM records behind the header.**
