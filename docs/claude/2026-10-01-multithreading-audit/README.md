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
   that further out. (Re-measured on a quiet machine in follow-up 2 below: 1.8M pairs/s for gzip, 2.7–3.4M for
   BGZF, binding from about 50 and 75 threads; where the input nears its limit, aligners spinning on the reader
   lock keep the inflating threads off the cores.)
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

These ceilings were measured while other sessions loaded the machine. Follow-up 2 below measured them again
with only this session's runs, and with workers that stand in for N aligners ("The input path, again").

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

## Follow-up 2: cheaper SAM parsing, parallel decompression, and the input path again

2026-10-01 and 02, branch `audit-fixes`, the same machine and data. Two commits, then the input path of the
alignment audited again: could more parallel input make protal faster?

- **`85e53e6`**: SAM lines parsed on views of the line.
  - `SamReader` splits a line into `std::string_view`s (`sam_detail::SplitFields`, the fields `LineSplitter`
    gives), copies only those a `SamEntry` keeps and reads the tags from the views. It also reads a text in
    memory, the profiler's chunks, without a stream.
  - `NormalizeCigar` only checks a CIGAR that needs no rewriting (`IsNormalCigar`: ops M, X, I, D and S, no two
    alike in a row, counts without leading zeros, query length equal to the sequence's), as protal writes them.
  - CIGAR counts (`CigarCount`) and an RNAME's taxon and gene (`NameNumber`) are parsed in place, with the
    values and exceptions of the `std::stoi` and `std::stoul` they replace.
- **`a212559`**: a `.sam.zst` decompressed on threads of its own, from 8 threads on.
  - `zstd::ParallelFrameStreambuf` reads a seekable zstd file in order while threads of its own decompress the
    frames ahead (a window of frames, each in a buffer of its own). A frame that fails stops them; the reader
    gets the frames before it, then the error, as from `IStreambuf`. `SamInput` takes the number of threads.
  - `ProfileSamParallel` uses none below 8 threads and `-t`/6 (2 to 8) from 8 on. Below 8, the reading thread
    keeps ahead of the parsers, and an extra thread only takes cores from them: the 5M-pair SAM on 6 threads
    took 3.0 s without and 3.25–3.44 s with one to three (`results/followup2/speed_decompress.tsv`).
  - `docs/running.md` says so in the `-t` row.

### Checks

- The committed source (`f359d49` + both, built from `git archive` and the patch):
  unit tests 269/269, e2e 121/121 (`results/followup2/final_check.txt`). `--profile_only` of three SAMs
  (paired-end 500k pairs, Nanopore 90 Mb, single-end 500k reads) at 1, 6 and 8 threads gives the same files as
  `f359d49`'s build, byte for byte; at 8 threads with two decompression threads.
- On `25d457e`, before the other sessions' last commits: the parsing commit alone 264/264 unit tests; both
  267/267 and e2e 121/121; 37 `--profile_only` runs byte-identical to `25d457e`'s build in every output file,
  logs the same apart from the output path: the eight SAMs of the full database and three of the training
  database, each alone at 1, 6 and 8 threads, each database's together at 6 and 16
  (`results/followup2/identity_25d457e.txt`).
- New unit tests: `NormalizeCigar` against its former code on 200,000 generated CIGARs; `SplitFields` against
  `LineSplitter`; the counts and ids against `stoi` and `stoul`, errors included; a text read as the same
  stream; the frames in order on 1–9 threads and windows of 1–64 frames, a reader that stops early, a corrupt
  frame; a `.sam.zst` of many frames profiled with 0, 1 and 3 decompression threads.
- ThreadSanitizer: no warnings on the 32 profiling, chunk, frame and reader tests (`results/followup2/tsan.txt`).

### Speed

Instructions (callgrind, inclusive), `--profile_only` of the 500k-pair SAM on one thread, `f359d49` against
`f359d49` + both, the same outputs (`results/followup2/callgrind_*.txt`):

| | `f359d49` | + both |
|---|---|---|
| whole run | 20.82G | 16.81G (−19%) |
| `SamReader::Advance` (zstd 2.33G in both) | 5.86G | 3.82G |
| `PrepareMAPQ` | 2.97G | 1.54G |
| `Taxon::AddSam` | 3.51G | 2.95G |
| `LineSplitter::Split` | 1.23G | not called |
| `NormalizeCigar` | 1.30G | not in the profile (the check is inlined) |

Parsing without zstd went from 3.5G to 1.5G. The largest piece left is `SequenceRangeHandler::CoveredPortion`
(2.38G, 14%: the genes' covered lengths, for the features and the output lines), unchanged.

The profiler's reader path alone (`readsam`: `SamInput` and the chunk reader, nothing done with the chunks), the
5M-pair SAM (1.71 GB of SAM text), medians of three, load 3.1–3.5 (`readsam.tsv`):

| decompression threads | 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| GB/s of SAM text | 0.75 | 1.06 | 2.08 | 2.07 | 1.69 |

With two, the reading thread hands out the SAM in 0.8 s instead of 2.3 s. Past two, cutting the text into chunks
on the reading thread is the limit; four compete for the 6 vCPUs.

Profiling the 5M-pair SAM (`--profile_only`), `f359d49` against `f359d49` + both, alternated, niced, load 3.2–5.7
before each run; medians of three (`results/followup2/speed.tsv`):

| `-t` | `f359d49` | + both | CPU seconds |
|---|---|---|---|
| 1 | 11.23 s | 9.43 s (−16%) | 11.3 → 9.5 |
| 6 | 3.63 s | 3.56 s (−2%; the ranges overlap) | 16.2 → 15.2 |
| 8, on the 6 vCPUs (two decompression threads) | 4.30 s | 3.62 s (−16%) | 18.9 → 15.6 |

At 6 threads the parsing is spread over the threads and the parallel phase as a whole sets the time, so the
cheaper parsing shows as CPU (−6%) more than as time. At 8 threads on 6 vCPUs the old build slows down; the new
one does not.

### The input path, again

The question: when protal aligns on many threads, does the input path (one inflating thread per FASTQ file, the
reader lock, batches of 32 pairs) hold it back, and would more parallel input make whole runs faster?

What was run (`scripts/followup2/`, results in `results/followup2/`):
- `protal-inexp`: `a212559` with `batch_experiment.patch`, which only makes the batch size of the readers
  settable (`PROTAL_EXPERIMENT_BATCH`, default 32 as in protal). `--no_profile`, so the alignment alone.
- The 5M-pair sample as BGZF (as the simulator writes it), as standard gzip (`pigz -6`, one member: the
  zlib-ng path) and uncompressed, against the full database (765 species) and the mini database of the tests
  (3 species: a pair costs k-mer extraction and lookups and little else, about 1.9 µs on one thread, so the input
  path can become the limit on 6 threads).
- `readbench` of the first audit, now built against the same tree, for the input path alone, and with workers
  that spin a fixed time per pair to stand in for aligners (below).
- The machine. The first runs (`input_audit.sh`, first repetition) ran while other sessions' simulations kept
  the load at 9–13 and the inputs partly out of the page cache; they are left out. Everything below ran with
  only this session's runs (load 3–6 before each, mostly the decay of the run before; each run waited up to a
  minute for less than 6). Even then the full-database runs at 6 threads vary by ±8% from one run to the next
  for the same input and build, more than some effects measured, so only medians of three or more, alternated,
  are used.

#### At 6 threads, input does not hold protal back

Full database, 6 threads, batch 32, medians of three (`runs2.tsv`):

| input | aligning | range | CPU seconds |
|---|---|---|---|
| uncompressed | 22.57 s | 22.52–23.49 | 137 |
| BGZF | 23.05 s (+2%) | 22.88–23.36 | 141 |
| gzip | 25.25 s (+12%) | 24.08–26.89 | 154 |

The aligners spend 0.7–1.1 s of 23–25 s in the reader stage (taking and parsing their batches, waits
included; `Sequence reader` in `_runtime.tsv`). In the three other series with this build at the default
(`waitpolicy.tsv`, `mutex.tsv`, `hybrid.tsv`), BGZF took 0–9% longer than uncompressed input and gzip was
within 3% of BGZF; the +12% here is mostly the spread of the gzip runs. The difference to uncompressed input
is the inflating threads' CPU on the 6 shared vCPUs (4–6 CPU seconds on a free core, `readbench` below), not a
limit on how many reads arrive: parallel inflation would spend the same CPU and not remove it.

Larger batches change nothing here: batches of 128 and 512 pairs were within the noise of 32 for every input on
both databases (`runs.tsv`, `runs2.tsv`).

#### The ceilings, on a quiet machine

`readbench` with one reader thread that does nothing with the reads, medians of three (`readbench.tsv`):

| input | pairs/s | CPU seconds |
|---|---|---|
| uncompressed | 3.71M | 1.7 |
| BGZF (libdeflate) | 3.37M | 4.1 |
| gzip (zlib-ng) | 1.99M | 6.1 |
| gzip, R2 twice (the larger file alone) | 1.83M | 7.1 |

These are higher than the first audit's (1.15M for gzip, 1.4–1.7M for BGZF), which ran at load about 6 from
other sessions: an inflating thread that has to share its core is what sets the ceiling. One zlib-ng thread
inflates R2 at about 1.83M pairs/s (each file's 1.68 GB in 2.7 s); libdeflate is at least 1.8 times as fast,
close to what the reader itself takes.

#### Can the input path feed N aligners?

`demand.sh`: `readbench` with 4 workers, which with the 2 inflating threads fill the 6 vCPUs. Each worker
spins 28.8 × 4 / N µs per pair, so together they ask for as many pairs per second as N aligners at 34.7k pairs/s
each (the full database's rate per thread at 6 threads). Two repetitions, batches of 32, 128 and 512
(`demand.tsv`); delivered pairs/s as a share of uncompressed input's at the same demand:

| aligners | demand | BGZF / uncompressed | gzip / uncompressed | uncompressed / demand |
|---|---|---|---|---|
| 6 | 0.21M | 0.97–1.02 | 0.99–1.03 | 0.93–0.97 |
| 12 | 0.42M | 1.00 | 0.99–1.00 | 0.94 |
| 24 | 0.83M | 0.95–0.99 | 0.99 | 0.91 |
| 48 | 1.67M | 0.97–1.02 | 0.91–0.98 | 0.81–0.84 |
| 96 | 3.33M | 0.94–0.96 | 0.70–0.73 | 0.70–0.73 |
| 192 | 6.66M | 0.70–0.74 | 0.45–0.49 | 0.54–0.57 |

- Up to 48 aligners, gzip and BGZF deliver what uncompressed input does: the input path keeps up.
- gzip tops out at 1.6–1.8M pairs/s, the single inflating thread of R2; it binds between 48 and 96 aligners.
- BGZF tops out at 2.7–2.8M; it binds between 96 and 192.
- Uncompressed input's share falls as a fixed cost per pair in the workers predicts (taking the records,
  about 0.5 µs: 2.4/2.9 = 0.83, 1.2/1.7 = 0.71, 0.6/1.1 = 0.55), so the reader lock did not limit at
  3.8M pairs/s. 4 workers make the lock's entries as often as N aligners would, but not N contenders; on a
  machine with 100 cores the hand-offs could cost more.
- The batch size changes nothing beyond the noise, except uncompressed input at 192 (+6% at 512).

#### Inside protal, the inflating threads compete with the aligners

`-t` counts the aligners, and each read file has its inflating thread besides. On an allocation of N cores with
`-t N` the inflating threads share the cores, and near the ceiling that is what limits. The mini database
(`runs2.tsv`, medians of three):

| input | 2 threads | 4 threads | 6 threads |
|---|---|---|---|
| uncompressed | 0.96M pairs/s | 1.45M | 1.64M |
| BGZF | 1.01M | 1.37M | 1.18M |
| gzip | 0.98M | 1.24M | 0.95M |

From 4 to 6 threads uncompressed input still gains, while BGZF and gzip lose; the aligners then spend 1.5 s
(BGZF) and 3.0 s (gzip) of 4.3–5.3 s in the reader stage. gzip delivers 0.95M pairs/s inside protal at 6
threads, against 1.7–1.8M with cores to spare.

Part of it is how the aligners wait. The reader lock is `#pragma omp critical(reader)`, and its holder reads the
next batch from the stream; when the inflated block is used up, it waits for the inflating thread inside the
lock. libgomp's waiters spin before they sleep, 300,000 rounds by default when there are no more OpenMP threads
than cores (OpenMP does not count the two `std::thread`s that inflate). So five aligners spin while the one
holding the lock waits for an inflating thread that the spinners keep from the cores.

`OMP_WAIT_POLICY=passive` (waiters sleep at once), alternated with the default, medians of three
(`waitpolicy.tsv`):

| run | default | passive |
|---|---|---|
| mini, 6 threads, BGZF | 4.01 s | 2.88 s (−28%) |
| mini, 6 threads, gzip | 4.62 s | 4.04 s (−13%) |
| mini, 6 threads, uncompressed | 3.23 s | 3.33 s (+3%) |
| mini, 4 threads, each input | | +2–4% |
| full, 6 threads, BGZF | 24.47 s | 23.74 s (−3%) |
| full, 6 threads, gzip | 23.81 s | 22.05 s (−7%; all three runs faster, ranges apart) |
| full, 6 threads, uncompressed | 22.45 s | 23.30 s (+4%) |

Sleeping helps when the cores are oversubscribed and the holder waits for inflated data (6 threads, compressed
input), and costs a little where the waits are short hand-offs (uncompressed input, 4 threads).

Passive waiting is a process-wide setting, so three ways to wait on this lock alone were measured, each
alternated with `critical(reader)` at libgomp's default in the same series, medians of three:
- a `std::mutex`, whose waiters sleep at once (`mutex.patch`, `mutex.tsv`);
- libgomp's spin count, `GOMP_SPINCOUNT` 300k (the default), 30k, 3k, 300 and 0, for gzip and uncompressed
  input (`spin.tsv`);
- a lock that tries `try_lock` with a pause N times and then sleeps on a `std::mutex`, N = 200, 2,000 and
  20,000 (`hybrid_exp.patch`, `hybrid.tsv`).

Time against `critical(reader)` in the same series:

| run | `std::mutex` | spin count 3k | 2,000 tries |
|---|---|---|---|
| mini, 6 threads, gzip | −19% | −9% | −23% |
| mini, 6 threads, BGZF | −29% | | −3% |
| mini, 6 threads, uncompressed | +4% | +6% | +2% |
| mini, 4 threads, gzip | +10% | −3% | +1% |
| mini, 4 threads, BGZF | +3% | | 0% |
| mini, 4 threads, uncompressed | −6% | +1% | −4% |
| full, 6 threads, gzip | −7% | −6% | −6% |
| full, 6 threads, BGZF | +2% | | +1% |
| full, 6 threads, uncompressed | +1% | −2% | −2% |

- Sleeping at once (`std::mutex`) gains where the holder waits for inflated data and loses up to 10% where the
  lock changes hands quickly (mini, 4 threads, gzip).
- A short try before sleeping keeps the gain and loses nothing measurable: 2,000 tries (tens to hundreds of
  microseconds) is never worse than 2% and is 6% faster for gzip on the full database. 20,000 tries is as slow
  as libgomp's default; 200 is close to 2,000.
- The full database with gzip input is the one full-database result that holds across the series: medians of
  23.4–23.8 s with `critical(reader)` in the four lock series (25.3 s in the format series above), 22.0–22.1 s
  in five of the seven variants that sleep soon (passive, `std::mutex`, spin counts 300 and 3k, 2,000 tries),
  22.8 s and 23.6 s in the other two. For BGZF and uncompressed input the differences are within the runs'
  spread.

`reader_lock.patch` (on `a212559`, not committed) is that lock: `ReaderLock` in `SeqReader.h` with 2,000 tries,
used by `SeqReader`, `SeqReaderSE` and `SeqReaderPE` in place of `critical(reader)`, and a unit test that
single-end readers on 8 threads take every read once (`ThreadsTakeEverySingleEndReadOnce`). On a clean tree it
passes the unit tests (270/270) and the e2e tests (121/121) (`reader_lock_check.txt`). ThreadSanitizer was not
run on it: the readers run on OpenMP threads, and libgomp is not instrumented.

#### What more parallel input would buy

A whole 5M-pair run on N threads, from the parts measured here and in the first audit: index load 0.8 s, the
SAM records behind the header 1.4 s (serial), profiling about max(0.8 s, 15 CPU s / N) + 0.3 s (the reader path
at 2.1 GB/s, the CPU of 6 threads, the serial steps), and aligning at 36k pairs/s per thread (BGZF at 6
threads), but not faster than the input ceiling (gzip 1.7M, BGZF 2.7M, the reader at least 3.8M pairs/s). That
per-thread rate is the laptop's at 6 threads; on a server with many cores it is likely lower, which moves the
ceilings to more threads.

| N | aligning, no ceiling | run, input without limit | gzip | BGZF |
|---|---|---|---|---|
| 16 | 8.7 s | 12.1 s | 12.1 s | 12.1 s |
| 32 | 4.3 s | 7.6 s | 7.6 s | 7.6 s |
| 48 | 2.9 s | 6.2 s | 6.2 s | 6.2 s |
| 64 | 2.2 s | 5.5 s | 6.2 s (+14%) | 5.5 s |
| 128 | 1.1 s (reader: 1.3 s) | 4.6 s | 6.2 s (+35%) | 5.2 s (+12%) |

These are the ceilings with cores to spare. Inside protal the inflating threads share the aligners' cores: on
these 6 vCPUs gzip then delivered 0.95M pairs/s (mini database, 6 threads), and 1.24M with 4 aligners or with
a reader lock that sleeps. Two extra threads are a third more than 6 cores but 6% more than 32, so on a large
node the loss should be smaller. Where it matters, the remedy is to keep the inflating threads running (the
lock above, or `-t` two below the cores), not to inflate in parallel.

So:
- **Below about 50 threads with gzip input, and 75 with BGZF, more parallel input would not make protal
  faster** on this database; on GTDB, where a pair costs more to align, it would take more threads still.
- **What does help now is how the aligners wait for the reader lock**: with gzip input, the usual case for
  real reads, a lock whose waiters sleep after a short try aligns about 6% faster at `-t` equal to the cores
  (about 5% of a whole run), and up to 23% where the input is the limit (`reader_lock.patch`, not committed).
- **From about 64 threads, a run with gzip input would be 11–26% shorter** if R2 were inflated faster.
  Standard gzip cannot be split; a faster inflater (ISA-L, declined earlier for portability) or speculative
  parallel decoding (as rapidgzip does) would be needed. BGZF can be inflated block by block on several
  threads (as the `.sam.zst` frames now are), which helps only above about 75 threads, up to the reader's
  3.8M pairs/s or more.
- **Larger batches** buy nothing measurable at any of the demands tested.
- At 32 threads and more, the serial copy of the SAM records behind the header (1.4 s, a fifth of a 32-thread
  run) costs more than the input path does.

### How it was run

```bash
bash scripts/followup2/final_check.sh        # a212559 from git archive: build, unit and e2e tests; identity against f359d49
bash scripts/followup2/compare_profiles3.sh  # the 37 identity runs against 25d457e's build
bash scripts/followup2/tsan2.sh PATCH        # ThreadSanitizer on the profiling, chunk, frame and reader tests
bash scripts/followup2/cg_commit.sh          # callgrind, f359d49 against a212559's src
bash scripts/followup2/speed_final.sh        # profiling the 5M-pair SAM at 1, 6 and 8 threads
bash scripts/followup2/readsam.sh            # the profiler's reader path alone, 0-4 decompression threads
bash scripts/followup2/input_audit.sh        # builds protal-inexp; inputs, threads, batches; readbench per file
bash scripts/followup2/input_audit2.sh       # the same, three alternated repetitions
bash scripts/followup2/demand.sh             # readbench standing in for 6-192 aligners
bash scripts/followup2/waitpolicy.sh         # OMP_WAIT_POLICY passive against the default
bash scripts/followup2/mutex_bench.sh        # the reader lock as a std::mutex (mutex.patch)
bash scripts/followup2/spin_sweep.sh         # GOMP_SPINCOUNT 300k to 0
bash scripts/followup2/hybrid_bench.sh       # tries, then sleep (hybrid_exp.patch)
bash scripts/followup2/check.sh lock scripts/followup2/reader_lock.patch a212559
```

`speed_decompress.sh` ran before `a212559` was committed, with an experiment build of it in which the number of
decompression threads was settable. The work folder is `~/mt-work` in WSL; the uncompressed and pigz copies of
the sample are the first audit's, in `~/mt-audit/plain`. `summary2.sh`, `policy_summary.sh`,
`demand_summary.sh` and `rb_summary.sh` make the `*_summary.tsv` files in `results/followup2` from the runs.
