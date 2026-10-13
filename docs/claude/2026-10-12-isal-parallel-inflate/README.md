# One gzip member on several threads with ISA-L alone: a feasibility prototype

2026-10-12. Question: can protal inflate a single-member gzip FASTQ (as `gzip` writes it) on several threads
with ISA-L's own decoder, faster than rapidgzip ([report](../2026-10-12-rapidgzip-input/README.md): ~0.3 GB/s a
thread on FASTQ, 0.87-1.08 GB/s on 4 threads, against ISA-L's ~1.7 GB/s on one)? Answer: yes. A prototype
outside protal decodes every piece **byte-identical** to the file. It runs at **2.69 GB/s on 4 threads**,
1.56× sequential ISA-L and 2.5-3× rapidgzip, for 1.42 CPU s/GB against ISA-L's 0.58.

Follow-up (2026-10-13): implemented in protal in place of rapidgzip, on the system's stock ISA-L (no fork, no
nasm), at 2.0 and 2.4 GB/s on 2 and 4 threads through protal's stream ([In protal](#in-protal-2026-10-13)).

## Setup

- **Code**: `scripts/isalpar.cpp` (the prototype; modes `bounds`, `seq`, `cost`, `par`) and `scripts/bfinder.cpp`
  (rapidgzip's block finder timed alone), against the ISA-L protal builds since `f9df09d` (`lib/isa-l`: 2.32.1 with
  rapidgzip's changes, of which the prototype uses the stopping points), built by `scripts/build.sh`
  (`g++ 13.3 -O2 -march=native`). Not in protal.
- **Data**: `~/gzpar/data/R1.member.gz`, the rapidgzip report's file: 1,681,851,361 bytes of simulated 150 bp
  Illumina FASTQ as one gzip member of 397,409,724 bytes (6,392 deflate blocks, ~257 KiB of output each). Also a
  `pigz -b 32` file and a `gzip --rsyncable` file (30,293 pieces), both generated and since deleted.
- **Machine**: WSL on an Intel Core Ultra 7 258V, `taskset -c 0-3 nice -n 5`, the file in the page cache. The
  reported timings (`results/timing3*.tsv`, `scripts/run_quiet.sh`) are the median of 4 alternated rounds, each
  started only after 5 s with no other busy process. Two earlier sets ran while other work kept cores 0-3 busy
  (sequential ISA-L fell to 0.55 GB/s) and were discarded.
- The piece starts come from one sequential pass (`isalpar bounds`, 1.07 s, as long as a full decode), used only
  to set up the prototype. A real decoder finds them with a block finder in parallel (see below).

## How it decodes

The file is cut into pieces of 1, 4 or 16 MiB of compressed data, each starting at a deflate block.

**Decoding a piece without its window.** A piece's decoder does not know the 32 KiB window before it, so it starts
with a *labelled* fake dictionary (`isal_inflate_set_dict` right after `isal_inflate_init`):
- **Start.** The input is primed mid-byte (`read_in`, `read_in_length`, as rapidgzip's `inflatePrime` does) and
  decoded as `ISAL_DEFLATE`.
- **End.** Input is fed up to 64 bytes before the piece's end without stops, then the decoder stops at every
  `END_OF_BLOCK` until its bit position is the piece's end.

**Recovering the window references.** A byte copied from window position p comes out as that position's label.
Two passes with different labels tell them apart:
- A(p) = p & 0xff, and B(p) = (A + 1 + (p >> 8)) & 0xff, so A ≠ B for every p.
- Wherever the two outputs differ, the byte is a *marker*, and p = A | (((B − A − 1) & 0xff) << 8).
- Where they agree, the byte is the real one.
- The labels the task brief proposed (three passes) were dropped: 82.8% of markers point at p ≥ 16,384, and
  0.54% of them would have been ambiguous after two passes.

**Resolving the markers.** Once the window before the piece is known (the end of the piece before), the markers
are filled in. The piece's own last 32 KiB is resolved first: that is the next piece's window, so the windows
chain in order (~0.005 s for the whole file). The rest is resolved in place in pass B's buffer (merge and resolve
fused, no marker list) in the worker pool. The CRC-32 and ISIZE are checked on the consumer (~0.04 s/GB).

**Two schedules.**
- **spec**: every thread decodes pieces speculatively, and the main thread only chains windows and consumes.
- **lead**: the main thread decodes in order with the true window (one pass, as sequential ISA-L) whenever it
  reaches a piece no thread has started. The other threads speculate at least 2 × threads + 4 pieces ahead.
  At 1 thread this is ISA-L in pieces (3% slower than plain sequential). At 4 threads with 1 MiB pieces the lead
  decoded 160 of the 376 pieces itself.

## Correctness

Every output was byte-identical to the file: CRC-32 and ISIZE checked, every piece's offset, and a chunked
`memcmp` against `R1.fq`.
- All 24 verify configurations (`results/verify.log`): 1, 4 and 16 MiB pieces, both schedules, both label
  schemes, 1 and 4 threads.
- All 112 timed runs (`results/timing3.tsv`, column `OK`).
- Every one of R1's 6,392 blocks as a piece of its own (`results/every_block_pieces.log`).
- The `gzip --rsyncable` file: 30,293 pieces, of which 16,066 are non-empty pieces under 32 KiB and 14,127 are
  empty stored blocks (`results/short_pieces_rsync.log`).
- The `pigz -b 32` file (`results/short_pieces_mixed.log`).

## Markers: why a window-less decode stays costly on FASTQ

| Pieces | Markers, share of the output | Pieces whose markers die out |
|---|---|---|
| 1 MiB | 47% | 0 of 375 |
| 4 MiB | 38% | 0 of 93 |
| 16 MiB | 19% | 0 of 23 |

- **The markers never die out.** In 4 MiB pieces their share falls from ~56% in the first MiB of output to ~27% in
  the last, and to ~7.5% after 70 MB. rapidgzip's switch to ISA-L (once a piece's window is fully resolved)
  therefore never triggers on this file, which is the rapidgzip report's finding.
- **By line type** (4 MiB pieces, `results/cost_4MiB_linetypes.log`): quality lines hold 78% of the markers (66% of
  quality bytes are markers), headers 19% (73% of header bytes), sequence lines 2% (1.7% of sequence bytes).
- 83% of the markers point into the newer half of the window, and runs of markers are short (mean length 1.3).

## CPU cost

One thread, s per GB of output (`results/cost_quiet.log`):

| Step | 1 MiB | 4 MiB | 16 MiB |
|---|---|---|---|
| ISA-L sequential (gzip, CRC) | 0.59 | 0.60 | 0.59 |
| one pass with the true window | 0.57 | 0.59 | 0.60 |
| one labelled pass | 0.58-0.61 | 0.58-0.59 | 0.58-0.61 |
| merge + resolve, fused (2-pass labels) | 0.385 | 0.348 | 0.234 |
| merge + resolve, fused (3-pass labels) | 0.40 | 0.36 | 0.25 |
| the brief's scheme: 3 passes, a marker list, its resolution | 2.91 | 2.72 | 2.45 |
| **a speculative piece, 2 passes, fused** | **1.54** | **1.50** | **1.41** |

A speculative piece costs ~2.6× a sequential decode: two ISA-L passes and the per-marker merge. ISA-L's decoder
writes bytes, and a byte cannot carry a marker flag and a 15-bit position. One pass would take a 16-bit-output
variant of ISA-L's assembly decoder.

## Throughput

GB/s of output, median of 4 rounds (`results/timing3_summary.tsv`):

| Configuration | 1 thread | 2 | 3 | 4 | CPU s/GB at 4 | max RSS MiB at 4 |
|---|---|---|---|---|---|---|
| ISA-L sequential | 1.73 | | | | 0.58 | 383 |
| N ISA-L processes side by side (the machine's capacity) | 1.73 | 3.34 | | 5.59 | | |
| **lead, 1 MiB** | 1.68 | 2.17 | 2.40 | **2.69** | 1.42 | 560 |
| lead, 4 MiB | 1.59 | 2.08 | 2.33 | 2.56 | 1.48 | 1,072 |
| lead, 16 MiB | 1.55 | 1.81 | 2.01 | 2.23 | 1.73 | 2,031 |
| spec, 1 MiB | 0.64 | 1.25 | 1.80 | 2.03 | 1.97 | 474 |
| spec, 4 MiB | 0.66 | 1.24 | 1.74 | 2.09 | 1.88 | 628 |
| spec, 4 MiB, the brief's 3-pass labels | 0.47 | 0.93 | 1.21 | 1.54 | 2.55 | 750 |
| rapidgzip (its report) | | 0.60-0.66 | 0.76-0.96 | 0.87-1.08 | | |

- **Memory.** Max RSS includes the 383 MiB mapped input. At 4 threads the extra over sequential is +91 MiB (spec,
  1 MiB pieces), +177 MiB (lead, 1 MiB) and +689 MiB (lead, 4 MiB).
- **lead against spec.** lead beats spec by 29% at 4 threads and 75% at 2: whatever the lead decodes, it decodes
  once.
- **The serial parts are negligible.** Chaining the windows and the CRC take little time; at 4 threads all four
  cores are busy decoding.

**Block finder.** rapidgzip's block finder (`results/bfinder.log`) scans 115-120 MB/s of compressed data on one
thread. It found 6,390 of 6,391 block starts with no false candidates. Finding a piece's start takes ~0.3 ms
(a 30-37 KiB scan), i.e. 0.02-0.06 s/GB, and runs in parallel.

**Estimate at 8 and 16 threads (extrapolated, not measured).** The fit is lead ≈ e × (1/0.57 + (W − 1)/1.54),
with e = 0.73 from 4 threads. It gives ~4.5 GB/s at 8 threads and ~8 GB/s at 16, 2.7× and 5× one ISA-L thread.
Beyond that, memory bandwidth (~5-6 bytes moved per output byte) and the one in-order consumer per file would
limit it. With the prototype's in-flight bound, lead with 1 MiB pieces would hold up to ~0.6 GB of buffers at 16
threads, so the bound needs tightening.

## What it means for protal

- **Replace rapidgzip.** A decoder on this scheme in `src/IO/ParallelGzip.cpp` (lead, 1 MiB pieces, two-pass
  labels, rapidgzip's block finder for the piece starts) would replace rapidgzip's decoding for single-member
  gzip. Per thread it is ~2.5× faster on FASTQ.
- **Where it pays.** It pays from 2 threads per file (+25% at +57% CPU) and 4 (+56% at 2.4× CPU), where rapidgzip
  never beats one ISA-L thread. `GzipThreadsFor` could then hand gzip 2-4 threads per file wherever the input
  limits the run. Paired-end runs on 32 threads hit that floor, and more so with the sparse-first gates (the
  alignment is 35-52% faster), so a paired-end gzip run on 16-32 threads would gain the most.
- **The checks stay.** `ParallelGzip`'s checks (pieces in file order, ISIZE, the last trailer at the file's end)
  still apply. A speculative piece is used only if the piece before ends exactly where it starts; otherwise it is
  discarded and decoded in order.
- **Not tried in the prototype:**
  - several members (gzip headers in mid-stream);
  - stored and fixed-Huffman blocks as piece starts found by the real finder;
  - rejecting a finder's false candidate once ISA-L fails on it.

  These, the in-flight bound and protal's tests would take an estimated 2-3 days.
- **BGZF stays the fastest input:** 1-4 ISA-L workers per file, each at sequential speed.

## In protal (2026-10-13)

At the user's request the decoder replaced rapidgzip completely, on the ISA-L that protal linked before
`f9df09d`, the system's: `libisal-dev` 2.31 or conda's `isa-l`, and, for `just static`, a downloaded 2.32.1 built
with nasm. The vendored `lib/isa-l`, `lib/rapidgzip` and `lib/rapidgzip-zlib` are gone, and so is nasm for an
ordinary build. `src/IO/ParallelGzip.cpp` has the same interface as before, and `ThreadedGzStreambuf` uses it for
a regular gzip file that is not BGZF and gets more than one thread. The tables below are for that commit,
measured in WSL on cores 0-3 as above.

**No fork of ISA-L.** The prototype stopped at a piece's end with the stop points of the ISA-L that rapidgzip
patched. Stock ISA-L has none, so the input given for a piece ends at the piece's end bit instead:
- **The masked last byte.** That byte is given with its bits from the end on zeroed. Zero bits begin a stored
  block, whose header then waits for more input.
- **Where the header began.** When stock ISA-L runs out of input inside a block header, it keeps the header's
  bytes and restores its bit buffer to where the header began (`read_header_stateful`: `tmp_in_size`,
  `read_in_length`). So the bit at which ISA-L waits for a header is
  (input consumed − `tmp_in_size`) × 8 − `read_in_length`.
- **The test.** A piece ended at a block boundary if ISA-L waits for a header that begins exactly at the piece's
  end. Otherwise it failed, or it ended in a final block before that, which is the member's end.

**Finding piece starts.** rapidgzip's finder went with rapidgzip, so the decoder has its own:
- **Candidates.** The first three bits (BFINAL 0, BTYPE 2) are tested for 56 positions at once.
- **The common refusal.** Each candidate's code-length code is checked for completeness by a branch-free Kraft
  sum, four lengths per table lookup.
- **Full check.** Survivors, 0.36% of candidates, decode their code lengths and check the literal/length and
  distance codes as zlib does.
- **Cost.** On R1 that is ~1 candidate per byte scanned, ~30 KiB scanned per piece, and 0.09-0.10 s of CPU for the
  whole file. A per-length loop took 0.27-0.30 s.
- **False starts.** A start found by chance is never used: a piece is taken only if the one before ends exactly at
  its start. When the reader's own in-order decode of a range fails, the range's end start is dropped and the
  range decoded again to the next start. After two such drops in a row, the range is decoded to the member's end
  and its bytes given as they come, so a real error shows where it is.

**Scheduling.** Workers take, in this order of priority:
1. resolving pieces the reader waits for;
2. finding piece starts ahead;
3. speculating from 2 × threads + 4 pieces ahead of the chain, up to 2 × threads further.

The reader's thread decodes in order whatever piece no worker took. While its next piece waits for a worker to
resolve it, it decodes further pieces in order rather than resolving that piece itself. Buffers are pooled; each
new allocation was faulted in and reallocated once per piece. Members, trailers (CRC-32 and ISIZE on the reader's
side), zero padding and data after the last member are handled as `ThreadedGzStreambuf`'s one-thread path handles
them. `GzipThreadsFor` gives gzip a thread per 4 of `-t` and file, 1 to 8: 2 per file for a paired-end run on 16
threads, 4 on 32.

**Tuning (R1, 4 threads).** Values are the CPU seconds of each step and the speed:

| Version | Finder | Speculation | Worker resolves | Reader's in-order decoding | Reader's copy and CRC | GB/s |
|---|---|---|---|---|---|---|
| first draft (no pool, reader resolves) | | | | | | 1.91-1.94 |
| pool, reader decodes while waiting | 0.29 | 1.25 | 0.36 | 0.46 | 0.15 | 2.37-2.43 |
| + table Kraft sum in the finder | 0.09-0.10 | 1.27-1.33 | 0.37-0.39 | 0.43 | 0.15 | 2.23-2.43 |

- **Piece size.** Pieces of 2 and 4 MiB were slower: 1.57-1.80 and 1.45-1.46 GB/s on 2 threads, 1.72-1.80 and
  1.84-1.85 on 4.
- **Against the prototype.** The prototype gave 2.17 and 2.69 GB/s on 2 and 4 threads. The difference is protal's
  stream: the reader's thread copies every byte into `ThreadedGzStreambuf`'s blocks, and the benchmark's own
  thread copies them out again, a fifth thread on 4 cores. The prototype only computed the CRC.

**Tests.**
- **Small pieces.** `ParallelGzip::SetPieceBytes(bytes, gap)` lets the unit tests use 16 KiB pieces, with the
  workers starting at the next piece, so a few MB of gzip are hundreds of pieces.
- **Pieces ahead.** `GzipPiecesDecodedAheadReadAsTheFile` covers levels 1, 6 and 9, FASTQ of repeats, and members of
  every level with zero padding. It runs on 2, 3 and 5 threads, at gap 1 and the default.
- **False starts.** `AFalseBlockStartIsNeverUsed` copies a real dynamic block header into stored data.
- **Damaged and cut files.** The damaged-file and cut-file tests run with small pieces too.
- **Counts.** `ThreadedGzStreambuf::parallel_gzip_counts()` reports the pieces decoded ahead, in order and to the
  member's end, and the starts dropped.

**Other gzip readers and writers in protal.**
- **Reads and SAM input.** Read files and SAM input (`SamFile.h`), as well as `LongReadSimulator` and the read-type
  detection, all read through `ThreadedGzStreambuf`, so they take the pieces path wherever a gzip file is not BGZF.
  A `.sam.gz` that protal wrote is BGZF and keeps the BGZF workers.
- **Genome files.** `GenomeStore`'s `Inflate` reads whole genome files of a few MB, and its callers already read
  many genomes on several threads. Pieces of 1 MiB would not pay there.
- **Writing.** protal writes gzip only as BGZF (`Bgzf.h`, ISA-L's deflate on several threads); this decoder does
  not apply to writing.
- **gzstream and zlib-ng.** `lib/gzstream` over zlib-ng is still built and its header included (`Options.h`,
  `Profiler.h`, `RunProtal.h`), but no class of it is used in `src/`. Only the tests use zlib-ng, to write gzip at
  given levels. Dropping it from `protal_lib` would be a cleanup of its own.

**Final speed** through `ThreadedGzIstream` (`protal/bench_gz2.cpp`, `protal/pgz_tune.sh`; R1, cores 0-3, 3
alternated rounds, `protal/bench_final.txt`):

| Threads | GB/s | CPU s | max RSS MB (with the 397 MB mapped file) |
|---|---|---|---|
| 1 (ISA-L in order, as before) | 1.67-1.70 | 1.05-1.07 | 10 |
| 2 | 1.99-2.05 | 1.70-1.74 | 438-442 |
| 3 | 2.34-2.37 | 2.17-2.20 | 509 |
| 4 | 2.35-2.48 | 2.39-2.43 | 557 |

rapidgzip on the same file and cores gave 0.56-0.66 GB/s on 2 threads and 0.87-1.08 on 4. A file without dynamic
blocks within the look-ahead, such as zlib-ng's level 1 (fixed codes only), reads in order, with its pieces looked
for only that far ahead (`protal/pgz_l1.sh`; GNU gzip's level 1 has dynamic blocks and reads in pieces). Tests and
commands: `protal/pgz_full2.sh`.
