# Reading gzip input on several threads with ISA-L

- **Date**: 2026-10-08.
- **Code**: branch `audit-fixes` at `dd57c6f` plus this change: `src/IO/ThreadedGzStream.h`, `src/RunProtal.h` (the
  threads per file from `-t`, a `--verbose` line per input), `tests/test_ReaderAndTimers.cpp`, `docs/running.md`.
- **Machine**: WSL2 Ubuntu 24.04 on the Core Ultra 7 258V, GCC 13.3, ISA-L 2.31.0 (Ubuntu), 4 cores
  (`taskset -c 0-3`, niced). Other sessions loaded the machine at times: the uncompressed reference fell from 6.0M to
  3.5M pairs/s between two rounds, so only numbers within one round are compared.
- **Data**: `rl150_p5000000_s_1` (5M pairs, 1.68 GB of FASTQ per mate) as the simulator wrote it (BGZF), as one gzip
  member (`gzip -6`) and uncompressed (`~/gzpar/data`); for whole runs `rl150_p500000_s_1` (BGZF) on the v0.7.5 world.
- **Scripts**: [`scripts/`](scripts/): `gz_setup.sh` (the data, `readbench` of the
  [multithreading audit](../2026-10-01-multithreading-audit/README.md) built against HEAD), `gz_base.sh` and
  `gz_over.sh` (HEAD's input path), `bench.sh` and `onefile.sh`/`onefile.cpp` (this change), `check.sh` (unit tests,
  whole runs), `tsan.sh`. They were run from a scratch directory, whose path is in them.
- **Question**: can protal read gzip on several threads with ISA-L, now that the r226 paired-end run waits for its
  input again ([AVX-512 report](../2026-10-08-avx512-kernels/README.md#on-the-cluster-f53db5e-on-a-zen-4-node))? If so,
  implement it.

## Summary

1. **For BGZF, yes, and it is implemented.** BGZF blocks are independent deflate streams of at most 64 KB with their
   size in the header. A dispatcher thread reads them in order into batches of up to 1 MB of output, worker threads
   inflate the batches with ISA-L, and the reader takes them in file order. One file's inflation scales almost
   linearly: 1.1-1.2 GB/s of FASTQ on one thread, 2.1-2.3 on two, 2.7-3.4 on three. protal uses 1 to 4 threads per
   file, one per 4 alignment threads and file (`-t 32` paired-end: 4 per file). Outputs are identical.
2. **For one gzip member, not with ISA-L alone.** A deflate stream's block boundaries are found only by decoding it
   from the start, and its matches reach 32 KB back. Parallel inflation of one member needs speculative decoding at
   guessed block starts with an unknown window, as rapidgzip does (MIT or Apache-2.0, a header-only C++17 library that
   hands the decoding to ISA-L once a chunk's window is known). That would be a new dependency, so it is not added
   here.
3. **Which one the r226 run reads is not known yet.** bcl2fastq writes BGZF by default, but the thirteenth run's
   reader sped up with the switch of gzip streams from zlib-ng to ISA-L (2.87 to 1.19 s per thread), which only the
   non-BGZF path had. A run with this build and `--verbose` prints, per input, "BGZF, inflated on N thread(s)" or
   "not BGZF". If the files are one-member gzip, `bgzip -@ 8` makes them BGZF once.
4. **Locally (4 cores) the paired-end input path does not gain**: with two or three readers and two files, the
   inflating threads already share the 4 cores. The gain needs spare cores, as on the cluster, where the sequence
   reader's time rose by 0.65 s per thread when the seeding got faster.

## How it works

`ThreadedGzStreambuf` opened a file and inflated it on one thread of its own, into 4 output blocks of 1 MB, which the
reader takes in order. Now, for a BGZF file with `SetBgzfThreads(n)`, n ≥ 2:

- **The dispatcher** (`Dispatch`, the file's own thread) takes a free output block and reads whole BGZF blocks into
  its batch (`ReadBgzfJob`, as `FillBgzf` read them, without inflating) while their sizes from the footers leave room
  for one more 64 KB block. It queues the batch for the workers and, in the same order, for the reader. If a gzip
  member of another kind follows the BGZF blocks (`cat a.bgzf.gz b.gz`), the dispatcher inflates it itself, as before,
  into the next blocks.
- **The workers** (`Work`, n of them) inflate a batch block by block with `bgzf::DecompressBlock` (ISA-L, a state per
  thread, the CRC checked) into its output block.
- **The reader** (`underflow`) takes the next block in file order once it is complete; a worker may finish later
  blocks first. Output blocks: 2 per worker + 2.
- **Errors where the file has them.** A corrupt block ends its batch with the bytes before it. The reader reads those,
  then stops with that block's error. What stops the dispatcher, a cut file or a missing end-of-file block, becomes the
  reader's error once it has read everything before it. The messages are those of one thread.
- **Settings**: `ThreadedGzStreambuf::BgzfThreadsFor(threads, files)` = `threads / (4 × files)`, 1 to 4, set by
  protal before each sample's alignment; `PROTAL_INFLATE_THREADS=N` overrides it. With one thread, the old code path
  runs unchanged. Other gzip files, zstd and plain files are read on one thread as before.

## Checks

| check | result |
|---|---|
| unit tests, whole suite | 460 passed, 3 skipped (as before) |
| new: `BgzfOnSeveralThreadsReadsAsOnOne` | empty, small, 12 MB, 60 BGZF files concatenated (end-of-file blocks inside), BGZF followed by two gzip members; 1, 2, 3, 4, 7 threads: the content |
| new: `ADamagedBgzfFileReadsAlikeOnAnyThreads` | 11 damaged files: a flipped byte at five places, a damaged block header, no end-of-file block, cut inside a block, at a block boundary, inside a header, a cut gzip member after the blocks; 2, 4, 7 threads give one thread's bytes, failure and message |
| new: `ThreadsTakeEveryPairOnceFromParallelBgzf` | 4 reader threads, 2 files on 3 threads each: every pair once |
| `ClosingBeforeTheEndStopsTheInflatingThread` | now also with 4 threads per BGZF file, closed while blocked or busy |
| protal, pe 500k pairs (BGZF), `-t 4`, 1 against 3 threads per file | SAM records identical (418,464), every other file identical; `--verbose` says "BGZF, inflated on 3 thread(s)" |
| ThreadSanitizer, `ThreadedGzStream.*` (Debug, `-DPROTAL_NO_CLONES`: TSan cannot run `target_clones`' resolvers before it starts) | 18 tests pass, the new ones included; no report touches `ThreadedGzStream.h`. The 35 reports are all in code that `#pragma omp critical` synchronises, which TSan cannot see without a TSan-built libgomp: `SeqReaderPE`/`SE` copied and merged inside `omp parallel`, and the older tests' own `omp parallel` loops |

## Measurements

One file, one reader thread that only takes the bytes (`onefile`, 3 rounds; MB/s of FASTQ):

| | round 1 | round 2 | round 3 |
|---|---|---|---|
| one gzip member (HEAD's path) | 1,463 | 1,291 | 1,216 |
| BGZF, 1 thread | 1,219 | 1,075 | 1,087 |
| BGZF, 2 threads | 2,346 | 2,096 | 2,113 |
| BGZF, 3 threads | 3,375 | 2,758 | 2,725 |

The paired-end input path (`readbench`, 2 files, readers that do nothing else, 4 cores; pairs/s, 3 rounds): with 2
readers, uncompressed 3.52-3.59M, one member 2.25-2.30M, BGZF on 1, 2, 3 threads per file 1.52-2.10M, 1.75-2.05M,
1.73-1.97M. Two readers plus two files' inflating threads already fill 4 cores, so more threads per file cannot help
here. With HEAD's reader, readers that also work (`gz_over.sh`) lost 10-25% against uncompressed input, about the CPU
the inflating takes, not stalls.

## What to run on the cluster

The thirteenth-run setup with this build and `--verbose` (measure_performance.sh passes it): the logs' "Input" lines
say whether the r226 sample's files are BGZF. If they are, the sequence reader's time per thread (1.67 s in the AVX-512
run) should fall back to about 1 s, and the aligning by up to the seeding's gain (1.2 s per thread). If they are not:
either `zcat R1.fastq.gz | bgzip -@ 8 > R1.bgzf.fastq.gz` once per file, or rapidgzip as a dependency.
