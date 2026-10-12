# gzip input on several threads with rapidgzip, and ISA-L built with protal

- **Date**: 2026-10-12.
- **Question (the user's)**: would rapidgzip be compatible with ISA-L, or at least similar in performance? In any
  case, implement it, so that it can be compiled statically. (The r226 runs of the sparse-first gate had shown
  that a single-stream gzip file caps a run at one inflating thread per file:
  [report](../2026-10-12-sparse-first/README.md), runs 4-7.) Then, asked which ISA-L rapidgzip should use: the
  patched ISA-L for every build.
- **Code**: the working tree on `85259b0` (uncommitted): `lib/isa-l` (ISA-L 2.32.1 with rapidgzip's changes,
  `lib/isa-l.cmake`), `lib/rapidgzip` (rapidgzip 0.16.0, `lib/rapidgzip.cmake`), `lib/rapidgzip-zlib/zlib.h`,
  `src/IO/ParallelGzip.{h,cpp}`, `src/IO/ThreadedGzStream.h`, `src/RunProtal.h`, the CMake files, CI, the conda
  recipe, the docs, `tests/test_ReaderAndTimers.cpp`.
- **Machine**: WSL2 on the Core Ultra 7 258V, cores 0-3 (`taskset -c 0-3`), Release build.
- **Data**: `~/gzpar/data/R1.member.gz`, 1,681,851,361 bytes of simulated 150 bp Illumina FASTQ (the 2026-10-08
  BGZF benchmark's) as one gzip member of 397,409,724 bytes (4.2x).

## Answer

- **Compatible with ISA-L**: rapidgzip uses ISA-L, but its own fork of it
  ([mxmlnkn/isa-l](https://github.com/mxmlnkn/isa-l), 15 commits on intel's master of August 2025), and needs
  three of the fork's changes: stopping points in `isal_inflate` (it hands the rest of a chunk to ISA-L at a
  deflate block once the chunk's 32 KB window is known), the inflate's Huffman-table functions public (its
  own decoder builds its tables with them), and a `const` dictionary in `isal_inflate_set_dict`. The fork's
  commits did not apply to 2.32.1 (upstream had reformatted `igzip_inflate.c`); the stopping points applied
  with fuzz and the rest were carried over by hand. Built without ISA-L, rapidgzip decodes everything with its
  own decoder and needs zlib's API, which protal's zlib-ng (native API) does not have.
- **Performance on FASTQ: no.** One rapidgzip thread inflates ~0.3 GB/s of this FASTQ, one ISA-L thread
  1.2-1.6 GB/s: rapidgzip spends about four times the CPU per byte. The reason is in its statistics: 98.9% of
  the bytes were decoded with markers (the 32 KB before a chunk unknown), so ISA-L almost never took over. A
  FASTQ record's name and qualities are copies of the record before, which were copies of the one before it,
  back to the chunk's start: the unknown window never drops out of the references. (rapidgzip's own README
  calls FASTQ "almost the worst case".) It is faster than one ISA-L thread only from about 6 threads per file.
- **Statically compiled: yes.** rapidgzip is header-only; it is built with ISA-L (`lib/isa-l`, now linked into
  every binary) and zlib-ng (its zlib calls are compiled, not called; `lib/rapidgzip-zlib/zlib.h` maps them to
  zlib-ng's native API). The static binary was built and tested (below).

## What was done

- **`lib/isa-l`**: ISA-L 2.32.1 (the release `just static` built so far) without `doc/`, `examples/`, `programs/`,
  `tests/`, `tools/`, with the three changes (marked "protal" in the code, listed in `lib/isa-l.cmake`):
  1. the fork's commit 18660e0 (stopping points), applied with fuzz and checked hunk by hunk; 2.32.1 had added a
     second check for corrupt data after the copy into `tmp_out`, which now ends through the fork's restore of the
     input (two copies of the check came out of the fuzzy apply, one dead);
  2. `set_codes`, `set_and_expand_lit_len_huffcode`, `make_inflate_huff_code_lit_len`, `make_inflate_huff_code_dist`
     public with const counts, `struct huff_code` moved to `igzip_lib.h` (the fork's 68027bb, part of 3d91c4d);
  3. `isal_inflate_set_dict(..., const uint8_t*)`.
  Deflate is untouched. `lib/isa-l.cmake` builds it with upstream's CMake for every binary (static, PIC; nasm from
  the PATH, else nasm 2.16.03 downloaded and built, as `lib/static-deps.cmake` did for the static build), so no
  system ISA-L is used: Ubuntu's `libisal-dev` is replaced by `nasm` in CI, `isa-l` by `nasm` in the conda recipe and
  `envs/protal-db-build.yaml`; `lib/static-deps.cmake` builds zstd only.
- **`lib/rapidgzip`**: `src/{core,filereader,huffman,indexed_bzip2,rapidgzip}` of the rapidgzip 0.16.0 source
  distribution (PyPI, sha256 `8b124f29…c14e`), unmodified, 1.5 MB; `lib/rapidgzip.cmake` makes it an INTERFACE target
  with `LIBRAPIDARCHIVE_WITH_ISAL`, without rpmalloc. `lib/rapidgzip-zlib/zlib.h` gives its (uncalled) zlib code
  zlib-ng's native functions, as inline functions (rapidgzip has a namespace `deflate` and calls `::deflate`).
- **`src/IO/ParallelGzip.{h,cpp}`**: the only file that includes rapidgzip (`ParallelGzipReader`, 1 MiB pieces,
  CRC-32 on, index windows not kept). rapidgzip reads a gzip file cut inside a member as a shorter one, often
  without a word, sometimes as an empty one (below), so `ParallelGzip` checks the pieces as they come (rapidgzip's
  read callback hands over each piece with its gzip footers): each piece must begin where the one before ended in
  the file, each member's bytes must number its ISIZE, and at the end the bytes must end with a member whose
  trailer ends the file. rapidgzip itself checks the CRC-32s (only while pieces follow each other, which this checks).
- **`ThreadedGzStreambuf`**: a regular file that starts as gzip but not BGZF, with more than one thread for gzip
  (`SetGzipThreads`, protal: `GzipThreadsFor`), is read through `ParallelGzip` on that many threads: the file's own
  thread copies rapidgzip's bytes into its 1 MB blocks (`FillParallel`), the rest is as before. A pipe, a file
  ending in four zero bytes (padding, which rapidgzip fails on; an empty last member), or a file rapidgzip cannot
  open is read by ISA-L on one thread. `--verbose` says "gzip (not BGZF), inflated on N thread(s) by rapidgzip".
  `PROTAL_INFLATE_THREADS` sets the threads of BGZF and other gzip files alike (for comparing runs).
- **The default** (`GzipThreadsFor`): a thread per 2 alignment threads and file, if that is 8 or more, at most 16;
  else 1 (ISA-L). So a paired-end run on 32 threads or more, or a single read file on 16 or more; every smaller run
  reads gzip as before.

## Correctness

- **Tests** (`tests/test_ReaderAndTimers.cpp`): rapidgzip on 2, 3, 4 and 7 threads gives ISA-L's bytes for one
  member of 12 MB of random FASTQ (~9 pieces) at levels 6 and 1, stored blocks, a header with every optional field,
  25 members of random sizes, gzip followed by BGZF, a tiny file; zero padding and an empty last member go to ISA-L
  and read the same. Cut in half, in the trailer, without trailer, corrupt at three places, a wrong CRC-32, a wrong
  ISIZE, a damaged second header, data or one byte after the last member: all fail on 2 and 4 threads, a cut file's
  bytes a prefix of the content. 60 cut files (100 KB to 5 MB of FASTQ, gzip levels 1, 6, 9, cut at 0-99.9%) all
  fail on 2 and 4 threads. A pipe stays on one thread; read pairs from two files on rapidgzip's threads each come
  once; closing a stream early stops rapidgzip's threads.
- **Before the checks** (`scripts/rgz_trunc.sh`, `bench_gz` on 3 threads,
  [results/trunc_before_checks.txt](results/trunc_before_checks.txt)): of 72 cut files (100 KB to 20 MB of the FASTQ
  above, levels 1, 6, 9, cut at 25%, 50%, 90%, 99.9%), rapidgzip read 62 as whole files without an error, 17 of them
  as empty files (all 12 cut files of 100 KB); ISA-L on one thread failed all 72. After the checks, none is read as
  whole. The end-to-end test of a cut `.sam.gz` (with `PROTAL_INFLATE_THREADS=3`) had found it.
- **The suites**: 537 unit tests; the 163 end-to-end tests on the mini database as usual, with
  `PROTAL_INFLATE_THREADS=3` (every gzip input, FASTQ and SAM, through rapidgzip, every BGZF on 3 workers), and with
  the static binary (`build-static/protal_0.8.2_static`, `-DPROTAL_STATIC_FETCH_DEPS=OFF` with Ubuntu's
  `libzstd.a`; "statically linked", `ldd`: not a dynamic executable) under `PROTAL_INFLATE_THREADS=3`: all pass
  (`scripts/rgz_full.sh`, `rgz_round3.sh`, `rgz_final.sh`). Before the checks, the end-to-end run with
  `PROTAL_INFLATE_THREADS=3` failed one test: a cut `.sam.gz` profiled with a SAM parse error instead of "the file
  is truncated or corrupt".
- **What differs from ISA-L on one thread**: messages ("corrupt gzip data (corrupt or truncated file?): ..." with
  rapidgzip's or `ParallelGzip`'s words); a corrupt member may hand over other bytes before the error than ISA-L
  does (both fail); rapidgzip prints a line of its own to stderr on some damaged files ("Unexpected end of file when
  getting block at ..."). The outputs of intact files are the same.

## Speed

Through `ThreadedGzIstream` (`scripts/bench_gz.cpp`: read the whole file in 1 MB pieces), the FASTQ above, cores
0-3, three alternated rounds ([results/threadedgz_bench.txt](results/threadedgz_bench.txt), 4 MiB pieces, before the
checks):

| inflating | GB/s of FASTQ | CPU s | peak memory |
|---|---|---|---|
| ISA-L, 1 thread | 1.38-1.60 | 1.1-1.3 | 10 MB |
| rapidgzip, 2 threads | 0.56-0.66 | 5.1-6.1 | 630 MB |
| rapidgzip, 3 threads | 0.83-0.94 | 5.4-6.0 | 670-700 MB |
| rapidgzip, 4 threads | 0.87-0.96 | 6.4-7.2 | 730-760 MB |

rapidgzip alone (`scripts/bench_rgz.cpp`, [results/rapidgzip_direct.txt](results/rapidgzip_direct.txt)): 1 thread
1.18-1.22 GB/s (it then decodes in order, with ISA-L), 3 threads 0.81-0.98, 4 threads 1.04-1.17; `-march=native`
the same as protal's `-march=x86-64`; reads of 1-16 MiB the same; pieces of 1 MiB 5% slower than 4 MiB at a quarter
of the memory (200 against 780 MB on 4 threads), 8 MiB no faster at 1.3 GB. Its statistics: 98.9% of the bytes were
decoded with markers, 1.1% without; ISA-L took over almost nothing. One rapidgzip thread does ~0.3 GB/s of this
FASTQ: it needs about 5-6 threads to match one ISA-L thread, and the 4 cores here cannot show more.

The final code (1 MiB pieces, `ParallelGzip`'s checks; `scripts/rgz_final.sh`,
[results/final_bench.txt](results/final_bench.txt)), three alternated rounds, the same bytes in every run:

| inflating | GB/s of FASTQ | CPU s | peak memory |
|---|---|---|---|
| ISA-L, 1 thread | 1.39-1.61 | 1.1-1.3 | 10 MB |
| rapidgzip, 2 threads | 0.60-0.66 | 5.1-5.6 | 171-176 MB |
| rapidgzip, 3 threads | 0.76-0.96 | 5.2-6.6 | 186-192 MB |
| rapidgzip, 4 threads | 0.87-1.08 | 5.4-7.0 | 206-219 MB |

The checks cost nothing measurable; the smaller pieces took the memory from ~750 to ~210 MB at 4 threads.

## What to make of it

- **For FASTQ, BGZF stays the answer**: protal inflates BGZF on 4 ISA-L threads per file at ~1 GB/s each, at
  ISA-L's CPU; rapidgzip needs about four times the CPU for the same bytes and pays only with many threads. Turning
  the reads into BGZF once (`zcat | bgzip -@ 8`) costs a compression of the file.
- **rapidgzip is the fallback for gzip files as they come**, on many cores: with the default it starts at 8 threads
  per file (a paired-end run on 32 threads), where ~0.3 GB/s a thread would give ~2.4 GB/s per file against ISA-L's
  ~1.2-1.6, at ~50 MB per thread. That is extrapolated from 4 cores: rapidgzip's own README says FASTQ fails to scale
  above 64 cores, and the cluster has not run it.
- **ISA-L built with protal** is the cost of rapidgzip at full speed (its fork's stopping points and Huffman
  functions): nasm instead of `libisal-dev`/`isa-l` for building from source, and the same ISA-L in every build.

## Next run

On the EPYC node, the paired-end gzip sample (not BGZF), 32 threads, gated and not: as it is (rapidgzip on 8 threads
per file), with `PROTAL_INFLATE_THREADS=1` (ISA-L on one), and with `PROTAL_INFLATE_THREADS=16`, beside the BGZF copy
(run 6: aligning 23.6 s plain, 11.6 s at C = 2):

```bash
sbatch --exclusive -c 32 --mem 100G --wrap "bash scripts/measure_performance.sh --test_sparse runs/rgz8 DB.protal pe:R1.fq.gz:R2.fq.gz"
sbatch --exclusive -c 32 --mem 100G --wrap "PROTAL_INFLATE_THREADS=1 bash scripts/measure_performance.sh --test_sparse runs/rgz1 DB.protal pe:R1.fq.gz:R2.fq.gz"
sbatch --exclusive -c 32 --mem 100G --wrap "PROTAL_INFLATE_THREADS=16 bash scripts/measure_performance.sh --test_sparse runs/rgz16 DB.protal pe:R1.fq.gz:R2.fq.gz"
```

## Files

- `scripts/`: `rgz_diff.sh`, `rgz_diff2.sh` (the fork against 2.32.1), `rgz_patch.sh`, `rgz_patch3.sh` (the fork's
  commits on 2.32.1), `rgz_vendor_isal.sh` (lib/isa-l), `rgz_build.sh`, `rgz_full.sh` (build and tests),
  `bench_gz.cpp`, `rgz_bench.sh` (through ThreadedGzIstream), `bench_rgz.cpp`, `rgz_bench2.sh`, `rgz_bench3.sh`
  (rapidgzip alone), `rgz_trunc.sh` (cut files), `rgz_round3.sh`, `rgz_final.sh` (the checks after the fixes).
- `results/`: `threadedgz_bench.txt` (4 MiB pieces, before the checks), `rapidgzip_direct.txt`,
  `trunc_before_checks.txt`, `final_bench.txt`.
