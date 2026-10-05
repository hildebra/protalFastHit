# zstd-compressed read files: protal reads .fq.zst

2026-10-05, branch `audit-fixes` at `f3e6e87` plus this change (committed with this report). The user asked for a zstd FASTQ reader after the
collector's report ([2026-10-05-collector-profiling](../2026-10-05-collector-profiling/README.md)) measured zstd
against the gzip the simulators write: at its fast levels ~7% smaller and several times faster to write and read.
Until now protal read plain, gzip and BGZF files; zstd needed a pipe (`-1 <(zstd -dc a_R1.fq.zst)`), and a `.fq.zst`
given as a file failed its sample as "not FASTQ".

## What changed

`src/IO/ThreadedGzStream.h`: the input stream every read path uses (paired-end, single-end and long reads,
`--benchmark_alignment`, the read-type detection of a single read file, the long-read check of single-end samples)
decompresses zstd in its own thread, as it inflates gzip:

- **Told by content**: in the state that looks at a file's first bytes (the one that tells gzip from plain), a zstd
  frame's magic number, or a skippable frame's (pzstd, the seekable format), switches to `FillZstd`. The name does not
  matter (a zstd file named `.fq` reads), and pipes and process substitution work as for gzip.
- **Streaming decompression** with `ZSTD_decompressStream` into the stream's 1 MB blocks, in the thread that
  inflated gzip: the alignment threads still only copy whole lines under the reader lock (`TakeLines`). Frames one
  after the other (`cat a.zst b.zst`) and skippable frames are read; the window limit is raised to libzstd's maximum,
  as the database reader's (`zstd::IStreambuf`), so files of `zstd --long=28` or more read (libzstd's default stops
  at 2^27).
- **Errors as for gzip**: a file that ends inside a frame ("the file ends inside zstd frame N (truncated file?)"),
  a corrupt frame (libzstd's reason "in zstd frame N (corrupt file?)", e.g. a checksum mismatch), or data after the
  frames that is not a frame read as a prefix of the reads and fail the sample ("The FASTQ files of sample ... are
  truncated or corrupt"). A file cut exactly between frames reads as complete, as one cut between gzip members.
- One file is decompressed by one thread. libzstd decompresses 1.2-1.4 GB/s of FASTQ on a core (below), well above
  what the alignment takes from a file, so the seekable format's frames are not decompressed in parallel here (the
  SAM reader does that: `zstd::ParallelFrameStreambuf`).

Also `docs/running.md` (the read formats), and the comments of `Options::LongestRead` and `SampleReads`.

## Tests

Built in WSL from `git archive f3e6e87` with the changed files (`~/zstdreader`, Release, GCC 13, system libzstd
1.5.5):

- **Unit tests** (`tests/test_ReaderAndTimers.cpp`, `tests/test_InputValidation.cpp`): new `ReadsAZstdFile` (one
  frame; several one after the other; skippable frames before, between and after; an empty frame; a frame of window
  2^28, which libzstd's defaults refuse with `frameParameter_windowTooLarge`, as the test checks; a zstd file named
  `.fq`), `ACutOrCorruptZstdFileIsReported` (cut inside a frame, the checksum cut, a corrupt frame, data after the
  frames, a cut between frames); zstd added to `ReadsAPipe` (a FIFO), to `SeqReader.FastqRecordsAreTheSameFromAnyStream`
  (records equal to those of a string stream) and to `ReadTypeDetection.LongReadsByTheirQuality`. **365 of 365 pass**
  (`ctest`).
- **End-to-end** (`tests/e2e/test_protal_e2e.py`, on a mini database built by that binary): new `test_zstd_reads`
  (mates as two zstd frames each, as files and from FIFOs: the same profile as the plain files);
  `test_corrupt_zstd_reads_fail_their_sample` (the zstd magic and random bytes, which the old test used as reads that
  are not FASTQ: now "truncated or corrupt ... zstd frame 1"); `test_reads_that_are_not_fastq_fail_their_sample` now
  with bzip2 bytes. **133 of 133 pass.**

## Speed (`bench_read_formats.cpp`, `run_bench.sh`, `results/`)

The first reads of 500,000 simulated 150 bp pairs (`simulate_metagenomes`, HSXt), three times over: 477 MB of FASTQ.
The benchmark takes the records as protal's readers do (`TakeLines`, 32 records at a time) while the stream's thread
decompresses; median of 5 runs, the files in the page cache, the laptop idle (load ~1):

| File | Size | Read in | FASTQ per second | CPU |
|---|---:|---:|---:|---:|
| plain | 477 MB | 0.057 s | 8.4 GB/s | 0.10 s |
| BGZF (the simulator's, libdeflate level 6) | 119 MB | 0.320 s | 1.49 GB/s | 0.36 s |
| gzip -6 (one member: zlib-ng) | 116 MB | 0.525 s | 0.91 GB/s | 0.56 s |
| zstd -1 | 124 MB | 0.400 s | 1.19 GB/s | 0.44 s |
| zstd -3 | 118 MB | 0.405 s | 1.18 GB/s | 0.44 s |
| zstd -9 | 113 MB | 0.345 s | 1.38 GB/s | 0.39 s |
| zstd -19 | 89 MB | 0.334 s | 1.43 GB/s | 0.38 s |

zstd reads as fast as BGZF (within 20%) and 30-55% faster than one gzip member; all of them hand over reads far
faster than protal's alignment threads take them, and in a thread of their own, so the format does not change
protal's speed. What zstd changes is the size at high levels (-19: 25% smaller than BGZF for
these reads; at -3 the same), and the writing speed (zstd -3 ~340 MB/s against libdeflate's BGZF; the collector
report's table).

## Not done

- bzip2 and xz still go through a pipe.

The website's page on running protal lists the read formats: it needs zstd added.

## Follow-up (2026-10-05): the database build's simulations write zstd by default

The user then asked for the build's simulations to write zstd by default.

- `simulate_metagenomes --reads_compression bgzf|zstd` (its own default stays `bgzf`, `.fq.gz`): with `zstd` each
  sample's reads go to `_R1.fq.zst`/`_R2.fq.zst` as one zstd frame at level 3 with a checksum, no long window
  (`zstd::OStream`, compressed as each genome's reads are appended as before: the same bytes for any number of
  threads); the protal map names them, and its prefixes strip `.zst`.
- `collect_training_data.py --read_compression zstd|gzip` (default `zstd`), and `build_gtdb_database.py
  --read-compression` (default `zstd`), which passes it on: the simulator's reads, the long and Ultima reads
  (`hifi_reads.py`, the pbsim3 renaming), the host's paired-end chunks (zstd frames appended to the sample's zstd
  file, as gzip members were to its BGZF file) and the chunks' templates are zstd (level 3, the templates level 1).
  `scripts/compressed.py` writes by the file name and reads by the first bytes, with Python 3.14's
  `compression.zstd` or else the `zstd` command (in the protal-db-build environment and on CI). A long-read point's
  key gets the compression, and a paired-end point's the simulator command with its option, so a rerun after the
  change simulates its points again.
- **The same reads**: the A/B of the collector report (`ab_collector.py`, now comparing files of either compression
  decompressed) gives the same 21 samples (Ultima, HiFi, Nanopore, chunked or not, a host share) as `22c445d`'s gzip
  collector, both with zstd (the default) and with `--read_compression gzip`; the zstd files took 197 MB against
  216 MB of gzip (8.5% less). In the end-to-end test, `simulate_metagenomes --reads_compression zstd` writes the same
  FASTQ as its BGZF files decompressed, the same bytes on 1 and 3 threads.

Tests (WSL, `git archive e7b391e` + this change): `ctest` 366 of 366; e2e `SimulatorTest` and `BadInputTest` 11 of 11
(new: zstd reads of `simulate_metagenomes` against its BGZF ones); `test_mini_db.py` 68 of 68 with Python 3.12 (the
`zstd` command), `GtdbBuildTest` included (the scenario build's host and Ultima reads are `.fq.zst`; the build
followed in blocks gives the same tables as the one protal run); the long-read and host tests also with Python 3.14
(`compression.zstd`). New in `test_mini_db.py`: the default writes zstd, gzip (`--read_compression gzip`) gives the
same reads, the host's chunks in zstd hold the reads of its gzip ones.
