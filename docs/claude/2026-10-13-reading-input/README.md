# Reading the input: a log line per sample, read-ahead, the piece decoder's memory, threads against CPUs

2026-10-13, on `0ab85f4` (the member reader of the [gzip-members report](../2026-10-13-gzip-members-regression/README.md)).
Four follow-ups to that report, items 1, 3, 4 and 5 of its last list:
1. **Say in the log how the input was read.**
3. **Read files ahead for cold reads.**
4. **Bound the piece decoder's memory.**
5. **Measure threads against CPUs.**

The code is done for all four. The speed measurements of sections 3-5 are rough: from about 18:30 other sessions'
builds and tests ran on the same four cores, and the timing was stopped there. Later that evening the tests were rerun
on a quiet machine, at the committed `a758a4a` ([Speed tests on a quiet machine](#speed-tests-on-a-quiet-machine-2026-10-13-evening-at-a758a4a)):
- **Item 4:** the shorter look-ahead is as fast (fastest runs within 2%) for 29-38% less memory.
- **Item 5:** the reading threads cost the alignment nothing measurable.
- **Item 3:** the read-ahead hints neither help nor hurt measurably on a local disk; a network file system is still to
  test.

## Setup

- **Machine.** WSL (Intel Core Ultra 7 258V), `taskset -c 0-3`, protal and the benchmark built from the working
  tree.
- **Data.**
  - The gzip-members report's 1.48 GB of Illumina-like FASTQ (`make_fastq.pl 4000000 7`), as one member (`one.gz`)
    and as 64 KiB members (`m64k.gz`).
  - For whole runs: 2M read pairs of the mini DB's genes (`mini_reads.pl`), as one member and as 64 KiB members.
- **Scripts** in [scripts/](scripts/), results in [results/](results/). The data and the WSL folders (`~/pgzfix`,
  `~/pgzread`) were removed afterwards.

## 1. The input line

After each sample's alignment, every run's log now has one line, e.g. from the mini-DB runs:

    Sample s input: members_R1.fq.gz gzip of 6606 members, inflated by members on 4 threads; members_R2.fq.gz gzip of
    6606 members, inflated by members on 4 threads; the 4 alignment threads spent 0.43 s each getting reads, 6% of the
    alignment's 7.55 s; 4 alignment and 10 reading threads on 4 CPUs, more threads than CPUs

**Each file** (`ThreadedGzStreambuf::files_read`, recorded as the file closes):
- **Format and count:** BGZF blocks, zstd frames, or gzip members, counted as read. For a file read in pieces, the
  count is ParallelGzip's.
- **Threads and path:** how many threads read it, and how: by members, in pieces, or on one thread.
- **Switches:** where it went on by members (from byte B), and where members fell back to one thread.

**The sample:**
- **Waiting:** the "Sequence reader" time (`Statistics::reader_seconds`) per alignment thread, as a share of the
  alignment's time. A quarter or more adds "the input held the alignment back".
- **Thread budget:** the alignment threads, the reading threads (per stream the most its files used, the file's own
  thread counted with its workers), and the CPUs the process may use (`AvailableCpus`, the affinity: a SLURM
  cpuset, taskset).

**measure_performance.sh** turns the line into two `runs.tsv` columns:
- `input`: e.g. `gzip-members:6606:4t,gzip-members:6606:4t`, with `+members@B` or `+one-thread@B` added where a file
  switched or fell back.
- `reader_share`: e.g. `0.06`.

A run with a share of 0.25 or more is flagged on its progress line. Had 0.8.3a written the line, the gzip-members
regression would have shown at once: pe1's line would have given B10's member count, its files read in pieces, and the
threads waiting 67% of the alignment.

Tests: `FilesReadSayHowEachWasRead` (a stream of five formats on 3 threads and on one), and the switch and fallback
checks in the members tests.

## 3. Read-ahead for files read cold

A sample's reads are usually read once, from a network file system, not from the page cache as repeated benchmarks
read them. In 0.8.3a's job the ONT file's first read waited 128 s a thread, about 50 MB/s.

**The hints:**
- **Mapped (ParallelGzip, in pieces):**
  - the mapping is marked sequential (`MADV_SEQUENTIAL`);
  - the pages up to 64 MB past the workers' furthest piece are asked for (`MADV_WILLNEED`), 32 MB at a time,
    beside the release of the pages behind.
- **read() (one thread, BGZF, zstd frames, members):**
  - `POSIX_FADV_SEQUENTIAL` when a regular file is opened;
  - the next 64 MB asked for (`POSIX_FADV_WILLNEED`) every 16 MB read.
- **The switch:** `PROTAL_READ_AHEAD=0` turns all of it off, so one binary can be compared with itself.

**Measuring:** `measure_performance.sh --cold_reads` drops the read files from the page cache before each timed run
(`dd iflag=nocache count=0`: posix_fadvise DONTNEED, no root). Every run then reads its input from the disk, as a
sample's first run does.

**Locally** ([results/cold_reads.txt](results/cold_reads.txt); `0ab85f4` against the hints, 3 rounds, rough):
- **Cold against warm:** a cold read took 0-20% longer than a warm one (the 64 KiB members 10-50%); the virtual disk
  is fast.
- **With and without the hints:** no difference beyond the noise. one.gz on 4 threads read cold in 0.63-0.71 s
  before and 0.66-0.79 s with them; m64k.gz in 0.42-0.58 s and 0.45-0.46 s.
- **What it cannot show:** a local disk shows no network latency, which is what read-ahead is for.

**To do on the cluster:** run `measure_performance.sh --cold_reads` twice, with and without `PROTAL_READ_AHEAD=0` in
the environment.

## 4. The piece decoder's memory

**What holds the memory.** At 4 threads, one file in pieces held ~470 MB once its mapped pages were released (the
gzip-members report). The decoder was instrumented ([scripts/piece_buffer_count.sh](scripts/piece_buffer_count.sh),
[results/piece_buffer_count.txt](results/piece_buffer_count.txt)):

| threads | look-ahead (slots) | buffers made | live at most | pooled at the end | max RSS MB |
|---|---|---|---|---|---|
| 2 | 12 | 28-29 | 14-15 | 13 (177 MB) | 215-244 |
| 4 | 20 | 52-58 | 38 | 33-36 (454-494 MB) | 430-453 |

- **The look-ahead sets it.** Each slot within it holds up to two output buffers of 4-16 MB: a piece decoded ahead,
  and its second pass until it is resolved.
- **Some buffers are made more than once.** Buffers that grew past the pool's 16 MiB cap are freed and made again,
  and the allocator keeps what was freed.

**Changed:** the workers now start `threads + 2` pieces ahead, up to `threads` further. Before, they started
`2 × threads + 4` ahead, up to `2 × threads` further.

Max RSS of one file read in pieces, 3 rounds, `piece_buffers.sh` and at 3 threads `piece_buffers_t3.sh`
([results/piece_buffers.txt](results/piece_buffers.txt), [results/piece_buffers_t3.txt](results/piece_buffers_t3.txt)):

| threads | before (MB) | after (MB) |
|---|---|---|
| 2 | 197-227 | 156-178 |
| 3 | 346-409 | 236-267 |
| 4 | 467-485 | 285-327 |

- **Memory:** -20 to -36% per file. A paired-end run at 32 threads (4 a file) holds ~0.35 GB less.
- **Speed, unproven.** At 2 and 3 threads the new look-ahead was as fast within the noise (3 threads: 0.67-0.86
  against 0.67-0.93 s).
- **4 threads on 4 cores** (oversubscribed with the reading thread) was slower in 2 of 3 rounds (0.71-0.98 against
  0.65-0.70 s).
- **What would prove it.** On the cluster the decoder's threads have cores of their own. A single-member paired-end
  run there is the check (item 2 of the gzip-members report). It reverts in one line (`ParallelGzip::Open`).

## 5. Threads against CPUs

The reading threads come on top of `-t`. At `-t 32` paired-end that is 10 more: 4 workers and the file's own thread
per file.

**The data.** protal on the mini DB at `-t 4` on 4 CPUs, with 1, 2 and 4 reading threads per file, 6 rounds
([scripts/protal_threads.sh](scripts/protal_threads.sh), [scripts/protal_threads_t4.sh](scripts/protal_threads_t4.sh),
[results/](results/)). The machine was loaded at times, so these are rough. Median aligning seconds:

| input | 1 per file | 2 per file | 4 per file |
|---|---|---|---|
| 64 KiB members | 7.63 | 7.66 | 7.76 |
| one member (pieces from 2) | 8.02 | 8.04 | 8.71 |

- **Members: no cost.** Their workers wait for room most of the time.
- **Pieces at 4 a file: ~9% slower.** On CPUs the alignment fills, the speculation takes CPU from the alignment, and
  the threads spent 10-26% of the alignment getting reads (4-6% on one thread).
- **A CPU left free is worse.** `-t 3` with 4 reading threads per file was slower than `-t 4` with either input
  (9.6-12.3 s).
- **No policy change.** The cost shows only where the alignment alone fills every CPU and the input is
  single-member gzip, and the data is rough. The input line now says it ("more threads than CPUs").

**To do on the cluster:** a run with `-t` equal to the job's CPUs, single-member gzip, with
`PROTAL_INFLATE_THREADS=1` against the default. That would decide whether pieces should get fewer threads when there
are no CPUs to spare.

## Speed tests on a quiet machine (2026-10-13, evening, at a758a4a)

The three questions again, at the committed `a758a4a` ([scripts/speed_tests.sh](scripts/speed_tests.sh),
[results/speed_tests.txt](results/speed_tests.txt) and the `speed_*.tsv` rows).
- **When.** Every timed run started only once no other protal or Python job ran and cores 0-3 had been under 10%
  busy for 5 s (`/proc/stat`). So they ran between other sessions' jobs (20:50 to 23:24).
- **Noise.** This machine still varies 10-50% between rounds (Windows moves its vCPUs between fast and slow
  cores), with an occasional 2-3× outlier. So the medians come with the fastest run of each setting, which is the
  steadier comparison here.

**Item 4, the look-ahead:** one member, 5 rounds.

| threads | before: median (fastest) s, RSS MB | a758a4a: median (fastest) s, RSS MB |
|---|---|---|
| 2 | 0.698 (0.693), 235 | 0.709 (0.691), 152 |
| 3 | 0.629 (0.625), 361 | 0.722 (0.637), 255 |
| 4 | 0.693 (0.606), 463 | 0.653 (0.599), 288 |

- **Speed:** the same. The fastest runs are within 2% at every thread count. The medians go both ways (3 threads
  +15%, 4 threads -6%), and the new look-ahead's 3-thread rounds included a 2.2 s outlier.
- **Memory:** -29 to -38%. The new look-ahead stays.

**Item 3, read-ahead:** 5 rounds, median (fastest) s.

| file, threads | cold, hints off | cold, hints on | warm, off | warm, on |
|---|---|---|---|---|
| one member, 4 (pieces, a mapping) | 0.628 (0.614) | 0.663 (0.606) | 0.604 (0.578) | 0.608 (0.590) |
| one member, 1 (read()) | 0.763 (0.756) | 0.795 (0.720) | 0.699 (0.628) | 0.719 (0.644) |
| 64 KiB members, 4 (read()) | 0.440 (0.422) | 0.478 (0.427) | 0.374 (0.358) | 0.375 (0.343) |

- **Cold reads:** 4-27% slower than warm on this disk (medians).
- **The hints:** no gain here. The cold medians are 4-9% higher with them; the fastest runs are the same or faster.
- **What remains:** a local virtual disk has no network latency to hide, which is what the hints are for. Whether
  they pay is the cluster's question (`--cold_reads`, with and without `PROTAL_READ_AHEAD=0`). They stay on, a
  variable away from off.

**Item 5, reading threads against CPUs:** protal at `-t 4` on 4 CPUs, 6 rounds, median aligning seconds and the
reader's seconds a thread.

| input | 1 per file | 2 per file | 4 per file |
|---|---|---|---|
| 64 KiB members | 6.15 (0.40) | 6.01 (0.36) | 6.15 (0.32) |
| one member (pieces from 2) | 6.37 (0.34) | 6.01 (0.75) | 6.35 (1.23) |

- **No cost to the alignment** from the extra reading threads, by members or in pieces, even with every CPU taken by
  the alignment.
- **The earlier ~9%** (section 5) came from the loaded machine.
- **In pieces,** the reader's time rises with the threads (0.34 → 1.23 s a thread), but the alignment does not
  wait for it here.
- **No policy change** is needed.

## Not done

- **The cluster runs above.** Item 3 still needs one (network file system, cold). Items 4 and 5 are settled by the
  quiet tests above as far as this machine can tell.
- **The pool's 16 MiB cap.** Buffers that outgrow it are freed and made again; keeping them, or capping their growth,
  would cut the churn. It needs a speed check first.
