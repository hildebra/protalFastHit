# Why 0.8.3 runs the paired-end sample slower: gzip files of many small members

2026-10-13. Question: the cluster benchmark of protal 0.8.3 (`local/Perf0.8.3a`, commit `9c27756`) ran the
paired-end sample much slower than 0.8.2 did the day before (`local/Perf0.8.2g`, `5222386`). What happened?

Answer: **the new gzip reader** (`5334f2c`, ParallelGzip, gzip that is not BGZF inflated in pieces on several
threads). On files made of many small gzip members it is slower than the single ISA-L thread it replaced. Its
block finder then scans whole 1 MiB slices of the file without finding a block to start a piece at, and the
reading thread does much of that scanning itself. Alignment, profiles and every other stage are unchanged.

## Setup

- **Cluster runs** (copied to [results/](results/)): `scripts/measure_performance.sh`, 32 threads, database
  `protal0.8.1_r226_v24`, samples pe1 (PEARL_AGE B10, two `.fastq.gz`), pb2 (a HiFi barcode) and ont3 (a soil ONT
  run), all gzip and not BGZF by their first bytes.
  - 0.8.2: `5222386`, node q512n10, 2026-10-12 16:25 ([results/cluster_082g](results/cluster_082g/)).
  - 0.8.3: `9c27756`, node q512n11, 2026-10-13 12:58 ([results/cluster_083a](results/cluster_083a/)).
  - Both nodes are AMD EPYC 9634 (84 cores). Another tool timed on the same files took the same time in this job and
    in 0.8.2's earlier 32-thread job on q512n19 (pe1 64.5-66.5 s).
- **Reproduction** in WSL (Intel Core Ultra 7 258V, `taskset -c 0-3 nice -n 5`):
  - `bench_pgz` ([scripts/bench_pgz.cpp](scripts/bench_pgz.cpp)) reads a file through `ThreadedGzIstream` with
    1, 2 or 4 gzip threads. It is built from `git archive 9c27756` against the system ISA-L 2.31
    ([scripts/setup.sh](scripts/setup.sh)).
  - The data is 1.48 GB of Illumina-like FASTQ: 4 M records of 151 bp with binned qualities
    ([scripts/make_fastq.pl](scripts/make_fastq.pl)).
  - It is compressed four ways, all with `gzip -6`: as one member, and cut by `split -b 8M / 1M / 64K` into
    members that are concatenated.
  - Timings: [scripts/run_bench.sh](scripts/run_bench.sh), 2 alternated rounds ([results/bench.txt](results/bench.txt);
    `bench_pgz.cpp` has since gained two columns for the fix below).
  - Instruction profile: [scripts/run_callgrind.sh](scripts/run_callgrind.sh) ([results/callgrind.txt](results/callgrind.txt)).
  - The data was removed afterwards; `setup.sh` remakes it.

## What changed between the two cluster runs

Warm runs (repeats 2 and 3; repeat 1 of 0.8.3 read the database and the ONT file cold from NFS, see below):

| | 0.8.2 (`5222386`) | 0.8.3 (`9c27756`) |
|---|---|---|
| pe1 wall s | 38.0-38.2 | **67.2-68.0** |
| pe1 aligning s | 23.6-24.5 | **58.2-58.3** |
| pe1 "Sequence reader", s summed over 32 threads | 39.8 | **1,210** |
| pe1 aligning cores busy | | 15.9 of 32 |
| pe1 user s / instructions | 944-947 / 6,352 G | 1,078-1,079 / 7,822-7,830 G |
| pe1 max RSS GB | 37.7 | 45.0 |
| pb2 wall s (aligning s) | 12.7-13.8 (7.6-7.7) | 11.6-11.8 (6.7) |
| pb2 max RSS GB | 37.6 | 39.4 |
| ont3 wall s (aligning s) | 51.4-52.0 (44.0-44.2) | 50.5-50.9 (44.2-44.6) |
| ont3 max RSS GB | 38.1 | 45.4-45.5 |

- **The output is the same.** The pe1 records written are identical (2,223,279), and so are the SAM statistics.
  Every profile in `sparse/`, plain and at each sparse-first level of all three samples, is byte-identical between
  the two versions.
- **Only the input stage moved.** In `stages.tsv` (repeat 2) the stages of pe1 are within 5% (seeding 369 / 370 s,
  K-mer screen 44.0 / 44.8 s, ...), the output handler aside (21.8 / 26.2 s), but "Sequence reader" went from 39.8
  to 1,210 thread-seconds. The 32 alignment threads spend most of the run waiting for reads.
- **The gated runs show the floor.** pe1 aligns in 58.0-60.6 s at every sparse-first level (C 1-8), against
  23.1-24.1 s with 0.8.2. A floor that does not move with the work is the input.
- **The long-read samples gained.** pb2 at C 4-32: 4.3-4.6 s aligning (0.8.2: 8.8-10.0 s); ont3: 27.3-29.4 s
  (36.5-39.7 s). On their files the new reader does what it was written for: it lifts the one-thread gzip floor
  that the sparse-first report found.

## Why: the reader on files of small members

The reproduction, seconds to read 1.48 GB of FASTQ (2 rounds; 1 thread is the 0.8.2 path, ISA-L on the file's own
thread):

| file | members | 1 thread | 2 threads | 4 threads | CPU s at 4 |
|---|---|---|---|---|---|
| `gzip -6` | 1 | 1.08-1.59 | 1.09-1.54 | 0.99-1.00 | 3.6 |
| members of 8 MiB of FASTQ | 177 | 1.06-1.08 | 0.93-0.94 | 0.86-0.90 | 3.1 |
| members of 1 MiB | 1,413 | 1.11-1.14 | 1.02-1.10 | 1.20-1.26 | 1.8-1.9 |
| members of 64 KiB | 22,606 | 1.23-1.28 | **4.19-4.20** | **4.39-4.72** | 8.0-8.6 |

With 64 KiB members the reader delivers 0.31-0.35 GB/s on 2 or 4 threads, 3.5× slower than one thread, for 6× the
CPU. On the cluster, 0.8.2's one thread a file gave pe1 its 24 s floor. 0.8.3 takes 58 s, 2.4× longer, with
+130 user s and +1.47 T instructions: the same pattern.

callgrind, the first 64 MiB of the FASTQ as 1,024 members of 64 KiB:
- **1 thread:** 0.60 G instructions, 64% of them ISA-L's inflate.
- **2 threads:** 2.35 G instructions.
  - 1.54 G (66%) are the block finder: `IsDynamicBlockHeader` 1.30 G, `FindBlock` 0.24 G. That is 131
    instructions per byte of compressed file.
  - The inflate is the same 0.38 G.
  - Every member was streamed (`streamed 1024`), none decoded ahead (`ahead 0`).

The mechanism, in [ParallelGzip.cpp](../../../src/IO/ParallelGzip.cpp):
- **Piece starts.** A piece starts at the first *non-final* dynamic block of its 1 MiB slice (`FindBlock`,
  `IsDynamicBlockHeader`: BFINAL 0, BTYPE 2).
- **Small members give none.** A member of 64 KiB of FASTQ is ~11.5 KB of deflate, usually one final block. A slice
  of ~90 such members has no start, so `FindBlock` scans the whole slice, finds nothing, and records `kNone`.
- **Every member goes back to the finder.** The reader's thread has no start within its look-ahead
  (`NextStartAfter` asks for 21 slices), so it streams the member to its end (`Advance`, `streaming`). At the
  member's end (`StreamSegment` → `MemberEnd`) it goes back to `Advance` and asks for the look-ahead again.
- **The reader scans for itself.** Each slice the reader moves into brings a slice at the far end of its look-ahead
  that has not been scanned. The workers look for starts at most 2 slices beyond that (`Work`, step 2), and a slice
  costs ~140 M instructions (~10-15 ms, estimated from the count). So the reader often scans the slice itself
  (`StartOf`) instead of decoding.
- **Nothing is decoded ahead.** The workers have no starts to speculate from, so the reader's thread does all the
  decoding and much of the scanning.
- **Larger members work.** At 1 MiB members, slices hold non-final blocks again and the reader is about as fast as
  one thread. At 8 MiB members it is faster.
- **The report's cost estimate does not hold here.** The [inflate report](../2026-10-12-isal-parallel-inflate/README.md)
  measured ~0.1 s of finder CPU for a whole one-member file, because there each slice's scan stops at its first
  block, ~30 KB in (~3% of the file scanned). On a file without starts it scans every byte, ~35× more, and some
  slices twice (the reader and a worker).

**Is B10 such a file?** Probably, but not checked. The benchmark script tells BGZF from gzip by the first bytes only,
and nothing in the logs counts members; the [sparse-first report](../2026-10-12-sparse-first/README.md)'s "one gzip
stream" means "not BGZF". The cluster numbers match the 64 KiB reproduction:
- the same rate, ~0.33 GB/s a file;
- the extra instructions: +1.47 T is about 11 GB of finder scanning at 131 instructions a byte, one to two passes
  over two files of the size the memory below suggests;
- the long-read files got faster, not slower, so the slowdown belongs to the B10 files, not to the node or the build.

One command on the cluster counts the member headers in the first 64 MiB (a one-member file gives 1):

    perl -e 'open(my $f, "<:raw", $ARGV[0]) or die; my $n = read($f, my $b, 64 << 20);
             my $c = () = $b =~ /\x1f\x8b\x08[\x00-\x1f]/g; printf "%d member headers in %.0f MiB\n", $c, $n / 2**20' \
        /hpc-home/hildebra/data/raw/Illumina/PEARL_AGE/00_fastq/B10_R1_001.fastq.gz

## The memory: the compressed files are mapped whole

Max RSS rose during aligning by ~7.3 GB for pe1 and ont3 and by ~1.8 GB for pb2. It fell back once the reads were
read: pe1's log says "Memory after aligning: 37.59 GB resident, peak 44.99 GB", where 0.8.2 had 37.67 / 37.67.
- **The mapped file.** `ParallelGzip::Open` maps the whole compressed file (`mmap`, `MAP_PRIVATE`). Every page read
  stays in the process's resident set until the file is closed.
- **The reproduction shows it.** Max RSS is ~260 MB for the 260 MB m64k file, against 9-10 MB with one thread. The
  output buffers add 100-350 MB on the other files.
- **The cost.** These are clean page-cache pages, which the kernel can drop. They still count in RSS, in protal's
  "peak", and in a job's cgroup memory. The 7.3 GB fits two compressed FASTQ files of ~3.6 GB (`ls -l` on the
  cluster checks it).

## Not caused by the version

- **pe1 repeat 1** (519 s): the index loaded in 453 s, read cold from NFS. The v24 database had not been used on
  q512n11 before. 0.8.2's first 32-thread job on q512n19 had a 935 s cold load.
- **ont3 repeat 1** (152 s, 62,471 major faults, 128 s a thread waiting for reads): the ONT file was read from NFS
  for the first time on that node. If it is ~7 GB, as the memory rise suggests, that is ~50 MB/s, close to the
  database's cold 66 MB/s in the same job. Repeats 2 and 3 are at 0.8.2's time.

## What to do

1. **Now, for such files.**
   - `PROTAL_INFLATE_THREADS=1` brings back the one-thread reader (and 0.8.2's ~24 s floor for pe1). It also puts BGZF
     and zstd files on one thread.
   - Better: rewrite the files as BGZF once (`bgzip -@ 8`). From BGZF, pe1 at 32 threads had no input floor
     (sparse-first runs 5-7).
2. **Fix the reader (small).**
   - The reader's own thread should never run `FindBlock`: a slice not yet scanned counts as "no start known".
   - A file whose first member ends within a few MiB should not be read in pieces. Its members are independent, so
     they need no window.
3. **Fix the reader (better).** Decode the members of such a file on several threads, as BGZF blocks are decoded:
   - the members are independent deflate streams, so there is no window to recover and one pass each;
   - member starts are found by their header (1f 8b 08 ...) and confirmed by the member before ending there.

   That would read multi-member files about as fast as BGZF.
4. **Memory.** Release the mapped pages behind the reader (`madvise(MADV_DONTNEED)` on the consumed range), or read
   the file with `pread` into a window instead of mapping it whole.
5. **Benchmarks.** `measure_performance.sh` should record each read file's member count beside "gzip"/"BGZF".

## The fix: members on several threads, mapped pages released (2026-10-13, later that day)

Fixes 3 and 4 above are implemented on `9c27756` (working tree, not committed yet). Fix 2 is not.

**Members as BGZF blocks** ([ThreadedGzStream.h](../../../src/IO/ThreadedGzStream.h)):
- **Which files.** A regular gzip file that is not BGZF and whose second member begins within `kMaxGzipMember`
  (16 MiB) of its start (`SecondMemberWithin`, a `pread` of the file's first bytes when it is opened). Such a file is
  read by the BGZF dispatcher and workers (1-4 threads per file, `InflateThreadsFor`). Other gzip files go in pieces
  to ParallelGzip, as before.
- **Cutting.** The dispatcher reads the file with `read()`, not a mapping, and cuts it where a member's header begins
  (`MemberStart`):
  - the magic and the method;
  - flags without reserved bits, extra flags 0, 2 or 4, and a known system;
  - optional fields within bounds, and the header CRC if there is one;
  - a possible first deflate block (not type 3; a stored block's LEN and NLEN complements; HLIT and HDIST at most 29);
  - the trailer's length before it no more than 1,032 times the bytes since the last cut.

  Whole members are batched while their trailers' lengths fit an output block (1 MiB; one larger member grows it, up
  to 64 MiB), as zstd frames are.
- **Inflating.** A worker inflates a batch's members one after the other from the batch's start
  (`InflateGzipMembers`; ISA-L checks each CRC-32 and length). The batch must end exactly where it was cut. A header
  the check misses only puts two members into one batch.
- **When a batch fails.** This covers a header found by chance that cut a member, and a damaged member. The reader
  stops the dispatcher and workers and seeks to the batch's start, which the batch before ended at exactly. It reads
  the rest of the file on one thread (`ReadOnOneThreadFrom`), so a damaged file gives one thread's bytes and error.
- **Long members.** A member longer than 16 MiB, or holding more than 64 MiB, sends the rest of the file to the
  dispatcher's own thread, as a zstd frame without a size does.
- **`--verbose`** says "gzip of several members (not BGZF), its members inflated on N thread(s)".

**Mapped pages released** ([ParallelGzip.cpp](../../../src/IO/ParallelGzip.cpp), `Release`): every 32 MiB, the
reader's thread releases the mapped pages before the point where it decodes next (less 1 MiB) with
`MADV_DONTNEED`. The mapping is private and never written, so a worker that looks back reads the page from the page
cache again.

**Tests.**
- **Unit tests:** 549 pass. New or changed:
  - `GzipMembersOnSeveralThreadsReadAsOnOne`: members of 1-64 KB at levels 0-9, full headers, empty members, BGZF
    after them, zero padding, a member over the limit, and a file whose second member comes too late, which is read
    in pieces.
  - `AGzipMemberHeaderFoundByChanceIsNeverUsed`: fake headers written into a stored member; at least one fallback to
    one thread, and the file's bytes.
  - `ADamagedGzipMembersFileReadsAlikeOnAnyThreads`: the same bytes and the same error as one thread for flipped
    bytes, a damaged header, a wrong CRC-32 or length, a cut, and data after the last member.
  - Pairs taken by four threads from two member files.
  - Closing a member file mid-read.
  - The pieces tests run with members turned off (`SetMaxGzipMemberBytes(0)`).
- **End-to-end:** 175 pass (mini DB rebuilt at the working tree).
- **protal itself** ([scripts/protal_check.sh](scripts/protal_check.sh)): 400,000 pairs on the mini DB at `-t 16`, as one
  member per file and as 64 KiB members. `--verbose` reports pieces and members respectively. The profile, the SAM
  counts and every side file are the same, apart from the sample name.

**Before and after** ([scripts/run_compare.sh](scripts/run_compare.sh) on the files of
[scripts/bench_data.sh](scripts/bench_data.sh), [results/compare.txt](results/compare.txt)):
- **Setup.** bench_pgz (since extended: members and restarts columns) built against `9c27756` and against the fix.
  3 alternated rounds, WSL, cores 0-3. The 1-thread runs of one file vary up to 2× on this machine.
- **Columns.** Seconds to read 1.48 GB of FASTQ (ranges over the rounds), and max RSS in MB at 4 threads.
- **"mixed"** is the file's first 400 MB as one member (67 MB of gzip), the rest in 64 KiB members.

| file | 1 thread (both) | 2 threads before | 2 after | 4 threads before | 4 after | RSS at 4, before → after |
|---|---|---|---|---|---|---|
| one member | 0.86-2.01 | 0.97-1.67 | 1.16-1.78 | 1.14-1.52 | 1.19-1.46 | 565-601 → 424-452 |
| members of 8 MiB | 1.07-1.43 | 1.00-1.28 | 0.64-1.23 | 0.98-1.32 | 0.47-0.77 | 504-532 → 126 |
| members of 1 MiB | 0.92-1.38 | 0.90-1.45 | 0.47-0.81 | 0.97-1.37 | 0.46-0.61 | 292 → 24 |
| members of 64 KiB | 1.06-1.82 | 3.83-6.18 | 0.59-1.45 | 4.14-6.50 | 0.47-1.43 | 258 → 22 |
| mixed | 1.02-2.34 | 3.96-4.80 | 3.82-5.09 | 3.56-5.92 | 3.68-5.77 | 509-538 → 369-394 |

- **Small members.** Files of small members are now read faster than on one thread, with 4-9× less time and CPU
  (1.3-3.5 s against 7.4-12.0 s at 4 threads for 64 KiB members) than in pieces. Their RSS is the dispatcher's
  buffers alone, with no mapping.
- **One member.** Speed is unchanged. RSS falls by ~150 MB of the 249 MB file; the rest is ParallelGzip's output
  buffers.
- **Mixed files.** A file whose first member is larger than 16 MiB is still read in pieces, so the 64 KiB members after
  it are as slow as before. Fix 2 (the reader's thread never scanning for block starts itself) is what would help
  there.
- **Not measured on the cluster.** The cluster run and B10's member count are still to come.

## Mixed files: from pieces to members at a member's end (2026-10-13, evening)

Fix 2 above, as written, would have made the reader's thread treat a slot not yet scanned as one without a start. In a
large member that the workers have not scanned far enough ahead, that thread would then stream the member to its end
alone, so one-member files would lose their pieces. It is done another way.

- **Where it stops.** ParallelGzip, at each member's end (`MemberEnd`), looks for the headers of the next two
  members within the member limit (`gzip_member::Next`). The check that tells member headers apart now lives in
  [GzipMember.h](../../../src/IO/GzipMember.h) for both readers. If both members are there, ParallelGzip stops
  before them (`StoppedAt`, `MembersRead`); every member before has been read and checked.
- **Where it goes on.** The file's own thread, which was copying pieces, takes over (`SwitchToMembers`). It ends
  ParallelGzip, so its threads end and its mapping goes, and seeks to the member. It starts the BGZF workers and
  becomes their dispatcher for the rest, as for a file of members from its start. Damage after that point, and
  members longer than the limit, are handled as there.
- **Blocks.** A file opened in pieces gets output blocks for that many workers (10 at 4).
- **Tests.**
  - `AGzipFileInPiecesGoesOnByMembers`: small members after one or two large members, and a large member again after
    them. Damage after the switch (a flipped byte, a damaged header, a cut, data after the last member) gives one
    thread's bytes and error.
  - Closing such a file before, during and after the switch.
  - Unit 550 and end-to-end 175 pass; the protal check gives the same outputs again.

Before (`9c27756`) and after, the same files and setup as above, 3 rounds ([results/compare_fix1.txt](results/compare_fix1.txt)):

| file | 1 thread (both) | 2 threads before | 2 after | 4 threads before | 4 after | RSS at 4, before → after |
|---|---|---|---|---|---|---|
| one member | 0.70-1.90 | 0.72-1.02 | 0.85-1.07 | 0.77-0.97 | 0.85-1.08 | 597-647 → 456-495 |
| members of 8 MiB | 0.85-1.46 | 0.90-1.06 | 0.66-0.76 | 0.88-0.96 | 0.46-0.49 | 501-508 → 126 |
| members of 1 MiB | 0.84-1.38 | 0.81-1.94 | 0.45-1.80 | 0.91-2.56 | 0.51-0.93 | 292-296 → 24 |
| members of 64 KiB | 1.03-2.13 | 3.35-4.46 | 0.55-0.60 | 3.26-4.15 | 0.38-0.68 | 258 → 22 |
| mixed (400 MB as one member, then 64 KiB members) | 0.78-1.25 | 2.33-2.77 | 0.72-0.88 | 2.34-5.89 | 0.57-2.04 | 495-511 → 363-398 |

- **Mixed files.** They now read faster than on one thread, not 2-5× slower. Their peak memory is the pieces' while the
  first member is read.
- **Noise.** One-member reads changed by no more than this machine's noise (its 1-thread runs vary 2× between rounds).
  The member reads are as in the section above.
