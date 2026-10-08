# A simulated sample's .fq.gz differed between 1 and 4 threads: ISA-L hashes a byte from the address of its state

2026-10-08, `audit-fixes` at `401c4f5` (the C++ of `ec7ddb8`) plus the change: `src/IO/Bgzf.h`,
`tests/test_RunStatusAndBgzf.cpp`, `docs/installation.md`. WSL Ubuntu 24.04 on the shared laptop (Core Ultra 7 258V),
gcc 13.3, the system ISA-L 2.31.0 (`libisal2 2.31.0-0.1build1`, which the test binaries link); ISA-L's 2.32.1 source
(`~/isal/isa-l-2.32.1`) and its GitHub master of 2026-10-08 read for the cause. Every build and run on cores 0-3,
nice 5.

**In short.** The two files held the same reads. One BGZF block was deflated differently, because ISA-L's level-1
compressor (its SSE4.2/AVX/AVX2 assembly, 2.31.0 to its current master) hashes the third byte of every stateless call
from a register it never loads there. That register holds the address of the caller's `isal_zstream`, so bits 16-47
of the address choose one hash bucket per block, and with it now and then a match (in the failing block, the bucket of
its first four bytes: their next occurrence lost its match). The one-thread run compressed the block on the main
thread, whose stream's address happened to hit that bucket (one bucket value in 4096 for this block); the four-thread
runs compressed it on pool threads. protal now maps each thread's stream at an address whose hash is the same for all (`Deflater`,
`detail::MapStream`), so a block gets the same bytes on any thread and in any run. That holds for the simulator's
`.fq.gz` and for protal's `.sam.gz` (both `bgzf::Compress`).

## The failure

`LongReadSimulation.SamplesTemplatesAndThreads` (`tests/test_LongReadSimulation.cpp`) failed once in 50 runs of the
unit suite under load (2026-10-07/08, `~/pipefix/runs/fixed_suite.39.log`, binary `~/pipefix/bin/tests_fixed` built
from `0fe84e3` plus the three files of `973e2fb`). Line 278 compares sample s2 (`s2.fq.gz`, BGZF, 150,000 bases,
genomes GB then GA) written with 1 thread against 4 threads; line 288 compares the two genome-store runs (4 threads)
against the 1-thread file. The three 4-thread files agreed, the 1-thread file was the odd one.

## What differed

gtest prints both strings whole (escaped). [`unescape.pl`](scripts/unescape.pl) turned them back into the two files
and [`compare_s2.sh`](scripts/compare_s2.sh) compared them ([results](results/probes.txt)):

- Both are valid gzip and hold the same FASTQ: 302,899 bytes, 151 reads, the same MD5. The reads did not differ.
- Of the nine BGZF blocks, blocks 0-5 are identical. Block 6 (3,580 bytes of content at offset 294,045: the whole
  piece of round 2, genome g0, reads 148 and 149) is 1,905 bytes in the 1-thread file and 1,901 in the 4-thread
  files; blocks 7 (a piece of 5,274 bytes) and 8 (end of file) are the same bytes, 4 bytes earlier.

[`deflate_dump.pl`](scripts/deflate_dump.pl) decoded block 6 of both into deflate tokens
([diff](results/block6_token_diff.txt)). The two streams make one different decision, at input position 1859, the
second read's name `@g0x_149`: the 4-thread stream matches the first read's name at position 0 (`M 7 1859`,
`@g0x_14`); the 1-thread stream misses it and writes `@`, `g`, `0` as literals, then a 4-byte match. Everything else
follows from that (the Huffman tables, 4 bytes).

Recompressing block 6 ([`probe_deflate.cpp`](scripts/probe_deflate.cpp)) always gave the 4-thread bytes (1,875 bytes
of deflate): with fresh Deflaters, after the sample's other blocks, with any bytes after the input (zeros, 0xff, the
file's next bytes, the input's own tail, random letters). Valgrind's memcheck found no read beyond the input and no
uninitialised value. So neither the input buffer's tail (the piece's `std::string` capacity) nor the level buffer's
history mattered.

## The cause: a register ISA-L does not load

ISA-L's level 1 (and 2) body is `igzip_icf_body_h1_gr_bt.asm` (`isal_deflate_icf_body_hash_hist_01/02/04`, chosen
for SSE4.2, AVX and AVX2; only the plain C `_base` version, for CPUs without SSE4.2, differs). It hashes two
positions per step: `hash` for `f_i` and `hash2` for `f_i + 1`. At the start of a stream without history, which
every `isal_deflate_stateless` call is (`reset_match_history` sets `IGZIP_NO_HIST`), it jumps to
`.write_first_byte`, which writes byte 0 as a literal and then prepares position 2's hash:

```
    mov    hash, hash2
    shr    tmp2, 16
    compute_hash  hash2, tmp2     ; crc32 of the 4 bytes at position 2 - if tmp2 held the input
```

On the way there `tmp2` (`rcx`) is never loaded with the input; the loop loads it (`mov tmp2, curr_data`) before the
same shift, and the level-0 body's `.write_first_byte` (`igzip_body.asm`) has its own `mov tmp6, curr_data`. On
Linux the function starts with `mov rcx, rdi`: `rcx` is the `isal_zstream` pointer. So position 2's "hash" is the
CRC-32C of bits 16-47 of the stream's address, masked by the block's hash mask (8191 for blocks of 4 KB and more;
`(1 << bsr(size)) - 1` below, 4095 for block 6). Position 2 is written into that bucket `G` (and its own bucket is not
updated). The output stays valid (every match is checked byte by byte), but it depends on `G`. The code is the same in
2.31.0, 2.32.1 (the latest release, 2026-07-01) and the GitHub master of 2026-10-08.

Two experiments confirm it on block 6:

- [`probe_address.cpp`](scripts/probe_address.cpp): the same call with the `isal_zstream` at 8,192 addresses 64 KB
  apart. 8,188 gave the 4-thread bytes, 4 gave exactly the 1-thread bytes.
- [`probe_buckets.cpp`](scripts/probe_buckets.cpp): one place per value of `G` (0-4095). Only `G` = 3586 gives the
  1-thread bytes, and 3586 is the bucket of the block's first four bytes (`@g0x`, position 0 and position 1859): `G`
  overwrote position 0 in its bucket, and position 1859 found position 2 there instead.

**Why the 1-thread run, and why the 4-thread runs agreed.** A Deflater is `thread_local`
(`bgzf::Compress`). With 1 thread the main thread compressed every block; its stream's address (fixed for the
process, chosen by ASLR) happened to give `G` = 3586 under mask 4095. In the 4-thread runs a pool thread compressed
block 6, and any of the 4095 other values gives the common bytes. A run of the test can fail on any of s2's eight
blocks whose first bytes come back within them; for each, one bucket value in 512 to 8192 (by the block's size) does
it, for the main thread's stream or for a pool thread's.

**Not involved:** the pipeline (the reads are the same; `Plan` runs once a round's items are all in; `Make` depends on
the item and the genome only), timing and load (the address alone decides; the load only accompanied the runs that
showed it), the level buffer's history, the bytes after the input.

## The fix

`protal::bgzf::Deflater` ([`Bgzf.h`](../../../src/IO/Bgzf.h)) no longer puts its `isal_zstream` on the heap.
`detail::MapStream` maps it (`mmap` with `MAP_FIXED_NOREPLACE`) at an address whose bits 16-47 have the CRC-32C
`kStreamHash` (0) under level 1's hash mask (`IGZIP_LVL1_HASH_SIZE - 1`, 8191; every block's mask is a part of it). It
looks 64 KB apart downward from where the kernel maps now: one place in 8192 qualifies, so one about every 512 MB of
address space, and a place in use is skipped (a freed one is used again). Every stream then hashes position 2 into
bucket 0, on any thread and in any process, so a block's bytes depend on its content alone. Where no place is found
within 64 GB, or `mmap` fails otherwise, the stream goes on the heap as before: valid output that may differ by
thread. `Crc32cWord` (a table CRC-32C, ISA-L's hash) matched the `crc32` instruction on 1,000,000 random words.

The bytes are not those of a fixed ISA-L (which would hash position 2's own bytes), and they differ from the old output
in the few blocks where bucket 0 or the old `G` mattered. The content is the same.

Cost ([`probe_fix.cpp`](scripts/probe_fix.cpp)): 0.5 ms to place each of 20 streams alive at once, 1.7 ms each for
100 (each search passes the places in use), about 0.1 ms when one is freed and the next created (the place is used
again). A Deflater is made once per thread.

The placement depends on ISA-L's internals (the register holds the stream's address; the hash is CRC-32C; the level-1
mask). If ISA-L fixes the bug, the placement does nothing. If it changed so that something else chose the bucket, the
new unit test would fail.

**The test.** `BgzfWriter.BlocksDoNotDependOnTheThread` ([`test_RunStatusAndBgzf.cpp`](../../../tests/test_RunStatusAndBgzf.cpp))
deflates 3,000 blocks of two FASTQ records named alike (`@XYZ_1`, `@XYZ_2`; random names, ~222 bytes, so mask 255 and
every bucket holds some block's first four bytes) with four Deflaters (streams at four addresses) and, through
`bgzf::Compress`, on the main thread and three other threads; all must be the same. On the old `Bgzf.h` it failed every
comparison (23, 20 and 25 blocks differed between Deflaters; all three threads differed); with the fix the outputs
are identical. It takes about 0.1 s (92 ms under load). `LongReadSimulation.SamplesTemplatesAndThreads` itself needed
no change: its reads never depended on the threads.

## Verification under load

[`loops_under_load.sh`](scripts/loops_under_load.sh), then [`loops_final.sh`](scripts/loops_final.sh) with the
final files (the comment in `Bgzf.h` reflowed; the tree checked equal to the working files). Binaries from
[`build_head.sh`](scripts/build_head.sh) (`401c4f5`) and [`build_variants.sh`](scripts/build_variants.sh) (the new
test on the old header; the fix), run on cores 0-3 at nice 5 while `protal_tests` was rebuilt from scratch (ccache
off, `ninja -j4`) on the same cores, round after round. Every run is a new process (new addresses), except the `--gtest_repeat` line. Load average 6-10 (other
sessions' niced analyses ran on the other cores).

| Binary | Test | Runs | Failed |
|---|---|---|---|
| new test, old `Bgzf.h` | `BgzfWriter.BlocksDoNotDependOnTheThread` | 20 | 20 |
| fix | `BgzfWriter.*` | 20 | 0 |
| `401c4f5` (old) | `LongReadSimulation.SamplesTemplatesAndThreads` | 500 | 2 (runs 53 and 316: s2 only, the 1-thread file against the three 4-thread ones, as in the original failure, other blocks) |
| fix | `LongReadSimulation.SamplesTemplatesAndThreads` | 500 | 0 |
| fix | the same, one process, `--gtest_repeat=500` | 500 | 0 |
| fix | the whole unit suite (452 tests, 2 skipped by design: ASan and debug builds only) | 5 | 0 |
| final files | `BgzfWriter.*` | 20 | 0 |
| final files | `LongReadSimulation.SamplesTemplatesAndThreads` | 2,000 | 0 |
| final files | the whole unit suite | 2 | 0 |

The old binary's rate (2 in 500 processes, 0.4%) fits the cause: s2 has eight data blocks, each changed by one bucket
value in 512 to 8192, for the main thread's stream or a pool thread's. At that rate 0 failures in 500 runs would still
happen by chance 13% of the time; 0 in the 2,500 runs of both rounds would happen about once in 22,000 (0.996^2500).
The direct evidence is the new test (fails every time on the old header) and the address and bucket experiments.

## What else this touches

- **protal's `.sam.gz`** goes through the same `bgzf::Compress` on its output threads: its bytes could also differ
  between runs and thread counts; now they do not.
- **Earlier files:** a `.fq.gz` or `.sam.gz` written before the fix can differ from one written now in a few blocks;
  their content is the same.
- **Docs:** [`installation.md`](../../installation.md) said the gzip is the same byte for byte for one ISA-L version on
  SSE4.2-AVX2 CPUs; it now adds "whichever thread compresses it". The website does not cover this.
- **Upstream:** the one-line fix in ISA-L is `mov tmp2, curr_data` before `shr tmp2, 16` in `.write_first_byte` of
  `igzip/igzip_icf_body_h1_gr_bt.asm`. Not reported to ISA-L (the user's decision).

## Commands and data

- The failing log: `~/pipefix/runs/fixed_suite.39.log` (WSL); lines 1004-1037 are the failure of line 278
  (`sed -n 1004,1037p` gives the excerpt the scripts read as `fail278.txt`).
- Rebuilt files and decoded blocks: `~/lrdet/fail39` (WSL); the probes' builds: `~/lrdet/fix`.
- The build tree and binaries: `~/bgzfdet` (`bin/tests_head`, `bin/tests_newtest_oldcode`, `bin/tests_fixed`,
  `bin/tests_final`); failing runs' logs in `~/bgzfdet/runs` and `runs_final`; `summary.txt`, `summary_final.txt`,
  `load_rounds.log`, `load_rounds_final.log`.
- The scripts ran from the session's scratch folder (`$SP` in them); copy them and `fail278.txt` together to run them
  again.
