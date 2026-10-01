# Performance, round 2: benchmark, hot paths, SIMD and other opportunities

2026-09-30. Branch `audit-fixes` at `0e379fc` (uncommitted files in the tree are not part of the build,
which was made from a copy). Follow-up to [the first profiling report](../2026-09-29-performance-profiling/README.md),
after its optimisations (threaded gzip reader, branch-free AVX2 syncmers, anchored alignment, zstd SAM,
zlib-ng, parallel index passes, 2-bit genes) were merged. Everything here was measured; estimates are
marked as such.

**Machine:** Intel Core Ultra 7 258V (4 performance and 4 low-power cores; WSL2 shows 6 vCPUs today), 23 GB,
Ubuntu 24.04 in WSL2, gcc 13, valgrind 3.22, on AC power. No `perf` (no hardware counters). Other Claude sessions
ran experiments on the same machine during the whole session (load average 1–6, one to three busy cores), so
**single wall-clock numbers swing by up to 25%**. Comparisons below are alternated A/B runs, pinned to one vCPU
where they say so, and instruction counts (callgrind), which do not depend on load.

**Data** (as in round 1, in WSL `~/protal-perf`): `db900n`, the 900-species world database rebuilt by this
build (26.86M index values, `database.protal` 125 MB; built in 96 s at 6 threads, 391 s at 8 threads in round 1);
`mix` (1M pairs, 5% from database species, the realistic case), `w900` (1M pairs, 53% align, stresses alignment).
"Previous binary" is the release build left in `~/protal-perf/build-rel` on 2026-09-29 20:25 (the `performance`
branch around the anchored-alignment commit; the exact commit was not recorded), run on the round-1 `db900`.

## Summary

_Follow-up of 2026-10-01 at the end of this file: the first three fixes are implemented, with a denser database._

1. **Big gains since round 1.** Alignment of `mix` at 1 thread takes 9.5 s instead of 18.0 s, at 6 threads
   3.2 s instead of 5.6 s; `w900` 21.9 s instead of 31.2 s at 1 thread. Instructions per pair in the alignment loop:
   `mix` 156k → 84k, `w900` 405k → 302k. Reading gzipped input is no longer the bottleneck (reader 2.3 s → 0.4 s
   at 1 thread; `.gz` costs 5–10% over plain FASTQ at 6 threads).
2. **The biggest remaining costs are not the ones the last round attacked.** On `mix` (per pair, callgrind,
   fixed costs removed by differencing): reverse complement and read copies 16% of the instructions, syncmer
   extraction 22%, WFA 23%, gzip inflate 13%; in wall time (gdb sampling) index lookups ≈ 35–39%, WFA 14%,
   syncmers 13–17%, stage timers 9–13%.
3. **Three cheap fixes remove 19% of the wall time on `mix` and 13% on `w900` (1 thread)**: no per-read clock reads
   (they cost 11–14%), a reverse complement that is not 1.2 µs per read (table instead of a switch, no copies, once
   per read), no per-anchor copy of the read. Measured in a prototype (`build_quick.sh`): identical output
   on 77k SAM records, −15.6k/−18.6k instructions per pair.
4. **Syncmer extraction is mostly scalar set-up, not the AVX2 part**: two per-base loops (~8.2k of the
   ~10k instructions per read). Working on a 2-bit packed read should save about 12k of `mix`'s 99k instructions per pair.
5. **Index loading is now the largest fixed cost** (11 G instructions, 0.6–0.9 s at 6 threads = 20–30% of a
   1M-pair `mix` run) and its decoder is scalar: 134M mostly empty blocks cost ~20 instructions each, and the per-key
   loops run 3 times. Vector fill and prefix sums should take the decode from 7.5 G to ~3 G instructions (estimate).
6. **Things that did not pay on this machine**: software prefetching of the key map (−4% of the seeding timer,
   −1% of the loop), WFA2 built for AVX2 (−3.6% of w900's loop instructions, no change in wall time). Profile-guided
   optimisation gave −3% (`mix`) and −8% (`w900`).
7. **Seeding is memory-bound and only half-understood.** 39% of the loop's wall time on `mix`, 3% of its instructions,
   60 simulated last-level misses per pair; the benchmark of random reads gets 8.5 ns per access with huge pages,
   protal about 40 ns. Prefetching did not fix it. Needs hardware counters (`perf stat`) on a cluster node.
8. **Gaps:** GTDB scale (flex blocks, `Entry::Get`, index load ×70) and thread counts above 6 are not measured here.

| # | Finding | Evidence | Opportunity | Gain (measured / est.) | Changes results? |
|---|---|---|---|---|---|
| 1 | Per-read stage timers: ~50 clock reads per pair (16 ns each) | timers compiled out: 9.3 → 8.3 s min, −11…−16% wall (`mix`, 1 thread) | time only with `--verbose`, or sample every 64th read | −11…−16% wall | no |
| 2 | `ReverseComplement(std::string)`: by value, appended per character through a switch; 4 per pair; the anchor finder and the alignment handler each compute it, and `AlignAnchor(..., std::string rev)` copies it per anchor | 16k of 99k instructions per pair on `mix`; 1200 ns vs 38 ns (table) per 150 bp in a micro-benchmark | table lookup or AVX2 shuffle, `const&`, once per read | −15.6k (`mix`), −18.6k (`w900`) instr/pair; −9…−12% wall | no |
| 3 | Syncmer extraction: the two scalar per-base loops of `ScanWindowsAvx2` | 8.2k of ~10k instr per read in these loops (line-level callgrind); 13–17% of wall in the sampler | 2-bit pack the read once (vector), cut k-mers and s-mers with shifts | est. −12k instr/pair on `mix` (−12%), −3% on `w900` | no (same k-mers, unit-testable) |
| 4 | Index decode `DecodeChunk`: 134M blocks, mostly empty, at ~20 instr each; per-key loops with three `(n+2)/3` | 7.5 G of 19 G instructions in a 50k-pair `mix` run; `Load Index` 0.57–0.90 s at 6 threads | fill empty runs by bitmap word with vector stores; AVX2 prefix sums and flex counts | est. −4 G instr, −0.3 s at 6 threads; more at GTDB size | no (byte-identical index) |
| 5 | Seeding (key map lookups) is memory-bound | 35–39% of wall samples, 60 LL misses/pair (simulated), 3% of instructions | find out why MLP is ~3 (TLB? VM?), then batch lookups across reads | unknown; prefetch prototype −1% | no |
| 6 | PGO (GCC) | −3.2% `mix`, −8.4% `w900` (1 thread, pinned) | PGO in the release/conda build | −3…−8% | no |
| 7 | WFA is 41% of `w900`'s instructions: ends-free bookkeeping and wf-adaptive ≈ 15% of the run | per-function callgrind; 2.04 WFA calls per read; 16% of flank calls have exactly one mismatch, 60% five or more | ungapped check for 1-mismatch flanks (provably optimal); fewer or cheaper calls for divergent flanks | est. −4…−5% of w900's loop for the 1-mismatch case | no (1-mismatch case); to validate otherwise |
| 8 | Flex-block scan in `GetFromLookup`: scalar similarity + a `vector<uint16_t>` | micro-benchmark: AVX2 is 3.7× faster at 64 cells, 3.8× at 1024; nothing at 3–8 | AVX2 (nibble popcount), two passes without the vector | small on `db900` (0.8k instr/pair); **scale-dependent, matters on GTDB** | no |
| 9 | Profiling runs after alignment on one thread | `w900` t6: 2.3 s of 10.4 s; `mix` 0.23 s | overlap with alignment, or parse fewer bytes | up to −20% on alignment-rich samples | no |
| 10 | Smaller: `Uppercase` (1.6k instr/pair), CIGAR run-length pass (3.3% of `w900`), malloc/free (4%), per-record `.Sequence()` decode in output (1%), `FastxReader.cpp` and other `protal_lib` sources are built for baseline x86-64 | callgrind | see below | 1–3% each | no |
| 11 | Start-up: model parse 1.14 G instr (0.41 G still in 24,674 exceptions), genomes 1.25 G | callgrind | the exceptions come from `MiningField`/`OpType` parsing | −0.4 G | no |

## 1. Benchmark: new against previous binary

`scripts/matrix.sh 3`: alignment stage only (`--no_profile`), alternated old/new runs, medians of 2–3 (1 thread: 2),
`--verbose` stage timers, `db900` (old) and `db900n` (new). Full table: `data/matrix_medians.txt`, raw: `data/matrix.tsv`.

| 1M pairs | previous wall | current wall | Δ | previous user CPU | current user CPU |
|---|---:|---:|---:|---:|---:|
| `mix`, `.gz`, 1 thread | 17.98 s | 9.53 s | −47% | 16.5 s | 9.8 s |
| `mix`, `.gz`, 6 threads | 5.55 s | 3.18 s | −43% | 28.9 s | 15.5 s |
| `mix`, plain, 6 threads | 4.71 s | 2.88 s | −39% | 23.6 s | 13.7 s |
| `w900`, `.gz`, 1 thread | 31.21 s | 21.91 s | −30% | 29.6 s | 22.1 s |
| `w900`, `.gz`, 6 threads | 8.62 s | 6.81 s | −21% | 47.3 s | 36.9 s |
| `w900`, plain, 6 threads | 9.23 s | 6.48 s | −30% | 50.1 s | 35.2 s |

Reader time (mean per thread, includes waiting for the lock), `.gz`: 1 thread 2.33 → 0.38 s (`mix`); 6 threads
1.17 → 0.27 s. The alignment loop, `mix` 1 thread: 14.96 → 7.45 s. Peak RSS 3.3–3.5 GB (was 3.4–3.5).

**Thread scaling** (`scripts/threads.sh`, `mix` plain, loop time, medians of 3, one other process busy):
1 thread 7.35 s, 2: 4.21 (1.75×), 4: 2.75 (2.7×), 6: 2.26 (3.25×); user CPU 8.8 → 10.2 → 12.8 → 13.65 s (+55%).
The CPU has 4 fast and 4 slow cores and was shared, so this says little about a cluster node; the extra
CPU time per pair is the thing to check there (memory contention in seeding, or frequency under load).

**Fixed and serial costs** (`mix`, plain, 6 threads, total 2.88 s): preload genomes 0.13 s, index load 0.6 s,
alignment loop 1.98 s, model + profile + rest ≈ 0.2 s. A 1M-pair sample at 16 or more threads will spend
more time in the serial start-up than in the loop.

**A full run** (with profiling, `--no_qcmsa`, 6 threads, loaded machine): `mix` 4.35 s, of which profiling
0.23 s; `w900` 10.44 s, of which loading 0.57 s, alignment 6.4 s, profiling 2.3 s (one thread, 894k SAM records
parsed from the compressed SAM).

## 2. Where the instructions go

callgrind, first 50,000 pairs, 1 thread, profiling build (`-g -fno-omit-frame-pointer`, same optimisation flags).
Per-pair costs are the difference between a 50k-pair and a 1k-pair run (`scripts/cg_diff.sh`), so index load,
model parsing and start-up cancel. Raw tables: `data/perpair_*.tsv`, `data/categories_*.txt`,
`data/callgrind_stages_*.txt`.

| Stage, instructions per pair | `mix` round 1 | `mix` now | `w900` round 1 | `w900` now |
|---|---:|---:|---:|---:|
| alignment loop | 156k | **84k** | 405k | **302k** |
| FASTQ reader incl. gzip inflate | 26k | 5k + 13k inflate (`inflate_fast_avx2`) | 26k | 5k + 12k |
| syncmer extraction | 63k | 22k | 63k | 26k |
| seeds + anchors (incl. reverse complement) | 18k | 17k | 40k | 42k |
| alignment handler | 47k | 37k | 255k | 200k |
| · of which WFA (`alignEndsFree`) | 35k | 23k | 214k | 145k |
| SAM output handler (incl. `ExtractSNPs`, `VariantHandler`, gene decode per record) | 1k | 1k | 16k | 26k |
| profiling, after alignment | 4k | 2k | 35k | 28k (SAM parse) |
| WFA calls per read | 0.16 | 0.20 | 1.46 | 2.04 |
| calls of `GetFromLookup` per read (round 1: "lookups turned into seeds") | 2.9 | 2.9 | 10.3 | 10.3 |

Exclusive cost by category (per pair, both datasets, `data/categories_*.txt`):

| `mix`, 99.3k instr/pair | | `w900`, 353k instr/pair | |
|---|---:|---|---:|
| strings: reverse complement and copies | 25.9k (26%) | WFA (`wavefront_*`) | 144.9k (41%) |
| WFA | 23.0k (23%) | strings | 74.1k (21%) |
| syncmer extraction | 22.1k (22%) | syncmer extraction | 26.2k (7%) |
| gzip inflate | 12.9k (13%) | profiling stage | 17.8k (5%) |
| malloc/free | 4.3k (4%) | seeds + anchors | 15.4k (4%) |
| memcpy/memset/memchr | 2.5k (3%) | malloc/free | 14.8k (4%) |
| seeds + anchors (instructions only) | 2.2k (2%) | alignment post-processing, SAM | 14.5k (4%) |
| stage timers (instructions only) | 1.4k (1%) | gzip inflate | 12.8k (4%) |

The "strings" category is every `std::string` operation callgrind sees; about 16k of `mix`'s 26k per pair are the
four reverse complements (2 in the anchor finder, 2 in the alignment handler, which also copies the read), the rest
FASTQ record assignment, read ids and CIGAR strings. Removing the reverse-complement cost and the copies in the
prototype removed 15.6k (`mix`) and 18.6k (`w900`) instructions per pair, which confirms it.

Fixed costs per run, instructions (1k-pair run): index load 11.0 G (zstd 4.5 G, `DecodeChunk` 6.2 G + 0.94 G
`memcpy`), model parse 1.14 G (0.41 G of it in 24,674 exceptions from `MiningField` and `OpType` parsing), genomes 1.25 G.

### Instructions are not time: where wall-clock goes

Instruction counts hide stalls. A gdb sampler (`scripts/sample.sh`: the profiling build under gdb, SIGINT every
20–50 ms, innermost inlined frames, pinned to one vCPU, `mix` plain, 160–300 samples each, so ±4 points):

| stage (first frame of the stack that belongs to it) | share of samples | share of instructions |
|---|---:|---:|
| index lookups / seeding (`Seedmap::Get` and callers) | 35–39% | ~3% |
| WFA | 14–16% | 23% |
| syncmer extraction | 13–19% | 22% |
| stage timers (`steady_clock::now`, `clock_gettime`) | 9–13% | 1.4% |
| reverse complement | 5–6% | 16% |
| FASTQ reader | 3–5% | 5% |
| gene decode | 2–3% | <1% |

The timers and the reverse complement shares agree with the A/B runs (−11…−16% and −9…−12%). The seeding share is
the odd one: 3% of the instructions, more than a third of the time.

## 3. SIMD opportunities

Ordered by expected gain. "Micro" = `scripts/bench_micro.cpp`, `scripts/bench_mem.cpp`.

**3.1 Reverse complement (and stop doing it four times).** Current code: 1,200 ns for 150 bp; a 256-entry
table 38 ns (0.25 ns per base); my AVX2 shuffle version 138 ns was slower than the table (too many blends;
the byte shuffle needs two lookups to keep non-ACGT bases as they are). A table is enough; the larger gain
is structural: one reverse complement per read, shared by the anchor finder (`ChainAnchorFinder::operator()`),
the alignment handler (`SimpleAlignmentHandler::operator()`) and the SAM writer (`ArtoSAM`), and no copy of the
read (`auto fwd = sequence`, `AlignAnchor(..., std::string rev)`). Measured in situ: −9…−12% wall, −15.6k
instructions per pair on `mix`.

**3.2 Syncmers from a 2-bit packed read.** `ScanWindowsAvx2` is AVX2 only for the last step (min over the
9 s-mers of each of 8 windows, ~130 instructions per 8 windows). Its two set-up loops are scalar and cost
~20 instructions per base: `FillSmers` (385M of the 935M instructions in a 50k-pair run) rolls the forward and
reverse s-mers, and the loop over the read (432M) rolls the forward and reverse k-mers and stores four arrays.
Plan: convert the read once to 2 bits per base (`(c >> 1) & 3` with the G/T swap, 32 bases per
`vpmovmskb`-style step; reads with a non-ACGT base keep the scalar path, since `N` codes as 0 on both strands),
build the reverse strand by complement and pair-reversal of the packed words, and take each window's k-mer,
core and s-mers with an unaligned 8-byte load and a shift (`vpsrlvq` for 4 windows). That is ~3 instructions per
window and strand. Estimate −12k of `mix`'s 99k instructions per pair, about −8…−10% wall; −3% on `w900`.
The unit tests that compare the scan with `WindowByWindow` cover it.

**3.3 Index decode.** `DecodeChunk` (`Hash/IndexCodec.h`) on `db900n`: the loop over blocks runs 134M times
(`4^15` key slots / 8 keys per block), 1.07 G instructions for the `bitmap` test alone and ~0.5 G for the store of
each empty block's 24-byte pattern. Test the bitmap a 64-bit word at a time; for an all-empty word write 64 × 24
bytes of the same pattern (the running value offset `v`) with AVX2 stores. For filled blocks, the counts → offsets
step is a prefix sum over 8 uint16 (one 128-bit vector), and `FlexCells(n) = (n+2)/3` and `n - FlexCells(n)` are
computed in three loops; computing them once per key in a vector (multiply-shift) removes two. `OrPlanes`
is already vectorised by the compiler (234M instructions for 27M values). Estimate: 7.5 G → ~3 G
instructions. At GTDB scale the dense part dominates (≈4 G values), where the per-value byte planes matter
more than the empty blocks; the same vector code applies.

**3.4 Flex-block scan.** `KmerLookupSM::GetFromLookup`/`Get` compute `Seedmap::Similarity` for every flex cell, push it to a
`std::vector<uint16_t>`, then scan the vector again. In the micro-benchmark, with eight 32-bit cells per
register (xor, or-shift, andnot, nibble-popcount with `vpshufb`, `vpmaddubsw`/`vpmaddwd`; two passes: maximum, then
a `movemask` of the cells at the maximum):

| cells | current | AVX2 |
|---:|---:|---:|
| 3 | 4.7 ns | 4.5 ns |
| 8 | 8.1 ns | 5.8 ns |
| 16 | 15.7 ns | 8.1 ns |
| 64 | 65 ns | 20 ns |
| 256 | 232 ns | 68 ns |
| 1024 | 948 ns | 255 ns |

On `db900` this function is 0.8k instructions per pair, too small to matter; on a GTDB index, where conserved marker
genes put hundreds of entries behind one core, it is one of the things that grows. The block-size histogram of
the real index decides; it is a 30-line function, so implement it when that is measured (`max_ubiquity`
limits what is kept after the scan, not what is scanned).

**3.5 Longest exact match.** `ExtendSeed` (left and right) in `ChainAnchorFinder`, `ExtendSeedLeft/Right` and the
link verification in `AnchoredAligner::Align` compare base by base. Micro: 30-base run 8 ns scalar, 1.6 ns
AVX2; 100-base run 27 ns, 1.5 ns. Low value per call (these are 2–5 short comparisons per aligned read) but
trivial. `IdenticalIgnoreAmbig`/`IdenticalIgnoreN` run per seed.

**3.6 CIGAR bookkeeping.** `AlignmentInfo::GetInstructionCountsAndCompress` (per-base loop with six counters and
a `std::to_string` per run) is 11.5k instructions per pair on `w900` (3.3%), and runs for every candidate alignment.
Six compares and `popcnt` per 32 operations, and run boundaries from `movemask` of `ops[i] != ops[i-1]`.
`IsAlignmentValid` and the `ops += 'M'` loops in `AnchoredAligner::Align` are the same kind.

**3.7 Smaller ones.** `Uppercase` (`FastxReader.cpp`: 1.6k instructions per pair; one 32-byte `or 0x20`-style
step, or skip when the first 32 bytes are upper case); the WFA2 kernels (see 4.3: no gain seen).

Not SIMD candidates: WFA's wavefront recurrences (already `-O3`, branchy, latency-bound: ~145k instructions per
pair on `w900`, but see 4.2 for the algorithmic side), zlib-ng inflate (already AVX2 with run-time dispatch).

## 4. Other opportunities

**4.1 Stage timers (11–16% of wall, no code change in the algorithm).** `Benchmark::Start/Stop` read
`steady_clock` (16.2 ns per call here, vDSO) about 50 times per `mix` pair: 16 in `RunPairedEnd`, 28 in
`ChainAnchorFinder::operator()` (14 clock reads per read, for its 7 timers), 2 per read in `RecoverAnchors`, 4 per aligned anchor. A build with both
functions empty ran 9.3 → 8.3 s at best (min of 6 alternated runs; medians 10.7 → 9.0 s). Options: take the
timestamps only with `--verbose`, or time one read in 64.

**4.2 WFA and alignment policy.** On `w900` WFA is 41% of the instructions: `wavefront_extend_matches_packed_endsfree`
27.6k, `wavefront_heuristic_wfadaptive` 16.9k, `compute_affine_idm` 16.1k, then the ends-free bookkeeping
(`compute_init_ends` 10.9k, `termination_endsfree` 10.4k, `compute_process_ends` 9.2k, `fetch_input` 7.2k: 10.7% of
the per-pair instructions together). Counting the flanks that `AnchoredAligner::Flank` hands to WFA (`scripts/build_flankcount.sh`, 100k pairs
of `w900`): 254,843 flank calls, mean length 52 bases; ungapped mismatches after the maximal exact match: 1 mismatch
16%, 2: 8%, 3: 5.5%, 4: 4.4%, 5 or more: 60% (the anchor is extended to the first mismatch, so the count is never
0). A flank with exactly one mismatch has penalty 4 against at least 8 for any alignment with a gap, so it needs no WFA;
that is 16% of the calls, less of the time (they are the cheap ones), est. 4–5% of `w900`'s loop. The 60% of flanks with five
or more mismatches are relatives' alignments; whether they should be aligned at all (they are what `-c`/`--align_top`
and the 90% identity cut-off decide) is a sensitivity question, not a speed one. `mix`: 33,330 flank calls per 100k pairs, 6%
with one mismatch, 81% with five or more, which is the hopeless-anchor case of round 1.

**4.3 WFA2-lib is built without AVX2 in `protal_avx2`.** `lib/wfa2-lib.cmake` builds one static `wfa_lib` with the
baseline flags for both binaries, so WFA2's own `#if __AVX2__` kernels (`wavefront_extend_matches_packed_endsfree_avx2`
and the `AVX512` ones) are in neither: `nm` finds no `wavefront*avx2` symbol in the release binary. Built for
`x86-64-v3` (`build_wfa_avx2.sh`, copy of the tree) the kernels are there: −0.54 G instructions of 31.8 G on `w900`
(extend −25%, `compute_affine_idm` −35%, `init_ends` +14%), **but no measurable wall-clock change** (21.70 →
21.85 s, pinned, 4 alternated runs; `mix` 9.03 → 9.19 s). Not worth a second library unless cluster nodes show otherwise.

**4.4 Profile-guided optimisation.** GCC 13 `-fprofile-generate`/`-fprofile-use` on the same source, trained on
300k pairs of `mix` and 150k of `w900` (`build_pgo.sh`), C and C++ both profiled: `mix` 9.11 → 8.81 s (−3%), `w900`
23.03 → 21.09 s (−8%), 1 thread, pinned, 6 and 4 alternated runs. The training reads are the first reads of the test
reads, so take the figure as an upper bound. The binary shrank from 3.3 to 2.6 MB. A PGO step needs a training
database and reads in the release build; not free, but a single-digit gain for no code change.

**4.5 Index lookups: memory behaviour (unexplained).**
- Sampler: 35–39% of the loop's wall time in `Seedmap::Get` and its callers (`mix`, 1 thread), for ~1.6k instructions per pair.
- Simulated caches (`scripts/cachesim.sh`, callgrind `--cache-sim` with a 12 MB LL, 20k against 2k pairs): 78 last-level data
  read misses per pair, **60 of them in the lookups** (~46 lookups per pair, counting round 1's 22.9 per read, since the k-mers are unchanged; ~1.3 lines per lookup), 8 in `memchr`
  (the FASTQ text) and 7 in `GetFromLookup` (values). No prefetcher or TLB is simulated.
- `bench_mem`, on this machine, a 3.4 GB table with 24-byte blocks: independent random reads cost 8.5 ns each
  with transparent huge pages (21.7 ns with 4 KB pages); batches of 22 with a branch on each: 10.9 ns; prefetching the next read's
  22: 8.7 ns. So 60 misses should cost ~0.6 µs per pair; the sampler says ~2.5 µs (39% of 6.5 µs), i.e. protal
  overlaps 2-3 misses where the benchmark overlaps ~10 (the sampled share also holds the lookups' own instructions, ~0.15 µs).
- The prototype (`build_prefetch.sh`: prefetch both lines of every lookup's control block before the `Get` loop, then the value
  and flex blocks of the 12 smallest lookups before `GetFromLookup`) gave: Seeding timer 1.86 → 1.78 s (−4%), loop −0.6%
  (pinned, 8 alternated runs). The huge-page madvise is active (`AnonHugePages` 3.35 GB of 3.43 GB RSS).
- What is left to try: prefetch across reads (issue read *n+1*'s lookups while read *n* aligns), which `bench_mem`
  shows is the only variant that helped (−20%); and hardware counters (`perf stat -e dTLB-load-misses,LLC-load-misses`) on a
  bare-metal node, because WSL2 adds a second level of page tables and the answer may differ there.

**4.6 Profiling stage.** One thread, after the last read: `w900` 2.3 s of a 10.4 s run, `mix` 0.23 s. 28k instructions per
pair on `w900` for parsing the SAM (`ReadSamGroups`, `strtol`, `SequenceRangeHandler::CoveredPortion`). Handing the alignments over in
memory, or parsing while aligning the next sample, removes it from the critical path (round-1 item 6, still open).

**4.7 Small things seen in the profile.**
- `Gene::Sequence()` now decodes the whole 2-bit gene on each call (`GeneSequence`, up to 4 KB on the stack). I expected this to
  hurt: it is called in `ExtractAndExtendAnchorFromSeed`, `ExtendAnchor`, `AlignAnchor` and once per SAM record in
  `ExtractSNPs`. Measured: `GeneSequence`'s constructor is about 3.7k instructions per pair on `w900` (1%), 0.4k on `mix`. Fine
  for now; a window decode (`UnpackRange`) only matters for long genes or long reads.
- malloc/free are 4% of the instructions (24 allocations per pair from the reverse complements alone).
- `KmerUtils::Complement` returns `N` for lower case; reads are upper-cased first, so the table version matches.
- `src/*.cpp` of `protal_lib` (FASTQ parsing, `BufferedOutput`, taxonomy) is compiled for baseline x86-64 also in `protal_avx2`.
  About 2–3k instructions per pair; ignore.
- cPMML still throws 24,674 exceptions while parsing the model (0.41 G instructions, ~0.1 s).

## 5. What I would do, in order

1. Stage timers off unless `--verbose` (11–16%), reverse complement once per read through a table, no read copies (9–12%). An hour of work; the prototype is in `scripts/build_quick.sh`.
2. Packed-read syncmer extraction (3.2): ~8–10% on realistic samples.
3. `DecodeChunk` vector fill and prefix sums (3.3), together with a measurement of the load on a bigger index.
4. PGO in the release build (3–8%).
5. On the cluster: `perf stat` on the real GTDB index for the lookups (4.5), the block-size histogram for 3.4, and thread scaling at 16 and 32.
6. Flank shortcut for one mismatch (4.2) and the profiling overlap (4.6) when those stages show on real samples.

## Gaps

- **GTDB scale is not measured.** `db900n` has 0.025 values per key slot; GTDB r226 has ~70 times the values. Everything
  that scales with values (flex scans, `Entry::Get`, `GetFromLookup`, index load, WFA calls through ties) is under-represented here.
- Single-vCPU pinned A/B runs are stable (±1–2%); the thread-scaling numbers are not (the machine was shared, 4 fast + 4
  slow cores, 6 vCPUs visible). The sampler has ±4 points at 160–300 samples.
- Simulated cache misses have no prefetcher and no TLB; real hardware counters were not available.
- Real reads (quality, adapters, host DNA), other read lengths, single-end and long reads, `qcmsa` were not run.
- The "previous binary" comparison mixes several changes (zlib-ng, libdeflate removal, gene store, parallel passes and whatever
  else was in the 2026-09-29 20:25 build); it shows where the project is, not the effect of one commit.

## Reproducing

Scripts are in [`scripts/`](scripts/) (`env.sh`: `PERF_DIR`, default `~/protal-perf`, on a Linux file system). Round-1 scripts
(`prep_db.sh`, `prep_reads.sh`) are in the [first report](../2026-09-29-performance-profiling/scripts/).

```bash
export PERF_DIR=~/protal-perf
bash scripts/build.sh                                   # build-rel and build-prof from the current tree
bash scripts/prep_db.sh db900n ~/tune/world             # the database, with this build
bash scripts/matrix.sh 3; bash scripts/summarize_matrix.sh            # old against new (set OLD_BIN, OLD_DB)
bash scripts/threads.sh $PERF_DIR/build-rel/protal_avx2 $PERF_DIR/db900n $PERF_DIR/reads/mix_plain 3 "1 2 4 6"
bash scripts/callgrind.sh r2_mix $PERF_DIR/db900n $PERF_DIR/reads/mix 50000      # and the same with 1000 pairs
bash scripts/cg_diff.sh r2_mix 50000 r2_mix1k 1000 50; bash scripts/cg_categories.sh r2_mix
PROF_BIN=$PERF_DIR/build-prof/protal_avx2 bash scripts/sample.sh mix3 $PERF_DIR/db900n $PERF_DIR/reads/mix3 400 5 0.02
PROF_BIN=... bash scripts/cachesim.sh cs_mix20k $PERF_DIR/db900n $PERF_DIR/reads/mix_plain 20000   # and 2000; cg_cachediff.sh
bash scripts/build_quick.sh; build_prefetch.sh; build_notimers.sh; build_wfa_avx2.sh; build_pgo.sh; build_flankcount.sh   # experiments, on copies of the tree
PIN=2 bash scripts/ab.sh A/protal_avx2 B/protal_avx2 $PERF_DIR/db900n $PERF_DIR/reads/mix_plain 1 6
bash scripts/same_output.sh A B $PERF_DIR/db900n $PERF_DIR/reads/w900_plain 100000
g++ -O3 -march=x86-64-v3 scripts/bench_micro.cpp -o bench_micro && ./bench_micro; g++ -O3 -march=x86-64-v3 scripts/bench_mem.cpp -o bench_mem && ./bench_mem
```

Note for WSL: background jobs started from a `wsl -e bash -c '...&'` die with SIGHUP when that call returns; use
`setsid nohup ... < /dev/null &`.

## Follow-up (2026-10-01): the cheap fixes, implemented, and a denser database

The first three rows of the ranked list (stage timers, reverse complement and copies) are in the tree, with two
more that a denser database made worth doing (the CIGAR counting and the exact-match appends in the anchored
aligner). Everything below is measured on the code as committed (HEAD plus only these changes, built from `git archive`
with the patch applied: `scripts/build_clean.sh`), against HEAD.

**What changed**

| | |
|---|---|
| `Benchmark` (`Utilities/Benchmark.h`) | the stages that run per read or per anchor (`Benchmark::kPerRead`) time every 61st call, the first included, and report the mean timed interval times the number of calls; the period is odd so that mate 1 and mate 2 are timed alike. Coarse stages (load, align, profile) are exact as before. `P_runtime.tsv` and `--verbose` keep their rows; the per-read ones are estimates (`docs/running.md`) |
| dead timers | the profiling stage's `bm_add_sam`, `m_add_sam`, `m_cigar_info`, `Strain`'s three and `VariantHandler`'s were started and stopped for every SAM record and never read; removed |
| `KmerUtils::ReverseComplement` | takes a `string_view`, 256-entry table; `ReverseComplementInto` fills a caller's buffer. The anchor finder computes a read's reverse complement once (into its buffer, `ReverseComplement()`) and the alignment handler takes it (`operator()` with a `reverse` argument; without it, as before); the SAM writer fills the entry's buffer; the handler no longer copies the read, and `AlignAnchor` no longer copies the reverse complement per anchor |
| `AlignmentInfo::GetInstructionCountsAndCompress` | counts run by run and writes the digits itself, no `std::to_string` per run (for every candidate alignment) |
| `AnchoredAligner::Align`/`Between` | a link that is exact is appended as one run of `M` (`memcmp`), the window bookkeeping counts with `std::count` |

Not changed, and why: the per-anchor whole-gene decode (`Gene::Sequence()`), which a dense database makes visible (below), needs
a windowed decode and positions relative to the window, not a cheap edit.

**Same outputs.** `scripts/compare_outputs.sh` runs HEAD and the new build on the same reads and compares every output
file but `*_runtime.tsv` (timings) byte by byte: 176 files for `mix`, 369 for `w900`, 241 for the dense world's `dense_w`, 140 for
`dense_mix`, with one thread (the zstd SAMs too); at 6 threads, and with `.gz` input, the SAMs equal when sorted (the order of
records differs from run to run) and every other file equal; single-end reads (SAM); the baseline-ISA binary as well as `protal_avx2`.
No differences.

**Tests**, on the clean tree: 215 unit tests (Release, and under ASan + UBSan), 24 mini-database generator tests, 113 end-to-end tests,
the 3 `GtdbBuildTest`s (synthetic release downloaded, built and trained), 7 model tests: all pass. New: the reverse complement
for every byte and random reads, the buffer reuse, the handler given a reverse complement against making its own, for both strands,
the CIGAR counts against the column loop (random strings, 16-bit wrap), and the benchmark's sampling (which calls are timed, both mates,
the estimate from below, join, print).

**The denser database.** `scripts/dense_world.sh`: 4 genera of 60 species (3 genomes each, species 1.5–6% apart, strains 0.2–2%), 240 species, 720
genomes, 7.7M index values; where `db900` has mostly one to three species per genus, here a gene has up to 60 relatives in its genus.
`dense_w` (1M pairs, 20 species) and `dense_mix` (5% from the database) are made as `w900` and `mix` were. Per read pair (6 threads, 100k pairs,
`--verbose`): seeds per read 61 (`dense_w`), 25 (`w900`), 7 (`mix`); anchors per pair 11.6, 3.1, 0.34; candidate alignments per pair 3.7, 1.1,
0.06; WFA calls per read 3.3, 2.0, 0.2. Real GTDB is denser still, but the direction is the one the fixes had to survive. `GetFromLookup` runs
2.9 times per read on `mix` and ~10 on `w900` and `dense_w` (the seed search stops once it has enough seeds), so what grows with the density is
everything after the lookups.

Instructions per pair (callgrind, 30k pairs against 1k, whole pipeline incl. profiling), HEAD → now:

| | HEAD | now | |
|---|---:|---:|---:|
| `mix` | 98.5k | 77.6k | −21% |
| `dense_mix` | 90.7k | 68.4k | −25% |
| `w900` | 347.8k | 303.1k | −13% |
| `dense_w` | 480.9k | 405.8k | −16% |

Alignment stage wall time, 1 thread pinned to a vCPU, alternated HEAD and new, minimum and median of 3–6 runs per experiment (the machine was loaded by other
sessions, so absolute times are 1.5× those of the first report; the ratios are what to read):

| | HEAD min / median | now min / median | Δ min / median |
|---|---:|---:|---:|
| `mix` (1M pairs), vCPU 2 / vCPU 4 | 16.6 / 17.9 s, 15.7 / 18.5 s | 12.6 / 13.6 s, 13.6 / 14.4 s | −24% / −24%, −13% / −22% |
| `w900` (1M), vCPU 2 / vCPU 4 | 37.9 / 44.3 s, 38.6 / 38.7 s | 30.7 / 32.1 s, 31.8 / 33.4 s | −19% / −27%, −18% / −14% |
| `dense_mix` (300k), two experiments | 5.3 / 5.4 s, 5.9 / 6.3 s | 4.5 / 4.8 s, 4.7 / 5.3 s | −16% / −11%, −20% / −15% |
| `dense_w` (300k), two experiments | 22.4 / 22.6 s, 24.4 / 25.3 s | 19.8 / 20.3 s, 21.4 / 22.2 s | −11% / −10%, −12% / −12% |

(On a quiet machine the first three fixes alone gave −19% on `mix` and −13% on `w900`, `data/ab_experiments.txt`.) The share of the saving falls
as the database gets denser, because the work that does not change grows: on `dense_w`, WFA is 173k of 481k instructions (36%; 43% now) and
3.3 calls per read.

What a dense database shows that the sparse ones hid, per pair on `dense_w` after the fixes (405.8k): WFA 172.8k (43%); strings and copies 69.5k
(17%; not broken down further); `ChainAnchorFinder::operator()` itself 20.4k and `GetFromLookup` 12.3k
(the flex scans; a block holds up to 180 entries here); `GeneSequence` 12.4k (3%: the whole gene is decoded again for the anchors of a read, 11.6 per pair; 0.4k on `mix`); `AlignAnchor` 13.5k; the syncmers 21.6k, no longer the largest. So the next steps, in the order they pay on dense
databases: fewer or cheaper WFA calls (item 4.2 of the first report: 3.3 calls per read; on `w900` 16% of the flank calls have exactly one mismatch, not measured here), a windowed
gene decode, moving candidate alignments instead of copying them, then the packed-read syncmers and the vector flex scan.

Scripts added: `build_pair.sh`, `build_clean.sh`, `compare_outputs.sh`, `run_all_tests.sh`, `dense_world.sh`, `threads.sh`, `cachesim.sh`.

## Follow-up (2026-10-01): windowed gene decode

Implemented in `f71e7d5` (on `062bdf9`/`b1027b3`, which add `GeneSequence` windows and `Gene::Window`).

`Gene::Window(begin, end)` decodes only the bases a read touches. The view it returns is as long as the gene and indexed by gene
position, so the code that takes a gene as a `string_view` works unchanged; the bases outside `[Begin(), End())` are not decoded.
Three places use it, each with the range it can read:

| where | window |
|---|---|
| `ChainAnchorFinder` (`GeneAround`, `ExtractAndExtendAnchorFromSeed`) | the stretch of the gene along each link's diagonal: `[genepos - readpos, genepos - readpos + read length)` per link, merged, with a margin of 8 on both sides. A single seed of unknown strand takes the union of its two diagonals. One decode per anchor now serves the checks of a short anchor and its extension (two before) |
| `SimpleAlignmentHandler::AlignAnchor` | the alignment window itself (`window.ref_start` to `ref_end`), which `AnchoredAligner` checks that every link and flank lies in |
| the SAM writers (`ReferenceOf`) | the reference the record's CIGAR spans, from its position, for `ExtractSNPs` and `PrintAlignment` |

Whole genes are still decoded where whole genes are needed (the strain MSAs, the build, long reads, the profiler's debugging check).

**An error that the review of the first version missed, and how it was found.** `ExtendSeed` extends a link to the right while
`query[qpos] == gene[rpos]`, and for the last link its limit is the gene's length, not the read's: it stops at the read's
terminating NUL only because no base is NUL, so it reads the gene base one past the read's diagonal. With a window that ended
there, that byte was whatever the stack held; when it was 0, the extension ran on and the link came out longer than the read.
One read pair of 20,000 gave a different CIGAR (`62M3D88M` for `62M3D85M3S`, which is the better alignment, by luck of the chain
lengths) in one build and not in another with the same source: a build that depended on stack residue. The
byte-identical comparison caught it (1 record of 1,027), but only because the stack held a 0 in that build; AddressSanitizer did not (it poisons
whole 8-byte granules and leaves the bytes of a partly poisoned one addressable, so a byte just past a window is missed more often than not);
`valgrind --tool=memcheck` did, deterministically (an uninitialised conditional jump, then an invalid read of the read's NUL + 1). The window has a margin of 8 on each side
now, and the tests make the class of error deterministic: `packed::OutsideFill(byte)` fills the undecoded part of every window with a
byte, and the new tests (`tests/test_GeneWindows.cpp`) compare the finder's extension on windows with the one on the whole
gene and the handler's alignments with different fills (0, `#`, `N`), over reads that lie inside the gene, run over either end,
have a gap between their links, on genes shorter than the read, of 1,945 bases (stack) and 6,000 (heap). With the margin set to 0 they
fail. A second trap: the `pair/new` build of that session was stale: the margin was added while a build was running and ninja judged
the objects up to date, so the first comparisons were of the old window. `build_pair.sh` now starts the build directory empty.

**Same outputs.** `compare_outputs.sh`, HEAD against the new build, as before: 176 (`mix`), 369 (`w900`), 241 (`dense_w`) and 140
(`dense_mix`) files byte-identical at one thread, SAMs equal sorted at 6 threads and with `.gz` input, single-end reads, and the
baseline-ISA binary. Under AddressSanitizer + UBSan (`protal_avx2` Debug build, the undecoded part of every window poisoned):
`mix` 40k pairs, `w900` 20k, `dense_w` 12k, `dense_mix` 40k: no error, SAMs and profiles identical to HEAD's. Under
`valgrind --tool=memcheck` (the check that found the error above), on the final build: 1,500 pairs of `dense_w` and 1,500 of `w900`: 0 errors.

**What it saves.** Instructions per pair (callgrind, 30k pairs against 1k), before this change → now: `dense_w` 405.8k → 398.3k (−1.8%),
`w900` 303.1k → 300.5k (−0.9%), `dense_mix` 68.4k → 67.9k, `mix` 77.6k → 77.3k. The constructor of the window is 200
instructions per decode against 603 for the whole gene, at 20.6 decodes per pair on `dense_w`: 12.4k → 4.1k instructions per pair.
Wall time (1 thread pinned, quiet machine, 3–6 alternated runs; `ab.sh`): `dense_w` 10.17 → 10.23 s (min), 11.14 → 10.71 s (median); `w900`
18.17 → 17.77 s (min), 19.48 → 18.88 s (median); `dense_mix` and `mix` unchanged. So 0–4%, within what the machine's noise
allows. The genes of these worlds (24 and 84 MB packed) are reused read after read, and stay in the caches. A cold gene decode is
where the window pays: `scripts/bench_window.cpp`, genes drawn at random from 256 MB of packed genes:

| gene | whole | window of 200 bases |
|---:|---:|---:|
| 1,000 bases | 170 ns | 82 ns |
| 2,000 | 234 ns | 93 ns |
| 4,000 | 278 ns | 85 ns |
| 6,000 (heap buffer) | 351 ns | 186 ns |
| 20,000 | 532 ns | 176 ns |

At ~20 decodes per pair that is up to 2 µs of a pair's ~30 µs on `dense_w` when the genes are cold, as they will be in a GTDB index of
5 GB of packed genes; not measured at that size. The 4 KB inline buffer still makes genes above 4,096 bases allocate their whole length on the
heap for a window (the 186 ns above); a thread-local scratch buffer would remove that.

Scripts added: `bench_window.cpp`, and `build_pair.sh` rebuilds from empty.

Note: `~/protal-perf`, the work folder of both rounds (databases, read sets, builds, callgrind outputs), was removed from the machine
after the checks above; what the report shows is in its `data/` folders, and the scripts regenerate the rest
(`prep_db.sh`, `prep_reads.sh`, `dense_world.sh`, `build_pair.sh`). The commit was checked again afterwards on fresh builds of
`git archive HEAD` (`f71e7d5`) and of `HEAD~1`: 230 unit tests and 116 end-to-end tests pass, and on a mini database (`build_mini_db.sh`)
with 200,000 simulated read pairs all 12 outputs are byte-identical to `HEAD~1`'s at one thread, the SAMs equal sorted at four.
