# Where AVX-512 could pay in protal

- **Date**: 2026-10-08.
- **Code**: branch `audit-fixes` at `7dba06d`, read only. Nothing was built or run.
- **Data**: the stage times of the thirteenth cluster run (`230bf64`, r226 v11, EPYC 9634, 32 threads; the
  [2026-10-06 profiling report](2026-10-06-performance-profiling/README.md)), its lookup benchmark, and the earlier
  SIMD reports ([CPU dispatch](2026-10-02-cpu-dispatch/README.md), [WFA2 AVX2](2026-10-02-wfa-avx2/README.md),
  [one binary](2026-10-02-one-binary/README.md), [GTDB scale](2026-10-04-performance-gtdb-scale/README.md)).
- **Question**: are there parts of protal that would benefit from an AVX-512 implementation?

## Summary

1. **One place clearly would: the flex-cell scoring of every seed lookup** (`flex_scan::BestAvx2` and `TiesAvx2` in
   `src/Hash/FlexScan.h`, called from `KmerLookupSM::GetFromLookup`). Seeding is 12.9 of 24.8 s per thread of the
   r226 paired-end aligning, and the last cut of this kernel's instructions (−43% per lookup, 2026-10-06 item 1) turned
   into −20% of the seeding on the cluster, so its arithmetic still costs time despite the memory waits.
   AVX-512 helps here through instructions AVX2 lacks, not mainly through width (below). Estimated: the scoring's
   instructions to about a fifth, a lookup's to about 70%, the seeding −1 to −2 s per thread (pe), about −5% of the
   paired-end aligning. HiFi's seeding (2.8 s per thread) likewise, smaller.
2. **Two places might, a little**: the closed-syncmer scan's window minima (`KmerIterator.h`, 16 windows per step instead
   of 8; the stage is 1.95 s per thread, and fewer passes, item 8 of 2026-10-06, comes first) and packing a read
   2 bits a base (`packed::Pack`, `AVX512_VBMI` byte permutes; a small share). Each is worth well under 0.5 s per thread.
3. **The rest would not**: WFA2 (its AVX2 kernels were slower than the portable ones, and 90-96% of candidates at r226
   never reach it), the k-mer screen's stamping (table writes and reads at k-mer codes: gathers and scatters, slow on
   Zen 4), the seed sort (now 1.06 s per thread on 18% of the seeds), the emission of seeds (one or two per lookup, a
   41-bit entry decode each), and the profiling stage (3.7 s, scattered map and EM code). zlib-ng, ISA-L and zstd already
   choose their own AVX-512 code at run time.
4. **One thing AVX-512 already changes, untested**: ISA-L's level-1 deflate has an AVX-512 path, which ISA-L picks on the
   EPYC nodes. Its output bytes were checked only on SSE4.2, AVX and AVX2 (`src/IO/Bgzf.h`), and the 2026-10-08 thread-bytes
   fix analysed only those assembly versions. `.sam.gz` and the simulator's `.fq.gz` written on the cluster may therefore
   differ in bytes from those written on the laptop. The content stays the same. Checking it needs an AVX-512 machine, or
   Intel SDE.

## The flex scan with AVX-512

Today, per step of 8 cells (`BestAvx2`, ~26 instructions): two loads, a shift pair and an OR to extract the cells at the
block's bit offset, two XORs, a shift and two ANDs for the equal base pairs, a nibble-table popcount (two shuffles, two
ANDs, a shift, an add, `maddubs`, `madd`), a compare and AND for the lanes past the block, the max, a shuffle and two
extracts and stores for the scores. Then `TiesAvx2` compares the scores 32 at a time. At r226's 70 cells a block that is
~9 steps, about 240 of a lookup's ~660 instructions.

With AVX-512 (`F`, `BW`, `VL`, `VPOPCNTDQ`, `VBMI2`; all on Zen 4 and on Intel since Ice Lake), per step of 16 cells,
~9 instructions:

| AVX2 today | AVX-512 |
|---|---|
| shift right, shift left, OR | `vpshrdvd` (funnel shift, VBMI2), one instruction |
| XOR, XOR with all-ones | `vpternlogd`, one |
| shift, AND, AND | shift, `vpternlogd` |
| nibble popcount: 8 instructions | `vpopcntd`, one |
| lane compare and AND for the last step | a `k` mask from `bzhi`, zero-masking for free |
| shuffle, two extracts, two stores of scores | `vpmovdb` to memory, masked, one |
| `TiesAvx2`: compare, movemask, tail mask, 32 a step | `vpcmpeqb` into a `k` register, 64 a step, tail by mask |

The scores are exact integers, so the outputs stay identical by construction; the test pattern of
`FlexScan.Avx2ScoresAndTiesAsTheScalarScan` (every size and bit shift against `ScoreScalar`) carries over.

**Width matters less than the instructions.** Zen 4 runs a 512-bit operation as two 256-bit halves, so going from 8 to
16 cells a step saves front-end slots but not execution throughput. Most of the gain is `vpopcntd`, `vpternlogd`, the
funnel shift and the masks, and those exist at 256 bits too (AVX-512VL). A 256-bit AVX-512VL version would get most of
it on Zen 4, avoid the clock drop older Intel server parts take on 512-bit code, and also suit AVX10/256 CPUs. Both
widths are worth benchmarking; I would start with 256 bits.

## How it would fit the code

- **Dispatch** as the AVX2 path is chosen today (`flex_scan::Avx2Enabled`), with a third level checked by
  `__builtin_cpu_supports` for `avx512vl`, `avx512bw`, `avx512vpopcntdq` and `avx512vbmi2`, and a
  `PROTAL_FLEX_SCAN=avx2` value to compare runs.
- `GetFromLookup` is a `PROTAL_CLONE_V3` clone. GCC's `arch=x86-64-v4` clone would not do: v4 has no `VPOPCNTDQ` or
  `VBMI2`. The AVX-512 functions get their own `__attribute__((target(...)))`, as `BestAvx2` has today; the masks they
  fill are the ones the shared emission loop already reads.
- **Padding**: the 512-bit version reads up to 64 bytes past a block's last cell. `Seedmap::kPackedPadding` is 64 today,
  so the `static_assert` holds; a masked load would remove the read-past entirely.
- **Testing**: the laptop (Core Ultra 7 258V) has no AVX-512. Correctness can be checked locally under Intel SDE
  (`sde64 -spr --`), speed only on the EPYC nodes, as `scripts/measure_performance.sh` runs there.

## What I would do

Item 1 (at 512 bits) and the syncmer scan of item 3 were implemented the same day: see
[AVX-512 kernels](2026-10-08-avx512-kernels/README.md).

| # | change | exact | expected at r226 |
|---|---|---|---|
| 1 | `BestAvx512` and `TiesAvx512` (256-bit VL first, 512 second), dispatched at run time, `lookup_bench.cpp` and a cluster run to compare | yes | seeding −1 to −2 s per thread (pe), the lookup's instructions to ~70% |
| 2 | check ISA-L's AVX-512 level-1 bytes against its AVX2 bytes (SDE, or a node), and the address-hash fix on that path | – | none in speed; settles whether written gzip is byte-stable on the cluster |
| 3 | the syncmer scan in fewer passes first (2026-10-06 item 8), then 16 windows a step if the scan is still 6% | yes | under −0.5 s per thread |

Not worth it on this evidence: WFA2, the screen's stamping, the seed sort, seed emission, the profiling stage.
