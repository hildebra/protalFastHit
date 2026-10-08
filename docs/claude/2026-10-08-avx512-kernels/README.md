# AVX-512 kernels for the seed lookup and the syncmer scan, chosen at run time

- **Date**: 2026-10-08.
- **Code**: branch `audit-fixes` at `7dba06d` plus this change (uncommitted when written): `src/Utilities/SimdLevel.h`
  (new), `src/Hash/FlexScan.h`, `src/Hash/KmerLookup.h`, `src/SequenceUtils/KmerIterator.h`,
  `src/SequenceUtils/PackedSequence.h`, the tests `test_FlexScan.cpp`, `test_PackedIndex.cpp`, `test_Syncmers.cpp`,
  `docs/installation.md`, `docs/development.md`, `CMakeLists.txt` (a comment).
- **Machine**: WSL2 Ubuntu 24.04 on the Core Ultra 7 258V, GCC 13.3, 4 cores (`taskset -c 0-3`, niced). This CPU has
  **no AVX-512**, so the AVX-512 code was compiled and inspected here, not run. No emulator is installed.
- **Data**: the v0.7.5 benchmark world (`~/bench071/V075/protal_db`), samples `rl150_p500000_s_1` (pe 500k pairs),
  `pb_b90000000_s_1`, `ont_b90000000_s_1`; callgrind on 100k pairs (`~/perf6/reads/pe100k`) and `pb_b3000000_s_1`.
- **Scripts**: [`scripts/build.sh`](scripts/build.sh) (HEAD and HEAD plus [`scripts/files.txt`](scripts/files.txt) in
  `~/avx512`), [`scripts/check_local.sh`](scripts/check_local.sh) (its output: [`check_local.log`](check_local.log)),
  [`scripts/check_avx512.sh`](scripts/check_avx512.sh) (for an AVX-512 node).
- **Question**: implement AVX-512 for the seed lookup and the syncmer scan, chosen at run time (following the
  [assessment](../2026-10-08-avx512-assessment.md)); which other advanced CPU instructions would help?

## Summary

1. **The level is chosen at run time, per CPU, with one switch.** `simd::Default()` takes the CPU's highest level
   (scalar, AVX2, or AVX-512 of Ice Lake / Zen 4: F, BW, VL, DQ, VBMI, VBMI2, VPOPCNTDQ), capped by `PROTAL_SIMD=auto|
   avx512|avx2|scalar` in the environment. It governs the flex scan, the syncmer scan and 2-bit packing.
   `PROTAL_FLEX_SCAN=scalar` still works. Skylake-X and Cascade Lake (AVX-512 without VBMI2 and VPOPCNTDQ) use AVX2.
2. **The flex scan with AVX-512** (`flex_scan::ScanAvx512`): scores, best score, tie masks and their count in one call,
   16 cells a step, through a funnel shift, two ternary-logic steps, a masked `vpopcntd`, a masked byte store, then the
   ties by mask compares 64 at a time. The loop is ~9 instructions per 16 cells where AVX2 takes ~26 per 8. Masked loads
   read nothing past a block's cells.
3. **The syncmer scan, reworked for both vector levels**: a window's cores are read from the 16-base codes and its
   k-mer is built only if the window is a syncmer, so the per-window k-mer and core arrays are gone. The AVX2 path
   takes **16% fewer instructions** (1.215 → 1.018 G on 100k pairs), the same k-mers. The AVX-512 path codes 64 bases a
   step (masked, no scalar tail), its doubling and s-mer loops are vectorised by the compiler for AVX-512, and it
   evaluates 16 windows a step with a masked last step.
4. **Locally everything that can run is identical.** 452 unit tests pass (3 skipped as before); paired-end, PacBio and
   ONT runs give the same SAM records and the same other files with the new build (AVX2) and with `PROTAL_SIMD=scalar`
   as with HEAD.
5. **The AVX-512 code is unverified on hardware.** Its tests skip here. Before relying on it, run
   `scripts/check_avx512.sh` on an AVX-512 node (the EPYC 9634 nodes are Zen 4), or the unit tests under Intel SDE.

## What changed

| part | before | now |
|---|---|---|
| choice | per kernel: `flex_scan::Avx2Enabled`, `SimpleKmerHandler::UsesAvx2`, `packed::Avx2Enabled` | `simd::Level` (scalar, avx2, avx512), `simd::Default()` capped by `PROTAL_SIMD`; `flex_scan::Kernel()`/`Use()`, `SimpleKmerHandler::Level()`/`Use()`; packing keeps its AVX2 switch, off under `PROTAL_SIMD=scalar` |
| flex scan | `BestAvx2` + `TiesAvx2` (AVX2), `ScoreScalar` | plus `ScanAvx512`, one call per lookup (not inlined into the `PROTAL_CLONE_V3` caller) |
| syncmer fill | codes, `Bases16` per strand into one array, k-mers and cores of every window into four arrays, s-mers | codes, `Bases16` into one array per strand, s-mers; cores from those arrays in the window loop, k-mers only for syncmers |
| syncmer windows | AVX2 8 a step, scalar tail | AVX2 the same (scalar tail from the arrays); AVX-512 16 a step, last step masked |

The shared fill steps (`Bases16`, `Bases16AndSmers`) are `always_inline` functions without a target of their own, so
each is compiled into `FillAvx2` for AVX2 and into `FillAvx512` for AVX-512 (`prefer-vector-width=512` under GCC).

## Checks

Unit tests (`check_local.log`): the whole suite, 452 passed and 3 skipped (the two opt-in benches and the
AddressSanitizer-only test, as before). The kernel suites again with `PROTAL_SIMD=scalar`: 24 passed. New or changed:
`FlexScan.Avx512ScoresAndTiesAsTheScalarScan` (sizes 1-140 and 10 larger, every bit shift, the block's cells copied to the
end of their buffer, no score or mask word written past the block; skipped here), `FlexScan.LevelsFollowTheCpu`,
`PackedIndex.LookupsGiveTheSameSeedsAtEveryVectorLevel` (each level the CPU has against scalar),
`Syncmers.ScanGivesTheDefinitionsKmersForOtherShapes` (eight k, core and s-mer shapes, from a 15-base k-mer with its core
at base 0 to 31-mers with 11-base cores; every level), the other `Syncmers` tests over every level the CPU has.

Whole runs, 4 threads, HEAD against this change (AVX2 here) and against this change with `PROTAL_SIMD=scalar`:

| sample | SAM records | identical, this change | identical, scalar | other files |
|---|---|---|---|---|
| pe 500k pairs | 418,464 | yes | yes | all identical |
| PacBio 90 Mb | 24,947 | yes | yes | all identical |
| ONT 90 Mb | 34,871 | yes | yes | all identical |

Callgrind, one thread (instructions):

| | HEAD | this change | |
|---|---|---|---|
| pe 100k pairs, whole run | 29.02 G | 28.83 G | −0.6% |
| – syncmer scan (`ScanWindowsAvx2`, inclusive) | 1.215 G | 1.018 G | −16% (of which `FillAvx2` 0.505 G) |
| – `GetFromLookup` | 0.515 G | 0.519 G | +0.7% (the level read; the AVX2 path otherwise unchanged) |
| PacBio 3 Mb, whole run | 15.40 G | 15.34 G | −0.3% |

The binary (`objdump`): `ScanAvx512` is 107 instructions with one `vpshrdvd`, two `vpternlogd` (0xC3, 0x80), a masked
`vpopcntd`, a masked `vpmovdb` store and a masked `vpcmpeqb`; `ScanWindowsAvx512` uses zmm throughout; `FillAvx512`
mixes zmm and ymm (the compiler's choice for the widening loops), all AVX-512VL encodings.

## What is expected at r226, and how to check it

The [assessment](../2026-10-08-avx512-assessment.md) estimated the seeding −1 to −2 s per thread of the paired-end
run's 12.9 s (twelfth to thirteenth cluster runs: the last cut of this kernel's instructions, −43%, gave −20% of the
seeding). The syncmer scan (1.95 s per thread) loses 16% of its instructions on AVX2 and more on AVX-512; its fill is
now half of it.

`scripts/check_avx512.sh OUT DB TYPE:R1[:R2] ...` on a node with AVX-512, from a Release build with tests:
1. prints the CPU's AVX-512 flags and the level the kernels use there;
2. runs the kernel tests, which compare each level the CPU has with the scalar code;
3. runs each sample with `PROTAL_SIMD=avx2` and with the default: the SAM records and every other file must be the same;
4. alternates `scripts/measure_performance.sh` at both levels (`ROUNDS`, `REPEATS`): compare "seeding" and "taking
   the k-mers" in the two `stages.tsv`.

Without a node: Intel's Software Development Emulator (a download from Intel, no root needed) runs the tests on this
laptop: `sde64 -icx -- build/tests/protal_tests --gtest_filter='FlexScan.*:Syncmers.*:PackedIndex.*'`.

If the node's checks fail, `PROTAL_SIMD=avx2` keeps runs on the verified AVX2 code until it is fixed.

## Other advanced CPU instructions

The cluster runs say where time goes now (thirteenth run, pe per thread: seeding 12.9 s, alignment handler 2.8 of which
the k-mer screen 1.4, extending anchors 2.1, k-mers 2.0, reader 1.2, seed sort 1.1). Ranked by what they could give:

| instructions | where | expected | why not more |
|---|---|---|---|
| AVX-512 or AVX2 sorting networks (e.g. Intel's x86-simd-sort, 64-bit key and index) | the shared-seed sort, after the 64-bit keys of 2026-10-06 item 6 | up to ~0.5 s per thread | only 18% of the seeds are sorted now; needs 64-bit keys first |
| `vpcompressd` (AVX-512) | turning tie masks into cell indices before the entries are decoded | small | 1-2 seeds per lookup; the bit walk is cheap |
| `vpmultishiftqb` (VBMI) | decoding the 42-bit packed entries of several seeds at once | small | same: few seeds per lookup, per-seed work is the push into the list |
| gathers (`vpgatherdd`) | the k-mer screen's read k-mers against the window's stamp table | uncertain, likely none on Zen 4 | gathers are slow on Zen 4 and the stamps are writes (scatters, slower) |
| GFNI (`vgf2p8affineqb`) | reverse complements of packed strands (bit pairs reversed in a byte) | tiny | reverse complements are ~1.5% of instructions |
| BMI2 `pext`/`pdep` | bit-field extraction | none | the hot code already uses shifts; slow (microcoded) on Zen 1-2 |
| AVX-512 VNNI, AMX, FP16 | none of protal's hot loops is a dot product or matrix product | none | the GBM is scored once per taxon |

The seeding is the largest stage, and the larger part of it is waiting on memory: 2.2G lookups of ~70 cells in a 27 GB
index. No instruction set removes those misses. More lookups in flight per thread (interleaving reads, as the
[GTDB-scale report](../2026-10-04-performance-gtdb-scale/README.md) proposed) and a smaller index entry are the
levers there.
