# AVX2 where the CPU has it, and not where it has not

- **Date**: 2026-10-02.
- **Code**: branch `audit-fixes` at `c40cec3`; `src/`, `lib/`, `CMakeLists.txt` and `protal_launcher` are the
  same as at `8798733`. The binary examined is the Release build of `5c7d64c` + `8798733`'s change (the same
  `src/`), built in WSL with gcc 13 (`~/mt-work/samroomhead/src/build/protal`).
- **Machine**: WSL2 Ubuntu 24.04 on an Intel Core Ultra 7 258V, which has every x86-64-v3 extension. No
  emulator (QEMU user mode, Intel SDE) is installed, so the binary was not run as a CPU without AVX2; see
  "Not checked".
- **Question**: does protal use AVX2 and the other x86-64-v3 instructions where the CPU has them, and leave
  them alone where it has not?

## Summary

Yes, by two mechanisms:

1. **The installed `protal` is a launcher** (`protal_launcher`). It runs `protal_avx2`, the whole program
   compiled for x86-64-v3, only if `/proc/cpuinfo` lists all of avx, avx2, bmi1, bmi2, fma, f16c, movbe and
   abm (LZCNT). Otherwise it runs `protal_baseline`, built for plain x86-64. The kernel leaves avx and avx2 out of
   those flags where the OS does not save their registers, and hypervisors that mask them do the same.
   Without `/proc/cpuinfo` (macOS) or on ARM it takes the baseline, and `PROTAL_NO_AVX2` forces it.
2. **The baseline build dispatches at run time itself.** Its only AVX2 code is protal's syncmer scan and
   sequence packing. These are compiled for AVX2 by function attribute and chosen with
   `__builtin_cpu_supports("avx2")` at first use, plus zlib-ng's per-CPU variants, chosen by its own detection.
   Everything else in it is plain x86-64.

So a CPU without AVX2 gets the baseline build through the launcher, and that build runs no AVX2 instruction
on such a CPU. Running `protal_avx2` directly on one would stop with an illegal instruction, as
`docs/installation.md` says ("built for x86-64-v3").

## How the builds are made

| target | compiled for | installed as |
|---|---|---|
| `protal` | `main.cpp` with `-march=x86-64 -mtune=generic` (`isa_baseline`) | `protal_baseline` |
| `protal_avx2` | `main.cpp` with `-march=x86-64-v3 -mtune=generic` (`isa_avx2`) | `protal_avx2` |
| `protal_static` | as `protal`, linked statically | `protal_<version>_static` |
| `protal_lib` (both) | `isa_baseline` | linked into both |
| `wfa_lib` (WFA2) | no `-march`: gcc's default, plain x86-64 (upstream's `-march=native` is not used) | linked into both |
| zlib-ng | `WITH_NATIVE_INSTRUCTIONS OFF`, `WITH_RUNTIME_CPU_DETECTION ON` | linked into both |
| libzstd, libdeflate | the system's (or conda-forge's) shared libraries, which choose their code at run time | |

`just install` and the conda recipe both install the launcher as `protal`, next to `protal_baseline` and
`protal_avx2`.

## The dispatch in the baseline build

- `SimpleKmerHandler` (`KmerIterator.h`): `ScanWindowsAvx2`, `FillAvx2` and `Bases16` carry
  `__attribute__((target("avx2")))`; the handler sets `m_avx2 = CpuHasAvx2()` in its constructor, at run time.
- `packed::Pack`/`Unpack` (`PackedSequence.h`): `PackAvx2`, `UnpackAvx2` and `Unpack128Avx2` likewise; the
  choice is a function-local static (`Avx2Enabled`), set at first use.
- Neither check runs during static initialisation, where `__builtin_cpu_supports` could come before libgcc's
  CPU detection. libgcc reports AVX2 only if the OS saves the YMM registers (XGETBV), not on the CPUID bit alone.
- Both have scalar paths and tests that compare them with the AVX2 ones (`tests/test_Syncmers.cpp`,
  `tests/test_PackedSequence.cpp`, which skip the AVX2 half on a CPU without it).

## The binary

`scripts/isa_check.sh` disassembles the baseline `protal` and lists, per function, the instructions beyond
plain x86-64: ymm/zmm registers, v-prefixed (VEX/EVEX) forms, BMI1/BMI2, LZCNT/TZCNT, MOVBE
(`results/isa.txt`). They occur only in:

- protal's `ScanWindowsAvx2`, `Bases16`, `PackAvx2` and `UnpackAvx2`, the dispatched functions above;
- zlib-ng's per-CPU implementations (`*_avx2`, `*_avx512`, `*_avx512_vnni`, `*_vpclmulqdq` and their static
  helpers `chunkcopy_safe`, `partial_fold` and `adler32_fold_copy_impl`), chosen by its run-time detection;
- TZCNT in plain functions (`std::from_chars`, WFA2's portable extend and CIGAR kernels, zlib-ng's `*_sse2`).
  That is gcc's encoding of `__builtin_ctz` for any x86-64: a CPU without BMI1 runs it as BSF, with the same
  result for the non-zero inputs `__builtin_ctz` allows. There is no LZCNT in the whole binary, which a
  pre-BMI CPU would run as BSR with a different result (`results/lzcnt_tzcnt.txt`).

The dynamic libraries are libzstd, libdeflate, libstdc++, libgomp, libm, libgcc_s and libc (`results/ldd.txt`).

## The launcher

`scripts/launcher_test.sh` runs a copy of `protal_launcher` that reads a fake `/proc/cpuinfo`, next to stub
binaries that print their names (`results/launcher.txt`):

| fake `/proc/cpuinfo` | runs |
|---|---|
| every x86-64-v3 flag | `protal_avx2` |
| any one of avx, avx2, bmi1, bmi2, fma, f16c, movbe, abm left out | `protal_baseline` |
| avx512f but no avx2 (the flags are matched as whole words) | `protal_baseline` |
| no file | `protal_baseline` |
| ARM (`Features`, no `flags` line) | `protal_baseline` |
| every x86-64-v3 flag, `PROTAL_NO_AVX2=1` | `protal_baseline` |

The e2e test `LauncherTest` checks that the launcher prefers the binaries next to it over those on `$PATH`.

## Not checked, and an observation

- **The baseline build on a CPU without AVX2.** That needs an emulator, for example
  `qemu-x86_64 -cpu Westmere build/protal ...` (or Intel SDE's `-snb`), and neither is installed here. The
  disassembly and the dispatch code say it would run. A run under an emulator would confirm it, and could
  join the CI.
- **WFA2's AVX2 extend kernels are never built.** `wavefront_extend_kernels_avx.c` compiles its AVX2 and
  AVX-512 kernels only under `#if __AVX2__`, and `wfa_lib` is built once, for plain x86-64, and linked into
  `protal_avx2` as well. So even the AVX2 build extends wavefronts with the portable kernel. That is safe but
  leaves speed unused. A second `wfa_lib` built with `isa_avx2` for `protal_avx2` would use them; whether that
  is measurable was not tried.
- Intel Macs have no `/proc/cpuinfo`, so the launcher runs the baseline there; protal's own AVX2 dispatch
  still applies.

## How it was run

```bash
bash scripts/run.sh ~/mt-work/samroomhead/src/build/protal   # isa.txt, lzcnt_tzcnt.txt, ldd.txt, launcher.txt
```
