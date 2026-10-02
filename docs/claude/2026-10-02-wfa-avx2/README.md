# AVX2 in WFA2, and one binary for every CPU

- **Date**: 2026-10-02.
- **Code**: branch `audit-fixes` at `f07acb3`; the experiment is `scripts/wfa_avx2.patch` on it (not committed).
- **Later the same day**: `scripts/wfa_avx2.patch` was removed, as WFA2 stays without AVX2 (the user's decision).
  It is in the history: `git show 1fced40:docs/claude/2026-10-02-wfa-avx2/scripts/wfa_avx2.patch`. `wfabench.cpp`,
  `wfabench_run.sh` and `wfa_real.sh` need it. The one binary was done with `target_clones` instead
  ([report](../2026-10-02-one-binary/README.md)).
- **Machine**: WSL2 Ubuntu 24.04 on an Intel Core Ultra 7 258V (4 fast and 4 low-power cores, 6 vCPUs), gcc 13.
  Other sessions kept the load at 6–17 for most of the time, which matters for the timings (below).
- **Questions**: (1) could every AVX2 use of protal be chosen at run time, as the syncmer scan is, so that one
  binary serves every CPU? (2) switch on AVX2 for WFA2-lib, at run time if possible.

## Summary

1. **WFA2 with AVX2 does not make protal faster on this CPU, so it is not switched on.**
   - The run-time choice works: `lib/wfa2-extend` builds the extend kernels twice and picks one when protal
     starts, and the outputs stay the same.
   - Upstream's AVX2 kernels were 1–14% *slower* than its portable ones, at every divergence from 1% to 30%.
     Their ends-free kernel also ended some alignments on another diagonal of the same score, which changed
     SAM records.
   - An AVX2 kernel of protal's own gives the same alignments by construction and runs as fast as the
     portable one (0.88–1.06; 1.00–1.02 on a mix like protal's).
   - The extend kernels are a third of WFA2's instructions in protal's use, and WFA2 is only part of the
     alignment. The performance report of 2026-09-30 (§4.3) found no wall-clock change either with WFA2
     built for x86-64-v3.
2. **One binary is possible, and the baseline build is close to being it already.**
   - The baseline `protal` chooses AVX2 at run time where it was found to matter: the syncmer scan, sequence
     packing, zlib-ng, libdeflate and zstd.
   - `protal_avx2` adds the compiler's x86-64-v3 code in the rest of `main.cpp`. On 2026-09-29 that was worth
     3–12% on the syncmer scan's scalar paths, and nothing once the scan used AVX2 at run time (590–620 ns
     per read in both builds).
   - Today's whole runs could not settle it: other sessions' load made even CPU seconds vary by 25%.
   - If a run on a quiet node (`scripts/isa_compare.sh`) shows 3% or less, ship the baseline build alone and
     drop `protal_avx2` and the launcher. If it shows more, mark the few hot functions that gain with
     `target_clones`. Compiling all of protal twice into one binary is not safe: the two copies share inline
     and template functions under the same symbols, and the linker keeps one of each.

## WFA2 with AVX2 at run time

### How upstream chooses

`wavefront_extend_kernels.c` calls its AVX2 kernels (`wavefront_extend_kernels_avx.c`) under `#if __AVX2__`, so
the choice is made when WFA2 is compiled. protal builds `wfa_lib` for plain x86-64, for `protal` and
`protal_avx2` alike, which leaves them out of both binaries. protal aligns with `WFAlignerGapAffine`, full
alignment, `MemoryHigh` (never the bialignment and its own AVX2 code) and ends-free, so the kernel that matters is
`wavefront_extend_matches_packed_endsfree`. The end-to-end ones are only used through the tests.

### The run-time choice (`scripts/wfa_avx2.patch`)

Upstream's files stay unmodified:

- `lib/wfa2-extend/extend_scalar.c` compiles upstream's `wavefront_extend_kernels.c` with its three packed extend
  functions renamed (`protal_wfa_scalar_extend_*`).
- `extend_avx2.c`, compiled with `-mavx2`, holds the AVX2 kernels (`protal_wfa_avx2_extend_*`).
- `extend_dispatch.c` defines upstream's three names and calls one kernel or the other. A constructor decides
  before `main` and before any thread exists, with `__builtin_cpu_init()` and `__builtin_cpu_supports("avx2")`.
  `protal_wfa_extend_use_avx2(bool)` lets tests compare both paths.
- `lib/wfa2-lib.cmake` does this on x86 with GCC or Clang; elsewhere it builds upstream's two files as before.
- In the baseline binary, WFA2's only AVX2 instructions are in those kernels (`objdump`).

### Upstream's AVX2 kernels

The new unit test (`WFA2Wrapper.TheAvx2ExtendKernelsAlignAsThePortableOnes`: 4,500 150 bp reads ends-free with
and without X-drop, 60 long reads end to end and ends-free, with and without wf-adaptive, each aligned with both
kernels) failed with upstream's kernels: reads running past their gene's end got other operations at the same
score (`results/upstream_endsfree_differs.txt`).

The cause is in the ends-free kernel. It extends 8 diagonals at a time and checks for the alignment's end on the
diagonals it extends further one by one *during* that loop, and on the others afterwards. The portable kernel
extends and checks each diagonal in turn from the lowest. Where several diagonals reach the end at the same
score, they stop on different ones. With the check moved after all diagonals are extended, the test passed. Whole
runs at `-t 1` (paired-end 500k and 1k pairs, Nanopore 3 Mb) then gave the same SAM text and outputs, and at
`-t 6` (5M pairs, Nanopore 90 Mb, three runs each) the same sorted records (`results/real_runs_upstream_kernels.txt`).

That version is slower. `scripts/wfabench.cpp` aligns the same cases with protal's wrapper, as protal aligns
reads. It switches the kernel every 1,000 short or 20 long alignments and adds up the thread's own CPU time, so
that both kernels see the same state of the machine. The rounds stayed within ±5% at a load of 11–13
(`results/wfabench_upstream_kernels.txt`):

| reads | 1% | 3% | 5% | 10% | 15% | 20% | 30% | mix (0–15%) |
|---|---|---|---|---|---|---|---|---|
| 150 bp, AVX2 / portable | 1.06 | 1.13 | 1.13 | 1.10 | 1.12 | 1.10 | 1.11 | 1.14 |
| 2–10 kb, AVX2 / portable | 1.01 | 1.08 | 1.06 | 1.05 | 1.03 | 1.06 | 1.08 | 1.10 |

The kernel gathers 4 characters of 8 diagonals and extends every diagonal whose 4 characters matched with the
portable kernel anyway. With reads close to their genes, most do, so the gathers are extra work.

### protal's AVX2 kernels

`extend_avx2.c` as in the patch: upstream's portable kernels, loop for loop (the same diagonals in the same
order, the same checks for the end, empty diagonals skipped as there), with the matches on a diagonal counted 32
characters per step (`_mm256_cmpeq_epi8`, `movemask`, `ctz`) instead of 8. That gives the same alignments by
construction, and the unit test passes.

Reading 32 characters is safe. WFA2-lib keeps 64 bytes after the pattern and after the text
(`WF_SEQUENCES_PADDING`), an extension starts at most at a sequence's end, and the end characters of pattern and
text differ (`!`, `?`), so a diagonal's matches stop there at the latest.

Speed, the same way as above (`results/wfabench_protal_kernels.txt`):

| reads | 1% | 3% | 5% | 10% | 15% | 20% | 30% | mix |
|---|---|---|---|---|---|---|---|---|
| 150 bp, AVX2 / portable | 1.01 | 0.96 | 0.97 | 1.03 | 1.04 | 1.02 | 1.06 | 1.02 |
| 2–10 kb, AVX2 / portable | 1.04 | 0.94 | 0.88 | 0.99 | 1.05 | 1.02 | 1.01 | 1.00 |

### Where WFA2's time goes

callgrind of `wfabench` (one round of the mix, both kernels in turn; `results/callgrind_wfabench.txt`):

| part | instructions |
|---|---|
| extend kernel (each of the two about half) | 32.5% |
| computing the next wavefront (`wavefront_compute_affine`) | 33.5%, of which the vectorisable kernel 11.3% |
| wf-adaptive heuristic (`wavefront_heuristic_cufoff`) | 18.0% |
| checks for the alignment's end (`wavefront_termination_endsfree`) | 9.5% |

The two extend kernels executed about the same number of instructions (1.50 G portable, 1.45 G AVX2). Even an
AVX2 compute kernel at twice the speed would save about 5% of WFA2's time. On 2026-09-30, WFA2 built for
x86-64-v3 (upstream's kernels and the compiler's AVX2 everywhere in it) saved 0.54 G of 31.8 G instructions and
no measurable time: 21.70 → 21.85 s and 9.03 → 9.19 s, pinned, four alternated runs
(`../2026-09-30-performance-round2/README.md`, §4.3).

### So

Not switched on. The patch stays here (`scripts/wfa_avx2.patch`, with protal's kernels and the unit test) in case
another CPU, for example a server's, shows a gain with `scripts/wfabench_run.sh`. The whole-run checks above were
made with the upstream-based version. Before protal's kernels were committed, they would need the same
`-t 1`/`-t 6` comparisons and the unit and e2e tests.

## One binary for every CPU

What decides it is how much `protal_avx2`'s x86-64-v3 code gains over the baseline build, whose run-time AVX2 kernels
are already in place:

- 2026-09-29 (`../2026-09-29-performance-profiling/README.md`, the syncmer scan, ns per read, 500k reads):

  | build | window by window | one by one | AVX2 at run time |
  |---|---|---|---|
  | `-march=x86-64` | 2,690–2,890 | 1,410–1,530 | 600–620 |
  | `-march=x86-64-v3` | 2,610–2,780 | 1,240–1,380 | 590 |

  So x86-64-v3 helped the scalar paths (3–12%) and left the run-time AVX2 path the same.

- Today (`scripts/isa_compare.sh`, both builds of `f07acb3`'s code, alternated; stopped after two of five
  repetitions as the load reached 15; `results/baseline_vs_avx2_partial.tsv`):

  | case | CPU seconds, baseline | CPU seconds, `protal_avx2` |
  |---|---|---|
  | 500k pairs, 1 thread | 19.56, 20.30 | 18.78, 19.89 |
  | 5M pairs, 6 threads | 172.1, 174.7 | 145.9, 182.5 |
  | Nanopore 90 Mb, 6 threads | 48.2, 51.5 | 50.9, 51.1 |

  The one-thread runs suggest 2–4% less CPU for `protal_avx2`. The six-thread ones swing by 25% from run to
  run (threads land on fast or low-power cores), so they show nothing.

Ways to one binary, in order of effort:

1. **The baseline build alone.** If `isa_compare.sh` on a quiet node shows 3% or less, install the baseline build
   as `protal` and drop `protal_avx2`, `protal_launcher` and their places in `justfile`, the conda recipe, CI,
   `docs/installation.md` and the e2e `LauncherTest`.
2. **`target_clones` on the functions that gain.** GCC builds each such function for `default` and
   `arch=x86-64-v3`, and the loader picks one (an ifunc). Profiling the two builds side by side would show which
   functions gain; they become separate calls, so only large hot functions qualify.
3. **All of `main.cpp` twice in one binary** (a `main` that chooses a copy) is not safe without care: the copies
   share the inline and template functions of protal, std and the header libraries under the same symbols, the
   linker keeps one of each, and x86-64-v3 code could then run on the baseline path. It would need one copy
   linked with its symbols made local.

## How it was run

```bash
git apply scripts/wfa_avx2.patch                   # on f07acb3, in a copy of the tree
bash ../2026-10-01-multithreading-audit/scripts/followup2/build_wt.sh wfa PATCH f07acb3 "protal_tests protal protal_avx2"
./tests/protal_tests --gtest_filter='WFA2Wrapper.*' # in its build folder
bash scripts/wfabench_run.sh                       # wfabench: the mix, then 1-30% divergence
bash scripts/cg.sh                                 # callgrind of wfabench
bash scripts/wfa_real.sh                           # whole runs, -t 1 and -t 6, old against new builds
bash scripts/isa_compare.sh                        # baseline against protal_avx2 (on a quiet machine)
```

`wfa_real.sh` ran with the upstream-based version of the patch (`results/real_runs_upstream_kernels.txt`).
