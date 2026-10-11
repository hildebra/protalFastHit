# Can WFA2 with Singletrack be made faster? SIMD and other exact changes, measured on protal's own calls

- **Date**: 2026-10-11.
- **Code**: `audit-fixes` at `18950c5` (WFA2-lib: the fork with Singletrack and protal's three patches,
  [2026-10-11-singletrack](../2026-10-11-singletrack/README.md)). In sections 1-5 nothing in protal was changed; the
  variants are scratch copies of the library. Section 6 is the implementation, on branch `wfa-speed` from `0a1beaa`.
- **Machine**: WSL2 on the Core Ultra 7 laptop (AVX2, no AVX-512), cores 0-3, load 1-4 from other sessions.
- **Data**:
  - Database: `~/realworld_s1/protal_db` (six GTDB r226 families, 530 genomes).
  - Samples simulated from 40 of its genomes: paired-end 1M pairs, PacBio HiFi 30 Mb, Nanopore 30 Mb.
  - A scratch protal build (`scripts/WFA2Wrapper2_dump.h`) wrote every WFA2 call of a `-t 1` run: free ends, budget,
    X-drop, read part, reference part.
  - That gave 40,026 calls (pe, all ends-free), 10,682 (PacBio: 10,104 end-to-end pieces between anchors, 578 ends-free
    flanks) and 28,069 (ONT: 25,596 and 2,473).
- **Method**: `scripts/wfareplay.cpp` replays the calls through protal's `WFA2Wrapper2`, with one reused aligner per
  kind and X-drop, as protal makes them. For each variant it gives:
  - thread CPU time, with the variants alternated round by round;
  - callgrind instructions, which are deterministic;
  - every call's outcome, score and CIGAR hash, compared with the unchanged library.
- **Question (the user)**: can WFA2, the new Singletrack version, be improved further, e.g. with SIMD/AVX2? The user
  added that it was probably only ever optimized for PacBio long reads.
- **Follow-up (the same day, the user)**: "implement all four and the slab fast path, also check the per-step
  bookkeeping". Done on branch `wfa-speed` (from `0a1beaa`), [section 6](#6-implemented-a-d-the-slab-and-the-bookkeeping-checked):
  every call, unit test and whole run gives the same results as before; WFA2 instructions −26% pe, −46 to −48% PacBio,
  −39 to −43% ONT. Then Singletrack became the default of `--build` (protal 0.8.2, end of section 6).

## Summary

1. **Singletrack's own parts are not where the time goes.**
   - Its M-only backtrace is 2-4% of WFA2's instructions.
   - The time is in the forward pass, which classic and Singletrack share. Only 11-15% of it is the loop compilers
     vectorize (`wavefront_compute_affine_idm`).
   - The rest is the extend kernel (16-33%), the wf-adaptive heuristic (9-20%), per-step bookkeeping (8-27%),
     per-diagonal termination checks (6-12%) and, for long reads, the slab reset at each alignment (2-5%).
2. **Why it looks PacBio-shaped.** The fork keeps every M wavefront and recycles only the gap wavefronts, so it saves
   memory traffic in wide wavefronts. Replayed, Singletrack against classic:

   | Calls | Singletrack against classic |
   |---|---|
   | Short reads | −5% |
   | PacBio | −24% |
   | ONT | −29% |

   Short reads' wavefronts already sit in cache. Most long-read calls are short end-to-end pieces between anchors (about
   70 bp), but the few ends-free flanks carry most of the instructions.
3. **Four exact improvements, all byte-identical** (outcome, score, CIGAR) on all 78,777 calls, in both modes:
   - **a. Vectorize the heuristic's loops** (one line per loop). GCC vectorized none of the wf-adaptive and X-drop
     loops at any target, not even x86-64-v3. The loops run `for (k=wavefront->lo; k<=wavefront->hi; ++k)` while
     storing into an array GCC cannot prove is apart from `*wavefront`, so the trip count "cannot be computed". With
     `lo`/`hi` copied into locals they vectorize.
   - **b. Pre-check before the termination call.** The ends-free extend kernel calls `wavefront_termination_endsfree`
     (`FORCE_NO_INLINE`) for every extended diagonal. It can only stop an alignment once the diagonal has reached the
     end of the text or the pattern, so a two-comparison pre-check skips nearly every call. A TODO in the kernel already
     suggests this.
   - **c. SSE4.1 or AVX2 for the vectorizable loops**, at run time with `target_clones` as protal does for its own hot
     functions. Baseline x86-64 (SSE2) has no 32-bit integer min/max, which these loops are made of.
     - Not by compiling WFA2 for x86-64-v3: `__AVX2__` switches it to its own AVX2 extend kernels, which change ends-free
       ties, as the [2026-10-02 report](../2026-10-02-wfa-avx2/README.md) found. Here too the results differed.
   - **d. Halve the penalties.** protal's 4/6/2 are all even, so every odd score step of WFA2 is empty. WFA2 does not
     normalize by the greatest common divisor. Aligning with 2/3/1, a budget of ⌈max/2⌉ and the score doubled back gives
     the same alignments.
4. **Gains.**

   Instructions inside WFA2 (Singletrack, callgrind):

   | Variant | pe | PacBio | ONT |
   |---|---|---|---|
   | a | −3% | −11% | −9% |
   | a + b | −9% | −24% | −17% |
   | a + b, x86-64-v2 | −13% | −33% | −24% |
   | a + b, x86-64-v3 without WFA2's AVX2 extend kernels | −13% | −38% | −28% |

   Thread CPU time against today's Singletrack (9 alternated rounds, ±5-10% noise):

   | Variant | pe | PacBio | ONT |
   |---|---|---|---|
   | a + b + d | −15% | −12% | −12% |
   | a + b + c (v2) + d | −20% | −25% | −12% |
   | a + b + c (v3) + d | −12% | −17% | −20% |

   In protal at r226 scale, where WFA2 was 31% of short-read and 91% of long-read alignment, that is about −6% for
   short-read alignment and about −20% for long-read alignment, if the r226 mix behaves like these samples.
5. **What does not pay: more SIMD in the extend kernel.** Its cost is per diagonal (an 8-byte compare that usually stops
   in the first block, across wide wavefronts), not long match runs.
   - A 32-byte compare per diagonal was neutral on 2026-10-02.
   - WFA2's gather-based AVX2 kernel was slower and changes ties.
   - AVX-512 (the cluster's Zen 4) could not be measured here.
6. **A side finding: WFA2's X-drop is broken for a match score of 0** (upstream too).
   - It converts scores with a "match bonus" of −1, so it ranks the least advanced diagonals best.
   - With protal's default `--x_drop 1000` it provably never prunes a short read: the drop it measures is at most
     (v+h+s)/2, under 1000 for anything below ~2,000 combined bases.
   - So `--x_drop 0` would give the same results, and save its 2.5-3% of short-read WFA2 time.
7. **Recommendation.**
   - Implement a, b and d in protal: a and b as patches 4 and 5 of `lib/wfa2-lib`, d in `WFA2Wrapper2`.
   - Add c as `target_clones` on the heuristic's distance loops and the compute kernel.
   - Check with this replay (identity on every call) and whole-run SAM identity, as for Singletrack.
   - Possibly next: the slab-reset fast path below (−5-8% on long reads), and `--x_drop 0` as default.
   - Done, with the slab fast path: [section 6](#6-implemented-a-d-the-slab-and-the-bookkeeping-checked). `--x_drop`
     stays 1000; the wrapper leaves X-drop out of every alignment it provably cannot prune.

## 1. Where Singletrack's instructions go

Self shares of `wavefront_align`'s inclusive instructions, the unchanged library, Singletrack (`results/profiles.txt`):

| Function | pe | PacBio | ONT |
|---|---|---|---|
| `wavefront_extend_matches_packed_endsfree` (+`_end2end`) | 16.1% | 29.0% (+4.4%) | 17.5% (+9.1%) |
| `wavefront_heuristic_wfadaptive` | 9.3% | 20.4% | 17.9% |
| `wavefront_compute_affine_idm` (the vectorizable kernel) | 10.7% | 14.5% | 13.9% |
| `wavefront_termination_endsfree` | 5.9% | 11.7% | 7.0% |
| compute bookkeeping (`init_ends`, `process_ends`, `fetch_input`, `allocate_output`) | 26.6% | ~8% | ~15% |
| `wavefront_slab_reap_repurpose` (per alignment) | 2.0% | 3.2% | 5.4% |
| `wavefront_backtrace_affine_m_only` (Singletrack's backtrace) | 2.4% | < 1.5% | < 1.5% |
| `wavefront_heuristic_xdrop` | 2.5% | – | – |

After a + b + c (v3) + d, the extend kernel is the largest part: pe 20%, PacBio 55%, ONT 39%. Then come the bookkeeping
(12-32%), the compute kernel (9-10%), the slab reset (PacBio 5.4%, ONT 8.0%) and wf-adaptive (6%).

## 2. The changes

| # | Where | Change | Exact because |
|---|---|---|---|
| a | `wavefront_heuristic.c`, the four `PRAGMA_LOOP_VECTORIZE` loops (wf-adaptive's distances end-to-end, weighted, ends-free; X-drop's scores) | `const int k_lo = wavefront->lo, k_hi = wavefront->hi;` before the loop, `for (k=k_lo;k<=k_hi;++k)` | the loop body never writes `lo`/`hi` |
| b | `wavefront_extend_kernels.c`, both calls of `wavefront_termination_endsfree` | call only if `WAVEFRONT_H(k,offset) >= text_length` or `WAVEFRONT_V(k,offset) >= pattern_length` | the function returns true only in those two cases |
| c | the vectorizable loops (a's, `wavefront_compute_affine_idm`, the init-ends loops) | `target_clones("default","arch=x86-64-v3")`, or v2; WFA2's `#if __AVX2__` kernels stay out, as a clone does not define `__AVX2__` | the same integer operations, wider registers |
| d | `WFA2Wrapper2` | penalties divided by their greatest common divisor (2: 2/3/1), the budget `⌈max/2⌉`, the score doubled | WFA2 stops when the score reaches the budget (`score >= max_alignment_steps`), and every score is even; wf-adaptive works in bases. X-drop's score conversion is not scale-free, but at protal's 1000 it cannot prune, at either scale (finding 6) |

GCC's verdict on the heuristic's loops (`scripts/vecinfo.sh`, protal's flags `-O3 -fopenmp`): unchanged, "not
vectorized: number of iterations cannot be computed" at x86-64, v2 and v3. With a, the end-to-end and ends-free
distance loops vectorize, with 16-byte vectors and with 32-byte ones at v3. The weighted (float) distance loop and
X-drop's score loop still do not.

## 3. Timings

Thread CPU seconds of a whole replay, median of 9 rounds, the four libraries alternated, each running the modes in one
process (`results/timing_9rounds.txt`; earlier 5-round run `results/variants_5rounds.txt`):

| Library | Mode | pe | PacBio | ONT |
|---|---|---|---|---|
| unchanged | classic | 0.091 | 0.156 | 0.241 |
| unchanged | Singletrack | 0.085 | 0.114 | 0.163 |
| unchanged | Singletrack, halved | 0.077 | 0.115 | 0.160 |
| a + b | Singletrack | 0.077 | 0.108 | 0.145 |
| a + b | Singletrack, halved | 0.072 | 0.100 | 0.144 |
| a + b, v2 | Singletrack, halved | **0.068** | **0.086** | 0.144 |
| a + b, v3 without AVX2 extend | Singletrack, halved | 0.075 | 0.095 | **0.130** |

Classic gains as well: a + b + v2 + halved is −21% pe, −19% PacBio and −17% ONT against unchanged classic.

## 4. Not tried, and next

- **Slab reset per alignment.** At each alignment `wavefront_slab_clear` walks every wavefront the aligner ever
  allocated (`wavefront_slab_reap_repurpose`), which after one gene-long flank means thousands of them, to free the few
  the last alignment used.
  - A fast path would free only the wavefronts handed out since the last clear, when no size change or deallocation
    happened. That is exact, since slab state does not affect alignments.
  - It needs care with the free stack, fresh allocations and Singletrack's frees during an alignment.
  - Worth −5-8% of long-read WFA2 time after the changes above.
  - Done: section 6.
- **Per-step bookkeeping** (init/process ends, allocation, input fetch) stays 12-32%. Halving removes half the steps; the rest is many small
  loops per step. Checked in section 6.
- **AVX-512**: a `target_clones` entry `arch=x86-64-v4` for a's loops would be cheap to add. It can only be measured on
  the cluster (`scripts/variants.sh NAME bounds,term,noavxext x86-64-v4` and the replay).
- **The X-drop bug** deserves an upstream note, with the two Singletrack issues
  ([upstream-issues.md](../2026-10-11-singletrack/upstream-issues.md)).

## 5. How it was run

```bash
bash scripts/setup.sh            # scratch build with the dump header, samples, dumping runs (~/wfasimd)
bash scripts/replay.sh 7         # replay: classic, singletrack, halved; identity between modes
bash scripts/profile.sh singletrack
bash scripts/vecinfo.sh          # GCC's vectorization of the hot files at x86-64, v2, v3
bash scripts/run_variants.sh     # variants a, a+b, v2, v3 built (variants.sh) and timed, 5 rounds
bash scripts/run_variants2.sh    # v3 without WFA2's AVX2 extend kernels; callgrind instructions; identity
bash scripts/timing.sh           # 9 alternated rounds
```

The WSL folder `~/wfasimd` was removed afterwards; `setup.sh` remakes it in about 5 minutes.

## 6. Implemented: a-d, the slab, and the bookkeeping checked

Branch `wfa-speed` from `0a1beaa` (audit-fixes after the Singletrack merge), checked against `0a1beaa` itself. The
calls are dumped again by a `0a1beaa` build with the dump header; the simulations are seeded, and they are again the
same 78,777 calls.

### What changed

| # | Where | Change |
|---|---|---|
| a | `wavefront_heuristic.c` | the four loops read `lo`/`hi` into locals (patch 4 in `lib/wfa2-lib.cmake`) |
| b | `wavefront_extend_kernels.c` | `wavefront_termination_endsfree` only for a diagonal at the end of the text or the pattern (patch 5) |
| c | `utils/commons.h`, `wavefront_heuristic.c`, `wavefront_compute_affine.c` | `WFA_CLONE_V3` (`target_clones("default","arch=x86-64-v3")`, off with `-DPROTAL_NO_CLONES`) on wf-adaptive's two distance loops and `wavefront_compute_affine_idm` (patch 6). v3 rather than v2, as protal's own clones; the weighted (float) distance loop and X-drop's loop do not vectorize and are not cloned |
| d | `src/Alignment/WFA2Wrapper2.h` | two aligners (below) |
| e | `wavefront_slab.c`, `.h` | the slab's reset at each alignment: only the wavefronts handed out since the last reset are marked free, and only the part of the free stack the alignment changed is rewritten (patch 7) |

**d in detail.** The wrapper's fast aligner has the penalties divided by their greatest common divisor (2/3/1) and no
X-drop; its budget is `max <= 0 ? max : (max-1)/g+1` (WFA2 stops at `score >= max`, and in the given scale every score
of the fast aligner is a multiple of g), and its score is multiplied back. A second aligner, made when first needed,
has the penalties as given and X-drop, for the alignments X-drop could prune: X-drop compares a state's
`(-(v+h)-score)/2` with the best seen (≤ 0), so it can prune only when read part + reference part + max(budget, 0)
reaches 2 × X-drop. Each has its own classic fallback for Singletrack. The accessors (`Mismatch()`, ...) still give
4/6/2, which the alignment screen's bounds use. In the dumped runs every short-read call is far below the bound (the
largest sum is 379, against 2,000 at `--x_drop 1000`) and long reads align without X-drop, so protal's default runs use
only the fast aligner.

### Checks

| Check | Result |
|---|---|
| Unit tests (`ctest`, 516) | all pass. Two are new: `WFA2Wrapper.HalvedPenaltiesGiveTheSameAlignments` (the wrapper against raw WFA2 with 4/6/2 and X-drop at budgets score−1 to score+2, 61, 0, 1 and none; ends-free and end to end; short and long reads; X-drop 0, 1000, 50, 20; both tracebacks) and `WFA2Wrapper.XDropAboveTheBoundPrunesNothing` (X-drop one above the bound changes nothing) |
| e2e (`tests/e2e`, 154, mini database) | all pass |
| AddressSanitizer + UBSan (Debug, `-DPROTAL_NO_CLONES`: no clone in `libwfa_lib.a`) | every unit test passes; the sanitized replay of all calls is clean and gives the base's results |
| Replay, every call (outcome, score, CIGAR hash) | identical to `0a1beaa` for all 78,777 calls, classic and Singletrack; no Singletrack fallback in either |
| Whole runs, pe, PacBio, ONT, 4 threads, classic and Singletrack | SAM records identical to `0a1beaa`'s (87,341 / 1,805 / 1,632 records); every other output file identical but `misc/cpu.tsv` (CPU times) |

### Gains

Instructions in `wavefront_align` (callgrind, one replay of all calls; the final code, `results/impl/queue4.log`):

| Calls | Mode | `0a1beaa` | `wfa-speed` | Change |
|---|---|---|---|---|
| pe | classic | 1,263 M | 931 M | −26.3% |
| pe | Singletrack | 1,229 M | 908 M | −26.1% |
| PacBio | classic | 2,054 M | 1,067 M | −48.1% |
| PacBio | Singletrack | 1,932 M | 1,049 M | −45.7% |
| ONT | classic | 2,764 M | 1,564 M | −43.4% |
| ONT | Singletrack | 2,505 M | 1,536 M | −38.7% |

Thread CPU seconds of a replay (median of 5 alternated invocations of 9 rounds each, cores 0-3, the final code). The
machine was busier than in queue3's round (the base 6-25% slower there), whose changes for the same code without the
inlined allocator were within 1-4 points of these:

| Calls | Mode | `0a1beaa` | `wfa-speed` | Change |
|---|---|---|---|---|
| pe | classic | 0.0597 | 0.0479 | −19.8% |
| pe | Singletrack | 0.0570 | 0.0441 | −22.6% |
| PacBio | classic | 0.1230 | 0.0841 | −31.6% |
| PacBio | Singletrack | 0.0972 | 0.0670 | −31.1% |
| ONT | classic | 0.1790 | 0.1002 | −44.0% |
| ONT | Singletrack | 0.1184 | 0.0880 | −25.7% |

Classic gains most on long reads: it keeps every wavefront, so its slab held thousands after a gene-long flank, and
the reset walked them all. Classic is what runs by default (`--wfa_traceback auto` follows the database, and databases
without `build_settings.tsv` were built classic). With WFA2 at 31% of short-read and 91% of long-read alignment at
r226 scale (summary item 4), that is about −6% of short-read alignment and −23 to −40% of long-read alignment, if the
r226 calls behave like these. The whole runs here are too short to time (0.5-4.5 s, mostly loading the index).

**Singletrack the default of `--build` since 0.8.2** (the user, the same day, after this section). After a-e
Singletrack still needs fewer instructions than classic (pe −2.4%, PacBio −1.6%, ONT −1.8%, the first table) and less
thread time (−8%, −20%, −12%, the second). `--build` now aligns with it unless `--wfa_traceback classic`, so a new
database records `singletrack` and its runs follow. A run on a database built before keeps that database's traceback
(classic) and gets no warning. Singletrack chooses otherwise among equally good alignments in 47-49% of ONT records
([2026-10-11-singletrack](../2026-10-11-singletrack/README.md)); a database built with it trains its models on its own
alignments, but no r226 build with Singletrack has been evaluated yet.

### The per-step bookkeeping

Self instructions after a-e (before the free-stack fix below), as a share of `wavefront_align` (Singletrack;
`results/impl/bookkeeping.txt`):

| Function | pe | PacBio | ONT |
|---|---|---|---|
| `wavefront_compute_init_ends` (fills the inputs' missing elements with NULL) | 11.4% | 4.9% | 7.8% |
| `wavefront_compute_process_ends` (trims the outputs' out-of-range ends) | 9.9% | 4.2% | 6.7% |
| `wavefront_compute_allocate_output` + `wavefront_slab_allocate` + `wavefront_init` + `resize_null__victim` | 15.0% | 6.0% | 10.0% |
| `wavefront_compute_fetch_input` + `wavefront_compute_affine` + `limits_input` | 9.8% | 4.2% | 6.3% |
| together | ~46% | ~19% | ~31% |

With wavefronts as narrow as wf-adaptive keeps short reads' (a few dozen diagonals), a step's fixed work is as large
as its kernels'. At the line level it is many small loops and calls: init_ends' two NULL-fill loops (a few elements per
input and side: WFA2 computes each output over the union of its inputs' ranges, so every input needs NULLs where it is
narrower), and the trims walking the NULL ends of the gap wavefronts. What was tried and found:

- **Inlining across WFA2's files** (all its sources with `-flto` into the replay, `scripts/impl/lto.sh`): the same
  results, but 1-2% more instructions and 0.4-1.8% more time in every sample and mode. Not adopted; a unity build
  would give the same inlining.
- **The slab reset's free stack**, found in this profile: the first version of e copied the whole wavefront list into
  the free stack at each reset, `memcpy` 3.4% (PacBio) and 4.6% (ONT) of classic's instructions. Only pops go below
  the stack's top during an alignment, so the stack's first `free_intact` entries are still the list's, and only the
  rest is copied now. Its allocator then had to be marked inline: GCC had stopped inlining it for the extra comparison
  (pe +1.4% instead of +0.4%). Against the build before: classic −2.7% (PacBio) and −3.7% (ONT) of WFA2's
  instructions, Singletrack −0.1 to −0.2%, pe +0.2-0.4% (the comparison per allocation).
- **Not done**: computing each output wavefront over its own inputs' range instead of the union would remove most NULL
  fills and trims, but that is a redesign of WFA2's compute step and its kernels, not a patch.

### How it was run

The scripts are in `scripts/impl/`, the logs in `results/impl/`. They read from one scratch folder: `base.tar` (`git
archive 0a1beaa`), `new.tar` (the worktree's tracked files), `scripts/WFA2Wrapper2_dump.h`,
[`compare_sams.py`](../2026-10-11-singletrack/scripts/compare_sams.py) and `scripts/impl/wfareplay.cpp` (section 5's,
with four decimals in its timings).

```bash
bash scripts/impl/setup.sh       # ~/wfaspeed: base (0a1beaa + dump header) and new (the worktree's files) built,
                                 # unit tests, samples, base runs dumping their calls
bash scripts/impl/queue1.sh      # the fixed test rebuilt; replay.sh (identity, timing, callgrind); lto.sh; identity.sh
                                 # (whole runs, e2e)
bash scripts/impl/queue2.sh      # lto.sh again (its source list fixed), asan.sh
bash scripts/impl/queue3.sh      # free_intact: rebuilt, unit tests, replay, against the build before, whole runs, ASan
bash scripts/impl/queue4.sh      # the allocator inlined again: the same
bash scripts/impl/bookkeeping.sh # self and line-level shares of the new build's callgrind profiles
bash scripts/impl/xdropsum.sh    # the dumped calls against the X-drop bound
```

The WSL folder `~/wfaspeed` was removed afterwards.
