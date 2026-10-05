# One binary for every CPU, with `target_clones`

- **Date**: 2026-10-02.
- **Code**: branch `audit-fixes` at `d11381f`, plus this change (`96d9c3d`). The builds
  compared are listed under [How it was run](#how-it-was-run).
- **Machine**: WSL2 Ubuntu 24.04 on an Intel Core Ultra 7 258V (4 fast and 4 low-power cores, 6 vCPUs), GCC 13.3.
- **Data**: the database of the v0.7.1 benchmark world (`~/bench071/V071/protal_db`, 288 MB). Instructions were
  counted (callgrind, one thread) on the first 100,000 pairs of `rl150_p500000_s_1`, on Nanopore reads of 3 Mb
  (`ont_b3000000_s_1`) and on profiling the 500k-pair SAM. Whole runs were timed on 500k pairs, Nanopore 90 Mb
  and the same profiling, one thread.
- **Questions**: (1) as proposed in [the previous report](../2026-10-02-wfa-avx2/README.md) (way 2): one binary
  that switches between AVX2 and plain code at run time, with `target_clones`, in as many places as it is not
  detrimental; (2) would x86-64-v2 as the baseline bring gains beyond excluding old CPUs?

## Summary

1. **protal is one binary now.** `protal_avx2` and the launcher are gone. Thirteen hot functions are marked
   `PROTAL_CLONE_V3`. GCC compiles each of them twice, for plain x86-64 and for x86-64-v3, and the loader picks
   the copy the CPU supports.
   - On this CPU the binary executes 3.5% fewer instructions than the plain x86-64 build when aligning short
     reads, 1.9% fewer on Nanopore reads and 0.3% fewer when profiling.
   - `protal_avx2` has all of `main.cpp` built for x86-64-v3. The clones save more than it does on alignment
     (3.2% and 1.5% for `protal_avx2`), and less on profiling (0.7%).
   - Every build gives the same outputs.
2. **Both paths were checked by running them.**
   - On this CPU, the x86-64-v3 copies run.
   - With libgcc's CPU feature bits cleared under gdb before the first function is resolved, only the plain
     copies run, and the outputs stay the same.
   - No plain copy contains an AVX instruction.
3. **Tests pass.** In Release, 274 unit tests and 121 e2e tests on the mini database pass. Under ASan +
   UBSan (CI's sanitizer job), 274 of 274 unit tests pass.
4. **x86-64-v2 as the baseline: not worth it.**
   - On its own it saves 1.9%, 0.4% and 1.0% of the instructions.
   - With the clones in place, almost all of its alignment gain is already in the cloned functions on any
     x86-64-v3 CPU. About 0.6% would remain for profiling, which spends that on popcount in the sparse maps.
   - The cost: protal would stop with "Illegal instruction" on CPUs without SSE4.2 and POPCNT. Those are Intel
     before Nehalem (2008), AMD before Bulldozer (2011), and virtual machines that show a generic CPU, such as
     QEMU's `qemu64` or `kvm64`.
5. **Time.** Compare the fastest of five alternated runs at one thread.
   - Aligning 500k pairs: 9.62 s for the clones against 10.07 s for the plain build (−4.5%) and 9.75 s for
     `protal_avx2`.
   - Nanopore 90 Mb: 23.78 s against 24.48 s (−2.9%).
   - Whole runs on this laptop vary by up to 50% with the core they land on, so these are indications. The
     instruction counts are the measure.

## What changed

- `src/Utilities/TargetClones.h` defines `PROTAL_CLONE_V3` as `__attribute__((target_clones("default",
  "arch=x86-64-v3")))`.
  - This needs GCC 12 or later on x86-64 Linux with glibc, which has the ifunc mechanism GCC uses.
  - Elsewhere, or with `-DPROTAL_NO_CLONES`, the macro is empty and every function is compiled once, for x86-64.
- `-ffp-contract=off` is added to `isa_baseline`, the ISA flags of every protal target.
  - GCC contracts `a*b+c` into an FMA instruction where the target has one: by default for C++ even with
    `-std=c++20`, unlike C.
  - Without the flag, the x86-64-v3 copies could round differently from the plain ones.
  - Plain x86-64 has no FMA, so the flag changes nothing in the plain code.
  - The outputs of `protal_avx2` had happened to be the same anyway.
- Removed:
  - the `protal_avx2` target, `isa_avx2`, `protal_launcher`, `PROTAL_NO_AVX2` and the e2e `LauncherTest`;
  - the launcher's names in `Options::ProgramName`;
  - `just avx2`.
- `just install`, the conda recipe, CI and `.vscode/tasks.json` build and install `protal` itself. `just install`
  removes an earlier install's `protal_baseline` and `protal_avx2`.
- Docs updated:
  - `docs/installation.md`: targets, install layout, and a new section "One binary for every CPU". It also
    fixes `just install DIR`: just passes recipe arguments by position, so the documented
    `just install prefix=DIR` installed into a folder named `prefix=DIR`;
  - `docs/running.md`: the environment variable is gone;
  - `docs/development.md`: what CI builds;
  - `docs/README.md`.
- The website (`main.php?site=documentation`) says nothing about AVX2 builds or the launcher, so it does not
  go out of date.

### The marked functions

These were chosen from callgrind's per-function counts of the plain build against `protal_avx2`. Each is a
function that saves instructions as x86-64-v3 code and that does enough work per call:

| function | what x86-64-v3 does there (pe100k unless noted) |
|---|---|
| `index_codec::detail::DecodeChunk` | decoding the index at start-up (23% of the short-read run's instructions here): −376 M |
| `KmerLookupSM::GetFromLookup` | `__builtin_popcount` of the flex k-mer similarity as one instruction instead of libgcc's software popcount: most of the −268 M of `__popcountdi2` |
| `GenomeLoader::GetGenome`, `HasGene`, `GeneLength` | the same for tsl's sparse maps (popcount on every lookup); mostly profiling |
| `AnchoredAligner::Align`, `Flank` | vectorised loops: −201 M |
| `SimpleAlignmentHandler::AlignAnchor` | −50 M (with what it inlines) |
| `CigarANI` | `count_if` over the CIGAR, 32 bytes at a time: −26 M |
| `ArtoSAM`, `ProtalPairedOutputHandler::operator()` | −16 M, −8 M |
| `Profiler::PrepareMAPQ` | −9 M when profiling |
| `ChainAnchorFinder::operator()` | x86-64-v3 itself saves little there (13 M in `protal_avx2`), but the read-order sort's comparator runs 89 M cheaper inside the x86-64-v3 copy: −48 M overall, −33 M on Nanopore reads, against not cloning it |

Per-function tables: `results/functions_pe100k.txt`, `functions_ont3M.txt`, `functions_prof500k.txt`. Clone
suffixes are merged there, so a function's two copies count as one.

### What a clone costs, and the one that was dropped

- A cloned function is never inlined into its callers. Each call goes through the PLT, an extra indirect jump.
- GCC inlines a function called from one place whatever its size. Cloning a caller makes two places, so large
  callees can stay separate and run as plain x86-64 code in both copies.
  - Cloning `ChainAnchorFinder::operator()` leaves `FindSeeds`, `ExtendAnchor` and `ReverseComplementInto`
    separate (+18 M).
  - It still pays, because of the sort: the build without that clone executed 29.04 G, 25.52 G and 16.94 G
    instructions (`cl3` in `results/callgrind_totals.tsv`), 0.1–0.2% more than with it.
- `RecordEvidenceCollector::NoteRecord` was dropped.
  - Its copy saved 28 M, but its `ForEachAlternative` was no longer inlined (+41 M).
  - The totals were the same with and without it (16.930 G and 16.932 G when profiling).
- Functions with function-local statics are not marked.

## Instructions

callgrind, one thread, `-t 1 --no_qcmsa`, on four builds of the same code (`results/callgrind_totals.tsv`):

| workload | plain x86-64 | x86-64-v2 | `protal_avx2` (x86-64-v3) | clones |
|---|---|---|---|---|
| 100k pairs | 30.05 G | 29.49 G (−1.85%) | 29.08 G (−3.20%) | 29.00 G (−3.49%) |
| Nanopore 3 Mb | 25.98 G | 25.89 G (−0.36%) | 25.61 G (−1.46%) | 25.49 G (−1.91%) |
| profiling 500k pairs | 16.98 G | 16.81 G (−1.03%) | 16.87 G (−0.67%) | 16.93 G (−0.29%) |

"clones" is the change as committed (`cl2` in the results), which is also the build timed below. `cl` also
cloned `NoteRecord`, and `cl3` did not clone `ChainAnchorFinder::operator()`.

- The clones get the gains of `protal_avx2` on alignment and a little more:
  - `DecodeChunk` saves 376 M as a clone against 255 M in `protal_avx2`;
  - the read-order sort inside `ChainAnchorFinder` saves 89 M;
  - the anchored aligner saves 201 M in both (`Align`, `Flank` and `Between`);
  - `AlignAnchor` saves 50 M against 58 M.
- `protal_avx2` does better than the clones when profiling (−0.67% against −0.29%). Profiling's popcounts are
  spread over many small sparse-map lookups (`Taxon::AddSam`, `MicrobialProfile::TaxonOf`, tsl's own
  functions), and x86-64-v2 does best there (−1.03%).
- The rest of the profiling is no better as x86-64-v3: `VariantHandler::AddAlignment` executes 61 M *more*
  instructions in `protal_avx2`.

The index decoding is a large share here because the inputs are small. It is a fixed cost per run and does
not grow with the reads.

## Time

Whole runs at one thread (`scripts/oneb_time.sh`, `results/timing_runs.tsv`). The four builds alternated, in
five rounds, niced. Times are protal's own stage times ("Aligning reads took", "Profiling took"), in seconds. The
load from other sessions was 1.0–2.8.

The same run of the plain build took 10.1 to 13.2 s, depending on whether it landed on a fast or a low-power
core. The fastest of the five runs is the closest to a like-for-like comparison:

| case | plain x86-64 | clones | `protal_avx2` | x86-64-v2 |
|---|---|---|---|---|
| 500k pairs, fastest | 10.07 | 9.62 (−4.5%) | 9.75 (−3.2%) | 9.83 (−2.4%) |
| 500k pairs, median | 10.91 | 10.33 | 10.21 | 10.08 |
| Nanopore 90 Mb, fastest | 24.48 | 23.78 (−2.9%) | 24.08 (−1.6%) | 24.29 (−0.8%) |
| Nanopore 90 Mb, median | 25.38 | 24.39 | 25.36 | 24.94 |
| profiling 500k pairs, fastest | 1.12 | 1.13 | 1.16 | 1.12 |

- The fastest runs are in the same order as the instruction counts, and of the same size.
- Differences of 1–2% are within this machine's noise. So are the medians, which mix fast and slow cores.
- Profiling takes 1.1 s here, too short to show 0.3–1%.

## x86-64-v2 as the baseline?

x86-64-v2 adds SSSE3, SSE4.1, SSE4.2, POPCNT and CMPXCHG16B to plain x86-64. On its own (table above) it
saves:

- 1.85% when aligning short reads. Most of that is popcount (−277 M) and the anchored aligner's loops (−157 M),
  plus `AlignAnchor`, `GetFromLookup`, `ArtoSAM` and `CigarANI`.
- 0.36% on Nanopore reads, almost all popcount.
- 1.03% when profiling: popcount (−99 M) and `SequenceRangeHandler::CoveredPortion` (−30 M).

With the clones in place, a CPU with x86-64-v3 already runs those functions as x86-64-v3 code:

- On alignment, what x86-64-v2 would add outside the clones is below 0.05%: `mate_guidance::BestDiagonal`
  −4 M, tsl's sparse set −1 M.
- On profiling it would add about 0.6%, the popcounts in the many small lookups.
- Only a CPU that has x86-64-v2 but not x86-64-v3 would get the full 0.4–1.9%: Intel before Haswell (2013),
  later Pentium, Celeron and Atom chips, and AMD before 2015.

The cost: the binary would no longer start, with "Illegal instruction", on CPUs without SSE4.2 or POPCNT:

- Intel Core 2 and older (before 2008-2009);
- AMD K10 (Phenom II, Opteron until 2011, which lack SSSE3 and SSE4.1);
- virtual machines with a generic CPU model, such as QEMU's `qemu64` and `kvm64`. That was Proxmox's default
  until version 8, and is still found on clusters.

So plain x86-64 stays the baseline. If profiling's popcounts ever matter, marking `Taxon::AddSam` and
`MicrobialProfile::TaxonOf` would be the targeted way.

## Checks

- **Symbols** (`results/clones.txt`):
  - the binary has a `.default`, an `.arch_x86_64_v3` and a `.resolver` symbol for each marked function: 14
    sets for 13 functions, because `ProtalPairedOutputHandler` has two instantiations (`<true>` and `<false>`);
  - no `.default` copy uses a `ymm` register;
  - `ymm` registers appear only in the x86-64-v3 copies and in the kernels that choose AVX2 themselves (syncmer
    scan, packing, zlib-ng, libdeflate, zstd).
- **The plain path, run** (`scripts/oneb_nocpu.sh`, `results/dispatch_gdb.txt`). There is no emulator on this
  machine, so the script simulates a CPU without x86-64-v3:
  - gdb stops at libgcc's `__cpu_indicator_init`. The first ifunc resolver calls it while the loader
    relocates protal, before `main` and before any constructor.
  - gdb lets it fill in the CPU's features, then clears the feature words of `__cpu_model` and
    `__cpu_features2`.
  - From then on, every resolver and protal's own `__builtin_cpu_supports("avx2")` see a CPU without AVX2,
    POPCNT or SSE4. zlib-ng, libdeflate and zstd ask the CPU themselves.
  - Breakpoints count the calls of four functions' copies on 1,000 pairs:

    | copy | `CigarANI` | `GetFromLookup` | `AlignAnchor` | `DecodeChunk` |
    |---|---|---|---|---|
    | plain | 354 | 14,506 | 1,751 | 52 |
    | x86-64-v3 | 0 | 0 | 0 | 0 |

  - The outputs equal the plain x86-64 build's, SAM text included.
- **The x86-64-v3 path** (`scripts/oneb_v3cpu.sh`): the same run without clearing the bits calls only the
  x86-64-v3 copies (the same counts), with the same outputs.
- **Outputs** (`results/outputs_identical.txt`):
  - every callgrind run (`-t 1`): the SAM text and every other output of the clones build, `protal_avx2` and
    the x86-64-v2 build equal the plain build's;
  - profiling: all 313 files;
  - the first round of the timed runs: 500k pairs, Nanopore 90 Mb and profiling, all the same, SAM text
    included.
- **Tests**:
  - Release, from a clean tree (`check.sh`): 274 unit tests and 121 e2e tests on a freshly built mini database
    pass (`results/final_checks.txt`);
  - CI's sanitizer job (Debug, ASan + UBSan, `detect_leaks=1`, `halt_on_error=1`): 274 of 274 unit tests
    pass, with 62 clones in the test binary (`results/asan.txt`);
  - `just install` into a prefix holding an earlier install's `protal_baseline` and `protal_avx2` removed them
    and installed `protal`, an ELF binary, with the other tools (`results/final_checks.txt`).

## How it was run

The scripts are in `scripts/`. `S` in them is the scratch folder they ran from, with the patches.
`build_wt.sh` and `check.sh` are those of
[the multithreading audit](../2026-10-01-multithreading-audit/scripts/followup2/).

```bash
bash build_wt.sh oneb-ref empty.patch d11381f "protal protal_avx2"   # plain x86-64 and protal_avx2 of d11381f
bash scripts/oneb_builds.sh      # the clones build, and an x86-64-v2 baseline (scripts/v2.patch)
bash scripts/oneb_nm.sh          # clone symbols, ymm registers
bash scripts/oneb_cg.sh          # callgrind of the four builds, three workloads; outputs compared
bash scripts/oneb_nocpu.sh       # the plain path under gdb; scripts/oneb_v3cpu.sh the x86-64-v3 path
bash scripts/oneb_rebuild.sh     # cl2: the clones build without NoteRecord's clone (the final code), measured again
bash scripts/oneb_time.sh        # whole runs of base, cl2, v3, v2, alternated, five rounds; scripts/oneb_timesum.sh
bash scripts/oneb_final.sh       # cl3 (ChainAnchorFinder not cloned): check.sh, callgrind, both paths
CG=0 bash scripts/oneb_final.sh  # the final code again: check.sh, both paths, just install
PATCH=oneb_final.patch bash scripts/oneb_asan.sh   # CI's sanitizer job on the final code
bash scripts/oneb_collect.sh     # results/
A=base B=cl2 bash scripts/oneb_diff.sh pe100k 25   # per-function differences of two builds
```
