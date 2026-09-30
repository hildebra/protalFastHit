# Performance profiling on simulated paired reads

2026-09-29. Branch `performance` at `995c4f1` (identical to `audit-fixes` at that commit; the
uncommitted work in the main checkout is not included). Binary: `protal_avx2`, Release (`-O3`,
`-march=x86-64-v3`); instruction profiles use the same flags plus `-g -fno-omit-frame-pointer`.

**Machine:** WSL2 Ubuntu 24.04 on an Intel Core Ultra 7 258V (4 performance and 4 low-power cores,
one thread each, 12 MB L3), 15 GB RAM, gcc 13, valgrind 3.22; no `perf` (no hardware counters).
Other sessions ran simulations and builds on the same machine throughout (load average 7–14 on 8
cores). Wall-clock times are therefore inflated and noisy; the comparisons below were run back to
back or alternated, and the instruction counts (callgrind) do not depend on load.

## Data

| Name | What | Size |
|---|---|---|
| `db64` | protal database of the 64-species synthetic GTDB world `~/audit4/gtdb64` (192 genomes) | 1.7M index values, `database.protal` 8.6 MB, built in 31 s |
| `db900` | protal database of the 900-species world `~/tune/world` (GTDB-like genus sizes from `gtdb_like_lineages.py`, 2700 genomes, 101,081 representative genes) | 26.9M index values, `database.protal` 125 MB, built in 391 s at 8 threads |
| `s64` | 1M pairs, ART HS25 2×150, fragment 350±50, 20 species of the 64-world with strains (`--strains_per_species 0.3,0.1`, seed 11) | 197 MB `.fq.gz` |
| `w900` | the same for 60 species of the 900-world; 53% of pairs align (genomes are ~270 kb, ~40% marker genes) | 197 MB |
| `mix` | every 20th `w900` pair interleaved with 950k ART pairs from 40 Mb of random sequence: 5% of pairs from database species, as in a real metagenome whose marker genes are a few % of each genome | 201 MB |

The synthetic genomes are mostly marker genes, so `s64` and `w900` stress alignment; `mix` has the
share of database reads of a real sample and is the realistic case.

## Summary

For the realistic `mix`, 8 threads align 1M pairs in 21 s with gzipped input and in 8.5 s with the
same reads uncompressed: the FASTQ reader inflates gzip inside `omp critical(reader)`, which caps
protal at about 200k pairs/s whatever the thread count. That is the largest gap found. After it
come syncmer extraction (40% of the per-read instructions on `mix`), WFA alignment (53% on
marker-rich reads, and on `mix` about 70% of it spent on reads that cannot align), TLB misses
in the 3 GB key map (huge pages cut seeding by 30–38% without a code change), and fixed start-up
costs: 3 GB and 2–3 s single-threaded for the key map of any database, and 229k C++ exceptions
while parsing the model. Nothing here was measured on the real GTDB database, where lookups and
alignments per read will be larger (see [Gaps](#gaps)).

| # | Finding | Evidence | Opportunity | Changes results? |
|---|---|---|---|---|
| 1 | gzip is inflated inside the reader's critical section | `mix`, 8 threads: 21.1 s `.gz` vs 8.5 s plain; `w900`: 23.0 vs 15.0 s. `zcat` of one input file takes 2.0–2.5 s | Inflate outside the lock (a reader thread per mate, or libdeflate/ISA-L), larger batches, a block reader instead of gzstream | no |
| 2 | Syncmer test recomputed from scratch for every window | 63k instructions per pair, ~260 per 31-mer window (~180 of them in the syncmer test), independent of the data: 40% of the loop on `mix`, 16–20% on `w900`/`s64` | Precompute the read's s-mers once per strand and slide the minimum; est. 2–3× on this part | no (must give the same syncmers; unit-testable) |
| 3 | WFA dominates marker-rich reads; alignments of hopeless anchors | 53% of the loop on `w900`, 73k instructions per call, 1.46 calls per read, ~39% of calls fail; on `mix` about 70% of the WFA work is for background reads | Ungapped fast path (≤1 mismatch is provably optimal), a pre-filter for anchors that cannot reach 90% identity, a cheaper WF-adaptive schedule | fast path: no; the rest: to validate |
| 4 | The 3 GB key map is accessed at random through 4 KB pages | `GLIBC_TUNABLES=glibc.malloc.hugetlb=1` (3.5 GB in huge pages): seeding 5.9→3.7 s and 4.4→3.1 s | `madvise(MADV_HUGEPAGE)` in `AllocateKeymap`/`AllocateValues`; document the tunable meanwhile | no |
| 5 | Fixed start-up costs | Key map: 3.1 GB RSS and 3.7 G instructions to decode for any database (2.2–2.7 s at 1 thread on `db64`). Model: 4.6 G instructions, 3.35 G of them in 229k exceptions | Two cPMML patches (below): −38% instructions for a 1000-pair run, profiles identical | no |
| 6 | Profiling runs one thread per sample | `w900`: 3.8 s after 19 s of 8-thread alignment; SAM text parsing is 72% of its instructions | Profile sample *i* while aligning sample *i+1*, or hand alignments over in memory | no |
| 7 | The built-in stage timers mislead | The k-mer extraction timer is measured but never printed; every interval is floored to whole µs; `<sample>_runtime.tsv` is in whole seconds and "Anchor recovery" is always 0 | Report k-mer extraction; accumulate in ns | no |
| 8 | Waiting threads spin | `mix`, 8 threads: 71.9 s user CPU with `.gz`, 35.3 s plain | Follows from 1; `OMP_WAIT_POLICY=passive` may reduce the burn on shared nodes (not tested) | no |

## Where the instructions go

callgrind, 50k pairs of each set, 1 thread (`scripts/callgrind.sh`, `scripts/cg_stage.sh`). The
alignment loop is `RunPairedEnd`; percentages are of that loop.

| Stage (instructions per pair) | `s64` / `db64` | `w900` / `db900` | `mix` / `db900` |
|---|---:|---:|---:|
| **alignment loop** | **310k** | **405k** | **156k** |
| FASTQ reader (gzip inflate 70% of it) | 26k (8.5%) | 26k (6.5%) | 26k (16.6%) |
| syncmer extraction | 63k (20.4%) | 63k (15.6%) | 63k (40.4%) |
| seeds and anchors | 32k (10.4%) | 40k (9.9%) | 18k (11.4%) |
| alignment handler | 168k (54.1%) | 255k (62.9%) | 47k (30.1%) |
| · of which WFA | 113k (36.4%) | 214k (52.9%) | 35k (22.1%) |
| · of which WF-adaptive cut-offs | 20k (6.5%) | 41k (10.2%) | 7k (4.3%) |
| SAM output | 16k (5.0%) | 16k (3.9%) | 1k (0.5%) |
| profiling, after alignment | 32k | 35k | 4k |
| WFA calls per read | 1.30 | 1.46 | 0.16 |
| index lookups (`Seedmap::Get`) per read | 22.9 | 22.9 | 22.9 |
| lookups turned into seeds per read | 9.3 | 10.3 | 2.9 |

Fixed costs per run, in instructions: index load 4.0 G (`db64`) and 11.0 G (`db900`), that is
~3.7 G for the key map plus ~280 per index value (zstd and column decoding); model parsing 4.6 G;
genomes, reference map, unique k-mers and SAM header 0.14 G (`db64`) to 2.4 G (`db900`).

In `mix`, its 2,500 database pairs would cost ~0.6 G in the alignment handler at `w900`'s rate; it
spends 2.35 G. Of its 16,012 WFA calls, about 7,300 are for the database reads (1.46 per read) and
8,700 for the random background reads; those take ~140k instructions each, 1.2 of the 1.7 G spent
in WFA, because they run to the score limit and fail.

Output SAM records of `w900` (`scripts/sam_stats.sh`, 894,815 records): 27% exact, 17% one
mismatch, 15% two, 29% three to five, 8% more, 4% with indels.

## Wall clock

1M pairs, `db900`, alignment stage only (`--no_profile`), "Aligning reads" time
(`scripts/scaling.sh`). The reader time is the mean per thread, including waiting for the lock.

| | `mix` | `w900` |
|---|---:|---:|
| 1 thread, `.fq.gz` | 35.8 s (reader 5.8 s) | 61.0 s (reader 5.0 s) |
| 8 threads, `.fq.gz` | 21.1 s (reader 11.7 s), 1.7× | 23.0 s (reader 7.0 s), 2.7× |
| 8 threads, plain `.fq` | 8.5 s (reader 1.7 s), 4.2× | 15.0 s (reader 1.7 s), 4.1× |
| user CPU, 8 threads `.gz` / plain | 71.9 / 35.3 s | 91.0 / 77.3 s |

Load average 11–14 during these runs, so the 8-thread speedups are lower than on a dedicated
node. The bound from the reader holds anywhere: inflating one pair takes ~4.5 µs on one core
(`zcat` alone: 2.0 s for R1, 2.3–2.5 s for R2), so gzipped input cannot go faster than ~200k
pairs/s, and 16- or 32-thread cluster jobs will not scale beyond it. Decompressing in separate
`zcat` processes (`-1 <(zcat R1)`) did not help: gzstream reads through a 303-byte buffer, and a
pipe costs a system call per buffer (reader time 17 s).

Full run of `w900` at 8 threads (`--verbose`): load index 1.2 s, alignment 19.2 s, profiling
4.7 s (3.8 s for the one sample, single-threaded), total 27 s, peak RSS 3.8 GB; the SAM is 350 MB.

### Huge pages

`mix` plain, 1 thread, `--no_profile`, alternated (`scripts/thp_test.sh`); THP is in `madvise`
mode on this system, so only the tunable puts the index in huge pages (AnonHugePages 0 vs 3.5 GB).

| Run | Seeding off → on | Alignment stage off → on |
|---|---:|---:|
| 1 | 5.94 → 3.67 s | 29.1 → 23.9 s |
| 2 | 4.39 → 3.08 s | 21.6 → 20.2 s |

Seeding is 30–38% faster; the load changed during the runs, so the effect on the whole stage
(6–18%) is less certain. The key map is the same size for every database, and a larger value
array adds TLB misses, so the gain on GTDB is likely at least as large.

## Memory

50k pairs of `s64`, `db64`, 8 threads (`scripts/mem_trace.sh`): 95 MB after loading the genomes, 3.14 GB during
the index load, 3.18 GB during alignment, 0.61 GB during profiling (the index is freed first).
With `db900`: peak 3.5 GB (1 thread) to 3.8 GB (8 threads). The key map is direct-indexed by the
15-mer core (4^15 slots), so it takes ~3 GB for any database; the values add 8 bytes each.

## Start-up: model parsing

cPMML throws and catches three exceptions per tree node (65k nodes in 256 trees): `Node` reads
`recordCount` with `std::stod("null")` (our PMML export does not write it), and `InternalScore`
converts every node's score with `stod` (internal nodes have none, leaves hold "TRUE"/"FALSE") and
throws again from the handler. [`patches/`](patches/) returns the same values without exceptions. 1000 pairs,
`db64`, 1 thread, three alternated runs each (`scripts/model_test.sh`):

| | as is | patched |
|---|---:|---:|
| instructions | 9.28 G | 5.73 G |
| in `__cxa_throw` | 3.35 G | 0.44 G |
| user CPU (mean) | 1.66 s | 1.40 s |

The profiles are byte-identical. The remaining 0.44 G are 24,674 `out_of_range` exceptions from
`MiningField`, `OutlierTreatmentMethod` and `OpType` parsing. The run takes 3.0–3.5 s, of which
2.2–2.7 s load the key map.

## Opportunities, in order

1. **Read input without serialising on gzip.** Inflate R1 and R2 in their own threads into blocks
   that workers parse, or use a faster inflater (libdeflate, ISA-L) with large reads instead of
   gzstream's `getline` into a `stringstream` and back. Gain on this machine up to the plain-input
   times: 2.5× for `mix` and 1.5× for `w900` at 8 threads, more at higher thread counts.
2. **Huge pages for the index:** `madvise(ptr, bytes, MADV_HUGEPAGE)` after the `calloc`/`malloc`
   in `Seedmap::AllocateKeymap` and `AllocateValues`. Until then, `GLIBC_TUNABLES=glibc.malloc.hugetlb=1`
   does it for users (glibc ≥ 2.35).
3. **Syncmer extraction** (`SimpleKmerHandler::operator()` with `ClosedSyncmer`): the 9 s-mer
   minimum of the canonical 15-mer core is recomputed for each of the ~120 windows of a read.
   Precomputing the read's 7-mer values once per strand and sliding the minimum should cut the
   ~180 instructions per window to well under 100. A unit test can check it returns exactly the
   current syncmers, which the index was built with.
4. **Fewer and cheaper WFA calls.**
   - An anchor chain on one diagonal that spans the read with at most one mismatch is optimal as
     an ungapped alignment: score 4 against at least 8 for any gap (mismatch 4, gap 6+2). 44% of
     the output records are such; they are the cheapest WFA calls, so the saving is well below
     that share (not measured). Reads that overhang a gene end still need WFA.
   - Background reads and relatives' anchors run WFA to the score limit and fail. A cheap bound
     before WFA (e.g. mismatches on the anchor diagonal, anchor coverage) could skip most of them;
     it must be checked against sensitivity, and the tie rule in `SimpleAlignmentHandler::operator()`
     (all anchors tied with the last of the top `-c`) should be watched on GTDB, where relatives tie.
   - WFA runs with its default WF-adaptive heuristic, cutting off after every step (10% of the
     loop on `w900`). `--x_drop` is parsed and stored in `WFA2Wrapper2` but never passed to WFA.
5. **Apply the cPMML patches** (−3.5 G instructions per run, same results) and, for small or custom
   databases, consider a key map that does not cost 3 GB and 2–3 s regardless of size.
6. **Profiling in parallel with alignment** for the next sample, or from alignments in memory
   instead of re-parsing the SAM (72% of profiling).
7. **Per-read allocations** (malloc and free are ~4.5% of the loop): `KmerUtils::ReverseComplement`
   takes its input by value and appends character by character; `SimpleAlignmentHandler::operator()`
   copies the read and computes the reverse complement the anchor finder already has;
   `AlignAnchor` takes `rev` by value and copies the reference window with `substr` for every call.
8. **Timers**: print "Retrieve k-mers", accumulate in nanoseconds, fix `_runtime.tsv`.

Smaller things seen on the way: `SimpleAlignmentHandler::ExtendSeed(ChainLink&, ...)` loops while
`qpos < qpos + s.length`, which is always true, and reads past the read; it is only reached from
`ExtendAllAnchors`, which nothing calls, so it is dead code to remove. `-t` defaults to 1.

## Gaps

- **GTDB scale was not profiled**, and it will change the balance. `db900` has 0.025 values per
  slot of the 4^15-slot key map. If the r226 database holds ~120 marker genes of ~1 kb for each
  of ~143k species (~15 Gbp), it has about 70 times `db900`'s values: nearly every lookup, also of off-target reads,
  returns a block whose flex keys are scanned, reads get more seeds and anchors, and conserved
  genes give ties among relatives, hence more WFA calls. Index load then scales with the values
  (~280 instructions each here). Run `scripts/` on the cluster with the real database: callgrind on
  50k pairs of a realistic mix, and the 1/8/32-thread runs with `.gz` and plain input.
- A quiet machine: absolute times, thread-scaling curves and the whole-stage THP gain need a
  dedicated node; this CPU also mixes fast and slow cores.
- No hardware counters (no `perf`): cache and TLB misses are inferred from the huge-page test.
- Real reads (quality trimming, adapters, host DNA, low complexity), other read lengths, and
  multi-sample maps (index loaded once, samples profiled in parallel) were not run; nor the strain
  MSA stage (it needs two or more samples) or the SAM compression step of map runs.
- The gains for syncmers and WFA are estimates from instruction counts, not prototypes.

## Follow-up: cPMML upstream, and the patches on `7c6b2f8`

Same day, after `performance` was fast-forwarded to `audit-fixes` at `7c6b2f8`.

**No newer cPMML.** [AmadeusITGroup/cPMML](https://github.com/AmadeusITGroup/cPMML) has no releases
or tags. Its `master` still ends at `2cd19f9` (2021-01-29), the commit protal imported (version
0.1); its only other branches are dependabot updates of the documentation's Python requirements.
The vendored copy differs from upstream only by protal's `Model::from_string` and its CMake
changes, so there is nothing to update to.

**The patches still apply, and still help.** The models now come in two kinds. The shipped
`scripts/random_forest.xml` (SoftwareAG PMML Generator, 65,280 nodes; also in databases built with
it) has no `recordCount`, so it needs both patches. Models written by `scripts/model_pmml.py`
(`random_forest_cmdline.py`, a few hundred to ~7,000 nodes) carry `recordCount` and a class label as
`score` on every node: the `recordCount` patch does nothing for them, and the score patch removes
their throws. 1000 pairs, `db64`, 1 thread, `--model` pointing at each model
(`patch_check` in the session scratchpad; the same measurement as `scripts/model_test.sh`):

| Model | Model parse, instructions | Exceptions |
|---|---:|---:|
| shipped, as is → patched | 4.71 G → 1.16 G | 229,000 → 24,674 |
| `~/tune/D/trained_model.xml` (6,818 nodes), as is → patched | 746 M → 249 M | 35,530 → 6,437 |

Profiles and their companion files are byte-identical with and without the patches for both
models. Both patches are now applied to `lib/cPMML` on `performance`. The exceptions left come
from `MiningField`, `OutlierTreatmentMethod` and `OpType` looking up absent attributes with
`.at()`.

## Follow-up: the reader, huge pages and the timers, implemented

Same day, on `performance` (`7c6b2f8` plus these changes and the cPMML patches), same machine and
data. The machine was quieter than in the first round (load 2–8, most of it these runs), but the
alternated runs still vary by up to 2×; the reader benchmark is the clean measurement.

**What changed**

- `src/IO/ThreadedGzStream.h`: an input stream that inflates the file with zlib in a thread of its
  own, into four 1 MB blocks ahead of the reader. The FASTQ inputs of paired, single-end and long
  reads use it instead of `igzstream` (`RunProtal.h`), so the reader lock only copies inflated
  bytes, and R1 and R2 inflate in parallel. Truncated or corrupt files are still reported
  (`read_failed()`, as gzstream). No new dependency: libdeflate or ISA-L would inflate faster
  still, but neither is on this machine or in the conda recipe.
- `Seedmap::AllocateKeymap`/`AllocateValues`: `madvise(MADV_HUGEPAGE)` on blocks of 64 MB or more,
  before they are touched. Without THP support it does nothing.
- `Benchmark`: sums nanoseconds on `steady_clock` (each interval was floored to whole
  microseconds on `high_resolution_clock`); printing no longer divides the stored sum by the
  threads, so `_runtime.tsv` no longer depended on `--verbose`. "Retrieve k-mers" is timed (also
  for single-end reads), joined, printed and written; "Anchor recovery" is joined (it was always 0);
  `RecoverAnchors` stops its timer on the early return; an unused timer per alignment in
  `AlignAnchor` is gone. `_runtime.tsv` has a header, `stage seconds threads seconds_per_thread`,
  with seconds to the microsecond instead of whole seconds.

**Reader alone**, 1M pairs of `mix` (`.gz`), threads that only take pairs
(`scripts/bench_reader.sh`), three runs each:

| | 1 thread | 8 threads |
|---|---:|---:|
| `igzstream` | 2.4–2.7 s (~400k pairs/s) | 3.1–3.5 s (~300k pairs/s) |
| `ThreadedGzIstream` | 0.82–0.95 s (~1.15M pairs/s) | 1.15–1.26 s (~830k pairs/s) |

The ceiling rises 2.7–2.9×. With 8 threads the new reader is slower than with one: the lock
changes hands every 32 records. Batches of 128 or 512 records made no consistent difference in
this benchmark (runs of one setting varied 2×), so the batch size is unchanged.

**Alignment stage**, 1M pairs, `db900`, `--no_profile`, alternated (`scripts/ab_alignment.sh`);
median of the runs, "old" is `7c6b2f8` as committed, "old + THP" the same with
`GLIBC_TUNABLES=glibc.malloc.hugetlb=1`:

| | old | old + THP | new |
|---|---:|---:|---:|
| `mix` `.gz`, 8 threads (5 runs) | 6.64 s | 6.27 s | 4.50 s |
| `w900` `.gz`, 8 threads (4 runs) | 10.2 s | 9.4 s | 8.1 s |
| `mix` `.gz`, 4 threads (3 runs) | 5.90 s | 5.75 s | 4.52 s |
| `mix` plain, 8 threads (3 runs) | 4.28 s | 3.27 s | 3.07 s |

With gzipped input the new build is 20–32% faster at 8 threads; most of it is the reader, as the
plain-input row (no inflating at all) shows. The gain depends on load: in the first round, on a
machine with load 11–14, the same old run took 21 s, because a thread preempted while holding the
reader lock stalls all others. The new build gets huge pages without the tunable (AnonHugePages
3.35 of 3.58 GB resident). Single-threaded runs varied too much (12.7–19.9 s) to show the huge-page
gain on the whole stage; the seeding gain measured before (30–38%) is the better number.

The old stage timers lost most of short intervals: in the same conditions the old build printed
0.47 s for seeding and the new one 0.81 s, and k-mer extraction, now printed, is the largest
stage on `mix` (0.73 of 1.82 s per thread at 4 threads). Stage times from the first round are
low for the short stages.

**Same results.** 200k pairs of `w900` and of `mix`, old against new: at 1 thread the SAM and
every output file are byte-identical; at 8 threads (and single-end at 4) the profiles, logs and
statistics are identical and the SAM holds the same records in another order. Unit tests 134/134
(8 new in `tests/test_ReaderAndTimers.cpp`), also under ASan/UBSan; the new stream's tests are
clean under TSan (5 repeats; the OpenMP test is left out, libgomp is not instrumented; TSan needs
`setarch -R` on this kernel); e2e 85/85; `examples/mini_db` passes.

**Also from the docs:** `docs/running.md` puts the full r226 database at about 59 GB in memory.
Even if half of that were reference sequence, the index would hold over 100 times `db900`'s
26.9M values (8 bytes each), more than the ~70 times estimated under [Gaps](#gaps); seeding and
huge pages matter more there than this small database shows.

## Follow-up: syncmers, where WFA spends its time, and --x_drop

Same day, on `performance` after `bdc5790`, same machine and data.

**Syncmers, prototype** (`scripts/bench_syncmer.cpp`, 500k reads of `mix`, one thread). It
computes each read's 7-mers once per strand instead of 9 times per 31-mer window and finds the first
minimum of a window branch-free (7-mer value packed with its index, both strands evaluated,
the canonical one selected). The k-mer lists are identical to `SimpleKmerHandler<ClosedSyncmer>`'s
for all 1.5M sequences tried (the reads, their reverse complements, and the reads with Ns):

| | ns per read (150 bp) |
|---|---:|
| `SimpleKmerHandler` (current) | 2,530–2,740 |
| 7-mers once, branching minimum | 1,890–2,060 (1.2–1.4×) |
| 7-mers once, branch-free | 1,150–1,360 (2.1–2.4×) |

The branching version gains little: the minimum and the canonical strand are data-dependent
branches. Vectorising 8 windows at a time (AVX2) would be the next step; not tried.

**Where WFA spends its time.** `scripts/wfa_calls.patch` logs every WFA call (outputs unchanged):
50k pairs, one thread, `scripts/wfa_calls.sh`.

| outcome | `w900` calls | ns/call | WFA time | `mix` calls | WFA time |
|---|---:|---:|---:|---:|---:|
| failed (score limit reached) | 46.7% | 5,534 | 66.3% | 75.3% | 86.4% |
| aligned, >5 mismatches | 19.3% | 4,464 | 22.1% | 9.0% | 9.0% |
| aligned, 2–5 mismatches | 15.9% | 1,286 | 5.2% | 7.4% | 2.1% |
| aligned, internal indel | 4.0% | 4,121 | 4.2% | 1.9% | 1.7% |
| aligned, 1 mismatch | 8.1% | 670 | 1.4% | 3.8% | 0.6% |
| aligned, exact | 6.0% | 469 | 0.7% | 2.7% | 0.3% |

68% (`w900`) and 85% (`mix`) of the anchors are one exact link, mostly 20–49 bp; 98% of the
failed calls have 16 or more mismatches along their anchor's diagonal (relatives' genes, and in
`mix` random reads). An ungapped fast path for reads with at most one mismatch would save about
2% of WFA time: those calls are already cheap. The time goes into proving that anchors fail.

**Anchored extension, prototype.** For single-link anchors, WFA extended from the two ends of the
link instead of aligning the whole read in a window with dovetails: right of the link, then left
(reversed), each with its far end free, the score limit shared between them, so a failing side
stops the other. Next to the current alignment, for the same calls:

| | `w900` | `mix` |
|---|---:|---:|
| single-link calls | 99,681 | 13,612 |
| WFA time, current → anchored | 0.452 → 0.296 s (1.52×) | 0.067 → 0.045 s (1.50×) |
| aligned by both: same score | 37,908 of 37,931 (99.94%) | 1,909 of 1,909 |
| anchored better / worse | 23 / 0 | 0 / 0 |
| aligned only by the current / only anchored | 0 / 7 | 0 / 2 |

A real version also needs the CIGAR (left reversed + the link + right), multi-link anchors (the
gaps between links aligned end to end) and the SAM position; its outputs would have to be checked
against today's, as the scores differ in 0.06% of alignments (always better).

**`--x_drop`** was parsed but never reached WFA2. It is now WFA2's X-drop on top of the
wf-adaptive heuristic that stays on (0 turns it off). WFA2 alone does not suit protal's scoring
(`scripts/xdrop_wfa2.cpp`: with a match scoring 0, X-drop alone aborted a 2 kb alignment at 2%
divergence with 50 and an 8 kb one at 1% with 200). Added to wf-adaptive: 200k pairs of `w900` and
`mix`, `-x 0`, the default 1000 and `-x 200` give byte-identical SAMs and profiles; `-x 50` reports
the same species with abundances up to 0.007% (`w900`) and 0.36% (`mix`) apart, and changes 10,093
SAM records of `w900` (45 fewer alignments, and MAPQs). The default costs 0.12% more instructions (callgrind, 20k pairs). With `-x 50`, WFA2
reported one alignment complete whose operations held 151 read bases for a 150 bp read, and protal
stopped (exit 90); the wrapper now counts such an alignment as dropped
(`tests/test_WFA2Wrapper.cpp` has that read). X-drop does not make protal faster: it only prunes
at steps where wf-adaptive did not.

Calling WFA2 from unit tests showed that UBSan reports its unaligned 8-byte loads and left shifts
of negative values, which stop CI's sanitizer job (`halt_on_error=1`); `lib/wfa2-lib.cmake` now
builds WFA2 without UBSan (ASan stays on).

## Follow-up: the syncmer scan, implemented

2026-09-30, on `performance` after `ff71972`. `SimpleKmerHandler<ClosedSyncmer>` now scans whole
sequences as the branch-free prototype does (`ScanClosedSyncmers`); the window-by-window loop stays
as `WindowByWindow` (the definition, and the path of any other minimizer). Reads and the database
build use the same code.

- The same k-mers: `tests/test_Syncmers.cpp` compares the scan with the loop and with a brute force
  from `ClosedSyncmer`'s definition, for both s-mer masks (index formats 1 and 2), on random,
  low-complexity and short sequences, with Ns, lower case and other symbols. The mini database built
  by the old and the new binary is byte-identical (`database.protal`), and every output of 200k
  pairs of `w900` and of `mix` (1 thread) is identical.
- 200k pairs of `mix`, 1 thread, alternated three times: "Retrieve k-mers" 1.08–1.21 s before,
  0.53–0.78 s after; the alignment stage 2.44–2.73 s before, 1.87–2.76 s after. Callgrind, 20k
  pairs: the handler 1.26 G → 0.80 G instructions, the alignment loop 2.75 G → 2.29 G (−17%).
- Unit tests 139/139, also under ASan/UBSan; e2e 85/85.

## Reproducing

The scripts are in [`scripts/`](scripts/) (settings in `env.sh`: `PERF_DIR`, default
`~/protal-perf`, on a Linux file system). Generated data stays in `PERF_DIR`, not in git.

```bash
export PERF_DIR=~/protal-perf
bash scripts/build.sh
bash scripts/prep_db.sh db64 ~/audit4/gtdb64
bash scripts/prep_db.sh db900 ~/tune/world          # the 900-species world, see below
bash scripts/prep_reads.sh s64 ~/audit4/gtdb64/simulation/genomes.tsv 20
bash scripts/prep_reads.sh w900 ~/tune/world/simulation/genomes.tsv 60
bash scripts/prep_reads.sh --mix mix w900
bash scripts/callgrind.sh w900 $PERF_DIR/db900 $PERF_DIR/reads/w900/reads 50000
bash scripts/scaling.sh $PERF_DIR/db900 $PERF_DIR/reads/mix mix 8
bash scripts/thp_test.sh $PERF_DIR/db900 $PERF_DIR/reads/mix_plain mix 1 2
bash scripts/mem_trace.sh mem_s64 $PERF_DIR/db64 $PERF_DIR/reads/s64/reads 8
bash scripts/model_test.sh $PERF_DIR/db64 $PERF_DIR/reads/s64/reads
# follow-ups: the reader alone, and two builds' alignment stage alternated
bash scripts/bench_reader.sh $PERF_DIR/reads/mix/mix_R1.fq.gz $PERF_DIR/reads/mix/mix_R2.fq.gz 3
bash scripts/ab_alignment.sh OLD/protal_avx2 NEW/protal_avx2 $PERF_DIR/db900 $PERF_DIR/reads/mix 8 5
```

The worlds: `~/audit4/gtdb64` is a 64-species `simulate_gtdb_release.py` release with 3 genomes
per species; `~/tune/world` was made with `gtdb_like_lineages.py --species 900` (then untracked in
the main checkout) and `simulate_gtdb_release.py --lineages`, with strain and species divergence
ranges.
