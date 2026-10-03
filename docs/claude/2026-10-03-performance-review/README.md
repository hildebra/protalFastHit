# Performance review, 2026-10-03: where a run's time goes now, and what would make it faster

- **Date**: 2026-10-03.
- **Code**: branch `audit-fixes` at `27423c6` (the merge of `denoise`), built from `git archive` in WSL
  (`scripts/build_ref.sh`, Release, one binary). Nothing in the tree was changed; the one experiment build
  (`scripts/tim_build.sh`) only prints the wall time of the profiler's passes.
- **Machine**: WSL2 Ubuntu 24.04 on an Intel Core Ultra 7 258V (4 fast and 4 low-power cores, 6 vCPUs), GCC 13.3,
  23 GB. The machine was idle (load 0.0 at the start; no other session ran). Wall times still swing with the
  core a thread lands on, so take few-percent differences as noise. `perf` is not available here; instruction
  counts are callgrind's, one thread.
- **Data**: the v0.7.3 benchmark world (`~/bench071/V073`, 765 species, 2,295 genomes, database with gene
  neighbours and congeners). Reads: `rl150_p500000_s_1` (500k pairs), `rl150_p5000000_s_1` (5M pairs, `samples_deep`),
  `ont_b90000000_s_1` and `pb_b90000000_s_1` (90 Mb of long reads), the first 100k pairs of the 500k sample for callgrind,
  and samples 1–4 of the 500k point as a four-sample map. `--no_qcmsa` everywhere (qcmsa is python, outside protal).
- **Earlier rounds**: [round 1](../2026-09-29-performance-profiling/README.md), [round 2](../2026-09-30-performance-round2/README.md),
  [round 3](../2026-10-02-performance-round3/README.md), [prefetch](../2026-10-02-prefetch/README.md),
  [multithreading audit](../2026-10-01-multithreading-audit/README.md). This review repeats their measurements on
  the code as it is now, with the features that came in since (gene neighbours, mate guidance, phasing, the
  sample-context features of `be35d15`), and ranks what is left.
- **Question**: where does protal spend its time today, and where are the gains?

## Summary

1. **The profiling stage is now 30–40% of a short-read run.** Round 3 measured it at ~10%. On 5M pairs at six
   threads the whole run takes 34.0 s: aligning 20.2 s, profiling 12.6 s. On 500k pairs at one thread: 17.0 s,
   of which profiling 5.1 s. In a four-sample cohort (4 × 500k pairs, six threads): 15.2 s, of which profiling 4.8 s
   and strain MSAs 0.6 s.
2. **What grew is the sample-context pass of `be35d15`** (`MicrobialProfile::ApplySampleContext`), and above all its
   abundance-weighted read assignment (`context::AbundanceWeightedShares`, the "read EM"):
   - instructions, 500k pairs: profiling 36.3 G in all; the EM 14.0 G (38%), the congener distances (sketching the
     references' genes) 4.7 G (13%), reading the SAM 9.7 G (27%), writing the profile 3.1 G (9%), parsing the model
     1.7 G (5%);
   - wall time, 5M pairs at six threads: the EM 7.5 s of the 12.0 s that profiling the sample takes, the distances
     1.1 s, parsing and adding the records 3.1 s, everything else 0.6 s;
   - **the EM runs on one thread**, so more threads do not help: the 5M SAM profiles in 12.2 s at 6 threads, 12.4 s at 8,
     11.5 s at 12. At one thread it takes 29–31 s, so the stage scales 2.4× on six threads where the alignment scales 4.7×.
   - The EM's cost is in how it is written, not in what it computes: 441k ambiguity classes (5M pairs) are visited up
     to 200 times, each visit doing a dozen hash-map lookups and two or three `pow` calls, and each iteration copies
     a map. Resolving the classes to indices once and taking `pow` out of the loop leaves a multiply-add per
     alternative per iteration: the same arithmetic in the same order, so the same output, for roughly 20× fewer
     instructions. **This is the one change that pays on every short-read run: about −20% of a deep sample's
     wall time at six threads, −50% of the profiling stage.**
3. **The alignment loop is where the earlier rounds left it.** 100k pairs, one thread: 17.1 G instructions in
   `RunPairedEnd` (round 3: 17.0 G on the v0.7.1 world); WFA2 6.7 G, almost all in the anchored aligner's flanks; the
   anchor finder 2.8 G; the output handler 1.7 G (zstd 0.8 G); syncmers 1.2 G. New since round 3 and visible: the
   gene-neighbour table lookups of mate guidance (0.43 G, 2.5% of the loop), `IsAlignmentValid` re-walking every
   alignment's CIGAR against the sequences as a safety net (0.46 G, 2.7%). Per pair the loop is 171k instructions.
   The stage timers at one thread: seeding 1.4 s of 9.8 s (the k-mer lookups, memory-bound), sorting seeds 0.7 s,
   alignment 3.3 s, joining pairs 1.0 s, output 0.8 s.
4. **Fixed costs per run**: loading the index 7.3 G instructions (1.6 s at one thread, 0.3–0.6 s at six), the genomes
   1.1 G, the model 1.7 G (0.4 s). A 500k-pair run at six threads takes 5.0 s, so these are a fifth of it; they do not
   matter for deep samples or cohorts.
5. **Long reads**: 90 Mb of Nanopore reads run in 6.7 s at six threads (aligning 4.7 s, of which the alignment handler
   is 4.3 s per thread; profiling 0.9 s, of which the congener distances 0.56 s); PacBio 5.0 s. The long-read
   alignment is what round 3 changed; its open item (the window's free ends from the last link) stands.
6. **A measurement artefact to know about**: the per-read stage timers time one call in 61. The output handler's
   calls are mostly tiny and now and then a buffer flush (compress and write a megabyte); whichever a sample hits
   decides the estimate. In one six-thread run it reported 1.25 s per thread, in three others 0.15–0.18 s; the
   latter is right (its instructions are 10% of the loop's, at one thread it is 0.8 s of 9.4 s). The alignment
   handler's and anchor finder's calls are uniform enough for the sampling.

## What I would do, in order

| # | change | where | expected | outputs | effort |
|---|---|---|---|---|---|
| 1 | the read EM on indices: resolve each class once to (own index, alternative indices, `pow(ratio, edits)` factors), dense weight vectors, no per-iteration map copy | `SampleContext.h`, `AbundanceWeightedShares` | EM 14.0 G → under 1 G per 500k pairs; 5M at 6 threads 34 → ~27 s; 500k at 1 thread 17 → ~15 s | identical (same operations, same order) | small |
| 2 | then: the EM's stopping rule. `kEmTolerance` is 1e-10 on the largest relative change, `kEmIterations` 200, and both samples here run 197–200 sweeps (the 500k one hits the cap unconverged); a tolerance in line with what the features resolve, or a parallel sweep with a fixed reduction order | same | iterations ÷ 2–5 | changes the `em_*` features in their last digits, so the profiles may change | small, needs a benchmark check |
| 3 | congener distances: `GeneSketch` keeps 64 hashes but `gene_conservation::Kmers` first builds and sorts every k-mer of the gene (3.1 of the 4.6 G); take the 64 smallest hashes with a bounded selection instead; and make the sketches of all the sample's taxa in one parallel pass before the pairs, rather than on first use under a mutex | `SampleContext.h` `GeneSketch`, `CongenerDistances` | 4.6 G → ~1.5 G; the 1.1 s (5M, six threads) and 0.5–0.6 s (long reads) → ~0.3 s | identical (the same 64 hashes) | small |
| 4 | cohorts: profile sample *i* while sample *i*+1 aligns (round 1's item 6). With 1–3 done the serial tail is short, so one or two threads suffice | `RunProtal.h` | hides most of the profiling stage in map runs (4.8 of 15.2 s here) | identical | medium (threads, memory per sample) |
| 5 | alignment loop, identical outputs, in the order they pay: (a) `IsAlignmentValid` on every alignment → one in 64, or debug builds only (2.7% of the loop); (b) memoise `gene_neighbours::Table::Partners/Assess` per thread by (taxon, gene, end) (2.5%); (c) the seeds sorted on one packed 64-bit key, radix or `std::sort` (the sort is 7% of the one-thread wall time for 1.6% of the instructions: memory- and branch-bound; a key sort halves it); (d) one pass over the CIGAR instead of `GetInstructionCountsAndCompress` + `NextCompressedCigar` + `AlignmentEdits` (~4%); (e) a small WFA2 aligner of their own for short-read flanks (`wavefront_slab_reap_repurpose` + `slab_allocate` 3.5%, round 3's item 4) | `AlignmentStrategy.h`, `AcrossGenes.h`/`MateGuidance.h`, `ChainAnchorFinder.h`, `AnchoredAlignment.h` | −8…−12% of the alignment loop together | identical | small each |
| 6 | PGO in the release build | CMake | −3…−8% (round 2) | identical | medium (training inputs in the build) |
| 7 | GTDB scale: nothing above was measured on the real database. Run `scripts/measure_performance.sh` on a cluster node (it records `perf stat` where allowed) on one deep short-read sample, one long-read sample and a four-sample map | cluster | tells whether the lookups (prefetch), the index load (value cells) or the sketches (more species) change the picture | — | an hour of cluster time |

Not worth it now: more decompression threads for the SAM (tested: 8 and 12 threads no faster than 6, the
stage is bound by the serial EM); the output handler (10% of the loop's instructions, half of it zstd at level 3,
and it scales with the threads); the index load for this database (0.3–0.6 s at six threads), unless GTDB-size
loads say otherwise; WFA2's AVX2 kernels (round 2: no wall gain).

## Measurements

### Whole runs (`results/whole_runs.txt`, `results/stage_times.txt`)

| run | threads | wall | load index | aligning | profiling | strains | max RSS |
|---|---|---|---|---|---|---|---|
| 500k pairs | 1 | 17.0 s | 1.64 s | 9.77 s | 5.14 s | – | 3.4 GB |
| 500k pairs | 6 | 5.0 s | 0.3 s | 2.22 s | 2.01 s | – | 3.7 GB |
| 5M pairs | 6 | 34.0 s | 0.30 s | 20.2 s | 12.6 s | – | 3.7 GB |
| Nanopore 90 Mb | 6 | 6.7 s | 0.62 s | 4.75 s | 0.90 s | – | 3.7 GB |
| PacBio 90 Mb | 6 | 5.0 s | 0.41 s | 3.42 s | 0.77 s | – | 3.6 GB |
| 4 × 500k pairs (map) | 6 | 15.2 s | 0.46 s | 8.88 s | 4.77 s | 0.56 s (53 species) | 3.7 GB |

The one-sample runs make no strain MSAs (a species needs two samples). User CPU: 500k at 1 thread 16.0 s, at 6
threads 19.3 s; 5M at 6 threads 148 s (so the 34 s wall is 4.4 cores busy on average, of 6).

Aligning 500k pairs alone (`--no_profile`): 1 thread 9.38 s, 3 threads 3.55 s, 6 threads 2.00 s (4.7×); plain
SAM instead of `.sam.zst` at six threads: 1.91 s (within noise).

### The profiling stage, pass by pass

Wall time printed by the experiment build (`scripts/tim_build.sh`; `results/profiling_split.txt`), `--profile_only` on the
runs' own SAMs:

| SAM | threads | parse + add records | EM (`AbundanceWeightedShares`) | congener distances | post-process SNPs, score, write | "Profile sample" | classes / pairs |
|---|---|---|---|---|---|---|---|
| 500k pairs | 1 | 0.84 s | 2.52 s | 1.31 s | 0.24 s | 4.75 s | 76,028 / 641 |
| 500k pairs | 6 | 0.35 s | 1.06 s | 0.45 s | 0.08 s | 1.88 s | 76,028 / 641 |
| 5M pairs | 1 | 8.24 s | 16.5 s | 3.26 s | 2.05 s | 28.8 s | 441,369 / 1,315 |
| 5M pairs | 6 | 3.13 s | 7.50 s | 1.06 s | 0.76 s | 12.0 s | 441,369 / 1,315 |
| Nanopore 90 Mb | 6 | 0.19 s | 0.03 s | 0.56 s | 0.11 s | 0.86 s | 2,936 / 451 |

- "classes" are the ambiguity classes the EM iterates over (distinct best-taxon / kept / alternatives-with-edits
  patterns, `RecordEvidenceCollector::NoteAmbiguity`); they grow with the reads. "pairs" are the congener pairs whose
  references are sketched and compared (`kMaxRelatives` = 4 congeners of more fragments per taxon).
- The EM is single-threaded (`ApplySampleContext(threads)` parallelises only the distances). Its one-thread times
  differ between the 500k and 5M runs by more than the class count (the 5M run at one thread also landed on slower
  cores for its whole 31 s; its user time equals its wall time).
- `FinishLink` is nothing; `PostProcessSNPs`, scoring and the profile writers are 0.3–2 s at one thread and
  parallel.
- Threads for the 5M SAM (`results/whole_runs.txt`, runs3): 6 → 12.2 s, 8 → 12.4 s, 12 → 11.5 s (8 and 12 add zstd
  decompression threads, `Profiler::DecompressThreads`; the machine has 6 vCPUs).

### Profiling instructions (callgrind, 500k pairs, `--profile_only`, one thread; `results/callgrind_profile_pe500k_*.txt`)

36.3 G in all, 418,466 SAM records.

| what | G | share | notes |
|---|---|---|---|
| `ApplySampleContext` | 18.6 | 51% | |
| – `AbundanceWeightedShares` | 14.0 | 38% | self 5.6 G; `pow` 5.1 G; `unordered_map` 2.6 G; `std::map` copies and increments |
| – `CongenerDistances::Compared` | 4.7 | 13% | `GeneSketch` 4.6 G, of which `gene_conservation::Kmers` 3.1 G (every k-mer collected, sorted, deduplicated) and the hash sort 1.1 G |
| `ReadSamGroups` / `CollectGroups` | 9.7 | 27% | `SamReader::Advance` 3.8 G, zstd 2.2 G, `getline` 2.4 G (inclusive, overlapping) |
| `Taxon::AddSam` | 3.1 | 8% | `VariantHandler::AddAlignment` 2.6 G |
| `PrepareMAPQ` | 1.9 | 5% | |
| `WriteSparseProfile` | 3.1 | 9% | `SequenceRangeHandler::CoveredPortion` 2.4 G, `CalculateCoverageVector2` 0.75 G (round 3's item 5) |
| model parse (`Node::Node`, cPMML) | 1.7 | 5% | fixed per run |

### Alignment instructions (callgrind, first 100k pairs, `--no_profile`, one thread; `results/callgrind_align_pe100k_*.txt`)

28.07 G in all.

| what | G | notes |
|---|---|---|
| `RunPairedEnd` (the loop) | 17.08 | 171k per pair |
| – `SimpleAlignmentHandler` | 8.79 | `AnchoredAligner::Align` 7.0 G; WFA2 (`wavefront_align`) 6.7 G, of which `Flank` 6.2 G; slab reset + allocate 0.59 G; `memset` 0.46 G (WFA2's `compute_process_ends`) |
| – `ChainAnchorFinder::operator()` | 2.79 | `GetFromLookup` 0.56 G, `ExtendAnchor` 0.44 G, `LookupKmers` 0.40 G, sort 0.16 G + comparator 0.12 G |
| – `ProtalPairedOutputHandler` | 1.68 | `SamOutput::Write` 0.82 G = `sam_zstd::Compress` 0.81 G; `ExtractSNPs` 0.20 G; `AlignmentEdits` 0.16 G |
| – syncmers (`ScanWindowsAvx2` + `Bases16`) | 1.18 | |
| – `IsAlignmentValid` | 0.46 | every alignment's CIGAR walked against read and gene (`AlignmentStrategy.h:580`) |
| – `AlignmentInfo::GetInstructionCountsAndCompress` | 0.39 | plus `NextCompressedCigar` 0.19 G |
| – `gene_neighbours::Table::Assess` | 0.43 | 419 M of it from `Table::Partners` (mate guidance's search on the neighbour genes), the rest from `PairAcrossNeighbours` |
| – `KmerUtils::ReverseComplementInto` | 0.25 | |
| – `mate_guidance::BestDiagonal` | 0.10 | |
| gzip inflate (reader thread) | 1.12 | |
| `Seedmap::Load` | 7.33 | `DecodeChunk` 4.92 G (self 4.43), zstd 3.4 G |
| `LoadAllGenomes` | 1.10 | `gene_neighbours::Table::Read` 0.20 G |

Stage timers, 500k pairs at one thread (`results/runtime_tsv.txt`): seeding 1.37 s, sorting seeds 0.70 s, pairing
0.38 s, extending anchors 0.65 s, sorting anchors 0.09 s, recovering anchors 0.19 s, reader 0.37 s, k-mers 0.58 s,
alignment handler 3.28 s, joining pairs 0.99 s, output 0.82 s; "Aligning reads" 9.77 s. Seeding is 14% of the time
for about 6% of the instructions (the lookups wait on memory; the prefetch of `87b857e` is in this build).

## The read EM, in detail

`context::AbundanceWeightedShares` (`SampleContext.h`):

- Input: `classes`, a `std::map<std::vector<uint32_t>, uint64_t>` from a key (own taxon, kept flag, then
  (alternative taxon, edits more) pairs) to the number of reads of that pattern; `counts` per taxon.
- Each iteration: `next = plain` (a map copy), then for every class a `posterior` that looks the own weight and
  every alternative's weight and ratio up in `unordered_map`s and computes `pow(ratio, edits)` for each alternative,
  then adds into `next` (hash lookups again), then a pass over `next` for the largest relative change.
- `ratio[taxid]` never changes between iterations, so `pow(ratio, edits)` per (class, alternative) is a constant;
  the own and alternative taxa of a class are constants too. Resolving a class once into indices into a dense
  weight vector plus its factors leaves, per iteration, one multiply per alternative, one division per class and
  the additions, in the same order as now, hence the same floating-point results.
- The map copy per iteration and the change pass over a map become vector operations.
- Estimated: from ~920 instructions per class and iteration (14.0 G / 76k classes / ~200 iterations) to 30–50.

**EM iterations** (`scripts/tim_iters.sh`, `results/profiling_split.txt`): the 500k sample runs all 200 iterations and
stops at the cap with the largest relative change still 4.7e-4 (a taxon of tiny weight keeps moving); the 5M
sample stops after 197 at 9.2e-11. So the pass costs its full 200 sweeps on every sample of any depth, which is
what item 2 is about; item 1 makes each sweep cheap regardless.

## Notes and gaps

- GTDB-scale behaviour remains unmeasured (as in every earlier round). What scales with the database: the lookups'
  memory stalls (the prefetch should pay more), the index load's value cells, the number of species with reads
  (sketches, EM classes), and the gene-neighbour table. `scripts/measure_performance.sh` exists for that
  (`def47d4`); its output is what a GTDB round should start from.
- The six-thread numbers come from a laptop whose six vCPUs are four fast and four slow cores shared by Windows.
  Scaling factors above (4.7× aligning, 2.4× profiling) are relative within this machine.
- Real reads (quality trimming, host DNA, adapters) were not run; neither was qcmsa.
- The strain stage ran only in the map run (0.56 s for 53 species of 4 shallow samples); it is not profiled here.
- No code was changed in the repository. The experiment build's prints (`scripts/tim_build.sh`, `tim_iters.sh`) are
  perl substitutions on a copy; they are not a patch to apply.

## How it was run

Paths are the session's WSL scratch (`~/mt-work/perf4`, `~/bench071`).

```bash
bash scripts/build_ref.sh        # HEAD's archive built in ~/mt-work/perf4/ref
bash scripts/runs1.sh            # whole runs: 500k (1 and 6 threads), 5M, Nanopore, PacBio; stage timers
bash scripts/cg_prof.sh          # callgrind of --profile_only on the 500k SAM (one thread)
bash scripts/cg_align.sh         # callgrind of aligning 100k pairs, --no_profile (one thread)
bash scripts/runs2.sh            # output-handler check (plain SAM, 1/3/6 threads), profiling of the 5M SAM at 6 and 1 threads
bash scripts/runs3.sh            # the 5M SAM at 8 and 12 threads
bash scripts/runs4.sh            # the four-sample map (strains on)
bash scripts/tim_build.sh        # the experiment build with wall-clock prints around the profiler's passes; its runs
bash scripts/tim_iters.sh        # the EM's iteration count
bash scripts/collect.sh          # results/
```

## Follow-up (2026-10-03): the read EM on indices, implemented

Item 1 of the list, at the user's request (`src/Profiling/SampleContext.h`, `AbundanceWeightedShares`; test
`SampleContext.SharesOnIndicesEqualTheReferenceToTheLastBit` in `tests/test_SampleContext.cpp`, which keeps the first
version as `ReferenceShares`).

- Every taxon of the counts or of a class gets an index once; a class is resolved once into its own index, its
  reads, its kept flag and its alternatives' indices with the constant factor `EditRatio^(edits more)` (0 for a
  taxon without records, whose weight is 0 throughout, as it had no map entry before). A sweep is a vector copy, a
  multiply per alternative, the same divisions and additions in the same order, and the change over the weight
  vector. The stopping rule is untouched (item 2).
- **Outputs**: `--profile_only` outputs of the reference build and the new one compared byte for byte
  (`scripts/em_build.sh`, `results/em_identical.txt`): identical on the 500k-pair SAM at 1 and 6 threads, the 5M-pair
  SAM at 6 threads, the Nanopore and PacBio 90 Mb SAMs at 6 threads (every file of the output folders, profiles,
  logs and per-species statistics). The test compares 40 random samples with ties, 0-4 alternatives at 0-5 edits,
  owners and alternatives without records: equal doubles everywhere.
- **Time** (`--profile_only`, "Profiling took", the same machine, one run each):

  | SAM | threads | before | after |
  |---|---|---|---|
  | 500k pairs | 1 | 5.7 s | 3.0 s |
  | 500k pairs | 6 | 2.0 s | 1.4 s |
  | 5M pairs | 6 | 13.2 s | 6.6 s |
  | Nanopore 90 Mb | 6 | 0.89 s | 0.85 s |
  | PacBio 90 Mb | 6 | 0.78 s | 0.74 s |

  On the 5M-pair run that is about 6.6 s of the 34 s whole run at six threads (−19%). The long-read samples have few
  classes (2,936) and gain nothing; their profiling time is the congener sketches (item 3).
- Instructions (callgrind, `--profile_only` on the 500k SAM, one thread, `scripts/cg_em.sh`): the EM 13.95 G → 2.97 G
  (−79%; 200 sweeps over 76k classes, about 190 instructions per class and sweep, most of them the four divisions
  per class and the posterior's small vector), the profiling stage 36.3 G → 25.3 G (−30%). What remains of the
  stage is the congener sketches (4.67 G, now 18% of it, item 3), reading the SAM (9.7 G), the EM (3.0 G) and
  writing the profile (3.1 G).
- Tests on the patched tree (`scripts/em_tests.sh`, `em_e2e_mini.sh`, `results/em_tests.txt`): 309 unit tests pass
  (`ctest`); the 131 end-to-end tests pass against a mini database. Against the V073 database the end-to-end suite
  cannot run as such: 13 of its 102 tests there fail before or beside the profiler (symlinks onto model files the
  database already has, a mock species' MSA, `unreported_species.tsv`, knob and depth checks tuned to the mini
  database), which the EM cannot touch; the reference build was not rerun there. Noted so nobody repeats it.

## Follow-up (2026-10-03): items 2-5 implemented, two of the alignment-loop changes dropped

At the user's request, after item 1. Scripts: `scripts/p5_build.sh` (items 3 and 5), `scripts/p4_build.sh` (item 4 on
top of them), `scripts/em_sweeps2.sh` with `sweeps.pl` (item 2's measurement); results in `results/items_2_to_5.txt`.

### Item 3: the congener sketches (identical outputs)

- `context::GeneSketch` hashes the k-mers as they come, repeats included, selects the 2 × 64 smallest hashes in linear
  time (`std::nth_element`) and sorts those; whenever they hold 64 distinct values these are the 64 smallest distinct
  hashes of the gene (every hash outside is at least as large as any inside), else the gene's hashes are all sorted
  (genes of many repeats). Before, every distinct k-mer was collected and sorted, then hashed, then partially sorted
  again. `ApplySampleContext` sketches every reference of the sample's congener pairs in one pass over the threads
  before the pairs are compared (`CongenerDistances::Sketch`), so no sketch is made twice by two pairs at once.
- Test `SampleContext.SketchesEqualTheReferencesToTheLastHash`: random genes, genes with Ns, short ones, repeats of 30
  k-mers (fewer distinct than the sketch holds), half-repeated genes, genes of one base: the same hashes as the first
  version, in the same order.
- Instructions (callgrind, `--profile_only` on the 500k SAM): `GeneSketch` 4.57 → 2.27 G, the profiling stage 25.3 → 23.0 G.
  `--profile_only` outputs identical on the five SAMs (`results/items_2_to_5.txt`).

### Item 5: the alignment loop (identical SAM)

Kept, all three verified by SAM text identical to the reference build's at one thread (500k pairs, the same reads as
single-end, Nanopore 90 Mb, PacBio 90 Mb) and sorted SAM text identical at six threads (500k pairs):

- (a) `IsAlignmentValid` on every alignment → the counts must cover the read and stay inside the gene (every alignment,
  O(1) from `AlignmentInfo`'s counts), and one alignment in 64 is walked base by base as before. 0.61 G of 17.08 G gone
  (with the `NextCompressedCigar` it called). The `info.Valid` pass that followed it (a count over the CIGAR) is
  covered by the same check.
- (b) `across_genes::Neighbours` keeps each thread's answers by (taxon, gene, end) for the table and database they came
  from (`gene_neighbours::Table::Generation`, which changes with every `Read` and `SetLineage`): `Table::Partners` and
  `Assess` (0.43 G, from mate guidance's search on the neighbour genes) gone.
- (d) `AlignmentEdits(AlignmentInfo const&)` from the counts `GetInstructionCountsAndCompress` took (the only way the
  compressed CIGAR is made), instead of parsing the compressed CIGAR again for every candidate's ZA entry.
- Together: `RunPairedEnd` 17.08 → 15.98 G per 100k pairs (−6.5%), the run 28.07 → 26.97 G.

Tried and dropped (both reverted before the commit, measured in `scripts/p5_build.sh`'s build):

- (c) the seeds sorted on one packed 64-bit key: the same order, but computing the key for both elements of every
  comparison cost more instructions than the three-field comparator's branches save; the anchor finder went 2.79 → 3.09 G
  and the "Sorting Seeds" timer moved from 0.70 to 0.64 s at one thread. Not worth +1.7% of the loop's instructions.
- (e) a WFA2 aligner of their own for short reads' flanks: `wavefront_slab_reap_repurpose` and `wavefront_slab_allocate`
  did not change at all (0.37 and 0.35 G). WFA2's per-call reset walks the wavefronts of the largest alignment the
  aligner made, and their number follows the score budget, not the sequence lengths, so a flank-sized aligner resets as
  many as the window's. Round 3's estimate for this item was wrong.

### Item 4: profiling a sample while the next one is aligned (identical outputs, no gain here, opt-in)

- Implemented as `ProfilingAhead` in `RunProtal.h`: `RunWrapper` hands over each sample whose SAM is complete (and each
  it skips because the SAM exists); a worker thread profiles them one after another on a quarter of the threads with
  `ProfileSample` (the body of the old `ProfileWrapper` loop, now a function any thread can call for a sample of its
  own, with `ProfilingContext` holding what the samples share); the profiling stage takes what the worker had not
  started, on all threads, and waits for it. The last sample is never handed over (no alignment follows it).
- Outputs: profiles, strain MSAs and per-species statistics identical between the overlapped and the sequential runs
  of the four-sample map (two rounds) and of a single sample; sorted SAM text identical for every sample. The
  eight-sample map differed in one value's sixth decimal (`ANISum` of one gene) between the two runs, which happens
  between any two six-thread runs of the same binary (the SAM's record order varies with the threads; round 3's note).
- Time, six threads, two rounds each (`results/items_2_to_5.txt`): four samples 13.7 s with the worker against 12.9 s
  without, then 14.8 against 15.9; eight samples 30.0 against 29.7. The "Processing all samples" stage grew by as
  much as "Profiling" shrank: on six vCPUs with the alignment on every one of them, the worker's CPU time comes out
  of the alignment's. After items 1-3 the profiling stage has little serial work left to hide, so there is nothing
  to gain unless the profiling stage leaves cores idle (many threads, a long tail of one deep sample). Hence
  **off by default**: `--profile_ahead` (dev option, in `docs/running.md`) turns it on, so that it can be measured on
  a cluster node with the same binary.

### Item 2: the EM's stopping rule (changes the last digits of the `em_*` features)

Measured first (`scripts/em_sweeps2.sh` inserts `scripts/sweeps_insert.cpp`; `results/em_sweeps.txt`): per sweep, the
largest relative weight change (the first rule's measure) and the largest change of any taxon's own share.

| SAM | sweeps | relative weight change at sweep 10 / 50 / 100 / last | own-share change at sweep 50 / 100 / last |
|---|---|---|---|
| 500k pairs | 200 (cap) | 0.046 / 0.0029 / 0.00046 / 0.00047 | 0.0029 / 0.00024 / 0.00023 |
| 5M pairs | 197 | 0.065 / 0.00039 / 2.2e-6 / 9.2e-11 | 0.00028 / 1.5e-6 / 6.6e-11 |
| Nanopore 90 Mb | 200 (cap) | 0.053 / 0.00027 / 0.00027 / 0.00026 | 0.00027 / 0.00026 / 0.00026 |
| PacBio 90 Mb | 200 (cap) | 0.036 / 0.00073 / 8.3e-5 / 1.2e-5 | 0.00062 / 8.3e-5 / 1.2e-5 |

- The 5M sample converges geometrically: the shares are stable to 1e-6 by sweep 100; the first rule (relative weight
  change below 1e-10) let it run to 197.
- The other three plateau: from sweep 50 on, some taxon's share moves by nearly the same 1e-5–2e-4 every sweep (the
  500k sample's change shrinks by 2.5% over sweeps 100–200: a geometric drain at a rate within 3e-4 per sweep of
  standing still). That is a taxon whose reads a congener explains about as well: the congener's weight times its
  edit factor is within a fraction of a percent of the taxon's reads, so each sweep moves only that fraction. No
  tolerance above the noise ends it within 200 sweeps; such taxa run to the cap under any rule, and their feature is
  whatever the cap leaves (as before).
- **The rule now**: stop when no taxon's own share (its reads left to it over its records, the feature the EM is for)
  changed by `kEmTolerance` = 1e-6 in the last sweep; `kEmIterations` = 200 stays. The `own` sum per sweep comes from
  the same products the weights are updated with, so the sweeps' arithmetic is unchanged. Samples that hit the cap are
  bit-identical to before (every sweep runs, as before); converging samples stop at about half the sweeps with shares
  within 1e-6 of the full run's. The exactness test now runs both versions with tolerance 0 (every sweep);
  `SampleContext.StopsWhenTheSharesAreStable` checks the early stop and the harmonic case.
- Effect on the outputs (`scripts/final_build.sh`, `results/items_2_to_5.txt`): none visible. The `--profile_only` outputs of the
  final build are identical to the EM build's on all five SAMs, the 5M-pair one included (its sweeps now end at about
  105 instead of 197: the shares differ by less than the features' printed precision); the four-sample map's profiles have
  the same calls as the reference build's, no taxon's `Predicted` or `Probability` differs, and the strain MSAs are identical.

### Times, final build against the references (alternated rounds; `scripts/fin2.sh`, `results/items_2_to_5.txt`)

The machine had slowed by 1.5–2× for these runs (the reference build's four-sample map took 30 s against 15.2 s in
the morning), so only the ratios within a round mean anything:

| run | reference | final | |
|---|---|---|---|
| `--profile_only`, 5M pairs, 6 threads, 3 rounds | 12.8 / 11.5 / 13.1 s (EM build) | 11.1 / 9.5 / 10.1 s | −14% (items 2 and 3 on top of item 1) |
| 500k pairs, whole run, 1 thread, 2 rounds | 22.9 / 23.0 s (27423c6) | 18.7 / 18.4 s | −19%; profiling 7.1 / 7.3 → 2.4 / 2.8 s, aligning the same |
| four-sample map, 6 threads, 2 rounds | 29.9 / 30.6 s (27423c6) | 29.1 / 24.5 s | profiling 8.7 / 9.1 → 6.8 / 3.7 s; the rest noise |

Against 27423c6 on the 500k-pair run at one thread the profiling stage is now a third of what it was (items 1–3); the
alignment stage gains −6.5% of its instructions (item 5), too little to see in wall time here.
