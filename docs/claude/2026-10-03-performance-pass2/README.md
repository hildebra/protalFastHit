# Performance pass 2 (2026-10-03): where the time goes at `0e4c5a4`, and four cheap exact savings

Based on `0e4c5a4` (audit-fixes), built in WSL from a copy on the Linux filesystem (`scripts/build_head.sh`), the
v0.7.3 benchmark world's database (`~/bench071/V073/protal_db`) and its simulated paired-end sample of 500,000
pairs (`rl150_p500000`), 1 and 6 threads, callgrind on 100,000 of those pairs (aligning) and on the 500k-pair SAM
(profiling). Instruction counts are exact; wall times here vary by 5-10% between identical runs (the machine's
vCPUs move between fast and slow cores), so only alternated rounds are compared. Scripts are in `scripts/`, raw
numbers in `results/`. The earlier rounds are in [round 3](../2026-10-02-performance-round3/README.md),
[prefetch](../2026-10-02-prefetch/README.md) and [the review before this one](../2026-10-03-performance-review/README.md).

## Where a run's time goes now

500k pairs, whole run (`results/stages_500k_t*.txt`):

| | 1 thread | 6 threads |
|---|---|---|
| index + genome preload | 1.8 + 0.13 s | 0.95 + 0.11 s |
| aligning the reads | 9.4 s (72%) | 3.0 s |
| profiling | 2.05 s (16%) | 1.13 s |
| whole run | 13.0 s | 4.5-5.6 s |

Six threads give 3.1× on the aligning, 1.8× on profiling, 2.4-2.9× on the run: the machine has 4 fast and 4 slow
cores shown as 6 vCPUs (CPU time at six threads is 19-23 s against 10-12 s at one), and the load is serial.

Aligning 100k pairs at one thread is 26.8 G instructions, of which the index load is 8.4 G (a fixed cost: 7.3 G
`Seedmap::Load` = `index_codec::DecodeChunk` 4.4 G + zstd 2.6 G, genomes 1.1 G) and `RunPairedEnd` 15.8-16.2 G:

| part of `RunPairedEnd` | G instructions | share |
|---|---|---|
| alignment handler, of which `AnchoredAligner::Flank` (WFA2 on the read's flanks) | 8.1, 6.2 | 50%, 38% |
| anchor finder (lookups, extending) | 2.8 | 17% |
| SAM output handler | 1.6 | 10% |
| syncmer scan (`ScanWindowsAvx2` 0.83) and k-mer handler | 1.2 | 7% |
| input inflate (second thread) | 1.1 | 7% |

## What the flanks cost (a new measurement)

`AnchoredAligner::Flank` aligns the read's two flanks around its anchor link with WFA2 unless the flank has at most
one mismatch on the link's diagonal. A counting build (`scripts/flank_stats.sh`, `results/flank_stats.txt`) shows,
for 100k pairs, about 180,000 flanks sent to WFA2 (6.1-6.7 G instructions, ~35k each):

- 47% have 2 or 3 mismatches (84k flanks) and cost about 7k each; 53% have 4 or more and cost ~90% of the WFA2 time.
- The ungapped alignment is WFA2's answer for 100% of the 2-mismatch flanks, 99.9% at 3, 99.5% at 4, 98.7% at 5,
  97.5% at 6, 96% at 7, 94% at 8, 92% at 9, 90% at 10, 87% at 11, 82% at 12.

## Changes (all exact: SAM text and profile outputs identical on the 500k-pair sample)

1. **`SequenceRange::CoverageVector` as a difference array** (`SequenceRange.h`). It added every read to every
   base it covers; `CoveredPortion` (9.4% of a profiling run's instructions on this sample: `CoveredBases` and the
   feature strings call it per taxon and gene) and `CalculateCoverageVector2` (3%) go through it. Now two updates
   per read and one prefix sum: the same counts (test `SequenceRange.CoverageVectorEqualsAddingEveryBase`, 200 random
   ranges, every strand and divergence filter).
2. **Flanks with 2 or 3 mismatches without WFA2** (`AnchoredAligner::NoCheaperGap`). The ungapped alignment of a
   flank costs 8 or 12; an alignment with two gap runs costs 16 or more, and one with a single run (of L bases:
   6 + 2L, plus 4 per mismatch before and after the change of diagonal, the cheapest change position found by one
   pass per shift) is checked for shifts of up to 3 bases. The flank is taken as ungapped only when every such
   alignment costs more; a tie goes to WFA2. Bases without a partner on the shifted diagonal are counted as free,
   which makes the check stricter, never looser, with free ends. Test
   `AnchoredAlignment.FlanksWithTwoOrThreeMismatchesAreWhatWFA2Gives`: 30,000 reads (genes of random bases and of
   short repeats, indels in the read, free ends of both sides, budgets down to the edge) give WFA2's status and
   operations, 42,585 flanks without WFA2.
3. **`AbundanceWeightedShares` without a heap vector per class visit** (`SampleContext.h`): the shares of a class go
   to a stack buffer; the arithmetic and its order are unchanged (the bit-exact test passes).
4. **`VariantHandler::AddAlignment`: an exact run of M without Ns is accepted with one `memcmp` and `memchr`**
   instead of a base-by-base loop.

Effect, one thread (`scripts/compare.sh`, `scripts/cg2.sh`):

| | before | after |
|---|---|---|
| aligning 100k pairs, instructions | 27.15 G | 26.76 G (`RunPairedEnd` 16.16 → 15.77 G, −2.4%; WFA2 6.73 → 6.11 G) |
| profiling the 500k-pair SAM, instructions | 24.05 G | 20.04 G (−16.7%) |
| `Raw alignment` timer, 500k pairs | 2.62-2.90 s (6 runs) | 2.53-2.68 s (5 runs) |
| profiling stage | 2.03-2.11 s (6 runs) | 1.92-1.99 s (4 runs; one 2.19) (−6%) |
| SAM, profile outputs | | identical |

The wall gain is small (item 2 is 1.4% of the run, items 1, 3, 4 6% of the profiling stage, ~1% of the run at one
thread; more on deep samples, where profiling is 30-40% of the run). All 329 unit tests pass (one skip, as before).

## What is left, ranked by what it could give

1. **Flanks with 4 or more mismatches: 38% of the aligning loop's instructions.** The ungapped alignment is optimal
   in 82-99.5% of them, but no cheap exact certificate exists once two gap runs can tie (16 or more): it needs a
   banded search, which is what WFA2 does. Reducing WFA2's per-call work (it resets and expands wavefronts for a
   score budget of up to the ANI floor) is the lever; limiting the budget changes which alignments exist.
2. **A per-read budget from the best candidate.** `SimpleAlignmentHandler::operator()` aligns up to `align_top`
   candidates each against the full ANI-floor budget (up to 60 for 150 bp). The profiler reads all of them
   (alternatives within 5 edits, the ambiguity classes of the EM), so a budget of best + margin would save WFA2
   work on far relatives at the price of fewer listed candidates: not exact, to be benchmarked on F1 first.
3. **Index load: 25-27% of a 100k-pair run, 1.8 s at one thread, 0.95 s at six, more at GTDB scale.**
   `DecodeChunk` is 4.4 G on its own (three passes over the keys, `count`/`FlexCells` per key each time) and zstd
   2.6 G. A single pass, or decoding the planes with SIMD, could take a third of it; it needs the GTDB-sized index
   to be worth the time.
4. **`NoteRecord` reads each record's CIGAR five times** (`DifferencesAndAligned`, `ReferenceBases`,
   `MismatchesByCodonPosition`, `ReadExcess`, `MeanErrorProbability`): 1.2 G, 6% of profiling.
5. **Congener sketches** (`MakeSketch` 2.3 G, 11% of profiling here): a sketch depends only on the database, so
   the build could store them (a database format change); or a table for the base codes in `BottomSketch` (−30%
   of its 1.4 G).
6. **`AbundanceWeightedShares` is still 2.2 G (11%)**: classes are visited 105 times; the division per share and
   the random access into the taxa's weights are what is left.
7. The SAM written by the aligning stage is compressed and read back by the profiling stage: 2.7 G of zstd
   decompression per 500k pairs, 13% of profiling.

Not measured: GTDB scale (the index load, the sketches and the EM grow with the database and the sample); it needs
`scripts/measure_performance.sh` on a cluster node.
