# Performance, round 3: what is left, and long reads through their chain

- **Date**: 2026-10-02.
- **Code**: branch `audit-fixes` at `aeb7bf4` (the one binary of `96d9c3d`).
  - Committed from this round, with identical outputs: `f343113`.
  - Not committed: long reads aligned through their chain, branch `long-reads-through-chain` (`scripts/long_reads_through_chain.patch`, on `f343113`). It changes long-read alignments, so it waits for a decision.
  - **Later the same day**: merged as `29369eb`, at the user's request once the other sessions had finished, after unit (295) and e2e (128) tests passed on it from a clean tree (as `35d46f1`, on 0.7.3).
- **Machine**: WSL2 Ubuntu 24.04 on an Intel Core Ultra 7 258V (4 fast and 4 low-power cores, 6 vCPUs), GCC 13.3.
  - Other sessions kept the load at 6–9 for most of the time.
  - Instruction counts (callgrind, one thread) are the measure. Whole runs were timed only for the long-read benchmark, with the two builds alternated.
- **Data**:
  - Short reads and 3 Mb long reads: the v0.7.1 benchmark world's database (`~/bench071/V071`). 100k pairs are the first of `rl150_p500000_s_1`, 500k pairs the whole sample; Nanopore and PacBio are the `ont_b3000000`/`pb_b3000000` samples (pbsim3).
  - HiFi: 3 Mb simulated with `scripts/hifi_reads.py` from 15 kb templates of community `rl150_p1000_s_1` (`scripts/templates.pl`, `hifi.sh`).
  - Long-read benchmark: the samples 0.7.2's collector made (`~/bench071/samples_lr072`), against the V072 database.
- **Question**: is anything else in protal worth optimising now?

## Summary

1. **Where a run's time goes now.**
   - Short reads (500k pairs, one thread, stage timers, wall time):
     - seed and anchor finding: 40%, including the k-mer lookups (19% of the time but ~6% of the instructions, so stalled on memory) and sorting the seeds (7%);
     - alignment: 33%;
     - joining pairs and output: 17%.
   - Per pair, alignment executes 170k instructions, 41% of them in WFA2.
   - Long reads: alignment is 88% of a Nanopore run (WFA2 11.7 G of 13.5 G instructions on 3 Mb).
   - Loading the index costs a fixed 9.6 G instructions per run (1.2 s at 6 threads on this database).
2. **Done, same outputs (`f343113`).**
   - Flanks with at most one mismatch are aligned without WFA2. That is provably what WFA2 returns (below), and 20,000 random flanks check it.
   - The index decoder fills runs of 64 empty blocks at once.
   - 100k pairs: 29.00 → 26.78 G instructions (−7.7%). Index decoding went 6.85 → 4.92 G, alignment 17.04 → 16.75 G.
   - The same SAM text and outputs on paired-end, single-end, Nanopore, PacBio and HiFi reads.
3. **Long reads through their chain (not committed, outputs change).**
   - Long reads are aligned as a whole window today, because their chains hold indels. They can be aligned from the anchor's exact matches as short reads are: through every link, with more links found where the read's seeds leave stretches, and with WFA2 end to end between links on an aligner of their own.
   - Instructions of the long-read alignment drop by 31% (Nanopore), 38% (PacBio, pbsim3) and 37% (HiFi).
   - On the long-read benchmark:
     - F1 is the same in every scenario;
     - the profiles differ by a Bray-Curtis of 0.003 or less from those of the whole-window alignment;
     - CPU time is 15–21% lower (90 Mb samples).
   - Alignment penalties are equal in total. 97.5–100% of records keep their MAPQ. 8–46% of CIGARs place an indel elsewhere among equally good places.
4. **Not done, by what they would give:**
   - prefetching the k-mer lookups across reads (19% of short-read time, memory-bound; not measurable on this machine);
   - PGO (−3…−8%, round 2);
   - the seed sort on one 64-bit key (same order, ~1–2%);
   - a WFA2 aligner of their own for short reads' flanks (the per-alignment reset is 3% of short-read alignment, identical outputs);
   - the profiling stage's `CoveredPortion` and `VariantHandler::AddAlignment` (13% and 8% of profiling, itself ~10% of a run);
   - index decoding at GTDB size (the value cells, which grow with the database).

## Where the time goes

**Stages** (`results/stage_times_1thread.txt`, the timed 500k-pair run of the one-binary report, plain build, one thread; seconds of wall time):

| stage | s | share |
|---|---|---|
| seed and anchor finding | 5.27 | 40% |
| – seeding (the k-mer lookups) | 2.48 | 19% |
| – sorting seeds | 0.91 | 7% |
| – extending anchors / pairing / sorting anchors | 0.90 / 0.48 / 0.11 | 11% |
| alignment handler | 4.33 | 33% |
| joining alignment pairs and sorting | 1.18 | 9% |
| output handler | 1.01 | 8% |
| retrieving k-mers | 0.74 | 6% |
| sequence reader | 0.53 | 4% |

- The seeding is 19% of the time for about 6% of the instructions (`FindSeeds` and `GetFromLookup`, 0.9 G of 17 G), so the k-mer lookups wait on memory.
- Round 2 (`../2026-09-30-performance-round2/README.md` §4.5) found the same and could not explain why only 2–3 lookups overlap. Its prefetch prototype gave −4% of the seeding timer.
- What it left untried is prefetching the next read's lookups while a read aligns. A few percent cannot be measured on this machine.

**Instructions**, 100k pairs (`results/callgrind_pe100k_head_self.txt`):

- Alignment (`RunPairedEnd`): 17.0 G, 170k per pair.
  - WFA2: 7.05 G (41%), most of it in the anchored aligner's flanks (6.5 G with what they call). Its per-call reset (`wavefront_slab_reap_repurpose`) is 0.49 G.
  - Chaining: 3.0 G.
  - Output: 1.7 G.
  - Syncmers: 1.2 G.
- Index load: 9.6 G. That is `DecodeChunk` 6.85 G and zstd 3.65 G.

**Nanopore** (3 Mb, `results/callgrind_ont3M_head.txt`):

- Aligning: 13.5 G, of which WFA2 11.7 G: extension 4.3 G, wf-adaptive 2.5 G, ends-free termination 1.7 G.
- 1,827 reads give 3,849 gene alignments (about 3 M instructions each), of which 1,564 are written.

## Done: the same outputs, fewer instructions (`f343113`)

**Flanks of one mismatch without WFA2** (`AnchoredAligner::Flank`).

When is it safe? A flank starts at its link, fixed. Suppose the link's diagonal, followed to the end of the read or of the reference, has at most one mismatch, and the bases left at the other end fit its free end. Then:

- That ungapped alignment costs 0 or 4.
- Any other alignment leaves the diagonal with a gap, which costs 8 or more, so the ungapped one is the only best.
- WFA2 ends an ends-free alignment only on the last row or column, within the free counts (`wavefront_termination_endsfree`). It writes the bases left over as trailing `I` (read) or `D` (reference) (`wavefront_backtrace_affine`).
- It stops when the score reaches `max_steps` (`wavefront_unialign_reached_limits`), so it fails where 4 × mismatches reaches the budget.

The shortcut gives exactly that.

- `AnchoredAlignment.UngappedFlanksAreWhatWFA2Gives` compares it with WFA2 on 20,000 random flanks with up to 4 changes and Ns, free ends on either side and budgets of 1–30: the same status and operations every time, 16,421 flanks without WFA2.
- On 100k pairs: alignment 17.04 → 16.75 G (−1.7%).

**Empty index blocks 64 at a time** (`index_codec::detail::DecodeChunk`).

- The key map's control blocks are mostly empty on a small database. A run of them is the same 24 bytes repeated.
- Where a bitmap word says 64 blocks are empty, they are now written as a 96-byte pattern, instead of being tested and stored one by one.
- `DecodeChunk`: 6.85 → 4.92 G per run, the same on every workload, since it is the index load.
- `IndexCodec.LongRunsOfEmptyBlocks` round-trips runs of 63 to 300 empty blocks at many offsets to the chunks.

**Checks** (`results/identical_changes_checks.txt`):

- SAM text identical to `aeb7bf4`: 100k pairs at `-t 1`, single-end 500k at `-t 1`, Nanopore, PacBio and HiFi 3 Mb at `-t 1`.
- 500k pairs at `-t 6`: the sorted SAM is identical. The profile is the reference's own on that SAM.
- Unit tests: 276. The commit also passes from a clean `git archive`.

## Long reads through their chain (not committed)

### Why the long reads take the whole window

The anchored aligner (round 1, `3a9cfaa`) takes only chains on one diagonal.

- A long read's indels shift its seeds' diagonal along a gene, so every long read went to the whole-window WFA2 alignment.
- `RunProtal` switched anchored alignment off for long reads altogether.

### What its chains hold

`results/long_read_chains.txt`, 3 Mb samples:

- The links of an aligned candidate span 38% (Nanopore) and 41% (PacBio) of its gene, and cover 24%.
- That leaves on average ~420 gene bases before the first link and ~440 after the last.
- Merging a read's anchors of one gene would not help. There are 1.07 anchors per read, gene and strand, and all of them together span 32–35% of the gene. Few of the index's k-mers survive these reads' error rates.
- 4% of the alignments repeat a read, gene and strand already aligned, from a second anchor.

### Steps, measured

Instructions of the long-read alignment (`RunLongReads`), G, one thread (`results/callgrind_parts.txt`):

| step | Nanopore 3 Mb | PacBio 3 Mb |
|---|---|---|
| reference: the whole window | 13.47 | 10.71 |
| through the chain: links in read order, overlaps cut, end-to-end WFA2 between links on different diagonals | 10.65 | 7.78 |
| + re-seeding: 12-mers unique in a band around the neighbouring link's diagonal | 10.83 | 7.81 |
| + the pieces between links on an aligner of their own, with WFA2's end-to-end kernels | 9.43 | 6.73 |
| + the middle before the flanks (the budget shared; dropped: no gain) | 9.69 | 6.87 |
| final, with `f343113`'s flank shortcut | **9.32** (−31%) | **6.67** (−38%) |

HiFi 3 Mb: 9.74 → 6.15 G (−37%).

Re-seeding alone did not pay: it moved the work from the flanks into many short pieces.

- WFA2 keeps the wavefronts of the largest alignment an aligner made, and resets every one of them for each alignment: `wavefront_slab_reap_repurpose`, about 55k instructions per call here (`results/wfa2_calls_and_costs_ont3M.txt`).
- The pieces shared the aligner with gene-long flanks, so each paid that.
- An aligner of their own keeps that small, and end-to-end kernels skip the ends-free bookkeeping.

After re-seeding, a Nanopore read's flank holds a median of 16 gene bases (75th percentile 40). One in ten holds more than 200, which are mostly relatives' genes where 12-mers rarely match. Those flanks are most of what is left (4.2 G).

**Reads left to the whole window.** The window's free ends come from the first link's diagonal
(`SimpleAlignmentHandler::AlignAnchor`, a dovetail of 9 bases).

- Where a read's net indels move the diagonal by more than that, more read bases lie past the gene's end than are free there.
- The anchored alignment then leaves the read to the whole-window alignment, which gives it what it gave before.
- In the unit test (`LongReadsThroughTheirChain`, Nanopore-like and HiFi-like reads) that was 9 of 60 reads, all at the right end.
- Taking the right end's free bases from the last link's diagonal would align those reads too, and clip their overhang instead of forcing it into the alignment. That changes their alignments, so it is left for the decision on this change.

### Alignments

Records of the same 3 Mb runs (`results/long_read_records.txt`), the prototype against the whole window:

| reads | records | same CIGAR | same MAPQ | same POS | penalty (total) |
|---|---|---|---|---|---|
| Nanopore | 1,564, one each way moves to another gene | 62.4% | 97.5% | 98.4% | −0.00% (2 lower, 2 higher) |
| PacBio (pbsim3) | 1,103, all the same genes | 53.6% | 99.5% | 99.8% | +0.00% |
| HiFi | 1,242, all the same genes | 92.5% | 100% | 100% | +0.00% |

The CIGARs that differ place indels elsewhere among equally good places: all but 4 of them (Nanopore) have the
same penalty.

### Profiles

`scripts/lr_bench.sh`, `lr_score.pl`; `results/long_read_benchmark.txt`, `long_read_benchmark_times.tsv`. Both builds ran on the 16 long-read samples against the full V072 database at `-t 6`, alternated per sample:

| scenario | F1 reference | F1 prototype | Bray-Curtis between them | profiles identical | user CPU s | alignment s (wall) |
|---|---|---|---|---|---|---|
| Nanopore 3 Mb | 0.669 | 0.669 | 0.003 | 0 of 4 | 10.7 → 8.3 | 1.0 → 1.1 |
| PacBio 3 Mb | 0.535 | 0.535 | 0.000 | 3 of 4 | 9.4 → 6.4 | 0.8 → 0.6 |
| Nanopore 90 Mb | 0.889 | 0.889 | 0.0005 | 0 of 4 | 146.0 → 124.7 | 39.1 → 37.3 |
| PacBio 90 Mb | 0.870 | 0.870 | 0.000 | 1 of 4 | 98.2 → 77.8 | 26.4 → 21.2 |

- F1 here is species-level on the truth lists, simpler than the benchmark's `score.py`, hence lower than its figures. 0.7.2's own runs score the same here as HEAD's: 0.669, 0.535, 0.889, 0.870.
- The 3 Mb CPU times include the index load, which `f343113` also shortened in the prototype's build.
- Wall times at six threads on this machine are indications only.

### Why it is not committed

- It changes long-read alignments (same penalties and calls, other indel placements, a few MAPQs).
- Another session is working on long reads now: phasing, models, the r226 build. A change under its benchmarks would confound them.
- The branch has the change cleanly, with unit tests:
  - `AnchoredAlignment.LongReadsThroughTheirChain`: 60 reads at two error rates against the whole-window alignment;
  - `AnchoredAlignment.OverlappingLinksOfALongReadAreCut`.
- `docs/running.md` there describes it.

## Not done

| # | what | evidence | expected | outputs |
|---|---|---|---|---|
| 1 | prefetch the next read's k-mer lookups while a read aligns | seeding 19% of short-read wall time for ~6% of the instructions; round 2's `bench_mem` gained 20% that way | up to ~10% of short-read time; needs a quiet machine or bare metal to measure | same |
| 2 | PGO in the release build | round 2: −3% (mix), −8% (w900) | −3…−8% | same |
| 3 | sort seeds by one 64-bit key (taxid, gene, read position packed, as `SortByReadComparator2`) | sorting seeds 7% of short-read time; the same comparison results give `std::sort` the same order | ~1–2% | same |
| 4 | an aligner of their own for short reads' flanks | `wavefront_slab_reap_repurpose` 0.49 G on 100k pairs (3% of alignment) | ~2% | same |
| 5 | profiling: `SequenceRangeHandler::CoveredPortion`, `VariantHandler::AddAlignment` | 13% and 8% of profiling instructions; profiling ~10% of a run | ~2% of a run | same if order-preserving |
| 6 | the value cells of the index decoder at GTDB size | up to 10 byte planes per value cell; ~4.4 billion cells at GTDB size (35 GB) | seconds per run at GTDB size | same |
| 7 | long reads: the window's free ends from the last link's diagonal (with the change above) | 9 of 60 test reads fall back to the whole window | fewer whole-window alignments; better clipping | changes those reads' alignments |

## Notes

- At more than one thread the SAM's records come in an order that varies from run to run. The profile then sums some values in another order, and the gene log's `ANISum` can differ in its sixth decimal between two runs of the same binary (seen once here, between the prototype and the reference). This was already so before this round.

## How it was run

`S` in the scripts is the scratch folder they ran from. `build_wt.sh` and `check.sh` are those of
[the multithreading audit](../2026-10-01-multithreading-audit/scripts/followup2/).

```bash
bash scripts/setup.sh                       # HEAD's archive, built, in ~/mt-work/perf3/ref
bash scripts/lr_run.sh                      # the first through-chain prototype; records compared (lrcmp.sh)
bash scripts/cg_lr.sh                       # callgrind of it; calls.sh, calls_all.sh: WFA2's calls and costs
bash scripts/lrdbg.sh                       # chains and anchors per gene (lrdbg.patch: a debug print, experiment only)
bash scripts/lr_run2.sh; bash scripts/lr_run3.sh   # re-seeding, the piece aligner, the middle first
bash scripts/flk.sh                         # flank sizes with and without re-seeding
bash scripts/hifi.sh                        # 3 Mb of HiFi reads (templates.pl, scripts/hifi_reads.py)
bash scripts/lr_run4.sh                     # the combined prototype: short reads identical, records, callgrind
bash scripts/id_check.sh                    # f343113's changes: unit tests, outputs identical, callgrind
bash scripts/lr_bench.sh; perl scripts/lr_score.pl; perl scripts/f1_v072.pl   # the long-read benchmark
bash scripts/lr2_check.sh                   # the long-read change as on its branch: check.sh, records
bash scripts/collect.sh                     # results/
```
