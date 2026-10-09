# ONT long reads: skipping the WFA2 alignments that cannot succeed

Date: 2026-10-09. Branch `audit-fixes` at `7ace9b8` (`src/Core/AlignmentStrategy.h`, `src/Core/LongReads.h`,
`src/Core/ChainAnchorFinder.h`, `src/Alignment/AlignmentScreen.h`, `src/Alignment/AnchoredAlignment.h`,
`src/Alignment/WFA2Wrapper2.h`, `lib/wfa2-lib`).

Data: the r226 v19 build's protal run logs (`local/v19/.../logs/protal_runs_*.log`, 84 threads; SLURM job 24099010).

Method: code reading only, by a read-only agent and checked here against the code; nothing was built or run. A
follow-up to the [build parallelism report](2026-10-09-build-parallelism/README.md), which found ONT alignment to be 70%
of the build's protal time.

**Question.** ONT alignment took 6,001 of the 8,587 s of aligning in the v19 build. What are the options for skipping
more of its WFA2 alignments?

## Summary

Nearly every ONT candidate goes through WFA2, and on host reads 99.8% of them fail. The k-mer screen
(`AlignmentScreen.h`) cannot refuse any of them at ONT's identity floor of 0.85: its bound is exact by the q-gram lemma,
and at 15% divergence an alignment within the budget may touch every k-mer of the window.

**The best option keeps every output the same:** a second exact bound before WFA2, from the indel distance. A failing
candidate is mostly a gene window against unrelated read sequence, and its indel distance is far above the budget. A
bit-parallel LCS of the read part and the window costs ~20-40 µs, against ~590 µs for the WFA2 alignment it would save.

Estimate: ONT alignment 60-80% faster, about 1 h less protal time per r226 build at 84 cores. Confidence medium: it
assumes the failing ONT candidates are as unrelated as PacBio's, of which the screen refuses 96%. That is not measured
yet.

| Rank | Option | Outputs | Saving on ONT alignment | Effort | Confidence |
|---|---|---|---|---|---|
| 1 | **Exact indel-distance bound** (bit-parallel LCS) before WFA2 | identical | −60-80% | medium, ~1 day with tests | medium |
| 2 | Probabilistic k-mer floor in the screen (long reads only) | change, rarely | about as 1 | small | medium |
| 3 | Minimum anchor evidence per candidate (seeds, or exact bases per gene length) | change (ZF; strains' FN risk) | large alone, ~3-5% after 1 | small | low |
| 4 | Pruning by chain score against the read's best, or a cap per segment (minimap2's secondary chains) | change (MAPQ, ZA, ZF) | up to 2/3 alone, small after 1 | small | low |
| 5 | `--long_read_budget` (exists, off) | change | < 5% at r226 | none | medium |
| 6 | Two-stage alignment (a short stretch around the best anchor first) | change | dominated by 1 and 2 | medium | low |
| 7 | WFA2 heuristics (z-drop or x-drop, a narrower wf-adaptive band) | change, succeeded alignments too | ~2× on failing candidates | small | low |
| 8 | Host or low-complexity filters | change | ~0 | n/a | medium-high |

## 1. Where the ONT time goes

**The cost per candidate is flat.** Each candidate takes 0.59-0.60 ms of thread time, whatever the sample:

| Sample | Candidates | Aligning | Thread time per candidate | Candidates per k-mer looked up |
|---|---|---|---|---|
| `sc_host_ont_b3000000000_s_1` (90% human reads) | 17,926,550 | 126.5 s | 0.59 ms | 0.0198 |
| `sc_moderate_ont_b3000000000_s_1` | 5,692,086 | 40.1 s | 0.59 ms | 0.0201 |
| `sc_host_pb_b3000000000_s_1` (PacBio, the control) | 17,846,585 | 16.7 s | 0.08 ms | 0.0197 |

**The candidates come with the reads, not the content.** There are about 6 candidates per kb of read, in human and
bacterial reads alike:
- An anchor needs two seeds within ±6 diagonals (`ChainAnchorFinder.h:153`, `:493`).
- Every place on the read aligns its `--align_top` 3 longest anchors, plus ties (`LongReads.h:557-572`).

**PacBio shows what is avoidable.** The PacBio host sample does the same seeding work as the ONT one:
- 908M against 904M k-mers looked up;
- 64.5M against 62.9M anchors;
- 17.85M against 17.93M candidates.

Its screen refuses 96.4% of them, and it aligns in 16.7 s. About 110 of the ONT sample's 126.5 s are spent aligning
candidates of the kind the PacBio screen refuses.

**A failing candidate pays the whole budget.**
- **The budget.** `MaxScore = ceil(0.15 × overlap) × 4 + 1` (`AlignmentStrategy.h:461`, applied at `:641`). That is
  601 for a 1 kb gene, and ~7.7k for the longest genes.
- **Why failures reach it.** An unrelated window costs ~1.4-2 per base, so WFA2 runs to the budget. It stops only at its
  step limit; the default wf-adaptive band trims but does not end it.
- **One budget per candidate.** The anchored aligner spends one budget across the left flank, every link and gap, and
  the right flank (`AnchoredAlignment.h:84-127`), so a failing candidate costs one full exploration.

## 2. Why the k-mer screen refuses nothing at 0.85

Take a 1 kb gene:
- the window is ~1,300 bases (margin 100 + L/20, `LongReads.h:247-268`), with ~159 free read bases at each end, so
  982 bases must align;
- `max_score` is 601;
- `MaxTouched(k = 7)` is max(600 × 7 / 4, 600 × 7 / 8, 600 / 2) = 1,050 k-mers;
- so `Required` = 976 − 1,050 < 0, and `Bound` returns "may align" before counting a k-mer (`AlignmentScreen.h:242-250`).

In general `Required` ≈ L(1 − 0.15k), negative for every k ≥ 7: 150 mismatches spaced every 6.7 bases destroy every
7-mer. Smaller k does not rescue it:

| k | Shared k-mers required | Shared by chance with an unrelated 1 kb window |
|---|---|---|
| 6 | 7.9% | ~22% |
| 5 | 25% | ~62% |

Bounds per link or per sub-window add only the clamped deficits, and chance sharing is even along the window, so they
gain nothing either. **No exact q-gram variant refuses at 0.85.**

## 3. Option 1: an exact indel-distance bound

**The bound.** With protal's penalties (mismatch 4, gap opening 6, extension 2), an alignment's score is at least twice
the indel distance of the parts it aligns:
- a mismatch (4) is a deletion and an insertion;
- a gap of g bases costs 6 + 2g ≥ 2g.

So, with n_req and m_req the read and window bases outside the free ends, any alignment within the budget needs
`2 × (n_req + m_req − 2 × LCS(read part, window)) < max_score`. A smaller sub-part can only have a smaller LCS, so the
bound holds for every alignment with those free ends. It also holds for the anchored path through the links, which is
one such alignment, constrained, under the same shared budget. A refused candidate is one that both the whole-window
WFA2 and the anchored aligner would fail: **the outputs are the same**.

**What it refuses.**

| Candidate | Lower bound | Budget | Refused? |
|---|---|---|---|
| unrelated 1 kb window (LCS of random sequence ~0.70-0.74 of the shorter) | ~1,000-1,160 | 601 | yes, margin ~1.7× |
| unrelated 300 bp gene | ~228 | 181 | yes, margin ~1.26× |
| congener at 85-90% ANI with ONT's indels | ~0.4-0.5 L | ~0.6 L | no (these fail cheaply in WFA2 anyway) |

**The cost.** A bit-parallel LCS (Allison-Dix / Hyyrö) over 1,300 window columns and 16 machine words is ~21k word
operations: ~150k instructions, 20-40 µs, against ~590 µs for the WFA2 alignment.
- The four base masks can come from the gene's packed bytes (`gene.Packed()`, as the k-mer screen's window does), so a
  refused window is never decoded.
- Ns, or the read's 2-bit packing, can only lower the LCS bound's left side. The bound stays valid.

**Where.**
- A function beside the k-mer screen in `AlignmentScreen.h`, called right after it in `AlignAnchor`
  (`AlignmentStrategy.h:~670`).
- Gated on long reads; harmless for PacBio, where the k-mer screen refuses most first.
- Its own counter in the per-sample counts line (`RunProtal.h` `PrintAlignmentCounts`), e.g. "refused by the
  indel bound".

**Outputs to keep the same.** The refusal happens inside `Align`, after the candidate is marked as tried
(`LongReads.h:564`, also `:395`, `:597`). So the following stay as they are:
- the failed taxa (`:481-492`) and the ZF tags (`:756`);
- the header's failed-candidate counts (`:704`);
- the profile's `FailedCandidateRate` and `failed_gene_share` (`Profiler.h:1487`, `:1436`).

Only the counts line moves: candidates move from "aligned from the anchor's exact matches" to "refused". This has to be
confirmed in the implementation, by the tests below.

**Saving.**
- ONT host sample: about 13 s of seeding (PacBio's non-alignment part), plus ~7% of the 110 s for the candidates not
  refused, plus 17.9M × ~25 µs / 84 for the bound itself ≈ 25 s against 126.5 s, about −80%.
- Gut-like samples should be similar or a little less.
- Over the build's 6,001 s of ONT aligning: −3,600 to −4,800 s, 45-55% of all aligning.

**Validation.**
- Unit:
  - extend `RefusesOnlyWhatWFA2Fails` (`tests/test_AlignmentScreen.cpp:82`) with ONT-like errors (indels,
    homopolymers, 5-25%) and random windows, for whole-window and anchored alignment;
  - extend `TheHandlerAlignsTheSameWithAndWithoutIt` (`:291`) to long-read windows.
- End to end: the ONT and PacBio bench samples with the bound on and off. The sorted SAM records and every profile
  output must be byte-identical, as was checked for the k-mer screen (2026-10-04 report, item 2).
- Cluster: the new counter on `sc_host_ont` and `sc_gut_ont` in the next build's protal runs, and the ONT samples' rows
  of `protal_cpu_<collection>.tsv`.

## 4. The other options

**2. A probabilistic floor in the k-mer screen.**
- For long reads: required = max(Required, f × k-mers) in `Bound`.
- A passing alignment keeps ≥ 0.85^7 ≈ 32% of its 7-mers with even errors, and more with ONT's clustered ones.
- Chance sharing grows with the window: ~6% at 1 kb, 26% at 5 kb, ~55% at 13 kb. So f must scale with the window's
  chance rate, and long windows need k = 8-10 (`kMaxK` is 10).
- Cheaper than option 1 (~33k instructions per candidate, as the HiFi screen measured), but inexact: an alignment WFA2
  accepts with few shared k-mers is lost.
- A quick experiment in shadow mode: count the candidates it would refuse that WFA2 aligns, on one build. Then an A/B of
  the ONT model.

**3. Minimum anchor evidence.**
- Require ≥ 3 seeds, or exact bases ≥ x% of the gene, in the candidate loop (`LongReads.h:560`).
- Mark the skipped candidates as tried, so ZF stays as it is when they would have failed.
- Do not drop them in `AddCandidate` (`:251`): the ZN crowding (`:540-554`) counts every one.
- The risk: divergent strains whose gene got two seeds would be lost, and the FN anatomy found the misses to be mostly
  strains with ≤ 2 fragments.
- After option 1 it saves only the bound's ~25 µs per refused candidate.

**4. Pruning against the read's best chain, or a cap per segment.**
- It does nothing on reads without a real gene (~75% of gut reads, ~90% of host reads), where every place is spurious.
- Lowering `--align_top` from 3 for ONT changes MAPQ (best against second), ZA and ZF.

**5. `--long_read_budget`.**
- It applies only after a success in the segment, and only 1-3% of ONT candidates succeed at r226.
- Locally it was −15% for ONT, in a world where 60% of candidates align, with 17,564 records changed (2026-10-04, item
  4d). At r226, under 5%.

**6. Two-stage alignment.**
- Inexact (ONT's error bursts), and dominated by options 1 and 2.
- It cannot be made exact: a sub-stretch's cost bounds the whole budget only weakly, and random 200-base stretches fit
  inside it.

**7. WFA2 heuristics.**
- x-drop was dropped for long reads because it lost the own species' alignments of genes a read ends in
  (`WFA2Wrapper2.h:22-24`), and z-drop works the same way.
- A narrower wf-adaptive band would roughly halve the cost of failing candidates, but changes the CIGARs of succeeded
  alignments everywhere.

**8. Host or low-complexity filters.** Candidates per k-mer are the same in human-dominated and gut samples, and only
4,711 lookups were dropped as too ubiquitous in the host sample. The candidates are not repeat-driven, so no gain is
expected.

## 5. Prior work this builds on

- 2026-10-04 (GTDB scale): the k-mer screen refused 20% of local PacBio candidates and none of Nanopore's.
  `--long_read_budget` stayed opt-in. WFA2's AVX2 kernels gave nothing on long reads. The window's right end comes from
  the last link.
- 2026-10-06 (profiling): the screen's read side and the packed window; for ONT "the screen never refuses".
- Round 3 (`29369eb`): long reads go through their chain, with re-seeding. 4% of alignments repeat a read, gene and
  strand from a second anchor, an exact but small deduplication.
- After option 1, the next ONT lever is seeding: ~13 s of the host sample's remaining ~25 s.

## 6. Recommendation

Implement option 1 with the tests above, and compare outputs byte for byte on the bench ONT and PacBio samples. Then
read its counter and the ONT samples' CPU rows in the next r226 build. Option 2 is worth a shadow-mode count only if the
bound refuses much less than PacBio's 96%. The others change outputs for less.

Nothing here changes how protal is run. The website is not affected until option 1 adds its counter to the per-sample
log line.

## 7. Implemented (follow-up, 2026-10-09)

The same day the user asked for option 1, and then for the bound as ONT's only screen.

**The bound** (`src/Alignment/AlignmentScreen.h`, class `IndelBound`):
- For ONT reads, in place of the k-mer screen, before the window is decoded (`AlignAnchor`; RunProtal sets it per
  read type). `--no_indel_bound` gives ONT reads the k-mer screen back, as before; `--no_alignment_screen` turns
  both off.
- Why in place of the screen: the screen refused none of the bench sample's 145,184 ONT candidates (1,191 of 17.9M on
  r226's host sample). A first version kept the screen and ran the bound only where the screen's bound was vacuous
  (every ONT window above 512 bases at 0.85). But in shorter windows the screen counts k-mers, and it still refused
  none there, so that gate lost 2,254 of the bound's 27,180 refusals on the bench sample. The user then asked for the
  bound as ONT's only screen.
- Not for PacBio reads, which keep the k-mer screen: the candidates it passes mostly align, and the bound on them
  cost 20% more instructions than it saved (below). Section 3's "harmless for PacBio" was wrong.
- With q = min(mismatch, 2 gap_extension) = 4, a candidate is refused when
  `q (n_req + m_req − 2 LCS(read, window)) > 2 (max_score − 1)`. n_req and m_req are the read and window bases outside
  the free ends, and the LCS is that of the whole read stretch and the whole window. This is section 3's bound in
  integers, with the penalties taken from the aligner.
- The LCS runs bit-parallel over the window's bases, read base by read base (V = (V + (V & M)) | (V & ~M), the LCS
  the zeros of V). The match masks come from the gene's packed bytes. A base other than A, C, G or T matches anything on
  either side, so the LCS can only grow.
- Exits every 64 read bases:
  - "may align" once the LCS so far is above the largest refused;
  - "refused" once no LCS can get there: the bit-parallel state holds the LCS of the read so far with every prefix of
    the window, and the r read bases left can only match the window's last bases, so the LCS is at most
    `ZerosBelow(m − r) + r`. This exit is tighter than "the LCS so far + r", though it cut the bound's instructions by
    only 8% on the bench sample.
- The step adds with carry (`_addcarry_u64` on x86-64), and the function is cloned for x86-64-v3 as the aligner's hot
  functions are.
- The lengths alone decide when the bound refuses nothing, or when even LCS = min(n, m) is refused.
- Windows of more than 2^27 read × window bases (genes over ~11 kb) are left to WFA2. The bound's cost grows with the
  product, WFA2's with the budget.

**Counts and timers.**
- The per-sample counts line now reads "T candidate alignments tried: S refused by the k-mer screen, I by the indel
  bound, F aligned from the anchor's exact matches and W as whole windows; ...". For ONT reads S is 0; for other reads I
  is.
- `scripts/measure_performance.sh` reads both forms; its `runs.tsv` has a last column, `indel_refused` (NA for an older
  protal).
- The long-read `<sample>_runtime.tsv` has an "Indel bound" row.
- The options summary has an "indel bound:" line for runs with long reads.

**Tests** (WSL on 4 cores):
- Unit (479 passed): `IndelBound.TheLcsAndTheAnswerAreTheDefinitions` (the LCS against dynamic programming at the word
  boundaries, with Ns, ambiguity codes and lower case; the answer, exits included, against its definition from packed
  bytes and from text, 1,500 candidates with free ends at both ends of both sequences, 417 refused);
  `IndelBound.RefusesOnlyWhatWFA2Fails` (ONT-like windows against WFA2 without X-drop at 0.85: of the gene itself, 0
  refused; of a relative 10-25% away, 183 of 300 refused, all failing; unrelated, 298 of 300 refused, where the k-mer
  screen refused 4); `IndelBound.TheHandlerAlignsLongReadsTheSameWithAndWithoutIt` (through `AlignAnchor` with the
  chain, either strand, anchored and whole-window, the bound in place of the screen against the screen: the same outcome,
  score, start and CIGAR for 439 candidates; the bound refused 209, the screen 0).
- End to end (149 OK): `OntTest.test_the_indel_bound_changes_no_alignment` (with `--no_indel_bound` every SAM record is
  the same, and the counts move from refused to aligned).
- The bench samples (`~/bench071`, V075 database, 90 Mb each, 4 threads): with the bound and with the k-mer screen
  instead (`--no_indel_bound`), the SAM records (ONT 34,871, PacBio 24,947) and every other output file are identical.
  The bound refused 27,180 of the ONT sample's 145,184 candidates, the k-mer screen 0.

**Cost and saving** (callgrind, the first 2,000 reads of each bench sample, 1 thread; instructions, not wall time, as
the machine was shared):

| | without the bound | with it | change |
|---|---|---|---|
| ONT: candidates refused | 0 (the k-mer screen) | 3,867 of 20,889 (18.5%) | |
| ONT: `AlignAnchor` | 35.51 G | 19.55 G | **−45%** |
| ONT: the whole run (index load included) | 53.67 G | 37.71 G | −30% |
| ONT: the bound itself (`LcsAbove`) | | 7.26 G | ~350k per candidate |
| PacBio, the bound after the screen on every candidate it passed (a first version) | 38.93 G | 46.71 G | +20% |
| PacBio, the k-mer screen only (as committed) | 38.93 G | 38.93 G | 0 |

- A refused ONT candidate saved ~6M instructions of WFA2, and the bound costs ~1/18 of that on every candidate it
  checks. So it pays wherever it refuses more than ~5% of the candidates.
- In the bench world 60% of ONT candidates align, and the bound refuses 18.5%. Wall and CPU times of the 90 Mb runs
  were too noisy to compare on the shared machine (between 21 and 37 s of aligning CPU with the bound, 29-37 s
  without); the instructions are not.
- At r226, 99.8% of the host sample's candidates fail and most windows are unrelated, which the bound refuses nearly
  always (298 of 300 in the unit test). Section 3's −60-80% of ONT aligning stands as the estimate. The refusal rate
  there is the number to read.
- The cost per candidate is ~2× section 3's estimate: ~11 instructions per word step, and the "may align" exit of a
  true candidate comes ~70% into the read.

What the next r226 build should show: the new counter on `sc_host_ont` and `sc_gut_ont` in
`logs/protal_runs_*.log` (section 3 expects most of the 17.9M host candidates refused), and the ONT samples' `aligning`
rows in `logs/protal_cpu_<collection>.tsv` against v19's 126.5 s for the host sample.
