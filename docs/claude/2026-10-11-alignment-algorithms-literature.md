# Faster optimal pairwise alignment than WFA2? A literature search

- **Date**: 2026-10-11.
- **Code**: branch `audit-fixes` at `7d5df4b` (`src/Alignment/WFA2Wrapper2.h`, `src/Alignment/AnchoredAlignment.h`,
  `lib/wfa2-lib` = WFA2-lib v2.3.6).
- **Method**: a web literature search (papers, preprints, repositories and their commit logs, up to 2026-10-11).
  Nothing was built or run. The speed figures below are the authors' own measurements, on their data and machines.
- **Question**: are there new optimal or near-optimal pairwise alignment algorithms that are as fast as WFA2, or
  faster?

## What protal asks of the aligner

The answer depends on the use case, and protal's is narrow:

- **Cost model**: gap-affine, match 0, mismatch 4, gap open 6, gap extension 2 (`AnchoredAlignment.h:141`).
- **Alignment type**: ends-free with the read's ends free (flanks, gene windows), plus end-to-end pieces between the
  anchors of long reads. Full CIGAR (`MemoryHigh`), wf-adaptive on, X-drop for short reads.
- **Budget**: a score limit from the ANI floor (`setMaxAlignmentSteps`). A candidate that would cost more is
  abandoned, so failing candidates cost a whole budget each.
- **Sequences**: short-read flanks of a few to ~150 bases at ≤ 10% divergence; long-read gene windows of 1-8 kb at
  ≤ 15% divergence (ONT).

**Where the WFA2 time goes** (earlier reports):

- **Short reads**: WFA2 is 4.93 of the 16.05 G instructions of `RunPairedEnd` (31%), almost all of it flanks
  ([2026-10-06](2026-10-06-performance-profiling/README.md)).
- **Long reads**: at GTDB scale WFA2 is 91% of the alignment stage
  ([2026-10-04](2026-10-04-performance-gtdb-scale/README.md)).
- **ONT in the build**: ONT is 70% of the build's aligning, at 0.59 ms per candidate. On host reads 99.8% of the
  candidates fail ([2026-10-09](2026-10-09-ont-wfa2-skipping.md)).

## Summary

1. **No newer algorithm replaces WFA2 for protal's case.**
   - The exact aligners that beat WFA2 (A\*PA2, QuickEd) do so on long, noisy, *global* alignments with *unit edit
     costs*.
   - On short or similar sequences their own papers report WFA as faster: QuickEd is 14× slower than WFA on 250 bp
     Illumina pairs and 6.8× slower on HiFi; A\*PA2 is slower than BiWFA below 1 kb or 2% divergence.
   - None of them has gap-affine costs or ends-free alignment.
2. **One new method speeds up WFA2 itself, and gives the same optimal scores: Singletrack** (López-Villellas,
   ..., Marco-Sola; bioRxiv, November 2025).
   - **What it does**: it stores only WFA's M wavefronts and recomputes the gap steps during the backtrace. The I and
     D wavefronts become a small recycled ring.
   - **Speed and memory**: WFA2-lib alignments get 1.2-2.1× faster, with 3× less memory for gap-affine.
   - **Benchmark penalties**: the paper benchmarks exactly protal's penalties (0, 4, 6, 2).
   - **Where to get it**: ends-free support (with adaptive, drop and band heuristics) is in a maintained WFA2-lib fork
     (`tmokveld/WFA2-lib`, mode `wavefront_memory_singletrack`).
   - **Catch**: its CIGARs can break ties between equally good alignments differently from `MemoryHigh`. That would
     change SAM records, so it needs an identity check, and a retrain if records change.
3. **For the largest WFA2 cost, failing ONT candidates, a faster aligner is the wrong lever.** Even a 2× aligner
   saves less than the exact indel-distance bound of the [ONT report](2026-10-09-ont-wfa2-skipping.md) (estimated
   −60-80%, outputs identical). The literature supports that route: bit-parallel LCS or edit distance as a
   lower bound, Sassy's SIMD semi-global search.
4. **The rest is not applicable here**:
   - heuristic aligners: Block Aligner, KSW2, miniwfa's heuristic mode;
   - SIMD variants of WFA or BiWFA, which do not beat WFA2's `MemoryHigh` in published numbers;
   - GPU libraries: WFA-GPU, TSUNAMI, Accelign, AGAThA;
   - FORAlign, which only wins at low similarity on 5-190 kb sequences;
   - sequence-to-graph aligners (Theseus).

| Rank | Option | Outputs | Expected effect in protal | Effort | Confidence |
|---|---|---|---|---|---|
| 1 | Exact indel-distance bound before WFA2 (already proposed in the [ONT report](2026-10-09-ont-wfa2-skipping.md)) | identical | ONT aligning −60-80% | ~1 day | medium |
| 2 | **Singletrack** (`tmokveld/WFA2-lib`, or port `wavefront_memory_singletrack` into `lib/wfa2-lib`) | same scores; CIGAR ties may differ | WFA2 time ×0.7-0.8 if the paper's gains carry over: short reads up to ~9% of `RunPairedEnd`, long-read aligning ~15-25% | ~1 day to build and measure | low-medium |
| 3 | Watch: A\*PA2 with affine or semi-global costs (listed as future work), FILTR (a DP-recurrence compiler, July 2026) | – | none yet | – | – |

## 1. Singletrack: the only drop-in speed-up of WFA2 found

**Paper**: L. López-Villellas, C. Iñiguez, A. Jiménez-Blanco, Q. Aguado-Puig, M. Moretó, J. Alastruey-Benedé,
P. Ibáñez, S. Marco-Sola. *Singletrack: An Algorithm for Improving Memory Consumption and Performance of Gap-Affine
Sequence Alignment.* bioRxiv 10.1101/2025.10.31.685625 (2 November 2025). Code: `github.com/LorienLV/singletrack`.

**Idea.** In the gap-affine recurrence each I or D cell traces back either to an M cell (a gap opened) or to the
same I or D matrix (a gap extended). So the backtrace can follow M alone:
- If the current M cell is not explained by a diagonal step, follow two tentative gap paths (I and D).
- Lengthen them until one meets an M value with the expected score. This costs O(n + m), at most about twice the
  classical backtrace.
- In WFA this means only the M̃ wavefronts are kept for the backtrace. The Ĩ and D̃ wavefronts live in a small scope
  that is recycled as the score advances.
- The time saved comes from better cache behaviour, not from fewer cells.

**Results** (single thread, global alignment, penalties a = 0, x = 4, o = 6, e = 2 as in protal; 2× Xeon Platinum
8480+). Times are µs per kbp of alignment, memory is peak RSS:

| Dataset | WFA (WFA2-lib) | WFA+Singletrack | Speed-up | BiWFA | KSW2 | WFA mem | WFA+ST mem |
|---|---|---|---|---|---|---|---|
| Illumina, 140-280 bp, ~1% error (GIAB, 100M pairs) | 4.8 | 3.4 | 1.41× | 8.9 | 65.8 | 5.3 MiB | 26.0 MiB |
| PacBio HiFi, 0.2-24.6 kb, 0.6% (10M pairs) | 13.3 | 10.9 | 1.22× | 18.6 | 2,980.6 | 428 MiB | 180 MiB |
| ONT PromethION, ≤ 50 kb, 12% (18 pairs) | 7,447.2 | 5,233.2 | 1.42× | 4,428.1 | 7,849.8 | 6,337 MiB | 2,124 MiB |

With dual gap-affine costs the gain is 2.1× on Illumina. Over many threads WFA+ST scales better than WFA, because it
is less bound by shared cache and memory: up to 3.2× over WFA at equal thread counts.

**What carries over to protal, and what may not:**
- **Penalties**: identical, so the score per base and the wavefront sizes are like the paper's.
- **Ends-free**: the paper benchmarks global alignment only. The `tmokveld/WFA2-lib` fork, which merges upstream
  regularly (last on 2026-10-05), supports Singletrack for end-to-end, ends-free and extension alignments.
- **Heuristics**: the fork documents Singletrack with the adaptive, drop and band heuristics. It does not document
  `max_alignment_steps` (protal's budget), so that needs a test.
- **Flank lengths**: protal's short-read flanks are shorter than the paper's 140-280 bp pairs, and wf-adaptive
  already keeps the wavefronts small. The gain from less memory traffic may be smaller than 1.41×.
- **Long failing candidates**: ONT's failing candidates run to scores of 600-7,700 under `MemoryHigh` (O(s²)
  wavefronts). There the 3× smaller memory could matter more. That is a guess, and rank 1 removes most of those
  candidates anyway.
- **Ties**: the fork's README says plainly that exact Singletrack "can resolve ties differently" from
  `wavefront_memory_high`. protal did not switch on upstream's AVX2 kernels for this same reason: alignments on
  another diagonal of the same score ([2026-10-02](2026-10-02-wfa-avx2/README.md)). The CIGAR feeds the SAM records,
  the strain alleles and the features. Changed ties are acceptable only with a retrain.
- **Memory per aligner**: on short reads Singletrack allocates more (26 vs 5.3 MiB in the paper; it reserves the
  worst-case scope). protal keeps one aligner per thread, so this is small.

**How to measure it.**
1. Swap `lib/wfa2-lib` for the fork, or port the memory mode.
2. Add a `MemorySingletrack` choice in `WFA2Wrapper2`.
3. Run the usual pe100k, pb and ont callgrind runs (`2026-10-06-performance-profiling/scripts/cc.sh`).
4. Compare the SAMs record by record with the current build.

If the records are identical, take it. If they differ, count how often and how (indel placement only, or other
mismatches too), then decide together with the next retrain.

## 2. Exact aligners that are faster than WFA2 elsewhere

**A\*PA2** (R. Groot Koerkamp, WABI 2024, LIPIcs 312:17). Edit distance only, global only, in Rust (C bindings in
`astarpa-c`, MPL-2.0).
- **Method**: Ukkonen band doubling with Myers' bit-parallel DP, with these additions:
  - 256-column blocks and 256-bit SIMD;
  - a 2-bit profile, so the alphabet is ACGT only;
  - incremental doubling;
  - an optimistic diagonal-transition traceback;
  - optionally A\*PA's gap-chaining seed heuristic.
- **Speed**: 19× faster than BiWFA and Edlib on > 500 kb ONT reads (6% divergence), 5.6× on 11 kb reads, and 0.81×
  on 800 bp reads (11%).
- **Where BiWFA wins**: the paper's own discussion says BiWFA is faster for sequences under 1 kb or under 2%
  divergence, which is protal's whole range.
- **Future work listed**: affine costs (which would need SIMD affine blocks in place of bitpacking) and semi-global
  alignment. The repository still describes A\*PA as edit distance only.

**QuickEd** (M. Doblas et al., Bioinformatics 2025, btaf112). Edit distance only, end-to-end, C library.
- **Method**: "bound-and-align". Fast heuristics give an upper bound ŝ, then a bit-parallel banded DP runs in
  O(nŝ).
- **Real datasets**, single thread:

| Dataset | WFA | Edlib | QuickEd | QuickEd vs WFA |
|---|---|---|---|---|
| Illumina 250 bp | 4.3 s | 186.0 s | 60.6 s | 14× slower |
| PacBio HiFi | 58.8 s | 3,405.6 s | 400.8 s | 6.8× slower |
| ONT PromethION | 353.4 s | 574.2 s | 96.6 s | 3.7× faster |
| ONT UltraLong | 2,824.8 s | 3,858.6 s | 835.2 s | 3.4× faster |

  The authors say WFA wins on short, similar sequences, because QuickEd's smallest band is 192 cells.

**FORAlign** (Wei et al., Briefings in Bioinformatics 2025). Gap-affine, Four Russians lookup blocks inside
Hirschberg's divide and conquer, multithreaded C++.
- **Data**: only 5-190 kb sequences.
- **High similarity**: WFA is far faster. The 190 kb MPox self-pair takes WFA 1.3-1.8 s against FORAlign's 126-217 s
  per batch.
- **Low similarity**: FORAlign wins (16.8× over the paper's WFA2 baseline at 128 threads, where WFA's high-memory
  mode runs out of memory).
- Not applicable.

## 3. Heuristic aligners (not optimal)

- **Block Aligner** (D. Liu and M. Steinegger, Bioinformatics 2023). Adaptive blocks that shift, grow and shrink;
  global, local and X-drop with traceback, no ends-free; Rust with C bindings, SSE2/AVX2/NEON.
  - **Speed**: about 6× faster than exact WFA on erroneous Nanopore reads, but WFA is "much faster" on Illumina
    reads.
  - **Suboptimal alignments**: 0% of Illumina reads, 0.2% of 1 kb Nanopore reads, 2.2-2.8% of longer Nanopore reads.
  - **Against WFA-adaptive** (what protal uses): about 5× fewer suboptimal results.
  - **In A\*PA2's benchmark**: about as fast as BiWFA on the small datasets.
  - Not a gain for short reads. Long reads would trade optimality for a speed protal can mostly get from rank 1.
- **KSW2 and mm2-fast**. KSW2 is minimap2's banded Suzuki-Kasahara DP, and mm2-fast its AVX-512 version (Kalikar et
  al., Nature Computational Science 2022; up to 1.8× end-to-end).
  - **Against WFA**: in the Singletrack table KSW2 is 14× slower than WFA on Illumina and 220× slower on HiFi.
  - **Elsewhere**: Li and Homer's minibwa (arXiv 2606.15357, June 2026) uses minimap2's base alignment, not a new
    algorithm.
- **miniwfa** (H. Li). A reimplementation of WFA for dual gap-affine costs with a low-memory mode and a chaining
  heuristic. In its exact mode its README says BiWFA is generally faster. Only its heuristic mode is faster, and it
  is not optimal.
- **Scrooge / GenASM** (Lindegger et al., Bioinformatics 2023). Bit-parallel edit distance (Bitap-based), windowed,
  so approximate. On the CPU 1.7× faster than Edlib and 20× faster than KSW2 on long reads. Edit distance only.

## 4. Faster WFA kernels, and other hardware

- **WFA-vect** (IEEE Xplore 10822352): hand-vectorised WFA with load-synchronous and mask-based branching. It reports
  3.87-3.98× over *scalar* WFA at 1% and 20% error. The full text was not accessible. WFA2-lib's kernels are already
  auto-vectorised and its baseline is unclear, so the numbers say little about WFA2.
- **Bi-WaveX** (Chan et al., WCTP 2025, published April 2026): AVX-512 BiWFA, up to 4.67× over BiWFA. BiWFA is the
  O(s)-memory mode, which is slower than `MemoryHigh` on short sequences (8.9 vs 4.8 µs/kbp in the Singletrack
  table). No evidence that it beats protal's mode.
- **protal's own test**: AVX2 extend kernels for WFA2 gave no gain on short or long reads
  ([2026-10-02](2026-10-02-wfa-avx2/README.md), [2026-10-04 §4](2026-10-04-performance-gtdb-scale/README.md)).
- **GPUs**: WFA-GPU (Aguado-Puig et al., Bioinformatics 2023), TSUNAMI (IEEE 2023), Accelign (Kallenborn, Dabbaghie,
  Steinegger, Schmidt; bioRxiv December 2025, 9-16 TCUPS on one Blackwell GPU), AGAThA (PPoPP 2024). protal runs on
  CPU nodes. Out of scope.
- **Compilers**: FILTR (Vinaithirthan, Sundram, Goenka, Kjolstad; arXiv 2607.06225, July 2026) compiles DP
  recurrences with pruning and scheduling to C++, "0.95× to 30×" against hand-tuned code. It is too new to tell
  whether it matches WFA2 on gap-affine ends-free alignment. Worth a look later.

## 5. Skipping alignments: what the literature adds to rank 1

The ONT report's bound works like this:
- Each base of an alignment costs at least 2 per indel-distance unit: a mismatch costs 4 and is 2 indels, a gap base
  costs ≥ 2.
- So the affine score is ≥ 2 × the (ends-free) indel distance, and a candidate whose indel distance exceeds half the
  budget cannot succeed.
- The same holds for unit edit distance (score ≥ 2 × edit distance), but that bound is weaker, because indel distance
  ≥ edit distance.

Tools that compute such distances fast:
- **Myers (1999) / Hyyrö**: bit-parallel edit distance and LCS, one 64-bit word per 64 rows and a handful of
  instructions per column. Edlib uses Myers' algorithm with a threshold k and stops once k is exceeded, which is the
  shape of a budgeted reject.
- **Sassy** (R. Beeloo and R. Groot Koerkamp, Bioinformatics May 2026; Sassy2, bioRxiv March 2026): semi-global search
  of a pattern of up to ~1,000 bp in a text with ≤ k edits. SIMD bitvectors in the text direction, 4-15× faster than
  Edlib, CIGARs out, C bindings.
  - **Drawbacks for protal**: it needs AVX2 and BMI2 at compile time (`x86-64-v3`), which clashes with protal's one
    binary with run-time dispatch. It is Rust, and it is unit edit distance.
  - **Use**: a reference design rather than a dependency. A plain bit-parallel LCS in C++ (~100 lines) fits protal
    better.
- **Pre-alignment filters** (SneakySnake, Shouji, GateKeeper; Alser et al., 2017-2020) decide "within e edits or not"
  for short reads. SneakySnake has no false rejects at a given edit threshold. They are built for Hamming-like
  short-read windows. protal's q-gram screen (`AlignmentScreen.h`) already fills that role for short reads.

## 6. Upstream WFA2-lib since v2.3.6

- **Tags**: v2.3.6 is the last tag (2026-03-18) and is what protal vendors.
- **Commits since** (to 2026-09-30): fixes only, no speed work.
  - static-band fixes, also for BiWFA;
  - a SAM-CIGAR fix (protal reads the operations, not WFA2's SAM CIGAR);
  - a CIGAR buffer length fix;
  - [#126](https://github.com/smarco/WFA2-lib/issues/126): an end-to-end alignment on a reused aligner inherited the
    previous ends-free span, which gave wrong scores (fix `90582eb`).
- **#126 and protal**: by code reading protal is not affected. Its end-to-end pieces use their own aligner, which never
  aligns ends-free (`m_pieces`, `AnchoredAlignment.h:290-293`). The window aligner always aligns ends-free with
  explicit spans.
- **Singletrack upstream**: Cristian Iñiguez (a Singletrack author) has merged several PRs upstream in 2026, but
  Singletrack itself is so far only in the fork and the authors' repository.

## Follow-up, later on 2026-10-11

- **Rank 1 was already in protal.** The exact indel-distance bound for ONT candidates has been protal's default since
  2026-10-09 (`AlignmentScreen.h` IndelBound; `--no_indel_bound` turns it off). This report called it only proposed.
- **Rank 2 was tried**: [Singletrack in protal](2026-10-11-singletrack/README.md) (branch `singletrack`).
  - The fork needed three fixes. As it came, its Singletrack backtrace ended a PacBio run (`exit(-1)`), and its classic
    traceback gave invalid CIGARs after wf-adaptive pruning.
  - With the fixes, classic is v2.3.6 byte for byte, and Singletrack gives the same scores.
  - Speed: PacBio aligning −13% (handler −23%), ONT −10%, short reads no change.
  - Other records: 0.01% of short reads, 0.2-0.6% of PacBio reads, 47-49% of ONT reads (tie-breaks).
  - It is available as `--wfa_traceback`, recorded per database in `build_settings.tsv`.

## Sources

- Singletrack: [bioRxiv 10.1101/2025.10.31.685625](https://www.biorxiv.org/content/10.1101/2025.10.31.685625.full.pdf),
  [LorienLV/singletrack](https://github.com/LorienLV/singletrack), [tmokveld/WFA2-lib](https://github.com/tmokveld/WFA2-lib)
- WFA2-lib: [smarco/WFA2-lib](https://github.com/smarco/WFA2-lib), [tags](https://github.com/smarco/WFA2-lib/tags),
  [issue #126](https://github.com/smarco/WFA2-lib/issues/126), [commit 90582eb](https://github.com/smarco/WFA2-lib/commit/90582eb)
- A\*PA2: [WABI 2024](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.WABI.2024.17),
  [astar-pairwise-aligner](https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner)
- QuickEd: [Bioinformatics 2025 (PMC11937955)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11937955/)
- FORAlign: [Briefings in Bioinformatics 26(1) bbaf061](https://academic.oup.com/bib/article/26/1/bbaf061/8030539)
- Block Aligner: [Bioinformatics 2023 (PMC10457662)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10457662/)
- BiWFA: [Bioinformatics 39(2) btad074](https://academic.oup.com/bioinformatics/article/39/2/btad074/7030690)
- Bi-WaveX: [WCTP 2025](https://www.atlantis-press.com/proceedings/wctp-25/126023801)
- WFA-vect: [IEEE Xplore 10822352](https://ieeexplore.ieee.org/document/10822352/) (abstract via search only)
- Sassy: [bioRxiv 2025.07.22.666207](https://www.biorxiv.org/content/10.1101/2025.07.22.666207.full.pdf),
  [docs.rs/sassy](https://docs.rs/crate/sassy/latest)
- Scrooge: [arXiv 2208.09985](https://arxiv.org/pdf/2208.09985)
- mm2-fast: [bioRxiv 2021.07.21.453294](https://www.biorxiv.org/content/10.1101/2021.07.21.453294v2.full.pdf);
  minibwa: [arXiv 2606.15357](https://arxiv.org/abs/2606.15357); miniwfa: [lh3/miniwfa](https://github.com/lh3/miniwfa)
- WFA-GPU: [Bioinformatics 39(12) btad701](https://academic.oup.com/bioinformatics/article/39/12/btad701/7425447);
  Accelign: [bioRxiv 2025.12.17.694868](https://www.biorxiv.org/content/10.64898/2025.12.17.694868.full.pdf);
  AGAThA: [arXiv 2403.06478](https://arxiv.org/abs/2403.06478)
- FILTR: [arXiv 2607.06225](https://arxiv.org/abs/2607.06225); Theseus: [bioRxiv 2026.02.12.705572](https://www.biorxiv.org/content/10.64898/2026.02.12.705572v1.full)
- SneakySnake: [arXiv 1910.09020](https://arxiv.org/pdf/1910.09020)
