# The r226 v7 and v8 builds: the divergence features, the suspect gene copies and the sample's depth at GTDB scale

2026-10-03. Data: two SLURM runs of `build_gtdb_database.py` on GTDB r226 with the binary and scripts of `a3e397d`
(protal v0.7.4: round one of the false-positive work merged, [report](../2026-10-03-false-positive-fixes/README.md)),
uploaded as `local/protal_r226_v7_{logs,training}.tgz` and `local/protal_r226_v8_{logs,training}.tgz` (git-ignored;
extracted to `local/v7` and `local/v8`). Both runs: `--inputs refGenos --scratch` on the node's SSD, 12 samples per
design point, read pairs 1,000 to 10M, abundances σ 1.3, 20% of the species (7,899) and the clades held out, seed 1,
so that both have the same samples (56,737 held-out rows, 23,629 test taxa in 78 paired-end samples; the test design
σ 2.0, 500 to 5M pairs). They differ in the features and the calls:

| run | features | calls | threads |
|---|---|---|---|
| v7 (job 23904738) | `normalized+adjacency+distance+divergence` | knob curve over the sample's depth | 64 |
| v8 (job 23904739) | the default of `a3e397d`: `normalized+adjacency+distance+depth+divergence` | knob 0.5 (no curve with the depth feature) | 52 |
| v5 (`27423c6`, [report](../2026-10-03-r226-v5-v6-training/README.md)), for reference | `normalized+adjacency+distance` | curve (fdr shipped) | 64 |

v7 against v5 (other samples: the simulator changed between the commits) shows the divergence features and the
suspect-copy scan; v8 against v7 (the same samples) is a clean ablation of the depth feature at GTDB scale. The
per-taxon comparison is [`compare_v7_v8.py`](compare_v7_v8.py) → [`compare_v7_v8_output.txt`](compare_v7_v8_output.txt),
the suspect copies [`incongruence_r226.py`](incongruence_r226.py) → [`incongruence_r226_output.txt`](incongruence_r226_output.txt);
the rest is read off `model_logs/summary.txt`, `build_metadata.tsv` and the training reports.

## Summary

Independent test set, F1 at knob 0.5 and, where a curve was fitted, at the knob curve (how that protal calls by
default); false positives per test sample; F1 of the species held out at knob 0.5:

| read type | v5 test 0.5 / curve | v7 test 0.5 / curve | v8 test 0.5 | FP per sample v7 → v8 | held out v7 → v8 |
|---|---|---|---|---|---|
| pe | 0.9615 / 0.9642 | 0.9634 / 0.9654 | **0.9671** | 2.50 → 1.64 | 0.9644 → 0.9747 |
| se | 0.9570 / 0.9593 | 0.9575 / 0.9606 | **0.9685** | 2.91 → 1.46 | 0.9587 → 0.9719 |
| pb | 0.9714 / 0.9761 | 0.9720 / 0.9764 | 0.9763 | 0.93 → 0.64 | 0.9728 → 0.9743 |
| ont | 0.9704 / 0.9696 | 0.9717 / 0.9737 | 0.9730 | 1.00 → 0.78 | 0.9723 → 0.9751 |

- **The depth feature holds at GTDB scale.** Against v7's knob curve, v8 gains 0.0017 (pe) and 0.0079 (se) of test F1
  and is level for the long reads (pb −0.0001, ont −0.0007, within the noise of 3,000-3,800 test taxa); against knob
  0.5 it gains 0.0037, 0.0110, 0.0043 and 0.0013. False positives per sample fall by a third to a half for every
  read type. The 0.7.4 default, kept on the benchmark world against a loss there, was the right call.
- **The divergence features and the suspect-copy scan are worth a little**: v7 against v5, +0.001 to +0.002 for
  every read type at both operating points, false positives down (pe test 218 → 195). Not separable from each other
  or from the samples' change here; `excess_scaled_median` is the most important feature of the pe model in both v7
  and v8.
- **The suspect-copy scan at r226**: 5,589 of 14.5M gene copies (0.04%) of 2,908 species, 1,530 of them identical to
  another genus's copy; two thirds of the species lose one gene, ten lose 21-32 (MAG contamination between gut
  genera, see below). 8-9 minutes of the build.
- **The error budget has tipped from false positives to misses.** v8's pe test set has 128 false positives and 173
  misses (v7: 195 and 143); the species held out 330 and 255. The remaining false positives are nearly all congeners
  of a species the database lacks (104 of 128); the misses are strains with 1-10 fragments (123 of 173), most of them
  in the test samples of 50,000 and 200,000 read pairs, where the depth prior learned on the training design (σ 1.3,
  20-200 species) under-calls the thin present taxa of the test design (σ 2.0, 10-300 species). Round two's mixed
  training σ and its unfiltered features aim at exactly these two populations; the next build measures them.

## The depth feature: what it trades

### The test set by depth (pe; se alike)

| read pairs | present | absent | FP v7 | FP v8 | FN v7 | FN v8 |
|---:|---:|---:|---:|---:|---:|---:|
| 500 | 53 | 13 | 0 | 1 | 5 | 3 |
| 2,000 | 177 | 56 | 6 | 11 | 6 | 4 |
| 10,000 | 558 | 229 | 10 | 17 | 28 | 17 |
| 50,000 | 1,164 | 804 | 48 | 40 | 61 | **80** |
| 200,000 | 800 | 3,032 | 31 | 15 | 20 | **43** |
| 1,000,000 | 1,160 | 6,168 | 66 | 33 | 17 | 20 |
| 5,000,000 | 679 | 8,736 | 34 | 11 | 6 | 6 |

The forest has learned the knob curve and more: at the shallow points it calls more (fewer misses, a few more false
positives, net positive), at the deep points it rejects the singletons (at 1M and 5M pairs the false positives fall
to half and a third). The cost sits at 50,000 and 200,000 pairs, where the misses rise from 81 to 123.

### Which calls flipped (pe test set, v7 → v8 at knob 0.5)

| flip | taxa | ≤ 2 fragments | strains (another genome than the representative) | read pairs |
|---|---:|---:|---:|---|
| false positive removed | 89 | 79 | 0 | 50k: 10, 200k: 18, 1M: 37, 5M: 24 |
| true positive lost | 51 | 44 | 38 | 50k: 24, 200k: 23, 1M: 4 |
| false positive added | 22 | 12 | 0 | 500 to 10k: 13, deeper: 9 |
| true positive gained | 21 | 17 | 16 | 500 to 10k: 15, deeper: 6 |

Both populations are taxa of one or two fragments at identity ~0.98: the depth feature cannot tell a present strain's
single read from an absent congener's single read, it only knows that at 200,000 pairs most singletons are absent.
That prior is right for the training design and too strict for the test design, whose σ 2.0 abundances give every
sample a longer tail of thin present species. On the species held out (the training design) the same feature loses
nothing: 245 false positives removed against 55 added, 76 true positives gained against 22 lost, and the misses of
one-fragment taxa fall from 153 to 117 of 1,101.

### The species held out: singletons by depth

Absent taxa with one fragment, and how many were called:

| read pairs | absent, 1 fragment | FP v7 | FP v8 |
|---:|---:|---:|---:|
| 1,000 | 73 | 5 | 10 |
| 5,000 | 328 | 23 | 39 |
| 20,000 | 985 | 28 | 33 |
| 100,000 | 3,138 | 77 | 27 |
| 500,000 | 9,354 | 99 | **1** |
| 2,000,000 | 5,179 | 34 | **0** |
| 10,000,000 | 6,851 | 9 | **0** |

From 500,000 pairs on, no single-fragment taxon is called any more, present or not; the 52% of r226 false positives
that were perfect singletons in deep samples ([anatomy](../2026-10-03-false-positive-anatomy/README.md)) are gone.
The price is the present singletons there, which no feature of the record itself can rescue; what could is what
round two added: the taxon's reads before the MAPQ filter and the EM's share of them (`em_fragments`), the reads
that seeded on it and failed (`failed_candidate_rate`), and GTDB's priors on the species.

### Calibration

Best threshold on the test set and the F1 there (at 0.5 in brackets):

| read type | v7 held out | v7 test | v8 held out | v8 test |
|---|---|---|---|---|
| pe | 0.634 | 0.519 (0.9637 vs 0.9634) | 0.545 | 0.468 (0.9682 vs 0.9671) |
| se | — | 0.563 (0.9584 vs 0.9575) | — | 0.421 (0.9696 vs 0.9685) |
| pb | 0.401 | 0.310 (0.9756 vs 0.9720) | 0.392 | 0.372 (0.9774 vs 0.9763) |
| ont | 0.449 | 0.557 (0.9719 vs 0.9717) | 0.450 | 0.356 (0.9753 vs 0.9730) |

The short-read models are within 0.001 of their best at knob 0.5 on both sets; v8's best test thresholds below 0.5
(0.47, 0.42) say the same as the depth table: the model is a little too strict for the test design. The PacBio model
is consistently best below 0.5 (0.39-0.40 held out, 0.31-0.37 on the test set, in both runs), for +0.001 held out
and +0.001 to +0.004 on the test set; ONT's best thresholds scatter (0.45, 0.56, 0.36). The trainer's rule for knob
points (a gain of 0.002 on the species held out) would not move either, and the long-read test sets have shown
before that their thresholds do not generalize ([v4/v2 reruns](../2026-10-03-r226-v4-v2-rerun/README.md)).

### Held out against test

v7's F1 held out and on the test set were 0.001 apart (0.9644, 0.9634); v8's are 0.0076 apart (0.9747, 0.9671).
The depth feature fits the training design's prior, and the test set, by design, has another one. This gap is the
number the mixed training σ of round two (`--abundance lognormal:1.3,2.0`, with 30% of the species held out) should
close; on the benchmark world it took the depth feature's cost from −0.005 to within noise.

## The suspect gene copies at r226

`protal --build` compared 17.4 billion sketch pairs in 8 minutes (9.5 at 52 threads) and found 1,400,361 pairs of
gene copies within 0.05 of each other across genera (1,387,438 within one family: genus boundaries in GTDB are not an
ANI threshold), 1,682 of them identical. Flagged, by the rule of the round-one report (within 0.02 of another genus's
copy and farther from one's own congeners' by 0.02): 5,589 copies of 2,908 species.

| | |
|---|---|
| by the rank the pair shares | family 2,731, domain 1,010, order 856, class 730, phylum 257, none 5 |
| by distance to the other genus's copy | identical 1,530, ≤ 0.01: 2,035, ≤ 0.02: 2,024 |
| species by genes flagged | one gene 1,938, two 492, 3-5: 336, 6-20: 132, 21-32: 10 |
| genes | 156 of 168 have a flagged copy; the most, 154 copies each (genes 2 and 78); at most 1% of a gene's copies |

The ten species with 21-32 genes flagged are MAG pairs across gut genera, most of their copies identical:
*Dorea_A sp019421265* ↔ *Ruminococcus_B gnavus* (32 and 23 genes, 21 identical), *Phocaeicola sp948566455* ↔
*Bacteroides intestinalis* (26), *Oliverpabstia tarda* ↔ *Blautia_A fusiformis* (26), *JALHMQ01 sp020248785* ↔
*Rubrivivax sp020248915* (24 each, 14 identical, flagged on both sides), *Caccovivens sp017464845* ↔
*Avigastranaerophilus sp017467305* (14 each, 12 identical). A fifth of their markers gone leaves them detectable; a
single read on such a gene no longer creates both species. The 1,010 copies whose pair shares only the domain (a
bacterial gene copy within 0.02 of an archaeal one, or the reverse) are contamination by definition. Nothing in the
scan's output suggests the thresholds are off; the near-pair table (`model_logs/gene_incongruence.tsv`, 1.4M rows)
is the place to look if a species of interest goes missing.

## Does this change the strategy?

No; it confirms it and sharpens what to measure in the next build (`b305b71`, round two: `fragments_all`,
`em_fragments`, `failed_candidate_rate`, the species priors, σ 1.3/2.0, 30% of the species held out, a 30M-pair
point):

1. **Keep the depth feature as the default**, and the calls at knob 0.5 with no curve. Its gain at r226 (pe +0.002 to
   +0.004, se +0.008 to +0.011, false positives −35 to −50%) is the largest single step since the distance features.
2. **The next build's numbers to watch**: the test set's misses at 50,000 and 200,000 pairs (80 and 43 for pe in v8),
   the gap between held-out and test F1 (0.0076), and the test set's best threshold (0.468, 0.421). The mixed σ should
   move all three toward v7's values while keeping v8's false positives. If the misses at those depths do not fall,
   the training design's species-per-sample range (20-200 against the test's 10-300) is the next prior to mix.
3. **The misses are now the larger error**, and they are strains with 1-10 fragments at identity 0.94-0.97, the
   population round two's unfiltered features address from the records (half the r226 misses had most of their records
   below MAPQ 10, [FN anatomy](../2026-10-03-false-positive-fixes/README.md)). Should they not move it, the next
   lever is the index itself: strain alleles of the marker genes beside the representative's, so that a strain's reads
   align at identity 1 and are not ambiguous (the F1 strategy list's item 1, ceiling +0.009 pe).
4. **The remaining false positives are congeners of missing species** (104 of 128 on the test set), the class no
   feature of the taxon can resolve: the reads are real, the species is not in the database. The priors' cluster radius
   can help where the genus is crowded; beyond that only a report of such calls at genus level would turn them from
   errors into information (tried as a second forest: +0.0002, [genus fallback](../2026-10-03-genus-fallback/README.md)).
5. **No change to the suspect-copy scan**: 0.04% of copies, the right species, 8 minutes. A flat knob below 0.5 for
   the PacBio model is the only setting the data argue for, and only weakly (+0.001 held out); not worth a special case
   until a second build repeats it.

The incongruence report printed identical copies' distance as `-0` (the sign of `-log(1)`); `SketchDistance` now
clamps at 0 (no value changes, [`GeneConservation.h`](../../../src/SequenceUtils/GeneConservation.h)).
