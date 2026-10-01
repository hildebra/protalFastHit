# Depth-scaled identity, genus-only unsure calls, and other alignment-level evidence

- **Date**: 2026-10-01.
- **Question**: do (a) a read-identity threshold scaled by depth and (b) reporting unsure calls (p < 0.8)
  only as their genus when their reads fit a congener as well raise the F1 of the presence calls; and
  which other alignment-level information helps. Follow-up of
  [False positives traced to their reads](../2026-10-01-false-positives-v2-test/README.md).
- **Data**: the `V2` build of [Models per read type and a test set](../2026-09-30-read-type-models.md)
  (`~/tune/V2` in WSL): training tables of 96 samples (36 pe, 36 se, 12 pb, 12 ont) and the independent
  test set of 48 (18, 18, 6, 6), profiled against the training database (species and whole clades held out).
- **Code**: the V2 profiles, dumps and SAMs (branch `audit-fixes` at `5385719`, binary
  `~/protal-head/build/protal`, unchanged since). The analysis ran at `34f2201`; nothing in protal changed.
- **Realignment**: every V2 sample aligned again with the same binary and `-m 3` (the best alignment and the
  two next candidates of protal's default `align top 3`, as secondary records), alignment only
  ([`realign.sh`](realign.sh), 2 min for the test set, 4 min for the training set, 6 threads, 3.6 GB). The
  primary records are those of the V2 run (none of 20,000 compared differ).
- **Model**: the trainer's forest, refitted per experiment (64 trees, at most 128 leaves, balanced classes,
  `max_features` sqrt, the 24 normalized features). With seed 1 it reproduces the V2 reports exactly
  (pe: cross-validation 24 FP / 45 FN, test 11 / 58; se 47 / 63 and 14 / 92; pb test 2 / 22; ont 5 / 53).
- **Scoring**: species held out in 5-fold cross-validation of the training table (the trainer's
  `species` scheme), and the test set scored by a forest fitted on all training rows; 5 forest seeds; the
  test F1 difference with a paired bootstrap over samples (500 resamples per seed; the intervals below are
  their mean over seeds). Rules are tuned on the cross-validated training calls only and applied
  unchanged to the test set.
- **Run** (WSL, `protal-db-build` environment):

      bash realign.sh                          # ~/fpexp/{test,training}/sam
      python3 prep_features.py test training   # SAM features of the V2 run
      python3 prep_features.py --secondary test training
      python3 run_features.py 5 --with-secondary > features_output.txt
      python3 run_rules.py 5 base > rules_base_output.txt
      python3 run_rules.py 5 all > rules_all_output.txt
      python3 unsure_calls.py > unsure_output.txt

## The candidates

- **(a) identity scaled by depth.** A present taxon's identity varies between taxa (its strain is 0.2-2%
  from the reference) and, with few reads, by chance. Fitted on the training table's present taxa:
  var(identity | n fragments) = sb² + sr²/n, with mean μ at 20 or more fragments (pe: μ 0.973, sb 0.013,
  sr 0.012; se 0.976, 0.012, 0.011; pb 0.967, 0.010, 0.023; ont 0.956, 0.012, 0.041). `identity_z` =
  (identity − μ) / sqrt(sb² + sr²/n) is how far below a present taxon's identity a taxon is, given its
  depth. Tried as a **veto** (drop a call with `identity_z` < −k, k tuned) and as a **feature**.
- **(b) genus only for unsure calls.** A call with 0.5 ≤ p < U (U = 0.8 as proposed, also tuned) is
  reported as its genus instead of its species when its reads fit a congener as well: the share of its
  reads with an alternative alignment to another species of the genus within δ edits of the best one
  (mismatches, indels and clipped bases, from the realignment) is at least f (as proposed: δ = 1, f = 0.5;
  also tuned). A cheaper variant needs no realignment: a congener is a confident call (p ≥ U) in the same
  sample. A demoted call is no species call (an absent one stops being a false positive, a present one
  becomes a false negative at species level); its genus report is right when a species of that genus was
  simulated in the sample.
- **(c) other alignment evidence**, per taxon, as features: `mean_mapq` (in the dump, not in the model) and
  the share of reads with MAPQ < 10 (protal's MAPQ is the score gap to the next candidate); for pe, the
  share of fragments with both mates on the taxon (`mate_concordance`); the share of reads clipped by 10 or
  more bases; and the congener-fit share of (b) (`congener_share_d1`, with `other_genus_share_d1`).

## (a) Depth-scaled identity: no gain

| | pe | se | pb | ont |
|---|---|---|---|---|
| veto, tuned k | 2.75 | 3.25 | 4.00 | 4.75 |
| test F1 change, veto | −0.0004 | +0.0002 | 0 | −0.0008 |
| CV F1 change, feature | −0.0007 | −0.0004 | +0.0008 | −0.0005 |
| test F1 change, feature (interval) | −0.0004 (−0.0035, 0.0026) | +0.0007 (−0.0020, 0.0032) | −0.0013 (−0.0061, 0.0028) | +0.0016 (−0.0025, 0.0080) |

Tuned on the training calls, the veto settles where it drops almost nothing. The forest already splits on
identity, top identity and fragments together, and the low end of present taxa (strains 2% from their
reference, sequencing errors) reaches the relatives' 0.93-0.94: only the three thick false positives of
the previous report sit there, beside a few present taxa.

## (b) Genus only for unsure calls: right genera, lower species F1

Base model, test set, mean of 5 seeds (all-features model: the same pattern, [`rules_all_output.txt`](rules_all_output.txt)):

| read type | rule | test F1 | FP | FN | demoted: absent / present | genus right | output precision | species sensitivity |
|---|---|---|---|---|---|---|---|---|
| pe | none | 0.9761 | 9.4 | 61.4 | | | 0.9935 | 0.9593 |
| pe | congener fit, as proposed | 0.9741 (−0.0021) | 7.4 | 69.2 | 2.0 / 7.8 | 98% | 0.9948 | 0.9541 |
| pe | confident congener, U 0.8 | 0.9691 (−0.0071) | 5.8 | 85.0 | 3.6 / 23.6 | 100% | 0.9960 | 0.9436 |
| se | none | 0.9642 | 14.0 | 88.6 | | | 0.9900 | 0.9398 |
| se | congener fit, as proposed | 0.9634 (−0.0008) | 9.6 | 95.0 | 4.4 / 6.4 | 100% | 0.9931 | 0.9354 |
| se | confident congener, U 0.8 | 0.9570 (−0.0073) | 8.2 | 113.8 | 5.8 / 25.2 | 100% | 0.9941 | 0.9226 |
| pb | none | 0.9758 | 3.4 | 25.2 | | | 0.9942 | 0.9582 |
| pb | congener fit, as proposed | 0.9737 (−0.0023) | 2.4 | 28.6 | 1.0 / 3.4 | 100% | 0.9959 | 0.9526 |
| ont | none | 0.9515 | 5.4 | 52.2 | | | 0.9905 | 0.9154 |
| ont | congener fit, as proposed | 0.9501 (−0.0015) | 5.4 | 53.8 | 0 / 1.6 | 100% | 0.9905 | 0.9128 |

Output precision: of everything reported (species and genera), the share that is right. Tuning U, δ and f
on the training calls does not help either: every tuned variant loses on the test set (−0.001 to −0.010).

Why: the unsure band holds mostly present species. On the test set, 60 pe calls have 0.5 ≤ p < 0.8, 8 of
them absent and 52 present (se 11 / 46, pb 2 / 69, ont 2 / 35), although that band holds most false
positives (73%, 79%, 100% and 40% of them). The congener-fit share does separate the two (median 0-0.67,
mostly 0.17-0.67, for absent unsure calls against 0 for present ones; AUC 0.61-0.82,
[`unsure_output.txt`](unsure_output.txt)),
but at F1 ≈ 0.97 removing a false positive gains about as much as losing a true positive costs, so a
demotion pays only where more than half of the demoted calls are absent. With 3-35% absent in the band,
no threshold gets there. The genera reported instead are right in 98-100% of cases: as a policy that prefers
no species name to a wrong one, the rule raises the share of right output by up to 0.3 points (pe 0.9935
to 0.9948, se 0.9900 to 0.9931, ONT unchanged) for 0.3-0.6 points of species sensitivity. On F1 it loses.

## (c) Alignment evidence as features: small gains, two of them clear

Change of F1 against the shipped feature set (cross-validation; test set with interval):

| features added | pe | se | pb | ont |
|---|---|---|---|---|
| `mean_mapq`, low-MAPQ share | +0.0006; +0.0010 (−0.0019, 0.0039) | +0.0019; +0.0057 (−0.0009, 0.0156) | +0.0069; −0.0043 (−0.0092, 0.0006) | +0.0006; +0.0030 (−0.0028, 0.0093) |
| mate concordance | +0.0008; −0.0014 (−0.0047, 0.0016) | | | |
| clipped share | +0.0001; +0.0001 | +0.0018; +0.0004 | +0.0023; +0.0011 | −0.0002; +0.0045 (−0.0007, 0.0124) |
| congener-fit shares | +0.0009; +0.0017 (−0.0005, 0.0041) | +0.0040; +0.0029 (−0.0008, 0.0069) | +0.0015; −0.0008 | −0.0001; +0.0021 |
| all of these and `identity_z` | +0.0008; +0.0011 (−0.0030, 0.0055) | **+0.0054; +0.0081 (0.0006, 0.0192)** | +0.0082; −0.0027 (−0.0082, 0.0028) | −0.0001; **+0.0061 (0.0011, 0.0135)** |

Only single-end reads gain in both cross-validation and the test set (se with all of them: F1 0.9797 to
0.9851 cross-validated, 0.9642 to 0.9721 on the test set, false negatives 89 to 64). The congener-fit
share is the most consistent single addition for short reads (cross-validated false positives: pe 25 to
23, se 50 to 34). For PacBio, cross-validation improves and the test set does not; the test set has 6
samples, so these long-read numbers are the least certain.

## Ideas that this world cannot test

- **Where the mismatches fall in the codon.** A relative's differences are mostly at third codon
  positions, sequencing errors anywhere. The tuning world mutates every position alike
  (`simulate_gtdb_release.py`), so it would show nothing here; it needs real genomes.
- **Which genes the reads hit.** A relative's reads land on the genes that differ least or most between
  species ([the gene-scaled margin report](../2026-10-01-gene-scaled-margin/README.md) found them on the fast
  genes); V2's world has one rate for all genes. The stress worlds with gene rates could test a feature
  comparing a taxon's hit genes with `gene_conservation.tsv`.
- **Reassigning ambiguous reads** (an EM over the reads that fit several taxa) before computing the
  features, so that a present congener with unique support keeps them. Not implemented.

## What follows

- Neither proposed rule raises F1 here: the depth-scaled identity is already in the forest, and the unsure
  calls are mostly real species whose reads also fit a congener. Genus-only reporting is a choice of
  precision over sensitivity, not an F1 gain.
- The congener-fit share and the MAPQ shares are worth adding to protal's features: protal aligns the
  candidates anyway (align top 3), so they cost no extra alignment, only a per-taxon count. The evidence
  here is modest (clear for single-end, positive but uncertain for paired-end); the GTDB r226 build, with
  its larger and real-genome test set, decides. The same scripts run on it (`V2=OUTDIR FPEXP=...`) once its
  points are realigned with `-m 3`.
