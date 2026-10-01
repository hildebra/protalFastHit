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

## In protal: the features, a read consensus, and mates that guide each other

Implemented afterwards (uncommitted at the time of writing, on top of `a1bf586`, in a tree that also held
another session's performance changes; built in WSL as `~/protal-cons/build/protal`, kept as
`~/protal-cons/protal-before-fix`, which the results below use before another session extended `MateGuidance.h`
with the database's gene neighbours; the corrected PacBio and ONT figures come from the shared tree with the
correction, rebuilt in the same place, on a database without gene neighbours):

- **The features.** Each read's best record carries `ZA:Z:<taxid>:<edits more>,...`: its other candidates in
  other taxa among those protal aligns (align top 3), with at most 5 edits more than the best, or `*`
  (`AlternativesTag`, AlignmentOutputHandler.h; long reads per gene). The profiler counts every read's best
  record for its taxon before its MAPQ and length filters drop any (a read that fits a congener as well has
  MAPQ near 0, and the first version, which counted after the filters, missed exactly those: its se gain was
  gone; [`new_build_output.txt`](new_build_output.txt)), and the dump gets `low_mapq_share`,
  `congener_fit_share`, `other_genus_fit_share` (within one edit) and `linked_share` (reads with two records on
  the taxon: both mates, or two genes of a long read). `scripts/model_features.py` adds `mean_mapq`,
  `low_mapq_share`, `congener_fit_share` and `other_genus_fit_share` to the normalized features (28), so that
  `build_gtdb_database.py` trains with them; `linked_share` stays out (below).
- **A read consensus** (ReadConsensus.h). A long read's genes, and the two mates of a pair that aligned to two
  genes, take one taxon: taxa are compared on the parts where both have a candidate (so that a gene the
  database lacks for the read's species, which aligns to a relative alone, outvotes nothing), and every part
  takes its alignment of the winner; a part without one, and a long read's gene whose best hit is clearly
  another taxon's (MAPQ 4 or more against the winner's hit: a chimeric read, or a homolog of a gene elsewhere on
  the read), keeps its best and is written with MAPQ 0 (long reads: `ZR:i:2`). This replaces the long reads'
  vote (two thirds of the confident genes, which settled only the ambiguous ones). The clearly-other rule came
  after the results were first taken (correction below).
- **Mate guidance** (MateGuidance.h, `--no_mate_guidance` turns it off). When a pair's best candidate has one
  mate only, that mate has MAPQ 20 or more and the other mate has no candidate of its taxon: the other mate's
  own anchor of the taxon is aligned (only a read's best anchors are), else the other mate is looked for on the
  first one's gene where the fragment can reach (1,000 bases on its strand's side), placed by the diagonal of
  its 12-mers and aligned as from any anchor, in part if it runs past the gene's end.

Tested on the tuning world again: `~/tune` was deleted to free space during this work, and
[`regen_v3.sh`](regen_v3.sh) rebuilt it from its seeds and ran the V2 build command with the feature build
(`V3`); its test sets have the same taxa as V2's (pe 3,672, se 3,394, pb 758, ont 958), so the samples are
V2's. [`profile_new_build.sh`](profile_new_build.sh) profiled V3's samples again with the build that adds the
consensus and mate guidance, [`tables_new_build.py`](tables_new_build.py) made the tables, and
[`run_new_build.py`](run_new_build.py) compares feature sets (5 seeds, as above; the shipped 24 features are
frozen in `exp_lib.py`). Outputs: [`v3_features_output.txt`](v3_features_output.txt),
[`cons_features_output.txt`](cons_features_output.txt), and for PacBio and ONT after the correction
[`fix_long_reads_output.txt`](fix_long_reads_output.txt) ([`fix_long_reads.sh`](fix_long_reads.sh)).

Test F1 (mean of 5 seeds; in brackets the change against the shipped features of the same build, with the
paired bootstrap interval):

| read type | shipped features, feature build | + MAPQ, congener fit | + all five | shipped features, + consensus and guidance | + MAPQ, congener fit | + all five |
|---|---|---|---|---|---|---|
| pe | 0.9765 | 0.9774 (+0.0007) | 0.9782 (+0.0016) | 0.9769 | **0.9795** (+0.0026; −0.0007, 0.0072) | 0.9785 (+0.0016) |
| se | 0.9638 | 0.9703 (+0.0067; 0.0000, 0.0165) | **0.9707** (+0.0072; 0.0008, 0.0163) | (single-end reads have neither) | | |
| pb | 0.9731 | 0.9712 (−0.0019) | 0.9731 (0) | **0.9757** | 0.9741 (−0.0018) | 0.9747 (−0.0011) |
| ont | 0.9542 | 0.9552 (+0.0010) | 0.9536 (−0.0005) | 0.9541 | **0.9596** (+0.0056; −0.0036, 0.0179) | 0.9561 (+0.0020) |

Test precision and sensitivity (mean of 5 seeds; FP and FN per seed over the read type's test samples):

| read type | shipped features, feature build | + MAPQ, congener fit, feature build | + MAPQ, congener fit, + consensus and guidance |
|---|---|---|---|
| pe | 0.9926 / 0.9610 (FP 10.8, FN 58.8) | 0.9923 / 0.9629 (11.2, 56.0) | 0.9940 / 0.9654 (8.8, 52.2) |
| se | 0.9881 / 0.9406 (16.6, 87.4) | 0.9875 / 0.9538 (17.8, 68.0) | (as the feature build) |
| pb | 0.9931 / 0.9539 (4.0, 27.8) | 0.9917 / 0.9516 (4.8, 29.2) | 0.9828 / 0.9655 (10.2, 20.8) |
| ont | 0.9892 / 0.9216 (6.2, 48.4) | 0.9906 / 0.9222 (5.4, 48.0) | 0.9931 / 0.9284 (4.0, 44.2) |

- Single-end gains clearly from the congener fit (alone +0.0040, interval 0.0003-0.0083): cross-validated false
  positives 51 to 36.
- `linked_share` alone lowers the paired-end F1 on both builds (−0.0020; on the consensus build the interval is
  −0.0040 to −0.0003) and helps nowhere clearly: it is left out of the model's features.
- The consensus changes the calls of long reads most: with the shipped features, PacBio sensitivity 0.954 to
  0.967 (FN 27.8 to 20.0; FP 4.0 to 9.0), ONT FP 6.2 to 3.4; the tables lose taxa that only a read's
  inconsistent part supported (training pe 7,366 to 7,151 rows, pb 1,125 to 1,049, ont 1,458 to 1,334). On the
  test samples 46,892 long-read gene hits took their read's taxon (another best hit or a higher MAPQ) and 1,453
  had no hit of it or a clearly better one of another taxon.
- With the consensus build and the four features: pe 0.9795, se 0.9703, pb 0.9741, ont 0.9596, against
  0.9765, 0.9638, 0.9731 and 0.9542 with the shipped features of the feature build (pb does best with the
  shipped features on the consensus build, 0.9757; its 6 test samples cannot tell these apart).

Correction. The first consensus gave every long-read gene its hit of the read's taxon, also when another taxon's
gene fit clearly better. protal's long-read e2e tests (reads of genes drawn from any taxon) then failed: a gene
became the read taxon's homolog, which could stand twice on one read while the gene itself was lost. Such a gene
now keeps its best hit, with MAPQ 0, as the old vote did. On the tuning world's long reads 310 test gene hits
changed so (47,202 settled and 1,143 inconsistent before), and the test F1 moved by at most 0.0012 (before: pb
0.9762 shipped, 0.9746 with the four features, 0.9746 with all five; ont 0.9534, 0.9598 and 0.9573); the PacBio
and ONT figures above are
those of the corrected build (`fix_long_reads.sh`; paired-end and single-end are unchanged by it).

Mate guidance on the 18 paired-end test samples: 267,075 fragments had a sure mate and another without a
candidate of its taxon; 9 were aligned from their own anchor (the existing anchor recovery already aligns the
other mate's anchors of the sure mate's gene), and 33,358 of 267,066 looked for were found on the gene
([`rescued_mates.py`](rescued_mates.py), on the SAMs with and without guidance): 29,913 mate records that had
none before (1.5% of 1.93 million; the other mates found replace a record on another taxon; none is lost), 90%
of them partial (clipped at the gene's end, median 116 bases clipped), at a median identity of 0.979 (10th
percentile 0.920), and 90% simulated from the genome of the taxon they are written on. These are bases near
gene ends that the strain MSAs had none of.

Cost. Wall-clock times were no use on this machine while other sessions ran their benchmarks (two rounds of
[`time_guidance.sh`](time_guidance.sh), 2 threads pinned at nice 19: user CPU 115.0 s with guidance and 112.6 s
without in the first, about +2%; the second round disagreed with itself by more than that). Instructions instead
([`count_guidance.sh`](count_guidance.sh): callgrind on the OpenMP body of `RunPairedEnd`, the first 40,000
pairs of `rl150_p500000_s_1`, one thread, the same binary with and without `--no_mate_guidance`):

| | instructions | per pair |
|---|---|---|
| without guidance | 8,722,949,244 | 218,000 |
| with guidance | 9,014,050,981 (+3.3%) | 225,000 |
| of which `GuideMate` | 280,962,955 | 82,000 per fragment looked for |
| of which the 12-mer diagonal (`BestDiagonal`) | 178,473,394 | 52,000 per fragment looked for |

3,423 of the 40,000 fragments were looked for and 423 found (none from an anchor). Most of the cost is the
diagonal search of mates that are not there (88% of those looked for): a prefilter that gives up after a few of
the mate's k-mers miss the gene would cut most of it. Index loading and the SAM output are not in these counts;
over a whole run guidance costs less than 3.3%.

The website's documentation of protal's SAM output and options would need the `ZA` tag, `ZR:i:2` and
`--no_mate_guidance` (docs/running.md has them).
