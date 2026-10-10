# The r226 v22 build: column weights, ancestry sites and deep alleles, and what moved the false positives and negatives

**Data.**
- **The build:** r226 v22 (SLURM 24127982, q512n11, 2026-10-10 09:52 to 12:05, 84 threads). protal, the simulator
  and the scripts at `82bcc15` (`build_metadata.tsv`), seed 1. New against v21 (`ea6de54`): the column weights and
  their five `weights` features (91 default features), the weighted site shift (on in `82bcc15`), the four changes
  of the ancestry true-positive test (`ancestry_fixed_agreement`, alleles chosen for coverage, the tolerant
  consensus, `--insilico-multi 0.5`) and the deep alleles' k-mers in the index (`--index_alleles 0.01`).
  The share archive `local/protal0.7.9_r226_v22_share.tar.gz` (534 MB) is unpacked to `local/v22`.
- **The comparison:** v21 ([2026-10-09-r226-v21](../2026-10-09-r226-v21/README.md), `local/v21`): the same seed,
  design, scenarios and hold-out (10,620 species left out).
- **The refits:** the models refitted on v22's own tables (`work/`) with the build's trainer (scripts of
  `82bcc15`) and options (`--model gbm --ntree 64 --maxnodes 63 --seed 1 --scenario-weight 0.25 --evaluation basic
  --depth-knobs`), only `--features` changed ([ablate.sh](ablate.sh), [ablate_all.sh](ablate_all.sh); WSL, 4
  threads). The full set reproduces the cluster's model (pe with species held out 0.9569 against 0.9573).
- **The comparisons:** [variants.py](variants.py) (F1 at knob 0.5, log loss, a paired bootstrap over samples, FN and
  FP by class, the calls that flip between variants; [variants_pe.txt](variants_pe.txt),
  [variants_se.txt](variants_se.txt), [variants_pb.txt](variants_pb.txt), [variants_ont.txt](variants_ont.txt));
  [builds_by_class.py](builds_by_class.py) (v21 against v22 by class, and the F1 at the other build's class mix;
  [builds_by_class.txt](builds_by_class.txt)); [feature_auc.py](feature_auc.py) (each feature alone on real strains
  against absent species beside a held-out congener; [feature_auc.txt](feature_auc.txt)).

**Questions (the user's, 2026-10-10).** Is there evidence that the ancestral states or the column weights are useful?
How were the FP and FN affected by the new features? (The deep alleles had already been judged not worth their cost
from v22's headline numbers and turned off, `84af9a4`.)

## Summary

- **v21 → v22 is mostly an easier mix again.** `--insilico-multi` raised the in-silico strains from 40% to 49% of
  the present taxa (real strains 28% → 21%), and in-silico strains are missed half as often. At the same class mix
  the code changes together are worth +0.0004 (pe), +0.0006 (se), −0.0011 (pb), −0.0013 (ont); the mix
  +0.0014 to +0.0021. The error rates per class barely moved (real strains with alleles: pe 5.51% → 5.42%, pb 5.52%
  → 6.19%, ont 6.50% → 7.14%; absent beside a held-out congener pe 1.63% → 1.59%).
- **The column weights add nothing measurable**, in any read type (−0.0001 to −0.0004 of F1 left out, all within
  noise), and on the hard rows none of them separates (AUC 0.47-0.50). They could not: two of the five were 0 in
  every row (a cap in the estimate), the weighted shift discounted every column (another cap), 42% of their
  alignments failed, and their options were never read from the command line (section 4). All four are fixed since
  (`0d961e1`, `c60c272`); the weights need a new build to be judged.
- **The ancestry sites pay, a little, against false positives.** Left out: pe −0.0013 [−0.0019, −0.0007], se
  −0.0007, ont −0.0008, pb −0.0003 (noise); with the weights left out too −0.0009 to −0.0018, all significant (the
  weighted agreement stands in for the group). Most of it is FP beside held-out congeners (real species): +41 to +150
  of them without the sites. The FN they prevent are in-silico strains (+16 to +57 without them); for real strains
  with alleles the sites save nothing (−13 to +14), which says the in-silico strains, mutated copies of the
  representative that keep all its derived states, flatter them.
- **On the hard rows** (strains with p < 0.9, absent with p > 0.1) the ancestry agreement points the wrong way (AUC
  0.43-0.47), the fixed-site features and the alleles are the only ones above 0.55 (0.56-0.60).
- **The deep alleles** (`--index_alleles 0.01`): 115 M seeds of 1.5 M alleles at r226, `database.protal` 24.1 →
  25.2 GB, the index 36.7 → 39.5 GB resident, for at most the code's +0.0005; off by default since `84af9a4`.

## 1. v21 against v22

F1 at knob 0.5 with species held out (`build_metadata.tsv`): pe 0.9548 → 0.9573, se 0.9602 → 0.9627, pb 0.9754 →
0.9757, ont 0.9700 → 0.9704. By class ([builds_by_class.txt](builds_by_class.txt)), each build's own model:

| pe, species held out | v21 share | v22 share | v21 rate | v22 rate |
|---|---|---|---|---|
| FN real strain, alleles | 16.6% | 12.4% | 5.51% | 5.42% |
| FN real strain, none | 11.6% | 8.9% | 5.11% | 5.35% |
| FN in-silico strain | 39.8% | 48.6% | 2.62% | 2.78% |
| FN representative | 32.1% | 30.1% | 1.84% | 1.74% |
| FP beside a held-out congener | 69.8% | 69.2% | 1.63% | 1.59% |
| FP beside present congeners only | 12.1% | 12.4% | 0.67% | 0.62% |
| FP no congener in the sample | 18.1% | 18.4% | 1.30% | 1.28% |

F1 from the class rates and shares (rates of the row's build at the column's mix):

| | pe v21 mix | pe v22 mix | se v21 | se v22 | pb v21 | pb v22 | ont v21 | ont v22 |
|---|---|---|---|---|---|---|---|---|
| v21 rates | 0.9548 | 0.9570 | 0.9602 | 0.9621 | 0.9754 | 0.9767 | 0.9700 | 0.9714 |
| v22 rates | 0.9552 | 0.9573 | 0.9608 | 0.9627 | 0.9743 | 0.9757 | 0.9687 | 0.9704 |

## 2. The feature groups, on v22's samples

Species held out, the paired bootstrap's mean [2.5%, 97.5%] of the variant's F1 minus the full set's:

| left out | pe | se | pb | ont |
|---|---|---|---|---|
| `weights` (5) | −0.0002 [−0.0007, +0.0004] | −0.0004 [−0.0009, +0.0002] | −0.0002 [−0.0007, +0.0003] | −0.0001 [−0.0007, +0.0005] |
| `ancestry` (5) | −0.0013 [−0.0019, −0.0007] | −0.0007 [−0.0012, −0.0001] | −0.0003 [−0.0009, +0.0003] | −0.0008 [−0.0014, −0.0002] |
| both (10) | −0.0018 [−0.0024, −0.0012] | −0.0009 [−0.0016, −0.0002] | −0.0012 [−0.0018, −0.0005] | −0.0018 [−0.0024, −0.0011] |

On the test sets everything is within noise except both groups for pe (−0.0031 [−0.0054, −0.0009]) and se (−0.0021
[−0.0041, −0.0001]). Log loss with species held out, full → without both: pe 0.0518 → 0.0535, se 0.0562 → 0.0575,
pb 0.0811 → 0.0823, ont 0.0860 → 0.0894.

**The calls that flip** (species held out, at knob 0.5; net = errors the variant adds minus those it removes):

| net change of errors | pe −anc | pe −both | se −anc | se −both | pb −anc | pb −both | ont −anc | ont −both |
|---|---|---|---|---|---|---|---|---|
| FP beside held-out congener | +150 | +258 | +70 | +84 | +41 | +102 | +78 | +137 |
| FN in-silico strain | +18 | +28 | +16 | +22 | +23 | +48 | +57 | +102 |
| FN real strain, alleles | +14 | +11 | −7 | −11 | −3 | −1 | −13 | −15 |
| FP no congener in the sample | +23 | +16 | +7 | +21 | +2 | +4 | +3 | +10 |

Without the weights alone the nets are small in every class (pe: FP beside held-out +21, FN real strains without
alleles −17).

**How the models use them** (`model_logs/trained_model*.varimp.tsv`, share of the split gain): `ancestry` 1.1% (pe),
1.3% (se), 2.7% (pb), 4.8% (ont); `weights` 0.4-0.5%, two of its features ranked 90th and 91st of 91;
`ancestry_fixed_*` 0.4-0.9%; `alleles` 0.3-0.4%.

**Each feature alone** ([feature_auc.txt](feature_auc.txt)), present real strains against absent species beside a
held-out congener, all rows / the hard rows:

| | pe | pb |
|---|---|---|
| identity | 0.965 / 0.537 | 0.965 / 0.485 |
| ancestry_agreement | 0.884 / 0.471 | 0.938 / 0.435 |
| ancestry_agreement_weighted | 0.866 / 0.466 | 0.927 / 0.421 |
| ancestry_fixed_agreement | 0.719 / 0.572 | 0.709 / 0.562 |
| allele_explained_share | 0.656 / 0.593 | 0.638 / 0.581 |
| conserved_mismatch_ratio | 0.679 / 0.483 | 0.425 / 0.494 |
| nonsynonymous_share | 0.711 / 0.497 | 0.390 / 0.492 |
| conserved_mismatch_rate, nonsynonymous_conserved_rate | 0.500 / 0.500 (constant) | 0.500 / 0.500 |

## 3. The alignment shifts

`logs/protal_runs_*.log`, "strain alleles:" summed: training pe 10.1 M unsure reads took their shifts (v21 11.5 M),
582 k moved to another species (693 k); test pe 3.9 M, 238 k (207 k); se 1.5 M, 62 k (64 k). The weighted shift
moved about as many reads as v21's shift did; with every column discounted (section 4) it carried no column
information.

## 4. What was wrong with the column weights

1. **Two features constant 0.** `conserved_mismatch_rate` and `nonsynonymous_conserved_rate` were 0 in all 379,029
   training rows (and in every row of the truth-test worlds). Each genus's share took its own pseudocounts
   (k + 1)/(n + 2) before the family's mean, with at most 24 species per genus: no column could exceed 25/26 (code 7),
   below the "conserved" code 9 (98.9%); the amino acids alike. Fixed in `0d961e1` (the genera's raw shares averaged,
   the pseudocounts over all the species compared).
2. **The weighted shift discounted every column.** The among estimate of at most 10 genus references peaks at 11/12
   (code 5), below the shift's full-difference threshold (code 6). Not fixed (the shift is off by default since
   `0b34254`).
3. **42% of the alignments failed** (4.8 of 11.4 M; the cause was not logged). On real GTDB families the bases'
   limits (0.2 within a genus, 0.4 among genera) are below the real divergence
   ([2026-10-10-real-ancestry](../2026-10-10-real-ancestry/README.md)). Since `c60c272` the copies are aligned as
   proteins (98.8% agreement with GTDB-Tk's alignments, seven times the coverage), and the build logs the failures by
   cause and the codes' distribution.
4. **The options were never read.** `--column_weights`, `--no_site_weights` and `--weighted_site_shift` were declared
   but never copied from the command line: the defaults held whatever it said. Fixed in `c60c272`.

## 5. Other lines of the build

- The ancestry report's AUCs (errors' reads, FN against FP) fell from v21: pe fixed sites 0.621 → 0.573, se 0.625 →
  0.595, pb 0.453 → 0.436, ont 0.510 → 0.485 (different taxa; not a like-for-like comparison).
- Composition, pe test median explained-share error +0.019 → +0.003, unknown share −0.025 → −0.005; pb's
  explained share still +0.199 (v21 +0.191).

## 6. What to do next

1. **Judge the column weights on a build with the fixes** (`0d961e1`, `c60c272`): the build log's "Column weights
   codes" line says whether conserved columns now exist; the training tables whether reads land on them.
2. **Measure the ancestry sites on real strains**, not in-silico ones: the FN they prevent here are in-silico
   strains', which keep the representative's derived states by construction. The real-data report measures the
   sites on real strains and congeners.
3. **Keep the deep alleles off** (`84af9a4`).
