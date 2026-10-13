# The r226 v21 build: what the strain alleles, the polymorphic sites and the ancestry sites gain

**Data.**
- **The build:** r226 v21 (SLURM 24112881, q1536n4, 2026-10-09 11:25 to 14:42, 120 threads).
  - Code: protal and the simulator at `ea6de54`, the scripts at `7ace9b8` (`build_metadata.tsv`), seed 1.
  - It is the first r226 build with the strain alleles (`504de80`) and the polymorphic sites (`ea6de54`).
  - Its share archive `local/protal0.7.9_r226_v21_share.tar.gz` (543 MB) is unpacked to `local/v21`.
- **The comparison:** v19 (`a81cee3`, [2026-10-09-r226-v19](../2026-10-09-r226-v19/README.md), `local/v19`).
  - The same seed, design, scenarios and hold-out (10,620 species left out, against 10,621).
  - v21 differs in its code and in `--allele-genome-share 0.5` (below).
- **The retraining:**
  - The models were refitted on v21's own training and test tables (`work/`) with `machine_learning_cmdline.py` of
    `7ace9b8`, using the build's options (`--model gbm --maxnodes 63 --seed 1 --scenario-weight 0.25 --evaluation basic
    --depth-knobs`).
  - Only `--features` changed, leaving groups out ([ablate.sh](ablate.sh), [ablate_all.sh](ablate_all.sh),
    [ablate_more.sh](ablate_more.sh)).
  - They ran in WSL on 4 threads, with scikit-learn 1.9.1, numpy 2.5.3 and pandas 3.0.6, as on the cluster.
  - The refit of the full set reproduces the cluster's model: F1 with species held out 0.9549 against 0.9548, and the
    test set 0.9703 against 0.9716. The paired difference is −0.0001 ± 0.0004 with species held out and +0.0014 ± 0.0024
    on the test set, which is the noise floor of a refit.
- **The comparisons:**
  - [variants.py](variants.py): F1 at knob 0.5, log loss, a paired bootstrap over samples (1000 draws), and errors by
    class.
  - [alleles_classes.py](alleles_classes.py): v19 against v21 by genome class and allele status.
  - [fp_sources.py](fp_sources.py): where the false positives' reads come from (`model_logs/error_reads/*/taxa.tsv.gz`).
  - [cross_score.py](cross_score.py) and [cross.sh](cross.sh): each build's test samples scored with models fitted on
    either build's training samples, both with v19's feature set.

  Their outputs are the `*.txt` files here. The paths in the scripts are this laptop's.

**Question (the user's, 2026-10-09).** The gains of the recent additions seem rather minor. How much did the ancestral
sites gain, and how much the strain alleles?

## Summary

- **Each addition is worth about +0.001-0.002 of F1, and together +0.003-0.004.** With species held out, on v21's own
  samples, leaving out a group of the classifier's features costs:

  | Left out | pe | se | pb | ont |
  |---|---|---|---|---|
  | polymorphic sites (3 features) | −0.0003 (noise) | | | |
  | strain alleles + polymorphic sites (5) | −0.0015 | −0.0010 | −0.0017 | −0.0014 |
  | ancestry sites (5) | −0.0018 | −0.0016 | −0.0010 | −0.0022 |
  | all three groups (10) | −0.0029 | −0.0029 | −0.0027 | −0.0040 |

  - Every difference except the polymorphic sites' is outside its 95% interval (±0.0004 to ±0.0009).
  - On the test sets the differences are within noise, except the strain alleles for pb (−0.0042 ± 0.0026) and all
    three groups for pe (−0.0026 ± 0.0022).
  - The two groups add up: together they gain about the sum of their separate gains (88-112%).
  - **The polymorphic sites add nothing measurable** on top of the strain alleles.
- **The strain alleles do what they were meant to, at a price.**
  - **Real strains of species with alleles** are missed 15-27% less:
    - pe 6.5% → 5.5%;
    - pb 7.7% → 5.6%;
    - ont 8.1% → 6.5%;
    - se 7.1% → 5.7%.
  - **Real strains of species without alleles** are missed 10-15% more (pe 4.7% → 5.2%, pb 5.4% → 6.2%).
  - **FP beside a held-out congener:** fewer where the database species has no alleles, more where it has them.
  - In the shipped database every species with several genomes has alleles. So the costly class is partly an artifact of
    the 50% genome split. The in-silico strains, which stand in for strains of one-genome species, are unaffected.
- **Most of v21's step up from v19 is not the additions but an easier simulation.**
  - With species held out, v21 is +0.009 (pe) and +0.007 (se, pb, ont) above v19.
  - `--allele-genome-share 0.5` keeps half the species' other genomes for the alleles. So real strains fell from 47% to
    28% of the present species, and in-silico strains rose from 24% to 40%. Real strains are missed 2-4 times as often.
    - That mix alone is worth +0.003-0.004 of v21's F1.
    - The same v19 model scores +0.0023 to +0.0027 better on v21's test samples than on v19's.
  - **The model trained on v21's samples calls more conservatively:** up to 17% fewer FP and as many or more FN. The FP
    it saves are mostly the same-genus ones from held-out congeners.
    - That trade pays on v21's mix.
    - On v19's samples it gains nothing (pe +0.0001, se +0.0004) or loses (pb −0.0021, ont −0.0006).
- **The alignment shifts cannot be measured from this build.**
  - In the short-read runs, 9-10% of alignments were unsure (a candidate of another species within 3 mismatches), and
    0.4-0.6% of alignments moved to another species.
  - Separating their effect needs the same samples profiled with `--no_allele_scores`.
  - The long reads, which get no shifts, rose as much from v19 as the short reads did. So the shifts are no large part of
    the step.
- **Recommendations:**
  1. Restore the real-strain share in the next build, so that its numbers compare with v19's and its models are not tuned
     to an easier mix.
  2. Profile the test set once with `--no_allele_scores` to measure the shifts.
  3. Take the polymorphic group out of the default set unless the next build shows a gain.

## 1. v21 against v19

F1 at knob 0.5 (`build_metadata.tsv`; v19's pe calls at its knob 0.74 are in its report):

| | pe | se | pb | ont |
|---|---|---|---|---|
| species held out, v19 → v21 | 0.9462 → 0.9548 | 0.9531 → 0.9602 | 0.9681 → 0.9754 | 0.9625 → 0.9700 |
| test set, v19 → v21 | 0.9641 → 0.9716 | 0.9599 → 0.9661 | 0.9816 → 0.9841 | 0.9794 → 0.9755 |
| soil hold-out, v19 → v21 | 0.9541 → 0.9625 | 0.9517 → 0.9599 | 0.9684 → 0.9733 | 0.9586 → 0.9635 |
| shallow soil hold-out, v19 → v21 | 0.9269 → 0.9313 | | 0.9631 → 0.9684 | 0.9525 → 0.9614 |
| FP / FN with species held out, v19 | 5615 / 2698 | 3079 / 1906 | 1353 / 1859 | 1791 / 2083 |
| FP / FN with species held out, v21 | 4665 / 2426 | 2558 / 1762 | 1108 / 1422 | 1416 / 1731 |

**What changed in the samples.**
- **The genome table:** 37,422 genomes instead of 49,595. 12,173 genomes are kept for the strain alleles and not
  simulated (`--allele-genome-share 0.5`).
- **The in-silico strains:** 14,746 instead of 9,030, now with 44% of their marker substitutions at congener sites
  (`insilico_strains.log`: "congener share 0.437", measured on the real strains). 5,071 of the 14,746 have a congener.
- **The present species with species held out:**

  | | real strains | in-silico strains | representatives |
  |---|---|---|---|
  | v19 | 47% | 24% | 29% |
  | v21 | 28% | 40% | 32% |

  It is the same in every read type.

The FN rate by class with species held out, at knob 0.5:

| | real strain | in-silico strain | representative |
|---|---|---|---|
| pe v19 → v21 | 5.31% → 5.35% | 2.14% → 2.62% | 1.94% → 1.84% |
| se | 5.54% → 5.49% | 2.10% → 2.69% | 1.77% → 1.99% |
| pb | 6.60% → 5.82% | 1.42% → 1.95% | 0.94% → 1.11% |
| ont | 6.84% → 6.33% | 1.66% → 2.49% | 1.58% → 1.66% |

- The in-silico strains became harder, as `ea6de54` intended. They are still missed half as often as real strains.
- With v21's rates and v19's mix, v21's F1 with species held out would be:
  - pe 0.9521 (mix +0.0027);
  - se 0.9574 (+0.0028);
  - pb 0.9716 (+0.0038);
  - ont 0.9662 (+0.0037).

  This counts only the FN; real strains' reads also spill onto congeners (below).

**Where the false positives come from** ([fp_sources.txt](fp_sources.txt)). Each FP is classed by the source of most of
the reads on it, with species held out:
- **Reads of a held-out species:** 82-92% of the FP in both builds.
- **Reads of a present species:**
  - per 1000 present real strains: 9-14 (pe), 12-14 (se), 2.5 (pb), 3.3-3.7 (ont);
  - per 1000 in-silico strains or representatives: 5-9 (pe and se), 1.4-2.7 (long reads).

From v19 to v21 the FP of held-out species' reads in the same genus fell by 20-24%, in every read type. Those in another
genus did not change:

| | same genus, v19 → v21 | other genus, v19 → v21 |
|---|---|---|
| se | 2075 → 1655 | 392 → 385 |
| pb | 1119 → 884 | 127 → 125 |
| ont | 1464 → 1107 | 167 → 162 |

(pe: 2204 → 3224, at v19's knob 0.74 against 0.5.) The held-out species per sample are the same (pe 255.5 against
257.5).

## 2. The feature groups, on v21's samples

Outputs: [variants_pe.txt](variants_pe.txt), [variants_se.txt](variants_se.txt), [variants_pb.txt](variants_pb.txt),
[variants_ont.txt](variants_ont.txt). Log loss with species held out:

| | pe | se | pb | ont |
|---|---|---|---|---|
| all features (85) | 0.0531 | 0.0571 | 0.0816 | 0.0844 |
| without the polymorphic sites | 0.0534 | | | |
| without strain alleles + polymorphic sites (v19's set) | 0.0550 | 0.0591 | 0.0859 | 0.0880 |
| without the ancestry sites | 0.0546 | 0.0590 | 0.0833 | 0.0887 |
| without all three | 0.0562 | 0.0612 | 0.0873 | 0.0926 |

The F1 differences against the full set (paired bootstrap over samples, mean [2.5%, 97.5%]):

| Left out | pe held out | pe test | se held out | pb held out | pb test | ont held out |
|---|---|---|---|---|---|---|
| polymorphic | −0.0003 [−0.0008, +0.0002] | +0.0017 [−0.0003, +0.0039] | | | | |
| alleles + polymorphic | −0.0015 [−0.0019, −0.0010] | −0.0001 [−0.0018, +0.0016] | −0.0010 [−0.0017, −0.0004] | −0.0017 [−0.0023, −0.0011] | −0.0042 [−0.0069, −0.0018] | −0.0014 [−0.0021, −0.0009] |
| ancestry | −0.0018 [−0.0024, −0.0012] | −0.0003 [−0.0020, +0.0013] | −0.0016 [−0.0021, −0.0009] | −0.0010 [−0.0015, −0.0006] | −0.0009 [−0.0027, +0.0009] | −0.0022 [−0.0028, −0.0015] |
| all three | −0.0029 [−0.0036, −0.0022] | −0.0026 [−0.0047, −0.0003] | −0.0029 [−0.0035, −0.0022] | −0.0027 [−0.0033, −0.0020] | −0.0022 [−0.0045, +0.0000] | −0.0040 [−0.0048, −0.0030] |

**The scenarios' hold-out samples** gain the most from all three groups:
- pe: soil +0.0041, shallow soil +0.0051, moderate +0.0024;
- ont: soil +0.0041.

**Where the strain alleles act** (species held out, FN rate of present species, FP rate of absent ones). A species has
alleles when its records lie on gene copies with alleles: `allele_copy_share` above 0.2 in any row of 3 or more fragments.

| Class | pe: v19's set → all | se | pb | ont |
|---|---|---|---|---|
| real strain, alleles | 6.49% → 5.50% | 7.09% → 5.69% | 7.71% → 5.61% | 8.13% → 6.52% |
| real strain, no alleles | 4.67% → 5.17% | 4.59% → 5.20% | 5.41% → 6.22% | 5.52% → 6.24% |
| in-silico strain | 2.55% → 2.60% | 2.60% → 2.72% | 1.93% → 1.91% | 2.49% → 2.52% |
| representative | 1.80% → 1.82% | 1.92% → 1.94% | 1.09% → 1.10% | 1.70% → 1.65% |
| absent beside a held-out congener, alleles | 2.58% → 2.79% | 2.30% → 2.63% | 3.16% → 3.38% | 3.56% → 3.96% |
| absent beside a held-out congener, no alleles | 1.41% → 1.24% | 1.41% → 1.22% | 2.03% → 1.80% | 2.02% → 1.76% |

The features help where the reads have alleles to compare with. Where they do not, an explained share of 0 is read as
evidence against presence. That costs the real strains of species without alleles and spares the absent species
without alleles.
- **In the training database:** a species with several genomes has no alleles when the hash split gave all of its other
  genomes to the simulation.
- **In the shipped database:** every species with several genomes has alleles (4.2 million alleles of 33,279 species,
  against 29,019 in the training database). The strains without alleles are then those of one-genome species, which
  the in-silico strains stand for, and those are unaffected.

**How well the features separate** ([alleles_classes.txt](alleles_classes.txt), v21's training rows on copies with
alleles). The AUC is of present real strains over absent species beside a held-out congener:

| | all rows (pe / pb) | the hard rows (pe / pb) |
|---|---|---|
| allele_explained_share | 0.72 / 0.70 | 0.59 / 0.57 |
| allele_identity_gain | 0.60 / 0.55 | 0.56 / 0.56 |
| polymorphic_known_share | 0.50 / 0.37 | 0.55 / 0.54 |
| polymorphic_novel_share (lower for strains) | 0.42 / 0.24 | 0.45 / 0.45 |
| ancestry_fixed_gain | 0.59 / 0.47 | 0.52 / 0.54 |
| ancestry_agreement | 0.86 / 0.91 | 0.52 / 0.50 |
| identity | 0.96 / 0.96 | 0.57 / 0.53 |

"The hard rows" are the strains with p < 0.9 and the absent species with p > 0.1, with species held out. On those, every
feature is weak (0.52-0.59). The strain alleles are about as good as identity, which is why they add a little.

**Their share of the split gain** (`model_logs/trained_model*.varimp.tsv`):

| | pe | se | pb | ont |
|---|---|---|---|---|
| strain alleles + polymorphic sites | 0.7% | 0.8% | 1.2% | 1.0% |
| ancestry sites | 1.0% | 1.2% | 2.5% | 3.8% |
| `identity` alone | 72% | | | |
| `excess_scaled_median` alone | | | 68% | |

The polymorphic novel share and the ancestry indel features rank 58th to 85th.

## 3. Model or samples: each build's model on the other build's samples

[cross_score.txt](cross_score.txt). The models were fitted with v19's feature set on each build's full training table. The
test tables of both builds were scored at knob 0.5. All rows (design and scenarios):

| | pe | se | pb | ont |
|---|---|---|---|---|
| v19 model, v19 samples | 0.9511 (1657/1212) | 0.9584 (821/816) | 0.9725 (423/689) | 0.9649 (604/846) |
| v21 model, v19 samples | 0.9512 (1652/1214) | 0.9588 (778/839) | 0.9704 (392/799) | 0.9643 (527/943) |
| v19 model, v21 samples | 0.9538 (1641/1075) | 0.9607 (823/747) | 0.9750 (450/567) | 0.9674 (598/748) |
| v21 model, v21 samples | 0.9571 (1438/1075) | 0.9628 (681/799) | 0.9757 (373/615) | 0.9689 (500/781) |

(F1, FP/FN.)
- **The same model scores better on v21's samples:**
  - F1 +0.0023 to +0.0027;
  - 8-18% fewer FN;
  - about the same FP (pb +6%).

  These are the fewer real strains.
- **The model trained on v21's samples:**
  - makes 0-17% fewer FP and as many or more FN, on either build's samples;
  - gains on v21's samples (pe +0.0033, se +0.0021, ont +0.0015, pb +0.0007);
  - does not gain on v19's (pe +0.0001, se +0.0004, ont −0.0006, pb −0.0021).

  Trained on fewer real strains, it calls less readily. On a mix with more real strains (v19's, and probably real
  samples') that costs as much as it gains.
- **For pe and se**, v21's samples also carry the alignment shifts, which v19's model never saw. Their +0.0033 and
  +0.0021 on v21's samples, against +0.0007 to +0.0015 for the long reads, may hold a small effect of the shifts. It is
  within the test sets' noise.

## 4. The alignment shifts

The protal runs' per-sample lines (`logs/protal_runs_*.log`, "strain alleles:"), summed:

| | samples | alignments made | unsure reads that took their shifts | moved to another species |
|---|---|---|---|---|
| training, pe | 389 | 116.8 M | 11.5 M (9.8%) | 693 k (0.59%) |
| training, se | 309 | 16.1 M | 1.5 M (9.4%) | 64 k (0.40%) |
| test, pe | 110 | 42.6 M | 3.9 M (9.1%) | 207 k (0.49%) |
| test, se | 78 | 3.3 M | 0.3 M (9.6%) | 13 k (0.39%) |

The training database has 3,666,134 alleles of 1,969,331 gene copies of 29,019 species (`training_db_index.log`). Long
reads get the features but no shifts.

Whether the moved reads went the right way cannot be read from the build:
- the per-read truth is in the SAMs on the node's disk;
- the error-read SAMs hold only error taxa.

The direct test is to profile the test set's samples again with `--no_allele_scores` and score both with the same
model.

## 5. The ancestry report

The fixed-site AUC fell from v19 to v21 (`console.log`). The FN taxa differ between the builds (fewer real strains, the
shifts), so the two are not a like-for-like comparison of the sites:

| | v19 | v21 |
|---|---|---|
| pe | 0.664 | 0.621 |
| se | 0.679 | 0.625 |
| pb | 0.583 | 0.453 |
| ont | 0.599 | 0.510 |

## 6. What to do next

1. **Restore the real-strain share.** The split is what keeps the alleles honest, so do not undo it. Instead:
   - simulate a real strain more often where a species has a simulation genome left (`--strains`);
   - or lower `--allele-genome-share`, so that fewer species lose their only other genome;
   - or, with a new rule, keep alleles only for species with three genomes or more.

   Then the next build's F1 compares with v19's, and its model is not tuned to a mix easier than real samples.
2. **Measure the shifts:** profile the test set with `--no_allele_scores` (protal) and score it with the same model.
3. **The polymorphic group:** take it out of the default set if the next build confirms no gain. Its three features
   gained −0.0003 to +0.0017, all within noise.
4. **The strain alleles and the ancestry sites stay.** Each gains about +0.0015 with species held out, together
   +0.003-0.004 and 6-9% of log loss. For pe that is about the depth feature's gain at r226 v7/v8 (+0.002-0.004).
