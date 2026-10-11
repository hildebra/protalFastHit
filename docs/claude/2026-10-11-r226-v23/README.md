# The r226 v23 build: what changed, what drives its errors, what each feature group is worth and what it costs

**Data.**
- **The build:** r226 v23 (SLURM 24140515, q512n11, 2026-10-10 22:23 to 2026-10-11 00:29, 84 threads). protal, the
  simulator and the scripts at `8d8c89b` (`build_metadata.tsv`), seed 1, `--index-alleles 0.01` given explicitly
  (the default is off since `84af9a4`). The share archive `local/protal0.7.9_r226_v23_share.tar.gz` (546 MB) is
  unpacked to `local/v23`.
- **The comparison:** v22 ([2026-10-10-r226-v22](../2026-10-10-r226-v22/README.md), `local/v22`, `82bcc15`): the
  same seed, design and scenarios, but other inputs and hold-out (section 1).
- **The refits** (section 4): v23's models refitted on v23's own tables (`work/`) with the build's trainer (the
  trainer and `model_features.py` are the same at `82bcc15` and `8d8c89b`) and options (`--model gbm --ntree 64
  --maxnodes 63 --seed 1 --scenario-weight 0.25 --evaluation basic --depth-knobs`), only `--features` changed:
  every group of the default set left out in turn and three combinations ([ablate.sh](ablate.sh),
  [ablate_all.sh](ablate_all.sh); WSL, 3 threads on cores 0-2), two more combinations ([ablate_extra.sh](ablate_extra.sh))
  and the full set on tables whose eight ancestry-site columns are blanked ([ablate_columns.sh](ablate_columns.sh)),
  both on core 3 with one thread (`CORES=3 THREADS=1`; the scores do not depend on the threads). The local copy of
  the trainer (`~/v23abl/scripts`) skips the rows and samples cross-validations (the species folds are `GroupKFold`,
  deterministic, so the species-held-out probabilities are the same); the full refits reproduce the cluster's models
  with species held out (pe 0.9596 against 0.9595, pb 0.9792 against 0.9791). Compared by
  [ablation_table.py](ablation_table.py) (F1 at knob 0.5, a paired bootstrap over samples, the net change of errors by
  class).
- **The costs** (section 5): the r226 runs' timers and CPU logs (`logs/protal_cpu_*.tsv`, `logs/protal_runs_*.log`;
  [profile_stages.py](profile_stages.py)), and local runs under callgrind (WSL, one thread on core 3) of protal
  `8d8c89b` built with `-g` (`~/b23`) and of a cost-probe build (`~/b23p`: `8d8c89b` + [cost_probe.pl](cost_probe.pl),
  a run-time bitmask `PROTAL_SKIP` that skips one feature group's evidence), profiling (`--profile_only`) the SAMs of
  the small real world's training samples (`~/realworld_sp_s1`, [2026-10-10-real-small-world](../2026-10-10-real-small-world/README.md):
  170 species, its training database) ([cost_runs.sh](cost_runs.sh)), and aligning one new sample of that world
  with and without `--no_allele_scores` ([align_cost_runs.sh](align_cost_runs.sh)).

**Questions (the user's, 2026-10-11).** A summary of the changes, the good and bad parts, what drives the false
positives and negatives; (a) how much the feature groups add to F1, derived quickly; (b) how much they cost the
profiler, with local test runs allowed.

## Summary

- **v23 is better than v22, and this time not by an easier mix.** Species held out: pe 0.9573 → 0.9595, se 0.9627
  → 0.9650, pb 0.9757 → 0.9791, ont 0.9704 → 0.9741. At v22's class mix v23's error rates are worth +0.0016 to
  +0.0024 in every read type. Nearly all of it is fewer false positives beside held-out congeners (−11% pe, −15% se,
  −17% pb and ont). The inputs changed too (9,981 more one-genome species, 12,142 species held out against 10,620).
- **The bad parts:** real strains with alleles are missed more often in every read type (pe 5.42% → 6.51%); the
  independent test sets move within noise (se −0.0010, pb −0.0039, ont −0.0013 on small sets); the composition
  estimates got worse (pe unknown share −0.022, average genome size +7.3%; PacBio's explained share still +0.20);
  only 13.3% of simulated species are real strains against 40% in-silico.
- **What drives the errors:** the 30 complex scenario samples (moderate, soil, soil_shallow) carry 86% of the pe FP
  and FN. FP are absent species near a held-out congener (27% called within 0.02, 5.7% at 0.02-0.05), whose reads
  they take, with plenty of fragments. FN are thin taxa (40% have one fragment) and strains beside congeners.
- **(a) What the features add:** alone, almost every group is replaceable (only divergence −0.0009 and distance
  −0.0005 for pe, none for PacBio). Everything since the v9 core set is worth +0.0055 (pe), +0.0051 (se), +0.0042 (pb) and
  +0.0058 (ONT). The
  three groups built on the ancestry sites (ancestry, polymorphic, weights) are worth pe −0.0022, se −0.0016, pb
  −0.0013, ont −0.0031 when all three go, about v23's whole gain; the eight columns that need the sites carry pe
  −0.0015, se −0.0009, ont −0.0019 of it and PacBio nothing measurable (−0.0001), mostly as false positives beside
  held-out congeners (+173 pe, +130 ont without them).
- **(b) What they cost:** profiling is ~26% of short-read CPU at r226 and 1-4% of long-read CPU. Inside it, building
  the ancestry sites is the largest single cost: 49% of the profiling instructions of a 1 M-pair pe sample (60% for
  PacBio, 38% over five samples in one run), ~2.5 M instructions per (species, gene) copy at r226, and complex r226
  samples touch 46,000-325,000 copies (cache 200,000). Congener distances 9-15%, column weights 2-4%, alleles ~1%,
  everything else under 0.3%; the alleles' scores in the alignment cost nothing measurable.
- **So the ancestry sites are worth ~0.001-0.0015 of short-read F1 at roughly 10-15% of a complex short-read
  sample's CPU; for ONT 0.0019 at little cost (profiling is ~4% of its CPU); for PacBio nothing measurable.**
  Computing them in `--build` and storing them would keep the gain without the run-time cost (section 7).
- **`--weighted_site_shift` has never taken effect since it became opt-in** (`0b34254`): the per-thread copy of the
  alignment handler drops the flag. The earlier tests of the shift itself were valid; it needs retesting with the
  fixed column weights.
- **Follow-up the same day (section 8):** the shift fixed and retested on the small world (pe +0.0008 with species
  held out, +0.0027 on the test set; that world has no conserved column, so the r226 test is still to come); the
  ancestry sites computed by `--build` and stored (`ancestry_sites.bin`: profiles byte-identical, the profiling −47% of
  its instructions); `--insilico-strains auto` and `download_gtdb.py --species 40000 --per_species 4` for real strains
  at ~54% of the simulated strains (v23: 25%).

## 1. What changed from v22

**The inputs.** 47,403 genomes of 34,956 species (Archaea 1,698, Bacteria 33,258) against v22's 37,422 of 24,975:
9,981 more species, all with one genome. The species with real strains (10,227), the genomes left to the strain
alleles (12,173) and the alleles themselves (3.70 M of 29,170 species in the training database, v22 3.67 M of 29,019)
are about the same; the in-silico strains of one-genome species grew from 14,746 to 24,727. Of the simulated
species, 13.3% are real strains (v22 18.6%), 40.0% in-silico strains (36.0%), the rest representatives. The hold-out
follows the table: 12,142 species left out (10,620), of them 1,755 in the 12 genus clades (3,276) and 10,387 alone
(7,344); the scenarios' genome tables grew alike (moderate 34,544 species, 10,363 of them not in the training
database; v22 24,410 and 7,323). So v22 and v23 differ in their samples, their hold-out and their in-silico share as
well as in the code; section 2 separates the class mix from the error rates.

**The code** (`82bcc15..8d8c89b`):
- **Column weights fixed and taken from GTDB's alignments** (`0d961e1` pooled genus shares, `c60c272` protein
  alignment and the options read, `6f8e6bb` GTDB's MSA columns and tree): "0 genes without an alignment" (v22: 42% of
  the alignments failed), within-genus codes up to 14 with 2.80% of columns at 9 or more ("conserved"), among-genus
  codes up to 12 with 2.46% at 6 or more (v22: capped at 7 and 5; two features constant 0).
- **Ancestry sites from the sites the species' known genomes share** (`371ac84`), not the representative's alone.
- **The weighted site shift opt-in** (`0b34254`): in v22 it ran on every unsure read (with the capped weights); v23
  ran without it (and `--weighted_site_shift` could not have turned it on, section 6).
- **Index alleles default off** (`84af9a4`), but the build was given `--index-alleles 0.01`: 101.6 M seeds of 1.33 M
  alleles in the training database (3.7% of its 2.77 billion index entries), 114.8 M in `protal_db`, as in v22.

**The build:** ready after 2:05:44 (v22 2:12:24) with 40% more species simulated; training database 25:50 (21:31),
`protal_db` 1:27:35 in the background (1:28:58); `database.protal` 25.3 GB (25.2); parity OK.

## 2. The results

F1 at the knob protal calls at, with species held out (knob 0.5) and on the independent test set (pe knob 0.74 in
v23, 0.75 in v22; se 0.5 in v23, 0.65 in v22; long reads 0.5), from `build_metadata.tsv`:

| | pe held out | pe test | se held out | se test | pb held out | pb test | ont held out | ont test |
|---|---|---|---|---|---|---|---|---|
| v22 | 0.9573 | 0.9700 | 0.9627 | 0.9672 | 0.9757 | 0.9867 | 0.9704 | 0.9833 |
| v23 | 0.9595 | 0.9709 | 0.9650 | 0.9662 | 0.9791 | 0.9828 | 0.9741 | 0.9820 |
| change | +0.0022 | +0.0009 | +0.0023 | −0.0010 | +0.0034 | −0.0039 | +0.0037 | −0.0013 |

The test sets are small (pb 2,057 taxa, 22 FP and 22 FN; ont 2,591): their changes are within noise (the v22 report's
bootstrap intervals on these sets were ±0.002-0.005). Per scenario ([scenarios_v22_v23.txt](scenarios_v22_v23.txt))
the long reads gained most on the complex ones with species held out (pb soil 0.9696 → 0.9736, soil_shallow 0.9670 →
0.9718; ont moderate 0.9783 → 0.9829, soil_shallow 0.9607 → 0.9647); the 4-sample hold-out sets move both ways.

**This time it is the error rates, not the mix** ([builds_by_class.txt](builds_by_class.txt), each build's own
model, species held out). Each build's class error rates applied to the other's class mix:

| F1 from the class rates | pe v22 mix | pe v23 mix | se v22 | se v23 | pb v22 | pb v23 | ont v22 | ont v23 |
|---|---|---|---|---|---|---|---|---|
| v22 rates | 0.9573 | 0.9572 | 0.9627 | 0.9630 | 0.9757 | 0.9772 | 0.9704 | 0.9717 |
| v23 rates | 0.9594 | 0.9595 | 0.9643 | 0.9650 | 0.9774 | 0.9791 | 0.9728 | 0.9741 |

At the same mix v23's rates are worth +0.0021 (pe), +0.0016-0.0020 (se), +0.0017-0.0019 (pb), +0.0024 (ont); the mix
adds nothing for pe and +0.0003-0.0017 for the others (fewer real strains). By class (pe; the other read types alike):

| pe, species held out | v22 share | v23 share | v22 rate | v23 rate |
|---|---|---|---|---|
| FN real strain, alleles | 12.4% | 8.6% | 5.42% | 6.51% |
| FN real strain, none | 8.9% | 6.5% | 5.35% | 4.88% |
| FN in-silico strain | 48.6% | 53.8% | 2.78% | 2.75% |
| FN representative | 30.1% | 31.1% | 1.74% | 1.84% |
| FP beside a held-out congener | 69.2% | 69.6% | 1.59% | 1.41% |
| FP beside present congeners only | 12.4% | 11.8% | 0.62% | 0.61% |
| FP no congener in the sample | 18.4% | 18.6% | 1.28% | 1.22% |

The gain is the false positives beside a held-out congener: −11% (pe 1.59% → 1.41%), −15% (se 1.55% → 1.32%), −17%
(pb 2.15% → 1.79%), −17% (ont 2.20% → 1.82%). The real strains with alleles are missed more often in every read type
(pe 5.42% → 6.51%, se 5.70% → 6.40%, pb 6.19% → 6.74%, ont 7.14% → 7.31%); the strains without alleles and the
in-silico ones about as before. Section 4 says which groups carry the FP gain.

**The features on the hard rows** ([feature_auc.txt](feature_auc.txt): present real strains against absent species
beside a held-out congener; hard: strains with p < 0.9, absent with p > 0.1): the ancestry agreement is no longer
inverted (pe 0.471 → 0.526, pb 0.435 → 0.469), as the shared-sites change meant; `ancestry_fixed_agreement` 0.572
→ 0.613 (pe), `allele_explained_share` 0.593 → 0.627. The column weights now separate the easy rows
(`conserved_mismatch_ratio` over all rows 0.679 → 0.821, `nonsynonymous_share` 0.711 → 0.803) but not the hard ones
(0.47-0.50). Their share of the models' split gain doubled ([varimp_groups.txt](varimp_groups.txt): pe 0.38% →
0.83%, pb 0.51% → 1.59%, ont 0.46% → 1.84%).

**Other lines.**
- Composition (test, median over samples without a host): the explained share +0.003 → +0.013 (pe), the unknown share
  −0.005 → −0.022 (pe), −0.005 → −0.056 (se), the average genome size +4.4% → +7.3% (pe; se +6.7%, pb +6.9%, ont
  +8.3%): worse than v22 in every read type, likely from the 10,000 one-genome species added (not examined here).
  PacBio's explained share is still +0.196 (v22 +0.199), the open depth bias.
- The ancestry report's AUCs (errors' reads, FN against FP) rose: pe fixed sites 0.573 → 0.593, identity 0.518 →
  0.552; se fixed sites 0.595 → 0.615; long reads still below 0.53 (other taxa than v22's; not like for like).

## 3. What drives the false positives and negatives

From the build's own models with species held out at knob 0.5 ([error_anatomy.txt](error_anatomy.txt), and
`logs/classifier_training*.log`).

**The complex scenarios carry most of both.** The 30 scenario samples of moderate, soil and soil_shallow
communities (8% of the pe training samples) hold 86% of the pe false positives (3,474 of 4,049) and 86% of the false
negatives (1,931 of 2,255); the 309 design samples (20-200 species) have 532 FP and 285 FN between them. pb: 89%
of the FP and 90% of the FN, ont 90% and 92%; se, without soil_shallow, 74% and 78%. Their F1 with species held out:
pe soil_shallow 0.948, soil 0.956, moderate 0.972, gut 0.991.

**False positives are congeners of species the database lacks.** 77% of the pe FP (3,104) are absent species with a
held-out congener in the sample, and the rate falls steeply with the distance to it: 27% of the absent species within
0.02 (Mash, marker genes) of a held-out congener are called, 5.7% at 0.02-0.05, 0.8% beyond (pb 18% / 3.3% / 0.7%,
ont 20% / 4.2% / 0.8%). The held-out congener's reads land on its nearest database species, so these FP do not lack
evidence: the FP rate rises with the taxon's fragments (1.1% of absent taxa with one fragment, 3.2% with 100 or more).
The rest are cross-family and deeper (family 421, order 85, class 64, phylum 20).

**False negatives are thin taxa and strains beside congeners.** 40% of the pe FN have one fragment (907; miss rate
11.5%), 62% three or fewer, against a miss rate of 0.75% with 100 fragments or more. 85% of them (1,918) share their
genus with another species of the sample, and a present species within 0.02-0.05 of a held-out congener is missed
23% of the time. By genome: real strains are missed 5.8% of the time (10.1% with 1-10 fragments), in-silico strains
2.7%, representatives 1.8%. The missed strains' identity is lower (median 0.966 against 0.994 for those found): they
look like the absent species that hold a missing congener's reads, which no feature tells apart on the hard rows
(section 2).

So the remaining errors are where the simulation is hardest and the evidence thinnest: complex communities with
many species the database lacks, where a missing species' reads and a thinly covered strain look alike.

## 4. What each feature group adds to F1

Species held out, F1 at knob 0.5, the paired bootstrap's mean [2.5%, 97.5%] of the variant's F1 minus the full
set's (full set: pe 0.9596, se 0.9650, pb 0.9792, ont 0.9741; [ablation_pe.txt](ablation_pe.txt),
[ablation_se.txt](ablation_se.txt), [ablation_pb.txt](ablation_pb.txt), [ablation_ont.txt](ablation_ont.txt); tables by
[group_tables.py](group_tables.py)). One group left out at a time:

ont no_sites exit 0 209 s
ont no_anc_dist exit 0 262 s
EXTRA DONE
pb core exit 0 93 s
queue stopped after pb Sun Oct 11 10:44:12 BST 2026
| left out | pe | se | pb | ont |
|---|---|---|---|---|
| adjacency | +0.0000 [−0.0003, +0.0003] | +0.0000 [+0.0000, +0.0000] | −0.0002 [−0.0006, +0.0002] | +0.0002 [−0.0003, +0.0007] |
| distance | −0.0005 [−0.0009, −0.0001] | −0.0006 [−0.0012, +0.0001] | −0.0000 [−0.0003, +0.0003] | −0.0001 [−0.0006, +0.0004] |
| depth | +0.0002 [−0.0002, +0.0006] | −0.0000 [−0.0005, +0.0005] | +0.0000 [−0.0003, +0.0003] | +0.0002 [−0.0002, +0.0006] |
| divergence | −0.0009 [−0.0014, −0.0004] | +0.0001 [−0.0004, +0.0006] | +0.0002 [−0.0002, +0.0006] | +0.0000 [−0.0005, +0.0006] |
| unfiltered | +0.0001 [−0.0003, +0.0006] | +0.0001 [−0.0004, +0.0005] | +0.0000 [−0.0003, +0.0004] | +0.0000 [−0.0004, +0.0004] |
| ref | −0.0001 [−0.0005, +0.0003] | −0.0002 [−0.0007, +0.0003] | −0.0001 [−0.0005, +0.0004] | −0.0003 [−0.0009, +0.0002] |
| complexity | −0.0003 [−0.0009, +0.0002] | −0.0003 [−0.0010, +0.0002] | +0.0001 [−0.0003, +0.0004] | −0.0001 [−0.0006, +0.0004] |
| consistency | −0.0002 [−0.0005, +0.0002] | −0.0003 [−0.0008, +0.0003] | −0.0003 [−0.0006, −0.0000] | −0.0002 [−0.0007, +0.0001] |
| shape | −0.0003 [−0.0007, +0.0001] | −0.0003 [−0.0007, +0.0002] | −0.0003 [−0.0007, +0.0002] | +0.0003 [−0.0003, +0.0007] |
| neighbourhood | +0.0001 [−0.0003, +0.0005] | −0.0005 [−0.0010, −0.0001] | −0.0002 [−0.0005, +0.0001] | +0.0000 [−0.0004, +0.0004] |
| ancestry | +0.0002 [−0.0002, +0.0006] | +0.0000 [−0.0004, +0.0005] | −0.0002 [−0.0005, +0.0003] | −0.0003 [−0.0007, +0.0002] |
| gaps | +0.0003 [−0.0001, +0.0007] | −0.0001 [−0.0005, +0.0004] | +0.0001 [−0.0001, +0.0004] | +0.0000 [−0.0004, +0.0004] |
| untried | −0.0000 [−0.0005, +0.0004] | −0.0000 [−0.0005, +0.0004] | +0.0000 [+0.0000, +0.0000] | +0.0000 [+0.0000, +0.0000] |
| alleles | −0.0001 [−0.0004, +0.0002] | −0.0003 [−0.0007, +0.0001] | +0.0001 [−0.0003, +0.0005] | −0.0002 [−0.0006, +0.0001] |
| polymorphic | −0.0002 [−0.0007, +0.0003] | −0.0000 [−0.0005, +0.0005] | −0.0000 [−0.0004, +0.0004] | −0.0006 [−0.0010, −0.0002] |
| weights | −0.0001 [−0.0004, +0.0003] | −0.0001 [−0.0005, +0.0004] | −0.0003 [−0.0006, −0.0000] | +0.0000 [−0.0004, +0.0004] |

**Alone, almost every group is replaceable.** Left out by themselves only divergence and distance cost pe F1,
neighbourhood se (−0.0005; distance −0.0006 at the edge), polymorphic ONT (−0.0006), and for PacBio none
(consistency and weights −0.0003 at the edge); the others are within ±0.0005 and often slightly
positive, as overlapping groups cover for each other. For PacBio even divergence, 71% of the model's split gain
(`excess_scaled_median`), is replaced by its unscaled twin in normalized. So what a group
is worth shows only when its stand-ins go with it:

| left out | pe held out | pe test | se held out | se test | pb held out | pb test | ont held out | ont test |
|---|---|---|---|---|---|---|---|---|
| ancestry | +0.0002 [−0.0002, +0.0006] | +0.0016 [−0.0002, +0.0032] | +0.0000 [−0.0004, +0.0005] | +0.0013 [−0.0006, +0.0033] | −0.0002 [−0.0005, +0.0003] | −0.0008 [−0.0025, +0.0012] | −0.0003 [−0.0007, +0.0002] | +0.0005 [−0.0023, +0.0037] |
| ancestry, distance | −0.0004 [−0.0009, +0.0000] | +0.0007 [−0.0014, +0.0026] | −0.0007 [−0.0013, −0.0001] | +0.0004 [−0.0021, +0.0031] | +0.0001 [−0.0003, +0.0005] | −0.0020 [−0.0057, +0.0010] | −0.0002 [−0.0008, +0.0002] | −0.0038 [−0.0077, +0.0000] |
| the 8 ancestry-site columns only | −0.0015 [−0.0021, −0.0010] | −0.0010 [−0.0033, +0.0014] | −0.0008 [−0.0013, −0.0002] | +0.0001 [−0.0018, +0.0021] | −0.0002 [−0.0006, +0.0002] | −0.0020 [−0.0042, +0.0000] | −0.0018 [−0.0024, −0.0013] | −0.0021 [−0.0055, +0.0018] |
| ancestry, polymorphic, weights | −0.0022 [−0.0029, −0.0016] | −0.0006 [−0.0029, +0.0017] | −0.0014 [−0.0020, −0.0008] | −0.0011 [−0.0036, +0.0013] | −0.0013 [−0.0018, −0.0007] | −0.0012 [−0.0041, +0.0012] | −0.0029 [−0.0037, −0.0022] | −0.0033 [−0.0067, +0.0002] |
| ancestry, alleles, polymorphic, weights | −0.0026 [−0.0032, −0.0020] | −0.0013 [−0.0035, +0.0010] | −0.0021 [−0.0027, −0.0015] | −0.0014 [−0.0037, +0.0009] | −0.0020 [−0.0026, −0.0014] | −0.0028 [−0.0057, −0.0003] | −0.0032 [−0.0040, −0.0024] | −0.0045 [−0.0080, −0.0010] |
| the same + gaps, untried | −0.0030 [−0.0036, −0.0024] | −0.0023 [−0.0045, −0.0001] | −0.0027 [−0.0033, −0.0019] | −0.0010 [−0.0035, +0.0014] | −0.0020 [−0.0026, −0.0014] | −0.0012 [−0.0048, +0.0023] | −0.0032 [−0.0038, −0.0026] | −0.0037 [−0.0069, −0.0005] |
| all but the core set (normalized..unfiltered) | −0.0055 [−0.0065, −0.0046] | −0.0050 [−0.0086, −0.0016] | −0.0051 [−0.0060, −0.0042] | −0.0023 [−0.0051, +0.0006] | −0.0042 [−0.0050, −0.0033] | −0.0016 [−0.0049, +0.0018] | −0.0058 [−0.0068, −0.0048] | −0.0051 [−0.0103, −0.0012] |

**The core set** (normalized, adjacency, distance, depth, divergence, unfiltered: the features up to v9) is 0.0055
below the full set for pe (test −0.0050 [−0.0086, −0.0016]), 0.0051 for se, 0.0042 for pb (test −0.0016, within
noise) and 0.0058 for ONT (test −0.0051 [−0.0103, −0.0012]): everything
added since is worth about half a point of F1, for pe mostly fewer false positives beside held-out congeners (+648
without), for pb fewer missed strains (in-silico +179, real with alleles +62 without) and false positives (+209).

**The three groups built on the ancestry sites** (ancestry, polymorphic, weights) stand in for each other: each alone
is worth nothing (pe +0.0002, −0.0002, −0.0001), the three together pe −0.0022, se −0.0016, pb −0.0012/−0.0013, ont −0.0031 (table below), about the whole step from
v22. Blanking only the eight columns that need the sites
([blank_site_columns.py](blank_site_columns.py), [ablate_columns.sh](ablate_columns.sh): the ancestry group,
`ancestry_fixed_gain`, `ancestry_fixed_agreement`, `ancestry_agreement_weighted`; the groups' other columns kept)
costs, against the build's own models ([site_columns.txt](site_columns.txt)):

| against the build's models | pe | se | pb | ont |
|---|---|---|---|---|
| ancestry, polymorphic, weights left out: F1 with species held out | −0.0022 [−0.0029, −0.0015] | −0.0016 [−0.0021, −0.0009] | −0.0012 [−0.0018, −0.0006] | −0.0031 [−0.0039, −0.0023] |
| ancestry and distance left out | −0.0004 [−0.0009, +0.0001] | −0.0009 [−0.0014, −0.0003] | +0.0002 [−0.0004, +0.0007] | −0.0004 [−0.0009, +0.0001] |
| eight site columns blanked: F1 with species held out | −0.0015 [−0.0020, −0.0009] | −0.0009 [−0.0015, −0.0003] | −0.0001 [−0.0006, +0.0003] | −0.0019 [−0.0025, −0.0014] |
| eight site columns blanked: test set | −0.0017 [−0.0038, +0.0007] | −0.0001 [−0.0023, +0.0018] | −0.0003 [−0.0021, +0.0019] | −0.0014 [−0.0042, +0.0011] |
| eight site columns blanked: net FP beside held-out congeners | +173 | +79 | +7 | +130 |
| eight site columns blanked: net FN in-silico / real strains | +54 / −19 | +25 / −16 | +36 / −21 | +81 / −17 |

The site columns carry two thirds of the three groups' worth for pe (−0.0015 of −0.0022), about 60% for se and ont,
nothing measurable for PacBio (there the allele and column-weight columns carry the groups' −0.0012). The local full refits
equal the build's models (pe +0.0001, pb +0.0001 against them).

**What the site features prevent, by class** (pe, net errors without the three groups): +277 FP beside held-out
congeners, +46 other FP, +71 FN of in-silico strains, but −25 FN of real strains (pb: +78, +0, +93, −39). As in v22,
part of the strain benefit flatters the in-silico strains, which keep the representative's derived states; the FP
benefit is on absent real GTDB species and should carry over to real samples.

**The strain alleles** (alleles group) add nothing measurable alone (pe −0.0001, pb +0.0001); without the three site
groups they add pe −0.0004 and pb −0.0007 (no_strain4 against no_sites; for pb the missed real strains with alleles go from −8 to +45 against the full
set), so the site groups cover for them. The alleles also act through the alignment scores and through polymorphic, so their total is not separable here.

## 5. What each feature group costs

**Where the features are computed.** protal computes every feature whatever the model reads: the evidence is gated
only by whether the database has the group's table (`--no_site_weights` and `--no_gene_neighbours` work by leaving
theirs empty), and `TaxonFeatures` (`Profiler.h:2507`) builds all ~110 columns for every taxon. Most of the work is
per SAM record, in `RecordEvidenceCollector::NoteRecord` (`Profiler.h:3152`), on the SAM-parsing threads, inside the
"reading the SAM" timer; it runs for every best record, the ones the MAPQ and length filters drop included. The
ancestry sites and the column weights of a (species, gene) copy are built lazily on its first record in the run
(`GenomeLoader::AncestrySitesOf`, `ancestry::Cache`: up to 17 `CompareCopies` against congeners' copies and a
consensus) and kept in run-wide caches of 200,000 copies that clear completely when full (`AncestrySites.h:72`;
r226 has ~4.2 M copies). The congener distances sketch every reference of a taxon with a congener in the sample
(`SampleContext.h`, cache 5,000 sketches).

**At r226** ([profile_stages.txt](profile_stages.txt)), the SAM parsing threads used 14,479 of the build's 155,600
protal CPU seconds (9.3%): for short reads 0.14-0.31 of the alignment's CPU (pe scenario samples 5,524 s against
17,584 s; 79 s per sample), for long reads 0.01-0.04. Profiling matters for short reads and hardly for long ones.

**Locally** ([costs.tsv](costs.tsv), [cost_runs.sh](cost_runs.sh)): instructions of the profiling (`ProfileSample`,
without the start-up) of one sample per run, and of five in one run (the caches shared), with each group's evidence
skipped by the cost probe:

| skipped | pe, 1 M pairs | PacBio, 100 Mb | pe, 5 × 1 M pairs in one run |
|---|---|---|---|
| nothing (G instructions) | 6.48 | 1.47 | 20.47 |
| ancestry sites (and with them the fixed-site part of polymorphic) | −3.19 (49.2%) | −0.88 (59.5%) | −7.76 (37.9%) |
| congener distances (distance) | −0.61 (9.4%) | −0.23 (15.3%) | |
| column weights (weights) | −0.25 (3.9%) | −0.03 (2.2%) | |
| strain alleles (alleles, polymorphic) | −0.08 (1.2%) | −0.01 (0.7%) | |
| congener gaps (gaps) | −0.014 (0.2%) | −0.0003 | |
| gene neighbours (adjacency; `--no_gene_neighbours`) | −0.004 (0.06%) | −0.002 (0.1%) | |
| untried, consistency's congener fit and gene overlap | < 0.002 | < 0.001 | |
| all seven | −4.19 (64.6%) | −1.17 (79.2%) | −10.24 (50.0%) |

The rest (2.3 G for the pe sample) is the SAM parsing, the records added to the taxa, the variants, the read EM, the
scoring and the profile. The groups not in the table (normalized, depth, divergence, unfiltered, ref, complexity,
shape, neighbourhood) are counters on the same parse or O(taxa) at the end and cannot be skipped apart from it; the
callgrind profile puts all of them within that 2.3 G.

**The ancestry sites cost per copy, and r226 samples touch many copies.** The pe sample's kept records touch 2,642
copies (2,835 built, counting the copies only dropped records reach) and the PacBio sample's 775: ~1.2 M
instructions per copy in both. From the callgrind profile, a build costs ~0.28 M instructions (decoding the copy, its
unique 12-mers, the consensus) plus ~0.56 M per congener compared (`CompareCopies`: decoding the congener's copy,
sorting its 12-mers, chaining, comparing; 3,994 comparisons, 1.4 per build: this world lists 1.35 neighbours per
species). r226's training database lists 3.9 (507,287 neighbours of 131,472 species), so a build there costs
~2.5 M instructions. Five samples in one run save half of it here (the same 170 species recur).

At r226 ([copies_per_sample.txt](copies_per_sample.txt)) a design sample's kept records touch a median of 903
copies (90%: 4,969), a complex scenario sample's 46,232 (90%: 105,642, largest 324,552, beyond the cache), a
PacBio scenario sample's 14,812. At ~2.5 M instructions per copy that is ~2 G instructions (about a second) for a
design sample, ~115 G (some 40-50 CPU seconds) for a median complex pe sample, and ~800 G for the largest, where the
cache also clears and starts over within the sample. The build's runs profile hundreds of samples per process and
share the cache, so their 79 CPU seconds of parsing per pe scenario sample include less of it; a run of one complex
sample pays it whole, against ~330 CPU seconds of aligning for a 9.4 M-pair soil sample at r226
(`sc_soil_pe_p20000000_s_2`, 84 threads): the ancestry sites are then the largest part of short-read profiling and some 10-15% of the run. An estimate: protal logs neither the builds nor
their time (section 7).

**In the alignment** ([align_cost_runs.sh](align_cost_runs.sh), [align_costs.tsv](align_costs.tsv)): the strain
alleles' scores (`ScoreAlleles`, `ShiftOf` on every candidate on a copy with alleles) cost nothing measurable. A new
pe sample of the small world (200,000 pairs; 2,271 of its 18,836 aligned candidates on copies with alleles, 12%, about
the share of an r226 soil sample: 299,686 of 2.95 M) took 16.69 G instructions in all with them and 16.74 G with
`--no_allele_scores`, the alignment 13.19 G both times (the profiling differs, as the records do).

**In the index:** the index alleles (`--index-alleles 0.01`, a build option, not a model feature) add 101.6 M of
2.77 billion entries (3.7%) to the training database's index; v22 measured `database.protal` +1.1 GB and +2.8 GB
resident for the index, and the training database's build grew 12:37 → 21:31 with them and the column weights. Their
cost per read was not measured here (it needs a second r226 index); at most a few percent of the seeding, plus the
candidates the allele seeds add for strains' reads.

## 6. `--weighted_site_shift` has had no effect since it was made opt-in

Found while mapping the features to the code. Every alignment thread works on a copy of the run's
`SimpleAlignmentHandler` (`Classify.h:89` and `:514`), and its copy constructor (`AlignmentStrategy.h:272-286`) copies
`m_allele_scores` but not `m_weighted_shift`, which `SetWeightedShift` sets on the original only (`RunProtal.h:763`).
So the threads always run with the flag false, whatever the command line says. History:

- **Until `0b34254`** (the truth-test worlds of [2026-10-09-site-weighted-evidence-plan](../2026-10-09-site-weighted-evidence-plan/README.md)
  section 9, and the v22 build at `82bcc15`): there was no flag; the shift discounted whenever `column_weights.tsv`
  was loaded, which a run that profiles does at start-up, before it aligns (`RunProtal.h:3215`). The "−0.006 F1 on the
  all-variable world" was measured with a working shift; v22 ran with it on, but with the among-genus codes capped at
  5, which discounted every column (the v22 report, section 4).
- **Since `0b34254`** (off by default, `--weighted_site_shift` to turn it on; the option read from the command line
  only since `c60c272`): the flag never reaches the threads, so the shift has not run since, and no test with the
  fixed weights (`0d961e1`, `c60c272`, `6f8e6bb`) has been possible. v23 ran without it, as intended.

`build_gtdb_database.py` has no option that passes `--weighted_site_shift` to its protal runs either. Retesting it
needs the copy constructor fixed (one line: `m_weighted_shift(other.m_weighted_shift)`), that option in the build
script, and a run with and without (section 7).

## 7. What to do next

1. **Make the ancestry sites cheap rather than drop them.** They carry most of what v23 gained (section 4) and half
   of the profiling (section 5). In order of effort: log the builds per sample (count and time, next to the SAM-read
   line) so that the next r226 build measures them; build sites only for copies with a kept record (`NoteRecord`
   builds them for every best record, filtered or not; here 7% more copies than the kept ones); keep each copy's
   sorted unique 12-mers (`CompareCopies` sorts the congener's again for every species that lists it); replace the
   clear-all caches (200,000 copies, 5,000 sketches) by LRU ones. The full answer is to compute the sites in
   `--build` and store them (~4.2 M copies, a few hundred bytes each: ~1-2 GB before compression; ~2.5 M
   instructions each, about 10^13 in all, minutes on the build's 84 threads), which also makes them the same in every
   run whatever the order of the samples.
2. **Retest the weighted site shift** (section 6): fix the copy constructor, add `--weighted-site-shift` to
   `build_gtdb_database.py`, and compare on the small real world (three seeds, minutes each) before an r226 build.
3. **The composition got worse** (section 2: unknown share −0.022, average genome size +7.3% for pe): look at the
   genome sizes of the 10,000 one-genome species added.
4. **Real strains with alleles are missed more** (pe 5.42% → 6.51%): fewer real strains in the training mix (13.3%
   of simulated species, against 40% in-silico) and site features that help in-silico strains more than real ones.
   More real strains (`download_gtdb.py --species`) or fewer in-silico ones (`--insilico-strains`), as the build's
   warning says.
5. **The website** is not affected: the follow-up's new options are the build's (`--no_ancestry_sites`,
   `build_gtdb_database.py --weighted-site-shift`, `--insilico-strains auto`, `--min-real-strains`) and the download's
   defaults, documented in `docs/databases.md`; a run's options and outputs are unchanged.

## 8. Follow-up the same day: the weighted shift fixed, the ancestry sites prebuilt, more real strains

The user asked for the first three of section 7's items, and for fewer in-silico strains by default with at least half
of the simulated strains real. All in the working tree on `8d8c89b` (WSL build `~/bnew`: unit tests 498 passed, 3 skipped as
before; e2e 149 on a mini database rebuilt with the table (222 of 348 copies with sites); the build pipeline 6, its
functions 27, the collector 18, the downloads 9, the in-silico strains 13, the mini DB 19).

**The weighted shift.** `SimpleAlignmentHandler`'s copy constructor copies `m_weighted_shift`
(`AlignmentStrategy.h`); `ColumnWeights.TheShiftDiscountsDifferencesAtVariableColumns` now also scores a copy of the
handler, as every alignment thread does, and fails with the fix reverted. `build_gtdb_database.py --weighted-site-shift`
passes `--weighted_site_shift` to the collections' protal runs (`collect_training_data.py --weighted_site_shift`, in
the profiles' cache key when on) and says so in `build_metadata.tsv`. The retest on the small world: below.

**The ancestry sites prebuilt** (`ancestry_sites.bin`). `protal --build` computes every copy's sites after the
column weights (`Build.h WriteAncestrySites`: the species neighbours, the congener gaps and the column weights it has
just written set on the genome loader, `GenomeLoader::ComputeAncestrySites` per copy, gene by gene on all threads, the
same code as a run's, `ancestry::ComputeSites` factored out of the cache) and stores them in the database
(`ancestry::Table`: sorted keys and offsets, then per copy varints of the positions' steps, the bases four to a byte,
the indels; only copies with sites). A run loads the table with the other tables (`RunProtal.h LoadAncestrySites`)
and `AncestrySitesOf` decodes a copy's sites on its first record through the same cache; without the table (older
databases, `--build --no_ancestry_sites`) or with `--no_site_weights` (the build polarised them with the weights) it
computes them as before. Unit tests: the encoding round trip (packed and unpacked bases, indels both ways, damaged
entries refused) and a table written from `ComputeSites` read back equal to it, damaged files refused.

- **Exact.** The small world's training database rebuilt (`~/rw_ws_off_s1`): 9,098 of 10,589 copies with sites, 68.3
  sites per copy, 1.03 MB (~113 bytes a copy), 4 s on 4 threads. Two pe samples and one PacBio sample of the earlier
  build profiled against the earlier training database (sites computed) and against the new one (decoded): all 15
  profile files byte-identical ([prebuilt_sites_costs.tsv](prebuilt_sites_costs.tsv)).
- **Cost.** The 1 M-pair pe sample under callgrind: the profiling 6.47 → 3.45 G instructions (−47%), the whole run
  8.33 → 5.33 G (−36%); decoding costs ~0.16 G against the 3.19 G of computing (the cost probe's skip, section 5). Wall
  for the three samples on one thread 2.79 → 1.81 s, resident memory +3.6 MB.
- **At r226** (expected, not run): ~113 bytes for each copy with sites, ~1 GB in memory if most of the ~14.5 M copies
  have congeners, and some minutes of the build on 84 threads (~1.5 ms per copy and thread here, more where species
  have more congeners).

**Fewer in-silico strains, more real ones.** [strain_supply.py](strain_supply.py) runs `download_gtdb.py`'s own
`pick()` on GTDB r226's metadata (`~/GTDB/r226`, no download) and the build's arithmetic of the mix (the simulator draws
a species, then one of its genomes; the allele split leaves k of a species' P strains simulated, k ~ Binomial(P, 0.5))
([strain_supply.txt](strain_supply.txt)). It reproduces v23 from its options (16,000 species × 2 strains, 19,000
representative-only, the in-silico defaults of then): real strains 13.2% of the simulated species, in-silico 40.0%
(v23: 13.3% and 40.0%). GTDB r226 has 36,515 species with strains passing the CheckM2 filters; v23 used 16,000.

| download (`--species`, `--per_species`, `--rep_only_species`) | in-silico (`--insilico-multi`, `--insilico-strains`) | real | in-silico | representative | real of the strains |
|---|---|---|---|---|---|
| 16000, 2, 19000 (v23's) | 0.5, 1 (v23's) | 13.2% | 40.0% | 46.8% | 24.7% |
| 16000, 2, 19000 | 0.25, auto (new defaults) | 14.4% | 14.4% | 71.2% | 50.0% |
| 40000, 4, 9000 (new defaults) | 0.5, 1 | 28.1% | 30.0% | 41.9% | 48.3% |
| 40000, 4, 9000 (new defaults) | 0.25, auto (new defaults) | 30.4% | 25.9% | 43.7% | 54.0% |
| 40000, 4, 19000 | 0.25, auto | 24.9% | 24.9% | 50.1% | 50.0% |

- `build_gtdb_database.py --insilico-strains auto` (the new default; a share as before): after the genome table, as
  many one-genome species get an in-silico strain as keep the real strains at `--min-real-strains` (0.5) of the
  simulated strains or more (`auto_insilico_share`; on the small world's table it gives the earlier mix exactly at
  the old settings, and 0.76 of the one-genome species at the new ones); `--insilico-multi` 0.5 → 0.25.
- `download_gtdb.py --species 40000 --per_species 4` (were 16000 and 2): every r226 species with strains, up to 4
  each, ~80,700 strains against ~24,700. A rerun into the existing folder keeps what is there, gives strains to its
  representative-only species first and downloads only the new genomes: ~56,000 strains, 194 Gb of sequence by
  GTDB's genome sizes, ~60 GB gzipped (an estimate at ~0.3 bytes a base), and the representatives of species new to
  the folder. On the cluster:
  `python3 scripts/download_gtdb.py -o /hpc-home/hildebra/DB/protal/refGenos/ --release 226` (the new defaults), then
  the build as before.

**The weighted shift on the small world** ([weighted_shift_runs.sh](weighted_shift_runs.sh),
[compare_ws.py](compare_ws.py), [weighted_shift.txt](weighted_shift.txt)). The builds of ~/realworld_sp_s1..3 again with
the new protal and scripts, with `--weighted-site-shift`, the in-silico strains as then (`--insilico-strains 1
--insilico-multi 0.5`). Without the shift, seed 1 was built again and its training and test tables are byte-identical
to the earlier build's (every feature of every row, so the prebuilt sites change nothing end to end); for seeds 2 and 3
the earlier builds serve as the arm without (`~/rw_ws_off_s2`, `_s3` link to them). The shift now acts: 6-9 times as
many candidates get a shift (seed 1: 4.57 M against 0.74 M), as every candidate on a copy with column weights is
discounted, and of the unsure reads (the same 0.67 M) 48,462 move to another species against 53,578 (seeds 2 and 3: 42,197
against 45,920, 24,595 against 23,996).

| F1 at knob 0.5 | pe, species held out | pe, test | pb |
|---|---|---|---|
| without the shift | 0.9665 | 0.9713 | 0.9751 / 0.9775 |
| with | 0.9672 | 0.9741 | the same |
| with - without, paired bootstrap over samples | +0.0008 [−0.0003, +0.0019] (900 samples) | +0.0027 [+0.0005, +0.0052] (135) | 0 (short reads only) |

So on this world the shift no longer costs (the truth-test worlds of section 6 lost 0.006 with it) and gains a little
on the test set. But this world cannot test what the shift is for: its families have about four genera each, no column
reaches the among-genus code of 6 (largest 4), so every column is discounted alike and the shift weighs nothing by
conservation. The test that counts is an r226 build with `--weighted-site-shift` (2.46% of r226 columns at 6 or more),
against v23 or the next build without it.
