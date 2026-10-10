# A small real world: do the ancestry features separate congeners?

**Data.**
- **The world:** the six real GTDB r226 families of [2026-10-10-real-ancestry](../2026-10-10-real-ancestry/README.md)
  (Cyclobacteriaceae, Chitinophagaceae, Ruminococcaceae, Azotimanducaceae, Acidimicrobiaceae and the archaeal
  Micrarchaeaceae): 24 genera (4 per family, 4-12 species each), 170 species (147 bacterial, 23 archaeal), and their
  whole genomes from NCBI: 170 representatives and 360 strains of 88 species (354 of them MAGs; most of these genera
  have no isolate). [make_inputs.sh](make_inputs.sh) assembles the inputs folder (WSL `~/realworld`) from the two
  subsets and fetches the genomes with `download_gtdb.py`'s own genome step ([fetch_genomes.py](fetch_genomes.py)).
- **The builds:** `build_gtdb_database.py` on it ([run_build.sh](run_build.sh)), protal, simulator and scripts at
  `6f8e6bb` (WSL `~/bidx`, 4 cores), three seeds ([run_all.sh](run_all.sh); WSL `~/realworld_s1..3`), each 13-15 min:
  - 30% of the species and 3 genera left out of the training database (64-70 species; other ones per seed);
  - strain alleles from half of each species' genomes (`--allele-genome-share 0.5`), in-silico strains for the
    one-genome species and for half of those with strains (135; 23% of the simulated species are real strains, 35%
    in-silico);
  - 300 paired-end samples (one setup, 2k to 1M pairs) and 72 HiFi samples (1-100 Mb), 15-60 species each, half of
    them in groups of 2-5 congeners, 2 archaeal species per sample; an independent test set of 45 pe and 24 HiFi
    samples; no scenarios;
  - the column weights from GTDB's alignments and tree (`--column_weights_alignment auto` takes them).

  Summaries: [builds.txt](builds.txt) ([collect_builds.sh](collect_builds.sh)).
- **The refits:** each build's pe and pb models refitted with the build's trainer and options without one group of
  features ([ablate.sh](ablate.sh), [ablate_all.sh](ablate_all.sh)): the ancestry sites (`ancestry`, 5), the
  polymorphic group with the fixed-site features (`polymorphic`, 4), the column weights (`weights`, 5), all three (14),
  and for comparison the strain alleles (`alleles`, 2).
- **The comparisons:** [variants.py](variants.py) pools the three builds (each build's samples apart): F1 at knob 0.5,
  a paired bootstrap over the samples, the errors by class and the calls that flip
  ([variants_pe.txt](variants_pe.txt), [variants_pb.txt](variants_pb.txt)); `2026-10-10-r226-v22/feature_auc.py` on
  each build's training rows ([feature_auc_pe_s1.txt](feature_auc_pe_s1.txt), `_s2`, `_s3`).

**Question (the user's, 2026-10-10).** In a small world like the true-positive test's, but of real genomes, genera
and families: how effective are the ancestry-related features at separating congeners, i.e. at reducing false
positives and false negatives?

## Summary

- **They are not, in this world.** Leaving out any group of the ancestral-state features, or all of them, changes F1
  by less than the noise in either read type (species held out: pe −0.0001 to +0.0006, pb +0.0003 to +0.0009; the
  95% intervals all span 0). The strain alleles are no different.
- **The errors they move cancel.** Without all three groups the pe model makes 12 fewer false positives beside
  held-out congeners (55 fixed, 43 new) and 12 more false negatives of real strains with alleles (8 fixed, 20 new);
  without the ancestry sites alone 12 fewer of those FP and as many FN. 554 FP and 606 FN in all.
- **Why:** on the hard rows (strains with p < 0.9, absent species beside a held-out congener with p > 0.1, 85-111
  strains and 200-289 absent per build) the ancestry agreement points the wrong way (AUC 0.34-0.43): a real strain
  that diverged before the representative's own mutations carries the congener's base at those "species" sites,
  while a close novel congener shares much of the species' lineage. The fixed-site agreement (0.52-0.69) and the
  alleles (0.56-0.63) are the only ones above 0.5, and weak. The same holds at r226 (v22: 0.43-0.47).
- **The column weights cannot be judged here:** with 4 genera per family no column reaches "conserved" (largest code
  6), so two of the five weights features are 0 everywhere, as before.

## 1. The builds

[builds.txt](builds.txt). With species held out, at knob 0.5:

| seed | species left out (genera) | pe FP / FN | pe F1 | pb FP / FN | pb F1 |
|---|---|---|---|---|---|
| 1 | 64 (3, 18 species) | 165 / 204 | 0.970 | 14 / 30 | 0.980 |
| 2 | 70 (3, 27 species) | 189 / 186 | 0.966 | 19 / 29 | 0.975 |
| 3 | 64 (3, 18 species; other genera than seed 1) | 200 / 216 | 0.965 | 28 / 31 | 0.972 |

A first build with 48 pe and 32 HiFi samples had 29 FP and 34 FN (pe): too few to see a group's effect, hence the
six-fold samples and three seeds.

## 2. The feature groups

Pooled over the three builds, species held out; the paired bootstrap's mean [2.5%, 97.5%] of the variant's F1 minus
the full set's:

| left out | pe | pb |
|---|---|---|
| `ancestry` | +0.0003 [−0.0004, +0.0009] | +0.0009 [−0.0003, +0.0022] |
| `polymorphic` (fixed sites) | +0.0003 [−0.0004, +0.0010] | +0.0005 [−0.0008, +0.0018] |
| `weights` | +0.0006 [−0.0001, +0.0013] | +0.0006 [−0.0006, +0.0019] |
| all three | −0.0001 [−0.0009, +0.0007] | +0.0003 [−0.0013, +0.0020] |
| `alleles` | −0.0002 [−0.0007, +0.0005] | +0.0001 [−0.0011, +0.0013] |

On the test sets all within noise (pe −0.0007 to +0.0008; pb −0.0023 to −0.0006, each interval reaching 0).

**The errors by class** (pe, species held out; error rate and count, full set → without all three groups):

| class | n | full | without all three |
|---|---|---|---|
| FN real strain, alleles | 3,552 | 4.50% (160) | 4.84% (172) |
| FN real strain, no alleles | 1,243 | 3.54% (44) | 3.62% (45) |
| FN in-silico strain | 9,146 | 3.39% (310) | 3.40% (311) |
| FN representative | 3,523 | 2.61% (92) | 2.72% (96) |
| FP beside a held-out congener | 7,381 | 6.12% (452) | 5.96% (440) |
| FP beside present congeners only | 3,670 | 2.70% (99) | 2.56% (94) |

**The calls that flip** (pe; fixed / new errors of the variant): without all three groups FP beside a held-out
congener 55 / 43, FN real strains with alleles 8 / 20; without the ancestry sites 40 / 28 and 14 / 12. pb: a few
calls either way.

## 3. Each feature alone

`feature_auc.py`, present real strains against absent species beside a held-out congener, all rows / the hard rows
(seeds 1, 2, 3):

| | seed 1 | seed 2 | seed 3 |
|---|---|---|---|
| identity | 0.977 / 0.382 | 0.971 / 0.400 | 0.963 / 0.442 |
| ancestry_agreement | 0.920 / 0.342 | 0.915 / 0.412 | 0.915 / 0.434 |
| ancestry_congener_share (lower for strains) | 0.114 / 0.638 | 0.109 / 0.571 | 0.111 / 0.574 |
| ancestry_fixed_agreement | 0.761 / 0.601 | 0.745 / 0.517 | 0.804 / 0.688 |
| ancestry_fixed_gain | 0.585 / 0.564 | 0.607 / 0.479 | 0.620 / 0.599 |
| allele_explained_share | 0.679 / 0.597 | 0.691 / 0.561 | 0.725 / 0.628 |
| conserved_mismatch_ratio | 0.840 / 0.456 | 0.824 / 0.439 | 0.834 / 0.523 |
| ancestry_agreement_weighted | 0.916 / 0.337 | 0.911 / 0.438 | 0.913 / 0.437 |

The ancestry features are known for 91-99% of the rows. Over all rows the agreement is strong (0.92) but below
identity (0.96-0.98), so the models gain nothing from it; on the hard rows, where the errors are, it is inverted.
Real genomes say why ([2026-10-10-real-ancestry](../2026-10-10-real-ancestry/README.md) section 5): the sites come from
one representative, a quarter of real strains carry the species' base at under 88.5% of them, and only the fixed
sites (shared by the species' other strains) separate cleanly.

## 4. What to do next

1. **Define the ancestry sites from the species' known genomes**, not its representative alone: a site only where the
   species' alleles (the other genomes of the full reference) share the derived base. The agreement then measures the
   species, not one genome. Done in `371ac84` and tested on this world the same evening (section 5): the feature
   improves a little, the calls do not change.
2. **Judge the column weights on whole families** (r226, or this world with every genus of the six families): here
   no column is conserved.
3. **Train on more real strains:** 23% of the simulated species are real strains here (35% in-silico), and the
   ancestry features' gain at r226 v22 came through in-silico strains, which keep the representative's sites.

## 5. Follow-up: the sites the species' genomes share (`371ac84`)

The same evening the ancestry features were changed to count only the sites the species' known genomes share
(`AncestrySites.h SharedBySpecies`: a site goes where one of the copy's strain alleles carries another base, or an
indel near; a copy without alleles keeps all its sites; the fixed-site features keep all the sites). The world was
built again at `371ac84` with the same three seeds ([run_all.sh](run_all.sh) `sp_`: WSL `~/realworld_sp_s1..3`), so the
samples, hold-outs and candidate taxa are the first run's and only the features differ; [compare_rerun.sh](compare_rerun.sh)
pairs the rerun's full models with the first run's ([rerun_vs_before_pe.txt](rerun_vs_before_pe.txt),
[rerun_vs_before_pb.txt](rerun_vs_before_pb.txt)), refits the variants ([rerun_variants_pe.txt](rerun_variants_pe.txt),
[rerun_variants_pb.txt](rerun_variants_pb.txt)) and takes each feature alone ([rerun_feature_auc_pe_s1.txt](rerun_feature_auc_pe_s1.txt),
`_s2`, `_s3`).

| species held out | pe F1 (FP / FN) | pb F1 (FP / FN) |
|---|---|---|
| first run, the representative's sites | 0.9667 (554 / 606) | 0.9757 (61 / 90) |
| rerun, the species' shared sites | 0.9665 (567 / 603) | 0.9751 (63 / 92) |
| paired difference | −0.0003 [−0.0010, +0.0005] | −0.0007 [−0.0020, +0.0008] |

- **The calls do not change** beyond noise; by class the counts move by a few (pe FP beside held-out congeners 452 →
  462, FN of real strains with alleles 160 → 160).
- **The feature improves a little and is still inverted on the hard rows:** `ancestry_agreement` all rows 0.920 /
  0.915 / 0.915 → 0.926 / 0.927 / 0.919, hard rows 0.342 / 0.412 / 0.434 → 0.392 / 0.467 / 0.476 (seeds 1-3).
- **Within the rerun the groups still add nothing:** without the ancestry features +0.0002 [−0.0005, +0.0008], without
  all three +0.0002 [−0.0007, +0.0011] (pe, species held out).
- **Why so little:** the filter needs the species' other genomes as alleles, and here 40 species of the training
  database have any (4,138 alleles; the other genomes are split between the alleles and the simulation, and most
  species have one genome). For the rest the sites are the representative's as before. A strain of a species without
  alleles, or one that diverged where no allele genome covers, still carries the congeners' base at the
  representative's private mutations; and the hard absent rows are close novel congeners, which share most of the
  species' lineage whatever the sites.

So, in this world, neither the representative's nor the species' shared sites make the ancestry features separate
congeners where the models need it. The fixed sites and the alleles carry what little there is (hard-row AUC 0.50-0.69).
Keeping the change costs nothing (the calls are the same) and makes the feature mean what it says where alleles exist.
