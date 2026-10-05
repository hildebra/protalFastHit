# r226 v10: the models without the species priors

Date: 2026-10-04. Build: `281a4ba` (v0.7.5, priors opt-in), `build_gtdb_database.py --inputs refGenos -t 64`,
seed 1, on the HPC (job 23918389). Data: `local/protal_r226_v10_{logs,training}.tgz`, unpacked in `local/v10`
(git-ignored). Compared with v9 (`26b065c`, same design with `+priors`, `local/v9`,
[2026-10-03-r226-v9-evaluation](../2026-10-03-r226-v9-evaluation/)).

## What v10 is

The same samples as v9: identical row counts (pe 70,732 / se 60,559 / pb 11,912 / ont 14,386 training rows,
pe 25,392 test rows), and the feature columns are identical apart from float noise (`excess_scaled_median`
differs in 61k rows, by more than 1% in 8 — the known nondeterminism of the conservation factors; `mean_ani`,
`depth_cv`, `identity` in a few hundred rows at the last digit). The only change is the feature set:
`normalized+adjacency+distance+depth+divergence+unfiltered` (49 features), without `+priors`. v10 is therefore a
clean ablation of the priors at GTDB scale, with the 0.7.5 binaries.

## Headline (knob 0.5)

| read type | v9 held out F1 | v10 held out F1 | v9 test F1 | v10 test F1 | v9 test FP/sample | v10 test FP/sample | v10 best test threshold |
|---|---|---|---|---|---|---|---|
| pe  | 0.9695 | 0.9643 | 0.9699 | 0.9630 | 1.73 | 2.15 | 0.577 (F1 0.9638) |
| se  | 0.9665 | 0.9616 | 0.9668 | 0.9596 | 1.86 | 2.26 | 0.562 (F1 0.9611) |
| pb  | 0.9730 | 0.9692 | 0.9727 | 0.9675 | 0.76 | 0.90 | 0.412 (F1 0.9710) |
| ont | 0.9724 | 0.9685 | 0.9741 | 0.9689 | 0.98 | 1.10 | 0.545 (F1 0.9699) |

The v9 ablation predicted ~0.963 pe / ~0.960 se on these samples; v10 reads 0.9630 / 0.9596. Prediction
confirmed; the in-report feature-set study on v10's own tables gives the same priors gain (pe 49 → 56
features: 0.9643 → 0.9691 species held out; se 0.9616 → 0.9667; pb 0.9692 → 0.9735; ont 0.9685 → 0.9722).

PMML parity holds for all four models (probabilities equal, features equal). Training took 41-133 s per model.

## Where the priors' loss lands (test set, pe; se alike)

Errors at 0.5 by GTDB cluster size of the taxon (singleton = one-genome species) and fragments:

| cluster | truth | fragments | taxa | errors v9 | errors v10 |
|---|---|---|---|---|---|
| singleton | absent | ≤2 | 8,840 | 18 | 37 |
| singleton | absent | 3-10 | 2,298 | 5 | 15 |
| singleton | absent | >10 | 1,051 | 1 | 18 |
| multi | absent | all | 9,123 | 111 | 98 |
| multi | present | ≤2 | 702 | 67 | 76 |
| multi | present | 3-10 | 671 | 15 | 20 |
| multi | present | >10 | 2,005 | 19 | 31 |
| singleton | present | all | 702 | 10 | 8 |

(se: singleton FP 33 → 74, multi-genome FN 108 → 139.) As in v9's ablation, the priors' gain is the cluster-size
rule: absent singletons are now called 2.9× as often, and present multi-genome strains lose their "this species
has divergent strains" prior. Present singletons are no worse without priors. Whether the rule is real biology or
a simulation artefact (strains are only simulated for multi-genome species, so present taxa are 17% singletons
vs 56% of absent ones) cannot be settled in simulation.

## Strengths

1. **No memorisation.** Rows vs species held out: F1 0.9642 vs 0.9643; held out by genus up to phylum:
   0.9635-0.9643, FP and FN rates flat across ranks. The features describe alignments, not taxa, so the model
   transfers to clades the training never saw — the property that matters for GTDB's 143k species.
2. **Held out and test agree** (pe 0.9643 vs 0.9630, ont 0.9685 vs 0.9689): the test design (σ 2.0, other
   depths) is not much harder for the model; no overfit to the training design.
3. **Deep samples are nearly solved.** pe test from 200k pairs: FP rate 0.78% → 0.14% at 5M pairs; FN 1.7-4%.
   With ≥100 fragments representatives are missed 0/481 (test), strains 1.3%.
4. **Archaea** sensitivity 0.983 held out, 0.978 test (pe); 1.0 for ont test. The domain skew of old training
   sets is gone.
5. **Feature additions pull their weight** (pe, species held out): depth +0.007, divergence +0.002, relatives
   (instead of distance) +0.004 on the reduced set. `unfiltered` is neutral for pe (+0.0002) but +0.002 se,
   +0.0016 pb, so keeping it is right.
6. **Forest size is right**: 64 → 256 trees +0.0006 F1; 512 leaves bind but unlimited leaves are worse.
7. **Clade-held-out FP are rare**: absent taxa near a species held out at genus/family/order level are called
   0.1-0.55%; the classic "relative of a missing genus" FP is under control.

## Weaknesses

1. **The congener of a missing species is the FP.** 359 of 430 pe FP (83%) are absent taxa whose closest
   simulated species is one the database lacks at species level (FP rate 0.96%, 5.8 FP per 100 such species);
   320 share the genus. Top-scored FP have `top_identity` 1.0 and up to 237 fragments (e.g. *Spiroplasma citri*,
   *Aristaeella* sp900320755 repeatedly) — near-identical marker genes between GTDB sibling species, likely not
   separable by these 168 genes at all. This is the precision ceiling of the marker set, not of the model.
2. **Divergent strains are the FN.** Strains (another genome than the representative) are missed 4.6% vs
   0.9% for representatives; with 1-10 fragments 9.3%. Missed strains' identity median 0.967 vs 0.987 for
   found ones — they sit at the species boundary (~96% ANI), where their reads look like a relative's
   (`conserved_fast_kept_ratio` −1.23 missed vs −0.19 found). 53% of present taxa are strains, so this is the
   largest FN pool. Without priors FN of multi-genome species with >10 fragments rose 19 → 31.
3. **Shallow samples are noisy, and the test depths fall between the training points.** pe test FP rate of
   absent taxa: 11% at 500 pairs, 15% at 2k, 7% at 10k and 50k (71 FP at 50k = 6 per sample, 42% of all test
   FP); FN rate peaks at 50k (5.2%). Training has points at 20k and 100k, none at 50k; 2k and 500 are also off
   grid. se is worse (21% FP rate at 2k). The model is fine where it was trained densely and weakest exactly in
   the gaps.
4. **Calibration drifted with the priors gone.** Best threshold moved 0.50 (v9) → 0.577 pe, 0.593 se held
   out (the same on test); calling at 0.58 would gain +0.0024 held out but only +0.0008 test — below the 0.002
   rule, so no knob change, but the default now over-calls short reads. pb under-calls for the fourth build
   running (best 0.468 held out, 0.412 test, +0.0035 test F1).
5. **Learning curve flat for pe/pb/ont**: 25% of the samples is within 0.001 F1 of 100%. More samples of the
   same design do not help; only se's log loss still falls (0.0423 → 0.0375 from half to all). The gains have
   to come from different data (strains of singletons, depths in the gaps), not more of it.
6. **Small long-read sets**: pb test 42 samples / 3,034 taxa, ont 41 / 3,848. The pb/ont differences between
   builds (±0.005) are within the noise of ~100 errors; the shallowest long-read point (150 kb) has 13-20
   present taxa. Long-read conclusions from one build are weak.
7. **Archaea FP at shallow depth** (held out): 9 FP of 335 absent at 500k, 7 of 61 at 20k pairs — higher rates
   than Bacteria at the same depths, from few rows.

## Recommendations

1. **Decide the priors on real data, not on another simulation.** Simulation cannot separate the rule from the
   artefact. Profile real gut and a MAG-heavy environment (soil/marine) with v9 (`+priors`) and v10 models and
   compare calls of singleton species: if v9 drops credible singleton MAG species that v10 keeps, stay with v10.
2. **Make the singleton rule learnable**: simulate in-silico strains of one-genome species (mutated
   representative, codon-aware, 95-99% ANI) so present singletons carry strain divergence too. Then the
   cluster-size prior either loses its value (it was an artefact) or keeps it for a real reason, and the strain
   FN (weakness 2) gets training rows at the species boundary.
3. **Fill the depth grid**: add 2k and 50k pairs to the training read-pair design (and 200k), at a fraction of
   the deep points' cost; the test's worst bins are exactly those.
4. **Threshold**: leave pe/se at 0.5 for now (gain < 0.002 on test); revisit after (1). Consider the pb knob at
   0.45 if the fifth build again shows < 0.5 with ≥ +0.002 test.
5. **Try `relatives` with the current set without priors** (`normalized+adjacency+relatives+depth+divergence+unfiltered`):
   on the reduced set relatives beats distance by +0.004 pe / +0.006 se, but the full combination was only
   evaluated with priors (0.9706 vs 0.9691 pe). A cheap retrain on the v10 tables answers it.
6. **The FP ceiling** of near-identical sibling species (weakness 1) needs either more genes for those genera
   (gene-subset work in reverse) or reporting such calls at genus level; not a model problem.

## Follow-up the same day: recommendations 2, 3 and 5

**5, relatives instead of distance, without priors** ([`scripts/relatives_without_priors.sh`](scripts/relatives_without_priors.sh): the trainer on the
v10 tables in WSL, the build's options, `--evaluation basic`; the distance runs reproduce the HPC models' numbers
exactly):

| model | held out F1 | FP | FN | log loss | test F1 | FP | FN | log loss |
|---|---|---|---|---|---|---|---|---|
| pe distance (default) | 0.9643 | 430 | 284 | 0.0309 | 0.9630 | 168 | 135 | 0.0358 |
| pe relatives | 0.9654 | 404 | 286 | 0.0298 | 0.9633 | 161 | 139 | 0.0355 |
| se distance | 0.9616 | 461 | 292 | 0.0375 | 0.9596 | 176 | 148 | 0.0436 |
| se relatives | 0.9636 | 434 | 280 | 0.0363 | 0.9603 | 172 | 146 | 0.0430 |
| pb distance | 0.9692 | 169 | 235 | 0.0948 | 0.9675 | 38 | 62 | 0.0895 |
| pb relatives | 0.9684 | 174 | 241 | 0.0946 | 0.9688 | 33 | 63 | 0.0887 |
| ont distance | 0.9685 | 206 | 219 | 0.0851 | 0.9689 | 45 | 53 | 0.0786 |
| ont relatives | 0.9694 | 203 | 211 | 0.0855 | 0.9689 | 45 | 53 | 0.0781 |

The reduced set's +0.004 does not survive the depth, divergence and unfiltered features: test F1 +0.0003 pe,
+0.0007 se, +0.0013 pb, 0 ont, all below the 0.002 rule; held out +0.001-0.002 with 5-6% fewer FP for short
reads, the same story as v5/v6 (lower log loss, no test gain). **The default stays distance**; the result is
noted in `scripts/model_features.py`.

**3, the depth grid**: `build_gtdb_database.py --read-pairs` default now
`1000,2000,5000,20000,50000,100000,200000,500000,2000000:4,10000000:2,30000000:1` (three more points of 12
samples per read setup, 108 shallow samples, ~25M read pairs in all, against ~175M of the deep points). The test
design is unchanged, so its 2k, 50k and 200k points are no longer between training points: the next build's
test numbers there are not comparable with v10's as a measure of interpolation, only of the errors.

**2, in-silico strains of one-genome species**: `scripts/insilico_strains.py`, run by `build_gtdb_database.py`
after the conversion (`--insilico-strains 1`, `--insilico-ani`; docs/databases.md):
- Every species of the genome table with one genome (2,017 of 7,998 in the r226 pool) gets a copy of its
  representative, `insilico_<accession>`, which the simulator draws as another genome; the training and test
  collections simulate from `genomes_simulated.tsv`.
- Divergence: not a fixed 95-99% ANI as proposed above, but drawn from the table's real strains (their median
  marker gene's divergence from the k-mer traces of `gene_positions.tsv`, 1 - share^(1/24)), so that a one-genome
  species' strain is as far from its reference as a multi-genome species' strains are from theirs. A 95-99% range
  would put the in-silico strains further out than real ones (the r226 strains' cluster mean ANI is 99.15%
  median, 97.3% at the 10th percentile) and teach the opposite artefact. The genome-wide divergence is the
  marker divergence / 0.45 (v10's present strains with ≥300 fragments: identity loss over 1 - cluster mean ANI,
  median 0.43), capped at 5% (95% ANI). `--insilico-ani 95-99` gives the uniform draw.
- Genes: each placed marker gene by its conservation factor, estimated from the same placements as protal
  estimates them (gene over the strain's median gene, median over strains, median gene 1, 0.25-4), so the
  divergence features (`excess_scaled_median`, `excess_conserved_fast_ratio`) see a strain's pattern, not a
  uniform one. Other open reading frames (stop to stop, ≥300 bases, longest first) at a gamma rate (shape 2),
  non-coding 1.3x.
- Codons: synonymous changes kept, non-synonymous with probability 0.15, stop codons never (also when two
  changes meet in one codon); ts/tv 3. On a synthetic genome 60%+ of the codon changes fall on third positions,
  no stop codon is made, a 3% target gives 2.98%. ~10 s per 4 Mb genome on one core (~5 min for 2,000 on 64).
- The collector marks such taxa (`meta_insilico_strain`), and the trainer's report lists "another genome,
  real" and "another genome, in silico" apart (both strain tables): the check that the mutations look like a
  strain's is that both are missed alike at the same fragments.
- Substitutions only: no indels, rearrangements or gene content differences, which real strains have; the
  marker genes, which carry the features, are what the model of mutation is fitted to.

**Fix after v11's first start (2026-10-05):** the strain step stopped with `IndexError: index 71 is out of bounds` (job
23954360). Two substitutions in one codon were grouped by the codon's first base only, and where a minus- and a
plus-strand frame meet, both have a codon starting at one base: their changes were summed into one codon index.
The synthetic test genome had no such overlap; real genomes do. Codons are now grouped by first base and strand.
A regression test (`test_codons_of_both_strands_at_one_base`) fails with the old code and passes. A stress
run on dense genes overlapping on both strands, with N runs and tiny contigs, at 95% ANI also passes. The
strains' contigs are now named after the strain (`<strain>_<contig>`), so their reads are told from the
representative's ([2026-10-05-trace-relatives.md](../2026-10-05-trace-relatives.md)).

What the next r226 build (v11) should show: present one-genome species in the test set simulated from in-silico
strains about half the time; their FN rate by fragments near that of real strains; and, the point of it, the
priors' gain (`+priors` in the feature-set study) shrinking from +0.005 held out toward 0 if the cluster-size
rule was the simulation's. Its test set is not comparable with v10's (other depths in training, in-silico strains
in the test).

## How the numbers were made

Reports: `local/v10/model_logs/summary.txt`, `trained_model*.report.txt`, `parity*.txt`, `build_metadata.tsv`.
v9 vs v10 table identity and the error table by cluster size: Perl one-offs joining
`v{9,10}/model_logs/trained_model{,_se}.test_predictions.tsv.gz` with `v10/test/training_data{,_se}.tsv` on
(`meta_sample`, `taxon`), scoring `p >= 0.5` against `truth`, grouped by `cluster_genomes_log10 == 0` and
`fragments` (≤2, 3-10, >10).
