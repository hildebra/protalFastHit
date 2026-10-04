# The r226 v9 build: round two at GTDB scale

2026-10-03. Data: a SLURM run of `build_gtdb_database.py` on GTDB r226 with the binary of `26b065c` (protal v0.7.4,
round two of the false-positive work, [report](../2026-10-03-false-positive-fixes/README.md)) and the scripts of
`ec1db6b` (the same build scripts), job 23905310, 52 threads, `--inputs refGenos --scratch` on the node's SSD and
otherwise the defaults; uploaded as `local/protal_r226_v9_{logs,training}.tgz` (git-ignored), extracted to
`local/v9`. Against v8 ([v7/v8 report](../2026-10-03-r226-v7-v8-evaluation/README.md)) the run changed five things at
once, all defaults of round two:

| | v8 (`a3e397d`) | v9 (`26b065c`) |
|---|---|---|
| features | `normalized+adjacency+distance+depth+divergence` (46 columns) | `+unfiltered+priors` (56): `fragments_all`, `em_fragments`, `failed_candidate_rate`; `rep_*`, `cluster_*` |
| training abundances | lognormal σ 1.3 | σ 1.3 and 2.0, the samples of a design point in turn |
| species held out of the training database | 20% (7,899) | 30% (8,668) |
| deepest paired-end point | 10M pairs, 2 samples per setup | + 30M pairs, 1 per setup (3 samples) |
| the database | | `species_priors.tsv` (the `sp_clusters_r226.tsv` fetched: cluster values for all 143,614 species, CheckM for all, 0 representatives with a duplicated marker), the same 5,589 suspect copies |

The test design is the one of v7 and v8 (σ 2.0, 500 to 5M pairs, 78 paired-end samples), but with 30% of the
species held out the test set has fewer present and more absent taxa, more of them relatives of a species the
database lacks (pe: 4,080 present and 21,312 absent against 4,591 and 19,038; 15,677 of the absent taxa beside a
missing species against 12,505). So the comparison is by rates, not by F1 alone, and the feature groups are judged
by an ablation on v9's own tables: the trainer of the same commit on `local/v9/training`, scored on `local/v9/test`
([`ablation_v9.sh`](ablation_v9.sh), [`ablation_v9_summary.py`](ablation_v9_summary.py) →
[`ablation_v9_summary_output.txt`](ablation_v9_summary_output.txt), [`ablation_by_cluster.py`](ablation_by_cluster.py) →
[`ablation_by_cluster_output.txt`](ablation_by_cluster_output.txt), with the two variants
[`priors_without_singleton_flag.sh`](priors_without_singleton_flag.sh) and [`priors_checkm_only.sh`](priors_checkm_only.sh);
the default set reproduces the pipeline's model exactly). The v8/v9 comparison by class is
[`compare_v8_v9.py`](compare_v8_v9.py) → [`compare_v8_v9_output.txt`](compare_v8_v9_output.txt), the cluster-size
analysis [`priors_bias.py`](priors_bias.py) → [`priors_bias_output.txt`](priors_bias_output.txt). Parity passes for all
four models; the run took no longer than v8 (from the end of the conversion to the summary 2:02 against 2:17; the
index build 33.5 min, of which the suspect-copy scan 11.7).

## Summary

Independent test set at knob 0.5, v8 → v9:

| read type | F1 | sensitivity | FP per 100 absent taxa | FP per 100 absent beside a missing species | misses / false positives | best test threshold | held out F1 |
|---|---|---|---|---|---|---|---|
| pe | 0.9671 → **0.9699** | 96.2 → 97.3% | 0.67 → 0.63 | 0.83 → 0.68 | 173 / 128 → 111 / 135 | 0.47 → **0.50** | 0.9747 → 0.9695 |
| se | 0.9685 → 0.9668 | 96.3 → 97.0% | 0.73 → 0.82 | 0.90 → 0.92 | 167 / 114 → 121 / 145 | 0.42 → 0.48 | 0.9719 → 0.9665 |
| pb | 0.9763 → 0.9727 | 96.9 → 96.7% | 2.07 → 2.16 | 2.87 → 2.65 | 54 / 27 → 52 / 32 | 0.37 → 0.40 | 0.9743 → 0.9730 |
| ont | 0.9730 → 0.9741 | 96.5 → 97.4% | 1.58 → 1.77 | 2.00 → 1.75 | 62 / 32 → 42 / 40 | 0.36 → 0.50 | 0.9751 → 0.9724 |

- **The paired-end model is the best so far**, on a harder test set: F1 +0.003 against v8 (+0.006 against v7's curve,
  +0.008 against v5's), misses down by a third, the log loss 0.0393 → 0.0313, and the two symptoms of v8 gone: the
  best test threshold is 0.50 (v8: 0.47) and the species held out and the test set agree (v8's gap of 0.0076 is
  −0.0004). The misses at 50,000 and 200,000 pairs, v8's cost, fell from 6.9% and 5.4% of the present taxa to 4.7%
  and 3.3%; the misses of one-fragment taxa from 18% to 12%, of strains with 1-10 fragments from 11.3% to 7.1%.
- **The single-end model traded misses for false positives** at about par (F1 −0.002 on the harder set): sensitivity
  +0.7, but false positives per absent taxon +12%, mostly "other" absent taxa (the closest species in the sample is in
  the database: 0.40 → 0.55%) at 50,000-200,000 pairs. The long-read models are within noise, with the same pattern
  (ONT misses −32%, false positives +25% per absent taxon).
- **The ablation says what did it**: the species priors, +0.007 to +0.009 of test F1 for pe and se on the same
  tables; the unfiltered group +0.002 (pe) and −0.001 (se) at knob 0.5 for a 4% lower log loss, although
  `failed_candidate_rate` is the most important feature of both models; the depth feature +0.001 to +0.002 on top of
  them. And the priors' gain is the cluster-size information alone: with CheckM quality and the ANI radius only, the
  group is worth nothing.
- **That gain rests on a rule the simulation cannot test.** Every present species of the simulation comes from a pool
  of species with several genomes (strains to simulate from), while 57% of the database's species, and of the absent
  taxa, are clusters of one genome. The model learns that a divergent read cloud on a one-genome species is a
  relative of a missing species, never a strain of that species: true in training (a singleton has no strain to
  simulate), true in environments GTDB has sampled densely, and untestable here for the environments where most
  species are singleton MAGs. The paired-end model with the priors calls 5 of 3,349 such clouds, the one without 32.

## v8 → v9 by depth, fragments and class

Paired-end test set, miss rate of the present taxa and false positives per 100 absent taxa
([`compare_v8_v9_output.txt`](compare_v8_v9_output.txt)):

| read pairs | FN rate v8 | FN rate v9 | FP / 100 absent v8 | FP / 100 absent v9 |
|---:|---:|---:|---:|---:|
| 500 | 5.7% | 4.2% | 7.7 | 22.2 (4 of 18) |
| 2,000 | 2.3% | 2.4% | 19.6 | 9.9 |
| 10,000 | 3.0% | 2.4% | 7.4 | 5.9 |
| 50,000 | **6.9%** | **4.7%** | 5.0 | 5.3 |
| 200,000 | **5.4%** | **3.3%** | 0.49 | 0.69 |
| 1,000,000 | 1.7% | 1.6% | 0.54 | 0.34 |
| 5,000,000 | 0.9% | 0.7% | 0.13 | 0.09 |

| present taxon's fragments | FN rate v8 | FN rate v9 | v9 strains |
|---|---:|---:|---:|
| 1 | 18.0% (108 of 599) | 11.8% (64 of 543) | 13.2% |
| 2 | 6.3% | 4.6% | 6.0% |
| 3-10 | 2.8% | 1.9% | 3.0% |
| 11-100 | 0.8% | 1.1% | 1.8% |
| > 100 | 0.6% | 0.4% | 0.6% |

False positives per 100 simulated missing species (the measure that does not depend on how many are missing): held
out 6.65 → 4.81, test set 4.53 → 3.12 (v7: 6.0). Of the test set's 135 paired-end false positives, 105 are beside
a species the database lacks and 2 come from held-out clades above the species rank. The 30M-pair point behaves
(3 samples: 186 of 188 present taxa found, 5 false positives among 12,051 absent). The single-end miss rates fall
alike (50,000 pairs 6.2 → 4.6%, 200,000 5.8 → 4.2%), the false positives per 100 absent taxa rise at 50,000 (5.6 →
6.4) and 200,000 (0.69 → 1.18).

## The ablation on v9's tables

Test F1 at knob 0.5, with the false positives split by whether the closest species in the sample is one the
database lacks, and the misses by fragments ([`ablation_v9_summary_output.txt`](ablation_v9_summary_output.txt)):

| features (pe) | F1 | best threshold | FP beside a missing species / other | FN (1 fragment) | log loss | AP |
|---|---:|---:|---|---|---:|---:|
| `normalized+adjacency+distance+depth+divergence` (v8's set) | 0.9614 | 0.58 | 147 / 35 | 135 (66) | 0.0379 | 0.9924 |
| + unfiltered | 0.9634 | 0.54 | 137 / 28 | 135 (65) | 0.0363 | 0.9926 |
| + priors | **0.9702** | 0.56 | 106 / 27 | 111 (62) | 0.0317 | 0.9946 |
| + unfiltered + priors (the default, v9's model) | 0.9699 | **0.50** | 107 / 28 | 111 (64) | **0.0313** | **0.9948** |
| the default without depth | 0.9689 | 0.48 | 107 / 28 | 119 (68) | 0.0321 | 0.9941 |
| the default, cluster size neutral and unknown ANI imputed | 0.9674 | 0.56 | 111 / 28 | 127 (67) | 0.0315 | 0.9945 |
| the default, priors = CheckM quality and ANI radius only | 0.9629 | 0.56 | 140 / 29 | 135 (70) | 0.0353 | 0.9933 |

| features (se) | F1 | best threshold | FP beside / other | FN (1 fragment) | log loss |
|---|---:|---:|---|---|---:|
| v8's set | 0.9606 | 0.56 | 151 / 28 | 137 (79) | 0.0460 |
| + unfiltered | 0.9598 | 0.56 | 146 / 29 | 147 (80) | 0.0442 |
| + priors | 0.9664 | 0.56 | 122 / 26 | 121 (71) | 0.0387 |
| the default | **0.9668** | 0.48 | 119 / 26 | 121 (71) | **0.0377** |
| the default without depth | 0.9648 | 0.50 | 125 / 32 | 125 (76) | 0.0397 |
| cluster size neutral, ANI imputed | 0.9641 | 0.52 | 120 / 29 | 138 (78) | 0.0397 |
| priors = CheckM and radius only | 0.9594 | 0.60 | 144 / 29 | 152 (83) | 0.0430 |

- **The design change alone (v8's set on v9's tables)** made the model more permissive than v8: best threshold 0.58,
  false positives per absent taxon 0.85% (v8 on its own test set 0.67%), misses 3.3% of the present (v8 3.8%). The
  mixed σ gave the model the test design's tails; the groups added on top pull the best threshold back to 0.50.
- **The unfiltered group ranks and does not call**, as on the benchmark world: log loss −4%, false positives beside
  missing species −10 (pe), F1 within ±0.002. `failed_candidate_rate` is the top feature (importance 0.125 pe, 0.162
  se; 10th for ONT, 18th for PacBio) because it separates cleanly: median 0.94 for absent taxa, 0.13 for present ones,
  0 for present taxa with one or two fragments; taxa above 0.8 are present in 0.27% of cases (16,454 rows), below
  0.0001 in 63%. But what it separates, the depth and distance features already mostly did at the knob.
  `fragments_all` and `em_fragments` have importances of 0.001: `fragments_all` is `fragments` plus ~15% (median 23
  against 20 for present strains), the MAPQ filter drops few whole reads of these simulated strains.
- **The priors are the cluster size.** `cluster_genomes_log10` and the unknown (−1) intra-species ANI of one-genome
  clusters carry the group (ranks 24-27 by importance, 0.004 each); `cluster_ani_radius` is 95 for nearly every
  species (importance 0.0005), `rep_completeness` and `rep_contamination` 0.001, `rep_duplicate_share` 0 for every
  species (GTDB's marker files hold one copy per genome; the docs now say so). Neutralising the size and imputing
  the unknown ANI with the median leaves +0.004 (pe), part of which is the imputed spike the forest can still find;
  CheckM and radius alone leave nothing.

### Where the cluster-size prior acts

Miss rate of present taxa and false-positive rate of absent taxa by the GTDB cluster's size and the taxon's fragments,
paired-end test set ([`ablation_by_cluster_output.txt`](ablation_by_cluster_output.txt)):

| features | FN, 1-genome species, ≤ 2 fragments | FN, ≥ 5 genomes, > 2 fragments | FP, 1-genome species, > 2 fragments | FP, 1-genome, ≤ 2 fragments | FP, ≥ 5 genomes, > 2 fragments |
|---|---|---|---|---|---|
| v8's set | 4.8% (6/124) | 2.6% (27/1,039) | 1.0% (32/3,349) | 0.5% (42/8,840) | 2.0% (23/1,164) |
| + unfiltered | 5.6% (7) | 2.6% (27) | 0.9% (31) | 0.4% (34) | 2.1% (25) |
| + priors | 9.7% (12) | 1.6% (17) | **0.1% (5)** | 0.2% (20) | 2.1% (24) |
| the default | 8.1% (10) | 1.3% (14) | 0.2% (6) | 0.2% (18) | 2.1% (24) |
| CheckM and radius only | 5.6% (7) | 2.4% (25) | 0.9% (31) | 0.4% (37) | 1.9% (22) |

Composition ([`priors_bias_output.txt`](priors_bias_output.txt)): of the test set's present taxa 17% are one-genome
species, 50% clusters of 2-4 genomes, 33% of 5 or more; of the absent taxa 57%, 27%, 16%; of the 37,978 species that
appear as taxa in the training dump 56% are one-genome clusters, as in GTDB. The present one-genome species are
simulated from their representative (there is no other genome), so their reads align at identity 1 and the model
finds them at more than two fragments without fail, prior or not. What the prior changes is the two populations
above: a divergent read cloud on a one-genome species is called a relative of a missing species (false positives
1.0% → 0.1%, 27 fewer on the test set, about a third of the priors' gain), and a present multi-genome species with
few reads is believed more readily (misses 2.6% → 1.6%, the other two thirds), at the price of more doubt about
present one-genome species with one or two reads (4.8% → 9.7%, 6 → 12 of 124).

In the simulation the rule is right by construction. In a sample from an environment GTDB has sampled densely it is
right too: the species present are those found before. In soil, sediment or sea water, where most GTDB species are
single MAGs and a present population is a strain of one of them at 96-99% identity, the rule rejects exactly those
calls, and the simulation has no way to show it (a one-genome species cannot be simulated as a strain). The model
without the priors, trained on the same tables, is in `local/v9/models_without_priors/` (`model_pe.xml`,
`model_se.xml`, with their reports): protal's `--model` and `--model_se` take them, so the two models can be
compared on real samples of both kinds, and the calls that differ are by construction the divergent clouds on
one-genome species.

## What worked, what did not

Worked:

1. **The mixed training σ with 30% held out**: calibration at 0.50, the held-out/test gap closed, the mid-depth
   misses halved, the 30M-pair point unproblematic; run time unchanged. The design is right.
2. **The cluster-size prior**: +0.007 to +0.009 of F1, the largest gain of the round, with the caveat above.
3. **The depth feature**, still +0.001 to +0.002 with everything else in, and the divergence features (part of the
   base here).
4. **`failed_candidate_rate`** as a ranking feature: the most important one, −4% log loss, the cleanest separation of
   present from absent singletons protal has; at the knob it adds little because the knob was already placed by the
   depth and distance features. It would matter more under `--fdr` or for the probabilities themselves.
5. **The infrastructure**: the species clusters file fetched and converted for all species, parity for all four
   models, the `ZF` records costing no time the logs show, the suspect-copy scan identical.

Did not:

6. **`fragments_all` and `em_fragments`**: nothing (importance 0.001). The misses the MAPQ filter was blamed for are
   not whole reads lost but ties; the EM's share (`em_own_share`) was already a feature. Harmless to keep.
7. **CheckM quality, ANI radius, duplicated markers**: nothing; the last cannot be non-zero for a GTDB release.
8. **The single-end and long-read false positives** rose by 12-25% per absent taxon while the misses fell more;
   F1 at par. The PacBio model's best test threshold is below 0.5 for the third run (0.40; held out 0.51 this time).

## Recommendations

1. **Ship v9's paired-end and single-end models as they are, and test the cluster-size prior on real data before
   relying on it for environmental samples**: a few gut and a few soil or marine samples, profiled with the
   database's model and with `--model local/v9/models_without_priors/model_pe.xml`; the calls that differ are the
   divergent clouds on one-genome species, and whether they look like strains (identity 96-99% across many genes,
   reads on both conserved and fast genes) or like relatives (conserved genes only) decides. If the prior is wrong
   there, the fix in the trainer is one line: `--features normalized+adjacency+distance+depth+divergence+unfiltered`,
   retrained from the uploaded tables; no rebuild.
2. **Make the prior learnable instead of a rule**: simulate strains of one-genome species in silico, the
   representative mutated at a divergence drawn from GTDB's intra-species ANI distribution with the codon-position
   bias of real strains, so that the training data hold present divergent clouds on singletons and the forest weighs
   the prior against the evidence. This is the F1 list's new item; it also enlarges the strain population the misses
   are made of.
3. **Leave the defaults otherwise.** The unfiltered group costs nothing and improves the probabilities; the two inert
   features of each group are constant or near-constant and do not bind the forest. A flat knob for the PacBio model
   is still not supported by the species held out.
4. **For the error budget**: the paired-end test set now has 111 misses (64 at one fragment, 87 strains with 1-10
   fragments) and 135 false positives (105 beside a missing species). Both populations are now at the level where the
   next gains need the index (strain alleles) or the simulation (item 2), not features.

## Follow-up (2026-10-04)

Asked which defaults to change before the next build, the choice fell on the prior: the `priors` group left the
default feature set (`DEFAULT_FEATURE_SET` is `normalized+adjacency+distance+depth+divergence+unfiltered`) and stays
available with `--features ...+priors` for a database meant for environments GTDB has sampled densely, where the
cluster-size rule is a fair bet. The reasoning: the gain is a base rate the simulation sets by accident (the species
pool has strains, so its present species have several genomes) and cannot contradict (a one-genome species has no
strain to simulate), so the model holds it as a rule rather than weighing it; a general-purpose database should be
neutral about one-genome species until the simulation can present them as strains (recommendation 2). On v9's
tables the default without the priors scores 0.9634 (pe) and 0.9598 (se) on the test set against 0.9699 and 0.9668
with them; the next build's summary will read lower than v9's for that reason and no other. The database still
carries `species_priors.tsv`, and the dumps the columns, so the opt-in model can be trained from any build's tables.
