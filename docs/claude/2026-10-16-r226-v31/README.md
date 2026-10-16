# The r226 v31 build (protal 0.8.5): level with v29 where the samples are the same, ONT better, and a blind spot for deep, species-rich, uneven samples

Date: 2026-10-16. Branch `audit-fixes` at `956e18b`, the commit the build ran. Nothing in protal or the scripts was
changed for this report.

**Data.**
- **v31:** SLURM 24292401 (a 1.5 TB node), 2026-10-16 14:42 to 19:10, 120 threads, `--scratch` on the node's SSD.
  protal, the simulator and the scripts at `956e18b` (0.8.5 with the F1 audit's fixes, the ONT identity floor 0.8,
  `--gene_end_clips credit`, `--sparse_first_min_picks 6`, spike mode `auto`, the hashed key map; see
  `protal_db/build_metadata.tsv`). The command: `build_gtdb_database.py --inputs .../refGenos/ --outdir ... --share-logs
  --protal ... --simulator ... -t 120 --scratch ...`. Everything else is a default of `956e18b`:
  - database gates `pe:4,se:4,pb:16,ont:8`, picked by stride; models also at `pe:1/8,se:1/8,pb:8/32,ont:1/16`;
  - GBM models as the database's default, forests beside them (`--model gbm,forest`);
  - the trainer's crossed folds (`--cv-scheme samples+species`) choose the knob; the test set no longer vetoes it;
  - paired-end training at 75, 100, 150 (HiSeq X and NovaSeq) and 250 bp; single-end at 75, 100 and 150 bp and Ultima;
  - a test seed drawn for the build (283979), so the test samples are new draws;
  - new steps: the final check on the finished database, and the test samples without the allele scores.
  - **No `--test-c-gates`:** the test samples were profiled only at the database's gate and the models' gates. There
    are no grid, hash, ungated, C 2 or `--flex_indels` rows this time.
- **The archive:** `local/protal0.8.5_r226_v31_share.tar.gz` (1.1 GB). Unpacked to `local/v31/` for this report and
  removed afterwards (not in git).
- **Compared with v29** ([2026-10-15-r226-v29](../2026-10-15-r226-v29/README.md); `local/protal0.8.4_r226_v29_share.tar.gz`,
  its logs and `model_logs` unpacked to `local/v29/` and removed afterwards). Same download, training seed (1),
  hold-out list (16,120 species) and genome table.
  - **The same simulations in both builds:** the 50 scenario hold-in samples per read type (30 for se), 154 of the
    paired-end and 76 of the single-end design training samples, but only 9 PacBio and 12 ONT design samples.
  - **Not the same:** every test sample (new test seed), and most design samples of se, pb and ont (new read setups).
- **Commands** (WSL, `taskset -c 0-1 nice -n 5`, Python of the `protal-db-build` environment):
  - unpack: `tar -xzf local/protal0.8.5_r226_v31_share.tar.gz -C local/v31`, and for v29
    `tar -xzf local/protal0.8.4_r226_v29_share.tar.gz -C local/v29 --wildcards '*/console.log' '*/logs/*'
    '*/model_logs/*.txt' '*/model_logs/*.tsv' '*/model_logs/*.json' '*/model_logs/*calls.tsv.gz' '*/protal_db/*'`;
  - [pairs.py](pairs.py) `local/v31/protal0.8.5_r226_v31 local/v29/protal0.8.4_r226_v29` → [pairs.txt](pairs.txt):
    1. v31's GBM against v29's GBM on the training samples simulated alike, with a paired bootstrap over samples;
    2. every moderate, soil and shallow-soil sample: present taxa, misses, and the misses among species with at most
       2 marker reads;
    3. the share of thin training rows (1-2 fragments) that are present, by the sample's depth and number of taxa.
  - awk and grep over `console.log`, the reports, `final_db_check.tsv`, `logs/protal_cpu_*.tsv`,
    `logs/protal_runs_*.log`, `composition_accuracy_*.tsv` and `c_gates.tsv` (sections 4-7);
  - section 8: [retrain_without_taxa.sh](retrain_without_taxa.sh) `REPO local/v31/protal0.8.5_r226_v31 ~/v31_no_taxa`
    (`taskset -c 0-3 nice -n 5`, 28 minutes), then [compare_no_taxa.py](compare_no_taxa.py)
    `~/v31_no_taxa/out local/v31/protal0.8.5_r226_v31/model_logs` → [no_taxa.txt](no_taxa.txt).

## Summary

- **The build ran clean.** All 13 steps finished in 4:27:50. Parity of all 24 models with protal holds on the training
  database. The only warnings are the final check's (section 3).
- **Where v29 and v31 profiled the same samples, the models are level, and ONT is better** (section 1; species held out,
  each model at its own knob, a present species without a row counted as missed):
  - pe: design +0.0012 [−0.0001, +0.0028], scenarios −0.0010 [−0.0017, −0.0003];
  - se: +0.0005 and +0.0006, pb: +0.0058 (9 samples) and +0.0004; all within noise except pb's small set;
  - **ONT scenarios +0.0056 [+0.0041, +0.0071].** Present species without a row fell from 849 to 40 of 46,887: the
    identity floor of 0.8 (`bb2ac82`) did what its report predicted. ONT's rowless share with species held out is
    now 2.4% (v29 3.7%), on the test set 2.8% (v29 5.2%).
- **Obvious weakness: deep, species-rich, uneven samples lose most of their thin species** (section 2). The soil
  hold-out F1 fell from 0.954 (v29) to 0.914 (pe), 0.953 to 0.916 (se), 0.970 to 0.960 (pb) and 0.960 to 0.946 (ont).
  This is not a regression in the models. The test seed drew new samples, and one of them is of a kind the training
  samples never contain:
  - Soil sample 4 has the most uneven abundances (lognormal sigma 2.5), the most species (10,521) and is deep
    (1.38× the scenario's depth). pe misses 894 of its 3,483 present species, 697 of them among the 823 with 1-2
    marker reads (85%). Eight of the ten training soil samples miss 22-40% of such species (the other two have 9 and
    16 such species and miss nearly all). Without that sample the pe soil hold-out F1 is 0.952, as in v29.
  - The rows' own evidence looks like that of the thin species the training samples call (identity, allele identity,
    fragments, held-out congeners). Their median score is 0.11 against 0.89. What differs is the sample: its
    `sample_log_fragments` and `sample_log_taxa` are both at the top of the training range.
  - In the training rows, a thin row (1-2 fragments) in a sample that deep and that rich is present 34 times in
    25,676 (0.1%). The deep, rich training samples all had even abundances (sigma 1.0), so their thin rows were
    the spill of held-out congeners. The model learned "thin row in a deep, rich sample = absent" and applies it to
    real tail species.
  - The same holds for se (soil sample 4: 559 of 691 thin species missed), shallow soil sample 2 (pe 460 of 598) and,
    weaker, ONT (soil sample 4: 229 of 1,029). PacBio is barely affected (139 of 920).
  - Real deep soil and sediment metagenomes look like this sample: thousands of species and a long tail. This is the
    most important finding of the build.
- **The final check works as designed and says the finished database is better, but its gain is partly a leak**
  (section 3). On the test samples' rows of species not held out, F1 rises by 0.015-0.033 from the training database
  to the finished one. Most of that is expected: the held-out species are in the finished database, so their reads
  no longer make false positives and confused rows (pe FP 665 → 302). In the samples without a held-out congener
  the gain is +0.004..+0.008, from false negatives halving (pe 556 → 289). The finished database takes its strain
  alleles from every genome (`--final-allele-genome-share 1`), including the real strains the test samples were
  simulated from. pe's model gives `allele_identity` 72% of its split gain, so part of that gain is the
  simulated strains finding themselves. The archive cannot split it.
- **18 pe rows (3 se, 1 pb, 1 ont) where protal's calls on the finished database differ from its model file's.** Of
  48,221 pe rows. The parity on the training database passed (equal probabilities). The rows are not in the
  archive; near-knob rounding is the likeliest cause, unverified.
- **The allele scores do nothing at r226** (the F1 audit's B1, measured now): pe −0.0008 with them (with strain alleles
  −0.0013), se +0.0001, forests ±0.0003. FN rise and FP fall slightly with them.
- **The C 4 grid with a sparse index was not used at all** (your question; section 5). The build had no
  `--test-c-gates`, so no grid rows. Every protal run loaded the whole index (23.6 GB and a 1.1 GB hashed key map;
  42-47 GB peak) and every gated sample was gated by stride (2,891 in the training runs, 894 in the test runs, at
  1/4, 1/8, 1/16 or 1/32). The last grid numbers are v29's emulated rows.
- **Other observations** (sections 4-7):
  - The extra gate models (`--model-gates`, kept despite v29's recommendation) are 70% of both collections' alignment
    CPU. Their held-out F1 is again within 0.001 of the database gate's model.
  - On the same scenario samples, alignment CPU is +20-23% for ONT (the lower floor; predicted +17%), +5-11% for pe
    (shallow soil +47%), level for PacBio and se.
  - Spike mode fired once, falsely, on a simulated sample (`sc_gut_pb_b6000000000_s_7`: 14 "spike genes", 89,618
    records left out; there are no PCR piles in simulations). That sample's F1 (0.982) is within its scenario's range.
  - PacBio's species depth bias is unchanged on the same samples (median depth ratio 1.159 in both builds). The test
    set's +35.3% (v29 +27.6%) is a different sample mix. Still open since v19.
  - The new single-end 75 and 100 bp test samples are the weakest se setups (F1 0.951 and 0.963); se's test F1 fell
    from 0.967 to 0.961 because they were added, not because the old setups got worse.

## 1. Headline numbers, and v31 against v29 on the same samples

F1 from the consoles, GBM models (v29's default was the forest; its GBM is the comparison). Species held out at knob
0.5; the test set and the hold-out scenario samples as protal calls. **Only the species-held-out column profiles the
same scenario samples in both builds; the test samples are new draws in v31.**

| read type | build | species held out (rows) | marker truth, held out | test set (rows) | marker truth, test | soil hold-out | soil_shallow hold-out |
|---|---|---:|---:|---:|---:|---:|---:|
| pe | v29 GBM | 0.9491 | 0.9519 | 0.9647 | 0.9546 | 0.9540 | 0.9311 |
| | v31 | 0.9507 | 0.9517 | 0.9671 | 0.9583 | **0.9143** | **0.9166** |
| se | v29 GBM | 0.9573 | 0.9525 | 0.9668 | 0.9558 | 0.9531 | — |
| | v31 | 0.9566 | 0.9522 | 0.9612 | 0.9496 | **0.9159** | — |
| pb | v29 GBM | 0.9718 | 0.9702 | 0.9847 | 0.9814 | 0.9703 | 0.9618 |
| | v31 | 0.9718 | 0.9704 | 0.9836 | 0.9779 | 0.9599 | 0.9551 |
| ont | v29 GBM | 0.9656 | 0.9525 | 0.9751 | 0.9542 | 0.9604 | 0.9517 |
| | v31 | 0.9628 | **0.9579** | 0.9747 | **0.9657** | 0.9461 | 0.9428 |

- ONT's row-based F1 fell while its marker-truth F1 rose: the lower floor gives more species a row (absent ones too:
  88,222 absent rows with species held out against 60,970, on partly different design samples). Only the marker
  truth compares floors
  ([ONT rowless species](../2026-10-15-ont-rowless-species/README.md)).
- Knobs: pe 0.69 (v29 0.74; the species-only scheme would again have chosen 0.74), se 0.59 (v29's was reset to 0.5 by
  the test veto, now gone), pb and ont 0.5. se's test F1 at its knob is 0.9612 against 0.9616 at 0.5.

**The same samples** ([pairs.txt](pairs.txt), section 1). Training rows scored with species held out, each model at
its own knob, a present taxon without a row in a build counted as missed there:

| read type | set | samples | present | rowless v29/v31 | F1 v29 | F1 v31 | v31 − v29 | 95% interval |
|---|---|---:|---:|---|---:|---:|---:|---|
| pe | design | 154 | 7,680 | 4/2 | 0.9676 | 0.9688 | +0.0012 | [−0.0001, +0.0028] |
| pe | scenarios (hold-in) | 50 | 61,217 | 18/53 | 0.9512 | 0.9502 | −0.0010 | [−0.0017, −0.0003] |
| se | design | 76 | 3,978 | 1/2 | 0.9671 | 0.9676 | +0.0005 | [−0.0023, +0.0033] |
| se | scenarios (hold-in) | 30 | 38,272 | 14/11 | 0.9539 | 0.9545 | +0.0006 | [−0.0002, +0.0015] |
| pb | design | 9 | 611 | 6/1 | 0.9844 | 0.9902 | +0.0058 | [+0.0013, +0.0093] |
| pb | scenarios (hold-in) | 50 | 45,335 | 1/2 | 0.9699 | 0.9703 | +0.0004 | [−0.0001, +0.0009] |
| ont | design | 12 | 735 | 4/5 | 0.9811 | 0.9837 | +0.0027 | [−0.0021, +0.0079] |
| ont | scenarios (hold-in) | 50 | 46,887 | **849/40** | 0.9547 | 0.9604 | **+0.0056** | [+0.0041, +0.0071] |

- pe's −0.0010 on the scenarios is the knob: at 0.69 instead of 0.74 v31 calls 374 more false positives and misses
  220 fewer. pe also lost the row of 35 more present species (53 against 18 of 61,217); `--gene_end_clips credit` and
  the gate's least picks are the 0.8.5 changes that touch rows. Small, not traced here.
- The design rows of se, pb and ont changed (new read setups), so their larger sets cannot be paired.

## 2. Deep, species-rich, uneven samples: the thin species are missed

The scenario samples' evenness comes from their number: lognormal sigma 1.0, 1.5, 2.0, 2.5 for samples 1, 2, 3, 4,
then again (`scenarios.py` `EVENNESS`). Species count and depth are drawn apart from it. The training draw (seed 1,
the same in v29 and v31) gave the soil scenario's deep, rich samples sigma 1.0 (samples 5 and 9: 3,153 and 3,464
present, 16 and 9 present species with 1-2 marker reads). v31's test draw gave soil sample 4 sigma 2.5, 10,521
species and 1.38× depth.

The soil and shallow-soil samples with the most misses ([pairs.txt](pairs.txt), section 2; "thin": present
species with at most 2 marker reads; p thin: their median score):

| read type | sample | set | sigma | sample_log_fragments | sample_log_taxa | present | FN | thin | thin FN | p thin |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| pe | soil 4 | hold-out | 2.5 | 5.925 | 4.326 | 3,483 | 894 | 823 | 697 (85%) | 0.11 |
| pe | soil_shallow 2 | hold-out | 1.5 | 5.315 | 4.211 | 3,746 | 662 | 598 | 460 (77%) | 0.16 |
| pe | soil 8 | hold-in | 2.5 | 5.898 | 4.189 | 2,128 | 222 | 427 | 159 (37%) | 0.89 |
| pe | soil 10 | hold-in | 1.5 | 5.244 | 4.208 | 3,967 | 316 | 727 | 189 (26%) | 0.95 |
| pe | soil 9 | hold-in | 1.0 | 5.902 | 4.487 | 3,464 | 120 | 9 | 9 (100%) | 0.07 |
| pe | soil 5 | hold-in | 1.0 | 5.695 | 4.375 | 3,153 | 111 | 16 | 15 (94%) | 0.12 |
| se | soil 4 | hold-out | 2.5 | 5.853 | 4.130 | 3,392 | 754 | 691 | 559 (81%) | 0.12 |
| ont | soil 4 | hold-out | 2.5 | 5.367 | 3.934 | 2,463 | 272 | 1,029 | 229 (22%) | 0.96 |
| pb | soil 4 | hold-out | 2.5 | 5.324 | 3.743 | 2,288 | 177 | 920 | 139 (15%) | 0.99 |

- The hold-in samples 9 and 5 already show it: their few thin present species are missed almost all. Those samples
  have so few that it cost nothing until a test sample had hundreds.
- The thin species of the two affected test samples look like those of the training samples by their own evidence.
  Medians, pe, soil 4 and shallow soil 2 against all training soil and shallow-soil thin species: identity 0.992 /
  0.993, `allele_identity` 0.993 / 0.993, fragments 2 / 1. A held-out congener is in the sample for 71% of soil 4's
  thin species against 62% of training soil 8's. Ranked over all 98 features, the farthest from the training
  samples are the sample-level ones (`sample_log_taxa` +1.1 IQR, `sample_log_fragments` +1.2 IQR), then
  `failed_candidate_rate` and the `relative_*` features of a crowded sample.
- **What the training rows teach** ([pairs.txt](pairs.txt), section 3; rows with 1-2 fragments of all training samples,
  every read type):

  | sample_log_fragments | sample_log_taxa | training thin rows | present | hold-out thin present rows | called |
  |---|---|---:|---:|---:|---:|
  | 4.0-5.0 | 3.6-4.0 | 50,995 | 18.5% | 4,195 | 85% |
  | 5.0-5.5 | 3.6-4.0 | 58,326 | 7.0% | 1,896 | 81% |
  | 5.0-5.5 | > 4.2 | 9,359 | 8.0% | 584 | 25% |
  | 5.8 and more | 4.0-4.2 | 15,118 | 3.8% | 665 | 24% |
  | 5.8 and more | > 4.2 | 25,676 | **0.1%** | 705 | **19%** |

  A thin row in a deep, rich training sample is almost always the spill of a congener the database lacks (soil leaves
  60% of its species out of the training database). The model has learnt that prior from 10 soil and 10 shallow-soil
  samples, and it holds only as long as deep, rich samples have even abundances.
- Without soil sample 4, pe's soil hold-out F1 is 0.952 (v29 0.954); without shallow-soil sample 2, pe's shallow
  soil is 0.931 (v29 0.931). The trainer flagged the soil scenario ("hold-out samples score F1 0.914 against 0.945 on
  its training samples ... too few training samples of it, or they differ").
- v29 had the same prior: the hold-in samples, identical in both builds, score their thin species alike (pe soil 8:
  159 of 427 missed in v31, 166 of 428 in v29). v29's test draw simply had no such sample.

## 3. The final check on the finished database

`model_logs/final_db_check.tsv`: the test samples written to disk (the streamed, deepest ones left out) profiled again
on `protal_db`, rows of species not held out, the training database's calls → the finished database's.

| read type | rows of | F1 training DB → finished DB | FP | FN |
|---|---|---|---|---|
| pe (GBM) | species not held out | 0.9488 → 0.9782 (+0.0294) | 665 → 302 | 1,622 → 675 |
| | of them without a held-out congener in the sample | 0.9752 → 0.9831 (+0.0078) | 159 → 204 | 556 → 289 |
| se | species not held out | 0.9615 → 0.9766 (+0.0151) | 340 → 150 | 680 → 467 |
| | without a held-out congener | 0.9747 → 0.9796 (+0.0049) | 121 → 122 | 397 → 297 |
| pb | species not held out | 0.9696 → 0.9946 (+0.0250) | 190 → 22 | 456 → 92 |
| | without a held-out congener | 0.9888 → 0.9968 (+0.0081) | 12 → 8 | 126 → 31 |
| ont | species not held out | 0.9603 → 0.9902 (+0.0299) | 356 → 61 | 535 → 155 |
| | without a held-out congener | 0.9887 → 0.9938 (+0.0052) | 45 → 33 | 100 → 46 |

- **The first row of each read type is the expected gain.** The held-out species are in the finished database. Their
  reads go to them and no longer make false positives of their congeners, nor confused rows of present congeners.
  The finished database has half the rows (pe 48,221 against 95,951).
- **The second row is partly a leak.** In samples without a held-out congener, false negatives halve for every read
  type. The finished database takes its strain alleles from every genome of a species (`--final-allele-genome-share 1`,
  `build_gtdb_database.py` lines 2564-2566); the training database leaves out the genomes the simulations draw strains
  from. 30.4% of the simulated species are real strains, whose own genomes are now alleles. pe's GBM puts 72% of its
  split gain on `allele_identity`. Real samples' strains are mostly not in GTDB, so the finished database's gain on
  real data lies between the two rows. The check does not split rows by genome source, and the share archive does not
  hold its rows, so the size of the leak is unknown.
- pe's false positives rise in those samples (159 → 204; 80 absent taxa newly called). The other read types' fall.
- **Calls unlike the model file's: pe 18, se 3, pb 1, ont 1 rows** (`database_check.py` `table_sides`: protal's
  `prediction` column against the PMML's calls on the same features and knob). Spike mode, which changes protal's
  evidence, entered no sample of the test collection's runs; the final check's own protal log is not in the archive.
  The parity of step 9 found equal probabilities for every model on 4,977 pe taxa.
  The likeliest cause is a feature printed with fewer digits than protal used, at rows scored within a hair of the
  knob. Unverified: the rows are not in the archive.

## 4. The allele scores (`--no_allele_scores`, the F1 audit's B1)

Short-read test samples written to disk, profiled again on the training database without the allele scores; the
same models (trained with them) score both.

| model | taxa | F1 without → with the allele scores | FP | FN |
|---|---|---|---|---|
| pe GBM, all | 52,834 | 0.9496 → 0.9488 (−0.0008) | 687 → 665 | 1,566 → 1,622 |
| pe GBM, species with strain alleles | 14,749 | 0.9511 → 0.9499 (−0.0013) | 365 → 350 | 867 → 913 |
| pe GBM, without | 38,085 | 0.9475 → 0.9473 (−0.0002) | 322 → 315 | 699 → 709 |
| se GBM, all | 32,474 | 0.9614 → 0.9615 (+0.0001) | 363 → 340 | 660 → 680 |
| pe forest, all | 52,834 | 0.9485 → 0.9483 (−0.0002) | 934 → 923 | 1,390 → 1,408 |
| se forest, all | 32,474 | 0.9565 → 0.9565 (0.0000) | 470 → 465 | 688 → 692 |

The allele scores trade a few false positives for slightly more misses, and the misses fall on species with strain
alleles, as B1 predicted (an allele-rich species takes reads from its neighbours). Net: no gain at r226. A model
trained without them might differ; this run could not test that.

## 5. The gates, and the grid's sparse index

`c_gates.tsv` (GBM; test samples' rows at each model gate, unseen species counted as missed):

| read type | gate | samples | unseen | F1 | Δ to the database's gate | aligning CPU s |
|---|---|---:|---:|---:|---:|---:|
| pe | 1 | 150 | 593 | 0.9333 | −0.0001 | 5,231 |
| | **4 stride** | | 593 | **0.9334** | | **3,687** |
| | 8 stride | | 593 | 0.9330 | −0.0004 | 3,323 |
| se | 1 | 142 | 449 | 0.9337 | −0.0001 | 3,131 |
| | **4 stride** | | 455 | **0.9338** | | **2,107** |
| | 8 stride | | 463 | 0.9350 | +0.0012 | 1,789 |
| pb | 8 stride | 61 | 272 | 0.9614 | +0.0002 | 5,363 |
| | **16 stride** | | 274 | **0.9611** | | **4,679** |
| | 32 stride | | 280 | 0.9614 | +0.0002 | 3,979 |
| ont | 1 | 62 | 459 | 0.9456 | −0.0007 | 9,713 |
| | **8 stride** | | 478 | **0.9463** | | **5,868** |
| | 16 stride | | 528 | 0.9454 | −0.0009 | 5,117 |

- The database's gates still hold. se C 8 is free again (+0.0012, −15% CPU), as in v29. ONT C 16 now loses 0.0009 and
  50 rowless species: with the lower floor, ONT keeps more reads that a coarser gate drops.
- **The grid's sparse index was not used.** Without `--test-c-gates` the build profiles no grid rows (and no hash,
  ungated, C 2 or `--flex_indels` rows). Every protal run of both collections loaded the whole index ("Index in memory
  (hashed key map ...): 23.6 GB; key map 1.1 GB"; peak 42.5-46.9 GB), and every gated sample's log line reads "sparse
  first (1/C of the cores by stride": 853, 1,327, 474 and 237 samples at 1/4, 1/8, 1/16 and 1/32 in the training runs, 292,
  416, 124 and 62 in the test runs. The database is stride-picked (`build_metadata.tsv`: `picking by
  pe:stride,se:stride,pb:stride,ont:stride`), so a grid run would need `--sparse_first_pick grid` per run anyway.
  The last grid numbers are v29's, emulated on the whole index: pe C 4 −0.0004 against the stride, pb C 16 −0.0020,
  ont C 16 −0.0040, pe C 8 −0.0046. A real sparse index's memory and speed at r226 have not been measured in a build.
- The models at other gates: held-out F1 within 0.001 of the database gate's model for every read type and kind,
  as in v29. Their rows took 70% of the alignment CPU (training 116,299 of 166,248 CPU s; test 37,647 of 53,988).

## 6. Speed and build time

Alignment CPU seconds summed over the scenario samples that both builds simulated alike (`protal_cpu_training.tsv`,
database gate only, 10 samples per row):

| read type | gut | moderate | soil | soil_shallow | host |
|---|---:|---:|---:|---:|---:|
| pe | +5% | +8% | +11% | +47% | −5% |
| se | — | −1% | +4% | — | −31% |
| pb | −1% | −2% | −1% | −2% | −32% |
| ont | +20% | +22% | +23% | +23% | +3% |

- ONT's +20-23% is the identity floor of 0.8 (the floor's report predicted +17% at 150 Mb).
- pe's +5-11% fits `--sparse_first_min_picks 6` (the B2 fix: up to +21% gate lookups of a 150 bp read). The shallow
  soil's +47% is 176 CPU s over 10 samples; the host rows' −31..−32% are small sums too. These are sums on a busy
  120-thread node: ±10% is within noise.
- The build took 4:27:50 (v29 4:16, with `--test-c-gates`): collections 2:18:44, training 9.6 minutes for both kinds,
  the final check 20:45, the allele-score run 5:25, the error-read and ancestry reports beside them.

## 7. Other observations

- **Spike mode fired on a simulation.** `sc_gut_pb_b6000000000_s_7` (training, and its C 8 and C 32 rows): "2 of its
  10 taxa with the most reads and 7 in all have a gene whose reads pile far beyond their other genes, 22.5% of its
  reads", 14 spike genes, 89,618 records left out. Simulations have no PCR piles, so this is a false trigger: a genome
  with a real high-copy gene, or long reads piling on one gene of a dominant species. That sample's F1 (0.982) is
  within the gut samples' 0.981-0.998. The `*.spike_genes.tsv` that names the genes is on the node, not in the archive.
- **Composition.** pe, se and ont: species' depth −1.8%, −2.6%, +1.2% (ONT was −6.1% in v29: the floor keeps more of
  its reads). PacBio: +35.3% on the test set (v29 +27.6%), explained share +0.257, unknown share −0.198. On the
  scenario samples both builds simulated alike, PacBio's median depth ratio is 1.159 in both, so the larger test figure
  is the new sample mix. The bias itself is still open (since v19).
- **Feature importance** is as concentrated as in v29: pe `allele_identity` 0.72, se `lu_per_kb` 0.66, ont
  `excess_scaled_median` 0.70; pb moved from `excess_scaled_median` (0.35) to `allele_top_identity_bound` (0.46),
  likely the F1 audit's quality fix (C3).
- **Single-end by read length** (test design samples, F1 as protal calls): 75 bp 0.9513, 100 bp 0.9626, 150 bp
  0.9627, 250 bp 0.9643, Ultima 0.9632. Paired-end: 75 bp 0.9635, 100 bp 0.9680, 150 bp 0.9695, 150 bp NovaSeq
  0.9645, 250 bp 0.9697.

## 8. Quick test: the pe and se models without `sample_log_taxa`

Asked after the report's first version. Does the sample's number of taxa carry the prior of section 2?

**How.** [retrain_without_taxa.sh](retrain_without_taxa.sh) copies the build's pe and se tables (`work/training`,
`work/test`) twice: as built, and with `sample_log_taxa` set to 0 in every row. A constant column can never be split
on, so no tree uses it. The trainer has no option to drop one feature of a group: `sample_log_taxa` is one of the
three "complexity" features. It then trains a GBM on each copy with the build's trainer options (report line 3),
except `--threads 4` and without the C 1 / C 8 extra test tables, which the archive lacks. It ran in WSL on 4 cores:
8:48 per pe model, 4:50 per se model. [compare_no_taxa.py](compare_no_taxa.py) compares the two retrained models on
the same rows → [no_taxa.txt](no_taxa.txt).

**The retraining is close to the build, not identical.** The as-built retraining differs from the build's calls at
1,939 of 566,352 pe rows (0.34%) and 1,755 se rows. se's knob came out 0.5 instead of 0.59: the build's bootstrap kept
0.59 at 97% support, just over the 95% rule. The likeliest cause is the thread count (4 against 10), which changes
the order of the boosting histograms' sums. So the comparison is the two local models, as built against without
the feature. A refit moves pe soil sample 4's misses by about 20 (873 local against 894 in the build).

**Results** (each model at its own knob, as protal calls; 95% interval over samples):

| read type | set | F1 as built | F1 without `sample_log_taxa` | difference | 95% interval |
|---|---|---:|---:|---:|---|
| pe | held out, design | 0.9706 | 0.9700 | −0.0006 | [−0.0014, +0.0002] |
| | held out, scenarios | 0.9506 | 0.9509 | +0.0003 | [−0.0002, +0.0007] |
| | test, design | 0.9674 | 0.9697 | +0.0023 | [+0.0007, +0.0039] |
| | soil hold-out | 0.9161 | 0.9197 | +0.0035 | [−0.0028, +0.0104] |
| | soil_shallow hold-out | 0.9157 | 0.9186 | +0.0029 | [−0.0084, +0.0105] |
| | gut hold-in / hold-out | 0.9875 / 0.9931 | 0.9863 / 0.9919 | −0.0012 / −0.0012 | [−0.0029, −0.0002] / [−0.0044, −0.0000] |
| se | held out, design | 0.9649 | 0.9651 | +0.0002 | [−0.0006, +0.0012] |
| | held out, scenarios | 0.9519 | 0.9522 | +0.0003 | [−0.0003, +0.0009] |
| | test, design | 0.9617 | 0.9613 | −0.0005 | [−0.0016, +0.0007] |
| | soil hold-out | 0.9148 | **0.9319** | **+0.0171** | [+0.0021, +0.0349] |

The affected samples (thin: present species with 1-2 marker reads):

| sample | thin | thin missed, as built → without | median p of thin | FP |
|---|---:|---|---|---|
| pe soil 4 | 823 | 675 → 621 | 0.15 → 0.30 | 47 → 60 |
| pe soil_shallow 2 | 598 | 454 → 400 | 0.18 → 0.41 | 48 → 63 |
| se soil 4 | 691 | **546 → 270** | 0.10 → 0.70 | 81 → 191 |

- **For se, `sample_log_taxa` is half of the problem.** Without it, se's soil sample 4 recovers half of its thin
  species, for 110 more false positives in that sample. se's soil hold-out F1 rises by 0.017, and nothing else moves
  beyond noise.
- **For pe, it is a small part.** Thin species' scores in the two affected samples only double (0.15 → 0.30) and
  their misses fall by 8-12%. pe's model still finds the prior in other features of a deep, crowded sample:
  `sample_log_fragments` and the row features that grow with the crowd (`failed_candidate_rate`, the `relative_*`
  features; section 2).
- **Dropping it costs nothing measurable elsewhere**, except pe's gut samples (−0.0012: 9 more false positives in
  the 10 hold-in samples, 3 more misses in the 4 hold-out ones), and pe's test set gains 0.0023. The crossed-fold
  estimate (the trainer's headline) is −0.0024 for pe (0.9516 → 0.9492, rows only); species held out alone is level
  (pe 0.9555 both, se 0.9563 → 0.9565).
- So removing one feature is not the fix. The training samples must contain deep, rich, uneven communities, so that
  the depth and crowding features stop meaning "a thin row is spill" there.

## Recommendations

1. **Teach the models deep, species-rich, uneven samples.** The scenarios tie evenness to the sample number and draw
   richness and depth apart, so the training draw can miss whole corners. At 10 samples per scenario, it did.
   Options, cheapest first:
   - draw each scenario sample's sigma at random rather than in turn, and add samples that are deep, rich and uneven
     on purpose (for example 2 of the 10 soil samples at sigma 2.5 with 8,000-11,000 species and 1-2× depth);
   - give the model a feature that describes the sample's tail directly, such as the share of the sample's rows
     with 1-2 fragments, so it does not have to infer the tail from depth × richness;
   - dropping `sample_log_taxa` alone (section 8) recovers half of se's lost thin species and little of pe's. It is
     neutral elsewhere, so it can go with either of the above, but it is not a fix by itself.
2. **Split the final check by genome source** (representative, real strain, in-silico strain) and keep its differing
   rows (sample, taxon, p, knob) in the archive. The real-strain rows show the leak; the representatives and in-silico
   strains are the honest gain; the 18 rows show whether the parity issue is rounding.
3. **Decide on the allele scores.** No F1 at r226 in either direction. Keep them if they cost little; measure their CPU
   first.
4. **Drop `--model-gates` from the default, or profile only the test collection at those gates** (v29's
   recommendation 2, still open). They are 70% of the alignment CPU and add nothing as models.
5. **If the grid's sparse index is to be judged**, run a build with `--test-c-gates` (it brings `--test-grid`) and a
   real sparse-index run (`measure_performance.sh --test_sparse`).
6. **Look at the spike-mode false trigger**: add `*.spike_genes.tsv` to `--share-logs`, and see which genes of the gut
   PacBio sample piled.
7. Still open: PacBio's depth bias; se C 8 as the se database gate (free in v29 and v31).

The website needs no change from this report: no default changed.
