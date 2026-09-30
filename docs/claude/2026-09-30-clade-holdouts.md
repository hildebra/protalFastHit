# Training the presence model with families, classes and phyla the database lacks

- **Date**: 2026-09-30.
- **Code**: branch `audit-fixes` at `56217e1` with the uncommitted changes of this report
  (`scripts/build_gtdb_database.py`, `collect_training_data.py`, `random_forest_cmdline.py`, the new
  `scripts/lineages.py`, defaults of `download_gtdb.py`) and the build speed-ups of
  [the index-build report](2026-09-30-index-build-gains/README.md) in the binaries.
- **Data**: the GTDB-like tuning world of [the tuning study](2026-09-29-model-training-tuning/README.md):
  765 species in the release the workflow sees (17 phyla, 31 classes, 53 orders, 128 families, 295
  genera; 3 genomes each, strains 0.4-4% from their representative), and its independent test sets
  (`test`, `test_congeners`: 120 samples each from all 900 species, 135 of which no database has).
- **Machine**: WSL Ubuntu 24.04, Intel Core Ultra 7 258V, 8 threads, 15 GB.
- **Run** (`G`, design C of the tuning study plus the new clade holdouts):

      python3 scripts/build_gtdb_database.py --gtdb ~/tune/release_p --outdir ~/tune/G \
          --protal PROTAL --simulator SIMULATE_METAGENOMES -t 8 --samples 8 \
          --read-pairs 1000,3000,10000,30000,100000 --species-per-sample 10-40 --archaea 2 --seed 1 \
          --extra-genomes ~/tune/world/simulation/genomes_nonreps

  with the new defaults `--holdout 0.1 --holdout-clades phylum:2,class:4,family:8
  --holdout-max-share 0.02 --novel-clades-per-sample 1`; then `tune_eval.py ~/tune PROTAL test C F G`
  (and `test_congeners`) from the tuning study's folder.

## What changed

- **Whole clades held out.** Before, the training database lacked a random 10% of the species. Now
  `build_gtdb_database.py` first leaves out whole clades (`--holdout-clades`, default 2 phyla, 4
  classes, 8 families), drawn from phylum down among clades with at least two species to simulate,
  at most 2% of the database's species, and not inside a clade drawn before; then 10% of the other
  species. `heldout_species.txt` has each species with its rank and clade; `model_logs/holdout.txt`
  sums them up.
- **In every sample.** Species are drawn uniformly, so the few of held-out clades would hardly be
  in a sample. `collect_training_data.py --novel_clades 1` draws one held-out clade per rank for each
  design point and asks the simulator for one of its species in every sample (`--taxon`).
- **Attributed.** For every absent taxon the collector finds the species simulated in its sample
  that is closest to it; when that is a species the database lacks (at least as close as any present
  species), `meta_novel_level` is the rank it was held out at, `meta_relative_rank` the rank they
  share. `meta_novel_levels` counts the sample's novel species by rank.
- **Trained on and reported.** The rows are training rows like any other (negatives). The trainer
  adds a section "Species the database lacks, by the rank they were held out at" (novel species
  simulated, absent taxa closest to them, called at the knob with species held out and by the
  collection model, per 100 novel species), and cross-validation that holds out whole genera,
  families, classes and phyla (with `--taxonomy`, which the workflow passes).
- **More species.** Defaults: `--samples 12` (was 8; the tuning study's learning curve was still
  falling at 120 samples), `--species-per-sample 20-50` (was 10-40), and `download_gtdb.py --species
  6000 --rep_only_species 2000` (was 4000 and 1000; about 20,000 genomes instead of 13,000). Run G
  kept the old sample settings so that it compares with design C.

## Results

Holdout: 2 phyla (14 species), 4 classes (37), 8 families (53) and 66 single species, 170 of 765.
Training table: 8,090 taxa in 120 samples (2,010 present).

Species the database lacks, by the rank they were held out at (trainer report, species held out,
knob 0.5):

| held out at | simulated | absent taxa closest to them | called | called by the collection model |
|---|---|---|---|---|
| species | 223 | 1,775 | 10 | 11 |
| family | 310 | 652 | 0 | 0 |
| class | 248 | 9 | 0 | 0 |
| phylum | 156 | 0 | 0 | 0 |

By the rank the absent taxa share with that species: genus 967 (8 called), family 576 (2), order
751, class 132, phylum 10 (none called).

Whole clades held out in cross-validation (F1, false positives per sample): species 0.995 / 0.08,
genus 0.994 / 0.09, family 0.995 / 0.07, class 0.995 / 0.06, phylum 0.996 / 0.05. The forest does
as well on clades it never saw.

Independent test sets (profiled against the full database, knob 0.5), against design C (the same
without clade holdouts) and F (C with congeners):

| test set | model | F1 | sensitivity | FP per sample | congeners of unknown species called |
|---|---|---|---|---|---|
| `test` | C | 0.9969 | 0.9971 | 0.067 | 0.34% |
| | F | 0.9965 | 0.9967 | 0.075 | 0.34% |
| | G | 0.9979 | 0.9988 | 0.058 | 0.25% |
| `test_congeners` | C | 0.9952 | 0.9956 | 0.108 | 0.76% |
| | F | 0.9952 | 0.9956 | 0.108 | 0.76% |
| | G | 0.9950 | 0.9960 | 0.125 | 0.88% |

Time: 535 s end to end on 8 threads (training database 25 s, collection 388 s, training 114 s);
the finished database was built in the background meanwhile (530 s, slowed by the collection; about
2 minutes alone).

## What it shows

- Clade holdouts cost nothing on species-level novelty: G is as good as C on both test sets (the
  differences are one or two taxa).
- In this world, reads of a species whose family the database lacks reach 652 taxa of other
  families (sharing order or class), and the model rejects all of them; reads of missing classes and
  phyla hardly reach any taxon, because they match nothing at protal's identity thresholds. The
  synthetic world's divergence between families may be larger than GTDB's at the most conserved
  marker genes, so the same table from the r226 run decides whether clade holdouts matter there.
- Samples need the forced species: with uniform draws, the 104 species of held-out clades would
  have come up about 3 times per sample in total, with no guarantee per rank.

## For the GTDB r226 run

Nothing extra to pass: the defaults hold out 2 phyla, 4 classes and 8 families (each at most 2% of
the species) besides 10% of the species, and put one species of each held-out rank in every
sample. In `model_logs/`: `holdout.txt` (which clades), the report's section on species the database
lacks (false positives per 100 novel species by rank, and by shared rank) and the clade rows of the
evaluation table; `build_metadata.tsv` records the clade counts.
