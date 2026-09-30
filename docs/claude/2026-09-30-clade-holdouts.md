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

## Follow-up: every rank held out, false positive and false negative rates by rank

Same day, after `0332a09`. Changes:

- Defaults `--holdout-clades phylum:2,class:4,order:6,family:8,genus:12` (orders and genera too) and
  `--holdout 0.2` (was 0.1): many more species the training database lacks, since false positives
  from organisms a database lacks are the main risk. Each design point takes the held-out clades of
  each rank in turn (every clade is used before one is used again) instead of at random.
- `collect_training_data.py` adds `meta_neighbour_rank`: for a present taxon, the deepest rank it
  shares with another species of its sample.
- The trainer's section "False positives and false negatives by taxonomic rank", rows from species
  to phylum: FP rates of absent taxa closest to species the database lacks, by the rank held out
  (and of the other absent taxa); FN rates of present taxa by `meta_neighbour_rank`; FN and FP rates
  with the taxon's species, genus, family, order, class or phylum held out of training
  (cross-validation, now with orders). The summary repeats the rates.
- **Out of bag.** scikit-learn 1.9 draws each tree's rows with probability proportional to their
  class weight (`class_weight="balanced"`). With many more absent than present taxa, a present row is
  drawn for a tree with 96-98% probability, and some are drawn for all 64 trees: they have no
  out-of-bag score, and scikit-learn reports 0 for them. In run G2, 128 of 1,436 present rows (9%)
  had out-of-bag probability 0 while their species-held-out one was 0.8-1.0, so the out-of-bag
  estimate said F1 0.951 against 0.994-0.996 for every cross-validation. The trainer now leaves those
  rows out of the out-of-bag estimate and says how many (retrained on G2: 123 rows left out, out-of-bag
  F1 0.995). Earlier reports'
  out-of-bag rows are affected where absent taxa far outnumber present ones; the cross-validations
  are not.

Runs (tuning world, 8 threads): **G2** as G (8 samples per point, 1,000-100,000 read pairs, 10-40
species) with the new holdout defaults; **H** with every default of the time (12 samples,
5,000-500,000 read pairs, 20-50 species). Both hold out 2 phyla (14 species), 4 classes (37), 6 orders
(33), 8 families (47), 12 genera (48) and 117 single species: 296 of the 765 species, a large share of
this small world (at GTDB scale a few percent).

H's false positives from what the training database lacks (FP rate, species held out, knob 0.5):

| held out at | simulated | absent taxa closest to them | FP | FP rate | collection model |
|---|---|---|---|---|---|
| species | 676 | 6,154 | 7 | 0.11% | 0.31% |
| genus | 414 | 2,839 | 0 | 0% | 0% |
| family | 437 | 1,322 | 0 | 0% | 0% |
| order | 333 | 126 | 0 | 0% | 0% |
| class | 351 | 12 | 0 | 0% | 0% |
| phylum | 205 | 0 | - | - | - |
| other absent taxa | - | 7,691 | 0 | 0% | 0% |

H's false negatives by the rank a present taxon shares with its closest other species in the
sample: genus 0.31% (2/641), family 0.29%, order 0.46%, class 0.35%, phylum and domain 0% (the
collection model: 3.0%, 9.8%, 6.0%, 2.1%, 12.1%, 11.2%). With the taxon's clade held out of
training, species to phylum: FN rate 0.19-0.34%, FP rate 0.03-0.04%, F1 0.997-0.998. G2's tables look
the same at about half the rows (FP rate 0.19% for species, 0% above; FN 0.4-1.4%).

Independent test sets (knob 0.5):

| test set | model | F1 | sensitivity | FP per sample | present found at 1,000 read pairs |
|---|---|---|---|---|---|
| `test` | C | 0.9969 | 0.9971 | 0.067 | 474/480 |
| | G | 0.9979 | 0.9988 | 0.058 | 478/480 |
| | G2 | 0.9971 | 0.9967 | 0.050 | 473/480 |
| | H | 0.9892 | 0.9815 | 0.058 | 440/480 |
| `test_congeners` | C | 0.9952 | 0.9956 | 0.108 | 499/509 |
| | G | 0.9950 | 0.9960 | 0.125 | 502/509 |
| | G2 | 0.9954 | 0.9952 | 0.092 | 498/509 |
| | H | 0.9868 | 0.9770 | 0.067 | 456/509 |

G2 (every rank held out, 20% of the species) is as good as C and G, with the fewest false positives
on the congener set. H lost its sensitivity at the shallowest depth only (from 10,000 read pairs
up, all models find the same): the default depths started at 5,000 read pairs, so the model had
never seen samples as shallow as the test set's 1,000 and 3,000, and its best threshold there moved
to 0.32. Its own estimate (F1 0.997) could not show this, since all its samples were as deep. The
default depths are now `1000,5000,20000,100000,500000` (both scripts).

**H2**, every default now (12 samples of 20-50 species at 1,000-500,000 read pairs, 180 samples;
2,813 s on 8 threads, collection 2,556 s): the shallow samples are found again.

| test set | model | F1 | sensitivity | FP per sample | present found at 1,000 read pairs |
|---|---|---|---|---|---|
| `test` | H2 | 0.9965 | 0.9975 | 0.092 | 475/480 |
| `test_congeners` | H2 | 0.9946 | 0.9952 | 0.125 | 499/509 |

Its own report, FP rate by the rank held out: species 0.20% (13/6,636), genus 0%, family 0.05%
(1/1,849), order, class and phylum 0%, other absent taxa 0%; FN rate by the closest other species in
the sample: genus 0.26%, family 0.63%, order 0.35%, class 0.55%, phylum 0.80%, domain 1.44%; with the
clade held out of training, species to phylum: FN rate 0.46-0.52%, FP rate 0.05-0.08%, F1
0.995-0.996. The collection (shipped) model: species 0.33%, FN rates 3-12%. Out of bag: 344 of the
present rows were drawn for every tree and are left out. The parity check passed.

What it shows: whole genera, families, orders, classes and phyla held out are, in this world, not a
source of false positives the model misses (one of 4,743 taxa closest to them was called), while
missing single species are (0.2% of their congeners, 1.6 per 100 such species); the model rejects
relatives of missing clades as well as a model trained without them, and does as well on the
independent test sets (H2 vs C: F1 0.9965 vs 0.9969 and 0.9946 vs 0.9952, a few taxa). Whether
GTDB's more conserved marker genes let reads of missing classes reach more taxa, the r226 run's
tables tell.

## For the GTDB r226 run (updated)

The defaults hold out 2 phyla, 4 classes, 6 orders, 8 families and 12 genera (each at most 2% of the
species) and 20% of the other species to simulate, put one species of each held-out rank into every
sample, and simulate 180 samples of 20-50 species at 1,000-500,000 read pairs (collection about
twice as long as with the old defaults). In `model_logs/`: `holdout.txt`, the report's section
"False positives and false negatives by taxonomic rank" and its summary lines, and the clade rows of
the evaluation table.
