# Training the presence model: the new trainer, and harder training data for GTDB

- **Date**: 2026-09-29.
- **Code**: branch `audit-fixes` (`995c4f1`: Java-free trainer, species-held-out evaluation;
  `822f1d4`: simulator quota fix), and branch `audit-fixes-zstd` (worktree `../protal-merge`: merge
  of `zstd-compression`, then the harder training data and placeholder models of this report).
- **Machine**: WSL Ubuntu 24.04, 8 threads, 15 GB; scikit-learn 1.9.0.
- **Data**: synthetic GTDB releases of `scripts/mini_db/simulate_gtdb_release.py`, reads by
  `simulate_metagenomes` (ART), profiled by protal. All of it under `~/tune` and `~/audit4` in WSL;
  regenerate with the scripts in this folder (their paths point at those folders).

## Questions

1. What do the trainer's changes (no grid search, 64 trees, Java-free PMML, evaluation on species
   held out) change, and do protal and the trainer agree on every probability?
2. What must the simulated training data hold for a model that works on real samples profiled
   against all of GTDB, and which settings should the full r226 build use?

## Answers

1. The new trainer calls as well as the old procedure (within 1-2 taxa on independent test
   samples), in 1 s instead of 65-112 s, with a model 5-18 times smaller; protal's probabilities
   equal scikit-learn's bit for bit (the sklearn2pmml models differed by up to 0.008). Its estimate
   on species held out matches independent test samples.
2. Training data need **both** strains other than the database's references **and** species the
   database lacks. With neither (the workflow before this work), the model called 27% of the
   relatives of unknown species present, 6 false positives per sample, and rated itself F1 0.9998
   while scoring 0.870. Species held out without strains made it miss 20% of strains. With both
   (10% of species held out), F1 0.997 and 0.07 false positives per sample, and its own estimate
   matched. Defaults for `build_gtdb_database.py`: `--holdout 0.1`, 8 samples per design point,
   strains from NCBI (`gtdb_strain_genomes.py`); see "Settings for GTDB r226".

## 1. The trainer (harder 64-species world)

World (`hard_world_lineages.py`): 64 species in 16 genera (4 archaeal), strains 0.4-3% from their
representative, congeneric species 3-8% apart at the markers; 16 species (one per genus) left out
of the database. Training: 180 samples from 8 genera; test: 75 samples from all 16 (half the genera
never in training), 518 present and 1,099 absent taxa. Every model scored by protal
(`--profile_only` on the test SAMs). `trainer_comparison.sh`, `trainer_comparison.py`.

| model | F1 | archaea found | FP (75 samples) | file | protal run | training |
|---|---|---|---|---|---|---|
| shipped | 0.853 | 56% | 55 | 10.6 MB | 0.42 s | - |
| old trainer, GTDB workflow settings (normalized, 512 trees) | 0.993 | 97.6% | 1 | 2.8 MB | 0.15 s | 112 s |
| old trainer, defaults (all features, 256 trees) | 0.987 | 97.6% | 8 | 1.4 MB | 0.09 s | 56 s |
| new trainer (normalized, 64 trees) | 0.994 | 98.1% | 2 | 0.59 MB | 0.05 s | 1.3 s |

The old procedure's grid search kept 8 of 24 features; its threshold chosen on its own test rows
(0.145) gave 33 false positives on species held out instead of 6. The export: protal (cPMML), the
Python scorer of `model_pmml.py` and scikit-learn agree exactly on every taxon, including values on
every threshold's float32 rounding boundary; the sklearn2pmml models differed from scikit-learn by
up to 0.0078 for 1-5 taxa per test set (cPMML reads their `float` fields as doubles).

Found on the way: `simulate_metagenomes --taxon d__Archaea:2` filled samples with archaea up to
`--species_per_sample` (fixed in `822f1d4`; every earlier training set is skewed).

## 2. Harder training data (GTDB-like 900-species world)

World (`setup.sh`): `gtdb_like_lineages.py --species 900 --archaea 0.08` (324 genera, 61% with one
species, the largest 30; 78% of species have a congener), 3 genomes per species, strains 0.4-4%
from their representative (`--strain_divergence 0.002-0.02`), congeneric species 3-12% apart at the
markers (`--species_divergence 0.015-0.06`). 135 species (15%) are unknown to every database
(`make_release_p.py`), like the species GTDB lacks; the release the workflow sees has 765.

Designs (`designs.sh`): the full `build_gtdb_database.py` on that release, 8 samples x 15 design
points (3 read setups x 1k-100k read pairs), 10-40 species and 2 archaea per sample.

| design | strains (`--extra-genomes`) | species held out (`--holdout`) | other |
|---|---|---|---|
| A (the workflow before) | no | 0 | |
| B | yes | 0 | |
| C | yes | 0.1 | |
| D | yes | 0.25 | |
| E | no | 0.1 | |
| F | yes | 0.1 | `--congeners 3` |

Test sets: 120 samples each from all 900 species (strains two thirds of the time), profiled against
the full database: `test` (uniform draws) and `test_congeners` (3 species of one genus in every
sample). `tune_eval.py` scores every model with protal on them.

**Table 3: `test`** (10,082 taxa; 1,609 present from strains, 828 from representatives; 2,126
absent taxa congeneric with an unknown species in the sample). Knob 0.5.

| model | F1 | sens. | prec. | FP/sample | strains found | archaea found | congeners of unknown called | own F1 estimate |
|---|---|---|---|---|---|---|---|---|
| shipped | 0.930 | 0.875 | 0.992 | 0.15 | 0.865 | 0.358 | 0.9% | - |
| A | 0.870 | 0.998 | 0.771 | 6.01 | 0.996 | 0.998 | 26.9% | 0.9998 |
| B | 0.953 | 0.999 | 0.911 | 1.98 | 0.998 | 0.998 | 10.1% | 0.999 |
| C | 0.997 | 0.997 | 0.997 | 0.07 | 0.996 | 0.995 | 0.3% | 0.996 |
| D | 0.998 | 0.996 | 1.000 | 0.008 | 0.994 | 0.995 | 0.05% | 0.993 |
| E | 0.929 | 0.868 | 1.000 | 0 | 0.800 | 0.902 | 0% | 0.9995 |
| F | 0.997 | 0.997 | 0.996 | 0.08 | 0.995 | 0.990 | 0.3% | 0.995 |

**Table 4: `test_congeners`** (9,255 taxa).

| model | F1 | sens. | FP/sample | strains found | archaea found | congeners of unknown called |
|---|---|---|---|---|---|---|
| shipped | 0.930 | 0.873 | 0.10 | 0.858 | 0.365 | 0.8% |
| A | 0.891 | 0.995 | 4.96 | 0.994 | 0.989 | 27.8% |
| B | 0.956 | 0.997 | 1.84 | 0.996 | 0.995 | 11.4% |
| C | 0.995 | 0.996 | 0.11 | 0.994 | 0.989 | 0.8% |
| D | 0.995 | 0.992 | 0.03 | 0.989 | 0.981 | 0.2% |
| E | 0.933 | 0.875 | 0 | 0.813 | 0.885 | 0% |
| F | 0.995 | 0.996 | 0.11 | 0.993 | 0.986 | 0.8% |

What they show:
- **Without strains or missing species there are almost no negatives.** A's first design points
  had 197-225 present and 0-1 absent taxa: reads of a species' own reference hardly land elsewhere.
  The model learnt "reads, so present", and its own estimate (0.9998) said nothing about real
  samples.
- **Strains alone** (B) give negatives (a strain's reads fit relatives nearly as well), but not the
  ones that matter most: reads of species the database lacks, which land on their congeners.
- **Species held out alone** (E) teach the model that reads matching the reference less than
  perfectly belong to a missing species: it rejects a fifth of the strains.
- **Both** (C, D): F1 0.995-0.998, and the trainer's estimate on species held out matches. 25%
  held out removes the last false positives but loses a little sensitivity (strains, archaea); 10%
  is the default, 0.25 for users who value precision.
- **Congeners in the training samples** (F: 3 species of one genus in each) changed nothing, not
  even on the congener test set: the relatives' cross-mapping that matters is already in C's data.
  `--congeners` stays off by default.
- In C's own report, log loss on species held out still fell from 60 to 120 training samples
  (0.0135, 0.0104): more samples help. Trees: 128 scored slightly better than 64 (log loss 0.0094
  vs 0.0104, F1 equal); the leaf limit did not bind (37 leaves per tree).

## Settings for GTDB r226

1. Strains: pick and download other genomes of GTDB species, and a pool of species to simulate:

       python3 scripts/gtdb_strain_genomes.py --gtdb GTDB_R226 -o strains --download

   4,000 species with up to 2 strains each and 1,000 simulated from their representative only (the
   simulator draws species uniformly, so without the pool the few species with strains would hardly
   be drawn among ~130,000); a simulated species is then another strain than the reference about 53%
   of the time (67% in the designs above). Needs the NCBI `datasets` CLI; about 8,000 genomes.
2. Build and train:

       python3 scripts/build_gtdb_database.py --gtdb GTDB_R226 --outdir OUT \
           --extra-genomes strains/genomes --simulate-species strains/simulation_species.txt \
           --protal build/protal --simulator build/simulate_metagenomes -t 32

   with the defaults `--holdout 0.1 --samples 8`. The training database (`OUT/training_db`) costs a
   second index build and its disk space.
3. Upload `OUT/model_logs/`. Check, in order: `genome_table.txt` (strain share), the report's
   per-domain counts (about 2 archaea per sample), "present, simulated from" and "absent taxa"
   (strains found; congeners of held-out species called), species-held-out F1 against the shipped
   model, the learning curve (more samples?), the trees table, and `parity.txt`.

## Placeholders for the other read types

A database holds `model_pe.xml`, `model_se.xml`, `model_PB.xml` and `model_ONT.xml`
(`zstd-compression`). Only paired-end reads can be simulated and aligned for training so far, so
`build_gtdb_database.py` packs placeholders (`scripts/placeholder_models.py`) for the other three:
they score every taxon 0, so no species is reported (`--knob 0` lists every taxon with reads, to
test the plumbing), and protal warns whenever it loads one (`<Annotation>protal:placeholder`).

## Next: training data for single-end, PacBio and ONT reads

Everything above is reusable per read type: genome tables and the strain pool, species held out
and the training database, and each design point's communities (`sim/manifest.tsv`) and truth
files. Only the reads and their alignment change:

- **se**: the R1 files of the paired-end samples (no new simulation), aligned by protal's
  single-end support (uncommitted in the main checkout as of this report).
- **pb**: each community replayed at the same total bases with a long-read simulator (PBSIM3; the
  PacBio branch's `simulate_hifi.py` as a start), HiFi and CLR; aligned by the native long-read
  aligner of branch `pacbio-long-reads`. FASTQ, since the profiler skips records without qualities.
- **ont**: replayed likewise at several accuracies (R9.4.1 to R10.4.1 SUP); protal cannot align
  them yet (chaining assumes few indels). minimap2 SAMs would do, but a model trained on them fits
  minimap2 alignments only.

A collector step `collect_training_data.py --replay PE_TRAINING --read_type T` would produce the
table, the trainer is unchanged, and `protal --add_model MODEL --read_type T` stores the model. Not
installed yet: PBSIM3, minimap2 (only `art_illumina` is). The read-type names differ between
branches (`zstd-compression`: `--read_type pe|se|pb|ont`, `model_PB.xml`; `pacbio-long-reads`:
`--read_type short|pacbio`, `model_pacbio.xml`; the single-end work: detected per sample,
`--model_se`) and must be reconciled before these branches merge.

## Follow-up (2026-09-29, later): one download stage

`gtdb_strain_genomes.py` became `scripts/download_gtdb.py`, which also fetches GTDB's own files
(r226 by default, `--release` 207 to 232) and the representatives of the simulated species from
NCBI, so the 127 GB archive of GTDB's representative genomes is not needed. The commands of
"Settings for GTDB r226" are now

    python3 scripts/download_gtdb.py -o /shared/protal_inputs/gtdb_r226            # node with internet
    python3 scripts/build_gtdb_database.py --inputs /shared/protal_inputs/gtdb_r226 --outdir OUT \
        --protal build/protal --simulator build/simulate_metagenomes -t 16

The read-type names were reconciled in `26f73db` (`pe|se|pb|ont`, `model_PB.xml`, `model_ONT.xml`).
