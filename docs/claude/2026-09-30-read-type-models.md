# One presence model per read type, an independent test set, and samples like real ones

- **Date**: 2026-09-30.
- **Code**: branch `audit-fixes` at `5385719` (the changes are in `3dec65d`, committed as "stash"
  before the merges of `strain-fixes` and `performance`), with the timing message of this report's
  commit. Binaries built from `5385719` (unit tests 183/183, mini database, e2e 103/103).
- **Data**: the GTDB-like tuning world of [the tuning study](2026-09-29-model-training-tuning/README.md)
  (765 species in the release, 3 genomes each) and pbsim3 3.0.x from the `protal-db-build`
  environment.
- **Machine**: WSL Ubuntu 24.04, Intel Core Ultra 7 258V, 8 threads, 23 GB.
- **Run** (`V2`, the defaults at a reduced size: 4 samples per design point, 3 depths):

      python3 scripts/build_gtdb_database.py --gtdb ~/tune/release_p --outdir ~/tune/V2 \
          --protal PROTAL --simulator SIMULATE_METAGENOMES -t 8 --seed 1 \
          --extra-genomes ~/tune/world/simulation/genomes_nonreps --samples 4 \
          --read-pairs 1000,20000,200000 --test-samples 2 --test-read-pairs 500,10000,500000 \
          --long-read-bases 300000,6000000,60000000 --test-long-read-bases 150000,3000000,90000000 \
          --pbsim ~/micromamba/envs/protal-db-build/bin/pbsim

## What changed

From the list of [the clade-holdout report](2026-09-30-clade-holdouts.md)'s follow-up:

1. **An independent test set in every build.** A second collection of another design and seed
   (`--test-samples 4` per point; 500 to 1,000,000 read pairs, 10-300 species, Poisson-lognormal
   sigma 2.0 instead of 1.3, strains 0.5/0.2), profiled against the training database. Each model
   scores it (`random_forest_cmdline.py --test-file`): the report's section "Independent test set"
   (metrics, FN and FP rates by depth and by rank, the threshold with the highest F1) and a warning
   when it scores clearly below cross-validation.
2. **Complex communities**: `--species-per-sample 20-200` (was 20-50).
3. **Co-existing strains**: `--strains-per-species 0.3,0.1` (a second strain in 30% of species, a
   third in 10%). `simulate_metagenomes` failed a sample when a rare species had fewer read pairs
   than strains; such a species now keeps as many strains as it has read pairs
   (`CommunityProfileDesigner.cpp`, `tests/test_CommunityDesign.cpp`).
4. **One model per read type**, trained in parallel (`--read-types pe,se,pb,ont`): single-end from
   the paired-end samples' first reads; PacBio and Nanopore from pbsim3 reads of the same communities
   (the paired-end point's manifest, each genome's share of `--long-read-bases` by abundance times
   length); all samples profiled in one protal run (the map's `READ_TYPE`); one table per read type;
   `check_model_parity.py --read_type`; each model stored with `--add_model --read_type`.
5. **Sequencer profiles**: HiSeq X (HSXt, ART's closest to NovaSeq) for 150 bp instead of HiSeq
   2500; a read setup `file=R1.txt+R2.txt` gives ART quality profiles of `art_profiler_illumina`
   (checked: ART then ignores the built-in profile and simulates 150 bp reads with a 100 bp one
   named).
6. **Provenance** in `build_metadata.tsv`: date, protal version and binary, the scripts' commit,
   the command, seed, genome table, holdout, training design, read types, and each model's F1 on
   species held out and on the test set.

## Results

1,097 s end to end: collection of 96 training samples (36 pe, 36 se, 12 pb, 12 ont) in 502 s, of
48 test samples in 368 s, the four models trained in parallel in 203 s, the finished database built
in the background (1,090 s alongside). The parity check passed for every read type (protal's
probabilities equal the model files', the features those of the collection), and the database holds
the four trained models.

| model | training taxa | F1, species held out | FP per sample | F1, test set | test FN rate at its shallowest depth |
|---|---|---|---|---|---|
| pe | 7,365 | 0.988 | 0.67 | 0.977 | 14.5% (500 read pairs) |
| se | 6,720 | 0.980 | 1.31 | 0.963 | 21.5% (500 reads) |
| pb | 1,126 | 0.979 | 0.92 | 0.980 | 9.4% (150 kb) |
| ont | 1,478 | 0.982 | 0.92 | 0.951 | 22.0% (150 kb) |

At the test set's other depths the FN rates were 0-6.8% and the FP rates 0.4-3%. With the clade
held out of training, species to phylum, F1 stayed within 0.003 of the species-held-out value for
every read type. The collection models of se, pb and ont are the placeholders (F1 0), which the
reports list for comparison.

## What it shows

- The whole workflow runs for all four read types, and the test set does what it is for: its
  shallowest depths are below the training design's (500 read pairs against 1,000; 150 kb against
  300 kb), the models miss taxa there, and the report warns. Cross-validation could not show it. At
  full size the defaults train from 1,000 read pairs and test from 500, so the warning will show how
  far below the training depths the models still hold.
- ONT is the weakest at shallow depths: on the test set its best threshold was 0.14, not 0.5, i.e.
  its probabilities are too low for taxa with little evidence. More long-read training samples at
  shallow depths (the defaults have 5 long-read depths and 12 samples each, 4 times this run) are
  the first thing to check in the r226 reports.
- Single-end reads trail paired-end ones (more FP, lower sensitivity), as expected without the mate.
- This run is small (4 samples per point, 12 long-read samples per read type), so its numbers say
  that the workflow works, not how good the models are at full size.

## For the GTDB r226 run

The defaults give 180 paired-end samples of 20-200 species (also profiled as single-end), 60 PacBio
and 60 Nanopore samples, and a test set of 72 paired-end samples with their single-end and long-read
counterparts. Collection takes several times this run's (4 times the samples, deeper points, larger
genomes); memory: the finished database's build in the background plus the collection's protal runs
(`--one-build-at-a-time` if the node is short). Check in `model_logs/`: each read type's report
(`trained_model*.report.txt`: the test-set section and its warning, FN rates by depth, the rank
tables), `parity*.txt`, and `build_metadata.tsv`.
