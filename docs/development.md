# Development: tests, mini databases and simulated metagenomes

```bash
just test          # C++ unit tests (GoogleTest)
just e2e           # builds the mini database, then the end-to-end tests against it
just example       # mini database + reads from a known mock community: is the profile right?
just mini-db-test  # Python tests: the mini database scripts, the GTDB build pipeline, in-silico strains, tracing, ...
just model-test    # the trainer, and its model's PMML export scores as scikit-learn does
```

The unit tests need GoogleTest (`libgtest-dev` on Ubuntu). The Python tests need numpy
(`python3-numpy`), `model-test` and the GTDB build test also pandas, joblib and scikit-learn, as
training does, and several tests the `zstd` CLI.
Build requirements are in [installation.md](installation.md#building-from-source).

## Unit tests

`tests/test_*.cpp`, one file per unit (the file's first lines say what it covers), build as one
binary, `protal_tests`: about 510 tests, 10 seconds on 4 cores. Most compare protal's code with an
independent oracle: libzstd and zlib-ng's readers, WFA2 alignments, brute-force definitions (syncmers,
flex neighbours, the alignment screen), the text tables against the binary gene table, serial against
parallel profiling. Their seeds are fixed. The helpers they share are in `tests/TestUtil.h` (scratch
directories, files, random and mutated sequences), `tests/TestReference.h` (taxa of genes written as
`reference.fna` and `reference.map`, and loaded) and `tests/TestSamSink.h`.
`tests/data/golden_model_rules.tsv` holds golden vectors of the rules protal and the trainer both
implement (the prior adjusted to a sample, the calls at a share of false calls, the knob by the
sample's depth): `test_GoldenModelRules.cpp` checks protal on them, `scripts/test_model_pmml.py` the
trainer.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON
cmake --build build --target protal_tests
ctest --test-dir build --output-on-failure -j4
```

Under AddressSanitizer and UndefinedBehaviorSanitizer, as CI runs them (a Debug build, so `assert()`
is on, about 3 minutes on 4 cores). `-DPROTAL_NO_CLONES` compiles the hot functions once, for plain
x86-64 (`src/Utilities/TargetClones.h`, and WFA2-lib's, `lib/wfa2-lib.cmake`), so that this run tests
the copies a CPU without AVX2 runs; the Release build on a machine with AVX2 runs the others:

```bash
FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer"
cmake -S . -B build-asan -G Ninja -DCMAKE_BUILD_TYPE=Debug -DPROTAL_BUILD_TESTS=ON \
  -DCMAKE_CXX_FLAGS="$FLAGS -DPROTAL_NO_CLONES" -DCMAKE_C_FLAGS="$FLAGS" -DCMAKE_EXE_LINKER_FLAGS="$FLAGS"
cmake --build build-asan --target protal_tests
ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1 ctest --test-dir build-asan --output-on-failure -j4
```

The AVX-512 kernels' tests (`FlexScan.Avx512ScoresAndTiesAsTheScalarScan`,
`PackedIndex.LookupsGiveTheSameSeedsAtEveryVectorLevel`, the `Syncmers` tests) run only on a CPU with the AVX-512
level of `src/Utilities/SimdLevel.h` (Intel Ice Lake, AMD Zen 4 and later); elsewhere they skip it or print
"no AVX-512 here". Intel's Software Development Emulator runs them on any x86-64 CPU:

```bash
sde64 -icx -- build/tests/protal_tests --gtest_filter='FlexScan.*:Syncmers.*:PackedIndex.*'
```

Benchmarks are not unit tests: the flex-cell scan's is
[`docs/claude/2026-10-06-performance-profiling/scripts/flex_scan_bench.cpp`](claude/2026-10-06-performance-profiling/scripts/flex_scan_bench.cpp).

WFA2-lib is built without UBSan (`lib/wfa2-lib.cmake`): its unaligned loads and shifts of negative
offsets would stop every test that aligns. zlib-ng is built with both sanitizers.

## End-to-end tests

`tests/e2e/test_protal_e2e.py` simulates reads from the mini database's reference genes and runs the
real `protal` and `simulate_metagenomes`: exit codes, output files, SAM records, strain MSAs, reruns,
failures reported and kept to their sample, a profile's truth counts and abundances, its composition
(the reads scanned, the genome sizes, the unknown share), a profile of the unknown alone (`? ? 1`)
for reads of nothing in the database, every default feature group the model may use, and a small
gradient-boosted model (`tests/e2e/data/model_gbm_small.xml`) scored end to end against probabilities
computed by hand. It covers paired-end, single-end, PacBio and ONT reads and phasing.

The tests are written for the mini database (`just mini-db`): they name its species, taxids and
genes, so another database fails them. They need Linux, numpy and the `zstd` CLI (or Python 3.14),
and take about a minute on 4 cores, 3.5 GB of memory and 3-4 GB of `/tmp`:

```bash
PROTAL_TEST_DB=data/mini_db/protal_db PROTAL=build/protal SIMULATE=build/simulate_metagenomes \
    PROTAL_TESTS_REQUIRED=1 python3 -m unittest -v tests/e2e/test_protal_e2e.py
```

A test whose prerequisite is missing (the database, a binary, the `zstd` CLI) is
skipped; `PROTAL_TESTS_REQUIRED=1` makes it an error instead, as CI runs them. The tests that need no
database (the version, the simulator, qcmsa's contract, small builds, gene neighbours, the WFA2
traceback) run without `PROTAL_TEST_DB`.

## Script tests

`scripts/mini_db/test_*.py` test the mini database scripts (`test_mini_db.py`), the GTDB downloads
(`test_downloads.py`), the training data collector and its scenarios (`test_collector.py`), the gene
neighbours (`test_gene_neighbours.py`), the build script's parts (`test_gtdb_build.py`), the GTDB
build end to end (`test_gtdb_pipeline.py`) and the ancestry true-positive test's worlds, oracle and
quick run (`test_ancestry_world.py`); `scripts/test_*.py` the trainer and the model's PMML export
(`test_model_pmml.py`), the in-silico strains, `trace_relatives.py`, `error_reads.py`,
`composition_accuracy.py`, the profile scripts, the strain test scripts, `ancestry_sites.py` and
`foreign_rates.py`; `scripts/mock_benchmark/test_mock_benchmark.py` the
[benchmark on public mock communities](#benchmark-on-public-mock-communities) (downloads from a local
stand-in of ENA, scores, run folders, and runs of stand-ins of protal and sylph). A test whose prerequisite is missing is skipped; `PROTAL_TESTS_REQUIRED=1` makes it
fail instead (`scripts/prerequisites.py`). With everything present they run in about 6 minutes on 4
cores, the GTDB build 3 of them:

```bash
PROTAL=$PWD/build/protal SIMULATE=$PWD/build/simulate_metagenomes PROTAL_TESTS_REQUIRED=1 \
    python3 -m unittest scripts/mini_db/test_*.py
python3 -m unittest scripts/test_insilico_strains.py scripts/test_trace_relatives.py scripts/test_error_reads.py \
    scripts/test_composition_accuracy.py scripts/test_profile_scripts.py scripts/test_strain_scripts.py \
    scripts/test_ancestry_sites.py scripts/test_foreign_rates.py scripts/test_model_pmml.py \
    scripts/mock_benchmark/test_mock_benchmark.py
```

**The GTDB build end to end** (`test_gtdb_pipeline.py`, `GtdbBuildTest`) needs `$PROTAL`, `$SIMULATE`,
and scikit-learn, joblib and pandas in `$PROTAL_TRAIN_PYTHON` (default: the Python
running the tests). On a synthetic release of 60 species, downloaded from stand-ins of GTDB's mirror
and NCBI, each build serves every check of what it does:
- a build trained for pe and se with its genes ranked from the training database (`--rank-genes`),
  its samples' composition checked against the truth (`model_logs/composition_accuracy*.tsv`: the
  species' genome sizes are their genomes'), and a rerun that builds nothing;
- a reduced database of the 3 best genes (`--n-genes`) ranked from a full build of the training
  database, the same ranking as `--rank-genes`'s, its models with the relatives features and calls at
  a target share of false calls;
- a reduced database ranked by a given table (`--gene-ranking`) whose background build fails, which
  stops the run;
- another seed, which rebuilds only the training database, and `SIGTERM`, after which no command is
  left running;
- the default scenarios, made small by a `--scenario-file` of their names (one with 90% host reads,
  soil scaled down to the genome table), with Illumina reads at a target quality and Ultima reads,
  their hold-in and hold-out samples scored in every report and in `summary.txt`, the feature sets the
  trainers chose and why, the reads behind each model's errors (`model_logs/error_reads/`), and the
  host scenario left out without a host genome. The other builds pass `--scenarios none` and
  `--error-reads none`;
- the full build's samples profiled as they are simulated (`--profile-blocks` of a few kB), and
  streamed into protal (`--stream-above`): the same tables. Both run in copies of the full build
  without its samples, whose databases they take as they are;
- no species held out: the models trained on a training database of their own.

The runs use 4 threads. `build_gtdb_releases.py` is tested with a stand-in of the build script
(`test_gtdb_build.py`).
The trainer's tests fit on one thread (`OMP_NUM_THREADS=1`): on every core of a busy machine,
gradient boosting made `test_model_pmml.py` take 20 minutes (scikit-learn 1.9.1). scikit-learn 1.8 and
later bins a weighted fit (the scenario rows weigh 0.25) by computing each of a feature's 254 weighted
percentiles on its own, ~10 times the time of an unweighted fit; the trainer puts
`scripts/hgb_binning.py` in its place, which computes them at once, bit for bit the same
(`HgbBinningTest`), and steps aside for any scikit-learn whose thresholds it does not reproduce. The
features of `TrainerScenariosTest`'s tables take 255 values or fewer, binned without percentiles.

## Mini databases

`scripts/mini_db/` builds a small database from a synthetic GTDB release (3 species, 3 genomes each
by default), to test without the full database. It is about 1.6 MB as `database.protal`
(`PROTAL_BUILD_ARGS=--no_bundle` or `--no_compress` for separate files; a raw `index.prx` is 3 GB
whatever the reference, since the k-mer key map has a fixed size):

```bash
just mini-db          # or: PROTAL=build/protal bash scripts/mini_db/build_mini_db.sh data/mini_db
protal --db data/mini_db/protal_db -1 r1.fq -2 r2.fq -o out/
```

`build_mini_db.sh` runs the same steps as a GTDB database ([databases.md](databases.md#step-by-step)):
the release, the converter, the gene neighbours, `protal --build` (on 4 threads at zstd level 3, ~4 s).
The scripts:

| Script | |
|---|---|
| `simulate_gtdb_release.py` | a directory laid out like an extracted GTDB release, plus `simulation/genomes.tsv` for `simulate_metagenomes` and `simulation/marker_positions.tsv` |
| `gtdb_like_lineages.py` | lineages shaped like GTDB's for `--lineages`: up to 999 species, most genera with one species, a share of archaea; the worlds the model's training was tuned on |
| `gtdb_to_protal_db.py` | the converter, for synthetic and real releases ([databases.md](databases.md#1-convert-the-release)) |
| `gene_neighbours.py` | the gene neighbours and positions ([databases.md](databases.md#gene-neighbours)) |
| `make_gene_rates.py` | the table of `--gene_rates r226` from a GTDB build's `gene_congeners.tsv`; rerun it on a newer release |
| `simulate_reads.py` | read pairs drawn straight from a mock community's genomes (substitutions only), with the truth |

`simulate_gtdb_release.py` evolves each marker from one random coding sequence down the lineage.
Its options:
- `--seed`, `--lineages FILE` (one GTDB lineage per line), `--genomes_per_species`,
  `--genome_length` (the background DNA, or a range `LOW-HIGH` from which each species' is drawn,
  log-uniformly, so that genome sizes differ), `--contigs`, `--marker_loss`;
- `--strain_divergence`, `--species_divergence`: a rate or a range `LOW-HIGH`, drawn per genome or
  species (`simulation/divergence.tsv`);
- `--gene_rates r226`: each marker evolves at the real r226 gene's speed (strains at its
  within-species factor, the lineage at its between-congener factor, each marker set's mean 1).
  protal's factors estimated on such a world match the real ones at Spearman +0.99 (archaea +0.97).
  Use it for training and benchmark worlds. `--gene_rates categories`, the older guess from gene
  names, reaches +0.54 (archaea +0.15)
  ([report](claude/2026-10-05-mini-database-genes/README.md));
- `--operons`: markers in clusters of up to 6 genes 0-150 bases apart, in the same order in every
  species (each family breaking some up, `--operon_breaks` 0.25), so that pairs and long reads span
  neighbouring genes. Without it, markers are shuffled and about 1.2 kb apart.

The mini genes are sized from the HMM model length, so the few Pfam markers that model one domain of a
longer protein are shorter than real ones.

`just example` runs [`examples/mini_db/run.sh`](../examples/mini_db/README.md): a seeded build of the
mini database, 30,000 read pairs from three genomes (two of them strains), a protal run, and a check
that every species is found with the right abundance (within 0.05).

## Simulating metagenomes

`simulate_metagenomes` draws mock communities from a table of genomes and simulates paired-end
Illumina reads for them with its own models of the instruments (no ART: [databases.md](databases.md#illumina-reads)),
and long reads of given communities (`--long_samples`). It makes the tests' and benchmarks' samples and the model's training data
([databases.md](databases.md#training-data)). Build it with the other binaries
(`cmake --build build --target simulate_metagenomes`, or `just simulate`).

**Input**: a tab-separated genome table without header: name (accession), GTDB taxonomy, FASTA path
(gzip or zstd work), and optionally the genome's length, which spares reading every genome at the start
(~5 ms per genome of GTDB's size). `build_gtdb_database.py` writes one with lengths
(`OUTDIR/work/genomes.tsv`), and synthetic releases have `simulation/genomes.tsv`. A first row is taken
for a header only if its fields name the columns (name, taxonomy, path, length) and none holds a
lineage (`;`), a path (`/`) or a number. A genome's FASTA is read whole and parsed in bulk each
time a sample needs it, unless a genome store (`--genome_store`) holds it.

```bash
./build/simulate_metagenomes --genome_table genomes.tsv --output_dir sims/ --samples 3 \
  --sample_prefix sim --total_read_pairs 100000 --species_per_sample 15 --distribution power_law \
  --strains_per_species "0.4,0.2" --seed 1 -t 8 --protal_metafile sims/protal
protal --db DB --map sims/protal.meta -t 8      # prints true and false positives per sample
```

| Output | |
|---|---|
| `reads/<sample>_R1.fq.gz`, `_R2.fq.gz` | the reads (BGZF, compressed as they are made, the genomes' in order; with `--reads_compression zstd` `.fq.zst`, zstd frames at level 3, as the database build writes them; with `--first_reads_only` no `_R2`, named nowhere: `-` in the protal map); files that are named pipes are written in place, a sample at a time, and a failed run ends them with a cut frame or record so that the reader fails too |
| `manifest.tsv`, `manifests/<sample>.tsv` | per (sample, genome): read pairs, relative abundance, vertical coverage (read bases over the genome's length; contigs shorter than a read get no reads, and the run notes a genome with 1% or more of its bases in them), the FASTA, the seed of its reads (`art_seed`) and its sample's seeds of the quality offset and the host reads (`run_seed`, `host_seed`) |
| `abundance_matrix.tsv` | relative abundance of each species in each sample |
| `run_params.tsv` | the command line, the seed (also when not given) and the read settings |
| `protal.meta`, `protal_goldstd/<sample>.profile_truth` | with `--protal_metafile DIR`: a protal map with a `PROFILE_TRUTH` column; protal then writes `<profile>.truth_annotated` |
| `plots/<sample>.png` | with `--plot_png` (needs `Rscript`) |

| Option | Default | |
|---|---|---|
| `-n, --samples`, `--sample_prefix` | 1, `sample` | samples, named `<prefix>_<n>` |
| `--total_read_pairs` | 100000 | read pairs per sample |
| `--species_per_sample` | 10 | a number, a range, e.g. `20-80` (each sample's drawn from it), or numbers given to the samples in turn, e.g. `20,80,45` (as `--total_read_pairs`) |
| `--distribution` | `lognormal` | `lognormal` (`--pln_mu`, `--pln_sigma`, several sigmas going to the samples in turn: a continuous long tail, no species below `--abundance_floor`, 0.001, of the median), `power_law` (`--alpha`), `negative_binomial` (`--nb_r`, `--nb_p`) or `poisson_lognormal` (Poisson counts + 1 of a lognormal mean, the default before 2026-10-07: at `--pln_mu` 0 about 40-45% of the species share the lowest weight). Every species gets a read pair or more |
| `--strains_per_species` | none | probabilities of a 2nd, 3rd, ... strain of a species, e.g. `0.4,0.2,0.1` |
| `--include_species`, `--genus`, `--taxon` | | species in every sample; `g__A:10,g__B:2` species from those genera, `d__Archaea:10` from any taxon (`--pick_random_demand_if_fail` caps instead of failing) |
| `--congener_groups` | none | `SHARE:MIN-MAX`, e.g. `0.25:2-5`: about SHARE of each sample's species in groups of MIN to MAX congeners |
| `--strain_sharing_file` | | strains shared across samples (below) |
| `--seed` | random | |
| `--read_length`, `--fragment_mean`, `--fragment_stdev`, `--sequencer` | 150, 350, 50, `HS25` | the reads, fragments and instrument (`HS20`, `HS25`, `HSXt`, `NovaSeq`, `MSv3`; `--illumina_report PAIRS` prints its qualities and errors) |
| `--mean_quality`, `--host_folder`, `--host_pairs`, `--first_reads_only` | | each read's qualities shifted to average this; a host genome (`scenarios.prepare_host`) and its read pairs per sample, after the community's; only the `_R1` files |
| `-t, --threads` | 1 | the reads are made in work items on all threads, several samples at once when one leaves threads idle; the same files for any number |
| `--test` | off | the design, manifests and truth, no reads |
| `--reads_compression` | bgzf | `bgzf` (`.fq.gz`, ISA-L at level 1) or `zstd` (`.fq.zst`: about 15% smaller, about half as fast to write; protal reads both) |
| `--plain_pipes` | off | outputs that are named pipes get plain FASTQ, whatever their names say (protal takes it as it is); regular files stay compressed |
| `--genome_store DIR` | none | genomes decoded once into `DIR` (2 bits a base, ~0.25 bytes a base, memory-mapped): the first run that needs a genome reads its FASTA and writes its file there, later samples, long-read rounds and runs map it; a FASTA newer than its file is read again; the same reads as without it |
| `--long_samples`, `--long_genomes`, `--long_setup`, `--long_model` | | long or Ultima reads of given communities ([databases.md](databases.md#one-model-per-read-type)); `--long_templates FASTA --long_out FILE`: one read of each sequence |
| `-v` | | the version and commit |

**Strains shared across samples.** `--strain_sharing_file` takes a tab-separated file, one species
per row: `SPECIES`, `SAMPLE_FRACTION` (share of samples with it), `N_STRAINS` (distinct strains
across samples), `MIN_OCCURRENCE` (samples each strain is in at least), and optionally `MIN_VCOV`
(minimum coverage per strain and sample) and `CONSPECIFIC_STRAINS` (as `--strains_per_species`).

**Replaying a dataset.** A manifest holds the genomes, read pairs and the seeds of their reads (`art_seed`)
and of each sample's quality offset and host reads (`run_seed`, `host_seed`), so `--from_manifest
sims/manifest.tsv --output_dir replay/` reproduces the reads byte for byte with the same simulator and
read settings, whatever `--seed`, and a per-sample manifest (`manifests/<sample>.tsv`) reproduces that
sample's. The read settings come from the command line (`--read_length`, `--sequencer`, ...,
`--host_folder`, `--host_pairs`); the replay warns about each that differs from the original
`run_params.tsv`. Sampling options are ignored. A manifest without `run_seed` and `host_seed` (before
2026-10-07) replays its genomes' fragments with this run's qualities, so other errors and reads; one
without `art_seed` (the ART days) replays the composition with fresh reads (pass `--genome_table`).

## Benchmark on public mock communities

`scripts/mock_benchmark/` scores protal on real reads of communities whose members are known. The datasets come from
the [survey of public mocks](claude/2026-10-13-public-mock-communities/README.md). This is unlike the builds' test
samples, which are simulated.

| File | |
|---|---|
| `samples.tsv` | 1,330 samples in 13 sets, one row per library: its runs, read type, truth table and replicate group (from `prepare_tables.py samples`, ENA's metadata) |
| `truth/*.tsv` | what each mock holds as GTDB r226 species, by kind (see below), with each member's strain (from `prepare_tables.py truth --gtdb ~/GTDB/r226`, which prints the evidence for each assignment) |
| `download_mocks.py` | fetches the runs from ENA |
| `score_mocks.py run`, `score` | runs protal (or sylph) per set, or per chunk of a set as a SLURM array, into a run folder per binary, and scores the profiles of every run |

**Sets.** `d6331_dna`, `d6331_host`, `meslier`, `treichel`, `d6331_hifi`, `atcc1003_hifi`, `d6331_ont`,
`d6322_basecalls`, `strain_mixes`, `nibsc_multilab`, `hcom2_invitro`, `hcom1a_dropouts`, `plant_syncoms`.

**Truth tables.** Their kinds:
- **`abundance`:** shares are known, as genome copies, since protal's abundances are shares of cells. A design by DNA mass is converted by genome size.
- **`closed`:** members are known but not their shares, so any other species called is false.
- **`dropout`:** strains left out by design must not be called.
- **`negative`:** a control without a community.

The shares add up to 1 together with a `?` row: the members protal cannot place, such as yeasts and strains the release lacks. The plant SynComs have no truth, only replicate groups.

**Downloading.** Run it on a node that reaches the internet. It resumes where it stopped.
- A sample keeps its first `--max-pairs` read pairs (default 20M) or, for long reads, its first `--max-bases` bases (default 20 Gb). A run is cut while it streams, so only what is kept is downloaded.
- A run within the limit is copied whole and checked against ENA's MD5.
- With these limits all sets take about 1.6 TB; `--dry-run` prints what each set would take.

**Running.** `score_mocks.py run` puts each run in a folder of its own, `$MOCKS/<run>/<set>/`, and skips samples that
already have a profile there.
- **protal** (the default): a map per set, run with `--no_strains`. The run folder is named from the binary's
  `--version`, `<version>_<commit>`, e.g. `0.8.3_1dcc325`, with `-modified` appended when it was built with uncommitted
  changes. Another build or commit therefore gets its own folder. `--run-name` names one, e.g. for `--protal-args`
  variants.
- **sylph** (`--tool sylph --sylph-db DB --gtdb DIR[,DIR]`): `sylph profile` per set, of the paired samples and of the
  single-read ones; a sample of several files is joined for it. Its run folder is `sylph_<version>_<database>`.
  - **Genomes to species:** each genome of its output gets its GTDB species by accession: first from the truth's
    release (`--gtdb` r226), then from later folders (e.g. r232 for a database of that release). A genome found in
    none of them is kept under its accession.
  - **Abundances:** each genome's `Taxonomic_abundance` (a share of cells) is summed per species into a profile in
    protal's format, without `?`.

**Scoring.** `score_mocks.py score` scores every run folder, or those given with `--runs`, into `$MOCKS/<run>/scores/`:
- `samples.tsv`: true and false positives, precision, recall and F1, the L1 distance and Bray-Curtis dissimilarity to
  the truth, and the profile's `?` against the expected share.
- `species.tsv`: each species of the truth or the profile, with both shares and its status.
- `replicates.tsv` and `groups.tsv`: Bray-Curtis of profiles and Jaccard of calls within each replicate group, e.g.
  one ONT flow cell basecalled several ways, or one library sequenced on several platforms.
- `summary.tsv`: medians per set, read type and truth kind.

All runs' rows also go into `$MOCKS/scores/summary_by_run.tsv` and `samples_by_run.tsv`, with a `run` column, so the
runs can be compared side by side.

```bash
python3 scripts/mock_benchmark/download_mocks.py -o $MOCKS --dry-run
python3 scripts/mock_benchmark/download_mocks.py -o $MOCKS --sets d6331_dna,meslier,d6331_hifi,d6331_ont --jobs 8
sbatch -c 32 --mem 64G --wrap "python3 scripts/mock_benchmark/score_mocks.py run --reads $MOCKS --db $DB -t 32 \
    --sets d6331_dna,meslier,d6331_hifi,d6331_ont"
sbatch -c 16 --mem 48G --array 0-9 --wrap "python3 scripts/mock_benchmark/score_mocks.py run --reads $MOCKS \
    --db $DB -t 16 --sets hcom2_invitro --chunks 10"
sbatch -c 32 --mem 64G --wrap "python3 scripts/mock_benchmark/score_mocks.py run --tool sylph --reads $MOCKS \
    --sylph-db gtdb-r232-c200-dbv2.syl2db --gtdb $GTDB/r226,$GTDB/r232 -t 32 --sets d6331_dna,meslier,d6331_hifi,d6331_ont"
python3 scripts/mock_benchmark/score_mocks.py score --reads $MOCKS
```

## Continuous integration

`.github/workflows/ci.yml` runs on every push and pull request (Ubuntu 24.04; `python3` with numpy,
pandas, scikit-learn and joblib from the distribution; the `zstd` CLI), with
`PROTAL_TESTS_REQUIRED=1`, so that a test whose prerequisite is missing fails instead of passing unseen
as a skip:
- a Release build with the unit tests, `just example`'s accuracy check (its mini database, built in
  `data/mini_db`, is the one the end-to-end tests use), the end-to-end tests, the script tests including
  the GTDB build end to end, and the trainer's tests;
- the unit tests of a Debug build under ASan and UBSan, with `-DPROTAL_NO_CLONES`.

## Strain test harness

The `strain-*` recipes of the [justfile](../justfile) run the strain pipeline on a simulated dataset
and build an HTML report. Set the database (`PROTAL_DB_PATH` or `strain_db=...`) and the dataset
(`strain_input=...`); the defaults are paths on the developers' cluster.

| Recipe | |
|---|---|
| `just strain-input strain_genomes=genomes.tsv` | simulate a dataset (`strain_sim_args`) and align it |
| `just strain-test` | profile with qcmsa, count each species' markers, write `strain_test_out/<variant>/report/report.html` |
| `just strain-test-raw` | the same with `--msa_min_hcov 0` (variant `test2`): all gene and sample filtering left to qcmsa |
| `just strain-refilter` | qcmsa again on a run's raw MSAs with other thresholds ([strains.md](strains.md#filtering-an-existing-run-again)) |
| `just strain-trees` | IQ-TREE trees from each MSA (`strain_tree_input=raw` for the raw MSAs) |
| `just strain-report`, `just strain-clean` | rebuild the report; remove the outputs |

## Other scripts

| Script | |
|---|---|
| `scripts/measure_performance.sh [--threads N] [--cold_reads] [--sylph SYLPH_DB] [--stage1_probe] [--test_sparse] OUT_DIR DB TYPE:R1[:R2] ...` | repeated protal runs on N threads (`--threads`, else `THREADS`, SLURM's CPUs or all; sylph's runs too): wall and CPU time, peak memory (and protal's peak after the preload, the index load, aligning and profiling), stage times, protal's read and alignment counts, the seeding's lookups and flex cells, the seeds sharing a gene, the lookups dropped as too ubiquitous and the anchors, the SAM header's genes and finishing time, and with `perf` instructions and cache misses; medians per sample. From protal's input line (since 2026-10-13) also how each read file was read (`input`, e.g. `gzip-members:22606:4t` for 22,606 gzip members read by members on 4 threads) and the share of the alignment spent getting reads (`reader_share`); a run at 0.25 or more is flagged on its progress line as held back by its input. `--cold_reads` drops the read files from the page cache before each timed run (`dd iflag=nocache`, no root needed), so that every run reads its input from the disk or network file system, as a sample's first run does. With two samples or more, then a cohort run of all of them from their SAMs (`cohort.tsv`): profiling, building the strain MSAs and qcMSA (`COHORT=0` leaves it out, `QCMSA=0` runs it without qcMSA). With `--sylph` (a sylph `.syldb` or `.syl2db`, best of the same GTDB release), last, [sylph](https://github.com/bluenote-1577/sylph) profiles each sample from the same reads as often, for comparison (`sylph.tsv`: wall and CPU time, peak memory, `perf` counters, genomes found; the profiles in `sylph/`), and the medians end with protal/sylph ratios; sylph from `$PATH` (its conda environment) or `SYLPH`, more options in `SYLPH_ARGS`. With `--stage1_probe`, after each sample's timed runs, one more run with protal's `--stage1_probe` (what the first stage of a two-stage protal would keep; not in `runs.tsv`): its tables in `stage1_probe/` (`<sample>_stage1_probe.tsv`, `_stage1_probe_taxa.tsv`, for short reads `_first_lookups_probe.tsv`), a line in `probe.tsv` (wall and CPU time, peak memory, aligning, the probe's reads and reads with records, and whether its profile is byte-identical to the last timed run's: a run whose profile differs is kept), and at the end each sample's probe line and the probe run's wall time against the timed runs' median. With `--test_sparse`, after each sample's timed runs, one run at each of four levels of protal's `--sparse_first` (C = 1, 2, 4, 8 for short reads, 4, 8, 16, 32 for long reads; `SPARSE_SHORT`, `SPARSE_LONG`, more options in `SPARSE_ARGS`; a level with an `s`, e.g. `8s`, adds `--sparse_seeding`, and `0` is an ungated run; not in `runs.tsv`): a line each in `sparse.tsv` (wall and CPU time, peak memory, aligning, `perf` counters, the reads given up and the sampled lookups, the records written against the timed runs', and the profile against the last timed run's: species called by both, lost, gained, Bray-Curtis), the profiles in `sparse/`, and at the end each level's wall and aligning time against the timed runs' medians ([report](claude/2026-10-12-sparse-first/README.md)). On a cluster, run it on a whole node |
| `scripts/db_compression_benchmark.sh` | compression ratio and speed per zstd level on a database, and load times |
| `scripts/protal_profile_utils merge`, `composition` | profiles into one abundance table (the unknown share `?` last); the samples' `<profile>.composition` files into one table ([running.md](running.md#what-the-called-species-explain-the-unknown-share)) |
| `scripts/protal_map_utils` | `generate`, `merge` (keeping the runs' SAMs), `flatten` and `validate` map files |
| `scripts/recurrent_calls.py`, `scripts/prevalence_calls.py` | a run's thin recurring calls, and calls adjusted by prevalence across samples ([running.md](running.md#where-the-outputs-go)) |
| `scripts/plot_abundances.R` | abundance plots for `simulate_metagenomes --plot_png` |

## Conventions

- Helper scripts and tests are Python 3 (argparse, a module docstring saying what the script does and
  how to call it), or plain shell for glue.
- Audits, benchmarks and reviews written with Claude Code go to [`docs/claude/`](claude/README.md).
