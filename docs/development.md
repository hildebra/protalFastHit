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
binary, `protal_tests`: about 400 tests, 10 seconds on 4 cores. Most compare protal's code with an
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
x86-64 (`src/Utilities/TargetClones.h`), so that this run tests the copies a CPU without AVX2 runs;
the Release build on a machine with AVX2 runs the others:

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
and take about 1.5 minutes on 4 cores, 3.5 GB of memory and 3-4 GB of `/tmp`:

```bash
PROTAL_TEST_DB=data/mini_db/protal_db PROTAL=build/protal SIMULATE=build/simulate_metagenomes \
    PROTAL_TESTS_REQUIRED=1 python3 -m unittest -v tests/e2e/test_protal_e2e.py
```

A test whose prerequisite is missing (the database, a binary, the `zstd` CLI) is
skipped; `PROTAL_TESTS_REQUIRED=1` makes it an error instead, as CI runs them. The tests that need no
database (the version, the simulator, qcmsa's contract, small builds, gene neighbours) run without
`PROTAL_TEST_DB`.

## Script tests

`scripts/mini_db/test_*.py` test the mini database scripts (`test_mini_db.py`), the GTDB downloads
(`test_downloads.py`), the training data collector and its scenarios (`test_collector.py`), the gene
neighbours (`test_gene_neighbours.py`), the build script's parts (`test_gtdb_build.py`) and the GTDB
build end to end (`test_gtdb_pipeline.py`); `scripts/test_*.py` the trainer and the model's PMML
export (`test_model_pmml.py`), the in-silico strains, `trace_relatives.py`, `error_reads.py`,
`composition_accuracy.py`, the profile scripts and the strain test scripts. A test whose prerequisite is missing is skipped;
`PROTAL_TESTS_REQUIRED=1` makes it fail instead (`scripts/prerequisites.py`). With everything present
they run in about 4 minutes on 4 cores, the GTDB build 3 of them:

```bash
PROTAL=$PWD/build/protal SIMULATE=$PWD/build/simulate_metagenomes PROTAL_TESTS_REQUIRED=1 \
    python3 -m unittest scripts/mini_db/test_*.py
python3 -m unittest scripts/test_insilico_strains.py scripts/test_trace_relatives.py scripts/test_error_reads.py \
    scripts/test_composition_accuracy.py scripts/test_profile_scripts.py scripts/test_strain_scripts.py \
    scripts/test_model_pmml.py
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
- the full build's samples profiled as they are simulated (`--profile-blocks` of a few kB): the same
  tables.

`build_gtdb_releases.py` is tested with a stand-in of the build script (`test_gtdb_build.py`).
The trainer's tests fit on one thread (`OMP_NUM_THREADS=1`): on every core of a busy machine,
gradient boosting made `test_model_pmml.py` take 20 minutes (scikit-learn 1.9.1).

## Mini databases

`scripts/mini_db/` builds a small database from a synthetic GTDB release (3 species, 3 genomes each
by default), to test without the full database. It is about 1 MB as `database.protal`
(`PROTAL_BUILD_ARGS=--no_bundle` or `--no_compress` for separate files; a raw `index.prx` is 3 GB
whatever the reference, since the k-mer key map has a fixed size):

```bash
just mini-db          # or: PROTAL=build/protal bash scripts/mini_db/build_mini_db.sh data/mini_db
protal --db data/mini_db/protal_db -1 r1.fq -2 r2.fq -o out/
```

`build_mini_db.sh` runs the same steps as a GTDB database ([databases.md](databases.md#step-by-step)):
the release, the converter, the gene neighbours, `protal --build`. The scripts:

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
(`OUTDIR/genomes.tsv`), and synthetic releases have `simulation/genomes.tsv`. A first row is taken
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

## Continuous integration

`.github/workflows/ci.yml` runs on every push and pull request (Ubuntu 24.04; `python3` with numpy,
pandas, scikit-learn and joblib from the distribution; the `zstd` CLI), with
`PROTAL_TESTS_REQUIRED=1`, so that a test whose prerequisite is missing fails instead of passing unseen
as a skip:
- a Release build with the unit tests, the mini database, the end-to-end tests, `just example`'s
  accuracy check, the script tests including the GTDB build end to end, and the trainer's tests;
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
| `scripts/measure_performance.sh OUT_DIR DB TYPE:R1[:R2] ...` | repeated protal runs: wall and CPU time, peak memory (and protal's peak after the preload, the index load, aligning and profiling), stage times, protal's read and alignment counts, the seeding's lookups and flex cells, the seeds sharing a gene, the lookups dropped as too ubiquitous and the anchors, the SAM header's genes and finishing time, and with `perf` instructions and cache misses; medians per sample. With two samples or more, then a cohort run of all of them from their SAMs (`cohort.tsv`): profiling, building the strain MSAs and qcMSA (`COHORT=0` leaves it out, `QCMSA=0` runs it without qcMSA). On a cluster, run it on a whole node |
| `scripts/db_compression_benchmark.sh` | compression ratio and speed per zstd level on a database, and load times |
| `scripts/protal_profile_utils merge`, `composition` | profiles into one abundance table (the unknown share `?` last); the samples' `<profile>.composition` files into one table ([running.md](running.md#what-the-called-species-explain-the-unknown-share)) |
| `scripts/protal_map_utils` | `generate`, `merge` (keeping the runs' SAMs), `flatten` and `validate` map files |
| `scripts/recurrent_calls.py`, `scripts/prevalence_calls.py` | a run's thin recurring calls, and calls adjusted by prevalence across samples ([running.md](running.md#where-the-outputs-go)) |
| `scripts/plot_abundances.R` | abundance plots for `simulate_metagenomes --plot_png` |

## Conventions

- Helper scripts and tests are Python 3 (argparse, a module docstring saying what the script does and
  how to call it), or plain shell for glue.
- Audits, benchmarks and reviews written with Claude Code go to [`docs/claude/`](claude/README.md).
