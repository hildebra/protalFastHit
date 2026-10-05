# Development: tests, mini databases and simulated metagenomes

```bash
just test          # C++ unit tests (GoogleTest)
just e2e           # builds the mini database, then the end-to-end tests against it
just example       # mini database + reads from a known mock community: is the profile right?
just mini-db-test  # Python tests: the mini database scripts, the GTDB build pipeline, in-silico strains, tracing
just model-test    # the presence model's PMML export scores as scikit-learn does
```

The unit tests need GoogleTest (`libgtest-dev` on Ubuntu). The Python tests need numpy
(`python3-numpy`), and `model-test` also pandas and scikit-learn, as training does. Build
requirements are in [installation.md](installation.md#building-from-source).

## Unit tests

`tests/test_*.cpp` cover SNP calling, the SAM round trip, index building and lookup, the index
column codec, zstd and the single-file database, the binary gene table, input validation, parsing,
strain output, the alignment screen, exact sums, and that a SAM's profile is the same on any number
of threads and in any order of its reads (`test_ProfileThreads.cpp`). They build as one binary:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON
cmake --build build --target protal_tests
ctest --test-dir build --output-on-failure
```

Under AddressSanitizer and UndefinedBehaviorSanitizer, as CI runs them (a Debug build, so `assert()`
is on):

```bash
FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer"
cmake -S . -B build-asan -G Ninja -DCMAKE_BUILD_TYPE=Debug -DPROTAL_BUILD_TESTS=ON \
  -DCMAKE_CXX_FLAGS="$FLAGS" -DCMAKE_C_FLAGS="$FLAGS" -DCMAKE_EXE_LINKER_FLAGS="$FLAGS"
cmake --build build-asan --target protal_tests
ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1 ctest --test-dir build-asan --output-on-failure
```

WFA2-lib is built without UBSan (`lib/wfa2-lib.cmake`): its unaligned loads and shifts of negative
offsets would stop every test that aligns. zlib-ng is built with both sanitizers.

## End-to-end tests

`tests/e2e/test_protal_e2e.py` simulates reads from a database's reference genes and runs the real
`protal` and `simulate_metagenomes`: exit codes, output files, SAM records, strain MSAs, reruns and
failure reporting, for paired-end and single-end reads. Point it at any database (a
`database.protal`, its folder, or separate files):

```bash
PROTAL_TEST_DB=data/mini_db/protal_db PROTAL=build/protal SIMULATE=build/simulate_metagenomes \
    python3 -m unittest -v tests/e2e/test_protal_e2e.py
```

The tests need the `zstd` CLI (or Python 3.14) to read `.sam.zst` outputs.

**The GTDB build pipeline.** `scripts/mini_db/test_mini_db.py` also runs `build_gtdb_database.py` end
to end (`GtdbBuildTest`, a few minutes) when `$PROTAL` and `$SIMULATE` name the binaries,
`art_illumina` is on `$PATH`, and `$PROTAL_TRAIN_PYTHON` (default: the Python running the tests) has
scikit-learn. On a synthetic release of 60 species, downloaded from a fake GTDB server and a fake
NCBI, it covers:
- a build trained for pe and se, and a rerun that skips the conversion and both builds;
- another seed that rebuilds only the training database;
- a background build that fails and stops the run;
- `SIGTERM`, after which no command is left running;
- a reduced database (`--n-genes`, then `--genes`) with a gene archaea have;
- `build_gtdb_releases.py` with its summary;
- the default scenarios, made small by a `--scenario-file` of their names (one with 90% host reads
  from a random host genome, soil scaled down to the genome table), with Illumina reads at a target
  quality and Ultima reads, their hold-in and hold-out samples scored in every report and in
  `summary.txt`, the feature sets the trainers chose and why, and the host scenario left out without
  a host genome. The other build tests pass `--scenarios none`: the presets have real depths.

```bash
PROTAL=build/protal SIMULATE=build/simulate_metagenomes PROTAL_TRAIN_PYTHON=~/protal-train/bin/python \
    python3 -m unittest -v scripts.mini_db.test_mini_db.GtdbBuildTest
```

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
| `simulate_reads.py` | read pairs drawn straight from a mock community's genomes (substitutions only, no ART), with the truth |

`simulate_gtdb_release.py` evolves each marker from one random coding sequence down the lineage.
Its options:
- `--seed`, `--lineages FILE` (one GTDB lineage per line), `--genomes_per_species`,
  `--genome_length`, `--contigs`, `--marker_loss`;
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
Illumina reads for them with [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art)
(`art_illumina`). It makes the tests' and benchmarks' samples and the model's training data
([databases.md](databases.md#training-data)). Build it with the other binaries
(`cmake --build build --target simulate_metagenomes`, or `just simulate`).

**Input**: a tab-separated genome table without header: name (accession), GTDB taxonomy, FASTA path
(`.gz` works), and optionally the genome's length, which spares reading every genome at the start
(~20 ms per genome of GTDB's size). `build_gtdb_database.py` writes one with lengths
(`OUTDIR/genomes.tsv`), and synthetic releases have `simulation/genomes.tsv`.

```bash
./build/simulate_metagenomes --genome_table genomes.tsv --output_dir sims/ --samples 3 \
  --sample_prefix sim --total_read_pairs 100000 --species_per_sample 15 --distribution power_law \
  --strains_per_species "0.4,0.2" --seed 1 -t 8 --protal_metafile sims/protal
protal --db DB --map sims/protal.meta -t 8      # prints true and false positives per sample
```

| Output | |
|---|---|
| `reads/<sample>_R1.fq.gz`, `_R2.fq.gz` | the reads (BGZF, compressed as each genome's reads are appended) |
| `manifest.tsv`, `manifests/<sample>.tsv` | per (sample, genome): read pairs, relative abundance, vertical coverage, the FASTA and ART's seed |
| `abundance_matrix.tsv` | relative abundance of each species in each sample |
| `run_params.tsv` | the command line, the seed (also when not given) and the ART settings |
| `protal.meta`, `protal_goldstd/<sample>.profile_truth` | with `--protal_metafile DIR`: a protal map with a `PROFILE_TRUTH` column; protal then writes `<profile>.truth_annotated` |
| `plots/<sample>.png` | with `--plot_png` (needs `Rscript`) |

| Option | Default | |
|---|---|---|
| `-n, --samples`, `--sample_prefix` | 1, `sample` | samples, named `<prefix>_<n>` |
| `--total_read_pairs` | 100000 | read pairs per sample |
| `--species_per_sample` | 10 | a number or a range, e.g. `20-80` |
| `--distribution` | `poisson_lognormal` | `power_law` (`--alpha`), `negative_binomial` (`--nb_r`, `--nb_p`) or `poisson_lognormal` (`--pln_mu`, `--pln_sigma`; several sigmas go to the samples in turn) |
| `--strains_per_species` | none | probabilities of a 2nd, 3rd, ... strain of a species, e.g. `0.4,0.2,0.1` |
| `--include_species`, `--genus`, `--taxon` | | species in every sample; `g__A:10,g__B:2` species from those genera, `d__Archaea:10` from any taxon (`--pick_random_demand_if_fail` caps instead of failing) |
| `--congener_groups` | none | `SHARE:MIN-MAX`, e.g. `0.25:2-5`: about SHARE of each sample's species in groups of MIN to MAX congeners |
| `--strain_sharing_file` | | strains shared across samples (below) |
| `--seed` | random | |
| `--read_length`, `--fragment_mean`, `--fragment_stdev`, `--sequencer` | 150, 350, 50, `HS25` | ART's read, fragment and error profile; `--extra_art_args` passes more |
| `-t, --threads` | 1 | samples at a time; more threads run ART on a sample's genomes side by side. The samples are the same for any number |
| `--test` | off | the design, manifests and truth, no reads |
| `--keep_tmp`, `--art_path`, `-v` | | keep each genome's reads; ART's path; the version and commit |

**Strains shared across samples.** `--strain_sharing_file` takes a tab-separated file, one species
per row: `SPECIES`, `SAMPLE_FRACTION` (share of samples with it), `N_STRAINS` (distinct strains
across samples), `MIN_OCCURRENCE` (samples each strain is in at least), and optionally `MIN_VCOV`
(minimum coverage per strain and sample) and `CONSPECIFIC_STRAINS` (as `--strains_per_species`).

**Replaying a dataset.** A manifest holds the genomes, read pairs and ART seeds, so
`--from_manifest sims/manifest.tsv --output_dir replay/` reproduces the reads byte for byte with the
same `art_illumina` and ART settings. The settings come from the command line; the replay warns about
each that differs from the original `run_params.tsv`. Sampling options are ignored, and a per-sample
manifest replays that sample. Manifests without seeds replay the composition with fresh reads (pass
`--genome_table`).

## Continuous integration

`.github/workflows/ci.yml` runs on every push and pull request (Ubuntu 24.04, `python3` with
`python3-numpy`): a Release build with the unit tests, the mini database tests and the in-silico
strain and trace tests; and the unit tests of a Debug build under ASan and UBSan. The end-to-end
tests, `GtdbBuildTest` and `just example` are not part of CI; run them before a merge.

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
| `scripts/measure_performance.sh OUT_DIR DB TYPE:R1[:R2] ...` | repeated protal runs: wall and CPU time, peak memory, stage times, protal's read and alignment counts, the seeding's lookups and flex cells, the SAM header's genes and finishing time, and with `perf` instructions and cache misses; medians per sample. With two samples or more, then a cohort run of all of them from their SAMs (`cohort.tsv`): profiling, building the strain MSAs and qcMSA (`COHORT=0` leaves it out, `QCMSA=0` runs it without qcMSA). On a cluster, run it on a whole node |
| `scripts/db_compression_benchmark.sh` | compression ratio and speed per zstd level on a database, and load times |
| `scripts/protal_profile_utils merge` | profiles into one abundance table |
| `scripts/protal_map_utils` | `generate`, `merge`, `flatten` and `validate` map files |
| `scripts/recurrent_calls.py`, `scripts/prevalence_calls.py` | a run's thin recurring calls, and calls adjusted by prevalence across samples ([running.md](running.md#where-the-outputs-go)) |
| `scripts/plot_abundances.R` | abundance plots for `simulate_metagenomes --plot_png` |

## Conventions

- Helper scripts and tests are Python 3 (argparse, a module docstring saying what the script does and
  how to call it), or plain shell for glue.
- Audits, benchmarks and reviews written with Claude Code go to [`docs/claude/`](claude/README.md).
