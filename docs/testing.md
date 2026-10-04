# Testing and development

```bash
just test        # C++ unit tests (GoogleTest): cmake -DPROTAL_BUILD_TESTS=ON, then ctest
just e2e         # builds the mini database, then runs tests/e2e/test_protal_e2e.py against it
just example     # mini database + reads from a known mock community: is the profile right?
just mini-db-test  # unit tests of the mini database generator and converter (no protal binary needed)
just model-test    # the presence model's PMML export scores as scikit-learn does (needs scikit-learn)
```

The unit tests need GoogleTest (`libgtest-dev` on Ubuntu); the Python tests use the standard
library only, except `model-test`, which needs numpy, pandas and scikit-learn, as training does.

`mini-db-test` also runs `build_gtdb_database.py` end to end (`GtdbBuildTest`, about 2 minutes) when
`$PROTAL` and `$SIMULATE` name the binaries and `art_illumina` is on `$PATH`, with a Python that has
scikit-learn (`$PROTAL_TRAIN_PYTHON`, default the one running the tests): a synthetic release of
60 species, downloaded from a fake GTDB server and a fake NCBI `datasets`, built and trained for
pe and se; a rerun that skips the conversion and both builds; another seed that rebuilds only the
training database and collects again; a build failing in the background that stops the run at
once; `SIGTERM`, after which no command of the run is left; a reduced database of the three
most distinctive genes (`--n-genes`, then `--genes`), one of them a gene archaea have; and
`build_gtdb_releases.py` building the full database (its genes ranked with `--rank-genes`) and
the reduced one from that ranking, with their summary. `GeneSubsetTest` covers the converter's
`--genes`, the neighbours counted over the genes kept and, with `$PROTAL`, the build of such a
folder, `rank_genes.py` on a full build (the domains' columns, a gene reserved for archaea) and
`--build_gene_subset`; `MiniDbTest.test_download_releases` the download phase
([building-a-database.md](building-a-database.md#reduced-marker-sets)).

```bash
PROTAL=build/protal SIMULATE=build/simulate_metagenomes PROTAL_TRAIN_PYTHON=~/protal-train/bin/python \
    python3 -m unittest -v scripts.mini_db.test_mini_db.GtdbBuildTest
``` Build requirements are in [installation.md](installation.md#building-from-source).

## Unit tests

`tests/test_*.cpp` cover SNP calling, the SAM round trip, index building and lookup, the index
column codec, zstd and the single-file database, input validation, parsing, and strain output.
They build as one binary:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON
cmake --build build --target protal_tests
ctest --test-dir build --output-on-failure
```

Under AddressSanitizer and UndefinedBehaviorSanitizer, as CI runs them (a Debug build, so
`assert()` is on):

```bash
FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer"
cmake -S . -B build-asan -G Ninja -DCMAKE_BUILD_TYPE=Debug -DPROTAL_BUILD_TESTS=ON \
  -DCMAKE_CXX_FLAGS="$FLAGS" -DCMAKE_C_FLAGS="$FLAGS" -DCMAKE_EXE_LINKER_FLAGS="$FLAGS"
cmake --build build-asan --target protal_tests
ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1 ctest --test-dir build-asan --output-on-failure
```

WFA2-lib is built without UBSan (`-fno-sanitize=undefined` in `lib/wfa2-lib.cmake`; ASan stays
on): its unaligned 8-byte loads and left shifts of negative offsets would otherwise stop every
test that aligns. zlib-ng (`lib/zlib-ng`) is built with both sanitizers.

## End-to-end tests

`tests/e2e/test_protal_e2e.py` simulates reads from a database's reference genes (sequencing
errors only, fixed seed) and runs the real `protal` and `simulate_metagenomes` binaries: exit
codes, output files, SAM records, strain MSAs, reruns and failure reporting, for paired-end and
single-end reads (`SingleEndTest` stands in the database's own model for `model_se.xml`). Point
them at any database:

```bash
PROTAL_TEST_DB=data/mini_db/protal_db PROTAL=build/protal SIMULATE=build/simulate_metagenomes \
    python3 -m unittest -v tests/e2e/test_protal_e2e.py
```

`PROTAL_TEST_DB` may be a `database.protal`, its folder, or separate raw or compressed files.
The tests need the `zstd` CLI (or Python 3.14) to read protal's default `.sam.zst` outputs, and
for compressed database files. `just e2e` builds the mini database first.

## Mini database

`scripts/mini_db/` builds a small protal database from a synthetic, sparse GTDB release (3
species, 3 genomes each by default), so you can test without the full database:

```bash
just mini-db          # or: PROTAL=build/protal bash scripts/mini_db/build_mini_db.sh data/mini_db
protal --db data/mini_db/protal_db -1 r1.fq -2 r2.fq -o out/
```

- `simulate_gtdb_release.py` writes a directory laid out like an extracted GTDB release
  (taxonomy, metadata, `*_marker_genes_{reps,all}` per-marker FASTAs, representative genomes)
  plus `simulation/genomes.tsv`, which you can pass to `simulate_metagenomes --genome_table`.
  Options: `--seed`, `--lineages FILE` (one GTDB lineage per line), `--genomes_per_species`,
  `--genome_length`, `--contigs`, `--marker_loss`, `--strain_divergence`, `--species_divergence`
  (these two take a rate or a range `LOW-HIGH` drawn per genome or species, written to
  `simulation/divergence.tsv`), `--gene_rates categories` (markers evolve at different speeds by what
  they do: ribosomal proteins 0.4, translation and transcription 0.8, tRNA synthetases and modification
  1.1, the rest 1.4, with noise, mean 1; `simulation/gene_rates.tsv`), `--operons` (markers in
  clusters of up to 6 genes 0-150 bases apart, in the same order in every species, each family
  breaking some up with `--operon_breaks`, 0.25: read pairs and long reads then span neighbouring
  genes, as in real genomes; without it each species' markers are shuffled and about 1.2 kb apart).
  `simulation/marker_positions.tsv` says where each gene lies.
- `gene_neighbours.py` finds the database's genes in every genome of their species (the
  representative's by their sequence, other strains' by their k-mer trace) and writes
  `gene_neighbours.tsv` (how often which genes lie next to which, per clade) and
  `gene_positions.tsv` (where each gene lies in each genome;
  [building-a-database.md](building-a-database.md#gene-neighbours)); `build_mini_db.sh` runs it.
- `gtdb_like_lineages.py` writes lineages shaped like GTDB's for `--lineages`: up to 999 species,
  most genera with one species and a few with many, unique names, a share of archaea. It makes the
  world the presence model's training was tuned on ([model-training.md](model-training.md)).
- `gtdb_to_protal_db.py` turns such a release (synthetic or a real, extracted one) into the input
  files of a database (`--exclude_species` leaves species out of a training database, keeping their
  taxids), and `build_mini_db.sh` runs `protal --build` on them
  ([building-a-database.md](building-a-database.md)).
- `simulate_reads.py` draws paired reads from a mock community of those genomes (no ART needed)
  and writes the truth table next to them.
- The build packs the database into `protal_db/database.protal` (see
  [database-files.md](database-files.md)), about 1 MB here, while a raw `index.prx` is about 3 GB
  even for a tiny reference, because the k-mer key map has a fixed size.
  `PROTAL_BUILD_ARGS=--no_bundle` builds separate compressed files, `PROTAL_BUILD_ARGS=--no_compress`
  a raw database.

`just example` runs [`examples/mini_db/run.sh`](../examples/mini_db/README.md): a seeded, fully
reproducible build of the mini database, 30,000 simulated read pairs from three genomes (two of
them strains that differ from the reference), a protal run, and a check that every species is
detected with the right abundance (within 0.05) and that aligned reads hit their own species.

## Continuous integration

`.github/workflows/ci.yml` runs on every push and pull request, on Ubuntu 24.04:

- a Release build of `protal`, `simulate_metagenomes` and the unit tests, the unit tests, and the
  mini database generator tests;
- the unit tests of a Debug build under AddressSanitizer and UndefinedBehaviorSanitizer.

The end-to-end tests and `just example` are not part of CI; run them before a merge.

## Strain test harness

The `strain-*` recipes of the [justfile](../justfile) run the strain pipeline on a simulated
dataset and build a self-contained HTML QC report. Set the database (`PROTAL_DB_PATH` or
`strain_db=...`) and the dataset (`strain_input=...`) on the command line; the defaults are paths
on the developers' cluster.

| Recipe | |
|---|---|
| `just strain-input strain_genomes=genomes.tsv` | simulate a dataset (`strain_sim_args`) and align it, in the layout the other recipes read |
| `just strain-test` | profile the existing alignments with the qcmsa post-filter, count each species' markers in the database, write `strain_test_out/<variant>/report/report.html` |
| `just strain-test-raw` | the same with `--msa_min_hcov 0` into variant `test2`: only the SNP filters, all gene and sample filtering left to qcmsa |
| `just strain-refilter` | re-run qcmsa on a run's raw MSAs with other thresholds ([qcmsa.md](qcmsa.md)) |
| `just strain-trees` | IQ-TREE trees from each MSA (`strain_tree_input=raw` for the raw MSAs) |
| `just strain-report`, `just strain-clean` | rebuild the report, remove the outputs |

The report comes from `scripts/strain_test/strain_report.py` (standard library only), the marker
counts from `scripts/strain_test/db_gene_counts.py`.

## Other scripts

| Script | |
|---|---|
| `scripts/db_compression_benchmark.sh` | compression ratio and speed per zstd level on a database, and protal's load times ([database-files.md](database-files.md#measuring)) |
| `scripts/measure_performance.sh OUT_DIR DB TYPE:R1[:R2] ...` | repeated protal runs on read files: wall and CPU time, peak memory, the index load, alignment, profiling and strain times, protal's counts of reads, anchors and candidate alignments (tried, refused by the k-mer screen, aligned, made, written), each run's stage timers (`misc/<prefix>_runtime.tsv`: the alignment stage's per thread, the profiling steps' wall times) and, where `perf` works, instructions, cycles and cache misses; medians per sample. For comparing machines, builds and databases; on a cluster, run it on a whole node (`sbatch --exclusive`) |
| `scripts/protal_profile_utils merge` | merge profiles into one abundance table |
| `scripts/protal_map_utils` | `generate`, `merge`, `flatten` and `validate` map files |
| `scripts/plot_abundances.R` | abundance bar plots for `simulate_metagenomes --plot_png` |

## Conventions

- Helper scripts and tests in this repository are Python 3 (argparse, a module docstring saying
  what the script does and how to call it), or plain shell for glue.
- Audits, benchmarks and reviews written with Claude Code go to [`docs/claude/`](claude/README.md).
