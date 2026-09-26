# Protal

Protal is a reference-based taxonomic profiler for bacterial communities and uses paired-end short reads from shotgun metagenomic sequencing as an input. The index is prebuilt and covers the whole taxonomic space from GTDB version r214. The index is available for download under https://protal.earlham.ac.uk/main.php?site=downloads

# Installation
Protal is in the final steps of development and is also available via conda. In the meantime, you can use a local build process via conda as described below.

## Requirements?
- git
- conda
- A linux machine (no support for mac or windows)

## Steps?
1. Install conda-build
2. clone git repository
3. build protal locally with conda-build
4. install in conda environment from local build

## 1. Install conda-build
This is needed to build a conda project from local files.
```bash
conda install conda-build
```

Alternatively, if you are using micromamba or mamba, you can also install conda-build with
```bash
micromamba install conda-build
# or
mamba install conda-build
```

## 2. Clone this repository
Clone this repository.
```bash
git clone git@github.com:4less/protal.git
```

## 3. build protal locally with conda-build
Compiles protal from the source files with instructions supplied in conda-recipe/meta.yml and conda-recipe/build.sh.
```bash
cd protal
mkdir conda-build
conda build conda-recipe -c conda-forge --output-folder conda-build

# If everything is successful, the local conda package is here
conda-build/linux-64/protal-<CURRENT_VERSION>.tar.bz2
```

## 4. Install in conda 

```bash
# Current directory is your local clone of this repository
conda create -n protal_env conda-build/linux-64/protal-<CURRENT_VERSION>.tar.bz2
#or
micromamba create -n protal_env conda-build/linux-64/protal-<CURRENT_VERSION>.tar.bz2
```

## Test the installation

```bash
conda activate protal_env
protal
```

## Compressed databases

protal reads its database files raw or compressed with [zstd](https://facebook.github.io/zstd/):
for `index.prx` and `reference.fna` it uses the raw file if it exists, else `index.prx.zst` /
`reference.fna.zst`. Compression roughly halves the bytes read (or more), which makes loading faster
wherever storage is slower than decompression, e.g. on network file systems.

protal writes these files in zstd's *seekable* format: independent frames of 64 MB each (plus a
seek table at the end), so that loading uses `-t` threads, each decompressing whole frames straight
into memory (~1-1.5 GB/s per thread). Raw files are read with `-t` threads too. A `.zst` with a
single frame (e.g. from the zstd CLI) still loads, but with one thread.

`index.prx.zst` holds the index in columns rather than as the raw bytes of `index.prx`: each
frame is a chunk of about 64 MB of the index, stored as the number of values per k-mer, then each
field of the values (taxon, gene, position, flags) in its own byte planes, only as wide as the
chunk needs; empty key blocks cost one bit. In tests this was 25-55% smaller than the raw bytes
compressed alike, and loaded about as fast (within ~10% on a synthetic index with 680 MB of values,
faster on sparse ones), since chunks decode in parallel. The writer decodes every chunk it writes and keeps the raw cells of any
that would not round-trip exactly. Because of the columns, plain `zstd -d` does not give an
`index.prx`; use `protal --decompress_db` (below). `reference.fna.zst` is plain seekable zstd.

`protal --build` writes a compressed database by default: `index.prx.zst`, and `reference.fna` is
replaced by `reference.fna.zst` (both verified; `reference.fna` is removed only afterwards).
Frames are compressed in parallel with `-t` threads.

| Build option | Default | |
|---|---|---|
| `--no_compress` | off | write `index.prx` and keep `reference.fna` raw |
| `--compress_level` | 19 | zstd level 1-22. Level 19 compresses ~3 MB/s per thread, level 12 ~40 MB/s; decompression speed barely depends on it |
| `--compress_frame_mb` | 64 | frame size; 0 writes a single frame (one loading thread) |
| `--compress_window_log` | 27 | long-distance matching window (2^27 = 128 MB, capped at the frame size); 0 turns it off |

An existing database (raw, e.g. a download, or compressed as a single frame) is converted in place,
without rebuilding it; the content stays byte-identical and each file is verified first:

```bash
protal --compress_db --db /path/to/protal-db -t 16      # -> index.prx.zst, reference.fna.zst
protal --decompress_db --db /path/to/protal-db -t 16    # back to raw, e.g. for protal versions before zstd support
```

`--decompress_db` writes `index.prx` and `reference.fna` and removes the `.zst` files; the index
is byte-identical to one built with `--no_compress`.

The order of the genes in `reference.fna` matters for its size: related species' copies of a
marker gene are similar, so a file ordered by gene (then by taxon, in taxonomic order) compresses
about 2x better than one ordered genome by genome, if related genomes lie further apart in the
file than zstd's window. protal finds genes through `reference.map`, so any order works;
`scripts/mini_db/gtdb_to_protal_db.py` writes gene order by default.

`--preload_genomes_off` (loading reference genes on demand) needs a raw `reference.fna`.
`scripts/db_compression_benchmark.sh` measures ratio and speed per level on your database and
compares protal's load times for raw and compressed copies and several thread counts.

## Metagenome simulation (C++)

Build the simulator helper binary:
```bash
cmake -S . -B cmake-build-release
cmake --build cmake-build-release --target simulate_metagenomes
```

Input TSV format (three columns): genome name, GTDB taxonomy string, path to genome FASTA (supports .gz). Example run:
```bash
./cmake-build-release/simulate_metagenomes \
  --genome_table genomes.tsv \
  --output_dir sims/ \
  --samples 3 \
  --sample_prefix sim \
  --total_read_pairs 100000 \
  --species_per_sample 15 \
  --distribution power_law \
  --strains_per_species "0.4,0.2"
```
Reads are simulated with `art_illumina`, concatenated per sample into `reads/<sample>_R1.fq.gz` and `reads/<sample>_R2.fq.gz` (compressed with `pigz`), and a `manifest.tsv` records the composition.

### Reproducing a simulated dataset

Every run writes its full provenance next to the reads:

- `manifest.tsv` — one row per (sample, genome), including the FASTA it came from and
  the `art_seed` ART used for it. `manifests/<sample>.tsv` holds the same rows split
  per sample.
- `run_params.tsv` — the command line, the RNG seed (resolved and recorded even when
  `--seed` was not given), and the ART settings.

A manifest fully describes a dataset's composition (genomes, read pairs, ART seeds), so
replaying one does not depend on reproducing the community-design RNG:

```bash
./cmake-build-release/simulate_metagenomes \
  --from_manifest sims/manifest.tsv \
  --output_dir sims_replay/
```

The replay reproduces the reads byte for byte, provided it uses the same `art_illumina`
build and the same ART settings. The manifest stores the ART seeds but not the settings,
so `--read_length`, `--fragment_mean`, `--fragment_stdev`, `--sequencer` and
`--extra_art_args` still come from the command line; the replay warns about each one that
differs from the original run's `run_params.tsv`, and also checks `--read_length` against
the depths the manifest implies. All sampling options (`--distribution`,
`--species_per_sample`, `--seed`, …) are ignored. A per-sample manifest replays just that
sample.

Manifests written before the `fasta_path` and `art_seed` columns existed still replay:
pass `--genome_table` so the genomes can be resolved by name, and the composition and
per-genome depth are reproduced exactly while the reads themselves are fresh
realizations. The replay's own manifest carries seeds, so it is exactly reproducible
from then on.

## Mini database for local testing

`scripts/mini_db/` builds a small protal database from a synthetic, sparse GTDB
release (3 species, 3 genomes each by default), so you can test without the full DB:

```bash
just mini-db          # or: PROTAL=build/protal bash scripts/mini_db/build_mini_db.sh data/mini_db
protal --db data/mini_db/protal_db -1 r1.fq -2 r2.fq -o out/
```

- `simulate_gtdb_release.py` writes a directory laid out like an extracted GTDB release
  (taxonomy, metadata, `*_marker_genes_{reps,all}` per-marker FASTAs, representative
  genomes) plus `simulation/genomes.tsv`, which you can pass to `simulate_metagenomes --genome_table`.
- `gtdb_to_protal_db.py` turns such a release (synthetic or a real, extracted one) into
  `reference.fna`, `reference.map`, `internal_taxonomy.dmp`, `full_reference.fna` and
  `model.xml`. `build_mini_db.sh` then runs `protal --build` on them.
- `simulate_reads.py` draws paired reads from a mock community of those genomes (no ART needed)
  and writes the truth table next to them.
- The build compresses the database (see [Compressed databases](#compressed-databases)):
  `index.prx.zst` is about 0.5 MB here, while a raw `index.prx` is about 3 GB even for a tiny
  reference, because the k-mer key map has a fixed size. `PROTAL_BUILD_ARGS=--no_compress`
  builds a raw one.

## Testing

```bash
just test    # C++ unit tests (GoogleTest, needs libgtest-dev): cmake -DPROTAL_BUILD_TESTS=ON, then ctest
just e2e     # builds the mini database, then runs tests/e2e/test_protal_e2e.py against it
just example # mini database + reads from a known mock community: is the profile right?
```

`just example` runs [`examples/mini_db/run.sh`](examples/mini_db/README.md): a seeded, fully
reproducible build of the mini database, 30,000 simulated read pairs from three genomes (two of
them strains that differ from the reference), a protal run, and a check that every species is
detected with the right abundance (±0.05) and that aligned reads hit their own species.

The end-to-end tests simulate reads from the database's reference genes and run the real
`protal` and `simulate_metagenomes` binaries: exit codes, output files, SAM records, strain
MSAs, reruns and failure reporting. Point them at any database with
`PROTAL_TEST_DB=<db> python3 -m unittest -v tests/e2e/test_protal_e2e.py` (`PROTAL` and
`SIMULATE` select the binaries). CI (`.github/workflows/ci.yml`) runs both, plus the unit tests
of a Debug build under AddressSanitizer and UndefinedBehaviorSanitizer.
