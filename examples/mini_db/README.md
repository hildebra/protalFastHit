# Mini database example: build a DB, profile simulated reads, check the result

A reproducible end-to-end test of protal that needs no downloads. It builds a tiny
database from a synthetic GTDB release, simulates reads from a mock community with a
known composition, runs protal on them and compares the profile with the truth.

```bash
just example                                  # builds build/protal first
# or, with any protal binary:
PROTAL=build/protal bash examples/mini_db/run.sh [--rebuild] [WORKDIR]
```

`WORKDIR` defaults to `examples/mini_db/work/` (git-ignored). `THREADS` (default 4) and
`PYTHON` (default `python3`, standard library only) can be set in the environment. The
script exits 0 if every check passes and 1 otherwise, so it can be used in CI.

## What it does

**a) Build the database** (`scripts/mini_db/build_mini_db.sh`)

1. `simulate_gtdb_release.py` (`BUILD_ARGS` in `run.sh`: `--seed 42`) writes `WORKDIR/gtdb_r226/`, a GTDB-layout release
   with 3 made-up species and 3 genomes each (one representative, two strains ~1% apart).
   *Mockella alpha* and *M. beta* share a genus (~92% identical marker genes); *Fakibacter
   gamma* is in another phylum.
2. `gtdb_to_protal_db.py` converts it into `WORKDIR/protal_db/`.
3. `protal --build` writes `index.prx` (~3 GB, whatever the reference size) and `unique_kmers.tsv`.

The DB is reused on later runs as long as the protal binary and the generator scripts are
unchanged (`WORKDIR/db.stamp`). `--rebuild` forces a new build.

**b) Test simulated reads**

1. `scripts/mini_db/simulate_reads.py --seed 7` draws 30,000 read pairs (2x150 bp, 0.2%
   substitution errors) from the genomes in [`community.tsv`](community.tsv). Abundances there
   are relative cell abundances, so each genome gets reads in proportion to abundance x length.
   Two of the three genomes are strains that are *not* the species' reference genome.
2. `protal --db WORKDIR/protal_db -1 ... -2 ... --prefix mini -o WORKDIR/run --no_qcmsa`
3. [`check_results.py`](check_results.py) compares `run/mini.profile` and `run/mini.sam` with
   `reads/mini.truth.tsv`. Per species it checks that the species is detected, that nothing
   else is, that the abundance is within 0.05, that at least 20% of the reads align (only
   marker genes, ~45% of a genome, are indexed), and that at least 95% of aligned reads hit
   a gene of their source species.

Expected output:

```
species              truth  profile  aligned_frac  assign_acc  status
s__Fakibacter gamma  0.200  0.197    0.432         1.0000      ok
s__Mockella alpha    0.500  0.502    0.449         0.9828      ok
s__Mockella beta     0.300  0.300    0.445         1.0000      ok

PASS: 3 species detected, abundances within 0.05, read assignment >= 0.95
```

The ~2% of *M. alpha* reads that land on *M. beta* come from marker regions the two
congeneric species share.

## Reproducibility

Both simulators are seeded and write byte-identical files for the same settings (also the
`.gz` files: no timestamps in the headers). `WORKDIR/checksums.md5` lists the generated DB
inputs and reads, so you can compare two runs or two machines. `WORKDIR/run_info.txt`
records the git commit, the protal binary and its md5, the Python version and the settings.
protal's own output may differ slightly with `THREADS`.

## Changing the test

- **Community:** edit `community.tsv`. Accessions are listed in
  `WORKDIR/gtdb_r226/simulation/genomes.tsv`; a species may appear with several strains.
- **Thresholds:** run `check_results.py` yourself with `--max_abundance_error`,
  `--min_assignment_accuracy` or `--min_aligned_fraction`.
- **Species or genomes:** add options for `simulate_gtdb_release.py` (e.g. `--lineages FILE`,
  one GTDB lineage per line, or `--genomes_per_species 5`) to `BUILD_ARGS` in `run.sh`. They
  are part of the DB stamp, so the next run rebuilds the DB; use `--rebuild` after editing a
  lineages file.

Profiling uses the shipped random forest (`scripts/random_forest.xml`, copied to `model.xml`).
