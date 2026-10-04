# Reduced marker sets: the gene subset audited, and chosen by rank in the build pipeline

Date: 2026-10-04. Code: `281a4ba` (0.7.5) plus this work. Data: the synthetic releases of
`scripts/mini_db/` (no GTDB-scale run yet, see "Open"). Commands: the tests at the end.

## What there was

`protal --build --build_gene_subset genes.txt` (commit `18470cf`, February 2026) takes a file of
gene ids and builds the index from those genes only; `scripts/subset_genes.py` (`7f38562`) cut a
`reference.fna` down to a gene list. The documentation named both and said that a model trained
on the full marker set sees fewer genes per species with a subset. Nothing chose the genes: the
reduced database on the downloads page was made by hand. The question asked was to audit this
code, bring it up to what 0.7.5 does at build time, add a `--n-genes N` option and a named-gene
option to `build_gtdb_database.py`, and document the whole.

## Audit

Where the subset reached, and where it did not, at `281a4ba`:

1. **Only the two index passes applied the subset** (`Build.h`, pass 1 and pass 2). The uniqueness
   check still read every record of the full reference (every genome's copy of every gene) and
   looked its k-mers up; with 12 of 120 genes that is ~90% of the longest build phase for nothing,
   as those k-mers are not in the index. The gene conservation estimate, the congeners report and
   the suspect-copy scan ran over every gene of `reference.map`, so `gene_conservation.tsv` and
   `suspect_copies.tsv` carried rows for genes nothing could hit. Fixed: the full-reference loop
   skips records of genes outside the subset, and the three later phases take only the subset's
   genes (`options.BuildGeneAllowed`). The unique k-mer table was right already: it writes a row
   per gene with values in the index, and since the hittable-gene logic of 0.7 (`GenomeLoader::
   LoadUniqueKmers`) a run counts only those genes, so the per-species denominators of a reduced
   database were correct.
2. **The gene neighbours were wrong for a subset.** `gene_neighbours.tsv` (0.7.1 onwards) is
   counted over every gene of every genome. Build a subset on top of it and a gene whose true
   neighbour was left out keeps that neighbour in the table; at run time the gene a long read
   meets next is the nearest gene of the subset, which the table calls an unlikely partner (the
   foreign-gene rule), and the long-read chain looks for partner genes that have no k-mers in the
   index. Fixed two ways: `protal --build` refuses `--build_gene_subset` on a folder that has the
   table (exit 8, with the way out in the message), and the converter's new `--genes` counts the
   table anew from `gene_positions.tsv` over the genes kept (`gene_neighbours.derive(..., genes=)`),
   so that in the derived folder a gene end faces the nearest kept gene, with the larger gap.
3. **The other genes stayed in the database.** `reference.fna` was indexed by the subset but packed
   whole into `database.protal` and loaded whole at run time (the 2-bit gene store), so the
   reference's share of memory and of the file did not shrink. The folder route (below) leaves
   them out.
4. **`subset_genes.py` did a quarter of the job.** It wrote a smaller `reference.fna` only: no
   `reference.map` (a second script), no full reference (so the uniqueness check would see copies
   of genes the database lacks, or none at all), no `gene2geneid.tsv`, no neighbours; its list
   format (`taxid_geneid`) was not the one `--build_gene_subset` reads; unknown ids were ignored
   without a word; a `.zst` reference could not be read. Removed; `gtdb_to_protal_db.py --genes`
   replaces it.
5. **Nothing ranked the genes.** The downloads-page subset's choice is not recorded anywhere in
   the repository.
6. Left as it was, on purpose: `species_priors.tsv`'s `markers` column counts the markers found in
   the representative over GTDB's whole set, also in a reduced database (the priors describe the
   genome's completeness, not the database; the prior features are opt-in since `281a4ba`).

## What changed

- `src/Build.h`: the subset applied to the uniqueness check, the conservation factors, the
  congeners and the suspect copies; `--build_gene_subset` refused with `gene_neighbours.tsv` in
  the folder. `src/Options.h`: the option's help says what it does and does not do.
- `scripts/mini_db/gtdb_to_protal_db.py`: `--genes LIST` with `--gtdb` and with `--from_db` (GTDB
  marker ids, with or without the Pfam version, or protal gene ids; comma-separated or a file's
  first column). The genes keep their ids; `reference.fna`, `reference.map`, the full reference and
  `gene2geneid.tsv` hold them alone; `--from_db` derives the neighbours over them. `exclude_from_db`
  became `derive_db(src, dst, names, genes)`; `--from_db` without either option is a plain copy.
  `read_gene_ids` and `read_gene_list` are shared with the pipeline.
- `scripts/mini_db/gene_neighbours.py`: `read_positions` and `derive` take a gene set.
- `scripts/rank_genes.py` (new): the ranking table and the gene list, from a build's separate files.
- `scripts/build_gtdb_database.py`: `--n-genes N`, `--genes LIST`, `--gene-ranking FILE`; the
  release converted whole into `OUTDIR/.converted` when a subset is asked for, a ranking build of
  the training database (every gene, species held out; `--no_bundle`, the training level, removed
  once ranked), `gene_ranking.tsv` and `gene_subset.txt` in OUTDIR and `model_logs/`, both database
  folders derived with `--genes`, the subset in the resume keys of both builds and of the ranking,
  `marker_genes` in `build_metadata.tsv`, a step of its own on the console.
- `docs/building-a-database.md` (reduced marker sets, the converter's option table),
  `docs/testing.md`; `scripts/subset_genes.py` removed.

## The score

A gene's score is prevalence times unique share. Prevalence is the share of the database's species
whose copy of the gene has k-mers in the index (one row of `unique_kmers.tsv` each): a gene that
most species lack cannot detect them. Unique share is the share of a copy's k-mers that are unique
to their species (short or long unique in `unique_kmers.tsv`), averaged over the species that have
the gene: it is what a read of that gene contributes to naming its species, and it is low for a
gene conserved between species or shared with another genus's copies. Both are what the build
already measures, at no extra cost, and the product is between 0 and 1.

Not in the score, but in the table so that a subset can be chosen otherwise: the copies' mean
length (a longer gene catches more reads per genome, at proportionally more index), the congeners'
between-species factor and identical share (`gene_congeners.tsv`: how a gene differs between
congeneric species against within species) and the number of suspect copies. On the test
release the three best of 168 genes scored 0.83 down to 0.77.

## The pipeline with `--n-genes`

1. The release is converted whole into `OUTDIR/.converted` (the full folder is what the two
   database folders, and the ranking build, are derived from; it is removed at the end, and a
   rerun that needs it converts again, as the non-subset flow does for a new seed).
2. The species to hold out are chosen as before.
3. A ranking copy (every gene, those species left out) is derived and built with `--no_bundle`
   at the training database's compression level, ranked by `rank_genes.rank`, and removed. The
   ranking is kept (`gene_ranking.tsv`; resume key: conversion, held-out species, protal, the
   ranker). `--gene-ranking FILE` takes a ranking of an earlier full build of the release
   instead; `--genes LIST` needs no ranking.
4. The N best go to `gene_subset.txt`; `protal_db` and the training database are derived with
   `--genes gene_subset.txt` (`protal_db_files.log`, `training_db.log`), built, collected on and
   trained as before.

The ranking build costs what the training database's build costs on the full gene set (about
15 min at r226 with 16 threads, from the 2026-10-03 build profiling) plus its folder for that
time; the two subset builds are correspondingly cheaper. A rerun with the same inputs keeps the
ranking and both databases; other genes (`--genes`) rebuild both databases from the samples
already simulated, since the collector reuses samples simulated with the same design.

## Tests

- `GeneSubsetTest` (`scripts/mini_db/test_mini_db.py`): a six-species operon release; a subset of
  every other gene along one representative's longest contig, named three ways. The copy holds
  those genes with their ids in every file, the direct conversion gives the same files, the
  neighbours of the copy are exactly what `neighbour_ends` gives on the representative's
  placements restricted to the subset, and the genes whose neighbour was left out face the next
  kept gene with a larger gap. With `$PROTAL`: the subset's build lists only its genes in
  `unique_kmers.tsv`, a full build is ranked (table shape, bounds, order, markers), the list is
  taken by the converter and by `--build_gene_subset`, which is refused with the neighbours table
  and accepted without it.
- `GtdbBuildTest.test_d`: `build_gtdb_database.py --n-genes 3` on the 60-species release: the
  ranking build, the derived folders, the models, the metadata, the rerun that builds nothing, and
  `--genes` with other genes rebuilding both databases without simulating again.
- The C++ tests of the index, the database, the neighbours and parsing (42) pass.

Commands (WSL, the checkout rsynced to `~/protal-pack/src`, built with Ninja):

```bash
cmake --build ~/protal-pack/build --target protal simulate_metagenomes protal_tests
PROTAL=~/protal-pack/build/protal SIMULATE=~/protal-pack/build/simulate_metagenomes \
  python3 -m unittest -v scripts.mini_db.test_mini_db.GeneSubsetTest scripts.mini_db.test_mini_db.MiniDbTest \
  scripts.mini_db.test_mini_db.GeneNeighboursTest
PROTAL=... SIMULATE=... PROTAL_TRAIN_PYTHON=~/micromamba/envs/protal-db-build/bin/python \
  python3 -m unittest -v scripts.mini_db.test_mini_db.GtdbBuildTest.test_d_a_reduced_database_of_the_most_distinctive_genes
```

Results (2026-10-04, WSL, 6 threads): `GeneSubsetTest` 4 tests OK in 17 s (three protal builds
of the mini database among them); `MiniDbTest`, `GeneNeighboursTest`, `SpeciesLinesTest`,
`CircularGeneNeighboursTest`, `BuildOptionsTest` OK (47 tests, 40 s); `GtdbBuildTest.test_d` OK
in 53 s (three pipeline runs: `--n-genes 3`, its rerun, `--genes`); the 42 C++ tests of
`Index|Database|GeneNeighbours|Operations|Parsing` OK. On the 60-species test release the
ranking build took 5 s and ranked 168 genes; the three best were PF00410.20, PF00466.21 and
TIGR00064 (scores 0.828 down to 0.772).

## Open

- No GTDB-scale run yet: the r226 ranking (which 12 genes, their scores) and the F1 of a
  12-gene database's models are to be measured on the HPC (`--n-genes 12`, or `rank_genes.py` on
  the unpacked v9 database with `--gene-ranking` to skip the ranking build).
- The website's documentation page describes the old two ways (`subset_genes.py` and
  `--build_gene_subset` without the neighbours caveat) and is out of date.
- The score takes no account of gene length. If the goal is the most sensitivity per gene rather
  than per k-mer, weight by length from the table's `mean_length` column and pass the list with
  `--genes`.

## Follow-up the same day: archaea in the subset, and the two-phase scripts

**Archaea were not covered.** GTDB's bacterial and archaeal marker sets (120 and 53 genes,
`scripts/mini_db/markers_r226.tsv`) share five genes; archaea are a few percent of the species
(r226: about 5,000 of 143,000). Prevalence over all species puts every gene of the archaeal set
alone below 0.05, so the twelve best by the overall score would have been bacterial genes and the
database would have had no gene of most archaeal genomes. On the test release (120 bacterial, 48
archaeal-only genes in the reference) the three best by the overall score are indeed bacterial
only (`GeneSubsetTest` checks it with `--per-domain 0`).

Fix in `rank_genes.py`: the table carries prevalence and score per domain (`bacteria_*`,
`archaea_*`, from `internal_taxonomy.dmp`), and `--top N` reserves `--per-domain M` of the N
genes for each domain's best (default a third of N, at least 1: 4 of 12), the rest by the
overall score, a gene good for both counting for both; the ranker reports how many of the genes
chosen half or more of each domain's species have. `build_gtdb_database.py --genes-per-domain M`
passes it on; `--rank-genes` makes a run with every gene write the ranking (the training database
unpacked once built), for the reduced run to take with `--gene-ranking` (no ranking build then).
The test release's `--n-genes 3` run now chooses the two best overall and the best archaeal gene;
`gene_subset.txt` says which was chosen for which domain, and `build_metadata.tsv` how many
genes each domain has.

**Two phases across releases.** `scripts/download_gtdb_releases.py` (a node with internet) runs
`download_gtdb.py` per release into `INPUTS/gtdb_r<release>` and writes `download_summary.tsv`:
version, the taxonomy's genomes and clades per rank and the species per domain, genomes
delivered and missing, species and strains to simulate from, marker files, bytes, time.
`scripts/build_gtdb_releases.py` (no internet) runs `build_gtdb_database.py` per release and
variant (`full` with `--rank-genes`, then `n12` with `--gene-ranking` from it) into
`OUT/r<release>_<variant>`, keeps what is built, and writes `build_summary.tsv` and `.txt`: protal
version, marker genes, taxa per rank and per domain (from the database's taxonomy), genomes
simulated from, species held out, `database.protal` size, the build's time and peak memory, the
profiling run's peak memory (collector and protal together, from the run's console log), the run's
time, and per read type the test-set F1, false positives per sample, sensitivity and precision
and the held-out F1 and false positives per sample (`model_logs/summary.txt`). Both pass the
options they do not know on to the script they run.

Tests (all on the WSL box, 2026-10-04): `GeneSubsetTest` with two archaeal species added (4
tests, 30 s), `MiniDbTest.test_download_releases` (the download phase on the fake mirror, with a
failing release reported), `GtdbBuildTest.test_d` (the archaeal gene among the three chosen) and
`test_e` (`build_gtdb_releases.py` with `n3,full`: the full database ranked, the reduced one from
that ranking, the summary's columns, a rerun that keeps both; 82 s). Still no GTDB-scale run.
