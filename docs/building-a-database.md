# Building a database

The database on the [downloads page](https://protal.earlham.ac.uk/main.php?site=downloads) is built
from GTDB r226. You can build one yourself, from another GTDB release or from part of one. A
database has three ingredients:

- the marker genes of each species' representative genome: the reference reads are aligned to;
- the marker genes of all genomes: used to decide which k-mers are unique to a species;
- a presence model trained on simulated samples.

Two routes lead there: step by step with the converter and `protal --build`
([below](#step-by-step)), or in one command that also trains the model
([build and train](#build-and-train-in-one-command)). The small databases used for testing are
made the same way from a synthetic release, see [testing.md](testing.md#mini-database).

## GTDB files

Download the release from https://data.gtdb.ecogenomic.org/releases/ (for r226:
`release226/226.0/`) and keep GTDB's layout:

| File | Needed for |
|---|---|
| `bac120_taxonomy_r226.tsv`, `ar53_taxonomy_r226.tsv` (or `.tsv.gz`) | lineages |
| `bac120_metadata_r226.tsv.gz`, `ar53_metadata_r226.tsv.gz` | species representatives (optional: without them, the genomes in the `*_marker_genes_reps_*` files are the representatives) |
| `genomic_files_reps/{bac120,ar53}_marker_genes_reps_r226.tar.gz` | the reference |
| `genomic_files_all/{bac120,ar53}_marker_genes_all_r226.tar.gz` | unique k-mers (optional, but without it uniqueness is judged on the representatives alone) |
| `genomic_files_reps/gtdb_genomes_reps_r226.tar.gz` | whole genomes, only to simulate training data |

Extract the tarballs where they are, so that for example
`genomic_files_reps/bac120_marker_genes_reps_r226/fna/*.fna` exists. The converter takes the marker
id (`PFxxxxx.x` or `TIGRxxxxx`) from each FASTA's file name and the genome accession from each
record's header.

## Step by step

### 1. Convert the release

```bash
python3 scripts/mini_db/gtdb_to_protal_db.py --gtdb /data/gtdb_r226 --outdir /data/protal_r226_db
```

| Option | Default | |
|---|---|---|
| `--release` | detected from the taxonomy file names | e.g. `226` |
| `--model` | `scripts/random_forest.xml` | PMML model copied to `model.xml` |
| `--order` | `gene` | order of `reference.fna`: by gene, then taxon (compresses about 2x better), or `genome` |
| `-t, --threads` | 1 | marker files read in parallel; the output is the same for any number |
| `--exclude_species` | | file of species whose marker genes are left out; the taxonomy keeps them with their taxids (a training database) |
| `--from_db` | | instead of `--gtdb`: copy a folder this script wrote, without the species of `--exclude_species` |

It writes `reference.fna`, `reference.map`, `internal_taxonomy.dmp`, `full_reference.fna` (only if
`genomic_files_all` is there), `model.xml`, and two tables for your own use: `gene2geneid.tsv`
(marker id to protal's gene id) and `genome2tiid.tsv` (accession, species taxid, representative,
lineage). The converter spools each marker's genes to a temporary folder in the output, so it
holds one gene's sequences per worker in memory, not the whole release.

The default `model.xml` is the model shipped with protal. It was trained on older databases and
does not call archaea reliably ([model-training.md](model-training.md)); for a database you will
use, train a model on it, or use the one-command route below.

### 2. Build the index

```bash
protal --build --no_profile -t 16 --db /data/protal_r226_db \
    --reference /data/protal_r226_db/reference.fna \
    --full_reference /data/protal_r226_db/full_reference.fna
```

This indexes `reference.fna`, checks every k-mer's uniqueness against `full_reference.fna`
(`--reference` is used when `--full_reference` is not given), writes `unique_kmers.tsv`, and packs
everything into `database.protal`, which it reads back and compares before it removes the separate
files. `full_reference.fna`, `gene2geneid.tsv` and `genome2tiid.tsv` stay next to it.
[database-files.md](database-files.md#build-options-for-the-format) lists the options for
separate or uncompressed files and for the compression level. Pass `--no_profile`: without it,
build mode goes on to profile an empty sample list.

### 3. Check it

Profile samples simulated from genomes whose species you know, and compare. With
[`simulate_metagenomes`](simulation.md) this is one step, since its map carries the truth:

```bash
simulate_metagenomes --genome_table genomes.tsv -o sim -n 4 --total_read_pairs 200000 \
    --species_per_sample 20 --seed 1 -t 8 --protal_metafile sim/protal
protal --db /data/protal_r226_db --map sim/protal.meta -t 16
```

For each sample protal prints the true and false positives and the false negatives, and writes
`<profile>.truth_annotated`.

## Reduced marker sets

The reduced database on the downloads page uses a subset of the marker genes, for about five
times less memory. Two ways to build one:

- `--build_gene_subset genes.txt` with `protal --build` indexes only the genes listed, one gene id
  per line (the ids of `gene2geneid.tsv`, `#` lines are comments). `reference.fna` keeps all genes,
  but only the listed ones are in the index and can be hit.
- `scripts/subset_genes.py --ref reference.fna --build_gene_subset genes.txt --out subset.fna`
  writes a smaller reference; `scripts/index_reference.py subset.fna > reference.map` writes its
  map. Build from those.

A model trained on the full marker set sees fewer genes per species with a subset. Train one on
the reduced database, or check its calls on simulated samples first.

## Build and train in one command

`scripts/build_gtdb_database.py` runs the converter, builds and packs the index, simulates training
data from whole genomes, trains a random forest on the normalised features, and packs the database
again with the new `model.xml`. Its inputs come from `scripts/download_gtdb.py`, the one step that
needs the internet, so it can run on a download node:

```bash
python3 scripts/download_gtdb.py -o /shared/protal_inputs/gtdb_r226          # on a node with internet
python3 scripts/build_gtdb_database.py --inputs /shared/protal_inputs/gtdb_r226 --outdir /data/protal_r226 \
    --protal build/protal --simulator build/simulate_metagenomes -t 16
```

`download_gtdb.py` fetches GTDB r226 (`--release 220`, or another release from 207 on, for older
ones; the newest point release unless one is named, e.g. `214.1`): taxonomy, metadata and the
marker genes of the representatives and of all genomes (17.7 GB for r226), checked against the
release's `MD5SUM.txt` and extracted. Then it downloads the genomes to simulate from with NCBI's
`datasets` CLI (below). The folder serves every later build of that release; a rerun downloads
only what is missing, and `--dry_run` lists what it would fetch.

| `download_gtdb.py` writes | |
|---|---|
| `release/` | GTDB's files as `--gtdb` reads them (archives removed once extracted; `--keep_archives`) |
| `genomes/` | `<accession>.fna.gz` from NCBI: other strains, and the representatives of the species to simulate |
| `genomes.tsv`, `missing.txt` | the genomes there (species, role, lineage, CheckM2 values), and those NCBI did not deliver |
| `simulation_species.txt` | the species to simulate from |
| `download.json` | the release, the options, the counts, the checksums |

| Option | Default | |
|---|---|---|
| `-o` | required | the folder |
| `--release` | 226 | GTDB release (207 or later) or point release |
| `--rep_genomes` | `ncbi` | representatives' genomes of the simulated species from NCBI, or `gtdb`: GTDB's archive of all of them (137 GB for r226) |
| `--species`, `--per_species`, `--rep_only_species` | 4000, 2, 1000 | species with strains, strains each, species simulated from their representative only |
| `--min_completeness`, `--max_contamination` | 90, 5 | CheckM2 filters for strains |
| `--no_genomes`, `--dry_run` | | GTDB's files only; list the files and their sizes |
| `--mirror`, `--datasets`, `--batch`, `-t` | | GTDB server, NCBI CLI, genomes per NCBI request (500), parallel downloads and compression (8) |

It needs, besides a built `protal` and `simulate_metagenomes`: `art_illumina` and `pigz` on
`$PATH` (for the simulations), Python 3 with numpy, pandas, joblib and scikit-learn for the
training (no Java: the trainer writes the PMML itself), and NCBI's `datasets` for strain genomes.
The conda environment [`envs/protal-db-build.yaml`](../envs/protal-db-build.yaml) has them all,
with the compilers to build protal ([installation.md](installation.md#tools-to-build-a-database)).

With `--gtdb` instead of `--inputs`, whole genomes are looked up under
`genomic_files_all/gtdb_genomes_all_r<R>` or `genomic_files_reps/gtdb_genomes_reps_r<R>` (GTDB
publishes only the representatives' genomes) and in `--extra-genomes`. If they are elsewhere, pass
`--genome-table` in the simulator's three-column format (accession, GTDB taxonomy, FASTA path).

### Training data like real samples

Two things make the simulated samples like real ones, and on a GTDB-like simulated world the model
needed both ([report](claude/2026-09-29-model-training-tuning/README.md)): without them it called
a quarter of the relatives of species the database lacks present and misjudged itself completely;
with only one of them it either kept 2 false positives per sample or missed a fifth of the strains.

- **Other strains.** With the representatives alone, every simulated species is the database's
  own reference, whose reads match it more closely than those of real strains. `download_gtdb.py`
  picks, per domain in GTDB's proportions, 4,000 species with non-representative genomes and up to
  2 of those each, plus 1,000 species simulated from their representative only (about 13,000
  genomes, ~50 GB as `.fna.gz`), and `--inputs` simulates from exactly these species. The pool
  matters: the simulator draws species uniformly, and among all GTDB species the few with strains
  would hardly be drawn; with it, a simulated species is another strain than its representative
  about 53% of the time. `genome_table.txt` says how often, and the training report warns when
  nearly every present taxon has reads identical to its reference.
- **Species the database lacks.** The samples are profiled against a training database
  (`training_db/`) that leaves `--holdout` of the species out; their reads land on relatives, as
  those of species GTDB lacks do in real samples. The finished database has all species, and the
  model trained so (the training database keeps their taxids). It costs a second index build.

| Option | Default | |
|---|---|---|
| `--inputs` or `--gtdb` | required | a folder of `download_gtdb.py`, or an extracted release |
| `--outdir` | required | the output root |
| `--release` | detected | GTDB release number |
| `--genome-table` | built from the release | genomes to simulate from |
| `--extra-genomes` | | folder of more genomes of GTDB species (by accession in the file names); repeatable |
| `--simulate-species` | all | file of the species to simulate from |
| `--holdout`, `--holdout-species` | 0.1 | share of species left out of the training database, or a file naming them; 0 trains on the database itself |
| `--training-db-level` | 3 | zstd level of the training database (only read while training; level 19 would take about half its build) |
| `--protal`, `--simulator` | `protal`, `simulate_metagenomes` | the binaries |
| `-t, --threads` | 8 | |
| `--samples` | 8 | samples per design point |
| `--congeners` | 0 | species of one genus in every sample of a design point |
| `--no-placeholder-models` | | leave out the placeholder models for se, pb and ont (below) |
| `--read-pairs` | `5000,20000,100000,500000` | depths, one design point each |
| `--read-setups` | `100:HS20:300:40,150:HS25:350:50,250:MSv3:550:50` | read length : ART profile : fragment mean : fragment SD, one design point each |
| `--species-per-sample` | `10-40` | |
| `--archaea` | 2 | archaeal species per sample |
| `--seed` | 1 | |
| `--ntree`, `--maxnodes` | 64, 128 | random forest size (more trees did not score better, see [model-training.md](model-training.md#training)) |
| `--evaluation` | `full` | how much the trainer evaluates: `full`, `basic` or `none` |

With the defaults that is 3 read setups x 4 depths x 8 samples = 96 simulated samples, each
profiled against the training database. The simulations take most of the compute and disk space;
training and its evaluation take minutes. Whether 96 samples are enough, the training report's
learning curve says.

The database gets one model per read type: the trained one as `model_pe.xml`, and placeholders
(`scripts/placeholder_models.py`) as `model_se.xml`, `model_PB.xml` and `model_ONT.xml`, since only
paired-end reads can be simulated and aligned for training so far. A placeholder scores every taxon
0, so no species is reported (`--knob 0` lists every taxon with reads), and protal warns whenever
it loads one; replace it with `protal --add_model MODEL --read_type se --db DB`.

The output root holds:

| Path | |
|---|---|
| `protal_db/database.protal` | the finished database, and `protal_db/build_metadata.tsv` (GTDB release, feature set, forest size) |
| `genomes.tsv`, `genome_table.txt` | the genome table used for the simulations, and what it holds (species by domain, how often a simulated species is not its representative) |
| `training_db/`, `heldout_species.txt` | the training database and the species it leaves out |
| `training/` | the simulated samples, their profiles and `training_data.tsv`; a rerun skips the design points already done |
| `trained_model.*` | the model and the trainer's outputs ([model-training.md](model-training.md#training)) |
| `model_logs/` | what tells whether the model is good, in one folder: the training report and its numbers (`trained_model.report.txt`, `.metrics.json`), per-taxon predictions, the threshold table, feature importances, the parity check with protal (`parity.txt`), `genome_table.txt`, the collection log and `build_metadata.tsv` |
| `*.log` | one log per stage: `convert`, `training_db`, `index_and_package`, `training_db_index`, `training_data`, `classifier_training`, `parity`, `final_package` |
