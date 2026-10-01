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
| `--model` | `scripts/random_forest.xml` | PMML model copied to `model_pe.xml` (the model of paired-end reads) |
| `--order` | `gene` | order of `reference.fna`: by gene, then taxon (compresses about 2x better), or `genome` |
| `-t, --threads` | 1 | marker files read in parallel; the output is the same for any number |
| `--exclude_species` | | file of species whose marker genes are left out; the taxonomy keeps them with their taxids (a training database) |
| `--from_db` | | instead of `--gtdb`: copy a folder this script wrote, without the species of `--exclude_species` |

It writes `reference.fna`, `reference.map`, `internal_taxonomy.dmp`, `full_reference.fna` (only if
`genomic_files_all` is there), `model_pe.xml`, and two tables for your own use: `gene2geneid.tsv`
(marker id to protal's gene id) and `genome2tiid.tsv` (accession, species taxid, representative,
lineage). The converter spools each marker's genes to a temporary folder in the output, so it
holds one gene's sequences per worker in memory, not the whole release.

The default `model_pe.xml` is the model shipped with protal. It was trained on older databases and
does not call archaea reliably ([model-training.md](model-training.md)); for a database you will
use, train a model on it, or use the one-command route below.

### 2. Build the index

```bash
protal --build --no_profile -t 16 --db /data/protal_r226_db \
    --reference /data/protal_r226_db/reference.fna \
    --full_reference /data/protal_r226_db/full_reference.fna
```

This indexes `reference.fna`, checks every k-mer's uniqueness against `full_reference.fna`
(`--reference` is used when `--full_reference` is not given), writes `unique_kmers.tsv` and
`gene_conservation.tsv` (below), and packs everything into `database.protal`, which it reads back
and compares before it removes the separate files. `full_reference.fna`, `gene2geneid.tsv` and
`genome2tiid.tsv` stay next to it.

`gene_conservation.tsv` says how fast each gene diverges within species compared with the species'
other genes. It is in every database built with other genomes' copies, so that the depth identity
margin can be scaled by it, wider on fast genes and narrower on conserved ones, without a rebuild
(`--gene_conservation db`; by default the margin is the same on every gene,
[running.md](running.md#options-the-website-does-not-list)). For
each species and gene, the build compares up to 16 other genomes' copies in `full_reference.fna`
with the representative's (k-mer distance, k = 12); divided by the distance of the species' median
gene, this cancels how far a species' strains are from its representative. A gene's factor is the
median of that ratio over the species with 10 genes or more and strains at least 0.2% from the
representative, scaled to 1 for the median gene and shrunk towards 1 when few species inform it.
The log reports the range (`Gene conservation: factors 0.26-3.3 for 168 genes, from 120 species`);
`protal --unpack_db` writes the table out, and `build_gtdb_database.py` records the range in
`build_metadata.tsv`. Without other genomes' copies (no `--full_reference`, or one genome per
species) the build writes no table. On two simulated worlds whose genes evolve at different rates
the factors correlated 0.985 and 0.999 with the true rates
([report](claude/2026-10-01-gene-scaled-margin/README.md)).
[database-files.md](database-files.md#build-options-for-the-format) lists the options for
separate or uncompressed files and for the compression level. Pass `--no_profile`: without it,
build mode goes on to profile an empty sample list.

Every phase uses `-t` threads: the two passes over `reference.fna` that fill the index (each
thread owns a range of k-mers, so the index is the same for any `-t`), the value pointers, the
uniqueness check, the gene conservation (a second pass over `full_reference.fna`), the unique k-mer
statistics and the compression. The log times each phase (`Pass 1 (count the k-mers) took ...`,
`Value pointers`, `Pass 2 (place the values)`, `Uniqueness check`, `Gene conservation`, `Unique k-mer
statistics`, `Write index`) and counts the uniqueness check's work (k-mers, flex parts compared per
k-mer, k-mers found under another taxon, entries read back from their genes): at GTDB scale that
tells where a build spends its time ([report](claude/2026-09-30-index-build-gains/README.md)).

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
data from whole genomes, trains a random forest per read type on the normalised features, and adds
the trained models to the database (`model_pe.xml`, `model_se.xml`, `model_PB.xml`, `model_ONT.xml`). Its inputs come from `scripts/download_gtdb.py`, the one step that
needs the internet, so it can run on a download node:

```bash
python3 scripts/download_gtdb.py -o /shared/protal_inputs/gtdb_r226          # on a node with internet
python3 scripts/build_gtdb_database.py --inputs /shared/protal_inputs/gtdb_r226 --outdir /data/protal_r226 \
    --protal build/protal --simulator build/simulate_metagenomes -t 16
```

`download_gtdb.py` fetches GTDB r226 (`--release 220`, or another release from 207 on, for older
ones; the newest point release unless one is named, e.g. `214.1`): taxonomy, metadata and the
marker genes of the representatives and of all genomes (17.7 GB for r226), checked against the
release's `MD5SUM.txt` and extracted. Then it downloads the genomes to simulate from NCBI. The
folder serves every later build of that release; a rerun downloads
only what is missing, and `--dry_run` lists what it would fetch. A download whose connection drops
goes on from where it stopped (a `.part` file, also on a rerun), and one that stops after the last
byte is checked and kept. With NCBI not answering, the run stops after a few failed requests in a
row, and it fails if no genome arrived.

The genomes come from NCBI's FTP server (`--ftp_url`), `--connections` (8) at a time: the address
of a genome is built from its accession and the assembly name in GTDB's metadata, the file is kept as
NCBI compresses it, and it is checked against its announced length and read through gzip. A request
that fails is repeated after growing waits (or the wait the server asks for); a missing file is
not. What this does not deliver (an assembly NCBI renamed or withdrew, or NCBI not answering) goes
through NCBI's `datasets` CLI, a batch of `--batch` genomes at a time, which asks NCBI for each file
and is slower. The log prints the rate (MB/s, genomes/s) as it goes: raise `--connections` if the
link allows, and lower it if the log says NCBI limited the requests (it answered HTTP 503 to 32
connections at once in a test). `--connections 0` uses `datasets` only. On one test link (about
3 MB/s in total, so bandwidth-bound), 4 connections fetched 48 genomes in a third of the time of the
`datasets` route and one connection in four fifths of it; a faster link gains more from more connections.

| `download_gtdb.py` writes | |
|---|---|
| `release/` | GTDB's files as `--gtdb` reads them (archives removed once extracted; `--keep_archives`) |
| `genomes/` | `<accession>.fna.gz` from NCBI: other strains, and the representatives of the species to simulate. A build reads all of it, so a rerun that chooses other genomes (other options, or the quality ranking below on inputs downloaded before it) moves the ones it no longer wants to `genomes_unused/` |
| `genomes.tsv`, `missing.txt` | the genomes there (species, role, lineage, CheckM2 values, genome category, assembly level, contig count and, for strains, the sequencing technology at NCBI), and those NCBI did not deliver |
| `ncbi_info.tsv` | the sequencing technology NCBI gave for the candidate strains (a rerun asks for none it has) |
| `simulation_species.txt` | the species to simulate from |
| `download.json` | the release, the options, the counts, the checksums |

| Option | Default | |
|---|---|---|
| `-o` | required | the folder |
| `--release` | 226 | GTDB release (207 or later) or point release |
| `--rep_genomes` | `ncbi` | representatives' genomes of the simulated species from NCBI, or `gtdb`: GTDB's archive of all of them (137 GB for r226) |
| `--species`, `--per_species`, `--rep_only_species` | 6000, 2, 2000 | species with strains, strains each, species simulated from their representative only |
| `--min_completeness`, `--max_contamination` | 90, 5 | CheckM2 filters for strains |
| `--tech_candidates` | 30 | candidate strains per species whose sequencing technology is asked of NCBI (the best by category and assembly level) |
| `--no_tech_lookup` | | do not ask NCBI for sequencing technologies: no preference for PacBio and Nanopore assemblies |
| `--progenomes` | | a proGenomes ANI-clustering table ([`pg4_ANI_clustering.tsv.gz`](https://progenomes.embl.de/download.cgi), 5 MB, downloaded by hand) or any list of accessions: strains are taken from these genomes only |
| `--no_genomes`, `--dry_run` | | GTDB's files only; list the files and their sizes |
| `--connections`, `--ftp_url` | 8, `https://ftp.ncbi.nlm.nih.gov/genomes/all` | genomes fetched at a time straight from NCBI's FTP server (0: only through `datasets`), and that server |
| `--mirror`, `--datasets`, `--batch`, `-t` | | GTDB server, NCBI CLI, genomes per `datasets` request (500), parallel compression of what `datasets` delivers (8) |

It needs, besides a built `protal` and `simulate_metagenomes`: `art_illumina` on
`$PATH` (for the simulations; pbsim3 for PacBio and Nanopore reads), Python 3 with numpy, pandas, joblib and scikit-learn for the
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
  picks, per domain in GTDB's proportions, 6,000 species with non-representative genomes and up to
  2 of those each, plus 2,000 species simulated from their representative only (about 20,000
  genomes, ~75 GB as `.fna.gz`), and `--inputs` simulates from exactly these species. The pool
  matters: the simulator draws species uniformly, and among all GTDB species the few with strains
  would hardly be drawn; with it, a simulated species is another strain than its representative
  about 50% of the time. `genome_table.txt` says how often, and the training report warns when
  nearly every present taxon has reads identical to its reference.

  Which strains of a species are taken is a choice by quality, not by chance. Among the strains
  that pass the CheckM2 filters, in this order: isolate genomes before single-cell and
  metagenome-assembled ones (GTDB's `ncbi_genome_category`, and an "uncultured" or "metagenome"
  organism name); complete genomes before chromosome-, scaffold- and contig-level assemblies;
  PacBio or Nanopore assemblies (NCBI's sequencing technology, asked for the best 30 candidates of
  each species, in at most a few hundred batched requests) before the others; fewer contigs.
  Genomes that tie are drawn at random, so species with thousands of complete genomes are not
  always represented by the same two. A species with MAGs only keeps its best MAGs. The
  representatives are GTDB's, the database's references, whatever their quality. The download
  prints how many strains are isolates, at which assembly level, and how many are long-read
  assemblies (`download.json` keeps the counts). Contigs under 100 bases are dropped from a genome
  before pbsim3 simulates long reads from it, since pbsim3 stops at them.
- **Species the database lacks.** The samples are profiled against a training database
  (`training_db/`) that leaves species out; their reads land on relatives, as those of species GTDB
  lacks do in real samples. The finished database has all species, and the model trained so (the
  training database keeps their taxids). It costs a second index build, which runs while the
  first one does (below).

Besides, samples are as complex as real ones: 20-200 species each (`--species-per-sample`; gut
samples hold 100-300 GTDB species, most of them rare), a second and third strain of a species in
30% and 10% of cases (`--strains-per-species 0.3,0.1`; mixed strains change the allele-frequency
features), and reads from 100 bp (HiSeq 2000) to 250 bp (MiSeq), with HiSeq X, the closest of ART's
built-in profiles to NovaSeq, for 150 bp. A read setup `150:file=R1.txt+R2.txt:350:50` uses quality
profiles that `art_profiler_illumina` made from real reads instead, e.g. from a NovaSeq run.

The training database lacks whole clades at every rank as well as many single species, since
false positives from organisms the database lacks are what a presence model must avoid. First
`--holdout-clades` (default 2 phyla, 4 classes, 6 orders, 8 families, 12 genera) are drawn from
phylum down: clades with at least two species to simulate, at most `--holdout-max-share` (2%) of the
database's species, and not inside a clade drawn before. Then `--holdout` (20%) of the species to
simulate that no clade took. Reads of an organism whose genus the database has land on congeners;
those of one whose genus, family, order, class or phylum it lacks land on ever more distant
relatives, or nowhere, and the model has to learn all of it. So that such species are in the
samples, each design point takes one held-out clade of each rank, in turn (every clade is used
before one is used again), and puts `--novel-clades-per-sample` (1) of its species into every
sample. `heldout_species.txt` lists the species with the rank they were held out at and the clade,
`model_logs/holdout.txt` sums them up, and the training report gives false positive and false
negative rates by rank (see [model-training.md](model-training.md#species-and-clades-the-database-lacks)).

### One model per read type

`--read-types` (default `pe,se,pb,ont`) trains a model for each kind of reads protal profiles,
from the same communities: paired-end reads (ART); their first reads alone, profiled as single-end
reads; and PacBio and Nanopore reads simulated with [pbsim3](https://github.com/yukiteruono/pbsim3)
(`--pb-setup`: HiFi-like reads of the Sequel error model at 99.9% accuracy, 15 kb; `--ont-setup`:
the high-quality ONT model at 97%, 8 kb), `--long-read-bases` per sample, about as many bases as the
paired-end depths. All samples are profiled in one protal run, each with its read type's settings,
and the models are trained in parallel. pbsim3 must be installed for pb and ont (it is in
`envs/protal-db-build.yaml`; `micromamba install -c conda-forge -c bioconda pbsim3`); without it,
leave them out of `--read-types`, and they keep placeholder models.

### An independent test set

Cross-validation on the training data judges the model only on what the training design has: a
model trained on samples of 5,000 read pairs and more missed 8% of the present taxa of 1,000-pair
test samples while its own estimate said F1 0.997
([report](claude/2026-09-30-clade-holdouts.md)). So every build also profiles an independent test
set of another design, `--test-samples` (4) per design point: depths `--test-read-pairs`
(500 to 1,000,000), 10-300 species (`--test-species-per-sample`), more uneven abundances
(`--test-abundance lognormal:2.0`), more mixed strains (`--test-strains-per-species 0.5,0.2`),
long-read depths `--test-long-read-bases`, another seed. Each model scores it (the trainer's
`--test-file`): the report's section "Independent test set" gives F1, false positives per sample,
FN and FP rates by depth and by rank, and the threshold with the highest F1 there, and warns when
the test set scores clearly worse than cross-validation. `--test-samples 0` skips it.

| Option | Default | |
|---|---|---|
| `--inputs` or `--gtdb` | required | a folder of `download_gtdb.py`, or an extracted release |
| `--outdir` | required | the output root |
| `--release` | detected | GTDB release number |
| `--genome-table` | built from the release | genomes to simulate from |
| `--extra-genomes` | | folder of more genomes of GTDB species (by accession in the file names); repeatable |
| `--simulate-species` | all | file of the species to simulate from |
| `--holdout` | 0.2 | share of the species to simulate (that no held-out clade took) left out of the training database |
| `--holdout-clades` | `phylum:2,class:4,order:6,family:8,genus:12` | whole clades left out of the training database, `RANK:COUNT` for phylum, class, order, family, genus; `none` for single species only (with `--holdout 0`: train on the database itself) |
| `--holdout-max-share` | 0.02 | largest share of the database's species a held-out clade may have |
| `--novel-clades-per-sample` | 1 | species of held-out clades in every sample, per rank |
| `--holdout-species` | | a file of the species to leave out instead (optionally with rank and clade, as `heldout_species.txt`) |
| `--one-build-at-a-time` | | build the finished database after training, not while the training data are collected |
| `--training-db-level` | 3 | zstd level of the training database (only read while training; level 19 would take about half its build) |
| `--protal`, `--simulator` | `protal`, `simulate_metagenomes` | the binaries |
| `-t, --threads` | 8 | |
| `--samples` | 12 | samples per design point |
| `--congeners` | 0 | species of one genus in every sample of a design point |
| `--no-placeholder-models` | | leave out the placeholder models of read types not trained (below) |
| `--read-pairs` | `1000,5000,20000,100000,500000` | depths, one design point each (without the shallowest, a model missed 8% of the present taxa of 1000-pair samples) |
| `--read-setups` | `100:HS20:300:40,150:HSXt:350:50,250:MSv3:550:50` | read length : ART profile (or `file=R1.txt+R2.txt`) : fragment mean : fragment SD, one design point each |
| `--species-per-sample` | `20-200` | drawn per sample |
| `--strains-per-species` | `0.3,0.1` | probabilities of a second, third, ... strain of a species |
| `--abundance` | the simulator's | `lognormal:SIGMA`, `powerlaw:ALPHA` or `negbin:R:P` |
| `--archaea` | 2 | archaeal species per sample |
| `--read-types` | `pe,se,pb,ont` | the read types to train a model for |
| `--long-read-bases` | `300000,1500000,6000000,30000000,150000000` | bases per pb and ont sample, one design point each |
| `--pb-setup`, `--ont-setup` | `errhmm:ERRHMM-SEQUEL:15000:3000:0.999`, `qshmm:QSHMM-ONT-HQ:8000:6000:0.97:39/24/36` | pbsim3 method : model : length mean : length SD : accuracy (: substitution/insertion/deletion mix, for qshmm) |
| `--pbsim`, `--pbsim-models` | `pbsim`, found next to it | pbsim3 and its models |
| `--test-samples` | 4 | samples per design point of the independent test set; 0 for none |
| `--test-read-pairs`, `--test-species-per-sample`, `--test-abundance`, `--test-strains-per-species`, `--test-long-read-bases` | `500,2000,10000,50000,200000,1000000`, `10-300`, `lognormal:2.0`, `0.5,0.2`, `150000,1000000,5000000,25000000,250000000` | the test set's design |
| `--seed` | 1 | |
| `--ntree`, `--maxnodes` | 64, 128 | random forest size (more trees did not score better, see [model-training.md](model-training.md#training)) |
| `--evaluation` | `full` | how much the trainer evaluates: `full`, `basic` or `none` |

With the defaults that is 3 read setups x 5 depths x 12 samples = 180 paired-end samples of 20-200
species, profiled as paired-end and as single-end reads, 5 x 12 PacBio and 5 x 12 Nanopore samples
of the same communities, and a test set of 3 x 6 x 4 = 72 paired-end samples (and their single-end
and long-read counterparts), all against the training database. The simulations take most of the
compute and disk space; training and its evaluation take minutes per read type. Whether 180
samples are enough, the training report's learning curve says.

The finished database is needed only at the end, to take the trained models (`--add_model`), so it
is built in the background from the start, while the training database is built and the training
data are collected; the script waits for it before packing the models. That needs the memory of two
builds at once (about 50-60 GB each at GTDB scale before the genes were held at two bits per
base, which takes ~14 GB off each; estimated), or of one build and the collection's protal
runs. `--one-build-at-a-time` builds it after the training instead.

### Stopping and rerunning

The script checks before it starts that protal, the simulator, `art_illumina` (and pbsim3 for pb
and ont) are there and that its Python can import what the trainer needs. Each command it runs is
stopped with the processes it started when the script stops, whether a command failed, a Python
error, or `SIGTERM`, `SIGINT` (Ctrl-C) or `SIGHUP` stopped it; a rerun never races a build left
running. The build in the background is looked at every few seconds: when it fails, the run stops
then, not after the collection and the training.

A rerun into the same `--outdir`, after a failure in training, say, resumes. The conversion and each
of the two index builds are skipped when their inputs are those of the run that completed them:
the release's files and the converter for the conversion; also protal for the finished database; and
the species left out and `--training-db-level` for the training database (`.stages/`). A conversion
that completed but whose files a stopped build did not pack is used as it is. The collector reuses a
design point's samples if the simulator, the genome table and the design (seed, depth, read setup,
species per sample, held-out clades, ...) are the same, and its profiles if the database and protal
are too; else it simulates or profiles them again and says so, rather than mixing samples of an
earlier design or database into the table. So another `--seed` or `--holdout` rebuilds the training
database (from the release converted again, since the finished database's build packed the
converted files) and collects again, but keeps the finished database. The training, the parity
checks and `--add_model` run on every rerun.

The database gets the trained model of each read type in `--read-types` (`model_pe.xml`,
`model_se.xml`, `model_PB.xml`, `model_ONT.xml`), and a placeholder (`scripts/placeholder_models.py`)
for each read type left out. A placeholder scores every taxon 0, so no species is reported
(`--knob 0` lists every taxon with reads), and protal warns whenever it loads one; replace it with
`protal --add_model MODEL --read_type se --db DB`.

The output root holds:

| Path | |
|---|---|
| `protal_db/database.protal` | the finished database, and `protal_db/build_metadata.tsv`: GTDB release, date, protal version and binary, the scripts' git commit, the command, seed, genome table, what the training database leaves out, the training design, and each model's F1 on species held out and on the test set |
| `genomes.tsv`, `genome_table.txt` | the genome table used for the simulations, and what it holds (species by domain, how often a simulated species is not its representative) |
| `training_db/`, `heldout_species.txt` | the training database and the species it leaves out (species, the rank they were held out at, the clade) |
| `training/`, `test/` | the simulated samples, their profiles and one table per read type (`training_data.tsv` for pe, `training_data_se.tsv`, `_pb`, `_ont`); a rerun reuses the design points simulated and profiled from the same inputs (below) |
| `.stages/` | the inputs of the conversion and the two builds that completed, for a rerun (below) | 
| `trained_model.*`, `trained_model_se.*`, `_pb.*`, `_ont.*` | the models and the trainer's outputs ([model-training.md](model-training.md#training)) |
| `model_logs/` | what tells whether the models are good, in one folder: `summary.txt` (per read type: TP, FP, TN, FN, sensitivity, specificity, precision, F1 and false positives per sample, with species held out and on the test set; also printed at the end), each read type's training report and its numbers (`trained_model*.report.txt`, `.metrics.json`), per-taxon predictions (also on the test set), the threshold table, feature importances, the parity checks with protal (`parity*.txt`), `genome_table.txt`, what the training database leaves out (`holdout.txt`, `heldout_species.txt`), the collection logs and `build_metadata.tsv` |
| `*.log` | one log per stage: `convert`, `training_db`, `index_and_package`, `training_db_index`, `training_data`, `test_data`, `classifier_training*`, `parity*`, `final_package*` |
