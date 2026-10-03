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

It writes `reference.fna`, `reference.map`, `internal_taxonomy.dmp`, `full_reference.fna.zst` (only
if `genomic_files_all` is there: the marker genes of every genome, 86 GB uncompressed at r226,
compressed with the `zstd` command at level 6; `full_reference.fna` if there is no `zstd`),
`model_pe.xml`, and two tables for your own use: `gene2geneid.tsv`
(marker id to protal's gene id) and `genome2tiid.tsv` (accession, species taxid, representative,
lineage). The converter spools each marker's genes to a temporary folder in the output, so it
holds one gene's sequences per worker in memory, not the whole release.

The default `model_pe.xml` is the model shipped with protal. It was trained on older databases and
does not call archaea reliably ([model-training.md](model-training.md)); for a database you will
use, train a model on it, or use the one-command route below.

### Species priors

The converter also writes `species_priors.tsv`, what GTDB knows of each species before any read, for
the model's prior features ([model-training.md](model-training.md#reads-before-the-filters-and-the-species-priors)):
per species its taxid and representative, the marker genes found in the representative and how many of
them twice (the converter keeps a genome's first copy of a single-copy marker; a second copy is CheckM's
contamination signature), the representative's CheckM completeness and contamination from the
metadata (`checkm2_*` where the release has them, else `checkm_*`), and from GTDB's species clusters
file, `auxillary_files/sp_clusters_r226.tsv` (GTDB's spelling; `download_gtdb.py` fetches it when the
release lists it), the cluster's ANI circumscription radius, mean and minimum intra-species ANI and
number of genomes. -1 stands for unknown: a release without the clusters file or the quality columns
gives those columns -1, and a database converted by an earlier protal has no table at all, which a run
reports (`Species priors: the database has no species_priors.tsv`). `--build` packs the table into
`database.protal` ([database-files.md](database-files.md)). The log line says how many representatives
have a duplicated marker and how many species have quality and cluster values
([report](claude/2026-10-03-false-positive-fixes/README.md)).

### Gene neighbours

Optional, between conversion and build: how often each marker gene end faces which other in the
genomes of each clade. GTDB's marker files hold the accession and nothing about where a gene lies,
so this needs whole genomes: all those downloaded to simulate training data, the representatives
and the other strains:

```bash
python3 scripts/mini_db/gene_neighbours.py --db /data/protal_r226_db --genome_table genomes.tsv -t 16
```

The genome table is the simulator's (accession, GTDB taxonomy, FASTA path; `build_gtdb_database.py`
writes one as `OUT_DIR/genomes.tsv`).

**Placing.** Each database gene is a gene call of its species' representative, so it is found there
by its exact sequence on either strand. In another strain it differs by the strain's mutations and is
found by its k-mer trace: the gene's 24-mers are looked up among the genome's 24-mers at every 16th
position; the band of near diagonals (32 bases) with the most hits places the gene if it holds 10% or
more of the gene's k-mers that can hit there, and 3 at least (genes up to about 9% different).
`--min_placed` (0.8) skips a genome in which fewer of its species' genes are placed (another
assembly version).

**Counting.** For each gene's two ends (5' and 3' in its coding orientation) the script records, in
each genome, the next placed marker within `--max_gap` (3000) bases, which end of it faces this one
and how far apart they are, or that there is none. An end within `--max_gap` of its contig's end says
nothing: a fragmented assembly is no evidence that a gene has no neighbour. A circular sequence has
no ends: its last gene faces its first across the origin, and every gene end on it counts. A sequence
is read as circular if its header marks it a whole replicon (NCBI's "complete genome" or "complete
sequence", `topology=circular`), or if it is the whole genome of a species representative (a
representative in one sequence is a closed chromosome). A species counts once, so that species with
many genomes do not outweigh the others: at each gene end, the partner most of its genomes show
there (the representative's on a tie), and the end counts if it does in any of its genomes, so that
another strain's assembly fills in where the representative's contigs end. Per clade (`--ranks`,
default family, order, class, phylum, domain) this gives, for each gene end and partner, in how many
of the clade's species the end faces that partner, of how many in which the end counts (it is
*informative*): the observed frequency.

**A species' own lines.** Where a species' gene order differs from its family's, the clades'
frequencies alone would judge its own neighbours rare or unlikely (below), and a run would take the
genes next to each other on its own genome's reads for genes from elsewhere. So the script also
writes lines of the species itself (clade: its taxid, each partner 1 of 1 species) for the partners
that any of its genomes shows at a gene end (also a strain's, where the clades count the partner most
of its genomes show) whose share, smoothed over its clades as a run smooths it, is below 20%, in
clades with enough species to judge (the top one informative in 5 or more); `--no_species_lines` leaves them out, and
`--from_positions` follows what the positions file's comment says. A run reads the species as the
nearest clade of its lineage, so its own neighbours there are expected (a share of at least a
quarter) and its family's stay so (three quarters of their share). On the synthetic operon world,
whose simulator rearranges each species' clusters at random, 60,943 of the 319,535 lines were such
species lines ([report](claude/2026-10-02-phasing-and-foreign-genes/README.md)); real genomes,
whose order changes less within families, should need fewer.

The script writes two files into the folder ([database-files.md](database-files.md)), which
`--build` checks and packs:
- `gene_neighbours.tsv`, the frequencies per clade, which a run loads;
- `gene_positions.tsv`, where each gene lies in each genome (contig, its length, whether circular,
  start, end, strand, placed exactly or by its trace, the share of its k-mers that hit), which a run
  does not load. The frequencies are derived from it: `--from_positions FILE` does so without
  counting the genomes again, and `gtdb_to_protal_db.py --from_db` does so for a training database,
  without the species it leaves out, as a database does not know the gene order of an organism it
  lacks.

**Use.** A clade with few species says little on its own, so for a species a clade's share of a
pairing leans on the clades above it, the more the fewer species it has: from the top clade with data
on the gene end down to the species' family, each clade's share is (species + 3 × its parent's share)
/ (informative + 3), and on to the species itself where it has lines of its own. A family with the end
informative in 3 species counts as much as its order, one with 1 a quarter, one with 30 nine tenths.
At the nearest clade, two ends that face each other in 20%
or more are expected, in 5% or less unlikely (never seen included), and rare in between. If even the
top clade has the end informative in fewer than 5 species, a pairing seen is expected and another
unknown. Mates are paired and looked for, and long reads followed, only over expected neighbours;
[running.md](running.md#options-the-website-does-not-list) says what a run does with them. Leaving a
sparse clade out altogether, and letting the next clade with 5 species decide, was tried and lost a
species' own gene order where its family is small: on a synthetic world a quarter of the true
neighbours of species in such families came out unlikely, and a species' own pairs were barely more
often expected than a relative's ([report](claude/2026-10-02-gene-neighbour-frequencies/README.md)).

The script prints a summary: the genes placed exactly and by their trace, the genomes read as
circular, the share of gene ends with a neighbour and the gaps, and for each rank how alike its
clades' species are (the share of species that have the most common partner of a gene end) and how
many of its gene ends are informative in fewer than 3 species, leaning mostly on the rank above. Gene order changes little within families, so a
family's frequencies hold for its species without genomes.

### 2. Build the index

```bash
protal --build --no_profile -t 16 --db /data/protal_r226_db \
    --reference /data/protal_r226_db/reference.fna \
    --full_reference /data/protal_r226_db/full_reference.fna
```

This indexes `reference.fna`, checks every k-mer's uniqueness against `full_reference.fna`
(`--reference` is used when `--full_reference` is not given; for either name protal reads the `.zst`
file when there is no plain one), writes `unique_kmers.tsv` and
`gene_conservation.tsv` (below), and packs everything into `database.protal`, which it reads back
and compares before it removes the separate files. `full_reference.fna`, `gene2geneid.tsv` and
`genome2tiid.tsv` stay next to it.

`gene_conservation.tsv` says how fast each gene diverges within species compared with the species'
other genes. It is in every database built with other genomes' copies. Queries read it for two of
the model's features, which genes a taxon's reads hit by their conservation
([model-training.md](model-training.md)), and the depth identity margin can be scaled by it, wider
on fast genes and narrower on conserved ones, without a rebuild (`--gene_conservation db`; by
default the margin is the same on every gene,
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

The factors, and the model's conservation features, assume that a gene conserved within species is
conserved between them, so that a relative the database lacks aligns best on the conserved genes.
To check this on real genomes, the build also compares each gene between the representatives of
congeneric species (the genus from `internal_taxonomy.dmp`; of a genus with more than 64 species, a
fixed sample of 64, each against the next 4 of them) and writes `gene_congeners.tsv` beside the
database: per gene its within-species factor, its between-species factor (the median over pairs of
species of the gene's k-mer distance over that of the pair's median gene, the median gene 1), the
pairs and species it is from, and the share of species whose nearest congener's copy is identical,
or less than 1% apart, so that most of its reads fit both. The log sums it up (`Gene congeners:`:
the pairs and genera compared, the Spearman correlation of the two factors over the genes, and for
the genes of factor below 1 and the others the median between-species factor and the identical and
near-identical shares). Queries do not read
the file, and `database.protal` does not hold it; `build_gtdb_database.py` keeps it in
`model_logs/` and the summary in `build_metadata.tsv`
([report](claude/2026-10-01-conservation-pattern/README.md)).

The build then looks for **suspect gene copies**: it sketches every species' copy of each gene
(the 128 smallest hashes of its 12-mers) and compares each copy with the other species' copies that
share one of its 8 smallest hashes. A copy within `--suspect_copy_distance` (0.02, about the share of
bases that differ) of a copy of a species of another genus, family, order, class, phylum or domain,
while its nearest congener's copy is 0.02 farther or it has none, is suspect: a contaminating contig in
a MAG, or a transferred gene. Every present organism with that gene puts a perfect read on such a copy,
and at GTDB r226 a fifth of the false species calls were reads of a present species of another genus at
identity 0.99, recurring per reference ([report](claude/2026-10-03-false-positive-anatomy/README.md)). The suspect copies go into the database as
`suspect_copies.tsv` (taxid, geneid, the partner, the rank shared, the two distances), and a run
leaves records on them out of the evidence as if the reads had not aligned (`--keep_suspect_copies`
counts them; the log says how many per sample). Every near pair across genera (within 0.05) is
reported in `gene_incongruence.tsv` beside the database, with the verdict on each copy, for the
record and for reporting the references; `build_gtdb_database.py` keeps it in `model_logs/` and the
summary in `build_metadata.tsv`. The log sums it up (`Suspect copies: N of M gene copies (x%) of S
species ...`). `--suspect_copy_distance 0` skips the scan. A copy of a conserved gene that a whole
clade shares nearly unchanged is not suspect: its congeners' copies are as near as the other genus's.
[database-files.md](database-files.md#build-options-for-the-format) lists the options for
separate or uncompressed files and for the compression level. Pass `--no_profile`: without it,
build mode goes on to profile an empty sample list.

Every phase uses `-t` threads: the two passes over `reference.fna` that fill the index (each
thread owns a range of k-mers, so the index is the same for any `-t`), the value pointers, the
uniqueness check, the gene conservation (a second pass over `full_reference.fna`), the suspect copies
(every gene copy sketched and compared), the unique k-mer statistics and the compression. The log times each phase (`Pass 1 (count the k-mers) took ...`,
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

`scripts/build_gtdb_database.py` runs the converter, records the gene neighbours from all its
genomes ([above](#gene-neighbours); `OUT_DIR/gene_neighbours.log`,
`--no-gene-neighbours` to leave them out), builds and packs the index, simulates training
data from whole genomes, trains a random forest per read type on the normalised features and the
gene neighbours' (`--features`), and adds
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
`$PATH` (for the simulations; pbsim3 for Nanopore reads), Python 3 with numpy, pandas, joblib and scikit-learn for the
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
  assemblies (`download.json` keeps the counts). Long reads come from the contigs of 100 bases or
  more only.
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
reads; PacBio HiFi reads (`--pb-setup`: 15 kb, their quality by their length, Q50 at 5 kb to Q30 at
25 kb and Q20 at 50 kb, by `scripts/hifi_reads.py`); and Nanopore reads simulated with
[pbsim3](https://github.com/yukiteruono/pbsim3) (`--ont-setup`: the high-quality ONT model at 97%,
8 kb), `--long-read-bases` per sample, about as many bases as the paired-end depths, and
`--long-read-samples` (36) per design point, which replay the communities of the paired-end points of
the same depth of every read setup (3 x 12 of them; at r226 the Nanopore model's learning curve still fell
at 24 per point). The collector
draws each long-read sample's reads (genome, length, start, strand;
[model-training.md](model-training.md#training-data)) and `hifi_reads.py` or pbsim3 adds the errors
and qualities, in one run per sample. (Until 2026-10-02 pbsim3 made the PacBio reads too, with
quality 0 at every base; [report](claude/2026-10-02-pacbio-hifi-reads/README.md).) All samples are profiled in one protal run, each with its read type's settings,
and the models are trained in parallel. Every model also gets a knob curve over the sample's depth
(the trainer's `--depth-knobs`, [model-training.md](model-training.md#knobs-by-sample-depth);
`--depth-knob-read-types`, default `pe,se,pb,ont`, `""` for none), which `build_metadata.tsv` records
(`classifier_depth_knobs`, `model_<type>_depth_knobs`): at GTDB r226 the best threshold went from
~0.1 for samples of 1,000-5,000 read pairs to ~0.9 for 500,000, and thresholds by depth raised every
read type's test F1 (by 0.006 to 0.033; [report](claude/2026-10-02-r226-build-evaluation/README.md)).
With `--call-mode fdr` (default `curve`) every model also gets calibrated calls at a target share of
false calls per sample (the trainer's `--fdr-calls`,
[model-training.md](model-training.md#calls-at-a-target-share-of-false-calls)), which protal uses ahead
of the curve; the trainer reports both on the test set.
The deepest design points (2M and 10M read pairs, 1.5 and 6 Gb of long reads, a few samples each)
are there because real samples are that deep: the absent taxa a sample holds grow with its depth,
and a knob is only known for depths the training covered. pbsim3
must be installed for ont (it is in
`envs/protal-db-build.yaml`; `micromamba install -c conda-forge -c bioconda pbsim3`); without it,
leave ont out of `--read-types`, and it keeps a placeholder model.

### An independent test set

Cross-validation on the training data judges the model only on what the training design has: a
model trained on samples of 5,000 read pairs and more missed 8% of the present taxa of 1,000-pair
test samples while its own estimate said F1 0.997
([report](claude/2026-09-30-clade-holdouts.md)). So every build also profiles an independent test
set of another design, `--test-samples` (4) per design point: depths `--test-read-pairs`
(500 to 5,000,000, the deepest with 2 samples), 10-300 species (`--test-species-per-sample`), more uneven abundances
(`--test-abundance lognormal:2.0`), more mixed strains (`--test-strains-per-species 0.5,0.2`),
long-read depths `--test-long-read-bases` with `--test-long-read-samples` (8) per point, another seed. Each model scores it (the trainer's
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
| `--holdout` | 0.3 | share of the species to simulate (that no held-out clade took) left out of the training database (0.2 until 2026-10-03: more missing species give the model more of the false positives it must reject) |
| `--holdout-clades` | `phylum:2,class:4,order:6,family:8,genus:12` | whole clades left out of the training database, `RANK:COUNT` for phylum, class, order, family, genus; `none` for single species only (with `--holdout 0`: train on the database itself) |
| `--holdout-max-share` | 0.02 | largest share of the database's species a held-out clade may have |
| `--novel-clades-per-sample` | 1 | species of held-out clades in every sample, per rank |
| `--holdout-species` | | a file of the species to leave out instead (optionally with rank and clade, as `heldout_species.txt`) |
| `--one-build-at-a-time` | | build the finished database after training, not while the training data are collected |
| `--training-db-level` | 3 | zstd level of the training database (only read while training; level 19 would take about half its build) |
| `--final-db-level` | 9 | zstd level of the finished database: at r226 level 19 made the index 2.7% smaller than level 3 for 21 more minutes |
| `--protal`, `--simulator` | `protal`, `simulate_metagenomes` | the binaries |
| `--no-binary-check` | | run the binaries even if they were not built from the script's source (version, commit) |
| `-t, --threads` | 8 | |
| `--samples` | 12 | samples per design point |
| `--congeners` | `0.25:2-5` | relatives that share a sample, in the training data and the test set: `SHARE:MIN-MAX`, about SHARE of each sample's species in groups of MIN to MAX species of one genus, the genera drawn per sample; `N`, N species of one genus per design point; `0`, none (uniform draws, which hardly ever put congeners together). The relatives features need them: trained without, a model learns that an abundant congener means absence ([model-training.md](model-training.md#features)); `build_metadata.tsv` records the setting in `classifier_training_design` |
| `--no-placeholder-models` | | leave out the placeholder models of read types not trained (below) |
| `--no-gene-neighbours` | | do not record the gene neighbours ([above](#gene-neighbours)); protal then pairs no mates over neighbouring genes |
| `--read-pairs` | `1000,5000,20000,100000,500000,2000000:4,10000000:2,30000000:1` | depths, one design point each, `DEPTH:SAMPLES` for other samples than `--samples` (the deepest point because a model with the sample's depth as a feature cannot extrapolate past the deepest sample it saw; without the shallowest, a model missed 8% of the present taxa of 1000-pair samples; the deepest are as deep as real samples) |
| `--read-setups` | `100:HS20:300:40,150:HSXt:350:50,250:MSv3:550:50` | read length : ART profile (or `file=R1.txt+R2.txt`) : fragment mean : fragment SD, one design point each |
| `--species-per-sample` | `20-200` | drawn per sample |
| `--strains-per-species` | `0.3,0.1` | probabilities of a second, third, ... strain of a species |
| `--abundance` | `lognormal:1.3,2.0` | `lognormal:SIGMA`, `powerlaw:ALPHA` or `negbin:R:P`; `lognormal:S1,S2,...` gives a design point's samples the sigmas in turn, half 1.3 and half the test set's 2.0 by default, so that the model's depth prior does not rest on one abundance distribution |
| `--archaea` | 2 | archaeal species per sample |
| `--read-types` | `pe,se,pb,ont` | the read types to train a model for |
| `--long-read-bases` | `300000,1500000,6000000,30000000,150000000,1500000000:4,6000000000:2` | bases per pb and ont sample, one design point each, `DEPTH:SAMPLES` as for `--read-pairs` |
| `--long-read-samples` | 36 | samples per long-read design point, at most the communities of the paired-end points of its depth |
| `--pb-setup` | `hifi:15000:3000:3` | `hifi` : length mean : length SD : SD of the reads' quality around their length's (`hifi_reads.py`, [model-training.md](model-training.md#training-data)); or a pbsim3 setup as `--ont-setup`'s |
| `--ont-setup` | `qshmm:QSHMM-ONT-HQ:8000:6000:0.97:39/24/36` | pbsim3 method : model : length mean : length SD : accuracy (: substitution/insertion/deletion mix, for qshmm) |
| `--pbsim`, `--pbsim-models` | `pbsim`, found next to it | pbsim3 and its models (for ont, and pb with a pbsim3 setup) |
| `--test-samples` | 4 | samples per design point of the independent test set; 0 for none |
| `--test-read-pairs`, `--test-species-per-sample`, `--test-abundance`, `--test-strains-per-species`, `--test-long-read-bases` | `500,2000,10000,50000,200000,1000000,5000000:2`, `10-300`, `lognormal:2.0`, `0.5,0.2`, `150000,1000000,5000000,25000000,250000000,3000000000:2` | the test set's design |
| `--test-long-read-samples` | 8 | samples per long-read design point of the test set (with 4, the 22 samples of a long-read type could not tell its knob curve from one knob) |
| `--seed` | 1 | |
| `--ntree` | 64 | trees (more did not score better, see [model-training.md](model-training.md#training)) |
| `--maxnodes` | `512,pb:128,ont:128` | leaves per tree at most, `N` for the read types not named and `TYPE:N`: at r226 512 leaves gave the short-read models a lower log loss and fewer false positives than 256, and 128 the long-read models a lower log loss at the same F1 ([report](claude/2026-10-02-r226-v3-training/README.md)); `build_metadata.tsv` records each model's |
| `--features` | `normalized+adjacency+distance+depth+divergence+unfiltered+priors` | the models' features ([model-training.md](model-training.md#features)), feature groups joined by `+`: the normalised and gene neighbour features, the four relatives features by the references' distance (at r226 +0.004 paired-end F1 at the knob curve over `normalized+adjacency`, the other read types within noise), the sample's depth (`depth`: on the r226 v5 tables false positives halved at knob 0.5, test F1 +0.004 pe / +0.009 se / +0.005 pb over the distance set; no knob curve is fitted with it), the divergence features (`divergence`: divergence by gene conservation and codon position, lost mates; [report](claude/2026-10-03-false-positive-anatomy/README.md)), the reads before the filters and the failed candidates (`unfiltered`) and what GTDB knows of the species (`priors`, [report](claude/2026-10-03-false-positive-fixes/README.md)); `normalized+adjacency+distance+depth+divergence` leaves the last two groups out, `normalized+adjacency+distance` the depth and divergence features too, `normalized+adjacency` the relatives too, `normalized+adjacency+relatives+depth+divergence+unfiltered+priors` takes all the relatives features (a better ranking, but no better at the knob curve at r226), `normalized` leaves the gene neighbour features out too, to test them on real data; `build_metadata.tsv` records the set |
| `--call-mode` | `curve` | `curve`: the models carry a knob curve over depth where one is fitted (none with the sample's depth as a feature, the default: protal then calls at 0.5 or `--knob`); `fdr`: the models also carry calibrated calls (`random_forest_cmdline.py --fdr-calls`, [model-training.md](model-training.md#calls-at-a-target-share-of-false-calls)), which protal uses only with `--fdr F`: they called 0.0005-0.004 F1 below the curve on the benchmark world and 0.001-0.007 at r226 ([report](claude/2026-10-03-r226-v5-v6-training/README.md)). `build_metadata.tsv` records it (`classifier_call_mode`, `model_<type>_false_calls`) |
| `--evaluation` | `full` | how much the trainer evaluates: `full`, `basic` or `none` |
| `--previous-procedure` | off | the trainer also compares each model with its previous procedure ([model-training.md](model-training.md#training)), for the first builds of a release; `build_metadata.tsv` says whether it did |
| `--progress-every` | 0 | seconds between status lines of the stages running, besides each step's start and end (below); 0 for none |
| `--scratch` | | a node's own disk for the simulated samples (below) |

With the defaults that is 3 read setups x (5 depths x 12 samples + 4 samples of 2M and 2 of 10M read
pairs) = 198 paired-end samples of 20-200 species, profiled as paired-end and as single-end reads,
5 x 36 + 4 + 2 = 186 PacBio and as many Nanopore samples of the same communities, and a test set of
3 x (6 x 4 + 2) = 78 paired-end samples (and their single-end counterparts, and 5 x 8 + 2 = 42 PacBio and
Nanopore samples), all against
the training database. The simulations take most of the compute and disk space, the deep points most
of the alignment; training and its evaluation take minutes per read type. Whether the samples are
enough, the training report's learning curve says. The simulator writes a design point's samples on
threads of their own (as many as the point's share of the work), so that the deep points do not
take as long as their samples one after the other.

The finished database is needed only at the end, to take the trained models (`--add_model`), so it
is built in the background from the start, while the training database is built and the training
data are collected; the script waits for it before packing the models. That needs the memory of two
builds at once (about 50-60 GB each at GTDB scale before the genes were held at two bits per
base, which takes ~14 GB off each; estimated), or of one build and the collection's protal
runs. `--one-build-at-a-time` builds it after the training instead.

The simulations need no database. Both collections simulate from the moment the species to leave
out are chosen, in the background and at a lower priority (`nice` 10) than the builds
(`collect_training_data.py --simulate_only`; `training_data_simulation.log`,
`test_data_simulation.log`), so that they use what the builds leave of the cores. Within a
collection, the paired-end points and the long-read samples share the cores in one queue, the long
reads starting from the paired-end points' communities (a design run of seconds) rather than after
their reads, the longest samples first and those above 250 Mb in chunks side by side
(`--long_read_chunk`): at r226 (v5) the long reads had waited for the paired-end points (15 min), and
the 1.5-6 Gb samples, single runs queued last, ran alone for 27 more minutes. Once the training
database is built, the test set's samples are profiled in the training data's protal run, which
loads the database once ([report](claude/2026-10-03-build-profiling-r226/README.md)). The genome
table has each genome's length as a fourth column, so the simulator does not read every genome for
it at each design point.

Both builds read a `full_reference.fna.zst` (the marker genes of every genome; the training
database's without the species it leaves out), and nothing reads it afterwards: the script removes
each once its build is done, and a rerun that finds a build done removes one an earlier version
left (an uncompressed `full_reference.fna`, too). The converter's workers compress their chunks of
it as they write them (each a zstd frame; joined, they are the file), so the ~86 GB of r226 are
never written uncompressed, and `convert.log` gives each of the converter's steps its time; with
`--scratch` it spools the marker genes there (`--tmp`). The training database's copy is written with
`-t` threads of zstd (up to 8). The finished database is compressed at `--final-db-level` (9), the
training database at `--training-db-level` (3).

### What it prints

Each stage writes its output to its own log (listed below). On the console, every line starts
with the time and how long the run has taken so far (`[14:03:22 +1:25:00]`). The run has eight
steps (seven without a test set): the genome table, the release, the databases, the training data,
the test set, the models, the parity check and packing the models. A step says when it starts,
numbered and with its log (`4/8 training data (training_data.log): ...`), and on indented lines how
it went: how long it took, the peak memory of the largest process it ran, and a few numbers of what
it made (for a build the size of `database.protal`; for a collection the present and absent taxa of
each read type's table; for the models their F1 with species held out and on the test set). The
finished database's build in the background gets its indented line when it ends, whatever step the
run is at.

With `--progress-every N` each stage still running also gets a status line every N seconds: how
long it has run, the memory it and the commands it started use now (Linux only), and the last line
of its log. The collector's log (`training_data.log`, `test_data.log`), which such a line shows,
says when each design point is simulated, how many are done, and which sample protal is aligning,
then how many it has profiled. All this reads `/proc` and the end of the logs once per status line
(and with `--scratch` the free space of its file system every 5 s), which costs nothing next to
the stages.

From a small build of the tuning world (765 species, one sample per design point, two depths,
`--scratch SCRATCH`; the summary table at the end left out):

```
[02:48:16 +0:00:01] A protal database of GTDB r226 in OUT, with 4 threads; each step logs to a file there; the simulated samples go to SCRATCH (883.6 GB free)
[02:48:18 +0:00:02] 1/8 genome table (genomes.tsv, genome_table.txt): 2295 genomes of 765 species (Archaea 61, Bacteria 704); a simulated species is another genome than its representative 66.7% of the time
[02:48:18 +0:00:02] 2/8 converting GTDB r226 (convert.log)
[02:48:21 +0:00:05]     converted in 0:00:01, peak memory 174 MB; the genes' neighbours found in 0:00:02, peak memory 203 MB (gene_neighbours.log)
[02:48:21 +0:00:05] 3/8 training database (training_db, training_db_index.log): 296 species left out, 2 phylum, 4 class, 6 order, 8 family, 12 genus clades (179 species) and 117 species alone (model_logs/holdout.txt)
[02:48:21 +0:00:06]     simulating the training data and the independent test set meanwhile, in the background (training_data_simulation.log, test_data_simulation.log)
[02:48:21 +0:00:06]     building protal_db meanwhile, in the background (index_and_package.log)
[02:48:35 +0:00:20]     built training_db in 0:00:14, peak memory 4.1 GB; database.protal 83 MB; full_reference.fna.zst removed (14 MB); its files written in 0:00:01, peak memory 174 MB (training_db.log)
[02:48:35 +0:00:20] 4/8 training data (training_data.log): 2 pe, 2 se, 2 pb, 2 ont samples, 1 per design point
[02:48:35 +0:00:20]     waiting for its simulations in the background (training_data_simulation.log)
[02:48:41 +0:00:25]     simulated the training data in the background in 0:00:19, peak memory 41 MB
[02:48:44 +0:00:29]     collected in 0:00:03, peak memory 3.6 GB; taxa present/absent: pe 139/108, se 139/87, pb 79/7, ont 110/20; 25 MB in SCRATCH/training; scratch: 404 MB in use by the run (at most 404 MB), 883.2 GB free
[02:48:44 +0:00:29] 5/8 independent test set (test_data.log): 1 pe, 1 se, 1 pb, 1 ont samples, 1 per design point
[02:48:44 +0:00:29]     waiting for its simulations in the background (test_data_simulation.log)
[02:48:45 +0:00:30]     simulated the independent test set in the background in 0:00:24, peak memory 35 MB
[02:48:46 +0:00:31]     collected in 0:00:02, peak memory 3.6 GB; taxa present/absent: pe 140/37, se 130/27, pb 31/5, ont 47/6; 6 MB in SCRATCH/test; scratch: 307 MB in use by the run (at most 404 MB), 883.3 GB free
[02:48:46 +0:00:31] 6/8 training the pe, se, pb, ont models (classifier_training*.log) in parallel, 1 thread each
[02:48:48 +0:00:33]     trained in 0:00:02 (pe 0:00:02, se 0:00:02, pb 0:00:02, ont 0:00:02); F1 with species held out/on the test set: pe 0.986/0.771, se 0.993/0.860, pb 0.981/0.951, ont 0.995/0.894
[02:48:48 +0:00:33] 7/8 checking that protal scores the models as the trainer does (parity*.log)
[02:48:52 +0:00:37]     pe, se, pb, ont: the same probabilities and features, checked in 0:00:04
[02:48:52 +0:00:37]     the held-out species' reads by gene conservation, in 0:00:00: model_logs/relatives_by_gene_conservation.txt
[02:48:52 +0:00:37] 8/8 adding the pe, se, pb, ont models to protal_db (final_package.log)
[02:48:52 +0:00:37]     waiting for protal_db's build in the background (index_and_package.log)
[02:49:56 +0:01:41]     built protal_db in the background in 0:01:35, peak memory 4.6 GB; database.protal 108 MB; full_reference.fna.zst removed (23 MB)
[02:49:57 +0:01:42]     added in 0:00:01; database.protal 108 MB
...
[02:49:57 +0:01:42] The run took at most 489 MB on SCRATCH; the simulated samples there (38 MB) are left for a rerun
[02:49:57 +0:01:42] Ready protal database: OUT/protal_db; database.protal 108 MB; model evaluation: OUT/model_logs (start with trained_model.report.txt, and trained_model_<read type>.report.txt)
```

(The scratch figures count everything written to that file system during the run, other jobs'
too.)

The conversion's line now also says how the genes were placed for the gene neighbours, which takes
longer since every genome counts, e.g. on an operon-like world of the same 2,295 genomes:
`converted in 0:00:03, peak memory 177 MB; the genes placed in 2295 genomes (86236 exactly, 168006 by
their k-mer trace) and their neighbours counted in 0:00:23, peak memory 176 MB (gene_neighbours.log)`
([report](claude/2026-10-02-gene-neighbour-frequencies/README.md)).

### Local scratch

The collection makes many files and deletes them again: the simulator writes a plain copy of
each genome of a sample and ART's reads of each genome before it compresses the sample, and for
long reads the collector writes each sample's reads as templates for pbsim3 (until 2026-10-02 it
copied every genome of every sample for pbsim3, which wrote three files per contig). With the
defaults that is about 240 GB of temporary writes in ~150,000 files (estimated; 470 GB in millions
of files before), which a network file system is slow at (the collector then waits in state D).
`--scratch DIR`
makes the training and test samples on a node's own disk, in `DIR/training` and `DIR/test`, and
copies their tables to `OUTDIR/training` and `OUTDIR/test`; the converter spools the release's
marker genes there too, before the samples (~35 GB at r226, estimated). The simulator compresses a
sample's reads as it appends each genome's, and the collector removes a long-read sample's
templates once its reads are written, so the samples' temporary files stay small. With the 0.7.2
design the r226 build took at most 23.7 GB there and left 15.5 GB; with the deep design points it took at
most 120 GB (24 long-read samples per point; [docs/claude/2026-10-01-scratch-space](claude/2026-10-01-scratch-space/README.md),
[2026-10-02-r226-v3-training](claude/2026-10-02-r226-v3-training/README.md)), and the 36 long-read samples of the
defaults add a few GB. The training database is built there too (`DIR/training_db`, ~22 GB at r226): only the
collections and the parity check read it, and writing its `database.protal` to the network file system was 6 of
the 15 minutes of its r226 build. Give it 175 GB, or fewer deep samples (`DEPTH:SAMPLES`). The
line after each collection (and the status lines, with `--progress-every`) says how much the run
takes on DIR (what its file system holds more than at the start), and the last line the most it
took. A rerun reuses the samples only from the same DIR; on a node-local disk that is gone after
the job, the collection starts again (the builds are still kept in OUTDIR).

### Stopping and rerunning

The script checks before it starts that protal, the simulator, `art_illumina` (and pbsim3 for ont,
and pb with a pbsim3 setup) are there and that its Python can import what the trainer needs. protal and the
simulator must be built from the source the script is at: of its version and, as their `--version` says since
0.7.3 (`protal v0.7.3 (commit ...)`), of its commit, with nothing changed since in `src/`, `lib/` or the build
files; otherwise it stops at once with what to rebuild (`--no-binary-check` runs them anyway). A binary built
outside a git checkout, or of uncommitted changes, is only noted. (A run of 2026-10-02 collected its training data
with an older protal and failed two hours in, when the trainer found its features missing.) Each command it runs is
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
| `genomes.tsv`, `genome_table.txt` | the genome table used for the simulations (accession, GTDB taxonomy, FASTA path, genome length; a rerun counts only the genomes changed since), and what it holds (species by domain, how often a simulated species is not its representative) |
| `training_db/`, `heldout_species.txt` | the training database and the species it leaves out (species, the rank they were held out at, the clade); the `full_reference.fna.zst` of both databases is removed once their builds are done |
| `training/`, `test/` | the simulated samples, their profiles and one table per read type (`training_data.tsv` for pe, `training_data_se.tsv`, `_pb`, `_ont`); a rerun reuses the design points simulated and profiled from the same inputs (below). With `--scratch`, only the tables; the samples are in the scratch folder |
| `.stages/` | the inputs of the conversion and the two builds that completed, for a rerun (below) | 
| `trained_model.*`, `trained_model_se.*`, `_pb.*`, `_ont.*` | the models and the trainer's outputs ([model-training.md](model-training.md#training)) |
| `model_logs/` | what tells whether the models are good, in one folder: `summary.txt` (per read type: TP, FP, TN, FN, sensitivity, specificity, precision, F1 and false positives per sample, with species held out and on the test set; also printed at the end), each read type's training report and its numbers (`trained_model*.report.txt`, `.metrics.json`), per-taxon predictions (also on the test set), the threshold table, feature importances, the parity checks with protal (`parity*.txt`), `genome_table.txt`, what the training database leaves out (`holdout.txt`, `heldout_species.txt`), the collection logs and `build_metadata.tsv`; and what the conservation features rest on, on the release's real genomes: `gene_congeners.tsv` and `gene_incongruence.tsv` (the finished database's, see above) and `relatives_by_gene_conservation.txt` and `.tsv` (`scripts/trace_relatives.py`: the paired-end reads of the species held out of the training database followed to the genes they align to, per unit coverage against a species' own reads, before and after the MAPQ filter, by the genes' factors; a failure there is reported and does not stop the build) |
| `*.log` | one log per stage: `convert`, `training_db`, `index_and_package`, `training_db_index`, `training_data_simulation`, `training_data`, `test_data_simulation`, `test_data`, `classifier_training*`, `parity*`, `final_package` |
