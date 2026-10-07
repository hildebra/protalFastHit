# Databases: files, building and the presence model

A protal database holds the marker genes of GTDB's species, a k-mer index of them, and a presence
model per read type that decides which species are in a sample. The database on the
[downloads page](https://protal.earlham.ac.uk/main.php?site=downloads) is built from GTDB r226.

| Part | What it covers |
|---|---|
| [The files of a database](#the-files-of-a-database) | what a database holds, `database.protal`, compression, converting and unpacking, memory |
| [Building a database](#building-a-database) | from GTDB in one command, for several releases, reduced marker sets, step by step |
| [The presence model](#the-presence-model) | how protal calls species, training data, training and checking a model, installing it |

The model's features are listed in [features.md](features.md); simulated test worlds and mini
databases are in [development.md](development.md).

## The files of a database

### What a database holds

| File | What it holds |
|---|---|
| `index.prx` | the k-mer index of the reference marker genes |
| `reference.fna`, `reference.map` | the reference genes (`>taxid_geneid`), and each gene's taxid, gene id and byte range |
| `internal_taxonomy.dmp` | the taxonomy: id, parent, external id, name, rank, level, representative genome |
| `unique_kmers.tsv` | per species and gene, the k-mers unique to it in the database |
| `gene_conservation.tsv` | optional: per gene, how fast it diverges within species against the species' other genes ([below](#2-build-the-index)); for the conservation and divergence features, and `--gene_conservation db` |
| `suspect_copies.tsv` | optional: gene copies near-identical to another genus's (contamination, transferred genes), whose reads a run leaves out of the evidence ([below](#2-build-the-index)) |
| `species_neighbours.tsv` | optional (since 2026-10-06): each species' nearest congeners in the database by the distance of their marker genes, for the database-neighbourhood features and `unexpected_congener_fit_share` ([below](#2-build-the-index), [features.md](features.md#against-false-positives-in-complex-communities-consistency-shape-neighbourhood-2026-10-06)) |
| `species_priors.tsv` | optional: what GTDB knows of each species before any read ([below](#species-priors)) |
| `gene_neighbours.tsv`, `gene_positions.tsv` | optional: which genes lie next to which, per clade, and where each gene lies in each genome ([below](#gene-neighbours)); a run loads only the first |
| `gene_table.bin` | only in `database.protal`: `reference.map` and `unique_kmers.tsv` in binary, loaded on all threads without parsing (r226-sized tables, six threads: 1.76 → 0.43 s) |
| `model_pe.xml`, `model_se.xml`, `model_PB.xml`, `model_ONT.xml` | the presence models of paired-end, single-end, PacBio and Nanopore reads ([below](#the-presence-model)); `model.xml` in databases before 0.7 |

A database without an optional table works; the features that need it are then fixed values (0, or
0.5). A model trained with a table should run with it.

Every file is checked when a run starts, and protal stops at the first problem, with the file and
line. Ids and names in the taxonomy are unique, every lineage ends at the root, and the tables name
only genes and clades the database has.

### One file: database.protal

`protal --build` packs the files into `database.protal`, compressed with
[zstd](https://facebook.github.io/zstd/): one file to copy, download, checksum or version, whose
parts cannot get out of step. `--db` (or `$PROTAL_DB_PATH`) takes that file, a folder that holds it,
or a folder of separate files (raw, or `index.prx.zst` and `reference.fna.zst`); a folder with both
uses the separate files. Databases of earlier versions (raw files, index format 1) load unchanged.
protal prints the index's features when it loads it (`Index features: ...`) and checks the index
against the reference.

### Compression

Compression saves disk space and makes loading faster wherever storage is slower than decompression,
as on network file systems. The mini database is 1.0 MB as `database.protal` and 3.2 GB raw (the
index's key map has a fixed size).

protal writes zstd's *seekable* format: independent frames of 64 MB, so loading decompresses frames
on `-t` threads (~1-1.5 GB/s each). In `database.protal` each part is a range of frames. The index is
stored in columns (per chunk, the number of values per k-mer, then each field in byte planes only as
wide as the chunk needs), 25-55% smaller than its raw bytes compressed alike and as fast to load.
Because of the columns, plain `zstd -d` does not give an `index.prx`; use `protal --unpack_db` or
`--decompress_db`. `scripts/db_compression_benchmark.sh` measures ratio and speed per level on your
database.

### Build options for the format

`protal --build` writes `database.protal` into the `--db` folder, reads it back, compares, and
removes the separate files it holds. The build inputs (`full_reference.fna`, ...) stay.

| Build option | Default | |
|---|---|---|
| `--no_bundle` | off | separate compressed files (`index.prx.zst`, `reference.fna.zst`, ...) |
| `--no_compress` | off | separate raw files |
| `--compress_level` | 19 | zstd level 1-22; 19 compresses ~3 MB/s per thread, 12 ~40 MB/s; decompression speed hardly depends on it |
| `--compress_frame_mb` | 64 | frame size; 0 writes a single frame (one loading thread; needs `--no_bundle`) |
| `--compress_window_log` | 27 | long-distance matching window (128 MB); 0 turns it off |

### Converting a database

A database is converted in place without rebuilding it; the content stays byte-identical and is
verified before the old files are removed:

```bash
protal --compress_db --db /path/to/protal-db -t 16      # separate files -> database.protal
protal --unpack_db --db /path/to/protal-db -t 16        # database.protal -> separate files beside it (it is kept)
protal --decompress_db --db /path/to/protal-db -t 16    # -> separate raw files
```

- `--unpack_db` writes `index.prx.zst`, an uncompressed `reference.fna` and the other files but
  `gene_table.bin` (beside `database.protal`, or into `--unpack_dir`).
- `--compress_db` on a `database.protal` without a current `gene_table.bin` (packed before 0.7.6)
  adds it, rewriting the file once with the other parts' frames copied.
- `--compress_db` on a folder of raw files writes the same `database.protal` that `--build` would.

**Models.** `protal --add_model MODEL --read_type pe --db DB` (or `se`, `pb`, `ont`) checks a model
and replaces the database's model of that read type. Several at once
(`--add_model pe.xml,se.xml --read_type pe,se`) are all checked before anything is written.
`--model` and `--model_se` use another model for a run without changing the database.

The models are the last members of `database.protal`, so `--add_model` replaces them in place: it
rewrites only the end of the file (the models and the seek table) and the directory frame at its
start, a few MB, and leaves the rest untouched. That takes seconds, mostly compressing the models
(about 8 s for a 12 MB model at level 19; the models are compressed in parallel). A
`database.protal` packed before this change, or one that gains a read type it had no model for, is
rewritten once instead, with the models put last. That copies the whole file (21 GB at r226; 2.5
minutes on a local SSD, ~20 minutes on a network file system). The GTDB build's final database
already holds a model or placeholder for every read type (unless `--no-placeholder-models`), so its
`--add_model` step runs in place. While the models are replaced, `database.protal.journal` holds the
bytes being overwritten. If `--add_model` stops halfway, running it again first puts those bytes
back; other protal commands say so instead of reading a broken file. Do not start protal runs on the
database while `--add_model` replaces its models.

### The database in memory

- **Genes** are held at two bits per base. `N` is stored as `A`, an IUPAC code as the first base it
  stands for (A C G T order), so alignments, SNP and MSA reference rows show the stored base. `--build`
  leaves k-mers with an ambiguous base out of the index. A SAM of an older protal with an `M` against
  such a base is set aside (`misc/<sample>.err`).
- **The index** is held packed: each value takes the bits the database needs, its taxon and gene as
  one number (41 instead of 64 at r226: 26.6 GB instead of 35 GB), by a run and, since 2026-10-05, by
  `--build` too. The files are unchanged.
  A run and a build print `Index in memory: ...`.
- A full r226 run peaked at 36.5 GB, about 34 GB expected now ([running.md](running.md#memory)).
- `--preload_genomes_off` reads genes on demand from a raw `reference.fna`: less memory, slower, and
  it needs the database as separate files (protal prints the `--unpack_db` command).

`reference.fna` is ordered by gene, then taxon (`gtdb_to_protal_db.py --order gene`), which compresses
about twice as well as genome by genome; protal finds genes through `reference.map` in any order.

## Building a database

You can build a database from any GTDB release from r207 on, with all marker genes or a subset of
them:

| Route | Use it for | Section |
|---|---|---|
| **One command**: download, then build and train | a database to profile with (recommended) | [Build and train in one command](#build-and-train-in-one-command) |
| **Several releases** | each release's full and reduced database, with a summary table | [Several releases](#several-releases-full-and-reduced) |
| **Step by step** | your own pipeline, or a database without training | [Step by step](#step-by-step) |

The mini databases used for testing are made the same way from a synthetic release
([development.md](development.md#mini-databases)).

### What you need

- `protal` and `simulate_metagenomes` built from this checkout. The build script checks that they
  match the scripts' version and commit ([installation.md](installation.md#building-from-source)).
- pbsim3's model files for Nanopore reads (`QSHMM-ONT-HQ.model`: install pbsim3 or pass
  `--pbsim-models`; pbsim3 itself is not run).
- Python 3 with numpy, pandas, joblib and scikit-learn for the training.
- NCBI's `datasets` CLI for the download.

The conda environment [`envs/protal-db-build.yaml`](../envs/protal-db-build.yaml) has all of them,
with the compilers to build protal ([installation.md](installation.md#tools-to-build-a-database)).

At GTDB r226 (143,614 species) the whole pipeline takes:

| | |
|---|---|
| download | 17.7 GB of GTDB files and about 185 GB of genomes to simulate from (estimated for the defaults since 2026-10-05; 75 GB with the earlier 8,000 species), on a node with internet |
| time | about 2.5-3 hours on a 52-64-thread node (conversion 5 min; the r226 build of 2026-10-03 then took 2 h, and 0.7.6's design simulates half as many short-read samples again), and several hours more for the [scenarios](#scenarios-kinds-of-studies) (an estimate; `--scenarios none` leaves them out); training gradient-boosted models (the default since 2026-10-06) takes a few minutes per model with the build's defaults (an estimate; `--features auto` and `--evaluation full` take hours) |
| memory | each index build ~35-37 GB with 64 threads (estimated; it was 64 GB before 2026-10-05; its log's `Memory after ...` lines say), up to ~75 GB while both run at once; `--one-build-at-a-time` needs about half |
| node-local disk | for the simulated samples (`--scratch`): with `--profile-blocks` (the default) each design point's reads are removed once profiled, and the simulations wait rather than leave less than `--keep-free` GB, so ~150-200 GB is enough; all at once (`--profile-blocks 0`) ~180 GB for the design and ~500 GB more for the default scenarios (6 + 3 samples each since 2026-10-06; ~260 GB for the 3 + 2 measured on 2026-10-05) |
| result | `database.protal`, ~27 GB |

### Build and train in one command

Two commands: the download needs the internet, the build a compute node.

```bash
# 1. on a node with internet: GTDB's files and the genomes to simulate from
python3 scripts/download_gtdb.py -o /shared/protal_inputs/gtdb_r226

# 2. on a compute node: convert, build, simulate, train, package
python3 scripts/build_gtdb_database.py --inputs /shared/protal_inputs/gtdb_r226 --outdir /data/protal_r226 \
    --protal build/protal --simulator build/simulate_metagenomes -t 64 --scratch /local/scratch
```

The result is `/data/protal_r226/protal_db/database.protal`, with a model for each read type
(paired-end, single-end, PacBio HiFi, Nanopore). Profile with `protal --db /data/protal_r226/protal_db`.

#### 1. Download

`download_gtdb.py` fetches the release (r226 by default; `--release 220`, or a point release such as
`214.1`):

- **GTDB's files**: taxonomy, metadata, the species clusters and the marker genes of the
  representatives and of all genomes. They are checked against the release's `MD5SUM.txt` and
  extracted.
- **Genomes to simulate from**, from NCBI: 16,000 species with up to 2 non-representative genomes
  each, and 9,000 species with their representative only, about 50,000 genomes (estimated; until
  2026-10-05 6,000 and 2,000). Species are picked per domain in GTDB's proportions, in a random
  order fixed by the seed. 25,000 species hold the soil scenarios at full size, and most simulated
  strains are real ones: with 6,000 and 19,000 (r226 v12) three quarters were in-silico, and the
  models missed real strains 2.4 times as often ([report](claude/2026-10-05-r226-v12-scenarios/README.md)).
- **The host genome** of the [host scenario](#scenarios-kinds-of-studies): the human T2T-CHM13v2.0
  assembly (`GCF_009914755.1`), kept gzipped as NCBI serves it (0.9 GB, 3.1 GB unpacked). One that
  cannot be fetched is only noted: the other builds do not need it.

The folder serves every later build of that release. A rerun fetches only what is missing, and an
interrupted download resumes. A rerun with other counts keeps what the folder has as far as they
allow: a larger `--species` gives strains to its representative-only species first, so only their
strains are fetched. `--dry_run` lists what it would fetch.

Strains are chosen by quality among those passing the CheckM2 filters, in this order:
1. isolates before single-cell genomes and MAGs;
2. complete genomes before chromosome-, scaffold- and contig-level assemblies;
3. PacBio or Nanopore assemblies before others (NCBI is asked for the best 30 candidates per
   species);
4. fewer contigs.

Ties are drawn at random. The representatives are GTDB's, whatever their quality.

Genomes come from NCBI's FTP server, `--connections` (8) at a time, each checked against its length
and gzip. What that cannot deliver (renamed or withdrawn assemblies) goes through NCBI's `datasets`
CLI. The log shows the rate. Raise `--connections` on a fast link, and lower it if NCBI answers with
HTTP 503.

| `download_gtdb.py` writes | |
|---|---|
| `release/` | GTDB's files as `--gtdb` reads them (archives removed once extracted; `--keep_archives` keeps them) |
| `genomes/` | `<accession>.fna.gz`. A rerun that chooses other genomes moves the ones no longer wanted to `genomes_unused/` |
| `genomes.tsv`, `missing.txt` | the genomes (species, role, lineage, CheckM2, category, assembly level, contigs, sequencing technology), and those NCBI did not deliver |
| `simulation_species.txt` | the species to simulate from |
| `host/` | `<accession>.fna.gz`, the host genome |
| `ncbi_info.tsv`, `download.json` | NCBI's sequencing technologies; the release, options, counts and checksums |

| Option | Default | |
|---|---|---|
| `-o` | required | the folder |
| `--release` | 226 | GTDB release (207 or later) or point release |
| `--species`, `--per_species`, `--rep_only_species` | 16000, 2, 9000 | species with strains, strains each, species simulated from their representative only (and an in-silico strain) |
| `--min_completeness`, `--max_contamination` | 90, 5 | CheckM2 filters for strains |
| `--rep_genomes` | `ncbi` | where the representatives' genomes come from; `gtdb`: GTDB's archive of all of them (137 GB for r226) |
| `--tech_candidates`, `--no_tech_lookup` | 30 | candidates per species whose sequencing technology is asked of NCBI; or ask for none |
| `--progenomes` | | a proGenomes ANI-clustering table ([`pg4_ANI_clustering.tsv.gz`](https://progenomes.embl.de/download.cgi)) or any list of accessions: strains only from these |
| `--no_genomes`, `--dry_run` | | GTDB's files only; list what would be fetched |
| `--connections`, `--ftp_url` | 8, NCBI's | genomes fetched at once from NCBI's FTP server (0: only through `datasets`) |
| `--host_genome` | `human` | the host genome: `human`, `ACCESSION_ASSEMBLYNAME` of another NCBI assembly (e.g. `GCF_000001405.40_GRCh38.p14`), or `none` |
| `--mirror`, `--datasets`, `--batch`, `-t` | | GTDB server, NCBI CLI, genomes per `datasets` request (500), compression threads (8) |

#### 2. Build and train

`build_gtdb_database.py` runs these steps. Each logs to its own file in `--outdir` and prints one
numbered line when it starts and indented lines when it ends.

1. **Genome table** (`genomes.tsv`, `genome_table.txt`): the genomes to simulate from and their
   lengths.
2. **Conversion** (`convert.log`): the release into a database folder. The genes are then placed in
   every genome to record which genes are neighbours (`gene_neighbours.log`,
   [details](#gene-neighbours)).
3. **In-silico strains** (`insilico_strains.log`): a mutated copy of the representative of every
   species with one genome, so that these species are simulated from a strain too
   ([details](#training-data-like-real-samples)).
4. **Training database** (`training_db_index.log`): the database without the species and clades held
   out (30% of the species, plus whole clades of every rank). Simulating the samples and building the
   finished database start in the background here.
5. **Training data** (`training_data.log`): 309 paired-end samples of 20-200 species each (3 read setups x 103), profiled
   as paired-end and as single-end reads, and 186 PacBio and 186 Nanopore samples of the same
   communities, and the hold-in samples of four scenarios of real studies (gut, soil, shallow soil,
   host-dominated; [below](#scenarios-kinds-of-studies)), all against the training database.
6. **Independent test set** (`test_data.log`): samples of another design, for an honest score, and
   the scenarios' hold-out samples.
7. **Models** (`classifier_training*.log`): one model of gradient-boosted trees per read type (`--model`;
   a random forest before 2026-10-06), trained in parallel on the default feature set and evaluated with
   rows, samples and species held out (`--features`, `--evaluation basic`, the build's defaults since
   2026-10-06); with `--features auto` each trainer chooses its set and the console says which won and why.
8. **Parity check** (`parity*.log`): protal scores each model exactly as the trainer does.
9. **Packaging** (`final_package.log`): the models go into the finished database.

A rerun into the same `--outdir` resumes. Completed steps are skipped when their inputs are
unchanged, and simulated samples are reused when their design is the same
([details](#stopping-and-rerunning)).

#### 3. Check the result

`OUTDIR/model_logs/summary.txt` (also printed at the end) gives each model's sensitivity, precision,
F1 and false positives per sample, with species held out in training and on the independent test
set. The r226 build of 2026-10-04 scored F1 0.963 (paired-end), 0.960 (single-end), 0.968 (PacBio)
and 0.969 (Nanopore) on its test set ([report](claude/2026-10-04-r226-v10-evaluation/README.md)).
Each model's full report is `trained_model*.report.txt`. Start with these sections:

- **Independent test set**: the score per depth, and the threshold with the highest F1. A test score
  clearly below the cross-validated one is flagged.
- **False positives and false negatives by taxonomic rank**: how the model handles species and
  clades the database lacks.
- **Strains**: how often strains are missed, for real and for in-silico strains.
- **Scenarios**: F1, FP and FN rates per scenario, on its hold-out and hold-in samples;
  `summary.txt` has a row for each.
- **Feature set chosen**: each candidate set's F1 with species held out and on every test set (the
  design's and each scenario's hold-out); `summary.txt` and the console give the winner and why.

[The presence model](#the-presence-model) explains the reports; [features.md](features.md) lists
what the models use.

#### What it prints

Every console line starts with the time and how long the run has taken. The first lines of the
r226 build of 2026-10-05:

```
[07:53:00 +0:00:44] A protal database of GTDB r226 in /hpc-home/hildebra/DB/protal/protal0.7.5_r226_v11, with 64 threads; each step logs to a file there; the simulated samples go to /nbi/local/ssd/23954360/protalDBbuild_1 (263.1 GB free)
[07:55:58 +0:03:41] 1/9 genome table (genomes.tsv, genome_table.txt): 17248 genomes of 7998 species (Archaea 388, Bacteria 7610); a simulated species is another genome than its representative 44.2% of the time
[07:55:58 +0:03:41] 2/9 converting GTDB r226 (convert.log)
[08:00:47 +0:08:31]     converted in 0:03:25, peak memory 2.7 GB; the genes placed in 17180 genomes (1261838 exactly, 563962 by their k-mer trace) and their neighbours counted in 0:01:25, peak memory 3.6 GB (gene_neighbours.log)
[08:00:47 +0:08:31] 3/9 in-silico strains of the species with one genome (insilico_strains.log, genomes_simulated.tsv)
```

An indented line says how long a step took, the peak memory of its largest process, and what it
made: a database's size, a collection's present and absent taxa per read type, a model's F1. The
background build gets its line when it ends, whatever step the run is at. With `--progress-every N`
each running stage also prints a status line every N seconds (its run time, memory and the last line
of its log). The run ends with `Ready protal database: ...` and the path of the model reports.

#### Outputs

| Path in `--outdir` | |
|---|---|
| `protal_db/database.protal` | the finished database with its models; `protal_db/build_metadata.tsv` records the release, protal version and commit, command, design, held-out species and each model's scores. The versions and the scripts' commit are those the run started with (since 2026-10-06; read at its end before), followed by "at the end of the run: ..." if a pull or a rebuild changed them meanwhile, which the console then warns of |
| `model_logs/` | everything to judge the models: `summary.txt`, each read type's report (`trained_model*.report.txt`, `.metrics.json`), predictions (`.scenario_predictions.tsv.gz`: the scenarios' hold-out samples; `.calls.tsv.gz`: every row's call), thresholds, feature importances, parity checks, `genome_table.txt`, `holdout.txt`, `build_metadata.tsv`; what the conservation features rest on: `gene_congeners.tsv`, `gene_incongruence.tsv`, `relatives_by_gene_conservation.txt`; and `error_reads/`, where the reads behind each model's errors in every sample went (with `--share-logs` their SAM records too; [below](#the-reads-behind-the-errors)) |
| `<name>_share.tar.gz` | with `--share-logs`: the logs, `model_logs/` and the tables, to copy off the cluster ([below](#logs-to-share)) |
| `trained_model*` | the models and the trainer's outputs ([the presence model](#training)) |
| `genomes.tsv`, `genome_table.txt` | the genomes simulated from (accession, taxonomy, FASTA, length), and a summary |
| `genomes_simulated.tsv`, `insilico_strains/` | the same with the in-silico strains, their FASTAs and `insilico_strains.tsv` (per strain: divergence drawn and reached, substitutions) |
| `heldout_species.txt` | the species the training database leaves out, with the rank they were held out at and the clade |
| `training/`, `test/`, `training_db/` | one table per read type (`training_data.tsv` for pe, `_se`, `_pb`, `_ont`); without `--scratch` also the samples, their profiles and the training database |
| `.stages/`, `*.log` | what a rerun may skip; one log per step, and `console.log`, the console's lines |

### Several releases, full and reduced

Two scripts run the pipeline for a list of GTDB releases, each as a full database and a reduced
one of 12 marker genes ([reduced marker sets](#reduced-marker-sets)). Only the first needs the
internet:

```bash
# on a node with internet: INPUTS/gtdb_r220, INPUTS/gtdb_r226
python3 scripts/download_gtdb_releases.py -o /shared/protal_inputs --releases 220,226 -t 8

# on a compute node: OUT/r220_full, OUT/r220_n12, OUT/r226_full, OUT/r226_n12
python3 scripts/build_gtdb_releases.py --inputs /shared/protal_inputs --outdir /data/protal_dbs \
    --protal build/protal --simulator build/simulate_metagenomes -t 64 --scratch /local/scratch
```

`download_gtdb_releases.py` runs `download_gtdb.py` for each release. Options it does not know pass
through (e.g. `--species 6000`). It writes `INPUTS/download_summary.tsv`: per release the taxa per
rank and domain, the genomes delivered and missing, the species to simulate from, the size on disk
and the time. A release that fails is reported, and the others go on.

`build_gtdb_releases.py` runs `build_gtdb_database.py` for each release (`--releases`, default all
in `INPUTS`) and variant (`--variants`, default `full,n12`; `--n-genes` sets the 12):
- The full database is built first, with `--rank-genes`.
- The reduced one is built from that ranking, so it needs no extra ranking build.

Its other options pass through (`--samples`, `--read-types`, `--features`, ...). With `--scratch`
each database gets its own folder there. Databases already built and trained are kept (`--rerun`
builds them again). At the end `OUT/build_summary.tsv` and `.txt` compare the databases: release,
protal version, genes, taxa per rank and domain, size, build time and peak memory, profiling memory,
and per read type the test-set F1, false positives per sample, sensitivity and precision.

### How the build works

#### Training data like real samples

The model has to learn what reads look like when they do not come from the database's own
references. On a GTDB-like simulated world it needed two things
([report](claude/2026-09-29-model-training-tuning/README.md)). Without them it called a quarter of
the relatives of missing species present. With only one, it either kept 2 false positives per
sample or missed a fifth of the strains.

- **Other strains.** Simulated only from the representatives, every species would be the database's
  own reference, closer to it than real strains are. With the downloaded pool and the in-silico
  strains below, a simulated species is another genome than its representative more than half the
  time, mostly a real strain; `genome_table.txt` and the console say how often, real and in-silico
  apart, and the build warns when in-silico strains outnumber real ones (the models miss real
  strains more often: download more species with strains, `download_gtdb.py --species`).
- **In-silico strains** (since 0.7.6). A species with one genome would always be simulated from its
  reference. A model given GTDB's cluster sizes (`+priors`) learned from that to reject a divergent
  read cloud on a one-genome species, which is a rule of the simulation, not of nature. So each such
  species gets a strain (`scripts/insilico_strains.py`; `--insilico-strains`, the share given one,
  default 1):
  - a copy of the representative with substitutions only;
  - as far from it as the table's real strains are from theirs (marker genes from the k-mer traces
    in `gene_positions.tsv`), each marker gene by its conservation factor;
  - the genome at that divergence / 0.45, at most 5% (95% ANI; `--insilico-ani MIN-MAX` draws the
    ANI instead);
  - codon-aware: amino-acid changes kept with probability 0.15, stop codons never made, so most
    changes fall on third codon positions, as a strain's do.

  The training table marks taxa simulated from one (`meta_insilico_strain`), and the report lists
  them apart.
- **Species the database lacks.** The samples are profiled against a training database that leaves
  species out. Their reads land on relatives, as reads of organisms GTDB lacks do in real samples.
  Two kinds are held out:
  - whole clades (`--holdout-clades`: 2 phyla, 4 classes, 6 orders, 8 families, 12 genera, each at
    most 2% of the species);
  - then 30% of the remaining species to simulate (`--holdout`).

  Every sample holds one species of a held-out clade of each rank. `heldout_species.txt` lists them,
  and the report gives false positive and false negative rates by rank
  ([below](#species-and-clades-the-database-lacks)). The finished
  database has every species.

The samples are as complex as real ones:
- 20-200 species each (gut samples hold 100-300 GTDB species, most of them rare);
- a second and third strain of a species in 30% and 10% of cases;
- about a quarter of the species in groups of 2-5 congeners (`--congeners`), which the relatives
  features need;
- abundances lognormal with sigma 1.3 and 2.0;
- 100 bp (HiSeq 2000), 150 bp (HiSeq X Ten) and 250 bp (MiSeq) reads, made by
  `simulate_metagenomes`'s own model of each instrument (no ART; `NovaSeq` and `HS25` are there too,
  [below](#illumina-reads)).
- depths from 1,000 to 30M read pairs (`--read-pairs`). A model with the sample's depth as a feature
  cannot extrapolate beyond its deepest training sample, and real samples are 5-50M pairs.

#### One model per read type

`--read-types` (default `pe,se,pb,ont`) trains a model for each kind of reads, from the same
communities:

| read type | how the reads are made |
|---|---|
| `pe` | Illumina paired-end reads by `simulate_metagenomes`'s instrument models ([below](#illumina-reads)) |
| `se` | the same reads' first mates, profiled as single-end |
| `pb` | PacBio HiFi reads by `hifi_reads.py`'s model (`--pb-setup`: 15 kb, quality by length, Q50 at 5 kb to Q20 at 50 kb) |
| `ont` | Nanopore reads by [pbsim3](https://github.com/yukiteruono/pbsim3)'s quality-score model (`--ont-setup`: the high-quality model at 97%, 8 kb) |

The long-read samples (`--long-read-bases`, 300 kb to 6 Gb; `--long-read-samples` 36 per point) replay
the paired-end communities of the same depth. `simulate_metagenomes --long_samples` makes them in
process, one run per design point: it draws the reads' templates and writes each read as it is made,
compressed, with no template or pbsim3 files on the disk. Once the training database is built, the samples are
profiled as they are simulated, in protal runs of at least `--profile-blocks` GB of reads (20; the
training data's and the test set's runs take turns), and each design point's reads are removed once
every read type that reads them is profiled (the SAMs, profiles and dumps stay); `--profile-blocks 0`
profiles both collections in one protal run once all is simulated, and keeps the reads. A design
point whose largest sample would take more than `--stream-above` GB of compressed reads (2) is not
written to the disk at all: a protal run of its own reads its samples from named pipes while
`simulate_metagenomes` makes them, one sample at a time, as plain FASTQ (`--compressed-pipes`:
compressed as before), its single-end samples from a second run that writes only the first reads
(read 2 draws from a random stream of its own and is not made then); protal aligns a sample whose
reads come through a pipe whatever SAM it has. Every genome simulated is read from its FASTA once and
memory-mapped by the later samples and simulations of both collections, from a genome store on the
samples' disk (`--genome-store`). The models are trained in parallel. Without pbsim3's model file, leave `ont` out; a read type
left out gets a placeholder model that reports nothing and makes protal warn. Replace it later with
`protal --add_model MODEL --read_type se --db DB`.

#### Illumina reads

`simulate_metagenomes` makes the paired-end reads itself (ART is no longer used): a fragment of normal
length (`--fragment_mean`, `--fragment_stdev`) from a genome's contigs, either strand; read 1 from its
start, read 2 from its reverse complement's. A base's quality follows the instrument's published mean
by cycle (a curve per read, read 2 below read 1), with a run's, a cluster's and a read's offset and
noise along the read; now and then a read falls into a low-quality state for a few cycles, more
often towards its end, and a few reads end in Q2. A base is miscalled with its quality's probability
(more after GG), as the instrument's substitutions do; rare insertions and deletions (more in
homopolymers) and Ns. Binned instruments write their bins. The instruments (`--sequencer`, the
second field of a read setup) and what their reads give:

| instrument | written qualities | substitutions R1 / R2 | from |
|---|---|---|---|
| `HS20` HiSeq 2000 | 2-41 | ~0.34% / 0.46% | HiSeq 2x126 curve, Schirmer 2016 |
| `HS25` HiSeq 2500 | 8 bins | ~0.30% / 0.44% | the same, Illumina's 8-level bins |
| `HSXt` HiSeq X Ten / 4000 | 7 bins | ~0.26% / 0.47% | NovaSeq 2x151 curve |
| `NovaSeq` NovaSeq 6000 | 2, 12, 23, 37 | ~0.20% / 0.43% | NovaSeq curve, Illumina's RTA3 note |
| `MSv3` MiSeq v3 (250/300 bp) | 2-41 | ~0.34% / 0.72% (250 bp) | MiSeq 2x301 curve, Schirmer 2015 |

The curves are InSilicoSeq's (Gourlé et al. 2019); the sources and the fit are in the
[report](claude/2026-10-07-illumina-model/README.md). `simulate_metagenomes --illumina_report PAIRS
--sequencer X --read_length L` prints an instrument's qualities by cycle and error rates;
`--mean_quality Q` shifts each read's qualities (and errors) to average Q, as the scenarios do.

#### An independent test set

Cross-validation only judges the model on what the training design contains: a model trained from
5,000 read pairs up missed 8% of the present taxa of 1,000-pair samples while its own estimate said
F1 0.997 ([report](claude/2026-09-30-clade-holdouts.md)). So every build also profiles a test set of
another design: depths of 500 to 5M pairs, 10-300 species, more uneven abundances (sigma 2.0), more
mixed strains, another seed (`--test-*` options; `--test-samples 0` skips it).

#### Scenarios: kinds of studies

The design spans depths, community sizes and read setups so that one model learns them all; it does
not say how a model does on one kind of study. So every build also simulates samples of five such
kinds (`scripts/scenarios.py`; `--scenarios` picks others, `none` leaves them out): each scenario's
hold-in samples (`--scenario-samples`, 10; 6 before 2026-10-07) join the training data and inform the
models as the design's samples do, at a quarter of their weight (`--scenario-weight` 0.25), its
hold-out samples (`--scenario-test-samples`, 4, another seed) the test set, and each read type's report scores both
(section "Scenarios"; `summary.txt` has a row per scenario and set). A scenario is a community
sequenced by several technologies at the same bases, each technology reading the same communities,
each sample at a depth of its own around the scenario's (below):

| Scenario | Community | Reads |
|---|---|---|
| `gut` | 150-1,000 species, 5% lacking from the database | Illumina PE 150 bp at Q35, 20M pairs; PacBio HiFi and Nanopore at 6 Gb |
| `moderate` | 1,000-5,000 species, 30% lacking (freshwater, marine, sludge) | Illumina PE 150 Q35 at 10M pairs; Ultima at 10M reads; PacBio and Nanopore at 3 Gb |
| `soil` | 3,000-11,000 species, 60% lacking | Ultima Genomics SE 300 bp at Q25, 20M reads; Illumina PE 150 Q35 at 20M pairs; PacBio and Nanopore at 6 Gb |
| `soil_shallow` | soil communities (3,000-11,000 species, 60% lacking) | Illumina PE 150 Q35 at 5M pairs; PacBio and Nanopore at 1.5 Gb |
| `host` | 90% of the reads human, 2-50 species of power-law abundances (alpha 1: rank-abundance slope -1 on log-log axes), 5% lacking | Illumina PE 150 Q35 at 10M pairs; Ultima at 10M reads; PacBio and Nanopore at 3 Gb |

All but `host` have lognormal abundances of sigma 1.0, 1.5, 2.0 and 2.5, the samples' in turn (sigma
1.5 and 2.0 before 2026-10-07, when gut had 350-450 species, soil 9,000-11,000 and there was no
`moderate`: no training sample lay between the design's samples, at most ~4,600 taxa with reads, and
the soils', 9,800 and more; [report](claude/2026-10-07-r226-v15/README.md)).

`--scenarios gut,host` picks some (`all`, the default: every one, `gut:5`: 5 hold-in samples of one);
`--scenario-file` (JSON of scenarios by name) adds scenarios or changes a preset's fields, e.g.
`{"soil": {"species": "4000-5000"}}`. A scenario's fields: `species` (N or MIN-MAX), `novel_share`,
`abundance`, `strains`, `congeners` (as the build's options), `host_share`, and `reads`, one entry per
read type with its `depth` (read pairs, reads, or bases for long reads): `pe` with `length`,
`profile`, `fragment_mean`, `fragment_sd` and `quality` (the mean base quality), `se` with an Ultima
`setup` (`ultima:LENGTH_MEAN:LENGTH_SD:Q_MEAN:Q_SD`, default `ultima:300:40:25:2`), `pb` and `ont` with
an optional `setup` (the build's `--pb-setup` and `--ont-setup` by default); and optionally
`depth_range` (`LOW-HIGH`, default `0.125-2`) or `depth_spread` (S: 1/S to S), not both (below).

How the parts are made:
- **The depths.** Each sample's depth is the scenario's times a factor from `depth_range`'s LOW to
  HIGH (by default from 1/8 to twice it: soil 2.5M to 40M pairs, shallow soil 0.6M to 10M; from half to
  twice it before 2026-10-07), log-uniform and stratified so that a few samples
  still spread over the range, drawn from the build's seed (the hold-out samples, of another seed,
  fall between the hold-in ones). Sample s of every technology has the same factor, so the
  technologies still read the same communities at the same bases; `meta_read_pairs` in the tables is
  the scenario's depth, and the collection's log lists each sample's factor and species (`scenario soil
  (10 samples): ...; its samples at 1.62, 0.15, ... times its depths; of 4211, 9873, ... species`).
  Before 2026-10-06 every sample of a scenario had its depth, so the sample's depth feature told a
  scenario's few samples apart: a model could learn each training sample's own offset, which
  cross-validation by species does not see, and gradient boosting fell from 0.937 on the shallow-soil
  training samples to 0.881 on the hold-out ones ([report](claude/2026-10-06-r226-v13-soil/README.md)).
  From half to twice the depth a scenario's samples were still a narrow cluster, and the one at its edge
  was scored like the design's samples ([report](claude/2026-10-07-r226-v15/README.md)); below 1x the
  samples are cheaper, so ten samples at 1/8-2x cost about what six at 0.5-2x did. `depth_spread` 1
  gives every sample the scenario's depth.
- **The species per sample.** Each sample's count is drawn the same way from the scenario's MIN to MAX
  (log-uniform, stratified), from a random stream of its own, so that a community's size does not
  follow its depth (`simulate_metagenomes --species_per_sample N1,N2,...`).
- **The share the database lacks.** A scenario draws its species from a genome table of its own
  (`scenarios/<name>/genomes.tsv` beside the samples): every species of the build's table on the side
  of the split (held out or not) that is short of the share, and a random part of the other side, so
  that a species drawn uniformly is held out with the scenario's share. A sample's share varies
  around it. The table needs more species than a sample takes, so `soil` (up to 11,000 species, 60%
  held out: 6,600 held-out species and 4,400 others at least) needs a pool of about 25,000 species to
  simulate from, the default download's since 2026-10-05, or a larger `--holdout`. A scenario the
  table cannot hold is scaled down to it, and the run warns (`WARNING: scenario soil scaled from
  9000-11000 to ...`, also in `build_metadata.tsv`): its largest sample takes two thirds of the table,
  its smallest in proportion. With a download of 8,000 species (the defaults before 2026-10-05, about
  2,600 of them held out) soil and shallow soil hold about 2,400-2,900 species per sample.
- **Illumina reads at a quality.** The instruments' qualities have their own mean (HiSeq X Ten, 150 bp:
  about Q36.5 for first reads, Q35 for second reads); `--mean_quality` shifts each read's qualities to
  the scenario's, and the errors follow them ([Illumina reads](#illumina-reads)).
- **Ultima reads**: single-end, length from a gamma distribution, errors mostly homopolymer length
  errors (homopolymers from two bases), base qualities averaging the read's quality, made by
  `hifi_reads.py`'s flow model (in `simulate_metagenomes`) from templates drawn like long reads'. They are profiled as
  single-end reads and train the `se` model.
- **A host.** The host genome (`--host-genome`, by default the download's) is written once as plain
  sequence beside the samples (3.1 GB for the human one) and read by memory map. Without one the
  default scenarios run without `host` and the run warns (rerunning `download_gtdb.py` on an older
  inputs folder fetches only the host genome); `--scenarios host` without one stops the build. It gives `host_share`
  of a sample's read pairs (Illumina: the community is simulated at the rest; the host's pairs follow,
  made in the same run from fragments drawn from it, named `h_<n>`) or of its bases (Ultima, PacBio, Nanopore:
  templates drawn among the community's by their share). Host reads are in no truth file; they reach
  the profile only through spurious alignments, and the sample's depth feature counts only what lands
  on taxa.

Scenario samples are deep: one sample of every preset takes ~52 GB of reads (gut 15, soil 22,
soil_shallow 4, host 11; measured per Gb on 2026-10-05; the varied depths average ~1.08 times the
preset's), ~500 GB for the default 6 + 3 samples if they were all on the disk at once, which
`--profile-blocks` avoids, and hours more of simulation and profiling on a 64-thread node (the r226 v13
build, with 3 + 2, took 3 h 05 min in all; [report](claude/2026-10-05-collector-profiling/README.md)). Their rows also weigh in the training: a
soil sample adds ~10,000 rows, as many as 70 design samples, and at r226 v12 the scenarios' rows were
52-80% of every training table. At full weight they made the forests call more freely on the design's
samples (test F1 -0.002 to -0.004; on gut -0.008 for long reads), while without them soil's F1 fell by
0.01-0.15; at a quarter of the weight (`--scenario-weight`) soil kept its gain within 0.002 and half to
two thirds of the cost went ([report](claude/2026-10-05-r226-v12-scenarios/README.md)). The section
"Scenarios" gives the design's samples with species held out beside them, to show what the scenarios
cost the rest.

#### The reads behind the errors

The build keeps what each model error rests on, read by read, for every sample of the training data and
the test set, of every read type (`--error-reads`, default `all`; `none` turns it off, and `pe`,
`pe:design` or `pe:soil` keep a read type's samples, its design's or one scenario's). protal writes these
samples' SAMs with an unmapped record (flag 4, no sequence, `ZF` naming the taxa the read seeded on) for
every read that seeded on taxa but aligned nowhere: the non-hits, which protal otherwise only counts in
the SAM header (the collector's `--unmapped_reads`, a map's `UNMAPPED_READS` column; the profiles are the
same either way). Once the models are trained, `scripts/error_reads.py` takes each trainer's calls of
every row (`trained_model*.calls.tsv.gz`: the training rows' with species held out, the test table's by
the final model, at the model's knob) and follows the reads (the fragments, of paired-end reads) that
align to a false positive (`FP:<taxid>`) or a false negative (`FN:`), that seeded on a false negative
or an unseen species but did not align to it (`seeded:`), or that come from a genome of a false negative
or unseen species, wherever they went (`source:`). It writes per sample, to
`model_logs/error_reads/<read type>/<training|test>/<design point>/`:

| File | |
|---|---|
| `<sample>.taxa.tsv` | each error taxon: FP, FN, or `unseen` (a species of the sample that the training database has but that protal profiled no reads to: never scored, so not in the models' counts, but in protal's; not the species held out of it, which its `genome2tiid.tsv` keeps without genes: the build passes `heldout_species.txt`, `--heldout`); its score and knob; for FN and unseen species their genomes (and read pairs simulated, paired-end), their fragments with a record and where their best records went (on the taxon, at MAPQ 4 or more as the profiler counts them, elsewhere and on which taxa, or nowhere) and those that seeded on the species but did not align to it (`own_seeded_not_aligned`); for all, the fragments on the taxon and their sources, and those that seeded on it but failed to align. It counts every fragment |
| `<sample>.FP.sam.zst`, `<sample>.FN.sam.zst` | with `--share-logs`: the records of the reads taken for the false positives, and for the false negatives (`FN:`, and `seeded:` and `source:` of a false negative), of at most 20 fragments per taxon and reason (the same ones at any cap: the lowest CRC-32 of the read name); each record with QUAL left out (`*`), its source genome and species (`xg:Z:`, `xs:Z:`, "(not in the database)" for a species the training database lacks) and why it was taken (`xe:Z:`, its reasons); the `@SQ` lines cut to the genes they name |

and `summary.tsv` (one line per sample; the console sums it per read type). Unseen species are counted
in the tables but their reads not kept: a soil sample has hundreds (v15, which counted the held-out
species too, thousands), and their records made up most of v15's 9 GB of error-read SAMs. `scripts/error_reads.py` itself keeps
them with `--sams FP,FN,unseen` (`<sample>.unseen.sam.zst`), every fragment with `--max-fragments 0`, and
QUAL with `--qualities`. A read's source is its name:
`simulate_metagenomes` names a paired-end read after its contig (`<contig>-<n>`; single-end samples are the paired-end reads' first), which the genomes'
FASTAs tell (read once per build, `genome_contigs.tsv.gz` beside the samples, shared with
`trace_relatives.py`), and the collector names a drawn read (PacBio, Nanopore, Ultima) `g<i>x_<n>`, `i` the
genome's place in its community's manifest; only the host's paired-end reads have no known source. A read
that seeded on no taxon has no record: most of a genome's reads lie outside its marker genes. The
simulations are seeded, so the collector replays a sample's reads byte for byte when they are needed.
The unmapped records make the SAMs on the samples' disk larger (at r226 a deep paired-end sample had 45M
of them; an estimate of 10-20 GB more for the whole build, mostly the deep scenario samples).

`error_reads.py` extracts `--threads` samples at once, the largest SAMs first, while their estimated
memory (300 MB and twice the SAM's size on disk) fits in `--memory` GB (default 60% of the least of the
machine's memory, `SLURM_MEM_PER_NODE` and the process's cgroup limit), and keeps only what the tables
look up per fragment (its reasons, and of the taxa it aligned or failed to align to only the sample's
error taxa). A worker killed from outside (out of memory) costs no other sample: those it ran beside run
again one at a time, a sample killed again is named and left out. At r226 v15 the paired-end extraction
was killed (84 workers started on the 18 soil samples at once, each tracking millions of fragments at
1.2-2.3 kB with every taxon kept and the held-out species' reads counted), and its read type lost every
sample's tables ([report](claude/2026-10-07-r226-v15/README.md)).

#### Logs to share

`--share-logs` keeps the error reads' SAMs ([above](#the-reads-behind-the-errors)) and ends the build by
packing `OUTDIR/<name>_share.tar.gz` (`<name>`: the last part of `--outdir`), to copy off the cluster and
unpack elsewhere (`tar xzf`; everything is under `<name>/`):

| Path | |
|---|---|
| `console.log` and the other `*.log` | the console's lines (every run of the build, after a line with its time and command) and each step's log |
| `model_logs/` | as [above](#outputs), with the error reads' tables and SAMs |
| `training/`, `test/` | the training and test tables of each read type, their numbers to 9 significant digits instead of Python's 17, at about 72% of the gzipped size (v15's pe table). That keeps every float32 value, so a forest (float32 comparisons) trains the same; boosting compares doubles, so its thresholds can move in the tenth digit, but values closer than that are rare and share a bin |

It holds each file once: the trainer's outputs beside `model_logs/` (predictions, thresholds) are copies
of what is in it, and `internal_taxonomy.dmp` and the models are left out. The console's last line gives
its size.

#### Local scratch

The simulations write and delete many files, which a network file system is slow at.
`--scratch DIR` makes the samples, the training database and the converter's temporary files on a
node's own disk, and copies the tables to `--outdir`. The samples are profiled as they are simulated
and their reads removed once profiled (`--profile-blocks`), and a simulation that would leave less
than `--keep-free` GB (30) waits for that, so the disk holds the training database (~24 GB at r226),
what is simulated but not yet profiled, the SAMs and profiles, and the genome store (`--genome-store`,
~50 GB at r226, kept for the next build); give it 200-250 GB. With
`--profile-blocks 0` every sample is on the disk at once: the r226 design took up to 120 GB without
the scenarios (give it 175 GB), the default scenarios ~500 GB more (estimated for 6 + 3 samples each). The console says how much the run
takes there. A rerun reuses the samples only from the same DIR, so on a disk that is cleared
after the job, a rerun simulates again (the builds are kept in `--outdir`).

#### Time and memory

The simulations need no database, so they start as soon as the held-out species are chosen. They run
in the background at a lower priority than the index builds, long reads and paired-end points in one
queue, the longest first. The finished database is built in the background too, beside the training
database; it is only needed at the end, to take the models. That needs the memory of two builds at
once. `--one-build-at-a-time` builds the finished database after the training instead.

The marker genes of every genome (`full_reference.fna.zst`, 86 GB uncompressed at r226) are written
compressed and removed once each build has used them. The finished database is compressed at
`--final-db-level` (9), the training database at `--training-db-level` (3).

#### Stopping and rerunning

Before it starts, the script checks the tools. protal and the simulator must be built from the
source it is at: of its version and commit, with nothing changed since in `src/`, `lib/` or the
build files. Otherwise it stops at once and says what to rebuild (`--no-binary-check` runs them
anyway). When the script stops (a failure, Ctrl-C, `SIGTERM`), it stops every command it started. A
background build that fails stops the run within seconds.

A rerun into the same `--outdir` resumes:
- The conversion and both index builds are skipped when their inputs are unchanged (the release,
  the converter, protal, the held-out species; recorded in `.stages/`).
- The in-silico strains are kept when the genome table and options are the same.
- The collector reuses samples and profiles of the same design, database and protal. It simulates
  again, and says so, rather than mixing designs.
- Training, parity checks and packaging run on every rerun.

Another `--seed` or `--holdout` therefore rebuilds the training database and collects again, but
keeps the finished database.

#### Options

| Option | Default | |
|---|---|---|
| `--inputs` or `--gtdb` | required | a `download_gtdb.py` folder, or an extracted release (then `--extra-genomes DIR`, `--simulate-species FILE` or `--genome-table FILE` say what to simulate from) |
| `--outdir` | required | the output root |
| `--protal`, `--simulator` | on `$PATH` | the binaries |
| `-t, --threads` | 8 | |
| `--scratch` | | a node-local folder for the samples ([above](#local-scratch)) |
| `--read-compression` | zstd | how the simulated reads are written: `zstd` (`.fq.zst`; `simulate_metagenomes --reads_compression zstd`, the long and Ultima reads (`--long_samples`), the host's reads) or `gzip` (`.fq.gz`); smaller (about 15% against the simulator's gzip, which ISA-L writes about twice as fast) and several times faster to write than Python's gzip ([report](claude/2026-10-05-zstd-reads/README.md)) |
| `--profile-blocks`, `--keep-free` | 20, 30 | profile the samples as they are simulated, in protal runs of at least this many GB of reads, removing the reads profiled (0: one protal run once all is simulated, reads kept); the GB a simulation leaves free on the samples' disk, or it waits |
| `--stream-above` | 2 | with `--profile-blocks`: a design point whose largest sample would take more than this many GB of compressed reads is streamed into protal through named pipes, never written (0: none) |
| `--compressed-pipes` | off | with `--stream-above`: what goes through the pipes compressed as the files' names say, not as plain FASTQ |
| `--genome-store` | auto | `simulate_metagenomes --genome_store` for both collections: each genome simulated read and parsed from its FASTA once, written at 2 bits a base and memory-mapped by every later sample and simulation (at r226 ~1.5M genome reads of ~54k genomes otherwise); `auto`: `SCRATCH/genome_store` (`OUTDIR/genome_store` without `--scratch`), kept for the next build; a folder; or `none`. ~0.25 bytes a base of the genomes simulated, ~50 GB at r226 besides the samples ([report](claude/2026-10-07-simulate-metagenomes-audit/README.md)) |
| `--seed` | 1 | |
| `--holdout`, `--holdout-clades`, `--holdout-max-share` | 0.3, `phylum:2,class:4,order:6,family:8,genus:12`, 0.02 | species and clades left out of the training database; `--holdout-clades none` for species only |
| `--holdout-species` | | a file of the species to leave out instead |
| `--novel-clades-per-sample` | 1 | species of held-out clades per sample and rank |
| `--insilico-strains`, `--insilico-ani` | 1, | share of one-genome species given an in-silico strain; or their ANI drawn from `MIN-MAX` |
| `--samples` | 12 | samples per design point |
| `--read-pairs` | `1000,2000,5000,20000,50000,100000,200000,500000,2000000:4,10000000:2,30000000:1` | depths, one design point each; `DEPTH:SAMPLES` for another number of samples |
| `--read-setups` | `100:HS20:300:40,150:HSXt:350:50,250:MSv3:550:50` | read length : instrument ([Illumina reads](#illumina-reads)) : fragment mean : SD |
| `--species-per-sample`, `--archaea` | `20-200`, 2 | species per sample; archaeal species per sample |
| `--strains-per-species` | `0.5,0.2` | probabilities of a second, third, ... strain (the test set's; `0.3,0.1` before 2026-10-07) |
| `--congeners` | `0.25:2-5` | about this share of a sample's species in groups of 2-5 congeners; `0` for none |
| `--abundance` | `lognormal:1.3,2.0` | `lognormal:SIGMA[,SIGMA...]` (the samples take the sigmas in turn), `powerlaw:ALPHA` or `negbin:R:P` |
| `--read-types` | `pe,se,pb,ont` | the read types to train |
| `--long-read-bases`, `--long-read-samples` | 300 kb to 6 Gb, 36 | long-read depths and samples per point |
| `--pb-setup`, `--ont-setup`, `--pbsim`, `--pbsim-models` | | how long reads are made ([above](#one-model-per-read-type)) |
| `--test-samples` | 4 | test-set samples per design point; 0 for none |
| `--test-read-pairs`, `--test-species-per-sample`, `--test-abundance`, `--test-strains-per-species`, `--test-long-read-bases`, `--test-long-read-samples` | `500,...,5000000:2`, `10-300`, `lognormal:2.0`, `0.5,0.2`, 150 kb to 3 Gb, 8 | the test set's design |
| `--scenarios`, `--scenario-file` | `all` | kinds of studies to train on and score ([above](#scenarios-kinds-of-studies)): `gut`, `moderate`, `soil`, `soil_shallow`, `host`, `all`, `NAME:N`, `none`; JSON of more or changed ones |
| `--scenario-samples`, `--scenario-test-samples` | 10, 4 | each scenario's hold-in samples (training data; 0: scored, not trained on) and hold-out samples (test set); 6 and 3 before 2026-10-07, 3 and 2 before 2026-10-06 |
| `--scenario-weight` | 0.25 | the weight of the scenarios' hold-in rows in the models, the design's 1 |
| `--host-genome` | the download's | the host genome of scenarios with host reads |
| `--error-reads` | `all` | the samples whose SAMs keep the non-hits, and whose reads behind each model's false positives and false negatives `model_logs/error_reads/` follows ([above](#the-reads-behind-the-errors)): `all`, `none`, or `READ_TYPE`, `READ_TYPE:design`, `READ_TYPE:SCENARIO` |
| `--share-logs` | off | keep the error reads' SAM records (`<sample>.FP.sam.zst`, `<sample>.FN.sam.zst`), and end by packing `<name>_share.tar.gz` ([above](#logs-to-share)) |
| `--features` | `normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood` | the models' features (default since 2026-10-06, `auto` before: the set every model of the r226 v12 and v13 builds chose, the reference's k-mer uniqueness, the sample's complexity, which needs a protal of 2026-10-06 or later for the training data, and since 2026-10-07 the three groups against the false positives of complex communities, untested at r226, which need a protal of 2026-10-07 or later: `normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity` trains without them); `auto`: each trainer chooses its set ([below](#training)), which doubles a boosted model's training with `--evaluation basic`; or feature groups ([features.md](features.md)); `+priors` adds GTDB's species constants (opt-in since 0.7.6) |
| `--model` | `gbm` | the models: `gbm`, gradient-boosted trees (the default since 2026-10-06), or `forest`, a random forest ([below](#training)) |
| `--rounds`, `--ntree`, `--maxnodes` | 250, 64, `63` (`512,pb:128,ont:128` for forests) | boosting's rounds, a forest's trees, and leaves per tree by read type (`N` or `TYPE:N` items) |
| `--call-mode` | `curve` | `fdr` also stores calibrated calls at a target share of false calls ([below](#calls-at-a-target-share-of-false-calls)) |
| `--evaluation`, `--previous-procedure` | `basic`, off | how much the trainer evaluates: `basic` (default since 2026-10-06, `full` before), the models with rows, samples and species held out, which the summary and the knob need; `full` adds the clades held out and the studies ([below](#training)), several times a boosted model's training |
| `--n-genes`, `--genes`, `--gene-ranking`, `--genes-per-domain`, `--rank-genes` | | a reduced database ([below](#reduced-marker-sets)) |
| `--no-gene-neighbours` | | skip the gene neighbours; protal then pairs no mates across neighbouring genes |
| `--one-build-at-a-time` | | build the finished database after the training |
| `--training-db-level`, `--final-db-level` | 3, 9 | zstd levels of the two databases |
| `--no-placeholder-models` | | no placeholder models for read types not trained |
| `--no-binary-check` | | run binaries of another source anyway |
| `--progress-every` | 0 | seconds between status lines; 0 for none |

### Reduced marker sets

A database of a subset of the marker genes takes less memory: 12 of the 120 bacterial markers give
about a tenth of the index (the reduced database on the downloads page is one). In one command:

```bash
python3 scripts/build_gtdb_database.py --inputs /shared/protal_inputs/gtdb_r226 --outdir /data/protal_r226_n12 \
    --n-genes 12 --protal build/protal --simulator build/simulate_metagenomes -t 64 --scratch /local/scratch
```

The genes are chosen by `scripts/rank_genes.py` from a full build of the training database, and
the models are trained on the reduced database, as they must be. `--gene-ranking FILE` takes the
ranking of an earlier full build instead (a full build writes one with `--rank-genes`).
`--genes LIST` names the genes. The choice is in `gene_subset.txt`, `gene_ranking.tsv` and
`build_metadata.tsv`.

#### Which genes

`rank_genes.py` scores each gene by prevalence (the share of species whose copy has k-mers in the
index) times unique share (the share of a copy's k-mers unique to its species). It reads a build's
separate files (`--no_bundle`, or `protal --unpack_db`):

```bash
python3 scripts/rank_genes.py --db /data/protal_r226_files -o gene_ranking.tsv --top 12 --subset genes.txt
```

GTDB's bacterial and archaeal marker sets (120 and 53 genes) share five genes, and archaea are a
few percent of the species. A subset by the overall score alone would hold no archaeal gene, so
`--top N` reserves `--per-domain M` genes for each domain's best (a third of N by default).
`gene_ranking.tsv` also gives each gene's mean length, its congeners' divergence and its suspect
copies, so that a subset can be chosen on other grounds. A gene list is GTDB marker ids (`PF00380.20`,
`TIGR00001`) or protal gene ids, comma-separated or one per line.

#### A folder of the subset

`gtdb_to_protal_db.py --genes LIST` writes a converted folder with those genes only, from the release
or from a converted folder (`--from_db`), with the gene neighbours counted anew over the genes kept.
Build it as any other:

```bash
python3 scripts/mini_db/gtdb_to_protal_db.py --from_db /data/protal_r226 --genes genes.txt --outdir /data/protal_r226_12
protal --build --no_profile --db /data/protal_r226_12 --reference /data/protal_r226_12/reference.fna \
    --full_reference /data/protal_r226_12/full_reference.fna -t 16
```

`protal --build --build_gene_subset genes.txt` (protal gene ids) is the quick version on a full folder:
the other genes stay in `reference.fna` but get no k-mers, so no read lands on them. It refuses a
folder with `gene_neighbours.tsv`, which was counted over every gene.

### Step by step

#### GTDB files

Download the release from https://data.gtdb.ecogenomic.org/releases/ (for r226: `release226/226.0/`),
keep GTDB's layout and extract the tarballs where they are (e.g.
`genomic_files_reps/bac120_marker_genes_reps_r226/fna/*.fna`):

| File | Needed for |
|---|---|
| `bac120_taxonomy_r226.tsv`, `ar53_taxonomy_r226.tsv` (or `.gz`) | lineages |
| `bac120_metadata_r226.tsv.gz`, `ar53_metadata_r226.tsv.gz` | species representatives and quality (optional) |
| `genomic_files_reps/{bac120,ar53}_marker_genes_reps_r226.tar.gz` | the reference |
| `genomic_files_all/{bac120,ar53}_marker_genes_all_r226.tar.gz` | unique k-mers and gene conservation (optional, but without it uniqueness is judged on the representatives alone) |
| `auxillary_files/sp_clusters_r226.tsv` | species priors (optional) |
| `genomic_files_reps/gtdb_genomes_reps_r226.tar.gz` | whole genomes, only to simulate training data |

#### 1. Convert the release

```bash
python3 scripts/mini_db/gtdb_to_protal_db.py --gtdb /data/gtdb_r226 --outdir /data/protal_r226_db -t 16
```

It writes:
- `reference.fna` and `reference.map`: the representatives' marker genes;
- `internal_taxonomy.dmp`;
- `full_reference.fna.zst`: every genome's marker genes, 86 GB uncompressed at r226; only if
  `genomic_files_all` is there;
- `species_priors.tsv` ([below](#species-priors));
- `model_pe.xml`;
- for your own use, `gene2geneid.tsv` (marker id to protal gene id) and `genome2tiid.tsv`
  (accession, species taxid, representative, lineage).

The default `model_pe.xml` was trained on older databases and does not call archaea reliably.
Train a model for a database you will use ([below](#the-presence-model)), or take the
one-command route.

| Option | Default | |
|---|---|---|
| `--release` | from the file names | e.g. `226` |
| `--model` | `scripts/random_forest.xml` | the model copied to `model_pe.xml` |
| `-t, --threads` | 1 | marker files read in parallel; the output is the same for any number |
| `--exclude_species` | | species whose genes are left out, keeping their taxids (a training database) |
| `--from_db` | | copy a converted folder instead of reading the release (with `--exclude_species` and `--genes`) |
| `--genes` | | only these marker genes ([reduced marker sets](#reduced-marker-sets)) |
| `--order` | `gene` | `reference.fna` by gene, then taxon (compresses ~2x better), or `genome` |

#### Species priors

`species_priors.tsv` holds what GTDB knows of each species before any read, for the model's
optional prior features ([features.md](features.md#the-species-priors-priors-075-opt-in)):
- the marker genes found in the representative, and how many twice (always 0 for a GTDB release,
  whose marker files hold one copy per genome);
- the representative's CheckM completeness and contamination;
- from the species clusters file, the cluster's ANI radius, mean and minimum intra-species ANI and
  number of genomes.

-1 means unknown. `--build` packs the table into `database.protal`.

#### Gene neighbours

Optional, between conversion and build: which marker genes lie next to each other, per clade. With
it, protal pairs mates across neighbouring genes and follows long reads from gene to gene. GTDB's
marker files do not say where a gene lies, so this needs whole genomes, those downloaded to simulate
from:

```bash
python3 scripts/mini_db/gene_neighbours.py --db /data/protal_r226_db --genome_table genomes.tsv -t 16
```

The genome table is the simulator's (accession, GTDB taxonomy, FASTA path).

- **Placing.** A database gene is found by its exact sequence in its representative, and in other
  strains by its k-mer trace: its 24-mers on a band of near diagonals, which places genes up to about
  9% different. A genome in which fewer than 80% of its species' genes are placed is skipped
  (`--min_placed`).
- **Counting.** For each gene end, the next marker within 3,000 bases (`--max_gap`) and how far it
  is, or none. An end near a contig's end says nothing; a circular sequence has no ends. Each species
  counts once, per clade (family to domain): in how many species an end faces a partner, of how many
  in which the end is informative.
- **A species' own lines.** Where a species' gene order differs from its family's, it gets lines of
  its own, so that its own genes are not taken for genes from elsewhere (`--no_species_lines` leaves
  them out).
- **Use.** A small clade's frequency leans on the clade above it. Two ends facing each other in 20%
  or more of the nearest clade are expected neighbours, in 5% or less unlikely. Mates and long reads
  are only followed over expected neighbours ([running.md](running.md#options-the-website-does-not-list)).

It writes `gene_neighbours.tsv` (the frequencies, which a run loads) and `gene_positions.tsv` (where
each gene lies in each genome, which runs do not read). `--build` checks and packs both.
`gtdb_to_protal_db.py --from_db` derives a training database's frequencies from the positions,
without the species it leaves out ([report](claude/2026-10-02-gene-neighbour-frequencies/README.md)).

#### 2. Build the index

```bash
protal --build --no_profile -t 16 --db /data/protal_r226_db \
    --reference /data/protal_r226_db/reference.fna \
    --full_reference /data/protal_r226_db/full_reference.fna
```

Pass `--no_profile`; without it, build mode goes on to profile an empty sample list. protal reads
the `.zst` file when the plain one is missing. The build:

1. Indexes `reference.fna`, checks every k-mer's uniqueness against the full reference, and writes
   `unique_kmers.tsv`.
2. Writes `gene_conservation.tsv`: how fast each gene diverges within species against the species'
   other genes (strains' copies compared with the representative's, k = 12; the median gene 1).
   Queries use it for the conservation and divergence features. It needs other genomes' copies.
3. Writes `gene_congeners.tsv` beside the database: how each gene differs between congeneric species
   against within species, to check that conserved genes are conserved between species too
   ([report](claude/2026-10-01-conservation-pattern/README.md)).
4. Looks for **suspect gene copies**: a species' copy within 0.02 (`--suspect_copy_distance`) of
   another genus's copy while its own congeners' copies are farther, a contaminating contig or a
   transferred gene. At r226 such copies drew a fifth of the false species calls. They go into the
   database as `suspect_copies.tsv`, and runs leave reads on them out of the evidence
   (`--keep_suspect_copies` keeps them). Every near pair across genera is listed in
   `gene_incongruence.tsv` ([report](claude/2026-10-03-false-positive-anatomy/README.md)).
5. Compares every two species of a genus by their marker genes (the median Mash distance of the
   genes both have, k = 12, as a run's `relative_distance`) and writes each species' nearest
   congeners, at most 16 within 0.15, to `species_neighbours.tsv`: how crowded the database is
   around a reference, and how far apart two congeners a read fits are (since 2026-10-06; genera in
   batches of about 10,000 species, so the sketches of the whole database are never held at once).
6. Packs everything into `database.protal`, with `gene_table.bin` for a fast load, reads it back,
   compares, and removes the separate files. `full_reference.fna`, `gene2geneid.tsv` and
   `genome2tiid.tsv` stay beside it.

Every phase uses `-t` threads, and the index is the same for any `-t`. The log times each phase.
[Build options for the format](#build-options-for-the-format) lists the options for separate
or uncompressed files and the compression level.

#### 3. Check it

Profile samples simulated from genomes whose species you know. With
[`simulate_metagenomes`](development.md#simulating-metagenomes) the map carries the truth:

```bash
simulate_metagenomes --genome_table genomes.tsv -o sim -n 4 --total_read_pairs 200000 \
    --species_per_sample 20 --seed 1 -t 8 --protal_metafile sim/protal
protal --db /data/protal_r226_db --map sim/protal.meta -t 16
```

protal prints each sample's true and false positives and false negatives, and writes
`<profile>.truth_annotated`.

## The presence model

### How protal calls species

protal scores every species with reads with a model of trees (PMML: gradient-boosted trees since
2026-10-06, a random forest before; protal reads either) and reports those whose probability of
`TRUE` reaches a threshold:

- **`--knob`** (default 0.5), for every sample;
- **a knob curve** over the sample's depth, if the model carries one and `--knob` is not given
  ([below](#knobs-by-sample-depth)). Models trained with the sample's depth as a feature (the
  default) have none, but may carry one knob for every sample, chosen in training;
- **a target share of false calls** with `--fdr F`, for a model with calibrated calls
  ([below](#calls-at-a-target-share-of-false-calls));
- **the singleton rule** with `--singleton_congener N`, which vetoes a single-fragment species
  beside an abundant congener ([below](#the-singleton-rule)). Off by default.

A model belongs to the kind of database and reads it was trained on: its features come from the
reads and the database, and a pair counts twice. Each read type has its own model in the database
(`model_pe.xml`, `model_se.xml`, `model_PB.xml`, `model_ONT.xml`; or `--model`, `--model_se`,
`--model_pb`, `--model_ont`).

### What protal expects of a model

At start, before any alignment, protal checks that the model takes only inputs protal computes (the
features of `TaxonFeatures` in `src/Profiling/Profiler.h`, the columns of the training dump), predicts
a field with the value `TRUE`, and scores a species. A model that fails stops the run with exit code
2, naming the problem.

### Training data

With a truth file per sample (`--profile_truth`, or a map's `PROFILE_TRUTH` column), protal writes
`<profile>.truth_annotated`: one row per species with reads, whether it was present (`truth`), the
current model's call and probability, and every feature exactly as the model is given it. A truth
file names one present species per line, by GTDB lineage (`d__...;s__Genus species`) or by internal
taxid in the first field; `simulate_metagenomes --protal_metafile` writes such files and a map that
uses them ([development.md](development.md#simulating-metagenomes)).

`scripts/collect_training_data.py` makes such data at scale: it simulates samples over a grid of read
setups and depths, profiles them against a database, and joins the dumps into one table per read type
(`training_data.tsv` for pe, `training_data_se.tsv`, `_pb`, `_ont`). `build_gtdb_database.py` runs it
with a realistic design ([above](#training-data-like-real-samples)); by hand:

```bash
python3 scripts/collect_training_data.py --db DB --genome_table genomes.tsv -o training \
    --archaea 2 --species_per_sample 10-40 -t 8
```

- **Read types** (`--read_types`): `se` profiles the paired-end samples' first reads alone; `pb` and
  `ont` replay each paired-end point's communities as long reads, drawn read by read (genome by
  abundance × length, length from a gamma distribution, start uniform, cut at the contig's end), by
  `simulate_metagenomes --long_samples`, one run per design point, each read written compressed as it
  is made. PacBio HiFi reads follow `scripts/hifi_reads.py`'s model (quality by length, Q50 up to 5 kb
  to Q20 at 50 kb, errors mostly as homopolymer indels, calibrated base qualities); Nanopore reads
  pbsim3's quality-score model in its template mode (its `.model` file; pbsim3 is not run), drawing
  only the accuracy levels the model has an HMM for: pbsim3 gave ~20% of its reads level 100, Q93
  throughout and no error. pbsim3
  makes no HiFi reads: its one-pass reads have quality 0, which broke the excess features until
  2026-10-02 ([report](claude/2026-10-02-pacbio-hifi-reads/README.md)). The C++ models were checked
  against `hifi_reads.py` and pbsim3 on the same templates
  ([report](claude/2026-10-07-long-read-simulator/README.md)); `simulate_metagenomes --long_templates`
  makes one read of each sequence of a FASTA, as they do.
- **Scheduling**: all simulations share `--jobs` cores in one queue: the paired-end points with
  threads by their work (their samples' genomes and read pairs), the long-read points on threads by
  theirs, from the communities' design (seconds) rather than after the paired-end reads; then all
  samples are profiled in one protal run. A rerun skips design points simulated and profiled from the
  same inputs. The host's paired-end fragments (Python) run in worker processes
  ([report](claude/2026-10-05-collector-profiling/README.md)).
- **Columns** starting with `meta_` say where a row comes from: design point, sample, read type,
  depth, the taxon's domain, the species the database lacks in the sample and whether the taxon
  shares a genus with one (`meta_novel_*`), whether a present species was simulated from its
  representative (`meta_rep_genome`) or an in-silico strain (`meta_insilico_strain`), and how close
  its relatives in the sample are (`meta_relative_rank`, `meta_neighbour_rank`), and which scenario a
  sample is of (`meta_scenario`, empty for the design's). The training report breaks its errors down
  by them.
- **Scenarios** (`--scenarios NAME[:SAMPLES],...`, `--scenario_samples`, `--scenario_file`,
  `--host_genome`): samples of kinds of studies beside the design's, in the same tables
  ([above](#scenarios-kinds-of-studies)); `--samples 0` collects them alone. Their folders are
  `points/sc_*`, which `check_model_parity.py` and `trace_relatives.py` leave out.

What the genome table must hold for a good model: species the database lacks (their reads land on
relatives, the false positives to learn; build the training database with
`gtdb_to_protal_db.py --exclude_species`, which keeps their taxids), archaea (`--archaea`), and other
strains than the references. Congeners must share samples (`--congeners`) for the relatives
features.

| Option | Default | |
|---|---|---|
| `--db`, `--genome_table`, `-o` | required | the database, the simulator's genome table, the output folder |
| `--protal`, `--simulator` | on `$PATH` | the binaries |
| `--samples`, `--read_pairs`, `--read_setups` | 4, `1000,5000,20000,100000,500000`, three setups | samples per point, depths, read length : instrument : fragment mean : SD |
| `--species_per_sample`, `--archaea`, `--strains_per_species`, `--abundance`, `--congeners` | `5-30`, 0, one, sigma 1.3, 0 | the communities |
| `--read_types`, `--long_read_bases`, `--pb_setup`, `--ont_setup`, `--pbsim`, `--pbsim_models` | `pe` | other read types and how they are made |
| `--novel_species`, `--novel_clades`, `--taxonomy` | | the species the database lacks (`heldout_species.txt`), held-out clades per sample, and the taxonomy for the `meta_*` ranks |
| `--scenarios`, `--scenario_samples`, `--scenario_file`, `--host_genome` | none, 6 | scenarios to collect besides the design, their samples (each at a depth of its own), more or changed scenarios, the host genome |
| `-t`, `--jobs`, `--seed` | 4, `-t`, 1 | threads, cores for the simulations, seed |
| `--simulate_only`, `--prepare_profiling`, `--also_profile MAP` | | simulate and stop; write the map to profile and stop; profile another collection in the same run |
| `--read_compression` | zstd | `zstd` (`.fq.zst`) or `gzip` (`.fq.gz`): how the simulated reads are written |
| `--follow`, `--profile_block`, `--protal_lock`, `--min_free` | , 20, , 0 | profile as a `--simulate_only` run of the same `-o` simulates, in protal runs of at least this many GB of reads, removing each point's reads once profiled, then write the tables (points it cannot profile once the simulations have ended it simulates itself); a lock file for the protal runs of two followers to take turns; GB a simulation leaves free on the disk, or it waits (for a follower to remove reads) |
| `--stream_above` | 0 | `--follow` (and its `--simulate_only` run, the same value): points whose largest sample would take more than this many GB are streamed into a protal run of their own through named pipes, not written |

### Training

```bash
python3 scripts/machine_learning_cmdline.py --truth-file training/training_data.tsv --output-prefix training/model
```

The trainer needs Python 3 with numpy, pandas, joblib and scikit-learn, no Java: it writes the PMML
itself, so that protal's probabilities equal scikit-learn's bit for bit, and checks this on every
training row.

**The model.** Gradient-boosted trees (`--model gbm`, the default since 2026-10-06; scikit-learn's
`HistGradientBoostingClassifier`: 250 rounds at a learning rate of 0.1 of trees of up to 63 leaves,
20 rows a leaf at least, an L2 penalty of 1, balanced classes, no early stopping) or a random forest
(`--model forest`, the model before). Refitted on the r226 v13 tables, boosting scored higher than the
forest on every test set of the long reads (soil +0.004 to +0.008, the design's test set +0.003, gut
+0.008 to +0.013) and of Ultima's (+0.005 in soil), and, with the scenarios' depths varied (rounded in
the experiment), of paired-end reads (soil +0.004, shallow soil +0.006, the design's test set +0.002;
[report](claude/2026-10-06-r226-v13-soil/README.md)). It is written as a chain of one regression tree
per round into a logit sum (`model_pmml.write_boosted`), which protal's cPMML scores bit for bit as
scikit-learn does; a boosted model file is ~3 MB, a forest's ~13 MB. 250 rounds at 0.1 (500 at 0.05
until 2026-10-06) score within 0.001 of F1 of 500 at 0.05 on the v13 paired-end table, at half the trees
([report](claude/2026-10-06-training-time/README.md)). It trains slower than a forest: a fit of 500 rounds
took 88 s on one thread on four fifths of the r226 paired-end table and 52 s on two or three, and gained
little from more threads, so the folds of a cross-validation are fitted side by side in worker processes of
a few threads each (`--fold-jobs`: by default as many as there are folds, up to `--threads` / 3, each on
`--threads` / jobs threads; the scores are the same, as neither model depends on its threads). Its
OpenMP threads are held at `--threads`. A boosted model has no out-of-bag estimate; its importances
(`.varimp.tsv`) are the splits' gains.

| Option | Default | |
|---|---|---|
| `--truth-file`, `--output-prefix` | required | the training table; the prefix of the outputs |
| `--features` | `auto` | `auto` chooses (below); or feature groups joined by `+` ([features.md](features.md)); `+priors` is opt-in; `all` takes every column. A table of an older protal lacks newer groups: `auto` leaves out the sets it lacks, a set named must leave them out |
| `--model` | `gbm` | `gbm`: gradient-boosted trees; `forest`: a random forest |
| `--rounds`, `--learning-rate`, `--l2` | 250, 0.1, 1 | boosting's rounds, learning rate and L2 penalty (500 and 0.05 before 2026-10-06) |
| `--fold-jobs` | 0 | folds fitted side by side (0: as many as there are folds, up to `--threads` / 3; 1: one after another) |
| `--maxnodes`, `--min-samples-leaf` | 63, 20 for boosting; 256, 1 for a forest | leaves per tree at most (0: no limit), rows a leaf at least (the GTDB build gives forests 512 leaves for short reads, 128 for long reads) |
| `--ntree`, `--max-features` | 64, `sqrt` | a forest's trees and features tried per split |
| `--knob` | 0.5 | the threshold protal will use; errors are counted at it |
| `--depth-knobs`, `--fdr-calls`, `--singleton-congener` | off, off, 0 | a knob curve (with the sample's depth a feature: one knob for every sample, if it pays); calibrated calls; the singleton rule the counts follow ([below](#knobs-by-sample-depth)) |
| `--scenario-weight` | 0.25 | the sample weight of the scenarios' rows (`meta_scenario`) in every model fitted, the design's 1 ([above](#scenarios-kinds-of-studies)) |
| `--taxonomy` | | the database's taxonomy: domains, and whole clades held out |
| `--test-file` | | an independent test table, scored by the final model; its scenarios' rows (`meta_scenario`) are their hold-out samples, scored apart |
| `--evaluation`, `--folds`, `--previous-procedure` | `full`, 5, off | `basic`: held-out scores only; `none`: fit and export only |
| `--seed`, `--threads` | 1, 4 | |

**How a model is judged.** Rows of one sample share its reads, and rows of one species its reference,
so random rows would score a model on what it was trained on. The trainer scores each row with models
that saw neither its sample nor its species ("by species", the estimate that matters: most species
protal meets were never in training), and with `--taxonomy` and `--evaluation full` with its whole genus,
family, order, class or phylum held out. `--evaluation full` adds studies: the other feature sets, the model's size (a
forest's leaves and trees, boosting's leaves and rounds), and a learning curve.

**Scenarios.** For a table with scenarios (`meta_scenario`), the section "Scenarios: hold-in and
hold-out samples" gives, per scenario, as protal calls by default: TP, FP, FN, sensitivity,
precision, F1, the FP rate (of the absent taxa, those called), the FN rate (of the present taxa, those
missed), FP per sample, the false positives whose closest species in the sample is one the database
lacks, and the highest F1 at any threshold. On the scenario's training samples (hold-in) three times:
scored by the final model, which was fitted on them (in sample: how well it fits), with their samples
held out, and with their species held out; and on its samples in `--test-file` (hold-out), which no
model saw. A row for the design's training samples with species held out comes first, for
comparison. The summary has a line per scenario, and warns when a scenario's hold-out F1 is more than
0.02 below its hold-in F1 with species held out. The independent test set's section leaves the
scenarios' rows out.

**`--features auto`** (the default) scores each candidate set with species held out and keeps the one
of highest F1 at the knob, but the default set (`normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood`) unless
another beats it by 0.002 (`AUTO_MIN_GAIN`, the gain below which the depth knobs changed between fits
at r226). The candidates are the named sets without the priors, and
`normalized+adjacency+relatives+depth+divergence+unfiltered` (`auto+priors` adds the sets with the
priors, which win on simulated data for a reason the simulation makes: see
[features.md](features.md)). The section "Feature set chosen" lists each candidate's F1, AP, log loss
and errors, and, scored by each candidate's model fitted on all rows, its F1 on every test set apart:
the independent test set and each scenario's hold-out samples. These are for comparison only, since a
set chosen on them would leave them no longer held out. A line then says why the winner won: its F1
against the other sets' mean and the best of them, the rule that kept or chose it, and its F1 on each
test set against the other sets' mean (also the summary's last line, `metrics.json`'s
`features_auto.why`, and the build's console and `summary.txt`). It costs 8 candidates x 5 fits on four fifths of the rows, and one fit
on all rows each for the held-out tables: for forests at the r226 tables' size (the "Feature sets" study, 11 sets
x 2 schemes, took 10-170 s per read type, a final fit 1-6 s) roughly 10-100 s per read type, against
1-8 min for a model's whole training and 2.5-3 h for a build; for boosting about 15-30 min for the
paired-end table (48 fits of ~20-35 s; an estimate). The evaluation and, with `--evaluation
full`, the "Feature sets" study reuse the candidates' scores with species held out.

| Output | |
|---|---|
| `<prefix>.xml` | the model |
| `<prefix>.report.txt` | the evaluation (also printed) |
| `<prefix>.metrics.json`, `.predictions.tsv.gz`, `.thresholds.tsv`, `.varimp.tsv`, `.joblib` | its numbers, every taxon's held-out probability, F1 by threshold, feature importances, the scikit-learn model |
| `<prefix>.test_predictions.tsv.gz`, `.scenario_predictions.tsv.gz` | each taxon's probability on the independent test set and on the scenarios' hold-out samples |
| `<prefix>.calls.tsv.gz` | every row's score, knob and call as protal calls by default: the training table's with species held out, the test table's by the final model (for `error_reads.py`) |

`scripts/check_model_parity.py --db DB --model training/model.xml --training training` re-profiles
saved training samples and checks that protal's probabilities are the model file's and its features
those of the training data (`--read_type se`, `pb`, `ont` for other models). It profiles into
`training/parity` (`parity_<read type>`, or `-o`), which every run rewrites. Differences in the last
digits (1e-12) are noted, not failed. A sample without taxa (a shallow one of a database of few genes)
passes if it had none during collection either; the check fails only if no sample has a taxon to
compare. `gradient_boosted_cmdline.py` and
`hist_gradient_boosted_cmdline.py` are older trainers of gradient-boosted trees (their PMML export
needs Java), which `machine_learning_cmdline.py --model gbm` replaces; the R scripts are the older caret
pipeline.

### Species and clades the database lacks

The report's section "False positives and false negatives by taxonomic rank" (and the same for the test
set) shows how the model handles organisms the database lacks:

- **False positives from what the training database lacks**: absent taxa whose closest species in
  the sample is one the database lacks, by the rank it was held out at, against the other absent taxa.
- **False negatives** of present taxa by the deepest rank they share with another species of the
  sample.
- **With the taxon's clade held out of training**: FN and FP rates and F1 when the model has seen
  nothing of that species, genus, ... phylum.

The section "Strains" compares species simulated from another genome than the representative (real
and in-silico strains apart) with those simulated from it, by fragments: how often each is missed,
and whether missed strains look like found ones or like a missing species' congeners. "The
conservation features by class of taxon" checks the conservation features on the training data.

### Using a new model

Try it without changing the database, on new samples or on existing alignments:

```bash
protal --db DB --model training/model.xml --map test.map -t 16
protal --profile_only a.sam.zst,b.sam.zst --db DB --model training/model.xml -o test/
```

Then make it the database's model of its read type with `protal --add_model training/model.xml
--read_type pe --db DB`. Choose `--knob` on held-out samples like the ones it will profile;
`<prefix>.thresholds.tsv` is a start.

### Knobs by sample depth

A deep sample holds many more absent taxa with a few reads (at r226 they grew as depth^0.84), so its
best threshold is higher. With `--depth-knobs` the trainer fits a knob curve over log10 of the
sample's fragments: a point per half decade with 6 or more samples, its knob the threshold with the
highest F1 on species held out within half a decade of it, kept at `--knob` unless it gains 0.002.
The curve goes into the model (`<Extension name="protal_depth_knob_curve" .../>`), and protal
interpolates it per sample, keeping the deepest knob for deeper samples; `--knob` overrides it. The log
lists each sample's knob. At r226 (0.7.0's features) curves raised the test F1 by 0.006-0.033 per read
type ([report](claude/2026-10-02-r226-build-evaluation/README.md)).

With the sample's depth among the features (the default since 0.7.5) no curve is fitted: the model
already knows the depth, and a curve on top corrects twice (it lost 0.010). Instead the trainer chooses
one knob for every sample: the median of the best thresholds on species held out (rows weighted as in
the fits) over 200 bootstrap resamples of the training samples, if it beats `--knob` in 95% of them
and does not lose F1 on the test set (its rows weighted alike); it goes into the model as a curve of
one point, which protal reads as that knob at every depth (`--knob` overrides it). The report gives
the resamples' range of best thresholds and the test set's F1 at both. Until 2026-10-07 it was the
best threshold if it gained 0.002: on a flat curve near-equal models chose differently (se 0.73 at
r226 v14, 0.5 at v15; [report](claude/2026-10-07-r226-v15/README.md)). At r226 v12 the paired-end forest's
best threshold was 0.70 both with species held out and on the test set (+0.003 of test F1); the other
read types' sat near 0.5 ([report](claude/2026-10-05-r226-v12-scenarios/README.md)). Either way a
model cannot extrapolate beyond its deepest training samples: train at the depths you profile. Models
of 0.7.2 carry knobs per decade instead, which protal still reads.

### Calls at a target share of false calls

With `--fdr-calls` the trainer calibrates the scores (an isotonic fit on species held out) and picks the
target share of false calls (0.005-0.5) with the highest F1. With `protal --fdr F`, each sample's scores
become probabilities, adjusted to the sample's own share of present candidates, and protal reports
the highest-scoring taxa while the mean of their 1 − probability stays at or below F. At r226 these
calls were 0.001-0.007 F1 below the knob curve for every read type
([report](claude/2026-10-03-r226-v5-v6-training/README.md)), so they are off by default;
`build_gtdb_database.py --call-mode fdr` trains them.

### The singleton rule

With `--singleton_congener N`, a species with a single fragment beside a congener of N or more
fragments is never reported if its read looks like the congener's: the abundance-weighted assignment
(`em_own_share`) leaves it less than half of it, or its identity is below 0.95. Off by default (0): at
r226 the model already called none of the rows the rule would veto, and on worlds with congener groups
a rule without the read test removed true minor congeners
([report](claude/2026-10-03-false-positive-anatomy/README.md)). The profiles, statistics and MSAs
follow it, and the trainer counts it with `--singleton-congener`.
