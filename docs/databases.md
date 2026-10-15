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
| `build_settings.tsv` | since 2026-10-11: the settings `--build` used that runs should match: the WFA2 traceback (`--wfa_traceback`, [below](#2-build-the-index); a database without it was built with the classic traceback), how the `fragments` feature counts (`--fragment_count`), and since 2026-10-12 each read type's sparse-first gate (`--sparse_first`, `sparse_first_pe` to `sparse_first_ont`, with `sparse_first_score` and `sparse_first_hits`), since 2026-10-13 with its pick (`--sparse_first_pick`, `sparse_first_pick_pe` to `sparse_first_pick_ont`: `hash`, `stride` or `grid`; `--build` records `stride` where not given) |
| `unique_kmers.tsv` | per species and gene, the k-mers unique to it in the database |
| `gene_conservation.tsv` | optional: per gene, how fast it diverges within species against the species' other genes ([below](#2-build-the-index)); for the conservation and divergence features, and `--gene_conservation db` |
| `suspect_copies.tsv` | optional: gene copies near-identical to another genus's (contamination, transferred genes), whose reads a run leaves out of the evidence ([below](#2-build-the-index)) |
| `species_neighbours.tsv` | optional (since 2026-10-06): each species' nearest congeners in the database by the distance of their marker genes, for the database-neighbourhood features and `unexpected_congener_fit_share` ([below](#2-build-the-index), [features.md](features.md#against-false-positives-in-complex-communities-consistency-shape-neighbourhood-2026-10-06)) |
| `congener_gaps.tsv` | optional (since 2026-10-07): per species' copy of each marker gene, the alignment distance to its nearest congener's copy and to a typical one, and which congener is the nearest, for the `gaps` features and the ancestry sites ([below](#2-build-the-index), [features.md](features.md#the-gene-copies-gaps-and-foreign-reads-and-the-untried-candidates-gaps-foreign-untried-2026-10-07)) |
| `column_weights.tsv` | optional (since 2026-10-10, `--build --column_weights`, default 10 genera per family; 0: none): per family and marker gene, how conserved each column of the family's reference copy is among the genus references and within the genera (averaged over genera), the amino acid per codon, and the family's consensus base, each weight in half nats on 4 bits; and every copy's mapping onto the columns. For the `weights` features and the polarised ancestry sites of small genera ([features.md](features.md#the-columns-conservation-weights-2026-10-10)); `--no_site_weights` leaves it unread |
| `strain_alleles.tsv` | optional (since 2026-10-08, `--build` with `--full_reference`): per species' copy of each marker gene, up to 16 alleles (4 until 2026-10-11) of the species' other genomes (`--strain_alleles`) as edits of the representative's copy, each nearer the representative than the gene's nearest congener's copy; for the alignment scores (`--no_allele_scores`) and the `alleles` features ([below](#2-build-the-index), [features.md](features.md#the-strain-alleles-alleles-2026-10-08)) |
| `ancestry_sites.bin` | optional (since 2026-10-11, `--build`; `--no_ancestry_sites`: none): every species' copy of each marker gene with its ancestry sites (where its copy differs from the consensus of its congeners' copies, and the indels they share), computed by the build from `species_neighbours.tsv`, `congener_gaps.tsv` and `column_weights.tsv` as a run would compute them on a copy's first record; a run decodes them instead (computing them was half of a sample's profiling at r226 v23: [report](claude/2026-10-11-r226-v23/README.md)). Binary: per copy the positions, bases and indels; since 2026-10-14 (version 2) the positions as Rice-coded gaps, the index in 4-byte keys and offsets and the identity in 16 bits, 0.65-0.69 of version 1's bytes (1.13 GB at r226 v27, ~0.73 GB in version 2; a run reads either version, [report](claude/2026-10-14-length-free-features/README.md)). Without it, or with `--no_site_weights`, a run computes them as before; the features are the same ([features.md](features.md#which-side-the-reads-take-ancestry-079)) |
| `unique_cores.bin` | optional (since 2026-10-14, `--build`; `--no_unique_cores`: none): for every gene copy, in the order of its cores (the syncmers a run finds in its 2-bit sequence), two bits: whether its index entry is unique and unique at distance two (~0.7 GB at r226: two bits for each of its 2.9 G entries). A run counts each record's `ZS`, `ZU` and `ZT` from the gene's cores its alignment covers with exact k-mers, which is what the whole index's seeds find, so a sparse index or `--sparse_seeding` gives the counts of the whole index. Without it (or with `--no_unique_cores` in a run) they are the seeds of the record's anchor, as before ([features.md](features.md#independent-of-the-read-length-and-the-seeding-length_free-read_length-2026-10-14)) |
| `foreign_rates.tsv` | optional (since 2026-10-07): per gene copy, the reads of a tiled scan of the full reference's marker genes (every species' alike, a few genomes each) that landed on it and how many came from other species and genera, every copy listed (`scripts/foreign_rates.py`, stored with `--add_tables`; `build_gtdb_database.py --foreign-rates` does both for each database), for the `foreign` features ([features.md](features.md#the-gene-copies-gaps-and-foreign-reads-and-the-untried-candidates-gaps-foreign-untried-2026-10-07)); until 2026-10-08 the scan read the genomes at hand, which leaked the simulation's species ([features.md](features.md#the-foreign-features-leak)) |
| `species_priors.tsv` | optional: what GTDB knows of each species before any read ([below](#species-priors)); since 2026-10-08 also its genome size and the share of it its marker genes cover, from which a run tells how much of a sample the called species explain and the profile's unknown share ([running.md](running.md#what-the-called-species-explain-the-unknown-share)) |
| `ncbi_taxids.tsv` | optional (since 2026-10-11): each species' NCBI species taxid, the one most of its GTDB genomes have, NCBI's name of it and whether that is a named species or one of NCBI's placeholders ([below](#ncbi-taxids)); read only by runs with `--ncbi_taxids` |
| `gene_neighbours.tsv`, `gene_positions.tsv` | optional: which genes lie next to which, per clade, and where each gene lies in each genome ([below](#gene-neighbours)); a run loads only the first |
| `gene_table.bin` | only in `database.protal`: `reference.map` and `unique_kmers.tsv` in binary, loaded on all threads without parsing (r226-sized tables, six threads: 1.76 → 0.43 s) |
| `column_weights.bin`, `strain_alleles.bin` | only in `database.protal` (since 2026-10-11), beside the text tables they are made from: the same tables as a run holds them (the column weights' codes packed, ~1.5 bytes a column), read straight into memory on several threads without parsing. At r226 the column weights' text took 12 s of every run on one thread and ~7.6 GB more while it was read ([report](claude/2026-10-11-performance-audit/README.md)). Each records the size of its text table; one made from another text (`--add_tables` of the text alone, by an older protal) is not used, and the text is read |
| `model_pe.xml`, `model_se.xml`, `model_PB.xml`, `model_ONT.xml` | the presence models of paired-end, single-end, PacBio and Nanopore reads ([below](#the-presence-model)); `model.xml` in databases before 0.7 |
| `model_pe_c1.xml`, `model_PB_c32.xml`, ... | since 2026-10-13: models trained at other sparse-first gates than the database's (`--add_model --model_gate`, [`--model-gates`](#choosing-the-sparse-first-gates)), which runs gated nearer that gate are profiled with |
| `model_pe_gbm.xml`, `model_pe_c1_gbm.xml` (or `..._forest.xml`), ... | since 2026-10-13: each of those models of the other kind (`..._gbm.xml` where the default is a forest, as builds make it since 2026-10-14; `..._forest.xml` beside boosted trees, as from 2026-10-13 to 2026-10-14; `--add_model --model_kind`), which runs given `--model_kind` of that kind are profiled with ([below](#kinds-of-models)) |
| `build_metadata.tsv`, `gene2geneid.tsv`, `genome2tiid.tsv` | optional (since 2026-10-13): notes no run reads, so that `database.protal` holds the whole database: what it was built from and with and its models' scores, the gene ids, the genomes by name ([below](#the-notes-a-database-carries)) |

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
removes the separate files it holds. The build inputs (`full_reference.fna`, ...) stay. The index
goes straight from memory into `database.protal` (since 2026-10-07; before, it was written as
`index.prx.zst`, read back, copied into `database.protal` and read again).

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
  `gene_table.bin`, `column_weights.bin` and `strain_alleles.bin` (beside `database.protal`, or into
  `--unpack_dir`); a folder reads the text tables.
- `--compress_db` on a `database.protal` without a current `gene_table.bin` (packed before 0.7.6), or
  without current binary forms of its column weights and strain alleles (packed before 2026-10-11, such as
  r226 v23), adds them, rewriting the file once with the other parts' frames copied. `--add_tables` of a
  `column_weights.tsv` or `strain_alleles.tsv` stores its binary form with it.
- `--compress_db` on a folder of raw files writes the same `database.protal` that `--build` would.

**Models.** `protal --add_model MODEL --read_type pe --db DB` (or `se`, `pb`, `ont`) checks a model
and replaces the database's model of that read type. Several at once
(`--add_model pe.xml,se.xml --read_type pe,se`) are all checked before anything is written.
`--model_gate C` (one per model, `-` for the database's own gate) stores a model trained on runs
gated at another gate than the database's as `model_<type>_c<C>.xml` (`--add_model pe.xml,pe_c1.xml
--read_type pe,pe --model_gate -,1`): runs gated at or nearest C are profiled with it
([below](#choosing-the-sparse-first-gates)).
`--model_kind KIND` (one per model, `-` for the read type's default file) stores a model of another kind
beside the default one ([below](#kinds-of-models)), and checks that the file is of that kind.
`--add_tables FILE` with `--add_model` stores the tables in the same write (the GTDB build stores its
models and `build_metadata.tsv` so).
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

### Kinds of models

Since 2026-10-13 a database may hold its models in two kinds: gradient-boosted trees (`gbm`) and a
random forest (`forest`). The read types' files (`model_pe.xml`, `model_pe_c1.xml`, ...) hold the
database's default kind; a model of the other kind sits beside each as `model_pe_forest.xml`,
`model_pe_c1_forest.xml`, ... (or `..._gbm.xml` where the default is a forest). A GTDB build trains
both at every gate it trains (`build_gtdb_database.py --model forest,gbm`, the default since 2026-10-14): the
random forest is the default again, the boosted trees go beside it (`gbm,forest`, boosting the default, from
2026-10-13 to 2026-10-14). A run picks the kind with `--model_kind`:

```bash
protal --db DB -1 r1.fq.gz -2 r2.fq.gz -o out --model_kind forest            # every read type
protal --db DB --profile_only 'out/*.sam.zst' -o out_forest --model_kind pe:forest   # again, pe by forests
```

Each sample takes the model of that kind trained at the gate nearest its own (`--sparse_first`), as
without `--model_kind` it takes the default files; a default file counts as a model of its own kind
(protal tells the kind from how the PMML combines its trees), so a database of forests alone, older
ones too, serves `--model_kind forest`. A run whose read type has no model of the kind stops (exit 2).
The model's log line names its kind, and without `--model_kind` the other kind the database also has.
`--profile_only` with `--model_kind` re-profiles a run's SAMs with the other kind without aligning
again. The build scores both kinds side by side in `model_logs/summary.txt` (`pe forest`, `pe c1
forest`), and `build_metadata.tsv` gives the forests' scores (`model_pe_forest`).

### The notes a database carries

`database.protal` holds every file a run reads, and since 2026-10-13 also the notes no run reads, so
that the one file is the whole database: `build_metadata.tsv` (what `build_gtdb_database.py` built it
from and with, its models' scores), and the converter's `gene2geneid.tsv` (marker ids to protal gene
ids) and `genome2tiid.tsv` (genomes by name). `--build` and `--compress_db` pack them where the folder
has them and, unlike the other files, leave them there, where scripts read them. `--add_tables` stores
them as they are (text, checked line by line); `build_metadata.tsv` is the last member, after the
models, so the build's last step (`--add_model ... --add_tables build_metadata.tsv`, one write) replaces
both in place. `--unpack_db` writes them back.

### The database in memory

- **Genes** are held at two bits per base. `N` is stored as `A`, an IUPAC code as the first base it
  stands for (A C G T order), so alignments, SNP and MSA reference rows show the stored base. `--build`
  leaves k-mers with an ambiguous base out of the index. A SAM of an older protal with an `M` against
  such a base is set aside (`misc/<sample>.err`).
- **The index** is held packed: each value takes the bits the database needs, its taxon and gene as
  one number (41 instead of 64 at r226: 26.6 GB instead of 35 GB), by a run and, since 2026-10-05, by
  `--build` too. The files are unchanged.
  A run and a build print `Index in memory: ...`. A run whose every sample is gated with the grid pick loads
  only the entries of its smallest gate's grid level (`--sparse_index auto`, since 2026-10-13; about 1/C of the
  values), and says so on that line ("a sparse index: the entries of grid level L or more, N of the file's M, on K
  of its J cores with entries").
  The file's key map gives each of the 2^30 cores an offset (3.2 GB whatever the reference; r226 uses 9.2% of its
  slots, [report](claude/2026-10-13-compact-key-map.md)). A run keeps it for the cores with values alone (the compact
  key map, the default since 2026-10-14; presence bits with ranks, 153 MB, and 3 bytes a used core: ~0.45 GB at
  r226, less on a sparse index), built as the index loads; the same lookups, a second cache line each
  (`--direct_keys`: the direct one, [running.md](running.md#developer-options)). A sparse or compact load decodes no
  key map of the file's: the decoder names each core with values (since 2026-10-14; a grid's sparse index loads
  2.5-3x faster than before, [report](claude/2026-10-14-sparse-load-key-map/README.md)).
- A full r226 run peaked at 36.5 GB, about 31 GB expected now ([running.md](running.md#memory)).
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
| memory | each index build ~35-37 GB with 64 threads (estimated; it was 64 GB before 2026-10-05; its log's `Memory after ...` lines say), up to ~75 GB while the finished database's build runs beside the profiling; `--one-build-at-a-time` needs about half |
| node-local disk | for the simulated samples (`--scratch`): the run estimates what it needs there and streams into protal only what would not fit ([local scratch](#local-scratch)); with streaming ~200-250 GB is enough, ~1 TB holds every read at once (~825 GB estimated with the default scenarios of 10 + 4 samples, ~670 GB real) |
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
- **Genomes to simulate from**, from NCBI: up to 40,000 species with up to 4 non-representative genomes
  each (since 2026-10-11: every species of GTDB r226 with genomes passing the filters, 36,515, about 80,700
  strains), and 9,000 species with their representative only (16,000 with 2 each until 2026-10-11, 6,000 and 2,000
  until 2026-10-05). The build leaves about half of each species' strains to its alleles, so with 4 a species keeps a
  strain to simulate from 15 times in 16: real strains are ~54% of the simulated strains and ~30% of the simulated
  species with the build's defaults, against 25% and 13% in r226 v23 ([report](claude/2026-10-11-r226-v23/README.md),
  `strain_supply.py` there estimates it from the metadata). Species are picked per domain in GTDB's proportions, in a random
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
| `--species`, `--per_species`, `--rep_only_species` | 40000, 4, 9000 | species with strains, strains each, species simulated from their representative only (and an in-silico strain) |
| `--min_completeness`, `--max_contamination` | 90, 5 | CheckM2 filters for strains |
| `--rep_genomes` | `ncbi` | where the representatives' genomes come from; `gtdb`: GTDB's archive of all of them (137 GB for r226) |
| `--tech_candidates`, `--no_tech_lookup` | 30 | candidates per species whose sequencing technology is asked of NCBI; or ask for none |
| `--progenomes` | | a proGenomes ANI-clustering table ([`pg4_ANI_clustering.tsv.gz`](https://progenomes.embl.de/download.cgi)) or any list of accessions: strains only from these |
| `--no_genomes`, `--dry_run` | | GTDB's files only; list what would be fetched |
| `--no_msa` | | leave out GTDB-Tk's alignments of the representatives' marker proteins and GTDB's trees (2.2 GB for r226, kept packed), from which `--build` takes the column weights' columns and ancestral sequences (since 2026-10-10); without them it aligns the proteins itself |
| `--connections`, `--ftp_url` | 8, NCBI's | genomes fetched at once from NCBI's FTP server (0: only through `datasets`) |
| `--host_genome` | `human` | the host genome: `human`, `ACCESSION_ASSEMBLYNAME` of another NCBI assembly (e.g. `GCF_000001405.40_GRCh38.p14`), or `none` |
| `--mirror`, `--datasets`, `--batch`, `-t` | | GTDB server, NCBI CLI, genomes per `datasets` request (500), compression threads (8) |

#### 2. Build and train

`build_gtdb_database.py` runs these steps. Each logs to its own file in `OUTDIR/logs/` and prints one
numbered line when it starts and indented lines when it ends.

1. **Genome table** (`work/genomes.tsv`, `model_logs/genome_table.txt`): the genomes to simulate from
   and their lengths. The genomes that give strain alleles are left out first, so that only the genomes
   simulated are read. Each genome's length and contig names are kept in a cache across builds
   (`--genome-cache`, by default `genome_cache.tsv.gz` in the `--inputs` folder), so that a rebuild from
   the same inputs reads none of them; reading ~50,000 genomes from a network file system took 8-22 min
   at r226. The contig names also go to the reports' cache, so that the reports read no genome either.
2. **Conversion** (`convert.log`): the release into a database folder, in the background while the
   genome table is made (since 2026-10-09; when nothing converted is there to reuse). Once it has
   written `reference.fna` and `reference.map`, the genes are placed in every genome to record which
   genes are neighbours (`gene_neighbours.log`, [details](#gene-neighbours)), while the converter
   writes the full reference. The gene neighbours have a stage of their own: a changed genome table
   finds them again without converting again.
3. **In-silico strains** (`insilico_strains.log`): a mutated copy of the representative of species with one genome,
   so that these species are simulated from a strain too (since 2026-10-11 only as many as keep the real strains at
   half of the simulated strains or more, `--insilico-strains auto`), and of a quarter of the species with real
   strains (half from 2026-10-09 to 2026-10-11), a strain their alleles do not reach (`--insilico-multi`;
   [details](#training-data-like-real-samples)). Since 2026-10-09 the species clouds come first (`species_clouds.log`:
   every species' nearest congeners by their references), which give each strain its congener sites and steer the
   hold-out of step 4. Only the simulations read the strains, so they are made in the background, at a lower
   priority, beside step 4; the simulations wait for them.
4. **Training database** (`training_db_index.log`): the database without the species and clades held
   out (30% of the species, plus whole clades of every rank), built alone. Its files are derived from the
   converted release on `-t` processes (`training_db.log`; one process took 5:17 at r226 v19):
   `reference.fna` in pieces of whole records, the full reference frame by frame through its seek table,
   and the gene neighbours beside them. Once it is built, the
   finished database's build starts in the background at idle priority, and the samples' simulations
   (what to stream into protal chosen from the room on the samples' disk). The two databases' strain
   alleles differ: the training database's come from half of each species' genomes, those no sample
   is simulated from (`--allele-genome-share 0.5`), the finished database's from all of them
   (`--final-allele-genome-share 1`, since 2026-10-11). No model is trained or evaluated on the
   finished database; with no species held out the build still makes a training database for that.
5. **Training data** (`training_data.log`): 309 paired-end samples of 20-200 species each (3 read setups x 103), profiled
   as paired-end and as single-end reads, and 186 PacBio and 186 Nanopore samples of the same
   communities, and the hold-in samples of four scenarios of real studies (gut, soil, shallow soil,
   host-dominated; [below](#scenarios-kinds-of-studies)), all against the training database.
6. **Independent test set** (`test_data.log`): samples of another design, for an honest score, and
   the scenarios' hold-out samples.
7. **Models** (`classifier_training*.log`): one random forest per read type (`--model`; gradient-boosted
   trees from 2026-10-06 to 2026-10-14), trained in parallel on the default feature set and evaluated with
   rows, samples and species held out (`--features`, `--evaluation basic`, the build's defaults since
   2026-10-06); with `--features auto` each trainer chooses its set and the console says which won and why.
   Since 2026-10-13 the same models again in the other kind once the default ones are trained (`--model
   forest,gbm`, the default since 2026-10-14; `classifier_training*_gbm.log`, `trained_model*_gbm.*`), at
   every gate, for `protal --model_kind gbm` ([kinds of models](#kinds-of-models)); `--model forest` trains
   one kind.
8. **Parity check** (`parity.log`): protal scores each model exactly as the trainer does; the read
   types' checks side by side (`protal --profile_only`, no index loaded), those of both kinds.
9. **Packaging** (`final_package.log`): the models of both kinds go into the finished database, with
   `build_metadata.tsv`, in one write in place (the database's build packed placeholders of them). The reports
   of what the models' errors rest on (`trace_relatives.log`, `error_reads_<read type>.log`) start
   right after the models are trained, beside steps 8 and 9 and the wait for the finished database's
   build (since 2026-10-09).

Where the run spends its cores is in `logs/` only: `cpu_jobs.tsv` has a row per command that ended (its
wall-clock and CPU seconds, and the cores it kept busy on average), `cpu_timeline.tsv` the cores the run
kept busy every 15 s (from its cgroup, the SLURM job's on a cluster) with the step and the commands
running, and `protal_cpu_<collection>.tsv` each protal run of a collection, as a whole and by protal's own
stages and samples (`misc/cpu.tsv`, [running.md](running.md)). The share archive carries them.

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
- **Marker truth** (since 2026-10-12): F1 as protal calls with only the species whose marker genes
  got reads counted as present. The simulator counts each genome's reads that overlap a marker gene
  (`simulate_metagenomes --gene_positions`, the release's `gene_positions.tsv`; the manifests'
  `marker_read_pairs` and `marker_first_reads`), the training tables carry them as
  `meta_marker_reads`, and a species of the training database with such reads but no row is in
  `<table>.unseen.tsv` and counts as missed. The species held out of the training database are not:
  protal cannot report them (until 2026-10-12 they counted, which put r226 v25's marker-truth F1 far
  too low, [report](claude/2026-10-12-r226-v25/README.md)). At r226 v24b 60-76% of the false negatives had at most
  two reads of their own on the taxon ([report](claude/2026-10-12-r226-v24b-errors/README.md)): this
  F1 tells the method's errors from what the reads cannot show. `summary.txt` has its rows
  (`marker truth: ...`).

[The presence model](#the-presence-model) explains the reports; [features.md](features.md) lists
what the models use.

Since 2026-10-08 the build also checks the samples' composition
([running.md](running.md#what-the-called-species-explain-the-unknown-share)) against the simulator's
truth, with the trained models' calls (`scripts/composition_accuracy.py`, a step after the parity
check). For every training and test sample it computes the composition again from the depth and
genome size of each taxon in the sample's `.profile.log` over the species the model calls
(`trained_model*.calls.tsv.gz`), and compares it with the community the sample was drawn from:

| Estimate | Against |
|---|---|
| the share of the reads explained | the share of the read pairs (drawn reads: bases) of the species called and present, and of all the database's species |
| the unknown share `?` | the share of the cells of the species not called (correctly), and of the species the database lacks |
| the average genome size | the cell-weighted size of the species called and present, and of all the community's cells |
| a species' depth and genome size | the depth of its simulated reads, and the length of its simulated genomes |
| the species missing at the median depth | the species present and not called |

`model_logs/composition_accuracy.tsv` (`_se`, `_pb`, `_ont`) has a line per sample;
`composition_accuracy.txt` and `summary.txt` the median and 10-90% range of each error, by set and
read type, samples with a host's reads apart (the unknown share counts them as genomes of the
called species' average size, so it is high there by design); the console a line of medians per
read type. A sample's host share is its reads beyond its community's (a drawn sample's that of the
paired-end sample it replays). The script runs on its own too, on a collection of
`collect_training_data.py`, with the profiles' own calls without `--calls`.

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
background build gets its line when it ends, whatever step the run is at. During the long steps,
indented lines also say each protal run of a collection as it starts (its design points and samples,
and with `--test-c-gates` how many samples at which other gates) and ends, and each trainer's main
evaluations as they begin (cross-validation, the test set, the tables at other gates, the scenarios,
the marker truth). With `--progress-every N`
each running stage also prints a status line every N seconds (its run time, memory and the last line
of its log). The run ends with `Ready protal database: ...` and the path of the model reports.

#### Outputs

`--outdir` holds four folders and the console's lines, each file once. Only `protal_db/` is needed to
use the database. Before 2026-10-08 the logs, the trainer's outputs, the genome tables and the stage
keys lay beside `protal_db/`, copies of many of them in `model_logs/`, and the in-silico strains and
one table per sample of the error reads made a build of about 12,000 files (r226 v17); now about 120
with `--scratch`, and about 1,400 with `--share-logs`, whose error-read SAMs are a file per sample
([report](claude/2026-10-08-build-outputs.md)).

| Path in `--outdir` | |
|---|---|
| `protal_db/` | the database: `database.protal` with its models of both kinds; `build_metadata.tsv`, the release, protal version and commit, command, design, held-out species and each model's scores (the versions and the scripts' commit are those the run started with, since 2026-10-06, followed by "at the end of the run: ..." if a pull or a rebuild changed them meanwhile, which the console then warns of); and the converter's `gene2geneid.tsv` and `genome2tiid.tsv`, which name its genes and genomes. Since 2026-10-13 `database.protal` holds these three too ([the notes](#the-notes-a-database-carries)): it alone is the database |
| `model_logs/` | everything to judge the models, and what the run chose: `summary.txt`; each read type's model and trainer's outputs, `trained_model*` ([the presence model](#training)): the model (`.xml`, as in the database; `.joblib`, scikit-learn's), the report (`.report.txt`, `.metrics.json`), predictions (`.scenario_predictions.tsv.gz`: the scenarios' hold-out samples; `.calls.tsv.gz`: every row's call), thresholds, feature importances; `parity.txt`, every read type's parity check; `composition_accuracy*.tsv` and `.txt` (the samples' composition against the truth, [above](#3-check-the-result)) and `composition_species*.tsv.gz` (every sample's species present or called: protal's depth and genome size beside the truth); `genome_table.txt`, the genomes simulated from, with the in-silico strains; `heldout_species.txt`, the species the training database leaves out, with the rank they were held out at and the clade, since 2026-10-08 also the distance to the nearest kept congener (`-`: none within 0.15) and the species complex held out with it (`c<number>`), and `holdout.txt`, their counts; `species_clouds.tsv`, every species' nearest congeners by their references' marker genes (`protal --write_species_neighbours` on the converted release, `species_neighbours.tsv`'s format), which steer the hold-out and give the training table's `meta_novel_distance`; with a gene subset `gene_subset.txt` and `gene_ranking.tsv`; what the conservation features rest on: `gene_congeners.tsv`, `gene_incongruence.tsv` (the build's reports, moved out of `protal_db/`), `relatives_by_gene_conservation.txt`; `error_reads/`, where the reads behind each model's errors in every sample went, a table per read type (with `--share-logs` their SAM records too; [below](#the-reads-behind-the-errors)), and with `--share-logs` `ancestry_sites/`, which side those reads take where the species differs from its congeners |
| `logs/` | one log per step (a step of several commands, such as the foreign scan and `--add_tables`, in one log, each command after a `---` line; the read types' parity checks, run side by side, each after a `==> parity_<read type> <==` line), and with `--scratch` every protal run's log of a collection one after the other (`protal_runs_training.log`, `protal_runs_test.log`); the run's CPU use: `cpu_jobs.tsv`, `cpu_timeline.tsv`, `protal_cpu_<collection>.tsv` ([above](#2-build-and-train)) |
| `console.log` | the console's lines, of every run into the folder |
| `work/` | what a rerun reuses, not needed to use the database: `stages/` (what a rerun may skip), `genomes.tsv` (the genomes simulated from: accession, taxonomy, FASTA, length) and `genomes_simulated.tsv` (with the in-silico strains), `internal_taxonomy.dmp`, `gene2geneid.tsv`, `foreign_rates/` (the finished database's scan, stored in it), and `training/`, `test/`: one table per read type (`training_data.tsv` for pe, `_se`, `_pb`, `_ont`), to retrain a model. Without `--scratch` also the samples, their profiles, the training database, the genome store and the in-silico strains (with `--scratch` these are on the scratch disk) |
| `<name>_share.tar.gz` | with `--share-logs`: the logs, `model_logs/` and the tables, to copy off the cluster ([below](#logs-to-share)) |

The in-silico strains (`insilico_strains/`: a FASTA per strain, 9,030 at r226 v17, and
`insilico_strains.tsv`, per strain the divergence drawn and reached and the substitutions) are made
on the samples' disk, `--scratch` or `work/`, beside the genome store; a new `--scratch` makes them
again (the same strains: they are seeded). An `--outdir` of the layout before 2026-10-08 is not
resumed: a rerun into it does its steps again.

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
  species gets a strain (`scripts/insilico_strains.py`; `--insilico-strains`, the share given one; default `auto`
  since 2026-10-11: as many as keep the real strains at `--min-real-strains`, 0.5, of the simulated strains or more;
  until then 1, and in r226 v23 the in-silico strains were 40% of the simulated species against 13.3% real ones):
  - a copy of the representative with substitutions only;
  - as far from it as the table's real strains are from theirs (marker genes from the k-mer traces
    in `gene_positions.tsv`), each marker gene by its conservation factor;
  - the genome at that divergence / 0.45, at most 5% (95% ANI; `--insilico-ani MIN-MAX` draws the
    ANI instead);
  - codon-aware: stop codons are never made, and an amino-acid change is kept with a probability
    (dN/dS) calibrated on the table's real strains, so that the in-silico strains put the same share
    of their differences on third codon positions as real strains do (0.15 without real strains;
    `--omega` in the script). The r226 v15 in-silico strains, made with 0.15, had 0.06-0.10 less of
    their differences on third positions than real strains at the same identity, a spectrum the
    presence models could learn ([report](claude/2026-10-07-error-read-signatures/README.md)).
  - congener sites (since 2026-10-09): as many of its marker substitutions at the sites where its species' nearest
    congener differs (`species_clouds.tsv`), with the congener's base, as the table's real strains have there, each
    in place of another of the gene's substitutions (`--congener-share` in the script). Drawn by gene and codon alone,
    the r226 v19 in-silico strains lacked their species' base at 3.5% of those sites against 10% for the real strains,
    so the models learnt from them that a strain never sides with a congener
    ([report](claude/2026-10-09-r226-v19/README.md), section 7).
  - strains the alleles do not reach (since 2026-10-09 evening): a quarter of the species with real strains in the
    table (half until 2026-10-11) get one as well (`--insilico-multi`, the share; `--multi-share` in the script). A species' other genomes are
    split by a hash into its strain alleles and the genomes simulated from, so every simulated real strain has a
    close relative among the alleles, and the models learnt that a strain the alleles explain is present and one
    they do not is not: in a world grown on known trees they called the strains of lineages without an allele
    genome at 0.3-0.6 against 0.98 for the others
    ([report](claude/2026-10-09-ancestry-true-positive-test/README.md)). A mutated copy of the representative
    shares none of the alleles' variants, as a strain of a lineage GTDB has not sampled does not.

  The training table marks taxa simulated from one (`meta_insilico_strain`), and the report lists
  them apart.
- **Species the database lacks.** The samples are profiled against a training database that leaves
  species out. Their reads land on relatives, as reads of organisms GTDB lacks do in real samples.
  Two kinds are held out:
  - whole clades (`--holdout-clades`: 2 phyla, 4 classes, 6 orders, 8 families, 12 genera, each at
    most 2% of the species);
  - then 30% of the remaining species to simulate (`--holdout`), species complexes whole: before
    anything is built, `protal --write_species_neighbours` compares every species' reference with
    its congeners' (`species_clouds.tsv`, what `--build` stores as `species_neighbours.tsv`), and
    congeners within `--holdout-complex-distance` (0.01, a strain's distance) are held out or kept
    together. Torn, a held-out twin's reads teach "absent" at the identity where a divergent strain
    teaches "present", and no feature can satisfy both
    ([report](claude/2026-10-08-r226-v17/README.md)).

  Every sample holds one species of a held-out clade of each rank. `model_logs/heldout_species.txt` lists them
  (with the distance to the nearest kept congener and the complex), `model_logs/holdout.txt` counts
  them by that distance, the training table carries it per row (`meta_novel_distance`), and the
  report gives false positive and false negative rates by rank and by that distance
  ([below](#species-and-clades-the-database-lacks)); `machine_learning_cmdline.py --twin-weight`
  can weigh down the absent rows beside a held-out twin (1 by default). The finished
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
| `pe` | Illumina paired-end reads by `simulate_metagenomes`'s instrument models ([below](#illumina-reads)): 75, 100, 150 and 250 bp (`--read-setups`; 100, 150 and 250 before 2026-10-14) |
| `se` | the first mates of the Illumina setups of `--se-read-setups` (150 and 250 bp), profiled as single-end, and Ultima reads (`--se-ultima-setup`, 300 bp, SD 50, cut at 320 bases; a unit a depth replaying its communities, since 2026-10-14; the scenarios' before) |
| `pb` | PacBio HiFi reads by `hifi_reads.py`'s model (`--pb-setup`: 15 kb, quality by length, Q50 at 5 kb to Q20 at 50 kb) |
| `ont` | Nanopore reads by [pbsim3](https://github.com/yukiteruono/pbsim3)'s quality-score model (`--ont-setup`: the high-quality model at 97%, 8 kb) |

The long-read samples (`--long-read-bases`, 300 kb to 6 Gb; `--long-read-samples` 36 per point) replay
the paired-end communities of the same depth. `simulate_metagenomes --long_samples` makes them in
process, one run per design point: it draws the reads' templates and writes each read as it is made,
compressed, with no template or pbsim3 files on the disk. Once the training database is built, the samples are
profiled as they are simulated, in protal runs of at least `--profile-blocks` GB of reads (20) and at
most `--profile-block-max` (200; the training data's and the test set's runs take turns), and each
design point's reads are removed once every read type that reads them is profiled (the SAMs, profiles
and dumps stay); `--profile-blocks 0` profiles both collections in one protal run once all is
simulated, and keeps the reads. A design point whose largest sample would take more than
`--stream-above` GB of compressed reads (by default chosen from the room on the samples' disk,
[local scratch](#local-scratch)) is not written to the disk at all: a protal run reads its samples,
with those of the other points streamed whose communities are there, from named pipes while
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

#### Choosing the sparse-first gates

protal's `--sparse_first C` gives a read up before seeding when none of the 1/C of its k-mer cores
in the sample finds a near-exact hit ([running.md](running.md#developer-options)); the gate's pick
(`--sparse_first_pick`) chooses those cores: `stride`, every C-th core of a read in read order, or
`hash`, those whose hash falls in 1/C of its range.
A read given up counts no failed candidates, so the models must be trained on runs gated as the
database's runs will be: the build records `--sparse-first` (default `pe:4,se:4,pb:16,ont:8`;
`pe:1,se:1,pb:8,ont:1` on 2026-10-12/13) and `--sparse-first-pick` (default `stride`; hash before
2026-10-13) in both databases' `build_settings.tsv`, and every run on them, the collections' too,
follows them. At r226 (v27) these gates cost no measurable F1 against C 1 (pb C 8) and saved 36-44%
of the alignment's CPU (pb 11%), and the stride lost 2.5-17 times fewer present species than the
hash ([report](claude/2026-10-13-r226-v27-c-gates.md)).

**Models at other gates.** A user may gate otherwise (`protal --sparse_first`): a sensitive gate,
or a faster one. `--model-gates` (default `pe:1/4/8,se:1/4/8,pb:8/16/32,ont:1/8/16`; `none` for the
database's gate alone) trains a model at each other gate of a read type too: both collections
profile all their samples at those gates as well (map rows `<sample>__c<C>`; a simulation streamed
into protal is streamed again for each, by a simulator of its own into pipes of its own, the same
reads), a trainer per gate takes the tables in `training/c_gates/c<C>/` and tests on
`test/c_gates/c<C>/` (`model_logs/trained_model[_<type>]_c<C>.*`, `classifier_training*_c<C>.log`),
and the models go into the finished database beside its gate's, as `model_pe_c1.xml`,
`model_PB_c32.xml`, ... (`protal --add_model --model_gate C`; its build packs placeholders of them,
which `--add_model` replaces in place). A run is profiled with the model of the gate nearest its own
(on the scale of log2 C, off next to C 1; between two, the lower). Each other gate costs the
collections' alignment once more, at that gate's speed, and a trainer; the console's
"models at other gates" line and `summary.txt` (rows `pe c1`, ...) give their F1.

**Comparing gates.** `--test-c-gates` profiles the test set (the independent test samples and the
scenarios' hold-out samples) at more gates, in the same protal runs (map rows `<sample>__c<C>` with
the map's `SPARSE_FIRST` column; samples streamed through pipes only at the database's gate and those
of `--model-gates`). Every level of the test set, those of `--model-gates` too, is scored as protal
would call it: by the model of the gate nearest it (`--extra-test-file` of that model's trainer: the
report's "Extra test tables" section and `extra_tests` in `<prefix>.metrics.json`), and by the
database gate's model beside the test file's rows of the same samples. The present species of a
table's samples that have no row in it (`<table>.unseen.tsv`, beside it: a gate that gives up every
read of a species leaves it no row) count as missed; before 2026-10-13 they were left out, so that
a gate losing more species looked better ([report](claude/2026-10-13-r226-v27-c-gates.md)). The build
writes `model_logs/c_gates.tsv` and a console line per read type:

| Column | Holds |
|---|---|
| `read_type`, `gate`, `pick`, `flex_indels`, `database_gate`, `model_gate` | the read type, the gate (0: off), its pick (`stride` or `hash`; `-` below 2, where both look up every core), 1 with `--flex_indels`, 1 for the database's gate (the test file itself), and the gate of the model that scored it |
| `samples`, `taxa`, `present`, `unseen`, `TP`, `FP`, `FN`, `precision`, `sensitivity`, `F1` | the model's calls on those samples at that gate, `unseen` present species without a row counted in `present` and `FN` |
| `aligning_cpu_s`, `aligning_wall_s` | the samples' alignment, summed over the protal runs' `aligning <sample>` rows (`protal_cpu_<collection>.tsv`; empty when a sample has none) |
| `baseline_F1`, `baseline_aligning_cpu_s` | the same samples at the database's gate without `--flex_indels`, scored by its model, beside every other row (for se's rows with `--flex_indels`, the Ultima samples alone) |

A gate without a model of its own is scored by a model trained at another: a gate whose F1 holds even
so is a safe choice, while one that loses F1 may recover with a model trained at it (add it to
`--model-gates`). The wall seconds of rows profiled side by side in one protal run say little; compare
the CPU seconds. The test set is small (its F1 moves by a few thousandths between seeds: compare the
FP and FN counts too).

With `--test-c-gates` the test samples of ONT, at every gate (the database's C 1 too), and the
Ultima (single-end) ones, at se's gates, are also profiled with protal's `--flex_indels`
(`--test-flex-indels`, default `ont:0/1/4/8,ultima:0/1/2/4`; `none` for none): map rows
`<sample>__c<C>f` with the map's `FLEX_INDELS` column, tables in `c_gates/c<C>f/`, rows of
`c_gates.tsv` with `flex_indels` 1. The test set's Ultima simulations are then written to the disk
rather than streamed (at r226 ~50 GB, all three above the stream threshold; if they do not fit they are
streamed, with a warning, and have no such rows). It needs a protal whose `--full_help` lists
`--flex_indels` and `FLEX_INDELS`; an older one would ignore the column, so the build leaves these
rows out, with a warning. The models are trained without `--flex_indels`, as with other gates.

With `--test-c-gates` every read type's test samples are also profiled at its database's gate (if 2
or more) with the pick the database does not use (`--test-other-pick`, `--test-stride` before
2026-10-13; default `auto`; given gates, or `none`): with the default stride pick, the hash's, the
cores its hash puts in the sample (map rows `<sample>__c<C>h` with the map's `SPARSE_FIRST_PICK`
column `hash`, tables in `c_gates/c<C>h/`); with `--sparse-first-pick hash`, the stride's (rows
`<sample>__c<C>s`). Their rows of `c_gates.tsv` follow the same gate's row with the database's
pick, whose F1 the console line gives beside them. It needs a protal whose `--full_help` lists
`--sparse_first_pick` and `SPARSE_FIRST_PICK`; an older one would ignore the column, so the build
leaves these rows out, with a warning.

With `--test-c-gates` every read type's test samples are also profiled with the grid pick
(`--test-grid`, since 2026-10-13; default `auto`: each gate of `--model-gates` and the database's of 2
or more; given gates, or `none`): map rows `<sample>__c<C>g` with `SPARSE_FIRST_PICK` `grid`, tables in
`c_gates/c<C>g/`, scored by the model of the gate nearest it (trained with the database's pick), and in
`c_gates.tsv` after the same gate's row with the database's pick. These rows are profiled in the
collection's protal runs, on the whole index, seeded as the sparse index seeds them (protal
`--sparse_index emulate`, which every collection run gets: it changes only grid-picked samples, so a
database picked by grid is trained on rows seeded so too), so their records and calls are a sparse
index's but their CPU seconds are not; `scripts/measure_performance.sh --test_sparse` with levels such
as `8g` times the sparse index itself. It needs a protal whose `--full_help` lists `--sparse_index` and
`emulate`; with an older one the build leaves the rows
out, with a warning.

The `--test-c-gates` default (`pe:0/2/16,se:0/2/16,pb:0/4,ont:0/4/32`, since 2026-10-13) adds to
the gates of `--model-gates` the ungated runs, gates between the ladder's and beyond its fast end.
Each sample written to the disk is aligned once per way; streamed samples only at the database's
gate and those of `--model-gates`. The estimate of the defaults before 2026-10-13, 33 ways, ~20 CPU
hours at r226: [c-gates cost](claude/2026-10-13-c-gates-cost.md).

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
`setup` (`ultima:LENGTH_MEAN:LENGTH_SD:Q_MEAN:Q_SD[:MAX]`, default `ultima:300:40:25:2`, the reads cut at MAX bases, 320 unless given), `pb` and `ont` with
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
- **Ultima reads**: single-end, length from a gamma distribution cut at 320 bases (a run's flows end its
  reads; `ultima:...:MAX` sets another), errors indels, mostly homopolymer length errors (homopolymers from
  two bases), with substitutions rare at Q35 to Q40 (drawn per read), base qualities averaging the read's
  quality and falling over its last 100 bases to 15 below it at its end (both since 2026-10-14), made by
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
or unseen species, wherever they went (`source:`). It writes to `model_logs/error_reads/<read type>/`:

| File | |
|---|---|
| `taxa.tsv.gz` | every sample's error taxa, a line each, with the sample, its set (`training` or `test`), design point and scenario: FP, FN, or `unseen` (a species of the sample that the training database has but that protal profiled no reads to: never scored, so not in the models' counts, but in protal's; not the species held out of it, which its `genome2tiid.tsv` keeps without genes: the build passes `heldout_species.txt`, `--heldout`); its score and knob; for FN and unseen species their genomes (and read pairs simulated, paired-end), their fragments with a record and where their best records went (on the taxon, at MAPQ 4 or more as the profiler counts them, elsewhere and on which taxa, or nowhere) and those that seeded on the species but did not align to it (`own_seeded_not_aligned`); for all, the fragments on the taxon and their sources, and those that seeded on it but failed to align. It counts every fragment. Before 2026-10-08 a file per sample, `<training\|test>/<design point>/<sample>.taxa.tsv` (about 1,500 at r226 v17) |
| `summary.tsv` | one line per sample (the console sums it per read type) |
| `<training\|test>/<design point>/<sample>.FP.sam.zst`, `.FN.sam.zst` | with `--share-logs`: the records of the reads taken for the false positives, and for the false negatives (`FN:`, and `seeded:` and `source:` of a false negative), of at most 20 fragments per taxon and reason (the same ones at any cap: the lowest CRC-32 of the read name); each record with QUAL left out (`*`), its source genome and species (`xg:Z:`, `xs:Z:`, "(not in the database)" for a species the training database lacks) and why it was taken (`xe:Z:`, its reasons); the `@SQ` lines cut to the genes they name |

With `--share-logs` the build
then runs `scripts/ancestry_sites.py` on every read type's SAMs, in one run that reads the references once and
analyses the read types side by side (`model_logs/ancestry_sites/<read type>.*`, `logs/ancestry_sites.log`; a log per
read type before 2026-10-09): for every false positive and false negative, the sites where the
species' gene copies differ from its nearest congener's in the training database and, from the species'
other genomes in the build's full reference (kept until then), the sites where its own strains vary, a read's
own source genome left out of them (a missed real strain was simulated from a GTDB genome that is usually among
them, and with it in, its every difference was polymorphic by construction: the r226 v18 report); per
counted read, the sites it covers and whether it has the species' base or the congener's there. The
`.summary.txt` and `.auc.tsv` say how well each signal tells the false negatives' own reads (a strain's)
from the false positives' (a novel congener's) per taxon, overall and within identity bands, against
plain identity; `.taxa.tsv.gz` and `.fragments.tsv.gz` hold the counts. The same sites are features of
the models since 0.7.9 (`ancestry_*`, [features.md](features.md)); this report adds the alleles'
half, which the database does not hold yet. Unseen species are counted
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

The build runs these reports once the models are trained, beside the parity check and the packaging; the
database does not wait for them. `error_reads.py` for each read type and `trace_relatives.py` run side by
side, each on all of `-t` (`logs/error_reads_<read type>.log`, `logs/trace_relatives.log`): they differ in
work (r226 v19 on 16 threads each: pe 186 s, se 150 s, PacBio and ONT ~80 s), and one that ends leaves its
cores to the others. The genomes' contig names come from the genome table's cache (step 1); only the
in-silico strains' are read, from the samples' disk.
`error_reads.py` extracts `--threads` samples at once, the largest SAMs first, while their estimated
memory (300 MB and twice the SAM's size on disk) fits in `--memory` GB (default 60% of the least of the
machine's memory, `SLURM_MEM_PER_NODE` and the process's cgroup limit; the build gives each read type's
run its share of that), and keeps only what the tables
look up per fragment (its reasons, and of the taxa it aligned or failed to align to only the sample's
error taxa). A worker killed from outside (out of memory) costs no other sample: those it ran beside run
again one at a time, a sample killed again is named and left out. At r226 v15 the paired-end extraction
was killed (84 workers started on the 18 soil samples at once, each tracking millions of fragments at
1.2-2.3 kB with every taxon kept and the held-out species' reads counted), and its read type lost every
sample's tables ([report](claude/2026-10-07-r226-v15/README.md)).

#### Logs to share

`--share-logs` keeps the error reads' SAMs ([above](#the-reads-behind-the-errors)) and ends the build by
packing `OUTDIR/<name>_share.tar.gz` (`<name>`: the last part of `--outdir`), to copy off the cluster and
unpack elsewhere (`tar xzf`; everything is under `<name>/`, where it is in `--outdir`):

| Path | |
|---|---|
| `console.log`, `logs/` | the console's lines (every run of the build, after a line with its time and command) and each step's log |
| `logs/simulations_<collection>.log` | from the samples' disk, which the end of a cluster job may clear: each design point's simulator logs (`simulate.log`, `design.log`, `stream*.log`) and parameters (each `run_params.tsv`), one after the other under a `==> <point>/<file> <==` line (a log over 1 MB keeps its first and last 256 kB); and `logs/protal_runs_<collection>.log` where the build did not gather the protal runs' logs there (without `--scratch`) |
| `model_logs/` | as [above](#outputs), with the error reads' tables and SAMs; not the models (`.xml`, `.joblib`), which are in the database; and `insilico_strains.tsv` from the samples' disk (per in-silico strain: its divergence drawn and reached, its substitutions, its congener and those placed at its sites) |
| `protal_db/build_metadata.tsv`, `work/internal_taxonomy.dmp` | what the database was built from and with, and its taxonomy |
| `work/training/`, `work/test/` | the training and test tables of each read type, their numbers to 9 significant digits instead of Python's 17, at about 72% of the gzipped size (v15's pe table). That keeps every float32 value, so a forest (float32 comparisons) trains the same; boosting compares doubles, so its thresholds can move in the tenth digit, but values closer than that are rare and share a bin |

The database and the rest of `work/` are left out. The console's last lines give its size.

#### Local scratch

The simulations write and delete many files, which a network file system is slow at.
`--scratch DIR` makes the samples, the training database and the converter's temporary files on a
node's own disk (with the in-silico strains), and copies the tables to `--outdir/work`. The samples are profiled as they are simulated
and their reads removed once profiled (`--profile-blocks`), and a simulation that would leave less
than `--keep-free` GB (30) waits for that, so the disk holds the training database (~24 GB at r226),
what is simulated but not yet profiled, the SAMs and profiles, and the genome store (`--genome-store`,
~50 GB at r226, kept for the next build).

Before its builds the run estimates what it needs there besides the reads: the genome store's growth
(0.25 bytes a base of the genomes simulated from, less what the store holds), each database still to
be built on that disk (1.5 times its `reference.fna`: r226's training database is 21.8 GB from 15.4
GB), the host genome as plain sequence, the samples' SAMs and profiles (a tenth of the reads'
estimate; r226 v15 kept ~62 GB for 591 GB of reads; with `--test-c-gates` a tenth again for each other
gate a written point's test samples are profiled at, ~74 GB at r226 with the defaults,
[estimate](claude/2026-10-13-c-gates-cost.md)) and `--keep-free`; and what the reads of each
design point take (the collector's estimate, ~1.25 times what r226 v15 wrote: ~825 GB for both
collections at the defaults since 2026-10-07). It says both, and stops at once if even streaming every
point would not fit. Once the training database is built it measures the free space again and, with
`--stream-above auto` (the default), streams into protal only the largest points, until the reads of
the others, with their test samples' SAMs and profiles at other gates, fit at once. At the defaults
(116 points) with 1 TB free that streams none, with 700 GB the
9 with a sample above ~10 GB, with 500 GB the 21 above ~6 GB
([report](claude/2026-10-07-build-ordering/README.md)). Streaming saves the space but not time: a
streamed point's samples are aligned as they are made, and the profiling stage at the end of a protal
run is shared only by the points in that run. The streamed points whose communities are there share
one protal run (up to `--profile-block-max` GB until 2026-10-09; their reads are on no disk, and protal reads one
sample at a time, so the simulators it has not reached yet wait at their pipes). The console line `--stream-above auto: ...` says what it
chose, `build_metadata.tsv` (`samples_streamed`) keeps it, and the run's end gives the most it took on
the disk beside the estimate. Each protal run's log is copied to `OUTDIR/logs/protal_runs_<collection>.log`.
With `--profile-blocks 0` nothing is streamed and every sample is on the disk at once. A rerun reuses the
samples only from the same DIR, so on a disk that is cleared after the job, a rerun simulates again
(the builds are kept in `--outdir`).

#### Time and memory

Since 2026-10-07 the profiling is the build's longest stage: the simulations take minutes
(`simulate_metagenomes` makes every read type itself, the genomes from the store), protal hours
([report](claude/2026-10-07-build-ordering/README.md)). So the training database, which the profiling
waits for, is built first and alone, with every core and no simulation writing to its disk. The
simulations start once it is built, in the background at a lower priority than the profiling, long
reads and paired-end points in one queue, the longest first. The finished database is built in the
background after the training database, at the idle scheduling class (`SCHED_IDLE`): on the cores the
profiling leaves; it is only needed at the end, to take the models, and is paused while they are
trained (boosting's threads wait on each other at every step, so nothing else runs beside them). That
needs the memory of a build beside protal's. `--one-build-at-a-time` builds the finished database after
the training instead. The models go into the finished database right after the parity check; the
reports of what their errors rest on (`trace_relatives.py`, `error_reads.py`) follow, side by side.

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
  the converter, protal, the held-out species; recorded in `work/stages/`). The gene neighbours have a
  stage of their own since 2026-10-09 (they depend on the genome table), so an `--outdir` of an earlier
  version converts and builds once again. The training database on
  `--scratch` keeps that record in its own folder too (`built_for.json`): one that another
  `--outdir`'s build left there is built again. The folder of a reduced database (`--n-genes`,
  `--genes`) is never taken for the whole release by a later run without them.
- With the finished database kept, the in-silico strains of another `--seed` take the release
  converted again for the real strains' gene divergence (`gene_positions.tsv`, which the finished
  database packed), not a uniform ANI.
- The in-silico strains are kept when the genome table and options are the same and they are still on the samples' disk.
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
| `--profile-blocks`, `--profile-block-max`, `--keep-free` | 20, 200, 30 | profile the samples as they are simulated, in protal runs of at least this many GB of reads and at most that many of written reads (0: no limit; the streamed ones share a run, uncapped since 2026-10-09), removing the reads profiled (0: one protal run once all is simulated, reads kept); the GB a simulation leaves free on the samples' disk, or it waits |
| `--stream-above` | auto | with `--profile-blocks`: a design point whose largest sample would take more than this many GB of compressed reads is streamed into protal through named pipes, never written (0: none); `auto` (2 before 2026-10-07) streams only what the room on the samples' disk requires ([above](#local-scratch)) |
| `--profile-ahead`, `--no-profile-ahead` | on (since 2026-10-09) | protal `--profile_ahead` in the collections' runs: a sample profiled while the next is aligned, on a quarter of the threads, so that a run's profiling stage after its last alignment shrinks to what was left (r226 v19, without: 15.5 min over six runs); the profiles are the same |
| `--genome-cache` | `auto` | each genome FASTA's length and contig names by path, size and modification time, kept across builds: `auto` puts it in the `--inputs` folder (`genome_cache.tsv.gz`, where writable), else `OUTDIR/work/`; a path; or `none` |
| `--compressed-pipes` | off | with `--stream-above`: what goes through the pipes compressed as the files' names say, not as plain FASTQ |
| `--genome-store` | auto | `simulate_metagenomes --genome_store` for both collections: each genome simulated read and parsed from its FASTA once, written at 2 bits a base and memory-mapped by every later sample and simulation (at r226 ~1.5M genome reads of ~54k genomes otherwise); `auto`: `SCRATCH/genome_store` (`OUTDIR/work/genome_store` without `--scratch`), kept for the next build; a folder; or `none`. ~0.25 bytes a base of the genomes simulated, ~50 GB at r226 besides the samples ([report](claude/2026-10-07-simulate-metagenomes-audit/README.md)) |
| `--seed` | 1 | |
| `--holdout`, `--holdout-clades`, `--holdout-max-share` | 0.3, `phylum:2,class:4,order:6,family:8,genus:12`, 0.02 | species and clades left out of the training database; `--holdout-clades none` for species only |
| `--holdout-species` | | a file of the species to leave out instead |
| `--holdout-complex-distance`, `--species-clouds` | 0.01, | species within this distance of a congener on their references' marker genes form a complex that the training database leaves out or keeps whole, never torn (0: species drawn one by one, as before 2026-10-08); the distances come from `species_clouds.tsv`, which the build makes with `protal --write_species_neighbours` on the converted release, or from the table given ([below](#species-and-clades-the-database-lacks)) |
| `--novel-clades-per-sample` | 1 | species of held-out clades per sample and rank |
| `--insilico-strains`, `--min-real-strains`, `--insilico-multi`, `--insilico-ani` | `auto`, 0.5, 0.25, | share of one-genome species given an in-silico strain (`auto`: as many as keep the real strains at `--min-real-strains` of the simulated strains or more); share of the species with real strains given one too; or their ANI drawn from `MIN-MAX` |
| `--weighted-site-shift` | off | the collections' protal runs with the weighted site shift (`protal --weighted_site_shift`), and the database's users should profile with it too (`build_metadata.tsv` says so); the default from 2026-10-11 to 2026-10-12, off again since it lowered r226 v24b's short-read depths by 1-2% without a measurable F1 gain |
| `--samples` | 9 | samples per design point (12 before 2026-10-14, with three paired-end setups) |
| `--read-pairs` | `1000,2000,5000,20000,50000,100000,200000,500000,2000000:3,10000000:1,30000000:1` | depths, one design point each; `DEPTH:SAMPLES` for another number of samples (`2000000:4,10000000:2` before 2026-10-14) |
| `--read-setups` | `75:HS25:250:40,100:HS20:300:40,150:HSXt:350:50,250:MSv3:550:50` | the paired-end reads: read length : instrument ([Illumina reads](#illumina-reads)) : fragment mean : SD; 75 bp since 2026-10-14 |
| `--se-read-setups` | `150:HSXt:350:50,250:MSv3:550:50` | the single-end Illumina reads (the first reads of these setups' simulations; a setup in both lists is simulated once); since 2026-10-14, `--read-setups`' before |
| `--se-ultima-setup` | `ultima:300:50:25:2` | single-end Ultima reads of the design besides the Illumina ones, a unit a depth of that many reads and the depth's samples, drawn from the depth's communities ([Ultima reads](#scenarios-kinds-of-studies): cut at 320 bases, the last 100 poor, indels); since 2026-10-14; empty for none |
| `--species-per-sample`, `--archaea` | `20-200`, 2 | species per sample; archaeal species per sample |
| `--strains-per-species` | `0.5,0.2` | probabilities of a second, third, ... strain (the test set's; `0.3,0.1` before 2026-10-07) |
| `--congeners` | `0.25:2-5` | about this share of a sample's species in groups of 2-5 congeners; `0` for none |
| `--abundance` | `lognormal:1.3,2.0` | `lognormal:SIGMA[,SIGMA...]` (the samples take the sigmas in turn; a continuous long tail, no species below 1/1000 of the median, every species a read pair or more), `poisson_lognormal:SIGMA` (Poisson counts + 1, the lognormal before 2026-10-07: ~40-45% of the species at the lowest weight), `powerlaw:ALPHA` or `negbin:R:P` |
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
| `--features` | `normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood+ancestry+gaps+untried+alleles+polymorphic+weights+strain_identity+length_free+read_length` | the models' features (default since 2026-10-06, `auto` before: the set every model of the r226 v12 and v13 builds chose, the reference's k-mer uniqueness, the sample's complexity, which needs a protal of 2026-10-06 or later for the training data, and since 0.7.9 the three groups against the false positives of complex communities, the ancestry sites, and the congener gaps and untried candidates (`gaps`, `untried`, since the congener-gaps merge; `foreign` there too until 2026-10-08), in 0.8.0 the strain alleles and polymorphic sites (`alleles`, `polymorphic`), in 0.8.1 the column weights (`weights`), and since 2026-10-13 the identities against the best strain allele (`strain_identity`), and since 2026-10-14 `length_free` and `read_length` (features that do not change with the read length or a sparse index in place of those that do, and the sample's read length), which need a protal of 2026-10-14 or later (2026-10-13 for `strain_identity`, 0.8.1 for `weights`, 0.8.0 for `alleles` and `polymorphic`, 0.7.9 for the groups up to `ancestry`): `normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity` trains without them); `auto`: each trainer chooses its set ([below](#training)), which doubles a boosted model's training with `--evaluation basic`; or feature groups ([features.md](features.md)); `+priors` adds GTDB's species constants (opt-in since 0.7.6) |
| `--model` | `forest,gbm` | the kinds of models, comma-separated: `forest`, a random forest (the default again since 2026-10-14), `gbm`, gradient-boosted trees (the default from 2026-10-06 to 2026-10-14; [below](#training)); the first is the database's default (`model_pe.xml`, ...), the others are trained after it and stored beside it (`model_pe_gbm.xml`, ...: [kinds of models](#kinds-of-models)); both since 2026-10-13 (`gbm,forest` to 2026-10-14), `gbm` alone before |
| `--rounds`, `--ntree`, `--maxnodes` | 250, 64, `63` (`512,pb:128,ont:128` for forests) | boosting's rounds, a forest's trees, and leaves per tree by read type (`N` or `TYPE:N` items; `--maxnodes` for `--model`'s first kind, the others their defaults) |
| `--call-mode` | `curve` | `fdr` also stores calibrated calls at a target share of false calls ([below](#calls-at-a-target-share-of-false-calls)) |
| `--evaluation`, `--previous-procedure` | `basic`, off | how much the trainer evaluates: `basic` (default since 2026-10-06, `full` before), the models with rows, samples and species held out, which the summary and the knob need; `full` adds the clades held out and the studies ([below](#training)), several times a boosted model's training |
| `--n-genes`, `--genes`, `--gene-ranking`, `--genes-per-domain`, `--rank-genes` | | a reduced database ([below](#reduced-marker-sets)) |
| `--no-gene-neighbours` | | skip the gene neighbours; protal then pairs no mates across neighbouring genes |
| `--strain-alleles` | 16 | strain alleles per species and gene in both databases (`protal --build --strain_alleles`; 0: none; 4 until 2026-10-11, [report](claude/2026-10-11-strain-false-negatives/README.md)) |
| `--allele-genome-share` | 0.5 | the share of each species' genomes that give the **training database's** strain alleles, by a hash of the accession (`protal --allele_genome_share`). The genome table loses those (but the representatives), so that no simulated strain is its species' own allele, a gain no real sample would see; the in-silico strains are made after it, for the species left with one genome. The ancestry report takes its alleles from the same genomes |
| `--final-allele-genome-share` | 1 | the share of each species' genomes that give the **finished database's** strain alleles (since 2026-10-11; `--allele-genome-share` before): all of them, so that the database knows every strain GTDB has. No model is trained or evaluated on the finished database, whose alleles include the strains the samples are simulated from: the collections, the parity check and the reports read the training database, and a run with no species held out (`--holdout 0 --holdout-clades none`) builds a training database of its own when the two shares differ |
| `--index-alleles` | 1 (off) | the k-mers of the strain alleles at least this far from the representative's copy (edits per base) go into both databases' indexes under the species as non-unique entries, so that a read of a deep strain finds its species (`protal --build --index_alleles`, [below](#2-build-the-index); 1 or more: none, the default since the r226 v22 build; 0.01 was the default before) |
| `--fragment-count` | `mapq` | how the models' `fragments` feature counts (`protal --build --fragment_count`): recorded in both databases' `build_settings.tsv`, and followed by every run on them, the collections' included; `soft` adds each fragment the MAPQ and length filters drop by its read-EM share, so that thin divergent strains are not counted as nothing ([report](claude/2026-10-12-r226-v24b-errors/README.md), section 3; untested at r226) |
| `--wfa-traceback` | protal's (`singletrack`; `classic` before 0.8.2) | the WFA2 traceback of both databases (`protal --build --wfa_traceback`): recorded in their `build_settings.tsv`, and followed by the collection's protal runs, which align as the training database was built. `singletrack` gives the same scores but may choose otherwise among equally good alignments ([below](#2-build-the-index)) |
| `--sparse-first` | `pe:4,se:4,pb:16,ont:8` | each read type's sparse-first gate (`protal --build --sparse_first`; `C` for all, `none` for no gates; since 2026-10-12, `pe:1,se:1,pb:8,ont:1` until 2026-10-13): recorded in both databases' `build_settings.tsv` and followed by every run on them, the collections' included, so that the models are trained on runs gated as the database's runs are (a read given up counts no failed candidates). At r226 v27 these gates cost no measurable F1 against C 1 (pb C 8) for 36-44% less alignment CPU (pb 11%) ([report](claude/2026-10-13-r226-v27-c-gates.md); [earlier](claude/2026-10-12-sparse-first/README.md)) |
| `--sparse-first-pick` | `stride` | which 1/C of a read's cores each gate looks up first (`protal --build --sparse_first_pick`, recorded with the gates): `stride`, every C-th core in read order, `hash`, those whose hash falls in 1/C of its range (the only pick before 2026-10-13), or `grid`, the database's grid of each gene's cores (since 2026-10-13: runs gated at C 2 or more load a sparse index, [running.md](running.md#developer-options); measure it with `--test-grid` first); one for all, or read types and theirs (`pb:hash`) |
| `--model-gates` | `pe:1/4/8,se:1/4/8,pb:8/16/32,ont:1/8/16` | the gates each read type also gets a model of: both collections profiled at them too (streamed simulations streamed again), a model trained at each and stored as `model_<type>_c<C>.xml`, which runs gated nearer that gate are profiled with; `none`: the database gate's alone ([above](#choosing-the-sparse-first-gates)) |
| `--test-c-gates [GATES]` | off; without a value `pe:0/2/16,se:0/2/16,pb:0/4,ont:0/4/32` | the test samples written to the disk also profiled at these gates (the database's left out: its row in `c_gates.tsv` is the test set itself, profiled at it; those of `--model-gates` profiled anyway), in the same protal runs as `<sample>__c<C>`; each table scored by the model of the gate nearest it beside the test file's rows of the same samples, and `model_logs/c_gates.tsv` gives per read type and gate the F1, precision, sensitivity, FP and FN (present species without a row counted as missed) and the samples' alignment CPU and wall seconds ([above](#choosing-the-sparse-first-gates)); `c_gates.tsv` has the gates of `--model-gates` without it |
| `--test-other-pick [GATES]` (`--test-stride`) | with `--test-c-gates`: `auto`, each read type's database gate of 2 or more | the gates at which the test samples are also profiled with the pick the database does not use (`hash` rows `<sample>__c<C>h` beside a stride database); `none`: none ([above](#choosing-the-sparse-first-gates)) |
| `--test-grid [GATES]` | with `--test-c-gates`: `auto`, each read type's gates with a model of 2 or more | the gates at which the test samples are also profiled with the grid pick (rows `<sample>__c<C>g`, on the whole index: their calls are a sparse index's, their times not); `none`: none ([above](#choosing-the-sparse-first-gates)) |
| `--test-flex-indels [GATES]` | with `--test-c-gates`: `ont:0/1/4/8,ultima:0/1/2/4` | the gates at which the ONT test samples (key `ont`) and the Ultima ones (key `ultima`) are also profiled with protal's `--flex_indels`; `none`: none ([above](#choosing-the-sparse-first-gates)) |
| `--foreign-rates`, `--no-foreign-rates` | off | scan each database's full reference for the gene copies other species' reads reach (`scripts/foreign_rates.py`: the training database's right after its build, the finished database's before its models go in) and store `foreign_rates.tsv` in it, for the `foreign` features ([features.md](features.md#the-gene-copies-gaps-and-foreign-reads-and-the-untried-candidates-gaps-foreign-untried-2026-10-07)); without `--foreign-rates` (the default since the r226 v18 build, which measured 0.001 of AUC for two scans of ~21 min: [report](claude/2026-10-08-r226-v18/README.md)), or with `--no-foreign-rates`, which overrides it, there is no scan and the features are unknown (-1). From 2026-10-07 to 2026-10-08 the scan read the genomes the samples are drawn from, which told the models which species those are ([features.md](features.md#the-foreign-features-leak)) |
| `--foreign-stride`, `--foreign-per-header` | 250, 10 | the scan's reads: 150 bases every this many bases of a gene copy, of at most this many copies of each species' gene (the full reference's first; 0: every genome's) |
| `--one-build-at-a-time` | | build the finished database after the training, not beside the profiling |
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
`--genes LIST` names the genes. The choice is in `model_logs/gene_subset.txt`, `model_logs/gene_ranking.tsv` and
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
- `full_reference.fna.zst`: every genome's marker genes (`>taxid_geneid accession`, the accession
  since 2026-10-08), 86 GB uncompressed at r226; only if
  `genomic_files_all` is there. A zstd frame per marker file (per gene in a `--from_db` copy) and a
  seek table at the end, so that `protal --build` decompresses it on several threads; the `zstd`
  command reads it as any other file;
- `species_priors.tsv` ([below](#species-priors));
- `ncbi_taxids.tsv` ([below](#ncbi-taxids));
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
| `--priors_only` | off | write `species_priors.tsv` and `ncbi_taxids.tsv` alone into `--outdir`, touching nothing else there ([below](#species-priors)) |

#### Species priors

`species_priors.tsv` holds what GTDB knows of each species before any read, for the model's
optional prior features ([features.md](features.md#the-species-priors-priors-075-opt-in)):
- the marker genes found in the representative, and how many twice (always 0 for a GTDB release,
  whose marker files hold one copy per genome);
- the representative's CheckM completeness and contamination;
- from the species clusters file, the cluster's ANI radius, mean and minimum intra-species ANI and
  number of genomes.

Since 2026-10-08 it also holds each species' genome size, which no model feature uses: runs take it
for the sample's composition (how much of the sample the called species explain, its average
genome size, the profile's unknown share; [running.md](running.md#what-the-called-species-explain-the-unknown-share)):

| Column | |
|---|---|
| `genome_size` | the mean of the species' genomes' sizes in GTDB's metadata, each corrected for CheckM: assembly size × (100 − contamination) / completeness (Rodríguez-Gijón et al. 2022). Over the species' high-quality genomes (completeness ≥ 90%, contamination ≤ 5%), else its medium-quality ones (≥ 50%, ≤ 10%), else its genomes' assembly sizes as they are |
| `sized_genomes` | the genomes averaged |
| `rep_genome_size` | the representative's assembly size |
| `marker_bases` | the bases of the representative's marker genes, a copy of each (of every marker found, as `markers`: a [reduced database](#reduced-marker-sets) holds fewer) |
| `marker_share` | `marker_bases` over `rep_genome_size`: the share of the genome the markers cover |

-1 means unknown. `--build` packs the table into `database.protal`. Runs read the size columns by
their names, and a table without them (converted before 2026-10-08) leaves the sizes unknown: the
profiles then have no unknown share, and the run says so when it loads the database. Such a
database gets the sizes without being rebuilt: the converter's `--priors_only` writes the table
alone from the same release (the taxids are the same for a release, with or without
`--exclude_species` and `--genes`), and `--add_tables` stores it. protal checks that every taxid of
the table is a species of the database with the same representative genome, so a table of another
release is refused:

```bash
python3 scripts/mini_db/gtdb_to_protal_db.py --gtdb /data/gtdb_r226 --priors_only --outdir /tmp/priors -t 16
protal --add_tables /tmp/priors/species_priors.tsv --db /data/protal_r226_db
```

The first ten columns come out as before for a database `build_gtdb_database.py` converted (it
passes `--release`), so a model trained with `+priors` scores the same. A conversion run without
`--release` had left CheckM's and the clusters' values unknown (fixed on 2026-10-08: the release
found from the file names serves them too), and its table now gains them.

#### NCBI taxids

GTDB's species have no numeric ids, but its metadata gives every genome NCBI's species taxid
(`ncbi_species_taxid`). Since 2026-10-11 the converter writes `ncbi_taxids.tsv`, which gives each
species the NCBI species taxid most of its genomes have. A tie goes to the representative's taxid,
else to the smallest. `protal --ncbi_taxids` adds that taxid and its name to the profiles.

| Column | |
|---|---|
| `taxid`, `rep_genome`, `species` | the species in the database, its representative and its GTDB name |
| `ncbi_taxid`, `ncbi_name` | the majority NCBI species taxid, and NCBI's name of it (the species of `ncbi_taxonomy_unfiltered`) |
| `named` | 1 if NCBI names a species by it; 0 if it is one of NCBI's placeholders ("bacterium", "uncultured Lachnospiraceae bacterium", "Streptomyces sp. IB2014 016-6"), which GTDB's filtered NCBI lineage (`ncbi_taxonomy`) leaves out |
| `ncbi_genomes`, `genomes` | the species' genomes with that taxid, and all its genomes |
| `shared_by` | the species of the release with the same taxid (1: this one alone) |

At GTDB r226, 17.5% of the 143,614 species get a named NCBI species and 82.5% a placeholder.
Of the placeholders, 16,440 species have one of their own and 102,010 share one, up to 2,311
species on one taxid ("uncultured Clostridia bacterium"). NCBI files a placeholder under the
deepest rank it can name, so NCBI's tools place such a species there, not at species level. NCBI
also lumps species that GTDB splits: *E. coli* covers eight GTDB species. And it takes over
mislabelled genomes: 94,012 species have a single genome, and 58 GTDB species, most of them MAGs
of other genera, have "Akkermansia muciniphila" as their majority. `ncbi_genomes` and `shared_by`
show how far a taxid can be trusted ([report](claude/2026-10-11-ncbi-taxids.md)). A database
converted before 2026-10-11 gets the table as it gets the priors, without a rebuild:

```bash
python3 scripts/mini_db/gtdb_to_protal_db.py --gtdb /data/gtdb_r226 --priors_only --outdir /tmp/priors -t 16
protal --add_tables /tmp/priors/ncbi_taxids.tsv --db /data/protal_r226_db
```

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
the `.zst` file when the plain one is missing.

The build aligns with the WFA2 traceback `--wfa_traceback` names (`singletrack`, the default since protal 0.8.2, or
`classic`; [running.md](running.md#options-the-website-does-not-list)) and records it in
`build_settings.tsv`. Runs on the database align with the same unless their `--wfa_traceback` says
otherwise, and are warned when it does: the two tracebacks give the same scores but may choose
otherwise among equally good alignments, in the strain alleles, congener gaps and column weights
below and in the runs a model is trained on ([report](claude/2026-10-11-singletrack/README.md)).
`--sparse_first` with `--build` records each read type's sparse-first gate there too (`pe:1,pb:8`,
or one gate for all; `--sparse_first_score` and `--sparse_first_hits` with them), which runs on the
database follow unless their `--sparse_first` says otherwise: the database's models must have been
trained on runs gated so ([above](#choosing-the-sparse-first-gates)). `--sparse_first_pick` with
`--build` records each gate's pick likewise (`pe:hash`, or one pick for all; `stride` where it says
nothing, since 2026-10-13), which runs follow unless their `--sparse_first_pick` says otherwise (a
database that records its gates but no pick, built on 2026-10-12/13, picks by hash). The build
itself is the same. `protal --add_model FILE --read_type pe --model_gate 1` stores a model trained
at another gate than the database's as `model_pe_c1.xml`; a run gated at or nearest that gate is
profiled with it (`build_gtdb_database.py --model-gates`).

The build first makes the tables that need the genes alone (before the index takes its memory),
then the index:

1. Writes `gene_conservation.tsv`: how fast each gene diverges within species against the species'
   other genes (strains' copies compared with the representative's, k = 12; the median gene 1).
   Queries use it for the conservation and divergence features. It needs other genomes' copies.
2. Writes `gene_congeners.tsv` beside the database: how each gene differs between congeneric species
   against within species, to check that conserved genes are conserved between species too
   ([report](claude/2026-10-01-conservation-pattern/README.md)).
3. Looks for **suspect gene copies**: a species' copy within 0.02 (`--suspect_copy_distance`) of
   another genus's copy while its own congeners' copies are farther, a contaminating contig or a
   transferred gene. At r226 such copies drew a fifth of the false species calls. They go into the
   database as `suspect_copies.tsv`, and runs leave reads on them out of the evidence
   (`--keep_suspect_copies` keeps them). Every near pair across genera is listed in
   `gene_incongruence.tsv` ([report](claude/2026-10-03-false-positive-anatomy/README.md)). One gene
   at a time, its copies on all threads; a pair whose sketches cannot be within 0.05 is skipped
   unmerged (at r226 most of them), which changes no result.
4. Compares every two species of a genus by their marker genes (the median Mash distance of all the
   marker genes both references have, k = 12) and writes each species' nearest congeners, at most 16
   within 0.15, to `species_neighbours.tsv`: how crowded the database is around a reference, and how
   far apart two congeners a read fits are (since 2026-10-06; genera in batches of about 10,000
   species, so the sketches of the whole database are never held at once). A run's
   `relative_distance` takes the median over the genes with unique k-mers only, which the build does
   not know yet here: for near-identical congeners the table's distance is the smaller one. Training
   and profiling read the same table.
   Aligns every species' copy of each marker gene against the copies of its genus (WFA2, the ends partly free:
   all of them up to 24 others, else the 4 nearest by k-mer sketch and 16 drawn by a hash of the pair) and writes each
   copy's distance to its nearest congener's copy and the median distance to `congener_gaps.tsv` (since 2026-10-07):
   where in that gap a taxon's reads lie tells a strain from a congener the database lacks (the `gaps` features).
   Then takes the species' strain alleles from the full reference (since 2026-10-08): every other genome's copy of a
   gene (of `--allele_genome_share` of the genomes, by a hash of the accession the converter writes after each record's
   name; `build_gtdb_database.py` gives the training database 0.5 and the finished database 1, all genomes, since
   2026-10-11) aligned against the representative's, its differences kept as edits; per copy up to `--strain_alleles` (16;
   4 until 2026-10-11) of the distinct alleles, chosen among 16 sampled by hash for the sample's coverage (greedily:
   each sampled allele stands for the strains near it, and the chosen ones together explain as much of each as a read
   of it would be scored with; farthest first until 2026-10-09, which kept the deep lineages and dropped a strain's
   near relatives, [report](claude/2026-10-09-ancestry-true-positive-test/README.md)), each nearer the representative
   than the gene's nearest congener's copy (`congener_gaps.tsv`) and within 0.1 of it, into `strain_alleles.tsv`. A
   run scores each candidate alignment of a short read with its species' best allele (a read of a known strain then
   scores on its species as on the strain's gene; the records keep their alignment to the reference) and counts the
   alleles features. One pass over the full reference, as for the conservation factors. `--strain_alleles 0`: none.
   Then the column weights (since 2026-10-10, `column_weights.tsv`). Where the database folder has GTDB's
   alignments and tree (`column_msa/`, `species_tree.nwk`: `gtdb_to_protal_db.py` writes them from the release's
   `*_msa_marker_genes_reps` archives and trees, which `download_gtdb.py` fetches unless `--no_msa`), nothing is
   aligned (`--column_weights_alignment msa`, which `auto`, the default, takes when they are there): each family's
   columns are those of GTDB-Tk's alignment of the marker, each copy's residues are placed on its row (the row is
   the protein less its insertions), every genus of the family votes on the `among` conservation through its
   ancestral sequence (Fitch parsimony on GTDB's tree pruned to the genus), every species votes within its genus,
   and the family's ancestral base (parsimony over all its species, where three genera or more carry a base)
   polarises small genera's ancestry sites. On a real archaeal family this mapped 99.9% of GTDB's aligned residue
   pairs, against 56% by the protein alignment and 8% by the bases'. Without them, gene by gene, per family one
   genus reference (chosen by hash) gives the columns, the other genus references are aligned to it (up to `--column_weights`
   genera, 10, vote on the consensus among genera), every species' copy is aligned to its genus reference for the
   conservation within the genus. The alignments (WFA2) are of the proteins, the copies translated in frame, each
   aligned amino-acid pair mapping its codon's three bases (`--column_weights_alignment aa`, the default since the
   evening of 2026-10-10; up to 0.3 of the amino acids different within a genus, 0.5 among genera): GTDB's genera are
   too wide for the bases (`nt`: up to 0.2 within, 0.4 among), which lost two thirds of the copies of a real family
   ([report](claude/2026-10-10-real-ancestry/README.md)). The conservation within is averaged over the family's genera (since the evening
   of 2026-10-10 with the pseudocounts over all the species compared, so that the evidence of the genera adds up:
   before, no column could reach "conserved", 98.9%, and two of the `weights` features were 0 everywhere), and the
   amino acids per codon alike; every copy's mapping onto the columns is kept. One alignment per copy and one per genus
   reference and gene (3 minutes at r226 on 84 threads). `--column_weights 0`: none. The build log gives the failed
   alignments by kind and cause and how the columns' codes fall ("Column weights alignments", "Column weights codes").
   Then every copy's ancestry sites (since 2026-10-11, `ancestry_sites.bin`; `--no_ancestry_sites`: none): with the
   species neighbours, the congener gaps and the column weights just written, each copy compared with its congeners'
   copies as a run would on the copy's first record ([features.md](features.md#which-side-the-reads-take-ancestry-079)),
   gene by gene on all threads, and stored for the copies that have sites (~110 bytes each on the small real world;
   "Ancestry sites" in the build log). A run decodes them instead of computing them, which halved a sample's
   profiling instructions at the same outputs ([report](claude/2026-10-11-r226-v23/README.md), section 8).
5. Indexes `reference.fna` (its records taken from the genes it loaded, not read again), checks every
   k-mer's uniqueness against the full reference, and writes `unique_kmers.tsv`. With `--index_alleles 0.01` (on
   by default from 2026-10-09 until the r226 v22 build, off since: its F1 gain was too small for the entries it
   adds) the index also holds the k-mers of the deep strain alleles, those at least 0.01 from the representative's
   copy: a read of a strain that far from the representative seeds on few of its k-mers and could find its species
   only through a congener, so the allele's k-mers go in under the species as non-unique entries (about three per
   edit; the read is then aligned and scored against the representative as before). They are left out of
   `unique_kmers.tsv`'s totals, so a species' unique shares are the representative's alone.
   Since 2026-10-13 every entry also carries its grid level (0 to 5, in three bits of its position; index feature
   `grid levels`), for the sparse-first gate's grid pick and the sparse index (`--sparse_first_pick grid`,
   `--sparse_index`, [running.md](running.md#developer-options)): each gene copy's cores (the
   syncmers it is indexed by), in gene order, go through a nested tournament, each window of 2^l consecutive cores
   keeping one at level l or more, so that a run gated at C = 2^l counts about every C-th core of each gene and a
   sparse index of that gate holds about 1/C of the values. A window keeps the core a read of the species most likely
   carries intact: the share of the species' strain alleles (`strain_alleles.tsv`) that leave its 15-base core
   unedited, its 31-mer without indels and its flanks with two substitutions at most, pooled with a prior from how
   conserved its columns are within the genera of its family (`column_weights.tsv`; the prior weighs as two alleles,
   0.5 without weights); less a penalty for the entries sharing its key (log2, so that a core hitting many genes
   loses) and a small one for its distance from the window's middle (ties go to the middle). A core whose key the
   index leaves out never wins; a gene shorter than a window keeps its winner up to level 5. The build log gives the
   cores, the share the alleles cover and how many each level keeps with their mean conservation ("Grid levels"). A
   database whose genes are too long for the bits (a position of 2^17 or more) is built without grid levels, as are
   databases built before 2026-10-13; the stride and hash picks need none
   ([proposal and implementation](claude/2026-10-13-gene-grid-index.md)). `--grid_choice` (a developer option) builds
   the levels otherwise for comparison: `even`, every C-th core of each gene whatever the scores (a stride on the gene),
   or `phase`, every C-th core with the phase chosen by the scores; at C 4 on the mini database they lost 3.9% and 4.1%
   of the paired-end records against 3.2% for the default `scored` (all three before single-seed anchors, which brought
   `scored` to 0.05%).
   Since 2026-10-14, once the flags are final (after `unique_kmers.tsv`) and before the genes are released, every copy's
   cores are found again in its 2-bit sequence as a run finds them and each core's entry is looked up for its two flags
   (`unique_cores.bin`, "Unique cores" in the build log with the share unique; `--no_unique_cores`: none). The genes
   were released before `unique_kmers.tsv` until then; the uniqueness check before it held them beside the index as
   long, so the peak is the same.
6. Packs everything into `database.protal`, with `gene_table.bin`, `column_weights.bin` and
   `strain_alleles.bin` for a fast load: the index straight
   from memory, the other files compressed or copied; it reads each member back, compares, and
   removes the separate files. `full_reference.fna`, `gene2geneid.tsv` and `genome2tiid.tsv` stay
   beside it.

Every phase uses `-t` threads, and the index is the same for any `-t`. A full reference in zstd
frames with a seek table, as the converter writes it, is decompressed on several threads, a frame
each; the gene conservation factors are then the same on any run too (each species' copies of a
gene lie in one frame), and otherwise depend on which copies the threads reach first. A full
reference that cannot be read to its end (cut short, corrupt) stops the build. The log times each
phase. [Build options for the format](#build-options-for-the-format) lists the options for separate
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
  abundance × length, length from a gamma distribution, start uniform over the places where it fits;
  cut at a contig's end only if longer than every contig, so a genome of short contigs gets its share), by
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
| `--se_read_setups`, `--se_ultima_setup` | `--read_setups`', none | the single-end units' Illumina setups (a setup of either list simulated once); the design's Ultima reads, a unit a depth (since 2026-10-14) |
| `--species_per_sample`, `--archaea`, `--strains_per_species`, `--abundance`, `--congeners` | `5-30`, 0, one, sigma 1.3, 0 | the communities |
| `--read_types`, `--long_read_bases`, `--pb_setup`, `--ont_setup`, `--pbsim`, `--pbsim_models` | `pe` | other read types and how they are made |
| `--novel_species`, `--novel_clades`, `--taxonomy` | | the species the database lacks (`heldout_species.txt`), held-out clades per sample, and the taxonomy for the `meta_*` ranks |
| `--scenarios`, `--scenario_samples`, `--scenario_file`, `--host_genome` | none, 6 | scenarios to collect besides the design, their samples (each at a depth of its own), more or changed scenarios, the host genome |
| `-t`, `--jobs`, `--seed` | 4, `-t`, 1 | threads, cores for the simulations, seed |
| `--simulate_only`, `--prepare_profiling`, `--also_profile MAP` | | simulate and stop; write the map to profile and stop; profile another collection in the same run |
| `--read_compression` | zstd | `zstd` (`.fq.zst`) or `gzip` (`.fq.gz`): how the simulated reads are written |
| `--follow`, `--profile_block`, `--protal_lock`, `--min_free` | , 20, , 0 | profile as a `--simulate_only` run of the same `-o` simulates, in protal runs of at least this many GB of reads, removing each point's reads once profiled, then write the tables (points it cannot profile once the simulations have ended it simulates itself); a lock file for the protal runs of two followers to take turns; GB a simulation leaves free on the disk, or it waits (for a follower to remove reads) |
| `--stream_above` | 0 | `--follow` (and its `--simulate_only` run, the same value): points whose largest sample would take more than this many GB are streamed into protal through named pipes, not written; those whose communities are there share one protal run |
| `--profile_block_max`, `--profile_ahead` | 0, off | `--follow`: the most GB of written reads one protal run takes, at least one point's (0: no limit; streamed points are not capped since 2026-10-09); protal `--profile_ahead` in the collection's runs |

### Training

```bash
python3 scripts/machine_learning_cmdline.py --truth-file training/training_data.tsv --output-prefix training/model
```

The trainer needs Python 3 with numpy, pandas, joblib and scikit-learn, no Java: it writes the PMML
itself, so that protal's probabilities equal scikit-learn's bit for bit, and checks this on every
training row.

**The model.** A random forest (`--model forest`, the default before 2026-10-06 and again since 2026-10-14)
or gradient-boosted trees (`--model gbm`, the default from 2026-10-06 to 2026-10-14; scikit-learn's
`HistGradientBoostingClassifier`: 250 rounds at a learning rate of 0.1 of trees of up to 63 leaves,
20 rows a leaf at least, an L2 penalty of 1, balanced classes, no early stopping). Refitted on the r226 v13 tables, boosting scored higher than the
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
| `--model` | `forest` | `forest`: a random forest (the default again since 2026-10-14); `gbm`: gradient-boosted trees |
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
of highest F1 at the knob, but the default set (`normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood+ancestry+gaps+untried`) unless
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
those of the training data (`--read_type se`, `pb`, `ont` for other models). It takes the first and
last (`--points`) of the design points that have samples of that read type: with `--se_read_setups` a
point may hold paired-end or single-end samples only. It profiles into
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

The per-copy tables a run reads go in the same way, by their file names: `protal --add_tables foreign_rates.tsv --db DB`
(also `congener_gaps.tsv`, `strain_alleles.tsv`, `species_priors.tsv` and `ncbi_taxids.tsv`: [species priors](#species-priors), [NCBI taxids](#ncbi-taxids)). `scripts/foreign_rates.py --db DB --full-reference DB/full_reference.fna.zst --out
foreign_rates.tsv` makes the first: error-free reads every `--stride` bases (250) of at most `--per-header` (10) copies
of each species' gene in the full reference (`simulate_metagenomes --tiles --tile_fasta --tile_per_header`; the reads
named `<taxid_gene>:<record>:<start>`, the headers' records and tiles in `<tiles>.sources.tsv`), aligned against the
database once, each read's best record (MAPQ 4 or more) counted for its gene copy as a read of the copy's own species,
of another species or of another genus; every copy of the full reference gets a row, reached by a read or not, its
reads the foreign ones plus one genome's worth of its own (own reads over the records tiled), so that the species'
genome count stays out of the share. It
needs the database's `internal_taxonomy.dmp` (`--taxonomy`: a built `database.protal` packs it) and the full reference
the database was built from (a training database's lacks its held-out species, so the table knows nothing of them).
Until 2026-10-08 the scan tiled the genomes the training samples are drawn from, and a copy was listed only when a read
reached it, which told the models which species the simulation could draw
([features.md](features.md#the-foreign-features-leak)): the `foreign` features are not in the default set until a build
says what the new scan is worth; `--features auto` tries them.

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
