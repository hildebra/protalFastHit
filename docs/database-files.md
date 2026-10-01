# Database files

A protal database consists of six parts, optional tables of the genes' conservation and of their
neighbours, and models of other kinds of reads:

| File | What it holds |
|---|---|
| `index.prx` | the k-mer index of the reference marker genes |
| `reference.fna` | the reference marker genes, one record per gene and species (`>taxid_geneid`) |
| `reference.map` | taxid, gene id, start and end byte of each gene in `reference.fna` |
| `internal_taxonomy.dmp` | the taxonomy: id, parent, external id, name, rank, level, representative genome |
| `unique_kmers.tsv` | per species and gene, counts of the k-mers unique to it in the database; the uniqueness features of the model |
| `gene_conservation.tsv` | optional: per gene id, how fast the gene diverges within species compared with the species' other genes (1 for a typical gene, below 1 for conserved ones), which `--build` estimates from `full_reference.fna` ([building-a-database.md](building-a-database.md#2-build-the-index)); every query reads it for the model's conservation features (`conserved_fast_depth_ratio`, `conserved_hit_share`), and `--gene_conservation db` also scales `--depth_identity_margin` per gene by it, which is not the default. A database of an earlier version, or built from one genome per species, has none (the features are then 0 and 0.5) |
| `gene_neighbours.tsv` | optional: per clade (family up to domain), which end of which gene faces which end of another gene within 3 kb in the representative genomes of its species, how often and how far apart; from whole genomes by `scripts/mini_db/gene_neighbours.py` ([building-a-database.md](building-a-database.md#gene-neighbours)). protal uses it to pair mates and follow long reads over neighbouring genes ([running.md](running.md#options-the-website-does-not-list)). A database of an earlier version, or built without genomes, has none |
| `model_pe.xml` | the presence model (a random forest in PMML) of paired-end reads, see [model-training.md](model-training.md); `model.xml` in databases of earlier versions |
| `model_se.xml` | optional: the presence model of single-end reads; without it, single-end samples need `--model_se` ([running.md](running.md#single-end-reads)) |
| `model_PB.xml`, `model_ONT.xml` | optional: the presence models of PacBio and Nanopore reads (`--model_pb`, `--model_ont` without them) |

`protal --build` packs them into one file, `database.protal`, compressed with
[zstd](https://facebook.github.io/zstd/): one file to copy, download, checksum or version, whose
parts cannot get out of step. `--db` (or `$PROTAL_DB_PATH`) takes that file, a folder that holds
it, or a folder with the separate files (raw, or `index.prx.zst` / `reference.fna.zst`). A folder
with both uses the separate files.

Databases from earlier protal versions (a folder of raw files, index format 1) load unchanged.
The index records which features it was built with; protal prints them when it loads it
(`Index features: ...`) and checks the index against `reference.map` and `reference.fna`.

In `internal_taxonomy.dmp`, ids and names are unique (truth files and `--msa_species` name
taxa), every parent is defined, and every lineage ends at the root, a taxon that is its own
parent. `reference.map` and `unique_kmers.tsv` list at least one gene. `gene_conservation.tsv` has a
line `geneid<TAB>factor<TAB>species` per gene (the header line and the species column may be left
out), factors above 0. `gene_neighbours.tsv` has ten numbers per line (`clade gene end partner
partner_end species informative gap_median gap_min gap_max`; ends 5 or 3, partner 0 for no gene
within 3 kb), every clade in the taxonomy and every gene in `reference.map`; `--build` checks it
before packing it. protal stops at the first problem, with the file and line.

## Compression

Compression saves disk space and makes loading faster wherever storage is slower than
decompression, e.g. on network file systems. The mini database takes 1.0 MB as `database.protal`
and 3.2 GB raw (the index's k-mer key map has a fixed size); in `database.protal` the text files
and the model are compressed too (`model.xml` about 20x).

protal writes zstd's *seekable* format: independent frames of 64 MB each (plus a seek table at the
end), so that loading uses `-t` threads, each decompressing whole frames straight into memory
(~1-1.5 GB/s per thread). In `database.protal` each part is a range of frames (after a directory
frame), read exactly as the separate `.zst` file would be, so a single file loads as fast as
separate files. Raw files are read with `-t` threads too. A `.zst` with a single frame (e.g. from
the zstd CLI) still loads, but with one thread.

`index.prx.zst` holds the index in columns rather than as the raw bytes of `index.prx`: each
frame is a chunk of about 64 MB of the index, stored as the number of values per k-mer, then each
field of the values (taxon, gene, position, flags) in its own byte planes, only as wide as the
chunk needs; empty key blocks cost one bit. In tests this was 25-55% smaller than the raw bytes
compressed alike, and loaded about as fast (within ~10% on a synthetic index with 680 MB of values,
faster on sparse ones), since chunks decode in parallel. The writer decodes every chunk it writes
and keeps the raw cells of any that would not round-trip exactly. Because of the columns, plain
`zstd -d` does not give an `index.prx`; use `protal --unpack_db` or `--decompress_db` (below).
`reference.fna.zst` is plain seekable zstd, and so is `database.protal`, but `zstd -d` on it gives
all parts one after another.

## Build options for the format

`protal --build` ([building-a-database.md](building-a-database.md)) writes `database.protal` into
the `--db` folder, reads it back and compares it, and then removes the separate files it holds
(`--unpack_db` writes them back). Frames are compressed in parallel with `-t` threads.
`full_reference.fna` and the other build inputs stay as they are.

| Build option | Default | |
|---|---|---|
| `--no_bundle` | off | keep separate compressed files (`index.prx.zst`, `reference.fna.zst`, ...) instead of `database.protal` |
| `--no_compress` | off | separate raw files (`index.prx`, `reference.fna`, ...) |
| `--compress_level` | 19 | zstd level 1-22. Level 19 compresses ~3 MB/s per thread, level 12 ~40 MB/s; decompression speed barely depends on it |
| `--compress_frame_mb` | 64 | frame size; 0 writes a single frame (one loading thread; needs `--no_bundle`) |
| `--compress_window_log` | 27 | long-distance matching window (2^27 = 128 MB, capped at the frame size); 0 turns it off |

## Converting a database

An existing database is converted in place, without rebuilding it; the content stays
byte-identical and everything is verified before the old files are removed:

```bash
protal --compress_db --db /path/to/protal-db -t 16      # separate files -> database.protal (--no_bundle: index.prx.zst, reference.fna.zst, ...)
protal --unpack_db --db /path/to/protal-db -t 16        # database.protal -> separate files next to it (it is kept)
protal --decompress_db --db /path/to/protal-db -t 16    # -> separate raw files, e.g. for protal versions before zstd support
```

`--unpack_db` writes `index.prx.zst` (its frames as they are in `database.protal`), an
uncompressed `reference.fna`, and the other files, into the folder `database.protal` is in, or into
`--unpack_dir`. `--decompress_db` writes all files raw and removes the compressed ones
(`index.prx.zst`, `reference.fna.zst` or `database.protal`); the index is byte-identical to one
built with `--no_compress`, and `--compress_db` on such a folder writes `database.protal`
byte-identical to the one `--build` wrote. `--compress_db` needs the index in memory.

`protal --add_model MODEL --read_type pe --db db/` (or `se`, `pb`, `ont`) checks a model and
replaces the database's model of that read type, as `scripts/build_gtdb_database.py` does. By hand:
unpack, replace the file and pack again; the separate files win over the old `database.protal`,
which `--compress_db` then replaces:

```bash
protal --unpack_db --db db/ -t 8
cp new_model.xml db/model_pe.xml
protal --compress_db --db db/ -t 8
```

A model of single-end reads is added the same way, as `db/model_se.xml`: `--build` and
`--compress_db` pack a `model_se.xml` they find in the folder. (`--model new_model.xml` and
`--model_se` use other models without changing the database at all.)

## Gene order in reference.fna

The order of the genes in `reference.fna` matters for its size: related species' copies of a
marker gene are similar, so a file ordered by gene (then by taxon, in taxonomic order) compresses
about 2x better than one ordered genome by genome, if related genomes lie further apart in the
file than zstd's window. protal finds genes through `reference.map`, so any order works;
`scripts/mini_db/gtdb_to_protal_db.py` writes gene order by default.

## Genes in memory

protal holds the reference genes in memory two bits per base, four bases to a byte (a quarter of
the byte per base that `reference.fna` uses), and decodes a gene where the alignment code reads it;
the decoding uses AVX2 where the CPU has it. Only A, C, G and T have codes of their own. Any other
character is stored as a base, so a gene that protal reads back differs from `reference.fna` at its
ambiguous bases: `N` is stored as `A`, and an IUPAC ambiguity code as the first base it stands for
in the order A C G T (`R`, `W`, `M`, `D`, `H`, `V` as `A`; `Y`, `S`, `B` as `C`; `K` as `G`).
Lowercase letters are read as uppercase, and any other character is stored as `A`. The files
(`reference.fna`, the index) are unchanged. What the coding changes:

- Alignments are made against the stored bases, so a read is scored against `A` or `C` where the
  reference has `N`, and the SNP and MSA reference rows show the stored base. A SAM file that an
  earlier protal made against a gene with an ambiguous base can hold an `M` there that no longer
  matches the stored base; the profiler sets such records aside (`<sam>.err`) with a warning.
- `--build` leaves every k-mer whose window holds an ambiguous base out of the index (and out of the
  uniqueness check): it is not counted, placed or looked up, so the index holds k-mers of A, C, G
  and T only, and the k-mer that `--build` reads back from a gene for the uniqueness check is always
  the one that was indexed. Before, such a k-mer was indexed with the ambiguous base read as `A`; a
  database built that way still works, it holds a few more k-mers. Reads are unchanged: a k-mer of
  a read with an `N` is looked up with the `N` as `A`.

## Loading genes on demand

`--preload_genomes_off` (loading reference genes on demand) reads single genes from a raw
`reference.fna`, so it needs the database as separate files. On `database.protal` protal stops and
prints the `--unpack_db` command for it and the command that reruns protal on the unpacked files.

## Measuring

`scripts/db_compression_benchmark.sh` measures ratio and speed per level on your database and
compares protal's load times for raw and compressed copies and several thread counts.
