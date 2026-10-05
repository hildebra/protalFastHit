# Why trace_relatives.py traced 0 genes at r226 v10

Date: 2026-10-05. Code at `00bebbf` (audit-fixes). Data: `local/v10/model_logs/relatives_by_gene_conservation.{txt,tsv}`
of the r226 v10 build. The training samples and SAMs it read were on the build node's scratch disk and are gone,
so the v10 table cannot be redone; v11 will produce it.

## Symptom

The v10 report read 201 paired-end samples and summed the coverage of both classes:

- species in the database: 17,836.3
- species held out at species rank: 6,169.1

Yet it listed 0 genes ("Spearman ... over 0 genes: +nan"), and the `.tsv` held only its header. Because coverage
was found, the script did not say "Not traced", and the build reported the file as written.

## Cause

The script finds the genome a read came from by its name: `acc = QNAME.split("_contig")[0]`. That assumed the
read name starts with the genome's accession. It does not:

- `simulate_metagenomes` runs ART on each genome FASTA without renaming anything.
- ART names a read after its contig: `<contig>-<n>`.
- The synthetic genomes of `simulate_gtdb_release.py` name their contigs `<accession>_contigN`, so in every mini
  world and in `GtdbBuildTest` the split gave the accession and the trace worked.
- GTDB's genomes carry NCBI's contig names (`NZ_CP012345.1`, `JAAAAA010000001.1`). The split returned the whole
  contig name, no record matched a source genome, and no gene got a record.

The gene factors were found: `gene_congeners.tsv` stays a loose file in the training database folder.

Reproduced with the old script on a small world (`scripts/test_trace_relatives.py`):

| contig names | old script | fixed script |
|---|---|---|
| NCBI-style | 0 genes | 3 genes |
| `<accession>_contigN` (mini worlds) | 3 genes | 3 genes |

## Fix

**`trace_relatives.py`**

- Maps contigs to genomes from the genomes' own FASTAs, which the samples' `manifest.tsv` lists (`fasta_path`).
  Only the headers' first words are read, once per genome, on `--threads` processes. Each sample's records are
  traced through that sample's map.
- A contig name that two species of a sample share traces no read; the report counts such records. One shared
  within a species is kept.
- Every genome of a sample is in the map, those of species held out at other ranks too, so their reads are skipped
  rather than counted as untraced.
- Old manifests without `fasta_path` keep the `_contig` rule.
- The report says how many records matched no contig. When no record is traced at all, the script now says
  "Not traced: ... none of their aligned records was traced", so the build's console line shows the failure.

**`build_gtdb_database.py`** passes its `--threads` to the trace.

**`insilico_strains.py`** names an in-silico strain's contigs `<strain>_<representative's contig>`. Its reads were
otherwise named like its representative's: harmless for this trace, which needs only the species, but wrong for
anything that maps reads to genomes by name.

**Tests**

- `scripts/test_trace_relatives.py` (standard library only, run in CI and by `just mini-db-test`): NCBI and
  synthetic contig names, an old manifest without FASTA paths, the "nothing traced" message, read names with
  dashes.
- `test_insilico_strains.py` checks the renamed contigs.

## Cost

One header pass over each genome simulated in the paired-end training samples. That is up to ~17,000 genomes at
r226, ~75 GB gzipped, read on the build's threads: a few minutes on the 64-thread node. The SAM pass, at most 5M
records per sample, is unchanged.
