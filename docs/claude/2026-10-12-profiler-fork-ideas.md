# Changes made to another k-mer profiler: which ones protal could use, and two-stage profiling

The question: a fork of another metagenomic profiler (whole-genome k-mer sketches) lists a set of
changes, mostly written by Claude. Which of them would also work for protal? A follow-up at the end
takes up the user's variant of two-stage profiling: candidate reads kept in a temporary file.

## Based on

- The fork's README (its list of changes) and `ARCHITECTURE.md`, read on 2026-10-12. I did not read
  its code or run it. Its numbers below are as its README states them; it gives no benchmark tables
  or data sets for them.
- protal at `714649b` (audit-fixes), read only. No protal code was changed for this report.
- The 16th cluster run (`local/Perf0.8.2`: protal 0.8.2 `4a42f3b`, r226 v24 database 30 GB on NFS,
  EPYC 9634, 8 threads; the real pe sample PEARL_AGE B10 with 50,607,792 pairs and a HiFi barcode).
  Log lines are quoted from `logs/pe1.2.log`.
- One check in WSL: `zstd -t` on a test `database.protal` (5.9 MB, `~/protal-composition/val/db_all`)
  and on a copy with one byte changed in the middle. The first passes; the copy fails with
  "Decoding error (36) : Data corruption detected".

## Summary

| Change in the fork | protal today | Worth porting? |
|---|---|---|
| Two-stage profiling: a sparse screen, then dense data only for the genomes it passes (the fork reports 7–27× faster and ~6× less RAM) | One index of all taxa, loaded whole: 26.6 GB values + 3.2 GB key map of the 36.3 GB after loading | **Yes, the only large one, but only as an approximation, and measured first** (below) |
| Growable databases: append genomes, rebuild only the screen, same result as a full conversion | Index flags, tables (column weights, neighbours, congener gaps, ancestry sites, alleles) and models depend on all genomes. `--add_model` already updates in place with a journal (`fd85db9`) | No. An appended genome changes other taxa's flags and tables, and the models were trained without it |
| Compressed sketches: sorted hash deltas, Golomb–Rice, zstd (−55% samples, −30% databases) | No sketches. The on-disk index is already in columns (byte planes per field) in seekable zstd frames, 25–55% smaller than plain zstd ([databases.md](../databases.md)) | Maybe, as a measurement: delta-coding the sorted taxon-gene codes within a key. Only cold loads and download size would gain (below) |
| Reference-delta and error-k-mer compression of samples | A sample's kept data is its `.sam.zst` (1.70 GB of SAM text for 50.6M pairs at r226) | Low. CRAM-style SAMs would be smaller, but the SAM is an intermediate file and small next to the FASTQ |
| Several read inputs or sketches merged into one sample | One read file (or pair) per sample. `protal_map_utils generate` makes `x_S1_L001` and `x_S1_L002` two samples | **Small, user-facing**: several files per sample (lanes) read as one stream, so depth-based features see the whole sample |
| Multi-threaded read sketching, identical output at any `-t` | Done: reader thread, BGZF on ISA-L workers (`822c4cd`), order-independent sums (`3c2d56c`) | Nothing to port |
| Adding the unknown estimate to an old profile | The `?` row and `.composition` come with every run (`8ccc31c`), and `--profile_only` re-profiles from the SAM (profiling the r226 pe sample takes 12 s) | Nothing to port |
| Empty read inputs tolerated | An empty read file gives a warning and the run goes on (e2e `test_empty_read_files_are_warned_about`) | Nothing to port |
| AVX-512 k-mer extraction | Done: syncmers and flex scan, `PROTAL_SIMD` (`f53db5e`) | Nothing to port |
| Seeded RNG, genome order kept, byte-identical at any `-t` | No random numbers in the run path (`grep` of `src/` outside the simulators finds none), order independence `3c2d56c` | Nothing to port |
| Paired read length: mean of both mates, mates shorter than k left out (for the unknown estimate) | protal counts the bases it scanned, overlapping mate bases once (`Composition.h`), and estimates no read length | Not applicable |
| Corruption checks: whole-file XXH64, zstd frame checks | Every zstd frame carries a content checksum (`Zstd.h:360`), checked on every load. `zstd -t database.protal` checks the whole file and finds one changed byte (checked above) | Nothing to port. A line in [databases.md](../databases.md) about `zstd -t` would help users check downloads |

## Two-stage profiling for protal

### Where protal's memory and time go

From the r226 pe run (8 threads, 115 s, peak 37.5 GB):

- "Index in memory: 2903557366 entries of 41 bits … 26.6 GB; key map 3.2 GB": 29.8 of the 36.3 GB
  after loading. Genes are 4.12 GB packed, and column weights 1.47 GB.
- Aligning is 92.6 s, and two-thirds of it is seeding: 2.24G k-mer lookups and 156G flex cells
  scanned (70 per block). Each block holds every taxon that shares the k-mer.
- Which taxa the reads reach:
  - The failed candidates of unaligned reads name all **143,614** species.
  - The read EM's ambiguity classes name **11,366** species (7.9%).
  - The SAM header lists **78,064** of 14,527,738 genes (0.54%).
  - The profile calls **590** species.

So almost all of protal's memory is spent on taxa that no read aligns to. Most of the seeding time
goes to scanning their entries in shared blocks. The win could be as large as the fork reports.

### Why it cannot be exact

A whole-genome sketch profiler decides per genome, from thousands of k-mers per genome. protal
decides per read, among competing candidates: the top 3 plus adaptive congeners. MAPQ and the
ambiguity classes come from those competitors. If a stage-2 index lacks a read's best competitor:

- the read's MAPQ rises;
- a cross-genus or relative read can become a confident record.

That is where protal's false positives already come from. At r226 v15, false-positive taxa were fed
by reads of held-out congeners at about 0.97 identity, 80–92% of them from the same genus
([error-read signatures](2026-10-07-error-read-signatures/README.md)). In the v5 model's paired-end
false positives, 52% were single reads at identity ≥ 0.97
([false-positive anatomy](2026-10-03-false-positive-anatomy/README.md)).

Exact equality would need every taxon that appears in some read's tried candidates. At r226 that is
all 143,614 species. A two-stage protal therefore changes outputs by design. It is only acceptable
if the changes stay within the build evaluation's noise.

### A screen that keeps the faint taxa

protal's misses are already thin and divergent strains. At r226 v15, false-negative taxa were
divergent strains of wide species ([error-read signatures](2026-10-07-error-read-signatures/README.md)).
They have fewer exact k-mers, so a screen loses their reads first. A sparse species-level index can
keep most of them:

- Keep 1/c of the keys, chosen by hash. The index holds every syncmer of every gene, so a read from
  a gene in the database finds all of its syncmers there, up to sequencing errors.
- The pe sample made about 44 lookups per pair (2.24G / 50.6M). Assume a marker pair makes about as
  many, so it expects 44/c sparse hits on its own taxon.
- If those hits follow a Poisson distribution, a pair has no sparse hit with probability:
  - c = 8: 0.4%;
  - c = 16: 6%;
  - c = 32: 25%.
- Divergent strains have fewer exact k-mers, so they would lose more.
- With c = 8 the screen holds about 26.6 / 8 = 3.3 GB of values, plus a key map sized to those keys
  (0.4 GB instead of 3.2 GB).

### What stage 2 would load

Stage 2 would load:

- the species the screen passes;
- their competitors. `species_neighbours.tsv` already lists 577,071 congeners within 0.15, the set
  `--adaptive_candidates` uses.

The size of that set is not measured. If it is 15–25% of species, stage 2 holds about 4–7 GB of
index and 0.6–1 GB of genes. With the screen, a run would need about 9–12 GB instead of 37 GB.
Seeding would scan correspondingly fewer flex cells per block, because blocks shrink with the taxa
left out.

### Costs

**A second pass over the reads.** Most mates hit the index somewhere: 79.5M of 101.2M mates have an
anchor. So keeping only "reads with a hit" saves little.

The pe sample is single-member gzip. That is one inflater at about 1.8M pairs/s, so a second read
costs about 28 s for this sample. BGZF or zstd input would cost less. The follow-up below avoids
most of it.

**A new database layout.** protal's index is keyed by k-mer, with all taxa in each block. Stage 2
needs one of these:

1. **Filter while loading.** Keep the entries of the passed taxa while the index decodes. This is
   simple, and the existing `ValueSink` packs entries key by key already. RAM drops, but the whole
   index member is still decoded (about 18 GB compressed in earlier builds): 9 s warm, 107 s cold
   over NFS in this run.
2. **Rebuild from the passed taxa's genes.** Extract their syncmers again from the 2-bit genes. This
   reads only those genes, but the build's per-entry flags (unique-distance, checked uniques) must
   be stored per gene, because they are global.
3. **Seekable per-clade index blocks.** This is the fork's layout, and a new database format.

**More work on every sample.** Thresholds for the screen have to be chosen, evaluated per read type,
and the models may need retraining on two-stage runs.

### Suggested first step: an oracle experiment (about 1–2 days)

1. Add a development option `--restrict_taxa FILE`. It drops other taxa's index entries while
   loading (option 1 above) and leaves their genes out of the preload.
2. Run the r226 pe and HiFi samples, and the build's held-out test samples, with these lists:
   1. the species a full run's ambiguity classes name;
   2. the same plus their congeners within 0.15;
   3. all species of those genera.
3. Compare the following against the full run:
   - calls and depths;
   - F1 from the build's evaluation;
   - SAM records;
   - RSS;
   - seeding time.

If even list 2, chosen with hindsight, moves F1 beyond noise, a two-stage protal is not worth
building. If it does not, stage 1 is next: a sparse screen with c around 8. A taxa list from a
previous run would also allow a "known community" mode for cohorts.

## The on-disk index (compressed-sketch idea)

The fork's gain comes from sorted uniform hashes, whose gaps follow a geometric distribution, so
Golomb–Rice codes them well. protal's equivalent is the list of entries within a key:

- `taxid × 169 + gene`, 25 bits;
- position, 14 bits;
- 2 flags.

Sorted by taxon-gene code, the deltas would be small where related species share a k-mer. 71% of
the entries are in keys with 33 or more values ([memory audit](2026-10-03-memory-audit/README.md)).

Earlier work on this:

- Delta coding was rejected for memory (−9 GB, but lookups 2–3× slower).
- On disk, entries are decoded once at load, where cost matters less.

What it would gain:

- **Cold loads.** In this run the first load took 107 s, and later ones 9 s. The cold load read
  about 30 GB in 107 s, about 280 MB/s, so every 3 GB less is about 10 s less on a cold node.
- **The database download.**

To measure first: on a world index, the per-key delta + Rice size against the current columns + zstd
(`scripts/db_compression_benchmark.sh` gives the latter).

## Several read files as one sample

The fork merges sketches by summing their counts. For protal, the useful form is reading a sample's
lane files one after another as one stream: `L001`..`L004`, the R1 and R2 of each lane paired.

- It must happen at the read level. The sample-depth features (`sample_log_fragments`, the knob by
  depth) and the read EM need the whole sample. Profiles of lanes cannot simply be summed.
- Today, users have to concatenate the files themselves. `protal_map_utils generate` makes each
  lane its own sample.

This needs:

- a list per map cell, or one row per lane with the same `#SAMPLEID`;
- `ThreadedGzStream` over a list of files.

It is small. The website's map-file page would need the new form.

**Implemented 2026-10-12** (working tree on `3931e48`, not committed): a map's `FIRST` and `SECOND` cells list a
sample's files comma-separated, and `--one_sample` makes the files of `-1`/`-2` one sample. They are read one after the
other as one stream, and a lane whose mate files differ in number fails the sample. `protal_map_utils generate` makes
Illumina's lanes one sample. ctest 522/522 and e2e 160/160 pass. See
[running.md](../running.md#several-files-of-one-sample).

## What was not done

- No code changed for this report, nothing committed.
- The fork was not run, and its code not read.
- The stage-2 set size (15–25% of species) and the memory estimate built on it are assumptions,
  until the oracle experiment measures them.

## Follow-up: stage 1 keeps its candidate reads in a temporary file

The user's question: instead of reading the input twice, could stage 1 write the reads it cannot rule
out to temporary disk (ISA-L gzip), and stage 2 read only those?

Yes. This is better than a second pass, provided stage 1 can rule out most reads. The numbers are
from the same r226 pe run as above.

### What it saves

The second pass over this input would cost about 28 s:

- 50.6M pairs of single-member gzip;
- one ISA-L inflater, about 1.8M pairs/s.

The pairs that end in the profile are only 1,400,393 fragments, 2.8%. If stage 1 keeps a few
percent of the pairs, stage 2 reads back seconds' worth instead.

The temporary file:

- **Format.** BGZF written with protal's ISA-L `Deflater` (level 1) by the alignment threads, batch by
  batch. Stage 2 reads it back with the parallel BGZF reader (1–4 workers per file, 2.7–3.4 GB/s).
- **Location.** Node-local disk (`$TMPDIR`, or a `--tmp_dir`), not the NFS input folder. The page
  cache keeps it in memory while memory is free.
- **Size.** At 5% of this sample, 2.5M pairs are about 1.8 GB as FASTQ and about 0.5 GB as BGZF.

### What decides it: how selective stage 1 can be

protal looks reads up by a 15-base syncmer core. The 16 bases around it (the flex cell) are scored by
the bases they share with the read. Cores alone are common: 79.5M of 101.2M mates found an anchor,
and 90.3% of the 176M candidate alignments were refused by the k-mer screen. So "a hit in the stage-1
index" would keep most reads.

A useful rule needs two parts:

- **Sparsity.** Keep only cores whose hash falls in 1/c. This is decided without a memory access, so
  stage 1 makes 1/c of the lookups, into an index about 1/c of the size.
- **Near-exact flanks.** Count a hit only if its flex score reaches a threshold T. Random flanks then
  rarely pass.

A pair makes about 44 lookups, so it expects 44/c sparse hits on its own taxon. Sequencing errors and
divergence lower that. Two numbers have to be traded off:

- **F:** the share of pairs written to the file. It decides the file size, and how much work stage 2
  does.
- **Recall:** the share of pairs that align in a full run and also pass stage 1. Every pair lost is
  lost depth, and taxa with one or two pairs may disappear. 60–76% of v24b's false negatives had two
  or fewer reads of their own ([v24b errors](2026-10-12-r226-v24b-errors/)).

If F is small, the gain is time as well as memory:

- Seeding is about two-thirds of aligning.
- Stage 1 would make 1/c of the lookups, into an index 1/c the size.
- Stage 2 would seed only F of the pairs, into a small index.

### Changes in what the models see

Pairs dropped by stage 1 never become failed candidates. Stage 2 also lacks the taxa it did not load,
so these change:

- the failed-candidate counts per taxon;
- the MAPQ competitors;
- the read EM classes.

The models would have to be trained on two-stage runs, which means a build whose training runs use
it.

### Measuring F and recall before building anything

A development option can measure both in an ordinary full run:

- Each fragment is also judged by the stage-1 rule, for several c (for example 2, 4, 8, 16) and T
  (for example 14–16 of 16).
- The log reports, per (c, T): F, recall on the fragments with SAM records, and the called taxa
  whose fragments mostly fail.
- No new database is needed: the sparse index is the full index with the other cores ignored.

This takes about half a day to write and one cluster run of the r226 pe and HiFi samples. It
complements the `--restrict_taxa` experiment above: one measures which reads stage 1 must keep, the
other which taxa stage 2 must load.

Implemented as `--stage1_probe`: see [2026-10-12-stage1-probe](2026-10-12-stage1-probe/README.md).

**Measured at r226** ([results](2026-10-12-stage1-probe/README.md#results-at-r226),
[lookups](2026-10-12-stage1-probe-r226.md)): with every lookup, T = 14 keeps 9.9% of the pairs and 99.99% of those
with records, so the temporary file is small. A sparse index costs recall faster than estimated above: c = 2 keeps
99.63%, c = 4 97.87%, c = 8 92.90% (the estimate said 99.6% at c = 8). For paired-end reads, the memory saving of a
sparse first stage is therefore small (c = 2: half the index); HiFi reads keep 99.17% at c = 8.
