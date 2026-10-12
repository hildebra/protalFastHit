# Performance audit of protal 0.8.1, 2026-10-11: one table costs 12 s per run, and what else 0.8.x costs

- **Date**: 2026-10-11.
- **Code**: 0.8.1 = `b3ac271` (the cluster run: `7d5df4b`, which adds docs only), against 0.7.8 = `5de3224` (the 13th
  cluster run's `230bf64` differs from it by `--add_model` and the version only). Built from `git archive` in WSL
  (Release, the one binary; [build_versions.sh](scripts/build_versions.sh); 0.8.0 = `a82df8f` was built too but not
  needed). One experiment build, `tim` = 0.8.1 with a timer and a peak-memory reset around each start-up table and an
  env-gated release of a sample's evidence ([audit_timers.patch](scripts/audit_timers.patch), [build_tim.sh](scripts/build_tim.sh)).
  Nothing in the tree was changed.
- **GTDB scale**: the user's 14th cluster run (SLURM 24161547, `scripts/measure_performance.sh` at `7d5df4b`, node
  `q512n10`, EPYC 9634, 32 threads, the r226 v23 finished database, the real gut pe sample of 50.6M pairs and the HiFi
  barcode of 497,656 reads; copied to [results/cluster_v14](results/cluster_v14/)) against the 13th (`230bf64`, the same
  node, samples and threads, the r226 v11 database; `local/Perf0.7.5_v13`). Besides, the r226 v19-v23 builds' own protal
  runs (their share archives in `local/`; [r226_log_stats.sh](scripts/r226_log_stats.sh)).
- **Local machine**: WSL2 Ubuntu 24.04 on the Core Ultra 7 258V, GCC 13.3, every job pinned to cores 0-3 (`taskset -c 0-3
  nice -n 5`). Another session compiled and tested in WSL meanwhile (load 4-12) and the cores change speed under
  Windows: wall times are medians of three alternated rounds, and the exact comparisons are callgrind's instruction
  counts (one thread).
- **Local world** ("arc"): every GTDB r226 archaeal species (6,968; the marker genes of all 17,245 genomes), built by
  0.8.1's own `build_gtdb_database.py` ([arc_build.sh](scripts/arc_build.sh)), the finished database with the pipeline's
  own `--build` arguments ([build_arcdb.sh](scripts/build_arcdb.sh): alleles from every genome, 16 per species and gene,
  none indexed): 303,776 gene copies, 425 MB, 1/48 of r226's copies and 1/28 of its column-weight columns. The pipeline
  was stopped once its databases were built (its training was not needed); the small real world's 0.8.x models
  (`~/realworld_s1`, pe and PacBio) were stored for all four read types. 0.7.8 refuses an index written by 0.8.1
  ("Unsupported index.prx version 0.8.1"), so it got its own build of the same genes ([db078.sh](scripts/db078.sh)).
  Samples from the 3,259 whole genomes ([sim_pe.sh](scripts/sim_pe.sh), [sim_long.sh](scripts/sim_long.sh)): 1M pairs
  of 150 bp from 120 species with strains and congeners (its R1 the single-end sample), HiFi and ONT reads of the same
  community (100 Mb each; 5 Mb for callgrind), its first 100,000 pairs for callgrind, and a cohort of 24 samples of
  50,000 pairs.
- **Earlier rounds**: [GTDB scale](../2026-10-04-performance-gtdb-scale/README.md), [profiling 2026-10-06](../2026-10-06-performance-profiling/README.md)
  (the 12th and 13th runs), [build parallelism](../2026-10-09-build-parallelism/README.md), [r226 v23](../2026-10-11-r226-v23/README.md)
  (section 5: the feature groups' profiling costs; section 8: the prebuilt ancestry sites).
- **Question**: audit protal's performance, with the latest version (0.8.1).
- **Follow-up** (the same day): items 1-5 implemented, on `2caffd7` and checked again on `0dea1ff`; tests, identity
  and measurements in [section 8](#8-follow-up-items-1-5-implemented). Uncommitted when written.

## Summary

1. **0.8.1 is slower than 0.7.8 at r226**: the real pe sample 33.4 → 51.0 s (+53%), the HiFi barcode 12.3 → 24.8 s
   (+102%), peak memory 34.2 → 44.6 GB (+10.4 GB). For pe that is +12.5 s of start-up, +7.4 s of profiling and −1.6 s
   of aligning; for HiFi +12.6 s of start-up and +0.8 s of profiling.
2. **The start-up's +12.5 s is one table**, `column_weights.tsv` (0.8.1's column weights: 602M columns at r226). Every
   run parses its ~2.2 GB of text on one thread before the alignment may start, at 140 instructions a column (exact,
   locally). Yet only profiling reads it, and only for the copies a sample's records touch: ~900 of 13.9M for a median
   sample of the build's design. The r226 builds' runs date the jump to the v22 build, the first with the table.
3. **The memory**: resident after loading +5.2 GB, peak +10.8 GB. The column weights hold ~2.5 GB, and loading them
   needs ~7.6 GB more for a moment: text, a full `istringstream` copy, rows, sorted copy. A profile-only run with no index
   shows this alone, its peak 7.7 → 15.9 GB. The rest is the strain alleles (~0.9 GB), the congener gaps (~0.5 GB) and
   the 111M index-allele seeds of the v23 build (~0.7 GB).
4. **The profiling's +7.4 s (pe) is the ancestry sites.** The v23 database predates `ancestry_sites.bin`, so 0.8.1
   builds every touched copy's sites while it reads the SAM (1.8 → 6.9 s). A database built by 0.8.1 stores them
   (−47% of a pe sample's profiling instructions, the v23 report). ~1.9 s more lie outside the sample's timer, not
   attributed here.
5. **The alignment stage held, and ONT got much faster.** On the cluster, aligning pe took −6% (seeding −9%, taking
   the k-mers −24%: the AVX-512 kernels) and HiFi the same. Locally, exact: pe −2.0%, PacBio −0.6%, ONT ≥ 10.5× fewer
   instructions (0.8.0's indel bound; 54.9 → 6.2 s). The new per-candidate work costs 0.14% of the pe aligning (allele
   scores, adaptive candidates, sites' settling), the gene neighbours 0.4%.
6. **A 0.8.1-built r226 database (v24) will weigh more** (estimates): the ancestry sites ~2 GB as loaded (1 GB if
   the string were trimmed), strain alleles from every genome at 16 per gene ~2.3 GB, and −0.7 GB without index
   alleles. That is ~41 GB resident and ~48 GB peak. The alleles must load before aligning (~5-7 s on one core), and
   the column weights keep the start-up at ~12-14 s.
7. **Cohorts keep every sample's evidence.** After a sample's outputs are written, its slot still holds a per-taxon
   evidence map: ~1 KB for every taxon with a record *or a failed candidate*, which at r226 is ~143,600 taxa per sample.
   Massif finds it the largest live item after a 24-sample cohort's profiling. Freeing it cut that cohort's kept
   memory 0.57 → 0.45 GB. The r226 builds' 1,032-sample run held +61 GB after aligning.
8. **The database build**: congener gaps are 19.5 of the r226 training database's 25.8 min of `--build` (on the
   build's path). The finished database's copy (1:15 h) now delays the build's end by 11 min. The cut proposed on
   2026-10-09 is still open. The collections' `.sam.gz` SAMs cost ~3 min of copying per r226 test collection.

## What I would do, in order

| # | change | where | exact | expected at r226 (pe 51.0 s, HiFi 24.8 s now) |
|---|---|---|---|---|
| 1 | **Load the tables that only profiling reads beside the alignment**: column weights, ancestry sites, congener gaps, species neighbours, priors, conservation, suspect copies, NCBI taxids. Join them when profiling starts. The taxonomy, models, strain alleles and gene neighbours stay before aligning, which needs them (`--weighted_site_shift` would need the column weights too) | `RunProtal.h:3259-3292`; the profiling stage's start | yes | start-up 14.7 → ~3 s on the v23 database (~5-7 s with 0.8.1's alleles until item 2): **pe and HiFi −10 to −12 s each**. ~10 loader threads share the first seconds of aligning with its 32 |
| 2 | **Strain alleles loaded fast**: a flat binary member, or parsed in parallel by species. They stay on the path before aligning | `StrainAlleles.h` (`Read`, `FromRows`), the bundle | yes | v23's 4.2M alleles ~2-3 s, 0.8.1's ~11M ~5-7 s → < 1 s; otherwise they cap item 1 |
| 3 | **Column weights compact and lazy**: a binary member, the codes two to a byte and the consensus in 2 bits, read without parsing. Better still, a family decoded only when a touched copy needs it (an offset index, the bytes kept compressed or mapped) | `ColumnWeights.h` (`Family`, `Read`, `Write`), `Build.h` | yes (the same codes) | −1.1 GB (packed) to −2.5 GB (lazy) resident, −7 GB of start-up peak, ~12 s of one core |
| 4 | **Every table read once**: `ReadAll` reserving the member's size (zstd's frame content size); parsed from the string, not through an `istringstream` copy; the string trimmed where a table keeps it (`ancestry::Table::Read`) | `Database.h:280-305`, the tables' `Read`, `AncestrySites.h:841-860` | yes | each table's transient about halved; the ancestry sites ~1 GB less |
| 5 | **A sample's record evidence freed once its outputs are written** (nothing after reads it) | `RunProtal.h:1525`, `Profiler.h` (`MicrobialProfile`, `RecordEvidenceCollector`) | yes by construction (outputs not compared here) | cohorts: ~140 MB less kept per r226 sample |
| 6 | **Rebuild r226 with 0.8.1** (`ancestry_sites.bin`, no `--index-alleles`) and run `measure_performance.sh` on it, after 1-4 | the build | – | profiling pe 10.6 → ~5-6 s; it shows what items 1-4 leave |
| 7 | **Congener gaps**: the sketch cut for far pairs (2026-10-09, item 8), and the training build's distances reused by the finished build | `CongenerGaps.h`, `Build.h:1394-1466` | cut: if safe; reuse: yes | the r226 build −10 to −25 min |
| 8 | The collections' SAMs as `.sam.zst` | `collect_training_data.py:1506,1517` | yes | ~−3 min per r226 test collection |
| 9 | Each table's time and size in the log, and the tables timed apart from the preload (the timer "Loading the taxonomy, models and tables" holds the preload) | `RunProtal.h` | – | the next such regression shows in the log |

Items 1 and 2 together bring a run back near 0.7.8's start-up; with item 6's database the pe sample would be ~35 s.
Items 3-5 are the memory. Items 1-5 change no output.

## 1. GTDB scale: 0.7.8 against 0.8.1 on the cluster

The 14th cluster run (0.8.1, the r226 v23 finished database) against the 13th (0.7.8's code, the r226 v11 database):
the same node, samples and 32 threads. Medians of runs 2 and 3 (run 1 reads the database from NFS)
([cluster_v14/runs.tsv](results/cluster_v14/runs.tsv), the logs' timers):

| | pe 0.7.8 | pe 0.8.1 | HiFi 0.7.8 | HiFi 0.8.1 |
|---|---|---|---|---|
| **wall** | **33.4 s** | **51.0 s** | **12.3 s** | **24.8 s** |
| gene tables | 0.46 | 0.47 | 0.46 | 0.47 |
| "Loading the taxonomy, models and tables" (the preload inside: 1.7-1.8 s) | 1.72 | **14.25** | 1.75 | **14.33** |
| the index (on a thread of its own, beside) | 2.68 | 2.81 | 2.68 | 2.81 |
| aligning | 25.48 | 23.88 | 7.92 | 7.86 |
| writing the SAM | 0.81 | 0.78 | 0.40 | 0.33 |
| profiling stage | 3.16 | **10.56** | 0.38 | 1.16 |
| of it the sample's timer | 3.06 | 8.61 | 0.36 | 0.84 |
| of that reading the SAM | 1.80 | **6.93** | 0.20 | 0.58 |
| user / system CPU | 955 / 24 s | 1030 / 94 s | 356 / 17 s | 385 / 21 s |
| resident / peak after loading the index | 33.27 / 33.87 GB | **38.49 / 44.66 GB** | 33.01 / 33.63 | 38.49 / 44.66 |
| max RSS | 34.2 GB | **44.6 GB** | 34.5 GB | **44.7 GB** |

The profile-only cohort run of the two samples (no index loaded) went 9.4 → 29.2 s. Its peak went 7.74 → 15.92 GB and
its resident after profiling 7.14 → 9.62 GB: that is the tables alone.

**The alignment stage held.** Per thread, pe: seeding 12.89 → 11.73 s, taking the k-mers 1.95 → 1.49 s (the AVX-512
flex and syncmer scans of `f53db5e` on Zen 4), extending 2.10 → 2.04 s. The alignment handler went 2.79 → 3.25 s, with
the allele scores, the polymorphic sites' settling, the untried and adaptive candidates. HiFi: seeding 2.83 → 2.59 s,
the handler 7.90 → 7.66 s. The work is the same: 176.0M → 176.9M candidates (+0.5%), 10.12M → 10.30M alignments made.

**The start-up is one table.** The tables load each on a thread of its own beside the preload and the index
(`RunProtal.h:3259-3283`), so the bracket is the slowest table. Take away the 1.80 s of preload and 12.4 s remain, which
the builds' runs below date to the column weights. This database has 586,368 family-gene rows of 602,094,588 columns and
13,944,532 mapped copies.

**The profiling is the ancestry sites.** The v23 database was built before `9d4cc3d` and has no `ancestry_sites.bin`
("Ancestry sites: the database has none ...: computed for each gene copy on its first record"). So 0.8.1 builds every
touched copy's sites while it reads the SAM, as 0.8.0 did. The reading thread then waits 4.85 s for room (0.61 s in
0.7.8) because the parsing threads are busy with sites. About 1.9 s of the profiling stage lie outside the sample's
timer (0.1 s in 0.7.8). After that timer come `ReleaseReadData` over the sample's 143,614 taxa, `malloc_trim(0)` and
the freeing of the caches; locally these are 2.4 M instructions, so this world does not show the 1.9 s. System time grew
24 → 94 s for pe (on 32 threads, ≤ 2.2 s of wall); fresh pages for the start-up's transient allocations and the run-time
sites are the likely reason (not measured).

### What the r226 builds' runs showed before

The builds' first test run each, on training databases with 84 threads ([r226_startup.txt](results/r226_startup.txt),
[r226_build_logs.txt](results/r226_build_logs.txt)), date the growth:

| build | code | tables bracket | of it beyond the preload | resident / peak after the index |
|---|---|---|---|---|
| v19 | 0.7.9 + congener gaps | 3.36 s | 0.3 s | 32.21 / 38.55 GB |
| v21 | + strain alleles (4 per gene) | 4.14 | 1.5 | 33.33 / 41.66 |
| v22 (`82bcc15`) | + column weights, index alleles | 17.38 | 10.4 | 36.71 / 42.44 |
| v23 (`8d8c89b`) | | 18.71 | 12.9 | 36.88 / 42.75 |

The 84-thread runs' peaks are above the 32-thread cluster runs' by the index load's per-thread buffers. The builds'
own collections pay the same: v23's three test runs spent 18.7, 62.3 and 29.5 s in the bracket.

## 2. The start-up's tables

**One by one** ([m3_tables.sh](scripts/m3_tables.sh), [tables.tsv](results/local/tables.tsv),
[callgrind_functions.txt](results/local/callgrind_functions.txt)). The experiment build ran with `--sequential_load`,
timed each table and reset the peak (`/proc/self/clear_refs`) before it, on the pe100k sample, profiled, with 4
threads; the times are medians of three runs. The instructions are from the full runs under callgrind. The r226 column
scales them by the units that grow (columns, copies or alleles), and the times also by 0.66, the cluster core's speed
against this one for the column weights (12.4 s measured, 0.680 × 27.8 = 18.9 s here):

| table | arc: rows | instructions | time | resident after | its peak | r226 v23 rows | r226: time | resident | peak |
|---|---|---|---|---|---|---|---|---|---|
| column weights | 21.66M columns, 303,767 copies | 3.04 G (140 a column) | 0.680 s | +90.4 MB | +273.6 MB | 602.1M columns, 13.94M copies | **12.4 s** | **2.5 GB** | **7.6 GB** |
| strain alleles | 131,236 alleles | 0.59 G (4.5k an allele) | 0.082 | +27.2 | +27.2 | 4.22M (0.8.1: ~10-12M) | 1.7-3 (~5-7) | 0.87 (~2.3) | same |
| ancestry sites (0.8.1) | 231,206 copies | 0.36 G | 0.115 | +44.1 | +64.1 | none (0.8.1: ~11M) | (~3.6) | (~2.1) | (~3.1) |
| congener gaps | 254,610 copies | 0.30 G | 0.046 | +9.6 | +9.6 | 12.99M | 1.6 | 0.49 | same |
| the four models | | | 0.149 | +16 | +16 | | 0.1 | | |
| taxonomy, gene neighbours, species neighbours, priors, conservation, suspect copies | | | < 0.07 each | | | | < 1 s | small | |

The r226 resident sum is 4.6 of the measured +5.2 GB: column weights 2.5, strain alleles 0.9, congener gaps 0.5, and the
111M index-allele entries of `--index-alleles 0.01` 0.7 GB. The column weights' 7.6 GB peak matches the profile-only
cohort's +8.2 GB. "Resident after" includes what glibc keeps of a table's freed temporaries.

**Why the column weights cost so much.** Every table loads the same way (`RunProtal.h:245-508`):

1. `DbFile::ReadAll` (`Database.h:280-305`) inflates the member into a `std::string` grown by `append` in 1 MB pieces,
   with no `reserve`. That can take twice the text, and the old buffer lives on while the string grows.
2. The whole text is copied into a `std::istringstream`.
3. `Read` takes it line by line with `std::getline` into a vector of rows.
4. `FromRows` sorts the rows (`std::stable_sort`, with its own buffer) and copies them into the flat table.

The column weights (`ColumnWeights.h:282-346`) have the most text. An F line per family and gene holds four strings of
one character per column: within, among, the consensus base, and the amino acid per codon. That is 3.33 characters a
column, ~2.2 GB at r226. A C line per copy holds its runs. Each line is split into a `std::vector<std::string_view>`,
the hex decoded digit by digit, and each copy's runs kept in a `std::vector<Run>` of its own (13.9M heap blocks) until
`FromRows`. In memory every code takes a byte (`Family`, `ColumnWeights.h:119-129`), though a code is 0-15 and a base
0-4: 3.33 bytes a column where ~1.5 would do. And all of it is read, although a run expands only the copies its records
touch (`GenomeLoader::ColumnWeightsOf`, `column_weights::Cache`). A median r226 design sample's kept records touch ~900
copies and a complex one ~46,000 (the v23 report, section 5), of 13.9 million. In this world's full pe100k run the
table's 3.04 G instructions are 10× the whole profiling of the sample (0.33 G) and 73% of its aligning (4.18 G).

**The ancestry sites (0.8.1)** move the string into the table (`LoadAncestrySites`, `AncestrySites.h:841-860`). The
table erases the index part but keeps the string at the capacity `ReadAll` grew it to: massif shows 33.5 MB allocated
for the 21.1 MB file here. At r226 that is up to ~2 GB for a ~1 GB file. A `shrink_to_fit` there, or `ReadAll`
reserving the member's size, saves the difference.

**The strain alleles** are the one large table the alignment needs (the allele scores), so they stay on the start-up's
path whatever happens to the others. At 4.5k instructions an allele, 0.8.1's finished r226 database (16 per species and
gene from every genome, ~10-12M alleles by the [strain report](../2026-10-11-strain-false-negatives/README.md),
section 6) would take ~5-7 s of one core before the first read is aligned.

## 3. The alignment stage: 0.7.8 against 0.8.1 on the same genes

Each version ran with its own database of the arc world's genes, aligning only (`--no_profile`). The exact counts are
from callgrind on one thread ([m2_callgrind.sh](scripts/m2_callgrind.sh), [callgrind_functions.txt](results/local/callgrind_functions.txt)):

| instructions (inclusive) | 0.7.8 | 0.8.1 | |
|---|---|---|---|
| **pe, 100,000 pairs: aligning** (`RunPairedEnd`) | **4.27 G** | **4.18 G** | **−2.0%** |
| of it seeds and anchors (`ChainAnchorFinder`) | 1.20 G | 1.24 G | +3% |
| the alignment handler | 0.536 G | 0.553 G | +3% |
| the allele scores (`ScoreAlleles`); `--no_allele_scores` on the whole run | – | 5.5 M; −5.3 M | 0.13% |
| the adaptive candidates (`TryCongeners`), the sites' settling (`SettleBySites`) | – | 0.19 M, 0.09 M | |
| the gene neighbours (`--no_gene_neighbours`: aligning) | | −18 M | 0.4% |
| **PacBio, 5 Mb: aligning** (`RunLongReads`) | **6.50 G** | **6.46 G** | **−0.6%** |
| **ONT, 5 Mb: aligning** | **≥ 64.1 G** (stopped) | **6.10 G** | **≥ 10.5×** |
| of it the candidates (`AlignAnchor`: the screen or the bound, then WFA2) | ≥ 63.2 G | 5.25 G | |
| of that the indel bound's LCS (`IndelBound::LcsAbove`); WFA2 | – | 4.30 G; 0.53 G | 70% of aligning |
| a whole pe run: loading the index (each version's file); the genome preload | 22.72 G; 3.40 G | 21.20 G; 3.84 G | |

Timed with 4 threads, medians of three alternated rounds ([m1_align.sh](scripts/m1_align.sh), [align_runs.tsv](results/local/align_runs.tsv)):

| sample | 0.7.8 aligning | 0.8.1 aligning | user CPU 0.7.8 / 0.8.1 | the work (candidates tried; alignments made) |
|---|---|---|---|---|
| pe, 1M pairs | 6.2 s | 6.4 s | 14.7 / 18.1 s | 880,077 → 880,212; 83,966 → 84,089 |
| se, 1M reads | 3.0 | 4.0 | 11.0 / 11.6 | 438,346 → 438,381 |
| PacBio HiFi, 100 Mb | 4.2 | 4.4 | 21.4 / 21.7 | 250,292, the same |
| **ONT, 100 Mb** | **54.9** | **6.2** | **168.6 / 17.0** | 239,459, the same; 227,422 refused by the indel bound |
| ONT, 0.8.1 with `--no_indel_bound` | | 65.7 | 142.9 | |

The short-read and HiFi differences are within this machine's noise (one configuration's user CPU ranged 14-26 s over
rounds); the instruction counts and the cluster run agree they are small. The ONT gain is 0.8.0's indel bound
(`0c94366`), which refuses 95% of the candidates before WFA2. 0.7.8 aligned 239,454 windows to make 4,945
alignments. The bound itself, a bit-parallel LCS with exits every 64 read bases (`AlignmentScreen.h`), is now 70% of
ONT's aligning (`IndelBound::LcsAbove`, 4.30 of 6.10 G): the next ONT lever, if one is wanted (how much of each long
read it scans against a gene window would be the first thing to look at).

## 4. Profiling

The cluster run is in section 1, and the v23 report's section 5 measured the feature groups (the ancestry sites half
of a pe sample's profiling instructions, before `ancestry_sites.bin`). Locally, with the sites stored, profiling is
small: 0.33 G instructions for the pe100k sample (whole run 37.3 G, of which the index 21.2 G and the column weights
3.0 G), 0.14 G for 5 Mb of HiFi. Decoding the stored sites (`AncestrySitesOf`) takes 1.5 M of them.

Profile-only runs of the four samples' SAMs, 4 threads, medians of three ([m5_profiling.sh](scripts/m5_profiling.sh),
[profiling.tsv](results/local/profiling.tsv)):

| | pe | se | PacBio | ONT |
|---|---|---|---|---|
| whole run, default | 2.55 s, 734 MB | 2.21 s, 719 MB | 2.11 s, 704 MB | 2.10 s, 739 MB |
| `--no_site_weights` (the column weights unread) | 2.35 s, 529 MB | 1.89, 530 | 1.21, 504 | 1.53, 505 |
| the sample's own profiling, default / `--no_site_weights` | 0.7 / 1.6 s | 0.5 / 1.0 | 0.2 / 0.5 | 0.3 / 0.5 |
| `--no_gene_neighbours`, `--no_allele_scores` | no change beyond the noise | | | |

Leaving the column weights unread takes their load off the run (−0.2 to −0.9 s, −205 MB of peak here), but the sample's
own profiling gets slower. The build polarised the ancestry sites with the weights, so without them the sites are
computed again: `--no_site_weights` is no switch for the cost.

## 5. Memory: one run, and a cohort

**One run** at r226 is in section 1 (+5.2 GB resident, +10.8 GB peak) and section 2 (which tables). Together with
section 2's estimates, a 0.8.1-built r226 database would be at ~41 GB resident and ~48 GB peak, where 0.7.8 was at
33.3 and 33.9 GB. Items 3 and 4 would take ~4 GB of resident and ~8 GB of peak off that.

**A cohort keeps every sample's evidence.** Each profiled sample stays in a slot until the run ends
(`ProfilingContext::profile_slots`). `ReleaseReadData` (`RunProtal.h:1525`) then erases its taxa, except those that
enter strain MSAs, and `malloc_trim` follows. But the profile's `RecordEvidenceCollector` (`Profiler.h:3595-3617`)
stays, with its per-taxon `RecordEvidence` map (`m_counts`), the per-taxid failed and untried candidate counts, the
read EM's ambiguity classes, and `shared_ptr`s to every touched copy's expanded column weights and ancestry sites.
Nothing after profiling reads it (only the SAM reader does). Measured here
([m4_cohort.sh](scripts/m4_cohort.sh), [cohort.tsv](results/local/cohort.tsv), [m4b_evidence.sh](scripts/m4b_evidence.sh),
[cohort_evidence.tsv](results/local/cohort_evidence.tsv), [massif_cohort24.txt](results/local/massif_cohort24.txt)):

- Resident after profiling: 0.28, 0.36, 0.57 GB for 1, 6, 24 samples (strains on or off alike), ~12 MB per sample.
- Massif on a profile-only run of the 24 samples: the heap grows from ~370 to ~550 MB over the profiling. Its largest
  live item at the end is the evidence maps, 167 MB (30.6%), allocated by `FoldFailedCandidates`. That is ~1 KB per
  taxon, and a taxon gets an entry for a record *or a failed candidate*.
- The experiment build freeing each sample's evidence after its outputs: resident after profiling 0.57 → 0.45 GB for
  24 samples, 0.28 → 0.27-0.28 for one. (The profiles were meant to be compared by checksum, but the script looked in
  the wrong folder; nothing after the outputs reads the evidence, so it cannot change them.)

At r226 the failed candidates of a sample name nearly every taxon: the pe sample's header lists 143,614. That makes
~140 MB kept per sample. The r226 v23 build's training run of 1,032 samples (with `--profile_ahead`) held 97.97 GB after
aligning, +61 GB over its index and tables, and 63.5 GB for 138 samples. A user's cohort of a few hundred r226 samples
would hit the same.

The column weights' run-wide cache (`column_weights::Cache`, up to 200,000 expanded copies of ~4 bytes a base) held
73.5 MB at the end of the 24-sample run, and can reach ~1 GB at r226. It is bounded, so it is no leak.

## 6. The database build

**Congener gaps are most of `--build` at r226** (`Build.h:1394-1466`, `CongenerGaps.h:112-227`). Every copy is aligned
by WFA2 against its congeners' copies: all of a genus up to 24 others, else the 4 nearest by sketch and 16 drawn. The
alignments run over whole genes with both ends partly free and a score budget of 0.4 divergence.

| build | "Congener gaps" | `Run build` | alignments | given up at 0.4 |
|---|---|---|---|---|
| r226 v21 training database | 8:14 | 12:37 | | |
| r226 v22 training database | 12:46 | 20:35 | | |
| r226 v23 training database (84 threads, alone on the node) | **19:29** | 25:23 | 84.4M | 6.36M (7.5%) |
| r226 v23 finished database (idle class, beside the profiling) | **1:15:36** | 1:27:13 | 94.4M | 7.05M |
| arc finished database (4 threads) | 2:36 | 4:16 | 1.28M | 60,911 (4.8%) |

The training database's build is on the r226 build's path, since the simulations wait for it: 19.5 of its 25.8
minutes. The finished database computes the gaps again over most of the same copies' pairs. At v23 the build then
waited 11:16 for it at the end (00:17:43 → 00:28:59 of the 2:05:44 run), so it is on the path too. The 2026-10-09
build-parallelism report proposed giving a pair far beyond 0.4 by its sketch distance its 0.4 without WFA2 (its item 8;
WFA2 runs to its whole budget on exactly those pairs). That is not done yet. A second, exact lever: the two builds align
the same reference copies, so the training build's distances, by gene and pair of copies, would spare the finished
build most of its alignments.

**The collections' SAMs are `.sam.gz`** (`collect_training_data.py:1506,1517`). Only zstd SAMs leave room for the header
(`SamFile.h:225-233`). A BGZF SAM's records go to a side file and are copied behind the header when the sample ends.
The r226 v23 test collection copied 12.8 GB in 184 s for 216 samples (70 MB/s on the node's busy SSD; v22 208 s), and
the training collection 30.9 GB in 29 s. `scripts/compressed.py` reads zstd already. Users get `.sam.zst` by default.

## 7. Confirming on the cluster

After items 1-4, or to see a v24 database as it is, the same command as the 13th and 14th runs (the user's job; the
paths are the cluster's):

```bash
sbatch --exclusive -c 32 --mem 100G --wrap "PROTAL=build/protal bash scripts/measure_performance.sh \
    runs/Perf0.8.x /hpc-home/hildebra/DB/protal/<r226 build>/protal_db/database.protal \
    pe:.../PEARL_AGE/00_fastq/B10_R1_001.fastq.gz:.../B10_R2_001.fastq.gz \
    pb:.../m64319e_230622_130841.bc1002--bc1002.hifi_reads.fastq.gz"
```

With `PROTAL_ARGS=--sequential_load` the log separates the preload from the tables.

## 8. Follow-up: items 1-5 implemented

On `2caffd7`, and checked again on `0dea1ff` after the Singletrack merge (tests below). What changed:

1. **The profiling's tables load beside the alignment** (`RunProtal.h`, `DeferredTables`). In a concurrent start-up
   that aligns, gene conservation, suspect copies, species priors, NCBI taxids, species neighbours, congener gaps,
   foreign rates and ancestry sites load on their threads while the index loads and the reads are aligned. The first
   sample to be profiled waits for them (`ProfileSample`, also with `--profile_ahead`), and then the log says
   "Loading the profiling's tables beside the alignment took X s; the profiling waited Y s for them" with their
   messages. The alignment still waits for the taxonomy, the models, the strain alleles and the gene neighbours. It
   also waits for the column weights when short reads are aligned with the weighted site shift, the default since
   `353b88b` (two hours before this follow-up). So for pe and se the column weights stay on the path, and item 3 is
   what makes them cheap there; long reads defer them. `--sequential_load`, `--profile_only` and `--build` load as
   before.
2. **Strain alleles fast**: `strain_alleles.bin`, a binary member of `database.protal` holding the table's four
   flat arrays (`BinaryTable.h`: magic, version, the text's size, the arrays at 8-byte offsets), read straight into
   them on several threads. The text is parsed straight into the flat arrays too (the row-by-row path stays for
   unsorted or duplicate lines).
3. **Column weights compact**: `column_weights.bin` the same way. In memory the codes are packed: one byte for
   within | among << 4, the consensus base in 2 bits, a no-base bit, the amino acid in 4 bits per codon, so ~1.54
   bytes a column where there were 3.33. A run decodes a copy's columns when a record touches it (`Expand`), as
   before. The text is parsed into the packed arrays directly (fixed fields, no vector per line, sorted once if
   needed). The table is held whole, not lazily: the audit's "better still" (families decoded from compressed bytes
   on demand) was not needed for the start-up.
4. **Every table read once**: `DbFile::ReadAll` reads into a buffer of the member's size (its frames' sizes, or a
   seekable file's seek table), frame by frame on up to 8 threads. Tables parse it through `db::TextStream`, a
   stream over the buffer, instead of an `istringstream` copy. The ancestry sites keep the file's bytes and read
   their index in place (no erase, no doubled capacity). A binary member's header is decoded from its first block
   only (`zstd::ReadFramesPrefix`): the first version decoded the member's whole first frame (up to 64 MB) for each
   of its two header reads, which callgrind caught (1.40 G instructions for the two members' headers; 2.4 M now).
5. **A sample's record evidence is dropped once its outputs are written** (`MicrobialProfile::ReleaseEvidence`, after
   `ReleaseReadData`).

Where the binary members come from: `--build` packs them beside their text tables. `--compress_db` adds them to an
existing single-file database (or remakes stale ones), rewriting `database.protal` once with the other members'
frames copied as they are. `--add_tables` writes them with the tables it stores, and `--unpack_db` leaves them out, as
a folder reads the text. A run uses a binary member only if its header records the text member's size as it is now;
otherwise it reads the text and says why in the table's log line. On the arc database `--compress_db` took 14.3 s (4
threads) and added 11.9 MB.

Found while measuring, and fixed:

- **`--no_profile` aligned without the column weights.** At HEAD, the column weights were read only when the run
  profiles, so a `--no_profile` run of short reads aligned without that part of the weighted site shift, and its SAM
  differed from a profiled run's: 2,258 of the pe100k SAM's lines (sorted records). They are read whenever short
  reads are aligned with the shift now. The SAM records of ref's profiled run, work's profiled run and work's
  `--no_profile` run are the same ([sam_checks.txt](results/impl/sam_checks.txt)). This changes `--no_profile` SAMs,
  for example `foreign_rates.py`'s tiles; a SAM aligned by itself and profiled later now matches a one-step run.
- **`--ncbi_taxids` without the table could exit 8 instead of 2.** Its `exit(2)` ran on a loader thread while the
  ancestry sites' loader was reading, which then failed on torn-down state and called `exit(8)`. One e2e run showed
  it (`test_ncbi_taxids`, which had passed before); HEAD has the same race. The check now runs on the main thread
  before any thread reads the database. Five reruns passed. The other loaders still `exit(8)` from their threads on a
  corrupt table, as before.
- The deferred tables' "took" first timed up to the join (the index load plus the aligning). It now times up to
  the last table loaded.

**Tests** (WSL builds of the changes on `2caffd7`; the last build has every change): unit tests 505 run, 502 passed,
3 skipped. e2e 150 OK, on the last build and on the one before. The mini-database and GTDB-build Python tests:
113 OK in 358 s, on the build before the `--ncbi_taxids` and header fixes (neither touches the build scripts' paths).
ASan + UBSan unit tests (CI's sanitizer job: Debug, `PROTAL_NO_CLONES`, leak checks): 505 run, all passed but
the AVX-512 one CI skips too, no sanitizer report. Then HEAD moved to `0dea1ff` (the Singletrack merge, which touches
`RunProtal.h`, `Build.h`, `Options.h` and `Database.h` too). Rebuilt there, HEAD and HEAD with the changes: unit tests
514 run, 511 passed, 3 skipped; a fresh mini database and e2e 154 OK; the mini-database and GTDB-build Python tests
114 OK; and the identity runs below the same again ([identity_merged_0dea1ff.txt](results/impl/identity_merged_0dea1ff.txt)). `measure_performance.sh` writes and reads the three new columns (44 columns; one
pe and one HiFi run).

**Identity.** Every output is byte-identical to HEAD's (the SAMs' records compared sorted, `misc/` left out):

- on `~/realworld_s1`'s database (the small real world: 170 species, 530 real genomes): pe, se, PacBio and ONT, from the text tables,
  from the binary members, and with `--sequential_load`; a 4-sample cohort map (205 files) with and without
  `--profile_ahead` ([identity.txt](results/impl/identity.txt));
- on the arc database: the pe1m sample's profile-only outputs, from text and binary; the pe100k SAM of a profiled run.

**Exact costs** of the loaders, profile-only runs of the pe1m SAM on one thread with `--sequential_load`, under callgrind
([callgrind_loaders.tsv](results/impl/callgrind_loaders.tsv); HEAD's and the text column from the build before the
header fix, [callgrind_loaders_before_header_fix.tsv](results/impl/callgrind_loaders_before_header_fix.tsv), whose text
path is the same):

| instructions | HEAD (text) | changes, text | changes, binary | per unit, binary |
|---|---|---|---|---|
| column weights (21.66M columns, 303,767 copies) | 3,018 M | 2,426 M (−20%) | **622 M (−79%)** | 29 a column (140 before) |
| strain alleles (131,236 alleles) | 587 M | 338 M (−42%) | **134 M (−77%)** | 1.0k an allele (4.5k) |
| ancestry sites | 364 M | 351 M | 351 M | |
| congener gaps | 296 M | 290 M | 290 M | |
| both binary headers | | | 2.4 M | |

Most of what is left of the binary loads is zstd decoding. The members are compressed at level 19 in 64 MB frames,
so the arc world's 44 MB of column weights is one frame on one thread, and r226's ~1.2 GB about 19 frames on 8.

**Memory.** The column weights hold 44.4 MB in memory here (2.05 bytes a column with the copies and runs). Profile-only
runs of the same SAM on 4 threads (tables only, no index; medians of three,
[profile_only.tsv](results/impl/profile_only.tsv)): max RSS 737 MB (HEAD) → 661 MB (text) → **569 MB** (binary),
−23%. The start-up bracket beyond the genome preload went 0.46 → 0.02 → 0.01 s. A 24-sample cohort's resident memory
after profiling went 0.61 → 0.41 GB, from the evidence (item 5) and the smaller tables.

**Timed runs** on the arc world: 1M pairs and 100 Mb of HiFi, 4 threads, medians of three alternated rounds
([runs2.tsv](results/impl/runs2.tsv), the quiet round, on the build before the header fix):

| | pe HEAD | pe text | pe binary | HiFi HEAD | HiFi text | HiFi binary |
|---|---|---|---|---|---|---|
| wall | 3.37 s | 3.04 | 3.42 | 3.36 | 2.99 | 3.03 |
| what the alignment waits for: "Loading the taxonomy, models and tables" | 0.77 | 0.43 | 0.53 | 0.96 | 0.44 | 0.37 |
| of it the genome preload | 0.47 | 0.39 | 0.48 | 0.44 | 0.43 | 0.36 |
| the index (beside it) | 0.99 | 0.90 | 1.03 | 1.19 | 0.88 | 0.88 |
| the deferred tables loaded after | – | 0.10 | 0.12 | – | 0.48 | 0.53 |
| the profiling waited for them | – | 0.00 | 0.00 | – | 0.00 | 0.00 |
| aligning | 1.97 | 1.88 | 2.03 | 2.05 | 2.00 | 2.05 |
| max RSS | 3,900 MB | 3,792 | 3,779 | 4,011 | 3,903 | 3,882 |

This world cannot show item 1 in wall time. Its tables were done before its index anyway, so in every variant the
alignment started once the index was loaded, and the walls are equal within the noise (±0.3 s between rounds). What
the runs do show: what the alignment waits for shrinks to about the preload; the deferred tables are loaded 0.1-0.5 s
into the run, and the profiling never waited for them; max RSS is down 110-130 MB. The round on the final build
([runs3.tsv](results/impl/runs3.tsv)) ran while another session compiled (load 4.4-5.1): walls 4.1-7.7 s, too noisy
to compare, but the same brackets (HEAD 1.40-1.51 s, the changes 0.58-0.72 s), deferred tables (0.16-0.62 s, never
waited for) and memory.

**What to expect at r226** (estimates from the counts above; the user's next cluster run will tell):

- Start-up: the alignment then waits for the preload, the index and, for short reads with the weighted site shift,
  the column weights. From `column_weights.bin`, that is ~17 G instructions of mostly zstd decoding, spread over 8
  threads: ~0.5-1 s, where the text took 12.4 s. The strain alleles take < 1 s. So pe and HiFi should start
  aligning after ~3 s instead of ~14.3 s: **about −11 s each** (pe ~51 → ~40 s, HiFi ~25 → ~14 s). The ancestry
  sites' +7.4 s of profiling stay until the database has `ancestry_sites.bin` (item 6).
- The v23 database has no binary members. `protal --compress_db --db <v23>/database.protal` adds them in place,
  copying the other members' frames. Without them, pe with the weighted site shift still parses the text (now ~10 s
  of one core) before aligning. The v23 database's models were trained without the shift, so `353b88b` says to
  profile it with `--no_weighted_site_shift`. Then the column weights load beside the alignment either way.
- Memory: the column weights ~1.2-1.5 GB instead of ~2.5 GB resident, and no 7.6 GB transient (the binary arrays are
  read in place; the text path still holds its text and the packed table, ~3.5 GB, for a moment). A cohort keeps
  ~140 MB less per sample.

## 9. Is every table read only where the run reads it?

Asked after section 8 was committed (`0a1beaa`): are tables loaded only when an option, an output or the model needs
them, the NCBI taxids for example only when the profile gets NCBI columns? Each table against what reads it, from the
code (`RunProtal.h` start-up, the loaders, and every accessor's callers):

| table | what reads it | loaded at `0a1beaa` | now |
|---|---|---|---|
| NCBI taxids | the profile's two NCBI columns | only with `--ncbi_taxids` | the same |
| gene neighbours | the alignment (mates past a gene's end, long reads' neighbouring genes), the adjacency features | aligning or profiling; not with `--no_gene_neighbours` | the same |
| ancestry sites | the ancestry features | profiling; not with `--no_site_weights` | the same |
| column weights | the weighted site shift (inside short reads' allele scores), the weights counts behind `conserved_mismatch_*`, `nonsynonymous_*`, `third_position_share`, `ancestry_agreement_weighted`, small genera's ancestry sites | profiling; before the alignment for short reads with the weighted shift, **also with `--no_allele_scores`**, which leaves the shift unused | before the alignment only when short reads are scored with the alleles, else beside it for profiling |
| strain alleles | short reads' allele scores, the alleles and polymorphic features, the species' shared ancestry sites | **before the alignment whenever a run aligns or profiles**: also for long reads and with `--no_allele_scores`, whose alignment never reads them, and in such runs with `--no_profile`, where nothing does | before the alignment only for short reads scored with them; else beside it when profiling; not at all otherwise |
| suspect copies | the evidence filter | profiling, **read whole even with `--keep_suspect_copies`**, then dropped | not read with `--keep_suspect_copies` |
| taxonomy (for the alignment) | the adaptive candidates' genera, short reads only | whenever `--adaptive_candidates` is on, **also for long reads alone** | only with short reads (profiling and the gene neighbours load it anyway) |
| gene conservation | the conservation features; `--gene_conservation` margins | profiling | the same |
| species priors | the composition's genome sizes (always), the prior features | profiling | the same |
| species neighbours, congener gaps | the `db_*` and `gap_*` features, ancestry sites computed at run time | profiling | the same |
| foreign rates | the foreign features | profiling (an opt-in table) | the same |
| models | profiling, one per read type of the run's samples | profiling | the same |
| index | the alignment | aligning | the same |

**Gating by the model's features was considered and not done.** A model reads only the fields of its PMML mining
schema, so a table whose features it does not use could stay unread. But the r226 v24 pe model
(`trained_model.varimp.tsv`, 92 features) uses features of every table above except the species priors, which the
composition needs anyway, and the foreign rates, which v24 does not have. So it would save nothing on a real database.
And `--profile_truth`'s `.truth_annotated` dump, which the build's collections train on, needs every feature.

**What a profile-only run cannot skip**: the genome preload. A single-file database loads every gene (`LoadAllGenomes`:
3.8 G instructions here, the largest part of a profile-only start-up). `--preload_genomes_off` loads genes on demand,
but only from an unpacked, raw `reference.fna`. Loading only the genes a SAM's records touch, from `database.protal`,
would be the next step for profile-only and cohort re-profiling runs; not done here.

The changes: `RunProtal.h` (`alleles_for_alignment`, `weights_for_alignment`, `LoadSuspectCopies`, the genera), a
log line for strain alleles left unread, and an e2e test (`PacBioTest.test_tables_only_what_reads_them_loads`).
Checked on `0a1beaa` against `0a1beaa` with the change, on the small real world's database (strain alleles, column
weights; it has no suspect copies): the outputs are the same (the SAMs' records sorted) for pe, pe with
`--no_allele_scores`, PacBio, PacBio with `--no_profile`, a PacBio + ONT cohort with `--profile_ahead`, and a
profile-only run with `--keep_suspect_copies` ([tables_only_where_read.txt](results/impl/tables_only_where_read.txt)).
The logs show the strain alleles and column weights loading beside the alignment for PacBio and for
`--no_allele_scores`, and left unread by PacBio with `--no_profile`. Unit tests 514 run, 511 passed, 3 skipped; e2e
155 OK.

## 10. The same run on 8 threads (the 15th cluster run)

The user's next run (SLURM 24167284, node `q512n11`, the same EPYC 9634, database and samples, 0.8.1 at `7d5df4b`,
so **before** sections 8 and 9) used `-t 8` instead of 32, with sylph 1.0.0 beside it
([cluster_v15_8threads](results/cluster_v15_8threads/)). PacBio took only twice as long with a quarter of the threads,
and the logs say why: the start-up does not use the threads. Runs 2 (the database in the page cache):

| | pe, 32 threads | pe, 8 | ratio | HiFi, 32 | HiFi, 8 | ratio |
|---|---|---|---|---|---|---|
| **wall** | **51.0 s** | **126.0** | 2.47 | **24.8** | **47.7** | 1.92 |
| "Loading the taxonomy, models and tables" (the preload inside: 1.8 / 3.1 s) | 14.25 | 13.91 | 0.98 | 14.33 | 13.89 | 0.97 |
| the index, beside it | 2.81 | 9.44 | 3.36 | 2.81 | 9.47 | 3.37 |
| aligning | 23.88 | 91.11 | 3.82 | 7.86 | 29.74 | 3.78 |
| profiling | 10.56 | 18.92 | 1.79 | 1.16 | 2.32 | 2.00 |

- **The start-up is ~14 s at any thread count**: the column weights parsed on one thread (section 2). At 32 threads
  it was 58% of the HiFi run, so only the other 42% could get slower. The aligning scaled almost linearly
  (3.8× for 4× fewer threads).
- **Profiling scales worse** (1.8× for pe): the run computes the ancestry sites (v23 has no `ancestry_sites.bin`), and
  the SAM's reading thread waited 13.0 s for room (4.9 s at 32 threads).
- **At 8 threads, the fixes gain less.** With the column weights off the path, the alignment waits for the index, and
  at 8 threads that takes 9.5 s. Expected with `0a1beaa` and the binary members: a start-up of ~9.5 s instead of
  13.9, so **pe ~121 s and HiFi ~43 s** at 8 threads, against ~40 s and ~14 s at 32 (section 8). After that, the
  index load is what an 8-thread run waits for. It scales with the threads (3.4×), and its 44 GB peak is the run's.
- The first pe run took 311 s: the index (198 s) and the preload (70 s) read the 29 GB database cold from NFS, at
  ~150 MB/s.
- Against sylph on the same 8 threads (medians of three): pe 126.0 s against 69.5 s (1.81× the wall, 2.6× the CPU,
  44.0 GB against 4.5 GB); HiFi 47.7 s against 12.5 s (3.81×, 18.6× the CPU, 44.0 against 3.8 GB). With the fixes,
  ~1.7× and ~3.4× the wall. For HiFi the rest is mostly aligning (29.7 s, the alignment handler 29.4 s per thread)
  and the index.

## 11. 0.8.2 on the v24 database, 8 threads (the 16th cluster run)

SLURM 24186577, node `q512n10`, 0.8.2 (`4a42f3b`: sections 8 and 9, the WFA2 speed patches `aa36fc3`) on the **r226
v24** database (built by 0.8.1 with Singletrack; it has `ancestry_sites.bin`, and the runs read its
`column_weights.bin` and `strain_alleles.bin`, added after the build), 8 threads, the same samples
([cluster_v16_082](results/cluster_v16_082/)). Against the 15th run (0.8.1, v23, 8 threads); runs 2. This is not a
pure A/B: the code and the database both changed.

| | pe 0.8.1 v23 | pe 0.8.2 v24 | HiFi 0.8.1 v23 | HiFi 0.8.2 v24 |
|---|---|---|---|---|
| **wall** | **126.0 s** | **115.5** (−8%) | **47.4** | **36.3** (−23%) |
| what the alignment waits for: the tables (the preload inside, 3.1 s) | 13.91 | **3.13** | 13.89 | **3.13** |
| the index, beside them | 9.44 | 9.03 | 9.47 | 9.05 |
| so the alignment starts after | 13.9 | **9.0** | 13.9 | **9.0** |
| aligning | 91.11 | 92.58 | 29.74 | 25.13 |
| of it, per thread: seed and anchor finding; the alignment handler | 62.4; 12.8 | 62.2; 13.8 | –; 29.4 | –; 24.9 |
| profiling | 18.92 | **12.23** | 2.32 | 1.06 |
| max RSS | 44.0 GB | **37.5** | 44.1 | **36.7** |
| user CPU; instructions | 981 s; 6,687 G | 925 s; 6,104 G | 363 s; 3,445 G | 307 s; 2,747 G |
| the profile-only cohort of both: wall; peak | 40.7 s; 15.2 GB | **21.8 s; 10.2 GB** | | |

- **The start-up did what section 10 expected at 8 threads.** The tables went 13.9 → 3.1 s (the column weights read
  from `column_weights.bin`, 1.47 GB in memory; the deferred tables loaded in 1.1 s, never waited for). But the
  alignment now waits for the index, 9.0 s at 8 threads. So the start-up gained 4.9 s, not the 11 s it gains at 32
  threads, where the index takes 2.8 s.
- **The profiling of pe: −6.7 s**, from v24's prebuilt ancestry sites (reading the SAM 14.4 → 7.8 s per thread).
- **HiFi's aligning: −4.6 s**, the WFA2 patches (the alignment handler 29.4 → 24.9 s per thread).
- **pe's aligning did not change** (+1.5 s). 62 s per thread of its 92.6 s is seed and anchor finding (seeding 44.2 s,
  the k-mers 6.0 s, the k-mer screen 5.4 s), which neither round touched; WFA2 is a small part of short reads. The
  alignment handler gained 1.0 s per thread: v24 has 9.2M strain alleles (v23 4.2M), and the weighted site shift is on.
- **Memory: the peak is 7 GB lower** (44.0 → 37.5 GB), the column weights' load transient gone.
- Against sylph on the same 8 threads: pe 115.5 s against 69.4 s (1.66×, was 1.81×), HiFi 36.3 against 12.5 s
  (2.91×, was 3.81×).

So at 8 threads, a pe run of this sample is now 80% aligning, and two thirds of that is finding seeds and anchors. That
is the lever for short reads now. The next one is the index load (9 s at 8 threads, 8% of pe, 25% of HiFi), which is
now on the start-up's path whenever few threads are used.

## How it was run

Everything in WSL under `~/perf9` (removed afterwards), from the scripts in [scripts/](scripts/), in this order:
`build_versions.sh` (three versions), `arc_build.sh` (the world; stopped after its training database, then
`build_arcdb.sh` for the finished database), `sim_pe.sh`, `sim_long.sh`, `build_tim.sh` (the experiment build), then
`queue1.sh`: the models stored, `db078.sh`, `m3_tables.sh`, `m1_align.sh`, `m5_profiling.sh` (rerun after a path fix),
`m4_cohort.sh` and `m2_callgrind.sh` (0.7.8's ONT run stopped once its long-read aligning passed 64 G instructions);
then `m4b_evidence.sh` and the massif run `m6_massif.sh`. `summarize.sh` and `collect_results.sh` made
[results/local/](results/local/). The r226 numbers come from [results/cluster_v14](results/cluster_v14/),
`local/Perf0.7.5_v13` and the build logs ([r226_log_stats.sh](scripts/r226_log_stats.sh)).

The follow-up (section 8), in WSL under `~/impl` (removed afterwards), from the `impl_*.sh` scripts in
[scripts/](scripts/), which name their scratch paths: `impl_build.sh` (HEAD's tree with the changed files copied in;
`REF=1` builds HEAD too) and `impl_unit.sh`, `impl_identity.sh` (the small real world), `impl_e2e.sh` (the mini
database and the e2e tests), then `impl_perf_setup.sh` (the arc world rebuilt by HEAD's pipeline, stopped after its
conversion; HEAD's finished database `db_ref`; a copy given the binary members by `--compress_db`, `db_bin`; the
samples) and `impl_perf_runs.sh` (the timed runs). Then, chained, `impl_after.sh`, `impl_tables.sh`, `impl_cg_align.sh`,
`impl_final.sh`, `impl_round2.sh`, `impl_post.sh`, `impl_asan.sh` and `impl_runs3.sh`; `impl_collect.sh` made
[results/impl/](results/impl/). Section 9: `verify_commit.sh` (`0a1beaa` from `git archive`: build, unit and e2e
tests; its build is the reference) and `tl_check.sh`, in `~/vc` and `~/tl` (removed afterwards).
