# Faster `protal --build`: what else costs time at GTDB scale

- **Date**: 2026-09-30.
- **Code**: branch `audit-fixes` at `56217e1` (after the uniqueness-check fix of
  [the index-build report](../2026-09-29-index-build-parallel.md)); no C++ in the branch changed.
  Measured on scratch copies of it: **base** (`56217e1` with timers and counters, `scripts/instrument.py`),
  **+THP** (base with `GLIBC_TUNABLES=glibc.malloc.hugetlb=1`, which puts the index in transparent huge
  pages), and **perf** (base with the `performance` branch's syncmer scan, `8ea37c6`, and its
  `madvise(MADV_HUGEPAGE)` for the index, the `Seedmap.h` part of `bdc5790`).
- **Data**: the 765-species GTDB-like tuning world (`~/tune/release_p` in WSL), converted by
  `gtdb_to_protal_db.py` with all species and with 50% and 25% of them (`--exclude_species`, a
  random half and three quarters left out, seed 1). All species: `reference.fna` 91 MB (85,929
  genes), `full_reference.fna` 272 MB (257,718 genes). GTDB r226 has 188 times the species
  (representatives) and 319 times the genomes.
- **Machine**: WSL Ubuntu 24.04, Intel Core Ultra 7 258V (4 fast + 4 low-power cores), 15 GB. During
  the runs another session ran three `iqtree2 -T 2` jobs (load 4-6), so single times vary up to 2x:
  compare runs next to each other, and trust the counters (exact, the same in every variant).
- **Method**: `protal --build --no_profile --no_compress -t 1` and `-t 8`, every log line
  timestamped on arrival (`scripts/stamp.py`), phases cut at the log lines that start them
  (`scripts/phases.py`); `scripts/scale.sh` runs it all. The prototype (`scripts/proto.py`,
  `scripts/proto.sh`) and the zstd level sweep (`scripts/levels.sh`) did not finish: WSL crashed
  (`Wsl/Service/E_UNEXPECTED`) during their first runs and stayed down.

Regenerate the inputs of the scripts with `git archive -o head56.tar 56217e1`,
`git diff 7c6b2f8 8ea37c6 -- src/SequenceUtils/KmerIterator.h src/SequenceUtils/Minimizer.h > syncmer.patch`
and `git show --format= bdc5790 -- src/Hash/Seedmap.h > huge.patch`, next to the scripts.

## Summary

The earlier estimate of 1-1.5 h per GTDB build (after the uniqueness fix) extrapolated each phase
linearly and missed three things: a whole extra pass (`Check`), the real cost of the unique k-mer
statistics, and two parts that grow faster than the database. A more honest range for the current
code is **~1.5-3 h per build**. The changes below keep the index and `unique_kmers.tsv` byte-identical
unless noted; together they should bring a build to **~45-80 min**, and the parallel passes of the
earlier report to ~30-60 min.

| # | Change | Exact? | Measured here (765 species, 1 thread) | At GTDB r226 (16 cores), per build | Work |
|---|---|---|---|---|---|
| 1 | Drop `build::Check`: a third pass over `reference.fna` that looks up every k-mer and discards the results | yes | 3.6-5.0 s of 36-39 s | -10-20 min (serial; grows like the uniqueness scan) | delete one call |
| 2 | Unique k-mer statistics (`CountUniqueKmers`): count per gene under an integer key, not in four string-keyed maps per value; parallel over key ranges | yes (same row order possible) | 4.2-5.9 s, 95% in the maps (~250 ns per value) | -25-40 min (serial now) | ~50 lines |
| 3 | Its all-pairs loop over each k-mer core's values: only unique values read the result, only "distance <= 1" matters | yes | 0.25 s, but x3.6 per doubling of species | minutes to over an hour now, ~1 min after | ~20 lines |
| 4 | Syncmer scan and huge pages from branch `performance` (`8ea37c6`, `bdc5790`'s `Seedmap.h` part) | yes (verified) | 37.5 s -> 24.5 s (-35%) | -15-25 min | merge |
| 5 | Uniqueness check: compare the flex parts for equality (all it uses) instead of scoring similarity into a vector | yes | not measured (WSL crash) | removes most of the part that grows faster than linear | ~20 lines |
| 6 | `full_reference.fna` without duplicate records (same species, same sequence) | yes | 0.2% duplicates in this synthetic world | unknown: measure on the r226 files (see below) | converter |
| 7 | Build the final database while the training data are collected (`build_gtdb_database.py`) | yes | - | hides one build of the build-and-train run | script |
| 8 | Phase timers and the counters in the build log | - | - | the first r226 build then shows where its time goes | ~20 lines |

Not changed by any of these: compression and packing (zstd 19, ~20-30 min at GTDB; the level sweep
did not run), and the two serial passes (the partitioned parallel build of the earlier report).

## How the phases grow with the database

Counters of the build (the same in every variant and thread count):

| species | values | cores with values | cores with >= 2 values | all-pairs comparisons (statistics) | most values of a core | uniqueness queries | of them single-entry |
|---|---|---|---|---|---|---|---|
| 191 | 4.35 M | 3.20 M | 0.52 M | 7.4 M | 34 | 13.0 M | 6.39 M |
| 382 | 8.72 M | 5.71 M | 1.03 M | 25.9 M | 61 | 26.1 M | 11.1 M |
| 765 | 17.5 M | 10.0 M | 2.12 M | 98.0 M | 118 | 52.5 M | 18.9 M |

Values and queries double with the species, as expected. The all-pairs comparisons grow 3.5-3.8x per
doubling: the more related species, the more values share a k-mer core, and the loop is quadratic in
them. The uniqueness check scans the same blocks once per query, so its work per query grows the same
way (the scanned-cell counter of this run overcounts misses, see `instrument.py`, but grew 3.5-3.7x per
doubling too). GTDB r226 is 7.5 doublings beyond this world. The index drops cores with 2,048 or more
values, which bounds both, but not tightly (at most 2,047 comparisons per value). How far the growth
goes depends on how clades are sampled, which a synthetic world does not reproduce; the first r226
build with phase timers (item 8) will tell.

Phases at 765 species, 1 thread, seconds, two runs each (runs of the three variants alternated):

| phase (log line that starts it) | base | +THP | perf |
|---|---|---|---|
| pass 1 (`Build: iterate records`) | 2.96 / 3.05 | 3.57 / 2.57 | 1.49 / 1.75 |
| value pointers (`minimizers:`) | 7.17 / 5.27 | 2.86 / 2.38 | 2.25 / 2.31 |
| pass 2 (`After first put`) | 4.79 / 4.63 | 4.39 / 4.04 | 3.24 / 2.85 |
| uniqueness check (`Check Uniqueness`) | 11.34 / 13.57 | 13.06 / 10.52 | 9.36 / 8.18 |
| unique k-mer statistics (`Save unique kmer info`) | 4.22 / 5.88 | 5.20 / 5.12 | 4.21 / 4.46 |
| write the raw index (`Write`) | 1.44 / 1.69 | 2.07 / 1.68 | 1.50 / 1.49 |
| `Check` | 3.58 / 4.99 | 4.57 / 4.00 | 2.95 / 2.61 |
| **total** | **35.7 / 39.3** | **35.8 / 30.4** | **25.1 / 23.8** |

Huge pages mainly cut the value-pointer phase (the page faults of the 1.5 GB key map); the syncmer scan
cuts every phase that extracts k-mers (passes, uniqueness check, `Check`) by 25-45%. `index.prx` and
`unique_kmers.tsv` were byte-identical across variants, sizes' repeats and 1 and 8 threads. On 8
threads the uniqueness check took 1.9-2.7 s (all variants; one perf run at 5.2 s under load).

The statistics phase at 765 species: the scan took 4.1-5.8 s, of which the all-pairs loop 0.23-0.26 s
and writing the table 0.07 s; the rest is building `std::string` keys and updating four
`tsl::sparse_map<std::string, uint32_t>` per value (~250 ns per value with 86k genes, more at GTDB's
16 M genes, whose maps do not fit in cache).

## The changes

1. **`Check` does nothing.** `RunProtal.h` calls `build::Check` after the index is written: one thread
   reads `reference.fna` again, extracts the k-mers and calls `GetFlex` for each, and discards the
   results; nothing is compared, counted or printed. At GTDB scale it costs as much as pass 2, and its
   lookups grow like the uniqueness check's. Remove the call (and the function).

2. **The statistics.** Keep one `robin_map<uint64_t, row>` from `(taxid << 32) | geneid` to a row of
   four counters. The row order of `unique_kmers.tsv` matters a little: `GenomeLoader` inserts the
   hittable genes into a `sparse_set` whose iteration order depends on insertion order. Inserting each
   gene's string key once, on first sight, into a `sparse_map<std::string, row>` and writing in its
   order gives the old file byte for byte, at one string per gene instead of four map updates per
   value (`scripts/proto.py` does this; its run did not finish). The scan over control blocks splits
   into key ranges for threads: going through the ranges in order and inserting each range's genes in
   its first-sight order, skipping those already seen, reproduces the serial order.

3. **The all-pairs loop** computes, for each value of a core, its closest other value, but the result
   is only read for values still flagged unique, and only as "distance above 1". Skipping non-unique
   values and stopping at the first neighbour within distance 1 gives the same flags. For blocks where
   many unique values have no close neighbour, an exact linear alternative: insert each flex part with
   one position masked out (16 variants) into a hash set per block; a neighbour within distance 1
   shares one masked variant.

4. **Syncmer scan and huge pages** are on branch `performance` (not merged): `8ea37c6` scans each
   s-mer once per strand without data-dependent branches (k-mers identical, tested against the
   definition); `bdc5790` asks for transparent huge pages for the key map and the values. Merging the
   two gives the "perf" column. The rest of `bdc5790` (threaded gzip reader for FASTQ, timers) is for
   profiling runs.

5. **Uniqueness check.** Of `GetFlex`'s result it uses only whether the best similarity is the whole
   flex part, i.e. equal flex parts. `GetFlex` scores every flex cell of the block and pushes each score
   into a vector; a `GetExact` that compares 32-bit flex parts for equality (vectorisable, no vector)
   gives the same entries (`scripts/proto.py` has it, not yet run). This is the part that grows
   faster than the database.

6. **Duplicate records in `full_reference.fna`.** The check's result depends only on each record's
   species and sequence (clearing a flag twice changes nothing), so identical copies of a gene in
   genomes of one species can be written once. Species with thousands of near-clonal genomes (E. coli,
   M. tuberculosis, S. aureus and others) may share many identical copies; the synthetic world here has
   only 0.2% duplicates, so this needs the real files. On the download node, after `download_gtdb.py`,
   count them in the converted database (unique records against all):
   `awk 'NR%2{h=$0;next}{print h"\t"$0}' full_reference.fna | sort -u -S 50% --parallel=16 | wc -l`
   against `grep -c '>' full_reference.fna`. In the converter it is a set of sequence hashes per marker
   and species in `_write_full_reference` (and `test_full_reference_covers_all_genomes` changes).

7. **Overlap in `build_gtdb_database.py`.** It builds the final database, then the training database,
   then collects the training data and trains. The final database is needed only for `--add_model`
   at the end. Built in the background while the training data are collected, it costs no wall-clock
   time as long as the node has the memory for both (the collection's protal runs hold the training
   index; a build at GTDB scale needs ~50-60 GB).

8. **Timers.** Each phase of `build::Run` as a `Benchmark` in the log, and the counters above
   (queries, cells scanned, pairs), cost nothing and make the first r226 build a measurement.

## GTDB r226 per build, 16 cores (estimates)

| phase | now (`56217e1`) | with 1-5 |
|---|---|---|
| passes 1 and 2 (serial) | 25-40 min | 15-25 min |
| value pointers, zeroing the values | ~1-2 min | ~1 min |
| uniqueness check | 10-40 min (linear part ~7 min, the rest the growing scan) | 5-15 min |
| unique k-mer statistics | 25-40 min plus the all-pairs loop (minutes to over an hour) | 1-3 min |
| `Check` | 10-20 min | 0 |
| compression and packing (zstd 19) | 20-30 min | 20-30 min |
| **one build** | **~1.5-3 h** | **~45-80 min** |

The ranges are wide because the parts that grow faster than the database cannot be extrapolated from
a synthetic world. The partitioned parallel passes of the earlier report would take the first row to
~3-5 min; compression is then the largest phase (the level sweep is still to be run).

## Follow-up: 1-5, 7 and 8 implemented (2026-09-30)

On `audit-fixes` after `56217e1` (uncommitted when measured):

- **1.** `build::Check` is gone (`RunProtal.h`, `Build.h`).
- **2, 3.** `Seedmap::CountUniqueKmers` scans the control blocks in `-t` threads and counts into one
  row per gene of the reference (`Build.h` `GeneRowsOf`: taxon, gene id from `reference.map`), with
  relaxed atomic increments; each gene's first value (the smallest value index, which is the order
  the serial scan met genes in) decides its place: the rows' string keys go into a
  `sparse_map<std::string, uint32_t>` in that order, which iterates as the four maps of before did.
  The distance-two flag is computed only for unique values, stops at the first value within flex
  distance 1, and for cores with 256 or more values sorts the flex parts with each position masked
  (`Seedmap::FlexNeighbours`, n log n; `tests/test_Index.cpp` checks it against all pairs).
- **4.** The syncmer scan of `8ea37c6` (`KmerIterator.h`, `Minimizer.h`, `tests/test_Syncmers.cpp`)
  and the huge pages of `bdc5790` (`Seedmap.h` only).
- **5.** `KmerLookupSM::GetExact` (equal flex parts) replaces `GetFlex` in the uniqueness check.
- **7.** `build_gtdb_database.py` builds the finished database in the background from the start and
  waits for it before `--add_model` (`--one-build-at-a-time` for the old order).
- **8.** Phase timers and the uniqueness check's counters in the build log (below).

**Same output.** `index.prx` and `unique_kmers.tsv` byte-identical to `56217e1`'s at 25, 50 and 100%
of the tuning world, on 1 and 8 threads (`scripts/val_build.sh`). Unit tests 129/129 (three new),
mini-database tests 15/15 (two new), end-to-end 85/85.

**Time** (`--no_compress`, so without zstd; `56217e1` and new alternated, one run each; the machine
was slower than during the runs above, so compare within a row):

| tuning world | `56217e1` | new |
|---|---|---|
| 25%, 8 threads | 21.2 s | 10.0 s |
| 50%, 8 threads | 24.9 s | 11.9 s |
| 100%, 8 threads | 41.9 s | 18.5 s |
| 100%, 1 thread | 70.4 s | 27.7 s |

The log of a new build at 100%, 8 threads:

    Pass 1 (count the k-mers) took 3s 431ms
    Value pointers took 3s 140ms
    Pass 2 (place the values) took 6s 285ms
    Uniqueness check took 2s 181ms
    Uniqueness check: 52497345 k-mers, 288201258 flex parts compared (5.5 per k-mer), 9407639 found under another taxon or more than once, 1314685 single entries read back from their genes
    Unique k-mer statistics took 332ms
    Distance-two flags: 31000477 flex parts compared
    Write index took 2s 847ms

The statistics went from 4-6 s to 0.33 s (1.4 s on one thread); the distance-two flags compared
31 M flex parts instead of 98 M pairs. The phase times of the first tables came from the log's
timestamps, which the 65,536-line histogram printed at the value pointers delays: there, value
pointers and pass 2 are right together, not apart (pass 2 alone is ~6 s, not ~3 s).

Against the "perf" variant of the first tables (syncmer scan and huge pages only), alternated, 100%,
1 thread, two runs each (perf's phases from timestamps, pass 2 and value pointers together):

| phase | perf | new |
|---|---|---|
| pass 1 | 2.9 / 2.4 s | 2.1 / 2.5 s |
| value pointers + pass 2 | 7.5 / 7.5 s | 7.3 / 8.6 s |
| uniqueness check | 11.6 / 13.9 s | 12.9 / 13.4 s |
| unique k-mer statistics | 7.0 / 12.4 s | 1.3 / 1.9 s |
| `Check` | 3.8 / 7.2 s | - |
| **total** | **36.1 / 47.6 s** | **26.8 / 31.0 s** |

`GetExact` (item 5) gains nothing measurable here: 5.5 flex parts per k-mer are too few for the scan
to matter next to the k-mer extraction and the random reads of the index. It is for GTDB's larger
blocks, where the log's "flex parts compared per k-mer" will show how much it matters.

At GTDB r226 the serial passes (pass 1 and 2, ~9 s here on one thread, times 188) are now the largest
part next to compression: ~30 min, which the partitioned parallel build of the earlier report
addresses. With 16 threads the uniqueness check and the statistics should take ~10 min and ~1-2
min, so a build should take ~55-75 min, and the first r226 build's log says exactly.
