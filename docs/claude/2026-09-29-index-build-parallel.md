# Parallelising `protal --build`: where the time goes, and what a parallel index build would gain

- **Date**: 2026-09-29.
- **Code**: branch `audit-fixes` at `7c6b2f8` (no C++ changed for this report).
- **Machine**: WSL Ubuntu 24.04, Intel Core Ultra 7 258V, 8 threads, 15 GB.
- **Data**: the 765-species release of the GTDB-like tuning world
  ([model-training tuning](2026-09-29-model-training-tuning/README.md)), converted by
  `gtdb_to_protal_db.py`: `reference.fna` 91 MB (85,929 genes), `full_reference.fna` 272 MB
  (257,718 genes of 2,295 genomes). Its marker genes have GTDB's lengths, so sizes scale with the
  species count: GTDB r226 has 188 times the representatives and 319 times the genomes.
- **Method**: `protal --build --no_profile -t 1` and `-t 8`, every log line timestamped as it was
  written (protal flushes each line).

## Where the time goes

| Phase (log line that starts it) | `-t 1` | `-t 8` | parallel today? |
|---|---|---|---|
| pass 1: read, minimizers, count (`Build: iterate records`) | 3.2 s | 3.8 s | no |
| value pointers: prefix sums over 134 M control blocks (`minimizers:`) | 5.0 s | 6.3 s | no |
| pass 2: read again, minimizers again, place values (`After first put`) | 5.6 s | 5.2 s | no |
| uniqueness check against all genomes (`Check Uniqueness`) | 14.1 s | 9.0 s | yes, 1.6x on 8 threads |
| unique k-mer statistics and table (`Save unique kmer info`) | 5.6 s | 6.3 s | no |
| write the index, zstd 19 (`Write ... index.prx.zst`) | 81.5 s | 28.8 s | yes, 2.8x |
| pack `database.protal`, zstd 19, verified | 74.8 s | 59.0 s | partly (64 MB frames: the 91 MB reference has 2) |
| **total** | **195 s** | **123 s** | |

Peak memory 3.7 GB (`-t 1`) and 5.0 GB (`-t 8`); the key map alone is 1.5 GB whatever the size.

At this size compression dominates. At GTDB scale it is the other way round, because the serial
passes grow with the reference and the uniqueness check with all genomes, while the key map is
fixed. Extrapolated linearly (random accesses into tens of GB will be slower still, so these are
lower bounds for the serial parts):

| Phase at GTDB r226, 16 threads | estimate |
|---|---|
| pass 1 + value pointers + pass 2 (serial) | ~30-40 min |
| uniqueness check (1.6x on 8 threads) | ~35-50 min |
| unique k-mer statistics (serial, grows with the values) | ~5-15 min |
| compression and packing (zstd 19; 17 GB reference, values of the index) | ~20-30 min |
| **one build** | **~1.5-2 h** |

If the values grow with the reference (22.6 M here, x188), the index holds ~4 G values (~34 GB),
and a build with the reference preloaded needs ~50-60 GB of memory. The index size of the current
r226 database tells directly.

## Why the index build is serial

`Build.h` forces one thread for both passes ("Building db with multiple threads is currently
broken"), and the code shows why the OpenMP loops around them cannot simply be switched on:

- **Reading.** Each thread's `SeqReader` takes a 1 KB block of the shared stream in an OpenMP
  critical section, so records reach threads in an arbitrary order.
- **Counting.** `FirstPut` calls `CountUpKey`, a plain saturating increment of a 16-bit cell:
  concurrent increments are lost. `CountUpKeyOMP` exists, but its `#pragma omp atomic` wraps
  `x += x < MAX`, which reads `x` in the update and so is not a valid atomic.
- **Placing.** `PutOMP` puts a value into the next empty slot of its key's range under
  `#pragma omp critical(put)`: correct, but all insertions are serialised, and the order of a key's
  values becomes the order threads arrive in.
- **Order matters.** Serially, a key's values are in reference order (by gene, then taxon). With
  threads the index would differ from run to run, and whether alignments are unchanged would have
  to be shown, since anchors are taken from these value lists.

## A parallel index build that gives the same index

Partition the key space by control blocks into P ranges (P = 4-16 x threads, for balance). One
reader cuts `reference.fna` into numbered batches of records. Workers compute the minimizers of a
batch (the expensive, independent part) and split them, in record order, into P buckets. The owner
of a range applies its buckets batch by batch in batch order: `CountUpKey` in pass 1, the value
placement in pass 2, both without atomics or locks, since no other thread touches the range. Each
key's values then arrive in reference order: **the index is byte-identical to the serial one**,
which a test can check (build the mini database both ways, `cmp` the index). The value pointers are
prefix sums: per range in parallel, then an exclusive scan over the P range totals, then the
offsets per range in parallel. The unique k-mer statistics can be counted per range and summed.

Expected gain: the extraction scales with the threads; the updates are random accesses and scale
with the memory system's parallelism (on a 16-core node about 6-12x). The serial ~40-55 min
(passes, pointers, statistics) would become ~5-10 min, per build. The work: a batched,
range-partitioned variant of `Build::Run` and of `Seedmap`'s count and put for one range, plus
the parallel prefix sums; a few hundred lines in `Build.h` and `Seedmap.h`, and the byte-identity
test. Memory: the batches in flight (minimizers of a few MB of reference each).

Two steps are cheaper and worth trying first:

1. **The uniqueness check**, the longest phase at GTDB scale, is already parallel but scaled 1.6x on
   8 threads. Its reader takes 1 KB per critical section, about one gene record, so threads queue
   on the lock; a block of a few MB (one constant in `SeqReader.h`) may give most of the missing
   scaling. Measure before and after; with 1.6x now and ~6x possible, ~20-30 min per build.
2. **Parallel extraction with one applying thread** (producer-consumer, in record order): same
   index, a small change, but it gains only the extraction's share of the two passes (≤2x).

Not worth it: keeping pass 1's minimizers for pass 2 (4 G minimizers x 12 bytes at GTDB scale), or
atomics with a final sort of each key's values (changes the value order and so, possibly, results).

## Already done at the script level (no C++)

- The training database of `build_gtdb_database.py` is packed at zstd level 3
  (`--training-db-level`): it is read only for the training samples and the parity check. On the
  64-species test, its build took 11.0 s instead of 22.6 s (packing 0.5 s instead of 11.1 s); at
  GTDB scale this saves most of its ~20-30 min of compression.
- The two builds still run one after the other (their passes would compete for memory bandwidth,
  and need ~2x the memory together).
