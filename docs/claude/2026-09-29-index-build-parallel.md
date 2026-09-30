# Parallelising `protal --build`: where the time goes, and what a parallel index build would gain

- **Date**: 2026-09-29.
- **Code**: branch `audit-fixes` at `7c6b2f8` (no C++ changed for this report; the follow-up on
  the uniqueness check changes `SeqReader.h`, `Seedmap.h` and `Build.h`).
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
   (Done; the lock on the unique flag, not the block size, was what held it back: see the
   follow-up below.)
2. **Parallel extraction with one applying thread** (producer-consumer, in record order): same
   index, a small change, but it gains only the extraction's share of the two passes (≤2x).

Not worth it: keeping pass 1's minimizers for pass 2 (4 G minimizers x 12 bytes at GTDB scale), or
atomics with a final sort of each key's values (changes the value order and so, possibly, results).

## Follow-up: the uniqueness check (done)

Tried step 1 above. The block size alone helped little: 1 MB instead of 1 KB took 8 threads from
9.6 s to 7.1 s and still left them slower than 4. The lock that mattered was the other one: each
k-mer found in another taxon clears its entry's unique flag inside `#pragma omp critical
(SetNonUnique)`, one global lock taken millions of times. `Entry::SetFlagNonUnique` now clears
the bit with an atomic `fetch_and` (the lock is gone from `Build.h`); clearing a bit commutes, so
the result does not depend on the order. The readers of `value` were already unlocked and read
bits that nothing changes during the check.

Uniqueness check, median of 3 repetitions on the same input (the laptop's times vary by ~10%):

| reader block | flag | `-t 4` | `-t 8` |
|---|---|---|---|
| 1 KB | lock (before) | 7.3 s | 9.6 s |
| 1 MB | lock | 5.4 s | 7.1 s |
| 1 KB | atomic | 5.1 s | 3.9 s |
| 1 MB | atomic (now) | 4.1 s | 2.9 s |

A final run, the commit before (`18783ac`) and a clean copy of it with the change, built the same
way and run back to back on `--build --no_profile --no_compress`:

| | `-t 1` | `-t 4` | `-t 8` |
|---|---|---|---|
| before | 10.3 s | 4.3 s | 9.2 s |
| after | 10.2 s | 3.0 s | 1.9 s |

On 8 threads the check is 5.3x faster than on one (before: slower on 8 threads than on 4). One
thread is unchanged. `unique_kmers.tsv` and `index.prx` are identical in every run (block sizes
1 KB-8 MB, both flag variants, 1-8 threads). With the change: 126 C++ tests, the mini-database
tests and the 85 end-to-end tests pass. The check's ~35-50 min at GTDB scale should become
~10-15 min per build on 16 threads.

## Already done at the script level (no C++)

- The training database of `build_gtdb_database.py` is packed at zstd level 3
  (`--training-db-level`): it is read only for the training samples and the parity check. On the
  64-species test, its build took 11.0 s instead of 22.6 s (packing 0.5 s instead of 11.1 s); at
  GTDB scale this saves most of its ~20-30 min of compression.
- The two builds still run one after the other (their passes would compete for memory bandwidth,
  and need ~2x the memory together).

## Follow-up: the partitioned passes (implemented)

2026-09-30, on `performance` after it merged `audit-fixes` (`c7cba5d`, with the build changes of
[the index build gains](2026-09-30-index-build-gains/README.md)). The design above, as described:

- **Both passes** go through `build::PartitionedPass` (`Build.h`). Each round, one thread reads
  4 x threads batches of whole records (`FastaBatches`, `SequenceUtils/FastaBatches.h`: a batch ends
  before the first record that starts 1 MB or more into it). The threads parse them
  (`ForEachFastaRecord`: the header, and the sequence's lines joined, as `SeqReader` gives them),
  extract the k-mers and sort each batch's items by key range, keeping their order (a counting
  sort). Then each of 64 x threads key ranges (whole control blocks: the main key's top bits) is
  applied by one thread, batch after batch in reference order. Pass 1's items are main keys,
  counted up; pass 2's are the k-mer and its value (the entry's bits as `PutOMP` packs them),
  placed by `Seedmap::PutOwned`: `PutOMP`'s placement without its lock. `PutOMP` now calls it under
  its lock; gone are the static counter it incremented without one and a `Get` whose result it
  dropped.
- **The value pointers** (`BuildValuePointers(threads)`): 16 parts of the control blocks per thread
  are laid out from 0 in parallel, the parts' sizes summed in order, and each part's start added to
  its blocks' control cells; the frequency histogram and counters are summed. The histogram's
  65,536 log lines are no longer flushed one by one.
- `--serial_index_passes` (a dev option) runs the passes and the value pointers on one thread as
  before, for comparison; `--index_batch_kb` sets the batch size (for tests). The warning
  "Building db with multiple threads is currently broken" and the forced single thread are gone.

No two threads write the same memory: a range's keys own whole control blocks, their values (laid
out in block order) a stretch of the value array, and the batches are written and read in phases
with a barrier between them.

**The same index.** On the 765-species tuning world (all species, converted by
`gtdb_to_protal_db.py`), `index.prx` and `unique_kmers.tsv` (md5 starting `fae024569c18` and
`23f64e6a581e`) are the same in all 22 builds: `audit-fixes` (`01fcc55`) and the merge, at 1 and 8
threads; the new code at 1, 4 and 8 threads (several times each), with 64 KB batches, and with
`--serial_index_passes` at 1 and 8 threads. The new end-to-end test `ParallelIndexBuildTest` builds
the mini database with `--serial_index_passes -t 1`, with the parallel passes at `-t 1`, and at
`-t 4` with 16 KB batches (many batches and rounds), and compares both files byte for byte. New
unit tests: `FastaBatches` against `SeqReader` for batches of 1 byte to 64 MB (300 records: several
lines, CRLF, trailing blanks, descriptions, one of 300 kb, no final line end), and a file that is
not FASTA. Unit tests 182/182, end-to-end 103/103, mini-database tests 15/15. Under ASan and UBSan
(CI's sanitizer configuration, Debug), the mini database built on 4 threads with 16 KB batches and
with `--serial_index_passes`: no finding, and the same index. ThreadSanitizer was not run: the build
touches the whole 3 GB key map, and its shadow memory would not fit in 15 GB.

**Time.** The same binary with and without `--serial_index_passes`, alternated, two runs each,
`--build --no_profile --no_compress`, load 0.6-4:

| | serial, `-t 1` | parallel, `-t 1` | serial, `-t 8` | parallel, `-t 4` | parallel, `-t 8` |
|---|---:|---:|---:|---:|---:|
| pass 1 | 1.8-2.3 s | 1.4-1.8 s | 1.3-2.1 s | 0.65-1.0 s | 0.64-0.79 s |
| value pointers | 3.4-3.9 s | 3.1-4.3 s | 3.2-3.8 s | 1.0-2.0 s | 1.1-1.6 s |
| pass 2 | 3.3 s | 2.2-2.7 s | 2.9-3.7 s | 0.73-0.89 s | 0.94-1.0 s |
| **together** | **9.0 s** | **6.6-8.8 s** | **7.4-9.6 s** | **2.4-3.9 s** | **2.7-3.4 s** |
| the whole build | 20.8-24.7 s | 18.2-21.6 s | 13.4-15.3 s | 9.4-10.1 s | 9.0-10.7 s |
| peak RSS | 3.45 GB | 3.48 GB | 3.49 GB | 3.58 GB | 3.75 GB |

On 8 threads the passes and pointers take a third of the time; this laptop (4 fast and 4 slow
cores, shared with other sessions) gains little from 4 to 8 threads. On one thread the partitioned
passes are not slower: a range's updates stay within a few MB of the key map and its stretch of
values. The value pointers are the same work for any database (2^27 control blocks). The batches in
flight cost 0.1 GB at 4 threads and 0.3 GB at 8, whatever the database's size.

At GTDB r226 (an extrapolation: 188 times the reference, and a value array of ~34 GB that makes the
serial passes' random accesses slower than here): passes 1 and 2 took 5.1-5.6 s on one thread here
(times 188: ~16-18 min; the earlier estimate was ~30 min) and 1.6-1.8 s on 8 threads (~5-6 min),
less with more cores; the value pointers stay at ~1-2 s. That is **~5-8 min instead of ~20-30 min
per build**, and the first r226 build's phase timers will tell.

Reproduce: the inputs as in [the index build gains](2026-09-30-index-build-gains/README.md)
(`gtdb_to_protal_db.py --gtdb ~/tune/release_p --outdir in100 -t 8`), then per run, on a copy
`D` of them:
`protal --build --no_profile --no_compress -t T [--serial_index_passes] --db D --reference D/reference.fna --full_reference D/full_reference.fna`,
and `md5sum D/index.prx D/unique_kmers.tsv`.
