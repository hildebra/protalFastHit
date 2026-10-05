# The training database build's memory: where its 64 GB went, and what came off

Date: 2026-10-05. Branch `audit-fixes`, baseline `18dceda`; the change is the commit that adds this report
(`src/Build.h`, `src/RunProtal.h`, `src/SequenceUtils/GenomeLoader.h`, `src/Hash/IndexCodec.h`,
`src/Utilities/Zstd.h`, tests in `tests/test_IndexCodec.cpp` and `tests/test_PackedSequence.cpp`).

**Question.** In the GTDB r226 build on 64 cores, `training_db` took 0:28:39 at a peak of **64.2 GB**
(`database.protal` 23.6 GB). Can that be reduced a little?

**Answer in short.** Yes, by about a third, with every output file byte-identical. Most of the 64 GB
was not the index. The index (36.1 GB in the build's 8-byte layout), the 2-bit genes (4.0 GB) and the
gene tables (0.6 GB) make 40.7 GB. On top of that came two things:

- The suspect-copy scan (`WriteSuspectCopies`, 7m49s at r226) ran with the index in memory. At r226
  scale with 64 threads it peaks at +19.7 GB. When it ends, **12.3 GB of its freed sketches stay
  resident** in glibc's per-thread arenas, and nothing later reuses them.
- The index write and its read-back hold 150-265 MB per worker: two extra copies of each chunk's
  values, and compressed-output buffers zeroed to the full `ZSTD_compressBound`. At 64 workers that is
  ~15 GB, on top of the retained 12 GB.

40.7 + ≤12.3 + ~15 ≈ 64-68 GB, which matches the measured 64.2. After the change:

- the gene-only phases run before the index is allocated, and their freed memory is returned (`malloc_trim`);
- the suspect sketches are moved, not copied;
- the genes are freed after the uniqueness check;
- the writer and its check no longer copy a chunk's values;
- compressed-output buffers are no longer zeroed;
- the index is freed before `database.protal` is packed.

The projected peak is **~42-44 GB**, close to the 40.7 GB floor of this layout. The build's log now
prints `Memory after ...` lines, so the next r226 run will show the real figure.

All sizes below are decimal (GB = 10^9 bytes, as the build script prints "peak memory"). protal's own
log prints binary units under the name GB: its "33.65 GB uncompressed" index is 36.1 GB here.

## 1. Setup

| | |
|---|---|
| The run asked about | the user's r226 build, 64 threads: "built training_db in 0:28:39, peak memory 64.2 GB; database.protal 23.6 GB" (`wait4`'s `ru_maxrss` of `protal --build`) |
| Its phases | `local/v10/training_db_index.log`, the r226 v10 training database: the same size (`database.protal` 21.85 GiB = 23.5 GB), 64 threads, zstd level 3, 64 MB frames. Index 33.65 GiB = 36.1 GB (values 30.65 GiB, key map 3.0 GiB), 13.63M gene copies of 168 genes. Phase times: pass 2 32 s, uniqueness 3m03s, gene conservation 1m34s, suspect copies 7m49s, write index 7m49s, write `database.protal` 5m28s |
| Machine | WSL Ubuntu 24.04, Core Ultra 7 258V (6 cores visible), 23 GB, idle. More threads than cores were run to measure per-worker memory: the buffers exist per worker whether or not it has its own core. Times are noisy (±20%) |
| Harnesses | `scripts/codec_memory.cpp`: `index_codec::Write` and `Verify` on an index of r226's composition in 2^26 keys (9.2% of keys filled, ~29 values per filled key, 47.2 MB of values per 64 MB chunk as at r226; 2.03 GB of values, 41 chunks), compiled against either tree's headers. `scripts/suspect_memory.cpp`: `WriteSuspectCopies`' collection loop and the real `gene_incongruence::Scan` on 81,130 species × 168 genes = 13.63M random sketches (the r226 v10 count). Both sample VmRSS every 5 ms per stage (`rss.h`) |
| Whole builds | the 900-species tuning world (`~/tune/world`, converted by `gtdb_to_protal_db.py`), `protal --build` with `--compress_level 3`, traced every 0.05 s by phase (`trace_build.sh`), baseline vs change (`compare_builds.sh`) |

## 2. Where the memory went (baseline)

**The index write and its read-back** (`codec_memory`, baseline headers). Memory above the index while
the stage runs:

| threads (workers) | write, per worker | read-back, per worker | left resident after the write |
|---|---:|---:|---:|
| 1 (1) | 264 MB | 195 MB | 0 |
| 6 (6) | 234 MB | 195 MB | 0.31 GB |
| 24 (24) | 195 MB | 179 MB | 1.09 GB |
| 64 (41 chunks) | 158 MB | 152 MB | 2.06 GB |

The per-worker figure falls with oversubscription, because workers sharing a core are not at their
peaks together. A 64-core node runs all 64 at once, so the 6-thread figure applies there: about 15 GB
for the write and 12.5 GB for the read-back, plus what the write leaves.

A worker writing a chunk held:

- the chunk's flex cells and entries, copied out of the index (~50 MB);
- the encoded chunk (~47 MB);
- a decoded copy for the writer's self-check (`km_check` and `vals_check`, ~55 MB);
- a `ZSTD_compressBound`-sized output buffer (~47 MB). `std::vector::resize` writes zeros into all of
  it, so the whole bound is resident although the frame needs about half of it.

The read-back decoded each chunk's values into a buffer (~50 MB) before comparing them with the index.

**The suspect-copy scan** (`suspect_memory`, 13.63M copies, as `18dceda` collects them):

| | 64 threads | 6 threads |
|---|---:|---:|
| sketches gathered (per-thread lists copied into the per-gene lists) | 14.4 GB | 15.1 GB |
| peak during `Scan` | 19.7 GB | 15.5 GB |
| resident after every sketch is freed | **12.3 GB** | 7.2 GB |
| after `malloc_trim(0)` | 0.08 GB | 0.00 GB |

A sketch is 128 hashes (512 bytes) in a vector of its own. Each thread's sketches were copied into the
per-gene lists inside a critical section, which doubled them, and the originals were then freed. The
millions of 512-byte blocks freed in 64 arenas stay resident: glibc returns only the top of an arena's
heap, not holes in it. `Scan` adds one sorted table of every (hash, copy) of a gene per thread, about
100 MB per gene in flight.

In the baseline this all ran after pass 2, with the 40.7 GB of index and genes underneath. So the
build peaked twice: in the scan (~60 GB), and in the index write (40.7 + ≤12.3 retained + ~15 ≈ 64-68
GB). Some of the retained memory is reused by the smaller allocations in between, which is why the
measured 64.2 sits at the low end. Packing `database.protal` came after the write, with the index and
genes still allocated; 64 workers compressing the reference held 64 MB of input and 67 MB of zeroed
output each (8.6 GB).

## 3. The change

| | what | effect |
|---|---|---|
| 1 | `build::Run`: the gene conservation factors, gene congeners, suspect copies and the gene-neighbour and gene-position checks run first, before the index passes (none of them reads the index), then `ReleaseFreeMemory` (`malloc_trim(0)` on glibc) | the scan's peak and its 12 GB of freed sketches no longer sit on the index. A bad `gene_neighbours.tsv` or `gene_positions.tsv` now stops the build before its passes |
| 2 | `WriteSuspectCopies`: the threads' sketches are moved into the per-gene lists (`std::make_move_iterator`), not copied | scan peak 19.7 → 13.5 GB and freed-but-resident 12.3 → 6.1 GB at 64 threads; collection 97 → 17 s here (the copies were made inside the critical section) |
| 3 | `GenomeLoader::ReleaseGeneSequences` after the uniqueness check, the last reader of the genes: the preload's arenas and genes loaded one by one are freed; genes keep ids, lengths and counts, and reading one again exits with a message | −4.0 GB for the k-mer statistics, the index write and the packing |
| 4 | `index_codec::EncodeChunk` writes the byte planes straight from the index's cells (`ForEachCell`; `PutPlanes` and the copies are gone), and its self-check compares the decoded values with the index (`DecodeChunk(..., expected)`) instead of decoding into a copy. `Verify` does the same | per worker: write 234 → 116 MB, read-back 195 → 119 MB (6 threads); left after the write at 64 threads 2.06 → 0.23 GB |
| 5 | `zstd::CompressedBuffer`: the compressed-output slots of `WriteSeekable` and `CompressFramesTo` are allocated without zeroing (`new char[]`), so only the bytes zstd writes become resident | a bound-sized buffer costs its compressed size: about half of it for the index, a quarter for the reference |
| 6 | `RunWrapper`: the index (`KmerPutterSM`) is freed when `Run` returns, before `BundleDatabase`/`CompressReference`, and the freed memory returned | packing `database.protal` holds its buffers only: on the 900-species world 4.16 → 0.42 GB |
| 7 | `build::PrintMemory`: "Memory after the gene tables / pass 2 / the uniqueness check / writing the index / writing database.protal: X GB resident, peak Y GB" (VmRSS, VmHWM; protal's binary GB) | the next r226 log shows the peak's phase |
| 8 | `WriteSuspectCopies` set `std::cout` to 2-digit precision and left it so. Harmless where the call was, but before the passes it shortened their lines ("values: 0.2 GB") | precision 6 restored at the end of its line |

The file format is unchanged and the chunks are byte for byte what the copying writer wrote (unit test
`IndexCodec.ChunksAreTheFormatsBytes` checks every chunk against a plain re-statement of the format).
Queries load through the same `DecodeChunk`. Write mode only adds one predictable branch per key: decode
at 1 thread took 2.29 s (baseline, median of 3) vs 2.19 s (change), and at 6 threads 0.79 vs 0.81 s.
Both differences are within this machine's noise.

## 4. Measured after the change

`codec_memory`, the change's headers (baseline in section 2):

| threads (workers) | write, per worker | read-back, per worker | left after the write | write time, baseline → change |
|---|---:|---:|---:|---:|
| 1 (1) | 147 MB | 142 MB | 0 | 12.1 → 10.1 s |
| 6 (6) | 116 MB | 119 MB | 0.05 GB | 4.3 → 4.4 s |
| 24 (24) | 99 MB | 102 MB | 0.09 GB | 5.4 → 5.1 s |
| 64 (41) | 92 MB | 92 MB | 0.23 GB | 7.7 → 5.9 s |

What a worker still holds: the encoded chunk (~47 MB), the compressed frame (random harness values
compress 1.3×, r226's index 1.8×, so ~37 MB here and ~26 MB at r226), the chunk's key map cells for the
check (~5 MB), and the counts. With r226's compression that is ~80-85 MB per worker.

`suspect_memory` with the sketches moved, 64 threads: collected 8.2 GB, `Scan` peak 13.5 GB, 6.1 GB
resident after freeing, 0.08 GB after `malloc_trim`. In the build this happens before the index exists.

The 900-species world at 6 threads (`compare_builds.sh 6`; the key map's fixed 3.2 GB dominates a world
this small):

| phase (`trace_build.sh`) | baseline | change |
|---|---:|---:|
| gene conservation, congeners, suspect copies | 3.63-3.67 GB (after pass 2) | 0.11-0.13 GB (before pass 1) |
| pass 2 | 3.71 GB | 3.73 GB |
| write index | 4.41 GB | 4.22 GB |
| write `database.protal` | 4.16 GB | 0.42 GB |
| VmHWM | 4.47 GB | 4.29 GB |

The log of the change:

    Memory after the gene tables: 0.06 GB resident, peak 0.14 GB
    Memory after pass 2: 3.29 GB resident, peak 3.47 GB
    Memory after the uniqueness check: 3.27 GB resident, peak 3.47 GB
    Memory after writing the index: 3.49 GB resident, peak 3.99 GB
    Memory after writing database.protal: 0.27 GB resident, peak 3.99 GB

## 5. Expected at r226 (64 threads)

From the parts above (index 36.1 GB, genes 4.0 GB, tables 0.6 GB; per-worker figures at 6 threads, as
on a node with a core per worker):

| phase | baseline | change |
|---|---|---|
| conservation, congeners, suspect copies | 40.7 + 19.7 ≈ **60 GB** | before the index: 4.6 + 13.5 ≈ 18 GB |
| pass 2 | 40.7 + batches (~1.3) ≈ 42 GB | ≈ 42 GB |
| uniqueness check | ≈ 41 GB | ≈ 41 GB, then the genes go |
| index write | 40.7 + ≤12.3 retained + 64 × 0.23 ≈ **64-68 GB** (measured peak 64.2) | 36.7 + 64 × 0.085-0.12 ≈ **42-44 GB** |
| index read-back | 40.7 + ≤12.3 + 2 + 64 × 0.195 ≈ 65 GB | ≈ 42-44 GB |
| packing `database.protal` | 40.7 + retained + 8.6 ≈ 50-60 GB | gene tables + 64 × ~0.08 ≈ 6-7 GB |
| **peak** | **64.2 GB measured** | **~42-44 GB** (−20 to −22 GB) |

The finished database is built the same way (`index_and_package.log`, in the background beside the
training database), so its peak should come down by about as much, and with it the node's total
while both builds run. The time should not change: the work is the same, only reordered, and merging
the sketches is cheaper.

## 6. What is left

The rest is the layout. The build holds the index in its 8-byte layout (36.1 GB) because `--build`
sets the uniqueness flags in place and writes from it. Query runs have held it packed since 0.7.5 (27
GB at r226, `Seedmap::PackedLayout`). Building into the packed layout would be the next lever (~−8
GB), but it is a larger change of `Seedmap`'s build path, not a slight one. The genes stay until the
uniqueness check, which reads 6.9M single entries back from them. Pass 2 keeps 4 batches per thread in
flight (~1.3 GB at 64 threads); halving that would save ~0.6 GB at some load balance. A smaller
`--compress_frame_mb` would shrink the per-worker buffers (and the load's) proportionally, but it
changes the files, so it was left alone.

## 7. Checks

- Every file the build leaves (`database.protal`, `gene_incongruence.tsv`, `gene_congeners.tsv`,
  `full_reference.fna.zst`, the tables) is byte-identical between baseline and change on the
  900-species world at 6 threads. Its log has the same lines, sorted, without timings and the new
  `Memory after` lines.
- Unit tests: 352 passed, 2 skipped (a benchmark and an AddressSanitizer-only test), with three new ones:
  - `IndexCodec.ChunksAreTheFormatsBytes`: every chunk against the format written plainly, three
    indexes × three chunk sizes;
  - `IndexCodec.VerifyFindsOneChangedValue`: a flex cell, an entry and a raw chunk's cell changed are
    each found;
  - `GenomeLoaderPacked.ReleasedGenesKeepTheirLengthsAndCannotBeReadAgain`.
- End-to-end tests on the mini database built with the change: 128 of 131 passed. The 3 failures fail
  the same way with the `18dceda` binary, so they predate this change, and they are in query runs, not
  the build:
  - `CompressedDatabaseTest.test_single_file_db` expects older `--compress_db` wording;
  - `DepthKnobsTest.test_the_samples_depth_knob_unless_knob_is_given` matches the alignment summary's
    "Sample sa: ";
  - `OutputFilesTest.test_statistics_for_a_single_sample`: per-taxon statistics are opt-in since 0.7.6.

  The next commit updates the three tests to the current behaviour; then all 131 pass.
- `GtdbBuildTest` (the whole pipeline on a 60-species release, with this protal): 4 of 5 passed (build
  and train, rerun, another seed, background failure, `SIGTERM`). The failing
  `test_d_a_reduced_database_of_the_most_distinctive_genes` expects the step numbers of 9 steps
  ("2/9 the release"); `build_gtdb_database.py` at `18dceda` has 10 ("2/10"). That is drift between the
  two scripts at `18dceda`, and uncommitted edits by another session in the checkout already update
  the test.

## 8. Reproducing

In WSL, from Linux-side copies (the `/mnt/c` checkout is case-insensitive):

    bash scripts/setup.sh                       # ~/buildmem/base: HEAD built; the 900-species world converted
    bash scripts/build_change.sh                # ~/buildmem/new: HEAD + the changed files, with unit tests
    bash scripts/build_harnesses.sh ~/buildmem/base ~/buildmem/h_base
    bash scripts/build_harnesses.sh ~/buildmem/new ~/buildmem/h_new
    ~/buildmem/h_base/codec_memory ~/buildmem/codec.zst 6 26 64     # THREADS KEY_BITS FRAME_MB
    ~/buildmem/h_base/suspect_memory 64 81130 168 copy              # ... move for the change
    bash scripts/compare_builds.sh 6            # both builds traced, outputs compared
    bash scripts/run_tests.sh                   # unit, e2e and GtdbBuildTest on the change

Results (section 7): unit tests 352 passed, 2 skipped; e2e 128 of 131 (the 3 failures predate the change);
`GtdbBuildTest` 4 of 5 (the failure is a stale step-number assertion). Work files in WSL `~/buildmem`
(`runs/` holds the traced builds).
