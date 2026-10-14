# Sparse index loads and the key map: where the time goes, and fixes

Date: 2026-10-14. Branch `audit-fixes` at `ce2940e`. The user asked why the sparse (compact) key map looks slightly
slower than the 3.2 GB direct one, and why loading a sparse index (keys and values) looks slower than loading the whole
one; to profile both on 1-2 cores with small examples and to suggest fixes for (a) the sparse load and (b) the step
from a k-mer's core to its slots in the value table. (The core is the k-mer's 15 middle bases, so the key map has
4^15 = 2^30 entries; it is not a 12-mer table.)

Everything was measured with an experiment harness, `keybench` ([scripts/keybench.cpp](scripts/keybench.cpp)), built
against `ce2940e` plus [scripts/experiment.patch](scripts/experiment.patch). The patch holds the prototype fixes
below, each with the old code path kept behind an environment variable for A/B runs in one binary
(`PROTAL_OLD_LOAD=1`, `KEYBENCH_OLD_ADVICE=1`), and phase timers (`PROTAL_LOAD_PHASES=1`). It is an experiment, not
a commit candidate as it stands. Nothing in `src/` was changed in the repository.

**Machine.** WSL2 Ubuntu 24.04 on a Core Ultra 7 258V, shared with other sessions (load average 1.5-13 during the
runs). Every run was pinned to cores 4 (one thread) or 4,5 (two threads). In this VM a page's first touch and a TLB
miss cost far more than on a bare-metal node, so absolute times (and the size of the huge-page effects) are
machine-specific; the comparisons below are always within one run, configurations alternating.

**Data.**

| | entries | cores with values | grid levels |
|---|---|---|---|
| strain world, `~/strainfn/db_grid` (whole) | 10.27 M | 5.55 M (0.52% of 2^30) | yes |
| the same, sparse level 2 (C 4) / level 4 (C 16) | 2.59 M / 0.67 M | 2.06 M / 0.61 M | |
| reduced r226 of the mock benchmark, `~/Mocks/reduced_db_r226` | 16.4 M | 10.25 M (0.95%) | no (whole loads only) |
| synthetic key map at r226's density (`keybench synth 0.092 44`) | 4.35 G slots | 98.8 M (9.2%) | |
| synthetic key map of a C 16 sparse r226 index (`synth 0.055 3`) | 0.18 G slots | 59.0 M (5.5%) | |

Reads: D6331 paired-end `SRR23926922_1` (300,000 reads, 22.4 cores each), D6331 HiFi `SRR13128012` (3,000 reads,
1,563 cores each), the strain world's `deep_*_1` (58,543 reads). The synthetic lookups are 99% cores with values, as at
r226 (99.2%), in "reads" of 27 cores (short) or 3,000 (long).

## Summary

**(a) Loads.** A sparse load used 2.0-2.8x the CPU time of a whole load of the same index (1 thread: 2.9-3.6 s
against 1.3-1.4 s) although it keeps 1/C of the values: it decoded the whole 2^30-key map of every chunk (into a
buffer it first zeroed, for `--compact_keys`), found each key by bisecting that map, and placed each kept entry with
`PutOwned`. Three changes (prototype, identical indexes) cut its instructions by 53-60% and its time by 55-66%; on 2
threads a sparse load now takes 0.44-0.71 s against 0.77-1.32 s for the whole direct load. The whole load with
`--compact_keys` gains as much (15.1 → 7.6 G instructions, −55% time).

**(b) Lookups.** In the key-map step alone, protal's compact map costs 2-15x the direct map's. Three of the causes
are fixable: its huge-page advice is refused (`madvise` on an unaligned address), its rank counts call libgcc's
`__popcountdi2` (protal's baseline build), and the starts of a read's first 8 cores are never prefetched. With all
three fixed, the compact map matches the direct one on the real reads (paired-end 6.6 against 7.0 ns, HiFi 13.3
against 15.3); on the synthetic worst case at r226's density (uniformly random cores, no cache reuse) it stays at
+12-13 ns a lookup (short 16.6 against 3.7, long 32 against 20): its second dependent random line. Separately, the
direct map can read one cache line per lookup instead of two (exact; −20 to −55% of its key-map step for short reads).

## (a) The sparse load

### Before: why it was slower than the whole load

Wall, user and system seconds of one load, medians of 3 alternating rounds
([results/load_base_medians.tsv](results/load_base_medians.tsv); busy machine):

| strain world, 1 thread | wall | user | sys | peak |
|---|---|---|---|---|
| whole, direct key map | 2.76 | 1.27 | 1.56 | 3.09 GB |
| whole, `--compact_keys` | 4.21 | 2.66 | 0.68 | 0.35 GB |
| sparse level 2 (C 4), direct | 5.71 | 3.39 | 2.96 | 3.08 GB |
| sparse level 2, compact | 5.06 | 3.59 | 0.92 | 0.23 GB |
| sparse level 4 (C 16), direct | 8.63 | 3.28 | 5.26 | 3.04 GB |
| sparse level 4, compact | 3.69 | 2.88 | 0.53 | 0.17 GB |

The system time is the 3.2 GB direct key map faulting in (VM page faults, very variable here: 0.5-7 s); the user
time is protal's. Where it went (callgrind, one thread, sparse level 2 with `--compact_keys`, 14.99 G instructions;
[results/callgrind_self.txt](results/callgrind_self.txt)):

1. **The whole key map was decoded for every chunk** (`index_codec::detail::DecodeChunk` writes 2^27 blocks, 3.2 GB,
   whatever the load keeps): into `m_keymap` for a direct load, which the pointer pass then rewrote, or, with
   `--compact_keys`, into a per-worker buffer that `std::fill_n` zeroed first (3.61 G instructions of `memset`) and
   `CollectUsed` read back.
2. **Each key with values was found by bisection** in that decoded map (`Cursor::BlockOf`, ~21 dependent reads in a
   66 MB buffer per key), because the decoder's sink was told a key's first slot but not which key it was.
3. **Each kept entry was placed with `PutOwned`**: `Locate`, a bisection over the key's filled entries, atomic byte
   reads and writes. 91 ns an entry with the direct map, 47 ns with the compact one (1.43 G instructions).
4. **The direct pointer pass wrote 2^27 blocks one at a time**, a `memcpy` and a `memset` each (5.0 G + 1.6 G
   instructions; 0.55-0.7 s on one thread).
5. Smaller: the compact map's presence lines (134 MB) zeroed and their ranks counted on one thread, the counts through
   libgcc's `__popcountdi2` (0.32 G + 0.35 G).

### The fixes (prototype in `experiment.patch`)

- **The decoder names the key** (`ValueSink::key_at`: the key's index among the chunk's keys, from the filled blocks'
  bitmap) and takes a null key map for a chunk in the split form, which it then only checks. The sparse load and the
  compact whole load decode no key map at all (a chunk stored raw, rare, still goes the old way); the sparse load's
  `Cursor` and the compact load's per-worker buffer and `CollectUsed` are gone.
- **Kept entries are packed straight into their slots**, chunk by chunk, as `PackKey` packs a whole load (their slots
  follow from the kept counts), instead of one `PutOwned` each: 9 ns an entry.
- **The direct sparse pointer pass fills the runs of empty blocks** with `FillEmptyBlocks`' pattern and writes only
  the blocks with kept keys.

**Checked:** every core's slots and every cell (`Locate` and `UnpackKey` over all 2^30 cores, hashed;
[scripts/check.sh](scripts/check.sh)) are identical between the old and the new path, whole and sparse (levels 0, 2, 4),
direct and compact, on 2 threads, and on the reduced r226 database.

### After

Medians of 3 alternating rounds, old → new ([results/load_ab2_medians.tsv](results/load_ab2_medians.tsv), the final
code; [results/load_ab1_medians.tsv](results/load_ab1_medians.tsv) is an earlier round without the run fill, on a
quieter machine, with the same picture):

| | 1 thread wall (user) | 2 threads wall (user) | peak |
|---|---|---|---|
| whole, direct (path unchanged: the noise) | 2.17 / 1.52 (1.01 / 0.92) | 1.32 / 0.77 | 3.09 GB |
| whole, compact | 2.08 → 0.92 (1.91 → 0.82) | 1.21 → 0.55 | 0.36 → 0.26 GB |
| sparse level 2, direct | 3.30 → 1.23 (2.24 → 0.89) | 2.05 → 0.71 | 3.07 GB |
| sparse level 2, compact | 2.48 → 0.94 (2.14 → 0.85) | 1.51 → 0.52 | 0.23 GB |
| sparse level 4, direct | 3.10 → 1.41 (2.12 → 0.76) | 1.54 → 0.61 | 3.03 GB |
| sparse level 4, compact | 2.27 → 0.78 (2.07 → 0.71) | 1.16 → 0.44 | 0.17-0.28 GB |
| reduced r226, whole, compact | 2.52 → 1.19 (2.29 → 1.05) | 1.57 → 0.78 | 0.44 → 0.33 GB |

Instructions (callgrind, one thread; the load is 99% of the program):

| | old | new |
|---|---|---|
| sparse level 2, compact | 14.99 G | 7.05 G (−53%) |
| sparse level 2, direct | 17.34 G | 6.98 G (−60%; 13.00 G before the run fill) |
| whole, compact | 15.07 G | 7.60 G (−50%) |
| whole, direct (unchanged, for scale) | 7.50 G | |

### What is left, and what would help next

Of the new sparse compact load's 7.05 G instructions: the entries composed and handed over by the decoder,
every level's (~2.8 G), zstd (~2.1 G), filtering and keeping them (~0.7 G), the presence lines zeroed and ranked
(~0.7 G, one thread, libgcc popcount), packing (~0.1 G). Not built:

- **Compose only the kept entries.** The grid level sits in the position field's top bits: read it from the position
  planes first and compose (taxon, gene, flags, flex) only the entries of the load's level or more. At C 4 three
  quarters of the composing goes; estimated 15-20% of the sparse load, more at higher C.
- **Presence lines zeroed in parallel and ranked with `popcnt`** (`CountRanks` in a `PROTAL_CLONE_V3` function):
  ~0.1 s off every compact load, whatever the thread count.
- **A sparse run should load the compact key map** once (b)'s fixes are in: it is the faster load (0.44-0.52 s against
  0.61-0.71 s on 2 threads here) and 2.9 GB smaller, and the direct sparse load must still write all 3.2 GB of its map.
- zstd: only a file format that groups the entries by grid level would let a sparse load read less; not worth it
  while decompression is ~30% of the load.

**At r226 (an estimate from the per-item costs, not measured):** `PutOwned` alone cost ~47-91 ns a kept entry, so for
~0.7 G kept entries at C 4 about 35-65 CPU-seconds against ~6 now; with the bisection per key (98.7 M keys, ~18 steps
in a chunk's map there) and the key map decoded and written twice, the old sparse load likely spent 40-80
CPU-seconds more than needed, 1.3-2.5 s of wall time on 32 threads.

## (b) The key map step: core → its slots

`keybench keys` and `keybench synth` time the key map step alone (`Locate` per core, with protal's prefetching), one
thread, layouts alternating per batch of 200 reads, every layout giving the same slots (checked on all 2^30 cores of
the real databases, every 7th of the synthetic ones). Short reads as `Classify` runs them: a read's lookups, then the
next read's keys prefetched and 2 µs of other work; long reads streamed. Layouts:

- **D**: protal's direct key map (24-byte blocks of 8 cores; a lookup reads its block and the next block's base,
  two cache lines for 3 blocks in 8).
- **D2**: the same map; a lookup of the block's cores 0-6 reads its end from its own block, only core 7 the next
  block's base, and prefetches only the line it needs.
- **C**: `CompactKeys` as protal runs it (`--compact_keys`): the presence-and-rank line, then the start (a 24-bit offset
  and a 64-bit base per 64 ranks); keys prefetched 16 cores ahead, starts 8.
- **Cp** / **Cf** / **Ca**: C with the starts of the read's first 8 cores / all its cores prefetched before its loop,
  or all its starts prefetched a read ahead (halfway through the gap, after the keys).
- **Cl**: C with twice the lead (32 / 16). **M** / **Ma**: C with each group's base in its offsets' line.
- **B**: one 64-byte line per 128 cores (presence bits, a 40-bit base, 20 16-bit ends; fuller lines overflow to a side
  table): one line a lookup, 0.54 GB whatever the reference.

### Three reasons the compact map is slower, two of them fixable bugs

1. **Its huge pages are never applied.** `CompactKeys::Advise` calls `madvise(data, ...)` with the pointer
   `aligned_alloc(64, ...)` returned: not page-aligned, so the kernel refuses it (EINVAL) and the presence lines
   (153 MB at r226) and the starts (300 MB) sit on 4 KB pages. With the old advice the process holds exactly the direct
   map on huge pages (`AnonHugePages` 3,147,776 kB, the 3.22 GB map); fixed (page-aligned range, 2 MB-aligned
   allocation, the bases too), the compact map joins it. `Seedmap::AdviseHugePages` rounds to pages and is fine.
2. **Its rank counts call libgcc.** protal is built for baseline x86-64, so `std::popcount` compiles to a call to
   `__popcountdi2` outside the `target_clones` (`PROTAL_CLONE_V3`) functions. In protal's binary the out-of-line
   `Seedmap::Locate` that every compact lookup goes through has 7 such call sites, `LongReadAligner::operator()` 15,
   the grid's `KmerLookupSM::ScoreAtLevel` 7 ([results/popcountdi2_callers.txt](results/popcountdi2_callers.txt)); a
   lookup counts its rank twice (`PrefetchStart`, `Locate`), over up to 7 words each. `keybench_v3`, the same source
   built for x86-64-v3 (`popcnt` inline), measures the cost.
3. **A read's first starts are never prefetched.** For short reads protal fetches the next pair's keys a read ahead,
   but the starts only from the 9th core of a read on (`LookupKmers`: `PrefetchStart(j + 8)`), so the first 8 lookups
   of each read wait for memory one after the other: 8 of 22-27.

And one that stays: a lookup needs two dependent random lines (presence, then start) where the direct map needs one.
A longer lead (Cl) and the bases in the offsets' line (M) change nothing.

### Measurements

ns per lookup, key map step only ([results/keys2.txt](results/keys2.txt), [results/keys3.txt](results/keys3.txt)). Each
row is one run; compare within a row (the machine's state moved between runs).

Short reads (keys fetched a read ahead):

| data | build, huge-page advice | D | D2 | C (protal) | Cp | Cf | Ca | M | B |
|---|---|---|---|---|---|---|---|---|---|
| synthetic r226 | baseline, protal's (refused) | 15.0 | 7.3 | 217 | 159 | 138 | 129 | 209 | 39.3 |
| | baseline, fixed | 12.4 | 6.9 | 133 | 94 | 82 | 77 | 131 | 28.2 |
| | baseline, fixed (later run) | 13.5 | 7.9 | 71.8 | 38.0 | 32.5 | 33.7 | 86.3 | 17.1 |
| | x86-64-v3, fixed | 3.7 | 3.4 | 35.6 | 18.9 | 18.3 | 16.6 | 42.7 | 12.4 |
| reduced r226, D6331 pe | baseline, protal's | 14.3 | 7.6 | 41.2 | 23.8 | 21.6 | 21.0 | 43.2 | 27.0 |
| | baseline, fixed | 16.5 | 7.6 | 36.4 | 23.1 | 21.2 | 19.9 | 38.2 | 45.9 |
| | baseline, fixed (later run) | 7.2 | 4.3 | 18.4 | 9.6 | 9.0 | 8.9 | 20.5 | 18.2 |
| | x86-64-v3, fixed | 7.0 | 4.3 | 13.7 | 7.3 | 6.8 | 6.6 | 14.9 | 16.3 |
| strain world | baseline, protal's | 6.4 | 4.9 | 89.9 | 55.6 | 47.6 | 47.1 | 91.6 | 28.7 |
| | baseline, fixed | 5.3 | 4.4 | 58.6 | 33.4 | 30.2 | 30.1 | 65.0 | 25.2 |
| synthetic C 16 sparse | baseline, protal's | 7.1 | 4.5 | 98.7 | 53.7 | 46.5 | 45.3 | 115 | 12.8 |
| | baseline, fixed | 12.4 | 5.4 | 91.1 | 55.8 | 49.3 | 44.8 | 97.5 | 16.3 |

Long reads (streamed):

| data | build, huge-page advice | D | D2 | C (protal) | Cl | M | B |
|---|---|---|---|---|---|---|---|
| synthetic r226 | baseline, protal's | 40.2 | 32.5 | 193 | 186 | 176 | 56.1 |
| | baseline, fixed | 28.4 | 27.1 | 107 | 107 | 106 | 28.6 |
| | baseline, fixed (later run) | 19.5 | 16.9 | 63.6 | 62.7 | 63.6 | 18.4 |
| | x86-64-v3, fixed | 20.2 | 17.3 | 32.2 | 31.5 | 32.5 | 15.7 |
| reduced r226, D6331 HiFi | baseline, protal's | 21.2 | 16.4 | 39.9 | 28.9 | 29.2 | 22.8 |
| | baseline, fixed | 20.7 | 16.0 | 24.6 | 30.9 | 24.0 | 19.6 |
| | baseline, fixed (later run) | 16.3 | 12.8 | 20.3 | 21.9 | 19.0 | 15.9 |
| | x86-64-v3, fixed | 15.3 | 16.2 | 13.3 | 11.4 | 12.9 | 13.4 |
| synthetic C 16 sparse | baseline, protal's | 36.5 | 27.0 | 119 | 116 | 121 | 33.6 |
| | baseline, fixed | 26.0 | 26.5 | 91.5 | 90.3 | 94.0 | 23.9 |

Reading the tables:

- **The huge-page fix** takes 30-45% off the compact step at r226's density (217 → 133, 193 → 107) and 38% on the
  real HiFi run (39.9 → 24.6).
- **Inline `popcnt`** takes 25-50% off what is left (synthetic r226: 71.8 → 35.6 short, 63.6 → 32.2 long; real pe
  18.4 → 13.7; real HiFi 20.3 → 13.3).
- **Prefetching the first starts** halves the short-read step again (35.6 → 18.9 primed, 16.6 a read ahead; real pe
  13.7 → 7.3 / 6.6, the direct map's 7.0).
- With all three, the compact map matches the direct one on the real short reads and beats it on the real HiFi run
  (13.3 against 15.3: its tables are a fraction of the size, and real reads come back to the same cores). On the
  synthetic map at r226's density, every lookup a random core, it keeps a gap of 12-13 ns (short 16.6 against 3.7,
  long 32 against 20): the second dependent line. Real r226 runs lie between the two.
- **D2** takes ~40% off the direct map's short-read step where it is not at its floor (13.5 → 7.9, 7.2 → 4.3; 17-56%
  over all short-read rows) and 5-23% off streamed lookups with the baseline build (within the noise with x86-64-v3
  on HiFi). It is exact.
- **B** is as fast as the direct map on the synthetic maps, but real cores cluster: k-mers of related genomes differ
  in their last bases and fall into the same 128-core line. 17.5% of the reduced r226's used cores (5.2% of the
  strain world's) are in lines of more than 20, against 1.2% for uniform cores at r226's density, and on the real
  short reads B is slower than the fixed compact map (16.3 against 6.6). It would need the values laid out in a
  scrambled core order at load time. Not worth it now.

How much of a lookup this is at r226: the key map step was 19.5 of ~175 ns per lookup there (the flex cells 105;
[seeding profile](../2026-10-12-seeding-profile/README.md)), with the direct map. With the three fixes the compact
map should cost 0-13 ns a lookup more than the direct one (the real reads' 0, the random cores' 13), ≲ 7% of the
seeding; D2 should save a few ns a lookup with the direct map. Both need an r226 measurement (`measure_performance.sh`, `PROTAL_ARGS=--compact_keys`,
before and after).

## Recommendations

In order of value for the effort:

1. **Load: decode with `key_at`, pack kept entries directly, fill empty runs** (the prototype). Sparse loads −55 to
   −66%, compact whole loads −55%, identical indexes. To land: drop the timers and `PROTAL_OLD_LOAD`, add a unit test
   that `DecodeChunk` with `key_at` and a null key map hands over the same keys and cells as with one (split and raw
   chunks), and keep the existing whole/sparse/compact equality tests.
2. **`CompactKeys::Advise`: a page-aligned range** (and 2 MB-aligned allocations for blocks of 2 MB or more, the bases
   included). A bug: the advice never took effect.
3. **`popcnt` inline in the compact lookups**: put `Rank`/`Locate`/`PrefetchStart` (or `ChainAnchorFinder::LookupKmers`,
   `PrepareLookups` and the long-read loop) into `PROTAL_CLONE_V3` code, or keep each word's count in the presence line
   (one popcount a lookup, the compact-key-map report's open idea); also `CountRanks` at load.
4. **Prefetch the starts of a read's first cores**: in `LookupKmers`, `PrefetchStart` for the first `kPrefetchStarts`
   cores before the loop when the keys were fetched a read earlier; better, a `PrefetchStarts(next pair)` in
   `Classify`'s loops after `FindSeeds`, once the next pair's keys have arrived.
5. **D2 for the direct map**: `Locate` reads the next block's base only for a block's 8th core, `PrefetchKey`
   fetches only the line it needs. Small and exact.
6. Later: compose only the kept entries in a sparse load; parallel presence init; the compact map by default for
   sparse runs once 2-4 are in.

Not recommended: longer prefetch leads or the bases in the offsets' line for the compact map (no gain), the one-line
128-core bucket layout (real cores cluster).

No user-facing option changes, so the website and `docs/` need nothing until these land (then `docs/development.md`
for the load).

## Reproducing

The WSL build tree (`~/sparseprof`) was removed after the runs; to rebuild: `git archive ce2940e`, extract, apply
[scripts/experiment.patch](scripts/experiment.patch), and run [scripts/build_bench.sh](scripts/build_bench.sh) (it
copies `keybench.cpp` to `misc/` and adds the `keybench` and `keybench_v3` targets). The scripts' `SP` paths point
at the session's scratch folder, where the patched headers and `keybench.cpp` were kept: set them to your copies. Runs:
[scripts/load_runs.sh](scripts/load_runs.sh) (baseline loads), [scripts/load_ab.sh](scripts/load_ab.sh) and
[scripts/load_ab2.sh](scripts/load_ab2.sh) (old against new), [scripts/check.sh](scripts/check.sh) (index checksums),
[scripts/cg_load2.sh](scripts/cg_load2.sh) (callgrind; `--toggle-collect` on `*LoadColumns*` is not usable here, as
nested matches switch collection off), [scripts/keys_runs.sh](scripts/keys_runs.sh),
[scripts/keys_runs2.sh](scripts/keys_runs2.sh), [scripts/keys_runs3.sh](scripts/keys_runs3.sh) (key map step),
[scripts/summ.sh](scripts/summ.sh), [scripts/summ_ab.sh](scripts/summ_ab.sh), [scripts/summ_ab2.sh](scripts/summ_ab2.sh)
(medians), [scripts/popcnt_callers.sh](scripts/popcnt_callers.sh), [scripts/collect.sh](scripts/collect.sh). The
databases and reads are the existing `~/strainfn/db_grid`, `~/Mocks/reduced_db_r226` and `~/Mocks/reads`.
`keybench keys DB READS N PAIRED ROUNDS` also prints the bucket layout's line occupancy, so a run on the cluster's r226
database would tell how clustered r226's cores are.

## Follow-up: recommendations 1-5 implemented (2026-10-14)

The user asked to implement 1-5 with tests, the tests kept fast. On `f98254a` (0.8.4 and its docs commit), files
`src/Hash/IndexCodec.h`, `Seedmap.h`, `CompactKeys.h`, `KmerLookup.h`, `src/Core/ChainAnchorFinder.h`, the tests below
and `docs/running.md`, `docs/databases.md`. Unlike the prototype there are no phase timers and no old paths to switch
back to.

1. **Loads.** `ValueSink::key_at` hands each key's cells over with the key's index among the chunk's keys, and
   `DecodeChunk` takes a null key map with it: a split chunk's key map is only checked (its filled blocks from the
   bitmap, a word at a time), a chunk stored raw is read from a copy of its own and its keys handed over one by one, so
   the loads have a single path (the prototype still kept the `Cursor` for raw chunks). The compact whole load marks
   the keys and keeps their starts as it packs them (no per-worker buffer, no `CollectUsed`, which is gone); the
   sparse load keeps each key's entries of its level as before, then writes the direct key map once (its empty runs
   with `FillEmptyBlocks`) or marks the compact one, and packs its entries straight into their slots (no `PutOwned`).
2. **`CompactKeys::AllocateAdvised`**: blocks of 2 MB or more 2 MB-aligned, advised over their whole pages
   (`AdvisedRange`) before they are touched: the presence lines, the offsets and, new, the bases.
3. **popcnt**: `ChainAnchorFinder::LookupKmers`, `SparsePrepare`, `SparseHits`, `KmerLookupSM::ScoreAtLevel` and
   `CompactKeys::CountRanks` are `PROTAL_CLONE_V3`; `KmerLookupSM::Get`, `PrefetchKey`, `PrefetchStart`,
   `Seedmap::Locate`, `GetPacked`, `PrefetchKey`, `PrefetchStart` and `CompactKeys`' lookup functions are inlined always,
   so that the clones count with `popcnt`. In the built binary the x86-64-v3 clones hold 7-21 `popcnt` each and no call
   of libgcc's `__popcountdi2` ([results/implementation_v3_clones.txt](results/implementation_v3_clones.txt); their
   default clones, for CPUs without x86-64-v3, still call it).
4. **Starts fetched with the keys**: `LookupKmers(kmer_list, keys_fetched)`, `true` from `PrepareLookups` (whose keys
   `PrefetchKeys` fetched a pair earlier), fetches every start of the read before its first lookup and none in the
   loop; `SparsePrepare` likewise for the sampled cores (their keys from `PrefetchSparseKeys`). Long reads keep the
   in-loop prefetches. `LookupKmers`' three loops (grid level, hash pick, all) are one, the same lookups in the same
   order.
5. **The direct key map's `Locate`** reads a core's end from its own block for a block's first seven cores and the next
   block's first slot only for the eighth; `PrefetchKey` fetches the line of the bytes it reads.

**Tests** (new or extended; times on 2 WSL cores):

- `IndexCodec.KeysAreNamedWithoutAKeyMap` (251 ms): four small indexes (empty runs, more cells than a batch, one with
  a chunk stored raw), chunks of 2 KB and 64 MB; every key the decoder names with `key_at` and no key map has the
  index, first slot and cells the index's key map and values say, in key order; without `key_at` a null key map is
  refused.
- `CompactKeys.TheDirectKeyMapGivesEveryCoreItsSlots` (27 ms): every core of a 2^20-core map built from known values
  (whole blocks, a block's last core alone, the map's last core) gets the slots its values take, in core order from
  slot 0, and the others none: an expectation independent of the key map's code.
- `CompactKeys.BlocksAreAdvisedAsHugePagesFromAPageBoundary` (0 ms): the advised range of unaligned and aligned blocks,
  none for a small one; 2 MB of presence lines start on a huge page.
- `CompactKeys.ALoadBuildsTheCompactKeyMapAsItDecodes` (359 ms): now grid levels 0, 2 and 5 on 1 and 3 threads, and
  the same packed entries. With `GridIndex.ASparseLoadHoldsTheEntriesOfItsLevel` (each key's entries, levels and flex
  cells against the file's, ~50 chunks, 3 threads) it covers both new loads.
- e2e `test_compact_keys_give_the_same_records` (4.2 s, +1.4 s): also six long reads of the operon release's genes,
  gated at C 8 (`SparseHits`, then the streamed `LookupKmers`), the same records with and without `--compact_keys`.

All suites from a tree of `f98254a` with these files ([scripts/implementation/build.sh](scripts/implementation/build.sh),
[tests.sh](scripts/implementation/tests.sh)): 585 unit tests passed (40 s), the mini database built, 182 end-to-end
tests passed (223 s).

**protal itself, before and after** (`f98254a` against it with the changes, one thread on one core, 3 alternating
rounds, medians; [scripts/implementation/perf.sh](scripts/implementation/perf.sh),
[results/implementation_perf_medians.tsv](results/implementation_perf_medians.tsv), every round in
[results/implementation_perf_rounds.tsv](results/implementation_perf_rounds.tsv)):

| | before | after |
|---|---|---|
| strain world, sparse index C 4 (grid), "Load Index took" (`--sequential_load`) | 2.36 s | 0.96 s |
| the same with `--compact_keys` | 1.89 s | 0.65 s |
| the same, whole index (path unchanged) | 1.33 s | 1.36 s |
| reduced r226, D6331 paired-end (100k pairs), seeding ns per k-mer | 77.9 | 55.5 |
| the same with `--compact_keys` | 86.0 | 56.0 |
| reduced r226, D6331 HiFi (1,000 reads), seeding ns per k-mer | 58.1 | 47.1 |
| the same with `--compact_keys` | 73.2 | 38.4 |

The seeding times are the `Seeding` row of `misc/<sample>_runtime.tsv` over the k-mers looked up; ungated
(`--sparse_first 0`), 23% of the k-mers in this database's index. The SAM records are the same before and after,
direct and compact, whole and sparse ([results/implementation_same_records.txt](results/implementation_same_records.txt)).
The direct map gains too (one line a lookup, the lookup loop compiled for x86-64-v3). With the fixes the compact key
map is as fast as the direct one for paired-end reads here and faster for HiFi, so it could become the default once an
r226 run confirms it (`measure_performance.sh` with `PROTAL_ARGS=--compact_keys` beside a run without).
