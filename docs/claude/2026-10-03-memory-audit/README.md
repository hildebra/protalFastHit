# Memory audit: what is left to save at GTDB r226 size

2026-10-03. Branch `audit-fixes`, commit `0e4c5a4`. Nothing in `src/` changed; this is measurement
and a ranked assessment. Follow-up to [Where protal's memory goes](../2026-09-30-memory-profiling/README.md)
(2026-09-30) and [The 2-bit gene store](../2026-09-30-gene-store/README.md), whose strategies 1 and 3
(genes at two bits per base, the 32-byte `Gene`) are committed since (`7cc34ed`, `f71e7d5`).

**Question**: audit protal's memory use again and say whether there are meaningful reductions left
for a full GTDB r226 run, without giving up results or much speed.

**Answer in short.** After the 2-bit genes a full r226 run needs about **43 GB** (estimated, not
measured at that size): **35 GB of it is the index's values array**, 3.2 GB the key map, 4.3 GB the
genes, under 1 GB everything else. The values array is the only lever left, and it has slack: every
value is stored in a full 8-byte slot (plus a 4-byte flex cell when a key has several values) although
its fields need 42 bits, and the values of one key are related. Measured on three 900-species worlds,
a **per-key bit-packed layout** (minima and bit widths per key, entries in 23-39 bits, flex cells in a
per-key dictionary) holds the same values in **54-65%** of the bytes, and a lookup microbenchmark shows
it **costs nothing measurable per lookup** because every entry stays randomly addressable. At r226 that
is an estimated **−12 to −17 GB** (43 → 26-31 GB) with identical results. Varint delta coding, the
earlier report's idea, saves as much but makes lookups on long blocks 2-3× slower, which GTDB-size
blocks would turn into a real seeding cost, so it is not recommended. A plain 6-byte entry layout is
the safe half-step (−6 to −7 GB, no speed cost, same engineering of the offsets). Everything else
(key map, genes, per-thread buffers) is 7% or less of the total, or changes results. The figures
are extrapolated from 900-species worlds; `scripts/index_layout_db` prints the exact ones from the
r226 `database.protal` in minutes and should be run before implementing (section 6).

## 1. Setup

| | |
|---|---|
| Machine | WSL Ubuntu 24.04, Core Ultra 7 258V (4 fast + 4 low-power cores), 24 GB; timings are noisy (±15% between repeats), byte counts are exact |
| Worlds | `db060_full` (the 0.6/0.7 benchmark's 900-species world, raw `index.prx`, format 1), `db900n` (the performance reports' 900-species world, `database.protal`), `dbdense` (240 species × 3 genomes, the dense world of the round-2 performance report, `database.protal`); the two single-file databases were decompressed to raw copies with `protal --decompress_db` (`~/r226-build/build/protal`, 2026-10-02) |
| Tools | `scripts/index_layout.cpp` (a raw index) and `scripts/index_layout_db.cpp` (the index inside `database.protal`, decoded chunk by chunk with `IndexCodec.h`, no index in memory) share `scripts/layout_stats.h`: per key the bytes under each layout, the field maxima, the distribution of entries per key, the distinct flex cells and the ubiquitous groups. `scripts/bench_lookup_layouts.cpp` copies an index's values into five layouts and times what `KmerLookupSM::GetFromLookup` does on each. `measure_layouts.sh`, `measure_db_tool.sh`, `bench_run.sh` run them |
| Check | `index_layout_db` on `db900n/database.protal` and `dbdense/database.protal` gives the same numbers, to the byte, as `index_layout` on their raw copies (0.5 s and 0.3 s, 640-700 MB) |
| Not here | a real GTDB database. The 900-species worlds have GTDB-like genus sizes but 160× fewer genes than r226, so their blocks are short (section 4) |

## 2. The budget at r226 size now

From the earlier report's per-gene costs, with the 2-bit genes:

| Part | How | GB | Share |
|---|---|---:|---:|
| index values | 16.6M genes × ~187 entries per gene (`db900n`: 20.5M entries / 110k genes) ≈ 3.1G entries in ~1.3-1.5 slots each (more multi-value keys than at 900 species) × 8 B | **~35** | 80% |
| key map | 2^30 keys × 3 B (16-bit offset + a 64-bit control word per 8 keys), fixed | 3.2 | 7% |
| reference genes | ~17 Gbp × 2 bits | ~4.3 | 10% |
| gene tables | 16.6M × 32 B (`Gene`) + the per-genome vectors | ~0.6 | 1% |
| taxonomy, gene neighbours (~35 MB), conservation factors (2 × 66 MB), species priors, suspect copies, models | | ~0.3 | <1% |
| per thread | ~20 MB (SAM buffers, read batches); plus, only while the index loads, ~75 MB of zstd frame buffers per thread (`zstd::ForEachFrame`: a 64 MB chunk's compressed input and decoded output per worker) | 0.02 × t (+0.075 × t at load) | |
| **total** | | **~43 + 0.1 × t** | |

Cross-check: the r226 index compresses about 2× (18.05 GB in the v1 build at zstd level 19, [r226
evaluation](../2026-10-02-r226-build-evaluation/README.md)), so the raw index (values + key map) is
36-39 GB, as the first two rows say. The exact peak of a real r226 run is in the build script's log
("peak memory" of the training-data protal run, `wait4`'s `ru_maxrss`), which was not at hand here;
the earlier report's 59 GB predates the 2-bit genes. The page cache of the 20 GB `database.protal`
comes on top of this while the file is read (reclaimable, but charged to a cgroup).

So: the key map is fixed, the genes are at their floor (loading them on demand was ruled out: alignment
touches 96% of their pages), the tables are small. **Only the values array matters.**

## 3. What the values array holds

`Seedmap.h`: a key with e values takes S = e + ⌈e/2⌉ slots of 8 bytes: ⌈S/3⌉ slots of 32-bit **flex
cells** (the 16 bases flanking the 15-mer core, one per value, compared with the read's by
`Similarity`), then e 8-byte **entries** (taxid 20 bits, gene 20, position 20, two flags; the fields of
r226 need 18 + 7 + 14 + 2 = 41). A key with one value is one 8-byte entry, no flex cell. So an entry
in a multi-value key costs 12 B (12.7 with the half-slot padding of odd e), a single 8 B. `--build`
drops keys with 2048 or more values (`max_key_multiplicity`), so no key holds more than 2047 entries.

The three worlds (`index_layout`):

| | `db060_full` | `db900n` | `dbdense` |
|---|---:|---:|---:|
| non-empty keys | 10.41M (1.0% of 2^30) | 11.45M (1.1%) | 2.57M (0.2%) |
| entries | 18.88M | 20.49M | 5.78M |
| entries per non-empty key | 1.8 | 1.8 | 2.2 |
| single-value keys (share of entries / of bytes) | 44.3% / 34.2% | 43.5% / 33.5% | 36.5% / 27.5% |
| entries per multi-value key | 5.1 | 4.6 | 7.9 |
| entries in keys of 17-64 values | 12% | 13% | 50% |
| field maxima (taxid / gene / position bits) | 10 / 8 / 12 | 10 / 8 / 12 | 8 / 7 / 12 |
| distinct flex cells among the multi-value keys' entries | 67.6% | 78.0% | 58.5% |
| bytes now | 195.6 MB (10.4 B/entry) | 213.2 MB (10.4) | 61.4 MB (10.6) |

The dense world (three genomes per species, so congeners' k-mers coincide) shows the direction GTDB
takes: fewer singles, longer blocks, repeated flex cells.

## 4. Layouts: bytes

Per entry and as a share of today's bytes, all keys (`index_layout`; singles 6 B in every alternative):

| Layout | `db060_full` | `db900n` | `dbdense` |
|---|---:|---:|---:|
| now: 8 B slots | 10.36 B, 100% | 10.40 B, 100% | 10.61 B, 100% |
| **6-byte entries** + 4 B flex cells, no padding | 8.23 B, **79.5%** | 8.26 B, 79.4% | 8.54 B, 80.5% |
| 5-byte entries + 4 B flex (fits these worlds' 32 bits, not r226's 41) | 7.23 B, 69.8% | 7.26 B, 69.8% | 7.54 B, 71.1% |
| **bit-packed per key**, flex 32 bits: 10 B header (minima, widths), entries in w = w_t + w_g + w_p + 2 bits | 6.91 B, **66.7%** | 7.07 B, 68.0% | 6.46 B, 60.9% |
| **bit-packed, flex dictionary**: + the key's distinct flex cells (4 B each), each entry carries its cell's index | 6.38 B, **61.6%** | 6.76 B, 65.0% | 5.69 B, **53.6%** |
| delta varint (taxid zigzag, gene if changed, position zigzag, flex XOR) | 6.48 B, 62.5% | 6.71 B, 64.5% | 6.24 B, 58.8% |

The bit widths are set per key by the spread of its values: taxids of one key's species, positions of
one k-mer in homologous genes. The dictionary pays where a key's flex cells repeat, i.e. where congeners
share the whole 31-mer; on the dense world it is already the best layout.

**At r226.** Entries per key rise (3.1G entries on at most 1.07G keys, so 3 or more per non-empty key
against 1.8 here), singles fall below a third of the entries, blocks grow to hundreds of entries for
conserved k-mers (capped at 2047), and a key's species are congeners whose internal taxids are
adjacent (the converter numbers species in taxonomy order), whose positions differ by indels only and
whose flex cells are often identical. So per entry: 6-byte layout 10 B (83% of today's 12), bit-packed
with 32-bit flex ~3 + 4 = 7 B (58%), with the dictionary ~3 + 1 B plus the dictionary (45-55%). On
35 GB: **−6 GB, −15 GB, −16 to −19 GB**. These are extrapolations; section 6 gives the tool that
replaces them with r226's own numbers.

## 5. Layouts: lookup speed

`bench_lookup_layouts`: the values of a raw index copied into each layout; 4M lookups of keys drawn in
proportion to their entries (as reads from database species hit them), each doing what
`GetFromLookup` does: the read's flex cell (the key's first cell with one base changed) against every
flex cell of the key (`Similarity`), the entries at the best similarity read (none if more than 256).
Entries in the bit-packed layouts are extracted by bit offset (random access); the delta layout
decodes the whole key sequentially. Nanoseconds per lookup by the key's number of entries (second of
two repeats; `-O2 -mpopcnt`, as the AVX2 clone of `GetFromLookup` has):

| Layout | 1 | 2-8 | 9-64 | 65+ | all (`db060_full`) | all (`dbdense`) |
|---|---:|---:|---:|---:|---:|---:|
| now | 63 | 243 | 507 | 503 | 235 | 227 |
| 6-byte | 31 | 228 | 496 | 560 | 217 | 229 |
| bit-packed, flex 32 bits | 35 | 254 | 519 | 576 | 233 | 227 |
| bit-packed, flex dictionary | 28 | 306 | 564 | 773 | 264 | 269 |
| delta varint | 28 | 405 | 633 | 1433 | 341 | 329 |

(The first repeat and `dbdense`'s bins are in `bench_run.sh`'s output; the spread between repeats is
up to 15%, and most of a lookup is the cache miss on the key and its block, which every layout pays.)

- **6-byte and bit-packed: no measurable difference** to today, in every bin. Extracting a field from
  a bit offset is a shift and a mask; the block is smaller, so it costs fewer cache lines.
- **The dictionary variant is 15-25% slower per lookup here**, where a key's cells are mostly distinct
  (section 3): it does a bit extraction per entry in the scan instead of a 4-byte load, for no fewer
  `Similarity` calls. At r226, where cells repeat, it does one `Similarity` per distinct cell and a
  byte load per entry, which is less work than today's scan; the r226 distinct share (section 6) decides
  between the two bit-packed variants.
- **Delta coding is 40-45% slower overall and 2-3× slower on keys with 65 or more entries**: each
  entry costs ~6 ns of branchy varint decoding, and the whole key must be decoded to reach the matching
  entries. At r226 most entries sit in long blocks, so seeding (14-19% of short-read alignment time,
  [round 3](../2026-10-02-performance-round3/README.md)) would roughly double. A SIMD codec would cut
  that, but not below the fixed-width extraction, which is already free. **Not recommended.**

So the bit-packed layouts take the saving of delta coding without its cost: on the whole run,
seeding being 14-19% of the time, the change is within the noise of a benchmark.

## 6. Measure r226 before building it

`scripts/index_layout_db` reads the index out of `database.protal` chunk by chunk and prints section 3
and 4 for the real database; nothing is held in memory but a few chunks per thread (~75 MB each). On a
cluster node with the r226 database (18 GB to read, a few minutes with 16 threads):

```bash
# in the checkout's root (on WSL from a Linux-side copy: compiling from /mnt/c resolves <zstd.h> to
# protal's Zstd.h); in a conda environment with zstd add -I $CONDA_PREFIX/include -L $CONDA_PREFIX/lib
g++ -O2 -std=c++20 -I src -I src/Utilities -I src/Hash \
    docs/claude/2026-10-03-memory-audit/scripts/index_layout_db.cpp -o index_layout_db -lzstd -pthread
./index_layout_db /path/to/r226/database.protal 16
```

What to read off:

1. **entries, slots, singles**: the exact values-array budget (slots × 8 B) and the 6-byte saving.
2. **bit-packed, flex dictionary** vs **flex 32 bits**: the saving each gives at r226, and from the
   **distinct flex cells** line whether the dictionary also makes the scan cheaper (a low share of
   distinct cells does).
3. **ubiquitous groups**: entries of a key whose flex cell occurs more than 256 times in that key
   (`--max_key_ubiquity`'s default). Such entries can never be returned (`GetFromLookup`: at the best
   similarity their count exceeds the limit, otherwise they are not at the best), and they only matter
   as the maximum, so one representative cell with a count ≥ 257 preserves every result. 0 at 900
   species (32 groups with 8,922 entries, 0.05%, in `db060_full`); at r226 with conserved markers across
   genera of hundreds of species it may be a few percent, worth having only if the tool says so. The
   build would have to keep the full counts for `unique_kmers.tsv` (collapse at load, not in the file).
4. **non-empty keys**: if far below 2^30, a two-level key map would save part of the 3.2 GB; at a
   lookup cost, so only worth a thought if the share is small (the 900-species worlds use 1%, r226 will
   use most of it).

## 7. Ranked options

Savings at r226 from the ~43 GB; results unchanged for all of them.

| # | Change | Saves | Speed | Effort | Notes |
|---|---|---:|---|---|---|
| 1 | **Bit-packed value blocks in memory** (per key: minima, widths, entries in w bits; flex cells as a per-key dictionary or 32 bits), built by `LoadColumns` from the decoded columns; file format unchanged | **−12 to −17 GB** (values 35 → 18-23) | none measured per lookup; the dictionary variant slower at 900 species, likely faster at r226 (section 5) | L | the addressing: per-key 16-bit offsets become byte offsets (a block of 8 keys is ≤ 65536 cells = 512 KB today; compact blocks over 64 KB are rare and need an escape to a side table, or a lower cap at build); `Get`, `GetFromLookup`, `GetExact`, `GetSingleEntry`, `RecoverFromLookup`, `PrefetchValues` adapt; `--build` keeps its 8-byte layout (it runs its own uniqueness check on it and writes the file), so a query run has one layout and a build another; verify with byte-identical SAM and profiles on the 900-species worlds and the e2e tests, speed on `dbdense` |
| 2 | **6-byte entries + 4-byte flex cells** (the half-step: same addressing rework, trivial decoding) | −6 to −7 GB | none | M-L | if 1 is too much at once; the same offsets work must be done, so 1 is the better use of it |
| 3 | Collapse **ubiquitous flex groups** at load (section 6.3), on top of 1 | 0 at 900 species, unknown at r226 | neutral or faster (shorter scans) | S-M on top of 1 | only if `index_layout_db` shows a share worth it |
| 4 | `posix_fadvise(POSIX_FADV_DONTNEED)` on the database file after each member is loaded (`Zstd.h` only advises `SEQUENTIAL`) | up to 20 GB of page cache | none | S | matters for cgroup-limited jobs and nodes shared with other jobs; not in RSS |
| 5 | Cap the index load's per-worker frame buffers (`ForEachFrame` holds a chunk's compressed input and decoded output per worker, ~75 MB; 64 threads ≈ 5 GB on top of the index while loading) | 0.075 GB per thread, at the load peak only | the load uses fewer workers than threads, or frees buffers as it goes | S | or build with `--compress_frame_mb 16` |
| 6 | Cohort (`--map`) retention per sample after profiling (~0.15 GB per dense sample for the strain MSAs, [earlier report](../2026-09-30-memory-profiling/README.md) §1), unchanged since | 0.15 GB × samples | | S-M | only for runs of hundreds of samples |
| 7 | Several protal processes on one node each load the index: run the node's samples as one `--map` run instead (one index, one sample on all threads) | N × 43 → 43 GB | | none (usage) | a shared `mmap` of a raw index would do the same with 4 KB pages, at the TLB cost the THP work removed |

Together 1 + 4 (+ 3 if r226 has the groups) bring a full r226 run to about **26-31 GB** of anonymous
memory, with the database's page cache released as it is read.

### Not worth it

- **Key map** (3.2 GB, 7%): a sparse or two-level map costs a dependent load on every lookup, which is
  what the THP and prefetch work spent its effort on; only if r226 leaves most of 2^30 keys empty.
- **Genes** (4.3 GB): at 2 bits per base; on demand not viable (96% of the pages are touched);
  decoding windows is done.
- **Per-thread buffers**: 20 MB per thread.
- **Fewer syncmers, clade shards, dropping values**: change which reads seed where, hence the results
  and the model's features.

## 8. Reproducing

```bash
bash scripts/measure_layouts.sh      # index_layout on db060_full and raw copies of db900n, dbdense (made with --decompress_db)
bash scripts/measure_db_tool.sh      # index_layout_db from a Linux-side copy of src/, checked against the raw copies
bash scripts/bench_run.sh            # bench_lookup_layouts under ASan/UBSan on dbdense, then -O2 on db060_full and dbdense
```

The scripts assume the WSL paths of the performance reports (`~/protal-perf`, `~/bench071`,
`~/r226-build/build/protal`) and write to `~/protal-mem2`; edit the variables at their top.
