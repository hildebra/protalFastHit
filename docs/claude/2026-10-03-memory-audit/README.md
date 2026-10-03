# Memory audit: what is left to save at GTDB r226 size

2026-10-03. Branch `audit-fixes`, commit `0e4c5a4`. Nothing in `src/` changed; this is measurement
and a ranked assessment. Follow-up to [Where protal's memory goes](../2026-09-30-memory-profiling/README.md)
(2026-09-30) and [The 2-bit gene store](../2026-09-30-gene-store/README.md), whose strategies 1 and 3
(genes at two bits per base, the 32-byte `Gene`) are committed since (`7cc34ed`, `f71e7d5`).

**Question**: audit protal's memory use again and say whether there are meaningful reductions left
for a full GTDB r226 run, without giving up results or much speed.

**Answer in short.** After the 2-bit genes a full r226 run needs about **43 GB**: **35.0 GB of it is
the index's values array** (measured on the r226 v7 database, section 9), 3.2 GB the key map, 4.3 GB
the genes, under 1 GB everything else. The values array is the only lever left, and its slack is one
thing: every entry is stored in a full 8-byte slot although its fields need 42 bits (taxid 18, gene 8,
position 14, two flags). The other 4 bytes per entry, the flex cell, are incompressible: 82% of the
cells are distinct within their key, so a per-key dictionary gains nothing and delta coding little.
**Packing the entries at 42 bits** (widths from the index header, a fixed bit stride, flex cells a plain
32-bit array, every entry still randomly addressable) holds the values in **26.8 GB, −8.2 GB (−23%)**,
with identical results and, by a lookup microbenchmark, no measurable cost per lookup; a 6-byte entry
layout is the simpler variant at −6.0 GB. Varint delta coding, the earlier report's idea, saves 9.2 GB
but makes lookups on long blocks 2-3× slower, and at r226 70% of the entries sit in keys of 33 or more
values, so it is out. Everything else (key map at 9% fill, ubiquitous groups at 0.7% of the entries,
genes, per-thread buffers) is worth 1.5 GB or less, or changes results. So the run goes from 43 to
**about 35 GB**, and that is the floor of this representation; the 900-species worlds (sections 3-5)
had suggested more because a 15-mer key there holds a few related values, while at r226 a key mixes 29
values of different genes and positions whose fields span nearly the whole widths.

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

**At r226** (section 9, measured): the picture is different. A key holds 29 entries on average, but
they are not 29 congeners' copies of one k-mer: a 15-mer core collects k-mers of different genes (gene
spread 8 bits within a key) at different positions (~12 bits), from species across the taxonomy
(taxid spread ~17 bits), with flex cells that are 82% distinct. Per-key widths are therefore nearly
the global ones (about 37 bits + 2 flags against 42), the dictionary does not pay (76.4% against 75.8%
without it), and delta coding cannot touch the flex cells (73.6%). What remains is the 22 bits of
padding in every 64-bit entry: **6-byte entries 82.9%, 42-bit entries 76%** of today's bytes, −6.0 and
−8.2 GB. The extrapolation this section first made from the 900-species worlds (45-65%) was wrong for
this reason.

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

## 6. Measuring r226 (done: section 9)

`scripts/index_layout_db` reads the index out of `database.protal` chunk by chunk and prints section 3
and 4 for the real database; nothing is held in memory but a few chunks per thread (~75 MB each). On a
cluster node with the r226 database (18 GB to read, a few minutes with 16 threads; the result for the
v7 database is in section 9):

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
| 1 | **42-bit entries in memory** (implemented, section 10): entries at a fixed bit stride of taxid + gene + position + 2 flag bits, the widths from the index's field maxima (in the header, or found while loading), flex cells a separate 32-bit array; built by `LoadColumns` from the decoded columns; file format unchanged | **−8.2 GB** (values 35.0 → 26.8; section 9) | none measured per lookup (section 5: an entry is a shift and a mask at its bit offset; the flex scan is unchanged) | M-L | the addressing: the key map's per-key 16-bit offsets count 8-byte cells today; with entries and flex cells in two arrays the key's entry index and flex index both follow from the cell offsets (entries e of a key of S cells: e = S − ⌈S/3⌉, as `Get` computes now), so the key map can stay as it is and only the two array indices are derived; `Get`, `GetFromLookup`, `GetExact`, `GetSingleEntry`, `RecoverFromLookup`, `PrefetchValues` adapt; `--build` keeps its 8-byte layout (it runs its own uniqueness check on it and writes the file), so a query run has one layout and a build another; verify with byte-identical SAM and profiles on the 900-species worlds and the e2e tests, speed on `dbdense` and r226 |
| 2 | **6-byte entries** instead of 42-bit ones (byte-aligned, an unaligned 8-byte load and a mask) | −6.0 GB | none | M | the same two-array rework with simpler indexing; 1 costs little more and saves 2.2 GB more |
| 3 | Per-key minima and widths (section 4's bit-packed layout) or a flex dictionary on top of 1 | −0.3 GB, or nothing | per-key headers on the lookup path | L | measured at r226: not worth it (section 9) |
| 3b | Collapse **ubiquitous flex groups** at load (section 6.3) | −0.25 GB (0.73% of the entries) | | S-M | not worth it on its own |
| 4 | `posix_fadvise(POSIX_FADV_DONTNEED)` on the database file after each member is loaded (`Zstd.h` only advises `SEQUENTIAL`) | up to 20 GB of page cache | none | S | matters for cgroup-limited jobs and nodes shared with other jobs; not in RSS |
| 5 | Cap the index load's per-worker frame buffers (`ForEachFrame` holds a chunk's compressed input and decoded output per worker, ~75 MB; 64 threads ≈ 5 GB on top of the index while loading) | 0.075 GB per thread, at the load peak only | the load uses fewer workers than threads, or frees buffers as it goes | S | or build with `--compress_frame_mb 16` |
| 6 | Cohort (`--map`) retention per sample after profiling (~0.15 GB per dense sample for the strain MSAs, [earlier report](../2026-09-30-memory-profiling/README.md) §1), unchanged since | 0.15 GB × samples | | S-M | only for runs of hundreds of samples |
| 7 | Several protal processes on one node each load the index: run the node's samples as one `--map` run instead (one index, one sample on all threads) | N × 43 → 43 GB | | none (usage) | a shared `mmap` of a raw index would do the same with 4 KB pages, at the TLB cost the THP work removed |

Together 1 + 4 bring a full r226 run from about 43 to about **35 GB** of anonymous memory, with the
database's page cache released as it is read. Below that the representation itself would have to
change (fewer values, or values that do not carry the species: both change results).

### Not worth it

- **Key map** (3.2 GB, 7%): 9.2% of the 2^30 keys are non-empty at r226, but the map's granularity is
  the block of 8 keys, 54% of which are non-empty at that fill, so a sparse map saves ~1.5 GB at the
  price of a dependent load (a bitmap and its rank) on every lookup, which is what the THP and prefetch
  work spent its effort on. Only if memory, not time, is the hard limit.
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

## 9. The r226 database, measured

`index_layout_db` on the r226 v7 database (`protal0.7.3_r226_v7/protal_db/database.protal`, built by
the [v7/v8 run](../2026-10-03-r226-v7-v8-evaluation/README.md); 700 chunks, run on the HPC login
node, 2026-10-03):

| | r226 v7 |
|---|---:|
| non-empty keys | 98.6M (9.2% of 2^30) |
| value slots | 4,375,430,329 = **35.0 GB** |
| entries | 2,903,557,366 (29.4 per non-empty key) |
| single-value keys | 5.2M (0.2% of entries) |
| entries per multi-value key | 31.0 (45.4M keys pad half a slot) |
| entries in keys of 33+ / 129+ / 513+ values | 71% / 36% / 11% |
| field maxima (bits) | taxid 143,614 (18), gene 168 (8), position 12,883 (14); 2 flags = **42 bits** |
| distinct flex cells among the multi-value keys' entries | **81.6%** |
| ubiquitous groups (> 256 equal cells in a key) | 46,763 groups, 21.2M entries (0.73%, 254 MB); none > 2048 |

| Layout | GB | B/entry | of now |
|---|---:|---:|---:|
| now (8 B slots) | 35.00 | 12.06 | 100% |
| 6-byte entries + 4 B flex | 29.01 | 9.99 | 82.9% |
| 42-bit entries (global widths, fixed stride) + 4 B flex, computed: 2.90G × (42/8 + 4) B | 26.8 | 9.25 | 76.6% |
| bit-packed per key, flex 32 bits | 26.54 | 9.14 | 75.8% |
| bit-packed per key, flex dictionary | 26.73 | 9.20 | 76.4% |
| delta varint | 25.76 | 8.87 | 73.6% |
| 5-byte entries (does not hold 42 bits; for reference) | 26.11 | 8.99 | 74.6% |

Reading: the flex cells are 4 × 2.90G = 11.6 GB and stay 4 bytes in every layout (distinct within
their key, so neither a dictionary nor XOR deltas shrink them); the entries' 42 bits are 15.2 GB
against 23.2 GB in 8-byte slots. Per-key widths gain 0.3 GB over the global ones, because a key's
values span different genes and positions and species across the taxonomy. The entries-per-key
distribution (71% of the entries in keys of 33 or more values) is also why delta coding's 2-3× on long
blocks (section 5) would be the common case here.

What this changes against the extrapolation in sections 4 and 7 as first written: the saving is
8.2 GB, not 12-17; the dictionary and the per-key headers are not worth having; the ubiquitous
groups are not either (0.73%). The 6-byte and 42-bit layouts, which keep random access and leave the
flex scan as it is, are the whole of what is left.

## 10. Implemented: the packed index (2026-10-03)

Option 1 of section 7, in the working tree on top of `6b432dd` and committed with this section.

| | |
|---|---|
| `Seedmap::PackedLayout` (`Seedmap.h`) | the widths: `taxid_bits + gene_bits + pos_bits + 2` per entry (`EntryBits`), and `SlotBits() = max(W, ceil(2 (32 + W) / 3))` bits of region per slot of the file, the least that holds every key (a key of S slots has e = S − ⌈S/3⌉ ≤ 2S/3 entries and, from S ≥ 2, e flex cells of 32 bits; S ≡ 1 mod 3 only for S = 1); 50 bits per slot at r226, 43 on the 900-species worlds |
| `Seedmap::Pack`, `PackBlocks`, `PackKey`, `PutBits` | the 8-byte layout converted in parallel by block ranges; a key's region starts at bit (first slot) × SlotBits, its flex cells first, then its entries; the bytes a range shares with its neighbours (regions end mid-byte) are OR'd atomically, the rest with plain 8-byte read-modify-writes into zeroed memory |
| `Seedmap::LoadColumns` with a layout | each chunk of the column format is decoded into the thread's own buffer and packed from there, so the 8-byte layout is never in memory; a raw `index.prx`, a seekable raw `.zst` or a single-frame `.zst` is read as before and then packed (both layouts for a moment) |
| `Seedmap::GetPacked`, `PackedBlock`, `FlexCell`, `EntryValue` | a key's values: where its entries and cells lie (byte address and bit shift) and how many; a cell is one 8-byte load and a shift, an entry a 16-byte load, a shift and a mask, then the fields moved into a `ValueEntry`'s 20/20/20 bits so the rest of the code is unchanged |
| `KmerLookupSM` (`KmerLookup.h`) | `LookupPointer` derives from `PackedBlock`; `Get`, `GetFromLookup`, `RecoverFromLookup`, `PrefetchValues` read through it; `GetSingleEntry` and `GetExact` (the build's uniqueness check) keep the 8-byte `Get`, which now stops if the index is packed; `GetFromLoopupSIMD` (unused) removed |
| `GenomeLoader::IndexFieldMaxima`, `RunProtal.h` | a query run takes the widths from the reference's largest taxid, gene id and gene length (+ flex_k for the core's offset) and loads the index packed; it prints `Index in memory: <entries> entries of <W> bits (...) and their 32-bit flex cells, <B> bits per slot of the file: <GB>; key map <GB>`. A value outside the widths stops the run: "holds a value outside the reference's ranges ...: it was built against a different reference" (the fingerprint check comes after the load) |
| `--build`, `--compress_db`, `--decompress_db`, tests of the raw format | unchanged: they load without a layout and keep the 8-byte layout; the files are unchanged |
| `tests/test_PackedIndex.cpp` | an index built as `--build` builds it (count, lay out, place; 3,000 cores, 1-300 values each, flags set), packed in 4 threads: every flex cell and entry equal to the 8-byte layout's, and `KmerLookupSM::Get` returns what the stored cells say; the same through `SaveCompressed` (1 MB chunks) and a packed load; a value outside the layout exits 8 with the message; `PutBits` on bytes shared between ranges; the widths and the fit of SlotBits for every S to 4000 |
| docs | `database-files.md` (The index in memory), `running.md` (Memory: ~35 GB at r226) |

### Verification

Built in WSL from Linux-side copies (`scripts/build_packed.sh`: `HEAD` as the baseline, `HEAD` with
this change's files as the new tree; Release, Ninja), `scripts/verify_packed.sh`, `scripts/ab_packed.sh`,
`scripts/ab_followup.sh`:

| | |
|---|---|
| Unit tests | 334 of 334 (6 new) |
| End to end | 131 of 131 on a fresh mini database (`CompressedDatabaseTest` first showed that a zstd read error of the packed load lacked the "(truncated or corrupt file?)" suffix; fixed) |
| Same outputs, 1 thread | `db900n`, `mix`: all 168 output files byte-identical (the timing tables left out) |
| Same outputs, 6 threads | `db900n` `mix` (176 files) and `w900` (373), `dbdense` `dense_mix` (132) and `dense_w` (240): every file identical, SAM records compared sorted (the thread order); in `mix` one float's last digit (0.987188 / 0.987187 in `*.profile.gene.log` and `*.profile.genes.log`), a sum in another thread order, as the gene-store report saw |
| Memory, `db900n` `mix`, 6 threads (`mem_trace.sh`) | RSS while aligning 3528 → **3428 MB**; while loading 3241 → 2761 MB (the 8-byte values are not allocated; the chunk buffers are freed after the load). The values: 213 MB → 143 MB (26.65M slots × 43 bits), the key map 3.2 GB as before |
| Time, 6 threads, 2 rounds | `w900`: load 0.3-0.5 s both, align 3.4 s both, run 5.1-5.3 s both; `dense_w`: align 6.3-6.5 s baseline, 6.4 s new (one run 8.9 s: the laptop); `mix`: align 2.0 / 2.1 s. Within the noise of this machine, as the microbenchmark predicted |

**At r226** (from section 9's counts): 4,375,430,329 slots × 50 bits = **27.3 GB** of values instead of
35.0, so a run of about 35 GB instead of 43 (key map 3.2, genes 4.3, tables < 1). The chunk buffers
of the load add ~50 MB per thread while the index loads. Not measured at that size: the `Index in
memory:` line of the next r226 run gives the exact figure, and its alignment time against the previous
run the cost, expected nil.
