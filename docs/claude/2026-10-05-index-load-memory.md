# The index load's tail: workers stop one by one, their buffers unmapped

2026-10-05, branch `audit-fixes` at `2d809cf` plus this change. Follows "The index load's buffers" in
[2026-10-04-performance-gtdb-scale](2026-10-04-performance-gtdb-scale/README.md) (`3565360`), which left ~2.3 GB of
frame buffers on 32 loading threads at r226 (per thread the compressed frame, ~25 MB, and the decompressed chunk,
~46 MB). The user asked whether the loading threads could be shut down one by one as the load nears the database's
theoretical resident memory, then to implement that.

## Why it works

The index's key map and packed values are `calloc`'d (`Seedmap::AllocateKeymap`, `AllocatePacked`): their pages become
resident as the loading threads write them, so the index's memory climbs to its full size during the load, with the
threads' buffers on top. The peak is at the end, and it was the full size plus every thread's buffers:
`zstd::ForEachFrame` kept all workers' buffers (`std::vector<State> states`) until the last worker had finished, and
glibc then kept the freed blocks under its mmap threshold (the ~25 MB compressed frames) resident in the threads'
arenas for the rest of the run.

## What changed

- **A worker's buffers are its own and go when it stops** (`ForEachFrame`'s worker-local `FrameBuffer`s), allocated by
  `zstd::MappedAllocator`: anonymous mappings, unmapped when freed, so the memory goes back to the system at once.
- **Near the end the workers stop one by one** (`zstd::LoadBudget`): a worker takes another frame only while the
  output of the frames not yet started (what they will write into the index, known from the container: each chunk's
  key map cells, 2 bytes each, and its values, 64 bits each unpacked or the slot's bits packed;
  `index_codec::ChunkOutput`) is at least the buffers of the workers then decoding, each counted at the largest
  frame's (compressed plus decompressed); otherwise it stops and its buffers are unmapped. What is written and what is
  buffered then stay within the index's full size, plus one worker's buffers: the last worker always goes on, and
  the last frame, after which nothing is left to write, starts once the others are done. Frames already in progress
  are counted as written (their unwritten part is not known), so the rule errs on the safe side.
- Used by both index loads into `calloc`'d memory: the packed load (`Seedmap::LoadColumns`) and the unpacked one
  (`index_codec::Decode`). `Verify` (the build's check) and the other frame readers are unchanged, but for the
  buffers' release.
- The run says when it happened: `Index: its last N chunks decoded on fewer threads, so that their buffers stay within
  the index's memory`, before `Memory after loading the index`.

The decoded index is the same: each chunk writes only its own region, whatever the order and the number of threads.

## Tests and measurements

Built in WSL (`~/zstdreader`, `2d809cf` + this change, Release).

- New `ZstdSeekable.ForEachFrameStopsWorkersNearTheEndOfALoad` (48 frames, 8 threads, frames of a frame's buffer of
  output, of a third, and of an ample 1 TB): every frame handled once with its content; whenever a frame starts, the
  frames running beside it times the largest buffer fit in the output still to come, or it runs alone; the tail runs
  on fewer workers (from 8 or 24 frames before the end); an ample budget throttles only the last frame; without a
  budget nothing changes. Passed 5 times in a row. `ctest` 367 of 367; the e2e suite 133 of 133 (the decoded index and
  all outputs unchanged).
- **Memory**, the 0.7.5 benchmark database (`~/bench071/V075`, 113 MB, 68 frames of up to 67 MB decompressed and 9 MB
  compressed; its key map 3 GB), a run of 2,000 random read pairs (the run's peak is the index load), `Memory after
  loading the index` of 3 alternated runs each:

| Threads | Peak before → after | Resident after the load, before → after | `Load Index` before → after |
|---|---|---|---|
| 6 | 3.33-3.35 → 3.22-3.24 GB | 3.28-3.31 → 3.21-3.22 GB | 404-591 → 416-447 ms (4 chunks on fewer threads) |
| 3 | 3.29-3.31 → 3.23-3.24 GB | 3.26-3.31 → 3.21-3.22 GB | 616-913 → 649-797 ms (2 chunks) |

  The peak drops by the buffers (~0.1 GB here: few chunks, most small), and so does the memory after the load, which
  the old code kept for the rest of the run (freed but resident in glibc's arenas). The load time is within the noise
  of this shared laptop. On the mini database (1 MB, chunks of almost nothing) nothing changes, as expected.

**Expected at r226** (32 loading threads, ~450 chunks): the peak from the index's size plus ~2.3 GB to its size plus
one thread's buffers (~70 MB), about −2.2 GB, and those ~2.3 GB no longer resident after the load; the tail runs on
fewer threads for the last few dozen chunks, a few tenths of a second (an estimate from ~0.2 s per chunk and thread).
The next cluster run's `Memory after loading the index` and the new `Index: its last N chunks ...` line will tell, as
will whether the index load is still the run's peak (`Memory after aligning`).

## What else could go

Asked after this change. A run at r226 then holds the index's values (27.3 GB at 50 bits per slot,
`Index in memory:` of the ninth cluster run), its key map (3.2 GB), the genes (4.1 GB at 2 bits a base,
`Preload:`) and ~1-2 GB of tables, threads and buffers: about 36 GB, nearly all of it the database's data.
The [memory audit](2026-10-03-memory-audit/README.md) (sections 7, 9) already ruled out per-key widths, a
flex dictionary, delta coding and genes on demand. What is left, by what it saves for what it costs:

| Option | Saves at r226 | Cost | Outputs |
|---|---|---|---|
| **Bits per slot as a fraction.** `PackedLayout::SlotBits()` is ⌈2(32+W)/3⌉ = 50 so that every key fits its S slots (S = 1: W bits; S ≥ 2: e = S − ⌈S/3⌉ entries of 32 + W bits, at most 2S/3 of them). The rounding is needed per key, not per slot: a fixed-point width of 1579/32 = 49.34 bits (≥ 148/3), a key's region at bit ⌊slot × 1579 / 32⌋, holds every S too | −0.36 GB (4.38G slots × 0.66 bit) | a multiply and shift per address; small | the same |
| **taxid and gene as one number**, taxid × genes + gene: 143,615 × 169 < 2^25, so W = 41 instead of 18 + 8 + 14 + 2 = 42, and 1558/32 = 48.69 bits per slot | another −0.36 GB | a multiply by a reciprocal per decoded entry; every reader of an entry adapts; medium | the same |
| **Runs from SAM files load only the genes their profiling reads**: the congener sketches (`SampleContext::MakeSketch`) and the strain MSAs' references, of the taxa the SAM headers name, instead of all 14.5M genes | −4 GB in those runs (the cluster's cohort run from two SAMs peaked at 7.7 GB) | a gene list before the preload; medium | the same |
| **Cohort retention**: ~0.15 GB per dense sample kept for the strain MSAs ([2026-09-30](2026-09-30-memory-profiling/README.md), 900 species; never measured at r226); drop the alleles under `--snp_min_cov`, or spill a sample's strain data to disk | 0.15 GB × samples, in `--map` runs | measure a large cohort first | the same if exact |
| **Each flex cell once per run of equal ones in a key**: 81.6% of them are distinct in their key, so 18.4% of the 11.6 GB of flex cells go, less a flag bit per entry; the entries sorted by flex | −1.8 GB | a key's flex cells no longer follow from its slots: the key map's words need entry and flex bases, and a flag bitmap's popcount tells a key's first cell, one more dependent load per lookup (the objection that ruled out the sparse key map); large | to verify (the order of entries in a key changes; the seeds are sorted on all their fields) |
| Sparse key map (audit) | −1.5 GB | one more dependent load per lookup | the same |
| Fewer syncmers | in proportion | different seeds, features and models | different |
| `posix_fadvise(POSIX_FADV_DONTNEED)` on `database.protal` once loaded (audit item 4) | up to 27 GB of page cache, not resident memory | every load cold (the cluster: 2.8 s warm, 539 s cold from NFS), so only as an option | the same |

The first two together take ~0.7 GB without changing a byte of output or the file format; the third
matters only for runs from SAM files; the flex runs and the sparse key map trade lookup speed for
memory and are worth it only where memory, not time, is the limit. The next cluster run's `Memory after
aligning` and `Memory after profiling` lines will tell whether any stage after the load comes near it.

## Implemented: bits per slot as a fraction, taxid and gene as one number

On `60fea63` (the load change above, committed), at the user's request.

- `Seedmap::PackedLayout` holds `genes` (the largest gene id + 1, at least 2), `taxon_gene_bits` (the bits of
  taxid × genes + gene for the largest of each) and `pos_bits`; an entry is that number, the position and the 2
  flags. `SlotBits32()` is the slot's width in 32nds of a bit, ⌈64 (32 + W) / 3⌉ (at least 32 W), and a key's
  region starts at bit `SlotBit(slot)` = ⌊slot × SlotBits32() / 32⌋: a region of S slots spans at least
  ⌊S × SlotBits32() / 32⌋ whole bits wherever it starts, enough for the e ≤ 2S/3 entries and flex cells of any S
  (the test checks every S below 4000 from all 32 starting phases, and that one 32nd less fails).
- `EntryValue` splits the number by `DivideByGenes`: the high 64 bits of taxid_gene × ⌈2^64 / genes⌉, exact for
  numbers below 2^40 and genes up to 2^20 (the excess is below 2^-24 < 1 / genes; the test checks the first and
  last 2^20 numbers and 10^6 random ones for 2, 256, 169, 1000 and 2^20 genes). `PackValue` refuses a gene id at
  or above `genes`, or a number wider than the layout, with the same message as before.
- `index_codec::ChunkOutput` takes the width in 32nds; the files, the key map and every reader of an entry
  (through `EntryValue`) are unchanged; `Index in memory:` now reads e.g. `entries of 41 bits (taxid * 169 + gene
  25, position 14, 2 flags) and their 32-bit flex cells, 48.69 bits per slot of the file`.

**At r226** (4,375,430,329 slots, widths of the ninth cluster run): entries 42 → 41 bits, 50 → 48.69 bits per slot,
values 27.35 → 26.63 GB, **−0.72 GB** (half from the fraction, half from the joined number). On the local 0.7.5
world (taxid 10 + gene 8 bits → 17) entries 32 → 31 bits and 43 → 42.00 bits per slot, ~2% of its 0.1 GB.

**Checks** (WSL, `~/zstdreader`, Release): `ctest` 368 of 368 (two tests new or rewritten: the widths and slot
phases, the joined number for six gene counts); the e2e suite 133 of 133. Outputs against `60fea63`'s build on the
0.7.5 world, a paired-end sample of 214k pairs (`~/bench071/strains/reads/str_s1`) and a PacBio sample of 215
reads (`samples_lr073/points/pb_b3000000/sim/reads/pb_b3000000_s_1`): at 1 thread every file identical, SAMs
included; at 4 threads all identical but the paired-end SAM's record order, which differs between two runs of the
old build too (the records sorted, and the headers, are identical). Cost: callgrind on 20,000 pairs at 1 thread,
`RunPairedEnd` 3,348.95 → 3,344.61 M instructions and `ChainAnchorFinder` (the lookups) 629.6 → 627.5 M: the
multiplication costs nothing measurable. Alternated timings (6 each, 2 threads, `Aligning reads`) say only that
this machine is noisy: a first version, which divided with a variable 128-bit shift, looked ~5% slower (medians
3.53 vs 3.31 s, load 7), the final one ~10% faster (2.06 vs 2.29 s, load 3).

## Reproducing

In WSL, with `protal_before` a build of `2d809cf` and `protal_after` one of this change, and `m_R1.fq`, `m_R2.fq` 2,000
random 150 bp pairs:

```
for t in 6 3; do for rep in 1 2 3; do for bin in before after; do
  ./protal_$bin --db ~/bench071/V075/protal_db -1 m_R1.fq -2 m_R2.fq --prefix m -o out_${bin}_t${t}_$rep -t $t \
      --no_qcmsa > out_${bin}_t${t}_$rep.log 2>&1
  grep -E "Memory after loading the index|Load Index took|^Index: its last" out_${bin}_t${t}_$rep.log
done; done; done
```
