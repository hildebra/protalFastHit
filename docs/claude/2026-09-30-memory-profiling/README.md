# Where protal's memory goes, and how to reduce it

2026-09-30. Branch `audit-fixes`, commit `4b21427` (Version 0.7.0), built in WSL from a copy on the
Linux file system (Release, `-g -fno-omit-frame-pointer`, baseline ISA). Nothing in the checkout
changed; everything here is measurement and proposals.

**Question**: profile protal's memory use and say whether a viable strategy exists to reduce it.

**Answer in short.** A query run's memory is the database and nothing else: at GTDB r226 size
(extrapolated, see the assumptions) **~59 GB = 35 GB index values + 17 GB reference genes + 3.2 GB
key map + 2 GB gene tables**. Reads, threads and depth do not matter (≤ 20 MB per thread, flat from
1M to 8M read pairs). The profiling stage takes ~0.1 GB per 1M alignments, and only a cohort run
retains that per sample. Three changes are viable and cut the 59 GB to about **38 GB (−36%)**
without changing any result: store the genes 2-bit packed (−13 GB), hold the index's multi-value
blocks delta-coded in memory (−7 GB or more), and slim the `Gene` records (−1.3 GB). Loading genes
on demand, which is the obvious idea, does **not** work: alignment touches 96% of the genes' pages.
Details, evidence and what each costs are below.

## Setup

| | |
|---|---|
| Machine | WSL Ubuntu 24.04, Intel Core Ultra 7 258V (4 fast + 4 low-power cores), 24 GB, THP `madvise` |
| Databases | `db900` (900-species GTDB-like world, 101k genes, 26.9M index slots, `database.protal` 125 MB) from `~/protal-perf`, the one the performance reports use. No real GTDB database here. |
| Reads | `mix` (1M pairs, ~5% from `db900` species, rest random sequence), `w900` (1M pairs from 60 species, 895k alignment records), `w900x8` (8M pairs), `cohort N` (N copies of `w900` as N samples). Simulated by [the performance report's scripts](../2026-09-29-performance-profiling/README.md). |
| Method | `/proc/PID/status` (`VmRSS`, `RssAnon`, `VmHWM`) every 0.1 s, stages cut at protal's own timing lines (`scripts/mem_trace.sh`, `stage_peaks.sh`); heap attribution with `valgrind --tool=massif` (`massif.sh`, `massif2.sh`); the index's layout read from a raw `index.prx` (`index_stats.cpp`, `block_delta.cpp`); a log of the gene pages alignment reads (`touch_patch.pl`, a patch for a copy of the sources); a window-decode benchmark (`bench_unpack.cpp`). |
| Caveat | Timings on this laptop swing with power state; memory does not, so the memory numbers are firm. The GTDB-size numbers are extrapolations from `db900`'s per-gene and per-slot costs, as in [the load report](../2026-09-30-load-and-output/README.md). |

## 1. Where the memory goes in a run

`db900`, `mix`, 8 threads, total run 8.6 s:

| Stage | RSS |
|---|---:|
| start, options, genes preloaded | 0.16 GB |
| model (cPMML) loaded | 0.22 GB |
| index loaded | 3.66 GB |
| aligning (plateau) | 3.73–3.75 GB (VmHWM 3.75 GB) |
| after the index is freed, profiling | 0.36–0.6 GB |

The plateau is reached at the end of loading; aligning adds 10–90 MB. Of the 3.73 GB, 2.28 GB are
in transparent huge pages (`AnonHugePages`), no bloat.

**Heap at the peak** (massif, 1 thread, 100k pairs; total heap 3.634 GB):

| Allocation | Bytes | Share | Scales with |
|---|---:|---:|---|
| `Seedmap::AllocateKeymap` | 3,221,225,480 | 88.6% | nothing: 2^30 keys × 2 B + control blocks |
| `Seedmap::AllocateValues` | 214,905,912 | 5.9% | 26.86M slots × 8 B |
| gene arena (`LoadAllGenomes`) | 110,881,098 | 3.05% | bases of `reference.fna`, 1 B each |
| the rest (zstd buffers, tables, taxonomy, model) | ~46 MB | 1.3% | |

### Threads, depth, cohort (`w900`, `db900`)

| | Peak RSS |
|---|---:|
| 1 / 2 / 4 / 8 / 16 / 32 threads, 1M pairs | 3.51 / 3.54 / 3.62 / 3.75 / 3.98 / 4.12 GB |
| 8 threads, 8M pairs (7.2M alignment records) | 3.73 GB while aligning |

About **20 MB per thread** (16 MB SAM buffer per output handler plus read batches), and flat in
depth: the read batches are bounded. Profiling the 8M-pair sample peaks at 0.86 GB (1M pairs:
0.36 GB), so about **80 B per alignment record** in the profiler, after the index is freed.

**Cohorts** (`--map` with N samples, each `w900`, 4 threads; RSS after profiling, i.e. what the
strain stage starts from):

| Samples | default | `--no_strains` |
|---:|---:|---:|
| 4 | 0.95 GB | |
| 12 | 2.18 GB | 1.07 GB |

That is **~154 MB per sample** retained for the strain MSAs (~60 MB per sample without strains),
on top of ~0.33 GB. Massif on a small cohort attributes it to the per-position `Variant` records
in `VariantHandler::m_variants` (a `robin_map` of `vector<Variant>`, ~100 B each plus two byte
vectors), to the `ReadInfo` of every read in the merged `SequenceRange`s (24 B per read and gene),
and to the reference alleles `PostProcessSNPBin` adds to each bin. It does not depend on the
database, and at ~400 samples of this kind it reaches the 59 GB of the GTDB run; the index is gone
by then, so it only matters for big `--map` runs. A real sample retains what its
species and depth give; `w900` is dense (60 species, ~25× each).

**Building** (`protal --build`, `db900`, `scripts/build_memory.sh`): peak RSS **3.75 GB** at `-t 1` (2 min 35 s)
and **5.24 GB** at `-t 8` (2 min 18 s), as in the [index-build reports](../2026-09-29-index-build-parallel.md)
(3.7 / 5.0 GB). The key map and the index values are the same 3.4 GB as in a query run; the 1.5 GB more
at 8 threads (~0.2 GB per added thread) are probably the zstd compressors of the packing step (level 19,
128 MB window; not attributed). At GTDB size a build holds the index (39 GB), the reference (17 GB) and
the tables, so ~58 GB plus that per-thread part, as the earlier estimate (50–60 GB) said; not measured
at that size.

## 2. GTDB r226, extrapolated

Assumptions: 16.6M genes (136,646 bacterial + 6,968 archaeal species × 119 / 52 marker genes,
[load report](../2026-09-30-load-and-output/README.md)), ~17 Gbp of reference, and `db900`'s
**266 index slots per gene** (measured: 26.86M slots / 100.8k genes; it counts the genes' syncmers,
so it should not change with database size).

| Part | Formula | GB | Share |
|---|---|---:|---:|
| index values | 16.6M × 266 × 8 B | 35.4 | 60% |
| reference genes | ~17 Gbp × 1 B | 17 | 29% |
| key map | fixed | 3.2 | 5.5% |
| gene tables | 16.6M × 112 B (`sizeof(Gene)`; the 16.6M-gene benchmark measured 1.8 GB) | 1.9 | 3.2% |
| taxonomy, model, buffers, threads | | ~0.3 | 0.5% |
| **total** | | **~58** | matches the 59 GB in [running.md](../../running.md) |

Reading the index file from disk leaves its compressed bytes in the page cache (33 GB for the full
database; `Zstd.h` advises `POSIX_FADV_SEQUENTIAL` but never `DONTNEED`). It is reclaimable, but a
job limited to 60 GB of cgroup memory shares it with the 58 GB of anonymous memory.

## 3. What the index holds

`index_stats` on `db900`'s raw index (1.07G keys, 26.86M slots):

| | keys | values' slots |
|---|---:|---:|
| non-empty keys | 11.3M (1.1% of the key space) | |
| keys with one value | 77.9% | 32.8% of the slots |
| keys with 2+ values (flex blocks) | 22.1% | 67.2%, of which 23.3% of all slots are flex cells |

A key with n ≥ 2 values stores n 4-byte flex cells in ⌈n/2⌉ 8-byte slots in front of its values, so
a value in a multi-value block costs **12 B** (8 B value + 4 B flex), one of a single-value key 8 B.
The value is 62 bits (20 taxid, 20 gene, 20 position, 2 flags) of which GTDB needs about
18 + 8 + 14 + 2 = 42. At GTDB size there are ~3 values per key of the whole key space on average
(3.4G values for 1.07G keys, more where the genes are conserved), so I expect most values to sit in
blocks of many related species' copies of a k-mer, which differ in the taxid by little and in the
position by a few bases; `db900` shows the direction but not the size.

`block_delta` codes each multi-value block as differences (taxid, gene, position as varints; flags
in the gene byte) and measures 4.47 B per value against 8 B, in the stored order (which is already
ordered by taxid, sorting does not help): 52.6 MB for 11.8M values instead of 94.2 MB. Applied to
`db900` that is the index 215 MB → 173 MB (−19%); `db900`'s blocks are short (4.7 values), GTDB's
will be much longer and their differences smaller, so I expect a larger saving there; but that is
an expectation, from a synthetic world whose relatedness structure is only GTDB-like.

## 4. Do the genes need to be resident?

The genes are 29% of the memory, and a sample reads a few species' genes, so loading them on demand
(`--preload_genomes_off`, or an `mmap` of a raw `reference.fna`) looks attractive. I logged which
4 KB pages of the gene arena `Gene::Sequence()` hands out during alignment (`touch_patch.pl`):

| Sample | `Sequence()` calls | Pages touched by alignment | per read pair |
|---|---:|---:|---:|
| `mix` (5% from species) | 2.26M | 26,225 of 27,071 (**96.9%**) | 2.3 |
| `w900` (60 species) | 24.5M | 25,924 of 27,071 (**95.8%**) | 24 |

(Profiling alone, from an existing SAM, touches 23–30% of the pages; it adds nothing to the union.)
The calls come from the chain finder verifying seed hits, which spurious 15-mer hits of any read
cause as well as true ones, so they land on genes all over the reference. The pages are counted per
call over the whole gene (genes are ~1 kb, a page holds ~4), which slightly overstates the touched
set. With 160× more genes in the index than `db900` has, spurious hits land in more places, not
fewer, so a page-granular on-demand store would end up reading almost the whole file in any real
sample, and pay page faults for it. Loading on demand saves memory only for tiny samples.
**Not viable.**

The genes are four letters: 2 bits per base are enough. Reading a window from a 2-bit store
(`bench_unpack.cpp`, 1 Gbase store, random windows, scalar 256-entry table):

| Window | bytes arena (today) | 2-bit store | added |
|---:|---:|---:|---:|
| 64 | 18 ns | 55 ns | 37 ns |
| 150 | 31 ns | 146 ns | 115 ns |
| 350 | 55 ns | 174 ns | 119 ns |
| 1000 | 155 ns | 218 ns | 63 ns |

At `w900`'s 24 calls per pair that is 1.5–3 µs per pair against ~64 µs of alignment CPU per pair
(8 s × 8 threads for 1M pairs): **2–5%**, and below 0.5% for `mix`. The numbers are a
micro-benchmark with independent accesses; SIMD decoding would lower them.

## 5. Strategies

Saving at GTDB r226 size, from the 58–59 GB; results are unchanged unless noted.

| # | Change | Saves | Cost | Effort | Risk |
|---|---|---:|---|---|---|
| 1 | **Genes 2-bit packed** in the arena (`Gene::Sequence()` returns a decoded window, not a `string_view` of the whole gene); non-ACGT bases (N, IUPAC) in a per-gene exception list | **13 GB (22%)** | 2–5% of alignment CPU at 24 calls per pair | M: ~40 call sites read `Sequence()` (alignment, chain finder, SAM output, profiler, MSA); most need a window, `WFA2` needs a contiguous reference | how many non-ACGT bases r226's `reference.fna` has (none in the synthetic worlds); check on the real file |
| 2 | **Multi-value blocks delta-coded in memory**, built when the index is loaded from the column format; the file format does not change | **~6.7 GB (19% of the index, at `db900`'s block lengths)**, probably more | one block decode per lookup that has a flex match (blocks are scanned for flex already) | L: `Seedmap::Get`, the lookup result path, the uniqueness flags `--build` clears in place (the build keeps the 8-byte layout), the codec's decoder | lookup speed on long blocks: needs a prototype on a dense index |
| 2b | alternative to 2: 6-byte values (18 + 8 + 18 + 2 bits, or widths taken from the header) | 6.8 GB | unaligned 8-byte loads | L: the 16-bit cells count 8-byte slots, which a 6-byte value breaks | group size limit of 64k slots |
| 3 | **`Gene` diet**: 112 B → 32 B (`uint32` length, id and the four unique counts, offset into the arena; the `std::string`, the file pointer and the start byte go) | **1.3 GB (2%)** | none | S | none |
| 4 | Cohort runs: after `PostProcess`, drop the alleles with fewer than `--snp_min_cov` observations (singleton errors) and `shrink_to_fit`; or write the per-sample strain data to disk and read it back per species | 0.1–0.15 GB per sample in a `--map` run | | S (compaction; needs byte-identical MSAs) to M (spill) | alleles below the coverage floor can never be called; the frequency filter works on subsets, so only the count floor is safe |
| 5 | `posix_fadvise(DONTNEED)` on the database after loading | up to 33 GB of page cache | | S | none |
| 6 | Several protal processes per node (array jobs): one index in shared memory (`/dev/shm` with `huge=within_size`, or `mmap` of a raw `index.prx`) | N × 39 GB → 39 GB | TLB misses without huge pages cost ~⅓ of seeding time ([performance report](../2026-09-29-performance-profiling/README.md)) | M–L | the raw index is 39 GB on disk (3.2 GB key map + 35 GB values) |

Together 1 + 2 + 3 take the 59 GB to **~38 GB**.

### Not recommended

- **Genes on demand** (section 4): 96% of the pages are read anyway.
- **Lower syncmer density** (fewer values, linear saving): changes seeding sensitivity and the
  model's features; needs an accuracy study and a retrained model.
- **Clade shards loaded after a screening pass**: changes which species a read can hit, hence the
  statistics; large.
- **A sparse key map**: the 3.2 GB are 5% at GTDB size; it would help only databases of a few GB
  (`db900` uses 1.1% of its key space), at the price of a slower lookup.
- **Fewer or smaller per-thread buffers**: 20 MB per thread.

## Reproducing

```bash
# build a copy on the Linux file system (see the performance report), then
OUT=~/protal-mem/runs BIN=~/protal-mem/build/protal bash scripts/mem_trace.sh NAME DB READS_DIR THREADS
bash scripts/stage_peaks.sh ~/protal-mem/runs/NAME                   # peak / last RSS per stage
MAP=cohort.map bash scripts/mem_trace.sh ...                          # a cohort (scripts/cohort.sh N THREADS)
bash scripts/thread_scaling.sh ; bash scripts/deep_sample.sh          # 1-32 threads; 8M pairs
bash scripts/massif.sh ; bash scripts/massif2.sh                      # heap attribution, single sample / cohort of 4
protal --decompress_db --db COPY_OF_DB -t 8                           # raw index.prx for the two tools below
g++ -O2 -o index_stats scripts/index_stats.cpp && ./index_stats COPY_OF_DB/index.prx
g++ -O2 -o block_delta scripts/block_delta.cpp && ./block_delta COPY_OF_DB/index.prx
g++ -O2 -o bench_unpack scripts/bench_unpack.cpp && ./bench_unpack
perl scripts/touch_patch.pl SOURCE_COPY                               # then PROTAL_TOUCH=1 on a build of the copy
```

The scripts assume `~/protal-mem` and `~/protal-perf` (the databases and reads of the performance
report); edit the paths at the top.
