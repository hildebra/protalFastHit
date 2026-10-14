# A compact key map: how far the 15-base core table can shrink, whole and in sparse grid runs

Date: 2026-10-13. Branch `audit-fixes` at `19ccc0e` plus the uncommitted whole-index change of the
[gene-grid report](2026-10-13-gene-grid-index.md) (follow-up 4) and a count of cores in the sparse load's log line
(`Seedmap::LoadColumnsAtLevel`: "on K of its N cores with entries"). An evaluation; the follow-up (2026-10-14) builds it
as `--compact_keys` and measures it.

## The questions

1. In a sparse grid run (say C 16), is the key map used, or is part of it wasted?
2. How could protal shrink the key map, in sparse grid runs in particular?

## How the key map works now

The index has two parts. The **values** are each core's entries (taxon and gene, position, flags, the grid level),
packed, each with a 32-bit flex cell of its flank bases. The **key map** says where a core's values are. It is a
directly addressed array over every possible 15-base core: 4^15 = 2^30 slots, each a 16-bit offset, with an 8-byte
control cell (the first value's slot) per 8 cores. That is 3 bytes a core, **3.2 GB whatever the reference holds**.
A lookup reads one key map block (one cache line), then the core's values.

A sparse load (`--sparse_index auto`) keeps the values of the gate's grid level alone, entries and flex cells both
(about 1/C of them). The key map keeps its size; a core whose entries were all left out gets an empty range.

## 1. Is the key map wasted in a sparse grid run?

Not in use, but in memory:

- **In use, nothing is lost.** Every core of every read is looked up, as a read cannot tell which of its cores are on
  the grid. A core without kept entries finds an empty range, a miss like a core the reference lacks. Once a read is
  anchored, the whole read is aligned against the gene's full sequence (the genes are always loaded whole), so its
  alignment and the calls do not depend on which values were kept. The left-out values cost only seeds, which the
  single-seed anchors make up for (gene-grid report, follow-ups 3 and 5).
- **In memory, most of it holds nothing, and not because of the gate.** At r226 the whole index has 2.90 G entries on
  **98.7 M cores, 9.2% of the 2^30 slots** (`local/v18/index_and_package.log`, "Index cores"). The other 91% are empty
  in every run. A grid gate empties more of them, but fewer than one might think: a core is kept if any of its entries
  is on a grid, and at r226 most cores have many entries.

### How many cores a gate keeps

Measured with the new count (one deep sample of the strain world, `protal --sparse_first C --sparse_first_pick grid
--no_profile`, the database's index line):

| database | cores with entries | C 2 | C 4 | C 8 | C 16 | C 32 |
|---|---|---|---|---|---|---|
| strain world (`~/strainfn/db_grid`, 10.3 M entries) | 5.55 M | 3.58 M (64%) | 2.06 M (37%) | 1.13 M (20%) | 0.61 M (11%) | 0.32 M (5.8%) |
| mini database (71 k entries) | 64,165 | 54% | 28% | 14% | 7.2% | 3.7% |

A core of e entries keeps one with probability 1 − (1 − 1/C)^e if its entries are on their genes' grids
independently. With the strain world's cores by entries (1: 3.82 M; 2-15: 1.70 M, mean 3.4; 16-255: 30 k, mean 21) that
gives 2.05 M at C 4 and 0.60 M at C 16: the measured 2.06 M and 0.61 M. Applied to r226's cores by entries (1: 5.2 M;
2-15: 52.6 M, mean 7.2; 16-255: 39.6 M, mean 47; 256-2047: 1.27 M) it estimates:

| r226 | whole | C 4 | C 8 | C 16 | C 32 |
|---|---|---|---|---|---|
| cores with entries (estimate) | 98.7 M | ~88 M | ~74 M | ~59 M | ~43 M |

This is an upper bound: within a band the spread of e lowers it (as the band means overstate small e), and the grid's
ubiquity penalty makes a many-entry core lose more often than chance. A sparse grid run at r226 now logs the real
number on its index line.

## 2. How to shrink it

### A rank-indexed key map (recommended)

Replace the direct array by three parts, built at load:

- a **presence bit per possible core** (2^30 bits, 134 MB), with a cumulative count every 512 bits (8 MB): a core's
  rank among the cores with values in one cache line (presence and count interleaved, 448 bits and a count per 64-byte
  line: 153 MB);
- the **offsets of the cores with values only**, in today's scheme (a 16-bit offset per core, an 8-byte control cell
  per 8 cores): 3 bytes per used core;
- the values as now.

| r226 (GB) | values | key map now | compact | run's other memory (genes at 2 bits a base, tables, buffers) | run now | run, compact |
|---|---|---|---|---|---|---|
| whole index | 26.6 | 3.2 | 0.45 | ~4.2 | ~34.0 | ~31.2 |
| grid C 4 (pe gate) | ~6.7 | 3.2 | ~0.42 | ~4.2 | ~14.1 | ~11.3 |
| grid C 8 | ~3.4 | 3.2 | ~0.38 | ~4.2 | ~10.8 | ~8.0 |
| grid C 16 (pb gate) | ~1.7 | 3.2 | ~0.33 | ~4.2 | ~9.1 | ~6.2 |
| grid C 32 | ~0.9 | 3.2 | ~0.28 | ~4.2 | ~8.3 | ~5.4 |

The values per gate take the measured entry shares (1/4.0, 1/7.9, 1/15.5, 1/29.8 of the entries on both databases
above); the other memory is r226's expected ~34 GB less the index. The compact map saves ~2.8 GB in every run, the
whole index's too; in a sparse grid run it is what lets the index shrink with C: at C 16 the key map is two thirds of
the index today and a sixth with the compact map. Of the ~2.8 GB, the gates' share is small (0.03-0.17 GB); the
emptiness is the reference's, not the grid's. After it, the run's floor is the genes and tables (~4.2 GB).

**Lookup cost.** A lookup now reads one key map line, then the values; it would read the presence-and-rank line, then
the offsets line, then the values: one more dependent cache miss per core found (a core absent stops at the first
line, as now). At r226 the key map takes 19.5 of a lookup's 175 ns, because the keys are prefetched a pair ahead
([seeding profile](2026-10-12-seeding-profile/README.md)); the flex cells take 105 ns. The extra line can be prefetched
in a second stage, when the pair's lookups start (its rank is known once the first line has arrived), so the expected
cost is a few ns a lookup, a few percent of the seeding at most. The smaller tables also need fewer page-table entries
(153 MB of presence and rank in 77 huge pages against 1,600 for the key map), which may offset it. This needs a
measurement (`seedbench` of the seeding profile, at r226) before it becomes the default.

**What changes in the code.** The build keeps its full key map (it counts into it); the file format stays as it is.
Only the loaded index changes: the whole load decodes each chunk's key map into a per-thread buffer and emits presence
bits and used-core offsets (chunks own disjoint core ranges, so this runs in parallel after a count of each chunk's
used cores); the sparse load emits them from its kept cores directly. The query side (`Seedmap::GetPacked`, `Locate`,
the prefetch functions) reads the compact map; the build and the uniqueness check keep the full one. The unit tests
would compare every lookup of a compact index with the full one's, whole and sparse.

### Alternatives considered

- **A hash table of the used cores** (open addressing, 4-byte core and 4-byte offset at 80% load): ~1.0 GB at r226 and
  a probe sequence per lookup. More memory and no simpler.
- **A shorter key** (14 bases: 2^28 slots, 0.8 GB): the core is the syncmer's 15 bases, so the 15th base would move into
  every value (2 bits × 2.9 G entries, 0.7 GB) and be checked there; ~1.5 GB, and the seeding changes. No.
- **A minimal perfect hash of the used cores**: about as small as the rank map, but it maps an absent core to some
  used one, so each core needs a fingerprint to be checked; more complex for no gain.

## Recommendation

Build the rank-indexed key map for the loaded index, whole and sparse, after `seedbench` shows its lookup cost at r226.
It is worth more than any further change to the grid itself: ~2.8 GB in every run, and with the grid at C 16 an r226
HiFi run in ~6 GB instead of ~34 GB today (or ~9 GB with the grid alone). The models for such runs must still be
trained on grid-picked runs (gene-grid report, follow-up 5).

## Follow-up: built and measured (2026-10-14)

The user asked to build it and measure its costs. Code on `b5a9aa9` plus the uncommitted `--compact_keys`:
`src/Hash/CompactKeys.h` (the structure), `Seedmap` (`Locate`, `PrefetchKey` and a new `PrefetchStart` read it; the
whole and the sparse load build it as they decode; `CompactKeymap` converts a loaded direct key map, for
`--sparse_seeding`, whose statistics read the direct one), `KmerLookupSM::PrefetchStart`, the lookup loops of
`ChainAnchorFinder` (the start of a key fetched half as far ahead as its presence line, and both twice as far ahead as
the direct key map's block), `--compact_keys` (a developer option, off by default). Tests: `tests/test_CompactKeys.cpp`
(the structure over line and group edges; an index's every core located the same after `CompactKeymap`; whole and
sparse loads that build it give every core the same slots and every entry the same value and flex cell as direct
loads) and an e2e test (the same records with and without, whole and on the grid's sparse index).

**As built.** 448 cores per 64-byte line of presence bits with the rank count before the line; a 64-bit base per 64
used cores and a 24-bit offset per used core. The loads decode each chunk's key map into a buffer of their worker's,
mark its used cores and keep their starts (4 bytes each, relative to the chunk), then place the ranks, bases and
offsets once all are known; the direct key map is never allocated.

**Data and commands.** The reduced r226 database of the mock benchmark (WSL `~/Mocks/reduced_db_r226`: 16.4 M entries,
10.25 M cores with values; the direct key map 3.2 GB, the compact one 0.2 GB) and two real mock samples: paired-end
`d6331_hku_SRR23926922` (44.8 M k-mers looked up) and the D6331 HiFi run `SRR13128012` (186 M). Ungated
(`--sparse_first 0`: every core looked up), 4 threads, `--no_qcmsa --no_strains`, runs alternating without and with
`--compact_keys`, 5 times each:

```bash
/usr/bin/time -v taskset -c 0-3 protal --db ~/Mocks/reduced_db_r226 -1 SRR23926922_1.fq.gz -2 SRR23926922_2.fq.gz \
    --read_type pe -t 4 --prefix pe -o run --no_qcmsa --no_strains --sparse_first 0 [--compact_keys]
```

(HiFi: `-1 SRR13128012.fq.gz --read_type pb`). A lookup's cost: the seeding stage's seconds per thread
(`misc/<prefix>_runtime.tsv`, row `Seeding`) times the 4 threads over the k-mers looked up (the log's seeding line);
the records compared sorted, without the SAM header. Three rounds: the conversion after the
load (stage 1), the build during the load (stage 2), and stage 2 with the longer prefetch lead (final).

| final round | without | with `--compact_keys` |
|---|---|---|
| SAM records (sorted, without the header) | | identical in all 10 pairs (every round) |
| peak RSS | 3.34-3.43 GB | **0.74-0.87 GB** |
| index load | 0.66-0.99 s | 0.76-0.87 s (stage 2) |
| paired-end: ns per lookup, median (paired mean difference) | 103.9 | 96.0 (−5.6 ns: noise) |
| HiFi: ns per lookup, median (paired mean difference) | 86.0 | 90.6 (+7.3 ns, +8.6%) |
| HiFi: wall, median | 13.30 s | 13.64 s (+2.6%) |

- **Memory: −2.6 GB a run (−76% here).** The run peaks at 0.8 GB instead of 3.4 GB on this database; the whole r226
  index would drop from 29.8 to ~27 GB (key map 3.2 → ~0.45 GB), a grid run at C 16 from ~4.9 to ~2 GB.
- **Lookups: a few nanoseconds.** Paired-end reads, whose keys are fetched a read ahead, show no difference beyond
  the noise (the medians moved −8% to +7% between rounds). HiFi reads, looked up chunk by chunk with a lead of a few
  k-mers, pay ~7 ns a lookup (+10 ns before the lead was doubled), +2.6% of the wall time. At r226 a lookup takes ~175
  ns (its flex cells most of it), so the same few ns are ~4% of the seeding; the far smaller tables (153 MB of presence
  against 3.2 GB) may need fewer page-table walks there, which this database, whose key map misses the caches either
  way, cannot show.
- **Load: +0.1 s** here, the used cores collected and their starts placed.

**Not yet measured:** an r226 run (`measure_performance.sh` with `PROTAL_ARGS=--compact_keys` beside a run without),
which would also settle the page-table question. A line layout with each word's count beside the presence bits
(fewer popcounts per rank, which `Locate` and `PrefetchStart` each compute) could take back part of the HiFi cost.
