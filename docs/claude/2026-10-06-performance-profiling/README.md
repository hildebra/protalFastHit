# Performance profiling, 2026-10-06: where 0.7.7 spends its time, and what to optimise next

- **Date**: 2026-10-06.
- **Code**: branch `audit-fixes` at `f30d8dc` (0.7.7 and the commits after it), built from `git archive` in WSL
  (`scripts/build_head.sh`, Release, the one binary), and once more with `-g` on top of the Release flags for
  callgrind's line attribution (`scripts/build_g.sh`: the same code, with line tables). Nothing in the tree was
  changed. Two prototypes were built in copies of the tree (`scripts/build_proto.sh`); their diffs are
  `scripts/flex_ties_avx2.patch` and `scripts/screen_packed.patch`.
- **Machine**: WSL2 Ubuntu 24.04 on the Core Ultra 7 258V (6 vCPUs over 4 fast and 4 low-power cores), GCC 13.3. Other
  sessions ran tests and benchmarks throughout (load 1.4-6.9), and the machine changed speed during the session: the same
  six-thread runs took 1.6-2.3× longer in rounds 2 and 3 than in round 1, user time included. Wall times here are
  only compared within one alternated round; instruction counts (callgrind, one thread) are exact and carry the
  comparisons.
- **Data, local**: the v0.7.5 benchmark world (`~/bench071/V075/protal_db`, 765 species) and its samples:
  `rl150_p500000_s_1` (500k pairs; its R1 alone as single-end), `pb_b90000000_s_1` and `ont_b90000000_s_1`
  (90 Mb), for callgrind the first 100k pairs and the 3 Mb samples `pb_b3000000_s_1`, `ont_b3000000_s_1`.
  `--no_qcmsa` everywhere.
- **Data, GTDB scale**: the eleventh cluster run of the [GTDB-scale report](../2026-10-04-performance-gtdb-scale/README.md)
  (`2093770`, r226 v11, EPYC 9634, 32 threads, a real paired-end gut sample of 50.6M pairs and a HiFi barcode of
  497,656 reads; `results_v11/`). Since then only 0.7.7's two changes touched the run's speed (the read EM's cap of 100
  sweeps, the index load's allowance); nothing in the alignment stage changed.
- **Earlier rounds**: [round 1](../2026-09-29-performance-profiling/README.md), [round 2](../2026-09-30-performance-round2/README.md),
  [round 3](../2026-10-02-performance-round3/README.md), [review](../2026-10-03-performance-review/README.md),
  [pass 2](../2026-10-03-performance-pass2/README.md), [GTDB scale](../2026-10-04-performance-gtdb-scale/README.md).
- **Question**: profile protal again; where could the code be optimised further?

## Summary

1. **At GTDB scale the run is the alignment stage, and the alignment stage is seed and anchor finding.** The r226
   paired-end run takes 42 s: start-up 3.3 s, aligning 32.5 s, profiling 3.7 s. Per thread, the aligning is the
   seeding 48%, the alignment handler 12%, extending the anchors 11%, sorting the seeds 10%, taking the k-mers 6%,
   reading the input 3%, the rest ≤ 2% each. The HiFi run (15.6 s) is likewise 10.5 s of aligning.
2. **The seed lookup's own arithmetic is a large part of the seeding at r226, and a third of it is a scalar
   pass over the scores.** At r226 a lookup scans 70 flex cells (local world: 9). A benchmark of lookups at r226's
   block sizes, in a 6 GB arena with huge pages (`scripts/lookup_bench.cpp`), counts 1,156 instructions per lookup
   for the tree's code: 2.2G lookups make 2.6 T instructions, ~29% of the paired-end run's 8.76 T. Collecting the best
   cells with 32-byte compares and bit masks, with the AVX2 scoring's last step masked instead of a scalar tail,
   takes 661 (−43%), the same seeds (checksum of every field), 111-166 instead of 153-208 ns per lookup in cache.
   The prototype in the tree (`flex_ties_avx2.patch`) writes identical SAMs (paired-end and PacBio) and saves 17% of
   `GetFromLookup` on the local world's small blocks.
3. **An anchor's gene is reached through a chain of dependent cache misses.** `GetGenome` probes the 143,614-genome
   `tsl::sparse_map` twice (`contains`, then `at`), then the genome object, its gene vector, the gene record, the
   packed bases. At r226 size (`scripts/genome_bench.cpp`, the real `Genome` and `Gene` classes) that costs 220-232 ns
   per gene reached against 54-62 ns through a flat taxid → `Gene*` table, and a chain taken one after another
   0.95-1.25 µs. Extending the anchors at r226 takes 1.16 µs per mate: about what the anchors' gene accesses cost if
   they queue. A flat table and prefetching the genes of a read's anchors before extending them should remove most of
   it (the anchors per mate at r226 are not counted yet), and the same chain recurs in the alignment handler, the mate
   guidance and the output handler.
4. **The k-mer screen is most of the alignment handler at r226 and can be cut by a third at once.** 90% (paired-end)
   and 96% (HiFi) of the candidates are refused by the screen; it costs ~5.5k instructions per short-read candidate and
   ~33k per HiFi one (local counts), an estimated 11% and 8% of the two runs' instructions. Stamping the window's k-mers straight from the
   gene's 2-bit packed bytes (packed bases are all A, C, G or T, as the decoded window is) skips the decoding of refused
   candidates and the per-base re-encoding: the prototype (`screen_packed.patch`) writes identical SAMs and takes 33%
   off the screen. The read side is now two thirds of what is left; its k-mer codes from SIMD passes, once per read,
   would bring the screen to about a third of today's.
5. **Smaller or measured-small**: 64-bit seed sort keys (the same order; −24% instructions, −5-10% time on the sort,
   10% of aligning); prefetching every cache line of a block (slower: 366-442 against 267-307 ns per lookup);
   the syncmer scan (6% at r226, ~39 instructions a base in five passes); the FASTQ reader (3%).
6. **Locally (765 species) a short run is a third start-up.** 100k pairs at one thread are 29.4 G instructions, of
   which loading the index 9.6 G, the genome preload 1.1 G, parsing the gene-neighbour text table 0.84 G and the
   database tables 0.5 G: ~12 G fixed per run, which every simulated sample of a database build's training pays.
   Profiling a 500k-pair SAM is 15.8 G: reading it 10.2 G, the congener sketches 2.3 G (computed per run, could be
   computed at build), the read EM 1.2 G.
7. **One stage timer misleads**: the one-thread "Output handler" reads 7.2 s of an 18 s paired-end alignment while
   callgrind gives the output handler 9.5% of the alignment loop; the per-read sampled timer catches the buffer flushes
   (as the review of 2026-10-03 saw at six threads).

Items 2-4 are exact (outputs unchanged by construction; the two prototypes checked on SAMs). Estimated together they
take a fifth to a quarter off the r226 alignment stage; how much of the instruction savings turns into time where the
seeding also waits on memory is for the next cluster run to show.

## What I would do, in order

Items 1, 2, 3 and 5 were implemented the same day, outputs unchanged: see
[the follow-up](#follow-up-the-same-day-items-1-3-and-the-counters-implemented). After the twelfth cluster run, items 4
and 7 and ISA-L for gzip followed: see [that section](#after-the-twelfth-run-shared-seeds-the-reads-packed-strands-isa-l).

| # | change | where | exact | evidence | expected at r226 |
|---|---|---|---|---|---|
| 1 | best cells by AVX2 masks, masked last scoring step, index padding 16 → 64 bytes | `FlexScan.h`, `KmerLookup.h::GetFromLookup`, `Seedmap::AllocatePacked` | yes | bench −43% instructions per lookup, same seeds; prototype: SAMs identical, `GetFromLookup` −17% locally | 1.1 T of 8.76 T instructions (pe), 0.23 T of 4.0 T (HiFi); seeding (15.6 s per thread) −2 to −5 s |
| 2 | flat taxid → genome (or → first `Gene`) table in `GenomeLoader`, `GetGenome` with one probe where the map stays; prefetch the gene records and packed bases of a read's anchors before `ExtendAnchor` | `GenomeLoader.h`, `ChainAnchorFinder::operator()` | yes | bench: 220-232 → 54-62 ns per gene reached (independent), ~1.0-1.25 → ~0.8 µs (dependent chain) | extending 3.67 s per thread → ~1-1.5 s; some of the alignment handler and output |
| 3 | the screen from the packed gene bytes, the window decoded only for candidates that pass | `AlignmentScreen.h`, `AlignAnchor` | yes | prototype: SAMs identical, screen −33% | ~4% of the pe run's instructions, ~3% of HiFi's |
| 4 | the screen's read side from k-mer codes made once per read (whole read and its reverse complement; long-read windows as offsets into them), exits tested every 16 k-mers | `AlignmentScreen.h`, `SimpleAlignmentHandler::operator()`, `LongReadAligner::Align` | yes (the exits are monotone) | per k-mer ~12 → ~5 instructions estimated from the line costs | with 3: the screen at about a third, ~7% of the pe run's instructions |
| 5 | counters for the next cluster run: anchors per mate, seeds in (taxid, gene) groups of 2+, lookups dropped as too ubiquitous, the screen's and the gene access's time | `ChainAnchorFinder`, `SimpleAlignmentHandler` | – | – | decides 6 and 7 |
| 6 | 64-bit seed sort keys where taxon-gene + position bits fit (39 at r226) | `ChainAnchorFinder::Sort` | yes | bench: same order, −24% instructions, ~−8% time | sort (3.4 s per thread) −0.3 s |
| 7 | sort and pair only the seeds of (taxid, gene) groups of 2+ (FindPairs skips the others) | `ChainAnchorFinder` | SAM and profile yes; the verbose seed statistics change | none yet: depends on how many seeds are singletons at r226 (5) | up to most of the sort, if most seeds are singletons |
| 8 | the syncmer scan in fewer passes | `KmerIterator.h` | yes | 5.9k instructions per 150-base read | k-mers (2.0 s per thread) −0.5-1 s |
| 9 | the congener sketches stored in the database | `SampleContext.h`, the build | yes | 2.3 G of a 500k-pair profiling's 15.8 G | small at r226 (0.23 s) |

Not worth it on this evidence: prefetching every cache line of a block (bench: slower), and the sort beyond 64-bit
keys. The non-exact levers are in [their own section](#non-exact-levers).

## Where the time goes

### At GTDB scale (the eleventh cluster run, run 1 of 3)

Paired-end, 50,607,792 pairs (101.2M mates), 32 threads: run 41.3 s (start-up 3.3 s, aligning 32.5 s, writing the
SAM 0.7 s, profiling 3.7 s). Per-thread seconds of the aligning (`results_v11/stages.tsv`):

| stage | s per thread | share of aligning | per unit |
|---|---|---|---|
| seeding (lookups and their scans) | 15.65 | 48% | 225 ns per lookup (2.22G lookups, 70.4 cells each) |
| alignment handler | 3.92 | 12% | 713 ns per candidate (176.0M, 90.4% refused by the screen) |
| extending anchors | 3.67 | 11% | 1.16 µs per mate |
| sorting seeds | 3.40 | 10% | 24.9 ns per seed (4.37G, 43 per mate) |
| taking the k-mers (syncmer scan) | 1.98 | 6% | 625 ns per mate |
| sequence reader | 1.02 | 3% | 323 ns per mate |
| output handler | 0.65 | 2% | |
| pairing (FindPairs) | 0.64 | 2% | |
| anchor recovery | 0.50 | 1.5% | |
| joining pairs, sorting | 0.48 | 1.5% | |
| sorting anchors, rest of the anchor finder | 0.80 | 2.5% | |

`perf stat` for the run: 8.76 T instructions, 4.46 T cycles (IPC 1.97), 1,218 s user (3.66 GHz).

HiFi, 497,656 reads, 32 threads: run 15.6 s (start-up 3.3 s, aligning 10.5 s, SAM 0.7 s, profiling 0.4 s). Long reads
are seeded inside the long-read aligner, so the "alignment handler" timer (10.41 s per thread) holds the seeding 3.26 s,
the seed sort 1.26 s, extending 0.97 s and pairing 0.19 s; the alignment proper is ~4.7 s per thread for 9.69M
candidates, of which 9.26M (95.6%) refused by the screen. 4.0 T instructions, IPC 2.52.

### Locally: whole runs (round 1, six threads; the fast phase of the machine)

| sample | wall | user | max RSS | index load | aligning | profiling |
|---|---|---|---|---|---|---|
| paired-end 500k pairs | 6.08 s | 28.0 s | 3.55 GB | 1.02 s | 3.88 s | 1.0 s |
| single-end 500k reads | 2.67 s | 11.9 s | 3.51 GB | 0.59 s | 1.39 s | 0.6 s |
| PacBio 90 Mb | 3.70 s | 19.5 s | 3.51 GB | 0.41 s | 2.73 s | 0.4 s |
| ONT 90 Mb | 8.66 s | 47.7 s | 3.63 GB | 0.38 s | 7.33 s | 0.7 s |

Rounds 2 and 3 (`results/runs.txt`) repeat these 1.6-2.3× slower with 1.7× the user time: the machine, not the code.
The local world differs from r226 where it matters most: 37% of the k-mers have a core in the index (r226: 99.2%),
9.3 flex cells per block (70.4), 5% of candidates refused by the screen (90%), the genes 0.02 GB packed and in cache
(4.1 GB).

### Locally: instructions (callgrind, one thread; `results/callgrind/`)

Aligning 100k pairs without profiling, 29.35 G:

| part | G | share |
|---|---|---|
| index load (`Seedmap::Load`: `DecodeChunk` 4.9, zstd 2.6+0.4, packing 1.4) | 9.65 | 33% |
| genome preload, gene-neighbour table (text), database tables | 1.10 + 0.84 + 0.50 | 8% |
| input inflate (its own thread) | 1.12 | 4% |
| `RunPairedEnd` | 16.05 | 55% |
| – alignment handler: WFA2 4.93 (flanks 4.61), the screen 1.34, post-processing 0.53, CIGAR counts 0.43 | 7.91 | |
| – anchor finder: lookups 0.77 + 0.37, pairing 0.57, sort 0.53, extending 0.46, anchors 0.34, reverse complement 0.25 | 2.98 | |
| – output handler (SAM zstd 0.82) | 1.53 | |
| – syncmer scan | 1.19 | |
| – FASTQ reader | 0.53 | |
| – mate guidance | 0.46 | |

Profiling the 500k-pair SAM, 19.20 G (15.82 G in `ProfileSample`): reading the SAM 10.18 G (its zstd 2.16, `PrepareMAPQ`
2.80, `Taxon::AddSam` 2.38, `NoteRecord` 1.63), record evidence 3.60 G, sample context 3.57 G (the congener sketches
2.28 G for 25,887 genes, 87k instructions each; the read EM 1.21 G), writing the profiles 1.42 G, SNP post-processing
0.58 G.

Long reads, whole runs of 3 Mb: PacBio 17.75 G (index 9.65, `RunLongReads` 3.23: the long-read aligner's alignments
1.95 of which WFA2 pieces 1.13 and the screen 0.20, seeding 0.53, neighbour genes 0.37; profiling 1.85); ONT 24.09 G
(`RunLongReads` 9.10, the screen never refuses at ONT's error rate).

### Since performance pass 2

Pass 2's build (`0e4c5a4`, still in `~/perf-pass2`) and `f30d8dc` on the same 100k pairs and the v0.7.3 database it
used (`scripts/drift.sh`, `results/callgrind/drift_*.incl.txt`):

| G instructions | `0e4c5a4` (pass 2) | `f30d8dc` | |
|---|---|---|---|
| whole run | 27.15 | 29.35 | +8.1% |
| index load | 7.33 | 9.65 | +2.32: the packed index built during the load (`3565360`, `959d46c`, `b454c97`; −8 GB at r226) |
| `RunPairedEnd` | 16.16 | 16.06 | −0.6% |
| – alignment handler | 8.15 | 7.91 | the screen (+1.34) pays for itself in WFA2 even here |
| – anchor finder | 2.79 | 2.98 | +7% (among the changes: tied seeds kept for their own anchors, `fdad375`; not broken down) |
| – output handler | 1.59 | 1.53 | |
| – syncmer scan, mate guidance | 1.18, 0.46 | 1.19, 0.46 | |

So the alignment loop costs what it cost three days ago on the small world, while it got much cheaper at r226 (pe 131
→ 42 s) through changes the small world barely exercises (the screen, the SAM read, the start-up). The local world
is no longer a guide to where r226's time goes: the r226 numbers above are.

## The findings in detail

### 1. The seed lookup: the best cells collected with masks

`KmerLookupSM::GetFromLookup` scores a block's flex cells with `flex_scan::ScoreAvx2` (8 cells a step, then up to 7
cells one by one, then a pass that counts the cells with the best score: 32 at a time, the rest one by one), and
then walks all cells again, one by one, for those with the best score. At r226's 70 cells per block the two scalar
parts and the walk are as large as the vector scoring.

`scripts/lookup_bench.cpp` builds blocks as the packed index holds them (32-bit flex cells at bit offsets, then
41-bit entries: r226's layout), sized by r226's shares of lookups per block-size group and the groups' mean sizes
(70.3 cells per block in a 6 GB arena), with 1-3 cells near the read's key so that ~1.7 seeds per lookup result, and
takes them as `FindSeeds` does (a mate's 22 lookups smallest first, values prefetched 16 ahead, a pair's 44 prefetched
first). Variants, all with the same seeds (a checksum of every field of every seed, in order):

| variant | instructions per lookup | ns per lookup, in cache | 6 GB arena, load 1.3 | 6 GB arena, load 4.3 |
|---|---|---|---|---|
| A: the tree's | 1,156 | 153-208 | 267-307 | 293-476 |
| B: A's scoring, the best cells by 32-byte compares and bit masks | 771 | 121-187 | 226-279 | 281-455 |
| C: B, and the scoring's last step masked (no scalar tail, no counting pass: B's masks count) | 661 | 111-166 | 229-251 | 226-522 |
| D: C, every cache line of the block prefetched (up to 16) | – | 113-133 | 366-442 | – |

Ranges are over 3 (`lookup_bench_abcd.txt`) and 6 (`lookup_bench_abc.txt`) alternated rounds; with other sessions'
work on the machine (load 4.3) the arena runs say little, the instruction counts (callgrind, `lookup_bench_callgrind.txt`)
are exact. The arena's lookups cost about what r226's do per thread (225 ns), which
suggests the bench's mix of memory and arithmetic is close to the cluster's.

Variant C needs reads up to 32 bytes past a block's last cell: the block's entries follow it, and at the end of the
index the padding grows from 16 to 64 bytes. The prototype (`scripts/flex_ties_avx2.patch`: `flex_scan::BestAvx2`,
`flex_scan::TiesAvx2`, the AVX2 branch of `GetFromLookup`, the padding) writes the same SAM records as `f30d8dc` on
100k pairs and on PacBio 3 Mb; locally (8 cells per block) `GetFromLookup` falls from 607 to 506 M instructions
(−17%), the anchor finder from 2.98 to 2.88 G. At r226 the bench's −495 instructions per lookup are 1.1 T of the
paired-end run's 8.76 T. A test like `FlexScan.Avx2ScoresAsTheScalarScan` for the masked scoring and the ties (every
size 1-80 and the large ones, every bit shift) would come with it.

Prefetching more of a block (D) made the arena runs slower: with 16 lookups ahead and ~5 lines each, the prefetches
outrun the core's miss buffers and push out what is about to be used. The first line of the cells and of the
entries, as now, is right.

### 2. Reaching an anchor's gene

Every anchor is extended (`ChainAnchorFinder::operator()`) before the anchors are ranked, so every anchor's gene is
read: `GenomeLoader::GetGenome` (`m_genomes.contains(key)`, then `m_genomes.at(key)`: two probes of a
`tsl::sparse_map`, whose buckets are a bitmap and a value array), `Genome::GetGeneOMP` (the loaded flag, the gene
vector), the `Gene` record (32 bytes; 14.5M of them, 465 MB at r226), then `Gene::Window` decodes the packed bases
(4.1 GB at r226). At r226 each of these is a cache miss, and each depends on the one before.

`scripts/genome_bench.cpp` holds 143,614 real `Genome` objects (232 bytes each) of 101 `Gene`s in the sparse map, the
packed bases in a 2 GB arena, and reads a gene's length and first packed word for 4M random (taxid, gene) pairs:

| path | independent lookups (ns each) | dependent chain (ns each) |
|---|---|---|
| `GetGenome` + `GetGeneOMP` (as the tree) | 220-232 | 950-1,250 |
| one `find` + `GetGeneOMP` | 149-198 | 880-1,120 |
| flat taxid → `Genome*` | 81-89 | 1,040-1,110 |
| flat taxid → `Gene*` (the genome's gene array) | 54-62 | 784-802 |

(The quiet run, load 1.45; a second run under load 6.9 keeps the order, `results/bench/genome_bench*.txt`.) Extending
the anchors at r226 costs 1.16 µs per mate. The extension loop takes a read's anchors one after another, and each
anchor's decoding and extension is longer than the core's out-of-order window, so the next anchor's misses start only
when the last one is done: the chains queue, as in the right-hand column (which is WSL's memory latency; the EPYC
node's is lower, its index and genes four times larger). The anchors per mate are not counted at r226 (item 5);
locally there are 3.4. That the extension is mostly these misses is the reading the numbers fit, not a measurement:
the next cluster run's timers around the gene access would confirm it.

The fix in two parts, both exact:

- a flat table from taxid to the genome's gene array (taxids up to 2^18 at r226: 2 MB of pointers), filled once the
  genome map stops changing (after the preload), used by `GetGenome`/`GetGeneOMP` and by a `GeneOf(taxid, gene)` for
  the hot paths; `GetGenome` with one `find` where the map stays;
- in `operator()`, after `FindPairs`: prefetch every anchor's `Gene` record, then (a second loop) the first lines of
  each anchor's window in the packed bases, then extend. The misses of a read's anchors then overlap instead of
  queueing.

The same chain recurs for each candidate in `AlignAnchor` (cheap there when the extension has just touched the gene),
in `GuideMate`, `RecoverAnchors` (`ExtractAndExtendAnchorFromSeed`), the long-read aligner's `AddCandidate`, and in the
output handler's `ReferenceOf` per record.

### 3. The k-mer screen from the packed bases

`AlignAnchor` decodes the gene window, then `AlignmentScreen::MayAlign` stamps the window's k-mers (6-mers, 7-mers for
windows over 512 bases) into a table, re-encoding every decoded base (`BaseToInt`, a validity count), then rolls along
the read's aligned part counting the k-mers the window has, with two exits. Line costs on 100k pairs (`-g` build,
`results/callgrind/lines_align_AlignmentScreen.txt`): the stamping ~0.69 G, the read pass ~0.62 G, of the screen's 1.34 G
(5.5k instructions per candidate; PacBio: 33.5k).

At r226 the screen refuses 90.4% (paired-end) and 95.6% (HiFi) of the candidates. With the local per-candidate costs:
176.0M × 5.5k ≈ 0.97 T instructions (11% of the paired-end run's) and 9.69M × 33.5k ≈ 0.32 T (8% of the HiFi run's).

The packed bases have no character other than A, C, G and T (`PackedSequence.h`: N and the ambiguity codes are stored as
bases), as the decoded window has none, so the window's k-mers can be taken from the packed bytes: four k-mers per
8-byte load (`(w >> 2j) & mask`, the first base lowest; the read's k-mers coded the same way), and the window decoded
only for a candidate that passes. The prototype (`scripts/screen_packed.patch`: `AlignmentScreen::MayAlignPacked`, the
call before `gene.Window` in `AlignAnchor`, the old path kept for a gene that is not loaded) writes the same SAM records
on 100k pairs and PacBio 3 Mb:

| | `f30d8dc` | prototype |
|---|---|---|
| paired-end 100k pairs: the screen | 1.342 G | 0.893 G (−33%) |
| – `RunPairedEnd` | 16.06 G | 15.60 G (−2.8%) |
| PacBio 3 Mb: `AlignAnchor` (the screen 0.20 G) | 1.868 G | 1.801 G |

What is left is two thirds the read pass, ~12 instructions per k-mer: a table lookup, a branch on the base, the
rolling code, the stamp, two exit tests. The read is the same for all candidates of a mate (2.2 per mate at r226), and
a long read's candidate windows are pieces of the read and of its reverse complement (19.5 candidates per HiFi read).
Its k-mer codes, made once per read for both strands by SIMD passes (as `SimpleKmerHandler::FillAvx2` makes the
syncmers' codes) with a code for "has a base other than ACGT", and the exits tested every 16 k-mers (shared k-mers
only grow and shared plus remaining only shrink, so the answer is the same), would leave ~5 instructions per k-mer.
With the packed stamping that is about a third of today's screen.

### 4. The seed sort

`ChainAnchorFinder::Sort` packs each seed into a 128-bit key, sorts, unpacks (since `22c445d`). At r226 a taxon and gene
fit 25 bits as the packed index's `taxid * genes + gene` (which orders as (taxid, gene) does), the gene position 14, the
read position 16, the flags 2: 57 bits, so a 64-bit key holds the same order. `scripts/sort_bench.cpp` (51 seeds per
mate, a third of the mates with an anchor's seeds among random ones): the same order, 24% fewer instructions, 29.1-32.9
against 31.2-33.2 ns per seed when the machine was quiet (no clear difference under load). Most of the sort's time is its
comparisons' branches on data, which the key width does not change. Worth it as a small exact change, with the 128-bit
key kept where the bits do not fit.

A larger lever is what is sorted: `FindPairs` uses only seeds that share their taxon and gene with another seed. If at
r226 most seeds are such singletons (a real sample's reads are mostly not from the database's genes, and a lookup
whose core matches but whose flex part shares 8 of 16 bases still gives seeds), counting (taxid, gene) in a small
per-read hash first and sorting the rest would remove most of the 3.4 s per thread. The seed counts in the verbose
statistics would change; the SAM would not. The next cluster run should count them first (item 5).

### 5. Smaller items

- **The syncmer scan** (`SimpleKmerHandler::ScanWindowsAvx2`) takes 5.9k instructions per 150-base read (39 per base):
  per-base codes, then `Bases16`'s four doubling loops per strand (widening `uint8`/`uint16`/`uint32` arrays,
  ~120 instructions each), the k-mers and cores per window, the s-mers, then the windows 8 at a time. Fusing the
  doubling into 32-bit lanes and the k-mer, core and s-mer loops into one pass per strand would roughly halve it: up to
  1 s of r226's 2.0 s per thread.
- **The FASTQ reader**: 2.6k instructions per mate locally (line splitting, the record's three strings copied), 1.0 s per
  thread at r226. Small.
- **The output handler's timer** (7.2 s of 18 s at one thread locally, against 9.5% of the loop's instructions): the
  sampled per-read timer counts a whole flush when it samples one. Timing the flushes apart would make it trustworthy.
- **Start-up for small samples**: ~12 G instructions per run on the local world (above). The gene-neighbour table is
  parsed from text in every run (0.84 G here; `gene_neighbours::Table::Read`); a binary member as `gene_table.bin` did
  for the gene tables would remove most of that. The database build runs protal on thousands of small simulated
  samples, where this is a real share.
- **Profiling**: the congener sketches (`CongenerDistances::MakeSketch`, 2.3 G of a 500k-pair profiling's 15.8 G
  locally, 0.23 s at r226) depend only on the reference and could be stored at build; `PrepareMAPQ` is 2.8 G and was not
  looked into.

### Non-exact levers

At r226 99.2% of the read k-mers have their 15-base core in the index, and a lookup keeps the entries of the best flex
score whatever it is, unless more than 256 tie. Reads that come from none of the database's genes (most of a real
sample) still make ~2 seeds per lookup, anchors on 78.6% of the mates, and 176M candidates of which 90% the screen
refuses. A minimum flex score for a lookup's seeds (say 12 of 16 bases), or dropping a read's lookups once its seeds
scatter over many genes, would cut the seeding's emits, the sort, the anchors and the screen together, but changes
which divergent strains are found: a question for the r226 evaluation, not for this report. More threads per node
(the EPYC node has 84 cores; the runs use 32) is a measurement, item 2 of the GTDB-scale report's seeding options.

## Follow-up the same day: items 1-3 and the counters implemented

At the user's request, on `1750475` (HEAD then; it changed none of the files below). Exact: every output the same.

**1. The seed lookup's best cells with AVX2 masks** (`FlexScan.h`, `KmerLookup.h`, `Seedmap.h`).
`flex_scan::ScoreAvx2` (scoring, a scalar tail, a counting pass) is replaced by `BestAvx2` (the scores and the best in one
pass, the last step masked: lanes past the block score 0 and do not raise the best) and `TiesAvx2` (the cells with the
best score as bit masks, 32 scores a step, and their count). `KmerLookupSM::GetFromLookup` takes the entries from the
masks in entry order; the scalar path (`ScoreScalar`, no AVX2 or `PROTAL_FLEX_SCAN=scalar`) builds the same masks, so one
emission loop serves both. `GetFromLookup` now returns false for a lookup dropped as too ubiquitous. The index's padding
after the packed values grows from 16 to 64 bytes (`Seedmap::kPackedPadding`; `BestAvx2` reads up to 32 bytes past a
block's last cell, a `static_assert` ties the two). A first version took the emission through a lambda, which GCC did
not inline into the target-cloned `GetFromLookup` (331 M instructions in the lambda alone, no gain); the shared masks
removed it.

**2. Reaching an anchor's gene** (`GenomeLoader.h`, `ChainAnchorFinder.h`, `AlignmentStrategy.h`, `LongReads.h`,
`MateGuidance.h`). The genome loader keeps a flat taxid → genome table (made when a load of the genes completes the
genome map; anything that adds to or clears the map empties it first) and, once `LoadAllGenomes` has loaded every genome,
a taxid → gene list table published through an atomic flag (the tables load beside the preload). `GetGenome` (both) uses
the first, else one `find` (it probed the map twice); `GetGeneOMP(taxid, gene)` the second, else the genome as before
(on-demand loading unchanged); `PreloadedGene` and `PrefetchGene` give a gene's place without reading it, and
`Gene::PrefetchBases` fetches up to four lines of its packed bases. The anchor finder's extension loop fetches the gene
records 8 anchors ahead and the bases around the anchors 4 ahead; the anchor finder, `AlignAnchor`, the handler's
`ExtendAnchor`, the long reads' candidates and neighbour genes and the mate guidance reach genes through
`GetGeneOMP(taxid, gene)`.

**3. The k-mer screen from the packed bases** (`AlignmentScreen.h`, `AlignAnchor`). `MayAlignPacked` stamps the window's
k-mers from the gene's packed bytes, four per 8-byte load; `AlignAnchor` screens before it decodes the window (a gene
without packed bytes goes through the decoded window as before). Both screens share the bound and the read pass, with
k-mers coded first base lowest on both sides (`Roll`).

**Counters** (`ChainAnchorFinder.h`, `AlignmentStrategy.h`, `Classify.h`, `measure_performance.sh`, `docs/running.md`,
`docs/development.md`). The seeding line ends "..., N seeds, P of them sharing their taxon and gene with another seed;
D lookups dropped as too ubiquitous; A anchors"; a "K-mer screen" timer (sampled per read, inside the alignment
handler's) is printed with `--verbose` and written to `misc/<sample>_runtime.tsv`, so `stages.tsv` has it;
`measure_performance.sh` adds the columns `paired_seeds`, `dropped_lookups` and `anchors` at the end of `runs.tsv` (NA
for an older protal). Locally (pe 500k): 21.8M seeds, 18.0M (83%) in shared genes, 965 lookups dropped, 3.41M anchors
(3.4 per mate). The gene accesses have no timer of their own: with the prefetch their misses move, so "Extending Anchors"
is the measure.

**Tests** (`results/implementation/`): `FlexScan.Avx2ScoresAndTiesAsTheScalarScan` (every size 1-80 and 13 larger ones,
every bit shift: scores, best, masks and count as `ScoreScalar`, garbage past the block masked out, no mask word written
past the block's), `ABlockWithoutASharedBaseHasEveryCellBest`, the bench (`PROTAL_FLEX_BENCH=1`, since 2026-10-07 [`scripts/flex_scan_bench.cpp`](scripts/flex_scan_bench.cpp): 0.41 against 2.0-2.2 ns
per cell for 70-cell blocks, scoring and the ties together, where `ScoreAvx2` alone took 1.07 on 2026-10-05);
`PackedIndex.LookupsGiveTheSameSeedsWithAndWithoutAvx2` (a packed index of 3,000 keys, keys whose flex part matches
exactly, in a few bases or not at all, ubiquity 256 and 3: the same seeds, flags and order, the same lookups dropped);
`AlignmentScreen.PackedWindowsGiveTheDecodedWindowsAnswer` (48,000 screens of genes with Ns and ambiguity codes,
windows at every offset up to the gene's last byte, empty or shorter than k, reads with Ns, free ends, budgets, k 4-10:
the same answer as on the decoded window; 9,410 refused); `GenomeLoaderPacked.TheFlatTablesReachTheGenomesAndGenesOfTheMap`
(preloaded and on demand, after the sequences are freed, gene 0 and past the list, unknown taxids).

| check (`scripts/build7.sh`, `scripts/check7.sh`) | result |
|---|---|
| unit tests | 369 passed, 3 skipped (the two opt-in benches, the AddressSanitizer-only test) |
| end-to-end tests on the mini database | 133 passed |
| pe 500k pairs, se 500k, PacBio and ONT 90 Mb at 6 threads, against `1750475` | SAM records identical (sorted), every other output file byte for byte; the same with `PROTAL_FLEX_SCAN=scalar` |
| `measure_performance.sh` | the three new columns filled, the screen's timer in `stages.tsv` |

Instructions (callgrind, one thread, `results/implementation/callgrind_*`):

| G instructions | `1750475` | this change | |
|---|---|---|---|
| aligning 100k pairs | 29.35 | 28.82 | −1.8% |
| – `RunPairedEnd` | 16.06 | 15.53 | −3.3% |
| – anchor finder | 2.98 | 2.89 | −3.0% (`GetFromLookup` 0.607 → 0.510, −16%; the prefetch loops added) |
| – the k-mer screen | 1.342 | 0.896 | −33% |
| PacBio 3 Mb, `RunLongReads` | 3.23 | 3.15 | −2.7% (the screen −0.067, the anchor finder −0.017) |

Locally the blocks have 9 cells and the genes are in cache, so this is the small end. At r226 (70 cells a block, genes
spread over 4.1 GB, 90-96% of the candidates screened) the bench numbers above apply: ~−43% of the lookup's instructions
(1.1 T of the paired-end run's 8.76 T), the gene chain at a quarter of its cost where it is waited for, the screen −33%.
Wall times here were not measurable (load 5-9 from other sessions); the next cluster run (`scripts/measure_performance.sh`
on the r226 database, as the eleventh) shows the time, the anchors per mate and the share of seeds in shared genes,
which decide items 6 and 7.

## The twelfth cluster run: `cdce3c0` at r226

SLURM job 24004543 (the user's), `scripts/measure_performance.sh` at `cdce3c0` on node `q512n17`, 32 threads, the r226 v11
database and the same two samples as the eleventh run; [`results/cluster_v12/`](results/cluster_v12/) (copied unchanged from
`local/Perf0.7.5v12/`). The first paired-end run loaded the index cold from NFS (251 s) and is left out of the paired-end
medians. Against the eleventh run (`2093770`), which predates 0.7.7: of the changes between them only `cdce3c0` touches the
alignment stage. 0.7.7's two (the index load's allowance, the read EM capped at 100 sweeps) and the gradient-boosted
models of `1750475` touch the load and the profiling stage.

| | eleventh (`2093770`) | twelfth (`cdce3c0`) | |
|---|---|---|---|
| pe wall (median) | 42.2 s | **36.6 s** (36.2, 37.0) | −13% |
| pe aligning | 32.5 s | 28.4-29.0 s | −11% |
| pe user time; instructions | 1,218 s; 8.76 T | 1,051 s; 7.33 T | −14%; −16% |
| HiFi wall | 15.6 s | **13.7 s** | −12% |
| HiFi aligning | 10.5 s | 9.4 s | −10% |
| HiFi instructions | 4.00 T | 3.60 T | −10% |
| index load | 3.26 s | 2.66 s | 0.7.7's allowance, as expected |
| pe profiling | 3.8 s | 3.1 s | read EM 1.18 → 0.64 s (0.7.7's cap) |
| peak memory | 34.2 GB | 34.2 GB | |

Per thread (medians over the runs, `stages.tsv`):

| stage, s per thread | pe eleventh | pe twelfth | | HiFi eleventh | HiFi twelfth | |
|---|---|---|---|---|---|---|
| seeding | 15.63 | 12.51 | −20% | 3.26 | 2.83 | −13% |
| extending anchors | 3.67 | 2.05 | −44% | 0.97 | 0.55 | −43% |
| alignment handler | 3.92 | 3.36 | −14% | 10.45 | 9.28 | −11% |
| – of which the k-mer screen (new timer) | | 1.91 | | | 0.99 | |
| sorting seeds | 3.40 | 3.39 | | 1.26 | 1.26 | |
| taking the k-mers | 1.98 | 1.94 | | | | |
| sequence reader | 1.08 | **2.87** | +166% | 0.48 | 0.51 | |
| output handler | 0.65 | 0.89 | +37% | 0.10 | 0.10 | |

(For HiFi the seeding, sorting and extending run inside the alignment handler's time.)

- **The predictions held.** The seeding fell by 3.1 s per thread (predicted 2-5 s), and the extension by 44%: the queued
  gene misses were most of it, as the numbers suggested. The k-mer screen, timed for the first time, is 57% of the
  paired-end alignment handler after its own −33%. The instructions fell by 1.43 T (pe), of which ~1.1 T were predicted
  for the lookups alone.
- **The paired-end run is now limited by its gzip input.** The sequence reader's time per thread rose by 1.8 s while
  everything around it fell: the aligners wait for reads. The run aligned 50.6M pairs in 28.4-29.0 s, 1.75-1.78M pairs/s,
  and each FASTQ file is inflated by one zlib-ng thread, whose ceiling the multithreading audit put at 1.8-2.0M pairs/s
  ([2026-10-01](../2026-10-01-multithreading-audit/README.md), follow-up 2). Without that wait the aligning would be about
  27.2 s. So for gzip input at 32 threads and more, further alignment savings show only once the input is faster:
  zstd input (`.fq.zst`, which protal reads since `e7b391e`; its frames could be inflated in parallel) or BGZF for
  measuring the alignment itself, and for gzip a faster single-stream inflater (ISA-L's igzip inflates about twice as
  fast as zlib-ng) or parallel inflation of a single gzip member (as rapidgzip does), a larger piece of work. The HiFi run
  is not input-bound (0.5 s per thread in the reader).
- **The counters at r226.** pe: 4.37G seeds, of which 550M (12.6%) share their taxon and gene with another seed of their
  read; 6.98M lookups dropped as too ubiquitous (0.31%); 211M anchors (2.09 per mate; an extension now costs ~310 ns of a
  thread per anchor). HiFi: 964M seeds, 155M (16.1%) shared; 124k lookups dropped; 46.6M anchors (93.6 per read). So 87%
  of the paired-end seeds (84% HiFi) are sorted only for `FindPairs` to skip them: item 7 (count each read's (taxid, gene)
  first, sort and pair only the shared ones) would take most of the 3.4 s per thread of sorting (1.3 s HiFi), exactly.
- **The strain stage** of the cohort run took 5.6 s against 3.0 s: building the MSAs (51 against 24 s summed) and qcMSA (a
  Python process per species: 92 against 44 s) both doubled, as on the tenth run's node: the node's NFS, not the code.

What is left, by what it would give at r226 now: the input path for gzip (pe at 32+ threads), sorting only the shared
seeds (pe −2-3 s, HiFi −1 s per thread, exact), the screen's read side from codes made once per read (−~1 s pe, −0.5 s
HiFi), then the seeding's remaining 12.5 s per thread (memory-bound lookups; interleaving reads, the options of the
GTDB-scale report).

## After the twelfth run: shared seeds, the read's packed strands, ISA-L

At the user's request, on `cdce3c0`: items 1 and 3 of the twelfth run's list (sort only the shared seeds; the screen's
read side once per read), and ISA-L for gzip in place of libdeflate, if protal can still be built statically (it can).
The first two are exact; ISA-L changes no content, only the bytes of the gzip protal writes.

### Sorting only the seeds FindPairs can pair

`SharedSeeds` (`ChainAnchorFinder.h`) counts a read's seeds per hash of their taxon and gene in byte counters (16-32
per seed, cleared per read) and sorts, as 128-bit keys as before, only the seeds whose counter reached two: every seed
whose gene another seed of the read has, and the 3-6% of the others whose counter another gene's seed shares.
`FindPairs` makes anchors of runs of two or more seeds of one gene only, and a gene's run is the same in both lists, so
it makes the same anchors in the same order; a seed of a gene of its own is a run of one and skipped. The caller's seed
list is left as found (only its length is used after the anchor finder). A first version counted taxon and gene exactly
in an open-addressing table; the byte counters cost the same at r226's mix and are simpler.

`scripts/shared_sort_bench.cpp` (`sort_bench.cpp`'s seed lists: 51 seeds per mate, 14.3% in shared genes, as the twelfth
run's 12.6%): the same groups for `FindPairs` (checksums equal), per seed (callgrind, the difference of one and two
rounds):

| | instructions per seed | seeds sorted |
|---|---|---|
| all seeds sorted (up to `cdce3c0`) | 163 | 100% |
| exact counting, the shared ones sorted | 66 | 14.3% |
| byte counters (implemented) | 65 | 17.6% |

The sort's time is mostly mispredicted branches, which fall with the seeds sorted, so at the twelfth run's 3.39 s per
thread (pe) about −2 s are expected, and about −0.8 s of HiFi's 1.26 s. Locally (83% of the seeds in shared genes: the
local world's reads come from its own genomes) the counting is a small extra cost: 528 → 599 M instructions on 100k
pairs (+13%; the exact counting 627 M), 121 → 138 M on PacBio 3 Mb, 0.2% and 0.1% of the runs.

### The screen's read side from the read's packed strands

`AlignmentScreen::ReadKmers` (`AlignmentScreen.h`) packs a strand of the read 2 bits a base (`packed::Pack`, AVX2, as the
genes are packed) with a bitmap of its bases other than A, C, G or T (AVX2 compares), when a candidate first needs it,
once per read for all its candidates and every k. The k-mer at base p is then the 2k bits at bit 2 (p % 4) of the 8 bytes
from byte p / 4, four k-mers to a load as the window's (`StampPacked`); a k-mer with another base counts as shared
(`Any(k)`, a stamp index always set). The exits are tested every 16 k-mers (shared k-mers only grow and shared plus
remaining only shrink, so the answer is the same). `SimpleAlignmentHandler::operator()` sets the read and its reverse
complement (the anchor finder's); `LongReadAligner` sets the read once and screens each candidate's window `[s, e)` as
the stretch from `s` of the forward strand or from `n - e` of the reverse one, made from the read when first needed
(`AlignmentScreen::ReadStretch`, a new last argument of `AlignAnchor`).

A first version made per-k code tables (32-bit codes rolled base by base, blocks of 64 made on demand): ~10 instructions
per base and k, which at r226 (2.2 candidates per mate, 90% refused after most of the read) would have cancelled the
gain. Callgrind (one thread):

| instructions | `cdce3c0` | code tables | packed strands (implemented) |
|---|---|---|---|
| pe 100k pairs: the screen (243,020 candidates, 5% refused) | 896 M | 748 M (−17%) | 612 M (−32%) |
| – of which the read side | ~600 M (two thirds, from the line costs) | 407 M | 271 M |
| PacBio 3 Mb: the screen (5,977 candidates, 20% refused) | 134 M | not measured | 99 M (−26%) |

What is left of the screen is now mostly the window's stamping. At r226, where most candidates are refused after most of
the read, the read side falls further than locally, where most pass early: about −0.6 s of the paired-end screen's 1.91 s
per thread and −0.5 s of HiFi's 0.99 s are expected.

### ISA-L for gzip, in place of libdeflate

**Static builds still work.** ISA-L builds as a static library (`libisal.a`, ~0.8 MB) with its own CMake build (since
2.32) and dispatches its SSE/AVX2/AVX-512 code at run time, so the static binaries still run on any x86-64; a fully
`-static` link works. Its x86-64 code is nasm assembly: `just static` (`lib/static-deps.cmake`) uses the nasm on the PATH
if it is 2.14.01 or later and otherwise downloads and builds nasm 2.16.03 too (a C compiler and make, ~30 s), as it does
zstd. Ubuntu's `libisal-dev` (2.31.0) has no `libisal.a`, so the static build always builds ISA-L 2.32.1 itself; without
an ISA-L on the system, the other binaries link that one too (no `libisal-dev` needed for `just static`). A normal build
needs ISA-L like zstd (Ubuntu `libisal-dev`, conda-forge `isa-l`, or a prefix of one's own through `CMAKE_PREFIX_PATH`).

**Where it is used** (`Bgzf.h`, `ThreadedGzStream.h`). Every gzip input: BGZF blocks whole (`isal_inflate_stateless`,
raw deflate, then `crc32_gzip_refl`), other gzip as a stream (`isal_inflate`), which zlib-ng did before. BGZF output
(`.sam.gz`, the simulator's `.fq.gz`, the in-place compressor) at ISA-L's level 1 with a fixed level buffer
(`isal_deflate_stateless`; a block that does not shrink is stored by ISA-L itself). zlib-ng stays for gzstream.

**Measured** (single thread, WSL under load 18-20, so the ratios matter: a 477 MB FASTQ, as one `gzip -6` member and as
libdeflate's BGZF; `~/isal`, 10 alternated runs):

| | ISA-L | before | |
|---|---|---|---|
| inflating one gzip member, MB/s (best / median) | 1,189 / 896 | zlib-ng 651 / 552 | 1.6-1.8× |
| inflating BGZF as a stream | 998 / 769 | zlib-ng 502 / 414 | ~1.9× |
| inflating BGZF block by block | 869 / 648 | libdeflate 1,048 / 653 | libdeflate 0-20% faster |
| deflating 64 KB blocks, MB/s | level 1: 282 / 275 | libdeflate level 6: 40 / 37 | ~7× |
| compressed size | level 1: 138.4 MB | libdeflate level 6: 119.0 MB | +16% |

ISA-L's levels 2 and 3 compress no better than level 1 on FASTQ, and level 3 is slower. Through protal's own writer and
reader (`scripts/read_bench.sh`, `results/isal/read_bench.log`; HEAD's build against this one, alternated, load 7-16):

| the 477 MB FASTQ | before | ISA-L |
|---|---|---|
| written as BGZF (`bgzf::CompressFile`), one thread: CPU s | 17.7 (libdeflate level 6) | 2.6 (level 1) |
| – six threads: wall s, CPU s | 19.3, 27.8 | 3.8, 3.5 |
| – size | 119.2 MB | 138.6 MB (+16%) |
| zstd -3, one thread, for the simulator's other format: user s, size | 5.3-8.7, 118.0 MB | |
| read by `ThreadedGzStream` (gzip member, BGZF, zstd), MB/s | 235-412, 289-456, 164-227 | 386-401, 289-515, 207-389 |

The reading rounds spread more than the formats differ (the old build read the same gzip member at 235 and 412 MB/s):
under this load the reader shows no difference, and the inflate benchmark above (ten alternated runs of the inflate
alone) is the measure. Written BGZF is read back to the same content (`zcat`, md5).

**Two things found in ISA-L 2.32.1.** `isal_inflate` reading a gzip header itself (`crc_flag = ISAL_GZIP`) keeps the
header's state in a local of one call: a header with a header CRC, or with two of extra field, name and comment, that is
split between two reads of the input fails (CRC errors or invalid blocks; plain headers, a name alone and BGZF's extra
field alone are fine). The reader keeps its own `isal_gzip_header`, reads the header with `isal_read_gzip_header` until it
is complete, and inflates with `ISAL_GZIP_NO_HDR_VER`; at a member's end ISA-L gives back what it read past the trailer, so
the next member starts where the last ended (the byte offsets in the reader's messages stay exact). And the compressed
bytes are not the same on every CPU at level 3 (its AVX2 path matches differently); levels 1 and 2 gave the same bytes
through ISA-L's SSE4.2, AVX and AVX2 code (AVX-512 was not available to check; CPUs without SSE4.2 hash differently). The
output also depends on the level buffer's size, so it is fixed. For one ISA-L version, the gzip protal writes is thus the
same for any thread count and on those CPUs; its content is the same everywhere.

**Checks** (`scripts/build8.sh ref|work|isal`, `scripts/check8.sh`, `scripts/check_isal.sh`, `scripts/check_sfetch.sh`;
`results/isal/`):

| check | result |
|---|---|
| unit tests, all three changes | 373 passed, 3 skipped (the two opt-in benches, the AddressSanitizer-only test) |
| end-to-end tests on the mini database | 133 passed (items 1 and 3; all three) |
| pe 500k pairs (BGZF), se 500k, PacBio and ONT 90 Mb at 6 threads: items 1 and 3 against `cdce3c0` | SAM records identical (sorted), every other output file byte for byte; seeds, shared seeds, dropped lookups and anchors the same |
| pe 100k pairs (one gzip member), pe and se 500k (BGZF), PacBio 90 Mb (one gzip member with a name): ISA-L against zlib-ng and libdeflate | the same |
| a `.sam.gz` written with ISA-L | `gzip -t` passes; the same records as the `.sam.zst`; the profile files the same |
| `-DPROTAL_STATIC_FETCH_DEPS=ON` without nasm on the PATH or ISA-L on the system (local tarballs) | nasm 2.16.03 and ISA-L 2.32.1 (`libisal.a` 0.8 MB) built, protal linked to that ISA-L; the same SAM records. (A first try failed: CMake gives extracted files the extraction's time, and nasm's Makefile then remakes `configure` with autotools; nasm now keeps its tarball's times.) |
| `protal_static` against `~/isal/prefix`'s `libisal.a` and the system's `libzstd.a` | "statically linked"; the same SAM records |

**All three, locally** (callgrind, one thread, `scripts/cg_final.sh`, `results/isal/cg_final.log`; `cdce3c0` against
this round's build; the paired-end input one gzip member):

| instructions | `cdce3c0` | this round | |
|---|---|---|---|
| aligning 100k pairs, the whole run | 28.82 G | 28.42 G | −1.4% |
| – `RunPairedEnd` | 15.53 G | 15.29 G | −1.5% |
| – the seed sort | 528 M | 599 M | +13% (83% of the seeds shared here) |
| – the k-mer screen | 896 M | 612 M | −32% |
| – inflating the input (the reader's thread) | 1,124 M (zlib-ng) | 980 M (ISA-L) | −13% (ISA-L's gain is mostly in instructions per cycle) |
| PacBio 3 Mb, the whole run | 15.39 G | 15.34 G | −0.3% |
| – the seed sort; the screen | 121 M; 134 M | 138 M; 99 M | +15%; −26% |

**What is new in the tests**: `SharedSeeds.FindPairsSeesTheRunsOfAllSeedsSorted` (4,000 reads of 0-5,000 seeds, few or
many taxa and genes, repeats, the largest ids: the runs `FindPairs` sees equal those of all seeds sorted, the sorted
seeds are sorted and seeds of the read, few singletons among them), `SharedSeeds.OneSeedOrNoneHasNoneShared`,
`AlignmentScreen.TheReadsSharedCodesGiveTheDefinitionsAnswer` (108,000 screens: stretches of both strands at any offset,
long-read windows as reverse-strand stretches, Ns, ambiguity codes and lower case in the read, every k, with and without
AVX2, the reverse strand given or made: the answer of the screen's definition without its exits, and that of the screen
packing the candidate's read alone), `ThreadedGzStream.AHeaderSplitBetweenReadsIsRead` (a header with every optional
field split between the reader's first two reads at all 43 places, and 300 such members through a pipe written a few
bytes at a time); `ReadsGzipOfAnotherWriter` and the pipe test now get their gzip from zlib-ng.

**For the next cluster run**: the pe run should show the sequence reader's wait gone or smaller (ISA-L's inflate at
~1.6-1.8× zlib-ng's), the seed sort at about a third, the screen at about two thirds; HiFi the sort and screen likewise.

## The thirteenth cluster run: `230bf64` at r226

SLURM job 24012922 (the user's), `scripts/measure_performance.sh` at `230bf64` (this round's three changes) on node
`q512n10`, 32 threads, the same r226 v11 database and samples as the eleventh and twelfth runs;
[`results/cluster_v13/`](results/cluster_v13/) (copied unchanged from `local/Perf0.7.5_v13/`), compared with
`scripts/cmp_runs.sh` (`results/cluster_v13/cmp_v12_v13.txt`). The first paired-end run loaded the index cold from NFS
(369 s, the database tables alone 2.5 min) and is left out of the paired-end numbers, as before. Between `cdce3c0` and
`230bf64` only this round's commits touch protal's code.

Every count is the twelfth run's: the reads, the mates with an anchor, the candidates tried, refused by the screen and
aligned, the alignments made and records written, the seeds, the shared seeds, the dropped lookups and the anchors.

| | twelfth (`cdce3c0`) | thirteenth (`230bf64`) | |
|---|---|---|---|
| pe wall (the two warm runs) | 36.2, 37.0 s | **33.1, 33.6 s** | −9% |
| pe aligning | 28.4-29.0 s | 24.8-25.5 s | −13% |
| pe user time; instructions; cycles | 1,051 s; 7.33 T; 3.85 T | 954 s; 6.58 T; 3.48 T | −9%; −10%; −10% |
| pe pairs aligned per second | 1.75-1.78M | 1.99-2.04M | |
| HiFi wall | 13.7 s | **12.3 s** | −10% |
| HiFi aligning | 9.4 s | 7.9 s | −16% |
| HiFi user time; instructions; IPC | 403 s; 3.60 T; 2.47 | 356 s; 3.28 T; 2.56 | −12%; −9% |
| index load, peak memory | 2.66 s, 34.2 GB | 2.68 s, 34.2-34.5 GB | |

Per thread (medians over the runs, `stages.tsv`):

| stage, s per thread | pe twelfth | pe thirteenth | | HiFi twelfth | HiFi thirteenth | |
|---|---|---|---|---|---|---|
| seeding | 12.51 | 12.90 | +3% | 2.83 | 2.83 | |
| sorting seeds | 3.39 | **1.06** | −69% | 1.26 | **0.32** | −74% |
| pairing | 0.64 | 0.51 | −21% | 0.14 | 0.11 | −23% |
| alignment handler | 3.36 | 2.79 | −17% | 9.28 | 7.80 | −16% |
| – the k-mer screen | 1.91 | **1.36** | −29% | 0.99 | **0.44** | −55% |
| extending anchors | 2.05 | 2.09 | | 0.55 | 0.55 | |
| taking the k-mers | 1.94 | 1.95 | | | | |
| sequence reader | 2.87 | **1.19** | −59% | 0.51 | 0.51 | |
| output handler | 0.89 | 0.77 | −13% | 0.10 | 0.10 | |

(For HiFi the seeding, sorting, pairing and extending run inside the alignment handler's time.)

- **Each change did what was predicted.** The seed sort fell by 2.3 s per thread (pe, predicted ~2 s) and 0.94 s (HiFi,
  predicted ~0.8 s); the screen by 0.55 s (pe, predicted ~0.6 s) and 0.55 s (HiFi, predicted ~0.5 s).
- **The paired-end run is no longer held by its input.** The sequence reader's time per thread fell from 2.87 to 1.19 s,
  near the eleventh run's 1.08 s from before the input became the limit, while the run aligned 1.99-2.04M pairs/s, above
  the 1.8-2.0M pairs/s at which one zlib-ng thread per file capped it. ISA-L's inflate (1.6-1.8× zlib-ng's) leaves room
  to about 3M pairs/s or more.
- **Seeding is now over half of the paired-end aligning** (12.9 of 24.8 s per thread; its +3% is the nodes' spread: no
  seeding code changed). The k-mer lookups wait on memory (2.2G lookups of 70 cells in a 27 GB index); the levers are the
  GTDB-scale report's (more reads in flight per thread, a smaller index entry). After it come the alignment handler
  (2.8 s, the screen half of it), extending the anchors (2.1 s) and taking the k-mers (2.0 s).
- **The strain stage** took 3.7 s (5.6 s on the twelfth run's node, 3.0 s on the eleventh's): the twelfth run's was the
  node's NFS, as thought.

Since the eleventh run (`2093770`, before this report's changes): pe 42.2 → 33.4 s (−21%), HiFi 15.6 → 12.3 s (−21%),
instructions 8.76 → 6.58 T and 4.00 → 3.28 T; since the first r226 run of the GTDB-scale report, pe 131 → 33 s and HiFi
104 → 12 s.

## How it was run

In WSL (`~/perf6`; the implementation in `~/perf7`: `build7.sh ref|work`, `check7.sh`), scripts in `scripts/`:

```bash
bash scripts/build_head.sh            # git archive HEAD -> ~/perf6/head, Release
bash scripts/build_g.sh               # the same with -g -> ~/perf6/headg
bash scripts/runs.sh 1                # timed runs, rounds 1-3 (THREADS=6 / 1)
bash scripts/cg.sh                    # callgrind: aligning 100k pairs, profiling, PacBio, ONT, single-end
bash scripts/cgg.sh                   # callgrind of the -g build; annot.sh / callers.sh read them
bash scripts/cc.sh scripts/lookup_bench.cpp lookup_bench && ./lookup_bench 6 150000 6 ABC
bash scripts/cc.sh scripts/sort_bench.cpp sort_bench && ./sort_bench 200000 4
bash scripts/cc.sh scripts/genome_bench.cpp genome_bench && ./genome_bench 4000000 3
bash scripts/build_proto.sh seed scripts/flex_ties_avx2.patch   # also: screen scripts/screen_packed.patch
bash scripts/cmp.sh seed              # callgrind against HEAD, SAM records compared
bash scripts/drift.sh                 # pass 2's build (0e4c5a4) and HEAD on the v0.7.3 world
bash scripts/collect.sh <report dir>  # results/
```

The r226 numbers are arithmetic on the eleventh run's `runs.tsv`, `stages.tsv` and logs in
`../2026-10-04-performance-gtdb-scale/results_v11/`.

After the twelfth run (in `~/perf8`, ISA-L 2.32.1 and nasm 2.16.03 built into `~/isal/prefix` from their release
tarballs, `results/isal/build_nasm.sh` and `build_isal_cmake.sh`; this machine has no `libisal-dev`):

```bash
bash scripts/build8.sh ref            # HEAD (cdce3c0); `work`: items 1 and 3; `isal`: all three, -DCMAKE_PREFIX_PATH=~/isal/prefix
bash scripts/check8.sh                # items 1 and 3 against HEAD: unit and end-to-end tests, outputs, callgrind
bash scripts/check_isal.sh            # all three against items 1 and 3: tests, outputs, .sam.gz, a static protal
bash scripts/check_sfetch.sh          # lib/static-deps.cmake's nasm and ISA-L from the local tarballs
bash scripts/read_bench.sh            # BGZF written and gzip read through protal's own code, HEAD against this round
bash scripts/cg_final.sh              # callgrind, HEAD against this round
bash scripts/cmp_runs.sh results/cluster_v12 results/cluster_v13 twelfth thirteenth   # the cluster runs side by side
H=~/perf8/benchtree bash scripts/cc8.sh scripts/shared_sort_bench.cpp shared_sort_bench && ./shared_sort_bench 50000 1 0
bash results/isal/build.sh && bash results/isal/run_inflate.sh   # ISA-L against zlib-ng, alternated (bench_inflate.c)
```
