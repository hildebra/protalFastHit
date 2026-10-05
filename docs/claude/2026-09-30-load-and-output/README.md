# Loading the database and reads, writing the output

2026-09-30. Branch `performance` at `3a9cfaa` (after the work in
[2026-09-29-performance-profiling](../2026-09-29-performance-profiling/README.md)). Binary:
`protal_avx2`, Release. Same machine as that report: WSL2 Ubuntu 24.04 on an Intel Core Ultra 7
258V (4 performance and 4 low-power cores), 15 GB RAM, gcc 13, zlib 1.3, libdeflate 1.19, zstd 1.5.5.

**Data:** `db64` and `db900` and the read sets `mix` and `w900` (1M pairs each, 2×150) of the
earlier report; 50k- and 200k-pair subsets of them; synthetic `reference.map` and
`unique_kmers.tsv` files at GTDB r226 size (below). No real GTDB database is on this machine, so
GTDB-scale numbers are measured on synthetic files where possible and otherwise extrapolated from
per-value and per-gene costs, and marked as estimates.

**Timing caveat:** the laptop ran on battery (Windows throttles it; the CPU reported 2.2 GHz) with
several other sessions busy. Each rerun was slower than the one before (the last by 35–40%, e.g.
zlib inflate 431 → 281 MB/s), uniformly across the alternatives, so the ratios between them are
stable; absolute times are given as the range of all runs. Instruction counts (callgrind) do not
depend on this.

## Summary

For a run on GTDB, start-up is dominated by per-value and per-gene work, and a few single-threaded
steps become seconds each: parsing `reference.map` and `unique_kmers.tsv` (7–12 s, measured at
GTDB size) and sizing 16.6M gene strings one by one (est. 6–14 s). Every SAM starts with a 398 MB
`@SQ` header of all 16.6M genes, about as large as a typical sample's alignments, which is then
compressed and parsed again. On the read side, gzip inflation (zlib) and the reader lock bound
throughput again now that alignment is faster: at 8 threads `mix` takes 39–45% longer from
`.fq.gz` than from plain FASTQ. On the output side, compressing the SAM with pigz afterwards takes as
long as aligning it (3.7 s against 4.2 s at 8 threads for `w900`), because zlib level 6 compresses
at ~20 MB/s per core; libdeflate gives the same ratio at 3.6× the speed, 6× at level 4.

| # | Where | Measured | At GTDB scale | Opportunity | Gain | Results |
|---|---|---|---|---|---|---|
| 1 | `.sam.gz` output | pigz after alignment: 3.7 s at 8 threads, 20.7 s at 1, for `w900`'s 350 MB (alignment stage 4.2 s) | proportional to the SAM, plus the header (#2) | compress each thread's buffer with libdeflate as gzip members (BGZF) while aligning | −3.7 s → +~0.6 s CPU spread over the threads (level 6, same ratio); no uncompressed temp file | same SAM text; the `.gz` bytes differ |
| 2 | `@SQ` header | 101k lines, 2.2 MB (`db900`) | 16.6M lines, **398 MB per SAM**, 1.3–2.2 s to write, ~20 s CPU to compress at zlib 6, parsed and checked again by the profiler | list only the genes that have alignments (valid SAM; the profiler's check already works so) | per sample at GTDB: ~400 MB, several seconds of writing, compression and profiling | header lines differ |
| 3 | gene tables | `reference.map` + `unique_kmers.tsv`: 0.93 G instructions for 101k genes | **7–12 s single-threaded** (measured at 16.6M genes) | a `from_chars` parser over the whole buffer, in parallel chunks | est. to ~1 s | none |
| 4 | gene sequences | sized one by one (`PrepareSequence`) before the parallel fill | **est. 6–14 s single-threaded** (4.2 GB take 1.3–3.4 s; 8 threads do not help) | one arena for all genes (0.09–0.14 s per 4.2 GB touched by 8 threads) | est. −6–14 s, and 16.6M fewer allocations | none |
| 5 | FASTQ reader lock | 4.5k instructions per pair inside the lock (5% of `mix`'s 86k per pair) | the same per pair | take raw blocks of whole records under the lock (memchr), parse outside with a string_view parser | lock time ~20× shorter; lifts the ~19-thread ceiling for `mix` | none |
| 6 | gzip inflate | zlib 281–431 MB/s, 12.8k instructions per pair; `mix` 2.2–3.3 s `.gz` vs 1.6–2.3 s plain (8 threads) | the same per pair | libdeflate (3.05–3.10× zlib) for multi-member/BGZF input; a streaming inflater (ISA-L, zlib-ng; not measured) for ordinary `.fq.gz` | up to the plain-input times | none |
| 7 | profiling the SAM | 2.3–2.4 s plain, 3.3–3.5 s `.sam.gz`, one thread per sample (`w900`, 895k records); parsing 8.3k instructions per record | plus the 398 MB header (#2) | faster SAM parsing (no `getline`/split into strings), inflate in its own thread (`ThreadedGzIstream`), parallel BGZF blocks | est. −30–45% of profiling | none |
| 8 | index decode | ~3.6 G instructions fixed (key map) + ~270 per value, in parallel over frames | est. 1.2 T instructions: 200–320 s on 1 thread, ~35–55 s at 8 | a leaner `DecodeChunk` (154 of the 270 per value: per-key `memcpy` calls, temporary vectors); for array jobs, a shared memory-mapped index | est. −25–35% of index CPU; mmap: one copy per node, no per-job decode | none |

Numbers 1 and 2 matter for every map run (`protal_map_utils generate` names SAMs `.sam.gz`);
3, 4 and 8 for every run on the full database; 5 and 6 above ~8 threads. None of this changes a
profile.

## Loading the database

1000 pairs, all three runs (`scripts/startup.sh`):

| | `db64` 1 thread | `db64` 8 threads | `db900` 1 thread | `db900` 8 threads |
|---|---:|---:|---:|---:|
| Load Index | 1.70–1.87 s | 0.19–0.39 s | 1.85–2.97 s | 0.31–0.70 s |
| Preload genomes | 17–30 ms | 12–35 ms | 218–272 ms | 123–193 ms |
| whole run | 1.90–2.19 s | 0.36–0.71 s | 2.34–3.65 s | 0.65–1.28 s |
| peak RSS | 3.2 GB | 3.2 GB | 3.6 GB | 3.8 GB |

Where the instructions go, `db900`, 1 thread, 15.0 G in all (`scripts/cg_startup.sh`):

| Step | Instructions | Scales with | Threads |
|---|---:|---|---|
| index: `DecodeChunk` | 7.55 G (`db64`: 3.68 G) | key map (~3.6 G, fixed) + ~154 per value | frames in parallel |
| index: zstd | ~3.25 G (`db64`: ~0.33 G) | ~116 per value | frames in parallel |
| model (cPMML, patched) | 1.26 G | fixed | 1 |
| genes: zstd and copy (`ParallelReadFrames`) | 1.15 G | reference bytes (~11.4k per gene) | frames in parallel |
| genes: sizing and sorting | ~0.08 G | genes (plus zero-filling their bytes) | 1 |
| `reference.map` (`LoadPositionMap`) | 0.42 G | genes (~4.2k each) | 1 |
| `unique_kmers.tsv` | 0.50 G | genes (~5.0k each) | 1 |
| SAM header (`WriteSamHeader`) | 0.18 G | genes (~1.7k each), per sample | 1 |
| taxonomy | 0.07 G | taxa | 1 |
| profiler's header check (`CheckReference`) | 0.06 G | genes (~600 each), per sample | 1 |

### At GTDB r226 size

GTDB r226 has 136,646 bacterial and 6,968 archaeal species (GTDB's release statistics, not
checked against the database); protal uses 119 and 52 marker genes
([model-training.md](../../databases.md#the-presence-model)), so up to **16.6M genes** (~17 Gbp at ~1 kb each).
`db900` has 266 index values per gene; at that density the full database has ~4.4 G values
(35 GB), and 3.2 GB key map + 35 GB values + ~17 GB genes + 1.8 GB gene tables ≈ 57 GB, close to
the 59 GB that [running.md](../../running.md) gives. The mini database has fewer genes per species.

**Measured at that size** (`scripts/benchmarks.sh`, `bench_genes.cpp`: protal's own `GenomeLoader`
on synthetic files with 16,623,210 genes, a sparse `reference.fna`; plain files, three runs):

| | `db900`-sized (100,800 genes) | GTDB-sized (16,623,210 genes) |
|---|---:|---:|
| `reference.map` (532 MB) | 0.03–0.08 s | 3.4–6.8 s |
| `unique_kmers.tsv` (796 MB) | 0.02–0.06 s | 3.4–5.6 s |
| SAM header, written to a file | 0.008–0.025 s, 2.2 MB | 1.3–2.2 s, **398 MB** |
| memory of the gene tables | 20 MB | 1.8 GB |

In `database.protal` both text files also pass through one zstd stream first.

**Sizing the gene sequences** (`bench_prepare.cpp`): `LoadAllGenomes` gives every gene its
`std::string` of the right length (`assign(length, '\0')`) on one thread before the threads fill
them. 4M genes of 900–1199 bytes (4.2 GB): 1.34–3.44 s serially, 1.60–3.90 s with 8 threads
(allocation and page faults contend), 0.09–0.14 s for one arena of the same size touched by 8
threads. Scaled to 17.5 GB: **6–14 s single-threaded**, ~0.5 s as an arena.

**Estimated** (from the per-value and per-gene instruction counts at `db900`'s 3.7–5.9 G
instructions per second and its 4–7× speed-up from 1 to 8 threads): the index 200–320 s on one
thread, **35–55 s at 8 threads**; filling the genes 30–50 s on one, 5–9 s at 8. With the text files
and the gene sizing, start-up at 8 threads is about 55–90 s on this laptop from the page cache, of
which 13–26 s single-threaded (#3, #4). From network storage, reading the ~30 GB of compressed
database (the download is 33 GB) comes on top and probably dominates the index;
`scripts/db_compression_benchmark.sh load` on the cluster measures that.

### Opportunities

- **Gene tables (#3).** Both files are parsed line by line with `getline`, `Utils::split` into a
  vector of strings, `std::all_of(isdigit)` and `std::stoull` per field, then two `sparse_map`
  lookups per gene. Parsing the decompressed buffer with `std::from_chars` (and splitting it into
  line-aligned chunks parsed in parallel, inserted serially) keeps the checks and messages and
  should take well under a second at GTDB size.
- **Gene arena (#4).** One allocation for all sequences, in `reference.fna` order, that the reader
  fills directly (the parallel `GeneSink` already copies by offset); each `Gene` then holds a
  pointer and length. `Gene::Sequence()` returns `std::string const&` today and is used widely, so
  the change is mechanical but touches many call sites (a `std::string_view` accessor).
- **Index decode (#8).** At 8+ threads on shared storage the load is probably I/O-bound. The CPU
  part: `DecodeChunk` calls `memcpy` twice per key for 1–3 values and builds batches in temporary
  vectors; writing values in place when a batch has no flex cells would cut most of the 154
  instructions per value. For array jobs that run many protal processes on one node, a raw index
  memory-mapped from local disk or `/dev/shm` would be decoded once and shared (with the gene arena,
  the genes could be too); that is a larger change of the loading model.
- The fixed 3.2 GB / ~3.6 G instruction key map (1.7–1.9 s at 1 thread on `db64`) matters only for
  small databases.

## Reading reads

1M pairs, `db900`, 8 threads, `--no_profile`, two alternated runs (`scripts/reads_output.sh`):

| | `mix` `.fq.gz` | `mix` plain | `w900` `.fq.gz` | `w900` plain |
|---|---:|---:|---:|---:|
| Aligning reads | 2.2–3.3 s | 1.6–2.3 s | 4.2–5.4 s | 4.1–5.8 s |
| user CPU | 17.3–27.8 s | 13.8–19.9 s | 34.6–43.9 s | 32.8–46.3 s |
| reader time per thread (incl. waiting for the lock) | 0.59–0.93 s | 0.22–0.36 s | 0.33–0.37 s | 0.28–0.33 s |

(For comparison, at `995c4f1` 1M pairs of `mix` took 21.1 s with `.gz` and 8.5 s plain, at a
load of 11–14.) `w900`
spends its time aligning; `mix`, the realistic case, is again bound by its input.

Per pair of `mix`, 50k pairs, 1 thread (`scripts/cg_reader.sh`): the whole loop 86k instructions;
the batch copy under `omp critical(reader)` (`BufferedFastxReader::LoadBatch`: 32 records per file,
`getline` per line, appended to a `stringstream`) **4.5k**; parsing the records from it outside the
lock 3.2k; inflating (the two `ThreadedGzStreambuf` threads) 12.8k, of it crc32 2.5k.

- **The lock.** 5% of each pair's work is serial, so `mix` cannot use more than ~19 threads even
  from plain files, and queueing shows well before (at 8 threads the threads together spend about
  as long in the reader as the stage lasts). Under the lock only the boundaries are needed: find the
  128th newline in each file's inflated block with `memchr` and move those bytes out in one copy;
  parse them outside with a string_view parser instead of `getline` into strings twice. FASTA and
  multi-line FASTA keep their own path.
- **Inflate.** One thread per file keeps up with ~13 threads' worth of `mix` alignment
  (6.4k of 86k instructions per pair). Inflating R1 alone: `zcat` 1.48–1.53 s, `pigz -dc` 1.0 s.
  In memory (`bench_compress.cpp`, R1 of `mix`, 327 MB): zlib `inflate()` 281–431 MB/s,
  libdeflate 858–1334 MB/s (**3.05–3.10×**), crc32 alone 4.5–5.8 GB/s (w900: 3.05×). libdeflate
  only inflates whole gzip members, which suits BGZF or multi-member files (bgzip, some pipelines),
  and those can then also be inflated in parallel; an ordinary single-member `.fq.gz` needs a
  streaming inflater, ISA-L's igzip or zlib-ng (neither installed here, not measured).
- The two inflating threads come on top of `-t`: `-t 8` on 8 cores runs 10 busy threads with
  `.gz` input. Worth documenting, or counting them in `-t`.

## Writing the output

**The SAM.** `w900`, 1M pairs: 350 MB, 894,815 records, 101,083 header lines. Formatting records
(`Output handler`) takes 0.28–0.43 s per thread of the 4.2–5.8 s stage (7–8%), for `mix` under 2%.
Each thread fills a 16 MB buffer and writes it under `omp critical(sam_output)`; that part is cheap.

**Compressing it (#1).** A `.sam.gz` is written uncompressed and then compressed by `pigz -p <-t>`
(zlib level 6): 20.7 s with 1 thread, **3.7 s with 8**, to 102 MB (ratio 3.42), against 4.2 s for
aligning it. The SAM cut into blocks, each compressed alone, one thread (`bench_compress.cpp`,
first run; the later runs were uniformly 35–40% slower):

| Codec | 64 KB blocks: ratio | compress | decompress | 1 MB blocks: ratio | compress | decompress |
|---|---:|---:|---:|---:|---:|---:|
| zlib 1 | 2.71 | 120 MB/s | 327 MB/s | 2.73 | 121 MB/s | 338 MB/s |
| zlib 6 (pigz's level) | 3.29 | 21 MB/s | 404 MB/s | 3.42 | 19 MB/s | 434 MB/s |
| libdeflate 1 | 2.96 | 284 MB/s | 1005 MB/s | 3.00 | 297 MB/s | 1092 MB/s |
| libdeflate 4 | 3.20 | 129 MB/s | 1202 MB/s | 3.27 | 133 MB/s | 1372 MB/s |
| libdeflate 6 | 3.34 | 76 MB/s | 1123 MB/s | 3.46 | 73 MB/s | 1300 MB/s |
| zstd 1 | 3.10 | 319 MB/s | 1065 MB/s | 3.42 | 366 MB/s | 1579 MB/s |
| zstd 3 | 3.18 | 251 MB/s | 1205 MB/s | 3.73 | 316 MB/s | 1529 MB/s |

Each alignment thread could compress its buffer into gzip members before writing it: a
multi-member gzip file, which `zcat`, `gzip -d`, zlib's `gzread` (protal's own reader) and htslib
read. With 64 KB blocks and the BGZF extra field it is BGZF, htslib's own format, whose blocks a
reader can inflate in parallel. At libdeflate level 6 that is 350 MB / 76 MB/s ≈ 4.6 s of CPU
spread over the alignment threads (about +0.6 s at 8 threads) instead of 20.7 s of CPU and 3.7 s
afterwards, at pigz's ratio; level 4 is 1.7× faster again for a 4% larger file. zstd compresses
better and faster, but `.sam.zst` is read by fewer tools. It also removes the uncompressed
temporary SAM and reading it back. libdeflate is small (MIT) and packaged; it would be a new
dependency, or vendored like WFA2.

**The header (#2).** `WriteSamHeader` lists every gene of the database. At GTDB size that is
16.6M lines and 398 MB in every SAM (measured, above), written in 1.3–2.2 s before the first
record; pigz then compresses it (at ~20 MB/s per core, ~20 s of CPU), and the profiler reads and
checks every line again. A sample with 10M pairs of which 3–5% align has SAM records of a similar
size (0.6–1M records of ~390 bytes). SAM needs `@SQ` lines only for references that records use, and
`Profiler::CheckReference` checks each `@SQ` line that is there against the database, so it
accepts such a header unchanged. The header has to follow the records then: collect the genes used
per thread, write the records to the `.partial` file, and at the end write the header and append
the records (`copy_file_range`, or as separate gzip members). Tools that combine SAMs of several
samples would see different `@SQ` sets per file; that was not tested.

**Reading it back (#7).** The profiler reads each SAM on one thread (samples in parallel):
2.29–2.41 s for `w900`'s plain SAM, 3.34–3.47 s for its `.sam.gz` (the extra second is zlib
inflating in the same thread, through gzstream). callgrind of `--profile_only` on 176k records
(`scripts/cg_profile.sh`): parsing (`SamReader::Next`) 8.3k instructions per record, of it
`LineSplitter::Split` into strings, `strtol`/`strtoul` and `SamEntry` construction; adding the
record to its taxon (`Taxon::AddSam`) 8.1k; writing the profile 1.29 G in all. Parsing with
string views and `from_chars` into reused entries, inflating in a thread of its own
(`ThreadedGzIstream`, as for reads) and, with BGZF, in parallel blocks would cut 30–45% of it. The
earlier report's larger options stay: profile a sample while the next one aligns, or hand the
alignments over in memory.

## Order

1. In-process libdeflate/BGZF compression of `.sam.gz` (#1): the largest gain for map runs today,
   self-contained in the output handler and `Compressor`.
2. `@SQ` lines for used genes only (#2): the largest per-sample saving at GTDB scale.
3. Gene tables with `from_chars` (#3) and the gene arena (#4): 13–26 s of single-threaded start-up
   at GTDB scale; #3 is small and local, #4 touches many call sites.
4. The reader lock (#5), then a faster inflater (#6): scaling beyond 8 threads.
5. SAM parsing in the profiler (#7).
6. `DecodeChunk` (#8) once the cluster shows whether loading is CPU- or I/O-bound.

## Gaps

- No real GTDB database: the index, gene-fill and memory numbers at that size are extrapolations;
  run `scripts/startup.sh` and `scripts/db_compression_benchmark.sh load` on the cluster.
- Laptop on battery with other sessions busy: absolute times are ranges, not benchmarks.
- ISA-L and zlib-ng were not available to measure; libdeflate was measured in memory, not inside
  protal.
- The effect of #1 and #2 on a full map run (many samples, strain MSAs) was not run.

## Follow-up: SAM formats and the trimmed header (#1, #2), implemented

2026-09-30, on `performance` after `3a9cfaa`. Asked whether zstd, which protal already links, would
do instead of gzip: yes, as an output format, so no library was added; gzip stays, for tools that
read only it.

- `src/IO/SamFile.h`. The SAM's name chooses the format: `.sam.zst` (zstd level 3, 1 MB frames,
  zstd's seekable format with a seek table at the end and a skippable marker frame at the start),
  `.sam.gz` (BGZF: 64 KB gzip members with the `BC` field and the end-of-file block, zlib level 6,
  as pigz), or plain. The output handlers hand each full 16 MB buffer, with the genes its records
  name, to `SamOutput`; the thread compresses it and only the append is serialised (replacing
  `omp critical(sam_output)`). `FinishSamFile` no longer runs pigz, and protal no longer needs pigz.
- The records collect in `<sam>.records.partial`; `Finish` writes `<sam>.partial` as the header
  (`@HD`, `@SQ` of the genes the records name, sorted, and the read-type `@CO`), then the records
  (`copy_file_range`), then the format's end, and the file is renamed as before.
  `--full_sam_header` (developer option) writes every gene first, as before, with no copy.
- Reading (`SamInput`): the profiler and `--profile_only`'s read-type check take any of the three
  by content. Besides zlib's and zstd's own checks, a BGZF file must end with its end-of-file block
  and a `.sam.zst` with protal's marker must end with a valid seek table, so a file cut exactly at
  a block or frame boundary is an error too, not a sample with fewer reads.
- `protal_map_utils generate` and `merge --use-sampleid` take `--zstd` (`.sam.zst` names) next to
  `--nogzip`.

**Same results** (`scripts/sam_validate.sh`, 200k pairs of `w900` and `mix`, against `3a9cfaa`):
with `--full_sam_header` every output is byte-identical at 1 thread; by default every output but
the SAM is, and the SAM's records are too. Its header lists exactly the genes that records name
(RNAME and RNEXT), all taken from the old header: 1,449 and 993 of 101,081 genes, 2.17 MB of header
→ 31 KB and 21 KB. `.sam.gz` and `.sam.zst` pass `gzip -t` and `zstd -t` and decompress to the plain
SAM byte for byte at 1 thread, to the same records at 8, with identical profiles; a rerun finds them
and profiles them alike; no temporary file is left. Unit tests 158/158, also under ASan/UBSan, and
TSan on the `SamFile` tests; e2e 92/92 on a fresh mini database, with new tests for the formats,
the header, reruns, `--full_sam_header` and `protal_map_utils --zstd`.

**Speed**, 1M pairs, `db900`, 8 threads, `--no_profile`, map runs, two alternated rounds
(`scripts/sam_speed.sh`; on battery, load 4–9 from other sessions):

| `w900` | wall | user CPU | aligning (incl. compressing) | SAM |
|---|---:|---:|---:|---:|
| `3a9cfaa`, `.sam.gz` (pigz afterwards) | 14.7–16.5 s | 111–117 s | 7.1–7.9 s | 102.2 MB |
| `.sam.gz` (BGZF, zlib 6, while aligning) | 14.5–15.9 s | 108–115 s | 13.2–14.3 s | 106.0 MB |
| `.sam.zst` | 9.9–11.3 s | 73–84 s | 8.7–9.9 s | 93.3 MB |
| `.sam` | 11.4–13.1 s | 82–96 s | 9.8–11.1 s | 347.8 MB |

`.sam.gz` did not get faster: zlib at level 6 costs as much CPU in the alignment threads as it did
in pigz, and at 8 busy threads there is nothing to overlap it with (the report's libdeflate estimate
for #1 assumed a 3.6× faster deflater). It no longer needs pigz, an uncompressed temporary SAM or a
second pass, and is 3.7% larger (64 KB blocks). `.sam.zst` takes 22–38% less wall time and ~30%
less CPU than `.sam.gz`, for a 12% smaller file (9% smaller than pigz's), and is even faster than
writing the 3.7× larger plain SAM here. For `mix` (a 17.5 MB SAM) the four were within noise
(4.6–6.5 s). Profiling the `w900` SAM (1 thread): plain 3.6–3.8 s, `.sam.zst` 3.8–4.1 s,
`.sam.gz` 5.3–5.7 s (old or new). The header saving is small at `db900` (2.2 MB per SAM; with 200k
pairs of `mix` the old `.sam.gz` was 27% larger than the new one because of it) and 398 MB per SAM
at GTDB size.

Not done: libdeflate for `.sam.gz` (3.6× less CPU at the same ratio, measured above) would make
gzip output about as cheap as zstd, at the cost of a new dependency.

## Follow-up: zstd by default, libdeflate for gzip, and the MSAs

2026-09-30, on `performance` after `24aa261`.

**zstd by default (`3249337`).** The SAM is an intermediate file, so the names protal picks (`-1/-2`,
or a map without a `SAM` column) now end in `.sam.zst`; `--sam_format gz|sam` picks another format.
`protal_map_utils` (`--gzip`, `--nogzip`) and `simulate_metagenomes --protal_metafile` follow.

**libdeflate as the gzip handler.** A new dependency (`libdeflate-dev`, conda `libdeflate`):
- Writing: `.sam.gz` BGZF blocks are deflated with libdeflate instead of zlib (`src/IO/Bgzf.h`),
  and the simulator's reads are compressed in process as BGZF (`bgzf::CompressFile`, parallel,
  the same bytes for any thread count) instead of by pigz, which protal no longer needs anywhere
  (`--pigz_path` is accepted and ignored).
- Reading: `ThreadedGzStreambuf` (reads, SAMs, genomes) inflates a BGZF file block by block with
  libdeflate and checks each block's CRC and the end-of-file block; any other gzip file, e.g. the
  single-member `.fq.gz` of sequencers, still with zlib, since libdeflate cannot stream (a member
  must be in memory whole). The profiler reads `.sam.gz` through it too (before: gzstream).
- Instructions, callgrind, 1 thread (`scripts/gzip_cg.sh`; wall-clock numbers were swamped by
  load 11 on battery):

| | total | in the inflating threads | in `SamOutput::Write` |
|---|---:|---:|---:|
| read 50k pairs of `mix`, single-member gzip (zlib) | 18.21 G | 640 M | |
| the same reads as BGZF (libdeflate) | 17.91 G | 337 M | |
| the same reads plain | 17.57 G | | |
| write 20k pairs of `w900` as `.sam.gz` (libdeflate) | 20.36 G | | 923 M |
| as `.sam.zst` | 19.68 G | | 239 M |
| as `.sam` | 19.44 G | | |

User CPU for 1M pairs of `w900` at 8 threads (noisy): `.sam.gz` 79–81 s now against 108–115 s with
zlib; `.sam.zst` 71–73 s. zstd stays the cheaper and smaller format (93 MB against 104 MB).

**The strain MSAs, `.gz` or `.zst`? (evaluated, not changed.)** Each MSA holds one row per sample,
the species' concatenated marker genes (~110–130 kb). Two runs (`scripts/msa_formats.sh`, 1M pairs
split into samples of the same species): 3 raw MSAs of 3 rows from `w900`, and 2 raw MSAs of 5 rows
plus one qcmsa-filtered MSA from `s64` in 6 samples:

| compression ratio | `w900` raw (1.17 MB) | `s64` raw (1.04 MB) | `s64` filtered (0.50 MB) |
|---|---:|---:|---:|
| gzip -6 | 4.0 | 3.7 | 3.2 |
| BGZF, libdeflate 6 | 4.3 | 3.8 | 3.4 |
| zstd -3 | 7.6 | 10.8 | 13.7 |
| zstd -19 | 9.3 | 12.3 | 15.8 |

gzip's 32 KB window is shorter than one row, so it compresses each row on its own (3–4×, whatever
the number of samples); zstd's window (2 MB at level 3) spans many rows and stores each further
sample's row mostly as matches to the others, so its ratio grows with the samples (7.6× at 3 rows,
10.8× at 5). Compressing is cheap either way (zstd -3 under 10 ms here, a few ms per MB), and
writing the MSAs is a small part of a run (0.4–0.9 s here). What reads them:
- `qcmsa.py` (protal runs it on each raw MSA) reads plain and `.gz` (Python's `gzip`), not `.zst`:
  the standard library has zstd only from Python 3.14 (`compression.zstd`); before, the `zstd` CLI
  (which the conda package brings) or the `zstandard` module would be needed. It read a raw MSA as
  fast from `.gz` as plain (0.44 and 0.48 s).
- `scripts/strain_test/strain_report.py` opens them as plain text.
- The filtered `.msa.fna` goes to phylogenetics tools (RAxML-NG, IQ-TREE, FastTree), which expect
  plain FASTA; some read gzip, few or none zstd.

So: keep the filtered `.msa.fna` plain; the raw MSA, an intermediate between protal and qcmsa, could
be `.raw.msa.fna.zst` (with qcmsa reading it through `compression.zstd` or `zstd -dc`, and
`strain_report.py` likewise). In absolute terms it matters for large cohorts only: 1,000 samples ×
200 species × 120 kb are ~24 GB of raw MSAs, ~6–8 GB as gzip, and likely a small fraction of that
as zstd (the ratio at 1,000 rows was not measured).

## Follow-up: #3–#5 implemented

2026-09-30, on `performance` after `401694b`.

**#3, the gene tables.** `reference.map` and `unique_kmers.tsv` are read in pieces of ~64 MB of
whole lines; each piece is cut into chunks parsed in parallel (`std::from_chars`, no `getline` or
strings per field), and the rows are added to the genomes in file order, so the first problem in
the file is still the one reported, with its line (tested with problems in late chunks, and a
duplicate gene against a malformed line in either order). `-t` threads (`ProtalDB`). A number too
large for 64 bits is now an error message instead of an uncaught exception. A first version read
each table whole and kept all rows: 2.5× less user CPU, but the 1.7 GB of short-lived buffers cost
4 s more in page faults (WSL2); pieces with reused buffers avoid that.

At GTDB size (`bench_gene_tables.cpp`, 16.6M synthetic genes, both tables; load 9–15 on battery,
so wall times are inflated):

| | wall | user + system CPU | peak RSS |
|---|---:|---:|---:|
| before (map, unique k-mers; the SAM header, ~2–3 s, also written) | 23–25 s | ~13 s without the header | 1.8 GB |
| after, 1 thread | 7.8–8.9 s | 6.2–7.4 s | 1.95 GB |
| after, 8 threads | 4.0–6.3 s | 7.2–7.5 s | 1.97 GB |

callgrind at 100k genes: 1.1k instructions per gene for `reference.map` (was 4.2k) and 1.5k for
`unique_kmers.tsv` (was 5.0k). What remains is the serial, memory-bound part: adding 16.6M genes
to 143k genomes (the gene objects alone are 1.6 GB). The report's estimate of "under a second" was
too optimistic for it.

**#4, the gene arena.** `LoadAllGenomes` no longer sizes and zero-fills one `std::string` per gene
before the parallel fill: it allocates one uninitialized buffer for all genes it loads, advised for
huge pages when over 64 MB, and points each gene at its bytes in `reference.fna` order; the reader
threads then write the bytes (and take the page faults) in parallel. `Gene::Sequence()` returns a
`std::string_view`, so the aligners, SNP and variant code, and the MSA take views (the variant
handler held a `std::string const&` to the gene; it now holds the view). A gene loaded on its own
(`Gene::Load`, a genome not preloaded) still gets its own string. Tested: a genome loaded gene by
gene, then `LoadAllGenomes` for the rest, at 1 and 3 threads.

`bench_gene_arena.cpp` at 2M synthetic genes (2.1 GB, the sparse reference in the page cache, so
this is the memory side only; load 11–12):

| `LoadAllGenomes` | wall | user | system | peak RSS |
|---|---:|---:|---:|---:|
| before, 1 thread | 4.1–4.8 s | 0.78 s | 1.8–1.9 s | 2.33 GB |
| after, 1 thread | 2.1–3.0 s | 0.55–0.60 s | 1.3–1.9 s | 2.34 GB |
| before, 8 threads | 2.8–3.5 s | 0.96–0.99 s | 2.2–2.7 s | 2.82 GB |
| after, 8 threads | 1.3–1.4 s | 0.66 s | 2.0–2.1 s | 2.80–2.84 GB |

The system time (page faults, copying from the page cache) stays; it is now spread over the reader
threads instead of taken by one thread first. Scaled linearly to GTDB's 17.5 GB (not measured; the
machine has 15 GB): 12–17 s less wall time at 8 threads. Memory is unchanged: the arena saves only
the per-string allocation overhead (~16 bytes a gene), and each gene keeps its empty `std::string`.
Every output of 200k pairs of `w900` and `mix`, at 1 and 8 threads, as before.

**#5, the reader lock.** For FASTQ, `LoadBatch` (under `omp critical(reader)`) now takes the batch's
4 × 32 lines as they are: from a `ThreadedGzStreambuf`, protal's read input for plain and gzipped
files alike, `TakeLines` finds them with `memchr` in the inflated block and copies them in one
`append` per block; from any other stream (tests, `std::istringstream`) line by line with `getline`.
`NextFastq` then parses the batch outside the lock with `std::string_view`s instead of `getline`
into a `std::stringstream`: the same records, ends (an empty header line, a record cut short) and
error message as `ReadNextSequence`, which FASTA and `SeqReader` (`LoadBlock`, used by the build)
keep. Tested: whole lines across 1 MB blocks (a line longer than a block, empty and CRLF lines, a
last line without `'\n'`) against `getline`; records with descriptions, tabs, CRLF, a 1.5 MB read
and a last record without `'\n'` from `std::istringstream`, a plain and a gzipped file, against
`SeqReader`'s old parser; a malformed header through both paths.

callgrind, 50k pairs of `mix` (`.fq.gz`), 1 thread, per pair:

| | before | after |
|---|---:|---:|
| `LoadBatch` (under the lock) | 4.5k | 1.08k |
| parsing (outside) | 3.2k | 1.5k |
| the whole run (start-up included) | 355k | 350k |

Of the 1.08k, 650 are the one copy, which glibc does with `rep movsb` and callgrind counts per byte
(so less in time than in instructions), 330 are `memchr` (one call per line) and 90 the rest.

The reader alone (`bench_reader.cpp`: threads take pairs and do nothing else; 1M pairs of `mix`;
load 10–17):

| pairs/s | 1 thread | 4 threads | 8 threads |
|---|---:|---:|---:|
| plain, before | 0.97–1.18M | 0.81–0.96M | 0.68–0.93M |
| plain, after | 2.07–2.24M | 1.27–1.87M | 1.18–1.64M |
| `.fq.gz`, before | 232–274k | 225–235k | 275k |
| `.fq.gz`, after | 227–266k | 257–262k | 307–324k |

From plain files the ceiling rises 1.7–2×. From gzip files the reader is now bound by zlib's
inflate (one thread per file, ~270k pairs/s here), which is #6: 1M pairs of `mix` align in 3–5 s
at 8 threads (200–330k pairs/s), so `mix` from an ordinary `.fq.gz` waits for its input; from BGZF
input (bgzip, protal's simulator) libdeflate inflates ~3× faster (above). The alignment stage on
1M pairs of `mix` at 8 threads was 4.9–5.0 s → 3.5–4.8 s (`.fq.gz`) and 4.8–5.4 s → 3.0–4.7 s
(plain), too noisy on this machine for more than "not slower". Every output of 200k pairs of `w900`
and `mix`, at 1 and 8 threads, as before.

## Follow-up: zlib-ng replaces zlib (#6)

2026-09-30, on `performance` after `d4105e8`. zlib-ng rather than ISA-L, for portability (ISA-L's
fast paths are assembly for x86-64 and aarch64, built with nasm) and the static build; the only
gzip path, no zlib fallback.

**What changed.** zlib-ng 2.3.3 is vendored in `lib/zlib-ng` (the release tarball, sha256
`f9c65aa9…07d1`, as conda-forge packages it, without `test/` and `doc/`: 228 files, 2.5 MB) and
built by `lib/zlib-ng.cmake` with upstream's own CMakeLists.txt: a static `libz-ng.a` with the
native API (`zlib-ng.h`, `zng_*`), runtime CPU detection (SSE2 to AVX-512 variants chosen at run
time, no `-march=native`), no tests (they would fetch googletest), no install rules. Upstream's CMake
target does not pass `WITH_GZFILEOP` on to its users, which `zlib-ng.h` needs for the `gzFile`
functions; `lib/zlib-ng.cmake` adds it. Everything that used zlib now uses zlib-ng: the gzip input
that is not BGZF (`ThreadedGzStreambuf`), `gzstream` (`igzstream`/`ogzstream`) and the tests. zlib
is gone from the build (`find_package(ZLIB)`, `-lz`), from CI (`zlib1g-dev`) and from the conda
recipe; zlib-ng needs nothing new, as it is built with protal and linked in. With the native API a
leftover `#include <zlib.h>` or `gzopen()` would fail to compile or link instead of silently using
a system zlib. Updating it: replace `lib/zlib-ng` with a new release (less `test/`, `doc/`) and
update the version and sha256 in `lib/zlib-ng.cmake`.

**Checks.** A fresh configure and build as CI does: no compile unit reads the system `zlib.h`
(`ninja -t deps`), no link command has `-lz`; `ldd` lists no `libz` for `protal`, `protal_avx2`,
`simulate_metagenomes` or the tests; the static binaries link, and the static `protal` gives every
output of 200k pairs of `w900` (`.fq.gz`, 8 threads) as before. Unit tests 167/167 (new: gzip that
libdeflate wrote, one member and several, read through zlib-ng), also in CI's sanitizer
configuration (Debug, ASan + UBSan, zlib-ng compiled with them too); e2e 93/93; every output of
200k pairs of `w900` and `mix`, at 1 and 8 threads, as before; the SAM records of 1M pairs of `mix`
identical.

**Speed.** This time the laptop was on mains power and idle (load 0.05–2.7): everything ran 3–6×
faster than in the sections above, so compare within this section only.

| inflate, 1 thread (`bench_inflate.cpp`) | `mix` R1 (99 → 327 MB) | `w900` R1 (97 → 336 MB) |
|---|---:|---:|
| zlib 1.3, `gzread` | 393–403 MB/s | 419–421 MB/s |
| zlib-ng 2.3.3, `zng_gzread` | 611–668 MB/s | 678–697 MB/s |
| libdeflate 1.19, the whole file in memory (BGZF's library) | 1247–1263 MB/s | 1259–1320 MB/s |

zlib-ng inflates 1.5–1.7× as fast as zlib; libdeflate, which protal uses for BGZF input, is still
~1.9× faster, so BGZF input keeps its advantage.

| the reader alone, `.fq.gz`, pairs/s (`bench_reader.cpp`) | 1 thread | 4 threads | 8 threads |
|---|---:|---:|---:|
| zlib (`d4105e8`) | 1.22M | 1.10–1.13M | 0.89–0.96M |
| zlib-ng | 1.93–1.94M | 1.70–1.81M | 1.39–1.46M |

From plain files the reader is unchanged (6.3–6.6M pairs/s at 1 thread).

The alignment stage, 1M pairs of `mix` from `.fq.gz`, 8 threads, `--no_profile`, six alternated
runs per build:

| | zlib (`d4105e8`) | zlib-ng |
|---|---:|---:|
| Aligning reads | 1.83–2.76 s | 1.49–2.22 s |
| wall (start-up included) | 2.39–3.60 s | 1.98–2.72 s |
| user CPU | 15.8–24.2 s | 13.4–18.7 s |
| reader time per thread (incl. waiting for the lock) | 0.42–0.67 s | 0.16–0.29 s |

Most of the gain is waiting that is gone: with zlib the lock holder waited for inflated bytes while
the others waited for the lock (likely spinning for part of it, as user CPU fell as well). An ordinary `.fq.gz` still inflates
at about half BGZF's speed, one thread per file; going further would need parallel inflating of
single-member gzip, which no library here does.

## Reproducing

The scripts are in [`scripts/`](scripts/) (settings in `env.sh`: `PERF_DIR`, `BIN`, `PROTAL_SRC`,
`OUT`); the data are those of the earlier report, made by its `prep_db.sh` and `prep_reads.sh`.

```bash
export PERF_DIR=~/protal-perf BIN=~/protal-perf/an/build/protal_avx2 PROTAL_SRC=~/protal-perf/an/src
bash scripts/startup.sh $PERF_DIR/reads/mix $PERF_DIR/db64 $PERF_DIR/db900
bash scripts/cg_startup.sh $PERF_DIR/db64 $PERF_DIR/db900
bash scripts/cg_reader.sh $PERF_DIR/db900 $PERF_DIR/reads/mix_plain
bash scripts/reads_output.sh $PERF_DIR/db900 mix $PERF_DIR/reads/mix $PERF_DIR/reads/mix_plain \
    w900 $PERF_DIR/reads/w900/reads $PERF_DIR/reads/w900_plain
bash scripts/cg_profile.sh $PERF_DIR/db900 SAM_OF_200K_PAIRS
bash scripts/benchmarks.sh $PERF_DIR/reads/mix/mix_R1.fq.gz $PERF_DIR/io/o_w900_gz/w900_gz.sam
# follow-up: the SAM formats, against the build before them (OLD) on 200k-pair sets, then speed
bash scripts/sam_validate.sh OLD/protal_avx2 NEW/protal_avx2 $PERF_DIR/db900 READS_200K_W900 w900 READS_200K_MIX mix
bash scripts/sam_speed.sh OLD/protal_avx2 NEW/protal_avx2 $PERF_DIR/db900 w900 W900_R1.fq.gz W900_R2.fq.gz \
    mix $PERF_DIR/reads/mix/mix_R1.fq.gz $PERF_DIR/reads/mix/mix_R2.fq.gz
# follow-up #4: LoadAllGenomes before and after the gene arena (two checkouts: e825cbf and after)
bash scripts/gene_arena.sh OLD_CHECKOUT NEW_CHECKOUT 2000000
# follow-up #5: the reader alone before and after (checkouts 94f6a12 and after), plain and gzipped
bash scripts/reader_lock.sh OLD_CHECKOUT NEW_CHECKOUT $PERF_DIR/reads/mix_plain/mix_R{1,2}.fq \
    $PERF_DIR/reads/mix/mix_R1.fq.gz $PERF_DIR/reads/mix/mix_R2.fq.gz
# follow-up #6: zlib, zlib-ng and libdeflate inflating; the reader before (d4105e8) and after zlib-ng
ZNG_BUILD=NEW_BUILD/zlib-ng bash scripts/inflate_libs.sh $PERF_DIR/reads/mix/mix_R1.fq.gz W900_R1.fq.gz
ZNG_BUILD=NEW_BUILD/zlib-ng bash scripts/reader_lock.sh OLD_CHECKOUT NEW_CHECKOUT \
    $PERF_DIR/reads/mix_plain/mix_R{1,2}.fq $PERF_DIR/reads/mix/mix_R1.fq.gz $PERF_DIR/reads/mix/mix_R2.fq.gz
```
