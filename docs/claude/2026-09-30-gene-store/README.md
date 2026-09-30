# The 2-bit gene store and the 32-byte `Gene`

2026-09-30. Branch `audit-fixes` on top of `4b21427` (Version 0.7.0), working tree, not committed.
Follow-up to [Where protal's memory goes](../2026-09-30-memory-profiling/README.md), strategies 1
and 3. Built and measured in WSL (Ubuntu 24.04, Core Ultra 7 258V, 24 GB, Release) from a copy on
the Linux file system; the old binary is `4b21427` built the same way.

## What changed

| | |
|---|---|
| `src/SequenceUtils/PackedSequence.h` (new) | the 2-bit coding; `packed::Pack` / `Unpack` with a table and with AVX2 (chosen at run time with `__builtin_cpu_supports`, like the syncmer scan); `packed::PackInto` for threads that fill pieces of one gene; `GeneSequence`, a decoded gene |
| `Gene` (`GenomeLoader.h`) | 112 → **32 bytes**: start byte or packed pointer in one 8-byte field, `uint32` id, 31-bit length and a loaded flag, four `uint32` counts; no `std::string`, no stream pointer, no start byte once loaded |
| `Genome` | holds the stream (`m_reader`) and the packed buffers of genes loaded one by one (`--preload_genomes_off`); `LoadGeneData` replaces `Gene::Load` |
| `GenomeLoader::LoadAllGenomes`, `GeneSink` | one zeroed arena (`calloc`, untouched until filled) of `Σ ⌈length/4⌉` bytes, every gene on a byte; the loading threads pack their pieces with `PackInto` (atomic or on the byte at either end of a piece, plain stores inside) |
| `Gene::Sequence()` | returns a `GeneSequence` (decoded on call, on the stack up to 4096 bases, else on the heap); ~20 call sites hold it in a variable instead of a `string_view`; `Sequence().length()` became `GetLength()` |
| `VariantHandler` | takes the `Gene` and decodes it once per alignment (`Reference()`), or a fixed `string_view` as in the tests; `StrainLevelContainer` passes its gene |
| `unique_kmers.tsv` | counts above 2^32 − 1 are refused as too large (they are stored in 32 bits; a gene has fewer than 2^20 bases) |
| tests | `tests/test_PackedSequence.cpp` (9 tests), the existing ones unchanged |
| docs | `docs/database-files.md` (Genes in memory), `docs/running.md` (Memory) |

**Coding** (your suggestion): A 0, C 1, G 2, T 3, four bases to a byte. `N` is stored as `A`, an
IUPAC code as the first base it stands for in the order A C G T (`R W M D H V` → `A`, `Y S B` → `C`,
`K` → `G`), lowercase as uppercase, anything else as `A`. No exception list.

**SIMD.** Packing takes 32 characters per step: a block of only A, C, G, T (either case) is coded by
one `pshufb` on the low nibble and two multiply-adds, a block with anything else goes through the
table; unpacking takes 32 bytes (128 bases) per step with four `pshufb` and the unpack/permute
instructions, the last <128 bases through a padded copy, so that nothing past the gene is read or
written. Both give the same bytes as the table (tested for every length to 700 and at 4095–65537,
all characters, both cases).

## Verification

| | |
|---|---|
| Unit tests | 199 of 199 (Release), and under ASan + UBSan (RelWithDebInfo) 199 of 199 |
| End to end | `tests/e2e/test_protal_e2e.py` on a fresh mini database built with the new binary: 112 of 112 (includes `--preload_genomes_off` with a raw reference, the compressed-reference refusal, and reruns) |
| Same results | `db900`, `mix` and `w900` (1M pairs each), 8 threads: the file lists of the outputs are equal, every output file but two is byte-identical (all profiles, statistics, strain tables); the SAM has the same 50,550 and 904,322 records in another order (threads), and `*.profile.gene.log` differs in the last printed digit of a float summed in another order (19.536563 / 19.536562) |

## Memory

`db900` (900 species, 101k genes), peak RSS of the whole run, 8 threads, 6 alternated runs each
(medians):

| | before | after |
|---|---:|---:|
| `mix` | 3.75 GB | 3.68 GB (−75 MB) |
| `w900` | 3.82 GB | 3.73 GB (−90 MB) |

The gene arena went from 111 MB to 28 MB and the gene tables from 11 MB to 3 MB; the key map
(3.2 GB) is unchanged, so the whole run changes by 2%. The part that scales with the database is
measured on a **synthetic reference of 8M genes** (8.4 Gbase, genes of 900–1199 bases, half of the
GTDB r226 size; `scripts/bench_load_genes.cpp` compiled against each version's `src/`, file in the
page cache, `reference.map` and `LoadAllGenomes` with 8 threads, a run that reads every gene once):

| | before | after |
|---|---:|---:|
| RSS after `reference.map` (the gene tables) | 0.94 GB | **0.35 GB** |
| RSS after `LoadAllGenomes` | 9.35 GB | **2.46 GB** |
| peak RSS | 10.0 GB | 3.1 GB |
| `reference.map`, 8 threads | 0.58–0.89 s | 0.38–0.39 s |
| `LoadAllGenomes`, 8 threads | 2.4–4.4 s | **0.95 s** |
| `LoadAllGenomes`, 1 thread | 3.5–12.0 s | 2.3 s |
| reading all genes once (`Sequence()` + one base) | 0.12–0.47 s (13–59 ns) | 0.27–0.37 s (34–46 ns) |

Doubled for 16.6M genes: the genes and their tables go from ~19.5 GB to **~5.1 GB**, so the full
r226 database from ~59 GB to **~45 GB** (−14 GB, −24%) with the key map (3.2 GB) and the index
values (35 GB) unchanged. Loading the genes is 2.5–4.6× faster at 8 threads: 4× fewer bytes to write, none to zero,
and packing runs at 13.6 GB/s of bases with AVX2 (3.5 GB/s with the table) against a copy.
Extrapolated, not measured at that size.

## Speed

Alignment stage (`Aligning reads took`) and whole run, `protal_avx2` of both versions alternated
6 times, 8 threads, `db900` (the laptop was on battery; single runs vary by ±10%):

| | before | after |
|---|---:|---:|
| `mix`, align / run (median) | 2.06 s / 3.10 s | 2.09 s / 3.10 s |
| `w900`, align / run (median) | 5.54 s / 8.89 s | 5.59 s / 8.98 s |

The alignment code decodes a gene per call (24 calls per pair in `w900`, 2.3 in `mix`, see the
memory report): **+1%, within the noise**. Decoding a whole gene, `scripts/bench_gene_decode.cpp`
(1 Gbase store, random genes, idle machine, ns per gene):

| Bases | byte copy from a byte arena (before) | table | AVX2 | `GeneSequence` (AVX2) |
|---:|---:|---:|---:|---:|
| 300 | 57 | 171 | 75 | 75 |
| 1000 | 162 | 245 | 172 | 171 |
| 3000 | 229 | 390 | 255 | 254 |

The packed store reads a quarter of the memory, which pays for the decoding: at 1000 bases it costs
10 ns more than copying from the big byte arena, the table path 80 ns more. So without AVX2 the
alignment stage would lose about 3% (estimated: 24 calls × 80 ns against ~64 µs of alignment CPU per pair in `w900`).

## Non-ACGT bases: what the coding changes

A real reference may well contain `N` or IUPAC codes; GTDB's marker genes were not counted
here. I tested with `db900` whose `reference.fna` got four ambiguity codes per gene (`N`, `Y`, `R`,
and a lowercase `k`, at random places, 1 per ~250 bases: far more than a real database), the index
rebuilt with each version, and `w900` mapped with each:

- **Index.** Key map and every value identical; 5,550 of 26.7M value words differ, all in the two
  uniqueness flag bits, and 4,582 of 101,081 rows of `unique_kmers.tsv` differ slightly. Cause: the
  build's uniqueness check for a k-mer whose core occurs once (`IndexedKmer`, `Build.h`) reads the
  k-mer back from the gene, which now gives the stored base where the index took `A` (forward) or
  the complement-coded 0 (reverse). `N` → `A` agrees with the index on the forward strand only.
  Only windows that overlap an ambiguous base are affected.
- **Alignments.** 195,042 of 893,137 SAM records differ at this density: with `N` in the reference
  the old aligner scored a mismatch (`X`) at the `N`, now the read is scored against `A` or `C`,
  matches a quarter of the time, and scores, CIGARs and MAPQ follow (for example `72M1X16M1X11M1X3M1X44M`
  became `72M1X16M1X15M1X44M`, MAPQ 110 → 111; some 111 → 68 where the true gene's score changed
  against its relatives'). No record went to `.err`. At a real density of a few ambiguous bases per
  thousand genes these differences shrink by the same factor.
- **Not byte-identical by design** for any gene with such bases: MSA and SNP reference rows show the
  stored base, and a SAM made by an older protal against such a gene can hold an `M` the stored
  base does not match (set aside into `.err`, with the existing warning).

If exact behaviour at ambiguous bases is wanted back, the smallest addition is a sorted list of
`(position, character)` per gene for the few genes that have them, applied in `GeneSequence` after
the unpack; the store stays 2-bit.

## Left as it was

- `Sequence()` decodes the whole gene on every call. A windowed call (`Sequence(pos, len)`) would
  cut it for the chain finder's per-seed checks; not needed at +1%.
- The build's `IndexedKmer` still reads genes from the store (see above).
- Not measured: the real r226 reference (ambiguity density), a cluster node, and the website's
  download sizes (memory figures on the website and in `docs/running.md` predate this).

## Reproducing

```bash
# build a copy on the Linux file system (see the performance report), then
ctest --test-dir build --output-on-failure                              # PackedSequence, GeneSequence, GenomeLoaderPacked
bash docs/claude/2026-09-30-memory-profiling/scripts/ab_genes.sh OLD/protal_avx2 NEW/protal_avx2 DB READS 8 6
g++ -O2 -std=c++20 -I src docs/claude/2026-09-30-gene-store/scripts/bench_gene_decode.cpp -o bench_gene_decode -pthread && ./bench_gene_decode
# bench_load_genes.cpp: compile against a version's src/ (include paths as in its build.ninja), then
./bench_load_X DIR 8000000 8 generate; ./bench_load_X DIR 8000000 8
```

## Follow-up: ambiguous bases are not indexed

The build discrepancy above (uniqueness flags of k-mers over an ambiguous base) is removed at its
source rather than patched: `--build` (`Build.h`, `DropAmbiguousKmers`) takes every k-mer whose
31-base window holds a base other than A, C, G, T out of a record's k-mers, in both passes that
count and place the values and in the uniqueness check. Nothing is indexed with an ambiguous base
read as `A` any more, so what `IndexedKmer` reads back from a gene (the stored bases) is exactly
the indexed k-mer, and a k-mer of `full_reference.fna` over an `N` no longer meets an `A` k-mer of
the index and flags it as shared. The two passes drop the same k-mers of the same record, which the
counting needs for the placing to fit. Reads are not changed (a read k-mer over an `N` is still
looked up with the `N` as `A`). The scan costs one pass over each record (vectorised), and only a
record that has an ambiguous base is scanned twice.

Effect on an index: a build of a reference without ambiguous bases is byte-identical to before (the
900-species world, `index.prx` and `unique_kmers.tsv`); with them the index loses the syncmers whose
window holds one (about 7 of the 31 windows over a base, at the syncmer density). Tests:
`AmbiguousKmers` (unit), `BuildAmbiguousBasesTest` (end to end: a gene with an `N`, `Y`, `R` or
lowercase `k` in a 3 kb gene loses between 1 and 31 k-mers, the other genes none).

After the change: 202 of 202 unit tests and 113 of 113 end-to-end tests. On the 900-species world
with four ambiguity codes per gene in `reference.fna`, the passes count 18.4M instead of 20.6M
syncmers (−10%), the index is the same for `-t 1`, `4`, `8` and `--serial_index_passes`, and the key
map and value slots are those of a build of the clean reference (byte-identical) only where no
window meets an ambiguous base.
