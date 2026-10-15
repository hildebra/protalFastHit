# The F1 audit's remaining findings: B2 measured and fixed, L2-L13 fixed or documented

2026-10-15. **The question** (the user's): fix B2 and the lower-severity findings L2-L13 of the
[F1 audit](../2026-10-15-f1-audit-since-0.6.0a.md), which its section 6 left open after 0.8.5.

**How it was done.** Branch `f1-lower-findings` from `audit-fixes` at `d6ad7b9` (protal 0.8.5). Each finding was
checked again against 0.8.5's code by a reviewer (three in parallel, read-only), then fixed, documented or left with
a reason; B2 was measured first. L2 was implemented on a branch of its own and measured there
([L2-gene-end-clips.md](L2-gene-end-clips.md)).

**Short answer.**

| Finding | Still true at 0.8.5 | Done |
|---|---|---|
| B2 the gates drop short divergent reads | yes; measured (below): at C 4, 9-12% of 75-100 bp single-end reads 5% from their reference, 3% of 150 bp, 0.1% of 2x150 | **fixed**: `--sparse_first_min_picks` (a short read's stride halved until it has 6 picks; at 5% 1.4-2.0% of those reads dropped instead of 9-12%, for 3x the gate's lookups of a 75 bp read), recorded by `--build`, off for databases built before |
| L2 clipped bases at a gene end score nothing against a congener's covering gene | yes | **option** `--gene_end_clips credit`, recorded with the database: the default of new databases since the integration, `off` for older ones: the clip beyond the least-clipped candidate's counts at the floor rate ([L2-gene-end-clips.md](L2-gene-end-clips.md)) |
| L3 adaptive candidates: anchor order, split genera never tried | partly: the order is by exact-match length (taxid only breaks ties); split genera (`Genus_A`) are other genera | **counted**: a log line per sample of the reads with more eligible congeners than `--adaptive_candidates` tried; split genera left (≤1.3% of FN taxa's own fragments landed outside their genus at r226 v15) |
| L4 exact ties go to the lower internal taxid | yes, for 16 candidates or fewer (libstdc++'s stable insertion sort), introsort's order beyond; MAPQ 0, but the tie's record counts in `fragments_all` and a few rates | **documented** (`Classify.h` SortAlignmentPairs3) |
| L5 a long read's `fragments` count once per gene segment, not once | yes | **documented** (features.md, the comments; training and runs count alike) |
| L6 a `unique_cores.bin` that does not fit falls back gene by gene, silently | yes | **fixed**: a table whose rows are not the database's genes is not used (whole run from the anchors, a warning); the genes that do not fit are counted over a sample's threads and warned about |
| L7 single-value cores of an index from before `78f7a92` count as unique | yes | **fixed**: they count only where the index checked them; the load warns for such an index |
| L8 a reduced (`--genes`) database leaves the other genes out of the uniqueness pass | yes | **documented** (databases.md; the comments' rationale corrected) |
| L9 the AVX-512 `--flex_indels` kernel never ran on AVX-512 hardware | yes (this machine has no AVX-512) | **documented**: `PROTAL_SIMD=avx2` with `--flex_indels` on such CPUs until checked |
| L10 an N in a read is looked up as A or T | yes | **fixed**: a read's k-mers over an ambiguous base are not looked up (as the index has none) |
| L11 mate files in another order only warn | yes | **fixed**: a sample whose pairs are mostly named apart fails; `r` and `r/2` now match |
| L12 `-march=native` turns on WFA2's AVX2 kernels | yes, with such flags | **fixed**: WFA2 is compiled at the baseline after any `-march` in the flags |
| L13 `--depth_identity_margin`'s help | yes | **fixed** (help, running.md) |

What changes for default runs on a database built before: L7 (only on an index from before 2026-09-26), L10 (reads
with an N: 1-2% of simulated 150 bp reads carry one), L11 (only misordered input, which now fails). B2 and L2 act only
on databases that record them (B2: every database `--build` makes from now on; L2: none unless asked).

## Setup

- **Code:** `d6ad7b9` (0.8.5) and the branch's commits; built in WSL (GCC 13, Release), tests as in
  [development.md](../../development.md).
- **Data (B2):** the reduced r226 database of the mock benchmark (`~/c4test/db/protal_db`, 761 species, built at
  `5369c99`); reads from 40 of its species' own marker genes (≥80 genes each), mutated at 0, 3, 5 and 10%
  substitutions, 500 fragments per species (20,000 reads or pairs a sample), Illumina-like substitution errors at 0.1%.
- **Machine:** WSL on 4 cores (`taskset -c 0-3`), shared.

| Script | What it does |
|---|---|
| [b2_reads.py](b2_reads.py) | reads of known divergence from a database's own marker genes |
| [b2_count.py](b2_count.py) | a sample's fragments aligned, and those on their own species, from its SAM |
| [b2_gate_drop.sh](b2_gate_drop.sh) | every divergence and read setup aligned ungated, at C 4, and at C 4 with 6 least picks |

## B2: what the gates drop, and the least picks

The gate looks a read up by every C-th of its cores (the stride pick, the default) and gives it up without a hit
among them. A read's cores are its syncmers: about 12 for a 75 bp read, 26 for a 150 bp one, so at C 4 a 75 bp read
has about 3 chances and a 150 bp one 7. The audit estimated the loss binomially; measured (`b2_gate_drop.sh`, the
fragments aligned to anything, as a share of those aligned ungated):

| Reads | 0% | 3% from the reference | 5% | 10% |
|---|---|---|---|---|
| single-end 75 bp, C 4 | 0.04% | 5.7% | 12.1% | 24.7% |
| … with 6 least picks | 0 | 0.3% | 1.4% | 6.0% |
| single-end 100 bp, C 4 | 0 | 2.9% | 9.0% | 24.6% |
| … with 6 least picks | 0 | 0.4% | 2.0% | 8.7% |
| single-end 150 bp, C 4 | 0 | 0.5% | 3.4% | 17.8% |
| … with 6 least picks | 0 | 0.2% | 2.0% | 13.7% |
| paired-end 2x75 bp, C 4 | 0 | 0.8% | 4.5% | 18.7% |
| … with 6 least picks | 0 | 0 | 0.2% | 3.5% |
| paired-end 2x150 bp, C 4 | 0 | 0 | 0.08% | 8.1% |
| … with 6 least picks | 0 | 0 | 0.04% | 5.6% |

(Of 20,000 reads or pairs; ungated, 32-78% of the 10% reads aligned at all, at the identity floor of 0.9, which is
why the audit's binomial estimate, 35-59% for the short reads at 10%, overstated the gate's share there.) The
fragments on their own species move with them (e.g. single-end 75 bp at 5%: 15,403 ungated, 13,536 at C 4, 15,188 with
the least picks). Without a mutation the gate loses nothing: the test reads come from the database's own genes.

**The fix:** `--sparse_first_min_picks N` samples a read (each mate, a long read's chunk) that would have fewer than N
picks at C at C halved until it has them (`stage1_probe::StrideFor`; with N 6 at C 4: a read of 11-20 cores at C 2,
of fewer than 11 at C 1). A read is still kept on one hit, so the rule changes only which reads are given up. It is a
setting of the database, as the gates are: `--build` records it (6 unless given), runs follow it, and a database built
before records none and keeps every read at its stride (0), as its models were trained. Tests:
`Stage1Probe.ShortReadsGetADenserStride` (every core count to 120, every C), the settings' round trip, e2e (the gates'
lines, a build's `build_settings.tsv`).

**The cost:** the gate's lookups (`sparse lookups` in each sample's log, unmutated reads):

| Reads | C 4 | with 6 least picks |
|---|---|---|
| single-end 75 bp | 2.45 a read | 7.55 |
| single-end 100 bp | 3.60 | 7.35 |
| single-end 150 bp | 5.91 | 7.17 |
| paired-end 2x75 bp | 4.90 a pair | 15.07 |
| paired-end 2x150 bp | 11.83 | 14.32 |

A read kept is seeded in full either way, so the cost falls on the reads given up (host DNA, organisms the database
lacks), where the gate saves its time: for those, a 75 bp read's gate costs three times as much, a 150 bp one's a
fifth more (a 150 bp read with fewer than 21 cores gets a denser stride too). Not measured on a host-rich sample; the
next build's `--test-c-gates` table gives the alignment time per gate with the least picks on.

## L2-L13

Each was checked again against 0.8.5's code (`d6ad7b9`); line numbers are of that commit.

**L2: clipped bases at a gene end.** Still true. A short read running past its gene's end is aligned with the gene
boundary fixed and its overhang free (`AlignmentStrategy.h:74-83, 664-665`); WFA2's free bases become soft clips
(`AlignmentUtils.h:727-728`), which the bitscore ranking every candidate and giving the MAPQ (`Score(2,3,1,2)`)
scores 0. A congener whose gene covers the read wins whenever its extra mismatches are fewer than 2/5 of the clipped
bases (30 aligned + 120 clipped scores 60 against a covering congener's 285 with 3 mismatches). Now
`--gene_end_clips credit` (recorded with the database; `--build` records `credit` unless given, a database built before
stays `off`): each candidate's bases clipped past its gene's
ends, beyond those of the read's least clipped candidate, score the identity floor's rate, every candidate of a taxon
as that taxon's least clipped one. On the reduced r226 database it removes 6% of the fragments on a wrong species and
gains a few right ones; crediting each candidate by its own clips had split a read's own candidates at gene junctions
(−2,645 right fragments). Details, data and the measurement: [L2-gene-end-clips.md](L2-gene-end-clips.md).

**L3: adaptive candidates.** Partly true. The anchors are tried in order of their exact-match length
(`ChainAnchorFinder.h:1291-1294`), the only score before alignment; the taxid order of the seeds only breaks ties.
"Genus" is the taxonomy's genus node, so a GTDB split genus (`Escherichia_A`) is another genus and never tried
(`Align.cpp:245-251`). Of the 75,195 own fragments of r226 v15's missed taxa that landed elsewhere, 98.7% landed in
the same genus ([error-read signatures](../2026-10-07-error-read-signatures/README.md)): split genera can explain at most
~1.3%, and merging genera by name would be arbitrary (the suffixed genera are separate clades). Whether the cap of 7
or their order loses own reads was never measured, so the log now counts, per sample, the reads that had more eligible
congeners than `--adaptive_candidates` tried (no output changes). If that count is large, the next step is to order
only those reads' eligible anchors (by exact-match length, then unique seeds) before the cap.

**L4: exact ties.** True in substance. Pairs of equal bitscore are sorted with `std::sort` (`Classify.h:507-511`),
which libstdc++ does by insertion, stably, for 16 or fewer (the seeds' taxid-first order: the lower internal taxid
first) and by introsort beyond (crowded reads): deterministic for a read, by taxid or not. Their MAPQ is 0, so no
taxon's `fragments` counts them, but the first one's record counts in `fragments_all`, `low_mapq_share`,
`congener_fit_share` and the denominators of the failed-candidate rates. Training and runs order alike; documented at
the sort. A stable sort with a taxid tie-break would move only the MAPQ-0 records of crowded reads: left for a
retraining.

**L5: a long read's `fragments`.** True. Each gene segment of a long read starts a group of its own
(`Profiler.h:5264-5276`, at every supplementary record) and gets its own read id, so the taxon's `fragments` count a
long read once for each segment with a kept record on it; `fragments_all` counts it once (by its link). features.md
said "a long read once"; it now says so, and that `--fragment_count soft` therefore adds mostly to taxa few of whose
reads keep a segment (absent ones: at r226 v25 the pb/ont rows of present taxa moved 3-4%, of absent ones 28-30%).
Training and runs count alike: no F1 effect.

**L6: `unique_cores.bin` that does not fit.** True. The fit test (`UniqueCores.h:267`) is the most exact the format
allows (the row's bytes hold its cores and up to three of padding), but a gene that does not fit, or one beyond the
table, fell back to the anchor's counts silently (`Mismatched()` was never read), so one sample could mix the two
schemes. Since the build writes the table the run reads, a default run never meets this; a copied table of another
database would. Now a run checks the table's rows per taxon against the database's genes before aligning
(`unique_cores::RowsMismatch`) and does not use one that differs (a warning; the whole run counts from the anchors,
as the SAM settings line then says), and the genes whose cores do not fit are counted over all of a sample's threads
(a count the counter's copies share) and warned about per sample. Tests: `UniqueCores.TheTableFitsTheDatabasesGenes...`,
the shared count.

**L7: single-value cores of old indexes.** True. A key of one value took its entry's unique flag whatever the index
said (`KmerLookup.h:390, 409`), though only an index built since `78f7a92` (2026-09-26) checked such values against
the full reference (feature bit 3, `Seedmap::ChecksSingleEntryUniques`). Such a database's models were trained when
single values never counted as unique (they do since `11cdf6b`, 2026-10-06). Now they count only where the index
checked them, and the load warns once for an index that did not. No change for databases built since 2026-09-26. Test:
`PackedIndex.SingleValuesCountAsUniqueOnlyWhereTheBuildCheckedThem` (the same index under an older header).

**L8: reduced databases' uniqueness.** True. With `--build_gene_subset` (and the converter's `--genes` folders) the
uniqueness pass leaves the other genes' copies out (`Build.h:2413, 3243`), so an exact k-mer that a kept gene shares
with another gene family of another taxon stays unique, and a read of that left-out gene can hit it. The comments
said such a read could not occur; they now give the real reason (the full reference's cost, and that such sharing is
rare: orthologs keep their gene id and are compared). Computing it against every gene would undo the 2026-10-04 speed
fix (the pass ~10x longer for 12 of 120 genes). Documented in databases.md; full databases are not affected.

**L9: the AVX-512 `--flex_indels` kernel.** True: its unit test (`FlexScan.Avx512ScoresAndTiesAsTheScalarScan`) skips
without AVX-512, and this machine has none; the non-indel kernels matched byte for byte on a Zen 4 node. Documented
(installation.md, running.md): use `PROTAL_SIMD=avx2` with `--flex_indels` on AVX-512 CPUs until a run there checks
it (`protal_tests --gtest_filter='FlexScan.*'` on such a node, or under Intel SDE). Default runs never use this kernel.

**L10: N in a read.** True. The k-mer scan codes anything but A, C, G, T as 0, A on one strand and T on the other
(`KmerIterator.h:215-226`, the AVX2 path too), while the index holds no k-mer over an N (`Build.h:60-82`): a read's
k-mers over an N were looked up as if it were A or T, which favours A/T references in seeding (the alignment counts
the N as a mismatch everywhere). Now a read's k-mers over an ambiguous base are not looked up, at its five k-mer sites
(single-end, both mates, a long read's chunks in the gate and in seeding); a read of A, C, G, T only is scanned once
and left as it was. About 1.6% of simulated 150 bp reads carry an N. Test:
`Syncmers.ReadKmersOverAnAmbiguousBaseAreLeftOut` (at every SIMD level).

**L11: mate files in another order.** True: with equal counts, mates named apart only warned
(`Align.cpp:602-606`), and every pair became two unrelated reads, which the mates' vote and the pair features then
judged. Now a sample whose pairs are more than half named apart fails (no SAM, exit 1 for the run), a few keep the
warning, and a mate named without its number matches one with it (`r` and `r/2`; swapped numbers still differ).
Tests: SeqReaderPE names, e2e `rotated` (the second file one read ahead: fails), `names` (one pair in ten renamed:
warned).

**L12: `-march=native`.** True, with such flags: WFA2 was compiled with `CMAKE_C_FLAGS` and without protal's baseline,
so `-march=native` turned on its AVX2 extend kernels (`wavefront_extend_kernels.c:106, 130, 160`), which choose
otherwise among ends-free alignments of the same score ([wfa-avx2](../2026-10-02-wfa-avx2/README.md)). WFA2 now takes
`isa_baseline` (`-march=x86-64` after the flags); its cloned loops still take x86-64-v3 at run time. Default builds
are unchanged. Checked: configured with `-march=native` in both flag sets, the kernels' compile line ends
`-march=native -march=x86-64 -mtune=generic`, and `__AVX2__` is not defined there.

**L13: `--depth_identity_margin`'s help.** True: it said the reads below the margin still count "for detection", but
the margin sets two of the model's features (`depth`, `low_identity_share`). The help and running.md now say the
margin changes the calls, which is why a run takes the database's.

## What this does not settle

- **Nothing here ran at GTDB scale.** B2's least picks change what the next build trains on for short reads (its
  `--test-c-gates` table will show their cost and gain); L2's credit is the default of new builds, so the next build trains with it;
  its effect at GTDB scale is unmeasured.
- **L10 acts on databases built before too:** reads with an N lose the k-mers over it, where those databases' models
  were trained on them looked up as A or T. Small (1-2% of simulated reads carry an N, and only the k-mers over it
  change), but not measured.
- **L3:** the new count says how often the cap is reached; whether reordering the capped reads' congeners pays is open.
- **L9** needs a run on an AVX-512 node.
- **L4's** tie-break and L5's soft counting for long reads are left as they are: changing either changes features and
  needs a retraining for little expected gain.
