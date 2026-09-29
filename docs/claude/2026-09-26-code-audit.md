# Protal code audit

2026-09-26 · Falk Hildebrand with Claude Code

> Exported on 2026-09-29 from the Claude Docs page "Protal code audit", where the audit was written
> during the work on branch `audit-fixes`. Commit hashes are on that branch. The last row (step H,
> `d35f158`) is the state when the page was last edited; later commits added the single-file
> database (`66016e9`) and the GTDB build-and-train workflow (`e917af3`, see
> [building-a-database.md](../building-a-database.md)).

Protal's read-mapping core is well designed and fast, but eight verified bugs in the SAM round-trip and the SNP/MSA path silently change profiles and strain output. Build, packaging and tests lag behind the code. Scope: branch `alpha` at commit `014f4a9`. A clean clone builds in about a minute with no errors; it was compiled with `-Wall -Wextra` under gcc 13 and clang 18, and the CLI was exercised on edge cases. There was no end-to-end profiling run: that needs a GTDB index, and the simulator needs `art_illumina`.

## What works well

The core design choices are sound and worth keeping; the problems below are mostly local and fixable.

- **Index design** (`src/Hash/Seedmap.h`): a direct-addressed table over the 15-mer core, so hash collisions are impossible. Control blocks of 8 keys cost about 2 bytes per key, and the file header carries a magic number and a version-range check.
- **Canonical k-mers**: orientation comes from an odd-length middle core, so there are no palindrome ties. Syncmers are computed on that core, so sampling is the same on both strands.
- **Hot path**: components are templates with no virtual dispatch. Each thread has its own aligner, lookup and buffers. R1/R2 are read in lockstep batches, and a read-count mismatch is detected. SAM output is buffered per thread and flushed in whole lines.
- **CLI validation**: every input is checked before the index loads, all missing files are listed at once, and bad input exits non-zero (verified: codes 2, 30, 31). Help text is detailed, and `--map_help` has a worked example.
- **Strain output**: MSA row and column order is deterministic, the IUPAC table is correct, each rejected site is attributed to a named filter, and the qcmsa call is properly shell-quoted.
- **Simulator**: zero compiler warnings under `-Wall -Wextra`. ART runs via `fork`/`execvp` with exit-status checks, read counts sum exactly to the target, and every run records its seed and a replayable manifest.
- **Documentation of recent work**: the strain M1–M5 commits have descriptive messages, and `scripts/qcmsa.md` is a good operator guide.

## Result-changing bugs

Eight defects alter what protal reports; each was confirmed by reading the code at the cited line, and the simulator hang was reproduced.

| # | Where | Defect | Effect | Fix step |
| --- | --- | --- | --- | --- |
| 1 | `IO/AlignmentOutputHandler.h:488` | Validity guard tests `ar1` twice instead of `ar2`; flags persist across candidates; `return` ends the loop | Pairs where only R1 aligns are never written | 2 |
| 2 | `SequenceUtils/VariantHandler.h:257`, `Profiling/Strain.h:246` | Synthetic REF allele gets forward-strand counts only; strand filter is on by default | Any site with one or more non-REF reads becomes `N` when REF is the consensus; N rate grows with depth | 1 |
| 3 | `SequenceUtils/VariantHandler.h:14` | Quality offset 36 (should be 33), stored in `uint8_t` | All base qualities understated by 3; Q0–2 wrap to 253–255, so NovaSeq's lowest bin counts as the best | 1 |
| 4 | `IO/SamHandler.h:187` | Read2 strand read from the mate-reverse bit (0x20) | Strand-bias filter tests fragment orientation, not read strand | 1 |
| 5 | `Profiling/Strain.h:726` | IUPAC candidates exclude the reference allele | REF/ALT mixtures (the usual two-strain case) never get an ambiguity code | 1 |
| 6 | `IO/AlignmentOutputHandler.h:103` | Read name always loses its last two characters | SRA-style IDs collide in blocks; uniqueness is judged by neighbouring names, so distinct reads look like multi-mappers | 2 |
| 7 | `Hash/Seedmap.h:1022` | Loop index `idx` never incremented while scoring flex keys | "Unique at distance ≥2" flag is wrong in built indexes, including the shipped DB; feeds the `ZT` tag | 5 |
| 8 | `SequenceUtils/Minimizer.h:38` | Syncmer mask stored in `uint8_t` | Syncmers compare 4-mers instead of 7-mers; build and query agree, so density changes rather than lookup breaking | 5 |

The reviewers also reported, not re-verified here: training and inference compute some random-forest features differently (`Profiler.h:1217` vs `:893`); `--no_strains` changes species calls; the ML scripts tune on the full data before the train/test split.

## Robustness and operations

Failures are hard to detect from outside: most errors print a warning and the process still exits 0.

| Area | Problem | Consequence | Fix step |
| --- | --- | --- | --- |
| Exit status | `main` always returns 0; FASTQ errors, pigz failures, missing SAMs and qcmsa failures only warn | Snakemake/Nextflow cannot see failed samples | 3 |
| Crash safety | SAMs are written straight to their final name, and existing SAMs are skipped on rerun | A truncated SAM from a crash is reused, then profiling calls `exit(120)` | 3 |
| Output dirs | In `-1/-2/-o` mode, `strains/` and `misc/` are never created, and those streams are never checked | Strain output (on by default) is silently lost | 3 |
| Launcher | `os.system(' '.join(argv))` re-splits arguments and discards the exit code; crashes if `protal_avx2` is missing | `--qcmsa_args "…"` breaks; the installed `protal` always exits 0 | 3 |
| Install | `just install` copies `scripts/protal_launcher` (does not exist) and installs `protal_baseline`, while the launcher looks for `protal_plain` | Install aborts; no fallback on CPUs without AVX2 | 3 |
| Simulator | Read-count correction loop cannot terminate when read pairs < species or < strains (`CommunityProfileDesigner.cpp:541`) | Hangs at 100% CPU with no message (reproduced) | 3 |
| Reproducibility | pigz runs without `-n`; the manifest omits fragment and ART settings | README's "byte for byte" replay holds for FASTQ content only, with default ART settings | 3 |
| Throughput | `IsAlignmentValid` wraps its whole body in a global `omp critical`, for every alignment (`SNPs/SNPUtils.h:110`) | Alignment threads serialise on one lock | 4 |
| Index loading | No checks after the header: sizes, stream state and derived fields are trusted | A truncated or corrupt index ends in out-of-bounds reads | 5 |
| Memory safety | `Variant` destructor test is inverted and uses `delete[]` on `new`; the struct is copied shallowly (`SequenceUtils/Variant.h:94`) | Every indel string leaks; fixing only the test would double-free | 1 |

## Engineering debt and repo hygiene

Nothing here changes results today, but it hides bugs like the ones above and slows every change.

- **Structure**: everything compiles as one translation unit of headers. The largest are `RunProtal.h` (1,734 lines), `Profiler.h` (2,098) and `Options.h` (1,653), and each option is defined in three places.
- **Asserts never fire**: `-DNDEBUG -O3` is hard-coded for every build type, including `CMAKE_CXX_FLAGS_DEBUG`.
- **Exits and output**: 132 hard-coded `exit(N)` calls, some inside OpenMP regions; 10 `Utils::Input()` stdin waits on production paths; 878 `std::cout` uses and no logger.
- **Dead or duplicate code**: `src/Utils/` (a drifted copy of `Utilities/`), `WFA2Wrapper.h` vs `WFA2Wrapper2.h` (identical), `AnchorFinder.h`, `SeedmapFlex.h`, `compact_map.*`, the fully commented-out `Coverage.h`, `Test.h`, `seeding_strategy_safe.txt`, `src/CMakeLists_bu.txt/`, the unused `lib/wfa2-lib-2.3.4`, and the htslib tarball.
- **Warnings**: 296 in own code under `-Wall -Wextra`, mostly `-Wreorder` and unused variables. One is `std:vector` at `RunProtal.h:925`, which compiles only as a goto label.
- **Tests**: three files, none wired into CMake, and no CI. `test_Coverage` cannot compile because `Coverage.h` is commented out, `test_CigarIterator` tests code nothing uses, and `/tests` is in `.gitignore`.
- **Committed build output**: `src/gzstream/*.o`, `libgzstream_lib.a`, `lib/gzstream/gzstream.o`. The build writes `src/protal_config.h` into the source tree, and neither it nor `build/` is ignored.
- **Packaging**: `LICENSE` is GPLv2 but `meta.yaml` says MIT with a placeholder summary. The recipe builds GitHub `alpha`, not the local checkout the README describes, and installs neither qcmsa nor the simulator.
- **Hard-coded personal paths** in `environment.sh`, the justfile and the R scripts. `strain-protal` in the justfile ends in `|| true`, so a failed run still produces a report.
- **README**: `{r bash}` code fences, GTDB r214 stated while the justfile uses an r226 DB, and no usage docs for profiling, the map file or strain mode.

## Fix plan

All six steps are done on branch `audit-fixes` (9 commits, not pushed). The Notes column records what landed and what still needs a human.

| Step | Scope | Status | Notes |
| --- | --- | --- | --- |
| 1 | SNP/MSA correctness: REF strand, quality offset, read2 strand, IUPAC with REF, `Variant` ownership | Done | Commit `2da2748`. 10 unit tests, clean under ASan/UBSan. Quality thresholds now compare true Phred scores (about 3 per base more permissive). New `refs_retained` column in `snp_stats.tsv`. |
| 2 | SAM round-trip: R1-only pairs, pair reading, read names; then revalidate the random forest | Done | Commit `fa0b5cd`. 19 unit tests, clean under ASan/UBSan. Single mates now also feed strain analysis. Still needs a human: revalidate the random-forest model, since R1-only pairs now reach the profiler. |
| 3 | Operations: exit codes, atomic SAM writes, output dirs, launcher and install, simulator hang, reproducibility | Done | Commit `ced6a2c`. Exit 1 plus a summary on any sample or output failure; SAMs written as `.partial` and then renamed; shell launcher. Checked end to end on the mini DB. Existing simulator seeds reproduce unchanged. Please confirm the licence choice in `meta.yaml` (GPL-2.0-only vs -or-later). |
| 4 | Remove the global lock in `IsAlignmentValid` | Done | Commit `fd5aa2d`. Lock-free check; out-of-bounds reads fixed. No measurable speed change on the mini DB (469 vs 453 ms, 208k pairs, 8 threads). SAM output unchanged. |
| 5 | Next index version: flex-key `idx` fix, syncmer mask, validation on load | Done | Commit `e02ed99`. Index format 2; format-1 DBs still load and give byte-identical output. Same alignment sensitivity on the mini DB (196,192 vs 196,196 of 208k pairs). Long super-uniques 2,870 → 5,059, so the random forest must be retrained before using a newly built DB. |
| 6 | Hygiene: dead code, build flags, gtest, CI, tiny fixture DB | Done | Commits `c559a84`, `443a1b3`, `740fb05`. 275 dead or generated files removed; Debug builds now enable asserts. Mini-DB generator committed. 11 end-to-end tests (Python) and a CI workflow; 30 unit tests pass under ASan/UBSan. |

## Measured on the mini DB

The fixes change results where the audit said they would, without costing sensitivity or speed. Before = `alpha` at `014f4a9`; after = `audit-fixes`. Data: 3 synthetic species, simulated reads with sequencing errors only.

| Check | Before | After |
| --- | --- | --- |
| Strain output in `-1/-2/-o` mode | none written, exit 0 | written |
| Pairs where only read1 aligns (4,905 pairs) | 0 SAM records | 4,905 |
| Variant positions written as N (*Mockella alpha*, one sample) | 1,747 of 1,747 | 521 of 1,739 |
| N in strain MSAs | 0.76% | 0.13% |
| Simulator with fewer read pairs than species | hangs | exits 1 in under 1 s |
| Truncated `index.prx` | not checked | exit 8 with a message |
| Alignment of 208k pairs, 8 threads (step 4 only) | 469 ms | 453 ms |
| Read pairs aligned, index format 1 vs 2 | 196,196 | 196,192 |
| Long super-unique k-mers in `unique_kmers.tsv` | 2,870 | 5,059 |

## Still open

- Retrain and revalidate the random forest. R1-only pairs now reach the profiler (step 2), and new indexes change the `lsu*` features (step 5). Existing DBs keep working with the shipped model.
- Confirm the licence in `conda-recipe/meta.yaml`: GPL-2.0-only or GPL-2.0-or-later.
- Re-check `--snp_min_phred_sum` and `--snp_min_mean_qual`: they now see true Phred scores, about 3 higher per base.
- CI has not run yet, because the branch is only local.

## Second round: integration findings

The components each work alone, but their hand-offs trust each other blindly. Bad or unusual inputs therefore give exit 0 with silent garbage, or crash after hours of alignment. Five reviewers tested the contracts with Release and ASan/UBSan builds, the mini DB, ART simulations, qcmsa and IQ-TREE; the rows marked below were re-checked in the code. Also landed since round one: licence GPL-2.0 (`75f8b5b`) and WFA2-lib v2.3.6 (`c9291e2`, byte-identical output).

| # | Severity | Contract | Finding | Verified |
| --- | --- | --- | --- | --- |
| 1 | High | Profiler ↔ model training | The training dump writes RAF/RA features from index 1 for every column, and `su_rate` is computed from `lu`; the dump's rate columns are also named differently from the model's. Retraining learns from the wrong quantities, and the Python trainer aborts on the missing `*_rate_ref`. | code |
| 2 | High | DB files ↔ query | A DB without `unique_kmers.tsv` silently loads without uniqueness data: every taxon is rejected and empty profiles are written with exit 0. | code |
| 3 | High | Reference ↔ build | Missing `--reference`/`--full_reference`, or a lowercase reference (encoded as A), yields an empty or unchecked 3 GB index with exit 0. | code |
| 4 | High | Index ↔ reference.fna/map | Nothing binds the two (no checksum, counts or lengths): a mismatched pair gives a wrong profile with exit 0, or out-of-bounds reads. | run |
| 5 | High | Truth file ↔ GetTruth | A truth file whose first line has no `d__`/`s__` field (including the documented integer format) reads `tokens[-1]` and crashes profiling for all samples; `--profile_truth` on the command line is ignored. | code |
| 6 | High | Map file ↔ parser | A row missing its PROFILE cell gets the next sample's profile written under its name, then exit 33 after alignment; an empty cell shifts all later columns (`LineSplitter` skips the character after a delimiter). | code |
| 7 | Medium | SAM ↔ strain MSA | Deletions get quality 0 and never pass, so the MSA shows N plus reference bases where the strain has a deletion. | code |
| 8 | Medium | Partition file ↔ IQ-TREE | Partitions are 0-based; `iqtree2 -p` fails with "Negative site ID" (qcmsa keeps the base). | run |
| 9 | Medium | MSA writer | With `--msa_min_hcov 0` (strain-test-raw), a species with no usable gene segfaults on `partitions.back()` of an empty list. | code |
| 10 | Medium | Coverage counters | 16-bit coverage wraps above 65,535 reads (depth reported as 464 for 66,000); a variant reaching 65,535 observations calls `exit(9)` silently and ends the whole run. | code |
| 11 | Medium | Profiler ↔ qcmsa | Samples absent from the MSA still count toward `--gene-min-samples`; "qcmsa wrote nothing" counts as success; the multi-allelic filter is inert below 4 samples; qcmsa's allele rule differs from the MSA's IUPAC rule. | run |
| 12 | Medium | Options ↔ random forest | `--no_strains` changes RF probabilities and can drop a taxon (allele features come from the strain containers). | run |
| 13 | Medium | Abundance model | Abundance is the median over genes with at least one read: 7.4× too high at 23 read pairs, 2.3× at 105, correct from \~1,000 pairs. | run |
| 14 | Medium | Reads with N ↔ SNP calling | N bases become X ops and 'N' alleles: 2% N gives 23,259 spurious variant positions that feed RF features and qcmsa. | run |
| 15 | Medium | DB side files | `reference.map` offsets and `internal_taxonomy.dmp` are not validated; one bad offset turns all later genes into NUL bytes; a dangling taxonomy parent waits on stdin (`Utils::Input`). The model is loaded only after alignment. | run |
| 16 | Medium | Harness and tools | Nothing produces the strain-test input layout, and its map assembly can realign everything into the repo; `protal_map_utils merge` and its validator resolve relative paths differently from protal; `--benchmark_alignment` aborts on simulator read names. | run |
| 17 | Medium | Marker completeness | Reads from a species' marker gene that is missing from its reference pile onto another species' copy: \~185 false ambiguous columns per sample in the mini DB. | run |
| 18 | Low–Medium | External SAMs | Tags are read by column position, unmapped `*` records crash, `=`/H CIGARs are handled inconsistently, and a taxon with zero valid reads divides by zero (SIGFPE). | run |

Smaller, all confirmed by runs: the last multi-mapped read of each SAM is dropped (`-m` ≥ 2 then depends on thread count); a blank line in a SAM ends reading; the launcher can pick `protal_avx2` from another install on `$PATH` (my step-3 bug); gzstream's putback `memcpy` overlaps (UB, aborts ASan); `.profile.genes.log` header is one column off; `protal_profile_utils` cannot read `.profile.gz`; a header-only SAM gets no profile and no error; ART aliases seeds ≥ 2³¹; strain-sharing specs are not exclusive.

## Contracts that hold

On well-formed input the pipeline is sound and reproducible.

- **Accuracy at normal coverage:** species detection was exact on 12 simulated samples across all distributions. Abundance was within 0.0065 of the simulator's truth, with no bias for multi-strain species.
- **Strain recovery:** identical strains are at p-distance ≤ 5e-5 and different strains at ≥ 0.0087; unpartitioned IQ-TREE recovers the clustering with 100% support.
- **Determinism:** index builds are byte-identical at 1 and 4 threads. Profiles, MSAs and statistics are identical across thread counts with `-m 1`, and between preloaded and lazily loaded genomes. AVX2 and baseline outputs are identical.
- **SAM writer invariants:** mates are adjacent and share a QNAME, orphans carry 0x8, POS is 1-based, and no alignment runs past its gene. Reverse-strand SEQ and QUAL orientation is handled correctly (3,177 of 3,177 SNP positions right).
- **File hand-offs:** protal, qcmsa, strain\_report and the justfile recipes agree on file names, sample names and meta columns. All 31 model inputs are provided under matching names. `.sam.gz` round-trips.
- **Simulator:** replays are byte-identical, and changed ART settings are warned about. Its `--protal_metafile` output runs directly in protal, with truth annotation.

## Proposed round-2 fixes

Most findings share one cause: inputs are trusted instead of checked. So the first step is failing fast, before any alignment.

| Step | Scope | Findings |
| --- | --- | --- |
| A | Fail fast on bad inputs: require `unique_kmers.tsv`, check build inputs exist and uppercase the reference, record reference size and gene count in the index header and check them at load, validate `reference.map`, taxonomy and truth files, load the model up front, and wire up `--profile_truth` | 2, 3, 4, 5, 15 |
| B | Parsing: fix `LineSplitter` empty fields, check map column lengths, fix the last multi-mapped read, SAM blank lines, unmapped/`*` records, tag names instead of positions | 6, 18, smaller |
| C | Strain output: call deletions, drop N alleles, 1-based partitions (with qcmsa), handle empty MSAs, 32-bit coverage and observation counters | 7, 8, 9, 10, 14 |
| D | Random-forest contract: one feature builder for inference and the training dump, fixed RAF/RA/`su_rate`, allele features independent of `--no_strains`; then retrain | 1, 12 |
| E | qcmsa and tools: count only MSA samples, report empty output, align allele rules; strain-test glue, `protal_map_utils` paths, launcher order, `.profile.gz`, benchmark read names | 11, 16, smaller |

Two items need a decision rather than a patch: the low-coverage abundance model (13) and marker completeness in the DB (17).

## Round-2 fixes: status

All five steps are on branch `audit-fixes`, one commit each, not pushed. 60 unit tests and 31 end-to-end tests pass; the unit tests also pass under ASan/UBSan. Each e2e test added for a fix was checked to fail on the commit before it.

| Step | Commit | What changed | Findings |
| --- | --- | --- | --- |
| A | `5c93647` | Every DB file and input is checked before alignment: `reference.map`, taxonomy, truth files, model, `unique_kmers.tsv`, `--reference` against the map. The index records a fingerprint of `reference.map` and the size of `reference.fna`. `--profile_truth` works. | 2, 3, 4, 5, 15 |
| B | `969b241` | Empty map fields are kept; every map row needs every declared column. New `SamReader`: blank lines, unmapped and `*` records, tags by name, `=` and H ops, per-sample errors. The last multi-mapped read is kept. | 6, 18, smaller |
| C | `ca2df11` | Deletions are called; N is no allele; 32-bit counters; 1-based partitions; no crash on empty MSAs. Abundance uses the blended estimator. Genes without long uniques are left out of the MSA for species with relatives in the DB. | 7-10, 13, 14, 17 |
| D | `050e6cd` | One feature builder for scoring and the training dump. Allele features are recorded with `--no_strains` too. | 1, 12 |
| E | `7418caa` | qcmsa counts only MSA samples and uses the MSA's allele rule; empty output is reported. Launcher lookup order, map-utils paths, `.profile.gz`, strain-test glue and a `strain-input` recipe, benchmark read names, ART seeds. | 11, 16, smaller |
| F | 78f7a92 | The DB build checks k-mers whose core occurs once in the index against the full reference; index feature bit 3. | found in C |

### Measured on the mini DB

- **Abundance at low coverage** (depth over the expected depth, 1.0 = right; old/new): 17 read pairs per species 5.5-7.3 / 0.76-1.41 (Poisson noise); 35 pairs 2.7-3.9 / 0.91-1.08; 83 pairs 1.5-2.0 / 0.90-1.08; from \~500 pairs both 0.97-1.12. The weight of the median ramps smoothly with the fraction of genes hit (0.80-0.95) or the median depth (0.5x-1x), so there is no step.
- **Marker artefact**: columns ambiguous in every single-strain alpha sample, 185 before, 57 after (all in gene 36, which keeps a few long uniques).
- **Partitions**: IQ-TREE `-p` failed with "Negative site ID" before, and now runs.
- **Determinism**: outputs are identical across `-m 1`/`-m 3` and `-t 1`/`-t 4`; before, `-m 3` depended on the thread count.
- **No change where none was intended**: at normal coverage, profiles and model probabilities are identical to before, and SAM records are identical as multisets.

### Found while fixing, needs a decision

- **Short unique k-mers were never checked (DB build): fixed in `78f7a92`.** A k-mer core that occurs once in the index has no flex keys, so the uniqueness pass never compared it (`KmerLookup::GetFlex` returned nothing), and every such k-mer stayed "unique". Every gene was therefore hittable, and `su_*` counted k-mers in small blocks. The pass now reads the entry's k-mer back from its gene and applies the flex-block rule. Index feature bit 3 marks fixed indices. On the mini DB, alpha and beta lose \~200 of \~16,800 short uniques and profiles stay the same; on a full reference the change will be larger. Rebuild the DB, then retrain the model on it.
- **Option 1 needed a guard.** `long_unique == 0` marks genes that relatives lack or share unchanged, but only for species with relatives in the DB. A species alone in its genus (gamma) has long uniques in only 24 of 117 genes, all its k-mers being unique at the 15-mer core, so literal option 1 would drop 93 of its genes. The filter therefore applies only when most of a species' genes have long uniques. Checked short uniques do not replace it: the four artefact genes keep \~90% of them, because at \~93% ANI few 31-mers are shared exactly.
- **Retrain the model** on a DB rebuilt with `78f7a92` or later, using dumps from `050e6cd` or later. The shipped model still runs, but it was trained on unchecked uniques.
- Not done: the strain-sharing spec is not exclusive, and MIN\_VCOV applies per species; the manifest records design counts rather than what ART wrote; TLEN is the alignment length; `protal_map_utils merge` moves SAM paths into the default `alignments/` directory; memory use with `--no_strains` now equals that of a default run.

## Round 3: the profiling stage

Four reviewers covered the profiling stage, each reviewing code and running experiments:

- read assignment
- model features
- the presence decision and output files
- robustness and performance

The experiments used a new synthetic world: 24 species in 11 genera, with congeners, archaea and singletons, and 3 genomes each. `db20` is built from 20 of them; the other 4 are simulated but not in the DB, to test false positives. In all, about 700 simulated samples. Read names give each read's true source. The features were recomputed independently and matched protal's exactly for 247 taxa, and a separate PMML evaluator reproduces the model's probabilities exactly.

**In short:** on bacteria, protal is accurate. At knob 0.5 it detected 2214 of 2217 present bacteria with 0 false positives in 494 samples, and the median abundance ratio of detected taxa is 0.96–1.01. The weak points are the model's features, which break on archaea and depend on depth and library; two SAM-writer bugs that drop reads; reads of relatives inflating abundances; and memory and time at scale.

| # | Severity | Area | Finding | Evidence |
| --- | --- | --- | --- | --- |
| 1 | High | Model | Archaea are never called confidently. Several features are absolute counts (`present_genes`, `hittable`, `expected_gene_presence`, `total_genome`, `su/lu/lsu_genome`): archaea have 52 marker genes and 7.7k k-mers, bacteria 119 and 24k. The model (2024) predates archaea support (2026-03). | Archaea detected in 161 of 374 samples (43%); p is pinned at 0.46–0.53 whatever the depth (Calidella tepida at 93x: p 0.484). Scaling the counts to 118 genes gives p 0.969. |
| 2 | High | SAM writer | Pairs whose best alignment is mate 2 alone are never written: `AlignmentOutputHandler.h:407` tests mate 1's CIGAR, which is empty, so its ANI is 0 < 0.8. Present since 2022. | 0 of 7310 such pairs reached the SAM. About 9% of hits are lost, which is why protal's depth runs at 0.88x the truth. A one-line fix raises sensitivity from 0.91 to 0.997 and depth matches the truth. |
| 3 | High on real genomes | SAM writer | A fragment whose mates land on two different genes gets MAPQ \~0 and is dropped: mates are paired only on the same gene (`Classify.h:209`), and MAPQ then compares the two single mates. | 744 of 744 constructed pairs dropped; same-gene controls kept. db20 cannot show it (all gene gaps are 1239 bp), but real bac120/ar53 markers sit in operons. |
| 4 | High | Model | About 30% of the forest's splits are on per-species database constants (`su/lu/lsu/total_genome`, `*_rate_ref`, `hittable`), many with thresholds outside db20's range. In training, `su_rate_ref` meant lu/all: the 2024 dump wrote `lu_ref/all` into that column. | `su_genome`: every split threshold is below db20's minimum. |
| 5 | High | Performance | The profiling stage does not scale: SAM loading is serialised by `critical(load_sam)`, the whole SAM is kept as string copies (\~1.1 KB per pair; `-m 5` keeps every alternative alignment although none is read), and every finished sample's profile, including per-taxon genome copies, stays in memory until the end. | `-t 4` is no faster than `-t 2`. Removing the lock: 23 s → 14 s with identical outputs. RSS: 9.0 GB with `-m 5` vs 3.7 GB; about 170 MB kept per 0.9M-pair sample. |
| 6 | Medium | Abundance | Reads of absent or dominant congeners inflate present species and can mask rare ones. Per-read identity separates them for depth (an identity cut at 0.96 removes 90% of foreign and 0.5% of own reads), but dropping those reads for detection creates false positives. | When over 50% of a taxon's bases are foreign, abundance is 1.87x (median). Fakibacter zeta at 2.6x; Mockella alpha's depth goes from 0.52 to 6.0 next to M. novus. |
| 7 | Medium | Strains | MSAs include samples in which the model rejected the species: `GetMSAForTaxon` gets no filter, and its parameter type cannot take the forest. | 84 of 175 MSA rows in one setting came from rejected, absent samples; some had \~90k called bases. |
| 8 | Medium | Model | Features depend on depth and library, not only on presence. `stddev`/`variance` scale with depth; mates count as independent draws in `ExpectedGenes`; `lu_per_read` scales with read length; `mean_mapq` with log bitscore; `mean_ani` counts indels as matches. | The same 30 fragments score 0.609 paired vs 0.482 single-end. Novel-congener false positives reach 0.461 while archaeal true positives score 0.46–0.53, so no knob separates them. |
| 9 | Medium | CLI | The `--knob` help claim ("0.4–0.6 does not change F1 much") does not hold, and the knob is not validated. | F1 is 0.995 at 0.45 (11 false positives), 0.957 at 0.5, and 0.921 at 0.55 (no archaea). `--knob 1.5` gives empty profiles and exit 0. |
| 10 | Medium | Model | The up-front model check only loads it. | Labels other than TRUE give all-zero scores and exit 0. A missing feature throws inside OpenMP and aborts (exit 134) after alignment. |
| 11 | Medium | MAPQ | MAPQ depends on `-c` (`--align_top`), and the profile's cutoffs (4, and 20 for "unique") assume the default. | `-c 1`: "unique" share goes from 74% to 99.8%, and false hits from held-out species rise 20–45%. Calibration: at MAPQ 0–4, 22% of reads are assigned correctly; at ≥10, ≥98%. |
| 12 | Medium | Performance | cPMML opens a 4-thread OpenMP team on every `Score()`, whatever `-t` is. `Score` runs 6–9 times per taxon. The model is parsed 3 times per run, with \~229k exceptions each. | 7.6 ms per score call at `-t 1`. Model loading is 54% of all instructions on a small sample; a 1-pair run takes 1.3 s. |
| 13 | Medium | Robustness | A truncated `.sam.gz` is profiled silently: gzstream reports it as a normal end of file. | 5 of 7 truncated copies gave exit 0 and a partial profile. |
| 14 | Medium | Robustness | A SAM aligned against another DB is not detected. | Records go silently to `.err` or are compared against the wrong gene (436 "FAULTY" dumps); exit 0 with an empty profile. |
| 15 | Medium | Performance | `MSA()` builds per-position debug strings that nothing reads. | Removing them cuts the strain stage from 6.5 s to 1.4 s, with identical output. |
| 16 | Low | Robustness | An empty hittable set means "every gene is hittable". Coverage ranges can overlap, which shifts coverage positions, and the overlap check is dead (its condition repeats the `if`). | Both shown by construction; neither occurred in the simulations. |
| 17 | Low | Outputs | Assorted output defects; see the list below. | Each confirmed in runs. |
| 18 | Low | Code | Assorted code hazards; see the list below. | From code review, mostly confirmed. |

The output defects (17):

- `.genes.log`: abundance is non-zero for rejected taxa and `inf` when no taxon passes; `UniqueTwoMersReads` is always 0; the MAPQ and ANI columns are sums.
- `.profile.log`: no header; its width varies between files; `VCOV -1` for rejected taxa; `mean_gene_cov_ratios` uses an int accumulator.
- `.gene.log`: its coverage sum is always 0.
- `statistics.tsv` is written only for taxa with reads in more than one sample.
- The printed FN count includes species that are not in the DB, and the lines carry no sample name.
- `#SAMPLEID` is ignored in favour of PREFIX.
- `--profile_only` writes next to the SAM and uses the full path as the sample name.

The code hazards (18):

- `profiler::Gene` has a copy assignment that suppresses its move constructor, so the gene map deep-copies genes, about 44 times per gene.
- `GetGeneList` returns by value.
- CIGARs are parsed about 6 times per record.
- `ProcessMAPQ` makes three score computations whose results are never used.
- `m_sams` keeps dangling pointers.
- Blocking `Utils::Input()` calls remain in `MSA` and `Gene::VerticalCoverage`.
- Several accumulators are `int`.
- The trainers compute sensitivity and precision with swapped labels.

### Verified to hold

- **Features and scores:** every feature matches the independent recomputation, and the probabilities match the independent PMML evaluator. No NaN or inf.
- **Bacteria:**
  - Detected almost perfectly, with 0 false positives at knob 0.5 in 494 samples, including novel species alone at up to 400k pairs.
  - Detection limit: about 20 read pairs (88% detected).
  - Median abundance ratio of detected taxa: 0.96–1.01; two-strain samples: 0.996.
- **Congeners:** at ratios from 1:1 to 100:1, read assignment is unaffected (35750 of 35751 correct).
- **Invariance:** `-m 5` and `--no_strains` give identical profiles.
- **Determinism:** byte-identical outputs across `-t 1/2/4` and sample order.
- **Sanitizers:** TSan clean in `ProfileWrapper`; ASan/UBSan clean on all edge cases.
- **Scaling:** time and memory are linear in the number of reads.

### Proposed round-3 fixes

| Step | Scope | Findings |
| --- | --- | --- |
| G | SAM writer: write lone-mate-2 pairs (test the mate that is set); give mates on two genes their own MAPQ and write both. Changes the SAMs, so it comes before retraining. | 2, 3 |
| H | Model contract and retraining: normalised features (fractions of hittable genes, per-read and per-kb rates, coefficient of variation instead of SD); drop per-database constants; count fragments, not mates; ANI with indels as differences; train over depths, read lengths, archaea and bacteria, with novel-species negatives. Also: validate the model up front (TRUE label, every input provided, a test score), check the knob range, fix the help text and the trainers' metrics. | 1, 4, 8, 9, 10, 11 |
| I | Abundance from a taxon's own reads (identity-mode depth or EM between congeners), reporting the low-identity share. MSAs only from samples where the species passes. | 6, 7 |
| J | Performance and memory: stream SAM read groups; drop unused alternatives; parallel loading bounded by a semaphore; free per-sample data the MSA does not need; reference the shared Genome; score each taxon once; cPMML with `THREADS=1`; load the model once; remove the MSA debug strings; give `Gene` a move constructor. | 5, 12, 15, 18 |
| K | Robustness and outputs: detect truncated gzip; check `@SQ` against the DB; fix the hittable and range-overlap edge cases; headers and correct columns in the logs; statistics for single samples; attributable TP/FP/FN; `#SAMPLEID`; `--profile_only` output directory. | 13, 14, 16, 17 |

Needs a decision: H needs training data generated from simulations and a new model, and it should follow G, whose changes to the SAMs shift the features.

### Round-3 fixes: status

Order agreed: G, I, J, K, then H. Each step is one commit on `audit-fixes` (not pushed); unit, end-to-end and ASan tests pass after each.

| Step | Commit | Result |
| --- | --- | --- |
| G | `c3c89ae` | Pairs where only mate 2 aligns, and fragments across two genes, are now written (0% before, over 90% now). Each mate of a split fragment gets its own MAPQ; the profiler applies the MAPQ threshold per mate. On db20, depth/truth went from 0.878 to 0.962; no call changed. |
| I | `f6b7ba6` | Depth counts only reads within 0.04 identity (`--depth_identity_margin`) of the taxon's best reads; `.profile.log` reports the low-identity share. The depth blend no longer trusts a high median depth from a few reads on a short gene. With reads from relatives, the p90 of depth/truth fell from up to 1.57 to at most 0.93, and species above 1.2x depth from 21% to at most 1%; detection is unchanged. MSAs use only samples where the species passes. |
| J | `1825800` | The SAM is streamed read by read; samples are profiled in parallel; each sample's read data is freed once its outputs are written, except the variants and read ranges of taxa that pass the model (none with `--no_strains`). Each taxon is scored once. The model is loaded once, and cPMML is single-threaded. `Gene` now moves without throwing, so `sparse_map` stops copying genes. Four 1M-pair samples, 4 threads: peak memory 3.6 GB to 0.74 GB, profiling 20 s to 8 s; with up to 5 alignments per read, 9.1 GB to 0.74 GB and 32 s to 12 s. Outputs are byte-identical, except that with multi-mapped reads the gene lines in `.gene.log` and `.genes.log` change order. The semaphore planned for loading turned out not to be needed, since no SAM is held in memory. |
| K | `0d04d2c` | A truncated `.sam.gz` or FASTQ `.gz` is now an error: all 7 cut copies fail, where before 2 gave exit 0. So is a SAM whose `@SQ` genes are missing from the database or have another length. Rejected reads are counted; the mismatch dumps are capped at 3. Other fixes: a genome listed with no unique k-mers has no hittable gene; a read bridging several ranges merges them all, with coverage positions correct; the MSA no longer waits for keyboard input; samples may not share output files. The logs get headers and fixed columns, with no -1 depths; abundance is 0 for rejected taxa, never inf; MAPQ and ANI are means, not sums. The dead counter and the zero coverage sum are fixed, and taxa and genes are listed in id order. `statistics.tsv` is written for single samples too. The truth counts name the sample and count species the database lacks separately. `#SAMPLEID` names the sample, and `--profile_only` honours `-o`. Profiles and MSAs are unchanged except for row order. |
| H | `d35f158` | Before aligning, protal now checks the model: its inputs must all be features protal computes, it must predict TRUE, and it must score a taxon. A failure stops the run with exit 2 and names the problem; before, the run aborted after the alignment or scored every taxon 0. `--knob` must be between 0 and 1. The dump gains normalised features next to the unchanged legacy ones: fragments, fraction of genes hit, expected-gene ratio, identity with indels, top identity, low-identity share and per-kb rates. On 36 simulated samples, `hit_gene_fraction` has AUC 0.96 (archaea 0.97, against 0.78 for `present_genes`), and its archaea/bacteria ratio is 1.0 where the counts give 0.44. The shipped model finds 38 of 108 archaea there. New: `collect_training_data.py`, trainers with `--features normalized` and fixed sensitivity/precision, and `scripts/model_training.md`. Profiles are unchanged. **Open:** no model has been trained yet; that needs scikit-learn installed. |
