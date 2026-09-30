# Audit round 5: the parts not audited before

Date: 2026-09-30. Branch `strain-fixes`, commit `39a8585`: `audit-fixes` with the strain fixes and the merged `performance`
branch, which `audit-fixes` fast-forwards to. Rounds 1–4 ([code audit](../2026-09-26-code-audit.md),
[strain MSAs](../2026-09-29-strain-audit/README.md)) covered read mapping and SAM records, the contracts between files,
profiling, and the strain MSAs. This round covers the rest, except the short-read alignment core, which the `performance`
work has just rewritten.

## Setup

- **Reviewers.** Six, one per area, each reading the code and testing it; their full tables are in `areas/`:

  | Area | Table |
  |---|---|
  | the single-file database and its compression | [database](areas/database.md) |
  | the I/O layer rewritten by `performance`: gzip/BGZF input, SAM output as sam, sam.gz and sam.zst, the FASTQ reader, loading | [io](areas/io.md) |
  | long reads (PacBio HiFi, ONT) | [longreads](areas/longreads.md) |
  | the GTDB download, build and training scripts | [gtdb](areas/gtdb.md) |
  | the simulator, taxonomy, profile output, `protal_map_utils` and `protal_profile_utils` | [simulator](areas/simulator.md) |
  | build, packaging, install, CI, install and testing docs | [build](areas/build.md) |

- **Machine.** WSL Ubuntu 24.04, 8 cores, GCC 13.3, CMake 3.28, shared with other sessions; every reviewer used at most 2
  threads. Binaries: one Release build of `39a8585`; ASan/UBSan and Debug builds of copies where needed.
- **Data.** The e2e mini DB, the audit and tuning worlds, reads simulated with ART or by the reviewers' own generators
  (with a truth map per base for long reads), a synthetic GTDB release served from a local mirror, crafted or corrupted
  files. Nothing was downloaded from the internet.
- **Scripts** are in `scripts/<area>/`. Data stays in WSL under `~/audit6/<area>/` and is not in git; the scripts
  regenerate it.
- Findings were confirmed by experiment unless marked "code". I re-read the code for the High findings and for F1, F2,
  LR1 and T1.

WSL crashed several times during the round (`Wsl/Service/E_UNEXPECTED`), with load around 13 from parallel `-j 8`
builds of other sessions. Once the kernel's OOM killer ended a protal process while about 23 protal processes ran. Each
reviewer re-ran what a crash interrupted.

## Result

The core formats are sound: the single-file database round-trips byte for byte and rejects forged or corrupt files, SAM
output is identical across thread counts and its three formats, crash safety holds, long reads are aligned with correct
coordinates and profile within about 0.02 of the truth, and the exported model scores exactly as scikit-learn.

The weak points are elsewhere:

- **Silent exit 0.** Several bad inputs give an empty or shortened result with exit 0 instead of an error: paired reads in
  a format protal cannot read, reads through a pipe (a regression from `performance`), unreadable files, a damaged later
  member of a multi-member gzip. A taxonomy cycle makes protal hang.
- **The GTDB build script** is fragile for a run of days: after any failure or kill its reruns fail, a failed background
  build is reported only at the end, a stop leaves child builds running, and nothing resumes.
- **Long reads** lose the right species' alignment on some ONT reads because of the new `--x_drop` default (a regression
  from `performance`).
- **The simulator's truth** records the read pairs it asked ART for, not those written, which is off for fragmented
  genomes.
- **Installation**: the documented `just install prefix=...` installs into a folder named `prefix=...`, and the
  installed binaries search build-tree folders for their libraries.

## Findings

2 High, 23 Medium, 60 Low and 1 informational (the UBSan report found by two reviewers counted once). The table lists
High and Medium; the Low ones are in `areas/`.

| # | Sev | Area | Finding | Evidence |
|---|---|---|---|---|
| H1 | High | GTDB build | After a failure or kill, every rerun into the same folder fails: `exclude_from_db` (`gtdb_to_protal_db.py:292-294`) copies every file of `protal_db` into `training_db`, including the partial `unique_kmers.tsv` of the background build, and the `training_db` build then exits 8 each time. | Killed run, then rerun: "Invalid unique k-mer file .../training_db/unique_kmers.tsv, line 1: gene 30_46 is not in reference.map", on every rerun until the folders are cleaned by hand. |
| B1 | High | install | `docs/installation.md:76` documents `just install prefix="$HOME/.local"`; `prefix` is a positional recipe parameter, so everything lands in `<checkout>/prefix=/.../bin`, and the final check (`$(prefix=/... --version)`) cannot fail: exit 0. | `just install prefix=$P` rc 0, files under `src/prefix=/home/.../bin`, log "--version: command not found". `just install "$P"` works. |
| F1 | Medium | I/O | Paired-end input neither FASTQ nor FASTA (`.fq.zst`, `.bz2`, `.xz`, BAM, text) prints "unrecognized file format" but exits 0 with a header-only SAM and an empty profile: `Next()` returns when both blocks are invalid without checking the error (`SeqReader.h:205-208`). Single-end fails correctly. | `-1 zst_R1.fq.zst -2 zst_R2.fq.zst`: rc 0, 0 records. |
| F2 | Medium | I/O | Regression (401694b): reads from a pipe or FIFO (`-1 <(zcat a.gz) -2 <(zcat b.gz)`) give an empty result: `open()` sniffs for BGZF with a separate stream, which eats the start of the pipe (`ThreadedGzStream.h:51`, `Bgzf.h:139-144`). Paired then hits F1 (exit 0); single-end exits 1. | process substitution and mkfifo, paired: rc 0, 0 records. |
| F3 | Medium | I/O | An unreadable input file (permission denied) gives exit 0, an empty SAM and profile, and no message: the stream only sets badbit, and nothing checks `is_open()`. | `chmod 000` on plain, gzip and BGZF inputs, paired and single: rc 0, 0 records, empty log. |
| F4 | Medium | I/O | A damaged header of a later gzip member is taken for trailing garbage: all reads after it are dropped without an error (zlib-ng `gz_look`, `ThreadedGzStream.h:180-192`). | Single-end: 1962 of 3924 reads, rc 0; both mates damaged alike: 3924 of 7848, rc 0. `gzip -t` exits 2. |
| DB1 | Medium | database | A frame without a stored content size (legal zstd) skips the size check, and the seek table's claimed size (up to ~4 GB per frame and per thread) is allocated before validation (`Zstd.h` ~730-735, ~829). protal's own frames store the size, so honest corruption is caught. | Crafted frame: `std::bad_alloc` under `ulimit -v 1200000`; with the size stored: rejected before allocating. |
| DB2 | Medium | database | Stray separate files next to `database.protal` (e.g. an interrupted `--unpack_db`, which writes the index first) make protal ignore the complete bundle, and the "does not exist" errors never mention it (`Database.h` ~291-294, `Options.h` ~1604-1628). | {`index.prx.zst`, `database.protal`}: exit 30, missing reference.fna, map, taxonomy, model, unique_kmers; no word of the bundle. |
| LR1 | Medium | long reads | Regression (ff71972, tested on 150 bp pairs only): with the default `--x_drop 1000`, WFA2 fails or returns an invalid CIGAR for the own species' alignment of a gene an ONT read ends in; the relative's weaker hit is then the only one and gets MAPQ 139-149, counted as unique. `--x_drop 0` fixes it; the help says the default "prunes nothing in practice". | 300 ONT reads: 3 wrong-species records and 2 invalid alignments with 1000, none with 0. HiFi and 5%-error ONT unaffected; no speed difference. |
| LR2 | Medium | long reads | The long-read check reads only the first 100 reads: long reads without `--read_type`, behind 100 short ones, are aligned as short reads (one gene per read, the rest soft-clipped; a 140 kb read gets no record), exit 0. | 100 × 150 bp then 20-140 kb reads: rc 0, 5 of 40 and 2 of 60 whole genes found. |
| G1 | Medium | GTDB build | A failure of the background `protal_db` build is noticed only in `finish()`, after collection, training and parity (`build_gtdb_database.py:254-280, 441-443, 479-480`); the likeliest one at GTDB scale is out of memory. | Background build exits 137; the script collects, trains, checks parity and reports it last. |
| G2 | Medium | GTDB build | SIGTERM or an uncaught exception leaves both child `--build`s running; a rerun then runs four, which delete each other's files. | After `kill -TERM`: both builds alive; rerun's training build: "reference.fna does not exist". |
| G3 | Medium | GTDB build | No resume: every rerun repeats the conversion and both index builds (~1.5-2 h each at GTDB scale); only the collector skips work. | Rerun of a finished run: both builds again, same md5. |
| G4 | Medium | GTDB build | The collector does not check which database and options made its dumps: after a change of `--seed`, holdouts, genome table or protal it reuses the old ones and trains on the mix; only the parity check fails, after training. | `--seed 2` rerun: "47 species left out" (was 45), 0 s collection, trained, parity fails on a foreign `@SQ` gene. |
| G5 | Medium | GTDB build | Missing protal, simulator, `art_illumina` or scikit-learn/joblib are found only when first used, e.g. joblib after the builds and the whole collection. | `/usr/bin/python3` without joblib: `ModuleNotFoundError` after the collection. |
| G6 | Medium | GTDB download | A complete `.part` gets HTTP 416, is retried 3×, and the run exits, on every rerun; the message does not say to delete it (`download_gtdb.py:150-169`). | Fake mirror: 416 × 3, exit 1. |
| G7 | Medium | GTDB download | No resume within a run: a dropped connection returns short data, the MD5 fails, the `.part` is deleted and the download restarts at byte 0; after 3 drops the run exits. | Server closing after 5000 B: 3 GETs without Range, exit 1. |
| S1 | Medium | simulator | The manifest, abundance matrix and gold-standard profiles record the requested read pairs and coverage, not what ART wrote. ART splits one coverage over the contigs in whole reads: fragmented genomes get from −100% to +80% of the requested reads, contigs shorter than a read none, while the truth still lists the species. | One genome, 2000 pairs requested: 1998 (3 contigs), 1644 (500 bp contigs), 0 (140 bp). Mini world (3 contigs): −0.07%. |
| T1 | Medium | taxonomy | `IntTaxonomy::Load` checks duplicate ids and dangling parents but not cycles or a root with a parent; `LineageStr` then loops forever (`Taxonomy.cpp:717`, stops only at a self-parented node) while the profile is written. | Family's parent set to its genus: killed by `timeout` after 120 s, 0-byte profile, no message. |
| T2 | Medium | taxonomy | Duplicate taxon names are accepted and the name lookup keeps the last: truth, `--msa_species` and statistics point at the wrong taxon. | "TP 2, FP 1", "Species unknown to taxonomy: s__Mockella beta"; truth labels of two taxa swapped. |
| U1 | Medium | map utilities | `protal_map_utils validate`, `flatten` and `merge` treat the single-end marker `-` in SECOND as a file: `flatten` writes `<dir>/-`, which protal rejects; `merge` refuses two single-end rows; `validate` rejects maps without SECOND, which protal accepts. | validate: "Missing read files: .../runA/-"; protal on the flattened map: exit 30. |
| M1 | Medium | simulator | Still open from round 2: the strain-sharing spec is not exclusive (the species is also drawn at random elsewhere, with any strain) and MIN_VCOV floors the species total, not each strain; `docs/simulation.md` says "per strain". | Spec for 3 of 6 samples and 3 strains: species in 6 of 6 samples with 5 strains; strains down to 0.34x under a floor of 5x. |
| B2 | Medium | build | The dynamic binaries, also once installed, carry a RUNPATH of 13 absolute build-tree folders (11 nonexistent, from `../IO`, `../Hash`, … relative to `src/`); the loader searches them first for libzstd, libdeflate, libgomp, libstdc++, and there is no `install()` to strip it (`src/CMakeLists.txt:141-172`). | `readelf -d prefix/bin/protal_baseline`: `RUNPATH [.../src/src/../IO:...]`. |
| B3 | Medium | build | `CMAKE_AR`/`CMAKE_RANLIB` are forced to `$GCC-ar`/`-ranlib` (`src/CMakeLists.txt:5-19`): any `GCC` other than a conda triplet breaks the first archive. | `GCC=gcc-13 cmake`: "Ranlib: gcc-13-ar", zlib-ng archive fails. |
| B4 | Medium | CI | Nothing automated runs `protal_avx2`, the binary the launcher picks on nearly every CPU; CI never runs the e2e tests, the example, the static build, the launcher or a conda build. The e2e suite cannot test `protal_avx2` as it stands (one test expects `$PROTAL` as the program name). | `just e2e` 102/102; with `PROTAL=build/protal_avx2` 101/102. |

Also noted, Low but cheap to fix or likely to confuse (details in `areas/`):

- **Inputs:** `StripString` reads before the start of an empty line (ASan heap-buffer-overflow on FASTA with a blank or
  `\r`-only line, also in `--build`); lowercase reads never align (0 records, exit 0); R1/R2 names are never compared; a
  cut last record of a plain single-end FASTQ is dropped; SEQ/QUAL lengths are not checked; a BGZF file followed by a
  plain gzip member is rejected; the gzip error's reason is never printed; empty inputs and a 0-byte `unique_kmers.tsv`
  pass without a warning.
- **Outputs:** two runs on one prefix overwrite each other's `.partial` (no lock); SAMs are not fsynced before the rename;
  `statistics.tsv` lacks a final newline; `.profile.log`'s `GeneCov0` is always 0.
- **Database:** folder-mode `--add_model` leaves a co-located bundle stale; `--add_model` skips the compression option
  checks; compression memory is threads × frame size.
- **Long reads:** a rerun reusing a SAM ignores its read-type header (profiled with the wrong model, no warning); the
  options summary prints the global `-a` and `snp_min_af`, not ONT's 0.85 and 0.2; single-hit segments at read ends get
  high MAPQ; genes over ~32 kb would break chunking unnoticed; no unit test of `LongReadAligner`.
- **GTDB scripts:** docs say `model.xml` (it is `model_pe.xml`); design points collide when setups share a read length;
  with NCBI down the download exits 0 with 0 genomes; a held-out clade can take all archaea; disk at r226 is undocumented
  (~170 GB kept, ~85 GB transient, besides the release); there is no end-to-end test of `build_gtdb_database.py`.
- **Simulator and taxonomy:** repeated `--taxon` keeps the last; `-rs` in `--extra_art_args` gives every genome the same
  seed; `-l`/`-f` overrides are not reflected in the manifest; duplicate genome names share one length; genome names are
  used as paths (`../` escapes); `--plot_png` runs a shell command built from the output paths (`$(...)` executes);
  headerless genome tables can lose their first row; unknown ranks leave lineage fields empty; `/` or `;` in names break
  file names or lineages; `protal_profile_utils` sums taxa with equal lineages and fails on a glob that catches
  `.profile.log`; `protal_map_utils generate` drops single-end and long-read files.
- **Build and docs:** the qcmsa lookup order is documented the wrong way round (`installation.md:98-100`,
  `running.md:170`; `PROTAL_QCMSA_SCRIPT` comes first); cPMML forces `-Ofast -flto` in every build type; a plain build runs
  cPMML's benchmarks and copies 65 MB into each build tree; global include directories cause the case-insensitive
  `options.h` clash; `--version` names no variant or commit; 456 warnings in own code under `-Wall -Wextra`.
- UBSan reports a `memcpy` from a null pointer on every index load (`IndexCodec.h:203`), which stops any sanitizer run with
  `halt_on_error=1`.

## Verified to hold

- **Database file:** compressing and decompressing round-trip byte for byte and are deterministic; every writer uses
  `.partial` + verify + rename; forged member names (`../`, absolute) are rejected; a flipped byte, a cut seek table, bad
  versions and overlapping members fail with clear messages; 4.2 M mutated index chunks decode under ASan/UBSan without a
  crash.
- **I/O:** SAM records are identical for 1 and 2 threads and for sam, sam.gz and sam.zst, and so are profiles and strain
  outputs; `@SQ` lists exactly the aligned genes; reruns reuse each format; `--profile_only` gives the same profile; cut,
  corrupt or empty SAMs of each format are rejected with a reason; `kill -9` during alignment and full disks leave nothing
  behind; the gene arena gives identical results from every database form; the parallel table parsers report the right
  line numbers; 46 I/O unit tests pass under ASan/UBSan.
- **Long reads:** about 33,000 records checked against the truth: clips, SEQ/QUAL orientation, flags and positions (≥98%
  of HiFi and ≥95% of ONT bases at their true position); every whole gene on reads of 65,000 bp to 200 kb found once;
  profiles within ~0.02; FASTA input identical to FASTQ; mixed pe/pb/ont maps use each sample's `-a` and model; strain
  MSAs from HiFi recover 938/942 true SNPs.
- **GTDB scripts:** a full build runs end to end offline; the PMML model equals scikit-learn 1.9.0 exactly on up to 63k
  rows, including values at every threshold; no label leakage; holdouts are deterministic and leak-free; the collector
  resumes correctly after a kill; the downloader refuses unsafe archive members.
- **Simulator:** reproducible across seeds and thread counts; read pairs sum exactly to the target; its map runs in protal
  (TP 3, FP 0, FN 0), with abundances within 0.012 of its truth; output files' header and row widths agree.
- **Build:** Release, Debug (asserts on) and static builds work; 180/180 unit tests in Release, Debug and ASan/UBSan;
  baseline, AVX2 and static binaries give byte-identical outputs; the launcher picks and falls back correctly and passes
  exit codes through; the installed layout runs the example (PASS).

## Proposed fixes

| Step | Scope | Findings |
|---|---|---|
| Q | Fail instead of exit 0: check the reader's error on paired input, fail unreadable inputs, open each input once (pipes), require a gzip header or the end after each member, reject taxonomy cycles and duplicate names, check read lengths beyond the first 100 reads (or on each read), warn on empty inputs; `StripString` on empty lines. | F1-F4, T1, T2, LR2, and the Low input items |
| R | Long reads: no x-drop in long-read windows (or `--x_drop 0` by default) and the help text; honour a reused SAM's read type; print the per-type values. | LR1, LR3, LR4 |
| S | GTDB build: copy only the files the converter wrote, poll the background build between stages, kill children on any exit, check tools and imports up front, skip builds whose outputs are complete, fingerprint the collector's dumps, resume downloads within a run and accept a complete `.part`; an end-to-end test on a synthetic release. | H1, G1-G7 |
| T | Database file: cap or reject frames without a stored size; name the bundle when separate files shadow it; `--add_model` in folder mode removes or warns about a stale bundle. | DB1, DB2, Low items |
| U | Simulator truth: record the reads ART wrote (per genome, from its output); make the strain spec exclusive and per strain, or document it; `protal_map_utils` single-end support. | S1, M1, U1, Low items |
| V | Build and CI: `just install <prefix>` in the docs, install rules without the build RUNPATH, CMake's own `ar`, CI running the e2e tests against `protal_avx2` and the installed launcher. | B1-B4, Low items |

## Website

The website is likely out of date or incomplete on: `just install` (B1, if it shows the `prefix=` form), which read
compressions and whether pipes are supported (F1, F2), the qcmsa lookup order ([build](areas/build.md) B10), the map and
profile utilities with single-end maps (U1; [simulator](areas/simulator.md) U3), and the definition of relative abundance
(vertical-coverage share).
