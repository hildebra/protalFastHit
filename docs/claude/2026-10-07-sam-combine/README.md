# Combining the SAMs of several runs into one strain analysis

2026-10-07, branch `audit-fixes` at `9ecf7ee` (protal 0.7.8). Binary: `~/testaudit/src/build/protal`, built from a
source tree identical to `git archive 9ecf7ee` (`diff -rq` of `src/`, `scripts/`, `CMakeLists.txt`), copied to
`~/samcombine/bin`. Database: the mini database (`~/ti/mini_db/protal_db/database.protal`, three species: Mockella alpha
and beta, Fakibacter gamma). Reads: the e2e tests' simulator (`scripts/simulate.py`), 12 pairs on every gene of
320 bp or more. WSL, `taskset -c 0-3`, `-t 4`, `--no_qcmsa`. Work folder `~/samcombine/work` (not kept in git).

**Question**: audit the path that combines `.sam.zst` files of several samples (from several runs) into one strain
analysis. Is it implemented effectively, are there bugs, can it be done dynamically within a folder or across folders
(wildcards), and is it documented?

## Summary

1. **The core works and is exact.** Strain MSAs are built over all samples of one protal invocation, and
   `--profile_only a.sam.zst,b.sam.zst,...` (or a map whose `SAM` column names existing files) makes such an
   invocation from SAMs of separate runs. The raw MSAs, partitions and `.meta.tsv` of two studies' SAMs combined so
   were byte-identical to those of one joint run over the same reads (path A, F), and so were the profiles. No index is
   loaded (`--profile_only`, or a map whose SAMs all exist); at r226 profiling was 3.7 s of a 42 s paired-end run
   ([2026-10-06](../2026-10-06-performance-profiling/README.md)), so re-profiling is about a tenth of a full run.
2. **Bug (crash): every run without `-o` aborts** (`std::filesystem::create_directories("")` throws, exit 134):
   `--profile_only` after profiling every sample, at the strain stage, although its help and `docs/running.md` say the
   outputs then go next to each SAM; `-1/-2` runs right after loading the index, also with `--no_strains`. With `-o .`
   the same run works. The three call sites are as old as `alpha`.
3. **`protal_map_utils merge` does not reuse the runs' SAMs.** It rewrites the `SAM` column to the bare file name under
   the merged `#OUTPUT_DIR/alignments`, so protal aligns every sample again (index loaded, new SAMs), and fails when the
   read files are gone. It is a "rerun all together" tool, not a "combine alignments" tool, and the docs do not say so.
4. **Not dynamic.** `--profile_only` takes a comma-separated list only: a folder fails (`does not end with .sam`,
   exit 35), a quoted glob fails (`sam file does not exist: run*/alignments/*.sam.zst`). The shell can expand one
   (`"$(ls run*/alignments/*.sam.zst | paste -sd,)"`), which works across folders. Sample IDs are the SAM file stems;
   two folders' `sa.sam.zst` are refused with a clear message, and `--prefix` (one per SAM, all or none) or a map's
   `#SAMPLEID` renames them.
5. **Documentation**: nowhere (docs/ or website) says that MSAs span only the samples of one invocation or how to
   combine runs; `--profile_only` is described as a quick way to try another knob or model. Details below.

## What was tried

Three studies were aligned and profiled on their own with maps from `protal_map_utils generate` (`run1`: sa, sb;
`run2`: sc, sd; `run3`: another sample named sa), and sa..sd once more in one joint run as the reference.
`scripts/run.sh` runs everything, `scripts/no_outdir.sh` the runs without `-o`; `results/run_summary.txt` and
`results/no_outdir.txt` are their outputs.

| | Path | Result |
|---|---|---|
| A | `--profile_only "$(ls run1/alignments/*.sam.zst run2/alignments/*.sam.zst \| paste -sd,)" -o comb_po` | exit 0; 3 MSAs with rows sa sb sc sd; all 9 strain files and the 4 profiles identical to the joint run's |
| A2 | `--profile_only` of the joint run's own SAMs | strain files identical |
| A3 | A again into the same `-o` | exit 0; `misc/sa_runtime.tsv` grows from 9 to 17 rows (the profiling rows appended again) |
| B | `--profile_only` without `-o` | exit 134 at the strain stage, after the profiles were written next to the SAMs (`run1/alignments/sa.profile`, ...) |
| C | `run1/.../sa.sam.zst,run3/.../sa.sam.zst,...` | exit 30: shared profile file and sample ID `sa`, both reported; with `--prefix run1_sa,run3_sa,sc` exit 0, rows run1_sa run3_sa sc |
| D | `--profile_only run1/alignments`; `--profile_only 'run*/alignments/*.sam.zst'` | exit 35 (`does not end with .sam, .sam.gz or .sam.zst`); exit 30 (`does not exist`) |
| E | `protal_map_utils merge --map run1.map run2.map --out comb_merge` | the merged map (`results/merged.map`) names `sa.sam.zst` under `comb_merge/alignments`: protal loaded the index and aligned all four again (MSAs identical to the joint run's); with the reads moved away, exit 30 (`-1 file does not exist`). Without `--out` the SAMs go to `work/alignments`, again not the runs' |
| F | a map of absolute SAM paths, `FIRST`/`SECOND` naming files that are gone (`results/sams.map`) | exit 0, "All alignments are present", 8 warnings (`-1 file does not exist ... ( But .sam file does )`); strain files identical to the joint run's |
| F2 | the same with `#SAMPLEID` study1_sa, study3_sa, study2_sc and `-` as `FIRST`/`SECOND` | exit 0; rows named by `#SAMPLEID`; three warnings that the SAMs hold paired-end reads "profiled as such, not as single-end reads" (`-` in `SECOND` means single-end) |
| G | a map of two existing SAMs and one new sample's reads (`results/grow.map`) | exit 0; the two SAMs skipped, the new sample aligned; rows sa sb sc |
| H | `--profile_only` of SAMs in a read-only folder | exit 0, no warning, no `.err` file written |
| L | `--map run1.map --profile_only anything` | exit 0: profiles the map's SAMs in place (two warnings that the read lists are cleared); undocumented |
| N | `protal_map_utils generate` on empty files | `x_R1.fq/x_R2.fq` → sample `x_R`; `x.R1.fastq.gz` → `x.R`; `x_S1_L001_R1_001.fastq.gz` → `x_S1_L001_R1_001`; `.fq.zst` files are not found |

`no_outdir.sh`: `--profile_only` of one SAM, `-1/-2` with one and two samples, and two samples with `--no_strains`
all abort with `filesystem error: cannot create directories: Invalid argument []`; with `-o .` exit 0.

## Findings

1. **Runs without `-o` abort** (B; high: a `--profile_only` run loses its strain stage after all the profiling).
   `create_directories` of an empty parent path throws in libstdc++: `src/RunProtal.h:531` (the SAM's folder:
   `--prefix sa` without `-o` gives `sa.sam.zst`), `src/RunProtal.h:1130` (the profile's folder; `exists("")` is false),
   `src/RunProtal.h:2555` (`GetOutputDir()`, empty for `-1/-2` and `--profile_only` without `-o`). Fix: treat an empty
   folder as the current one at all three (or default `output_dir` to `.` where the prefix logic does not need it
   empty), and add e2e runs without `-o`. Every e2e `--profile_only` call passes `-o`, and all but the SAM-input
   tests (`tests/e2e/test_protal_e2e.py:2430`, one usable SAM each) pass `--no_strains` (`:331`), so no test checks an
   MSA built by `--profile_only`.
2. **`merge` cannot combine existing alignments** (`scripts/protal_map_utils:532` replaces each SAM by its base name,
   `--use-sampleid` by `<id>.sam.zst`). Either keep a row's SAM as an absolute path when the file exists (then the
   merged map reuses it, path F shows protal handles that), or say in its help and the docs that a merged map aligns
   everything again.
3. **The `.err` file goes next to the input SAM** (`src/RunProtal.h:994`), not into `-o`: combining SAMs rewrites the
   original runs' `.err` files (A: its mtime changed), two concurrent combinations over the same SAMs write the same
   file, and a read-only folder silently gets none (H). Write it next to the profile or into `misc/`, and warn if it
   cannot be written.
4. **A map must name read files even when only its SAMs are used** (`src/Options.h:1700` needs `FIRST` and `PREFIX`);
   placeholders work but print two warnings per sample, and `-` as `SECOND` turns the sample into single-end until the
   SAM header overrides it (F2). Allow `FIRST` to be `-` (or the column to be absent) when the `SAM` exists.
5. **Repeated `--profile_only` into one `-o` appends** the profiling rows to `misc/<sample>_runtime.tsv` again
   (`src/RunProtal.h:915`, A3). Minor.
6. **Memory grows with the cohort.** Every sample's read ranges (`ReadInfo`, 24 bytes per read and gene) and variants of
   the species that enter MSAs stay in memory until the strain stage (`src/RunProtal.h:1181`,
   `Profiler.h:1313`). Not measured here (the mini database is too small to show it); an estimate for 1,000 samples
   with a million reads each on MSA species is tens of GB. Adding one sample to a study also re-profiles every SAM.
   Both would go away with a per-sample strain evidence file written at profiling (variants and read ranges of the
   MSA species), which a later run merges without the SAMs.
7. **`protal_map_utils generate` names samples badly** for the commonest read names (N): `remove_pair`
   (`scripts/protal_map_utils:253`) strips only a trailing `1`/`2` and the separators before it, so `_R1` leaves `_R`
   and Illumina's `_R1_001` keeps the whole name; `READ_EXTENSIONS` (`:31`) has no `.fq.zst`/`.fastq.zst`, which
   protal reads since `e7b391e`. These names become `#SAMPLEID`, the MSA row names.

## Doing it dynamically

What works today, within one folder or across folders:

```bash
protal --db DB -t 16 -o combined --profile_only "$(ls study*/alignments/*.sam.zst | paste -sd,)"
```

The sample IDs are the file stems; when two folders hold the same name, give `--prefix` for every SAM (e.g. built
from the folder names), or write a map with `#SAMPLEID`, absolute `SAM` paths and placeholder `FIRST`/`SECOND`
(path F2). Paths with commas cannot be listed.

What would make it dynamic in protal itself (not implemented):
- `--profile_only` accepting folders (their `*.sam`, `*.sam.gz`, `*.sam.zst`, without `.partial` files) and glob
  patterns (POSIX `glob(3)`), sorted for a fixed sample order, plus `@file` for a list of SAMs with an optional
  sample-ID column (no comma limit, no ID collisions);
- on an ID collision, an option to name samples `<folder>_<stem>` instead of failing;
- or a `protal_map_utils` subcommand that writes a map of existing SAMs from globs, with IDs from the folder where
  names collide (`generate --id-from-folder` does that for reads).

## Documentation

- `docs/strains.md` says a species needs two samples "so a run of one sample writes no MSAs", but not that MSAs span
  only one invocation's samples, nor how to combine runs. It needs a section: profile the runs' SAMs together with
  `--profile_only ... -o` (same database: protal checks the `@SQ` genes) or a map of their SAMs; sample IDs from file
  stems, `--prefix` for collisions; per-sample profiles are rewritten into `-o`; no index is loaded.
- `docs/running.md:170` and the `--profile_only` help (`src/Options.h:93`) promise outputs "next to each SAM" without
  `-o`, which aborts (finding 1); neither says that the listed SAMs form one strain analysis, where `.err` goes, or
  that `--map` with `--profile_only` re-profiles a map's SAMs.
- `protal_map_utils merge` is described only as "merge maps" (`docs/installation.md:121`, `docs/development.md:259`,
  `docs/running.md:37`): that the merged map aligns again into its new folder is not said.
- `--map_help` (`src/Options.h:1458`) still says `SAM` defaults to `<PREFIX>.sam`; three lines above it says
  `.sam.zst`, which is right.
- The website (not in this repository): its "Strain resolution" section does not say which samples enter an MSA or that
  MSAs are per run, says nothing about combining runs, lists `--profile_only` only as profiling without re-aligning,
  and says the `SAM` column defaults to `{OUTPUT_DIR}/{PREFIX}.sam` (now `.sam.zst`). It does not mention
  `protal_map_utils`.

## Reproducing

```bash
# in WSL; SRC: git archive of the commit, PROTAL: its protal, DB: the mini database
SRC=~/samcombine/head PROTAL=~/samcombine/bin/protal DB=~/samcombine/db/database.protal W=~/samcombine/work \
    bash docs/claude/2026-10-07-sam-combine/scripts/run.sh
bash docs/claude/2026-10-07-sam-combine/scripts/no_outdir.sh     # needs run.sh's work folder
bash docs/claude/2026-10-07-sam-combine/scripts/collect.sh       # copies the summaries into results/
```

## Follow-up: the fixes (2026-10-07)

On `9ecf7ee` plus the uncommitted changes below, built in WSL from a copy of the working tree
(`scripts/build.sh`, Release, unit tests on). The user asked to fix the `-o` crash and `merge`, document
combining runs, add wildcards for existing SAMs, and refuse map sample names that cannot be written as
file names on Linux and macOS.

- **Runs without `-o`** (finding 1): `MakeFolder` (`src/RunProtal.h`) takes an empty folder for the
  current one at the three call sites; the strain stage makes the strain output folder it writes to.
- **Wildcards** (`src/Options.h`, `ExpandSamFiles`): a `--profile_only` item with `*`, `?` or `[...]`
  that names no file stands for the `.sam`, `.sam.gz` and `.sam.zst` files it matches (`glob(3)`,
  sorted, other files passed over, with a note of the count); a pattern without a SAM and a folder are
  errors, reported with the other input errors. Arguments that follow no option are taken as SAM files
  after `--profile_only` (an unquoted pattern, which protal before silently cut to its first file, see
  `scripts/unquoted_glob.sh`), and stop protal anywhere else (exit 2; an unquoted `-1 *_1.fq` aligned
  the first file only), except one `--build` reference.
- **Sample names** (`SamSampleNames`): file names, or where they repeat, the folders in which the
  paths differ (`S1`, `S2` for a folder per sample; `study1_sa`, `study2_sa` for two studies), with a
  note. This replaces the error of finding C, which `--prefix` still overrides.
- **Unsafe sample IDs** (`SampleIdProblem`): `/`, `:`, control characters, a leading `-`, `.` and
  `..`, invalid UTF-8 (macOS file names must be UTF-8) and more than 200 bytes are refused: in a map as
  it is read, every such row with its line (exit 9), and for `--prefix` and `--profile_only` names
  with the other checks. Whitespace stays refused only with strain MSAs, as before.
- **`protal_map_utils merge`** (finding 2): a row whose run wrote its SAM (the `SAM` column's file, or
  `<OUTPUT_DIR>/<PREFIX>.sam*` without the column, which is then added) keeps it by its absolute path;
  `--new-sams` restores the old behaviour. `merge`, `validate` and `flatten` keep `-` in `SECOND`
  (single-end samples, which `merge` turned into a path and refused twice as a duplicate) and accept
  maps without a `SECOND` column; `validate` does not ask for the reads of a row whose SAM exists.
- **Docs**: `docs/strains.md` "Strain MSAs over several runs"; `docs/running.md` (outputs without `-o`,
  `#SAMPLEID` rules, `merge`, `--profile_only`); `--profile_only`, `-o` and `--map_help` texts (the
  stale `<PREFIX>.sam` default fixed); the tool tables in `installation.md` and `development.md`.

Checks (`scripts/check_fixes.sh` on run.sh's work folder, `results/fix_summary.txt`):

| | Run | Result |
|---|---|---|
| 1 | `-1/-2`, two samples, no `-o` | exit 0, outputs in the current folder; strain files identical to run1's |
| 2 | `--profile_only` of four SAMs, no `-o` | exit 0, profiles next to the SAMs, `strains/` in the current folder, identical to the joint run's |
| 3, 4 | `'run[12]/alignments/*.sam.zst'` quoted; the same unquoted | exit 0, rows sa sb sc sd, identical to the joint run's |
| 5 | `'run*/alignments/*.sam.zst'` (sa in run1 and run3) | exit 0, rows run1_sa run1_sb run2_sc run2_sd run3_sa, with a note |
| 6 | `'fix/persample/*/aln.sam.zst'` | exit 0, rows sa sb sc sd, identical to the joint run's |
| 7 | a pattern without SAMs, a folder, `*.err`, a missing file | exit 30, all four reported at once |
| 8 | `-1 a_1.fq b_1.fq` | exit 2, `unexpected argument(s): reads/run1/sb_1.fq` |
| 9 | `merge` of run1.map and run2.map, then `protal --map` with the reads moved away | SAMs by absolute path; `validate` exit 0; "All alignments are present", no index loaded, identical to the joint run's |
| 10 | `merge` of single-end maps | `-` kept; without `SECOND` exit 0 |
| 11 | a map with `a/b`, `c:d`, `-e`, `..`, a control character, invalid UTF-8, 201 bytes | exit 9, all seven listed by line |

Tests: `tests/test_Parsing.cpp` (`Options.SampleIdsThatCannotNameFiles`,
`SamPatternsExpandToTheSamFilesTheyMatch`, `SamSampleNamesFromFilesOrFolders`); `tests/e2e`:
`CombiningRunsTest` (six tests), `MapUtilsTest` (merge keeps SAMs and a merged map aligns nothing,
single-end maps, unsafe sample IDs), and `MateAssignmentTest` now runs without `-o`. Unit suite 399
passed (2 skipped, as before: `scripts/unit.sh`); e2e 134 passed in 111 s (`scripts/e2e.sh`); with the
binary of `9ecf7ee`, 9 of the 12 new and changed e2e tests fail (the other three test the script, or
what worked).

Still open: the `.err` file next to the input SAM (finding 3), repeated `--profile_only` runs appending
to `misc/<sample>_runtime.tsv` (5), a map needing `FIRST` when its SAM exists (4: protal warns,
`validate` accepts it), memory over large cohorts (6), and `generate`'s sample names (7). The website
does not describe combining runs, patterns, runs without `-o` or the sample-ID rules.

## Follow-up 2: reruns, generate's names, and the memory a cohort keeps (2026-10-07)

On `12b059b` (the fixes above, committed and verified from `git archive`: unit tests pass, e2e 134 OK)
plus the changes below, built in a clean tree of that commit with this session's files only
(`scripts/build_on_commit.sh`). Asked: stop `--profile_only` where results exist (the duplicate runtime
rows), name `generate`'s samples without the R1, and find out why memory grows with the samples and
whether it can be limited.

**Reruns.** `--profile_only` now stops before it starts where it would overwrite an earlier run's
profiles or the strain folder's `species.tsv` (`Options::PrepareAndCheckValidity`); `--force` writes
them again. `WriteProfilingTimes` replaces a sample's earlier `Profiling: ...` rows of
`misc/<sample>_runtime.tsv` instead of appending a second set, which also fixes ordinary reruns that
profile an existing SAM again.

**`generate`'s names.** `remove_pair` drops the mate number, an `R` before it, the separator before
that and Illumina's chunk number after it: `x_R1.fq` is `x` (was `x_R`), `x.R1.fastq.gz` is `x` (was
`x.R`), `x_S1_L001_R1_001.fastq.gz` is `x_S1_L001` (was the whole name), `x_R1_trimmed.fq` is
`x_trimmed`; `.fq.zst` and `.fastq.zst` files are found.

### Why memory grows with the samples

Data: `db900n` (the 900-species world of the performance reports), `w900` (1M pairs from 60 species,
~25x each) aligned once, then `--profile_only` over 1-16 hard-linked copies of its SAM, 4 threads
(`scripts/memory.sh`); heap profiles of 2 copies on one thread with valgrind massif (`scripts/massif.sh`).

Each sample's profile stays in memory until the strain MSAs are written. What it kept up to `12b059b`:

- **The strain evidence of the species that enter MSAs**, as the reads left it: a `Variant` of 88 bytes
  per allele at every position with a variant, each with a quality byte and a strand/divergence byte
  per read in two vectors, in a `robin_map` of vectors; a 24-byte `ReadInfo` per read and gene; the
  positions without a base. For `w900`: 59 species, 6,490 genes, 7.1M bases, 349,000 variant positions
  (the simulated strains differ from their references).
- **Every other taxon with reads**, its genes and counters, which only `--taxon_statistics` reads.

Heap (massif, 2 samples): base 103 MB; **159 MB kept per sample** (421 MB after profiling). The RSS
(16 samples, 4 threads): 2.82 GB after profiling, about 163 MB more per sample once more samples than
threads are done; with `--no_strains` 1.12 GB.

Nothing of the strain stage needs the reads themselves: a sample's MSA rows are made from its own records
filtered at its own identity threshold (`MSAItem`), and the other samples set only which genes and
columns an MSA has. Every reader of the kept data reads aggregates of an allele: its observations per
strand, their quality sum, its flags.

### What changed

1. **Kept packed, after the sample's outputs are written** (`Gene::KeepForMSAs`, `PackedMSAItem` in
   `Strain.h`): for each gene of a species that enters MSAs, its `MSAItem` (the parameters it is built
   with are fixed per sample: the identity threshold of the taxon's own reads, the SNP filters of the
   read type), each allele in 32 bytes (kind, position, bases, observations per strand, quality sum,
   flags, INDEL bases apart), the coverage in 1, 2 or 4 bytes per base by its largest value; the three
   SNP counts and the covered bases the strain stage reads. The records are freed. Asking a kept gene for
   anything else throws. `Variant::Restore` makes the alleles again (compact: no per-read bytes; its
   unused `variant_id` holds the quality sum, so a `Variant` does not grow).
2. **Taxa that enter no MSA are dropped** (`MicrobialProfile::ReleaseReadData`), unless
   `--taxon_statistics`.
3. **`malloc_trim(0)`** after each sample (glibc): the freed records' pages go back to the system.

| 16 samples, 4 threads | `12b059b` | 2, items kept as compact `Variant`s | 1 + 2 | 1 + 2 + 3 |
|---|---:|---:|---:|---:|
| RSS after profiling, strains | 2.82 GB | 2.14 GB | 1.44 GB | 1.11 GB |
| peak RSS, strains | 2.84 GB | 2.15 GB | 1.45 GB | 1.42 GB |
| RSS after profiling, `--no_strains` | 1.12 GB | 1.09 GB | 1.09 GB | 0.72 GB |
| heap kept per sample (massif, strains) | 159 MB | | 37 MB | |

(`results/memory/rss_*.tsv`, `heap_*_2.txt`; wall times vary up to 2x on this machine and are not
compared.) The outputs are byte-identical: profiles, raw MSAs, partitions, meta, SNP and coverage
tables of 4 and 16 samples, with and without strains, against `12b059b` (`scripts/compare_runs.sh`,
537 files for 16 samples), and the e2e suite's MSA and long-read phasing tests pass.

The peak is now about the samples profiled at once (each ~150 MB of records while it is read) plus 37 MB
per sample done. Without strains a sample keeps less than massif's 0.5% threshold; the RSS that still
grows a little there is the allocator's, not data (`heap_nostrains_2.txt`: after the first sample's
release the heap holds the database, the model, the shared congener sketches and the next sample's
SAM buffer).

**Further limits, not done.** The packed evidence is plain arrays, so it could be written to a file per
sample and read back per species in the strain stage: memory would then be one species over all samples,
at a disk cost of about what it holds now (37 MB per sample of this kind). Coverage could be run-length
coded (7 MB of the 37 here). At 1,000 dense samples the kept evidence is ~37 GB instead of ~159 GB.

Tests: `PackedMSAItem.UnpacksWhatTheMSAReads` (alleles of every kind, coverage in every width),
`Abundance.ReleasingReadDataKeepsWhatLaterStagesRead` (the kept item and counts equal those from the
records; other parameters throw), `Abundance.ReleasingReadDataDropsTheTaxaNoMsaTakes`; e2e
`CombiningRunsTest.test_earlier_results_stop_protal_unless_forced`,
`MapUtilsTest.test_generate_names_samples_without_mate_numbers`. Unit 401 passed (2 skipped), e2e 136
passed in 58 s.
