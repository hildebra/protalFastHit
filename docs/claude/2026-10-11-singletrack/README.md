# Singletrack in protal: the option, the database's record, and a SAM identity check

- **Date**: 2026-10-11.
- **Code**: branch `singletrack` (worktree `../protal-singletrack`) on `7d5df4b`, uncommitted when written. The
  reference build is `7d5df4b` itself (WFA2-lib v2.3.6).
- **Machine**: WSL2 on the Core Ultra 7 laptop. Every job ran on cores 0-3, niced; other sessions kept the load at 4-9.
  Wall-clock numbers are medians of alternated runs, instruction and cache counts from callgrind.
- **Data**:
  - Database: `~/realworld_s1/protal_db` (six real GTDB r226 families, 530 genomes; built 2026-10-10 at `6f8e6bb`, so
    without `build_settings.tsv`).
  - Samples: simulated by `simulate_metagenomes` from 40 of its genomes (representatives and strains):
    - paired-end, 300k and 3M pairs;
    - PacBio HiFi (`hifi:15000:3000:3`), 30 and 300 Mb;
    - Nanopore (`qshmm:QSHMM-ONT-HQ`), 30 and 150 Mb.
- **Scripts and outputs**: `scripts/`, `results/`. The builds, samples and runs were in WSL `~/stk` (removed after).
- **Follow-up of**: the [literature search](../2026-10-11-alignment-algorithms-literature.md), which found Singletrack
  the one drop-in speed-up of WFA2.
- **Task (the user)**: try Singletrack with a SAM identity check; add it as an option to protal; warn users if the
  database was built with another algorithm.

## Summary

1. **The option.** `--wfa_traceback auto|classic|singletrack`, default `auto`.
   - `--build` records its traceback (`auto`: classic) in a new database member, `build_settings.tsv`. It also aligns
     its strain alleles, congener gaps and column weights with it.
   - A run with `auto` aligns as its database was built. A database without the file was built classic.
   - A run told otherwise is warned before it aligns.
   - `build_gtdb_database.py --wfa-traceback` passes it to both builds. The collection's runs follow the training
     database (auto), so the models learn from alignments made as the finished database's runs make them.
2. **The library.** `lib/wfa2-lib` is now T. Mokveld's WFA2-lib fork (`85d063f`: v2.3.6, upstream's later fixes,
   Singletrack for ends-free alignment and heuristics), with three changes of protal's. As it came, the fork had two
   bugs, both triggered by WFA2's wf-adaptive heuristic, which protal always uses:
   - **Singletrack ended the process.** The backtrace called `exit(-1)` on a PacBio flank with the read's end free (the
     first PacBio run died).
   - **Classic CIGARs could be wrong.** The classic traceback gave operations that are no alignment of the sequences (an
     M of C against A).
   - **Scores changed after pruning.** A forward-pass change made heuristic-pruned alignments score differently from
     v2.3.6.

   Section 2 has the changes.
3. **Classic is v2.3.6, byte for byte.**
   - All six samples gave the same SAM records and every other output file (but `misc/cpu.tsv`, CPU times).
   - So did 33,500 synthetic protal-like alignments (status, score, operations).
   - So the default changes nothing, and the models stay valid.
4. **Singletrack gives the same scores, in a few alignments other operations.**
   - Same outcome and score in all 33,500 synthetic alignments.
   - Its backtrace never failed in the runs (no fallback to classic).
   - Reads with other records:

     | Sample | Reads with other records | Of them, another MAPQ |
     |---|---|---|
     | Paired-end | 0.009-0.013% | under 0.01% |
     | PacBio | 0.2-0.6% | |
     | Nanopore | 47-49% | 4-5% |

   - Profiles: identical for the paired-end and PacBio samples. ONT was not profiled (the database has no ONT model).
5. **Speed.**
   - **PacBio**: aligning 13% faster (handler −23%), every one of five rounds.
   - **Nanopore**: aligning 10% faster.
   - **Paired-end**: no measurable change; WFA2 is 0.8% of a short-read run's instructions here.
   - **Inside WFA2 on long reads**: −7% instructions and −70 to −87% L1 data misses.
   - **Peak RSS**: 30-40 MB lower.
6. **Recommendation.**
   - Keep `classic` the default.
   - Singletrack pays for long reads, but changes half of the Nanopore records. Use it with models trained on a database
     built with it (`build_gtdb_database.py --wfa-traceback singletrack`).
   - Measure it at r226 first, where WFA2 is 91% of long-read aligning
     ([2026-10-04](../2026-10-04-performance-gtdb-scale/README.md)).
   - The two fork bugs are worth reporting to `tmokveld/WFA2-lib` (the reproducers are in `scripts/`); that is for you
     to decide.

## 1. What protal does now

| Piece | Where | What |
|---|---|---|
| traceback choice | `src/Alignment/WfaTraceback.h` | `WfaTraceback {Classic, Singletrack}`, names, `ResolveWfaTraceback(given, database)` (given, else the database's, else classic) |
| aligner | `src/Alignment/WFA2Wrapper2.h` | a traceback per aligner (copies keep it); Singletrack = WFA2's `MemorySingletrack`. An alignment whose Singletrack backtrace finds no path is aligned again classic, same penalties and X-drop (`SingletrackFallbacks()` counts them; the run log prints the count) |
| record | `src/Utilities/BuildSettings.h` | `build_settings.tsv`: `name<TAB>value` lines, `#` comments; unknown names skipped (later versions), a bad value refused |
| option | `src/Options.h` | `--wfa_traceback` (exit 2 on another value), `WfaTracebackGiven()`, `GetWfaTraceback()`, `BuildSettingsDbFile()`; the file in `--no_bundle`/`--unpack_db`'s lists |
| build | `src/Build.h` | `WriteBuildSettings` first in the build; the congener gaps, strain alleles and column weights align with the build's traceback; `build_settings.tsv` packed into `database.protal` |
| run | `src/RunProtal.h` | `SetWfaTraceback` before the gene tables: reads the record, resolves, logs "WFA2 traceback: X; the database was built with Y", warns if X ≠ Y (not with `--profile_only`) |
| long reads | `src/Alignment/AnchoredAlignment.h` | the pieces' aligner traces back as the window's |
| GTDB builds | `scripts/build_gtdb_database.py` | `--wfa-traceback {classic,singletrack}` into both `protal --build` commands and, when given, both build keys (earlier builds' keys unchanged) |

The warning, for a database built classic and a run with `--wfa_traceback singletrack`:

```
Warning: this run aligns with the singletrack WFA2 traceback (--wfa_traceback singletrack), but the database was built
with classic (it has no build_settings.tsv: built before 2026-10-11). The two give the same alignment scores, but among
equally good alignments they may choose differently, so the SAM records, and the features the database's models were
trained on, can differ slightly from a run that aligns as the database was built. Leave out --wfa_traceback (auto) to
align as it was built.
```

## 2. The library: the fork, and three changes

`lib/wfa2-lib` is `github.com/tmokveld/WFA2-lib` at `85d063f` (2026-10-05), without `img/`, `tests/`, `tools/`,
`examples/` and `.github/`. The fork holds v2.3.6, upstream's fixes since (static bands, BiWFA, CIGAR buffers, #126) and
Singletrack (`wavefront_memory_singletrack`) for end-to-end, ends-free and extension alignments with heuristics. The
source list is unchanged. `lib/wfa2-lib.cmake` records the source and the changes.

**The mechanism behind all three.** wf-adaptive narrows a wavefront's `lo`/`hi` after it is computed, but leaves the
pruned offsets in the array. The compute kernels of a later wavefront read `[lo-1, hi+1]` of their inputs, and an
ends-free wavefront's limits can be wider than its inputs', so they read past `hi+1` into the pruned offsets
(`wavefront_compute_init_ends` nulls only beyond the initialized range). A state can thus descend from a pruned one.
v2.3.6 is consistent about it: its classic backtrace looks up M offsets in the whole initialized range, as the kernels
read them.

| # | File | Fork | protal's change | Why |
|---|---|---|---|---|
| 1 | `wavefront_heuristic.c` | `wf_heuristic_nullify_boundaries`: nulls `lo-1`/`hi+1` of a pruned wavefront | off (`#if 0`) | with it, heuristic-pruned alignments differed from v2.3.6: 180 of 33,500 synthetic ones, 71 scoring worse, all flanks given a budget of 2000 (none within the 0.15 budget) |
| 2 | `wavefront_backtrace.c` | classic M lookups (`wavefront_backtrace_m_offset`) only within `lo`/`hi` ("Use active M-wavefront bounds", `caecb60`) | the whole initialized range, as v2.3.6 | the classic traceback took a wrong predecessor for a state from a pruned one: operations that are no alignment (repro3: 152 × 218 bp flank, an M of C against A at op 151) |
| 3 | `wavefront_backtrace.c`, `.h`, `wavefront_unialign.c` | Singletrack's M-only backtrace looks only within `lo`/`hi`; no path: `exit(-1)` | a cell without a path within `lo`/`hi` is tried again in the initialized range; a backtrace still without a path returns false, and the alignment fails as `WF_STATUS_UNATTAINABLE` | the first PacBio Singletrack run died: a 309 × 544 bp flank, the read's end free and the reference's start not, a 3-base deletion opened from a pruned M state (repro: 36 of 280 configurations of free ends, X-drop and budget exited) |

The retry in change 3 changes nothing that worked before: a cell is looked up in the padding only after the active
bounds gave no path. `WFA2Wrapper2` turns a failed backtrace into a classic alignment of the same pair. After changes 1-3
no run and no synthetic alignment needed that.

## 3. Identity: classic against v2.3.6

`stk_check.sh identity`, 4 threads, the reference build against the work build (`auto` = classic on this database):

| Sample | SAM records | SAM | Other files |
|---|---|---|---|
| paired-end, 300k pairs | 26,568 | identical | identical |
| PacBio, 30 Mb | 1,805 | identical | identical |
| Nanopore, 30 Mb | 1,656 | identical | identical |
| paired-end, 3M pairs | 264,584 | identical | identical |
| PacBio, 300 Mb | 18,385 | identical | identical |
| Nanopore, 150 Mb | 8,458 | identical | identical |

"Other files" means every output file but `misc/cpu.tsv` (CPU times): profiles, strain tables, logs of evidence.

`wfadiff.sh` gave 0 differing lines (status, score, operations) between the reference library and the final one, over
33,500 synthetic alignments (`results/wfadiff.txt`):
- 20,000 short reads in gene windows with dovetails, X-drop 1000, budget 61;
- 1,500 long-read gene windows, a fifth unrelated;
- 6,000 long-read flanks with the read's end free, budgets 2000 or 0.15 of the read;
- 6,000 end-to-end pieces.

## 4. Singletrack against classic

**Synthetic alignments.** The same outcome and score in all 33,500 (`wfadiff.txt`), every operation string an alignment
of the sequences. Other operations in:

| Kind | Completed | Other operations |
|---|---|---|
| short reads | 12,810 | 561 (4.4%) |
| long-read windows | 1,008 | 568 (56%) |
| flanks | 5,266 | 2,786 (53%) |
| end-to-end pieces | 6,000 | 932 (16%) |

**Whole runs** (`results/identity_final.out`):

| Sample | Reads whose records differ | What differs |
|---|---|---|
| paired-end, 300k pairs | 2 of 15,283 (0.013%) | CIGAR |
| paired-end, 3M pairs | 14 of 152,201 (0.009%) | 9 CIGAR, 5 MAPQ |
| PacBio, 30 Mb | 2 of 893 (0.22%) | CIGAR |
| PacBio, 300 Mb | 53 of 9,096 (0.58%) | 48 CIGAR, 4 CIGAR + `ZA`, 1 MAPQ |
| Nanopore, 30 Mb | 447 of 910 (49%) | 389 CIGAR, 44 with MAPQ, 13 with `ZT`/`ZU`/SEQ or `ZA` |
| Nanopore, 150 Mb | 2,190 of 4,683 (47%) | 1,925 CIGAR, ~200 with MAPQ, a few POS, `ZR`, `ZT`/`ZU`, one RNAME, one FLAG |

**Profiles.** The profiles and strain outputs of the paired-end and PacBio samples were identical. The Nanopore runs had
`--no_profile`: this database has no ONT model.

**What the differences are.** They are tie-breaks:
- Where a gap sits within a repeat, or a mismatch beside a gap: equally good alignments, picked otherwise.
- Nanopore reads have the most indels, so they have the most ties.
- A tie-break changes a MAPQ when the two alignments have another number of mismatches against gap bases, which the
  identity and the allele scores count.

**Why the warning matters.** A model trained on runs of one traceback sees slightly other features in runs of the
other: hardly on short reads, noticeably on Nanopore.

## 5. Speed

**Wall clock** (`stk_check.sh timing`): 4 threads, `--no_profile`, classic and Singletrack alternated, 5 rounds,
medians.

| Sample | Aligning (stage), classic → Singletrack | Alignment handler (mean per thread) | Peak RSS |
|---|---|---|---|
| paired-end, 3M pairs | 5.47 → 5.05 s (noise: rounds 4.5-6.4 s) | 0.39 → 0.46 s (noise) | same |
| PacBio, 300 Mb | 2.13 → 1.85 s (−13%) | 1.96 → 1.52 s (−23%) | 3.31 → 3.27 GB |
| Nanopore, 150 Mb | 1.49 → 1.34 s (−10%) | 1.36 → 1.32 s (−3%) | 3.28 → 3.25 GB |

PacBio was faster with Singletrack in every round, and so was Nanopore's stage. On this small database the paired-end
handler is under a tenth of the stage, so its rounds are noise.

**Callgrind** (`stk_check.sh callgrind`): one thread, cache simulation, subsets of 50k pairs, 250 HiFi reads and 400
Nanopore reads. WFA2 inclusive (`wavefront_align`), classic → Singletrack:

| Subset | Instructions | L1 data read misses | L1 data write misses | LL write misses | WFA2 share of the run |
|---|---|---|---|---|---|
| paired-end | 48.1 → 46.9 M (−2.5%) | 109k → 50k (−54%) | 101k → 54k (−47%) | 1.2k → 0.6k | 0.8% |
| PacBio | 139.4 → 129.0 M (−7.5%) | 1.10M → 0.15M (−87%) | 498k → 143k (−71%) | 84k → 30k | 3.0% |
| Nanopore | 281.6 → 263.4 M (−6.5%) | 1.82M → 0.29M (−84%) | 959k → 304k (−68%) | 93k → 31k | 5.4% |

**What it means for protal.**
- Singletrack touches far less memory: the indel wavefronts live in a small recycled scope. That is where the
  long-read wall-clock gain comes from, more than from the instructions.
- These subsets understate WFA2's weight. At r226 scale WFA2 is 31% of a short-read run's alignment and 91% of a long
  read's ([2026-10-06](../2026-10-06-performance-profiling/README.md),
  [2026-10-04](../2026-10-04-performance-gtdb-scale/README.md)).
- So the PacBio handler's −23% is the better guide for r226 long reads. For short reads expect little: the WFA2
  instructions fall by 2.5%, and the flanks are short.

## 6. Tests

- **Unit tests**: 507 pass, 3 skipped (the AVX-512, debug-build and gene-window ones, skipped in the reference build
  too).
  - New: `WfaTraceback.*` (scores against classic on short reads, long windows, read-end-free flanks and end-to-end
    pieces; the failing PacBio flank; copies; names and resolution) and `BuildSettings.*`.
  - `WFA2Wrapper.AnAlignmentXDropCutsShortFails` runs again: with the unpatched fork it was skipped, as that library
    no longer cut the alignment short.
- **End-to-end**: `WfaTracebackTest`, 4 tests, all pass.
  - A synthetic release built with `--wfa_traceback singletrack`: the record is unpacked and checked, a run with
    `auto` follows the database without a warning, and a run with `classic` is warned.
  - A database whose `build_settings.tsv` was removed counts as classic.
  - An unknown value exits 2.
- **The whole e2e suite**, on a fresh mini database built by the work build (`scripts/stk_e2e_fresh.sh`, as
  `just mini-db`):
  - the work build: 154 tests OK;
  - the reference build: its 150 tests OK on the same database, so an older protal reads a database with
    `build_settings.tsv` (an unknown member to it) as before.

  On a borrowed older mini database that lacks `model_pe.xml` and `species_priors.tsv`, both builds failed the same
  29 tests (`stk_e2e_ref.sh`).
- **Build script**: `scripts/mini_db/test_gtdb_build.py` `BuildOptionsTest` passes, with the new
  `test_the_wfa_traceback_reaches_both_builds`.

## 7. Files

- **Changed**: `lib/wfa2-lib/` (the fork, changes 1-3), `lib/wfa2-lib.cmake`, `src/Alignment/WFA2Wrapper2.h`,
  `src/Alignment/AnchoredAlignment.h`, `src/SequenceUtils/CongenerGaps.h`, `src/SequenceUtils/ColumnWeightsBuild.h`,
  `src/Build.h`, `src/Options.h`, `src/RunProtal.h`, `src/Utilities/Database.h` (comment), `src/CMakeLists.txt`,
  `scripts/build_gtdb_database.py`.
- **New**: `src/Alignment/WfaTraceback.h`, `src/Utilities/BuildSettings.h`.
- **Tests**: `tests/test_WFA2Wrapper.cpp`, `tests/e2e/test_protal_e2e.py`, `scripts/mini_db/test_gtdb_build.py`.
- **Docs**: `docs/running.md` (option), `docs/databases.md` (file table, build, build script option),
  `docs/installation.md` (the fork).
- **Website**: [protal.earlham.ac.uk](https://protal.earlham.ac.uk/main.php?site=documentation) does not list
  `--wfa_traceback` or `build_settings.tsv`.

**How to rerun.**
1. `scripts/stk_build.sh` builds `~/stk/ref` (the archive of `7d5df4b`) and `~/stk/work` (with the worktree's diff as
   `st.patch`).
2. `scripts/stk_sim.sh` and `scripts/stk_sim_big.sh` make the samples.
3. `scripts/stk_check.sh unit|identity|timing|callgrind` runs the checks (`SETS="pe3m pbbig ontbig"` for the big
   samples).
4. `scripts/wfadiff.sh` runs the synthetic comparison; `repro*.sh` and `nonull.sh` hold the cases above.
