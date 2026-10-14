# Why `main.cpp.o` takes minutes to compile, and how splitting it would help

Date: 2026-10-14. Branch `audit-fixes` at `f98254a` (built from `git archive`, without the working tree's uncommitted
edits). The user asked why compiling `main.cpp.o` is so slow, and whether parts of it could be compiled separately
and reused while they are unchanged. The study changed nothing in `src/`; its split is a prototype in a WSL build
tree. The plan was then implemented the same day: [Implementation](#implementation) (on `5369c99`).

**Machine.** WSL2 Ubuntu 24.04 (GCC 13.3, Ninja, CMake 3.28) on a Core Ultra 7 258V, shared with other sessions:
load average 4-8 during the first runs, ~2 during the single-core runs. Every job was pinned to cores 0-3 (one core
each where timed alone). Wall-clock times on a loaded machine are up to 3x the quiet ones, so the tables say which
condition they come from, and the share of the work each part takes is measured with GCC's own counters
(`-ftime-report`'s memory allocated in "phase opt and generate"), which do not change with load.

**Build.** `cmake -DCMAKE_BUILD_TYPE=Release` (`-O3 -DNDEBUG`, `-march=x86-64 -mtune=generic -ffp-contract=off`,
`-fopenmp`, C++20), no compiler cache. Scripts: [scripts/](scripts/) (each says what it measures; `setup.sh` makes
the tree in `~/ctime` from a `git archive` tar).

## Summary

- **`main.cpp.o` is all of protal in one translation unit.** `main.cpp` includes `RunProtal.h`, which includes every
  other header: 296,908 preprocessed lines (9.1 MB), 4,452 functions, 3.9 MB of code. One compiler process does all
  of it, so `-j` cannot spread it over cores, and an edit to any header recompiles all of it (a compiler cache does
  not help either). The libraries next to it (`protal_lib`'s eight `.cpp` files, cPMML, WFA2, zlib-ng; 49 s on 4
  cores under load) compile in parallel with it and are rebuilt only when their own sources change.
- **The time is optimisation, spread over the whole program.** Quiet machine, one core: 114 s wall (103 s CPU, 3.0
  GB). Parsing and template instantiation take 14 s (12%); optimisation and code generation 89 s, spread over GCC's
  passes (none above ~5%). It is the amount of code at `-O3` that costs, not one bad construct. Under the load of the
  first runs the same compile took 236-347 s.
- **Not the cause:** the x86-64-v3 clones (`PROTAL_CLONE_V3`; 0.9% of the work), `--benchmark_alignment`'s second
  copy of the alignment templates (2% of the code), `-O3` against `-O2` (`-O2` does 10% less of the optimisation work,
  8% less CPU time), the standard library and intrinsics headers (2 s to parse).
- **Where the work is.** Cutting each stage out of `protal::Run` (GCC removes what nothing calls): alignment with the
  index build 30% of the optimisation work, profiling 6%, strains 6%, the database maintenance modes 1.4%; what is left
  with all four cut (`Run` itself, option parsing, loading the database and tables) 48%. Option parsing alone is ~10%
  of the code (`CxxOptions`, `PrepareAndCheckValidity`, `OptionsFromArguments`, `ToString`).
- **A split by stage works** (prototype: seven translation units, each including all of `RunProtal.h`, nothing else
  changed; it links and runs). But each unit takes 41-70 s alone, 369 s in all (3.2 times today's 114 s), because
  each one parses all of protal again and compiles again the header code it shares with the others. On 4 cores a full
  rebuild would take ~103 s instead of 114 s; with 7+ cores 70 s. An edit to one stage's own `.cpp` would recompile
  41-70 s instead of 114 s; an edit to any header still recompiles all seven.
- **Lean units are what pays.** The strain stage as a unit that includes only what it uses (`Options.h`,
  `Profiler.h`): **20 s** (43 s as a unit of the naive split, and today hidden in the 114 s). In the naive split
  12-22% of a unit's code (in the two units checked) is code it does not use at all: functions with `target_clones`
  (`PrepareMAPQ`, `DecodeChunk`, `GenomeLoader::GetGenome` and the rest) are compiled, with everything they call,
  wherever their callers are visible,
  and `RunProtal.h`'s unused `static` functions make them visible everywhere. In a unit that does not see those
  callers this cost disappears (the lean strains unit compiles the same with and without the clones).
- **Linking: 4.4-5.2 s of LTO for cPMML.** `lib/cPMML/CMakeLists.txt` builds cPMML with `-Ofast -flto`, so every link
  of `protal`, `simulate_metagenomes` and `protal_tests` re-optimises cPMML's code, serially. With cPMML built without
  `-flto` the link takes 0.12-0.15 s; with `-flto=auto` on the link 1.7-2.0 s.

**Recommendation.** (1) Build cPMML without `-flto` (−4.5 s on every link; check that profiling's model scoring is as
fast). (2) Split `RunProtal.h` into stage `.cpp` files that include only what they use, with option parsing and the
loaders as units of their own (plan below): estimated 45-55 s for a full rebuild on 4 cores (~30 s on 8) and 20-35 s
after an edit within one stage, against 114 s + 5 s link today. A split without the header work gains little on 4
cores and is only the first step. Edits to `Options.h` (the `Options` class, which every unit uses) would still
recompile everything, in parallel.

## How the time was measured

### The compile as it is

`scripts/phases.sh` (loaded machine, load 6-8):

| step | wall | peak RSS |
|---|---|---|
| front end only (`-fsyntax-only`) | 35.4 s | 1.6 GB |
| `-O0` | 91.1 s | 2.6 GB |
| `-O3` (the build's command) | 237.8 s | 3.0 GB |

`-ftime-report` of the `-O3` compile: parsing 16.2 s, deferred template instantiation 13.1 s, "opt and generate"
206.3 s (88%); within it "callgraph functions expansion" (per-function optimisation and code generation) 166 s and the
IPA passes 34 s. The largest single passes: operand scan 8.2 s, register allocation 8.0 s, alias walking 7.4 s,
dominator optimisation 6.8 s. The same compile alone on a quiet core (`scripts/sequential.sh`): 113.8 s wall, 102.6 s
CPU, parsing 6.1 s, instantiation 8.0 s, optimisation and code generation 88.5 s.

Lines of the preprocessed unit by file (top): the AVX-512 intrinsics headers ~47,000 (cheap: declarations),
`Profiler.h` 5,830, `RunProtal.h` 4,153, `Options.h` 4,122, `Build.h` 3,228, `Seedmap.h` 2,681, `cxxopts.hpp`
2,089. Front end of the main headers alone (`scripts/frontend.sh`, `-fsyntax-only`, four at a time, load ~2): the standard library with
`<immintrin.h>`, `<omp.h>` and `cxxopts.hpp` 2.2 s, `Options.h` 4.7 s, `Profiler.h` 7.0 s, `Classify.h` 7.9 s,
`Build.h` 8.9 s, `Strain.h` 5.0 s, all of `RunProtal.h` 13.6 s. (`GenomeLoader.h` and `Seedmap.h` do not compile on
their own: `Constants.h` uses `std::cout` and `uint32_t` without including `<iostream>` and `<cstdint>`.)

Code in the `-O3` object by part of protal (`scripts/symbols.sh`, function sizes by name; rough): alignment 17%,
profiling 14%, the build 12%, options 9%, strains 6%, standard library instantiations 18%, the rest 25%.

### Each stage's share (`scripts/variants.sh`)

`RunProtal.h` patched so that each stage's call in `protal::Run` can be left out (`-DNO_ALIGN` and so on; the four
functions of `RunProtal.h` that are not `inline` or `static` made `static`, so that GCC drops them when nothing calls
them). Four compiles at a time, one core each, so each batch ran under the same load:

| variant | CPU s | opt-phase memory | share of the full program's |
|---|---|---|---|
| full `main.cpp` | 261 | 3,873 MB | 100% |
| without alignment and index build | 189 | 2,704 MB | −30% |
| without profiling | 246 | 3,631 MB | −6% |
| without strains | 251 | 3,650 MB | −6% |
| without the maintenance modes (`--compress_db`, `--add_model`, ...) | 206* | 3,818 MB | −1.4% |
| all four cut: `Run`, options, loading | 120* | 1,865 MB | 48% left |
| `-DPROTAL_NO_CLONES` | 208* | 3,840 MB | −0.9% |
| `-O2` | 192* | 3,493 MB | −10% |

\* second batch, a lighter load: CPU seconds comparable within a batch only (the second batch's full compile would
have taken ~208 s, as `-DPROTAL_NO_CLONES` did).

The stage shares add up to 43%, the core to 48%; the other 8% is code that two or more stages use but `Run` does not.
In the core: options 23% of its code, loading the genomes, the database bundle and zstd 18%, the index 8%, the tables
9%, profiler code that `Run` calls itself (model loading, `--taxon_statistics`) 9%, `Run` 5.5%.

### Prototype: a translation unit per stage (`scripts/split.sh`)

`main.cpp` keeps `protal::Run`; each stage moves behind a function declared in [scripts/SplitDecls.h](scripts/SplitDecls.h)
and defined in a unit of its own that includes the patched `RunProtal.h`: alignment (`RunWrapper`), index build
(`build::Run` and the bundle), option parsing (`OptionsFromArguments`, `ToString`, `PrintHelp`), profiling
(`ProfileWrapper`, `ProfilingAhead`), strains (`StrainWrapper2`), maintenance modes. The seven objects link with the
build's own link command (no duplicate symbols: everything else in the headers is `inline`, a template or `static`)
and the binary runs. Not checked: whether its outputs equal the current binary's (nothing in the code changed, only
where it is compiled).

Each unit alone on one quiet core (`scripts/sequential.sh`):

| unit | wall | CPU | front end | opt and generate | text |
|---|---|---|---|---|---|
| today: `main.cpp` | 113.8 s | 102.6 s | 14.2 s | 88.5 s | 4.0 MB |
| `main` (`Run`, loading) | 55.6 s | 48.1 s | 13.6 s | 34.5 s | 1.2 MB |
| index build | 70.4 s | 61.6 s | 16.2 s | 45.5 s | 1.3 MB |
| alignment | 55.3 s | 47.8 s | 14.5 s | 33.3 s | 1.4 MB |
| options | 53.1 s | 44.7 s | 20.5 s | 24.3 s | 1.0 MB |
| maintenance modes | 49.7 s | 42.8 s | 14.7 s | 28.1 s | 1.1 MB |
| strains | 43.3 s | 36.9 s | 13.0 s | 24.0 s | 0.7 MB |
| profiling | 41.6 s | 35.3 s | 14.8 s | 20.5 s | 1.0 MB |
| sum of the seven | 369.0 s | 317.3 s | | | 7.6 MB |

On 4 cores the seven units take ~103 s (longest first: 70; 56 + 42; 55 + 43; 53 + 50), on 7 or more 70 s; the link
is 9.1 s for seven objects against 5.6 s for one (most of both is cPMML's LTO, below). The linked program has 17% more
code (4.89 to 5.71 MB): `RunProtal.h`'s `static` helpers are copied into every unit that uses them, and inlining
differs.

Why every unit is at least ~40 s, when the strain stage is 6% of the work (`scripts/duplicates.sh`,
`scripts/clones_per_unit.sh`):

- each unit parses all of protal again (13-20 s);
- 311 functions (1.5 MB of the 7.4 MB of functions) are compiled in all seven units, and 1,443 more in two to six:
  the header code the stages share (`GenomeLoader`, zstd, the tables, `Options`' inline members, `std::regex`) and
  code that only one unit uses. The second kind comes from `target_clones`: in the maintenance unit nothing calls
  `Profiler::PrepareMAPQ`, yet it is compiled there with `NoteRecord`, the ancestry-site decoder, `GetGenome` and the
  rest of its callees. Its callers (the profiler's SAM parsing, reached from `RunProtal.h`'s `ProfileSample`) are in
  every unit of the naive split, unused; GCC drops them but keeps the cloned function, apparently because the
  dispatcher it makes for a clone stays once a call to it was seen (in GCC's call graph of the unit the dispatcher
  is the root that holds the clones). Compiled with `-DPROTAL_NO_CLONES`, the strain
  unit loses `PrepareMAPQ`, `DecodeChunk`, `NoteRecord` and `ancestry::Decode` and 22% of its code (706 to 551 KB),
  the maintenance unit `PrepareMAPQ` and its callees (it keeps the `DecodeChunk` it uses). A unit in which no unused
  caller is visible does not have the problem, so it is not a reason to drop the clones (a small test with unused
  cloned functions, `scripts/clone_test.sh`, has GCC drop them all).

### A lean unit (`scripts/lean_strains.sh`)

The strain stage (`RunProtal.h` lines 2403-3734, `SharedSNP` to `StrainWrapper2`, with `MakeFolder`, `MSAItemFor` and
`MSAReleaseFor`) in a unit that includes `Options.h`, `Profiler.h` and `ProgressBar.h` instead of `RunProtal.h`:

| | wall | front end | code |
|---|---|---|---|
| strains in the naive split | 43.3 s | 13.0 s | 706 KB |
| strains, lean | 20.0 s | 7.8 s | 549 KB |
| strains, lean, `-DPROTAL_NO_CLONES` | 19.6 s | | 548 KB |

It compiled with the four missing helpers added and nothing else changed: the strain stage depends on little outside
its own lines.

### Linking (`scripts/link.sh`, `scripts/link_no_lto.sh`)

`lib/cPMML/CMakeLists.txt` sets `-std=c++11 -Wall -Ofast -flto` for cPMML, so `libcPMML.a` holds GCC's intermediate
code (2,188 LTO sections) and every link runs the link-time optimiser over it, serially ("using serial compilation of
5 LTRANS jobs"). protal's own objects are not LTO, so this optimises cPMML only within itself.

| link of protal (one `main.cpp.o`, 4 cores, three runs) | wall |
|---|---|
| as now | 4.4-5.2 s |
| `-flto=auto` added to the link (LTRANS jobs in parallel) | 1.7-2.0 s |
| cPMML compiled with `-fno-lto` instead | 0.12-0.15 s |

The same applies to every `simulate_metagenomes` and `protal_tests` link. Without `-flto` cPMML's model scoring may get
slower; it was not measured.

### Smaller points

- `protal` and `protal_static` (and the two `simulate_metagenomes`) compile the same source each: a build of all
  targets compiles `main.cpp` twice. The recipes in the `justfile` and the conda build name their targets (`protal`,
  `simulate_metagenomes`), so only a plain `cmake --build` or `just static` does that. An object library shared by
  both would avoid it, but with `PROTAL_STATIC_FETCH_DEPS` the static binary compiles against the downloaded zstd and
  ISA-L headers, so the two compiles are not always the same.
- No compiler cache is set up by CMake. With one translation unit it would only help an unchanged tree; with several
  units it would also reuse the unchanged ones across branches and the many build trees on this machine.
- A precompiled header could save only the front end (12% today); with lean units the front end is 5-9 s per unit.

## Plan for the split

1. **cPMML without `-flto`** (one line in `lib/cPMML/CMakeLists.txt`), and a check that the profiling stage's time
   is unchanged. Independent of the rest.
2. **Small fixes the split needs:** `Constants.h` includes `<iostream>` and `<cstdint>`; the four functions of
   `RunProtal.h` defined without `inline` or `static` (`ProfileWrapper` twice, `ExtractTaxa`, `GetTotalSystemMemory`)
   get one of them.
3. **Stage files.** A small header (e.g. `RunStages.h`) with `ProtalDB`, the type aliases (`Profiles`,
   `ReadTypeModels`) and one entry function per stage; `RunProtal.h`'s sections move into `.cpp` files that include
   only what they use:
   - `Run.cpp`: `protal::Run`, the start-up; `main.cpp` stays as it is or merges into it;
   - `Load.cpp`: the index and table loaders (`LoadIndex`, `Load*`, `ApplyBuildSettings`), the second-largest part of
     today's core;
   - `OptionsParse.cpp`: the bodies of `CxxOptions`, `OptionsFromArguments`, `PrepareAndCheckValidity`, `ToString`,
     `LoadFromMap`, `PrintHelp` out of `Options.h` (the class keeps the declarations and the accessors);
   - `Align.cpp`: `RunWrapper`'s alignment (it includes `Classify.h`, not `Profiler.h` or `Build.h`);
   - `BuildIndex.cpp`: `RunWrapper`'s build branch and `Build.h`;
   - `Profile.cpp`: `ProfileSample`, `ProfileWrapper`, `ProfilingAhead`, `--taxon_statistics`, model loading
     (`LoadModel`; or `TaxonFilterForest` out of `Profiler.h` into a header of its own, so that `Run.cpp` need not
     include `Profiler.h`);
   - `Strains.cpp`: the strain section (as the lean unit above);
   - `DbMaintenance.cpp`: `--compress_db`, `--decompress_db`, `--unpack_db`, `--add_model`, `--add_tables`.
4. **Header trimming that decides how much each edit recompiles:** `Options.h` includes `StrainAlleles.h`,
   `ColumnWeights.h`, `AncestrySites.h`, `GeneConservation.h` and seven other table headers mostly for their file
   names (`kFileName`): a header of file names instead keeps them out of every unit.

Expected (estimates from the measurements above, except the strain unit): every unit 15-35 s, the largest `Load.cpp`
or `Align.cpp`; a full rebuild 45-55 s on 4 cores, ~30 s on 8, plus the link; an edit within one stage (its `.cpp`, or
a header only it includes, e.g. `Profiler.h`, `Strain.h` or the alignment headers) 20-35 s. An edit to the `Options`
class (a new option) recompiles every unit, in parallel. Since 2026-09-26, 159 of 231 commits touched `RunProtal.h`,
`Options.h`, `main.cpp` or a CMake file; after the split most of `RunProtal.h`'s edits would fall within one stage
file.

Run time is not expected to change: the hot code (seeding, alignment, SAM parsing) stays inline in the headers of
the unit that uses it, and the stages call each other once per run. Points to check when doing it: outputs identical
(the test suites and an end-to-end comparison), callgrind instruction counts of the alignment and profiling benches,
and `static` functions with function-local state (only `RunQCMSA`'s `reported`, which stays in the strain unit). The
move rewrites most of `RunProtal.h` and part of `Options.h`, which other sessions edit daily: it is best done in one
commit while no other session has uncommitted changes in those files.

## Implementation

Done the same day on branch `split-units` (a worktree, while another session had uncommitted changes in
`RunProtal.h` and `Options.h`), first on `ab6bf25`, then rebased onto `5369c99` once that session had committed: its
changes to `RunProtal.h` (`LoadUniqueCores`, the unique cores in the alignment and the start-up) were ported into the
stage files by hand, and `Options.h` was split again from its new version by the same script. Scripts:
[scripts/implementation/](scripts/implementation/).

**What changed.**

- `src/Run/`: `RunProtal.h`'s sections moved into `Run.cpp`, `Load.cpp`, `Align.cpp`, `BuildStage.cpp`, `Profile.cpp`
  and `Strains.cpp`, verbatim by line range (`make_split.sh` checks each range's first and last line before it cuts).
  What changed in the code: the calls between the stages (`RunAlignment`, `RunBuild`, `RunDbMaintenance`,
  `StartIndexLoading`, `MakeProfilingAhead`, `SubmitProfilingAhead`, `TaxonomyGenera`, `ModelKindOfFile`), the
  functions another unit calls declared there instead of `static`, with their default arguments in the declarations,
  and `ProtalDB`'s three index members defined in `Load.cpp`, so that only the units that load or use the index
  include `Seedmap.h`. `RunStages.h` and `ProfileStage.h` declare what the units share; `ProtalDB.h` holds the database
  without the index's header. `RunProtal.h` declares `protal::Run`; `main.cpp` is unchanged.
- `src/Options.cpp` (in `protal_lib`, which the tests link): `CxxOptions`, `OptionsFromArguments`,
  `PrepareAndCheckValidity`, `LoadFromMap`, `ToString` and `PrintHelp` (`make_options2.sh`: found by their signatures,
  moved four spaces less indented; `PrintMapHelp` stays, its raw string would change). `Options.h` no longer includes
  `cxxopts.hpp` (4,226 to 2,401 lines).
- `src/Utilities/MemoryUse.h`: `PrintMemory`, `ReleaseFreeMemory` and `HumanBytes`, from `Build.h`, so that a run's
  units need not include `Build.h`. `Constants.h` includes `<iostream>` and `<cstdint>`.
- `lib/cPMML/CMakeLists.txt`: `-Ofast` without `-flto`.
- `docs/development.md`: a section on the source layout, which unit holds what.

What each unit includes of protal's large headers (`unit_includes.sh`): `Run.cpp`, `Profile.cpp` and `Strains.cpp`
`Profiler.h` and `Options.h`, not `Seedmap.h`, `Classify.h` or `Build.h`; `Load.cpp` `Seedmap.h`; `Align.cpp`
`Classify.h`; `BuildStage.cpp` `Build.h`; `Options.cpp` `Options.h` and cxxopts; `main.cpp` only `<iostream>`.

**Checks.** That the moved code is all there: the lines of the old `RunProtal.h`, `Options.h` and `Build.h` against
those of the new files, as multisets with whitespace stripped, differ only in the edited signatures, the changed
calls and three includes no longer needed. Built at `5369c99` (`final_check.sh`): the unit tests pass (587); the mini
database (`build_mini_db.sh`'s steps, `--build` by each binary from the same converted folder) is the same byte for
byte, `database.protal` included; the end-to-end tests pass with each binary on its database (182), and their kept
outputs (`compare_outputs.sh`, 44 folders, 3,511 files) are the same but for: SAMs whose records come in another order
(119; the alignment threads write them as they finish), run-time tables, paths, and the files `SamInputTest` cuts from
a SAM (`broken.sam`, `misfit.sam` and the profile of the latter: other records, as the order differs). The baseline
run twice against itself differs in the same files (on `ab6bf25`, where the split was compared the same way).

**Run time.** Callgrind instruction counts of one run of each binary (`runtime_check.sh`: one thread, a map of the
e2e tests' paired-end `sa`, PacBio `la` and ONT `oa` samples on their database with every read type's model, profiled,
no strains); the profiles, their logs and compositions identical:

| part | old | split | |
|---|---|---|---|
| alignment (`RunWrapper`, all three samples) | 1,857.4 M | 1,857.9 M | +0.03% |
| of it `RunPairedEnd` / `RunLongReads` | 465.5 M / 1,382.7 M | 467.0 M / 1,381.8 M | +0.3% / −0.07% |
| profiling (`ProfileWrapper`) | 512.5 M | 508.3 M | −0.8% |
| loading the index (`LoadIndex`) | 438.2 M | 432.4 M | −1.3% |
| loading the models (`LoadModel`: cPMML parsing the PMML) | 3,890.1 M | 3,650.8 M | −6.2% |
| the whole run | 6,996.0 M | 6,746.0 M | −3.6% |

protal's own code runs as before: the hot code is compiled within the unit that uses it, and its instruction counts
stay within 1%. cPMML without `-flto` parses the mini database's models with 6% fewer instructions; model scoring is
part of the profiling's line.

**Compile times.** Each unit alone on one core (`unit_times.sh`, load 1-3, two rounds; seconds wall):

| unit | round 1 | round 2 | front end | code |
|---|---|---|---|---|
| before: `main.cpp` | 125.1 | 96.8 | 15.2 | 4.07 MB |
| `Run/BuildStage.cpp` | 36.7 | 36.2 | 6.2 | 1.39 MB |
| `Run/Align.cpp` | 29.6 | 28.3 | 6.2 | 1.09 MB |
| `Run/Run.cpp` | 28.5 | 24.9 | 8.5 | 0.76 MB |
| `Run/Profile.cpp` | 28.0 | 28.1 | 6.4 | 0.92 MB |
| `Run/Strains.cpp` | 20.0 | 25.7 | 6.5 | 0.48 MB |
| `Options.cpp` | 17.8 | 20.7 | 3.6 | 0.78 MB |
| `Run/Load.cpp` | 16.1 | 15.4 | 5.1 | 0.43 MB |
| `main.cpp` | 0.3 | 0.4 | 0.2 | 0 |
| the eight together | 197.0 | 195.8 | | 5.85 MB |

(The same compile varies by up to 30% here between runs: Windows moves the VM's cores between fast and low-power
ones.) The units together take 1.6 times the old one's work (the naive split 3.2 times): what they share of the
headers is parsed and compiled in each, 44% more code that the linker keeps once. The largest unit takes 37 s; on 8
cores or more protal's own code compiles in about that, plus the link.

**Rebuilds.** `ninja -j4 protal` on cores 0-3 after an edit (the file touched), the old tree (`5369c99`) and the
split taking turns (`rebuild_times.sh`), two rounds; seconds wall (in brackets the CPU seconds of round 1). Other
sessions loaded the machine to 2-12 during these runs (the table's last column of `~/split-times.tsv`), so single
numbers vary by up to 2×; the pairs ran next to each other.

| edit | units compiled again | old | split |
|---|---|---|---|
| none: from scratch, the libraries too | all | 105, 157 (156) | 86, 90 (306) |
| `Profiling/Profiler.h` | Run, Profile, Strains | 107, 144 (104) | 36, 62 (90) |
| `Core/ChainAnchorFinder.h` (alignment) | Align | 133, 137 (129) | 36, 29 (35) |
| `Build.h` | BuildStage | 168, 103 (152) | 108, 36 (68) |
| strain stage code (old: `RunProtal.h`) | Strains | 187, 150 (170) | 29, 28 (27) |
| start-up code (old: `RunProtal.h`) | Run | 148, 157 (141) | 24, 33 (23) |
| `main.cpp` (the version history) | main | 112, 166 (109) | 0.6, 0.7 (0.3) |
| `Options.h` | all, and the library's files that include it | 108, 149 (104) | 62, 111 (215) |
| `SequenceUtils/GenomeLoader.h` | all but main and Options | 96, 124 (94) | 102, 53 (200) |

An edit within one stage, its code or a header only it includes, rebuilds in 24-36 s on a quiet machine instead of
100-190 s. An edit that every unit sees (`Options.h`, `GenomeLoader.h`, `Zstd.h`, ...) costs about twice the CPU of
the old build, spread over the four cores: 53-111 s here against 96-149 s, and on more cores about the 37 s of the
largest unit. From scratch on 4 cores the libraries compete with the units: 86-90 s against 105-157 s. The link is
part of each number: under half a second now (the `main.cpp` row: its compile and the link), 4.4-5.2 s before
(cPMML's LTO).

The binary's code is 12% larger (4.97 to 5.55 MB of text): what the linker cannot merge, each unit's own copies of
functions with internal linkage (the x86-64-v3 clones of `DecodeChunk` in three units, header functions declared
`static`, sort loops: 0.16 MB more than in the one unit), and each unit's own inlining and specialised copies.
