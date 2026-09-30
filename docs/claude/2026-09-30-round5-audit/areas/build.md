# Reviewer report: build, packaging, install, CI (round 5, strain-fixes @ 39a8585)

Saved from the reviewer's hand-back (it could not write findings.md). Source copy WSL ~/audit6/build/src (diff -rq empty);
Ubuntu 24.04, GCC 13.3.0, CMake 3.28.3, just 1.21.0, Intel Core Ultra 7 258V (x86-64-v3); -j 2 / taskset -c 0,1; fake nproc.
Scripts ../scripts/build/ (00_env.sh ... 11_asan.sh, run_detached.sh); logs ~/audit6/build/. WSL crashed twice.

| id | sev | file:line | finding | evidence |
|---|---|---|---|---|
| B1 | High | docs/installation.md:76, justfile:261-271 | documented `just install prefix="$HOME/.local"` does not set the prefix: prefix is a positional parameter, gets the literal "prefix=/..."; everything lands in <checkout>/prefix=/.../bin, exit 0; the final check cannot fail ($(prefix=/.../protal --version) = assignment + a command named --version) | 05_install.sh: rc 0, files in ~/audit6/build/src/prefix=/home/falk/audit6/build/prefix/bin; log "bash: line 1: --version: command not found", "Installed to prefix=...: "; `just install $P` works (7 files) |
| B2 | Medium | src/CMakeLists.txt:141-172 (+ lib/gzstream target_link_directories), justfile:262-270 | dynamic binaries, also installed ones, carry a RUNPATH of 13 absolute build-tree dirs (11 nonexistent, from ../IO ../Hash ... ../cPMML relative to src/); just install copies with cp, no install(), RPATH never stripped; loader searches them first for libzstd, libdeflate, libgomp, libstdc++ (planted library would load); not relocatable | readelf -d prefix/bin/protal_baseline: RUNPATH [.../src/src/../IO:...:.../src/../cPMML:.../lib/gzstream]; 11 missing |
| B3 | Medium | src/CMakeLists.txt:5-19 | CMAKE_AR/RANLIB forced to $ENV{GCC}-ar/-ranlib: GCC set to anything but a conda triplet breaks the first archive; without GCC the plain ar works only because Ubuntu auto-loads the LTO plugin cPMML's -flto objects need | GCC=gcc-13 cmake: "Ranlib: gcc-13-ar"; gzstream_lib rc 2 "Error running link command", "[zlib-ng/libz-ng.a] Error 2" |
| B4 | Medium | .github/workflows/ci.yml:22,25; tests/CMakeLists.txt:22; justfile:220-222; tests/e2e/test_protal_e2e.py:1100 vs src/Options.h:799-803 | nothing automated runs protal_avx2 (the launcher's choice on nearly every CPU): CI builds it, runs only baseline-linked unit tests and mini-DB generator tests; CI never runs e2e, example, static build, launcher/installed layout, conda build; e2e cannot validate protal_avx2 as is: protal names itself "protal" when its file is protal_avx2/_baseline, one test expects $PROTAL (so also fails with PROTAL=<prefix>/bin/protal) | 08_e2e.sh: just e2e 102/102; PROTAL=build/protal_avx2 -> 1 failure test_no_preload_on_a_single_file_says_how_to_unpack: 'protal' != '.../protal_avx2'; 101 pass |
| B5 | Low | lib/cPMML/CMakeLists.txt:15-29 | cPMML forces its own build type and -std=c++11 -Wall -Ofast -flto in every build type: Debug and the sanitizer job compile the scorer fast-math, LTO, no -g; every link an LTO link ("using serial compilation of N LTRANS jobs") | Debug flags.make for cPMML -Ofast -flto; ASan tree same plus -fsanitize |
| B6 | Low | lib/cPMML/CMakeLists.txt:177-193, lib/cPMML/benchmark/CMakeLists.txt:1-3 | plain cmake --build builds cPMML's model_benchmark.exe and runs 11 benchmarks; configure copies 65 MB of benchmark data into each build tree; cmake --install installs only libcPMML.a and cPMML.h; protal has no install rules | 04_static.sh: 11 "Benchmarking", "Scoring TPS: -2,147,483,648"; du 65M |
| B7 | Low | CMakeLists.txt:78 | global include_directories(src/ lib/) puts protal's src/ ahead of cPMML's src/: root cause of the documented /mnt/c failure (options.h vs Options.h); target-scoped includes would fix it | ln -s Options.h src/options.h -> cPMML build fails "src/options.h:12:10: fatal error: LineSplitter.h" |
| B8 | Low | CMakeLists.txt:30-31 | CMAKE_CXX_FLAGS_RELEASE/DEBUG overwritten; user's -DCMAKE_CXX_FLAGS_RELEASE ignored | cache -O1 -g, targets -O3 -DNDEBUG |
| B9 | Low | CMakeLists.txt:12-13, 21-24, 34-47; src/CMakeLists.txt:101-139, 258-299; lib/gzstream/CMakeLists.txt:4-5 | leftovers: unused __CHAR_UNSIGNED___EXITCODE; macOS OpenMP paths though Linux-only; configure deletes src/protal_config.h from the source tree; protal_lib_static / gzstream_lib_static duplicates compiled twice; stdc++fs linked; ~10 debug message() lines per configure; gzstream keeps -fopenmp with OpenMP disabled; MY_FLAGS undocumented | 07 (a): "Removed stale generated src/protal_config.h"; 07 (d) gzstream -fopenmp |
| B10 | Low | docs/installation.md:98-100, docs/running.md:170 | both pages give the qcmsa lookup order the wrong way round: code (RunProtal.h:1336-1344, comment correct) checks PROTAL_QCMSA_SCRIPT before the qcmsa next to the binary and $PATH | installed layout + PROTAL_QCMSA_SCRIPT=fake: the fake ran for all 3 species |
| B11 | Low | conda-recipe/meta.yaml:35-44, conda-recipe/build.sh:7-15, 35-38 | recipe tests only presence/--version/--help (the mini-DB generator allows a 30 s functional test); mini_db/*.sh installed mode 644; strain_test/ and db_compression_benchmark.sh not shipped; ${CMAKE_ARGS} not passed; version synced by hand (0.6.0 both today); source path: .. copies the whole tree and the documented --output-folder conda-build lies inside it, not gitignored | code |
| B12 | Low | src/RunProtal.h:1855-1860 | --version prints only "protal v0.6.0": no variant, ISA or commit; docs call bioconda's package 0.6.0a | launcher (AVX2), forced baseline, static binary, run_info.txt all "protal v0.6.0" |
| B13 | Low | justfile:216, 239-256, 262-270 | -j$(nproc) hard-coded; {{prefix}} unquoted (a prefix with a space fails rc 1 after creating stray dirs); just --list shows the wrong comment line for multi-line comments; install skips share/protal/scripts though installation.md:79 says "the same layout as the conda package"; strain recipes default to cluster paths | just install "$A/with space": rc 1 "cp: target 'space/bin/protal_baseline'" |
| B14 | Low | justfile:210-211 | without scikit-learn, just model-test skips all 5 tests and exits 0; CI does not run it | rc 0, skipped x5 |
| B15 | Low | docs/testing.md:162-163 | coverage list names 8 areas; 16 test files built (misses long reads, reader/timers, WFA2 wrapper, syncmers/AVX2, anchored alignment, SamFile) | code |
| B16 | Low | .github/workflows/ci.yml:11, 33 | actions/checkout@v4 pinned by tag, no permissions: block, no timeout-minutes | code |
| B17 | Low | environment.sh:3-4, varkit_protal_speed.txt, .idea/, .vscode/, .gitmodules | tracked leftovers: environment.sh hard-codes /usr/users/QIB_fr017/..., unreferenced; a benchmark note at the root; .idea/ and .vscode/ tracked though gitignored; empty .gitmodules | git ls-files |
| B18 | Info | CMakeLists.txt:105-107, src/CMakeLists.txt:227 | protal_avx2 links main.cpp (x86-64-v3) with libprotal_lib.a (x86-64): 95 inline/COMDAT functions in both, main.cpp.o first so its copies win; benign (identical outputs, AVX2 binary only runs on v3); the 8 library .cpp files stay baseline in the AVX2 binary (speed only) | 09_odr.sh: weak symbols main.cpp.o 1776, lib 287, both 95 |

Warnings under -Wall -Wextra (Debug -O0): 456 unique in own code (450 src/, 6 tests/): -Wsign-compare 169,
-Wunused-function 66, -Wreorder 64, -Wunused-variable 62, -Wunused-parameter 41, -Wunused-but-set-variable 20,
-Wdeprecated-copy 13; worst RunProtal.h 41, Utilities.h 40, Options.h 38, AlignmentUtils.h 36; lib/ 25 more
(~/audit6/build/dbg_warnings_own.txt). Default Release: none in own code (4 cPMML ptr_fun deprecations, LTO note).

Verified to hold: documented Release configure/build (protal, protal_avx2, simulate_metagenomes) 12m18s at -j 2 under load;
Debug -Wall -Wextra builds, protal gets -O0 -g without NDEBUG (round 1's "asserts never fire" fixed); static targets 4m06s,
protal_0.6.0_static 7.1 MB, "not a dynamic executable", only glibc's dlopen notice; dynamic binaries need no system zlib
(NEEDED: zstd, deflate, stdc++, m, gomp, gcc_s, c); ctest 180/180 Release (incl. AVX2 syncmer tests) and Debug; CI sanitizer
job locally (Ninja, Debug, ASan+UBSan) 180/180 in 57 s; mini-DB generator step / just mini-db-test 15/15; just e2e 102/102;
--help/--full_help/--map_help exit 0; baseline, AVX2, static and prebuilt Ninja binaries give byte-identical outputs on the mini
DB (3 samples x 20k pairs, MSAs + qcmsa, 64 files, .zst compared decompressed) at -t 1 except misc/*_runtime.tsv; at -t 2
only SAM record order differs; launcher under dash: picks protal_avx2, PROTAL_NO_AVX2=1 or 0 selects baseline, a failing CPU
check (no avx2, avx2 without abm, grep failing) selects baseline, hidden AVX2 binary -> baseline, none -> message and 127,
symlinked launcher works, exit codes pass through; just install <prefix> installs 7 executables mode 755; a 3-sample run
through the installed launcher finds qcmsa next to it and writes 3 filtered MSAs; examples/mini_db/run.sh via the installed
launcher PASS in 29 s (0.199/0.499/0.302 vs 0.2/0.5/0.3); installation.md build commands, names and the /mnt/c warning
accurate; meta.yaml and CMake both 0.6.0.

Not verified: conda build and environments; a real non-AVX2 CPU (static binary has 5,504 ymm instructions, mostly glibc and
zlib-ng ifunc variants).
