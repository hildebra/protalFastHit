#!/bin/bash
# Sync the instrumented Seedmap.h and keybench.cpp into ~/sparseprof and build protal + keybench on cores 4,5.
set -e
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/cf462e00-4c32-4801-b980-f96668bbdc3a/scratchpad
T=~/sparseprof
cd $T
[ -f bin/protal_base ] || { mkdir -p bin; cp build/protal bin/protal_base; }
rsync -a --checksum --no-times $SP/mod/src/ src/src/
mkdir -p src/misc
rsync -a --checksum --no-times $SP/keybench.cpp src/misc/keybench.cpp
if ! grep -q "add_executable(keybench" src/CMakeLists.txt; then
  cat >> src/CMakeLists.txt <<'EOF2'

# keybench: experiment only (docs/claude, sparse load and key map profiling), not part of protal.
add_executable(keybench misc/keybench.cpp)
target_compile_options(keybench PRIVATE ${PROTAL_COMMON_COMPILE_OPTIONS})
target_link_libraries(keybench PRIVATE protal_lib dl isa_baseline)
if (TARGET OpenMP::OpenMP_CXX)
    target_link_libraries(keybench PRIVATE OpenMP::OpenMP_CXX)
endif()
EOF2
fi
if ! grep -q "add_executable(keybench_v3" src/CMakeLists.txt; then
  cat >> src/CMakeLists.txt <<'EOF2'

# keybench_v3: the same compiled for x86-64-v3 (popcnt inline), as protal's target_clones paths are.
add_executable(keybench_v3 misc/keybench.cpp)
target_compile_options(keybench_v3 PRIVATE ${PROTAL_COMMON_COMPILE_OPTIONS} -march=x86-64-v3 -mtune=generic -ffp-contract=off)
target_link_libraries(keybench_v3 PRIVATE protal_lib dl)
if (TARGET OpenMP::OpenMP_CXX)
    target_link_libraries(keybench_v3 PRIVATE OpenMP::OpenMP_CXX)
endif()
EOF2
fi
taskset -c 4,5 nice -n 5 ninja -C build -j2 protal keybench keybench_v3 > ninja_bench.log 2>&1 || { grep -B2 -A12 "error" ninja_bench.log | head -80; exit 1; }
grep -c "march=x86-64-v3" build/build.ninja
tail -2 ninja_bench.log
ls -la build/protal build/keybench bin/protal_base
