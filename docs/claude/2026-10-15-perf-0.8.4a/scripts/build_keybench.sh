#!/usr/bin/env bash
# Build keybench (with the hashed-lines layout) against HEAD's protal_lib in ~/hkeys, on cores 0-3.
set -uo pipefail
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/95e4f938-fd80-46a7-b960-30fecd6a2871/scratchpad/hkeys
W=~/hkeys
mkdir -p $W/src
cd $W
if [ ! -f src/src/Hash/Seedmap.h ]; then tar -xf $SP/src_head.tar -C $W/src; fi
mkdir -p src/src/misc
rsync -a --checksum --no-times $SP/keybench.cpp src/src/misc/keybench.cpp
if ! grep -q "add_executable(keybench " src/src/CMakeLists.txt; then
  cat >> src/src/CMakeLists.txt <<'EOF2'

# keybench: experiment only (docs/claude, key map layouts), not part of protal.
add_executable(keybench misc/keybench.cpp)
target_compile_options(keybench PRIVATE ${PROTAL_COMMON_COMPILE_OPTIONS})
target_link_libraries(keybench PRIVATE protal_lib dl isa_baseline)
if (TARGET OpenMP::OpenMP_CXX)
    target_link_libraries(keybench PRIVATE OpenMP::OpenMP_CXX)
endif()
add_executable(keybench_v3 misc/keybench.cpp)
target_compile_options(keybench_v3 PRIVATE ${PROTAL_COMMON_COMPILE_OPTIONS} -march=x86-64-v3 -mtune=generic -ffp-contract=off)
target_link_libraries(keybench_v3 PRIVATE protal_lib dl)
if (TARGET OpenMP::OpenMP_CXX)
    target_link_libraries(keybench_v3 PRIVATE OpenMP::OpenMP_CXX)
endif()
EOF2
fi
if [ ! -f $W/build/build.ninja ]; then
  (cd $W/src && cmake -S . -B $W/build -G Ninja -DCMAKE_BUILD_TYPE=Release > $W/configure.log 2>&1) || { tail -30 $W/configure.log; exit 1; }
fi
taskset -c 0-3 nice -n 5 ninja -C $W/build -j4 keybench keybench_v3 > $W/ninja_bench.log 2>&1
rc=$?
echo "ninja exit $rc"
if [ $rc -ne 0 ]; then grep -B2 -A12 "error" $W/ninja_bench.log | head -60; exit 1; fi
tail -2 $W/ninja_bench.log
find $W/build -maxdepth 2 -name "keybench*" -type f -executable | xargs ls -la
echo BUILD_DONE
