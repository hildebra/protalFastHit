#!/bin/bash
# CI's sanitizer job, step by step (Ninja, Debug, ASan+UBSan), -j 2 on 2 cores.
set -u
A=~/audit6/build
cd $A/src
rm -rf build-asan
FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer"
cmake -S . -B build-asan -G Ninja -DCMAKE_BUILD_TYPE=Debug -DPROTAL_BUILD_TESTS=ON \
  -DCMAKE_CXX_FLAGS="$FLAGS" -DCMAKE_C_FLAGS="$FLAGS" -DCMAKE_EXE_LINKER_FLAGS="$FLAGS" > $A/asan_configure.log 2>&1; echo "configure rc=$?"
{ time taskset -c 0,1 cmake --build build-asan --target protal_tests -j 2 ; } > $A/asan_build.log 2>&1; echo "build rc=$?"; grep -E '^real' $A/asan_build.log
{ time ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1 taskset -c 0,1 ctest --test-dir build-asan --output-on-failure ; } > $A/asan_ctest.log 2>&1; echo "ctest rc=$?"
grep -E 'tests passed|tests failed|^real' $A/asan_ctest.log; grep -E 'runtime error|ERROR: AddressSanitizer|LeakSanitizer' $A/asan_ctest.log | sort | uniq -c | head
grep -E '^CXX_FLAGS|FLAGS = ' build-asan/build.ninja 2>/dev/null | grep -m2 cPMML
grep -A3 'build src/cPMML/CMakeFiles/cPMML.dir/src/api/model.cc.o' build-asan/build.ninja | grep FLAGS
echo DONE
