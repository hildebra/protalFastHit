#!/bin/bash
# Clean Debug configure + build with -Wall -Wextra (C++ only), -j 2; unit tests; the flags each
# library directory really gets (cPMML and gzstream set their own per-build-type flags).
set -u
A=~/audit6/build
cd $A/src
rm -rf build-debug
{ time cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug -DPROTAL_BUILD_TESTS=ON -DCMAKE_CXX_FLAGS="-Wall -Wextra" ; } > $A/dbg_configure.log 2>&1
echo "configure rc=$?"
{ time taskset -c 0,1 cmake --build build-debug --target protal protal_avx2 simulate_metagenomes protal_tests -j 2 ; } > $A/dbg_build.log 2>&1
echo "build rc=$?"
{ time taskset -c 0,1 ctest --test-dir build-debug --output-on-failure ; } > $A/dbg_ctest.log 2>&1
echo "ctest rc=$?"; tail -8 $A/dbg_ctest.log
echo "== per-directory flags (Debug tree)"
for f in build-debug/CMakeFiles/protal.dir/flags.make build-debug/CMakeFiles/protal_avx2.dir/flags.make \
         build-debug/src/CMakeFiles/protal_lib.dir/flags.make build-debug/src/cPMML/CMakeFiles/cPMML.dir/flags.make \
         build-debug/src/gzstream/CMakeFiles/gzstream_lib.dir/flags.make build-debug/src/CMakeFiles/wfa_lib.dir/flags.make; do
  echo "-- $f"; grep -E '^(CXX_FLAGS|C_FLAGS|CXX_DEFINES)' $f
done
echo "== Release tree flags"
for f in build/src/cPMML/CMakeFiles/cPMML.dir/flags.make build/src/gzstream/CMakeFiles/gzstream_lib.dir/flags.make build/CMakeFiles/protal_avx2.dir/flags.make build/src/CMakeFiles/protal_lib.dir/flags.make; do
  echo "-- $f"; grep -E '^(CXX_FLAGS|C_FLAGS)' $f
done
echo "== warnings, own code (src/, tests/), unique file:line:col [flag]"
grep -E "^$A/src/(src|tests)/[^:]+:[0-9]+:[0-9]+: warning:" $A/dbg_build.log \
  | sed -E "s#^$A/src/##" | sed -E 's/^([^:]+:[0-9]+:[0-9]+): warning: .*\[(-W[^]]+)\]$/\1 \2/' | sort -u > $A/dbg_warnings_own.txt
wc -l < $A/dbg_warnings_own.txt
echo "-- by flag"; awk '{print $2}' $A/dbg_warnings_own.txt | sort | uniq -c | sort -rn
echo "-- by file"; awk -F: '{print $1}' $A/dbg_warnings_own.txt | sort | uniq -c | sort -rn | head -40
echo "== warnings in lib/ (unique)"; grep -E "^$A/src/lib/[^:]+:[0-9]+:[0-9]+: warning:" $A/dbg_build.log | sort -u | wc -l
echo "== other warning lines (linker, lto)"; grep -iE 'warning' $A/dbg_build.log | grep -vE '^/' | sort | uniq -c | head
echo DONE
