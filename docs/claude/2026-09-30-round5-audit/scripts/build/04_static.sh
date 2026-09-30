#!/bin/bash
# Static targets in the Release tree (as `just static`), -j 2; what they link against.
set -u
A=~/audit6/build
cd $A/src
{ time taskset -c 0,1 cmake --build build --target protal_static simulate_metagenomes_static -j 2 ; } > $A/static_build.log 2>&1
echo "static build rc=$?"
grep -iE 'warning|error' $A/static_build.log | sort | uniq -c | head -20
ls -la build/protal_0.6.0_static build/simulate_metagenomes_static
file build/protal_0.6.0_static build/simulate_metagenomes_static
ldd build/protal_0.6.0_static 2>&1 | head -3
./build/protal_0.6.0_static --version; echo "rc=$?"
./build/simulate_metagenomes_static --help > /dev/null 2>&1; echo "simulate --help rc=$?"
echo "== link line of protal_static"
cat build/CMakeFiles/protal_static.dir/link.txt
echo "== ISA of static binary: any AVX (ymm) instructions outside zlib-ng/ScanWindowsAvx2?"
objdump -d --no-show-raw-insn build/protal_0.6.0_static | grep -c 'ymm'
echo "== the 'all' target (cmake --build build with no --target): does it build?"
{ time taskset -c 0,1 cmake --build build -j 2 ; } > $A/all_build.log 2>&1; echo "all rc=$?"; tail -5 $A/all_build.log
grep -E 'Built target|Benchmarking' $A/all_build.log | sort | head -40
echo "== cmake --install (protal has no install rules?)"
rm -rf $A/cminstall; cmake --install build --prefix $A/cminstall > $A/cminstall.log 2>&1; echo "install rc=$?"
cat $A/cminstall.log; find $A/cminstall -type f | head
echo DONE
