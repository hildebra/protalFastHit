#!/bin/bash
# Clean Release configure + build as docs/installation.md describes (default generator), -j 2.
# Then the unit tests in the same tree, as `just test` does.
set -u
A=~/audit6/build
cd $A/src
rm -rf build
{ time cmake -S . -B build -DCMAKE_BUILD_TYPE=Release ; } > $A/rel_configure.log 2>&1
echo "configure rc=$?" | tee -a $A/rel_configure.log
{ time cmake --build build --target protal protal_avx2 simulate_metagenomes -j 2 ; } > $A/rel_build.log 2>&1
echo "build rc=$?" | tee -a $A/rel_build.log
ls -la build/protal build/protal_avx2 build/simulate_metagenomes | tee -a $A/rel_build.log
./build/protal --version 2>&1 | tail -3 | tee -a $A/rel_build.log
echo "== warnings in release build"; grep -c 'warning:' $A/rel_build.log
# unit tests (as just test / docs/testing.md)
{ time cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON ; } > $A/rel_test_configure.log 2>&1
echo "test configure rc=$?" | tee -a $A/rel_test_configure.log
{ time cmake --build build --target protal_tests -j 2 ; } > $A/rel_test_build.log 2>&1
echo "test build rc=$?" | tee -a $A/rel_test_build.log
{ time ctest --test-dir build --output-on-failure ; } > $A/rel_ctest.log 2>&1
echo "ctest rc=$?" | tee -a $A/rel_ctest.log
tail -15 $A/rel_ctest.log
echo DONE
