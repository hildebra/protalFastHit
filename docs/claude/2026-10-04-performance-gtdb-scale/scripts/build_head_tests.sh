#!/usr/bin/env bash
# HEAD's (c76e838) tests with the working tree's test_InputValidation.cpp (its bench test), for a before/after bench of the gene table loader.
set -uo pipefail
E=$HOME/perf-gtdb/head
cp /mnt/c/Users/hildebra/Documents/locDev/protal/tests/test_InputValidation.cpp $E/tests/test_InputValidation.cpp
cd $E && nice cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON > $HOME/perf-gtdb/headtests.configure.log 2>&1 && \
  nice cmake --build build --target protal_tests -j 6 > $HOME/perf-gtdb/headtests.build.log 2>&1 && echo "OK head tests" || { echo FAIL; grep -E 'error' -A3 $HOME/perf-gtdb/headtests.build.log | head -20; exit 1; }
echo "== bench on HEAD (serial adds), 30000 taxa x 168 genes"
PROTAL_GENE_TABLE_TAXA=30000 $E/build/tests/protal_tests --gtest_filter='GeneTables.BenchLoadOfLargeTables' 2>&1 | grep -E 'genes in|thread\(s\)|PASSED|FAILED'
