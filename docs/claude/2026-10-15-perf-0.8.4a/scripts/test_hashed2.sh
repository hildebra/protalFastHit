#!/usr/bin/env bash
# The working tree's hashed key map (per-partition tallies, prefetched build, next-line prefetch): a fresh build of
# protal and its tests from the working tree's src/tests (cores 0-3), the mini database, the key map unit tests, the key
# maps' e2e test (next line on, the default, and off), then the whole unit and e2e suites.
set -uo pipefail
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/95e4f938-fd80-46a7-b960-30fecd6a2871/scratchpad/hkeys
W=~/hkeys2
mkdir -p $W/src
cd $W/src
tar -xf $SP/tree_now.tar
if [ ! -f $W/build/build.ninja ]; then
  cmake -S . -B $W/build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON > $W/configure.log 2>&1 || { tail -30 $W/configure.log; exit 1; }
fi
taskset -c 0-3 nice -n 5 ninja -C $W/build -j4 protal protal_tests simulate_metagenomes > $W/ninja.log 2>&1
rc=$?
echo "ninja exit $rc"
if [ $rc -ne 0 ]; then grep -B2 -A14 "error" $W/ninja.log | head -60; exit 1; fi
if [ ! -f data/mini_db/protal_db/database.protal ]; then
  PROTAL=$W/build/protal taskset -c 0-3 nice -n 5 bash scripts/mini_db/build_mini_db.sh data/mini_db > $W/mini_db.log 2>&1
  echo "mini db exit $?"
fi
cd $W/build/tests && taskset -c 0-3 ./protal_tests --gtest_filter='HashedKeys.*:CompactKeys.*' 2>&1 | tail -3
cd $W/src
for nl in 1 0; do
  echo "== PROTAL_HASHED_NEXT_LINE=$nl"
  PROTAL_HASHED_NEXT_LINE=$nl PROTAL_TEST_DB=data/mini_db/protal_db PROTAL=$W/build/protal SIMULATE=$W/build/simulate_metagenomes PROTAL_TESTS_REQUIRED=1 \
    taskset -c 0-3 nice -n 5 python3 -m unittest tests/e2e/test_protal_e2e.py -k test_compact_keys_give_the_same_records 2>&1 | tail -3
done
cd $W/build/tests && taskset -c 0-3 ./protal_tests 2>&1 | grep -E '^\[  PASSED|^\[  FAILED|tests from .* ran'
cd $W/src && PROTAL_TEST_DB=data/mini_db/protal_db PROTAL=$W/build/protal SIMULATE=$W/build/simulate_metagenomes PROTAL_TESTS_REQUIRED=1 \
  taskset -c 0-3 nice -n 5 python3 -m unittest tests/e2e/test_protal_e2e.py 2>&1 | tail -3
echo TESTS2_DONE
