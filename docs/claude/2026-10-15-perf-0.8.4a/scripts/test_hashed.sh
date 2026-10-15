#!/usr/bin/env bash
# The hashed key map's checks: rebuild with the scratch copy's files, the mini database (once), the key map unit tests and
# the e2e test of the key maps' records, on cores 0,1,3 (a bench may hold core 2).
set -uo pipefail
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/95e4f938-fd80-46a7-b960-30fecd6a2871/scratchpad/hkeys
W=~/hkeys
CORES=0,1,3
rsync -a --checksum --no-times $SP/mod/src/ $W/src/src/
rsync -a --checksum --no-times $SP/mod/tests/ $W/src/tests/
taskset -c $CORES nice -n 5 ninja -C $W/build -j3 protal protal_tests simulate_metagenomes > $W/ninja_protal2.log 2>&1
rc=$?
echo "ninja exit $rc"
if [ $rc -ne 0 ]; then grep -B2 -A14 "error" $W/ninja_protal2.log | head -60; exit 1; fi
cd $W/src
if [ ! -f data/mini_db/protal_db/database.protal ]; then
  PROTAL=$W/build/protal taskset -c $CORES nice -n 5 bash scripts/mini_db/build_mini_db.sh data/mini_db > $W/mini_db.log 2>&1
  echo "mini db exit $?"; tail -3 $W/mini_db.log
fi
cd $W/build/tests && taskset -c $CORES ./protal_tests --gtest_filter='HashedKeys.*:CompactKeys.*' 2>&1 | tail -4
cd $W/src
PROTAL_TEST_DB=data/mini_db/protal_db PROTAL=$W/build/protal SIMULATE=$W/build/simulate_metagenomes PROTAL_TESTS_REQUIRED=1 \
  taskset -c $CORES nice -n 5 python3 -m unittest -v tests/e2e/test_protal_e2e.py -k test_compact_keys_give_the_same_records 2>&1 | tail -15
echo TESTS_DONE
