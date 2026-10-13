#!/bin/bash
# Full check, as CI runs it: build (latest tree), the mini database built by examples/mini_db/run.sh, ThreadedGzStream
# tests, ctest, e2e as is and with every gzip input in pieces, the static binary's e2e in pieces, the scripts' suites.
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/0ff4a1d7-3599-4615-a178-e4ca253109af/scratchpad
W=~/pgzbuild
bash $SP/pgz_build.sh 2>&1 | grep -E "exit|FAILED|PASSED|Failure" | head -20
cd $W/src
rm -rf data/mini_db
taskset -c 0-3 nice -n 5 env PROTAL=$W/build/protal THREADS=4 bash examples/mini_db/run.sh data/mini_db > $W/mini_db.log 2>&1
echo "mini_db exit $?"; tail -2 $W/mini_db.log; ls data/mini_db/protal_db | head -20 | tr '\n' ' '; echo
cd $W
taskset -c 0-3 nice -n 5 ctest --test-dir $W/build -j4 > $W/ctest.log 2>&1; echo "ctest exit $?"; grep -E "tests passed|tests failed|\*\*\*Failed" $W/ctest.log | head
cd $W/src
export PROTAL_TEST_DB=$W/src/data/mini_db/protal_db SIMULATE=$W/build/simulate_metagenomes PROTAL_TESTS_REQUIRED=1
for mode in default inflate3; do
  extra=""; [ $mode = inflate3 ] && extra="PROTAL_INFLATE_THREADS=3"
  taskset -c 0-3 nice -n 5 env $extra PROTAL=$W/build/protal python3 -m unittest tests/e2e/test_protal_e2e.py > $W/e2e_$mode.log 2>&1
  echo "e2e $mode exit $?"; tail -1 $W/e2e_$mode.log; grep -A12 "^FAIL:\|^ERROR:" $W/e2e_$mode.log | cut -c1-300 | head -30
done
cd $W
taskset -c 0-3 nice -n 5 cmake --build build-static --target protal_static -- -j4 > static_build.log 2>&1
echo "static build exit $?"; grep -E "error|undefined" static_build.log | head
S=$(ls build-static/protal_*_static | head -1)
cd $W/src
taskset -c 0-3 nice -n 5 env PROTAL_INFLATE_THREADS=3 PROTAL=$W/$S python3 -m unittest tests/e2e/test_protal_e2e.py > $W/e2e_static.log 2>&1
echo "e2e static inflate3 exit $?"; tail -1 $W/e2e_static.log; grep -A12 "^FAIL:\|^ERROR:" $W/e2e_static.log | cut -c1-300 | head -30
export PROTAL=$W/build/protal PROTAL_TRAIN_PYTHON=$HOME/micromamba/envs/protal-db-build/bin/python
taskset -c 0-3 nice -n 5 python3 -m unittest scripts/mini_db/test_*.py > $W/scripts_mini.log 2>&1
echo "scripts mini_db exit $?"; tail -1 $W/scripts_mini.log; grep -A12 "^FAIL:\|^ERROR:" $W/scripts_mini.log | cut -c1-300 | head -30
