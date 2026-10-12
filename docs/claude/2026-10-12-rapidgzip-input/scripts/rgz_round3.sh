#!/bin/bash
# Rebuild; reader tests; the truncation sweep (bench_gz); all unit tests; e2e with every gzip on several threads.
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/0ff4a1d7-3599-4615-a178-e4ca253109af/scratchpad
W=~/rgzbuild
bash $SP/rgz_build.sh 2>&1 | grep -E "exit|error" | head -20
cd $W/build
taskset -c 0-3 nice -n 5 ./tests/protal_tests --gtest_filter='ThreadedGzStream.*:SeqReader*' > $W/unit_reader.log 2>&1
echo "reader tests exit $?"
grep -E "^\[  (FAILED|PASSED) |Failure|trace:|threads$|Which is|Actual|bytes read" $W/unit_reader.log | cut -c1-300 | head -60
bash $SP/rgz_trunc.sh 2>&1 | grep -c "three \[parallel 1 failed 0" | sed 's/^/truncated files read as whole on 3 threads: /'
bash $SP/rgz_trunc.sh 2>&1 | grep "three \[parallel 0\|failed 0" | head
taskset -c 0-3 nice -n 5 ctest --test-dir $W/build -j4 > $W/ctest.log 2>&1; echo "ctest exit $?"; grep -E "tests passed|tests failed|\*\*\*Failed" $W/ctest.log | head
cd $W/src
taskset -c 0-3 nice -n 5 env PROTAL_INFLATE_THREADS=3 PROTAL_TEST_DB=data/mini_db/protal_db PROTAL=$W/build/protal \
    SIMULATE=$W/build/simulate_metagenomes PROTAL_TESTS_REQUIRED=1 python3 -m unittest tests/e2e/test_protal_e2e.py > $W/e2e_inflate3.log 2>&1
echo "e2e inflate3 exit $?"; tail -3 $W/e2e_inflate3.log; grep -A12 "^FAIL:\|^ERROR:" $W/e2e_inflate3.log | cut -c1-300 | head -30
