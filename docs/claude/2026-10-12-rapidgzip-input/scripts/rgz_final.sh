#!/bin/bash
# The final code: inflate speed through ThreadedGzIstream (1-4 threads, 3 alternated rounds), e2e as usual, the static
# binary and its e2e with every gzip input on several threads.
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/0ff4a1d7-3599-4615-a178-e4ca253109af/scratchpad
W=~/rgzbuild
cd $W
F=~/gzpar/data/R1.member.gz
cat $F > /dev/null
for round in 1 2 3; do
  for t in 1 2 3 4; do
    taskset -c 0-3 ./bench_gz $F $t
  done
done
cd $W/src
taskset -c 0-3 nice -n 5 env PROTAL_TEST_DB=data/mini_db/protal_db PROTAL=$W/build/protal SIMULATE=$W/build/simulate_metagenomes \
    PROTAL_TESTS_REQUIRED=1 python3 -m unittest tests/e2e/test_protal_e2e.py > $W/e2e_default.log 2>&1
echo "e2e default exit $?"; tail -1 $W/e2e_default.log
cd $W
taskset -c 0-3 nice -n 5 cmake --build build-static --target protal_static -- -j4 > static_build.log 2>&1
echo "static build exit $?"; grep -E "error|undefined" static_build.log | head
S=$(ls build-static/protal_*_static)
file $S | cut -c1-120
cd $W/src
taskset -c 0-3 nice -n 5 env PROTAL_INFLATE_THREADS=3 PROTAL_TEST_DB=data/mini_db/protal_db PROTAL=$W/$S \
    SIMULATE=$W/build/simulate_metagenomes PROTAL_TESTS_REQUIRED=1 python3 -m unittest tests/e2e/test_protal_e2e.py > $W/e2e_static.log 2>&1
echo "e2e static inflate3 exit $?"; tail -1 $W/e2e_static.log
