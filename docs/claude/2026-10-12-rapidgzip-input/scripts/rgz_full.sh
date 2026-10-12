#!/bin/bash
# Full check of the working tree's rapidgzip change (cores 0-3): build, unit tests, mini DB, e2e as usual and with
# every gzip input on several threads (PROTAL_INFLATE_THREADS=3: rapidgzip / BGZF workers), then the static binary
# (system libzstd.a) and the e2e suite with it.
set -uo pipefail
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/0ff4a1d7-3599-4615-a178-e4ca253109af/scratchpad
W=~/rgzbuild
bash $SP/rgz_build.sh 2>&1 | grep -E "exit|error" | head -20
cd $W/build
taskset -c 0-3 nice -n 5 ctest -j4 > $W/ctest.log 2>&1; echo "ctest exit $?"; grep -E "tests passed|tests failed|\*\*\*Failed" $W/ctest.log | head
cd $W/src
if [ ! -d data/mini_db/protal_db ]; then
  taskset -c 0-3 nice -n 5 env PROTAL=$W/build/protal bash scripts/mini_db/build_mini_db.sh data/mini_db > $W/minidb.log 2>&1; echo "mini db exit $?"
fi
e2e() {  # label, protal binary, extra env
  local label=$1 bin=$2; shift 2
  taskset -c 0-3 nice -n 5 env "$@" PROTAL_TEST_DB=data/mini_db/protal_db PROTAL=$bin SIMULATE=$W/build/simulate_metagenomes \
      PROTAL_TESTS_REQUIRED=1 python3 -m unittest tests/e2e/test_protal_e2e.py > $W/e2e_$label.log 2>&1
  echo "e2e $label exit $?"; tail -3 $W/e2e_$label.log
}
e2e default $W/build/protal
e2e inflate3 $W/build/protal PROTAL_INFLATE_THREADS=3
ls -la /usr/lib/x86_64-linux-gnu/libzstd.a
cd $W
if [ ! -f build-static/build.ninja ]; then
  taskset -c 0-3 nice -n 5 cmake -S src -B build-static -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_STATIC_FETCH_DEPS=OFF > static_configure.log 2>&1
  echo "static configure exit $?"
fi
taskset -c 0-3 nice -n 5 cmake --build build-static --target protal_static -- -j4 > static_build.log 2>&1
echo "static build exit $?"; grep -E "error|undefined|warning: Using" static_build.log | head -20
S=$(ls build-static/protal_*_static)
file $S
ldd $S 2>&1 | head -3
$S --version | head -3
cd $W/src
e2e static_inflate3 $W/$S PROTAL_INFLATE_THREADS=3
