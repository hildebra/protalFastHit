#!/bin/bash
# The final working files (Bgzf.h with its comment reflowed) built as bin/tests_final, checked equal to the working
# files, then run in loops on cores 0-3 (nice 5) while protal_tests is rebuilt from scratch (ccache off) there.
set -u
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
W=~/bgzfdet
cd $W
mkdir -p runs_final
rm -f runs_final/* summary_final.txt stop_load load_rounds_final.log
cp $REPO/src/IO/Bgzf.h src/src/IO/Bgzf.h
cp $REPO/tests/test_RunStatusAndBgzf.cpp src/tests/test_RunStatusAndBgzf.cpp
taskset -c 0-3 nice -n 5 ninja -C build -j4 protal_tests > ninja_final.log 2>&1 || { tail -20 ninja_final.log; exit 1; }
cp build/tests/protal_tests bin/tests_final
# the tree against the commit and the working files: only the two changed files differ from 401c4f5
for f in src/IO/Bgzf.h tests/test_RunStatusAndBgzf.cpp; do
    cmp -s $REPO/$f src/$f && echo "$f: tree = working file" | tee -a summary_final.txt
done
(
    n=0
    while [ ! -e stop_load ]; do
        n=$((n + 1))
        rm -rf loadbuild
        CCACHE_DISABLE=1 taskset -c 0-3 nice -n 5 cmake -S src -B loadbuild -G Ninja -DCMAKE_BUILD_TYPE=Release \
            -DPROTAL_BUILD_TESTS=ON > load_cmake.log 2>&1
        echo "load build $n start $(date +%T) load $(cut -d' ' -f1-3 /proc/loadavg)" >> load_rounds_final.log
        CCACHE_DISABLE=1 taskset -c 0-3 nice -n 5 ninja -C loadbuild -j4 protal_tests > load_ninja.log 2>&1
        echo "load build $n end $(date +%T) rc $?" >> load_rounds_final.log
    done
) &
LOAD=$!
sleep 60

loop() {  # binary filter runs [extra gtest args]
    local bin=$1 filter=$2 runs=$3
    shift 3
    local fail=0
    for i in $(seq 1 $runs); do
        if ! taskset -c 0-3 nice -n 5 bin/$bin --gtest_filter="$filter" "$@" > runs_final/current.log 2>&1; then
            fail=$((fail + 1))
            mv runs_final/current.log "runs_final/${bin}_$(echo $filter | tr '.*' '__')_$i.log"
        fi
    done
    echo "$bin $filter x$runs $*: $fail failed ($(date +%T), load $(cut -d' ' -f1-3 /proc/loadavg))" | tee -a summary_final.txt
}
echo "start $(date +%T)" | tee -a summary_final.txt
loop tests_final 'BgzfWriter.*' 20
loop tests_final 'LongReadSimulation.SamplesTemplatesAndThreads' 2000
loop tests_final '*' 2
grep -E "tests from|PASSED|FAILED|BlocksDoNotDependOnTheThread \(" runs_final/current.log | tail -4 | tee -a summary_final.txt
touch stop_load
wait $LOAD
echo "done $(date +%T)" | tee -a summary_final.txt
cat load_rounds_final.log
