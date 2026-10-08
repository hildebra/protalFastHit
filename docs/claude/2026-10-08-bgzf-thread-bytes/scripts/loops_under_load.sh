#!/bin/bash
# The tests in loops on cores 0-3 (nice 5) while protal_tests is rebuilt from scratch (ccache off) on the same cores.
# Each loop line in summary.txt: binary, test, runs, failures. Failing logs are kept in runs/.
set -u
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
W=~/bgzfdet
cd $W
mkdir -p runs
rm -f runs/* summary.txt stop_load load_rounds.log
# the fixed binary from the working files as they are now
cp $REPO/src/IO/Bgzf.h src/src/IO/Bgzf.h
cp $REPO/tests/test_RunStatusAndBgzf.cpp src/tests/test_RunStatusAndBgzf.cpp
taskset -c 0-3 nice -n 5 ninja -C build -j4 protal_tests > ninja_fixed2.log 2>&1 || { tail -20 ninja_fixed2.log; exit 1; }
cp build/tests/protal_tests bin/tests_fixed
echo "tests_fixed from the working files: $(md5sum < src/src/IO/Bgzf.h | cut -c1-8) $(md5sum < src/tests/test_RunStatusAndBgzf.cpp | cut -c1-8)" | tee summary.txt
# load: clean builds of protal_tests, one after the other, until stop_load appears
(
    n=0
    while [ ! -e stop_load ]; do
        n=$((n + 1))
        rm -rf loadbuild
        CCACHE_DISABLE=1 taskset -c 0-3 nice -n 5 cmake -S src -B loadbuild -G Ninja -DCMAKE_BUILD_TYPE=Release \
            -DPROTAL_BUILD_TESTS=ON > load_cmake.log 2>&1
        echo "load build $n start $(date +%T) load $(cut -d' ' -f1-3 /proc/loadavg)" >> load_rounds.log
        CCACHE_DISABLE=1 taskset -c 0-3 nice -n 5 ninja -C loadbuild -j4 protal_tests > load_ninja.log 2>&1
        echo "load build $n end $(date +%T) rc $?" >> load_rounds.log
    done
) &
LOAD=$!
sleep 60  # let the build get going

loop() {  # binary filter runs [extra gtest args]
    local bin=$1 filter=$2 runs=$3
    shift 3
    local fail=0
    for i in $(seq 1 $runs); do
        if ! taskset -c 0-3 nice -n 5 bin/$bin --gtest_filter="$filter" "$@" > runs/current.log 2>&1; then
            fail=$((fail + 1))
            mv runs/current.log "runs/${bin}_$(echo $filter | tr '.*' '__')_$i.log"
        fi
    done
    echo "$bin $filter x$runs $*: $fail failed ($(date +%T), load $(cut -d' ' -f1-3 /proc/loadavg))" | tee -a summary.txt
}
echo "start $(date +%T)" | tee -a summary.txt

loop tests_newtest_oldcode 'BgzfWriter.BlocksDoNotDependOnTheThread' 20
loop tests_fixed 'BgzfWriter.*' 20
loop tests_head 'LongReadSimulation.SamplesTemplatesAndThreads' 500
loop tests_fixed 'LongReadSimulation.SamplesTemplatesAndThreads' 500
loop tests_fixed 'LongReadSimulation.SamplesTemplatesAndThreads' 1 --gtest_repeat=500
loop tests_fixed '*' 5
touch stop_load
wait $LOAD
echo "done $(date +%T)" | tee -a summary.txt
cat load_rounds.log
