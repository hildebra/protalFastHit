#!/usr/bin/env bash
# The EM on indices: HEAD's archive + em.patch built in ~/mt-work/perf4/em (protal and protal_tests); the
# SampleContext tests; --profile_only outputs compared byte for byte with the reference build's; timings.
set -uo pipefail
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
SCR=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/85c3e3d5-6df2-4cf2-a5bc-63aafb61f818/scratchpad
W=$HOME/mt-work/perf4; E=$W/em
rm -rf $E; mkdir -p $E
git -C $REPO archive HEAD | tar -x -C $E
(cd $E && patch -p1 -s < $SCR/em.patch) || { echo "FAIL patch"; exit 1; }
cd $E && nice cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON > $W/em.configure.log 2>&1 && \
  nice cmake --build build --target protal protal_tests -j 5 > $W/em.build.log 2>&1 && echo "OK em build" || { echo "FAIL em build"; grep -E 'error' -A3 $W/em.build.log | head -40; exit 1; }
$E/build/tests/protal_tests --gtest_filter='SampleContext.*' 2>&1 | tail -5
B=$E/build/protal; R=$W/ref/build/protal; DB=$HOME/bench071/V073/protal_db
cmpp() { # name threads sam
  local n=$1 t=$2 sam=$3
  for b in ref em; do
    local bin=$R; [ $b = em ] && bin=$B
    rm -rf $W/c.$b.$n
    /usr/bin/time -f "%e s wall, %U s user" $bin --db $DB --profile_only $sam --prefix s -o $W/c.$b.$n -t $t --no_qcmsa --verbose > $W/c.$b.$n.log 2> $W/c.$b.$n.time
    echo "   $b $n t=$t: $(tail -1 $W/c.$b.$n.time); $(grep -E '^Profiling took' $W/c.$b.$n.log)"
  done
  if diff -r $W/c.ref.$n $W/c.em.$n > $W/c.$n.diff; then echo "== $n: outputs identical ($(find $W/c.em.$n -type f | wc -l) files)"; else echo "== $n: OUTPUTS DIFFER"; head -20 $W/c.$n.diff; fi
}
cmpp pe500k_t1 1 $W/o.pe500k_t1/s.sam.zst
cmpp pe500k_t6 6 $W/o.pe500k_t1/s.sam.zst
cmpp pe5M_t6 6 $W/o.pe5M_t6/s.sam.zst
cmpp ont90M_t6 6 $W/o.ont90M_t6/s.sam.zst
cmpp pb90M_t6 6 $W/o.pb90M_t6/s.sam.zst
