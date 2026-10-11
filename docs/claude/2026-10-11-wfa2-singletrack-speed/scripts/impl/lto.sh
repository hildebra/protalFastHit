#!/usr/bin/env bash
# lto.sh [ROUNDS] [INVOCATIONS]: the new WFA2 sources (wfa-speed) compiled with -flto into wfareplay (so that GCC
# inlines across WFA2's files: the per-step helpers of wavefront_compute.c into the drivers of wavefront_compute_affine.c,
# wavefront_termination.c into the extend kernels, ...), against the same sources without LTO (wfareplay_new of
# replay.sh): per-call results, callgrind instructions, alternated timings. On cores 0-3.
set -uo pipefail
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/6a0fc5a0-51ae-4ff5-a1a7-da994d3168e8/scratchpad/wfaspeed
W=$HOME/wfaspeed; C=$W/calls; R=$W/replay; H=$W/new; L=$R/lto
ROUNDS=${1:-9}; INV=${2:-5}
rm -rf $L; mkdir -p $L
WD=$H/lib/wfa2-lib
SRCS=$(sed -n '/^set(WFA2_SOURCES/,/)/p' $H/lib/wfa2-lib.cmake | grep -oE '[A-Za-z_/0-9]+\.(c|cpp)')
for f in $SRCS; do
  o=$L/$(echo $f | tr / _).o
  case $f in
    *.cpp) taskset -c 0-3 g++ -O3 -DNDEBUG -g -std=c++20 -flto -D_FILE_OFFSET_BITS=64 -I$WD -I$WD/wavefront -I$WD/utils -c $WD/$f -o $o ;;
    *.c) taskset -c 0-3 gcc -O3 -DNDEBUG -g -flto -D_FILE_OFFSET_BITS=64 -I$WD -I$WD/wavefront -I$WD/utils -c $WD/$f -o $o ;;
  esac || { echo "FAIL $f"; exit 1; }
done
taskset -c 0-3 g++ -O3 -g -flto -march=x86-64 -mtune=generic -std=c++20 -I$H/src -I$H/src/IO -I$H/src/Hash -I$H/src/Core -I$H/src/SequenceUtils -I$H/src/Alignment \
  -I$H/src/Utilities -I$H/src/Taxonomy -I$H/lib -I$H/lib/robin -I$H/lib/tsl -I$WD -I$H/lib/gzstream -I$H/build/zlib-ng \
  -I$H/lib/zlib-ng -I$H/build/generated $S/wfareplay.cpp $L/*.o -o $R/wfareplay_lto 2>&1 | grep -E 'error|warning' | head
[ -x $R/wfareplay_lto ] && echo "built wfareplay_lto ($(echo $SRCS | wc -w) sources)" || { echo FAIL link; exit 1; }
for s in pe pb ont; do
  taskset -c 0-3 $R/wfareplay_lto $C/$s.tsv 1 $R/$s.lto classic singletrack > /dev/null 2>&1
  for t in classic singletrack; do
    cmp -s $R/$s.base.$t $R/$s.lto.$t && echo "   $s $t: lto results identical to base" || echo "   $s $t: lto results DIFFER from base"
  done
done
for s in pe pb ont; do
  : > $R/tl.$s.txt
  for i in $(seq $INV); do
    for b in new lto; do
      taskset -c 0-3 $R/wfareplay_$b $C/$s.tsv $ROUNDS /tmp/wfaspeed_t classic singletrack 2>/dev/null | awk -v b=$b '/median/ {print b, $1, $3}' >> $R/tl.$s.txt
    done
  done
  for t in classic singletrack; do
    bm=$(awk -v t=$t '$1=="new" && $2==t {print $3}' $R/tl.$s.txt | sort -g | awk '{a[NR]=$1} END {print a[int((NR+1)/2)]}')
    nm=$(awk -v t=$t '$1=="lto" && $2==t {print $3}' $R/tl.$s.txt | sort -g | awk '{a[NR]=$1} END {print a[int((NR+1)/2)]}')
    awk -v s=$s -v t=$t -v b=$bm -v n=$nm 'BEGIN {printf "   %-4s %-12s new %.4f s, lto %.4f s, %+.1f%%\n", s, t, b, n, 100*(n-b)/b}'
  done
done
tot() { callgrind_annotate $1 2>/dev/null | grep -m1 "PROGRAM TOTALS" | sed -E 's/^ *([0-9,]+) .*/\1/' | tr -d ,; }
for s in pe pb ont; do
  for t in classic singletrack; do
    o=$R/cg.$s.$t.lto.out
    taskset -c 0-3 valgrind --tool=callgrind --dump-instr=yes --callgrind-out-file=$o $R/wfareplay_lto $C/$s.tsv 1 /tmp/wfaspeed_cg $t > /dev/null 2>&1
    # with LTO wavefront_align may be inlined into its caller: program totals (the calls' parsing is the same code)
    lto=$(tot $o); newt=$(tot $R/cg.$s.$t.new.out)
    wfa=$(awk -v s=$s -v t=$t '$1==s && $2==t && $3=="new" {print $4}' $R/cg.tsv)
    awk -v s=$s -v t=$t -v b=$newt -v n=$lto -v w=$wfa 'BEGIN {printf "   %-4s %-12s program new %12d, lto %12d: %+.1f%% of new wavefront_align (%d)\n", s, t, b, n, 100*(n-b)/w, w}'
  done
done
echo "LTODONE $(date +%T)"
