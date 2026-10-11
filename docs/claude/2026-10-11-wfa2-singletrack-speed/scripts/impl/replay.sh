#!/usr/bin/env bash
# replay.sh [ROUNDS] [INVOCATIONS]: wfareplay built against ~/wfaspeed/base (0a1beaa: WFA2 as before, the wrapper with
# the penalties as given) and ~/wfaspeed/new (the wfa-speed changes); every dumped call replayed by both, classic and
# Singletrack, and the per-call results (outcome, score, operations' hash) compared; then thread CPU seconds of
# INVOCATIONS alternated invocations of ROUNDS rounds each (median of the invocations' medians); callgrind instructions of
# one round per sample, traceback and build. On cores 0-3.
set -uo pipefail
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/6a0fc5a0-51ae-4ff5-a1a7-da994d3168e8/scratchpad/wfaspeed
W=$HOME/wfaspeed; C=$W/calls; R=$W/replay
ROUNDS=${1:-9}; INV=${2:-5}
mkdir -p $R
for s in pe pb ont; do
  awk -F'\t' -v s=$s '{n[$1]++; q+=length($8); r+=length($9); if ($7 < 2147483647) m+=$7} END {for (k in n) printf "%s: %s %d calls; ", s, k, n[k]; printf "mean read part %.0f, reference %.0f, budget %.0f\n", q/NR, r/NR, m/NR}' $C/$s.tsv
done
for b in base new; do
  H=$W/$b; LIB=$(find $H/build -name libwfa_lib.a | head -1)
  g++ -O3 -g -march=x86-64 -mtune=generic -std=c++20 -I$H/src -I$H/src/IO -I$H/src/Hash -I$H/src/Core -I$H/src/SequenceUtils -I$H/src/Alignment \
    -I$H/src/Utilities -I$H/src/Taxonomy -I$H/lib -I$H/lib/robin -I$H/lib/tsl -I$H/lib/wfa2-lib -I$H/lib/gzstream -I$H/build/zlib-ng \
    -I$H/lib/zlib-ng -I$H/build/generated $S/wfareplay.cpp $LIB -o $R/wfareplay_$b 2>&1 | grep -E ' error' | head
  [ -x $R/wfareplay_$b ] && echo "built wfareplay_$b ($LIB)" || { echo "FAIL wfareplay_$b"; exit 1; }
done
echo "== identity ($(date +%T))"
for s in pe pb ont; do
  for b in base new; do
    taskset -c 0-3 $R/wfareplay_$b $C/$s.tsv 1 $R/$s.$b classic singletrack > $R/$s.$b.id.log 2>&1
  done
  for t in classic singletrack; do
    n=$(wc -l < $R/$s.base.$t)
    cmp -s $R/$s.base.$t $R/$s.new.$t && echo "   $s $t: $n calls, results identical" || \
      echo "   $s $t: $n calls, DIFFER in $(paste -d' ' $R/$s.base.$t $R/$s.new.$t | awk '$1!=$4 || $2!=$5 || $3!=$6' | wc -l) (outcome or score in $(paste -d' ' $R/$s.base.$t $R/$s.new.$t | awk '$1!=$4 || $2!=$5' | wc -l))"
  done
  grep -h "fallbacks" $R/$s.base.id.log $R/$s.new.id.log | tr '\n' ' ' | sed 's/^/   fallbacks base, new: /'; echo
done
echo "== timing: $INV alternated invocations of $ROUNDS rounds ($(date +%T))"
for s in pe pb ont; do
  : > $R/t.$s.txt
  for i in $(seq $INV); do
    for b in base new; do
      taskset -c 0-3 $R/wfareplay_$b $C/$s.tsv $ROUNDS /tmp/wfaspeed_t classic singletrack 2>/dev/null | awk -v b=$b '/median/ {print b, $1, $3}' >> $R/t.$s.txt
    done
  done
  for t in classic singletrack; do
    bm=$(awk -v t=$t '$1=="base" && $2==t {print $3}' $R/t.$s.txt | sort -g | awk '{a[NR]=$1} END {print a[int((NR+1)/2)]}')
    nm=$(awk -v t=$t '$1=="new" && $2==t {print $3}' $R/t.$s.txt | sort -g | awk '{a[NR]=$1} END {print a[int((NR+1)/2)]}')
    awk -v s=$s -v t=$t -v b=$bm -v n=$nm 'BEGIN {printf "   %-4s %-12s base %.4f s, new %.4f s, %+.1f%%\n", s, t, b, n, 100*(n-b)/b}'
  done
done
echo "== callgrind, one round ($(date +%T))"
for s in pe pb ont; do
  for t in classic singletrack; do
    for b in base new; do
      o=$R/cg.$s.$t.$b.out
      taskset -c 0-3 valgrind --tool=callgrind --dump-instr=yes --callgrind-out-file=$o $R/wfareplay_$b $C/$s.tsv 1 /tmp/wfaspeed_cg $t > /dev/null 2>&1
      tot=$(callgrind_annotate --inclusive=yes --threshold=100 $o 2>/dev/null | grep -E ':wavefront_align ' | head -1 | sed -E 's/^ *([0-9,]+) .*/\1/' | tr -d ,)
      echo "$s $t $b $tot"
    done
  done
done > $R/cg.tsv
awk '{k=$1" "$2; if ($3=="base") b[k]=$4; else printf "   %-4s %-12s base %12d, new %12d instructions in wavefront_align, %+.1f%%\n", $1, $2, b[k], $4, 100*($4-b[k])/b[k]}' $R/cg.tsv
echo "REPLAYDONE $(date +%T)"
