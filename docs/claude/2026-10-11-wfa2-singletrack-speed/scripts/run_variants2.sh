#!/usr/bin/env bash
# run_variants2.sh: the v3 variant without WFA2's AVX2 extend kernels; vectorization of the heuristic's loops per
# variant; callgrind instructions of one Singletrack round per sample and variant (deterministic, unlike the timings);
# identity of every variant's per-call results with the base build.
set -uo pipefail
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/6a0fc5a0-51ae-4ff5-a1a7-da994d3168e8/scratchpad/wfasimd
W=$HOME/wfasimd; C=$W/calls; O=$W/varout
bash $S/variants.sh bt_v3x bounds,term,noavxext x86-64-v3
for n in bounds bounds_term_v2 bt_v3x; do
  case $n in bounds) a=x86-64 ;; bounds_term_v2) a=x86-64-v2 ;; bt_v3x) a=x86-64-v3 ;; esac
  (cd $W/libs/$n/wfa2-lib && gcc -D_FILE_OFFSET_BITS=64 -I. -Iwavefront -Iutils -fopenmp -O3 -DNDEBUG -march=$a -fopt-info-vec-optimized \
     -c wavefront/wavefront_heuristic.c -o /tmp/vec.o 2>&1 | grep -oE 'heuristic.c:[0-9]+:[0-9]+: optimized: loop vectorized using [0-9]+ byte' | sort -u | tr '\n' ' ')
  echo "<- $n heuristic loops vectorized"
done
for s in pe pb ont; do
  for n in base bounds bounds_term bounds_term_v2 bt_v3x; do
    o=$O/cg.$s.$n
    taskset -c 0-3 valgrind --tool=callgrind --callgrind-out-file=$o.out $W/libs/$n/wfareplay $C/$s.tsv 1 $O/cgres.$s.$n singletrack > /dev/null 2>&1
    tot=$(callgrind_annotate --inclusive=yes --threshold=100 $o.out 2>/dev/null | sed -nE 's/^ *([0-9,]+) .*:wavefront_align \[.*/\1/p' | head -1 | tr -d ,)
    echo "$s $n $tot"
  done
done > $O/cg.tsv
awk '{if ($2=="base") b[$1]=$3; printf "%-4s %-16s %12d  %+6.1f%%\n", $1, $2, $3, 100*($3-b[$1])/b[$1]}' $O/cg.tsv
for s in pe pb ont; do
  cmp -s $O/$s.base.singletrack $O/cgres.$s.bt_v3x.singletrack && echo "$s bt_v3x: results identical to base" || echo "$s bt_v3x: results DIFFER"
done
echo "VAR2DONE $(date +%T)"
