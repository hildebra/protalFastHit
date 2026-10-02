#!/usr/bin/env bash
# Profiling the 5M-pair SAM (--profile_only): f359d49 (fbase2) against f359d49 + both commits (final2, = a212559's
# src), at 1, 6 and 8 threads (8: two decompression threads, on 6 vCPUs), alternated, niced, CPU seconds too.
set -uo pipefail
W=$HOME/mt-work/speedf; mkdir -p $W
SAM=$HOME/mt-audit/runs/b6/s1.sam.zst
DB=$HOME/bench071/V071/protal_db
declare -A BIN=( [base]=$HOME/mt-work/fbase2/src/build/protal [new]=$HOME/mt-work/final2/src/build/protal )
[ -f $W/runs.tsv ] || printf "binary\tthreads\trep\tprofiling_s\twall_s\tuser_s\tsys_s\tload\n" > $W/runs.tsv
seconds() { awk '{ s = 0; for (i = 1; i <= NF; i++) { v = $i; if (v ~ /ms$/) { sub(/ms/, "", v); s += v / 1000 } else if (v ~ /m$/) { sub(/m/, "", v); s += v * 60 } else if (v ~ /s$/) { sub(/s/, "", v); s += v } } printf "%.3f", s }'; }
quiet() { local i; for i in $(seq 6); do awk -v l=$(cut -d" " -f1 /proc/loadavg) "BEGIN { exit !(l < 6) }" && return; sleep 10; done; }
run() {  # run BIN THREADS REP
  quiet
  local bin=$1 t=$2 rep=$3 o=$W/out
  local load=$(cut -d' ' -f1 /proc/loadavg)
  rm -rf $o
  nice -n 5 /usr/bin/time -f "%e %U %S" -o $W/time ${BIN[$bin]} --db $DB --profile_only $SAM --prefix pe5M -o $o -t $t --no_qcmsa > $W/log 2>&1
  local p=$(grep 'Profiling took' $W/log | sed 's/.*took //' | seconds)
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" $bin $t $rep $p "$(tail -1 $W/time | tr ' ' '\t')" $load | tee -a $W/runs.tsv
}
for rep in 1 2 3; do
  for t in 1 6 8; do run base $t $rep; run new $t $rep; done
done
echo SPEED DONE
