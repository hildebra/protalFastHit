#!/usr/bin/env bash
# Input-path audit: how long should a waiter for the reader lock spin before it sleeps? protal-inexp (a212559,
# batch 32: critical(reader)) with GOMP_SPINCOUNT 300000 (libgomp's default), 30000, 3000, 300 and 0 (as
# OMP_WAIT_POLICY=passive), alternated, niced: mini database at 4 and 6 threads, full database at 6 threads,
# gzip and uncompressed. Rows into ~/mt-work/input/spin.tsv. Runs after the mutex benchmark.
set -uo pipefail
SP=$(cd "$(dirname "$0")" && pwd)  # this folder
W=$HOME/mt-work/input
B=$HOME/bench071
declare -A R1=( [gzip]=$HOME/mt-audit/plain/R1.std.gz [plain]=$HOME/mt-audit/plain/R1.fq )
declare -A R2=( [gzip]=$HOME/mt-audit/plain/R2.std.gz [plain]=$HOME/mt-audit/plain/R2.fq )
declare -A DB=( [mini]=$HOME/mt-work/rbfull/mini_db/protal_db [full]=$B/V071/protal_db )
[ -f $W/spin.tsv ] || printf "db\tinput\tthreads\tspincount\trep\taligning_s\tpairs_per_s\twall_s\tuser_s\tsys_s\tload\treader_s_per_thread\n" > $W/spin.tsv
seconds() { awk '{ s = 0; for (i = 1; i <= NF; i++) { v = $i; if (v ~ /ms$/) { sub(/ms/, "", v); s += v / 1000 } else if (v ~ /m$/) { sub(/m/, "", v); s += v * 60 } else if (v ~ /s$/) { sub(/s/, "", v); s += v } } printf "%.3f", s }'; }
quiet() { local i; for i in $(seq 6); do awk -v l=$(cut -d" " -f1 /proc/loadavg) "BEGIN { exit !(l < 6) }" && return; sleep 10; done; }
run() {  # run DB INPUT THREADS SPINCOUNT REP
  quiet
  local db=$1 in=$2 t=$3 spin=$4 rep=$5 d=$W/run
  rm -rf $d; mkdir -p $d
  local load=$(cut -d' ' -f1 /proc/loadavg)
  PROTAL_EXPERIMENT_BATCH=32 GOMP_SPINCOUNT=$spin nice -n 5 /usr/bin/time -f "%e %U %S" -o $d.time $HOME/mt-work/bin/protal-inexp --db ${DB[$db]} \
      -1 ${R1[$in]} -2 ${R2[$in]} --prefix s -o $d -t $t --no_profile > $d.log 2>&1
  local a=$(grep 'Aligning reads took' $d.log | sed 's/.*took //' | seconds)
  local reader=$(awk -F'\t' '$1 == "Sequence reader" { print $4 }' $d/misc/s_runtime.tsv)
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $db $in $t $spin $rep $a $(awk -v a=$a 'BEGIN { printf "%.0f", 4995812 / a }') \
      $(tail -1 $d.time) $load $reader | tee -a $W/spin.tsv
}
for rep in 1 2 3; do
  cat ${R1[@]} ${R2[@]} > /dev/null
  for t in 6 4; do for in in gzip plain; do for s in 300000 30000 3000 300 0; do run mini $in $t $s $rep; done; done; done
  for in in gzip plain; do for s in 300000 30000 3000 300 0; do run full $in 6 $s $rep; done; done
done
echo SPIN DONE
