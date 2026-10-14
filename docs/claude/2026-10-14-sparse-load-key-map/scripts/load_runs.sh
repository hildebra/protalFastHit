#!/bin/bash
# Index loads of the strain world's grid database (and the reduced r226 one, whole only), 1 and 2 threads, pinned to
# core 4 (1 thread) or cores 4,5; 3 rounds alternating the configurations. Phases to stderr.
K=~/sparseprof/build/keybench
OUT=~/sparseprof/runs/load_${1:-base}
mkdir -p $OUT
GRID=~/strainfn/db_grid/database.protal
R226=~/Mocks/reduced_db_r226/database.protal
for round in 1 2 3; do
  for threads in 1 2; do
    cpus=4; [ $threads = 2 ] && cpus=4,5
    for cfg in "0 0 0" "0 1 0" "0 0 1" "0 1 1" "2 0 1" "2 1 1" "4 0 1" "4 1 1"; do
      set -- $cfg
      echo "== grid level $1 compact $2 gridbits $3 threads $threads round $round" >> $OUT/phases.txt
      PROTAL_LOAD_PHASES=1 taskset -c $cpus $K load $GRID $1 $2 $threads $3 >> $OUT/loads.tsv 2>> $OUT/phases.txt
    done
    for cfg in "0 0" "0 1"; do
      set -- $cfg
      echo "== r226reduced level $1 compact $2 threads $threads round $round" >> $OUT/phases.txt
      PROTAL_LOAD_PHASES=1 taskset -c $cpus $K load $R226 $1 $2 $threads 0 | sed 's/^/r226reduced\t/' >> $OUT/loads.tsv 2>> $OUT/phases.txt
    done
  done
done
echo done >> $OUT/phases.txt
