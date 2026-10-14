#!/bin/bash
# Old (PROTAL_OLD_LOAD=1) against new load paths, alternating, 1 and 2 threads (cores 4 / 4,5), 3 rounds.
K=~/sparseprof/build/keybench
OUT=~/sparseprof/runs/load_ab2
mkdir -p $OUT
GRID=~/strainfn/db_grid/database.protal
R226=~/Mocks/reduced_db_r226/database.protal
for round in 1 2 3; do
  for threads in 1 2; do
    cpus=4; [ $threads = 2 ] && cpus=4,5
    for cfg in "grid 0 0 0" "grid 0 1 1" "grid 2 0 1" "grid 2 1 1" "grid 4 0 1" "grid 4 1 1" "r226 0 1 0"; do
      set -- $cfg
      db=$GRID; [ $1 = r226 ] && db=$R226
      for old in 1 0; do
        if [ $old = 1 ]; then e="PROTAL_OLD_LOAD=1"; else e="KEYBENCH_NEW=1"; fi
        echo "== $1 level $2 compact $3 bits $4 threads $threads old $old round $round" >> $OUT/phases.txt
        env $e PROTAL_LOAD_PHASES=1 taskset -c $cpus $K load $db $2 $3 $threads $4 2>> $OUT/phases.txt | sed "s/^/$1\told=$old\t/" >> $OUT/loads.tsv
      done
    done
  done
done
echo done >> $OUT/phases.txt
