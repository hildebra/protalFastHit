#!/bin/bash
# The loaded index's checksum (every core's slots and cells) of the old and the new load paths, 2 threads.
K=~/sparseprof/build/keybench
GRID=~/strainfn/db_grid/database.protal
R226=~/Mocks/reduced_db_r226/database.protal
for cfg in "0 0 1" "0 1 1" "2 0 1" "2 1 1" "4 0 1" "4 1 1"; do
  set -- $cfg
  for old in 1 0; do
    if [ $old = 1 ]; then e="PROTAL_OLD_LOAD=1"; else e="KEYBENCH_NEW=1"; fi
    env $e KEYBENCH_CHECKSUM=1 taskset -c 4,5 $K load $GRID $1 $2 2 $3 2>/dev/null | grep checksum | sed "s/^/grid old=$old\t/"
  done
done
for c in 0 1; do
  for old in 1 0; do
    if [ $old = 1 ]; then e="PROTAL_OLD_LOAD=1"; else e="KEYBENCH_NEW=1"; fi
    env $e KEYBENCH_CHECKSUM=1 taskset -c 4,5 $K load $R226 0 $c 2 0 2>/dev/null | grep checksum | sed "s/^/r226red old=$old\t/"
  done
done
