#!/bin/bash
# Callgrind of whole programs (no toggles: --toggle-collect on *LoadColumns* switched collection off inside nested
# matches), one load each, 1 thread, core 5; then the inclusive cost of the load's parts.
K=~/sparseprof/build/keybench
OUT=~/sparseprof/runs/cg2
mkdir -p $OUT
GRID=~/strainfn/db_grid/database.protal
for cfg in "2 1 1 sparse2c" "2 0 1 sparse2d" "0 1 1 wholec" "0 0 0 wholed"; do
  set -- $cfg
  for old in 1 0; do
    [ $4 = wholed ] && [ $old = 1 ] && continue
    if [ $old = 1 ]; then e="PROTAL_OLD_LOAD=1"; else e="KEYBENCH_NEW=1"; fi
    env $e taskset -c 5 nice -n 5 valgrind --tool=callgrind --callgrind-out-file=$OUT/$4_old$old.out $K load $GRID $1 $2 1 $3 > $OUT/$4_old$old.log 2>&1
    callgrind_annotate --inclusive=yes $OUT/$4_old$old.out 2>/dev/null | sed -n '/^Ir/,$p' | head -70 | cut -c1-260 > $OUT/$4_old$old.incl.txt
    callgrind_annotate $OUT/$4_old$old.out 2>/dev/null | sed -n '/^Ir/,$p' | head -40 | cut -c1-260 > $OUT/$4_old$old.self.txt
  done
done
echo done > $OUT/done
