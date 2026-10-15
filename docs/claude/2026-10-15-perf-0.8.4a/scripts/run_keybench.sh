#!/usr/bin/env bash
# The key map layouts (D2 direct, C/Cf compact, H hashed lines at loads 0.75 and 0.5, B 128-core buckets): synthetic
# r226 density (short reads paired with a 2 us gap, long reads streamed) and the reduced r226 database's real cores
# (D6331 pe, HiFi), each without and with an 8 GB values array read as FindSeeds reads blocks. One thread, core 2.
set -uo pipefail
W=~/hkeys
B=$W/build/keybench_v3
[ -x $B ] || B=$(find $W/build -maxdepth 2 -name keybench_v3 -type f | head -1)
OUT=$W/runs
mkdir -p $OUT
R226=~/Mocks/reduced_db_r226/database.protal
PE=~/Mocks/reads/d6331_pe_labs/d6331_hku_SRR23926922/SRR23926922_1.fq.gz
HIFI=~/Mocks/reads/d6331_hifi/d6331_hifi_SRR13128012/SRR13128012.fq.gz
run() {
  name=$1; shift
  for values in 0 8; do
    echo "== $name values $values GB: $*" >> $OUT/log.txt
    if [ $values = 0 ]; then unset KEYBENCH_VALUES_GB; else export KEYBENCH_VALUES_GB=$values; fi
    taskset -c 2 nice -n 5 $B "$@" > $OUT/${name}_v$values.tsv 2>> $OUT/log.txt
    echo "exit $?" >> $OUT/log.txt
    cat $OUT/${name}_v$values.tsv >> $OUT/log.txt
  done
  unset KEYBENCH_VALUES_GB
}
: > $OUT/log.txt
run synth_pe synth 0.092 44 8000000 3 27 1 2000
run synth_long synth 0.092 44 8000000 3 3000 0
run r226red_pe keys $R226 $PE 300000 1 3 0 2000
run r226red_hifi keys $R226 $HIFI 3000 0 3 0
echo ALL_DONE >> $OUT/log.txt
