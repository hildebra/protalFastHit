#!/usr/bin/env bash
# Hp (Hb plus the next line prefetched) against Hb; the key map layouts on a quiet machine, without the values array: D2 direct, Cf/C compact, H hashed (branchy scan),
# Hb hashed at the same load with the branchless scan (HashedKeys.h's), B buckets. Synthetic r226 density (short reads
# a read ahead, long reads streamed) and the reduced r226 database's real cores; one thread, core 2; 5 rounds.
set -uo pipefail
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/95e4f938-fd80-46a7-b960-30fecd6a2871/scratchpad/hkeys
W=~/hkeys
rsync -a --checksum --no-times $SP/keybench.cpp $W/src/src/misc/keybench.cpp
taskset -c 0-3 nice -n 5 ninja -C $W/build -j4 keybench_v3 > $W/ninja_bench2.log 2>&1 || { tail -30 $W/ninja_bench2.log; exit 1; }
B=$W/build/src/keybench_v3
OUT=$W/runs3
mkdir -p $OUT
R226=~/Mocks/reduced_db_r226/database.protal
PE=~/Mocks/reads/d6331_pe_labs/d6331_hku_SRR23926922/SRR23926922_1.fq.gz
HIFI=~/Mocks/reads/d6331_hifi/d6331_hifi_SRR13128012/SRR13128012.fq.gz
export KEYBENCH_H2_BRANCHLESS=1 KEYBENCH_H_PAIR=1
unset KEYBENCH_VALUES_GB
: > $OUT/log.txt
run() {
  name=$1; shift
  echo "== $name: $*" >> $OUT/log.txt
  taskset -c 2 $B "$@" > $OUT/$name.tsv 2>> $OUT/log.txt
  echo "exit $?" >> $OUT/log.txt
  cat $OUT/$name.tsv >> $OUT/log.txt
}
uptime >> $OUT/log.txt
run synth_pe synth 0.092 44 8000000 5 27 1 2000
run synth_long synth 0.092 44 8000000 5 3000 0
run r226red_pe keys $R226 $PE 300000 1 5 0 2000
run r226red_hifi keys $R226 $HIFI 3000 0 5 0
uptime >> $OUT/log.txt
echo ALL_DONE >> $OUT/log.txt
