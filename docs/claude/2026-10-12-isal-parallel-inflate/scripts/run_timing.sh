#!/bin/bash
# usage: run_timing.sh <outfile> <rounds>: alternated rounds over all configurations, with a calibration of
# N concurrent independent sequential ISA-L processes (N = 2..4) in each round
D=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/0ff4a1d7-3599-4615-a178-e4ca253109af/scratchpad/isal_parallel
X="taskset -c 0-3 nice -n 5 $D/bin/isalpar"
G=~/gzpar/data/R1.member.gz
B=$D/logs/R1.bounds.tsv
OUT=$1; R=${2:-3}
cat $G > /dev/null
uptime >> $OUT.uptime
for r in $(seq 1 $R); do
  $X seq $G | sed "s/^/round\t$r\t/" >> $OUT
  for N in 2 3 4; do
    for i in $(seq 1 $N); do $X seq $G > /tmp/isalpar_calib_$i.txt & done
    wait
    cat /tmp/isalpar_calib_*.txt | awk -F'\t' -v r=$r -v N=$N '$2=="gzip_crc"{s+=$12} END{printf "round\t%s\tcalib\tconcurrent_seq\tN\t%d\tGBps\t%.3f\n", r, N, s}' >> $OUT
    rm -f /tmp/isalpar_calib_*.txt
  done
  for t in 1 2 3 4; do
    for cfg in "lead 2 4" "spec 2 4" "lead 2 1" "spec 2 1" "lead 2 16" "spec 2 16" "lead 3 4" "spec 3 4"; do
      set -- $cfg
      $X par $G $B $3 $t $1 $2 0 - | sed "s/^/round\t$r\t/" >> $OUT
    done
  done
  uptime >> $OUT.uptime
done
