#!/bin/bash
# Runs timing rounds only when the machine is quiet (no other runnable processes for 5 s).
# usage: run_quiet.sh <outfile> <rounds> <max_wait_s_per_round>
D=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/0ff4a1d7-3599-4615-a178-e4ca253109af/scratchpad/isal_parallel
X="taskset -c 0-3 nice -n 5 $D/bin/isalpar"
G=~/gzpar/data/R1.member.gz
B=$D/logs/R1.bounds.tsv
OUT=$1; R=${2:-3}; MAXW=${3:-900}
quiet() {
  local ok=0 waited=0
  while [ $ok -lt 5 ]; do
    run=$(awk '{split($4,a,"/"); print a[1]}' /proc/loadavg)
    busy=$(ps -eo pcpu=,comm= | awk '$1 > 20 && $2 != "isalpar" {n++} END{print n+0}')
    if [ "$run" -le 2 ] && [ "$busy" -eq 0 ]; then ok=$((ok+1)); else ok=0; fi
    sleep 1; waited=$((waited+1))
    [ $waited -ge $MAXW ] && { echo "gave up waiting after $waited s" >> $OUT.uptime; return 1; }
  done
  echo "quiet after $waited s: $(cat /proc/loadavg)" >> $OUT.uptime
}
cat $G > /dev/null
for r in $(seq 1 $R); do
  quiet || continue
  $X seq $G | sed "s/^/round\t$r\t/" >> $OUT
  for N in 2 4; do
    for i in $(seq 1 $N); do $X seq $G > /tmp/isalpar_calib_$i.txt & done
    wait
    cat /tmp/isalpar_calib_*.txt | awk -F'\t' -v r=$r -v N=$N '$2=="gzip_crc"{s+=$12} END{printf "round\t%s\tcalib\tconcurrent_seq\tN\t%d\tGBps\t%.3f\n", r, N, s}' >> $OUT
    rm -f /tmp/isalpar_calib_*.txt
  done
  for t in 1 2 3 4; do
    for cfg in "lead 2 4" "spec 2 4" "lead 2 1" "spec 2 1" "lead 2 16" "spec 3 4" "lead 3 4"; do
      set -- $cfg
      $X par $G $B $3 $t $1 $2 0 - | sed "s/^/round\t$r\t/" >> $OUT
    done
  done
  $X seq $G | sed "s/^/round\t$r\tend\t/" >> $OUT
  echo "round $r done: $(cat /proc/loadavg)" >> $OUT.uptime
done
echo finished >> $OUT.uptime
