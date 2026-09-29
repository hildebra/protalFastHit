#!/bin/bash
# Resident memory over one run, by stage: mem_trace.sh LABEL DB READS_DIR THREADS [protal args]
# Samples VmRSS every 0.1 s; a sample belongs to the last stage message printed before it.
set -u
source "$(dirname "$0")/env.sh"
label=$1; db=$2; reads=$3; t=$4; shift 4
out=$PERF_DIR/runs/$label; rm -rf $out; mkdir -p $out
r1=$(ls $reads/*_R1.fq* | head -1); r2=$(ls $reads/*_R2.fq* | head -1)
t0=$(date +%s.%N)
( $BIN --db $db -1 $r1 -2 $r2 -o $out/out -t $t --no_qcmsa "$@" 2>$out/stderr.log | while IFS= read -r line; do
    printf '%.2f\t%s\n' "$(echo "$(date +%s.%N) - $t0" | bc)" "$line"; done > $out/stdout.ts ) &
sleep 0.2
pid=$(pgrep -n -f "$BIN --db $db")
while kill -0 $pid 2>/dev/null; do
  rss=$(awk '/VmRSS/{print $2}' /proc/$pid/status 2>/dev/null)
  [ -n "$rss" ] && printf '%.2f\t%s\n' "$(echo "$(date +%s.%N) - $t0" | bc)" "$rss" >> $out/rss.tsv
  sleep 0.1
done
wait
awk -F'\t' 'NR==FNR { ts[NR]=$1; msg[NR]=$2; n=NR; next }
  { while (i < n && ts[i+1] <= $1) i++; m = (i ? msg[i] : "(start)"); if ($2 > peak[m]) peak[m] = $2; if (!(m in first)) { first[m] = $1; order[++k] = m } }
  END { for (j = 1; j <= k; j++) printf "%7.2f s  %7.0f MB  %s\n", first[order[j]], peak[order[j]]/1024, substr(order[j], 1, 90) }' \
  $out/stdout.ts $out/rss.tsv
