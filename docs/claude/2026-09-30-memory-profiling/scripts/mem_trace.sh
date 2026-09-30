#!/bin/bash
# Memory over one protal run, by stage: mem_trace.sh LABEL DB READS_DIR THREADS [protal args]
# MAP=file runs a map file (several samples) instead of READS_DIR.
# Every 0.1 s: VmRSS, anonymous vs file-backed resident memory, huge pages and VmHWM from /proc;
# a sample belongs to the last stage message protal printed before it. Output: $OUT/LABEL/
# (stdout.ts, mem.tsv = seconds, rss_kb, anon_kb, file_kb, hwm_kb; stage summary on stdout).
set -u
BIN=${BIN:-$HOME/protal-mem/build/protal}
OUT=${OUT:-$HOME/protal-mem/runs}
label=$1; db=$2; reads=$3; t=$4; shift 4
out=$OUT/$label; rm -rf $out; mkdir -p $out
if [ -n "${MAP:-}" ]; then inputs="--map $MAP"; else
  r1=$(ls $reads/*_R1.fq* | head -1); r2=$(ls $reads/*_R2.fq* | head -1); inputs="-1 $r1 -2 $r2"; fi
t0=$(date +%s.%N)
( $BIN --db $db $inputs -o $out/out -t $t --no_qcmsa "$@" 2>$out/stderr.log | while IFS= read -r line; do
    printf '%.2f\t%s\n' "$(echo "$(date +%s.%N) - $t0" | bc)" "$line"; done > $out/stdout.ts ) &
sleep 0.2
pid=$(pgrep -n -f "$BIN --db $db")
while kill -0 $pid 2>/dev/null; do
  s=$(awk '/^VmRSS/{r=$2} /^RssAnon/{a=$2} /^RssFile/{f=$2} /^VmHWM/{h=$2} END{print r"\t"a"\t"f"\t"h}' /proc/$pid/status 2>/dev/null)
  [ -n "$s" ] && printf '%.2f\t%s\n' "$(echo "$(date +%s.%N) - $t0" | bc)" "$s" >> $out/mem.tsv
  sleep 0.1
done
wait
awk -F'\t' 'NR==FNR { ts[NR]=$1; msg[NR]=$2; n=NR; next }
  { while (i < n && ts[i+1] <= $1) i++; m = (i ? msg[i] : "(start)"); if ($2 > peak[m]) peak[m] = $2; if (!(m in first)) { first[m] = $1; order[++k] = m } if ($5>hwm) hwm=$5 }
  END { for (j = 1; j <= k; j++) printf "%7.2f s  %7.0f MB  %s\n", first[order[j]], peak[order[j]]/1024, substr(order[j], 1, 90); printf "VmHWM %.0f MB\n", hwm/1024 }' \
  $out/stdout.ts $out/mem.tsv
