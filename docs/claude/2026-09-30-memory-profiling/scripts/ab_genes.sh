#!/bin/bash
# Two protal builds alternated on one sample: ab_genes.sh OLD_BIN NEW_BIN DB READS_DIR THREADS ROUNDS
# Per run: seconds of the alignment stage and of the whole run (from protal's own timers) and the peak RSS.
set -u
old=$1; new=$2; db=$3; reads=$4; t=$5; n=$6
r1=$(ls $reads/*_R1.fq* | head -1); r2=$(ls $reads/*_R2.fq* | head -1)
out=${OUT:-$HOME/protal-mem/runs}/ab_genes; mkdir -p $out
secs() { awk '{ t = 0; for (i = 1; i <= NF; i++) { v = $i + 0; if ($i ~ /ms$/) t += v / 1000; else if ($i ~ /m$/) t += v * 60; else if ($i ~ /s$/) t += v } printf "%.3f", t }'; }
for i in $(seq 1 $n); do
  for v in old new; do
    bin=$old; [ $v = new ] && bin=$new
    rm -rf $out/$v; /usr/bin/time -f "%M" -o $out/rss $bin --db $db -1 $r1 -2 $r2 -o $out/$v -t $t --no_qcmsa > $out/stdout 2> $out/stderr
    a=$(grep "^Aligning reads took" $out/stdout | sed 's/Aligning reads took //' | secs)
    r=$(grep "^Run protal took" $out/stdout | sed 's/Run protal took //' | secs)
    printf '%s\t%s\talign %s s\trun %s s\tmaxrss %s KB\n' $i $v "$a" "$r" "$(tail -1 $out/rss)"
  done
done
