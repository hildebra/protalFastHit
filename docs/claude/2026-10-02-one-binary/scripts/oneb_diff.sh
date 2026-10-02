#!/usr/bin/env bash
# oneb_diff.sh WORKLOAD [N], with A=BUILD B=BUILD (default base, cl): per function, self instructions of two builds
# (callgrind_annotate output of oneb_cg.sh, clone suffixes merged), the N largest differences (A minus B).
W=$HOME/mt-work/onebcg; w=$1; n=${2:-25}
awk '
  FNR == 1 { file++ }
  /^ *[0-9,]+ \(/ {
    ir = $1; gsub(",", "", ir)
    line = $0; sub(/^ *[0-9,]+ \( *[0-9.]+%\) +/, "", line); sub(/ \[[^]]*\]$/, "", line); sub(/^[^:]*:/, "", line); gsub(/ \[clone \.[a-z_0-9.]+\]/, "", line)
    if (line ~ /PROGRAM TOTALS/) { total[file] = ir; next }
    v[file, line] += ir; names[line] = 1
  }
  END {
    for (f in names) { d = v[1, f] - v[2, f]; printf "%.0f\t%.0f\t%.0f\t%.2f\t%s\n", d, v[1, f], v[2, f], 100 * v[1, f] / total[1], substr(f, 1, 140) }
    printf "TOTAL\t%.0f\t%.0f\t%.2f\n", total[1], total[2], 100 * (total[1] - total[2]) / total[1] > "/dev/stderr"
  }' $W/$w.${A:-base}.txt $W/$w.${B:-cl}.txt | sort -t$'\t' -k1,1nr | head -$n | awk -F'\t' '{ printf "%8.1fM  A %7.1fM  B %7.1fM  %5.2f%%  %s\n", $1/1e6, $2/1e6, $3/1e6, $4, $5 }'
