#!/usr/bin/env bash
# Fastest of several alternated runs (the same run varies 2x on this laptop: fast or low-power cores, throttling):
# 500k pairs at -t 1 (8 rounds), 5M pairs at -t 6 (6 rounds); alignment stage time ("Aligning reads took").
W=$HOME/mt-work/pfrun
declare -A BIN=( [ref]=$HOME/mt-work/pfref/src/build/protal [pf]=$HOME/mt-work/pf/src/build/protal )
DB=$HOME/bench071/V071/protal_db; P=$HOME/bench071/samples/points; D=$HOME/bench071/samples_deep/points
R5=$P/rl150_p500000/sim/reads/rl150_p500000_s_1; R50=$D/rl150_p5000000/sim/reads/rl150_p5000000_s_1
seconds() { awk '{ s = 0; for (i = 1; i <= NF; i++) { v = $i; if (v ~ /ms$/) { sub(/ms/, "", v); s += v / 1000 } else if (v ~ /m$/) { sub(/m/, "", v); s += v * 60 } else if (v ~ /s$/) { sub(/s/, "", v); s += v } } printf "%.3f", s }'; }
printf "case\tbuild\trep\talign_s\tuser_s\n" > $W/runs2.tsv
run() { local c=$1 b=$2 rep=$3 t=1 args; d=$W/run2; rm -rf $d
  [ $c = pe500k ] && args="-1 ${R5}_R1.fq.gz -2 ${R5}_R2.fq.gz" || { args="-1 ${R50}_R1.fq.gz -2 ${R50}_R2.fq.gz"; t=6; }
  /usr/bin/time -f "%U" -o $W/time2 ${BIN[$b]} --db $DB $args --no_profile --prefix s -o $d -t $t --no_qcmsa > $W/log2 2>&1 || echo "FAIL $c $b"
  printf "%s\t%s\t%s\t%s\t%s\n" $c $b $rep $(grep 'Aligning reads took' $W/log2 | sed 's/.*took //' | seconds) $(tail -1 $W/time2) >> $W/runs2.tsv
}
for rep in 1 2 3 4 5 6 7 8; do o="ref pf"; [ $((rep % 2)) = 0 ] && o="pf ref"; for b in $o; do run pe500k $b $rep; done; done
for rep in 1 2 3 4 5 6; do o="ref pf"; [ $((rep % 2)) = 0 ] && o="pf ref"; for b in $o; do run pe5M $b $rep; done; done
awk -F'\t' 'NR > 1 { k = $1 "\t" $2; n[k]++; a[k, n[k]] = $4; u[k, n[k]] = $5 }
  END { for (k in n) { m = n[k]; for (i = 1; i <= m; i++) { x[i] = a[k, i]; y[i] = u[k, i] } asort(x); asort(y)
        printf "%s\talign fastest %.2f second %.2f median %.2f\tuser fastest %.1f median %.1f\t(n=%d)\n", k, x[1], x[2], x[int((m + 1) / 2)], y[1], y[int((m + 1) / 2)], m; delete x; delete y } }' $W/runs2.tsv | sort
echo PF2 DONE
