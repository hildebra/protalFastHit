#!/usr/bin/env bash
# The long-read change on the long-read benchmark (the samples 0.7.2's collector made, $B/samples_lr072: PacBio and
# Nanopore, 3 and 90 Mb, 4 samples each) against the full V072 database, -t 6: the reference (HEAD) and the
# prototype, alternated per sample; runs into ~/mt-work/perf3/lrbench/<build>.<scenario>_s_<n>, times in times.tsv.
B=$HOME/bench071; W=$HOME/mt-work/perf3/lrbench; mkdir -p $W
declare -A BIN=( [ref]=$HOME/mt-work/perf3/ref/build/protal [lr]=$HOME/mt-work/perf3-lr/src/build/protal )
[ -f $W/times.tsv ] || printf "build\tsample\twall_s\tuser_s\talign_s\tload\n" > $W/times.tsv
seconds() { awk '{ s = 0; for (i = 1; i <= NF; i++) { v = $i; if (v ~ /ms$/) { sub(/ms/, "", v); s += v / 1000 } else if (v ~ /m$/) { sub(/m/, "", v); s += v * 60 } else if (v ~ /s$/) { sub(/s/, "", v); s += v } } printf "%.3f", s }'; }
for scen in ont_b3000000 pb_b3000000 ont_b90000000 pb_b90000000; do
  t=ont; [[ $scen == pb_* ]] && t=pb
  for n in 1 2 3 4; do
    order="ref lr"; [ $((n % 2)) = 0 ] && order="lr ref"
    for b in $order; do
      name=$b.${scen}_s_$n; [ -f $W/$name/s.profile ] && continue
      rm -rf $W/$name; load=$(cut -d' ' -f1 /proc/loadavg)
      nice /usr/bin/time -f "%e %U" -o $W/$name.time ${BIN[$b]} --db $B/V072/protal_db -1 $B/samples_lr072/points/$scen/sim/reads/${scen}_s_$n.fq.gz \
        --read_type $t --prefix s -o $W/$name -t 6 --no_qcmsa > $W/$name.log 2>&1 || echo "FAIL $name"
      al=$(grep 'Aligning reads took' $W/$name.log | sed 's/.*took //' | seconds)
      printf "%s\t%s\t%s\t%s\t%s\n" $b ${scen}_s_$n "$(tail -1 $W/$name.time | tr ' ' '\t')" $al $load >> $W/times.tsv
    done
  done
done
echo LRBENCH DONE
