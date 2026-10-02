#!/usr/bin/env bash
# Whole runs of four builds of d11381f's code, alternated, one thread, niced: base = plain x86-64, cl = the
# target_clones build, v3 = protal_avx2 (all of main.cpp for x86-64-v3), v2 = x86-64-v2 baseline. Aligning 500k
# pairs and Nanopore 90 Mb (--no_profile), profiling the 500k-pair SAM. Five rounds. Rows into
# ~/mt-work/onebtime/runs.tsv; the first round's outputs are kept and compared with base's.
set -uo pipefail
W=$HOME/mt-work/onebtime; rm -rf $W; mkdir -p $W
declare -A BIN=( [base]=$HOME/mt-work/oneb-ref/src/build/protal [v3]=$HOME/mt-work/oneb-ref/src/build/protal_avx2
                 [cl]=$HOME/mt-work/oneb/src/build/protal [v2]=$HOME/mt-work/oneb-v2/src/build/protal )
DB=$HOME/bench071/V071/protal_db; P=$HOME/bench071/samples/points
declare -A ARGS=(
  [pe500k]="-1 $P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz -2 $P/rl150_p500000/sim/reads/rl150_p500000_s_1_R2.fq.gz --no_profile"
  [ont90M]="-1 $P/ont_b90000000/sim/reads/ont_b90000000_s_1.fq.gz --read_type ont --no_profile"
  [prof500k]="--profile_only $HOME/bench071/runs/v071.full.pe.rl150_p500000_s_1/rl150_p500000_s_1.sam.zst" )
seconds() { awk '{ s = 0; for (i = 1; i <= NF; i++) { v = $i; if (v ~ /ms$/) { sub(/ms/, "", v); s += v / 1000 } else if (v ~ /m$/) { sub(/m/, "", v); s += v * 60 } else if (v ~ /s$/) { sub(/s/, "", v); s += v } } printf "%.3f", s }'; }
printf "case\tbinary\trep\tstage_s\twall_s\tuser_s\tload\n" > $W/runs.tsv
run() {  # run CASE BINARY REP
  local c=$1 b=$2 rep=$3 d=$W/run.$1.$2; rm -rf $d
  local load=$(cut -d' ' -f1 /proc/loadavg)
  nice -n 5 /usr/bin/time -f "%e %U" -o $W/time ${BIN[$b]} --db $DB ${ARGS[$c]} --prefix s -o $d -t 1 --no_qcmsa > $W/log 2>&1 || echo "FAIL $c $b"
  local stage=$(grep -E 'Aligning reads took|Profiling took' $W/log | head -1 | sed 's/.*took //' | seconds)
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" $c $b $rep $stage "$(tail -1 $W/time | tr ' ' '\t')" $load | tee -a $W/runs.tsv
  [ $rep = 1 ] && { rm -rf $W/out.$c.$b; mv $d $W/out.$c.$b; }
}
for rep in 1 2 3 4 5; do
  for c in pe500k ont90M prof500k; do for b in base cl v3 v2; do run $c $b $rep; done; done
done
for c in pe500k ont90M prof500k; do
  for b in cl v3 v2; do
    r=$(diff -r -q -x '*_runtime.tsv' -x '*.sam.zst' $W/out.$c.base $W/out.$c.$b > /dev/null && echo same || echo DIFFER)
    s=""; [ -f $W/out.$c.base/s.sam.zst ] && { cmp -s <(zstdcat $W/out.$c.base/s.sam.zst) <(zstdcat $W/out.$c.$b/s.sam.zst) && s="SAM text same" || s="SAM text DIFFER"; }
    echo "$c $b vs base: outputs $r $s"
  done
done
# medians per case and build
awk -F'\t' 'NR > 1 { k = $1 "\t" $2; n[k]++; st[k, n[k]] = $4; us[k, n[k]] = $6 }
  END { for (k in n) { m = n[k]; for (i = 1; i <= m; i++) { a[i] = st[k, i]; b[i] = us[k, i] }
        asort(a); asort(b); printf "%s\tstage median %.2f\tuser median %.2f\t(n=%d)\n", k, a[int((m + 1) / 2)], b[int((m + 1) / 2)], m } }' $W/runs.tsv | sort
echo ONEBTIME DONE
