#!/usr/bin/env bash
# Seeding one read ahead with prefetching (pf.patch on HEAD): outputs against HEAD's build, then alternated whole
# runs: 500k pairs and 500k single-end reads at -t 1, 5M pairs at -t 6; five rounds; stage times.
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/1680d2f3-dfad-4cd7-9113-3d5aefceb3b5/scratchpad
BASE=$(git -C /mnt/c/Users/hildebra/Documents/locDev/protal rev-parse --short HEAD)
bash $S/build_wt.sh pfref $S/empty.patch $BASE protal || exit 1
bash $S/build_wt.sh pf $S/perf4/pf.patch $BASE protal || exit 1
W=$HOME/mt-work/pfrun; mkdir -p $W
declare -A BIN=( [ref]=$HOME/mt-work/pfref/src/build/protal [pf]=$HOME/mt-work/pf/src/build/protal )
DB=$HOME/bench071/V071/protal_db; P=$HOME/bench071/samples/points; D=$HOME/bench071/samples_deep/points; I=$HOME/mt-work/isacg
R5=$P/rl150_p500000/sim/reads/rl150_p500000_s_1; R50=$D/rl150_p5000000/sim/reads/rl150_p5000000_s_1
same() { cmp -s <(zstdcat $1) <(zstdcat $2) && echo "SAM text same" || echo "SAM DIFFER"; }
for b in ref pf; do
  rm -rf $W/i.pe.$b $W/i.se.$b $W/i.ont.$b $W/i.pb.$b $W/i.pe6.$b
  ${BIN[$b]} --db $DB -1 $I/r1.fq -2 $I/r2.fq --no_profile --prefix s -o $W/i.pe.$b -t 1 --no_qcmsa > /dev/null 2>&1
  ${BIN[$b]} --db $DB -1 ${R5}_R1.fq.gz --read_type se --no_profile --prefix s -o $W/i.se.$b -t 1 --no_qcmsa > /dev/null 2>&1
  ${BIN[$b]} --db $DB -1 $P/ont_b3000000/sim/reads/ont_b3000000_s_1.fq.gz --read_type ont --no_profile --prefix s -o $W/i.ont.$b -t 1 --no_qcmsa > /dev/null 2>&1
  ${BIN[$b]} --db $DB -1 ${R5}_R1.fq.gz -2 ${R5}_R2.fq.gz --prefix s -o $W/i.pe6.$b -t 6 --no_qcmsa > /dev/null 2>&1
done
echo "100k pairs -t 1: $(same $W/i.pe.ref/s.sam.zst $W/i.pe.pf/s.sam.zst)"
echo "single-end 500k -t 1: $(same $W/i.se.ref/s.sam.zst $W/i.se.pf/s.sam.zst)"
echo "Nanopore 3 Mb -t 1: $(same $W/i.ont.ref/s.sam.zst $W/i.ont.pf/s.sam.zst)"
cmp -s <(zstdcat $W/i.pe6.ref/s.sam.zst | sort) <(zstdcat $W/i.pe6.pf/s.sam.zst | sort) && echo "500k pairs -t 6: sorted SAM same" || echo "500k pairs -t 6: SAM DIFFER"
diff -q $W/i.pe6.ref/s.profile $W/i.pe6.pf/s.profile > /dev/null && echo "500k pairs -t 6: profile same" || echo "500k pairs -t 6: profile DIFFERS"
seconds() { awk '{ s = 0; for (i = 1; i <= NF; i++) { v = $i; if (v ~ /ms$/) { sub(/ms/, "", v); s += v / 1000 } else if (v ~ /m$/) { sub(/m/, "", v); s += v * 60 } else if (v ~ /s$/) { sub(/s/, "", v); s += v } } printf "%.3f", s }'; }
printf "case\tbuild\trep\talign_s\tseeding_s\twall_s\tuser_s\tload\n" > $W/runs.tsv
run() { local c=$1 b=$2 rep=$3 t=1 args; d=$W/run; rm -rf $d
  case $c in pe500k) args="-1 ${R5}_R1.fq.gz -2 ${R5}_R2.fq.gz";; se500k) args="-1 ${R5}_R1.fq.gz --read_type se";; pe5M) args="-1 ${R50}_R1.fq.gz -2 ${R50}_R2.fq.gz"; t=6;; esac
  local load=$(cut -d' ' -f1 /proc/loadavg)
  /usr/bin/time -f "%e %U" -o $W/time ${BIN[$b]} --db $DB $args --no_profile --prefix s -o $d -t $t --no_qcmsa > $W/log 2>&1 || echo "FAIL $c $b"
  local al=$(grep 'Aligning reads took' $W/log | sed 's/.*took //' | seconds)
  local sd=$(awk -F'\t' '$1 == "Seeding" { print $2 }' $d/misc/s_runtime.tsv)
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $c $b $rep $al $sd "$(tail -1 $W/time | tr ' ' '\t')" $load | tee -a $W/runs.tsv
}
for rep in 1 2 3 4 5; do
  for c in pe500k se500k pe5M; do
    order="ref pf"; [ $((rep % 2)) = 0 ] && order="pf ref"
    for b in $order; do run $c $b $rep; done
  done
done
awk -F'\t' 'NR > 1 { k = $1 "\t" $2; n[k]++; a[k, n[k]] = $4; s[k, n[k]] = $5 }
  END { for (k in n) { m = n[k]; for (i = 1; i <= m; i++) { x[i] = a[k, i]; y[i] = s[k, i] } asort(x); asort(y)
        printf "%s\talign fastest %.2f median %.2f\tseeding fastest %.2f median %.2f\n", k, x[1], x[int((m + 1) / 2)], y[1], y[int((m + 1) / 2)]; delete x; delete y } }' $W/runs.tsv | sort
echo PF DONE
