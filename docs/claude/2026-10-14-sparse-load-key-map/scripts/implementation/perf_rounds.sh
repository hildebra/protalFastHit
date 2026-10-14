#!/bin/bash
# Every round's value of perf.sh's runs (load seconds, seeding ns per k-mer), to see the spread.
cd ~/keysfix/perf
for c in whole sparse sparsec; do
  for b in base new; do
    printf "load_%s_%s" $c $b; for r in 1 2 3; do printf "\t%s" "$(grep -m1 '^Load Index took' load_${b}_${c}_$r.log | sed 's/Load Index took //')"; done; echo
  done
done
for t in pe hifi; do for keys in direct compact; do for b in base new; do
  printf "seed_%s_%s_%s" $t $keys $b
  for r in 1 2 3; do
    n=seed_${b}_${t}_${keys}_$r
    k=$(grep -m1 "seeding: " $n.log | sed 's/.*seeding: \([0-9]*\) k-mers looked up.*/\1/')
    s=$(awk -F'\t' '$1 == "Seeding" {print $2}' $n/misc/${t}_runtime.tsv | head -1)
    printf "\t%s" "$(awk -v s=$s -v k=$k 'BEGIN {printf "%.1f", 1e9 * s / k}')"
  done; echo
done; done; done
for b in base new; do for keys in direct compact; do printf "records_pe_%s_%s\t" $b $keys; grep -vc "^@" seed_${b}_pe_${keys}_1/pe.sam 2>/dev/null || find seed_${b}_pe_${keys}_1 -name "*.sam*" | head -2; done; done
