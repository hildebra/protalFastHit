#!/bin/bash
# Medians of perf.sh's runs: the index load's seconds ("Load Index took"), and per seeding run the Seeding row of
# misc/<prefix>_runtime.tsv over the k-mers looked up (the log's seeding line): ns per lookup.
OUT=~/keysfix/perf
cd $OUT
secs() {  # "Load Index took 1m 2s 345ms" -> seconds
  grep -m1 "^Load Index took" "$1" | sed 's/^Load Index took //' | awk '{t = 0; for (i = 1; i <= NF; i++) { v = $i; if (v ~ /ms$/) t += substr(v, 1, length(v) - 2) / 1000; else if (v ~ /h$/) t += substr(v, 1, length(v) - 1) * 3600; else if (v ~ /m$/) t += substr(v, 1, length(v) - 1) * 60; else if (v ~ /s$/) t += substr(v, 1, length(v) - 1) } print t}'
}
med() { sort -g | awk '{a[NR] = $1} END {if (NR) print a[int((NR + 1) / 2)]; else print "NA"}'; }
echo -e "case\tbase\tnew"
for c in whole sparse sparsec; do
  b=$(for r in 1 2 3; do secs load_base_${c}_$r.log; done | med)
  n=$(for r in 1 2 3; do secs load_new_${c}_$r.log; done | med)
  echo -e "load_${c}_s\t$b\t$n"
done
ns() {  # name prefix -> ns per lookup
  local kmers seeding
  kmers=$(grep -m1 "seeding: " "$1.log" | sed 's/.*seeding: \([0-9]*\) k-mers looked up.*/\1/')
  seeding=$(awk -F'\t' '$1 == "Seeding" {print $2}' "$1/misc/$2_runtime.tsv" 2>/dev/null | head -1)
  awk -v s="$seeding" -v k="$kmers" 'BEGIN {if (k > 0 && s != "") printf "%.2f\n", 1e9 * s / k; else print "NA"}'
}
for t in pe hifi; do
  for keys in direct compact; do
    b=$(for r in 1 2 3; do ns seed_base_${t}_${keys}_$r $t; done | med)
    n=$(for r in 1 2 3; do ns seed_new_${t}_${keys}_$r $t; done | med)
    echo -e "seed_${t}_${keys}_ns_per_lookup\t$b\t$n"
  done
done
grep -h "seeding: " seed_new_pe_direct_1.log seed_new_hifi_direct_1.log | cut -c1-160
grep -c "exit=0" runs.txt; grep -v "exit=0" runs.txt
