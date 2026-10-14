#!/bin/bash
# Functions emitted by several units of the split prototype (inline code each unit compiles again): by name, in how
# many units, and their code size in each.
D=~/ctime/split
cd $D
for u in main split_align split_buildindex split_options split_profile split_strains split_buildmaint; do
    nm -S -C $u.o | awk -v u=$u '$3 ~ /[tTwW]/ {sz=strtonum("0x"$2); $1=$2=$3=""; print u "\t" sz "\t" substr($0,4)}'
done > dup_all.tsv
awk -F'\t' '{n[$3]++; s[$3]+=$2} END {for (k in n) print n[k] "\t" s[k] "\t" k}' dup_all.tsv > dup_by_name.tsv
echo "functions in 1..7 units (count, total bytes over all units):"
awk -F'\t' '{c[$1]++; b[$1]+=$2} END {for (k in c) print k "\t" c[k] "\t" b[k]}' dup_by_name.tsv | sort -n
echo "total bytes: $(awk -F'\t' '{t+=$2} END {print t}' dup_all.tsv)"
echo "largest functions present in all 7 units (bytes summed over units):"
awk -F'\t' '$1 == 7' dup_by_name.tsv | sort -t$'\t' -k2 -rn | head -25 | cut -c1-180
echo "largest functions in 4-6 units:"
awk -F'\t' '$1 >= 4 && $1 < 7' dup_by_name.tsv | sort -t$'\t' -k2 -rn | head -25 | cut -c1-180
