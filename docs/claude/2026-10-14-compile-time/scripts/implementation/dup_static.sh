#!/bin/bash
# Local (static) functions present more than once in the split binary: copies of header functions declared static at
# namespace scope, one per unit that uses them. Their extra copies' bytes, and the largest.
B=~/split-build/build/protal
nm -S -C $B | awk '$3 == "t" {sz=strtonum("0x"$2); $1=$2=$3=""; print sz "\t" substr($0,4)}' > /tmp/locals.tsv
awk -F'\t' '{n[$2]++; s[$2]+=$1; m[$2]=$1} END {for (k in n) if (n[k] > 1) {extra += s[k] - s[k]/n[k]; print (s[k] - s[k]/n[k]) "\t" n[k] "\t" k} ; print extra "\tTOTAL" > "/dev/stderr"}' /tmp/locals.tsv | sort -rn | head -15 | cut -c1-160
echo "base binary local functions in duplicate:"
nm -S -C ~/split-base/build/protal | awk '$3 == "t" {sz=strtonum("0x"$2); $1=$2=$3=""; print sz "\t" substr($0,4)}' | awk -F'\t' '{n[$2]++; s[$2]+=$1} END {for (k in n) if (n[k] > 1) extra += s[k] - s[k]/n[k]; print extra}'
rm -f /tmp/locals.tsv
