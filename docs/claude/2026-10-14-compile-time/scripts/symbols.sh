#!/bin/bash
# Code size per function in the -O3 object (text symbols), demangled.
R=~/ctime/res
nm -S --size-sort -C $R/main_O3.o | awk '$3 ~ /[tTwW]/' > $R/nm_O3.txt
wc -l < $R/nm_O3.txt
# total text bytes
awk '{s+=strtonum("0x"$2)} END {print "text bytes:", s}' $R/nm_O3.txt
size $R/main_O3.o
tail -60 $R/nm_O3.txt | awk '{sz=strtonum("0x"$2); $1=""; $2=""; $3=""; print sz "\t" substr($0,4,200)}'
