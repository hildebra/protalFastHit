#!/usr/bin/env bash
# top.sh WORKLOAD BUILD N: the N functions with most self instructions, and inclusive for protal's functions.
W=$HOME/mt-work/onebcg; w=$1; b=$2; n=${3:-40}
echo "== $w.$b self"
grep -E '^ *[0-9,]+ \(' $W/$w.$b.txt | head -$n | sed -E 's/\[\/home[^]]*\]//; s/\(([^()]|\([^()]*\))*\)//g' | cut -c1-150
echo "== $w.$b inclusive"
callgrind_annotate --inclusive=yes --threshold=99 $W/cg.$w.$b 2>/dev/null | c++filt | grep -E '^ *[0-9,]+ \(' | head -$n | sed -E 's/\[\/home[^]]*\]//; s/\(([^()]|\([^()]*\))*\)//g' | cut -c1-150
