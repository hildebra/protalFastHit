#!/usr/bin/env bash
# xdropsum.sh: per sample, the calls with X-drop and the largest read part + reference part + budget among them
for s in pe pb ont; do
  awk -F'\t' -v s=$s '$2>0 {n++; b=($7>0?$7:0); t=length($8)+length($9)+b; if (t>m) m=t; if (t>=2*$2) p++} END {printf "%s: %d calls with X-drop, max read+reference+budget %d, at or above 2*x_drop: %d\n", s, n, m, p}' ~/wfaspeed/calls/$s.tsv
done
