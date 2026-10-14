#!/bin/bash
# Medians of wall, user, sys and peak RSS per load configuration in runs/load_<tag>/loads.tsv.
F=~/sparseprof/runs/load_${1:-base}/loads.tsv
awk -F'\t' '
function med(s,   a, n, i, j, t) { n = split(s, a, " "); for (i = 1; i <= n; i++) for (j = i + 1; j <= n; j++) if (a[j] < a[i]) { t = a[i]; a[i] = a[j]; a[j] = t }; return a[int((n + 1) / 2)] }
{
  db = "grid"; off = 0
  if ($1 == "r226reduced") { db = "r226red"; off = 1 }
  if ($(1 + off) != "load") next
  level = $(2 + off); compact = $(3 + off); threads = $(4 + off)
  wall = $(5 + off); user = $(6 + off); sys = $(7 + off); hwm = $(10 + off)
  sub(/wall /, "", wall); sub(/ s/, "", wall); sub(/user /, "", user); sub(/ s/, "", user); sub(/sys /, "", sys); sub(/ s/, "", sys)
  gsub(/[^0-9]/, "", hwm)
  bits = ($(11 + off) ~ /grid level 3\)/) ? "lvbits" : "nobits"
  key = db "\t" level "\t" compact "\t" bits "\t" threads
  W[key] = W[key] " " wall; U[key] = U[key] " " user; S[key] = S[key] " " sys; H[key] = H[key] " " hwm / 1048576
  if (!(key in seen)) { seen[key] = 1; order[++n] = key }
}
END {
  print "db\tlevel\tcompact\tbits\tthreads\twall_med\tuser_med\tsys_med\tpeak_GB_med"
  for (i = 1; i <= n; i++) { k = order[i]; printf "%s\t%.2f\t%.2f\t%.2f\t%.2f\n", k, med(W[k]), med(U[k]), med(S[k]), med(H[k]) }
}' $F
