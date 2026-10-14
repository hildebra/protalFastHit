#!/bin/bash
# Medians per configuration and path (old/new) of runs/load_ab/loads.tsv: wall, user, sys, peak.
F=~/sparseprof/runs/load_ab/loads.tsv
awk -F'\t' '
function med(s,   a, n, i, j, t) { n = split(s, a, " "); for (i = 1; i <= n; i++) for (j = i + 1; j <= n; j++) if (a[j] + 0 < a[i] + 0) { t = a[i]; a[i] = a[j]; a[j] = t }; return a[int((n + 1) / 2)] }
$3 == "load" {
  db = $1; old = $2; level = $4; compact = $5; threads = $6
  wall = $7; user = $8; sys = $9; hwm = $12
  sub(/wall /, "", wall); sub(/ s/, "", wall); sub(/user /, "", user); sub(/ s/, "", user); sub(/sys /, "", sys); sub(/ s/, "", sys)
  gsub(/[^0-9]/, "", hwm); sub(/level /, "", level); sub(/compact /, "", compact); sub(/threads /, "", threads)
  key = db "\tL" level "\tc" compact "\tt" threads "\t" old
  W[key] = W[key] " " wall; U[key] = U[key] " " user; S[key] = S[key] " " sys; H[key] = H[key] " " hwm / 1048576
  if (!(key in seen)) { seen[key] = 1; order[++n] = key }
}
END {
  print "db\tlevel\tcompact\tthreads\tpath\twall\tuser\tsys\tpeak_GB"
  for (i = 1; i <= n; i++) { k = order[i]; printf "%s\t%.2f\t%.2f\t%.2f\t%.2f\n", k, med(W[k]), med(U[k]), med(S[k]), med(H[k]) }
}' $F
