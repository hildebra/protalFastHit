#!/bin/bash
# Median wall and user CPU (s) per configuration of matrix.tsv: summarize_matrix.sh [matrix.tsv]
f=${1:-$HOME/protal-perf/runs2/matrix.tsv}
awk -F'\t' '
function secs(s,   a, n) { n = split(s, a, ":"); return n == 2 ? a[1]*60 + a[2] : s + 0 }
function st(s,   v) { v = 0; if (match(s, /[0-9]+m /)) { v += substr(s, RSTART, RLENGTH-2) * 60 }
  if (match(s, /[0-9]+s/)) { v += substr(s, RSTART, RLENGTH-1) } if (match(s, /[0-9]+ms/)) { v += substr(s, RSTART, RLENGTH-2) / 1000 } return v }
{ label = $1; sub(/_[0-9]+$/, "", label); n[label]++; w[label, n[label]] = secs($2); u[label, n[label]] = $3 + 0; r[label, n[label]] = $5 / 1048576
  for (i = 7; i <= NF; i++) { split($i, kv, "="); if (kv[1] == "loop" || kv[1] == "reader" || kv[1] == "seed" || kv[1] == "align" || kv[1] == "loadidx") s[kv[1], label, n[label]] = st(kv[2]) } 
  if (!(label in seen)) { seen[label] = 1; order[++m] = label } }
function med(arr, label, k,   i, j, t, v) { for (i = 1; i <= k; i++) v[i] = arr[label, i]; for (i = 1; i <= k; i++) for (j = i + 1; j <= k; j++) if (v[j] < v[i]) { t = v[i]; v[i] = v[j]; v[j] = t } return v[int((k + 1) / 2)] }
function smed(key, label, k,   i, j, t, v) { for (i = 1; i <= k; i++) v[i] = s[key, label, i]; for (i = 1; i <= k; i++) for (j = i + 1; j <= k; j++) if (v[j] < v[i]) { t = v[i]; v[i] = v[j]; v[j] = t } return v[int((k + 1) / 2)] }
END { printf "%-22s %3s %8s %8s %8s %8s %8s %8s %8s %8s\n", "config", "n", "wall", "user", "loop", "reader", "seed", "align", "loadidx", "rss GB"
  for (i = 1; i <= m; i++) { l = order[i]; k = n[l]; printf "%-22s %3d %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n", l, k, med(w, l, k), med(u, l, k), smed("loop", l, k), smed("reader", l, k), smed("seed", l, k), smed("align", l, k), smed("loadidx", l, k), med(r, l, k) } }' $f
