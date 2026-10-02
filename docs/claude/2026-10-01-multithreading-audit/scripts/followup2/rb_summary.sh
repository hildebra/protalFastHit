#!/usr/bin/env bash
# rb_summary.sh READBENCH.TSV: per input, the median pairs/s, wall and CPU seconds of one reader thread.
awk -F'\t' 'NR > 1 { k = $1; n[k]++; for (i = 7; i < NF; i++) if ($i == "pairs_per_s") p[k, n[k]] = $(i + 1); w[k, n[k]] = $3; c[k, n[k]] = $4 + $5 }
function med(arr, k, m,   i, j, t, x) { for (i = 1; i <= m; i++) x[i] = arr[k, i]; for (i = 1; i <= m; i++) for (j = i + 1; j <= m; j++) if (x[j] < x[i]) { t = x[i]; x[i] = x[j]; x[j] = t }
  return m % 2 ? x[(m + 1) / 2] : (x[m / 2] + x[m / 2 + 1]) / 2 }
END { printf "input\truns\tpairs_per_s\twall_s\tcpu_s\n"; for (k in n) printf "%s\t%d\t%.2fM\t%.2f\t%.2f\n", k, n[k], med(p, k, n[k]) / 1e6, med(w, k, n[k]), med(c, k, n[k]) }' "$1" | sort
