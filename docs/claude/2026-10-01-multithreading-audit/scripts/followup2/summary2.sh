#!/usr/bin/env bash
# summary2.sh RUNS2.TSV: per db, input, threads, batch: runs, median and range of aligning seconds, median pairs/s,
# median CPU seconds (user + sys), median reader seconds per thread.
awk -F'\t' 'NR > 1 { k = $1 "\t" $2 "\t" $3 "\t" $4; n[k]++; a[k, n[k]] = $6; c[k, n[k]] = $9 + $10; r[k, n[k]] = $12 }
function med(arr, k, m,   i, j, t, x) { for (i = 1; i <= m; i++) x[i] = arr[k, i]; for (i = 1; i <= m; i++) for (j = i + 1; j <= m; j++) if (x[j] < x[i]) { t = x[i]; x[i] = x[j]; x[j] = t }
  lo = x[1]; hi = x[m]; return m % 2 ? x[(m + 1) / 2] : (x[m / 2] + x[m / 2 + 1]) / 2 }
END {
  printf "db\tinput\tthreads\tbatch\truns\taligning_s\trange\tpairs_per_s\tcpu_s\treader_s\n"
  for (k in n) { m = n[k]; am = med(a, k, m); alo = lo; ahi = hi; cm = med(c, k, m); rm = med(r, k, m)
    printf "%s\t%d\t%.2f\t%.2f-%.2f\t%.0f\t%.1f\t%.2f\n", k, m, am, alo, ahi, 4995812 / am, cm, rm }
}' "$1" | (read -r h; echo "$h"; sort -t$'\t' -k1,1r -k3,3n -k4,4n -k2,2)
