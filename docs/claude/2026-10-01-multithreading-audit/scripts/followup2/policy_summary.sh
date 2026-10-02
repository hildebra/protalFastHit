#!/usr/bin/env bash
# policy_summary.sh TSV: per db, input, threads and the 4th column (policy or lock): runs, median aligning seconds
# and range, median CPU seconds, median reader seconds per thread.
awk -F'\t' 'NR > 1 { k = $1 "\t" $2 "\t" $3 "\t" $4; n[k]++; a[k, n[k]] = $6; c[k, n[k]] = $9 + $10; r[k, n[k]] = $12 }
function med(arr, k, m,   i, j, t, x) { for (i = 1; i <= m; i++) x[i] = arr[k, i]; for (i = 1; i <= m; i++) for (j = i + 1; j <= m; j++) if (x[j] < x[i]) { t = x[i]; x[i] = x[j]; x[j] = t }
  lo = x[1]; hi = x[m]; return m % 2 ? x[(m + 1) / 2] : (x[m / 2] + x[m / 2 + 1]) / 2 }
END {
  printf "db\tinput\tthreads\tvariant\truns\taligning_s\trange\tcpu_s\treader_s\n"
  for (k in n) { m = n[k]; am = med(a, k, m); alo = lo; ahi = hi
    printf "%s\t%d\t%.2f\t%.2f-%.2f\t%.1f\t%.2f\n", k, m, am, alo, ahi, med(c, k, m), med(r, k, m) }
}' "$1" | (read -r h; echo "$h"; sort -t$'\t' -k1,1r -k3,3nr -k2,2 -k4,4)
