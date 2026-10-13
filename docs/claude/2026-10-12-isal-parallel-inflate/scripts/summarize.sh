#!/bin/bash
# median / min / max GB/s per configuration of a timing log
awk -F'\t' '
function med(s,   a, n, i, j, t) { n = split(s, a, " "); for (i = 1; i <= n; i++) for (j = i + 1; j <= n; j++) if (a[j] + 0 < a[i] + 0) { t = a[i]; a[i] = a[j]; a[j] = t }; return (n % 2) ? a[(n + 1) / 2] : (a[n / 2] + a[n / 2 + 1]) / 2 }
function mn(s,  a, n, i, m) { n = split(s, a, " "); m = a[1]; for (i = 2; i <= n; i++) if (a[i] + 0 < m + 0) m = a[i]; return m }
function mx(s,  a, n, i, m) { n = split(s, a, " "); m = a[1]; for (i = 2; i <= n; i++) if (a[i] + 0 > m + 0) m = a[i]; return m }
{
  delete kv
  if ($3 == "end") { key = "seq_end_" $5; start = 6 } else if ($3 == "seq") { key = "seq_" $4; start = 5 } else if ($3 == "calib") { start = 5 } else { start = 4 }
  for (i = start; i < NF; i += 2) kv[$i] = $(i + 1)
  if ($3 == "calib") key = "calib_N" kv["N"]; if ($3 == "par") key = kv["sched"] "_v" kv["variant"] "_" kv["piece_MiB"] "MiB_t" kv["threads"]
  if (!(key in seen)) { seen[key] = 1; order[++no] = key }
  g[key] = g[key] " " kv["GBps"]; c[key] = c[key] " " kv["cpu_s_per_GB"]; m[key] = m[key] " " kv["maxrss_MiB"]
  if ($3 == "par") { ok[key] = ok[key] kv["OK"]; sp[key] = kv["spec_pieces"]; ld[key] = kv["lead_pieces"]; w[key] = w[key] " " kv["main_wait_s"] }
}
END {
  printf "config\tGBps_median\tGBps_min\tGBps_max\tcpu_s_per_GB_median\tmaxrss_MiB_max\tspec_pieces\tlead_pieces\tmain_wait_s_median\tOK\n"
  for (i = 1; i <= no; i++) { k = order[i]; printf "%s\t%.3f\t%.3f\t%.3f\t%.3f\t%s\t%s\t%s\t%s\t%s\n", k, med(g[k]), mn(g[k]), mx(g[k]), med(c[k]), mx(m[k]), sp[k], ld[k], (k in w) ? med(w[k]) : "", ok[k] }
}' "$1"
