#!/bin/bash
# The timers of a collection's protal run (profile_all/protal.log): index load, alignment per sample with
# its reads, the SAM header, profiling.
# usage: protal_log_summary.sh PROTAL_LOG
log=$1
awk '
  function secs(text,   n, a, i, v, s) { s = 0; n = split(text, a, " ")
    for (i = 1; i <= n; i++) { v = a[i]
      if (v ~ /ms$/) { sub("ms", "", v); s += v / 1000 } else if (v ~ /[0-9]s$/) { sub("s", "", v); s += v }
      else if (v ~ /[0-9]m$/) { sub("m", "", v); s += 60 * v } else if (v ~ /[0-9]h$/) { sub("h", "", v); s += 3600 * v } }
    return s }
  /^Align the / { type = $3; sample = $8 }
  /^Aligning reads took/ { t = $0; sub("Aligning reads took ", "", t); align[type] += secs(t); n[type]++
                           printf "  align %-12s %-26s %8.2f s\n", type, sample, secs(t) }
  /^Writing the SAM header and file took/ { t = $0; sub(/.*took /, "", t); sam += secs(t) }
  /^Load Index took/ || /^Preload genomes took/ || /^Processing all samples took/ || /^Profiling took/ || /^Run protal took/ {
    t = $0; name = $0; sub(/ took.*/, "", name); sub(/.* took /, "", t); printf "%-24s %8.2f s\n", name, secs(t) }
  END { for (k in align) printf "aligning %-12s %3d samples %8.2f s\n", k, n[k], align[k]; printf "SAM headers and files    %8.2f s\n", sam }
' "$log"
