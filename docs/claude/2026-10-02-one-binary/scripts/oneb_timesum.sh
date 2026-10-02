#!/usr/bin/env bash
# Per case and build of oneb_time.sh: the fastest and the median stage time, the five stage times, the load range.
T=$HOME/mt-work/onebtime/runs.tsv
cp $T /mnt/c/Users/hildebra/Documents/locDev/protal-clones/docs/claude/2026-10-02-one-binary/results/timing_runs.tsv
awk -F'\t' 'NR > 1 { k = $1 "\t" $2; n[k]++; v[k, n[k]] = $4; s[k] = s[k] " " $4 }
  END { for (k in n) { m = n[k]; for (i = 1; i <= m; i++) a[i] = v[k, i]; asort(a); printf "%s\tmin %.2f\tmedian %.2f\t%s\n", k, a[1], a[int((m + 1) / 2)], s[k]; delete a } }' $T | sort
awk -F'\t' 'NR > 1 { print $7 }' $T | sort -n | sed -n '1p;$p' | tr '\n' ' '; echo
