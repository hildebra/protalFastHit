#!/bin/bash
# A training table K times as large, for timing the trainer: K copies of every row, each copy's samples
# renamed (meta_design, meta_sample get _copyK) and its numeric features from column 19 on multiplied by
# 1 +- 2% noise (awk's srand: the same table on every run), so that the forests do not see exact duplicates.
# usage: enlarge_table.sh TABLE K OUT
set -euo pipefail
awk -F'\t' -v OFS='\t' -v k="$2" 'BEGIN { srand(11) }
  NR == 1 { print; next }
  { for (c = 1; c <= k; c++) {
      line = $0; split(line, f, "\t"); n = length(f)
      f[1] = f[1] "_copy" c; f[2] = f[2] "_copy" c
      for (i = 19; i <= n; i++) if (f[i] ~ /^-?[0-9.]+(e-?[0-9]+)?$/ && f[i] != 0) f[i] = f[i] * (1 + 0.04 * (rand() - 0.5))
      out = f[1]; for (i = 2; i <= n; i++) out = out OFS f[i]; print out } }' "$1" > "$3"
wc -l "$3"
