#!/usr/bin/env bash
# demand_summary.sh DEMAND.TSV: per emulated aligner count and batch, the mean delivered pairs/s of each input over
# the repetitions, and BGZF and gzip relative to uncompressed.
awk -F'\t' 'NR > 1 { k = $3 "\t" $2; s[k, $1] += $7; n[k, $1]++; d[k] = $6; keys[k] = 1 }
END {
  printf "aligners\tbatch\tdemand\tbgzf\tgzip\tplain\tbgzf/plain\tgzip/plain\tplain/demand\n"
  for (k in keys) {
    b = s[k, "bgzf"] / n[k, "bgzf"]; g = s[k, "gzip"] / n[k, "gzip"]; p = s[k, "plain"] / n[k, "plain"]
    printf "%s\t%s\t%.2fM\t%.2fM\t%.2fM\t%.2f\t%.2f\t%.2f\n", k, d[k] / 1e6, b / 1e6, g / 1e6, p / 1e6, b / p, g / p, p / d[k]
  }
}' "$1" | (read -r h; echo "$h"; sort -t$'\t' -k1,1n -k2,2n)
