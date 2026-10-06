#!/usr/bin/env bash
# cmp_runs.sh A B [label_a label_b]: two measure_performance.sh result directories side by side: the per-thread stage
# medians over each one's runs (stages.tsv), and every run's totals (runs.tsv). Used for the twelfth against the
# eleventh run (results/cluster_v12 against ../2026-10-04-performance-gtdb-scale/results_v11) and the thirteenth against
# the twelfth (results/cluster_v13 against results/cluster_v12).
set -uo pipefail
A=$1; B=$2; LA=${3:-a}; LB=${4:-b}
med() {  # stages.tsv -> sample \t stage \t median seconds per thread
  awk -F'\t' 'NR > 1 { k = $1 "\t" $3; v[k, ++n[k]] = $6; if (!(k in seen)) { seen[k] = 1; order[++m] = k } }
    END { for (j = 1; j <= m; j++) { k = order[j]; c = n[k]; for (i = 1; i <= c; i++) a[i] = v[k, i]
            for (i = 2; i <= c; i++) { x = a[i]; for (l = i - 1; l >= 1 && a[l] > x; l--) a[l + 1] = a[l]; a[l + 1] = x }
            printf "%s\t%.3f\n", k, c % 2 ? a[(c + 1) / 2] : (a[c / 2] + a[c / 2 + 1]) / 2 } }' "$1"
}
ta=$(mktemp); tb=$(mktemp)
med $A/stages.tsv > $ta
med $B/stages.tsv > $tb
printf '%-6s %-45s %8s %8s %8s\n' sample stage "$LA" "$LB" change
awk -F'\t' 'NR == FNR { a[$1 "\t" $2] = $3; next } { k = $1 "\t" $2; o = (k in a) ? a[k] : "NA";
    d = (o != "NA" && o > 0) ? sprintf("%+.0f%%", 100 * ($3 - o) / o) : "";
    printf "%-6s %-45s %8s %8.3f %8s\n", $1, $2, o, $3, d }' $ta $tb | grep -v "^cohort"
rm -f $ta $tb
echo
for f in $A/runs.tsv $B/runs.tsv; do
  echo "== $f"
  awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
    { printf "%s rep %s: wall %s user %s rss %s load_index %s aligning %s profiling %s instructions %.3g cycles %.3g ipc %s\n",
             $1, $3, $h["wall_s"], $h["user_s"], $h["max_rss_gb"], $h["load_index_s"], $h["aligning_s"], $h["profiling_s"],
             $h["instructions"], $h["cycles"], $h["ipc"] }' $f
done
