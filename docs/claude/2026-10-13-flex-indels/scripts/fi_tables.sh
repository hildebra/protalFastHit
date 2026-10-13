#!/usr/bin/env bash
# The report's tables (docs/claude/2026-10-13-flex-indels) from results/: the reads with a record of their own species
# by arm and error rate, what the gate and sparse seeding lose against the ungated run and the share --flex_indels wins
# back, and the gate's score as the alternative. Markdown to stdout.
#   bash fi_tables.sh RESULTS
set -uo pipefail
res=$1
accs=() files=()
[ -s "$res/sweep3_0.99/reads.tsv" ] && { accs+=(0.99); files+=("$res/sweep3_0.99/reads.tsv"); }
for a in 0.97 0.93 0.88; do accs+=("$a"); files+=("$res/sweep1/reads_$a.tsv"); done
rate() {  # the measured error rate of an accuracy setting (results/error_rates.tsv)
  awk -F'\t' -v a="$1" '$1 == a { printf "%.1f%%", 100 * $2 }' "$res/error_rates.tsv"
}
header="| arm |"; rule="|---|"
for a in "${accs[@]}"; do header+=" $(rate "$a") errors ($a) |"; rule+="---|"; done

echo "$header"; echo "$rule"
for arm in plain flex c1 c1f c8 c8f c16 c16f c32 c32f c8s c8sf c16s c16sf c32s c32sf; do
  row="| $arm |"
  for f in "${files[@]}"; do row+=" $(awk -F'\t' -v a="$arm" '$1 == a { printf "%d", $8 }' "$f") |"; done
  echo "$row"
done
echo
header="| reads lost against no gate: without / with the option (won back) |"; rule="|---|"
for a in "${accs[@]}"; do header+=" $(rate "$a") |"; rule+="---|"; done
echo "$header"; echo "$rule"
for arm in c1 c8 c16 c32 c8s c16s c32s; do
  case $arm in c1) name="gate, C = 1" ;; c*s) name="sparse seeding, C = ${arm:1:-1}" ;; *) name="gate, C = ${arm:1}" ;; esac
  row="| $name |"
  for f in "${files[@]}"; do
    row+=" $(awk -F'\t' -v a="$arm" '$1 == "plain" { p = $8 } $1 == a { x = $8 } $1 == a "f" { y = $8 }
             END { l = p - x; m = p - y; if (l >= 20) printf "%d / %d (%.0f%%)", l, m, 100 * (l - m) / l; else printf "%d / %d", l, m }' "$f") |"
  done
  echo "$row"
done
echo
echo "| arm | reads kept by the gate | with an own-species record | more kept than at score 14 | more with a record | per read kept more |"
echo "|---|---|---|---|---|---|"
for c in 16 32; do
  for arm in c$c c${c}f c${c}t13 c${c}t12 c${c}t11; do
    awk -F'\t' -v a="$arm" -v base="c$c" '
      FILENAME == ARGV[1] { if ($1 == a) { own = $8; reads = $5 } if ($1 == base) base_own = $8; next }
      $1 == a { kept = reads - $5 } $1 == base { base_kept = reads - $5 }
      END { label = a; sub(/f$/, " (--flex_indels)", label); sub(/t/, ", score ", label)
            if (a == base) printf "| %s (score 14) | %d | %d | | | |\n", label, kept, own
            else printf "| %s | %d | %d | %+d | %+d | %.2f |\n", label, kept, own, kept - base_kept, own - base_own, (own - base_own) / (kept - base_kept) }' \
      "$res/sweep2_score/reads_0.88.tsv" "$res/sweep2_score/runs_0.88.tsv"
  done
done
