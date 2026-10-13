#!/usr/bin/env bash
# The true species one arm calls and another misses (docs/claude/2026-10-13-sparse-seeding), from ss_accuracy.sh's
# truth-annotated profiles (WORK/arms/<arm>/<sample>.profile.truth_annotated): per such species, the model's
# probability, the genes present, the hits and the failed candidate rate in both arms (a species with no row in the
# second arm had no hit there).
#   bash ss_lost_species.sh WORK ARM_A ARM_B
set -uo pipefail
work=$1 a=$2 b=$3
printf 'sample\tspecies\tprob_%s\tprob_%s\tgenes_%s\tgenes_%s\thits_%s\thits_%s\tfailed_rate_%s\tfailed_rate_%s\n' "$a" "$b" "$a" "$b" "$a" "$b" "$a" "$b"
for fa in "$work/arms/$a"/*.profile.truth_annotated; do
  s=$(basename "$fa" .profile.truth_annotated); fb=$work/arms/$b/$s.profile.truth_annotated
  awk -F'\t' -v s="$s" '
    FNR == 1 { f++; for (i = 1; i <= NF; i++) col[f, $i] = i; next }
    { k = $(col[f, "taxon"]); t[f, k] = $(col[f, "truth"]); p[f, k] = $(col[f, "prediction"]); pr[f, k] = $(col[f, "probability"])
      g[f, k] = $(col[f, "present_genes"]); h[f, k] = $(col[f, "total_hits"]); n[k] = $(col[f, "taxon_name"])
      fr[f, k] = ((f, "failed_candidate_rate") in col) ? $(col[f, "failed_candidate_rate"]) : "NA"; seen[k] }
    END { for (k in seen) if (t[1, k] == 1 && p[1, k] == 1 && !((2, k) in p && p[2, k] == 1)) {
            inb = ((2, k) in p)
            printf "%s\t%s\t%.3f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", s, n[k], pr[1, k], inb ? sprintf("%.3f", pr[2, k]) : "-",
                   g[1, k], inb ? g[2, k] : 0, h[1, k], inb ? h[2, k] : 0, fr[1, k], inb ? fr[2, k] : "-" } }' "$fa" "$fb"
done
