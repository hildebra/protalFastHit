#!/bin/bash
# Genomes and species per simulated sample of collections (their points' manifest.tsv), and the contigs of
# those genomes (pad_genomes.sh's plan): ART and pbsim3 run once per genome of a sample.
# usage: genomes_per_sample.sh PLAN COLLECTION_DIR...
plan=$1; shift
for dir in "$@"; do
  cat "$dir"/points/rl*/sim/manifest.tsv | awk -F'\t' -v d="$dir" '
    NR == FNR { contigs[$1] = $4; next }
    $1 == "sample" { next }
    { g[$1]++; sp[$1 SUBSEP $3] = 1; acc = $2; sub(/_genomic.*/, "", acc); c[$1] += contigs[acc] }
    END { for (s in g) { n++; G += g[s]; C += c[s] }
          for (k in sp) S++
          printf "%s: %d samples, %.1f genomes and %.1f species per sample, %.0f contigs per genome\n", d, n, G / n, S / n, C / G }' "$plan" -
done
