#!/bin/bash
# Sums the timed ART and pbsim calls (timed_tool.sh) per design point: calls, wall and CPU seconds, the reads
# or bases asked for, and the cost per call and per million reads (a least-squares fit wall = a + b * reads
# over the point's calls).
# usage: tool_summary.sh TOOL_TIMES PLAN   (PLAN: pad_genomes.sh's plan.tsv, for the genome lengths)
set -euo pipefail
awk -F'\t' '
  NR == FNR { len[$1] = $3; next }
  {
    tool = $1; wall = $2; cpu = $3 + $4; n = split($5, a, " "); genome = ""; cov = 0; point = ""
    for (i = 1; i < n; i++) {
      if (a[i] == "-i" || a[i] == "--genome") genome = a[i + 1]
      if (a[i] == "-f" || a[i] == "--depth") cov = a[i + 1]
      if (a[i] == "-l") rl = a[i + 1]
    }
    if (match(genome, /points\/[^\/]+\//)) point = substr(genome, RSTART + 7, RLENGTH - 8)
    acc = genome; sub(/.*\//, "", acc); sub(/_genomic.*/, "", acc)
    glen = len[acc] + 0
    if (glen == 0 && match(genome, /tmp\/[^\/]+\/[0-9]+\//)) glen = 3.58e6  # pbsim: genome.fna copies
    units = tool == "pbsim" ? cov * glen / 1e6 : cov * glen / (2 * rl) / 1e6  # Mb of long reads; M read pairs
    k = tool " " point
    calls[k]++; W[k] += wall; C[k] += cpu; U[k] += units
    sx[k] += units; sy[k] += wall; sxx[k] += units * units; sxy[k] += units * wall
  }
  END {
    printf "%-26s %6s %9s %9s %11s %9s %12s\n", "tool point", "calls", "wall_s", "cpu_s", "units", "s/call", "s/unit(fit)"
    for (k in calls) {
      d = calls[k] * sxx[k] - sx[k] * sx[k]
      b = d > 0 ? (calls[k] * sxy[k] - sx[k] * sy[k]) / d : 0
      icpt = (sy[k] - b * sx[k]) / calls[k]
      printf "%-26s %6d %9.1f %9.1f %11.3f %9.3f %12.2f  (fit: %.3f s per call + %.2f s per unit)\n", k, calls[k], W[k], C[k], U[k], W[k] / calls[k], b, icpt, b
    }
  }' "$2" "$1" | sort
echo "units: M read pairs (art_illumina), Mb of reads (pbsim)"
