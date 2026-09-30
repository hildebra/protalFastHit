#!/bin/bash
# Trees of the low-depth runs with step O's binary and a Q30 allele floor (label Oq30).
set -u
S=$(cd "$(dirname "$0")" && pwd)
while pgrep -f "eval_step.sh P" > /dev/null; do sleep 30; done
RUNS="d2 d2r2 d2r3 d3 d3r2 d5 uneven" EXTRA="--snp_min_mean_qual 30" bash $S/eval_step.sh Oq30 phylo
echo "ALL_Oq30_PHYLO_DONE"
