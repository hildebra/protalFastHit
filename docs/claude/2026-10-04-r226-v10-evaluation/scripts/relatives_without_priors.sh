#!/usr/bin/env bash
# Recommendation 5 of the v10 evaluation: the relatives features in place of the distance features, without priors,
# on the v10 tables (local/v10), against the v10 default set; same trainer options as the build.
# Usage (WSL): relatives_without_priors.sh REPO OUTDIR
set -euo pipefail
repo=${1:-/mnt/c/Users/hildebra/Documents/locDev/protal}
out=${2:-$HOME/v10_relatives}
py=${PY:-$HOME/micromamba/envs/protal-db-build/bin/python}
v10=$repo/local/v10
mkdir -p "$out"
for t in "" _se _pb _ont; do
  leaves=512; [[ $t == _pb || $t == _ont ]] && leaves=128
  for set in normalized+adjacency+distance+depth+divergence+unfiltered normalized+adjacency+relatives+depth+divergence+unfiltered; do
    name=$( [[ $set == *relatives* ]] && echo relatives || echo distance )
    nice -n 10 "$py" "$repo/scripts/random_forest_cmdline.py" --truth-file "$v10/training/training_data$t.tsv" \
      --output-prefix "$out/${name}$t" --features "$set" --ntree 64 --maxnodes $leaves --seed 1 --threads 4 \
      --taxonomy "$v10/internal_taxonomy.dmp" --evaluation basic --depth-knobs \
      --test-file "$v10/test/training_data$t.tsv" > "$out/${name}$t.log" 2>&1
  done
done
echo done > "$out/DONE"
