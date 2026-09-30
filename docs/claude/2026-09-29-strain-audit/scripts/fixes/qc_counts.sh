#!/bin/bash
# qcmsa removals per label on the accuracy runs: samples, MRate2 genes and outlier cells.
A=$HOME/audit5/accuracy
for L in "$@"; do
  for r in A B Cs Cl; do
    d=$A/prot_${r}_$L/strains
    [ -d $d ] || continue
    s=$(cat $d/*.qcmsa_summary.tsv | awk -F'\t' '$2=="samples_filtered"{s+=$3} $2=="genes_filtered_mrate2"{g+=$3} $2=="outlier_cells"{c+=$3} $2=="samples_in"{n+=$3} END{printf "samples %d/%d, mrate2 genes %d, outlier cells %d", s, n, g, c}')
    echo "$L $r: $s"
  done
done
