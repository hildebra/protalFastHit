#!/bin/bash
# Genes in each species' raw MSA partition, per label (run A).
A=$HOME/audit5/accuracy
for L in "$@"; do
  echo "== $L"
  for p in $A/prot_A_$L/strains/*.raw.partition.txt; do
    echo "  $(basename $p .raw.partition.txt) $(wc -l < $p)"
  done
done
