#!/bin/bash
# Genotype accuracy per depth for the given labels (run A: bacteria except Dummya; raw and qcmsa).
# Usage: acc_cmp.sh RUN LABEL...   e.g. acc_cmp.sh A L M N   (labels become A_L, A_M, ...; "base" = A)
cd $HOME/audit5/accuracy
run=$1; shift
labels=""
for l in "$@"; do [ "$l" = base ] && labels="$labels $run" || labels="$labels ${run}_$l"; done
python3 compare.py --bacteria --nodummya $labels
python3 compare.py --bacteria --nodummya --stage qc $labels
