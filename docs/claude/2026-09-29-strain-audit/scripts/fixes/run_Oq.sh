#!/bin/bash
# Step O's binary with a higher base-quality floor for alleles (a lone read's mismatch needs Q30):
# accuracy runs only, label Oq30.
set -u
S=$(cd "$(dirname "$0")" && pwd)
cp $HOME/audit5/bin/protal_O $HOME/audit5/bin/protal_Oq30
cp $HOME/audit5/bin/qcmsa_O.py $HOME/audit5/bin/qcmsa_Oq30.py
EXTRA="--snp_min_mean_qual 30" bash $S/eval_step.sh Oq30 acc
echo "ALL_Oq30_DONE"
