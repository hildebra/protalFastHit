#!/bin/bash
# Profile differences between two labels on the accuracy runs. Usage: profile_diff.sh L M
S=$(cd "$(dirname "$0")" && pwd)
A=$HOME/audit5/accuracy
for run in A B Cs Cl; do
  [ -d $A/prot_${run}_$1 ] && [ -d $A/prot_${run}_$2 ] || continue
  echo "== run $run: $1 -> $2"
  python3 $S/profile_diff.py $A/prot_${run}_$1 $A/prot_${run}_$2
done
