#!/bin/bash
# Keep the worktree build as step P (with qcmsa's 0.002 floor) and evaluate it.
set -u
L=${1:-P}
S=$(cd "$(dirname "$0")" && pwd)
cp $HOME/strain-build/bin/protal $HOME/audit5/bin/protal_$L || exit 1
cp $HOME/strain-build/src/scripts/qcmsa.py $HOME/audit5/bin/qcmsa_$L.py || exit 1
echo "saved $L"
while pgrep -f "eval_step.sh [LMNO]" > /dev/null; do sleep 30; done
bash $S/eval_step.sh $L
echo "ALL_${L}_DONE"
