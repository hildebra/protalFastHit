#!/bin/bash
# Keep the worktree build as step N, wait for the L/M evaluation to finish, then evaluate N.
# Usage: run_N.sh LABEL  (default N)
set -u
L=${1:-N}
S=$(cd "$(dirname "$0")" && pwd)
cp $HOME/strain-build/bin/protal $HOME/audit5/bin/protal_$L || exit 1
cp $HOME/strain-build/src/scripts/qcmsa.py $HOME/audit5/bin/qcmsa_$L.py || exit 1
echo "saved $L"
while pgrep -f "eval_step.sh [LM]" > /dev/null; do sleep 30; done
bash $S/eval_step.sh $L
echo "ALL_${L}_DONE"
