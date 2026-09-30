#!/bin/bash
# Build step L (d403dac), take step M from the worktree build, and evaluate both.
set -u
S=$(cd "$(dirname "$0")" && pwd)
bash $S/build_commit.sh d403dac L || exit 1
cp $HOME/strain-build/bin/protal $HOME/audit5/bin/protal_M
cp $HOME/strain-build/src/scripts/qcmsa.py $HOME/audit5/bin/qcmsa_M.py
bash $S/eval_step.sh L
bash $S/eval_step.sh M
echo ALL_LM_DONE
