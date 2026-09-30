#!/bin/bash
# Re-run evaluate.py for existing protal outputs of a label (no reprofiling), with that label's
# instrumented qcmsa. Waits until no other accuracy evaluation runs (they share qcmsa_colmap.py).
# Usage: reeval_acc.sh LABEL RUN...
set -u
L=$1; shift
A=$HOME/audit5/accuracy
while pgrep -f "eval_step.sh [A-Za-z0-9]+ acc|evaluate.py" > /dev/null; do sleep 30; done
sed "s|~/audit5/src/scripts/qcmsa.py|$HOME/audit5/bin/qcmsa_$L.py|" $A/make_qcmsa_colmap.py > $A/make_qcmsa_colmap_$L.py
python3 $A/make_qcmsa_colmap_$L.py
cd $A
for run in "$@"; do
  man=$run; [ $run = Cs ] || [ $run = Cl ] && man=C
  python3 $A/evaluate.py $man $A/prot_${run}_$L ${run}_$L $A/logs/protal_${run}_$L.log > $A/logs/eval_${run}_$L.log 2>&1 \
    && echo "evaluate $run ok" || { echo "evaluate $run failed"; tail -3 $A/logs/eval_${run}_$L.log; }
done
echo "REEVAL_${L}_DONE"
