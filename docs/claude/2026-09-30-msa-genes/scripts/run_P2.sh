#!/bin/bash
# The worktree build with the 90% long-unique threshold (label P2): accuracy runs only.
set -u
S=$(cd "$(dirname "$0")" && pwd)
cp $HOME/strain-build/bin/protal $HOME/audit5/bin/protal_P2 || exit 1
cp $HOME/strain-build/src/scripts/qcmsa.py $HOME/audit5/bin/qcmsa_P2.py || exit 1
bash $S/eval_step.sh P2 acc
bash $S/gene_counts.sh P P2
cd $HOME/audit5/accuracy
python3 - <<'EOF'
import sys, os
sys.path.insert(0, os.path.expanduser("~/audit5/accuracy"))
import summarize as S
for lab in ("A_P", "A_P2"):
    df = S.load(lab)
    for stage in ("raw", "qc"):
        d = df[(df.stage == stage) & (df.kind == "strain") & (df.species == "Dummya solo")]
        t = S.table(d, ["depth"])
        print(lab, stage)
        print(S.fmt(t[["depth", "called", "acc", "sensCov", "N_snp", "fSNPppm", "iupac"]]))
EOF
echo "ALL_P2_DONE"
