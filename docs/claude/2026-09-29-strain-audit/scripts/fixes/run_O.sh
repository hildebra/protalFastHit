#!/bin/bash
# Keep the worktree build as step O, wait for earlier evaluations, evaluate O, then add filtered
# trees with qcmsa --gene-min-mean-depth 1 (variant filt_depth1) to the O results.
set -u
L=${1:-O}
S=$(cd "$(dirname "$0")" && pwd)
cp $HOME/strain-build/bin/protal $HOME/audit5/bin/protal_$L || exit 1
cp $HOME/strain-build/src/scripts/qcmsa.py $HOME/audit5/bin/qcmsa_$L.py || exit 1
echo "saved $L"
while pgrep -f "eval_step.sh [LMN]" > /dev/null; do sleep 30; done
bash $S/eval_step.sh $L

P=$HOME/audit5/phylo
S2=$P/scripts_$L
grep -q '"filt_depth1"' $S2/analyze.py || \
  sed -i 's|^    "filt": (\[\], "GTR+G", False),|&\n    "filt_depth1": (["--gene-min-mean-depth", "1"], "GTR+G", False),|' $S2/analyze.py
one() {
  r=$1; D=$P/runs/$r; mix=""; sp="Malpha,Tone"
  [ $r = mixed ] && mix="--mix mix=t03"
  { [ $r = congener ] || [ $r = congener5 ]; } && sp="Malpha"
  python3 $S2/analyze.py $D --species $sp $mix --protal_subdir protal_$L --tag $L --variants filt_depth1 >> $D/eval_$L.log 2>&1
  echo "depth1 $r done"
}
export -f one; export P S2 L
printf '%s\n' base20 d2 d2r2 d2r3 d3 d3r2 d5 d10 d50 uneven mixed congener congener5 | xargs -P 3 -I{} bash -c 'one {}'
echo "ALL_${L}_DONE"
