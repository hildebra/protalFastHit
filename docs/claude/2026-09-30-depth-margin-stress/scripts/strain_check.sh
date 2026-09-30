#!/usr/bin/env bash
# The strain MSAs under the new default: the strain audit's runs A (one strain per species, 1-50x) and B
# (two strains per species, minor 2-50%) aligned and profiled by 0.7 with --depth_identity_margin 0.04 and
# 0.08 (qcmsa on), and scored against the true genotypes by the benchmark's copy of the audit's evaluate.py.
set -uo pipefail
B=${BENCH:-$HOME/bench07}
SRC7=${SRC7:-$HOME/fix-build/src}
P7=${P7:-$HOME/fix-build/bin/protal}
ACC=${ACC:-$HOME/audit5/accuracy}
T=${T:-6}
HERE=$(cd $(dirname $0) && pwd)
EVAL=$HERE/../../2026-09-30-v07-vs-v06/scripts
AUDIT=$HERE/../../2026-09-29-strain-audit/scripts/accuracy
R=$B/runs
python3 $EVAL/make_colmap.py $SRC7/scripts/qcmsa.py $B/qcmsa_colmap_new.py
for run in A B; do
  for m in 0.04 0.08; do
    name=strain$run.m$m
    if [ ! -f $R/$name.done ]; then
      rm -rf $R/$name; mkdir -p $R/$name
      echo "$(date +%T) $name"
      /usr/bin/time -v -o $R/$name.time $P7 --db $B/strainw/db07 --map $ACC/sim_$run/protal.meta -o $R/$name -t $T \
        --qcmsa_script $SRC7/scripts/qcmsa.py --depth_identity_margin $m > $R/$name.log 2>&1 && touch $R/$name.done ||
        echo "  $name failed"
    fi
    QCMSA_COLMAP=$B/qcmsa_colmap_new.py python3 $EVAL/evaluate.py $run $R/$name margin_${run}_$m $R/$name.log > /dev/null
  done
done
python3 $AUDIT/summarize.py margin_A_0.04 margin_A_0.08 margin_B_0.04 margin_B_0.08 > $B/strain_margin_summary.txt
echo "summary: $B/strain_margin_summary.txt"
