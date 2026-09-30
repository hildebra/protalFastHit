#!/bin/bash
# Scenario queue 1: baseline (rest), depth series, mixed, congener. Sequential (2 threads).
set -uo pipefail
cd ~/audit5/phylo
S=scripts
DV=raw,filt,filt_keepsing,filt_nomr2,filt_nocov,filt_none
python3 $S/analyze.py runs/base20 --species Tone,Cferv
bash $S/run_protal.sh base20 protal_k04 --knob 0.4
python3 $S/analyze.py runs/base20 --species Cferv --protal_subdir protal_k04 --tag k04
for d in 5 2 10 50; do
  bash $S/run_sim.sh d$d --scenario depth --depth $d
  python3 $S/analyze.py runs/d$d --species Malpha,Tone,Cferv --variants $DV
  bash $S/run_protal.sh d$d protal_k04 --knob 0.4
  python3 $S/analyze.py runs/d$d --species Cferv --protal_subdir protal_k04 --tag k04 --variants $DV
done
for d in 2 5; do
  bash $S/run_protal.sh d$d protal_h1 --knob 0.4 --msa_min_hcov 1 --no_qcmsa
  python3 $S/analyze.py runs/d$d --protal_subdir protal_h1 --tag h1 --reapply_hcov 1 --variants raw,filt,filt_keepsing,filt_none
done
bash $S/run_sim.sh mixed --scenario mixed --depth 20
python3 $S/analyze.py runs/mixed --mix mix=t03 --variants raw,filt,filt_keepsing,filt_nomr2,filt_none
bash $S/run_protal.sh mixed protal_k04 --knob 0.4
python3 $S/analyze.py runs/mixed --species Cferv --mix mix=t03 --protal_subdir protal_k04 --tag k04 --variants raw,filt,filt_keepsing,filt_nomr2,filt_none
bash $S/run_sim.sh congener --scenario congener --depth 20
python3 $S/analyze.py runs/congener --species Malpha --variants raw,filt,filt_keepsing,filt_nomr2,filt_none
echo QUEUE1 DONE
