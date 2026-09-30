#!/bin/bash
# Queue 8 (replaces queues 4-7, which never started): uneven depth, low-depth replicates,
# low-depth caller settings, congeners over a 5x target, d5 MRate2-off redo.
set -uo pipefail
cd ~/audit5/phylo
S=scripts
DV=raw,filt,filt_keepsing,filt_nomr2,filt_nocov,filt_none
V=raw,filt,filt_keepsing,filt_nomr2,filt_none
bash $S/run_sim.sh uneven --scenario depth --tip_depths 2,100,3,60,5,40,8,25,1.5,80,12,4
python3 $S/analyze.py runs/uneven --species Malpha,Tone,Cferv --variants $DV,filt_zeros,filt_sabs2
bash $S/run_protal.sh uneven protal_k04 --knob 0.4
python3 $S/analyze.py runs/uneven --species Cferv --protal_subdir protal_k04 --tag k04 --variants $DV
for spec in "d2r2 2 101 8" "d2r3 2 202 9"; do
  set -- $spec
  SIMSEED=$3 bash $S/run_sim.sh $1 --scenario depth --depth $2 --seed $4
  python3 $S/analyze.py runs/$1 --species Malpha,Tone --variants $V
  bash $S/run_protal.sh $1 protal_k04 --knob 0.4
  python3 $S/analyze.py runs/$1 --species Cferv --protal_subdir protal_k04 --tag k04 --variants $V
done
bash $S/run_protal.sh d2 protal_mc1 --knob 0.4 --snp_min_cov 1 --snp_no_strand
python3 $S/analyze.py runs/d2 --protal_subdir protal_mc1 --tag mc1 --variants raw,filt,filt_keepsing
bash $S/run_sim.sh congener5 --scenario congener --depth 5
python3 $S/analyze.py runs/congener5 --species Malpha --variants raw,filt,filt_keepsing,filt_nomr2,filt_none
python3 $S/redo_variants.py runs/d5 "" filt_nomr2,filt_none
python3 $S/analyze.py runs/d5 --species Malpha,Tone --variants filt_nomr2,filt_none
for spec in "d3 3 303 10" "d3r2 3 404 11"; do
  set -- $spec
  SIMSEED=$3 bash $S/run_sim.sh $1 --scenario depth --depth $2 --seed $4
  python3 $S/analyze.py runs/$1 --species Malpha,Tone --variants $V
done
echo QUEUE8 DONE
