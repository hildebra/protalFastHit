#!/bin/bash
# Queue 3: redo MRate2-off variants that ran with the ineffective --iqr-mult switch.
set -uo pipefail
cd ~/audit5/phylo
S=scripts
until grep -q "QUEUE2 DONE" queue2.log; do sleep 20; done
python3 $S/redo_variants.py runs/base20 "" filt_nomr2,filt_none
python3 $S/analyze.py runs/base20 --species Malpha,Tone --variants filt_nomr2,filt_none
python3 $S/redo_variants.py runs/base20 k04 filt_nomr2,filt_none
python3 $S/analyze.py runs/base20 --species Cferv --protal_subdir protal_k04 --tag k04 --variants filt_nomr2,filt_none
echo QUEUE3 DONE
