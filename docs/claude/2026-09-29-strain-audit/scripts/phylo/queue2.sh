#!/bin/bash
# Scenario queue 2: reference row in/out, SNP-caller settings (IUPAC noise), after queue1.
set -uo pipefail
cd ~/audit5/phylo
S=scripts
until grep -q "QUEUE1 DONE" queue1.log; do sleep 20; done
python3 $S/analyze.py runs/base20 --species Malpha,Tone --variants raw_noref,filt_noref
python3 $S/analyze.py runs/base20 --species Cferv --protal_subdir protal_k04 --tag k04 --variants raw_noref,filt_noref
for r in base20 d50; do
  bash $S/run_protal.sh $r protal_af10 --knob 0.4 --snp_min_af 0.1
  python3 $S/analyze.py runs/$r --protal_subdir protal_af10 --tag af10 --variants raw,filt,filt_keepsing
  bash $S/run_protal.sh $r protal_ma1 --knob 0.4 --snp_max_alleles 1
  python3 $S/analyze.py runs/$r --protal_subdir protal_ma1 --tag ma1 --variants raw,filt,filt_keepsing
done
python3 $S/analyze.py runs/mixed --species Malpha,Tone --mix mix=t03 --variants filt_strict,filt_sabs2,filt_zeros
echo QUEUE2 DONE
