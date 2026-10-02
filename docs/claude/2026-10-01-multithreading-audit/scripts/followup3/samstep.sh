#!/usr/bin/env bash
# Where "Writing the SAM header and file" goes: 88f20d4 with samstep.patch (SAMSTEP lines: ms per step), the
# 5M-pair sample (BGZF) against the full database at 6 threads, --no_profile, .sam.zst three times and .sam.gz
# twice, niced. SAMSTEP lines and the step's time into ~/mt-work/samstep/steps.tsv.
set -uo pipefail
SP=$(cd "$(dirname "$0")" && pwd)  # this folder
W=$HOME/mt-work/samstep; mkdir -p $W
bash $SP/../followup2/build_wt.sh samstep $SP/samstep.patch 88f20d4 protal > $W/build.out 2>&1 || { cat $W/build.out; exit 1; }
B=$HOME/bench071; R=$B/samples_deep/points/rl150_p5000000/sim/reads/rl150_p5000000_s_1
cat ${R}_R1.fq.gz ${R}_R2.fq.gz > /dev/null
: > $W/steps.tsv
for run in zst1 gz1 zst2 gz2 zst3; do
  fmt=${run%?}; d=$W/run; rm -rf $d; mkdir -p $d
  nice -n 5 $HOME/mt-work/samstep/src/build/protal --db $B/V071/protal_db -1 ${R}_R1.fq.gz -2 ${R}_R2.fq.gz --prefix s -o $d \
      -t 6 --no_profile --sam_format $fmt > $W/log 2>&1
  grep SAMSTEP $W/log | sed "s/^/$run\t/" | tee -a $W/steps.tsv
  echo "$run	$(grep 'Writing the SAM header' $W/log)	$(ls -la $d/s.sam.* | awk '{print $5}')" | tee -a $W/steps.tsv
done
echo SAMSTEP DONE
