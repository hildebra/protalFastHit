#!/usr/bin/env bash
# 0.7.0 with 0.7.1's depth margin: 0.7.0's SAMs of the paired-end samples profiled again (--profile_only) with
# --depth_identity_margin 0.08, against both databases (runs/v070m008.<db>.pe.<sample>), to tell how much of
# 0.7.1's abundance gain is the margin.
set -uo pipefail
B=${BENCH:-$HOME/bench071}
R=$B/runs
for done in $R/v070.*.pe.*.done; do
  run=$(basename $done .done)
  db=$(echo $run | cut -d. -f2) sample=${run#v070.$db.pe.}
  name=v070m008.$db.pe.$sample
  [ -f $R/$name.done ] && continue
  dbdir=$B/V070/$([ $db = full ] && echo protal_db || echo training_db)
  sam=$(ls $R/$run/*.sam.zst 2>/dev/null | head -1)
  [ -n "$sam" ] || { echo "no SAM in $R/$run"; continue; }
  mkdir -p $R/$name
  cat /proc/loadavg > $R/$name.load
  /usr/bin/time -v -o $R/$name.time $B/bin/protal-0.7.0 --db $dbdir --profile_only $sam --prefix $sample -o $R/$name \
    -t 6 --no_qcmsa --depth_identity_margin 0.08 > $R/$name.log 2>&1 && touch $R/$name.done || echo "$name failed"
done
