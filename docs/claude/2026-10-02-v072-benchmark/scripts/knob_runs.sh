#!/usr/bin/env bash
# 0.7.2's long-read runs again with --knob 0.5 (v072k.*, into $B/runs_v072 beside profile.sh's): its models' depth
# knobs off, to tell their part in 0.7.2's long-read results from the rest.
set -uo pipefail
B=${BENCH:-$HOME/bench071}
T=${T:-6}
R=$B/runs_v072
declare -A DB=([full]=$B/V072/protal_db [missing]=$B/V072/training_db)
run() {
  local name=$1; shift
  [ -f $R/$name.done ] && return
  rm -rf $R/$name
  mkdir -p $R/$name
  cat /proc/loadavg > $R/$name.load
  if /usr/bin/time -v -o $R/$name.time "$@" -o $R/$name -t $T > $R/$name.log 2>&1; then
    touch $R/$name.done
  else
    echo "  $name failed ($?), see $R/$name.log"
  fi
  rm -f $R/$name/*.sam.zst $R/$name/*.sam.gz $R/$name/*.sam
}
for reads in $B/samples/points/{pb,ont}_b*/sim/reads/*.fq.gz; do
  id=$(basename $reads .fq.gz)
  type=${id%%_*}
  for db in full missing; do
    run v072k.$db.$type.$id $B/bin/protal-0.7.2 --db ${DB[$db]} -1 $reads --read_type $type --prefix $id --no_qcmsa --knob 0.5
  done
done
echo "$(date +%T) done"
