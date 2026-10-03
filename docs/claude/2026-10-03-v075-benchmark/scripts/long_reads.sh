#!/usr/bin/env bash
# 0.7.5 and 0.7.3 on the long-read samples of 0.7.3's collector ($B/samples_lr073, ../../2026-10-02-v073-benchmark/
# scripts/long_reads.sh: PacBio HiFi reads with qualities by length, Nanopore by pbsim3; 3 and 90 Mb, 4 samples
# each), against both databases of each. Runs into $B/runs_lr075, as profile.sh runs them.
set -uo pipefail
B=${BENCH:-$HOME/bench071}
T=${T:-6}
O=$B/samples_lr073
R=$B/runs_lr075
mkdir -p $R
[ -f $O/.simulated ] || { echo "no $O (the v0.7.3 benchmark's long_reads.sh makes it)"; exit 1; }
declare -A DB=([v073.full]=$B/V073/protal_db [v073.missing]=$B/V073/training_db
               [v075.full]=$B/V075/protal_db [v075.missing]=$B/V075/training_db)
binary() { echo $B/bin/protal-0.${1:2:1}.${1:3:1}; }
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
for reads in $O/points/{pb,ont}_b*/sim/reads/*.fq.gz; do
  id=$(basename $reads .fq.gz)
  type=${id%%_*}
  echo "$(date +%T) $id"
  for db in full missing; do
    args=(-1 $reads --read_type $type --prefix $id --no_qcmsa)
    for v in v073 v075; do
      run $v.$db.$type.$id $(binary $v) --db ${DB[$v.$db]} "${args[@]}"
    done
  done
done
echo "$(date +%T) done"
