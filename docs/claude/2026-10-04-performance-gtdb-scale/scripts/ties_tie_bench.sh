#!/usr/bin/env bash
# Seeds tied with the last one taken: dropped (before, 186c8ed) or kept for their own diagonal's anchor (after), on the
# 0.7.5 benchmark's samples and databases (~/bench071, docs/claude/2026-10-03-v075-benchmark): paired-end rl100/rl150 at
# 1k/10k/500k pairs and PacBio/Nanopore at 3/90 Mb, 4 samples each, full and missing databases, -t 6. Runs as the
# benchmark's, in a folder of their own ($F) where "v073" is before and "v075" after, so that its score.py pairs them.
set -uo pipefail
B=${BENCH:-$HOME/bench071}; F=${TIES:-$HOME/tiebench}; T=${T:-6}
declare -A BIN=([v073]=$F/bin/protal-before [v075]=$F/bin/protal-after)
declare -A DB=([full]=$B/V075/protal_db [missing]=$B/V075/training_db)
mkdir -p $F/V073 $F/runs_v075 $F/runs_lr075
ln -sfn $B/V075 $F/V075; for d in world samples samples_deep samples_lr073; do ln -sfn $B/$d $F/$d; done
ln -sf $B/V075/heldout_species.txt $F/V073/heldout_species.txt  # both are 0.7.5's databases: its held-out species
run() {
  local R=$1 name=$2 keep=$3; shift 3
  [ -f $R/$name.done ] && return
  rm -rf $R/$name; mkdir -p $R/$name
  cat /proc/loadavg > $R/$name.load
  if /usr/bin/time -v -o $R/$name.time "$@" -o $R/$name -t $T > $R/$name.log 2>&1; then touch $R/$name.done
  else echo "  $name failed ($?), see $R/$name.log"; fi
  [ $keep = 1 ] || rm -f $R/$name/*.sam.zst
}
for reads in $B/samples/points/rl*_p{1000,10000,500000}/sim/reads/*_R1.fq.gz; do
  id=$(basename $reads _R1.fq.gz); r2=${reads%_R1.fq.gz}_R2.fq.gz
  for db in full missing; do for v in v073 v075; do
    keep=0; case $id in *_p500000_*) [ $db = full ] && keep=1;; esac
    run $F/runs_v075 $v.$db.pe.$id $keep ${BIN[$v]} --db ${DB[$db]} -1 $reads -2 $r2 --prefix $id --no_qcmsa
  done; done
done
echo "$(date +%T) paired-end done"
for reads in $B/samples_lr073/points/{pb,ont}_b*/sim/reads/*.fq.gz; do
  id=$(basename $reads .fq.gz); type=${id%%_*}
  for db in full missing; do for v in v073 v075; do
    run $F/runs_lr075 $v.$db.$type.$id 0 ${BIN[$v]} --db ${DB[$db]} -1 $reads --read_type $type --prefix $id --no_qcmsa
  done; done
done
echo "$(date +%T) long reads done"
