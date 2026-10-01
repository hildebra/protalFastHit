#!/usr/bin/env bash
# The four model features and the depth knobs on the v0.7.1 benchmark world (docs/claude/2026-10-01-v071-benchmark,
# data in $B = ~/bench071):
#   1. protal at $REV (git archive, Release) as $B/bin/protal-$VERSION; its scripts in $B/src/$VERSION
#   2. its own build_gtdb_database.py pipeline into $B/V<VERSION without dots>: the 0.7.1 pipeline's design and
#      seed (../../2026-10-01-v071-benchmark/scripts/pipelines.sh), 0.7.1's simulate_metagenomes, so the same samples
#   3. the benchmark's samples against both databases (full, missing), as profile.sh ran 0.7.1:
#        $TAG.<db>.<reads>.<sample>    by default: the pb and ont models' depth knobs apply
#        ${TAG}k.<db>.<pb|ont>.<sample>  with --knob 0.5: the features without the depth knobs
# VERSION=0.7.2dev TAG=v072 REV=25d457e by default (the four features and the depth knobs); VERSION=0.7.3dev TAG=v073
# REV=<commit> for the next one. Runs that completed are not run again (runs/<run>.done).
set -uo pipefail
B=${BENCH:-$HOME/bench071}
REPO=${REPO:-/mnt/c/Users/hildebra/Documents/locDev/protal}
REV=${REV:-25d457e}
T=${T:-6}
v=${VERSION:-0.7.2dev}
tag=${TAG:-v072}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p $B/bin $B/src $B/logs $B/runs

if [ ! -x $B/bin/protal-$v ]; then
  S=$B/src/$v
  rm -rf $S; mkdir -p $S
  git -C $REPO archive $REV | tar -x -C $S || exit 1
  (cd $S && cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release > $B/logs/configure_$v.log 2>&1 &&
    cmake --build build --target protal -j 4 > $B/logs/build_$v.log 2>&1) || { echo "build of $v failed"; exit 1; }
  cp $S/build/protal $B/bin/protal-$v
  echo "$(date +%T) built $v ($REV)"
fi

VERSIONS=$v bash $HERE/../../2026-10-01-v071-benchmark/scripts/pipelines.sh || exit 1
P=$B/V${v//./}
cmp -s $B/V071/heldout_species.txt $P/heldout_species.txt || { echo "$P held out other species than V071"; exit 1; }

marker=$P/training_db/.models_added
if [ ! -f $marker ]; then
  for type in pe se pb ont; do
    model=$P/trained_model$([ $type = pe ] || echo _$type).xml
    $B/bin/protal-$v --add_model $model --read_type $type --db $P/training_db > $B/logs/add_model_$tag.$type.log 2>&1 ||
      { echo "--add_model $model failed"; exit 1; }
  done
  touch $marker
fi

R=$B/runs
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
}
declare -A DB=([full]=$P/protal_db [missing]=$P/training_db)
bin=$B/bin/protal-$v
for reads in $B/samples/points/rl*/sim/reads/*_R1.fq.gz $B/samples_deep/points/rl*/sim/reads/*_R1.fq.gz; do
  id=$(basename $reads _R1.fq.gz)
  r2=${reads%_R1.fq.gz}_R2.fq.gz
  echo "$(date +%T) $id"
  for db in full missing; do
    run $tag.$db.pe.$id $bin --db ${DB[$db]} -1 $reads -2 $r2 --prefix $id --no_qcmsa
    case $id in rl*_p5000000_*) continue;; esac  # the deep samples: paired-end only
    run $tag.$db.se.$id $bin --db ${DB[$db]} -1 $reads --read_type se --prefix $id --no_qcmsa
  done
done
for reads in $B/samples/points/{pb,ont}_b*/sim/reads/*.fq.gz; do
  id=$(basename $reads .fq.gz)
  type=${id%%_*}
  echo "$(date +%T) $id"
  for db in full missing; do
    run $tag.$db.$type.$id $bin --db ${DB[$db]} -1 $reads --read_type $type --prefix $id --no_qcmsa
    run ${tag}k.$db.$type.$id $bin --db ${DB[$db]} -1 $reads --read_type $type --prefix $id --no_qcmsa --knob 0.5
  done
done
echo "$(date +%T) done"
