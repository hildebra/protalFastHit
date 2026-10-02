#!/usr/bin/env bash
# 0.6.0a, 0.7.0, 0.7.1 and 0.7.2 on the v0.7.1 benchmark's samples (docs/claude/2026-10-01-v071-benchmark, data in
# $B = ~/bench071), all four run again now, one sample after the other with every version in turn, so that the
# machine's load hits them alike:
#   1. protal 0.7.2 (Version 0.7.2, $REV) from git (git archive, Release) as $B/bin/protal-0.7.2, its scripts in
#      $B/src/0.7.2; the other versions' binaries and databases from the v0.7.1 benchmark
#   2. 0.7.2's own build_gtdb_database.py pipeline into $B/V072 (the v0.7.1 benchmark's pipelines.sh: the design and
#      seed of 0.7.0's and 0.7.1's pipelines, 0.7.1's simulate_metagenomes, so the same training samples)
#   3. runs into $B/runs_v072 (time -v, the load average before each), -t 6, no qcmsa:
#        v060.<db>.pe.<sample>        0.6.0a, its own databases (db060_full, db060_missing), the model it ships
#        v07N.<db>.<reads>.<sample>   0.7.0, 0.7.1, 0.7.2, the databases and four models their pipelines built
#      <db>: full (765 species) or missing (the training database: 290 species left out; its pipeline's trained
#      models added first with --add_model). Paired-end for all four; single-end, PacBio and ONT for the 0.7 versions.
# A run that completed is not run again (<run>.done).
set -uo pipefail
B=${BENCH:-$HOME/bench071}
REPO=${REPO:-/mnt/c/Users/hildebra/Documents/locDev/protal}
REV=${REV:-fd1e719}
T=${T:-6}
R=$B/runs_v072
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p $B/bin $B/src $B/logs $R

if [ ! -x $B/bin/protal-0.7.2 ]; then
  S=$B/src/0.7.2
  rm -rf $S; mkdir -p $S
  git -C $REPO archive $REV | tar -x -C $S || exit 1
  (cd $S && cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release > $B/logs/configure_0.7.2.log 2>&1 &&
    cmake --build build --target protal -j 4 > $B/logs/build_0.7.2.log 2>&1) || { echo "build of 0.7.2 failed"; exit 1; }
  cp $S/build/protal $B/bin/protal-0.7.2
  echo "$(date +%T) built 0.7.2 ($REV): $($B/bin/protal-0.7.2 --version 2>&1 | tail -1)"
fi

VERSIONS=0.7.2 bash $HERE/../../2026-10-01-v071-benchmark/scripts/pipelines.sh || exit 1
cmp -s $B/V071/heldout_species.txt $B/V072/heldout_species.txt || { echo "V072 held out other species than V071"; exit 1; }

declare -A DB=([v060.full]=$B/db060_full [v060.missing]=$B/db060_missing
               [v070.full]=$B/V070/protal_db [v070.missing]=$B/V070/training_db
               [v071.full]=$B/V071/protal_db [v071.missing]=$B/V071/training_db
               [v072.full]=$B/V072/protal_db [v072.missing]=$B/V072/training_db)
binary() { echo $B/bin/protal-0.${1:2:1}.${1:3:1}; }  # v072 -> protal-0.7.2
binary060=$B/bin/protal-0.6.0a

# The trained models into each training database (the finished database has them from its pipeline; 0.7.0's and
# 0.7.1's were added by the v0.7.1 benchmark).
for v in 070 071 072; do
  marker=$B/V$v/training_db/.models_added
  [ -f $marker ] && continue
  for type in pe se pb ont; do
    model=$B/V$v/trained_model$([ $type = pe ] || echo _$type).xml
    $(binary v$v) --add_model $model --read_type $type --db $B/V$v/training_db > $B/logs/add_model_$v.$type.log 2>&1 ||
      { echo "--add_model $model failed"; exit 1; }
  done
  touch $marker
done

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
  rm -f $R/$name/*.sam.zst $R/$name/*.sam.gz $R/$name/*.sam  # the alignments are not scored; they take the space
}

for reads in $B/samples/points/rl*/sim/reads/*_R1.fq.gz $B/samples_deep/points/rl*/sim/reads/*_R1.fq.gz; do
  id=$(basename $reads _R1.fq.gz)
  r2=${reads%_R1.fq.gz}_R2.fq.gz
  echo "$(date +%T) $id"
  for db in full missing; do
    pair=(-1 $reads -2 $r2 --prefix $id --no_qcmsa)
    run v060.$db.pe.$id $binary060 --db ${DB[v060.$db]} "${pair[@]}"
    for v in v070 v071 v072; do
      run $v.$db.pe.$id $(binary $v) --db ${DB[$v.$db]} "${pair[@]}"
    done
    case $id in rl*_p5000000_*) continue;; esac  # the deep samples: paired-end only
    for v in v070 v071 v072; do
      run $v.$db.se.$id $(binary $v) --db ${DB[$v.$db]} -1 $reads --read_type se --prefix $id --no_qcmsa
    done
  done
done
for reads in $B/samples/points/{pb,ont}_b*/sim/reads/*.fq.gz; do
  id=$(basename $reads .fq.gz)
  type=${id%%_*}
  echo "$(date +%T) $id"
  for db in full missing; do
    for v in v070 v071 v072; do
      run $v.$db.$type.$id $(binary $v) --db ${DB[$v.$db]} -1 $reads --read_type $type --prefix $id --no_qcmsa
    done
  done
done
echo "$(date +%T) done"
