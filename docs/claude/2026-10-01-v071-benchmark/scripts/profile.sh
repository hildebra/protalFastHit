#!/usr/bin/env bash
# Every version on the benchmark samples (samples.sh), one timed run per sample and database
# (/usr/bin/time -v in runs/<run>.time, the load average before it in runs/<run>.load), -t 6, no qcmsa:
#   v060.<db>.pe.<sample>      0.6.0a, its own databases (db060_full, db060_missing), the model it ships
#   v070.<db>.<reads>.<sample> 0.7.0, the databases and the four models its pipeline built (V070)
#   v071.<db>.<reads>.<sample> 0.7.1, likewise (V071: with gene_neighbours.tsv; the margin 0.08)
#   v071m070.<db>.pe.<sample>  0.7.1 with 0.7.0's paired-end model: the binary's changes apart from the model's
# <db>: full (765 species) or missing (the training database: held-out species and clades left out; the four
# trained models stored in it first with --add_model). <reads>: pe, se (the first reads), pb, ont; 0.6.0a
# profiles paired-end reads only. A run that completed is not run again (runs/<run>.done).
set -uo pipefail
B=${BENCH:-$HOME/bench071}
T=${T:-6}
R=$B/runs
mkdir -p $R
declare -A DB=([v060.full]=$B/db060_full [v060.missing]=$B/db060_missing
               [v070.full]=$B/V070/protal_db [v070.missing]=$B/V070/training_db
               [v071.full]=$B/V071/protal_db [v071.missing]=$B/V071/training_db)

binary() { echo $B/bin/protal-0.${1:2:1}.${1:3:1}; }  # v071 -> protal-0.7.1
binary060=$B/bin/protal-0.6.0a

# The trained models into each training database (the finished database has them from its pipeline).
for v in 070 071; do
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
}

for reads in $B/samples/points/rl*/sim/reads/*_R1.fq.gz $B/samples_deep/points/rl*/sim/reads/*_R1.fq.gz; do
  id=$(basename $reads _R1.fq.gz)
  r2=${reads%_R1.fq.gz}_R2.fq.gz
  echo "$(date +%T) $id"
  for db in full missing; do
    pair=(-1 $reads -2 $r2 --prefix $id --no_qcmsa)
    run v060.$db.pe.$id $binary060 --db ${DB[v060.$db]} "${pair[@]}"
    run v070.$db.pe.$id $(binary v070) --db ${DB[v070.$db]} "${pair[@]}"
    run v071.$db.pe.$id $(binary v071) --db ${DB[v071.$db]} "${pair[@]}"
    run v071m070.$db.pe.$id $(binary v071) --db ${DB[v071.$db]} --model $B/V070/trained_model.xml "${pair[@]}"
    case $id in rl*_p5000000_*) continue;; esac  # the deep samples: paired-end only
    for v in v070 v071; do
      run $v.$db.se.$id $(binary $v) --db ${DB[$v.$db]} -1 $reads --read_type se --prefix $id --no_qcmsa
    done
  done
done
for reads in $B/samples/points/{pb,ont}_b*/sim/reads/*.fq.gz; do
  id=$(basename $reads .fq.gz)
  type=${id%%_*}
  echo "$(date +%T) $id"
  for db in full missing; do
    for v in v070 v071; do
      run $v.$db.$type.$id $(binary $v) --db ${DB[$v.$db]} -1 $reads --read_type $type --prefix $id --no_qcmsa
    done
  done
done
echo "$(date +%T) done"
