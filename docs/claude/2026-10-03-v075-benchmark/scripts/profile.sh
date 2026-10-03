#!/usr/bin/env bash
# 0.7.5 against 0.7.3 and 0.6.0a on the v0.7.1 benchmark's samples (docs/claude/2026-10-01-v071-benchmark, data in
# $B = ~/bench071), run now one sample after the other with every version in turn, so that the machine's load hits
# them alike. The 0.7.0-0.7.2 databases are gone from $B (their numbers: the v0.7.3 benchmark's tables, which that
# benchmark found reproducible to the fourth digit); 0.6.0a's databases and 0.7.3's pipeline (V073) are still there.
#   1. protal 0.7.5 ($REV, the checkout's HEAD) from git (git archive, Release) as $B/bin/protal-0.7.5, its scripts
#      in $B/src/0.7.5
#   2. 0.7.5's own build_gtdb_database.py pipeline into $B/V075: the design and seed of the earlier pipelines
#      (pipelines.sh of the v0.7.1 benchmark; --long-read-samples 4) and otherwise 0.7.5's defaults: the depth and
#      divergence features (no knob curve), suspect gene copies, species priors, the mixed training design (sigma
#      1.3 and 2.0, congeners), HiFi PacBio reads, 0.7.5's own simulate_metagenomes (0.7.1's, which the earlier
#      pipelines used for the same training communities, lacks what 0.7.5's collector asks of it), so the training
#      communities differ from V073's; 349 species held out instead of 290 (0.7.5's holdout design), so the
#      "missing" databases differ too. --no-binary-check: the binaries from git archive carry no commit.
#   3. runs into $B/runs_v075 (time -v, the load average before each), -t 6, no qcmsa:
#        v060.<db>.pe.<sample>        0.6.0a, its own databases (db060_full, db060_missing), the model it ships
#        v073.<db>.<reads>.<sample>   0.7.3, the databases and four models its pipeline built (V073)
#        v075.<db>.<reads>.<sample>   0.7.5, the databases and four models its pipeline built (V075)
#      <db>: full (765 species) or missing (the training database: 290 species left out; the pipeline's trained
#      models added first with --add_model). Paired-end for all; single-end for the 0.7 versions.
# A run that completed is not run again (<run>.done).
set -uo pipefail
B=${BENCH:-$HOME/bench071}
REPO=${REPO:-/mnt/c/Users/hildebra/Documents/locDev/protal}
REV=${REV:-$(git -C $REPO rev-parse --short HEAD)}
T=${T:-6}
W=$B/world
PY=${PY:-$HOME/micromamba/envs/protal-db-build/bin/python}
PBSIM=${PBSIM:-$HOME/micromamba/envs/protal-db-build/bin/pbsim}
R=$B/runs_v075
mkdir -p $B/bin $B/src $B/logs $R

if [ ! -x $B/bin/protal-0.7.5 ] || [ ! -x $B/bin/simulate_metagenomes-0.7.5 ]; then
  S=$B/src/0.7.5
  rm -rf $S; mkdir -p $S
  git -C $REPO archive $REV | tar -x -C $S || exit 1
  (cd $S && cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release > $B/logs/configure_0.7.5.log 2>&1 &&
    cmake --build build --target protal simulate_metagenomes -j 6 > $B/logs/build_0.7.5.log 2>&1) || { echo "build of 0.7.5 failed"; exit 1; }
  cp $S/build/protal $B/bin/protal-0.7.5
  cp $S/build/simulate_metagenomes $B/bin/simulate_metagenomes-0.7.5
  echo "$(date +%T) built 0.7.5 ($REV): $($B/bin/protal-0.7.5 --version 2>&1 | tail -1)"
fi

out=$B/V075
if [ ! -f $out/protal_db/database.protal ] || [ ! -s $out/model_logs/summary.txt ]; then
  echo "$(date +%T) pipeline of 0.7.5 into $out"
  rm -rf $out  # a failed attempt is not resumed
  /usr/bin/time -v -o $B/logs/pipeline_0.7.5.time $PY $B/src/0.7.5/scripts/build_gtdb_database.py --gtdb $W/release_p \
    --outdir $out --protal $B/bin/protal-0.7.5 --simulator $B/bin/simulate_metagenomes-0.7.5 -t $T --seed 1 \
    --extra-genomes $W/full/simulation/genomes_nonreps --samples 4 --read-pairs 1000,20000,200000 \
    --test-samples 2 --test-read-pairs 500,10000,500000 --long-read-bases 300000,6000000,60000000 \
    --test-long-read-bases 150000,3000000,90000000 --long-read-samples 4 --pbsim $PBSIM --no-binary-check \
    > $B/logs/pipeline_0.7.5.log 2>&1 ||
    { echo "pipeline of 0.7.5 failed, see $B/logs/pipeline_0.7.5.log"; tail -20 $B/logs/pipeline_0.7.5.log; exit 1; }
  grep -E "Elapsed|Maximum resident" $B/logs/pipeline_0.7.5.time
  cat $out/model_logs/summary.txt
fi
cmp -s $B/V073/heldout_species.txt $B/V075/heldout_species.txt || echo "note: V075 holds out other species than V073 (0.7.5's holdout design): the missing databases differ"

declare -A DB=([v060.full]=$B/db060_full [v060.missing]=$B/db060_missing
               [v073.full]=$B/V073/protal_db [v073.missing]=$B/V073/training_db
               [v075.full]=$B/V075/protal_db [v075.missing]=$B/V075/training_db)
binary() { echo $B/bin/protal-0.${1:2:1}.${1:3:1}; }  # v075 -> protal-0.7.5
binary060=$B/bin/protal-0.6.0a

# 0.7.5's trained models into its training database (the finished database has them from its pipeline).
marker=$B/V075/training_db/.models_added
if [ ! -f $marker ]; then
  for type in pe se pb ont; do
    model=$B/V075/trained_model$([ $type = pe ] || echo _$type).xml
    $B/bin/protal-0.7.5 --add_model $model --read_type $type --db $B/V075/training_db > $B/logs/add_model_075.$type.log 2>&1 ||
      { echo "--add_model $model failed"; exit 1; }
  done
  touch $marker
fi

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
    for v in v073 v075; do
      run $v.$db.pe.$id $(binary $v) --db ${DB[$v.$db]} "${pair[@]}"
    done
    case $id in rl*_p5000000_*) continue;; esac  # the deep samples: paired-end only
    single=(-1 $reads --read_type se --prefix $id --no_qcmsa)
    for v in v073 v075; do
      run $v.$db.se.$id $(binary $v) --db ${DB[$v.$db]} "${single[@]}"
    done
  done
done
echo "$(date +%T) done"
