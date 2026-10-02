#!/usr/bin/env bash
# 0.6.0a, 0.7.0, 0.7.1, 0.7.2 and 0.7.3 on the v0.7.1 benchmark's samples (docs/claude/2026-10-01-v071-benchmark, data
# in $B = ~/bench071), all five run again now, one sample after the other with every version in turn, so that the
# machine's load hits them alike:
#   1. protal 0.7.3 (Version 0.7.3, $REV) from git (git archive, Release) as $B/bin/protal-0.7.3, its scripts in
#      $B/src/0.7.3; the other versions' binaries and databases from the v0.7.1 and v0.7.2 benchmarks
#   2. 0.7.3's own build_gtdb_database.py pipeline into $B/V073: the design and seed of the other versions' pipelines
#      (docs/claude/2026-10-01-v071-benchmark/scripts/pipelines.sh) and 0.7.1's simulate_metagenomes, so the same
#      training communities; 0.7.3's other defaults (knob curves over depth for all four read types, HiFi PacBio
#      reads, 256 leaves per tree). --long-read-samples 4, the number of long-read samples per point the earlier
#      pipelines made (0.7.3's default is 24); --no-binary-check, as 0.7.1's simulator has no --version.
#   3. runs into $B/runs_v073 (time -v, the load average before each), -t 6, no qcmsa:
#        v060.<db>.pe.<sample>        0.6.0a, its own databases (db060_full, db060_missing), the model it ships
#        v07N.<db>.<reads>.<sample>   0.7.0 to 0.7.3, the databases and four models their pipelines built
#        v073k.<db>.<reads>.<sample>  0.7.3 at --knob 0.5: its knob curves over depth off
#      <db>: full (765 species) or missing (the training database: 290 species left out; its pipeline's trained
#      models added first with --add_model). Paired-end for all; single-end for the 0.7 versions. The long reads in
#      long_reads.sh (the benchmark's old long-read samples are half reads under 1 kb, ../2026-10-02-v072-benchmark).
# A run that completed is not run again (<run>.done).
set -uo pipefail
B=${BENCH:-$HOME/bench071}
REPO=${REPO:-/mnt/c/Users/hildebra/Documents/locDev/protal}
REV=${REV:-b146951}
T=${T:-6}
W=$B/world
PY=${PY:-$HOME/micromamba/envs/protal-db-build/bin/python}
PBSIM=${PBSIM:-$HOME/micromamba/envs/protal-db-build/bin/pbsim}
R=$B/runs_v073
mkdir -p $B/bin $B/src $B/logs $R

if [ ! -x $B/bin/protal-0.7.3 ]; then
  S=$B/src/0.7.3
  rm -rf $S; mkdir -p $S
  git -C $REPO archive $REV | tar -x -C $S || exit 1
  (cd $S && cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release > $B/logs/configure_0.7.3.log 2>&1 &&
    cmake --build build --target protal -j 4 > $B/logs/build_0.7.3.log 2>&1) || { echo "build of 0.7.3 failed"; exit 1; }
  cp $S/build/protal $B/bin/protal-0.7.3
  echo "$(date +%T) built 0.7.3 ($REV): $($B/bin/protal-0.7.3 --version 2>&1 | tail -1)"
fi

# 0.7.3's parity check (step 7/8) wants each feature exactly as during the collection, but adjacent_support is a sum
# over threads whose order depends on how protal shares its threads among a run's samples: the collector profiles a
# point's samples together, the check each alone, and the sums were 1.4e-15 (relative) apart, which stopped the
# pipeline. In this copy of the scripts only, a relative difference up to 1e-12 counts as equal.
sed -i 's/            if (d > 0).any():/            if (d > 1e-12).any():/' $B/src/0.7.3/scripts/check_model_parity.py
grep -q "d > 1e-12" $B/src/0.7.3/scripts/check_model_parity.py || { echo "check_model_parity.py not patched"; exit 1; }

out=$B/V073
if [ ! -f $out/protal_db/database.protal ] || [ ! -s $out/model_logs/summary.txt ]; then
  echo "$(date +%T) pipeline of 0.7.3 into $out"
  /usr/bin/time -v -o $B/logs/pipeline_0.7.3.time $PY $B/src/0.7.3/scripts/build_gtdb_database.py --gtdb $W/release_p \
    --outdir $out --protal $B/bin/protal-0.7.3 --simulator $B/bin/simulate_metagenomes -t $T --seed 1 \
    --extra-genomes $W/full/simulation/genomes_nonreps --samples 4 --read-pairs 1000,20000,200000 \
    --test-samples 2 --test-read-pairs 500,10000,500000 --long-read-bases 300000,6000000,60000000 \
    --test-long-read-bases 150000,3000000,90000000 --long-read-samples 4 --pbsim $PBSIM --no-binary-check \
    > $B/logs/pipeline_0.7.3.log 2>&1 ||
    { echo "pipeline of 0.7.3 failed, see $B/logs/pipeline_0.7.3.log"; tail -20 $B/logs/pipeline_0.7.3.log; exit 1; }
  grep -E "Elapsed|Maximum resident" $B/logs/pipeline_0.7.3.time
  cat $out/model_logs/summary.txt
fi
cmp -s $B/V071/heldout_species.txt $B/V073/heldout_species.txt || { echo "V073 held out other species than V071"; exit 1; }

declare -A DB=([v060.full]=$B/db060_full [v060.missing]=$B/db060_missing
               [v070.full]=$B/V070/protal_db [v070.missing]=$B/V070/training_db
               [v071.full]=$B/V071/protal_db [v071.missing]=$B/V071/training_db
               [v072.full]=$B/V072/protal_db [v072.missing]=$B/V072/training_db
               [v073.full]=$B/V073/protal_db [v073.missing]=$B/V073/training_db)
binary() { echo $B/bin/protal-0.${1:2:1}.${1:3:1}; }  # v073 -> protal-0.7.3
binary060=$B/bin/protal-0.6.0a

# 0.7.3's trained models into its training database (the finished database has them from its pipeline).
marker=$B/V073/training_db/.models_added
if [ ! -f $marker ]; then
  for type in pe se pb ont; do
    model=$B/V073/trained_model$([ $type = pe ] || echo _$type).xml
    $B/bin/protal-0.7.3 --add_model $model --read_type $type --db $B/V073/training_db > $B/logs/add_model_073.$type.log 2>&1 ||
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
    for v in v070 v071 v072 v073; do
      run $v.$db.pe.$id $(binary $v) --db ${DB[$v.$db]} "${pair[@]}"
    done
    run v073k.$db.pe.$id $B/bin/protal-0.7.3 --db ${DB[v073.$db]} "${pair[@]}" --knob 0.5
    case $id in rl*_p5000000_*) continue;; esac  # the deep samples: paired-end only
    single=(-1 $reads --read_type se --prefix $id --no_qcmsa)
    for v in v070 v071 v072 v073; do
      run $v.$db.se.$id $(binary $v) --db ${DB[$v.$db]} "${single[@]}"
    done
    run v073k.$db.se.$id $B/bin/protal-0.7.3 --db ${DB[v073.$db]} "${single[@]}" --knob 0.5
  done
done
echo "$(date +%T) done"
