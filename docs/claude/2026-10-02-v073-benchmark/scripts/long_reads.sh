#!/usr/bin/env bash
# The 0.7 versions' long reads, on two sets of samples of the same design and seed (samples.sh of the v0.7.1
# benchmark: PacBio and Nanopore, 3 and 90 Mb, 4 samples each, replaying the paired-end communities: 10-300 species,
# lognormal sigma 2, strains 0.5/0.2, 6 archaea):
#   $B/samples_lr073  made now by 0.7.3's collector (--simulate_only): PacBio HiFi reads by hifi_reads.py, their
#                     quality by their length (9037c8a), Nanopore by pbsim3 as before
#   $B/samples_lr072  made by 0.7.2's collector (../../2026-10-02-v072-benchmark/scripts/long_reads.sh): PacBio by
#                     pbsim3 errhmm, every base at quality 0
# 0.7.0's to 0.7.2's PacBio models were trained on quality-0 reads, 0.7.3's on HiFi reads; real HiFi reads are Q30
# and better, so samples_lr073 is the one like real data. Every 0.7 version on samples_lr073 against both databases,
# 0.7.2 and 0.7.3 also at --knob 0.5 (their depth knobs off); 0.7.3 and 0.7.3 at 0.5 also on samples_lr072, beside
# the v0.7.2 benchmark's runs there. Runs into $B/runs_lr073 and $B/runs_lr072, as profile.sh runs them.
set -uo pipefail
B=${BENCH:-$HOME/bench071}
T=${T:-6}
W=$B/world
PY=${PY:-$HOME/micromamba/envs/protal-db-build/bin/python}
PBSIM=${PBSIM:-$HOME/micromamba/envs/protal-db-build/bin/pbsim}
C=$B/src/0.7.3/scripts/collect_training_data.py
O=$B/samples_lr073
mkdir -p $B/runs_lr073 $B/runs_lr072
if [ ! -s $O/.simulated ]; then
  $PY $C --db $B/V073/protal_db --genome_table $W/full/simulation/genomes.tsv --protal $B/bin/protal-0.7.3 \
    --simulator $B/bin/simulate_metagenomes --species_per_sample 10-300 --abundance lognormal:2.0 \
    --strains_per_species 0.5,0.2 --archaea 6 --pbsim $PBSIM -t $T -o $O --samples 4 --read_pairs 1000,10000,500000 \
    --read_setups 150:HSXt:350:50,100:HS20:300:40 --read_types pb,ont --long_read_bases 3000000,90000000 --seed 501 \
    --simulate_only > $B/logs/samples_lr073.log 2>&1 || { echo "simulation failed"; tail -20 $B/logs/samples_lr073.log; exit 1; }
  touch $O/.simulated
fi
ls $O/points

declare -A DB=([v070.full]=$B/V070/protal_db [v070.missing]=$B/V070/training_db
               [v071.full]=$B/V071/protal_db [v071.missing]=$B/V071/training_db
               [v072.full]=$B/V072/protal_db [v072.missing]=$B/V072/training_db
               [v073.full]=$B/V073/protal_db [v073.missing]=$B/V073/training_db)
binary() { echo $B/bin/protal-0.${1:2:1}.${1:3:1}; }
run() {
  local R=$1 name=$2; shift 2
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
    for v in v070 v071 v072 v073; do
      run $B/runs_lr073 $v.$db.$type.$id $(binary $v) --db ${DB[$v.$db]} "${args[@]}"
    done
    for v in v072 v073; do
      run $B/runs_lr073 ${v}k.$db.$type.$id $(binary $v) --db ${DB[$v.$db]} "${args[@]}" --knob 0.5
    done
  done
done
for reads in $B/samples_lr072/points/{pb,ont}_b*/sim/reads/*.fq.gz; do
  id=$(basename $reads .fq.gz)
  type=${id%%_*}
  echo "$(date +%T) $id (0.7.2's samples)"
  for db in full missing; do
    args=(-1 $reads --read_type $type --prefix $id --no_qcmsa)
    run $B/runs_lr072 v073.$db.$type.$id $B/bin/protal-0.7.3 --db ${DB[v073.$db]} "${args[@]}"
    run $B/runs_lr072 v073k.$db.$type.$id $B/bin/protal-0.7.3 --db ${DB[v073.$db]} "${args[@]}" --knob 0.5
  done
done
echo "$(date +%T) done"
