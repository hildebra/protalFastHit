#!/usr/bin/env bash
# The long-read benchmark samples made again by 0.7.2's collector, which draws each long read's genome, length and
# start itself (pbsim3 --strategy templ): 0.7.1's collector had pbsim3 sample once per genome, which cut every
# contig's last read to the contig's share of the bases, so that its 3 Mb samples were half reads under 1 kb (mean
# 2.3 kb for PacBio, 1.7 kb for Nanopore, against 15 and 8 kb meant). 0.7.0's and 0.7.1's models were trained on such
# reads, 0.7.2's on the drawn ones, so on the old samples the comparison of the 0.7 versions' long-read models is
# not fair. The design and seed of samples.sh (the v0.7.1 benchmark): PacBio and Nanopore, 3 and 90 Mb, 4 samples
# each, replaying the paired-end communities (10-300 species, lognormal sigma 2, strains 0.5/0.2, 6 archaea); into
# $B/samples_lr072 (--simulate_only). Then every 0.7 version on them against both databases, and 0.7.2 again at
# --knob 0.5 (its depth knobs off): runs into $B/runs_lr072, as profile.sh runs them.
set -uo pipefail
B=${BENCH:-$HOME/bench071}
T=${T:-6}
W=$B/world
PY=${PY:-$HOME/micromamba/envs/protal-db-build/bin/python}
PBSIM=${PBSIM:-$HOME/micromamba/envs/protal-db-build/bin/pbsim}
C=$B/src/0.7.2/scripts/collect_training_data.py
O=$B/samples_lr072
R=$B/runs_lr072
mkdir -p $R
if [ ! -s $O/.simulated ]; then
  $PY $C --db $B/V072/protal_db --genome_table $W/full/simulation/genomes.tsv --protal $B/bin/protal-0.7.2 \
    --simulator $B/bin/simulate_metagenomes --species_per_sample 10-300 --abundance lognormal:2.0 \
    --strains_per_species 0.5,0.2 --archaea 6 --pbsim $PBSIM -t $T -o $O --samples 4 --read_pairs 1000,10000,500000 \
    --read_setups 150:HSXt:350:50,100:HS20:300:40 --read_types pb,ont --long_read_bases 3000000,90000000 --seed 501 \
    --simulate_only > $B/logs/samples_lr072.log 2>&1 || { echo "simulation failed"; tail -20 $B/logs/samples_lr072.log; exit 1; }
  touch $O/.simulated
fi
ls $O/points

declare -A DB=([v070.full]=$B/V070/protal_db [v070.missing]=$B/V070/training_db
               [v071.full]=$B/V071/protal_db [v071.missing]=$B/V071/training_db
               [v072.full]=$B/V072/protal_db [v072.missing]=$B/V072/training_db)
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
    for v in v070 v071 v072; do
      run $v.$db.$type.$id $(binary $v) --db ${DB[$v.$db]} -1 $reads --read_type $type --prefix $id --no_qcmsa
    done
    run v072k.$db.$type.$id $B/bin/protal-0.7.2 --db ${DB[v072.$db]} -1 $reads --read_type $type --prefix $id \
      --no_qcmsa --knob 0.5
  done
done
echo "$(date +%T) done"
