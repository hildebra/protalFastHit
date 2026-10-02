#!/usr/bin/env bash
# The strain benchmark (make_strains.py: 8 samples, one strain each of 12 species): every version profiles the 8
# samples in one run, with its strain MSAs and its own qcmsa post-filter (the env's python3 on PATH), against its
# finished database (all 765 species), -t 6, timed:
#   v060.pe 0.6.0a, v070.pe 0.7.0, v071.pe 0.7.1, v072.pe 0.7.2   paired-end
#   v07N.pb, v07N.ont                                              PacBio and Nanopore, the 0.7 versions only
# Runs into $B/strain_runs/<run>; one that completed is not run again (<run>.done).
set -uo pipefail
B=${BENCH:-$HOME/bench071}
T=${T:-6}
S=$B/strains
R=$B/strain_runs
ENV=${ENV:-$HOME/micromamba/envs/protal-db-build}
export PATH=$ENV/bin:$PATH
mkdir -p $R
declare -A DB=([v060]=$B/db060_full [v070]=$B/V070/protal_db [v071]=$B/V071/protal_db [v072]=$B/V072/protal_db)
declare -A VER=([v060]=0.6.0a [v070]=0.7.0 [v071]=0.7.1 [v072]=0.7.2)
samples=$(ls $S/reads/*_R1.fq.gz | xargs -n1 basename | sed 's/_R1.fq.gz//' | sort -V | paste -sd,)
list() { local suffix=$1; echo $samples | tr , '\n' | sed "s|^|$S/reads/|; s|\$|$suffix|" | paste -sd,; }

run() {
  local name=$1; shift
  [ -f $R/$name.done ] && return
  rm -rf $R/$name
  mkdir -p $R/$name/strains  # 0.6.0a writes no MSA unless its strains folder exists
  cat /proc/loadavg > $R/$name.load
  echo "$(date +%T) $name"
  if /usr/bin/time -v -o $R/$name.time "$@" -o $R/$name -t $T > $R/$name.log 2>&1; then
    touch $R/$name.done
  else
    echo "  $name failed ($?), see $R/$name.log"
  fi
  rm -f $R/$name/*.sam.zst $R/$name/*.sam.gz $R/$name/*.sam
}

for v in v060 v070 v071 v072; do
  bin=$B/bin/protal-${VER[$v]}
  qc=(--qcmsa_script $B/src/${VER[$v]}/scripts/qcmsa.py)
  run $v.pe $bin --db ${DB[$v]} -1 $(list _R1.fq.gz) -2 $(list _R2.fq.gz) --prefix $samples "${qc[@]}"
  [ $v = v060 ] && continue
  for type in pb ont; do
    run $v.$type $bin --db ${DB[$v]} -1 $(list _$type.fq.gz) --read_type $type --prefix $samples "${qc[@]}"
  done
done
echo "$(date +%T) done"
