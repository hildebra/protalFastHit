#!/usr/bin/env bash
# The phasing benchmark's runs: the 16 samples of make_mixtures.py (pure1..8, mix1..8) per read type in one protal
# run against the strain benchmark's 0.7.2 database (all 765 species), with phasing (the default) and with
# --no_phasing, with qcmsa (the env's python3 on PATH), -t 6, niced. Into $B/phasing/runs/<pb|ont>.<phase|nophase>;
# a run that completed is not run again (<run>.done).
#   PROTAL: the binary (default ~/protal-hap/build/protal), QCMSA: its qcmsa.py
set -uo pipefail
B=${BENCH:-$HOME/bench071}
P=${PROTAL:-$HOME/protal-hap/build/protal}
QC=${QCMSA:-$HOME/protal-hap/src/scripts/qcmsa.py}
DB=${DB:-$B/V072/protal_db}
ENV=${ENV:-$HOME/micromamba/envs/protal-db-build}
export PATH=$ENV/bin:$PATH
R=$B/phasing/runs
mkdir -p $R
samples=$(printf 'pure%d,' 1 2 3 4 5 6 7 8; printf 'mix%d,' 1 2 3 4 5 6 7 8)
samples=${samples%,}
for type in pb ont; do
  files=$(echo $samples | tr , '\n' | sed "s|.*|$B/phasing/reads/&_$type.fq.gz|" | paste -sd,)
  for arm in phase nophase; do
    name=$type.$arm
    [ -f $R/$name.done ] && continue
    extra=()
    [ $arm = nophase ] && extra=(--no_phasing)
    rm -rf $R/$name
    echo "$(date +%T) $name"
    if /usr/bin/time -v -o $R/$name.time nice $P --db $DB -1 $files --read_type $type --prefix $samples \
        --qcmsa_script $QC -o $R/$name -t 6 "${extra[@]}" > $R/$name.log 2>&1; then
      touch $R/$name.done
    else
      echo "  $name failed, see $R/$name.log"
    fi
    rm -f $R/$name/*.sam.zst $R/$name/*.sam.gz $R/$name/*.sam
  done
done
echo "$(date +%T) done"
