#!/usr/bin/env bash
# 0.7.5 and 0.7.3 on the v0.7.2 benchmark's strain samples (../../2026-10-02-v072-benchmark/scripts/make_strains.py:
# 8 samples, one strain each of 12 species), as that benchmark's strain_runs.sh runs them: the 8 samples in one run,
# with its strain MSAs and its own qcmsa (the env's python3 on PATH), against each version's finished database, -t 6,
# timed; paired-end, PacBio and Nanopore, the long reads also with --no_phasing. Runs into $B/strain_runs_v075.
set -uo pipefail
B=${BENCH:-$HOME/bench071}
T=${T:-6}
S=$B/strains
R=$B/strain_runs_v075
ENV=${ENV:-$HOME/micromamba/envs/protal-db-build}
export PATH=$ENV/bin:$PATH
mkdir -p $R
samples=$(ls $S/reads/*_R1.fq.gz | xargs -n1 basename | sed 's/_R1.fq.gz//' | sort -V | paste -sd,)
list() { local suffix=$1; echo $samples | tr , '\n' | sed "s|^|$S/reads/|; s|\$|$suffix|" | paste -sd,; }

run() {
  local name=$1; shift
  [ -f $R/$name.done ] && return
  rm -rf $R/$name
  mkdir -p $R/$name/strains
  cat /proc/loadavg > $R/$name.load
  echo "$(date +%T) $name"
  if /usr/bin/time -v -o $R/$name.time "$@" -o $R/$name -t $T > $R/$name.log 2>&1; then
    touch $R/$name.done
  else
    echo "  $name failed ($?), see $R/$name.log"
  fi
  rm -f $R/$name/*.sam.zst $R/$name/*.sam.gz $R/$name/*.sam
}

for v in 0.7.3 0.7.5; do
  tag=v${v//./}
  bin=$B/bin/protal-$v
  common=(--db $B/V${v//./}/protal_db --prefix $samples --qcmsa_script $B/src/$v/scripts/qcmsa.py)
  run $tag.pe $bin "${common[@]}" -1 $(list _R1.fq.gz) -2 $(list _R2.fq.gz)
  for type in pb ont; do
    run $tag.$type $bin "${common[@]}" -1 $(list _$type.fq.gz) --read_type $type
    run ${tag}np.$type $bin "${common[@]}" -1 $(list _$type.fq.gz) --read_type $type --no_phasing
  done
done
echo "$(date +%T) done"
