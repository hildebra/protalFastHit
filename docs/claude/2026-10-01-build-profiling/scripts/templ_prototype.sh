#!/bin/bash
# One pbsim3 call per long-read sample (--strategy templ, templates from long_read_templates.py) against the
# collector's call per genome: the same communities and depths as a collect_real_size.sh run's pb and ont
# points, timed; reads and bases of what came out.
# usage: templ_prototype.sh COLLECT_DIR OUT
set -euo pipefail
collect=$1 out=$2
here=$(cd "$(dirname "$0")" && pwd)
pbsim=${PBSIM:-$HOME/micromamba/envs/protal-db-build/bin/pbsim}
models=$HOME/micromamba/envs/protal-db-build/data
mkdir -p "$out"
printf 'sample\tbases_asked\ttemplates_s\tpbsim_s\treads\tbases\tmean_length\n' > "$out/templ.tsv"
for samples in "$collect"/points/{pb,ont}_b*/sim/samples.tsv; do
  point=$(basename "$(dirname "$(dirname "$samples")")")
  type=${point%%_*} bases=${point#*_b}
  if [ "$type" = pb ]; then mean=15000 sd=3000 args=(--method errhmm --errhmm "$models/ERRHMM-SEQUEL.model" --accuracy-mean 0.999)
  else mean=8000 sd=6000 args=(--method qshmm --qshmm "$models/QSHMM-ONT-HQ.model" --accuracy-mean 0.97 --difference-ratio 39:24:36); fi
  tail -n +2 "$samples" | while IFS=$'\t' read -r sample reads truth community; do
    cpoint=${community%_s_*}
    manifest=$collect/points/$cpoint/sim/manifest.tsv
    d=$out/$sample; rm -rf "$d"; mkdir -p "$d"
    t0=$(date +%s.%N)
    python3 "$here/long_read_templates.py" --manifest "$manifest" --sample "$community" --bases "$bases" \
      --length-mean $mean --length-sd $sd --seed 7 -o "$d/templates.fa" > "$d/templates.log"
    t1=$(date +%s.%N)
    "$pbsim" --strategy templ "${args[@]}" --template "$d/templates.fa" --seed 7 --prefix "$d/r" --id-prefix x > "$d/pbsim.log" 2>&1
    t2=$(date +%s.%N)
    stats=$(zcat "$d"/r*.fq.gz | awk 'NR % 4 == 2 { n++; b += length($0) } END { printf "%d\t%d\t%.0f", n, b, b / n }')
    printf '%s\t%s\t%.2f\t%.2f\t%s\n' "$sample" "$bases" "$(echo "$t1 - $t0" | bc)" "$(echo "$t2 - $t1" | bc)" "$stats" >> "$out/templ.tsv"
    ls "$d" | wc -l > "$d/files"
  done
done
column -t "$out/templ.tsv"
