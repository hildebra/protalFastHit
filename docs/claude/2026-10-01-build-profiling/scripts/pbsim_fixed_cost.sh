#!/bin/bash
# Where pbsim3's cost per call goes: one 3.6 Mb genome at a depth of a shallow sample's share (0.001x),
# in 3, 40 and 300 contigs, with each model of the collector's setups; three runs each, the median.
# usage: pbsim_fixed_cost.sh OUT
set -euo pipefail
out=$1
pbsim=${PBSIM:-$HOME/micromamba/envs/protal-db-build/bin/pbsim}
models=$HOME/micromamba/envs/protal-db-build/data
mkdir -p "$out"
# One genome's sequence, cut into N contigs.
genomes=(~/bprof/world/genomes/*.fna.gz)
src=${genomes[0]}
zcat "$src" | grep -v '^>' | tr -d '\n' > "$out/seq.txt"
size=$(wc -c < "$out/seq.txt")
for n in 3 40 300; do
  per=$(( size / n + 1 ))
  fold -w "$per" "$out/seq.txt" | awk '{ print ">c" NR; print }' > "$out/g$n.fna"
done
ls -la "$out"/g*.fna "$models"/ERRHMM-SEQUEL.model "$models"/QSHMM-ONT-HQ.model | awk '{print $5, $9}'
printf 'setup\tcontigs\twall_s\tuser_s\tsys_s\tfiles\n' > "$out/pbsim_fixed.tsv"
for setup in pb ont; do
  if [ $setup = pb ]; then args=(--method errhmm --errhmm "$models/ERRHMM-SEQUEL.model" --length-mean 15000 --length-sd 3000 --accuracy-mean 0.999)
  else args=(--method qshmm --qshmm "$models/QSHMM-ONT-HQ.model" --length-mean 8000 --length-sd 6000 --accuracy-mean 0.97 --difference-ratio 39:24:36); fi
  for n in 3 40 300; do
    for rep in 1 2 3; do
      d=$out/run_${setup}_${n}_$rep; rm -rf "$d"; mkdir -p "$d"
      /usr/bin/time -f "%e\t%U\t%S" -o "$d/time" "$pbsim" --strategy wgs "${args[@]}" --genome "$out/g$n.fna" \
        --depth 0.001 --seed $rep --prefix "$d/r" --id-prefix g > "$d/log" 2>&1
      printf '%s\t%s\t%s\t%s\n' $setup $n "$(tail -1 "$d/time")" "$(ls "$d" | wc -l)" >> "$out/pbsim_fixed.tsv"
    done
  done
done
cat "$out/pbsim_fixed.tsv"
