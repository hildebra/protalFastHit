#!/usr/bin/env bash
# Reproduces the measurements of this report in WSL. Usage: bash run_all.sh WORK_DIR PROTAL_CHECKOUT [COMMIT]
# Needs python3 with numpy, art_illumina is not needed; pbsim3 and its models (the protal-db-build environment).
set -euo pipefail
work=${1:?work dir}; repo=${2:?protal checkout}; commit=${3:-22c445d}
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$work" && cd "$work"
rm -rf scripts && (cd "$repo" && git archive "$commit" scripts) | tar x -C "$work"
cp "$here"/*.py "$work"/
export PATH=$HOME/micromamba/envs/protal-db-build/bin:$PATH
models=$HOME/micromamba/envs/protal-db-build/data
P=/usr/bin/python3
G=$work/genomes/genomes.tsv
[ -f "$G" ] || $P make_genomes.py "$work/genomes" --genomes 400 --threads 5
$P unit_costs.py --scripts "$work/scripts" --genomes "$G" --pbsim pbsim --pbsim_models "$models" --n 40 --community 200 | tee unit_costs_py312.txt
# the prototypes make the same reads as the collector
for setup in ultima ont; do
  for mode in threads onepass; do
    $P gil_scaling.py --scripts "$work/scripts" --genomes "$G" --out "$work/check" --mode $mode --slots 6 --setup $setup \
      --samples 1 --bases 12000000 --chunk 3000000 --community 120 --pbsim_models "$models" --keep
  done
  $P gil_scaling.py --scripts "$work/scripts" --genomes "$G" --out "$work/check" --slots 6 --setup $setup --check
done
# scaling: 1 slot, then the four modes on 6 slots, twice, alternated
rm -f "$work/scale/results.jsonl"
$P gil_scaling.py --scripts "$work/scripts" --genomes "$G" --out "$work/scale" --mode threads --slots 1 --setup ultima --pbsim_models "$models"
for rep in 1 2; do
  for mode in threads processes onepass onepass+processes; do
    $P gil_scaling.py --scripts "$work/scripts" --genomes "$G" --out "$work/scale" --mode $mode --slots 6 --setup ultima --pbsim_models "$models"
  done
done
for mode in threads processes onepass onepass+processes; do
  $P gil_scaling.py --scripts "$work/scripts" --genomes "$G" --out "$work/scale" --mode $mode --slots 6 --setup ont --samples 1 --pbsim_models "$models"
done
$P hifi_scaling.py --scripts "$work/scripts" --genomes "$G" --out "$work/hifi" --n 6 | tee hifi_scaling.txt
$P scenario_load.py --scripts "$work/scripts" --json scenario_load.json | tee scenario_load.txt
