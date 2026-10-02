#!/bin/bash
# The trainer's time by evaluation level and table size: the tuning world's pe table (V3) as it is and K
# times as large (enlarge_table.sh), one thread, as build_gtdb_database.py runs it.
# usage: trainer_scaling.sh OUT "K1 K2 ..." [THREADS]
set -euo pipefail
out=$1 ks=$2 threads=${3:-1}
here=$(cd "$(dirname "$0")" && pwd)
# The trainer of a32b60f, which trained V3: later ones want the features e680d91 added, which V3's table lacks.
trainer=${TRAINER:-$HOME/bprof/trainer_src/scripts/random_forest_cmdline.py}
if [ ! -f "$trainer" ]; then
  mkdir -p ~/bprof/trainer_src
  git -C "$here/../../../.." archive a32b60f scripts | tar -x -C ~/bprof/trainer_src
fi
python=${PYTHON:-$HOME/micromamba/envs/protal-db-build/bin/python3}  # with scikit-learn
mkdir -p "$out"
for k in $ks; do
  table=$out/train_x$k.tsv test=$out/test_x$k.tsv
  [ -s "$table" ] || bash "$here/enlarge_table.sh" ~/tune/V3/training/training_data.tsv "$k" "$table"
  [ -s "$test" ] || bash "$here/enlarge_table.sh" ~/tune/V3/test/training_data.tsv "$k" "$test"
  for level in ${LEVELS:-full basic}; do
    log=$out/x${k}_$level.log
    [ -s "$log" ] && grep -q '^time:' "$log" && continue
    /usr/bin/time -f 'wall %e s, user %U s, sys %S s, max %M kB' -o "$out/x${k}_$level.time" \
      nice -n "${NICE:-10}" "$python" "$trainer" --truth-file "$table" --output-prefix "$out/x${k}_$level" --features normalized \
      --ntree 64 --maxnodes 128 --seed 1 --threads "$threads" --taxonomy ~/tune/V3/internal_taxonomy.dmp \
      --evaluation "$level" --test-file "$test" > "$log" 2>&1
    echo "x$k $level: $(grep '^time:' "$log") | $(cat "$out/x${k}_$level.time")"
  done
done
echo TRAINER_DONE
