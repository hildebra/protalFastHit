#!/bin/bash
# Each read type's model of a build_gtdb_database.py run trained again on the run's own training table, with the
# normalized features and with the gene neighbour features as well (random_forest_cmdline.py --features
# normalized+adjacency), the build's settings otherwise, for each trainer seed of $SEEDS (default 1); all scored on
# the run's independent test set. Prints the means over the seeds (± their standard deviation).
# usage: [SEEDS="1 2 3"] compare_features.sh OUT [THREADS]
#   -> OUT/feature_comparison/<type>_<set>_s<seed>.{report.txt,metrics.json,log}
set -euo pipefail
out=$1 threads=${2:-5}
R=$HOME/bprof/gnb/src env=$HOME/micromamba/envs/protal-db-build
cmp=$out/feature_comparison
mkdir -p "$cmp"
for seed in ${SEEDS:-1}; do
  for t in pe se pb ont; do
    table=training_data$([ $t = pe ] || echo "_$t").tsv
    [ -f "$out/training/$table" ] || continue
    knobs=$([ $t = pb ] || [ $t = ont ] && echo --depth-knobs || true)
    for fs in normalized normalized+adjacency; do
      prefix=$cmp/${t}_${fs/+/_}_s$seed
      [ -f "$prefix.metrics.json" ] && continue
      "$env/bin/python3" "$R/scripts/random_forest_cmdline.py" --truth-file "$out/training/$table" \
        --test-file "$out/test/$table" --output-prefix "$prefix" --features "$fs" --ntree 64 --maxnodes 128 \
        --seed "$seed" --threads "$threads" --taxonomy "$out/internal_taxonomy.dmp" --evaluation basic $knobs \
        > "$prefix.log" 2>&1
    done
  done
done
"$env/bin/python3" - "$cmp" <<'PY'
import glob, json, os, statistics, sys
cmp = sys.argv[1]
print("| read type | features | seeds | species-held-out F1 | test F1 | test AP | test log loss | FP per test sample |")
print("|---|---|---|---|---|---|---|---|")
for t in ("pe", "se", "pb", "ont"):
    for fs in ("normalized", "normalized_adjacency"):
        runs = [json.load(open(p)) for p in sorted(glob.glob(os.path.join(cmp, f"{t}_{fs}_s*.metrics.json")))]
        if not runs:
            continue
        def cell(get, digits=4):
            values = [v for v in (get(m) for m in runs) if v is not None]
            if not values:
                return ""
            sd = f" ± {statistics.stdev(values):.{digits}f}" if len(values) > 1 else ""
            return f"{statistics.mean(values):.{digits}f}{sd}"
        test = lambda key: (lambda m: m.get("test", {}).get("this one", {}).get(key))
        held = lambda m: m.get("evaluation", {}).get("species", {}).get("F1")
        print(f"| {t} | {fs.replace('_', '+')} | {len(runs)} | {cell(held)} | {cell(test('F1'))} | {cell(test('AP'))} | "
              f"{cell(test('log_loss'))} | {cell(test('FP_per_sample'), 2)} |")
PY
