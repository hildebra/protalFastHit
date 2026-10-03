#!/usr/bin/env bash
# The priors without the singleton flag: copies of v9's tables where cluster_genomes_log10 is 0 for every row and an
# unknown mean or minimum intra-species ANI (-1, a cluster of one genome) is the median of the known values, so that
# the forest cannot tell singleton clusters apart, then the default feature set trained on them. pe and se.
set -euo pipefail
PY=$HOME/micromamba/envs/protal-db-build/bin/python
RF=$HOME/protal-fp/src/scripts/random_forest_cmdline.py
T=$HOME/v9_tables
I=$HOME/v9_tables_imputed
OUT=$HOME/v9_eval
mkdir -p $I/training $I/test
$PY - "$T" "$I" <<'EOF'
import csv, statistics, sys
src, dst = sys.argv[1], sys.argv[2]
for rt in ("", "_se"):
    known = {"cluster_mean_ani": [], "cluster_min_ani": []}
    with open(f"{src}/training/training_data{rt}.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for k in known:
                v = float(r[k])
                if v >= 0:
                    known[k].append(v)
    med = {k: statistics.median(v) for k, v in known.items()}
    print(rt or "pe", "medians of the known values:", med)
    for part in ("training", "test"):
        with open(f"{src}/{part}/training_data{rt}.tsv", encoding="utf-8") as fh, \
             open(f"{dst}/{part}/training_data{rt}.tsv", "w", encoding="utf-8", newline="") as out:
            rd = csv.DictReader(fh, delimiter="\t")
            wr = csv.DictWriter(out, fieldnames=rd.fieldnames, delimiter="\t", lineterminator="\n")
            wr.writeheader()
            n = 0
            for r in rd:
                r["cluster_genomes_log10"] = "0"
                for k, m in med.items():
                    if float(r[k]) < 0:
                        r[k] = f"{m:g}"
                        n += 1
                wr.writerow(r)
            print(f"  {part}: {n} values imputed")
EOF
for rt in "" _se; do
  $PY $RF --truth-file $I/training/training_data$rt.tsv --output-prefix $OUT/noflag$rt \
    --features normalized+adjacency+distance+depth+divergence+unfiltered+priors --ntree 64 --maxnodes 512 --seed 1 \
    --threads 6 --evaluation basic --depth-knobs --test-file $I/test/training_data$rt.tsv > $OUT/noflag$rt.log 2>&1
  echo "done noflag$rt $(date +%T)"
done
echo NOFLAGDONE
