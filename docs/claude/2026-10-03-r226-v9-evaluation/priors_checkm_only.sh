#!/usr/bin/env bash
# The priors without any cluster-size information: copies of v9's tables where cluster_genomes_log10, cluster_mean_ani
# and cluster_min_ani are 0 for every row (the ANI radius and the CheckM columns stay), then the default feature set.
set -euo pipefail
PY=$HOME/micromamba/envs/protal-db-build/bin/python
RF=$HOME/protal-fp/src/scripts/random_forest_cmdline.py
T=$HOME/v9_tables
I=$HOME/v9_tables_checkm
OUT=$HOME/v9_eval
mkdir -p $I/training $I/test
$PY - "$T" "$I" <<'EOF'
import csv, sys
src, dst = sys.argv[1], sys.argv[2]
for rt in ("", "_se"):
    for part in ("training", "test"):
        with open(f"{src}/{part}/training_data{rt}.tsv", encoding="utf-8") as fh, \
             open(f"{dst}/{part}/training_data{rt}.tsv", "w", encoding="utf-8", newline="") as out:
            rd = csv.DictReader(fh, delimiter="\t")
            wr = csv.DictWriter(out, fieldnames=rd.fieldnames, delimiter="\t", lineterminator="\n")
            wr.writeheader()
            for r in rd:
                for k in ("cluster_genomes_log10", "cluster_mean_ani", "cluster_min_ani"):
                    r[k] = "0"
                wr.writerow(r)
EOF
for rt in "" _se; do
  $PY $RF --truth-file $I/training/training_data$rt.tsv --output-prefix $OUT/checkm$rt \
    --features normalized+adjacency+distance+depth+divergence+unfiltered+priors --ntree 64 --maxnodes 512 --seed 1 \
    --threads 6 --evaluation basic --depth-knobs --test-file $I/test/training_data$rt.tsv > $OUT/checkm$rt.log 2>&1
  echo "done checkm$rt $(date +%T)"
done
cd /mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/3611c5a5-4e2b-44f3-a684-44ae56dc8a7c/scratchpad
sed -i 's/TAGS = \["v8set", "no_unfiltered", "no_priors", "default", "no_depth", "noflag"\]/TAGS = ["v8set", "no_unfiltered", "no_priors", "default", "no_depth", "noflag", "checkm"]/; s/"noflag": "default, priors without the singleton flag"}/"noflag": "default, priors without the singleton flag", "checkm": "default, priors = CheckM + radius only"}/' ablation_v9_summary.py
sed -i 's/("noflag", "default, priors without the singleton flag")\]/("noflag", "default, priors without the singleton flag"), ("checkm", "default, priors = CheckM + radius only")]/' ablation_by_cluster.py
python3 ablation_v9_summary.py ~/v9_eval ~/v9_tables > ablation_v9_summary_output.txt 2>&1
python3 ablation_by_cluster.py ~/v9_eval ~/v9_tables > ablation_by_cluster_output.txt 2>&1
grep -h 'CheckM + radius' ablation_v9_summary_output.txt | cut -c1-300
echo ---
grep -h 'CheckM + radius' ablation_by_cluster_output.txt | cut -c1-420
echo CHECKMDONE
