#!/bin/bash
# Retrain v31's pe and se GBM models on the build's own tables, as built (baseline) and with sample_log_taxa held
# constant (0) in the training and test tables, so that no tree can split on it (the trainer has no option to drop
# one feature of a group; sample_log_taxa is one of "complexity"). The trainer's options are the build's
# (trained_model*.report.txt line 3) except --threads 4 and without the C 1 / C 8 extra test tables, which the share
# archive lacks.
#
# Usage (WSL): retrain_without_taxa.sh REPO V31_FOLDER WORK_DIR
#   REPO: the protal checkout (its scripts/ at the build's commit, 956e18b)
#   V31_FOLDER: the unpacked share archive (protal0.8.5_r226_v31)
#   WORK_DIR: a new folder for the tables and the models (on the Linux file system)
set -euo pipefail
REPO=$1
A=$2
W=$3
FEATURES=normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood+ancestry+gaps+untried+alleles+polymorphic+weights+strain_identity+length_free+read_length

mkdir -p "$W"/as_built/training "$W"/as_built/test "$W"/no_taxa/training "$W"/no_taxa/test "$W"/out
for s in training test; do
    for t in "" _se; do
        cp "$A/work/$s/training_data$t.tsv" "$W/as_built/$s/"
        awk -F'\t' -v OFS='\t' 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "sample_log_taxa") c = i
                                          if (!c) { print "no sample_log_taxa column" > "/dev/stderr"; exit 1 }
                                          print; next }
                                {$c = 0; print}' "$A/work/$s/training_data$t.tsv" > "$W/no_taxa/$s/training_data$t.tsv"
    done
done

cd "$W/out"
for t in "" _se; do
    for v in as_built no_taxa; do
        name=${v}${t:-_pe}
        echo "=== $name $(date +%T)"
        /usr/bin/time -v python "$REPO/scripts/machine_learning_cmdline.py" \
            --truth-file "$W/$v/training/training_data$t.tsv" --output-prefix "$W/out/$name" \
            --features "$FEATURES" --model gbm --ntree 64 --maxnodes 63 --seed 1 --threads 4 --scenario-weight 0.25 \
            --taxonomy "$A/work/internal_taxonomy.dmp" --evaluation basic \
            --feature-versions "$A/model_logs/feature_versions.tsv" --depth-knobs \
            --test-file "$W/$v/test/training_data$t.tsv" > "$W/out/$name.log" 2> "$W/out/$name.err"
        grep -E "Maximum resident|Elapsed" "$W/out/$name.err" || true
    done
done
echo "=== done $(date +%T)"
