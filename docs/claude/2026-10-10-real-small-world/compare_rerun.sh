#!/bin/bash
# The rerun with the species' shared sites (371ac84, ~/realworld_sp_s1..3) against the first run (6f8e6bb,
# ~/realworld_s1..3): the first run's full models' calls linked into the rerun's refits as "before" (the same samples
# and taxa, so variants.py pairs them), the variants pooled, and each rerun build's features alone.
HERE=$(cd "$(dirname "$0")" && pwd)
PY=~/micromamba/envs/protal-db-build/bin/python
for s in 1 2 3; do
  for t in pe pb; do ln -sf ~/realworld_s${s}_abl/${t}_full.calls.tsv.gz ~/realworld_sp_s${s}_abl/${t}_before.calls.tsv.gz; done
done
B=~/realworld_sp_s1,~/realworld_sp_s2,~/realworld_sp_s3
for t in pe pb; do
  taskset -c 0-3 nice -n 5 $PY -I $HERE/variants.py $t $B before full > $HERE/rerun_vs_before_$t.txt 2>&1
  taskset -c 0-3 nice -n 5 $PY -I $HERE/variants.py $t $B full no_ancestry no_polymorphic no_weights no_ancestral no_alleles > $HERE/rerun_variants_$t.txt 2>&1
done
for s in 1 2 3; do
  taskset -c 0-3 $PY -I $HERE/../2026-10-10-r226-v22/feature_auc.py ~/realworld_sp_s$s pe > $HERE/rerun_feature_auc_pe_s$s.txt 2>&1
done
echo "== COMPARE DONE"
