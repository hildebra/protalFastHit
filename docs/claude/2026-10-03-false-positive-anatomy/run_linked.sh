#!/usr/bin/env bash
# linked_share (both mates with records on the taxon; computed by protal, not in any feature set) added to the sets.
set -euo pipefail
PY=/home/falk/micromamba/envs/protal-db-build/bin/python
W=$HOME/fpdepth
cd $W/rf
grep -q LINKED_FEATURES model_features.py || $PY - <<'PYEOF'
import re
p = 'model_features.py'; s = open(p).read()
s = s.replace('GENUS_FEATURES = ["genus_species"]', 'GENUS_FEATURES = ["genus_species"]\nLINKED_FEATURES = ["linked_share"]')
s = s.replace('"normalized+adjacency+distance+depth+genus", "all")', '"normalized+adjacency+distance+depth+genus", "normalized+adjacency+distance+linked", "normalized+adjacency+distance+depth+linked", "all")')
s = s.replace('"normalized+adjacency+distance+depth+genus": NORMALIZED_FEATURES + ADJACENCY_FEATURES + DISTANCE_FEATURES + SAMPLE_FEATURES + GENUS_FEATURES}',
 '"normalized+adjacency+distance+depth+genus": NORMALIZED_FEATURES + ADJACENCY_FEATURES + DISTANCE_FEATURES + SAMPLE_FEATURES + GENUS_FEATURES,\n              "normalized+adjacency+distance+linked": NORMALIZED_FEATURES + ADJACENCY_FEATURES + DISTANCE_FEATURES + LINKED_FEATURES,\n              "normalized+adjacency+distance+depth+linked": NORMALIZED_FEATURES + ADJACENCY_FEATURES + DISTANCE_FEATURES + SAMPLE_FEATURES + LINKED_FEATURES}')
open(p, 'w').write(s)
PYEOF
grep -c "linked" model_features.py
train() {  # read type, feature set, tag
  local rt=$1 feats=$2 tag=$3 leaves=512
  $PY $W/rf/random_forest_cmdline.py --truth-file $W/tables/training$rt.tsv --output-prefix $W/out/$tag$rt \
    --features $feats --ntree 64 --maxnodes $leaves --seed 1 --threads 6 --evaluation basic \
    --depth-knobs --test-file $W/tables/test$rt.tsv > $W/out/$tag$rt.log 2>&1
  echo "done $tag$rt $(date +%T)"
}
for rt in "" _se; do
  train "$rt" normalized+adjacency+distance+linked nad_linked
  train "$rt" normalized+adjacency+distance+depth+linked nad_depth_linked
done
echo LINKEDDONE
