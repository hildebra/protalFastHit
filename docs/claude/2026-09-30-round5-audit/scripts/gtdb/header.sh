#!/bin/bash
# columns of a training dump, and those feature_columns(all) would use
f=$(find ~/audit6/gtdb/${1:-b1}/training/points -name '*.truth_annotated' | head -1)
echo "$f"
head -1 "$f" | tr '\t' ' '
cd ~/audit6/gtdb/src/scripts
~/protal-train/bin/python - "$f" <<'EOF'
import sys, pandas as pd
from model_features import feature_columns, NON_FEATURE_COLUMNS
df = pd.read_csv(sys.argv[1], sep="\t")
cols = feature_columns(df.columns, "all")
print("feature set 'all':", len(cols), cols)
t = (df["truth"].astype(str).str.lower().isin(["1", "true"])).astype(int)
for c in cols:
    if pd.api.types.is_numeric_dtype(df[c]):
        if ((df[c] > 0).astype(int) == t).all() or (df[c].astype(float) == t).all():
            print("column equal to the label:", c)
EOF
