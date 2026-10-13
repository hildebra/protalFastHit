"""Score two builds' test tables with models fitted on either build's training table (the same feature set), at knob
0.5: whether v21's fewer false positives come from its model or from its samples.

usage: cross_score.py <read type> <name=joblib prefix>... -- <name=test table>...
"""
import json
import sys
import joblib
import numpy as np
import pandas as pd

T = sys.argv[1]
args = sys.argv[2:]
cut = args.index("--")
models = dict(a.split("=", 1) for a in args[:cut])
tables = dict(a.split("=", 1) for a in args[cut + 1:])

fitted = {}
for name, prefix in models.items():
    with open(prefix + ".metrics.json") as fh:
        cols = json.load(fh)["model"]["features"]
    fitted[name] = (joblib.load(prefix + ".joblib"), cols)

rows = []
for tname, path in tables.items():
    df = pd.read_csv(path, sep="\t", low_memory=False)
    y = df["truth"].to_numpy() == 1
    sc = df["meta_scenario"].fillna("").astype(str).to_numpy() if "meta_scenario" in df.columns else np.full(len(df), "")
    cong = pd.to_numeric(df["meta_novel_congener"], errors="coerce").to_numpy() == 1
    for mname, (model, cols) in fitted.items():
        p = model.predict_proba(df[cols].to_numpy(dtype=np.float64))[:, 1]
        call = p >= 0.5
        for part, m in (("design", sc == ""), ("scenarios", sc != ""), ("all", np.full(len(df), True))):
            tp = int((y & call & m).sum()); fp = int((~y & call & m).sum()); fn = int((y & ~call & m).sum())
            fpc = int((~y & call & m & cong).sum())
            rows.append({"test table": tname, "model": mname, "rows": part, "F1": round(2 * tp / (2 * tp + fp + fn), 4),
                         "FP": fp, "FP beside held-out": fpc, "FN": fn})
print(f"== {T}")
print(pd.DataFrame(rows).to_string(index=False))
