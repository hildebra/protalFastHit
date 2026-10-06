#!/usr/bin/env python3
"""soil_experiments.py's runs.tsv as two tables: F1 per read type and variant, and the change against "v13".

    python3 summarize_runs.py experiments/runs.tsv
"""
import sys

import pandas as pd

COLS = ["soil hold-out", "soil_shallow hold-out", "soil hold-in, species held out",
        "soil_shallow hold-in, species held out", "design test", "gut hold-out"]
r = pd.read_csv(sys.argv[1], sep="\t")
w = r.pivot_table(index=["read type", "variant"], columns="set", values="F1", sort=False)
cols = [c for c in COLS if c in w]
w = w[cols]
d = w.copy()
for rt in w.index.get_level_values(0).unique():
    d.loc[rt] = (w.loc[rt] - w.loc[(rt, "v13")]).to_numpy()
pd.set_option("display.width", 250)
print(w.round(4).to_string())
print("\nchange against v13:")
print(d.round(4).to_string())
