#!/usr/bin/env python3
"""cross_fit.py's sets.tsv as one table per read type: F1 at each refit's own knob (rows: model and training table;
columns: the design's test set and both builds' hold-out samples per scenario), and the same at 0.5.

    python3 summarize_cross.py cross/sets.tsv
"""
import sys

import pandas as pd

ORDER = ["design test", "soil hold-out v13", "soil hold-out v14", "soil_shallow hold-out v13", "soil_shallow hold-out v14",
         "gut hold-out v13", "gut hold-out v14", "host hold-out v13", "host hold-out v14"]


def main():
    df = pd.read_csv(sys.argv[1], sep="\t")
    pd.set_option("display.width", 250)
    for value in ("F1", "F1 at 0.5"):
        for rt, g in df.groupby("read type", sort=False):
            g = g[g["set"].isin(ORDER)]
            wide = g.pivot_table(index=["model", "trained on", "knob"], columns="set", values=value, sort=False)
            wide = wide[[c for c in ORDER if c in wide.columns]]
            wide.columns = [c.replace(" hold-out", "").replace("soil_shallow", "shallow") for c in wide.columns]
            print(f"== {rt}: {value}{' at the refit knob' if value == 'F1' else ''}")
            print(wide.to_string(float_format=lambda v: f"{v:.4f}"))
            print()


if __name__ == "__main__":
    main()
