#!/usr/bin/env python3
"""ablate.py's variants on the test rows split by the species' familiarity: present in the training table, only
absent there, or not in it at all. A final model can learn which species tend to be present (features constant
within a species can name it); the rows of species not in training show what carries over to species it never saw.

Per read type, variant and part (the design's test set; the scenarios' hold-out samples, pooled): F1 at the variant's
knob (ablate.py's knobs.tsv), at 0.5, and the best F1 with its threshold.

    python3 familiarity_variants.py --dir ablate > familiarity_variants.txt
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd


def f1(y, c):
    tp, fp, fn = int((c & y).sum()), int((c & ~y).sum()), int((~c & y).sum())
    return (2 * tp / (2 * tp + fp + fn) if tp else 0.0), fp, fn


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--dir", required=True)
    opts = ap.parse_args()
    knobs = pd.read_csv(os.path.join(opts.dir, "knobs.tsv"), sep="\t").drop_duplicates(["read type", "variant"], keep="last")
    knob = {(r["read type"], r["variant"]): float(r["knob"]) for _, r in knobs.iterrows()}
    grid = np.round(np.arange(0.05, 0.96, 0.01), 2)
    rows = []
    for path in sorted(glob.glob(os.path.join(opts.dir, "*_*.tsv.gz"))):
        rt, variant = os.path.basename(path)[:-len(".tsv.gz")].split("_", 1)
        d = pd.read_csv(path, sep="\t", usecols=["part", "meta_scenario", "taxon", "truth", "p"],
                        dtype={"meta_scenario": str}, keep_default_na=False, na_values=[""])
        tr, te = d[d["part"] == "training"], d[d["part"] == "test"]
        present, seen = set(tr.loc[tr["truth"] == 1, "taxon"]), set(tr["taxon"])
        fam = np.where(te["taxon"].isin(present), "present in training",
                       np.where(te["taxon"].isin(seen), "only absent in training", "not in training"))
        part = np.where(te["meta_scenario"].fillna("design").to_numpy() == "design", "design test", "scenario hold-out")
        y, p = te["truth"].to_numpy() == 1, te["p"].to_numpy()
        k = knob[(rt, variant)]
        for pt in ("design test", "scenario hold-out"):
            for fm in ("present in training", "only absent in training", "not in training", "all"):
                m = (part == pt) & ((fam == fm) if fm != "all" else True)
                if not m.any():
                    continue
                fk, fpk, fnk = f1(y[m], p[m] >= k)
                f5 = f1(y[m], p[m] >= 0.5)[0]
                best = max((f1(y[m], p[m] >= t)[0], t) for t in grid)
                rows.append({"read type": rt, "variant": variant, "part": pt, "species": fm, "taxa": int(m.sum()),
                             "present": int(y[m].sum()), "knob": k, "F1 at knob": fk, "FP": fpk, "FN": fnk,
                             "F1 at 0.5": f5, "best F1": best[0], "at": best[1]})
    out = pd.DataFrame(rows).sort_values(["read type", "part", "species", "variant"])
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 1000)
    print(out.to_string(index=False, float_format=lambda v: f"{v:.4f}"))


if __name__ == "__main__":
    main()
