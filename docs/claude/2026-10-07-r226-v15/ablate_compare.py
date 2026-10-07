#!/usr/bin/env python3
"""ablate.py's variants side by side, with a paired bootstrap over samples for each difference in F1.

Each variant is scored at its own knob (OUT/knobs.tsv, chosen as the trainer chooses its global one) on: the design's
test set and the scenarios' hold-out samples (final model, p), and the training rows per scenario with species held
out (p_species) and with samples held out (p_samples). Differences are B - A for each --pairs A:B; the 95% interval
resamples the set's samples with replacement (2000 times). With --build, the v15 variant's test scores are checked
against the build's own (*.test_predictions.tsv.gz).

    python3 ablate_compare.py --dir ablate --pairs no_complexity:v15,v15:no_depth --build local/v15 > ablate_compare.txt
"""
import argparse
import os

import numpy as np
import pandas as pd

SUFFIX = {"pe": "", "se": "_se", "pb": "_pb", "ont": "_ont"}


def counts(samples, y, p, knob):
    c = p >= knob
    g = pd.DataFrame({"s": samples, "tp": c & y, "fp": c & ~y, "fn": ~c & y})
    return g.groupby("s")[["tp", "fp", "fn"]].sum().sort_index().to_numpy(dtype=float)


def f1(m):
    tp, fp, fn = m[..., 0], m[..., 1], m[..., 2]
    return np.where(tp > 0, 2 * tp / np.maximum(2 * tp + fp + fn, 1), 0.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--dir", required=True)
    ap.add_argument("--pairs", default="no_complexity:v15,v15:no_depth")
    ap.add_argument("--read-types", default="pe,se,pb,ont")
    ap.add_argument("--build", help="the build whose test predictions the v15 variant should reproduce")
    ap.add_argument("--reps", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--fixed-knob", type=float, help="score every variant at this knob instead of its own")
    opts = ap.parse_args()
    rng = np.random.default_rng(opts.seed)
    knobs = pd.read_csv(os.path.join(opts.dir, "knobs.tsv"), sep="\t").drop_duplicates(["read type", "variant"], keep="last")
    knob = {(r["read type"], r["variant"]): opts.fixed_knob or float(r["knob"]) for _, r in knobs.iterrows()}
    rows = []
    for rt in opts.read_types.split(","):
        variants = {a for pair in opts.pairs.split(",") for a in pair.split(":")}
        if not all(os.path.exists(os.path.join(opts.dir, f"{rt}_{v}.tsv.gz")) for v in variants):
            continue
        data = {v: pd.read_csv(os.path.join(opts.dir, f"{rt}_{v}.tsv.gz"), sep="\t", dtype={"meta_scenario": str},
                               keep_default_na=False, na_values=[""]) for v in variants}
        if opts.build and "v15" in data:
            own = pd.read_csv(os.path.join(opts.build, "model_logs", f"trained_model{SUFFIX[rt]}.test_predictions.tsv.gz"),
                              sep="\t", usecols=["meta_sample", "taxon", "p"])
            mine = data["v15"][(data["v15"]["part"] == "test")]
            mine = mine[mine["meta_scenario"].fillna("") == ""]
            same = (own["taxon"].to_numpy() == mine["taxon"].to_numpy()).all()
            print(f"# {rt}: v15 variant against the build's design test scores: rows aligned {same}, max |p diff| "
                  f"{np.abs(own['p'].to_numpy() - mine['p'].to_numpy()).max():.2g}; knob {knob[(rt, 'v15')]}")
        base = data[next(iter(variants))]
        part, scen = base["part"].to_numpy(), base["meta_scenario"].fillna("design").to_numpy()
        sets = [("design test", (part == "test") & (scen == "design"), "p")]
        sets += [(f"{s} hold-out", (part == "test") & (scen == s), "p") for s in sorted(set(scen[part == "test"]) - {"design"})]
        for col, label in (("p_species", "species held out"), ("p_samples", "samples held out")):
            sets += [(f"{s} training, {label}", (part == "training") & (scen == s), col)
                     for s in ["design"] + sorted(set(scen[part == "training"]) - {"design"})]
        for name, m, col in sets:
            samples = base["meta_sample"].to_numpy()[m]
            y = base["truth"].to_numpy()[m] == 1
            n = len(np.unique(samples))
            draws = rng.multinomial(n, np.full(n, 1 / n), size=opts.reps).astype(float)
            for pair in opts.pairs.split(","):
                a, b = pair.split(":")
                ca = counts(samples, y, data[a][col].to_numpy()[m], knob[(rt, a)])
                cb = counts(samples, y, data[b][col].to_numpy()[m], knob[(rt, b)])
                d = f1(draws @ cb) - f1(draws @ ca)
                lo, hi = np.percentile(d, [2.5, 97.5])
                rows.append({"read type": rt, "set": name, "A": a, "B": b, "knob A": knob[(rt, a)], "knob B": knob[(rt, b)],
                             "samples": n, "F1 A": float(f1(ca.sum(0))), "F1 B": float(f1(cb.sum(0))),
                             "B - A": float(f1(cb.sum(0)) - f1(ca.sum(0))), "95% low": lo, "95% high": hi,
                             "FP A": int(ca[:, 1].sum()), "FP B": int(cb[:, 1].sum()),
                             "FN A": int(ca[:, 2].sum()), "FN B": int(cb[:, 2].sum())})
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 1000)
    print(pd.DataFrame(rows).to_string(index=False, float_format=lambda v: f"{v:.4f}"))


if __name__ == "__main__":
    main()
