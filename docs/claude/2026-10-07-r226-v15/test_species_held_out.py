#!/usr/bin/env python3
"""How much of the test set's F1 comes from the final model having seen the test rows' species? The test rows scored
twice: by the final model (all training rows, as the build scores them), and by the species fold models (5 folds of
the training rows grouped by taxon, as the trainer's "species held out"): a test row whose species has training rows
is scored by the fold model that did not see that species; one whose species has none, by the mean of the five.
A real sample's species were not drawn the way the simulations drew them, so a gain of the final model over the fold
models on these rows is a gain real samples do not get.

Variants as in ablate.py (v15, no_ref, no_complexity, ...); calls at the variant's knob (ablate.py's knobs.tsv) and
at 0.5, on the design's test set and the scenarios' hold-out samples, by the species' familiarity.

    python3 test_species_held_out.py --build local/v15 --ablate ablate --variants v15,no_ref --out tsho -t 4
"""
import argparse
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold
from threadpoolctl import threadpool_limits

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "scripts"))
import machine_learning_cmdline as trainer  # noqa: E402
from ablate import TABLES, columns_of, model  # noqa: E402


def f1(y, c):
    tp, fp, fn = int((c & y).sum()), int((c & ~y).sum()), int((~c & y).sum())
    return (2 * tp / (2 * tp + fp + fn) if tp else 0.0), fp, fn


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--build", required=True)
    ap.add_argument("--ablate", required=True, help="ablate.py's output (knobs.tsv, <read type>_<variant>.tsv.gz)")
    ap.add_argument("--variants", default="v15,no_ref")
    ap.add_argument("--read-types", default="pe,ont,se,pb")
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("-t", "--threads", type=int, default=4)
    opts = ap.parse_args()
    warnings.simplefilter("ignore", pd.errors.PerformanceWarning)
    os.makedirs(opts.out, exist_ok=True)
    knobs = pd.read_csv(os.path.join(opts.ablate, "knobs.tsv"), sep="\t").drop_duplicates(["read type", "variant"],
                                                                                          keep="last")
    knob = {(r["read type"], r["variant"]): float(r["knob"]) for _, r in knobs.iterrows()}
    results = os.path.join(opts.out, "test_species_held_out.tsv")
    for rt in opts.read_types.split(","):
        train = trainer.load_table(os.path.join(opts.build, "training", TABLES[rt]))
        test = trainer.load_table(os.path.join(opts.build, "test", TABLES[rt]))
        y = train["truth"].to_numpy()
        w = np.where(trainer.scenario_of(train) == "", 1.0, 0.25)
        taxa = train["taxon"].astype(str).to_numpy()
        folds = list(GroupKFold(opts.folds).split(train, y, taxa))
        fold_of = {}
        for i, (_, te) in enumerate(folds):
            fold_of.update({t: i for t in np.unique(taxa[te])})
        test_taxa = test["taxon"].astype(str).to_numpy()
        test_fold = np.array([fold_of.get(t, -1) for t in test_taxa])
        present = set(taxa[y == 1])
        fam = np.where(np.isin(test_taxa, list(present)), "present in training",
                       np.where(test_fold >= 0, "only absent in training", "not in training"))
        part = np.where(trainer.scenario_of(test) == "", "design test", "scenario hold-out")
        yt = test["truth"].to_numpy() == 1
        for variant in opts.variants.split(","):
            began = time.time()
            cols = columns_of(variant, train.columns)
            X, Xt = train[cols].to_numpy(dtype=np.float64), test[cols].to_numpy(dtype=np.float64)
            final_p = pd.read_csv(os.path.join(opts.ablate, f"{rt}_{variant}.tsv.gz"), sep="\t", usecols=["part", "p"])
            final_p = final_p.loc[final_p["part"] == "test", "p"].to_numpy()
            held = np.zeros(len(test))
            unseen = test_fold < 0
            with threadpool_limits(limits=opts.threads, user_api="openmp"):
                for i, (tr, _) in enumerate(folds):
                    m = model(opts.seed).fit(X[tr], y[tr], sample_weight=w[tr])
                    mine = test_fold == i
                    if mine.any():
                        held[mine] = m.predict_proba(Xt[mine])[:, 1]
                    if unseen.any():
                        held[unseen] += m.predict_proba(Xt[unseen])[:, 1] / len(folds)
            k = knob[(rt, variant)]
            rows = []
            for pt in ("design test", "scenario hold-out"):
                for fm in ("all", "present in training", "only absent in training", "not in training"):
                    m = (part == pt) & ((fam == fm) if fm != "all" else True)
                    if not m.any():
                        continue
                    for scoring, p in (("final model", final_p), ("species held out", held)):
                        fk, fpk, fnk = f1(yt[m], p[m] >= k)
                        rows.append({"read type": rt, "variant": variant, "part": pt, "species": fm,
                                     "taxa": int(m.sum()), "present": int(yt[m].sum()), "scored by": scoring, "knob": k,
                                     "F1 at knob": fk, "FP": fpk, "FN": fnk, "F1 at 0.5": f1(yt[m], p[m] >= 0.5)[0]})
            pd.DataFrame(rows).to_csv(results, sep="\t", index=False, mode="a", header=not os.path.exists(results),
                                      float_format="%.5f")
            print(f"{rt} {variant}: {time.time() - began:.0f} s", flush=True)


if __name__ == "__main__":
    main()
