"""The weighted site shift on the small real world (weighted_shift_runs.sh): per seed, (1) whether the build without the
shift reproduces the earlier one (~/realworld_sp_s<seed>, protal 371ac84 without the prebuilt ancestry sites): the
training and test tables byte for byte; (2) F1 with species held out and on the test set, without and with the shift,
and a paired bootstrap over the samples of all seeds (per sample TP, FP, FN of each arm's own calls, knob 0.5).

usage: compare_ws.py <home>   (the folders rw_ws_off_s1..3, rw_ws_on_s1..3, realworld_sp_s1..3 in it)
"""
import filecmp
import os
import sys
import numpy as np
import pandas as pd

H = sys.argv[1]
SEEDS = (1, 2, 3)
TYPES = ("pe", "pb")
rng = np.random.default_rng(1)


def suffix(t):
    return "" if t == "pe" else "_" + t


def calls(folder, t):
    c = pd.read_csv(f"{folder}/model_logs/trained_model{suffix(t)}.calls.tsv.gz", sep="\t",
                    usecols=["meta_sample", "meta_scenario", "taxon", "truth", "set", "p"], dtype={"meta_scenario": str})
    return c


def per_sample(c):
    call = c["p"] >= 0.5
    t = c["truth"] == 1
    g = pd.DataFrame({"s": c["set"] + ":" + c["meta_sample"], "tp": t & call, "fp": ~t & call, "fn": t & ~call})
    return g.groupby("s")[["tp", "fp", "fn"]].sum()


def f1(tp, fp, fn):
    return 2 * tp / (2 * tp + fp + fn) if tp else 0.0


print("== does the build without the shift reproduce the earlier one (realworld_sp, 371ac84)?")
for seed in SEEDS:
    old, new = f"{H}/realworld_sp_s{seed}", f"{H}/rw_ws_off_s{seed}"
    if os.path.realpath(old) == os.path.realpath(new):
        print(f"   seed {seed}: the earlier build itself (rw_ws_off_s{seed} links to it)")
        continue
    for d in ("training", "test"):
        for t in TYPES:
            a, b = f"{old}/work/{d}/training_data{suffix(t)}.tsv", f"{new}/work/{d}/training_data{suffix(t)}.tsv"
            if os.path.exists(a) and os.path.exists(b):
                same = filecmp.cmp(a, b, shallow=False)
                print(f"   seed {seed} {d} {t}: {'identical' if same else 'DIFFERENT'}")

print("== F1 at knob 0.5 by arm (species held out / test set), and the paired bootstrap of on - off over samples")
for t in TYPES:
    rows = []
    per = {"off": {}, "on": {}}
    for seed in SEEDS:
        r = {"seed": seed}
        for arm in ("off", "on"):
            c = calls(f"{H}/rw_ws_{arm}_s{seed}", t)
            for s, d in (("held out", c[c["set"] == "training"]), ("test", c[c["set"] == "test"])):
                ps = per_sample(d)
                tp, fp, fn = ps.sum()
                r[f"{arm} {s}"] = f"{f1(tp, fp, fn):.4f} ({fp}/{fn})"
                per[arm].update({(seed, s, k): v for k, v in ps.iterrows()})
        rows.append(r)
    print(f"-- {t}")
    print(pd.DataFrame(rows).to_string(index=False))
    for s in ("held out", "test"):
        keys = sorted(k for k in per["off"] if k[1] == s and k in per["on"])
        a = np.array([per["off"][k].to_numpy() for k in keys])
        b = np.array([per["on"][k].to_numpy() for k in keys])
        draws = rng.integers(0, len(keys), size=(1000, len(keys)))
        d = np.array([f1(*b[i].sum(0)) - f1(*a[i].sum(0)) for i in draws])
        print(f"   {s}: F1 off {f1(*a.sum(0)):.4f}, on {f1(*b.sum(0)):.4f}; on - off {d.mean():+.4f} "
              f"[{np.quantile(d, 0.025):+.4f}, {np.quantile(d, 0.975):+.4f}] over {len(keys)} samples")
