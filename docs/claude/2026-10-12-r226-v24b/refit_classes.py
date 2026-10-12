"""The refits of v24b's pe model (refit.sh) by class, species held out, knob 0.5: each refit's class error rates on the
rows it was fitted on, and F1 from those rates at the class mix of v24b, v24 and v23 (builds_by_class.py's classes
and arithmetic), beside the builds' own models.

usage: refit_classes.py <local dir> <refit out dir> <refit name> [<refit name> ...]
"""
import glob
import sys
import pandas as pd

L, OUT, names = sys.argv[1], sys.argv[2], sys.argv[3:]
KNOB = 0.5
CLASSES = ["FN real strain, alleles", "FN real strain, none", "FN in-silico strain", "FN representative",
           "FP beside held-out congener", "FP beside present congeners only", "FP no congener in sample"]
COLS = ["meta_sample", "meta_scenario", "meta_rep_genome", "meta_insilico_strain", "meta_novel_congener",
        "meta_relative_rank", "taxon_name", "truth", "set", "p"]


def classes(c, has):
    present = c["truth"] == 1
    rep = pd.to_numeric(c["meta_rep_genome"], errors="coerce") == 1
    ins = pd.to_numeric(c["meta_insilico_strain"], errors="coerce") == 1
    cong = pd.to_numeric(c["meta_novel_congener"], errors="coerce") == 1
    genus_rel = c["meta_relative_rank"].fillna("") == "genus"
    alle = c["taxon_name"].map(has).fillna(0) > 0.2
    return dict(zip(CLASSES, [present & ~rep & ~ins & alle, present & ~rep & ~ins & ~alle, present & ins,
                              present & rep, ~present & cong, ~present & ~cong & genus_rel,
                              ~present & ~cong & ~genus_rel]))


def counts(c, has, knob=KNOB):
    """{class: (rows, errors)} of the species-held-out rows (set training) at the knob."""
    c = c[c["set"] == "training"]
    call = c["p"].to_numpy() >= knob
    res = {}
    for name, m in classes(c, has).items():
        m = m.to_numpy()
        err = (~call[m]).sum() if name.startswith("FN") else call[m].sum()
        res[name] = (int(m.sum()), int(err))
    return res


def build_dir(v):
    return glob.glob(f"{L}/{v}/protal*_r226_v*")[0]


def allele_share(d):
    tab = pd.concat([pd.read_csv(f"{d}/work/{s}/training_data.tsv", sep="\t",
                                 usecols=["taxon_name", "fragments", "allele_copy_share"]) for s in ("training", "test")])
    return tab[tab["fragments"] >= 3].groupby("taxon_name")["allele_copy_share"].max()


def f1(rates, mix):
    tot = {p: sum(mix[k][0] for k in CLASSES if k[:2] == p) for p in ("FN", "FP")}
    fn = sum(mix[k][0] * rates[k][1] / max(rates[k][0], 1) for k in CLASSES if k.startswith("FN"))
    fp = sum(mix[k][0] * rates[k][1] / max(rates[k][0], 1) for k in CLASSES if k.startswith("FP"))
    tp = tot["FN"] - fn
    return 2 * tp / (2 * tp + fp + fn)


builds = {}
frames = {}
for v in ("v24b", "v24", "v23"):
    d = build_dir(v)
    has = allele_share(d)
    c = pd.read_csv(f"{d}/model_logs/trained_model.calls.tsv.gz", sep="\t", usecols=COLS,
                    dtype={"meta_scenario": str, "meta_relative_rank": str})
    builds[v] = (has, counts(c, has))
    frames[f"build {v}"] = (c, has)
has24b = builds["v24b"][0]
models = {f"build {v}": builds[v][1] for v in builds}
for n in names:
    c = pd.read_csv(f"{OUT}/pe_{n}.calls.tsv.gz", sep="\t", usecols=COLS,
                    dtype={"meta_scenario": str, "meta_relative_rank": str})
    models[f"refit {n}"] = counts(c, has24b)
    frames[f"refit {n}"] = (c, has24b)

print("== pe, species held out, knob 0.5: class error rates (rows) of each model")
rows = []
for m, r in models.items():
    row = {"model": m}
    for k in CLASSES:
        row[k] = f"{r[k][1] / max(r[k][0], 1):.4f} ({r[k][0]})"
    rows.append(row)
print(pd.DataFrame(rows).set_index("model").T.to_string())
print("== pe, species held out: F1 from each model's class rates at each build's class mix")
for m, r in models.items():
    print(f"   {m:<20} " + "  ".join(f"mix {v} {f1(r, builds[v][1]):.4f}" for v in builds))
# Thinning the positives a model is fitted on lowers its probabilities, as a higher knob would; protal calls at the
# knob the build chose (pe 0.72-0.74). So each model at the knob best for each mix (0.30-0.90 in steps of 0.01).
print("== pe, species held out: the best F1 over knobs 0.30-0.90 from each model's class rates at each build's mix (knob)")
KNOBS = [round(0.30 + 0.01 * i, 2) for i in range(61)]
for m, (c, has) in frames.items():
    by_knob = {k: counts(c, has, k) for k in KNOBS}
    out = []
    for v in builds:
        best = max(KNOBS, key=lambda k: f1(by_knob[k], builds[v][1]))
        out.append(f"mix {v} {f1(by_knob[best], builds[v][1]):.4f} ({best:.2f})")
    print(f"   {m:<20} " + "  ".join(out))
    if m in ("build v24b", "refit full", "refit thin_real", "refit thin_ins", "build v24"):
        k = max(KNOBS, key=lambda k: f1(by_knob[k], builds["v24b"][1]))
        r = by_knob[k]
        print("      at its best knob for v24b's mix: " + ", ".join(f"{n} {r[n][1] / max(r[n][0], 1):.4f}"
                                                                for n in CLASSES))
