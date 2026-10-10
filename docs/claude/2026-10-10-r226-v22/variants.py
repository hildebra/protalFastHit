"""Compare v22 models retrained with feature groups left out (ablate.sh): F1 with species held out and on the test set,
log loss, the scenarios' hold-out samples, a paired bootstrap over samples of each variant's F1 against the base, the
FN and FP by class, and the calls that flip between the base and each variant, by class.

Classes. Present taxa: the genome simulated (representative, in-silico strain, real strain) and, for real strains,
whether the species has alleles where its reads lie (allele_copy_share above 0.2 in a row of 3 or more fragments).
Absent taxa: beside a held-out congener (meta_novel_congener), beside present congeners only (closest present relative
in the genus), or with no congener in the sample (the closest relative beyond the genus).

usage: variants.py <ablation out dir> <v22 dir> <read type> <variant> [<variant> ...]   (the first variant is the base)
"""
import sys
import numpy as np
import pandas as pd

OUT, V22, T = sys.argv[1], sys.argv[2], sys.argv[3]
NAMES = sys.argv[4:]
SUF = "" if T == "pe" else "_" + T
KNOB = 0.5
rng = np.random.default_rng(1)

cols = ["taxon_name", "fragments", "allele_copy_share"]
tab = pd.concat([pd.read_csv(f"{V22}/work/{d}/training_data{SUF}.tsv", sep="\t", usecols=cols)
                 for d in ("training", "test")])
has = tab[tab["fragments"] >= 3].groupby("taxon_name")["allele_copy_share"].max()


def load(name):
    c = pd.read_csv(f"{OUT}/{T}_{name}.calls.tsv.gz", sep="\t",
                    usecols=["meta_sample", "meta_scenario", "meta_rep_genome", "meta_insilico_strain",
                             "meta_novel_congener", "meta_relative_rank", "taxon", "taxon_name", "truth", "set", "p"],
                    dtype={"meta_scenario": str, "meta_relative_rank": str})
    return c.sort_values(["set", "meta_sample", "taxon"]).reset_index(drop=True)


def f1(tp, fp, fn):
    return 2 * tp / (2 * tp + fp + fn) if tp else 0.0


def counts(d, call):
    t = d["truth"].to_numpy() == 1
    return int((t & call).sum()), int((~t & call).sum()), int((t & ~call).sum())


def per_sample(d, call):
    t = d["truth"].to_numpy() == 1
    g = pd.DataFrame({"s": d["meta_sample"].to_numpy(), "tp": t & call, "fp": ~t & call, "fn": t & ~call})
    return g.groupby("s")[["tp", "fp", "fn"]].sum()


def logloss(d):
    p = np.clip(d["p"].to_numpy(), 1e-6, 1 - 1e-6)
    y = d["truth"].to_numpy()
    return float(-np.mean(y * np.log(p) + (1 - y) * np.log(1 - p)))


data = {n: load(n) for n in NAMES}
base = data[NAMES[0]]
for n in NAMES[1:]:
    assert (data[n]["taxon"].to_numpy() == base["taxon"].to_numpy()).all(), n

sc = base["meta_scenario"].fillna("")
rep = pd.to_numeric(base["meta_rep_genome"], errors="coerce") == 1
ins = pd.to_numeric(base["meta_insilico_strain"], errors="coerce") == 1
cong = pd.to_numeric(base["meta_novel_congener"], errors="coerce") == 1
genus_rel = base["meta_relative_rank"].fillna("") == "genus"
alle = base["taxon_name"].map(has).fillna(0) > 0.2
present = base["truth"] == 1
groups = {
    "held out": base["set"] == "training",
    "test": (base["set"] == "test") & (sc == ""),
}
for s in sorted(set(sc[(base["set"] == "test") & (sc != "")])):
    groups[f"{s} hold-out"] = (base["set"] == "test") & (sc == s)

print(f"== {T}: F1 at knob {KNOB} (FP/FN); log loss with species held out")
rows = []
for n in NAMES:
    d = data[n]
    call = d["p"].to_numpy() >= KNOB
    r = {"variant": n, "logloss held out": round(logloss(d[groups["held out"]]), 5)}
    for g, m in groups.items():
        m = m.to_numpy()
        tp, fp, fn = counts(d[m], call[m])
        r[g] = f"{f1(tp, fp, fn):.4f} ({fp}/{fn})"
    rows.append(r)
print(pd.DataFrame(rows).to_string(index=False))

print(f"-- {T}: paired bootstrap over samples (1000), F1 of the variant - {NAMES[0]}: mean [2.5%, 97.5%]")
for g in ("held out", "test"):
    m = groups[g].to_numpy()
    ps = {n: per_sample(data[n][m], data[n]["p"].to_numpy()[m] >= KNOB) for n in NAMES}
    idx = ps[NAMES[0]].index
    draws = rng.integers(0, len(idx), size=(1000, len(idx)))
    b = ps[NAMES[0]].loc[idx].to_numpy()
    fb = np.array([f1(*b[k].sum(0)) for k in draws])
    for n in NAMES[1:]:
        v = ps[n].loc[idx].to_numpy()
        fv = np.array([f1(*v[k].sum(0)) for k in draws])
        dlt = fv - fb
        print(f"   {g:9s} {n:16s} {dlt.mean():+.4f} [{np.quantile(dlt, 0.025):+.4f}, {np.quantile(dlt, 0.975):+.4f}]")

classes = {
    "FN real strain, alleles": present & ~rep & ~ins & alle,
    "FN real strain, none": present & ~rep & ~ins & ~alle,
    "FN in-silico strain": present & ins,
    "FN representative": present & rep,
    "FP beside held-out congener": ~present & cong,
    "FP beside present congeners only": ~present & ~cong & genus_rel,
    "FP no congener in sample": ~present & ~cong & ~genus_rel,
}
for g in ("held out", "test"):
    print(f"-- {T}: {g}, by class: FN rate of present taxa, FP rate of absent ones (count)")
    rows = []
    for cname, m in classes.items():
        m = (m & groups[g]).to_numpy()
        r = {"class": cname, "n": int(m.sum())}
        for n in NAMES:
            call = data[n]["p"].to_numpy()[m] >= KNOB
            err = (~call).sum() if cname.startswith("FN") else call.sum()
            r[n] = f"{err / max(m.sum(), 1):.4f} ({err})"
        rows.append(r)
    print(pd.DataFrame(rows).to_string(index=False))

print(f"-- {T}: species held out, the calls that flip from {NAMES[0]} to each variant, by class "
      "(fixed: an error of the base the variant gets right; new: an error the variant adds)")
ho = groups["held out"].to_numpy()
cb = base["p"].to_numpy() >= KNOB
for n in NAMES[1:]:
    cv = data[n]["p"].to_numpy() >= KNOB
    rows = []
    for cname, m in classes.items():
        m = m.to_numpy() & ho
        if cname.startswith("FN"):
            fixed, new = (m & ~cb & cv).sum(), (m & cb & ~cv).sum()
        else:
            fixed, new = (m & cb & ~cv).sum(), (m & ~cb & cv).sum()
        rows.append({"class": cname, "fixed": int(fixed), "new": int(new), "net (variant - base errors)": int(new - fixed)})
    print(f"   {n}")
    print(pd.DataFrame(rows).to_string(index=False))
