"""The v23 models refitted with feature groups left out (ablate.sh, ablate_all.sh) against the full set, one line per
variant: F1 at knob 0.5 with species held out and on the test set, a paired bootstrap over samples of the difference
(1000 draws), log loss, and the net change of errors by class (as variants.py of the v22 report: present taxa by the
genome simulated, real strains split by whether their species has alleles where its reads lie; absent taxa by their
neighbourhood).

usage: ablation_table.py <ablation out dir> <build dir> <read type> [<variant> ...]  (default: every <type>_*.calls.tsv.gz,
full first)
"""
import glob
import os
import sys
import numpy as np
import pandas as pd

OUT, BUILD, T = sys.argv[1], sys.argv[2], sys.argv[3]
SUF = "" if T == "pe" else "_" + T
KNOB = 0.5
names = sys.argv[4:] or ["full"] + sorted(os.path.basename(f)[len(T) + 1:-len(".calls.tsv.gz")]
                                         for f in glob.glob(f"{OUT}/{T}_*.calls.tsv.gz")
                                         if not f.endswith(f"/{T}_full.calls.tsv.gz"))
rng = np.random.default_rng(1)

tab = pd.concat([pd.read_csv(f"{BUILD}/work/{d}/training_data{SUF}.tsv", sep="\t",
                             usecols=["taxon_name", "fragments", "allele_copy_share"]) for d in ("training", "test")])
has = tab[tab["fragments"] >= 3].groupby("taxon_name")["allele_copy_share"].max()


def load(name):
    c = pd.read_csv(f"{OUT}/{T}_{name}.calls.tsv.gz", sep="\t",
                    usecols=["meta_sample", "meta_scenario", "meta_rep_genome", "meta_insilico_strain",
                             "meta_novel_congener", "meta_relative_rank", "taxon", "taxon_name", "truth", "set", "p"],
                    dtype={"meta_scenario": str, "meta_relative_rank": str})
    return c.sort_values(["set", "meta_sample", "taxon"]).reset_index(drop=True)


def f1(tp, fp, fn):
    return 2 * tp / (2 * tp + fp + fn) if tp else 0.0


def per_sample(d, call):
    t = d["truth"].to_numpy() == 1
    g = pd.DataFrame({"s": d["meta_sample"].to_numpy(), "tp": t & call, "fp": ~t & call, "fn": t & ~call})
    return g.groupby("s")[["tp", "fp", "fn"]].sum()


def logloss(d):
    p = np.clip(d["p"].to_numpy(), 1e-6, 1 - 1e-6)
    y = d["truth"].to_numpy()
    return float(-np.mean(y * np.log(p) + (1 - y) * np.log(1 - p)))


base = load(names[0])
sc = base["meta_scenario"].fillna("")
present = base["truth"] == 1
rep = pd.to_numeric(base["meta_rep_genome"], errors="coerce") == 1
ins = pd.to_numeric(base["meta_insilico_strain"], errors="coerce") == 1
cong = pd.to_numeric(base["meta_novel_congener"], errors="coerce") == 1
genus_rel = base["meta_relative_rank"].fillna("") == "genus"
alle = base["taxon_name"].map(has).fillna(0) > 0.2
groups = {"held out": (base["set"] == "training").to_numpy(),
          "test": ((base["set"] == "test") & (sc == "")).to_numpy()}
classes = {
    "FP congener": ~present & cong,          # absent beside a held-out congener
    "FP other": ~present & ~cong,            # absent, no held-out congener in the sample
    "FN real+al": present & ~rep & ~ins & alle,
    "FN real": present & ~rep & ~ins & ~alle,
    "FN insilico": present & ins,
    "FN rep": present & rep,
}
draws = {}
base_ps = {}
for g, m in groups.items():
    base_ps[g] = per_sample(base[m], base["p"].to_numpy()[m] >= KNOB)
    draws[g] = rng.integers(0, len(base_ps[g]), size=(1000, len(base_ps[g])))

rows = []
cb = base["p"].to_numpy() >= KNOB
for n in names:
    d = base if n == names[0] else load(n)
    assert (d["taxon"].to_numpy() == base["taxon"].to_numpy()).all(), n
    call = d["p"].to_numpy() >= KNOB
    r = {"variant": n}
    for g, m in groups.items():
        ps = per_sample(d[m], call[m]).loc[base_ps[g].index].to_numpy()
        b = base_ps[g].to_numpy()
        r[f"F1 {g}"] = round(f1(*ps.sum(0)), 4)
        if n != names[0]:
            dl = np.array([f1(*ps[k].sum(0)) - f1(*b[k].sum(0)) for k in draws[g]])
            r[f"d {g} [95%]"] = f"{dl.mean():+.4f} [{np.quantile(dl, 0.025):+.4f}, {np.quantile(dl, 0.975):+.4f}]"
        else:
            r[f"d {g} [95%]"] = ""
    r["logloss"] = round(logloss(d[groups["held out"]]), 4)
    ho = groups["held out"]
    for cname, m in classes.items():
        m = m.to_numpy() & ho
        err_v = (~call[m]).sum() if cname.startswith("FN") else call[m].sum()
        err_b = (~cb[m]).sum() if cname.startswith("FN") else cb[m].sum()
        r[cname] = int(err_v) if n == names[0] else f"{int(err_v - err_b):+d}"
    rows.append(r)
pd.set_option("display.width", 300)
print(f"== {T}: v23 refits, knob {KNOB}; class columns: errors with species held out (full) and the variant's net change")
print(pd.DataFrame(rows).to_string(index=False))
