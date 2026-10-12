"""Builds (and refits) at one class mix and their best knob: each model's error rates by a finer class (the present taxa
by the genome simulated and their fragments, the absent ones by the distance of the held-out congener beside them,
meta_novel_distance, or by their neighbourhood), with species held out, applied to every build's class mix, at the knob
best for that mix (0.30-0.90 in steps of 0.01) and at 0.5. The finer classes keep a difference of depth or of
congener distance between the builds' samples from passing for a difference of the models.

usage: builds_best_knob.py <local dir> <read type> <build> [<build> ...] [--refits <dir> <name> [<name> ...]]
  a refit (pe only) is <dir>/pe_<name>.calls.tsv.gz, fitted on the first build's tables (refit.sh)
"""
import glob
import sys
import numpy as np
import pandas as pd

args = sys.argv[1:]
refits, rdir = [], None
if "--refits" in args:
    i = args.index("--refits")
    rdir, refits = args[i + 1], args[i + 2:]
    args = args[:i]
L, T, builds = args[0], args[1], args[2:]
SUF = "" if T == "pe" else "_" + T
KNOBS = np.round(np.arange(0.30, 0.9001, 0.01), 2)
COLS = ["meta_sample", "meta_rep_genome", "meta_insilico_strain", "meta_novel_congener", "meta_novel_distance",
        "meta_relative_rank", "taxon_name", "truth", "set", "p"]
FRAG = [("1-10", 0, 10), ("11-100", 11, 100), (">100", 101, 1e18)]
DIST = [("<0.02", -1, 0.02), ("0.02-0.05", 0.02, 0.05), ("0.05-0.15", 0.05, 0.15), ("farther", 0.15, 1e9)]
NAMES = ([f"FN {c} {f}" for c in ("real strain, alleles", "real strain, none", "in-silico strain", "representative")
          for f, _, _ in FRAG]
         + [f"FP beside held-out congener {d}" for d, _, _ in DIST]
         + ["FP beside present congeners only", "FP no congener in sample"])
IS_FN = np.array([n.startswith("FN") for n in NAMES])


def build_dir(v):
    return glob.glob(f"{L}/{v}/protal*_r226_v*")[0]


def tables(d):
    """The training rows' fragments by (sample, taxon), and the taxa's allele share (as builds_by_class.py)."""
    tab = {s: pd.read_csv(f"{d}/work/{s}/training_data{SUF}.tsv", sep="\t",
                          usecols=["meta_sample", "taxon_name", "fragments", "allele_copy_share"])
           for s in ("training", "test")}
    both = pd.concat(tab.values())
    has = both[both["fragments"] >= 3].groupby("taxon_name")["allele_copy_share"].max()
    frag = tab["training"].drop_duplicates(["meta_sample", "taxon_name"])[["meta_sample", "taxon_name", "fragments"]]
    return frag, has


def labelled(calls, frag, has):
    """The species-held-out rows with a class label (index into NAMES) and their probability."""
    c = calls[calls["set"] == "training"].merge(frag, on=["meta_sample", "taxon_name"], how="left")
    present = (c["truth"] == 1).to_numpy()
    rep = (pd.to_numeric(c["meta_rep_genome"], errors="coerce") == 1).to_numpy()
    ins = (pd.to_numeric(c["meta_insilico_strain"], errors="coerce") == 1).to_numpy()
    cong = (pd.to_numeric(c["meta_novel_congener"], errors="coerce") == 1).to_numpy()
    dist = pd.to_numeric(c["meta_novel_distance"], errors="coerce").fillna(1.0).to_numpy()
    genus_rel = (c["meta_relative_rank"].fillna("") == "genus").to_numpy()
    alle = (c["taxon_name"].map(has).fillna(0) > 0.2).to_numpy()
    fr = c["fragments"].fillna(0).to_numpy()
    lab = np.full(len(c), -1)
    kind = np.where(rep, 3, np.where(ins, 2, np.where(alle, 0, 1)))
    fbin = np.where(fr <= 10, 0, np.where(fr <= 100, 1, 2))
    lab[present] = (kind * 3 + fbin)[present]
    dbin = np.where(dist < 0.02, 0, np.where(dist < 0.05, 1, np.where(dist < 0.15, 2, 3)))
    a = ~present
    lab[a & cong] = 12 + dbin[a & cong]
    lab[a & ~cong & genus_rel] = 16
    lab[a & ~cong & ~genus_rel] = 17
    assert (lab >= 0).all()
    return lab, c["p"].to_numpy(dtype=float)


def rates_by_knob(lab, p):
    """[knob][class] (rows, errors)."""
    n = np.bincount(lab, minlength=len(NAMES))
    out = []
    for k in list(KNOBS) + [0.5]:
        call = p >= k
        err = np.where(IS_FN[lab], ~call, call)
        out.append((n, np.bincount(lab, weights=err, minlength=len(NAMES))))
    return out


def f1(rate, mix):
    n, e = rate
    r = e / np.maximum(n, 1)
    fn = (mix * r)[IS_FN].sum()
    fp = (mix * r)[~IS_FN].sum()
    tp = mix[IS_FN].sum() - fn
    return 2 * tp / (2 * tp + fp + fn)


models, mixes = {}, {}
first = None
for v in builds:
    d = build_dir(v)
    frag, has = tables(d)
    if first is None:
        first = (frag, has)
    calls = pd.read_csv(f"{d}/model_logs/trained_model{SUF}.calls.tsv.gz", sep="\t", usecols=COLS,
                        dtype={"meta_relative_rank": str})
    lab, p = labelled(calls, frag, has)
    models[f"build {v}"] = rates_by_knob(lab, p)
    mixes[v] = np.bincount(lab, minlength=len(NAMES)).astype(float)
for n in refits:
    calls = pd.read_csv(f"{rdir}/pe_{n}.calls.tsv.gz", sep="\t", usecols=COLS, dtype={"meta_relative_rank": str})
    lab, p = labelled(calls, *first)
    models[f"refit {n}"] = rates_by_knob(lab, p)

print(f"== {T}, species held out: class sizes (share of present or absent) by build")
rows = []
for i, name in enumerate(NAMES):
    r = {"class": name}
    for v, m in mixes.items():
        tot = m[IS_FN].sum() if IS_FN[i] else m[~IS_FN].sum()
        r[v] = f"{int(m[i])} ({m[i] / tot:.1%})"
    rows.append(r)
print(pd.DataFrame(rows).to_string(index=False))
print(f"== {T}, species held out: F1 from each model's class rates at each build's mix: best knob (knob) / at 0.5")
for m, rk in models.items():
    out = []
    for v, mix in mixes.items():
        vals = [f1(rk[i], mix) for i in range(len(KNOBS))]
        b = int(np.argmax(vals))
        out.append(f"mix {v} {vals[b]:.4f} ({KNOBS[b]:.2f}) / {f1(rk[-1], mix):.4f}")
    print(f"   {m:<22} " + "   ".join(out))


# The noise: each build's samples resampled (independently, the builds' samples are not paired), the F1 of its rates at
# the first build's mix at its best knob, and the difference to the first build's (BOOT replicates; 0 skips it).
BOOT = 300
if BOOT and len(builds) > 1:
    rng = np.random.default_rng(1)
    per = {}
    for v in builds:
        d = build_dir(v)
        frag, has = tables(d)
        calls = pd.read_csv(f"{d}/model_logs/trained_model{SUF}.calls.tsv.gz", sep="\t", usecols=COLS,
                            dtype={"meta_relative_rank": str})
        c = calls[calls["set"] == "training"]
        lab, p = labelled(calls, frag, has)
        s_idx = pd.factorize(c["meta_sample"])[0]
        ns, nc = s_idx.max() + 1, len(NAMES)
        n = np.bincount(s_idx * nc + lab, minlength=ns * nc).reshape(ns, nc)
        e = np.stack([np.bincount(s_idx * nc + lab, weights=np.where(IS_FN[lab], p < k, p >= k),
                                  minlength=ns * nc).reshape(ns, nc) for k in KNOBS])
        per[v] = (n, e)

    def best_f1(n, e, mix):
        return max(f1((n, e[i]), mix) for i in range(len(KNOBS)))

    mix0 = mixes[builds[0]]
    reps = {v: [] for v in builds}
    for _ in range(BOOT):
        for v, (n, e) in per.items():
            pick = rng.integers(0, n.shape[0], n.shape[0])
            reps[v].append(best_f1(n[pick].sum(0), e[:, pick].sum(1), mix0))
    print(f"== {T}, species held out, at {builds[0]}'s mix and the best knob: {BOOT} bootstraps over samples")
    base = np.array(reps[builds[0]])
    for v in builds[1:]:
        dlt = np.array(reps[v]) - base
        print(f"   {v} - {builds[0]}: {np.median(dlt):+.4f} [{np.quantile(dlt, 0.025):+.4f}, {np.quantile(dlt, 0.975):+.4f}],"
              f" {v} better in {np.mean(dlt > 0):.0%}")
