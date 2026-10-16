#!/usr/bin/env python3
"""Analyses for the r226 v31 report (docs/claude/2026-10-16-r226-v31).

Reads the trainers' calls files (PREFIX.calls.tsv.gz: every row's score, knob and call as protal calls by default;
training rows scored with species held out) and the collections' tables (work/{training,test}/training_data*.tsv)
of two builds' folders (the unpacked share archives), and prints:

1. v31's gbm against v29's gbm on the training samples both builds simulated alike (seed 1: the 50 scenario hold-in
   samples per read type and the design samples that did not change), with a paired bootstrap over samples;
2. the scenario samples one by one (v31; hold-in and hold-out): present taxa, false negatives, and those among the
   species with at most 2 reads on a marker gene, with the sample-level features;
3. what the training rows teach about a thin row (1-2 fragments): the share present by the sample's depth
   (sample_log_fragments) and taxa (sample_log_taxa), training rows of every read type's design and scenarios, and
   where the test samples' thin present rows fall.

Usage: pairs.py V31_FOLDER V29_FOLDER > pairs.txt
"""
import argparse
import os

import numpy as np
import pandas as pd

READ_TYPES = {"pe": "", "se": "_se", "pb": "_pb", "ont": "_ont"}
COLUMNS = ["meta_sample", "meta_scenario", "meta_marker_reads", "taxon", "truth", "set", "p", "knob", "call"]
RESAMPLES = 1000
SAME_SIMULATION = 0.95  # a shared sample name is the same simulation if this share of its present taxa is in both
THIN = 2  # marker reads (section 2) or fragments (section 3) of a thin species


def calls(folder, read_type, kind):
    """The calls file of a read type's model at the database's gate; kind "" (the default kind) or "_gbm"."""
    path = os.path.join(folder, "model_logs", f"trained_model{READ_TYPES[read_type]}{kind}.calls.tsv.gz")
    df = pd.read_csv(path, sep="\t", usecols=COLUMNS, dtype={"meta_scenario": str, "meta_sample": str},
                     keep_default_na=False)
    df["p"] = df["p"].astype(float)
    df["marker_reads"] = pd.to_numeric(df["meta_marker_reads"], errors="coerce")
    return df


def per_sample(df, call):
    """[samples x (TP, FP, FN)] of the rows given, with the sample names."""
    y = df["truth"].to_numpy() == 1
    t = pd.DataFrame({"s": df["meta_sample"].to_numpy(), "tp": call & y, "fp": call & ~y, "fn": ~call & y})
    g = t.groupby("s", sort=True)[["tp", "fp", "fn"]].sum()
    return g.index.to_numpy(), g.to_numpy(dtype=float, copy=True)


def f1(c):
    tp, fp, fn = c[..., 0], c[..., 1], c[..., 2]
    return 2 * tp / np.maximum(2 * tp + fp + fn, 1)


def paired(a, b, seed=1):
    """F1 of a and b (sample x count arrays of the same samples), b - a, its 95% interval over bootstrap resamples of the
    samples, and the share of resamples where b is better."""
    rng = np.random.default_rng(seed)
    n = len(a)
    w = rng.multinomial(n, np.full(n, 1 / n), size=RESAMPLES).astype(float)
    d = f1(w @ b) - f1(w @ a)
    return f1(a.sum(0)), f1(b.sum(0)), f1(b.sum(0)) - f1(a.sum(0)), np.quantile(d, 0.025), np.quantile(d, 0.975), \
        float((d > 0).mean())


def against_v29(v31, v29):
    print("## 1. v31 gbm against v29 gbm on the training samples both builds simulated alike")
    print(f"alike: a training sample name of both builds whose present taxa overlap by {SAME_SIMULATION:.0%} of their "
          "union or more. Rows scored with species held out (the calls files' training rows), each model at its own "
          "knob, as protal calls; a present taxon without a row in a build is a false negative there.")
    print()
    head = ["read type", "set", "samples (alike / shared names)", "present (union)", "rowless v29/v31", "knobs v29/v31",
            "v29 TP/FP/FN", "v31 TP/FP/FN", "F1 v29", "F1 v31", "v31 - v29", "95% interval", "P(v31 better)"]
    print("\t".join(head))
    for rt in READ_TYPES:
        old, new = calls(v29, rt, "_gbm"), calls(v31, rt, "")
        old, new = old[old["set"] == "training"], new[new["set"] == "training"]

        def present(df):
            d = df[df["truth"] == 1]
            return {s: set(t.tolist()) for s, t in d.groupby("meta_sample")["taxon"]}

        po, pn = present(old), present(new)
        scen_of = dict(zip(new["meta_sample"], new["meta_scenario"] != ""))
        shared = sorted(set(po) & set(pn))
        alike = [s for s in shared if len(po[s] & pn[s]) >= SAME_SIMULATION * len(po[s] | pn[s])]
        for label, scen in (("training, design", False), ("training, scenarios (hold-in)", True)):
            keep = [s for s in alike if scen_of[s] == scen]
            n_shared = sum(1 for s in shared if scen_of[s] == scen)
            if not keep:
                print("\t".join([rt, label, f"0 / {n_shared}"] + ["-"] * (len(head) - 3)))
                continue
            counts, rowless, knobs = [], [], []
            for df, own in ((old, po), (new, pn)):
                d = df[df["meta_sample"].isin(keep)]
                s, c = per_sample(d, d["call"].to_numpy() == 1)
                assert list(s) == sorted(keep)
                lost = np.array([len((po[x] | pn[x]) - own[x]) for x in s], dtype=float)
                c[:, 2] += lost
                counts.append(c)
                rowless.append(int(lost.sum()))
                knobs.append(f"{d['knob'].iloc[0]:g}")
            a, b = counts
            fa, fb, diff, lo, hi, better = paired(a, b)
            tot = [c.sum(0).astype(int) for c in counts]
            print("\t".join([rt, label, f"{len(keep)} / {n_shared}", str(int(tot[0][0] + tot[0][2])),
                             f"{rowless[0]}/{rowless[1]}", "/".join(knobs)] + ["/".join(map(str, t)) for t in tot] +
                            [f"{fa:.4f}", f"{fb:.4f}", f"{diff:+.4f}", f"[{lo:+.4f}, {hi:+.4f}]", f"{better:.3f}"]))
    print()


def sample_features(folder, read_type, columns):
    """{(set, sample): first row's values of the sample-level columns} from the collections' tables."""
    out = {}
    for s in ("training", "test"):
        path = os.path.join(folder, "work", s, f"training_data{READ_TYPES[read_type]}.tsv")
        d = pd.read_csv(path, sep="\t", usecols=["meta_sample"] + columns, dtype={"meta_sample": str})
        for sample, row in d.groupby("meta_sample")[columns].first().iterrows():
            out[(s, sample)] = row
    return out


def scenario_samples(v31):
    print("## 2. the scenario samples one by one (v31, gbm at its knob, as protal calls)")
    print(f"sigma: the sample's lognormal evenness (scenarios.py EVENNESS: 1.0, 1.5, 2.0, 2.5 by sample number in turn). "
          f"thin: present species with at most {THIN} reads on a marker gene (meta_marker_reads). p thin: the median "
          "score of the thin present rows. Rows of hold-in samples scored with species held out.")
    print()
    cols = ["sample_log_fragments", "sample_log_taxa"]
    head = ["read type", "scenario", "set", "sample", "sigma", "sample_log_fragments", "sample_log_taxa", "present",
            "FN", "FN rate", "thin", "thin FN", "thin FN rate", "p thin", "FP"]
    print("\t".join(head))
    sigmas = (1.0, 1.5, 2.0, 2.5)
    for rt in READ_TYPES:
        df = calls(v31, rt, "")
        df = df[df["meta_scenario"].isin(["moderate", "soil", "soil_shallow"])]
        feats = sample_features(v31, rt, cols)
        for (scen, s, sample), g in df.groupby(["meta_scenario", "set", "meta_sample"]):
            p = g[g["truth"] == 1]
            fn = p["call"] == 0
            thin = p["marker_reads"] <= THIN
            n = int(sample.rsplit("_s_", 1)[1])
            f = feats.get((s, sample))
            print("\t".join([rt, scen, "hold-in" if s == "training" else "hold-out", sample, f"{sigmas[(n - 1) % 4]:g}",
                             f"{f[cols[0]]:.3f}" if f is not None else "-", f"{f[cols[1]]:.3f}" if f is not None else "-",
                             str(len(p)), str(int(fn.sum())), f"{fn.mean():.3f}", str(int(thin.sum())),
                             str(int((fn & thin).sum())), f"{(fn & thin).sum() / max(thin.sum(), 1):.3f}",
                             f"{p.loc[thin, 'p'].median():.3f}" if thin.any() else "-",
                             str(int(((g['truth'] == 0) & (g['call'] == 1)).sum()))]))
    print()


def thin_prior(v31):
    print("## 3. what the training rows teach about a thin row")
    print(f"rows with at most {THIN} fragments (the feature), present or absent, of the training samples (design and "
          "scenarios; pe, se, pb and ont pooled), by the sample's sample_log_fragments and sample_log_taxa; share "
          "present, and the thin present rows of the hold-out samples in the same cells (as protal calls).")
    print()
    cols = ["meta_sample", "truth", "fragments", "sample_log_fragments", "sample_log_taxa"]
    train, test = [], []
    for rt in READ_TYPES:
        for s, out in (("training", train), ("test", test)):
            path = os.path.join(v31, "work", s, f"training_data{READ_TYPES[rt]}.tsv")
            d = pd.read_csv(path, sep="\t", usecols=cols, dtype={"meta_sample": str})
            d = d[d["fragments"] <= THIN].copy()
            d["read_type"] = rt
            out.append(d)
    train, test = pd.concat(train), pd.concat(test)
    # the hold-out samples' thin present rows with their calls
    c = pd.concat([calls(v31, rt, "").assign(read_type=rt) for rt in READ_TYPES])
    c = c[(c["set"] == "test") & (c["truth"] == 1)][["read_type", "meta_sample", "taxon", "call"]]
    tests = []
    for rt in READ_TYPES:
        path = os.path.join(v31, "work", "test", f"training_data{READ_TYPES[rt]}.tsv")
        d = pd.read_csv(path, sep="\t", usecols=["meta_sample", "taxon", "truth", "fragments", "sample_log_fragments",
                                                 "sample_log_taxa"], dtype={"meta_sample": str})
        d = d[(d["fragments"] <= THIN) & (d["truth"] == 1)].assign(read_type=rt)
        tests.append(d.merge(c, on=["read_type", "meta_sample", "taxon"]))
    tests = pd.concat(tests)
    depth_bins = [-np.inf, 4.0, 5.0, 5.5, 5.8, np.inf]
    taxa_bins = [-np.inf, 3.0, 3.6, 4.0, 4.2, np.inf]
    for d in (train, tests):
        d["depth"] = pd.cut(d["sample_log_fragments"], depth_bins)
        d["taxa"] = pd.cut(d["sample_log_taxa"], taxa_bins)
    g = train.groupby(["depth", "taxa"], observed=True)["truth"].agg(["size", "sum"])
    h = tests.groupby(["depth", "taxa"], observed=True)["call"].agg(["size", "sum"])
    print("\t".join(["sample_log_fragments", "sample_log_taxa", "training thin rows", "of them present",
                     "share present", "hold-out thin present rows", "of them called"]))
    for (dep, tax), r in g.iterrows():
        hh = h.loc[(dep, tax)] if (dep, tax) in h.index else None
        print("\t".join([str(dep), str(tax), str(int(r["size"])), str(int(r["sum"])), f"{r['sum'] / r['size']:.3f}",
                         str(int(hh["size"])) if hh is not None else "0",
                         str(int(hh["sum"])) if hh is not None else "-"]))
    for (dep, tax), hh in h.iterrows():
        if (dep, tax) not in g.index:
            print("\t".join([str(dep), str(tax), "0", "0", "-", str(int(hh["size"])), str(int(hh["sum"]))]))
    print()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("v31", help="v31's unpacked share folder (protal0.8.5_r226_v31)")
    ap.add_argument("v29", help="v29's unpacked share folder (protal0.8.4_r226_v29)")
    opts = ap.parse_args()
    against_v29(opts.v31, opts.v29)
    scenario_samples(opts.v31)
    thin_prior(opts.v31)


if __name__ == "__main__":
    main()
