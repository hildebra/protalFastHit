#!/usr/bin/env python3
"""Paired comparisons for the r226 v29 report (docs/claude/2026-10-15-r226-v29).

Reads the trainers' calls files (PREFIX.calls.tsv.gz: every row's score, knob and call as protal calls by default;
training rows scored with species held out) of one or two builds' model_logs, and prints:

1. forest against gbm in v29 (the same rows), per read type and set: F1 as protal calls (each kind at its own knob),
   the difference with a paired bootstrap over samples, and the true positives of each kind at the other's false
   positive count (matched FP);
2. v29 (forest and gbm) against v27 (gbm) on the test samples both builds simulated alike (the same present taxa in
   each, up to a few without a row), with the same bootstrap;
3. the paired-end test samples by read length.

The singleton rule is off in these builds (SINGLETON_CONGENER = 0), so a call is p >= knob.

Usage: pairs.py V29_MODEL_LOGS V27_MODEL_LOGS > pairs.txt
"""
import argparse
import os

import numpy as np
import pandas as pd

READ_TYPES = {"pe": "", "se": "_se", "pb": "_pb", "ont": "_ont"}
COLUMNS = ["meta_sample", "meta_read_length", "meta_scenario", "taxon", "truth", "set", "p", "knob", "call"]
RESAMPLES = 1000


def calls(folder, read_type, kind):
    """The calls file of a read type's model at the database's gate; kind "" (forest in v29, gbm in v27) or "_gbm"."""
    path = os.path.join(folder, f"trained_model{READ_TYPES[read_type]}{kind}.calls.tsv.gz")
    df = pd.read_csv(path, sep="\t", usecols=COLUMNS, dtype={"meta_scenario": str}, keep_default_na=False)
    df["p"] = df["p"].astype(float)
    return df


def groups(df):
    """(label, row mask) of the sets reported: species held out (training rows), the test set's design samples and
    its scenario samples (hold-out), and the soil scenarios' hold-out samples alone."""
    train, test, scen = df["set"] == "training", df["set"] == "test", df["meta_scenario"] != ""
    out = [("held out, design", train & ~scen), ("held out, scenarios (hold-in)", train & scen),
           ("test, design", test & ~scen), ("test, scenarios (hold-out)", test & scen)]
    for s in ("soil", "soil_shallow", "moderate", "gut"):
        m = test & (df["meta_scenario"] == s)
        if m.any():
            out.append((f"test, {s}", m))
    return out


def per_sample(df, call):
    """[samples x (TP, FP, FN)] of the rows given, with the sample names."""
    y = df["truth"].to_numpy() == 1
    t = pd.DataFrame({"s": df["meta_sample"].to_numpy(), "tp": call & y, "fp": call & ~y, "fn": ~call & y})
    g = t.groupby("s", sort=True)[["tp", "fp", "fn"]].sum()
    return g.index.to_numpy(), g.to_numpy(dtype=float)


def f1(c):
    """F1 of summed counts (..., 3)."""
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


def tp_at_fp(y, p, fp_max):
    """The most true positives at a threshold with at most fp_max false positives (rows called at p >= threshold)."""
    order = np.argsort(-p, kind="stable")
    ys, ps = y[order], p[order]
    tp, fp = np.cumsum(ys), np.cumsum(~ys)
    # a threshold calls every row with that score or more: only the last row of a run of equal scores is a cut
    last = np.append(ps[1:] != ps[:-1], True)
    ok = last & (fp <= fp_max)
    return int(tp[ok].max()) if ok.any() else 0


def fmt(x, nd=4):
    return f"{x:.{nd}f}"


def forest_against_gbm(v29):
    print("## 1. forest against gbm (v29, the same rows; each kind at its own knob, as protal calls)")
    print("bootstrap: 1,000 resamples of the samples; P(gbm better): the share of resamples where gbm's F1 is higher.")
    print("matched FP: the true positives each kind reaches at the other's false-positive count (its best threshold).")
    print()
    head = ["read type", "set", "samples", "present", "knobs f/g", "forest TP/FP/FN", "gbm TP/FP/FN", "F1 forest",
            "F1 gbm", "gbm - forest", "95% interval", "P(gbm better)", "TP at gbm's FP: f/g", "TP at forest's FP: f/g"]
    print("\t".join(head))
    for rt in READ_TYPES:
        fo, gb = calls(v29, rt, ""), calls(v29, rt, "_gbm")
        key = ["meta_sample", "taxon", "set"]
        m = fo.merge(gb[key + ["p", "knob", "call"]], on=key, suffixes=("_f", "_g"), validate="one_to_one")
        assert len(m) == len(fo) == len(gb), (rt, len(m), len(fo), len(gb))
        for label, mask in groups(m):
            d = m[mask]
            y = d["truth"].to_numpy() == 1
            cf, cg = d["call_f"].to_numpy() == 1, d["call_g"].to_numpy() == 1
            sa, a = per_sample(d, cf)
            sb, b = per_sample(d, cg)
            assert (sa == sb).all()
            fa, fb, diff, lo, hi, better = paired(a, b)
            ta, tb = a.sum(0).astype(int), b.sum(0).astype(int)
            pf, pg = d["p_f"].to_numpy(), d["p_g"].to_numpy()
            at_g = f"{tp_at_fp(y, pf, tb[1])}/{tb[0]}"
            at_f = f"{ta[0]}/{tp_at_fp(y, pg, ta[1])}"
            print("\t".join([rt, label, str(len(sa)), str(int(y.sum())),
                             f"{d['knob_f'].iloc[0]:g}/{d['knob_g'].iloc[0]:g}", "/".join(map(str, ta)),
                             "/".join(map(str, tb)), fmt(fa), fmt(fb), f"{diff:+.4f}", f"[{lo:+.4f}, {hi:+.4f}]",
                             f"{better:.3f}", at_g, at_f]))
    print()


SAME_SIMULATION = 0.95  # a shared sample name is the same simulation if this share of its present taxa is in both


def against_v27(v29, v27):
    print("## 2. v29 against v27 on the test samples both builds simulated alike")
    print(f"alike: a sample name of both builds whose present taxa (rows with truth 1) overlap by {SAME_SIMULATION:.0%} "
          "of their union or more. Its present taxa are the union: one without a row in a build is a false negative "
          "there (a gate's lost species). Each model at its own knob, as protal calls. v27's models are gbm (its "
          "default then); v29's forest is its database's default.")
    print()
    head = ["read type", "set", "samples (alike / shared names)", "present (union)", "rowless v27/v29",
            "v27 gbm TP/FP/FN", "v29 gbm TP/FP/FN", "v29 forest TP/FP/FN", "F1 v27 gbm", "F1 v29 gbm",
            "F1 v29 forest", "v29 gbm - v27", "95% interval", "P(v29 gbm better)", "v29 forest - v27", "95% interval",
            "P(v29 forest better)"]
    print("\t".join(head))
    for rt in READ_TYPES:
        old = calls(v27, rt, "")
        old = old[old["set"] == "test"]
        new_g, new_f = calls(v29, rt, "_gbm"), calls(v29, rt, "")
        new_g, new_f = new_g[new_g["set"] == "test"], new_f[new_f["set"] == "test"]

        def present(df):
            d = df[df["truth"] == 1]
            return {s: set(t.tolist()) for s, t in d.groupby("meta_sample")["taxon"]}

        po, pn = present(old), present(new_g)
        scen_of = dict(zip(old["meta_sample"], old["meta_scenario"] != ""))
        shared = sorted(set(po) & set(pn))
        alike = [s for s in shared if len(po[s] & pn[s]) >= SAME_SIMULATION * len(po[s] | pn[s])]
        for label, scen in (("test, design", False), ("test, scenarios (hold-out)", True)):
            keep = [s for s in alike if scen_of[s] == scen]
            n_shared = sum(1 for s in shared if scen_of[s] == scen)
            if not keep:
                print("\t".join([rt, label, f"0 / {n_shared}"] + ["-"] * (len(head) - 3)))
                continue
            counts, rowless = [], []
            for df, own in ((old, po), (new_g, pn), (new_f, pn)):
                d = df[df["meta_sample"].isin(keep)]
                s, c = per_sample(d, d["call"].to_numpy() == 1)
                assert list(s) == sorted(keep)
                lost = np.array([len((po[x] | pn[x]) - own[x]) for x in s], dtype=float)
                c[:, 2] += lost
                counts.append(c)
                rowless.append(int(lost.sum()))
            a, g, f = counts
            fa, fg, dg, lg, hg, bg = paired(a, g)
            _, ff, df_, lf, hf, bf = paired(a, f)
            tot = [c.sum(0).astype(int) for c in counts]
            print("\t".join([rt, label, f"{len(keep)} / {n_shared}", str(int(tot[0][0] + tot[0][2])),
                             f"{rowless[0]}/{rowless[1]}"] +
                            ["/".join(map(str, t)) for t in tot] +
                            [fmt(fa), fmt(fg), fmt(ff), f"{dg:+.4f}", f"[{lg:+.4f}, {hg:+.4f}]", f"{bg:.3f}",
                             f"{df_:+.4f}", f"[{lf:+.4f}, {hf:+.4f}]", f"{bf:.3f}"]))
    print()


def by_read_length(v29):
    print("## 3. the paired-end test set's design samples by read length (v29, as protal calls)")
    print()
    print("\t".join(["read length", "samples", "present", "forest TP/FP/FN", "gbm TP/FP/FN", "F1 forest", "F1 gbm",
                     "FP/sample forest", "FP/sample gbm"]))
    fo, gb = calls(v29, "pe", ""), calls(v29, "pe", "_gbm")
    for rl in sorted(fo["meta_read_length"].unique()):
        row = [str(rl)]
        cs = []
        for df in (fo, gb):
            d = df[(df["set"] == "test") & (df["meta_scenario"] == "") & (df["meta_read_length"] == rl)]
            s, c = per_sample(d, d["call"].to_numpy() == 1)
            cs.append((len(s), c.sum(0).astype(int)))
        (n, a), (_, b) = cs
        if n == 0:
            continue
        row += [str(n), str(a[0] + a[2]), "/".join(map(str, a)), "/".join(map(str, b)), fmt(f1(a.astype(float))),
                fmt(f1(b.astype(float))), f"{a[1] / n:.2f}", f"{b[1] / n:.2f}"]
        print("\t".join(row))
    print()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("v29", help="v29's model_logs")
    ap.add_argument("v27", help="v27's model_logs")
    opts = ap.parse_args()
    forest_against_gbm(opts.v29)
    against_v27(opts.v29, opts.v27)
    by_read_length(opts.v29)


if __name__ == "__main__":
    main()
