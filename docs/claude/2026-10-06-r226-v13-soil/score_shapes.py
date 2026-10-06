#!/usr/bin/env python3
"""The presence models' scores per environment, and cuts other than the knob (soil_errors.py --save rows).

Figures: score_histograms.png, the scores of present and absent taxa (counts on a log axis) per read type and
environment (design test set, soil and shallow soil, hold-in with species held out and hold-out together), the knob
marked; score_vs_distance.png, the score against the reads' distance from the reference (excess_scaled_median) for pe
and pb in soil and shallow soil.

Cuts compared on each environment's hold-out samples (soil, soil_shallow) and the design test set:
  knob              the model's knob (as protal calls)
  valley            per sample, unsupervised: the lowest point of the smoothed density of logit(score) between its
                    low and high modes (where the sorted scores' S curve is steepest)
  mixture           per sample, unsupervised: a two-component Gaussian mixture on logit(score), cut where the two
                    posteriors are equal
  oracle            the best single threshold for the set itself (a ceiling for any one cut per environment)
  by distance       one threshold per bin of excess_scaled_median, learned on the soil and shallow-soil hold-in
                    samples (species held out) by coordinate ascent on their F1, applied to the hold-out
  by fragments      the same per bin of fragments
  by distance x fragments   the same per cell of both
  ... oracle        the per-bin thresholds learned on the set itself (ceilings)

    python3 score_shapes.py --rows v13_rows.tsv.gz --out .
"""
import argparse
import os

import numpy as np
import pandas as pd

BLUE, ORANGE, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#1f1f1e", "#6b6a63", "#e4e3dc"
DIST_BINS = [-1, 0.0025, 0.005, 0.01, 0.015, 0.02, 0.03, 1]
FRAG_BINS = [0, 2, 3, 5, 10, 30, 1e12]
GRID_T = np.round(np.arange(0.05, 0.96, 0.025), 3)


def logit(p):
    p = np.clip(p, 1e-4, 1 - 1e-4)
    return np.log(p / (1 - p))


def f1(y, call):
    tp, fp, fn = (call & (y == 1)).sum(), (call & (y == 0)).sum(), (~call & (y == 1)).sum()
    return 2 * tp / (2 * tp + fp + fn) if tp else 0.0


def valley_cut(p):
    """The density minimum of logit(p) between its low and high modes (smoothed histogram)."""
    x = logit(p)
    h, e = np.histogram(x, bins=np.linspace(-9.3, 9.3, 94))
    k = np.exp(-0.5 * (np.arange(-6, 7) / 2.0) ** 2)
    s = np.convolve(np.log1p(h), k / k.sum(), mode="same")
    mid = (e[:-1] + e[1:]) / 2
    lo = np.argmax(np.where(mid < 0, s, -np.inf))
    hi = np.argmax(np.where(mid > 0, s, -np.inf))
    if hi <= lo + 1:
        return 0.5
    v = lo + np.argmin(s[lo:hi + 1])
    return 1 / (1 + np.exp(-mid[v]))


def mixture_cut(p):
    """Equal posteriors of a two-component Gaussian mixture on logit(p) (EM, 200 steps)."""
    x = logit(p)
    m = np.array([np.percentile(x, 25), np.percentile(x, 90)])
    s = np.array([x.std(), x.std()]) / 2 + 1e-3
    w = np.array([0.5, 0.5])
    for _ in range(200):
        d = np.stack([w[j] / s[j] * np.exp(-0.5 * ((x - m[j]) / s[j]) ** 2) for j in range(2)])
        r = d / d.sum(0).clip(min=1e-300)
        n = r.sum(1)
        w, m = n / len(x), (r * x).sum(1) / n
        s = np.sqrt((r * (x - m[:, None]) ** 2).sum(1) / n) + 1e-3
    grid = np.linspace(m.min(), m.max(), 400)
    d0 = w[0] / s[0] * np.exp(-0.5 * ((grid - m[0]) / s[0]) ** 2)
    d1 = w[1] / s[1] * np.exp(-0.5 * ((grid - m[1]) / s[1]) ** 2)
    hi = int(np.argmax(m))
    above = (d1 if hi == 1 else d0) >= (d0 if hi == 1 else d1)
    t = grid[np.argmax(above)] if above.any() else 0.0
    return 1 / (1 + np.exp(-t))


def bin_thresholds(y, p, b, n_bins, start=0.5, passes=4):
    """One threshold per bin maximizing the F1 of all rows together (coordinate ascent)."""
    t = np.full(n_bins, start)
    for _ in range(passes):
        for k in range(n_bins):
            at = b == k
            if not at.any():
                continue
            best, best_f = t[k], -1
            for g in GRID_T:
                t[k] = g
                f = f1(y, p >= t[b])
                if f > best_f:
                    best, best_f = g, f
            t[k] = best
    return t


def bins_of(df, kind):
    d = np.clip(np.digitize(df.excess_scaled_median.fillna(0).to_numpy(), DIST_BINS[1:-1]), 0, len(DIST_BINS) - 2)
    f = np.clip(np.digitize(df.fragments.fillna(0).to_numpy(), FRAG_BINS[1:-1]), 0, len(FRAG_BINS) - 2)
    nd, nf = len(DIST_BINS) - 1, len(FRAG_BINS) - 1
    return {"distance": (d, nd), "fragments": (f, nf), "distance x fragments": (d * nf + f, nd * nf)}[kind]


def figures(df, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                         "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False})
    envs = [("(design)", "design test set"), ("soil", "soil"), ("soil_shallow", "shallow soil")]
    rts = ["pe", "se", "pb", "ont"]
    fig, axes = plt.subplots(len(rts), len(envs), figsize=(11, 10), sharex=True)
    bins = np.linspace(0, 1, 41)
    for i, rt in enumerate(rts):
        for j, (env, title) in enumerate(envs):
            ax = axes[i, j]
            g = df[(df.read_type == rt) & (df.scenario == env)]
            ax.grid(axis="y", color=GRID, linewidth=0.6)
            if not len(g):
                ax.text(0.5, 0.5, "no samples", ha="center", va="center", color=MUTED, transform=ax.transAxes)
                ax.set_yticks([])
                continue
            for truth, colour, label in ((0, ORANGE, "absent"), (1, BLUE, "present")):
                h, _ = np.histogram(g.p[g.truth == truth], bins=bins)
                ax.stairs(np.maximum(h, 0.8), bins, color=colour, linewidth=2, label=label, baseline=0.8)
            k = g.knob.iloc[0]
            ax.axvline(k, color=INK, linewidth=1, linestyle=(0, (3, 2)))
            ax.set_yscale("log")
            ax.set_ylim(0.8, None)
            yy, cc = g.truth.to_numpy(), (g.p >= k).to_numpy()
            ax.set_title(f"{rt}, {title}: F1 {f1(yy, cc):.3f} at knob {k:g}", color=INK, fontsize=9, loc="left")
            if i == len(rts) - 1:
                ax.set_xlabel("model score")
            if j == 0:
                ax.set_ylabel("taxa")
    axes[0, 0].legend(frameon=False, loc="upper center")
    fig.suptitle("Scores of present and absent taxa (hold-in with species held out and hold-out together; dashed: knob)",
                 color=INK, fontsize=10, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "score_histograms.png"), dpi=130, facecolor="white")
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True, sharey=True)
    rng = np.random.default_rng(1)
    for i, rt in enumerate(["pe", "pb"]):
        for j, (env, title) in enumerate(envs[1:]):
            ax = axes[i, j]
            g = df[(df.read_type == rt) & (df.scenario == env) & (df.fragments >= 3)]
            x = g.excess_scaled_median.clip(-0.005, 0.06).to_numpy()
            for truth, colour, label, size in ((0, ORANGE, "absent", 3), (1, BLUE, "present", 3)):
                at = np.flatnonzero(g.truth.to_numpy() == truth)
                at = rng.choice(at, min(len(at), 6000), replace=False)
                ax.scatter(x[at], g.p.to_numpy()[at], s=size, color=colour, alpha=0.35, linewidths=0, label=label)
            ax.axhline(g.knob.iloc[0], color=INK, linewidth=1, linestyle=(0, (3, 2)))
            ax.grid(color=GRID, linewidth=0.6)
            ax.set_title(f"{rt}, {title} (>= 3 fragments; at most 6,000 points per class)", color=INK, fontsize=9, loc="left")
            if i == 1:
                ax.set_xlabel("distance from the reference (excess_scaled_median, clipped at 0.06)")
            if j == 0:
                ax.set_ylabel("model score")
    leg = axes[0, 0].legend(frameon=False, loc="center right", markerscale=4)
    for h in leg.legend_handles:
        h.set_alpha(1)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "score_vs_distance.png"), dpi=130, facecolor="white")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rows", required=True)
    ap.add_argument("--out", required=True)
    opts = ap.parse_args()
    df = pd.read_csv(opts.rows, sep="\t", low_memory=False)
    figures(df, opts.out)
    out = []
    for rt in ["pe", "se", "pb", "ont"]:
        r = df[df.read_type == rt]
        learn = r[r.scenario.isin(["soil", "soil_shallow"]) & (r.set == "hold-in species held out")]
        learned = {}
        for kind in ("distance", "fragments", "distance x fragments"):
            b, n = bins_of(learn, kind)
            learned[kind] = bin_thresholds(learn.truth.to_numpy(), learn.p.to_numpy(), b, n, start=learn.knob.iloc[0])
        for env, st in (("soil", "hold-out"), ("soil_shallow", "hold-out"), ("(design)", "design test")):
            g = r[(r.scenario == env) & (r.set == st)]
            if not len(g):
                continue
            y, p = g.truth.to_numpy(), g.p.to_numpy()
            row = {"read type": rt, "set": f"{env} {st}", "knob": f1(y, p >= g.knob.iloc[0])}
            for name, cut in (("valley", valley_cut), ("mixture", mixture_cut)):
                call = np.zeros(len(g), dtype=bool)
                cuts = []
                for s, idx in g.groupby("meta_sample").indices.items():
                    c = cut(p[idx])
                    cuts.append(c)
                    call[idx] = p[idx] >= c
                row[name] = f1(y, call)
                row[name + " cuts"] = ",".join(f"{c:.2f}" for c in cuts[:4])
            row["oracle"] = max(f1(y, p >= t) for t in GRID_T)
            for kind, t in learned.items():
                b, n = bins_of(g, kind)
                row["by " + kind] = f1(y, p >= t[b])
                row["by " + kind + " oracle"] = f1(y, p >= bin_thresholds(y, p, b, n, start=g.knob.iloc[0])[b])
            out.append(row)
        print(f"{rt}: distance-bin thresholds learned on soil hold-in: "
              + ", ".join(f"{lo:g}-{hi:g} {t:.3f}" for lo, hi, t in zip(DIST_BINS[:-1], DIST_BINS[1:], learned["distance"])),
              flush=True)
        print(f"{rt}: fragment-bin thresholds: "
              + ", ".join(f"{lo:g}-{hi:g} {t:.3f}" for lo, hi, t in zip(FRAG_BINS[:-1], FRAG_BINS[1:], learned["fragments"])),
              flush=True)
    out = pd.DataFrame(out)
    pd.set_option("display.width", 300)
    pd.set_option("display.max_columns", 30)
    print(out.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    out.to_csv(os.path.join(opts.out, "score_cuts.tsv"), sep="\t", index=False, float_format="%.5g")


if __name__ == "__main__":
    main()
