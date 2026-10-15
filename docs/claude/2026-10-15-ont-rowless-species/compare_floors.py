#!/usr/bin/env python3
"""Paired comparison of small builds that differ in protal's ONT identity floor (docs/claude/2026-10-15-ont-rowless-
species): for each build, the trainer's calls (model_logs/trained_model_ont.calls.tsv.gz: training rows with species
held out, test rows by the final model, each at the model's knob) and the collector's unseen tables
(<scratch>/<training|test>/training_data_ont.unseen.tsv: present species with marker reads and no row).

Per set (training: species held out; test), F1 as protal calls (truth: every simulated species with a row) and with the
marker truth (present: the species with marker reads; a present species without a row counts as missed), each against
the first build with a bootstrap over samples; then the marker truth by sample depth (meta_design), without intervals.

Usage: compare_floors.py LABEL=BUILD_DIR:SCRATCH_DIR ... > compare.txt
"""
import os
import sys

import numpy as np
import pandas as pd

RESAMPLES = 1000


def load(build, scratch):
    c = pd.read_csv(os.path.join(build, "model_logs", "trained_model_ont.calls.tsv.gz"), sep="\t",
                    usecols=["meta_design", "meta_sample", "meta_marker_reads", "taxon", "truth", "set", "p", "knob", "call"],
                    dtype={"meta_marker_reads": str})
    c["markers"] = pd.to_numeric(c["meta_marker_reads"], errors="coerce")
    unseen = {}
    for s in ("training", "test"):
        u = pd.read_csv(os.path.join(scratch, s, "training_data_ont.unseen.tsv"), sep="\t")
        unseen[s] = u.groupby("meta_sample").size()
    return c, unseen


def counts(c, unseen, which, samples):
    """[samples x (TP, FP, FN)] as protal calls and with the marker truth."""
    d = c[c["set"] == which]
    y, call = d["truth"].to_numpy() == 1, d["call"].to_numpy() == 1
    # present with marker reads; a genome without gene positions (no count) stays present, as the trainer keeps it
    m = (d["markers"].isna() | (d["markers"] > 0)).to_numpy()
    t = pd.DataFrame({"s": d["meta_sample"], "tp": call & y, "fp": call & ~y, "fn": ~call & y,
                      "mtp": call & y & m, "mfn": ~call & y & m})
    g = t.groupby("s")[["tp", "fp", "fn", "mtp", "mfn"]].sum().reindex(samples, fill_value=0)
    u = unseen[which].reindex(samples, fill_value=0).to_numpy()
    plain = g[["tp", "fp", "fn"]].to_numpy(dtype=float)
    marker = np.stack([g["mtp"], g["fp"], g["mfn"] + u], axis=1).astype(float)
    return plain, marker, int(u.sum())


def f1(c):
    return 2 * c[..., 0] / np.maximum(2 * c[..., 0] + c[..., 1] + c[..., 2], 1)


def main():
    builds = []
    for arg in sys.argv[1:]:
        label, dirs = arg.split("=", 1)
        build, scratch = dirs.split(":", 1)
        builds.append((label, *load(build, scratch)))
    for which in ("training", "test"):
        samples = sorted(set.intersection(*(set(c.loc[c["set"] == which, "meta_sample"]) for _, c, _ in builds)))
        print(f"## {which} ({'species held out' if which == 'training' else 'final model'}): {len(samples)} samples")
        print("\t".join(["build", "knob", "TP/FP/FN as protal calls", "F1", "d F1", "95% interval", "marker truth TP/FP/FN",
                         "unseen", "F1 markers", "d F1 markers", "95% interval"]))
        base = None
        rng = np.random.default_rng(1)
        w = rng.multinomial(len(samples), np.full(len(samples), 1 / len(samples)), size=RESAMPLES).astype(float)
        for label, c, unseen in builds:
            plain, marker, n_unseen = counts(c, unseen, which, samples)
            knob = c.loc[c["set"] == which, "knob"].iloc[0]
            row = [label, f"{knob:g}", "/".join(str(int(x)) for x in plain.sum(0)), f"{f1(plain.sum(0)):.4f}"]
            if base is None:
                base = (plain, marker)
                row += ["", "", "/".join(str(int(x)) for x in marker.sum(0)), str(n_unseen), f"{f1(marker.sum(0)):.4f}",
                        "", ""]
            else:
                dp = f1(w @ plain) - f1(w @ base[0])
                dm = f1(w @ marker) - f1(w @ base[1])
                row += [f"{f1(plain.sum(0)) - f1(base[0].sum(0)):+.4f}",
                        f"[{np.quantile(dp, 0.025):+.4f}, {np.quantile(dp, 0.975):+.4f}]",
                        "/".join(str(int(x)) for x in marker.sum(0)), str(n_unseen), f"{f1(marker.sum(0)):.4f}",
                        f"{f1(marker.sum(0)) - f1(base[1].sum(0)):+.4f}",
                        f"[{np.quantile(dm, 0.025):+.4f}, {np.quantile(dm, 0.975):+.4f}]"]
            print("\t".join(row))
        print()
        design = builds[0][1].drop_duplicates("meta_sample").set_index("meta_sample")["meta_design"].reindex(samples)
        depth = design.str.extract(r"_b(\d+)$", expand=False).astype(float).to_numpy()
        print(f"## {which} by depth: marker truth TP/FP/FN (unseen), F1")
        print("\t".join(["depth", "samples"] + [label for label, _, _ in builds]))
        for d in sorted(set(depth)):
            pick = depth == d
            row = [f"{d / 1e6:g} Mb", str(int(pick.sum()))]
            for label, c, unseen in builds:
                _, marker, _ = counts(c, unseen, which, samples)
                u = unseen[which].reindex(samples, fill_value=0).to_numpy()[pick].sum()
                m = marker[pick].sum(0)
                row.append(f"{'/'.join(str(int(x)) for x in m)} ({u}) {f1(m):.4f}")
            print("\t".join(row))
        print()


if __name__ == "__main__":
    main()
