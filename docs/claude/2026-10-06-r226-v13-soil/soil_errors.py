#!/usr/bin/env python3
"""Anatomy of the presence models' errors in the soil scenarios of an r226 build (soil, soil_shallow).

Reads a build's training and test tables (OUT/training, OUT/test) and its model_logs predictions: the hold-in samples
scored with species held out (trained_model*.predictions.tsv.gz, p_species, row for row with the training table) and
the hold-out samples scored by the final model (trained_model*.scenario_predictions.tsv.gz, row for row with the test
table's scenario rows). Calls at the model's knob (metrics.json global_knob, else 0.5), as protal calls.

Prints, per read type and scenario, the false negatives and false positives by class: fragments, strain kind
(representative, real strain, in-silico strain), the closest other species in the sample, the rank the source of an
absent taxon's reads was held out at, identity, low-MAPQ share and the congener features; the design's test set beside
it for contrast. --save writes the joined rows (meta columns, features used here, p, call) for later scripts.
"""

import argparse
import gzip
import json
import os

import numpy as np
import pandas as pd

READ_TYPES = {"pe": "", "se": "_se", "pb": "_pb", "ont": "_ont"}
FEATURES = ["fragments", "fragments_all", "em_fragments", "identity", "top_identity", "low_identity_share",
            "low_mapq_share", "congener_fit_share", "other_genus_fit_share", "em_own_share", "em_kept_own_share",
            "hit_gene_fraction", "gene_presence_ratio", "depth", "excess_median", "excess_high_share",
            "excess_scaled_median", "relative_skew", "relative_close_share", "genus_share", "genus_spill",
            "genus_top_fragments", "variant_sites_per_kb", "multiallelic_sites_per_kb", "sample_log_fragments",
            "failed_candidate_rate", "cluster_ani_radius", "cluster_min_ani", "cluster_genomes_log10", "lu_per_kb",
            "lsu_per_kb", "adjacent_support", "adjacent_unlikely_share", "conserved_fast_depth_ratio"]
META = ["meta_design", "meta_sample", "meta_read_pairs", "meta_rep_genome", "meta_insilico_strain", "meta_novel_level",
        "meta_relative_rank", "meta_neighbour_rank", "meta_novel_congener", "meta_scenario", "taxon", "taxon_name",
        "truth"]
FRAG_BINS = [0, 1, 2, 3, 5, 10, 30, 100, 1000, np.inf]
FRAG_LABELS = ["<1", "1", "2", "3-4", "5-9", "10-29", "30-99", "100-999", ">=1000"]


def read_table(path):
    head = pd.read_csv(path, sep="\t", nrows=0).columns
    cols = [c for c in META + FEATURES if c in head]
    return pd.read_csv(path, sep="\t", usecols=cols, low_memory=False,
                       dtype={"meta_scenario": str, "meta_novel_level": str, "meta_relative_rank": str,
                              "meta_neighbour_rank": str})


def knob_of(build, suffix):
    with open(os.path.join(build, "model_logs", f"trained_model{suffix}.metrics.json")) as f:
        m = json.load(f)
    k = (m.get("depth_knobs") or {}).get("global_knob")
    return float(k) if k is not None else 0.5


def load(build, rt):
    suffix = READ_TYPES[rt]
    name = "training_data.tsv" if rt == "pe" else f"training_data_{rt}.tsv"
    train = read_table(os.path.join(build, "training", name))
    pred = pd.read_csv(os.path.join(build, "model_logs", f"trained_model{suffix}.predictions.tsv.gz"), sep="\t",
                       usecols=["meta_sample", "taxon", "p_species"], low_memory=False)
    assert len(pred) == len(train) and (pred.taxon.to_numpy() == train.taxon.to_numpy()).all(), "predictions misaligned"
    train["p"] = pred.p_species.to_numpy()
    train["set"] = np.where(train.meta_scenario.fillna("") == "", "design species held out", "hold-in species held out")
    test = read_table(os.path.join(build, "test", name))
    sc = test.meta_scenario.fillna("") != ""
    design = test[~sc].copy()
    tp = pd.read_csv(os.path.join(build, "model_logs", f"trained_model{suffix}.test_predictions.tsv.gz"), sep="\t",
                     usecols=["taxon", "p"])
    assert len(tp) == len(design) and (tp.taxon.to_numpy() == design.taxon.to_numpy()).all(), "test misaligned"
    design["p"] = tp.p.to_numpy()
    design["set"] = "design test"
    held = test[sc].copy()
    sp = pd.read_csv(os.path.join(build, "model_logs", f"trained_model{suffix}.scenario_predictions.tsv.gz"), sep="\t",
                     usecols=["taxon", "p"])
    assert len(sp) == len(held) and (sp.taxon.to_numpy() == held.taxon.to_numpy()).all(), "scenarios misaligned"
    held["p"] = sp.p.to_numpy()
    held["set"] = "hold-out"
    df = pd.concat([train, design, held], ignore_index=True)
    df["scenario"] = df.meta_scenario.fillna("").replace("", "(design)")
    df["knob"] = knob_of(build, suffix)
    df["call"] = df.p >= df.knob
    df["outcome"] = np.select([df.truth.eq(1) & df.call, df.truth.eq(0) & df.call, df.truth.eq(1) & ~df.call],
                              ["TP", "FP", "FN"], "TN")
    df["frag_bin"] = pd.cut(df.fragments.fillna(0), FRAG_BINS, right=False, labels=FRAG_LABELS)
    rep = pd.to_numeric(df.meta_rep_genome, errors="coerce")
    ins = pd.to_numeric(df.meta_insilico_strain, errors="coerce").fillna(0)
    df["strain"] = np.where(df.truth.eq(1), np.where(rep.eq(1), "representative",
                                                     np.where(ins.eq(1), "in-silico strain", "real strain")), "")
    df["read_type"] = rt
    return df


def f1(g):
    tp, fp, fn = (g.outcome == "TP").sum(), (g.outcome == "FP").sum(), (g.outcome == "FN").sum()
    return 2 * tp / (2 * tp + fp + fn) if tp else 0.0


def rates(g, by, kind):
    """FN (present rows) or FP (absent rows) by class: rows, errors, rate %, share of the errors %."""
    if kind == "FN":
        g = g[g.truth == 1]
    else:
        g = g[g.truth == 0]
    err = g.outcome == kind
    t = g.groupby(by, observed=True).agg(rows=("outcome", "size"), errors=("outcome", lambda s: (s == kind).sum()))
    t["rate %"] = (100 * t.errors / t.rows).round(2)
    t["share %"] = (100 * t.errors / max(1, err.sum())).round(1)
    return t


def describe(df, label):
    print(f"\n=== {label}: {df.meta_sample.nunique()} samples, F1 {f1(df):.4f}, "
          + ", ".join(f"{k} {(df.outcome == k).sum()}" for k in ("TP", "FP", "FN", "TN")))
    print("\nFN by fragments:")
    print(rates(df, "frag_bin", "FN").to_string())
    print("\nFN by strain kind:")
    print(rates(df, "strain", "FN").to_string())
    print("\nFN by strain kind and fragments (rate %):")
    pres = df[df.truth == 1]
    print(pres.pivot_table(index="frag_bin", columns="strain", values="outcome", observed=True,
                           aggfunc=lambda s: round(100 * (s == "FN").mean(), 1)).to_string())
    print("\nFN by the closest other species in the sample (meta_neighbour_rank):")
    print(rates(df, df.meta_neighbour_rank.fillna("?"), "FN").to_string())
    print("\nFP by fragments:")
    print(rates(df, "frag_bin", "FP").to_string())
    print("\nFP by the rank their reads' source was held out at (meta_novel_level; '' = closest species in the DB):")
    print(rates(df, df.meta_novel_level.fillna(""), "FP").to_string())
    print("\nFP by the deepest rank shared with a simulated species (meta_relative_rank):")
    print(rates(df, df.meta_relative_rank.fillna("?"), "FP").to_string())
    cols = ["fragments", "identity", "top_identity", "low_mapq_share", "congener_fit_share", "em_own_share",
            "hit_gene_fraction", "excess_median", "relative_close_share", "genus_share", "cluster_ani_radius", "p"]
    cols = [c for c in cols if c in df]
    print("\nmedians by outcome:")
    print(df.groupby("outcome")[cols].median().round(4).to_string())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build", required=True, help="the build's local copy (training/, test/, model_logs/)")
    ap.add_argument("--read-types", default="pe,se,pb,ont")
    ap.add_argument("--save", help="write the joined rows of all read types here (tsv.gz)")
    opts = ap.parse_args()
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    frames = []
    for rt in opts.read_types.split(","):
        df = load(opts.build, rt)
        frames.append(df)
        print(f"\n\n######## {rt} (knob {df.knob.iloc[0]})")
        for scenario in ("soil", "soil_shallow"):
            for s in ("hold-in species held out", "hold-out"):
                g = df[(df.scenario == scenario) & (df.set == s)]
                if len(g):
                    describe(g, f"{rt} {scenario} {s}")
        describe(df[df.set == "design test"], f"{rt} design test")
    if opts.save:
        out = pd.concat(frames, ignore_index=True)
        out = out[out.scenario.isin(["soil", "soil_shallow"]) | (out.set == "design test")]
        out.to_csv(opts.save, sep="\t", index=False, float_format="%.6g")


if __name__ == "__main__":
    main()
