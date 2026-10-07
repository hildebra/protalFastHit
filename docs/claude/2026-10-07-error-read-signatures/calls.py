#!/usr/bin/env python3
"""The v15 build's rows with the scores protal calls on: the training tables with the species-held-out scores
(trained_model<sfx>.predictions.tsv.gz, p_species) and the test tables with the final model's (test_predictions, and
model_logs/...scenario_predictions for the scenarios' hold-out samples where the archive has them), keyed
<set>:<sample> as fragments.model_tables does; the knob from the error records' taxa table (the build's calls)."""
import glob
import os

import numpy as np
import pandas as pd

import fragments

SUFFIX = {"pe": "", "se": "_se", "pb": "_pb", "ont": "_ont"}
META = ["meta_sample", "meta_scenario", "meta_novel_level", "meta_relative_rank", "meta_novel_congener",
        "meta_rep_genome", "meta_insilico_strain", "meta_neighbour_rank", "meta_lineage_genus"]
NOT_FEATURES = {"truth", "prediction", "probability", "taxon", "taxon_name", "domain", "p", "set", "knob", "call",
                "p_rows", "p_samples", "p_species", "p_collection_model", "meta_read_pairs", "meta_read_length",
                "meta_novel_species"}


def knob_of(records, rt):
    t = pd.read_pickle(os.path.join(records, f"{rt}.taxa.pkl.gz"))
    return float(pd.to_numeric(t["knob"]).iloc[0])


def load(build, rt):
    """(rows, feature columns): every row of the training and test tables with its score p and its call."""
    sfx = SUFFIX[rt]
    t = fragments.model_tables(build, rt, None)
    t["taxon"] = t["taxon"].astype(str)
    scores = []
    p = pd.read_csv(os.path.join(build, f"trained_model{sfx}.predictions.tsv.gz"), sep="\t",
                    usecols=["meta_sample", "taxon", "p_species"])
    p["meta_sample"] = "training:" + p["meta_sample"].astype(str)
    scores.append(p.rename(columns={"p_species": "p"}))
    for path in [os.path.join(build, f"trained_model{sfx}.test_predictions.tsv.gz")] + glob.glob(
            os.path.join(build, "model_logs", f"trained_model{sfx}.scenario_predictions.tsv.gz")):
        if os.path.isfile(path):
            q = pd.read_csv(path, sep="\t", usecols=["meta_sample", "taxon", "p"])
            q["meta_sample"] = "test:" + q["meta_sample"].astype(str)
            scores.append(q)
    s = pd.concat(scores, ignore_index=True)
    s["taxon"] = s["taxon"].astype(str)
    s = s.drop_duplicates(["meta_sample", "taxon"])
    t = t.merge(s, on=["meta_sample", "taxon"], how="inner")
    t["truth"] = pd.to_numeric(t["truth"], errors="coerce").fillna(0).astype(int)
    features = [c for c in t.columns if not c.startswith("meta_") and c not in NOT_FEATURES]
    for c in features:
        t[c] = pd.to_numeric(t[c], errors="coerce")
    features = [c for c in features if t[c].notna().any()]
    return t, features


def classes(t):
    """present rep / present strain / present insilico / absent novel congener (meta_novel_level species, relative rank
    genus: a species of its genus in the sample that the database lacks) / absent near novel (other ranks) / absent."""
    present = t["truth"] == 1
    insilico = pd.to_numeric(t["meta_insilico_strain"], errors="coerce").fillna(0) == 1
    rep = pd.to_numeric(t["meta_rep_genome"], errors="coerce").fillna(0) == 1
    level = t["meta_novel_level"].fillna("").astype(str)
    rank = t["meta_relative_rank"].fillna("").astype(str)
    return np.select([present & rep, present & insilico, present, (level == "species") & (rank == "genus"), level != ""],
                     ["present rep", "present insilico", "present strain", "absent novel congener",
                      "absent near novel"], "absent")


def f1(tp, fp, fn):
    return 2 * tp / max(1, 2 * tp + fp + fn)
