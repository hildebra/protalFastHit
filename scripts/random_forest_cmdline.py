#!/usr/bin/env python3
"""Train protal's presence model and report how well it does on samples and species it has not seen.

The model is a random forest that gives each taxon with reads the probability that it is present;
protal reports taxa whose probability is at least --knob. Training needs no Java: the forest is
written as PMML by model_pmml.py, and protal scores exactly what scikit-learn scores (checked here on
every training row).

    python3 scripts/random_forest_cmdline.py --truth-file training/training_data.tsv \\
        --output-prefix training/model

The training table comes from collect_training_data.py: protal's training dumps
(<profile>.truth_annotated, one row per taxon with reads) joined, with meta_* columns naming the
sample (meta_sample) and, if known, the taxon's domain (meta_domain).

Evaluation. Rows of one sample share its reads, and rows of one species share its reference, so a
random split of rows scores a model on samples and species it was trained on. Here each row is
scored by forests that saw neither its sample ("by sample") nor its species ("by species"), and,
with --taxonomy, by forests that saw no taxon of its genus, family, class or phylum. By species is
what matters for a large database: of GTDB's ~130,000 species, a training set holds a few thousand,
so most species protal meets in real samples were never in training; the clades tell how far that
holds. The out-of-bag estimate (each tree scores the rows it was not grown on) comes free with the
fit. When the training database lacked species or whole clades (build_gtdb_database.py), the report
also counts the false positives their reads cause, by the rank they were held out at.

Written to PREFIX.*:
  xml                 the model (protal --model FILE, or protal --add_model FILE --read_type pe: model_pe.xml of a database)
  report.txt          the evaluation (also printed); metrics.json has its numbers
  predictions.tsv.gz  each taxon's probabilities out of fold, with its main features
  thresholds.tsv      precision and sensitivity by threshold, from species held out
  varimp.tsv          feature importances
  joblib              the fitted scikit-learn forest

--evaluation full adds studies that tell whether the training set and the settings suffice: other
feature sets, the procedure this script used before (grid search over max_features, then only the
top features, 512 trees), forest sizes and tree counts, and fewer training samples.

--depth-knobs also chooses a knob per depth of the sample: for each bin of the samples' fragments over
all their taxa (the digits of their number less one, 2 to 6), the threshold with the highest F1 on
species held out, where the bin has more than 50 taxa and 10 present. They go into the model's header,
and protal applies the knob of each sample's bin unless --knob is given (other bins: --knob's default).
On the 0.7.1 pipeline's own long-read test sets this raised F1 by 0.007 (PacBio) and 0.015 (ONT)
(docs/claude/2026-10-01-f1-opportunities), but on the v0.7.1 benchmark's samples it cost PacBio up to 0.016
and did not help Nanopore: a bin's knob rests on few samples, and samples near a bin edge switch knobs
(docs/claude/2026-10-01-features-depth-knobs). build_gtdb_database.py passes it for PacBio and ONT, to be
tested again with the new database build.
"""

from __future__ import annotations

import argparse
import datetime
import json
import os
import platform
import sys
import time

import joblib
import numpy as np
import pandas as pd
import sklearn
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import average_precision_score, brier_score_loss, log_loss, precision_recall_curve, roc_auc_score
from sklearn.model_selection import GridSearchCV, GroupKFold, StratifiedKFold

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lineages  # noqa: E402
from model_features import feature_columns  # noqa: E402
from model_pmml import PmmlForest, format_depth_knobs, read_depth_knobs, write_forest  # noqa: E402

# Clades held out in cross-validation (with --taxonomy): each row is scored by forests that saw no taxon of
# its genus, family, order, class or phylum.
CLADE_SCHEMES = ("genus", "family", "order", "class", "phylum")
NOVEL_RANKS = ("species", "genus", "family", "order", "class", "phylum")

# Main features written next to the predictions, when the table has them.
DIAGNOSTIC_FEATURES = ["fragments", "depth", "hit_gene_fraction", "gene_presence_ratio", "identity", "top_identity",
                       "low_identity_share", "lsu_per_kb", "lu_per_kb", "uniqueness"]
FRAGMENT_BINS = [0, 10, 100, 1000, np.inf]
# --depth-knobs: the bins of a sample's fragments (profiler::DepthKnobBin), the taxa a bin needs for a knob of its own,
# and the thresholds tried.
DEPTH_KNOB_BINS = range(2, 7)
DEPTH_KNOB_MIN_TAXA = 50
DEPTH_KNOB_MIN_PRESENT = 10
DEPTH_KNOB_GRID = np.round(np.arange(0.05, 0.955, 0.01), 2)
# The features the report shows by class of taxon (study_feature_classes): what the conservation of the genes a taxon's
# reads hit, before and after the MAPQ filter, and the divergence beyond the base qualities say of a relative's reads.
CLASS_FEATURES = ["conserved_fast_record_ratio", "conserved_fast_kept_ratio", "conserved_fast_depth_ratio",
                  "conserved_hit_share", "excess_median"]
TAXON_CLASSES = ["present, no congener held out", "present, beside a held-out congener",
                 "absent, congener of a held-out species", "absent, other"]


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--truth-file", required=True, help="training table (collect_training_data.py) or one training dump")
    p.add_argument("--output-prefix", required=True)
    p.add_argument("--features", choices=["normalized", "all"], default="normalized",
                   help="normalized: the features that do not depend on database, domain, depth and read length "
                        "(model_features.py; default); all: every feature column of the table")
    p.add_argument("--reference-pmml", help="train on the inputs of this PMML model instead of --features")
    p.add_argument("--ntree", type=int, default=64, help="trees (default 64)")
    p.add_argument("--maxnodes", type=int, default=128, help="leaves per tree at most, 0 for no limit (default 128)")
    p.add_argument("--min-samples-leaf", type=int, default=1)
    p.add_argument("--max-features", default="sqrt", help="features tried per split: sqrt, log2, a count or a fraction")
    p.add_argument("--knob", type=float, default=0.5, help="the threshold protal will use (its --knob, default 0.5)")
    p.add_argument("--depth-knobs", action="store_true",
                   help="also choose a knob per depth bin of the sample, on species held out, and store them in the "
                        "model, which protal then applies unless --knob is given (see above)")
    p.add_argument("--folds", type=int, default=5)
    p.add_argument("--evaluation", choices=["full", "basic", "none"], default="full",
                   help="full: also the studies (see above); basic: out of bag, by sample and by species; "
                        "none: fit and export only")
    p.add_argument("--taxonomy", help="internal_taxonomy.dmp of the database, for the taxa's domains when the "
                                      "table has no meta_domain column")
    p.add_argument("--test-file", help="an independent test table (collect_training_data.py with another design and "
                                       "seed): scored with the fitted forest and reported, by depth and by rank")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--threads", type=int, default=4)
    return p.parse_args(argv)


def max_features_value(text):
    if text in ("sqrt", "log2"):
        return text
    value = float(text)
    return int(value) if value >= 1 and value.is_integer() else value


# ---- data -------------------------------------------------------------------------------------------------

def as_truth(values):
    if values.dtype == bool:
        return values.astype(int)
    if np.issubdtype(values.dtype, np.number):
        return (values.fillna(0) > 0).astype(int)
    return values.fillna("").astype(str).str.lower().isin({"true", "t", "1", "yes", "y"}).astype(int)


def domains_from_taxonomy(path, taxa):
    tax = pd.read_csv(path, sep="\t", dtype=str)
    parent = dict(zip(tax["id"], tax["parent_id"]))
    name = dict(zip(tax["id"], tax["name"]))

    def domain(taxon):
        node = str(taxon)
        for _ in range(64):
            if str(name.get(node, "")).startswith("d__"):
                return name[node][3:]
            if node not in parent or parent[node] == node:
                break
            node = parent[node]
        return "unknown"

    return taxa.map({t: domain(t) for t in taxa.unique()})


def load_table(path, taxonomy=None):
    # round_trip: read each value back as the double protal gave the model (the default parser can be off by an ulp).
    df = pd.read_csv(path, sep="\t", float_precision="round_trip", low_memory=False)
    if "truth" not in df.columns:
        sys.exit(f"{path} has no truth column: is it a training dump (<profile>.truth_annotated)?")
    df["truth"] = as_truth(df["truth"])
    df["domain"] = df["meta_domain"].fillna("unknown").astype(str) if "meta_domain" in df.columns else "unknown"
    if taxonomy:
        unknown = df["domain"] == "unknown"
        df.loc[unknown, "domain"] = domains_from_taxonomy(taxonomy, df.loc[unknown, "taxon"])
        # meta_lineage_<rank>: the taxon's clades, for holding out whole clades (not features: meta_*).
        by_id, _ = lineages.from_taxonomy(taxonomy)
        taxa = df["taxon"].astype(str)
        for rank in CLADE_SCHEMES:
            names = {t: by_id.get(t, {}).get(rank, "unknown") for t in taxa.unique()}
            df["meta_lineage_" + rank] = taxa.map(names)
    return df


def check_features(df, cols):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        sys.exit("the training table lacks features: " + ", ".join(missing))
    text = [c for c in cols if not np.issubdtype(df[c].dtype, np.number)]
    if text:
        sys.exit("features that are not numbers: " + ", ".join(text))
    bad = {c: int((~np.isfinite(df[c])).sum()) for c in cols if (~np.isfinite(df[c])).any()}
    if bad:
        # protal gives the model the values as they are in the dump, so a model trained on replaced values
        # would see other values in protal than in training.
        sys.exit("features with values that are not finite numbers (the dump's values are what protal gives the "
                 "model, so they cannot be replaced for training): " + ", ".join(f"{c} ({n} rows)" for c, n in bad.items()))


# ---- models and folds -------------------------------------------------------------------------------------

def forest_params(opts, **overrides):
    params = dict(n_estimators=opts.ntree, max_leaf_nodes=opts.maxnodes or None, min_samples_leaf=opts.min_samples_leaf,
                  max_features=max_features_value(opts.max_features), class_weight="balanced",
                  random_state=opts.seed, n_jobs=opts.threads)
    params.update(overrides)
    return params


def folds(df, y, scheme, opts):
    """(train, test) index pairs; None if the table cannot be split that way."""
    if scheme == "rows":
        n = min(opts.folds, int(np.bincount(y).min()))
        return list(StratifiedKFold(n, shuffle=True, random_state=opts.seed).split(df, y)) if n >= 2 else None
    column = {"samples": "meta_sample", "species": "taxon"}.get(scheme, "meta_lineage_" + scheme)
    if column not in df.columns:
        return None
    groups = df[column].astype(str).to_numpy()
    n = min(opts.folds, len(np.unique(groups)))
    if n < 2:
        return None
    return list(GroupKFold(n).split(df, y, groups))


def predict_out_of_fold(X, y, splits, params, fit_rows=None, leaves=None):
    """Each row's probability from the forest that did not see its fold. fit_rows(train) may thin a fold's
    training rows; leaves, a list, gets each fold forest's mean leaves per tree."""
    p = np.full(len(y), np.nan)
    for train, test in splits:
        if fit_rows is not None:
            train = fit_rows(train)
        if len(np.unique(y[train])) < 2:
            continue
        rf = RandomForestClassifier(**params).fit(X[train], y[train])
        p[test] = rf.predict_proba(X[test])[:, 1]
        if leaves is not None:
            leaves.append(np.mean([e.tree_.n_leaves for e in rf.estimators_]))
    return p


# ---- metrics ----------------------------------------------------------------------------------------------

def metrics(y, p, df, knob):
    ok = ~np.isnan(p)
    y, p = y[ok], p[ok]
    sub = df[ok]
    call = p >= knob
    tp, fp = int((call & (y == 1)).sum()), int((call & (y == 0)).sum())
    fn = int((~call & (y == 1)).sum())
    both = len(np.unique(y)) == 2
    m = {
        "taxa": int(len(y)), "present": int(y.sum()),
        "AUC": float(roc_auc_score(y, p)) if both else None,
        "AP": float(average_precision_score(y, p)) if both else None,
        "log_loss": float(log_loss(y, np.clip(p, 1e-6, 1 - 1e-6), labels=[0, 1])),
        "brier": float(brier_score_loss(y, p)),
        "sensitivity": tp / (tp + fn) if tp + fn else None,
        "precision": tp / (tp + fp) if tp + fp else None,
        "F1": 2 * tp / (2 * tp + fp + fn) if tp + fp + fn else None,
        "FP": fp, "FN": fn,
    }
    if "meta_sample" in sub.columns:
        m["FP_per_sample"] = fp / max(1, sub["meta_sample"].nunique())
    for domain in sorted(sub["domain"].unique()):
        present = (sub["domain"].to_numpy() == domain) & (y == 1)
        if present.any():
            m[f"sensitivity_{domain}"] = float(call[present].mean())
    return m


def best_threshold(y, p):
    ok = ~np.isnan(p)
    precision, recall, thresholds = precision_recall_curve(y[ok], p[ok])
    precision, recall = precision[:-1], recall[:-1]
    f1 = np.where(precision + recall > 0, 2 * precision * recall / np.maximum(precision + recall, 1e-300), 0)
    i = int(np.argmax(f1))
    return float(thresholds[i]), float(precision[i]), float(recall[i]), float(f1[i])


def threshold_table(y, p):
    ok = ~np.isnan(p)
    y, p = y[ok], p[ok]
    rows = []
    for t in np.round(np.arange(0.0, 1.0001, 0.01), 2):
        call = p >= t
        tp, fp, fn = int((call & (y == 1)).sum()), int((call & (y == 0)).sum()), int((~call & (y == 1)).sum())
        rows.append({"threshold": t, "precision": tp / (tp + fp) if tp + fp else 1.0,
                     "sensitivity": tp / (tp + fn) if tp + fn else 0.0,
                     "F1": 2 * tp / (2 * tp + fp + fn) if tp + fp + fn else 0.0, "FP": fp, "FN": fn})
    return pd.DataFrame(rows)


# ---- report -----------------------------------------------------------------------------------------------

class Report:
    def __init__(self):
        self.lines, self.data = [], {}

    def section(self, title):
        self.add("")
        self.add("## " + title)

    def add(self, text=""):
        self.lines.append(text)
        print(text, flush=True)

    def table(self, frame):
        self.add(frame.to_string(index=False, float_format=lambda v: f"{v:.4f}", na_rep="-"))

    def write(self, prefix):
        with open(prefix + ".report.txt", "w") as fh:
            fh.write("\n".join(self.lines) + "\n")
        with open(prefix + ".metrics.json", "w") as fh:
            json.dump(self.data, fh, indent=1, default=_json_value)


def _json_value(value):
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return None if np.isnan(value) else float(value)
    return str(value)


def metrics_table(rows):
    """rows: list of (label dict, metrics dict) -> DataFrame with the usual columns."""
    keep = ["taxa", "present", "AUC", "AP", "log_loss", "F1", "sensitivity", "precision", "FP", "FN", "FP_per_sample"]
    out = []
    for label, m in rows:
        row = dict(label)
        row.update({k: m.get(k) for k in keep if k in m})
        row.update({k: v for k, v in m.items() if k.startswith("sensitivity_")})
        out.append(row)
    frame = pd.DataFrame(out)
    return frame.rename(columns={c: "sens_" + c[len("sensitivity_"):] for c in frame.columns if c.startswith("sensitivity_")})


def fmt(value, digits=3):
    return "-" if value is None or (isinstance(value, float) and np.isnan(value)) else f"{value:.{digits}f}"


# ---- the studies ------------------------------------------------------------------------------------------

def study_data(report, df, y, cols):
    report.section("Training data")
    samples = df["meta_sample"].nunique() if "meta_sample" in df.columns else 1
    present, absent = df[y == 1], df[y == 0]
    report.add(f"{len(df)} taxa in {samples} samples: {len(present)} present, {len(absent)} absent "
               f"({len(absent) / max(1, len(present)):.1f} absent per present); "
               f"{df['taxon'].nunique()} species ({present['taxon'].nunique()} ever present, "
               f"{absent['taxon'].nunique()} ever absent)")
    absent_count = ("truth", lambda s: int((s == 0).sum()))
    by_domain = df.groupby("domain").agg(present=("truth", "sum"), absent=absent_count, species=("taxon", "nunique"))
    report.table(by_domain.reset_index())
    if "meta_read_pairs" in df.columns and "meta_sample" in df.columns:
        by_depth = df.groupby("meta_read_pairs").agg(samples=("meta_sample", "nunique"), present=("truth", "sum"),
                                                     absent=absent_count)
        if "fragments" in df.columns:
            by_depth["median_fragments_present"] = present.groupby("meta_read_pairs")["fragments"].median()
            by_depth["median_fragments_absent"] = absent.groupby("meta_read_pairs")["fragments"].median()
        report.table(by_depth.reset_index())
    report.data["data"] = {"taxa": len(df), "samples": samples, "present": len(present), "absent": len(absent),
                           "species": int(df["taxon"].nunique()), "by_domain": by_domain.to_dict("index")}

    # Present taxa whose reads all match their reference closely: strains simulated from the database's own
    # genomes. Real strains differ from the representative by up to ~5%.
    quantiles = {}
    for c in ("identity", "top_identity"):
        if c in df.columns:
            q = [0.01, 0.05, 0.25, 0.5]
            quantiles[c] = {"present": present[c].quantile(q).tolist(), "absent": absent[c].quantile(q).tolist()}
            report.add(f"{c} quantiles 1/5/25/50%: present " + " ".join(fmt(v, 4) for v in quantiles[c]["present"])
                       + " | absent " + " ".join(fmt(v, 4) for v in quantiles[c]["absent"]))
    report.data["data"]["identity_quantiles"] = quantiles
    warnings = []
    # identity also counts relatives' reads, so it is lowered even where the strains are the references; a best
    # read identical to the reference in (nearly) every present taxon is what simulating the references gives.
    if "top_identity" in quantiles and quantiles["top_identity"]["present"][0] >= 0.999:
        warnings.append("99% of the present taxa have a read identical to their reference: the simulated genomes are "
                        "the references themselves (only representative genomes?). Real strains differ by up to a few "
                        "%, and the model may call them absent. Simulate from non-representative genomes as well "
                        "(for GTDB: from NCBI, by the accessions in GTDB's metadata).")
    # Features that differ between species but not within one (a feature constant everywhere names nothing).
    per_species = [c for c in cols if df[c].nunique() > 1 and df.groupby("taxon")[c].nunique().max() == 1]
    if per_species:
        report.add(f"features constant within every species ({len(per_species)}): {', '.join(per_species)}. "
                   "Together they can name a species, and the model can learn which species tend to be present.")
    report.data["data"]["features_constant_per_species"] = per_species
    if (df["domain"] == "unknown").all():
        report.add("no domains known (no meta_domain column, no --taxonomy): no sensitivity by domain")
    elif "Archaea" in df["domain"].values and (present["domain"] == "Archaea").sum() < 50:
        warnings.append(f"only {(present['domain'] == 'Archaea').sum()} present archaeal taxa: too few to judge archaea")
    return warnings


def study_evaluation(report, df, X, y, opts, oob):
    """Out of fold by rows, samples and species; returns the probabilities by scheme."""
    report.section(f"How well the model does on data it was not trained on (knob {opts.knob})")
    report.add("rows: random rows held out (the samples and species of the held-out rows are in training); "
               "samples: whole samples held out; species: whole species held out (as most GTDB species are when "
               "profiling real samples); genus, family, order, class, phylum (with --taxonomy): whole clades held out, "
               "so the forest saw no taxon of the row's clade (as for taxa of clades the training data barely "
               "cover). out of bag: each row scored by the trees not grown on it.")
    params = forest_params(opts)
    p = {"out of bag": oob}
    for scheme in ("rows", "samples", "species", *CLADE_SCHEMES):
        splits = folds(df, y, scheme, opts)
        if splits is not None:
            p[scheme] = predict_out_of_fold(X, y, splits, params)
    rows = [({"held out": k}, metrics(y, v, df, opts.knob)) for k, v in p.items()]
    if "probability" in df.columns:
        p["collection model"] = df["probability"].to_numpy(dtype=float)
        rows.append(({"held out": "collection model"}, metrics(y, p["collection model"], df, opts.knob)))
    report.table(metrics_table(rows))
    report.add("collection model: the model that profiled the training samples (their dumps' probability), "
               "e.g. the model shipped with protal.")
    report.data["evaluation"] = {label["held out"]: m for label, m in rows}
    return p


def study_breakdown(report, df, y, p, opts):
    report.section("Where the errors are (species held out; collection model for comparison)")
    new = p.get("species", p["out of bag"])
    old = p.get("collection model")
    call_new = new >= opts.knob
    call_old = old >= opts.knob if old is not None else None
    frame = df.assign(call_new=call_new, call_old=call_old if call_old is not None else np.nan)
    keys = ["domain"] + (["meta_read_pairs"] if "meta_read_pairs" in df.columns else [])
    rows = []
    for key, g in frame.groupby(keys):
        key = key if isinstance(key, tuple) else (key,)
        pr, ab = g[g.truth == 1], g[g.truth == 0]
        row = dict(zip(keys, key))
        row.update({"present": len(pr), "found": int(pr.call_new.sum()), "absent": len(ab), "FP": int(ab.call_new.sum())})
        if old is not None:
            row.update({"found_collection": int(pr.call_old.sum()), "FP_collection": int(ab.call_old.sum())})
        rows.append(row)
    report.table(pd.DataFrame(rows))
    if "fragments" in df.columns:
        report.add("by the taxon's fragments (read pairs) in the sample:")
        bins = pd.cut(df["fragments"], FRAGMENT_BINS, right=False, labels=["<10", "10-99", "100-999", ">=1000"])
        rows = []
        for b, g in frame.groupby(bins, observed=True):
            pr, ab = g[g.truth == 1], g[g.truth == 0]
            row = {"fragments": b, "present": len(pr), "found": int(pr.call_new.sum()), "absent": len(ab), "FP": int(ab.call_new.sum())}
            if old is not None:
                row.update({"found_collection": int(pr.call_old.sum()), "FP_collection": int(ab.call_old.sum())})
            rows.append(row)
        report.table(pd.DataFrame(rows))
    # The two things a harder simulation adds (collect_training_data.py --taxonomy, --novel_species): present
    # species simulated from other strains than the reference, and absent taxa that a species the database
    # lacks passes its reads to.
    groups = {}
    if "meta_rep_genome" in df.columns and df.loc[df.truth == 1, "meta_rep_genome"].notna().any():
        rep = df["meta_rep_genome"]
        groups["present, simulated from"] = [("the representative", (df.truth == 1) & (rep == 1)),
                                             ("another genome", (df.truth == 1) & (rep == 0))]
    if "meta_novel_congener" in df.columns and (df.loc[df.truth == 0, "meta_novel_congener"] == 1).any():
        congener = df["meta_novel_congener"] == 1
        groups["absent taxa"] = [("congeners of a species the database lacks", (df.truth == 0) & congener),
                                 ("others", (df.truth == 0) & ~congener)]
    for title, parts in groups.items():
        rows = []
        for label, mask in parts:
            g = frame[mask.to_numpy()]
            row = {title: label, "taxa": len(g), "called": int(g.call_new.sum())}
            if old is not None:
                row["called_collection"] = int(g.call_old.sum())
            rows.append(row)
        report.add(f"{title}:")
        report.table(pd.DataFrame(rows))
        report.data.setdefault("breakdown", {})[title] = rows
    shown = [c for c in ["meta_sample", "taxon_name", "domain", "fragments", "identity", "top_identity",
                         "hit_gene_fraction", "low_identity_share", "lsu_per_kb"] if c in df.columns]
    hard = df.assign(p=new)
    report.add("absent taxa scored highest:")
    report.table(hard[hard.truth == 0].nlargest(12, "p")[shown + ["p"]])
    report.add("present taxa scored lowest:")
    report.table(hard[hard.truth == 1].nsmallest(12, "p")[shown + ["p"]])


def rate(n, d):
    """n of d in %, None without d."""
    return 100 * n / d if d else None


def study_by_rank(report, df, y, p, opts, title="False positives and false negatives by taxonomic rank", key="by_rank",
                  clade_rows=True, scored="scored with species held out"):
    """False positive and false negative rates by taxonomic rank (collect_training_data.py: meta_novel_level,
    meta_novel_levels, meta_relative_rank, meta_neighbour_rank; the clades held out in cross-validation)."""
    ranks_in = lambda column: column in df.columns and (df[column].fillna("").astype(str) != "").any()
    if not (ranks_in("meta_novel_level") or ranks_in("meta_neighbour_rank") or any(k in p for k in CLADE_SCHEMES)):
        return
    report.section(title)
    report.add(f"Rates in % at the knob, {scored} (and by the collection model, _collection). "
               "FP rate: of the absent taxa, those called; FN rate: of the present taxa, those not called.")
    new = p["species"] if "species" in p else p["out of bag"]
    old = p.get("collection model")
    absent, present = y == 0, y == 1
    called_new = new >= opts.knob
    called_old = old >= opts.knob if old is not None else None
    data = {}

    def counts(mask, positive):
        """{taxa, errors, rate} for the rows of mask, errors: called (absent) or missed (present)."""
        err = lambda called: int((called[mask] if not positive else ~called[mask]).sum())
        row = {"taxa": int(mask.sum()), "FN" if positive else "FP": err(called_new)}
        row["FN rate" if positive else "FP rate"] = rate(row["FN" if positive else "FP"], row["taxa"])
        if called_old is not None:
            row[("FN" if positive else "FP") + "_collection"] = err(called_old)
            row[("FN" if positive else "FP") + " rate_collection"] = rate(err(called_old), row["taxa"])
        return row

    if ranks_in("meta_novel_level"):
        level = df["meta_novel_level"].fillna("").astype(str).to_numpy()
        simulated = {}
        if "meta_novel_levels" in df.columns and "meta_sample" in df.columns:
            for text in df.drop_duplicates("meta_sample")["meta_novel_levels"].fillna("").astype(str):
                for part in filter(None, text.split(",")):
                    rank, _, n = part.partition(":")
                    simulated[rank] = simulated.get(rank, 0) + int(n or 0)
        report.add("")
        report.add("False positives from what the training database lacks (build_gtdb_database.py --holdout, "
                   "--holdout-clades): absent taxa whose closest species in the sample is one the database lacks, by "
                   "the rank it was held out at (its reads are the likely source of the taxon's); simulated: such "
                   "species in the samples. Other absent taxa: closest to a species the database has.")
        rows = []
        for rank in NOVEL_RANKS:
            mask = absent & (level == rank)
            if mask.any() or simulated.get(rank):
                row = {"held out at": rank, "simulated": simulated.get(rank, 0), **counts(mask, False)}
                row["FP per 100 simulated"] = 100 * row["FP"] / simulated[rank] if simulated.get(rank) else None
                rows.append(row)
        rows.append({"held out at": "other absent taxa", "simulated": None, **counts(absent & (level == ""), False)})
        frame = pd.DataFrame(rows)
        frame["simulated"] = frame["simulated"].map(lambda v: "-" if v is None or pd.isna(v) else int(v))
        report.table(frame)
        data["false_positives_by_rank_held_out"] = rows
        if "meta_relative_rank" in df.columns:
            relative = df["meta_relative_rank"].fillna("").astype(str).to_numpy()
            rows = [{"shared rank": rank, **counts(absent & (level != "") & (relative == rank), False)}
                    for rank in (*NOVEL_RANKS[1:], "domain", "none") if (absent & (level != "") & (relative == rank)).any()]
            report.add("the same false positives by the deepest rank the taxon shares with that species:")
            report.table(pd.DataFrame(rows))
            data["false_positives_by_shared_rank"] = rows

    if ranks_in("meta_neighbour_rank"):
        neighbour = df["meta_neighbour_rank"].fillna("").astype(str).to_numpy()
        rows = [{"closest other species shares": rank, **counts(present & (neighbour == rank), True)}
                for rank in (*NOVEL_RANKS[1:], "domain", "none") if (present & (neighbour == rank)).any()]
        report.add("")
        report.add("False negatives by the deepest rank a present taxon shares with another species in its sample (a "
                   "congener's reads fit it nearly as well):")
        report.table(pd.DataFrame(rows))
        data["false_negatives_by_neighbour_rank"] = rows

    schemes = [k for k in ("species", *CLADE_SCHEMES) if k in p] if clade_rows else []
    if schemes:
        rows = []
        for k in schemes:
            call = p[k] >= opts.knob
            ok = ~np.isnan(p[k])
            fn, fp = int((present & ok & ~call).sum()), int((absent & ok & call).sum())
            m = report.data.get("evaluation", {}).get(k, {})
            rows.append({"held out in training": k, "present": int((present & ok).sum()), "FN": fn,
                         "FN rate": rate(fn, int((present & ok).sum())), "absent": int((absent & ok).sum()), "FP": fp,
                         "FP rate": rate(fp, int((absent & ok).sum())), "FP per sample": m.get("FP_per_sample"),
                         "F1": m.get("F1")})
        report.add("")
        report.add("With the taxon's whole species, genus, family, order, class or phylum held out of training "
                   "(cross-validation): how the model does on parts of the tree its training data lack:")
        report.table(pd.DataFrame(rows))
        data["clades_held_out_in_training"] = rows
    report.data[key] = data


def taxon_classes(df):
    """Each row's class (TAXON_CLASSES) by the collector's meta columns: present taxa beside a congener the database
    lacks in the sample (meta_novel_congener) or not, absent taxa whose closest species in the sample is one the
    database lacks at species rank and of their genus (meta_novel_level, meta_relative_rank: they hold its reads), and
    the other absent taxa; None without those columns."""
    if not {"meta_novel_congener", "meta_novel_level", "meta_relative_rank"} <= set(df.columns):
        return None
    congener = pd.to_numeric(df["meta_novel_congener"], errors="coerce").fillna(0).to_numpy() > 0
    near = ((df["meta_novel_level"].astype(str) == "species") & (df["meta_relative_rank"].astype(str) == "genus")).to_numpy()
    present = df["truth"].to_numpy() == 1
    return np.where(present, np.where(congener, TAXON_CLASSES[1], TAXON_CLASSES[0]),
                    np.where(near, TAXON_CLASSES[2], TAXON_CLASSES[3]))


def study_feature_classes(report, df, title="The conservation features by class of taxon", key="feature_classes"):
    """The medians and quartiles of CLASS_FEATURES by class of taxon: whether the features carry what a relative the
    database lacks does to its congeners' genes, in this training data (on real genomes in a GTDB build)."""
    features = [f for f in CLASS_FEATURES if f in df.columns]
    classes = taxon_classes(df)
    if not features or classes is None:
        return
    report.section(title)
    report.add("median (quartiles). conserved_fast_record_ratio: log2 of the depth of all best records (before the "
               "MAPQ filter) on genes of conservation factor below 1 over that on the others; conserved_fast_kept_ratio: "
               "the same of the records the filters keep; conserved_fast_depth_ratio "
               "and conserved_hit_share: the same of the hit genes' depths after the filters, and the conserved share of "
               "the hit genes. A species' own reads cover both kinds of genes alike; a relative's align best on the "
               "conserved ones (docs/claude/2026-10-01-conservation-pattern).")
    rows, data = [], {}
    for cls in TAXON_CLASSES:
        sel = classes == cls
        if not sel.any():
            continue
        row = {"taxa": cls, "rows": int(sel.sum())}
        data[cls] = {"rows": int(sel.sum())}
        for f in features:
            values = df.loc[sel, f].astype(float).to_numpy()
            values = values[np.isfinite(values)]
            if not len(values):
                continue
            q = np.quantile(values, [0.25, 0.5, 0.75])
            row[f] = f"{q[1]:+.3f} ({q[0]:+.3f}, {q[2]:+.3f})"
            data[cls][f] = [float(v) for v in q]
        rows.append(row)
    report.table(pd.DataFrame(rows))
    report.data[key] = data


def depth_bins(df):
    """Each row's depth bin as protal computes it (profiler::DepthKnobBin): the digits of its sample's fragments over
    all the sample's taxa (rows), less one, 2 to 6."""
    samples = df["meta_sample"].astype(str) if "meta_sample" in df.columns else pd.Series("", index=df.index)
    totals = df["fragments"].astype(float).groupby(samples).transform("sum").to_numpy()
    return np.clip([len(str(int(max(t, 1)))) - 1 for t in totals], DEPTH_KNOB_BINS[0], DEPTH_KNOB_BINS[-1])


def depth_knob_calls(p, bins, knobs, knob):
    """Calls at each row's depth knob, else at knob."""
    return p >= np.array([knobs.get(int(b), knob) for b in bins])


def f1_of(y, call):
    tp, fp, fn = int((call & (y == 1)).sum()), int((call & (y == 0)).sum()), int((~call & (y == 1)).sum())
    return 2 * tp / (2 * tp + fp + fn) if tp + fp + fn else None


def bin_fragments(b):
    """The fragments of depth bin b, as text."""
    if b == DEPTH_KNOB_BINS[0]:
        return f"< 10^{b + 1}"
    if b == DEPTH_KNOB_BINS[-1]:
        return f">= 10^{b}"
    return f"10^{b} - 10^{b + 1}"


def study_depth_knobs(report, df, X, y, p, opts):
    """--depth-knobs: the F1-optimal threshold of each depth bin on species held out; returns {bin: knob}."""
    report.section("Knobs by sample depth (species held out)")
    scores = p.get("species")
    if scores is None:
        splits = folds(df, y, "species", opts)
        if splits is None:
            report.add("the table cannot hold out species: no depth knobs")
            return {}
        scores = predict_out_of_fold(X, y, splits, forest_params(opts))
    bins = depth_bins(df)
    ok = ~np.isnan(scores)
    samples = df["meta_sample"].astype(str).to_numpy() if "meta_sample" in df.columns else np.full(len(df), "")
    knobs, rows = {}, []
    for b in DEPTH_KNOB_BINS:
        sel = ok & (bins == b)
        row = {"bin": b, "fragments": bin_fragments(b), "samples": len(np.unique(samples[sel])), "taxa": int(sel.sum()),
               "present": int(y[sel].sum()), "knob": None, f"F1 at {opts.knob}": f1_of(y[sel], scores[sel] >= opts.knob),
               "F1 at the knob": None}
        if row["taxa"] > DEPTH_KNOB_MIN_TAXA and row["present"] >= DEPTH_KNOB_MIN_PRESENT:
            f1s = [f1_of(y[sel], scores[sel] >= t) or 0.0 for t in DEPTH_KNOB_GRID]
            knobs[b] = float(DEPTH_KNOB_GRID[int(np.argmax(f1s))])
            row["knob"], row["F1 at the knob"] = knobs[b], max(f1s)
        rows.append(row)
    report.table(pd.DataFrame(rows))
    at_knob = f1_of(y[ok], scores[ok] >= opts.knob)
    with_knobs = f1_of(y[ok], depth_knob_calls(scores[ok], bins[ok], knobs, opts.knob))
    report.add(f"F1 of species held out: {fmt(at_knob, 4)} at {opts.knob}, {fmt(with_knobs, 4)} at the depth knobs "
               f"({format_depth_knobs(knobs) or 'none'}; bins without one: {opts.knob}). The knobs are chosen on these "
               "calls, so the second is optimistic; the test set (--test-file) tells.")
    report.data["depth_knobs"] = {"knobs": {str(b): k for b, k in knobs.items()}, "bins": rows,
                                  "F1_at_knob": at_knob, "F1_at_depth_knobs": with_knobs}
    return knobs


def study_test(report, rf, cols, opts, prefix, depth_knobs=None):
    """The fitted forest on an independent test table: samples of another design (depths, community sizes,
    abundance model, strains), which cross-validation on the training data cannot judge. Its probabilities are
    protal's (the PMML is checked to score as the forest does)."""
    report.section(f"Independent test set ({opts.test_file})")
    test = load_table(opts.test_file, opts.taxonomy)
    check_features(test, cols)
    y = test["truth"].to_numpy()
    p = rf.predict_proba(test[cols].to_numpy(dtype=np.float64))[:, 1]
    rows = [({"model": "this one"}, metrics(y, p, test, opts.knob))]
    if "probability" in test.columns:
        rows.append(({"model": "collection model"}, metrics(y, test["probability"].to_numpy(dtype=float), test, opts.knob)))
    report.add(f"{len(test)} taxa in {test['meta_sample'].nunique() if 'meta_sample' in test else 1} samples, "
               f"{int(y.sum())} present")
    report.table(metrics_table(rows))
    report.data["test"] = {label["model"]: m for label, m in rows}
    call = p >= opts.knob
    if "meta_read_pairs" in test.columns:
        depth_rows = []
        for depth, g in test.assign(call=call).groupby("meta_read_pairs"):
            pr, ab = g[g.truth == 1], g[g.truth == 0]
            depth_rows.append({"depth": depth, "samples": g["meta_sample"].nunique() if "meta_sample" in g else 1,
                               "present": len(pr), "FN": int((~pr.call).sum()), "FN rate": rate(int((~pr.call).sum()), len(pr)),
                               "absent": len(ab), "FP": int(ab.call.sum()), "FP rate": rate(int(ab.call.sum()), len(ab))})
        report.add("by depth (read pairs, or bases for long reads):")
        report.table(pd.DataFrame(depth_rows))
        report.data["test_by_depth"] = depth_rows
    if depth_knobs:
        knob_calls = depth_knob_calls(p, depth_bins(test), depth_knobs, opts.knob)
        fp, fn = int((knob_calls & (y == 0)).sum()), int((~knob_calls & (y == 1)).sum())
        report.add(f"at the depth knobs ({format_depth_knobs(depth_knobs)}), as protal calls by default: F1 "
                   f"{fmt(f1_of(y, knob_calls), 4)}, {fp} false positives, {fn} false negatives (at knob {opts.knob}: "
                   f"F1 {fmt(rows[0][1]['F1'], 4)}, {rows[0][1]['FP']} and {rows[0][1]['FN']})")
        report.data["test_depth_knobs"] = {"F1": f1_of(y, knob_calls), "FP": fp, "FN": fn}
    t, prec, rec, f1 = best_threshold(y, p)
    report.add(f"highest F1 on the test set at threshold {t:.3f} (F1 {f1:.4f}; at knob {opts.knob}: "
               f"{fmt(rows[0][1]['F1'], 4)})")
    report.data["test_best_threshold"] = {"threshold": t, "precision": prec, "sensitivity": rec, "F1": f1}
    study_feature_classes(report, test, "Independent test set: the conservation features by class of taxon",
                          "test_feature_classes")
    collection = {"collection model": test["probability"].to_numpy(dtype=float)} if "probability" in test.columns else {}
    study_by_rank(report, test, y, {"species": p, **collection}, opts,
                  title="Independent test set: false positives and false negatives by taxonomic rank",
                  key="test_by_rank", clade_rows=False, scored="scored by the forest fitted on all training rows")
    out = test[[c for c in test.columns if c.startswith("meta_")] + [c for c in ("taxon", "taxon_name", "domain", "truth")
                                                                         if c in test]].copy()
    out["p"] = p
    out.to_csv(prefix + ".test_predictions.tsv.gz", sep="\t", index=False, float_format="%.6g")


def study_threshold(report, y, p, opts, prefix):
    report.section("Threshold (species held out)")
    scores = p.get("species", p["out of bag"])
    t, prec, rec, f1 = best_threshold(y, scores)
    table = threshold_table(y, scores)
    table.to_csv(prefix + ".thresholds.tsv", sep="\t", index=False, float_format="%.6f")
    at = table.iloc[(table.threshold - opts.knob).abs().argmin()]
    report.add(f"at knob {opts.knob}: precision {at.precision:.4f}, sensitivity {at.sensitivity:.4f}, F1 {at.F1:.4f}")
    report.add(f"highest F1 at {t:.3f}: precision {prec:.4f}, sensitivity {rec:.4f}, F1 {f1:.4f}")
    report.add(f"table: {prefix}.thresholds.tsv")
    report.data["threshold"] = {"knob": opts.knob, "at_knob": at.to_dict(), "best_F1": {"threshold": t, "precision": prec,
                                                                                     "sensitivity": rec, "F1": f1}}


def study_features(report, df, y, opts, cols):
    report.section("Feature sets (held out by rows and by species)")
    report.add("A feature set that does much better on rows than on species has learned the training species.")
    sets = {"normalized": feature_columns(df.columns, "normalized"), "all": feature_columns(df.columns, "all")}
    sets = {k: [c for c in v if c not in ("domain",)] for k, v in sets.items()}
    rows = []
    for name, set_cols in sets.items():
        check_features(df, set_cols)
        X = df[set_cols].to_numpy(dtype=np.float64)
        for scheme in ("rows", "species"):
            splits = folds(df, y, scheme, opts)
            if splits is not None:
                rows.append(({"features": f"{name} ({len(set_cols)})", "held out": scheme},
                             metrics(y, predict_out_of_fold(X, y, splits, forest_params(opts)), df, opts.knob)))
    report.table(metrics_table(rows))
    report.data["feature_sets"] = [dict(**label, **m) for label, m in rows]


def study_old_procedure(report, df, y, opts, cols):
    """The procedure of this script before: a grid search over max_features on all rows, then a forest of 512
    trees on the top max_features features by importance, judged on a random 20% of the rows."""
    report.section("The previous procedure (grid search, top features only, 512 trees) against this one")
    X = df[cols].to_numpy(dtype=np.float64)
    t0 = time.time()
    end = min(25, len(cols))
    start = 1 if end < 6 else 6
    grid = GridSearchCV(RandomForestClassifier(**forest_params(opts, n_estimators=128, n_jobs=1,
                                                               class_weight="balanced_subsample")),
                        {"max_features": list(range(start, end + 1))}, cv=5, n_jobs=opts.threads)
    grid.fit(X, y)
    mtry = int(grid.best_params_["max_features"])
    importance = pd.Series(grid.best_estimator_.feature_importances_, index=cols).sort_values(ascending=False)
    top = importance.index[:mtry].tolist()
    grid_s = time.time() - t0
    Xt = df[top].to_numpy(dtype=np.float64)
    old = forest_params(opts, n_estimators=512, max_features="sqrt", max_leaf_nodes=128, min_samples_leaf=1)

    # Its own estimate: a random 20% of the rows.
    rng = np.random.RandomState(opts.seed)
    test = rng.choice(len(df), size=len(df) // 5, replace=False)
    train = np.setdiff1d(np.arange(len(df)), test)
    t0 = time.time()
    rf = RandomForestClassifier(**old).fit(Xt[train], y[train])
    fit_s = time.time() - t0
    p_own = np.full(len(df), np.nan)
    p_own[test] = rf.predict_proba(Xt[test])[:, 1]
    own_t = best_threshold(y[test], p_own[test])[0]
    nodes = sum(e.tree_.node_count for e in rf.estimators_)

    rows = [({"procedure": "previous", "judged on": "random 20% of rows (its report)"}, metrics(y, p_own, df, opts.knob))]
    splits = folds(df, y, "species", opts)
    if splits is not None:
        p_old = predict_out_of_fold(Xt, y, splits, old)
        rows.append(({"procedure": "previous", "judged on": "species held out"}, metrics(y, p_old, df, opts.knob)))
        rows.append(({"procedure": f"previous, its knob {own_t:.3f}", "judged on": "species held out"},
                     metrics(y, p_old, df, own_t)))
        p_new = predict_out_of_fold(X, y, splits, forest_params(opts))
        rows.append(({"procedure": "this one", "judged on": "species held out"}, metrics(y, p_new, df, opts.knob)))
    report.add(f"grid search: {end - start + 1} values of max_features x 5 folds in {grid_s:.1f} s (128 trees each); "
               f"chose {mtry}, so the forest used only these {mtry} of {len(cols)} features: {', '.join(top)}")
    report.add(f"its forest: {nodes} nodes, fitted in {fit_s:.1f} s")
    report.table(metrics_table(rows))
    report.data["previous_procedure"] = {"grid_seconds": grid_s, "mtry": mtry, "features": top, "nodes": nodes,
                                         "knob_from_its_test_rows": own_t,
                                         "results": [dict(**label, **m) for label, m in rows]}


def study_capacity(report, df, X, y, opts):
    report.section("Forest size (species held out)")
    report.add("leaves_per_tree below max_leaves: the limit does not bind. A lower log loss with more leaves: "
               "the data support larger trees (raise --maxnodes).")
    splits = folds(df, y, "species", opts) or folds(df, y, "rows", opts)
    rows = []
    for leaves, min_leaf in ((32, 1), (128, 1), (512, 1), (0, 1), (0, 5)):
        params = forest_params(opts, max_leaf_nodes=leaves or None, min_samples_leaf=min_leaf)
        t0 = time.time()
        grown = []
        p = predict_out_of_fold(X, y, splits, params, leaves=grown)
        m = metrics(y, p, df, opts.knob)
        rows.append(({"max_leaves": leaves or "none", "min_leaf": min_leaf, "leaves_per_tree": float(np.mean(grown)),
                      "fit_s": round((time.time() - t0) / len(splits), 2)}, m))
    report.table(metrics_table(rows)[["max_leaves", "min_leaf", "leaves_per_tree", "fit_s", "AP", "log_loss", "F1",
                                      "sensitivity", "precision", "FP"]])
    # Trees: grow 256 per fold once and score with the first k.
    ks = [k for k in (8, 16, 32, 64, 128, 256)]
    sums = {k: np.full(len(y), np.nan) for k in ks}
    for train, test in splits:
        rf = RandomForestClassifier(**forest_params(opts, n_estimators=max(ks))).fit(X[train], y[train])
        total = np.zeros(len(test))
        for i, tree in enumerate(rf.estimators_, 1):
            total += tree.predict_proba(X[test])[:, 1]
            if i in sums:
                sums[i][test] = total / i
    trows = [({"trees": k}, metrics(y, sums[k], df, opts.knob)) for k in ks]
    report.table(metrics_table(trows)[["trees", "AP", "log_loss", "F1", "sensitivity", "precision", "FP"]])
    report.data["capacity"] = {"leaves": [dict(**l, **m) for l, m in rows], "trees": [dict(**l, **m) for l, m in trows]}


def study_learning_curve(report, df, X, y, opts):
    report.section("More training samples? (species held out, training samples thinned)")
    splits = folds(df, y, "species", opts)
    if splits is None or "meta_sample" not in df.columns:
        report.add("needs meta_sample and several species")
        return
    samples = df["meta_sample"].astype(str).to_numpy()
    rows = []
    for fraction in (0.25, 0.5, 1.0):
        rng = np.random.RandomState(opts.seed)
        keep = set(rng.choice(np.unique(samples), size=max(2, int(round(fraction * len(np.unique(samples))))), replace=False))
        thin = (lambda train: train[np.isin(samples[train], list(keep))]) if fraction < 1 else None
        p = predict_out_of_fold(X, y, splits, forest_params(opts), fit_rows=thin)
        rows.append(({"samples": len(keep), "fraction": fraction}, metrics(y, p, df, opts.knob)))
    report.table(metrics_table(rows)[["samples", "fraction", "AP", "log_loss", "F1", "sensitivity", "precision", "FP"]
                                     + [c for c in metrics_table(rows).columns if c.startswith("sens_")]])
    report.data["learning_curve"] = [dict(**l, **m) for l, m in rows]
    first, last = rows[1][1]["log_loss"], rows[2][1]["log_loss"]
    report.add(f"log loss with half the samples {first:.4f}, with all {last:.4f}: "
               + ("still falling, more samples should help" if last < 0.9 * first else "little change, more of the same samples will not help much"))


# ---- main -------------------------------------------------------------------------------------------------

def main(argv=None):
    opts = parse_args(argv)
    prefix = opts.output_prefix
    os.makedirs(os.path.dirname(os.path.abspath(prefix)), exist_ok=True)
    report = Report()
    started = time.time()
    timing = {}

    report.add(f"# protal presence model: {os.path.basename(prefix)}")
    report.add(f"{datetime.datetime.now().isoformat(timespec='seconds')}  python {platform.python_version()}, "
               f"scikit-learn {sklearn.__version__}, numpy {np.__version__}, pandas {pd.__version__}, "
               f"{os.cpu_count()} CPUs")
    report.add("command: " + " ".join(sys.argv))
    report.data["run"] = {"args": vars(opts), "sklearn": sklearn.__version__, "python": platform.python_version()}

    t0 = time.time()
    df = load_table(opts.truth_file, opts.taxonomy)
    if opts.reference_pmml:
        cols = PmmlForest(opts.reference_pmml).features
        source = f"the inputs of {opts.reference_pmml}"
    else:
        cols = feature_columns(df.columns, opts.features)
        source = f"--features {opts.features}"
    cols = [c for c in cols if c != "domain"]
    check_features(df, cols)
    if df["truth"].nunique() < 2:
        sys.exit("the training table has only present or only absent taxa")
    X = df[cols].to_numpy(dtype=np.float64)
    y = df["truth"].to_numpy()
    timing["load"] = time.time() - t0
    report.add(f"{len(cols)} features ({source}): {', '.join(cols)}")

    warnings = study_data(report, df, y, cols) if opts.evaluation != "none" else []

    t0 = time.time()
    rf = RandomForestClassifier(oob_score=opts.evaluation != "none", **forest_params(opts)).fit(X, y)
    timing["fit"] = time.time() - t0
    nodes = sum(e.tree_.node_count for e in rf.estimators_)
    report.add("")
    report.add(f"forest: {opts.ntree} trees, max {opts.maxnodes or 'unlimited'} leaves, {nodes} nodes, "
               f"fitted in {timing['fit']:.2f} s")

    p = {}
    if opts.evaluation != "none":
        oob = rf.oob_decision_function_[:, 1].copy() if hasattr(rf, "oob_decision_function_") else np.full(len(y), np.nan)
        # scikit-learn (1.9) draws each tree's rows with probability proportional to their class weight, so with
        # many more absent than present taxa some present rows are in every tree's sample: they have no
        # out-of-bag score, and scikit-learn gives them 0. They are left out of the out-of-bag estimate.
        if hasattr(rf, "estimators_samples_"):
            in_bag = np.zeros(len(y), dtype=int)
            for rows in rf.estimators_samples_:
                drawn = np.zeros(len(y), dtype=bool)
                drawn[rows] = True
                in_bag += drawn
            never = in_bag == len(rf.estimators_)
            oob[never] = np.nan
            if never.any():
                report.add(f"out of bag: {int(never.sum())} rows ({int(y[never].sum())} present) were drawn for every tree "
                           "(rows are drawn by their class weight) and have no out-of-bag score; left out")
        t0 = time.time()
        p = study_evaluation(report, df, X, y, opts, oob)
        study_threshold(report, y, p, opts, prefix)
        study_breakdown(report, df, y, p, opts)
        study_by_rank(report, df, y, p, opts)
        study_feature_classes(report, df)
        timing["evaluation"] = time.time() - t0
    if opts.evaluation == "full":
        for name, study in (("feature_sets", lambda: study_features(report, df, y, opts, cols)),
                            ("previous_procedure", lambda: study_old_procedure(report, df, y, opts, cols)),
                            ("capacity", lambda: study_capacity(report, df, X, y, opts)),
                            ("learning_curve", lambda: study_learning_curve(report, df, X, y, opts))):
            t0 = time.time()
            study()
            timing[name] = time.time() - t0
    depth_knobs = {}
    if opts.depth_knobs:
        t0 = time.time()
        depth_knobs = study_depth_knobs(report, df, X, y, p, opts)
        timing["depth_knobs"] = time.time() - t0
    if opts.test_file:
        t0 = time.time()
        rf.n_jobs = 1  # sum the trees in file order, as protal does
        study_test(report, rf, cols, opts, prefix, depth_knobs)
        rf.n_jobs = opts.threads
        timing["test"] = time.time() - t0

    # Export, and check that the file scores as the forest does.
    report.section("Model file")
    t0 = time.time()
    notes = [f"trained {datetime.date.today().isoformat()} on {os.path.abspath(opts.truth_file)}",
             f"{len(df)} taxa, {int(y.sum())} present, {df['meta_sample'].nunique() if 'meta_sample' in df else 1} samples",
             f"{opts.ntree} trees, max leaves {opts.maxnodes or 'unlimited'}, min leaf {opts.min_samples_leaf}, "
             f"max features {opts.max_features}, seed {opts.seed}, scikit-learn {sklearn.__version__}"]
    species = report.data.get("evaluation", {}).get("species")
    if species:
        notes.append(f"species held out: AP {fmt(species['AP'])}, F1 {fmt(species['F1'])} at knob {opts.knob}")
    if depth_knobs:
        notes.append(f"knobs by sample depth (bin of the fragments: knob): {format_depth_knobs(depth_knobs)}")
    write_forest(rf, cols, prefix + ".xml", notes, depth_knobs)
    if read_depth_knobs(prefix + ".xml") != depth_knobs:
        sys.exit(f"{prefix}.xml: the depth knobs read back differ from those written")
    timing["export"] = time.time() - t0
    rf.n_jobs = 1  # sum the trees in file order, as protal does
    sk = rf.predict_proba(X)[:, 1]
    pmml = PmmlForest(prefix + ".xml").predict(X)
    diff = float(np.abs(sk - pmml).max())
    flips = int(((sk >= opts.knob) != (pmml >= opts.knob)).sum())
    size = os.path.getsize(prefix + ".xml")
    report.add(f"{prefix}.xml: {size / 1e6:.2f} MB, {nodes} nodes; written in {timing['export']:.2f} s")
    report.add(f"PMML scored as protal scores it vs scikit-learn, {len(X)} rows: max difference {diff:.3g}, "
               f"{flips} calls differ at knob {opts.knob}")
    report.data["model"] = {"bytes": size, "nodes": nodes, "trees": opts.ntree, "features": cols,
                            "pmml_vs_sklearn_max_diff": diff, "pmml_vs_sklearn_call_differences": flips}
    rf.n_jobs = opts.threads
    joblib.dump(rf, prefix + ".joblib")
    pd.DataFrame({"feature": cols, "importance": rf.feature_importances_}).sort_values(
        "importance", ascending=False).to_csv(prefix + ".varimp.tsv", sep="\t", index=False, float_format="%.6f")
    top = pd.Series(rf.feature_importances_, index=cols).sort_values(ascending=False).head(8)
    report.add("importance: " + ", ".join(f"{k} {v:.3f}" for k, v in top.items()))

    if p:
        keep = [c for c in df.columns if c.startswith("meta_")] + [c for c in ("taxon", "taxon_name", "domain", "truth") if c in df]
        out = df[keep].copy()
        for k, v in p.items():
            out["p_" + k.replace(" ", "_")] = v
        out = pd.concat([out, df[[c for c in DIAGNOSTIC_FEATURES if c in df.columns]]], axis=1)
        out.to_csv(prefix + ".predictions.tsv.gz", sep="\t", index=False, float_format="%.6g")

    report.section("Summary")
    if species:
        report.add(f"species held out: AP {fmt(species['AP'])}, F1 {fmt(species['F1'])}, sensitivity "
                   f"{fmt(species['sensitivity'])}, precision {fmt(species['precision'])}, "
                   f"{fmt(species.get('FP_per_sample'), 2)} false positives per sample at knob {opts.knob}; "
                   + ", ".join(f"{k[len('sensitivity_'):]} sensitivity {fmt(v)}" for k, v in species.items() if k.startswith("sensitivity_")))
        collection = report.data["evaluation"].get("collection model")
        if collection:
            report.add(f"collection model: AP {fmt(collection['AP'])}, F1 {fmt(collection['F1'])}, sensitivity "
                       f"{fmt(collection['sensitivity'])}, precision {fmt(collection['precision'])}; "
                       + ", ".join(f"{k[len('sensitivity_'):]} sensitivity {fmt(v)}" for k, v in collection.items() if k.startswith("sensitivity_")))
        rows = report.data["evaluation"].get("rows")
        if rows and rows.get("F1") is not None and species.get("F1") is not None:
            report.add(f"optimism of random row splits: F1 {fmt(rows['F1'])} on rows vs {fmt(species['F1'])} on species held out")
        by_rank = report.data.get("by_rank", {})
        pct = lambda v: "-" if v is None else f"{v:.2f}%"
        fps = [r for r in by_rank.get("false_positives_by_rank_held_out", []) if r["taxa"]]
        if fps:
            report.add("FP rate of absent taxa closest to species the database lacks, by the rank held out: "
                       + ", ".join(f"{r['held out at']} {pct(r['FP rate'])} ({r['FP']}/{r['taxa']})" for r in fps))
        fns = by_rank.get("false_negatives_by_neighbour_rank", [])
        if fns:
            report.add("FN rate of present taxa by the rank shared with the closest other species in the sample: "
                       + ", ".join(f"{r['closest other species shares']} {pct(r['FN rate'])} ({r['FN']}/{r['taxa']})" for r in fns))
        clades = by_rank.get("clades_held_out_in_training", [])
        if clades:
            report.add("with the clade held out of training: " + "; ".join(
                f"{r['held out in training']} FN rate {pct(r['FN rate'])}, FP rate {pct(r['FP rate'])}, F1 {fmt(r['F1'])}"
                for r in clades))
        previous = {r["judged on"]: r for r in report.data.get("previous_procedure", {}).get("results", [])
                    if r["procedure"] == "previous"}
        if "species held out" in previous:
            own, held = previous.get("random 20% of rows (its report)", {}), previous["species held out"]
            report.add(f"previous procedure: it reported F1 {fmt(own.get('F1'))}; on species held out F1 "
                       f"{fmt(held['F1'])}, AP {fmt(held['AP'])}, log loss {fmt(held['log_loss'], 4)} (this one: F1 "
                       f"{fmt(species['F1'])}, AP {fmt(species['AP'])}, log loss {fmt(species['log_loss'], 4)}); "
                       f"{report.data['previous_procedure']['nodes']} nodes vs {nodes}")
        if species.get("F1") is not None and species["F1"] < 0.9:
            warnings.append(f"F1 on species held out is {species['F1']:.3f}")
    test = report.data.get("test", {}).get("this one")
    if test:
        report.add(f"independent test set: F1 {fmt(test['F1'])}, sensitivity {fmt(test['sensitivity'])}, precision "
                   f"{fmt(test['precision'])}, {fmt(test.get('FP_per_sample'), 2)} false positives per sample at knob "
                   f"{opts.knob}; highest F1 at threshold {report.data['test_best_threshold']['threshold']:.3f}")
        worst = [r for r in report.data.get("test_by_depth", []) if r["FN rate"] is not None]
        if worst:
            w = max(worst, key=lambda r: r["FN rate"])
            report.add(f"independent test set, highest FN rate by depth: {w['FN rate']:.2f}% at {w['depth']}")
        if species and species.get("F1") is not None and test.get("F1") is not None and test["F1"] < species["F1"] - 0.01:
            warnings.append(f"the independent test set scores F1 {test['F1']:.3f} against {species['F1']:.3f} with species "
                            "held out in training: the training design misses what the test set has (see its depth "
                            "table)")
    if diff > 0 or flips:
        warnings.append(f"the PMML file does not score as scikit-learn does (max difference {diff:.3g}): do not use it")
    for w in warnings:
        report.add("WARNING: " + w)
    timing["total"] = time.time() - started
    report.add("time: " + ", ".join(f"{k} {v:.1f} s" for k, v in timing.items()))
    report.data["timing_seconds"] = timing
    report.data["warnings"] = warnings
    report.write(prefix)
    return 1 if diff > 0 or flips else 0


if __name__ == "__main__":
    sys.exit(main())
