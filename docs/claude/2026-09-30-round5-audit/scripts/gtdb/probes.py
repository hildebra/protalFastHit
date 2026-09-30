#!/usr/bin/env python3
"""Unit-level probes of the GTDB build scripts (run with the venv python from ~/audit6/gtdb/src/scripts)."""
import os
import sys
import tempfile
import traceback

S = os.path.expanduser("~/audit6/gtdb/src/scripts")
sys.path.insert(0, S)
sys.path.insert(0, os.path.join(S, "mini_db"))


def section(title):
    print(f"\n=== {title}", flush=True)


section("P1 collect_training_data.design_points: two read setups of one length")
import collect_training_data as collect
opts = collect.parse_args(["--db", "x", "--genome_table", "x", "-o", "x", "--read_setups",
                           "150:HS25:350:50,150:MSv3:550:50", "--read_pairs", "1000,1000"])
names = [p["name"] for p in collect.design_points(opts)]
print("points:", names, "distinct:", len(set(names)), "of", len(names))

section("P2 lineages.from_taxonomy on a taxonomy with a parent cycle")
import lineages
with tempfile.NamedTemporaryFile("w", suffix=".dmp", delete=False) as fh:
    fh.write("id\tparent_id\texternal_id\tname\trank\tlevel\trep_genome\n")
    fh.write("1\t2\t0\ts__A a\tspecies\t7\tGCF_1\n")
    fh.write("2\t3\t0\tg__A\tgenus\t6\t\n")
    fh.write("3\t2\t0\tf__A\tfamily\t5\t\n")
try:
    lineages.from_taxonomy(fh.name)
    print("no error")
except RecursionError as e:
    print("RecursionError:", str(e)[:80])

section("P3 collect: sample name from a dump path")
for path in ["/x/rl100_p1000_s1.profile.truth_annotated", "/x/my.profile_run_s1.profile.truth_annotated"]:
    print(path, "->", os.path.basename(path).split(".profile")[0])

section("P4 model_pmml: parity with scikit-learn on random data, several settings")
import numpy as np
import sklearn
from sklearn.ensemble import RandomForestClassifier
import model_pmml
print("scikit-learn", sklearn.__version__)
rng = np.random.RandomState(0)
n = 20000
X = np.column_stack([rng.lognormal(0, 3, n), rng.rand(n), rng.randint(0, 50, n).astype(float),
                     rng.normal(0, 1e-3, n), np.round(rng.rand(n), 3), rng.lognormal(0, 1, n) * 1e7])
y = ((X[:, 1] + 0.1 * rng.randn(n) > 0.8) | (X[:, 2] > 45)).astype(int)
feats = [f"f{i}" for i in range(X.shape[1])]
path = os.path.join(tempfile.mkdtemp(), "m.xml")
for params in (dict(n_estimators=16, max_leaf_nodes=128, class_weight="balanced"),
               dict(n_estimators=8, max_leaf_nodes=None, min_samples_leaf=1),
               dict(n_estimators=8, max_leaf_nodes=None, class_weight="balanced_subsample", bootstrap=False),
               dict(n_estimators=8, max_leaf_nodes=64, max_features=None, min_weight_fraction_leaf=0.001)):
    rf = RandomForestClassifier(random_state=1, n_jobs=1, **params).fit(X, y)
    model_pmml.write_forest(rf, feats, path)
    Xt = np.vstack([X, rng.lognormal(0, 3, (2000, X.shape[1]))])
    # values exactly at and next to the thresholds
    th = np.concatenate([e.tree_.threshold[e.tree_.feature >= 0] for e in rf.estimators_])
    fi = np.concatenate([e.tree_.feature[e.tree_.feature >= 0] for e in rf.estimators_])
    edge = np.tile(np.median(X, axis=0), (3 * len(th), 1))
    for k, (t, f) in enumerate(zip(th, fi)):
        edge[3 * k, f] = t
        edge[3 * k + 1, f] = np.nextafter(np.float32(t), np.float32(np.inf))
        edge[3 * k + 2, f] = float(np.float32(t)) + 1e-300
    Xt = np.vstack([Xt, edge])
    sk = rf.predict_proba(Xt)[:, 1]
    pm = model_pmml.PmmlForest(path).predict(Xt)
    print(params, "rows", len(Xt), "max |diff|", float(np.abs(sk - pm).max()))

section("P5 model_pmml.leaf_probabilities: tree_.value layout of this scikit-learn")
e = rf.estimators_[0].tree_
print("value shape", e.value.shape, "row sums of first nodes", e.value[:3, 0, :].sum(axis=1))

section("P6 model_pmml: a forest whose leaves are all pure 0 (one class in a bootstrap)")
Xs = rng.rand(50, 2)
ys = np.array([1] + [0] * 49)
rf2 = RandomForestClassifier(n_estimators=20, random_state=3, n_jobs=1).fit(Xs, ys)
pure = sum(1 for est in rf2.estimators_ if est.tree_.node_count == 1)
model_pmml.write_forest(rf2, ["a", "b"], path)
print("single-leaf trees", pure, "max |diff|", float(np.abs(rf2.predict_proba(Xs)[:, 1] - model_pmml.PmmlForest(path).predict(Xs)).max()))
try:
    print("per-tree classes of a single-class bootstrap tree:",
          [est.classes_.tolist() for est in rf2.estimators_ if est.tree_.node_count == 1][:1])
except Exception:
    traceback.print_exc()

section("P7 build_gtdb_database.choose_holdout: clades nested, species per domain")
import build_gtdb_database as build
print("parse_clades('phylum:2,genus:0'):", build.parse_clades("phylum:2,genus:0"))
try:
    build.parse_clades("species:3")
except SystemExit as e:
    print("parse_clades('species:3') ->", e)
