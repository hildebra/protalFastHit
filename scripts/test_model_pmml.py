#!/usr/bin/env python3
"""Checks for model_pmml.py: a forest or a gradient-boosted model written as PMML scores as scikit-learn
scores it, bit for bit.

The PMML file is scored with model_pmml.PmmlForest, which compares doubles and averages the trees in
file order as protal's cPMML does, or PmmlBoosted, which chains the trees into a logit sum as cPMML
does; check_model_parity.py compares them with protal itself. Also the trainer (machine_learning_cmdline.py): the knobs
by sample depth that --depth-knobs chooses and writes, the calls at a target share of false calls, the scenarios and
--features auto; and the rules the trainer shares with protal on the golden vectors of tests/data/golden_model_rules.tsv,
which tests/test_GoldenModelRules.cpp checks protal on.

Needs numpy; all but Float32SplitTest need scikit-learn, pandas and joblib (skipped without them, failed with
PROTAL_TESTS_REQUIRED=1). The fits run on one thread (OMP_NUM_THREADS, threadpoolctl), the trainers with --threads 1.

  python3 -m unittest scripts/test_model_pmml.py
"""

import json
import os
import re
import subprocess
import sys
import tempfile
import types
import unittest

# One OpenMP and BLAS thread for the fits here and in the trainers these tests run (they inherit the environment):
# gradient boosting takes every core by default, which on a machine of many cores, or a busy one, made this suite take
# 20 minutes instead of a few (BoostedExportTest 189 s against 1 s).
for _variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_variable, "1")

import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import prerequisites  # noqa: E402
from model_pmml import (PmmlBoosted, PmmlForest, float32_split, format_depth_knob_curve, format_depth_knobs,  # noqa: E402
                        load_model, read_depth_knob_curve, read_depth_knobs, read_false_calls, write_forest, write_model)

try:
    import joblib  # noqa: F401  (the trainer's)
    import pandas as pd
    from sklearn.ensemble import RandomForestClassifier
    from threadpoolctl import threadpool_limits
    threadpool_limits(1)  # also when numpy or scikit-learn were loaded before the variables above were set
    HAVE_SKLEARN = True
except ImportError:  # the checks need scikit-learn, as training does
    RandomForestClassifier = None
    HAVE_SKLEARN = False
NEEDS_SKLEARN = "needs scikit-learn, pandas and joblib"

GOLDEN = os.path.join(HERE, "..", "tests", "data", "golden_model_rules.tsv")


def golden_cases(rule):
    """The cases of a rule in tests/data/golden_model_rules.tsv: [(name, {field: text})]."""
    cases = []
    with open(GOLDEN) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0] == rule:
                cases.append((fields[1], dict(f.split("=", 1) for f in fields[2:])))
    return cases


def golden_numbers(text):
    """"0.1,0.2*3" -> [0.1, 0.2, 0.2, 0.2]"""
    values = []
    for item in filter(None, text.split(",")):
        value, _, copies = item.partition("*")
        values += [float(value)] * (int(copies) if copies else 1)
    return values


def golden_curve(text):
    """"2:0.2,4:0.8" -> [(2.0, 0.2), (4.0, 0.8)]"""
    return [tuple(float(v) for v in point.split(":")) for point in filter(None, text.split(","))]


class Float32SplitTest(unittest.TestCase):
    def test_matches_float32_comparison(self):
        rng = np.random.default_rng(1)
        thresholds = np.concatenate([rng.normal(0, 1, 2000), rng.uniform(0, 1e4, 2000), [0.0, -0.0, 0.5, 1.0, -1.0]])
        for t in thresholds:
            split = float32_split(t)
            candidates = [split, np.nextafter(split, np.inf), np.nextafter(split, -np.inf), t,
                          float(np.float32(t)), np.nextafter(float(np.float32(t)), np.inf)]
            for x in candidates:
                self.assertEqual(x <= split, float(np.float32(x)) <= t, f"x={x!r} t={t!r} split={split!r}")

    def test_ties_round_to_even(self):
        # 1 + 2^-24 is halfway between the float32s 1 and 1 + 2^-23 and rounds to 1 (even significand).
        self.assertEqual(float32_split(1.0), 1.0 + 2.0 ** -24)
        # halfway above 1 + 2^-23 (odd) rounds up, so the split lies just below the midpoint.
        odd = 1.0 + 2.0 ** -23
        self.assertEqual(float32_split(odd), np.nextafter(odd + 2.0 ** -24, -np.inf))


@prerequisites.requires(HAVE_SKLEARN, NEEDS_SKLEARN)
class ForestExportTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rng = np.random.default_rng(7)
        n = 3000
        X = np.column_stack([rng.random(n), rng.normal(0, 3, n), rng.integers(0, 50, n).astype(float),
                             rng.random(n) ** 8, rng.lognormal(0, 2, n)])
        y = ((X[:, 0] + 0.2 * X[:, 1] + rng.normal(0, 0.3, n)) > 0.6).astype(int)
        cls.features = ["a", "b", "count", "small", "wide"]
        cls.X, cls.y = X, y
        cls.tmp = tempfile.TemporaryDirectory()

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def export(self, **params):
        rf = RandomForestClassifier(random_state=3, n_jobs=1, class_weight="balanced", **params).fit(self.X, self.y)
        path = os.path.join(self.tmp.name, "model.xml")
        write_forest(rf, self.features, path, ["test & <annotation>"])
        return rf, PmmlForest(path)

    def test_same_probabilities(self):
        for params in ({"n_estimators": 16, "max_leaf_nodes": 32}, {"n_estimators": 8}, {"n_estimators": 4, "max_depth": 1}):
            rf, pmml = self.export(**params)
            self.assertEqual(pmml.features, self.features)
            self.assertEqual(len(pmml.trees), params["n_estimators"])
            np.testing.assert_array_equal(pmml.predict(self.X), rf.predict_proba(self.X)[:, 1])

    def test_values_on_split_boundaries(self):
        # Rows whose value is at, or one double either side of, a threshold and its float32 rounding points:
        # where comparing doubles with sklearn's thresholds unchanged would send rows the other way.
        rf, pmml = self.export(n_estimators=8, max_leaf_nodes=64)
        base = np.median(self.X, axis=0)
        rows, naive_disagree = [], 0
        for tree in rf.estimators_:
            t = tree.tree_
            for f, thr in zip(t.feature[t.feature >= 0], t.threshold[t.feature >= 0]):
                split = float32_split(thr)
                for x in (split, np.nextafter(split, np.inf), np.nextafter(split, -np.inf), thr, float(np.float32(thr))):
                    row = base.copy()
                    row[f] = x
                    rows.append(row)
                    naive_disagree += (x <= thr) != (float(np.float32(x)) <= thr)
        rows = np.array(rows)
        self.assertGreater(naive_disagree, 0)
        np.testing.assert_array_equal(pmml.predict(rows), rf.predict_proba(rows)[:, 1])

    def test_dataframe_input_in_any_column_order(self):
        rf, pmml = self.export(n_estimators=8)
        frame = pd.DataFrame(self.X, columns=self.features)[self.features[::-1]]
        np.testing.assert_array_equal(pmml.predict(frame), rf.predict_proba(self.X)[:, 1])

    def test_model_contract(self):
        self.export(n_estimators=2)
        with open(os.path.join(self.tmp.name, "model.xml")) as fh:
            xml = fh.read()
        self.assertIn('<MiningField name="truth" usageType="predicted"/>', xml)
        self.assertIn('<Value value="TRUE"/>', xml)
        self.assertIn("test &amp; &lt;annotation&gt;", xml)

    def test_depth_knob_curve(self):
        # Written as an Extension in the header, before the Application (PMML's order), and read back.
        rf = RandomForestClassifier(n_estimators=2, random_state=1).fit(self.X, self.y)
        path = os.path.join(self.tmp.name, "knobs.xml")
        write_forest(rf, self.features, path, depth_knob_curve=[(1.3, 0.12), (4.321, 0.9)])
        with open(path) as fh:
            xml = fh.read()
        self.assertIn('<Extension name="protal_depth_knob_curve" value="1.300:0.12,4.321:0.9"/>\n  <Application ', xml)
        self.assertEqual(read_depth_knob_curve(path), [(1.3, 0.12), (4.321, 0.9)])
        self.assertEqual(read_depth_knobs(path), {}, "no knobs by bin in a model with a curve")
        np.testing.assert_array_equal(PmmlForest(path).predict(self.X), rf.predict_proba(self.X)[:, 1])
        write_forest(rf, self.features, path)
        self.assertEqual(read_depth_knob_curve(path), [])
        self.assertEqual(format_depth_knob_curve([(2, 0.5), (3.25, 0.05)]), "2.000:0.5,3.250:0.05")
        # An older model's knobs by decade are still read.
        with open(path, "w") as fh:
            fh.write(xml.replace('name="protal_depth_knob_curve" value="1.300:0.12,4.321:0.9"',
                                 'name="protal_depth_knobs" value="2:0.6,4:0.37"'))
        self.assertEqual(read_depth_knobs(path), {2: 0.6, 4: 0.37})
        self.assertEqual(read_depth_knob_curve(path), [])
        self.assertEqual(format_depth_knobs({6: 0.05, 3: 0.5}), "3:0.5,6:0.05")

    def test_false_calls(self):
        # The calibration, prior and target as Extensions (profiler::ParseFalseCalls), each number as its double's
        # shortest text, read back exactly; none without them.
        rf = RandomForestClassifier(n_estimators=2, random_state=1).fit(self.X, self.y)
        path = os.path.join(self.tmp.name, "calls.xml")
        calls = {"curve": [(0.0, 0.0001), (0.1 + 0.2, 0.25), (1.0, 0.99)], "prior": 0.2137, "fdr": 0.05}
        write_forest(rf, self.features, path, false_calls=calls, depth_knob_curve=[(2.0, 0.5)])
        with open(path) as fh:
            xml = fh.read()
        self.assertIn('<Extension name="protal_calibration" value="0.0:0.0001,0.30000000000000004:0.25,1.0:0.99"/>', xml)
        self.assertIn('<Extension name="protal_prior" value="0.2137"/>', xml)
        self.assertIn('<Extension name="protal_fdr" value="0.05"/>', xml)
        self.assertEqual(read_false_calls(path), calls)
        self.assertEqual(read_depth_knob_curve(path), [(2.0, 0.5)])
        write_forest(rf, self.features, path)
        self.assertIsNone(read_false_calls(path))

    def test_rejects_other_classes(self):
        rf = RandomForestClassifier(n_estimators=2, random_state=1).fit(self.X, np.where(self.y == 1, "yes", "no"))
        with self.assertRaises(ValueError):
            write_forest(rf, self.features, os.path.join(self.tmp.name, "bad.xml"))


@prerequisites.requires(HAVE_SKLEARN, NEEDS_SKLEARN)
class BoostedExportTest(unittest.TestCase):
    """A HistGradientBoostingClassifier written as a chain of trees into a logit RegressionModel (write_boosted) scores
    as scikit-learn scores it, bit for bit (PmmlBoosted scores it as cPMML does)."""

    @classmethod
    def setUpClass(cls):
        from sklearn.ensemble import HistGradientBoostingClassifier
        cls.Boosted = HistGradientBoostingClassifier
        rng = np.random.default_rng(11)
        n = 4000
        X = np.column_stack([rng.random(n), rng.normal(0, 3, n), rng.integers(0, 50, n).astype(float),
                             rng.random(n) ** 8, rng.lognormal(0, 2, n)])
        y = ((X[:, 0] + 0.2 * X[:, 1] - 0.01 * X[:, 2] + rng.normal(0, 0.4, n)) > 0.5).astype(int)
        cls.features = ["a", "b", "count", "small", "wide"]
        cls.X, cls.y = X, y
        cls.weights = np.where(rng.random(n) < 0.3, 0.25, 1.0)
        cls.tmp = tempfile.TemporaryDirectory()

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def export(self, X=None, **params):
        X = self.X if X is None else X
        settings = dict(max_iter=40, learning_rate=0.1, max_leaf_nodes=15, min_samples_leaf=20, l2_regularization=1.0,
                        class_weight="balanced", early_stopping=False, random_state=3)
        model = self.Boosted(**{**settings, **params}).fit(X, self.y, sample_weight=self.weights)
        path = os.path.join(self.tmp.name, "boosted.xml")
        write_model(model, self.features, path, ["boosted & <annotation>"])
        return model, load_model(path)

    def test_same_probabilities(self):
        for params in ({}, {"max_leaf_nodes": 63, "max_iter": 25}, {"max_depth": 1, "max_iter": 10}):
            model, pmml = self.export(**params)
            self.assertIsInstance(pmml, PmmlBoosted)
            self.assertEqual(pmml.features, self.features)
            self.assertEqual(len(pmml.trees), model.n_iter_)
            np.testing.assert_array_equal(pmml.raw(self.X), model.decision_function(self.X))
            np.testing.assert_array_equal(pmml.predict(self.X), model.predict_proba(self.X)[:, 1])

    def test_values_on_split_boundaries(self):
        # Rows at a threshold and one double either side: the thresholds are written exactly.
        model, pmml = self.export()
        base = np.median(self.X, axis=0)
        rows = []
        for trees in model._predictors:
            nodes = trees[0].nodes
            for n in nodes[nodes["is_leaf"] == 0]:
                t = float(n["num_threshold"])
                for x in (t, np.nextafter(t, np.inf), np.nextafter(t, -np.inf)):
                    row = base.copy()
                    row[int(n["feature_idx"])] = x
                    rows.append(row)
        rows = np.array(rows)
        np.testing.assert_array_equal(pmml.predict(rows), model.predict_proba(rows)[:, 1])

    def test_missing_values_go_where_sklearn_sends_them(self):
        X = self.X.copy()
        X[::7, 1] = np.nan  # missing in training: each split learns a side for them
        model, pmml = self.export(X=X)
        test = self.X[:500].copy()
        test[::3, 1] = np.nan
        test[1::5, 3] = np.nan  # never missing in training
        np.testing.assert_array_equal(pmml.predict(test), model.predict_proba(test)[:, 1])

    def test_a_tree_of_one_leaf(self):
        # A round whose best split gains nothing is one leaf: written as the root's only child (cPMML scores a root's
        # children).
        model, pmml = self.export(max_iter=5, min_samples_leaf=len(self.y))
        self.assertTrue(all(trees[0].nodes[0]["is_leaf"] for trees in model._predictors))
        np.testing.assert_array_equal(pmml.predict(self.X), model.predict_proba(self.X)[:, 1])

    def test_model_contract(self):
        # protal reads every MiningField of the file: the target once, as predicted, the rest features; the chained
        # trees' outputs are no MiningFields.
        self.export(max_iter=3)
        with open(os.path.join(self.tmp.name, "boosted.xml")) as fh:
            xml = fh.read()
        import re
        fields = re.findall(r'<MiningField name="([^"]+)"( usageType="predicted")?/>', xml)
        self.assertEqual({name for name, predicted in fields if predicted}, {"truth"})
        self.assertTrue({name for name, predicted in fields if not predicted} <= set(self.features))
        self.assertIn('<Value value="TRUE"/>', xml)
        self.assertIn('multipleModelMethod="modelChain"', xml)
        self.assertIn('normalizationMethod="logit"', xml)
        self.assertIn('<OutputField name="tree_3" optype="continuous" dataType="double" feature="predictedValue"/>', xml)
        self.assertNotIn("Interval", xml)
        self.assertIn("boosted &amp; &lt;annotation&gt;", xml)

    def test_header_extensions(self):
        model = self.Boosted(max_iter=3, random_state=1).fit(self.X, self.y)
        path = os.path.join(self.tmp.name, "boosted_knobs.xml")
        calls = {"curve": [(0.0, 0.0001), (1.0, 0.99)], "prior": 0.2, "fdr": 0.05}
        write_model(model, self.features, path, depth_knob_curve=[(5.5, 0.7)], false_calls=calls)
        self.assertEqual(read_depth_knob_curve(path), [(5.5, 0.7)])
        self.assertEqual(read_false_calls(path), calls)

    def test_rejects_other_classes(self):
        model = self.Boosted(max_iter=2).fit(self.X, np.where(self.y == 1, "yes", "no"))
        with self.assertRaises(ValueError):
            write_model(model, self.features, os.path.join(self.tmp.name, "bad.xml"))

    def test_loads_either_model(self):
        rf = RandomForestClassifier(n_estimators=2, random_state=1).fit(self.X, self.y)
        path = os.path.join(self.tmp.name, "either.xml")
        write_model(rf, self.features, path)
        self.assertIsInstance(load_model(path), PmmlForest)


@prerequisites.requires(HAVE_SKLEARN, NEEDS_SKLEARN)
class FeatureSetsTest(unittest.TestCase):
    """The feature sets' names: a set joins its groups with "+" and takes their features in the groups' order, the
    default set too; normalized+adjacency leaves the rest out, normalized the gene neighbours' too, and
    normalized+adjacency+relatives adds all the relatives'; a table without a set's features says which it lacks."""

    def test_default_set_has_the_gene_neighbour_features(self):
        import model_features as mf
        import machine_learning_cmdline
        columns = (["truth", "taxon", "meta_sample"] + mf.NORMALIZED_FEATURES + mf.ADJACENCY_FEATURES +
                   mf.RELATIVE_FEATURES + mf.SAMPLE_FEATURES + mf.DIVERGENCE_FEATURES + mf.UNFILTERED_FEATURES +
                   mf.REF_FEATURES + mf.COMPLEXITY_FEATURES + mf.CONSISTENCY_FEATURES + mf.SHAPE_FEATURES +
                   mf.NEIGHBOURHOOD_FEATURES + mf.ANCESTRY_FEATURES + mf.GAP_FEATURES + mf.FOREIGN_FEATURES + mf.UNTRIED_FEATURES +
                   mf.PRIORS_FEATURES + ["genus_top_fragments", "other"])
        # The priors are opt-in: their gain at r226 is the cluster-size rule the simulation cannot test. The reference's
        # k-mer uniqueness (ref) and the sample's complexity are in the default set since 2026-10-06, the groups against
        # the false positives of complex communities (consistency, shape, neighbourhood) since 2026-10-07, the ancestry
        # sites since 0.7.9 (2026-10-07), and the per-copy tables and the untried candidates (gaps, foreign, untried) since
        # the merge of congener-gaps (2026-10-07).
        self.assertEqual(mf.DEFAULT_FEATURE_SET, "normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+"
                                                 "consistency+shape+neighbourhood+ancestry+gaps+foreign+untried")
        copies = mf.GAP_FEATURES + mf.FOREIGN_FEATURES + mf.UNTRIED_FEATURES
        self.assertEqual(len(copies), 8)
        self.assertEqual(len(set(copies)), 8)
        self.assertEqual(mf.REF_FEATURES, ["su_rate_ref", "lu_rate_ref", "lsu_rate_ref"])
        self.assertEqual(mf.COMPLEXITY_FEATURES, ["sample_log_taxa", "sample_low_identity", "sample_identity"])
        self.assertEqual(mf.ANCESTRY_FEATURES, ["ancestry_sites_per_record", "ancestry_agreement", "ancestry_congener_share"])
        new = mf.CONSISTENCY_FEATURES + mf.SHAPE_FEATURES + mf.NEIGHBOURHOOD_FEATURES + mf.ANCESTRY_FEATURES
        self.assertEqual(len(new), 18)
        self.assertEqual(len(set(new)), 18)
        self.assertEqual(mf.feature_columns(columns, mf.DEFAULT_FEATURE_SET),
                         mf.NORMALIZED_FEATURES + mf.ADJACENCY_FEATURES + mf.DISTANCE_FEATURES + mf.SAMPLE_FEATURES +
                         mf.DIVERGENCE_FEATURES + mf.UNFILTERED_FEATURES + mf.REF_FEATURES + mf.COMPLEXITY_FEATURES + new + copies)
        self.assertEqual(mf.feature_columns(columns, mf.DEFAULT_FEATURE_SET + "+priors"),
                         mf.NORMALIZED_FEATURES + mf.ADJACENCY_FEATURES + mf.DISTANCE_FEATURES + mf.SAMPLE_FEATURES +
                         mf.DIVERGENCE_FEATURES + mf.UNFILTERED_FEATURES + mf.REF_FEATURES + mf.COMPLEXITY_FEATURES +
                         new + copies + mf.PRIORS_FEATURES)
        # A table of a protal before the per-copy tables: a clear error with the default, the set without them works.
        before = [c for c in columns if c not in copies]
        with self.assertRaisesRegex(RuntimeError, "gap_informative_share"):
            mf.feature_columns(before, mf.DEFAULT_FEATURE_SET)
        self.assertEqual(mf.feature_columns(before, "normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+"
                                                    "consistency+shape+neighbourhood+ancestry"),
                         mf.NORMALIZED_FEATURES + mf.ADJACENCY_FEATURES + mf.DISTANCE_FEATURES + mf.SAMPLE_FEATURES +
                         mf.DIVERGENCE_FEATURES + mf.UNFILTERED_FEATURES + mf.REF_FEATURES + mf.COMPLEXITY_FEATURES + new)
        self.assertIn("normalized+adjacency+distance+depth+divergence+unfiltered", mf.AUTO_CANDIDATES)  # an old default
        # A table of the r226 v15 build (before the false-positive groups): a clear error with the default, its old set works.
        v15 = [c for c in columns if c not in new and c not in copies]
        with self.assertRaisesRegex(RuntimeError, "read_consensus_share"):
            mf.feature_columns(v15, mf.DEFAULT_FEATURE_SET)
        self.assertIn("normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity", mf.AUTO_CANDIDATES)
        self.assertEqual(mf.feature_columns(v15, "normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity"),
                         mf.NORMALIZED_FEATURES + mf.ADJACENCY_FEATURES + mf.DISTANCE_FEATURES + mf.SAMPLE_FEATURES +
                         mf.DIVERGENCE_FEATURES + mf.UNFILTERED_FEATURES + mf.REF_FEATURES + mf.COMPLEXITY_FEATURES)
        # A table of the r226 v14 build (before the sample's complexity): a clear error with the default, its old set works.
        v14 = [c for c in v15 if c not in mf.COMPLEXITY_FEATURES]
        with self.assertRaisesRegex(RuntimeError, "sample_log_taxa"):
            mf.feature_columns(v14, "normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity")
        self.assertEqual(mf.feature_columns(v14, "normalized+adjacency+distance+depth+divergence+unfiltered+ref"),
                         mf.NORMALIZED_FEATURES + mf.ADJACENCY_FEATURES + mf.DISTANCE_FEATURES + mf.SAMPLE_FEATURES +
                         mf.DIVERGENCE_FEATURES + mf.UNFILTERED_FEATURES + mf.REF_FEATURES)
        self.assertIn(mf.DEFAULT_FEATURE_SET, mf.FEATURE_SETS)
        self.assertIn(mf.DEFAULT_FEATURE_SET, mf.AUTO_CANDIDATES)
        self.assertEqual(mf.feature_columns(columns, "normalized+adjacency+distance+depth+divergence"),
                         mf.NORMALIZED_FEATURES + mf.ADJACENCY_FEATURES + mf.DISTANCE_FEATURES + mf.SAMPLE_FEATURES + mf.DIVERGENCE_FEATURES)
        self.assertEqual(mf.feature_columns(columns, "normalized+adjacency+distance"),
                         mf.NORMALIZED_FEATURES + mf.ADJACENCY_FEATURES + mf.DISTANCE_FEATURES)
        # Groups in any order, each once, normalized among them; the sample's depth is recognised.
        self.assertEqual(mf.feature_columns(columns, "depth+normalized"), mf.NORMALIZED_FEATURES + mf.SAMPLE_FEATURES)
        self.assertTrue(mf.has_sample_depth(mf.feature_columns(columns, "normalized+depth")))
        self.assertFalse(mf.has_sample_depth(mf.feature_columns(columns, "normalized+adjacency+distance")))
        for bad in ("normalized+adjacency+depth+depth", "adjacency", "normalized+ani", ""):
            with self.assertRaises(ValueError):
                mf.feature_set_columns(bad)
        self.assertEqual(mf.feature_set_name("normalized+divergence"), "normalized+divergence")
        self.assertEqual(mf.feature_set_name("all"), "all")
        # A table of an older protal lacks the new columns: a clear error, and the older set still works.
        older = ["truth", "taxon"] + mf.NORMALIZED_FEATURES + mf.ADJACENCY_FEATURES + mf.RELATIVE_FEATURES
        with self.assertRaises(RuntimeError):
            mf.feature_columns(older, "normalized+adjacency+depth")
        self.assertEqual(mf.feature_columns(older, "normalized+adjacency"), mf.NORMALIZED_FEATURES + mf.ADJACENCY_FEATURES)
        relatives = mf.feature_columns(columns, "normalized+adjacency+relatives")
        self.assertEqual(relatives, mf.NORMALIZED_FEATURES + mf.ADJACENCY_FEATURES + mf.RELATIVE_FEATURES)
        self.assertEqual(mf.feature_columns(columns, "normalized"), mf.NORMALIZED_FEATURES)
        self.assertTrue(set(mf.DISTANCE_FEATURES) <= set(mf.RELATIVE_FEATURES))
        self.assertNotIn("genus_top_fragments", relatives)  # it counts reads
        with self.assertRaisesRegex(RuntimeError, "adjacent_support"):  # a table of an older protal
            mf.feature_columns([c for c in columns if c != "adjacent_support"], "normalized+adjacency")
        with self.assertRaisesRegex(RuntimeError, "relative_spill"):  # a table of protal before the relatives features
            mf.feature_columns([c for c in columns if c != "relative_spill"], "normalized+distance")
        with self.assertRaisesRegex(RuntimeError, "em_own_share"):  # a table of protal before the relatives features
            mf.feature_columns([c for c in columns if c != "em_own_share"], "normalized+adjacency+relatives")
        # The trainer chooses its set by default, the default set among the candidates.
        opts = machine_learning_cmdline.parse_args(["--truth-file", "t.tsv", "--output-prefix", "p"])
        self.assertEqual(opts.features, "auto")
        self.assertIn(mf.DEFAULT_FEATURE_SET, mf.auto_candidates("auto"))


@prerequisites.requires(HAVE_SKLEARN, NEEDS_SKLEARN)
class ParityFeaturesTest(unittest.TestCase):
    """check_model_parity.py: a feature that differs in its last digits only (a sum added up in another order) is
    rounding, which the check notes; one that differs more, protal computes differently."""

    def test_last_digits_are_rounding(self):
        import check_model_parity as cmp
        joined = pd.DataFrame({"same": [0.0, 0.25], "same_collected": [0.0, 0.25],
                               "summed": [0.1 + 0.2 + 0.3, 7.0], "summed_collected": [0.1 + (0.2 + 0.3), 7.0],
                               "changed": [0.5, 1.0], "changed_collected": [0.5, 1.01]})
        found = cmp.feature_differences(joined, ["same", "summed", "changed"])
        self.assertEqual(set(found), {"summed", "changed"})
        self.assertGreater(found["summed"], 0)
        differs, rounded = cmp.split_rounding(found)
        self.assertEqual(set(differs), {"changed"})
        self.assertAlmostEqual(differs["changed"], 0.01 / 1.01)
        self.assertEqual(set(rounded), {"summed"})

    def test_sample_without_taxa(self):
        """A shallow sample of a database of few genes can have no taxon: its dump is the header alone, which is
        no problem when the collection's had none either."""
        import check_model_parity as cmp

        class Model:
            features = ["f"]

            def predict(self, X):
                return np.full(len(X), 0.25)

        columns = ["truth", "prediction", "probability", "taxon", "f"]
        none = pd.DataFrame({c: pd.Series(dtype=float) for c in columns})
        one = pd.DataFrame([[1, 1, 0.25, 7, 2.0]], columns=columns)
        self.assertEqual(cmp.compare_sample(Model(), "s", none, none), ([], {}))
        problems, _ = cmp.compare_sample(Model(), "s", none, one)
        self.assertEqual(problems, ["s: 0 taxa now, 1 during collection"])
        self.assertEqual(cmp.compare_sample(Model(), "s", one, one), ([], {}))
        problems, _ = cmp.compare_sample(Model(), "s", one.assign(probability=0.5), one.assign(probability=0.5))
        self.assertEqual(len(problems), 1)
        self.assertIn("differs from the model file's by up to 0.25", problems[0])


@prerequisites.requires(HAVE_SKLEARN, NEEDS_SKLEARN)
class FoldJobsTest(unittest.TestCase):
    """The trainer's cross-validation folds fitted side by side in worker processes (FOLD_JOBS, --fold-jobs) score as
    when fitted one after another: neither model depends on its threads."""

    def test_side_by_side_as_one_after_another(self):
        import machine_learning_cmdline as trainer
        from sklearn.model_selection import GroupKFold
        rng = np.random.default_rng(3)
        n = 1500
        X = rng.normal(size=(n, 6))
        y = ((X[:, 0] + X[:, 1] * X[:, 2] + rng.normal(0, 0.5, n)) > 0.2).astype(int)
        splits = list(GroupKFold(5).split(X, y, rng.integers(0, 40, n)))
        saved = trainer.FOLD_JOBS, trainer.FOLD_JOB_THREADS, trainer.ROW_WEIGHTS
        try:
            trainer.ROW_WEIGHTS = np.where(rng.random(n) < 0.3, 0.25, 1.0)  # the weights travel with the rows
            for model in ("gbm", "forest"):
                opts = trainer.parse_args(["--truth-file", "t.tsv", "--output-prefix", "p", "--model", model,
                                           "--rounds", "20", "--ntree", "8", "--threads", "2"])
                params = trainer.model_params(opts)
                scores, leaves = {}, {}
                for jobs, threads in ((1, 2), (3, 1)):
                    trainer.FOLD_JOBS, trainer.FOLD_JOB_THREADS = jobs, threads
                    leaves[jobs] = []
                    scores[jobs] = trainer.predict_out_of_fold(X, y, splits, params, leaves=leaves[jobs])
                np.testing.assert_array_equal(scores[1], scores[3], model)
                self.assertEqual(leaves[1], leaves[3], model)
                self.assertFalse(np.isnan(scores[1]).any())
        finally:
            trainer.FOLD_JOBS, trainer.FOLD_JOB_THREADS, trainer.ROW_WEIGHTS = saved


def metrics_knobs(test, name):
    """The knob points of the model `name` trained in a TrainerDepthKnobsTest."""
    with open(os.path.join(test.tmp.name, name + ".metrics.json")) as fh:
        return json.load(fh)["depth_knobs"]["points"]


@prerequisites.requires(HAVE_SKLEARN, NEEDS_SKLEARN)
class TrainerDepthKnobsTest(unittest.TestCase):
    """machine_learning_cmdline.py --depth-knobs on a table of shallow samples (hundreds of fragments, log10 ~2.7) and
    deep ones (tens of thousands, ~4.7), where a present taxon's evidence grows with depth; half the present taxa
    simulated from another genome than the representative. Each training runs once for the class (train)."""

    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, HERE)
        import machine_learning_cmdline
        cls.trainer = machine_learning_cmdline
        cls.tmp = tempfile.TemporaryDirectory()
        cls.trained = {}
        rng = np.random.default_rng(5)
        rows = []
        for sample in range(40):
            deep = sample % 2 == 1
            for taxon in range(30):
                present = taxon < 8
                fragments = rng.integers(500, 3000) if deep else rng.integers(5, 30)
                signal = (1.5 if deep else 0.6) * present
                near = not present and taxon < 12  # absent taxa that hold a held-out congener's reads
                rows.append({"truth": int(present), "taxon": 100 + (taxon + sample) % 60, "taxon_name": "t",
                             "meta_sample": f"s{sample}", "fragments": float(fragments),
                             "x": signal + rng.normal(0, 0.5), "y": rng.normal(0, 1),
                             "conserved_fast_record_ratio": (1.5 if near else 0.0) + rng.normal(0, 0.2),
                             "meta_novel_congener": int(present and taxon < 2),
                             "meta_novel_level": "species" if near else "",
                             "meta_relative_rank": "genus" if near else "",
                             "meta_rep_genome": (taxon % 2) if present else "",
                             "identity": (0.99 if taxon % 2 else 0.96) if present else 0.93})
        cls.table = os.path.join(cls.tmp.name, "training.tsv")
        pd.DataFrame(rows).to_csv(cls.table, sep="\t", index=False)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def train(self, name, *extra):
        """The trainer on the table with the options `extra`, once per name -> (model file, metrics)."""
        if name not in self.trained:
            prefix = os.path.join(self.tmp.name, name)
            result = subprocess.run([sys.executable, os.path.join(HERE, "machine_learning_cmdline.py"), "--truth-file",
                                     self.table, "--output-prefix", prefix, "--features", "all", "--ntree", "16", "--rounds",
                                     "30", "--evaluation", "basic", "--threads", "1", *extra], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout[-3000:] + result.stderr[-3000:])
            with open(prefix + ".metrics.json") as fh:
                self.trained[name] = (prefix + ".xml", json.load(fh))
        return self.trained[name]

    def test_depths_and_curve_as_protal(self):
        # profiler::DepthKnobAt: log10 of the sample's fragments over all its taxa (at least 1); the knobs along the
        # curve are checked on the golden vectors (GoldenModelRulesTest).
        frame = pd.DataFrame({"meta_sample": ["a", "a", "b", "c", "e"], "fragments": [60.0, 40.0, 1000.0, 3e7, 0.0]})
        np.testing.assert_allclose(self.trainer.sample_depths(frame), [2, 2, 3, np.log10(3e7), 0])
        curve = [(2.0, 0.2), (4.0, 0.8)]
        calls = self.trainer.depth_knob_calls(np.array([0.3, 0.3]), np.array([2.0, 4.0]), curve, 0.5)
        self.assertEqual(list(calls), [True, False])
        self.assertEqual(list(self.trainer.depth_knob_calls(np.array([0.3, 0.6]), np.array([2.0, 4.0]), [], 0.5)),
                         [False, True])

    def test_sparse_depths_join_their_neighbours(self):
        # A point needs DEPTH_KNOB_MIN_SAMPLES samples in its window (6 here, whatever the trainer's): a bin with fewer
        # joins the next deeper one, and bins left at the deep end the point before, so that a few deep samples set no
        # knob of their own.
        self.addCleanup(setattr, self.trainer, "DEPTH_KNOB_MIN_SAMPLES", self.trainer.DEPTH_KNOB_MIN_SAMPLES)
        self.trainer.DEPTH_KNOB_MIN_SAMPLES = 6

        def points(sample_depths):
            samples = np.array([f"s{i}" for i, _ in enumerate(sample_depths) for _ in range(3)])
            depths = np.repeat(np.array(sample_depths, dtype=float), 3)
            windows = self.trainer.depth_knob_windows(depths, samples, np.ones(len(depths), dtype=bool))
            return [(x, len(set(samples[in_group])), len(set(samples[window]))) for x, in_group, window in windows]

        self.assertEqual(points([2.2] * 12 + [3.2] * 12), [(2.2, 12, 12), (3.2, 12, 12)])
        # 3 samples at 1.2 join the 12 at 2.2; 3 at 5.2 and 2 at 5.8 (5 together) join the point at 3.2.
        self.assertEqual(points([1.2] * 3 + [2.2] * 12 + [3.2] * 12 + [5.2] * 3 + [5.8] * 2),
                         [(2.2, 15, 15), (3.2, 17, 17)])
        # 6 deep samples are enough for a point of their own (r226's 10M read pairs: 2 samples of 3 read setups).
        self.assertEqual(points([2.2] * 12 + [5.2] * 6), [(2.2, 12, 12), (5.2, 6, 6)])
        # Too few samples in all: one point.
        self.assertEqual(points([2.2] * 3 + [4.2] * 2), [(2.2, 5, 5)])

    def test_depth_knob_curve_in_the_model(self):
        model, metrics = self.train("knobs", "--depth-knobs")
        curve = read_depth_knob_curve(model)
        # One point per half decade with samples: the shallow ones (~525 fragments) and the deep ones (~52,000).
        self.assertEqual(len(curve), 2, curve)
        self.assertTrue(2.5 < curve[0][0] < 3.0 and 4.5 < curve[1][0] < 5.0, curve)
        self.assertTrue(all(0.05 <= k <= 0.95 for _, k in curve), curve)
        self.assertEqual([tuple(p) for p in metrics["depth_knobs"]["curve"]], curve)
        # The knobs are chosen on the species held out: their F1 there is at least that at the one knob.
        self.assertGreaterEqual(metrics["depth_knobs"]["F1_at_depth_knobs"], metrics["depth_knobs"]["F1_at_knob"])
        # Without --depth-knobs no curve, and without --fdr-calls no calibrated calls.
        plain, metrics = self.train("plain")
        self.assertEqual(read_depth_knob_curve(plain), [])
        self.assertNotIn("depth_knobs", metrics)
        self.assertNotIn("false_calls", metrics)
        self.assertIsNone(read_false_calls(plain))
        # A point keeps the knob unless its best knob gains DEPTH_KNOB_MIN_GAIN on its window.
        gain = self.trainer.DEPTH_KNOB_MIN_GAIN
        for point in metrics_knobs(self, "knobs"):
            if point["knob"] is not None:
                self.assertEqual(point["knob"], point["best knob"] if point["gain"] >= gain else 0.5, point)

    def test_no_depth_knob_curve_with_the_samples_depth_as_a_feature(self):
        # With sample_log_fragments among the features the forest sees the depth itself: no curve is fitted even with
        # --depth-knobs, and the report says why; one knob for every sample instead, if it gains DEPTH_KNOB_MIN_GAIN with
        # species held out, in the model as a curve of one point.
        table = pd.read_csv(self.table, sep="\t")
        totals = table.groupby("meta_sample")["fragments"].transform("sum")
        table["sample_log_fragments"] = np.log10(np.maximum(totals, 1.0))
        path = os.path.join(self.tmp.name, "with_depth.tsv")
        table.to_csv(path, sep="\t", index=False)
        prefix = os.path.join(self.tmp.name, "with_depth")
        result = subprocess.run([sys.executable, os.path.join(HERE, "machine_learning_cmdline.py"), "--truth-file", path,
                                 "--output-prefix", prefix, "--features", "all", "--ntree", "16", "--rounds", "30", "--evaluation", "basic",
                                 "--threads", "1", "--depth-knobs"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:] + result.stderr[-3000:])
        with open(prefix + ".metrics.json") as fh:
            metrics = json.load(fh)
        self.assertIn("depth is a feature", metrics["depth_knobs"]["skipped"])
        self.assertIn("No knob curve: the sample's depth is a feature", result.stdout)
        chosen = metrics["global_knob"]
        curve = read_depth_knob_curve(prefix + ".xml")
        if chosen["knob"] is not None:
            self.assertEqual([k for _, k in curve], [chosen["knob"]])
            self.assertEqual(metrics["depth_knobs"]["global_knob"], chosen["knob"])
            self.assertGreaterEqual(chosen["support"], self.trainer.KNOB_SUPPORT)
        else:
            self.assertEqual(curve, [])
            self.assertIsNone(metrics["depth_knobs"]["global_knob"])
        self.assertEqual([[x, k] for x, k in curve], metrics["depth_knobs"]["curve"])

    def test_one_knob_for_every_sample(self):
        # choose_global_knob: the median of the best thresholds with species held out over bootstrap resamples of the
        # training samples, kept if it beats 0.5 in KNOB_SUPPORT of them and does not lose F1 on the test set, the rows
        # weighted as the models weigh them (the scenarios' by --scenario-weight).
        rng = np.random.default_rng(2)
        n = 4000
        y = (rng.random(n) < 0.3).astype(int)
        scores = np.clip(np.where(y == 1, rng.normal(0.85, 0.08, n), rng.normal(0.45, 0.12, n)), 0, 1)
        frame = pd.DataFrame({"truth": y, "meta_scenario": "", "meta_sample": [f"s{i % 40}" for i in range(n)]})
        opts = types.SimpleNamespace(knob=0.5, seed=1)
        report = self.trainer.Report()
        knob = self.trainer.choose_global_knob(report, frame, y, {"species": scores}, opts)
        self.assertIsNotNone(knob)
        self.assertGreater(knob, 0.55)  # absent taxa score up to ~0.7: a knob above 0.5 calls fewer of them
        chosen = report.data["global_knob"]
        self.assertEqual((chosen["knob"], chosen["median"]), (knob, knob))
        self.assertGreaterEqual(chosen["support"], self.trainer.KNOB_SUPPORT)
        self.assertGreater(chosen["gain"], 0)
        self.assertLessEqual(chosen["range"][0], knob)
        self.assertGreaterEqual(chosen["range"][1], knob)
        self.assertEqual(self.trainer.choose_global_knob(self.trainer.Report(), frame, y, {"species": scores}, opts), knob)
        # A test set it loses on (its present taxa scored between 0.5 and the knob): none, protal calls at --knob.
        test_y = np.ones(50, dtype=int)
        test = (test_y, np.full(50, 0.52), np.ones(50))
        report = self.trainer.Report()
        self.assertIsNone(self.trainer.choose_global_knob(report, frame, y, {"species": scores}, opts, test))
        self.assertEqual(report.data["global_knob"]["test"]["F1_at_default"], 1.0)
        self.assertIn("loses on the test set", "\n".join(report.lines))
        # One it does not lose on: the knob stands.
        self.assertEqual(self.trainer.choose_global_knob(self.trainer.Report(), frame, y, {"species": scores}, opts,
                                                         (test_y, np.full(50, 0.95), np.ones(50))), knob)
        # Scores whose every threshold from 0.1 to 0.9 is as good: no knob beats 0.5, protal calls at --knob.
        even = np.where(y == 1, 0.9, 0.1)
        self.assertIsNone(self.trainer.choose_global_knob(self.trainer.Report(), frame, y, {"species": even}, opts))
        # Without scores with species held out (--evaluation none): none.
        self.assertIsNone(self.trainer.choose_global_knob(self.trainer.Report(), frame, y, {}, opts))
        self.assertEqual(self.trainer.knob_label([(4.0, 0.7)]), "knob 0.7")
        self.assertEqual(self.trainer.knob_label([]), "knob 0.5")
        self.assertTrue(self.trainer.knob_label([(2.0, 0.3), (4.0, 0.7)]).startswith("knob curve (2.000:0.3,"))

    def test_strains(self):
        _, metrics = self.train("plain")
        rows = {(r["simulated from"], r["fragments"]): r for r in metrics["strains"]["by_fragments"]}
        self.assertEqual(rows[("another genome", "all")]["present"], 40 * 4)
        self.assertEqual(rows[("the representative", "all")]["present"], 40 * 4)
        self.assertIn("strains missed", [r["taxa"] for r in metrics["strains"]["features"]])

    def test_conservation_features_by_class(self):
        _, metrics = self.train("plain")
        classes = metrics["feature_classes"]
        self.assertEqual(set(classes), set(self.trainer.TAXON_CLASSES))
        self.assertEqual(classes["absent, congener of a held-out species"]["rows"], 40 * 4)
        self.assertEqual(classes["present, beside a held-out congener"]["rows"], 40 * 2)
        near = classes["absent, congener of a held-out species"]["conserved_fast_record_ratio"][1]
        other = classes["absent, other"]["conserved_fast_record_ratio"][1]
        self.assertGreater(near, other + 1)

    def test_previous_procedure_only_when_asked(self):
        # Off by default. Asked for, its grid takes every second value of max_features (and the last), then the two
        # next to the best, 3 folds and up to PREVIOUS_GRID_ROWS rows; this procedure's side is the evaluation's own
        # forests with species held out.
        _, metrics = self.train("plain")
        self.assertNotIn("previous_procedure", metrics)
        args = ["--truth-file", "t.tsv", "--output-prefix", "p"]
        self.assertFalse(self.trainer.parse_args(args).previous_procedure)
        self.assertFalse(self.trainer.parse_args(args + ["--no-previous-procedure"]).previous_procedure)
        _, metrics = self.train("previous", "--previous-procedure")
        previous = metrics["previous_procedure"]
        self.assertEqual(previous["grid_folds"], 3)
        self.assertEqual(previous["grid_rows"], 40 * 30)
        values = previous["grid_values"]
        coarse = sorted(set(range(values[0], values[-1] + 1, 2)) | {values[-1]})
        self.assertTrue(set(coarse) <= set(values))
        self.assertLessEqual(len(values) - len(coarse), 2)
        self.assertIn(previous["mtry"], values)
        this = next(r for r in previous["results"] if r["procedure"] == "this one")
        self.assertEqual(this["judged on"], "species held out")
        self.assertAlmostEqual(this["AP"], metrics["evaluation"]["species"]["AP"])

    def test_previous_procedure_grid_on_whole_samples(self):
        frame = pd.DataFrame({"meta_sample": [f"s{i // 30}" for i in range(1200)]})
        opts = types.SimpleNamespace(seed=1)
        self.assertEqual(len(self.trainer.grid_rows(frame, opts)), 1200)
        limit = self.trainer.PREVIOUS_GRID_ROWS
        try:
            self.trainer.PREVIOUS_GRID_ROWS = 100
            rows = self.trainer.grid_rows(frame, opts)
        finally:
            self.trainer.PREVIOUS_GRID_ROWS = limit
        self.assertEqual(sorted(frame["meta_sample"].iloc[rows].value_counts().tolist()), [30, 30, 30])  # a 4th: 120


@prerequisites.requires(HAVE_SKLEARN, NEEDS_SKLEARN)
class TrainerFalseCallsTest(unittest.TestCase):
    """machine_learning_cmdline.py --fdr-calls and the singleton rule: the calibration, the prior adjusted to each sample
    and the calls at a target share of false calls, as protal makes them (context::FalseCallKnob; the rule itself on
    the golden vectors, GoldenModelRulesTest)."""

    @classmethod
    def setUpClass(cls):
        import machine_learning_cmdline
        cls.trainer = machine_learning_cmdline
        cls.tmp = tempfile.TemporaryDirectory()
        rng = np.random.default_rng(11)
        rows = []
        for sample in range(40):
            deep = sample % 2 == 1
            absent = 60 if deep else 10  # a deep sample has many more absent candidates
            for taxon in range(8 + absent):
                present = taxon < 8
                fragments = float(rng.integers(500, 3000) if deep else rng.integers(5, 30)) if present else \
                    float(1 if taxon % 3 == 0 else rng.integers(1, 20))
                rows.append({"truth": int(present), "taxon": 100 + (taxon + sample) % 90, "taxon_name": "t",
                             "meta_sample": f"s{sample}", "meta_read_pairs": 100000 if deep else 1000,
                             "fragments": fragments, "x": 1.2 * present + rng.normal(0, 0.6), "y": rng.normal(0, 1),
                             "genus_top_fragments": 500.0 if (not present and taxon % 3 == 0) else 0.0,
                             "em_own_share": 0.3 if (not present and taxon % 3 == 0) else 1.0,
                             "identity": 0.93 if not present else 0.99})
        cls.table = os.path.join(cls.tmp.name, "training.tsv")
        pd.DataFrame(rows).to_csv(cls.table, sep="\t", index=False)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def train(self, name, *extra):
        prefix = os.path.join(self.tmp.name, name)
        result = subprocess.run([sys.executable, os.path.join(HERE, "machine_learning_cmdline.py"), "--truth-file", self.table,
                                 "--output-prefix", prefix, "--features", "all", "--ntree", "16", "--rounds", "30", "--evaluation", "basic",
                                 "--threads", "1", "--test-file", self.table, *extra], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:] + result.stderr[-3000:])
        with open(prefix + ".metrics.json") as fh:
            return prefix + ".xml", json.load(fh), result.stdout

    def test_calls_at_a_target_as_protal(self):
        self.trainer.SINGLETON_CONGENER = 100  # the rule, off by default, vetoes a row here
        self.addCleanup(setattr, self.trainer, "SINGLETON_CONGENER", 0)
        # The calls over a table of two samples: each sample's as when it is alone (the rule, context::FalseCallKnob, on
        # the golden vectors); the singleton rule's taxon never; the counts those of the calls; tied scores together.
        frame = pd.DataFrame({"meta_sample": ["a"] * 6 + ["b"] * 3,
                              "fragments": [5, 5, 5, 5, 1, 5, 5, 5, 5], "genus_top_fragments": [0, 0, 0, 0, 500, 0, 0, 0, 0],
                              "em_own_share": [1, 1, 1, 1, 0.2, 1, 1, 1, 1], "identity": [0.99] * 9})
        p = np.array([0.99, 0.95, 0.95, 0.4, 0.98, 0.1, 0.9, 0.2, 0.05])
        y = np.array([1, 1, 0, 0, 0, 0, 1, 0, 0])
        curve = [(0.0, 0.0), (1.0, 1.0)]
        prepared = self.trainer.FalseCallSamples(frame, y, p, curve, 0.5)
        called = set()
        for fdr in (0.001, 0.03, 0.1, 0.3, 0.6):
            calls = prepared.calls(fdr)
            self.assertFalse(calls[4], "the singleton rule's taxon")
            for sample in ("a", "b"):
                rows = (frame["meta_sample"] == sample).to_numpy()
                alone = self.trainer.FalseCallSamples(frame[rows], y[rows], p[rows], curve, 0.5).calls(fdr)
                np.testing.assert_array_equal(calls[rows], alone, f"fdr {fdr}, sample {sample}")
            tp, fp = prepared.counts(fdr)
            self.assertEqual((tp, fp), (int((calls & (y == 1)).sum()), int((calls & (y == 0)).sum())))
            called.add(int(calls.sum()))
        self.assertGreater(len(called), 2, "the targets call different numbers of taxa")
        self.assertEqual(prepared.calls(0.3)[1], prepared.calls(0.3)[2], "tied scores are called together")

    def test_calibration(self):
        rng = np.random.default_rng(3)
        p = rng.random(5000)
        y = (rng.random(5000) < p ** 2).astype(int)
        curve, prior = self.trainer.fit_calibration(y, p)
        xs, ys = zip(*curve)
        self.assertEqual(xs[0], 0.0)
        self.assertEqual(xs[-1], 1.0)
        self.assertTrue(all(b > a for a, b in zip(xs, xs[1:])))
        self.assertTrue(all(b >= a for a, b in zip(ys, ys[1:])))
        self.assertLessEqual(len(curve), self.trainer.CALIBRATION_POINTS + 2)
        self.assertAlmostEqual(prior, y.mean())
        self.assertAlmostEqual(float(np.interp(0.5, xs, ys)), 0.25, delta=0.05)

    def test_singleton_rule(self):
        # Off by default (protal's --singleton_congener 0). With 100: vetoed is one fragment beside a congener of 100
        # or more, and a read that looks like the congener's (EM share below 0.5 or identity below 0.95); a minor
        # congener's own read (both high) is kept.
        self.assertEqual(self.trainer.SINGLETON_CONGENER, 0)
        frame = pd.DataFrame({"fragments": [1, 1, 2, 1, 1, 1], "genus_top_fragments": [100, 99, 500, 0, 200, 200],
                              "em_own_share": [0.2, 0.2, 0.2, 0.2, 0.98, 0.98], "identity": [0.99, 0.99, 0.99, 0.99, 0.97, 0.93]})
        try:
            self.trainer.SINGLETON_CONGENER = 100
            np.testing.assert_array_equal(self.trainer.vetoed(frame), [True, False, False, False, False, True])
            np.testing.assert_array_equal(self.trainer.call_scores(frame, np.full(6, 0.9)) >= 0.5,
                                          [False, True, True, True, True, False])
            self.trainer.SINGLETON_CONGENER = 0
            self.assertFalse(self.trainer.vetoed(frame).any())
            self.trainer.SINGLETON_CONGENER = 100
            self.assertFalse(self.trainer.vetoed(frame.drop(columns="genus_top_fragments")).any(), "an older dump")
        finally:
            self.trainer.SINGLETON_CONGENER = 0

    def test_fdr_calls_in_the_model(self):
        # (Without --fdr-calls, none: TrainerDepthKnobsTest's plain model.)
        model, metrics, stdout = self.train("fdr", "--fdr-calls", "--depth-knobs", "--singleton-congener", "100")
        calls = read_false_calls(model)
        self.assertIsNotNone(calls)
        self.assertEqual(calls["fdr"], metrics["false_calls"]["fdr"])
        self.assertIn(calls["fdr"], [float(v) for v in self.trainer.FDR_GRID])
        self.assertEqual([tuple(c) for c in metrics["false_calls"]["curve"]], calls["curve"])
        self.assertIn("test_false_calls", metrics)
        self.assertIn("singleton rule (one fragment beside a congener of 100 or more", stdout)
        self.assertIn("FP fdr", metrics["test_by_depth"][0])


@prerequisites.requires(HAVE_SKLEARN, NEEDS_SKLEARN)
class TrainerScenariosTest(unittest.TestCase):
    """machine_learning_cmdline.py on tables with scenarios (meta_scenario): the hold-in rows of the training table and
    the hold-out rows of the test table reported per scenario, apart from the independent test set; and --features
    auto, which chooses among the named sets with species held out. Each training runs once for the class (train)."""

    @classmethod
    def setUpClass(cls):
        import model_features
        cls.features = model_features
        cls.tmp = tempfile.TemporaryDirectory()
        cls.trained = {}
        columns = [c for name in model_features.AUTO_CANDIDATES for c in model_features.feature_set_columns(name)]
        columns = list(dict.fromkeys(columns))

        def table(path, samples, seed, scenario_of):
            rng = np.random.default_rng(seed)
            rows = []
            for sample in range(samples):
                scenario = scenario_of(sample)
                for taxon in range(20):
                    present = taxon < 6
                    row = {c: rng.normal(0, 1) for c in columns}
                    # the signal in a normalized feature and in the depth feature's group, so that the sets differ
                    row["hit_gene_fraction"] = 1.2 * present + rng.normal(0, 0.6)
                    row["sample_log_fragments"] = 0.8 * present + rng.normal(0, 0.6)
                    row.update({"truth": int(present), "taxon": 100 + (taxon + 3 * sample) % 50, "taxon_name": "t",
                                "meta_sample": f"{scenario or 'design'}{seed}_{sample}", "fragments": 10.0,
                                "meta_scenario": scenario, "meta_novel_level": "species" if taxon in (7, 8) else "",
                                "identity": 0.98})
                    rows.append(row)
            pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
        cls.training = os.path.join(cls.tmp.name, "training.tsv")
        cls.test = os.path.join(cls.tmp.name, "test.tsv")
        # Training: 30 design samples, 6 of scenario gut. Test: 8 design samples, 4 of gut, 3 of host (hold-out only).
        table(cls.training, 36, 1, lambda s: "gut" if s >= 30 else "")
        table(cls.test, 15, 2, lambda s: "gut" if 8 <= s < 12 else "host" if s >= 12 else "")

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def train(self, name, features, test=None, *extra):
        """The trainer on the training table, with --features `features` (None: its default), once per name ->
        (prefix, metrics, console)."""
        if name not in self.trained:
            prefix = os.path.join(self.tmp.name, name)
            result = subprocess.run([sys.executable, os.path.join(HERE, "machine_learning_cmdline.py"), "--truth-file",
                                     self.training, "--output-prefix", prefix, *(["--features", features] if features else []),
                                     "--ntree", "16", "--rounds", "30", "--evaluation", "basic", "--threads", "1",
                                     "--test-file", test or self.test, *extra], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout[-3000:] + result.stderr[-3000:])
            with open(prefix + ".metrics.json") as fh:
                self.trained[name] = (prefix, json.load(fh), result.stdout)
        return self.trained[name]

    def test_scenarios_hold_in_and_hold_out(self):
        prefix, metrics, stdout = self.train("normalized", "normalized")
        rows = {(r["scenario"], r["set"]): r for r in metrics["scenarios"]}
        self.assertEqual(set(rows), {("(design)", "training, species held out"), ("gut", "hold-in, in sample"),
                                     ("gut", "hold-in, samples held out"), ("gut", "hold-in, species held out"),
                                     ("gut", "hold-out"), ("host", "hold-out")})
        self.assertEqual((rows[("gut", "hold-in, in sample")]["samples"], rows[("gut", "hold-out")]["samples"],
                          rows[("host", "hold-out")]["samples"]), (6, 4, 3))
        for r in metrics["scenarios"]:
            self.assertEqual(r["TP"] + r["FN"], r["present"])
            self.assertAlmostEqual(r["FN rate %"], 100 * r["FN"] / r["present"])
            self.assertAlmostEqual(r["FP rate %"], 100 * r["FP"] / r["absent"])
            self.assertLessEqual(r["FP near novel"], r["FP"])
        self.assertEqual(rows[("gut", "hold-out")]["present"], 4 * 6)
        # The independent test set is the test table's design samples alone; the scenarios' rows are scored apart.
        self.assertEqual(metrics["test"]["this one"]["taxa"], 8 * 20)
        self.assertIn("## Scenarios: hold-in and hold-out samples", stdout)
        self.assertRegex(stdout, r"scenario gut: hold-out F1 [0-9.]+ \(FP rate [0-9.]+%, FN rate [0-9.]+%, 4 samples\); "
                                 r"hold-in with species held out F1")
        self.assertRegex(stdout, r"scenario host: hold-out F1 [0-9.]+")
        predictions = pd.read_csv(prefix + ".scenario_predictions.tsv.gz", sep="\t")
        self.assertEqual(sorted(predictions["meta_scenario"].unique()), ["gut", "host"])
        self.assertEqual(len(predictions), 7 * 20)
        # Every row's call at the model's knob (error_reads.py): the training rows' with species held out, the test
        # table's by the final model; the scenario table's errors are theirs.
        calls = pd.read_csv(prefix + ".calls.tsv.gz", sep="\t")
        self.assertEqual(calls.groupby("set").size().to_dict(), {"training": 36 * 20, "test": 15 * 20})
        self.assertTrue(((calls["p"] >= calls["knob"]) == (calls["call"] == 1)).all())
        held = calls[(calls["set"] == "test") & calls["meta_scenario"].notna()]
        self.assertTrue(np.allclose(held["p"].to_numpy(), predictions["p"].to_numpy(), atol=1e-5))
        tested = pd.read_csv(prefix + ".test_predictions.tsv.gz", sep="\t")
        design = calls[(calls["set"] == "test") & calls["meta_scenario"].isna()]
        self.assertTrue(np.allclose(design["p"].to_numpy(), tested["p"].to_numpy(), atol=1e-5))
        for (scenario, kind), group in calls[calls["meta_scenario"].notna()].groupby(["meta_scenario", "set"]):
            row = rows[(scenario, "hold-in, species held out" if kind == "training" else "hold-out")]
            self.assertEqual((int(((group["call"] == 1) & (group["truth"] == 0)).sum()),
                              int(((group["call"] == 0) & (group["truth"] == 1)).sum())), (row["FP"], row["FN"]))
        # The scenarios' rows weigh 0.25 in every fit by default (--scenario-weight); at 1 the model is another.
        self.assertEqual(metrics["scenario_weight"], {"weight": 0.25, "rows": 6 * 20})
        self.assertIn("the scenarios' 120 of 720 rows weigh 0.25 in every model fitted", stdout)
        full, _, _ = self.train("scenarios_full_weight", "normalized", None, "--scenario-weight", "1")
        again = pd.read_csv(full + ".scenario_predictions.tsv.gz", sep="\t")
        self.assertFalse(np.allclose(again["p"], predictions["p"]))
        # A test table of the scenarios' hold-out samples alone: no independent test set, the scenarios reported.
        only = os.path.join(self.tmp.name, "only_scenarios.tsv")
        frame = pd.read_csv(self.test, sep="\t")
        frame[frame["meta_scenario"].fillna("") != ""].to_csv(only, sep="\t", index=False)
        _, metrics, stdout = self.train("only", "normalized", only)
        self.assertNotIn("test", metrics)
        self.assertIn("holds the scenarios' hold-out samples only", stdout)
        self.assertIn(("host", "hold-out"), {(r["scenario"], r["set"]) for r in metrics["scenarios"]})

    def test_boosting_by_default_and_the_forest_on_request(self):
        # The default model is gradient-boosted trees (--model gbm), written as a chain; --model forest the forest,
        # averaged; each scores in PMML as scikit-learn does (the trainer's own check), with its importances.
        prefix, metrics, stdout = self.train("normalized", "normalized")
        self.assertEqual(metrics["model"]["kind"], "gbm")
        self.assertEqual(metrics["model"]["trees"], 30)
        self.assertEqual((metrics["model"]["pmml_vs_sklearn_max_diff"], metrics["model"]["pmml_vs_sklearn_call_differences"]),
                         (0.0, 0))
        self.assertIn("gradient-boosted trees: 30 rounds at a learning rate of 0.1, max 63 leaves, min leaf 20", stdout)
        self.assertIsInstance(load_model(prefix + ".xml"), PmmlBoosted)
        self.assertNotIn("out of bag", metrics["evaluation"])
        self.assertNotIn("genus", metrics["evaluation"])  # the clades held out only with --evaluation full
        self.assertIn("cross-validation folds fitted 1 at a time, on 1 threads each", stdout)
        importance = pd.read_csv(prefix + ".varimp.tsv", sep="\t")
        self.assertAlmostEqual(importance["importance"].sum(), 1.0, places=4)  # written with 6 decimals
        # The folds side by side in worker processes (--fold-jobs): the same scores as one after another.
        side, side_metrics, side_stdout = self.train("boosted_side_by_side", "normalized", None, "--fold-jobs", "2")
        self.assertIn("cross-validation folds fitted 2 at a time, on 1 threads each", side_stdout)
        self.assertEqual(side_metrics["evaluation"], metrics["evaluation"])
        pd.testing.assert_frame_equal(pd.read_csv(side + ".predictions.tsv.gz", sep="\t"),
                                      pd.read_csv(prefix + ".predictions.tsv.gz", sep="\t"))
        prefix, metrics, stdout = self.train("forest", "normalized", None, "--model", "forest")
        self.assertEqual((metrics["model"]["kind"], metrics["model"]["trees"]), ("forest", 16))
        self.assertEqual(metrics["model"]["pmml_vs_sklearn_max_diff"], 0.0)
        self.assertIn("random forest: 16 trees, max 256 leaves, min leaf 1", stdout)
        self.assertIsInstance(load_model(prefix + ".xml"), PmmlForest)
        self.assertIn("out of bag", metrics["evaluation"])

    def test_features_auto(self):
        # (auto is the trainer's default: FeatureSetsTest.)
        prefix, metrics, stdout = self.train("auto", "auto")
        auto = metrics["features_auto"]
        scored = [r["features"] for r in auto["candidates"]]
        self.assertEqual(scored, list(self.features.AUTO_CANDIDATES))  # the table has every candidate's features
        self.assertFalse(any("priors" in name for name in scored))
        f1 = {r["features"]: r["F1"] for r in auto["candidates"]}
        best, default = max(f1, key=f1.get), self.features.DEFAULT_FEATURE_SET
        self.assertEqual(auto["chosen"], best if f1[best] >= f1[default] + self.features.AUTO_MIN_GAIN else default)
        # Each candidate's forest also scored the samples never trained on, every test set apart, for comparison.
        self.assertTrue(all({"F1 test set", "F1 gut hold-out", "F1 host hold-out"} <= set(r) for r in auto["candidates"]))
        others = [r for r in auto["candidates"] if r["features"] != auto["chosen"]]
        self.assertAlmostEqual(auto["others_mean_F1"], float(np.mean([r["F1"] for r in others])))
        self.assertEqual(set(auto["held_out"]), {"test set", "gut hold-out", "host hold-out"})
        self.assertAlmostEqual(auto["held_out"]["host hold-out"]["others_mean_F1"],
                               float(np.mean([r["F1 host hold-out"] for r in others])))
        # Why, in a line: its F1 against the other sets', the rule, the test sets.
        self.assertRegex(auto["why"], rf"^F1 [0-9.]+ with species held out, the other {len(others)} sets [0-9.]+ on average "
                                      r"\(the best of them \S+ [0-9.]+\); (the highest|the default set, as no other is "
                                      rf"{re.escape(str(self.features.AUTO_MIN_GAIN))} better)[^;]*; on samples never trained "
                                      r"on \(not used to choose\): test set [0-9.]+ \(the others [0-9.]+\), gut hold-out "
                                      r"[0-9.]+ \(the others [0-9.]+\), host hold-out")
        # The model takes the chosen set's features, and the evaluation scored that set.
        self.assertEqual(load_model(prefix + ".xml").features, self.features.feature_set_columns(auto["chosen"]))
        self.assertIn("## Feature set chosen (--features auto, species held out)", stdout)
        self.assertIn(f"feature set (--features auto): {auto['chosen']}: {auto['why']}", stdout)
        self.assertAlmostEqual(metrics["evaluation"]["species"]["F1"], f1[auto["chosen"]])


@prerequisites.requires(HAVE_SKLEARN, NEEDS_SKLEARN)
class GoldenModelRulesTest(unittest.TestCase):
    """The trainer's side of the rules it shares with protal, on the golden vectors of tests/data/golden_model_rules.tsv
    (its header says how they were made); tests/test_GoldenModelRules.cpp checks protal's side on the same lines: the
    prior adjusted to a sample (sample_adjusted, context::SampleAdjusted), the calls at a target share of false calls
    (FalseCallSamples, context::FalseCallKnob) and the knob curve over the sample's depth (knob_at at sample_depths,
    profiler::DepthKnobAt)."""

    @classmethod
    def setUpClass(cls):
        import machine_learning_cmdline
        cls.trainer = machine_learning_cmdline

    def test_the_file_holds_every_rule(self):
        for rule, least in (("SampleAdjusted", 5), ("FalseCallKnob", 8), ("DepthKnobAt", 8)):
            self.assertGreaterEqual(len(golden_cases(rule)), least, rule)

    def test_sample_adjusted(self):
        for name, f in golden_cases("SampleAdjusted"):
            adjusted, rate = self.trainer.sample_adjusted(golden_numbers(f["q"]), float(f["prior"]))
            self.assertAlmostEqual(rate, float(f["rate"]), delta=1e-9, msg=name)
            np.testing.assert_allclose(adjusted, golden_numbers(f["adjusted"]), rtol=0, atol=1e-9, err_msg=name)

    def test_false_call_knob(self):
        for name, f in golden_cases("FalseCallKnob"):
            scores = np.array(golden_numbers(f["scores"]))
            prepared = self.trainer.FalseCallSamples(pd.DataFrame({"meta_sample": ["s"] * len(scores)}),
                                                     np.zeros(len(scores), dtype=int), scores, golden_curve(f["curve"]),
                                                     float(f["prior"]))
            calls = prepared.calls(float(f["fdr"]))
            self.assertEqual(int(calls.sum()), int(f["called"]), name)
            if f["knob"] == "none":
                self.assertFalse(calls.any(), name)
            else:  # the lowest score called, every taxon at or above it called
                self.assertEqual(float(scores[calls].min()), float(f["knob"]), name)
                np.testing.assert_array_equal(calls, scores >= float(f["knob"]), name)
            if len(scores):
                self.assertAlmostEqual(prepared.samples[0]["prior"], float(f["rate"]), delta=1e-9, msg=name)

    def test_depth_knob_at(self):
        for name, f in golden_cases("DepthKnobAt"):
            curve, knob = golden_curve(f["curve"]), float(f["knob"])
            depth = self.trainer.sample_depths(pd.DataFrame({"meta_sample": ["s"], "fragments": [float(f["fragments"])]}))
            if curve:
                self.assertAlmostEqual(float(self.trainer.knob_at(curve, depth)[0]), knob, delta=1e-12, msg=name)
            else:  # without a curve, protal calls at --knob's 0.5, as the trainer does
                calls = self.trainer.depth_knob_calls(np.array([knob, np.nextafter(knob, 0)]), np.repeat(depth, 2), curve, 0.5)
                self.assertEqual(list(calls), [True, False], name)


if __name__ == "__main__":
    unittest.main()
