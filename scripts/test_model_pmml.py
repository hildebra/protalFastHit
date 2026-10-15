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

import inspect
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
from model_pmml import (PmmlBoosted, PmmlForest, feature_ranges_of, float32_split, format_depth_knob_curve,  # noqa: E402
                        format_depth_knobs, format_feature_ranges, format_feature_versions, load_model, model_kind,
                        read_depth_knob_curve, read_depth_knobs, read_false_calls, read_feature_ranges,
                        read_feature_versions, read_feature_versions_file, write_forest, write_model, write_placeholder)

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

    def test_feature_versions_and_ranges(self):
        # Since 2026-10-15 (the audit's C1, C4): the features' versions (protal --feature_versions) and their ranges over
        # the rows the model was fitted on, as header Extensions before the Application, each number "%.9g"; read back;
        # the scores as without them.
        rf = RandomForestClassifier(n_estimators=2, random_state=1).fit(self.X, self.y)
        path = os.path.join(self.tmp.name, "features.xml")
        ranges = feature_ranges_of(self.X, self.features)
        write_forest(rf, self.features, path, feature_versions={"b": "3", "a": "1.2", "other": "9"}, feature_ranges=ranges)
        with open(path) as fh:
            xml = fh.read()
        self.assertIn('  <Extension name="protal_feature_versions" value="a:1.2,b:3"/>\n', xml)  # the model's, in its order
        value = re.search(r'<Extension name="protal_feature_ranges" value="([^"]*)"/>\n  <Application ', xml).group(1)
        items = value.split(",")
        self.assertEqual([item.split(":")[0] for item in items], self.features)
        for item, column in zip(items, self.X.T):
            expected = (column.min(), np.quantile(column, 0.01), np.quantile(column, 0.99), column.max())
            self.assertEqual(item.split(":")[1:], ["%.9g" % v for v in expected])
        self.assertEqual(read_feature_versions(path), {"a": "1.2", "b": "3"})
        read = read_feature_ranges(path)
        self.assertEqual(list(read), self.features)
        for f in self.features:
            np.testing.assert_allclose(read[f], ranges[f], rtol=1e-8)
        np.testing.assert_array_equal(PmmlForest(path).predict(self.X), rf.predict_proba(self.X)[:, 1])
        self.assertEqual(model_kind(path), "forest")
        # Without them none; with versions of none of the model's features no versions.
        write_forest(rf, self.features, path, feature_versions={"other": "1"})
        self.assertEqual((read_feature_versions(path), read_feature_ranges(path)), ({}, {}))
        with open(path) as fh:
            self.assertNotIn("protal_feature", fh.read())

    def test_rejects_other_classes(self):
        rf = RandomForestClassifier(n_estimators=2, random_state=1).fit(self.X, np.where(self.y == 1, "yes", "no"))
        with self.assertRaises(ValueError):
            write_forest(rf, self.features, os.path.join(self.tmp.name, "bad.xml"))


class FeatureExtensionsTest(unittest.TestCase):
    """The feature extensions' text (model_pmml.py): the ranges as C's "%.9g", names of letters, digits and "_" only, the
    table of protal --feature_versions read, and placeholder models without them."""

    def test_ranges_text(self):
        self.assertEqual(format_feature_ranges({"x": (-0.0, 1e-5, 0.1 + 0.2, 123456789012.0), "y_2": (1, 2, 3, 4)},
                                               ["y_2", "x"]),
                         "y_2:1:2:3:4,x:-0:1e-05:0.3:1.23456789e+11")
        for ranges, features in (({"a:b": (0, 0, 0, 0)}, ["a:b"]), ({"a b": (0, 0, 0, 0)}, ["a b"]),
                                 ({"a": (0, 0, 0)}, ["a"]), ({"a": (0, 0, 0, np.inf)}, ["a"]), ({"a": (0, 0, 0, np.nan)}, ["a"]),
                                 ({}, ["a"])):
            with self.assertRaises(ValueError, msg=str(ranges)):
                format_feature_ranges(ranges, features)
        self.assertEqual(format_feature_versions({"b": "2", "a": "x1"}, ["a", "c", "b"]), "a:x1,b:2")
        with self.assertRaises(ValueError):
            format_feature_versions({"a": "1,2"}, ["a"])
        with self.assertRaises(ValueError):
            feature_ranges_of(np.zeros((3, 2)), ["a"])

    def test_versions_file(self):
        with tempfile.TemporaryDirectory() as root:
            path = os.path.join(root, "versions.tsv")

            def read(text):
                with open(path, "w") as fh:
                    fh.write(text)
                return read_feature_versions_file(path)
            self.assertEqual(read("feature\tversion\tchanged\tnote\nfragments\t1\t2026-08-10\tfirst\n\n"
                                  "identity\t2\t2026-10-15\ta note, with commas: and colons\n# a comment\n"),
                             {"fragments": "1", "identity": "2"})
            self.assertEqual(read("#feature\tversion\tchanged\tnote\r\ndepth\t3\t\t\r\n"), {"depth": "3"})
            self.assertEqual(read("note\tversion\tfeature\nx\t4\tdepth\n"), {"depth": "4"})  # by the header's names
            for bad in ("", "name\tversion\nx\t1\n", "feature\tversion\nx\t1,2\n", "feature\tversion\na-b\t1\n",
                        "feature\tversion\nx\t1\nx\t2\n", "feature\tversion\nx\n"):
                with self.assertRaises(ValueError, msg=repr(bad)):
                    read(bad)
            self.assertEqual(read("feature\tversion\nx\t1\nx\t1\n"), {"x": "1"})

    def test_placeholder_models_have_none(self):
        with tempfile.TemporaryDirectory() as root:
            path = os.path.join(root, "model_pe.xml")
            write_placeholder(path, "pe")
            self.assertEqual(model_kind(path), "placeholder")
            self.assertEqual((read_feature_ranges(path), read_feature_versions(path)), ({}, {}))


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
        ranges = feature_ranges_of(self.X, self.features)
        write_model(model, self.features, path, depth_knob_curve=[(5.5, 0.7)], false_calls=calls,
                    feature_versions={f: "1" for f in self.features}, feature_ranges=ranges)
        self.assertEqual(read_depth_knob_curve(path), [(5.5, 0.7)])
        self.assertEqual(read_false_calls(path), calls)
        self.assertEqual(read_feature_versions(path), {f: "1" for f in self.features})
        self.assertEqual(list(read_feature_ranges(path)), self.features)
        np.testing.assert_array_equal(load_model(path).predict(self.X), model.predict_proba(self.X)[:, 1])
        self.assertEqual(model_kind(path), "gbm")

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
                   mf.ALLELE_FEATURES + mf.POLYMORPHIC_FEATURES + mf.WEIGHT_FEATURES + mf.STRAIN_IDENTITY_FEATURES + mf.PRIORS_FEATURES +
                   mf.LENGTH_FREE_FEATURES + mf.READ_LENGTH_FEATURES + ["genus_top_fragments", "other"])
        # The priors are opt-in: their gain at r226 is the cluster-size rule the simulation cannot test. The reference's
        # k-mer uniqueness (ref) and the sample's complexity are in the default set since 2026-10-06, the groups against
        # the false positives of complex communities (consistency, shape, neighbourhood) since 2026-10-07, the ancestry
        # sites since 0.7.9 (2026-10-07), and the congener gaps and the untried candidates (gaps, untried) since the merge of
        # congener-gaps (2026-10-07), the strain alleles since 2026-10-08, the polymorphic sites since 2026-10-09, the column
        # weights since 2026-10-10 and the identities against the best strain allele since 2026-10-13. The foreign
        # rates are not in the default set: their first scan read the genomes
        # the samples are drawn from, which they told the models (r226 v17, 2026-10-08); the scan of the full
        # reference since then is offered to --features auto by one named set.
        # The default of 2026-10-13, and since 2026-10-14 with the features that do not change with the read length
        # (length_free, beside or in place of those that do) and the sample's read length.
        previous = ("normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood+"
                    "ancestry+gaps+untried+alleles+polymorphic+weights+strain_identity")
        self.assertEqual(mf.DEFAULT_FEATURE_SET, previous + "+length_free+read_length")
        self.assertIn(previous, mf.AUTO_CANDIDATES)
        self.assertEqual([named for named in mf.FEATURE_SETS if "foreign" in named.split("+")],
                         [mf.DEFAULT_FEATURE_SET + "+foreign"])
        self.assertIn(mf.DEFAULT_FEATURE_SET + "+foreign", mf.AUTO_CANDIDATES)
        self.assertEqual(len(mf.FOREIGN_FEATURES), 3)
        # allele_copy_share is a column of the dump in no group: it tells a species with other genomes in GTDB from one
        # without, the cluster size the priors carry.
        self.assertEqual(mf.ALLELE_FEATURES, ["allele_explained_share", "allele_identity_gain"])
        self.assertNotIn("allele_copy_share", [c for group in mf.FEATURE_GROUPS.values() for c in group])
        # So are the polymorphic sites per kb and the share of the ancestry sites fixed within the species.
        self.assertEqual(mf.POLYMORPHIC_FEATURES, ["polymorphic_known_share", "polymorphic_novel_share", "ancestry_fixed_gain",
                                                   "ancestry_fixed_agreement"])
        for column in ("allele_sites_per_kb", "ancestry_fixed_share"):
            self.assertNotIn(column, [c for group in mf.FEATURE_GROUPS.values() for c in group])
        # The column weights (2026-10-10): five features; column_weight_coverage is in no group.
        self.assertEqual(mf.WEIGHT_FEATURES, ["conserved_mismatch_ratio", "conserved_mismatch_rate", "ancestry_agreement_weighted",
                                              "nonsynonymous_share", "nonsynonymous_conserved_rate"])
        self.assertNotIn("column_weight_coverage", [c for group in mf.FEATURE_GROUPS.values() for c in group])
        # The identities against the best strain allele (2026-10-13): two features.
        self.assertEqual(mf.STRAIN_IDENTITY_FEATURES, ["allele_identity", "allele_top_identity"])
        copies = (mf.GAP_FEATURES + mf.UNTRIED_FEATURES + mf.ALLELE_FEATURES + mf.POLYMORPHIC_FEATURES + mf.WEIGHT_FEATURES +
                  mf.STRAIN_IDENTITY_FEATURES)
        self.assertEqual(len(copies), 18)
        self.assertEqual(len(set(copies)), 18)
        self.assertEqual(mf.REF_FEATURES, ["su_rate_ref", "lu_rate_ref", "lsu_rate_ref"])
        self.assertEqual(mf.COMPLEXITY_FEATURES, ["sample_log_taxa", "sample_low_identity", "sample_identity"])
        self.assertEqual(mf.ANCESTRY_FEATURES, ["ancestry_sites_per_record", "ancestry_agreement", "ancestry_congener_share",
                                                "ancestry_indel_sites_per_record", "ancestry_indel_congener_share"])
        new = mf.CONSISTENCY_FEATURES + mf.SHAPE_FEATURES + mf.NEIGHBOURHOOD_FEATURES + mf.ANCESTRY_FEATURES
        self.assertEqual(len(new), 20)
        self.assertEqual(len(set(new)), 20)
        old_columns = (mf.NORMALIZED_FEATURES + mf.ADJACENCY_FEATURES + mf.DISTANCE_FEATURES + mf.SAMPLE_FEATURES +
                       mf.DIVERGENCE_FEATURES + mf.UNFILTERED_FEATURES + mf.REF_FEATURES + mf.COMPLEXITY_FEATURES + new + copies)
        self.assertEqual(mf.feature_columns(columns, previous), old_columns)
        # length_free puts the unique seed enrichments after the last of the eight unique k-mer rates, which stay (since
        # the counts come from the alignment, 2026-10-14), and its other features where the first of those they replace
        # were, those gone; read_length adds its three.
        replaced = {"top_identity": ["top_identity_bound"], "allele_top_identity": ["allele_top_identity_bound"],
                    "ancestry_sites_per_record": ["ancestry_sites_per_kb", "ancestry_indel_sites_per_kb"]}
        gone = {c for olds, _ in mf.LENGTH_FREE_REPLACEMENTS for c in olds}
        rates, enrichments = mf.LENGTH_FREE_ADDITIONS[0]
        self.assertEqual(enrichments, ("unique_seed_enrichment", "unique2_seed_enrichment"))

        def with_length_free(base):
            last = max((i for i, c in enumerate(base) if c in rates), default=None)
            return [n for i, c in enumerate(base)
                    for n in (replaced.get(c, []) if c in gone else [c] + (list(enrichments) if i == last else []))]

        expected = with_length_free(old_columns) + mf.READ_LENGTH_FEATURES
        self.assertEqual(mf.feature_columns(columns, mf.DEFAULT_FEATURE_SET), expected)
        self.assertEqual(len(gone), 4)
        self.assertEqual(len(expected), 98)
        self.assertEqual(expected.index("unique_seed_enrichment"), expected.index("lsu_gene_rate3") + 1)
        for column in ("top_identity", "ancestry_indel_sites_per_record", "allele_top_identity"):
            self.assertNotIn(column, expected)
        for column in ("lu_gene_rate3", "lsu_per_kb", "unique2_seed_enrichment", "top_identity_bound"):
            self.assertIn(column, expected)
        # The shares of the anchors' seeds that are unique are columns of the dump in no group.
        for column in ("unique_seed_share", "unique2_seed_share"):
            self.assertNotIn(column, [c for group in mf.FEATURE_GROUPS.values() for c in group])
        # length_free alone adds nothing the set lacks the old features of.
        self.assertEqual(mf.feature_columns(columns, "normalized+length_free"), with_length_free(mf.NORMALIZED_FEATURES))
        # A table of a protal before 2026-10-14: a clear error with the default, the previous one works.
        before_length_free = [c for c in columns if c not in mf.LENGTH_FREE_FEATURES + mf.READ_LENGTH_FEATURES]
        with self.assertRaisesRegex(RuntimeError, "unique_seed_enrichment"):
            mf.feature_columns(before_length_free, mf.DEFAULT_FEATURE_SET)
        self.assertEqual(mf.feature_columns(before_length_free, previous), old_columns)
        self.assertEqual(mf.feature_columns(columns, previous + "+priors"),
                         mf.NORMALIZED_FEATURES + mf.ADJACENCY_FEATURES + mf.DISTANCE_FEATURES + mf.SAMPLE_FEATURES +
                         mf.DIVERGENCE_FEATURES + mf.UNFILTERED_FEATURES + mf.REF_FEATURES + mf.COMPLEXITY_FEATURES +
                         new + copies + mf.PRIORS_FEATURES)
        # A table of a protal before the per-copy tables: a clear error with the default, the set without them works.
        before = [c for c in columns if c not in copies]
        with self.assertRaisesRegex(RuntimeError, "gap_informative_share"):
            mf.feature_columns(before, previous)
        self.assertEqual(mf.feature_columns(before, "normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+"
                                                    "consistency+shape+neighbourhood+ancestry"),
                         mf.NORMALIZED_FEATURES + mf.ADJACENCY_FEATURES + mf.DISTANCE_FEATURES + mf.SAMPLE_FEATURES +
                         mf.DIVERGENCE_FEATURES + mf.UNFILTERED_FEATURES + mf.REF_FEATURES + mf.COMPLEXITY_FEATURES + new)
        self.assertIn("normalized+adjacency+distance+depth+divergence+unfiltered", mf.AUTO_CANDIDATES)  # an old default
        # A table of a protal of 2026-10-08 (the strain alleles, before the polymorphic sites): a clear error with the default,
        # the set without them works.
        without_strain_identity = previous.removesuffix("+strain_identity")
        alleles_only = [c for c in columns if c not in mf.POLYMORPHIC_FEATURES + mf.WEIGHT_FEATURES + mf.STRAIN_IDENTITY_FEATURES]
        with self.assertRaisesRegex(RuntimeError, "polymorphic_known_share"):
            mf.feature_columns(alleles_only, previous)
        self.assertIn(without_strain_identity.removesuffix("+weights").removesuffix("+polymorphic"), mf.AUTO_CANDIDATES)
        # A table of a protal before the column weights (2026-10-10): the set without them is named too.
        before_weights = [c for c in columns if c not in mf.WEIGHT_FEATURES + mf.STRAIN_IDENTITY_FEATURES]
        with self.assertRaisesRegex(RuntimeError, "conserved_mismatch_ratio"):
            mf.feature_columns(before_weights, previous)
        self.assertIn(without_strain_identity.removesuffix("+weights"), mf.AUTO_CANDIDATES)
        # A table of a protal before the identities against the best strain allele (2026-10-13, the r226 builds up to v27):
        # a clear error with the default, the set without them works.
        before_strain_identity = [c for c in columns if c not in mf.STRAIN_IDENTITY_FEATURES]
        with self.assertRaisesRegex(RuntimeError, "allele_top_identity"):
            mf.feature_columns(before_strain_identity, previous)
        self.assertIn(without_strain_identity, mf.AUTO_CANDIDATES)
        self.assertEqual(mf.feature_columns(before_strain_identity, without_strain_identity),
                         mf.feature_columns(columns, previous)[:-len(mf.STRAIN_IDENTITY_FEATURES)])
        # A table of the r226 v18 build (before the strain alleles): a clear error with the default, the set without them works.
        v18 = [c for c in columns if c not in mf.ALLELE_FEATURES + mf.POLYMORPHIC_FEATURES + mf.WEIGHT_FEATURES +
               mf.STRAIN_IDENTITY_FEATURES]
        with self.assertRaisesRegex(RuntimeError, "allele_explained_share"):
            mf.feature_columns(v18, previous)
        self.assertIn("normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood+"
                      "ancestry+gaps+untried", mf.AUTO_CANDIDATES)
        self.assertEqual(mf.feature_columns(v18, "normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+"
                                                 "consistency+shape+neighbourhood+ancestry+gaps+untried"),
                         mf.NORMALIZED_FEATURES + mf.ADJACENCY_FEATURES + mf.DISTANCE_FEATURES + mf.SAMPLE_FEATURES +
                         mf.DIVERGENCE_FEATURES + mf.UNFILTERED_FEATURES + mf.REF_FEATURES + mf.COMPLEXITY_FEATURES + new +
                         mf.GAP_FEATURES + mf.UNTRIED_FEATURES)
        # A table of the r226 v15 build (before the false-positive groups): a clear error with the default, its old set works.
        v15 = [c for c in columns if c not in new and c not in copies]
        with self.assertRaisesRegex(RuntimeError, "read_consensus_share"):
            mf.feature_columns(v15, previous)
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

    def test_points_of_the_read_type(self):
        """With single-end setups of their own (--se_read_setups) a design point holds paired-end samples only,
        single-end ones only, or both: the check takes its points among those with samples of its reads (r226 v28
        took rl100 and rl75 for se, which had none)."""
        import check_model_parity as cmp
        with tempfile.TemporaryDirectory() as training:
            def point(name, pe, se):
                base = os.path.join(training, "points", name)
                os.makedirs(os.path.join(base, "sim"))
                sample = name + "_s_1"
                with open(os.path.join(base, "sim", "protal.meta"), "w") as fh:
                    fh.write(f"#OUTPUT_DIR\t{os.path.join(base, 'protal')}\n#SAMPLEID\tSAM\tPROFILE_TRUTH\n"
                             f"{sample}\t{sample}.sam.gz\t/truth\n")
                for wanted, folder, name_of in ((pe, "protal", sample), (se, "protal_se", sample + "_se")):
                    if wanted:
                        for sub, ending in (("alignments", ".sam.gz"), ("profiles", ".profile.truth_annotated")):
                            os.makedirs(os.path.join(base, folder, sub), exist_ok=True)
                            open(os.path.join(base, folder, sub, name_of + ending), "w").close()
            point("rl100_p10", pe=True, se=False)
            point("rl150_p10", pe=True, se=True)
            point("rl250_p10", pe=False, se=True)
            point("sc_gut_community", pe=True, se=True)
            names = {t: [(os.path.basename(p), [s[0] for s in samples]) for p, samples in cmp.points_with_samples(training, t)]
                     for t in ("pe", "se")}
            self.assertEqual(names["pe"], [("rl100_p10", ["rl100_p10_s_1"]), ("rl150_p10", ["rl150_p10_s_1"])])
            self.assertEqual(names["se"], [("rl150_p10", ["rl150_p10_s_1_se"]), ("rl250_p10", ["rl250_p10_s_1_se"])])
            self.assertEqual(cmp.points_with_samples(training, "pb"), [])

    @staticmethod
    def parity_layout(training, simulated, drawn, levels=()):
        """collect_training_data.py's layout of design points of one sample each, with its SAM and training dump:
        simulated {name: (pe, se)}, a simulator's map (the paired-end sample's SAM <sample>.sam.zst in protal/, its
        first reads' <sample>_se.sam.gz in protal_se/); drawn, the names of drawn units (sim/samples.tsv, their samples'
        SAMs in protal/); every sample also at the levels given (<folder>_c<level>/, <sample>__c<level>, its SAM's
        extension kept), as collect_training_data.unit_map_rows names them."""
        def files(folder, sample, ending):
            for level in (None, *levels):
                at, name = (folder, sample) if level is None else (f"{folder}_c{level}", f"{sample}__c{level}")
                for sub, end in (("alignments", ending), ("profiles", ".profile.truth_annotated")):
                    os.makedirs(os.path.join(at, sub), exist_ok=True)
                    open(os.path.join(at, sub, name + end), "w").close()
        for name, (pe, se) in simulated.items():
            base, sample = os.path.join(training, "points", name), name + "_s_1"
            os.makedirs(os.path.join(base, "sim"))
            with open(os.path.join(base, "sim", "protal.meta"), "w") as fh:
                fh.write(f"#OUTPUT_DIR\t{os.path.join(base, 'protal')}\n#SAMPLEID\tSAM\tPROFILE_TRUTH\n"
                         f"{sample}\t{sample}.sam.zst\t/truth\n")
            if pe:
                files(os.path.join(base, "protal"), sample, ".sam.zst")
            if se:
                files(os.path.join(base, "protal_se"), sample + "_se", ".sam.gz")
        for name in drawn:
            base, sample = os.path.join(training, "points", name), name + "_s_1"
            os.makedirs(os.path.join(base, "sim"))
            with open(os.path.join(base, "sim", "samples.tsv"), "w") as fh:
                fh.write(f"sample\treads\ttruth\tcommunity\n{sample}\t{sample}.fq.zst\t/truth_{name}\trl150_p10_s_1\n")
            files(os.path.join(base, "protal"), sample, ".sam.gz")

    def test_ultima_units_and_the_choice(self):
        """The design's Ultima units (se_ultima_r<reads>: drawn like long reads, sim/samples.tsv, SAMs in protal/) are
        single-end points too (since 2026-10-15), the scenarios' not; the points chosen take each kind's share, so
        that the default 2 re-profile an Illumina point and an Ultima one."""
        import check_model_parity as cmp
        with tempfile.TemporaryDirectory() as training:
            self.parity_layout(training, {"rl100_p10": (True, False), "rl150_p10": (True, True),
                                          "rl150_p20": (True, True), "rl250_p10": (True, True)},
                               ["se_ultima_r10", "se_ultima_r20", "pb_b10", "sc_soil_se_ultima_r10"])
            se = cmp.points_with_samples(training, "se")
            self.assertEqual([os.path.basename(p) for p, _ in se],
                             ["rl150_p10", "rl150_p20", "rl250_p10", "se_ultima_r10", "se_ultima_r20"])
            ultima = dict(se)[os.path.join(training, "points", "se_ultima_r10")]
            self.assertEqual(ultima, [("se_ultima_r10_s_1", os.path.join(training, "points", "se_ultima_r10", "protal",
                                                                         "alignments", "se_ultima_r10_s_1.sam.gz"),
                                       "/truth_se_ultima_r10",
                                       os.path.join(training, "points", "se_ultima_r10", "protal", "profiles",
                                                    "se_ultima_r10_s_1.profile.truth_annotated"))])

            def chosen(points, n):
                return [os.path.basename(p) for p, _ in cmp.choose(points, n)]
            self.assertEqual(chosen(se, 2), ["rl150_p10", "se_ultima_r10"])
            self.assertEqual(chosen(se, 3), ["rl150_p10", "rl250_p10", "se_ultima_r10"])
            self.assertEqual(chosen(se, 4), ["rl150_p10", "rl250_p10", "se_ultima_r10", "se_ultima_r20"])
            self.assertEqual(chosen(se, 9), [os.path.basename(p) for p, _ in se])
            # One kind of points: the first and the last, as before.
            pe = cmp.points_with_samples(training, "pe")
            self.assertEqual(chosen(pe, 2), ["rl100_p10", "rl250_p10"])
            self.assertEqual([os.path.basename(p) for p, _ in cmp.points_with_samples(training, "pb")], ["pb_b10"])
            self.assertEqual(cmp.points_with_samples(training, "ont"), [])

    def test_samples_at_a_level(self):
        """--level C: the samples as the collection profiled them at another sparse-first gate (--model_levels:
        <sample>__c<C> in <profile folder>_c<C>/, the SAM's extension that of the sample's own), for a model of that
        gate; a point without them is left out; the outputs in TRAINING/parity[_<type>]_c<C>."""
        import contextlib
        import io
        import shutil
        import check_model_parity as cmp
        with tempfile.TemporaryDirectory() as training:
            self.parity_layout(training, {"rl150_p10": (True, True), "rl250_p10": (True, True)}, ["se_ultima_r10"],
                               levels=("1", "8s"))
            shutil.rmtree(os.path.join(training, "points", "rl250_p10", "protal_c1"))
            points = os.path.join(training, "points")

            def at(t, level):
                return [s for _, samples in cmp.points_with_samples(training, t, level) for s in samples]
            self.assertEqual(at("pe", "1"), [("rl150_p10_s_1__c1",
                                              os.path.join(points, "rl150_p10", "protal_c1", "alignments",
                                                           "rl150_p10_s_1__c1.sam.zst"), "/truth",
                                              os.path.join(points, "rl150_p10", "protal_c1", "profiles",
                                                           "rl150_p10_s_1__c1.profile.truth_annotated"))])
            self.assertEqual([(s[0], os.path.relpath(s[1], points)) for s in at("se", "1")],
                             [("rl150_p10_s_1_se__c1", os.path.join("rl150_p10", "protal_se_c1", "alignments",
                                                                     "rl150_p10_s_1_se__c1.sam.gz")),
                              ("rl250_p10_s_1_se__c1", os.path.join("rl250_p10", "protal_se_c1", "alignments",
                                                                     "rl250_p10_s_1_se__c1.sam.gz")),
                              ("se_ultima_r10_s_1__c1", os.path.join("se_ultima_r10", "protal_c1", "alignments",
                                                                      "se_ultima_r10_s_1__c1.sam.gz"))])
            self.assertEqual([s[0] for s in at("pe", "8s")], ["rl150_p10_s_1__c8s", "rl250_p10_s_1__c8s"])
            self.assertEqual(at("pe", "4"), [])
            self.assertEqual([s[0] for s in at("pe", None)], ["rl150_p10_s_1", "rl250_p10_s_1"])
        self.assertEqual(cmp.default_out("T", "pe"), os.path.join("T", "parity"))
        self.assertEqual(cmp.default_out("T", "se", "1"), os.path.join("T", "parity_se_c1"))
        self.assertEqual(cmp.level_text("1f"), "sparse-first gate 1, --flex_indels")
        self.assertEqual(cmp.parse_args(["--db", "d", "--model", "m", "--training", "t", "--level", "8"]).level, "8")
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            cmp.parse_args(["--db", "d", "--model", "m", "--training", "t", "--level", "c8"])


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


@prerequisites.requires(HAVE_SKLEARN, NEEDS_SKLEARN)
class CrossedFoldsTest(unittest.TestCase):
    """The trainer's folds that hold out a row's sample and species together (--cv-scheme samples+species, CrossedFolds;
    the audit's E4): a fold's test rows share neither their sample nor their species with its training rows, which are
    all the rows of neither group; rotations score rotations/k of the rows, each once; the groups are drawn by the seed;
    an unseen table's species take the group of their name's taxon, or a seeded hash."""

    @classmethod
    def setUpClass(cls):
        import machine_learning_cmdline
        cls.trainer = machine_learning_cmdline
        rng = np.random.default_rng(4)
        rows = []
        for sample in range(23):
            for taxon in rng.choice(70, size=int(rng.integers(5, 30)), replace=False):
                rows.append({"meta_sample": f"s{sample}", "taxon": int(taxon), "taxon_name": f"species {taxon}",
                             "truth": int(taxon < 20)})
        cls.df = pd.DataFrame(rows)

    def opts(self, **kw):
        return types.SimpleNamespace(**{"folds": 5, "seed": 3, "cv_rotations": 2, **kw})

    def groups(self, crossed, df=None):
        df = self.df if df is None else df
        return (df["meta_sample"].map(crossed.sample_group).to_numpy(),
                df["taxon"].astype(str).map(crossed.species_group).to_numpy())

    def test_folds(self):
        crossed = self.trainer.crossed_folds(self.df, self.opts())
        self.assertEqual((crossed.k, crossed.rotations, len(crossed.splits)), (5, 2, 10))
        y = self.df["truth"].to_numpy()
        self.assertEqual([(list(a), list(b)) for a, b in self.trainer.folds(self.df, y, self.trainer.CROSSED, self.opts())],
                         [(list(a), list(b)) for a, b in crossed.splits])
        sg, tg = self.groups(crossed)
        samples, species = self.df["meta_sample"].to_numpy(), self.df["taxon"].to_numpy()
        times = np.zeros(len(self.df), dtype=int)
        for train, test in crossed.splits:
            i, j = sg[test[0]], tg[test[0]]
            np.testing.assert_array_equal(test, np.flatnonzero((sg == i) & (tg == j)))
            np.testing.assert_array_equal(train, np.flatnonzero((sg != i) & (tg != j)))
            self.assertFalse(set(samples[train]) & set(samples[test]))
            self.assertFalse(set(species[train]) & set(species[test]))
            times[test] += 1
        self.assertEqual(times.max(), 1)
        self.assertEqual(crossed.scored, int(times.sum()))
        self.assertAlmostEqual(times.mean(), 2 / 5, delta=0.08)
        np.testing.assert_array_equal(crossed.scored_cells(sg, tg), times == 1)
        self.assertIn(f"2 rotations of 5 folds (10 fits), each fitted on the rows of neither its sample group nor its "
                      f"species group; {crossed.scored} of {len(self.df)} rows scored", crossed.describe())
        # The groups hold about as many rows each.
        for g in (sg, tg):
            sizes = np.bincount(g, minlength=5)
            self.assertLess(sizes.max() - sizes.min(), 0.2 * sizes.mean(), sizes)
        # Every rotation (at most k): every row once.
        every = self.trainer.crossed_folds(self.df, self.opts(cv_rotations=9))
        self.assertEqual(every.rotations, 5)
        self.assertEqual(every.scored, len(self.df))
        # The same seed the same groups, another seed others.
        self.assertEqual(self.trainer.crossed_folds(self.df, self.opts()).sample_group, crossed.sample_group)
        self.assertNotEqual(self.trainer.crossed_folds(self.df, self.opts(seed=4)).sample_group, crossed.sample_group)
        # Without samples, or with one: none (the trainer falls back to species held out).
        self.assertIsNone(self.trainer.crossed_folds(self.df.drop(columns="meta_sample"), self.opts()))
        self.assertIsNone(self.trainer.crossed_folds(self.df[self.df["meta_sample"] == "s0"], self.opts()))

    def test_unseen_species_in_scored_cells(self):
        crossed = self.trainer.crossed_folds(self.df, self.opts())
        hashed = self.trainer.hashed_group
        unseen = pd.DataFrame({"meta_sample": ["s0", "s1", "s99", "s2"],
                               "taxon_name": ["species 3", "species 999", "species 4", "species 5"]})
        expected = [crossed.scored_cells(crossed.sample_group["s0"], crossed.species_group["3"]),
                    crossed.scored_cells(crossed.sample_group["s1"], hashed("species 999", 5, 3)),
                    crossed.scored_cells(hashed("s99", 5, 3), crossed.species_group["4"]),
                    crossed.scored_cells(crossed.sample_group["s2"], crossed.species_group["5"])]
        np.testing.assert_array_equal(crossed.unseen_scored(unseen), expected)
        self.assertEqual(len(crossed.unseen_scored(unseen.iloc[:0])), 0)
        self.assertEqual(hashed("species 999", 5, 3), hashed("species 999", 5, 3))
        # Over many unseen species, about the share of the rows scored.
        many = pd.DataFrame({"meta_sample": [f"s{i % 23}" for i in range(2000)],
                             "taxon_name": [f"new species {i}" for i in range(2000)]})
        self.assertAlmostEqual(crossed.unseen_scored(many).mean(), 2 / 5, delta=0.05)


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
        # choose_global_knob: the median of the best thresholds held out (--cv-scheme: samples and species by default)
        # over bootstrap resamples of the training samples, kept if it beats 0.5 in KNOB_SUPPORT of them, the rows
        # weighted as the models weigh them (the scenarios' by --scenario-weight). The test set does not choose (until
        # 2026-10-15 a knob that lost there was dropped: the audit's E2).
        crossed = self.trainer.CROSSED
        rng = np.random.default_rng(2)
        n = 4000
        y = (rng.random(n) < 0.3).astype(int)
        scores = np.clip(np.where(y == 1, rng.normal(0.85, 0.08, n), rng.normal(0.45, 0.12, n)), 0, 1)
        frame = pd.DataFrame({"truth": y, "meta_scenario": "", "meta_sample": [f"s{i % 40}" for i in range(n)]})
        opts = types.SimpleNamespace(knob=0.5, seed=1)
        report = self.trainer.Report()
        knob = self.trainer.choose_global_knob(report, frame, y, {crossed: scores}, opts)
        self.assertIsNotNone(knob)
        self.assertGreater(knob, 0.55)  # absent taxa score up to ~0.7: a knob above 0.5 calls fewer of them
        chosen = report.data["global_knob"]
        self.assertEqual((chosen["knob"], chosen["median"], chosen["scheme"]), (knob, knob, crossed))
        self.assertGreaterEqual(chosen["support"], self.trainer.KNOB_SUPPORT)
        self.assertGreater(chosen["gain"], 0)
        self.assertLessEqual(chosen["range"][0], knob)
        self.assertGreaterEqual(chosen["range"][1], knob)
        self.assertIsNone(chosen["species_only"])
        self.assertEqual(self.trainer.choose_global_knob(self.trainer.Report(), frame, y, {crossed: scores}, opts), knob)
        # A test set it loses on (its present taxa scored between 0.5 and the knob): the knob stands, the test set's F1
        # at it and at 0.5 reported.
        test_y = np.ones(50, dtype=int)
        report = self.trainer.Report()
        self.assertEqual(self.trainer.choose_global_knob(report, frame, y, {crossed: scores}, opts,
                                                         (test_y, np.full(50, 0.52), np.ones(50))), knob)
        tested = report.data["global_knob"]["test"]
        self.assertEqual((tested["knob"], tested["F1_at_knob"], tested["F1_at_default"], tested["F1_at_median"], tested["rows"]),
                         (knob, 0.0, 1.0, 0.0, 50))
        self.assertIn(f"on the test set (50 rows, weighted; for information, it does not choose): F1 0.0000 at the knob "
                      f"{knob:g}, 1.0000 at 0.5", "\n".join(report.lines))
        # Only the decisions' scheme chooses; the knob of species held out alone is reported beside it.
        report = self.trainer.Report()
        even = np.where(y == 1, 0.9, 0.1)  # every threshold from 0.1 to 0.9 as good: no knob beats 0.5
        self.assertEqual(self.trainer.choose_global_knob(report, frame, y, {crossed: scores, "species": even}, opts), knob)
        self.assertIsNone(report.data["global_knob"]["species_only"]["knob"])
        self.assertIn("with species held out alone (--cv-scheme species", "\n".join(report.lines))
        species = types.SimpleNamespace(knob=0.5, seed=1, cv_scheme="species")
        self.assertIsNone(self.trainer.choose_global_knob(self.trainer.Report(), frame, y, {crossed: scores, "species": even},
                                                          species))
        self.assertEqual(self.trainer.choose_global_knob(self.trainer.Report(), frame, y, {"species": scores}, species), knob)
        self.assertIsNone(self.trainer.choose_global_knob(self.trainer.Report(), frame, y, {crossed: even}, opts))
        # Without held-out scores (--evaluation none): none.
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

        # The features in steps of 0.1: a column of 255 values or fewer is binned at the midpoints between them. With
        # more, scikit-learn 1.8 and later bins a weighted fit (the scenario rows weigh 0.25, --scenario-weight) by
        # computing its 254 weighted percentiles one at a time, which made this class take 5 minutes instead of about one
        # (docs/claude/2026-10-12-test-suite-speed). The trainer computes them at once since (hgb_binning.py); the steps
        # keep the class fast where hgb_binning steps aside (another scikit-learn).
        def table(path, samples, seed, scenario_of):
            rng = np.random.default_rng(seed)
            rows = []
            for sample in range(samples):
                scenario = scenario_of(sample)
                for taxon in range(20):
                    present = taxon < 6
                    row = {c: round(rng.normal(0, 1), 1) for c in columns}
                    # the signal in a normalized feature and in the depth feature's group, so that the sets differ
                    row["hit_gene_fraction"] = round(1.2 * present + rng.normal(0, 0.6), 1)
                    row["sample_log_fragments"] = round(0.8 * present + rng.normal(0, 0.6), 1)
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
                                     ("gut", "hold-in, samples and species held out"),
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

    def test_the_forest_by_default_and_boosting_on_request(self):
        # The default model is a random forest again since 2026-10-14 (--model forest), averaged; --model gbm gradient-boosted
        # trees, written as a chain; each scores in PMML as scikit-learn does (the trainer's own check), with its importances.
        prefix, metrics, stdout = self.train("forest_default", "normalized")
        self.assertEqual((metrics["model"]["kind"], metrics["model"]["trees"]), ("forest", 16))
        self.assertIsInstance(load_model(prefix + ".xml"), PmmlForest)
        prefix, metrics, stdout = self.train("boosted", "normalized", None, "--model", "gbm")
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
        side, side_metrics, side_stdout = self.train("boosted_side_by_side", "normalized", None, "--fold-jobs", "2", "--model", "gbm")
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
        self.assertRegex(auto["why"], rf"^F1 [0-9.]+ with samples and species held out \(rows weighted as in the fits\), "
                                      rf"the other {len(others)} sets [0-9.]+ on average "
                                      r"\(the best of them \S+ [0-9.]+\); (the highest|the default set, as no other is "
                                      rf"{re.escape(str(self.features.AUTO_MIN_GAIN))} better)[^;]*; on samples never trained "
                                      r"on \(not used to choose\): test set [0-9.]+ \(the others [0-9.]+\), gut hold-out "
                                      r"[0-9.]+ \(the others [0-9.]+\), host hold-out")
        # The model takes the chosen set's features, and the evaluation scored that set with the same folds.
        self.assertEqual(load_model(prefix + ".xml").features, self.features.feature_set_columns(auto["chosen"]))
        self.assertIn("## Feature set chosen (--features auto, samples and species held out)", stdout)
        self.assertIn(f"feature set (--features auto): {auto['chosen']}: {auto['why']}", stdout)
        unweighted = {r["features"]: r["F1 unweighted"] for r in auto["candidates"]}
        self.assertAlmostEqual(metrics["evaluation"]["samples+species"]["F1"], unweighted[auto["chosen"]])
        self.assertEqual((auto["scheme"], auto["weighted"]), ("samples+species", True))
        # The F1 that chooses weighs the rows as the fits do: the gut scenario's rows 0.25 (--scenario-weight).
        predictions = pd.read_csv(prefix + ".predictions.tsv.gz", sep="\t")
        scored = predictions["p_samples+species"].notna().to_numpy()
        call = predictions["p_samples+species"].to_numpy()[scored] >= 0.5
        truth = predictions["truth"].to_numpy()[scored]
        weight = np.where(predictions["meta_scenario"].fillna("").to_numpy()[scored] != "", 0.25, 1.0)
        tp, fp, fn = (weight * (call & (truth == 1))).sum(), (weight * (call & (truth == 0))).sum(), \
            (weight * (~call & (truth == 1))).sum()
        self.assertAlmostEqual(f1[auto["chosen"]], 2 * tp / (2 * tp + fp + fn))
        self.assertNotAlmostEqual(f1[auto["chosen"]], unweighted[auto["chosen"]])


@prerequisites.requires(HAVE_SKLEARN, NEEDS_SKLEARN)
class TrainerHeadlineTest(unittest.TestCase):
    """The trainer's headline (the audit's E3): the marker-truth F1, first in report.txt and metrics.json, held out (with
    samples and species, and with species alone) and on the test set, the present species without a row counted as
    missed where their cells were scored, and said so where an unseen table is missing; the features' ranges in every
    model and their versions with --feature-versions; --cv-scheme species as before."""

    @classmethod
    def setUpClass(cls):
        import machine_learning_cmdline
        cls.trainer = machine_learning_cmdline
        cls.tmp = tempfile.TemporaryDirectory()
        cls.trained = {}

        def table(path, samples, seed):
            rng = np.random.default_rng(seed)
            rows = []
            for sample in range(samples):
                for i in range(25):
                    present, taxon = i < 10, (i + 2 * sample) % 70
                    rows.append({"truth": int(present), "taxon": 200 + taxon, "taxon_name": f"s__T {taxon}",
                                 "meta_sample": f"s{seed}_{sample}", "meta_scenario": "",
                                 "meta_marker_reads": int(rng.integers(1, 9)) if present else "",
                                 "fragments": float(rng.integers(1, 200)), "x": round(1.0 * present + rng.normal(0, 0.7), 1),
                                 "y": round(rng.normal(0, 1), 1)})
            frame = pd.DataFrame(rows)  # the sample's depth among the features: one knob for every sample (--depth-knobs)
            frame["sample_log_fragments"] = np.round(np.log10(frame.groupby("meta_sample")["fragments"].transform("sum")), 2)
            frame.to_csv(path, sep="\t", index=False)
        cls.training = os.path.join(cls.tmp.name, "training.tsv")
        cls.test = os.path.join(cls.tmp.name, "test.tsv")
        table(cls.training, 30, 1)
        table(cls.test, 8, 2)
        # The training table's unseen species (collect_training_data.py's <table>.unseen.tsv): of the table's species and
        # one it has no row of; the test table has none (missing).
        cls.unseen = pd.DataFrame({"meta_design": "d", "meta_sample": [f"s1_{i}" for i in range(12)], "meta_scenario": "",
                                   "taxon_name": [f"s__T {i * 5}" for i in range(11)] + ["s__T 999"], "marker_reads": 2})
        cls.unseen.to_csv(os.path.join(cls.tmp.name, "training.unseen.tsv"), sep="\t", index=False)
        cls.versions = os.path.join(cls.tmp.name, "versions.tsv")
        with open(cls.versions, "w") as fh:
            fh.write("feature\tversion\tchanged\tnote\nx\t2\t2026-10-15\tx's code\nfragments\t1\t2026-08-10\t\n"
                     "other\t5\t\t\n")

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def train(self, name, *extra):
        """The trainer on the training table with the options `extra`, once per name -> (prefix, metrics, console)."""
        if name not in self.trained:
            prefix = os.path.join(self.tmp.name, name)
            result = subprocess.run([sys.executable, os.path.join(HERE, "machine_learning_cmdline.py"), "--truth-file",
                                     self.training, "--output-prefix", prefix, "--features", "all", "--ntree", "16",
                                     "--evaluation", "basic", "--threads", "1", "--depth-knobs", "--test-file", self.test,
                                     *extra], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout[-3000:] + result.stderr[-3000:])
            with open(prefix + ".metrics.json") as fh:
                self.trained[name] = (prefix, json.load(fh), result.stdout)
        return self.trained[name]

    def test_headline(self):
        prefix, metrics, stdout = self.train("crossed")
        self.assertEqual(next(iter(metrics)), "headline")
        headline = metrics["headline"]
        self.assertEqual((headline["measure"], headline["scheme"]), ("marker truth", "samples+species"))
        self.assertEqual(set(headline) - {"knob", "measure", "scheme", "notes"}, {"held_out", "held_out_species", "test"})
        # Every present row has a marker read here: the marker truth is the rows' calls plus the unseen species as FN.
        for key, h in ((k, headline[k]) for k in ("held_out", "held_out_species", "test")):
            self.assertEqual((h["marker_truth_tp"], h["marker_truth_fp"]), (h["rows_tp"], h["rows_fp"]), key)
            self.assertEqual(h["marker_truth_fn"], h["rows_fn"] + h["unseen"], key)
            tp, fp, fn = h["marker_truth_tp"], h["marker_truth_fp"], h["marker_truth_fn"]
            self.assertAlmostEqual(h["marker_truth_f1"], 2 * tp / (2 * tp + fp + fn), msg=key)
        # Held out with samples and species: the unseen species of the scored cells; with species alone, all of them;
        # the test table has no unseen table, which the headline says.
        frame = self.trainer.load_table(self.training)
        crossed = self.trainer.crossed_folds(frame, types.SimpleNamespace(folds=5, seed=1, cv_rotations=2))
        self.assertEqual(headline["held_out"]["unseen"], int(crossed.unseen_scored(self.unseen).sum()))
        self.assertEqual(headline["held_out"]["rows_scored"], crossed.scored)
        self.assertEqual(headline["held_out_species"]["unseen"], len(self.unseen))
        self.assertEqual(headline["held_out_species"]["rows_scored"], len(frame))
        self.assertEqual(headline["test"]["unseen"], 0)
        self.assertEqual(len(headline["notes"]), 1)
        self.assertIn("independent test set: no ", headline["notes"][0])
        self.assertIn("test.unseen.tsv", headline["notes"][0])
        # First in report.txt, after the run's description; and first in the summary.
        with open(prefix + ".report.txt") as fh:
            report = fh.read()
        self.assertLess(report.index("cross-validation folds fitted"), report.index("## Headline"))
        self.assertEqual(report.index("\n## ") + 1, report.index("## Headline"), "the first section")
        self.assertRegex(report, r"\nsamples and species held out \(\d+ of 750 rows scored\): marker-truth F1 [0-9.]+ "
                                 r"\(TP \d+, FP \d+, FN \d+, \d+ of them unseen\); row-based F1 [0-9.]+ \(TP \d+, FP \d+, "
                                 r"FN \d+\)\n")
        summary = stdout[stdout.index("## Summary"):]
        self.assertIn("\nspecies held out alone, as the decisions before 2026-10-15 (750 rows): marker-truth F1", summary)
        self.assertIn("\nindependent test set (200 rows): marker-truth F1", summary)
        self.assertIn("Note: independent test set: no ", summary)
        # The evaluation: both schemes, the crossed one on its share of the rows; the knob chosen on it.
        self.assertEqual(metrics["evaluation"]["samples+species"]["taxa"], crossed.scored)
        self.assertEqual(metrics["evaluation"]["samples+species"]["fits"], 10)
        self.assertEqual(metrics["evaluation"]["species"]["taxa"], len(frame))
        self.assertEqual(metrics["cv_scheme"]["fits"], 10)
        self.assertEqual(metrics["global_knob"]["scheme"], "samples+species")
        self.assertEqual(metrics["global_knob"]["scored"], crossed.scored)
        self.assertIsNotNone(metrics["global_knob"]["species_only"])
        self.assertIn(("samples and species held out", "(all)"),
                      {(r["set"], r["scenario"]) for r in metrics["marker_truth"]["rows"]})
        # Every model's ranges over the training rows; no versions without --feature-versions.
        ranges = read_feature_ranges(prefix + ".xml")
        self.assertEqual(list(ranges), ["fragments", "x", "y", "sample_log_fragments"])
        self.assertEqual(ranges["x"][0], frame["x"].min())
        self.assertEqual(ranges["x"][3], frame["x"].max())
        self.assertEqual(read_feature_versions(prefix + ".xml"), {})

    def test_species_scheme_and_feature_versions(self):
        prefix, metrics, stdout = self.train("species", "--cv-scheme", "species", "--feature-versions", self.versions)
        headline = metrics["headline"]
        self.assertEqual(headline["scheme"], "species")
        self.assertEqual(set(headline) - {"knob", "measure", "scheme", "notes"}, {"held_out", "test"})
        self.assertEqual(headline["held_out"]["estimate"], "species held out")
        self.assertEqual(headline["held_out"]["unseen"], len(self.unseen))
        self.assertNotIn("samples+species", metrics["evaluation"])
        self.assertEqual(metrics["global_knob"]["scheme"], "species")
        self.assertIsNone(metrics["global_knob"]["species_only"])
        # The same species estimate as the crossed run's.
        _, crossed, _ = self.train("crossed")
        self.assertEqual(metrics["evaluation"]["species"], crossed["evaluation"]["species"])
        self.assertEqual(headline["held_out"]["rows_scored"], crossed["headline"]["held_out_species"]["rows_scored"])
        # --feature-versions: the model's features the file has, the others named.
        self.assertEqual(read_feature_versions(prefix + ".xml"), {"fragments": "1", "x": "2"})
        self.assertIn("feature versions (", stdout)
        self.assertIn("2 of the 4 features in the model's header (protal_feature_versions); not in the file, so left "
                      "out: y, sample_log_fragments", stdout)
        self.assertEqual(metrics["model"]["feature_versions"], {"fragments": "1", "x": "2"})

    def test_a_bad_versions_file_stops_at_once(self):
        bad = os.path.join(self.tmp.name, "bad_versions.tsv")
        with open(bad, "w") as fh:
            fh.write("feature\tversion\nx\t1,2\n")
        result = subprocess.run([sys.executable, os.path.join(HERE, "machine_learning_cmdline.py"), "--truth-file",
                                 self.training, "--output-prefix", os.path.join(self.tmp.name, "bad"), "--feature-versions",
                                 bad], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--feature-versions:", result.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.tmp.name, "bad.metrics.json")))


@prerequisites.requires(HAVE_SKLEARN, NEEDS_SKLEARN)
class HgbBinningTest(unittest.TestCase):
    """hgb_binning.py, which the trainer installs: the bin thresholds of a weighted gradient-boosting fit in one pass per
    feature (scikit-learn 1.8 and later computes each of a feature's 254 weighted percentiles on its own) are
    scikit-learn's, bit for bit, and so is the model; a scikit-learn that bins without the weights (before 1.8) is left
    as it is."""

    def setUp(self):
        import hgb_binning
        from sklearn.ensemble._hist_gradient_boosting import binning
        self.hgb, self.binning = hgb_binning, binning
        # scikit-learn's own, also when a trainer test of this process installed hgb_binning before
        self.original = hgb_binning._original or binning._find_binning_thresholds
        self.weighted = "sample_weight" in inspect.signature(self.original).parameters
        self.tearDown()

    def tearDown(self):
        self.binning._find_binning_thresholds = self.original
        self.hgb._installed = self.hgb._original = None

    def test_install(self):
        said = self.hgb.install()
        if not self.weighted:
            self.assertIn("bins without the sample weights", said)
            self.assertIs(self.binning._find_binning_thresholds, self.original)
            return
        self.assertIn("one pass per feature", said)
        self.assertIs(self.binning._find_binning_thresholds, self.hgb._find_binning_thresholds)
        self.assertIs(self.hgb._original, self.original)
        self.assertEqual(self.hgb.install(), said)  # once per process

    def test_the_thresholds_are_scikit_learns(self):
        if not self.weighted:
            self.skipTest("this scikit-learn bins without the sample weights")
        rng = np.random.default_rng(3)
        columns = self.hgb.check_columns(n=3000, seed=11) + \
            [("half the rows at 0.25", rng.normal(size=5000), np.where(rng.random(5000) < 0.5, 0.25, 1.0))]
        for name, column, weights in columns:
            for max_bins in (255, 64):
                with self.subTest(name, max_bins=max_bins):
                    theirs = self.original(column.copy(), max_bins, sample_weight=None if weights is None else weights.copy())
                    ours = self.hgb._find_binning_thresholds(column.copy(), max_bins,
                                                             sample_weight=None if weights is None else weights.copy())
                    self.assertEqual(ours.dtype, theirs.dtype)
                    np.testing.assert_array_equal(ours, theirs)

    def test_a_weighted_fit_is_the_same_model(self):
        if not self.weighted:
            self.skipTest("this scikit-learn bins without the sample weights")
        from sklearn.ensemble import HistGradientBoostingClassifier
        rng = np.random.default_rng(5)
        n = 2000
        X = rng.normal(size=(n, 6))
        X[:, 5] = np.round(X[:, 5], 1)  # a feature of few values
        y = (X[:, 0] + 0.5 * X[:, 1] + rng.normal(0, 0.7, n) > 0).astype(int)
        X[rng.random(n) < 0.05, 2] = np.nan
        weights = np.where(rng.random(n) < 0.3, 0.25, 1.0)

        def fit():
            model = HistGradientBoostingClassifier(max_iter=20, max_leaf_nodes=31, early_stopping=False, random_state=1)
            model.fit(X, y, sample_weight=weights)
            return model._bin_mapper.bin_thresholds_, model.predict_proba(X)[:, 1]
        thresholds, p = fit()
        self.assertIn("one pass per feature", self.hgb.install())
        installed_thresholds, installed_p = fit()
        for theirs, ours in zip(thresholds, installed_thresholds, strict=True):
            np.testing.assert_array_equal(ours, theirs)
        np.testing.assert_array_equal(installed_p, p)


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


class ExtraTestsUnseenTest(unittest.TestCase):
    """--extra-test-file (a test set profiled otherwise, e.g. at another sparse-first gate): the present species of its
    samples without a row in its table (its <table>.unseen.tsv) count as missed, so that a gate giving up all of a
    species' reads does not look better for leaving it out (2026-10-13)."""

    def test_unseen_species_count_as_missed(self):
        import unittest.mock
        import pandas as pd
        import machine_learning_cmdline as trainer

        class Forest:
            def predict_proba(self, x):
                return np.column_stack([1 - x[:, 0], x[:, 0]])
        table = pd.DataFrame({"meta_sample": ["a", "a", "b"], "truth": [1, 0, 1], "f": [0.9, 0.8, 0.1]})
        opts = types.SimpleNamespace(knob=0.5)
        with unittest.mock.patch.object(trainer, "call_scores", lambda frame, p: p):
            plain = trainer.extra_test_scores(Forest(), ["f"], table, opts, None)
            seen = trainer.extra_test_scores(Forest(), ["f"], table, opts, None, unseen=2)
        self.assertEqual({k: plain[k] for k in ("samples", "TP", "FP", "FN", "present", "unseen")},
                         {"samples": 2, "TP": 1, "FP": 1, "FN": 1, "present": 2, "unseen": 0})
        self.assertEqual((seen["TP"], seen["FN"], seen["present"], seen["unseen"]), (1, 3, 4, 2))
        self.assertAlmostEqual(seen["F1"], 2 / (2 + 1 + 3))
        self.assertAlmostEqual(seen["sensitivity"], 1 / 4)
        unseen = pd.DataFrame({"meta_design": ["d"] * 3, "meta_sample": ["a", "b", "c"], "meta_scenario": [""] * 3,
                               "taxon_name": ["x", "y", "z"], "marker_reads": [3, 1, 2]})
        self.assertEqual(trainer.unseen_in(unseen, {"a", "b"}), 2)
        self.assertEqual(trainer.unseen_in(None, {"a"}), 0)

    def test_over_the_samples_every_table_has(self):
        """--extra-test-samples (since 2026-10-15): each extra table scored over those samples only, beside the test
        file's rows of the same samples; a table of more samples (a model's gate, with the streamed ones) compares with
        the others, and one with none of them has no scores."""
        import unittest.mock
        import pandas as pd
        import machine_learning_cmdline as trainer

        class Forest:
            def predict_proba(self, x):
                return np.column_stack([1 - x[:, 0], x[:, 0]])

        class Report:
            def __init__(self):
                self.data = {}

            def section(self, *_):
                pass

            def add(self, *_):
                pass

            def table(self, *_):
                pass
        tables = {"c1": pd.DataFrame({"meta_sample": ["a", "b", "s"], "truth": [1, 1, 1], "f": [0.9, 0.1, 0.9]}),
                  "c8": pd.DataFrame({"meta_sample": ["s"], "truth": [1], "f": [0.9]})}
        test = pd.DataFrame({"meta_sample": ["a", "b", "s"], "truth": [1, 1, 1], "f": [0.9, 0.9, 0.9]})
        with tempfile.TemporaryDirectory() as root:
            keep = os.path.join(root, "compare_samples_pe.txt")
            with open(keep, "w") as fh:
                fh.write("a\nb\n")
            opts = types.SimpleNamespace(knob=0.5, taxonomy=None, test_unseen_file="none", test_file=None,
                                         extra_test_file=[(label, label) for label in tables], extra_test_samples=keep)
            report = Report()
            with unittest.mock.patch.object(trainer, "call_scores", lambda frame, p: p), \
                    unittest.mock.patch.object(trainer, "load_table", lambda path, taxonomy: tables[path]):
                trainer.study_extra_tests(report, Forest(), ["f"], opts, None, test)
        c1, c8 = report.data["extra_tests"]["c1"], report.data["extra_tests"]["c8"]
        self.assertEqual((c1["scores"]["samples"], c1["scores"]["TP"], c1["scores"]["FN"]), (2, 1, 1))
        self.assertEqual((c1["test_file_same_samples"]["samples"], c1["test_file_same_samples"]["TP"]), (2, 2))
        self.assertIsNone(c8["scores"])
        self.assertIsNone(c8["test_file_same_samples"])


if __name__ == "__main__":
    unittest.main()
