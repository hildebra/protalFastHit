#!/usr/bin/env python3
"""Checks for model_pmml.py: a forest written as PMML scores as scikit-learn scores it, bit for bit.

The PMML file is scored with model_pmml.PmmlForest, which compares doubles and averages the trees in
file order as protal's cPMML does; check_model_parity.py compares PmmlForest with protal itself.
Also the knobs by sample depth that random_forest_cmdline.py --depth-knobs chooses and writes.

  python3 -m unittest scripts/test_model_pmml.py
"""

import json
import os
import subprocess
import sys
import tempfile
import types
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from model_pmml import (PmmlForest, float32_split, format_depth_knob_curve, format_depth_knobs,  # noqa: E402
                        read_depth_knob_curve, read_depth_knobs, read_false_calls, write_forest)

try:
    import pandas as pd
    from sklearn.ensemble import RandomForestClassifier
except ImportError:  # the checks need scikit-learn, as training does
    RandomForestClassifier = None


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


@unittest.skipIf(RandomForestClassifier is None, "needs numpy, pandas and scikit-learn")
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


@unittest.skipIf(RandomForestClassifier is None, "needs scikit-learn")
class FeatureSetsTest(unittest.TestCase):
    """The trainer's default features are the normalised ones, the gene neighbours', the four relatives features by the
    references' distance, the sample's depth and the divergence features; a set's name joins its groups with "+";
    normalized+adjacency leaves the rest out, normalized the gene neighbours' too, and
    normalized+adjacency+relatives adds all the relatives'."""

    def test_default_set_has_the_gene_neighbour_features(self):
        import model_features as mf
        import random_forest_cmdline
        columns = (["truth", "taxon", "meta_sample"] + mf.NORMALIZED_FEATURES + mf.ADJACENCY_FEATURES +
                   mf.RELATIVE_FEATURES + mf.SAMPLE_FEATURES + mf.DIVERGENCE_FEATURES + ["genus_top_fragments", "other"])
        self.assertEqual(mf.DEFAULT_FEATURE_SET, "normalized+adjacency+distance+depth+divergence")
        self.assertEqual(mf.feature_columns(columns, mf.DEFAULT_FEATURE_SET),
                         mf.NORMALIZED_FEATURES + mf.ADJACENCY_FEATURES + mf.DISTANCE_FEATURES + mf.SAMPLE_FEATURES + mf.DIVERGENCE_FEATURES)
        self.assertEqual(mf.feature_columns(columns, "normalized+adjacency+distance"),
                         mf.NORMALIZED_FEATURES + mf.ADJACENCY_FEATURES + mf.DISTANCE_FEATURES)
        # Groups in any order, each once, normalized among them; the sample's depth is recognised.
        self.assertEqual(mf.feature_columns(columns, "depth+normalized"), mf.NORMALIZED_FEATURES + mf.SAMPLE_FEATURES)
        self.assertTrue(mf.has_sample_depth(mf.feature_columns(columns, mf.DEFAULT_FEATURE_SET)))
        self.assertFalse(mf.has_sample_depth(mf.feature_columns(columns, "normalized+adjacency+distance")))
        for bad in ("normalized+adjacency+depth+depth", "adjacency", "normalized+ani", ""):
            with self.assertRaises(ValueError):
                mf.feature_set_columns(bad)
        self.assertEqual(mf.feature_set_name("normalized+divergence"), "normalized+divergence")
        self.assertEqual(mf.feature_set_name("all"), "all")
        # A table of an older protal lacks the new columns: a clear error, and the older set still works.
        older = ["truth", "taxon"] + mf.NORMALIZED_FEATURES + mf.ADJACENCY_FEATURES + mf.RELATIVE_FEATURES
        with self.assertRaises(RuntimeError):
            mf.feature_columns(older, mf.DEFAULT_FEATURE_SET)
        self.assertEqual(mf.feature_columns(older, "normalized+adjacency"), mf.NORMALIZED_FEATURES + mf.ADJACENCY_FEATURES)
        relatives = mf.feature_columns(columns, "normalized+adjacency+relatives")
        self.assertEqual(relatives, mf.NORMALIZED_FEATURES + mf.ADJACENCY_FEATURES + mf.RELATIVE_FEATURES)
        self.assertEqual(mf.feature_columns(columns, "normalized"), mf.NORMALIZED_FEATURES)
        self.assertTrue(set(mf.DISTANCE_FEATURES) <= set(mf.RELATIVE_FEATURES))
        self.assertNotIn("genus_top_fragments", relatives)  # it counts reads
        with self.assertRaisesRegex(RuntimeError, "adjacent_support"):  # a table of an older protal
            mf.feature_columns([c for c in columns if c != "adjacent_support"], mf.DEFAULT_FEATURE_SET)
        with self.assertRaisesRegex(RuntimeError, "relative_spill"):  # a table of protal before the relatives features
            mf.feature_columns([c for c in columns if c != "relative_spill"], mf.DEFAULT_FEATURE_SET)
        with self.assertRaisesRegex(RuntimeError, "em_own_share"):  # a table of protal before the relatives features
            mf.feature_columns([c for c in columns if c != "em_own_share"], "normalized+adjacency+relatives")
        opts = random_forest_cmdline.parse_args(["--truth-file", "t.tsv", "--output-prefix", "p"])
        self.assertEqual(opts.features, mf.DEFAULT_FEATURE_SET)


@unittest.skipIf(RandomForestClassifier is None, "needs pandas")
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


def metrics_knobs(test, name):
    """The knob points of the model `name` trained in a TrainerDepthKnobsTest."""
    with open(os.path.join(test.tmp.name, name + ".metrics.json")) as fh:
        return json.load(fh)["depth_knobs"]["points"]


class TrainerDepthKnobsTest(unittest.TestCase):
    """random_forest_cmdline.py --depth-knobs on a table of shallow samples (hundreds of fragments, log10 ~2.7) and
    deep ones (tens of thousands, ~4.7), where a present taxon's evidence grows with depth; half the present taxa
    simulated from another genome than the representative."""

    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, HERE)
        import random_forest_cmdline
        cls.trainer = random_forest_cmdline
        cls.tmp = tempfile.TemporaryDirectory()
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
        prefix = os.path.join(self.tmp.name, name)
        result = subprocess.run([sys.executable, os.path.join(HERE, "random_forest_cmdline.py"), "--truth-file", self.table,
                                 "--output-prefix", prefix, "--features", "all", "--ntree", "16", "--evaluation", "basic",
                                 "--threads", "1", *extra], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:] + result.stderr[-3000:])
        with open(prefix + ".metrics.json") as fh:
            return prefix + ".xml", json.load(fh)

    def test_depths_and_curve_as_protal(self):
        # profiler::DepthKnobAt: log10 of the sample's fragments over all its taxa (at least 1), linear between the
        # curve's points, the ends' beyond them.
        frame = pd.DataFrame({"meta_sample": ["a", "a", "b", "c", "e"], "fragments": [60.0, 40.0, 1000.0, 3e7, 0.0]})
        np.testing.assert_allclose(self.trainer.sample_depths(frame), [2, 2, 3, np.log10(3e7), 0])
        curve = [(2.0, 0.2), (4.0, 0.8)]
        np.testing.assert_allclose(self.trainer.knob_at(curve, np.array([0.0, 2.0, 3.0, 3.5, 4.0, 7.5])),
                                   [0.2, 0.2, 0.5, 0.65, 0.8, 0.8])
        calls = self.trainer.depth_knob_calls(np.array([0.3, 0.3]), np.array([2.0, 4.0]), curve, 0.5)
        self.assertEqual(list(calls), [True, False])
        self.assertEqual(list(self.trainer.depth_knob_calls(np.array([0.3, 0.6]), np.array([2.0, 4.0]), [], 0.5)),
                         [False, True])

    def test_sparse_depths_join_their_neighbours(self):
        # A point needs DEPTH_KNOB_MIN_SAMPLES samples in its window: a bin with fewer joins the next deeper one, and
        # bins left at the deep end the point before, so that a few deep samples set no knob of their own.
        def points(sample_depths):
            samples = np.array([f"s{i}" for i, _ in enumerate(sample_depths) for _ in range(3)])
            depths = np.repeat(np.array(sample_depths, dtype=float), 3)
            windows = self.trainer.depth_knob_windows(depths, samples, np.ones(len(depths), dtype=bool))
            return [(x, len(set(samples[in_group])), len(set(samples[window]))) for x, in_group, window in windows]

        self.assertEqual(self.trainer.DEPTH_KNOB_MIN_SAMPLES, 6)
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
        plain, metrics = self.train("plain")
        self.assertEqual(read_depth_knob_curve(plain), [])
        self.assertNotIn("depth_knobs", metrics)
        # A point keeps the knob unless its best knob gains DEPTH_KNOB_MIN_GAIN on its window.
        self.assertEqual(self.trainer.DEPTH_KNOB_MIN_GAIN, 0.002)
        for point in metrics_knobs(self, "knobs"):
            if point["knob"] is not None:
                self.assertEqual(point["knob"], point["best knob"] if point["gain"] >= 0.002 else 0.5, point)

    def test_no_depth_knob_curve_with_the_samples_depth_as_a_feature(self):
        # With sample_log_fragments among the features the forest sees the depth itself: no curve is fitted even with
        # --depth-knobs, and the report says why.
        table = pd.read_csv(self.table, sep="\t")
        totals = table.groupby("meta_sample")["fragments"].transform("sum")
        table["sample_log_fragments"] = np.log10(np.maximum(totals, 1.0))
        path = os.path.join(self.tmp.name, "with_depth.tsv")
        table.to_csv(path, sep="\t", index=False)
        prefix = os.path.join(self.tmp.name, "with_depth")
        result = subprocess.run([sys.executable, os.path.join(HERE, "random_forest_cmdline.py"), "--truth-file", path,
                                 "--output-prefix", prefix, "--features", "all", "--ntree", "16", "--evaluation", "basic",
                                 "--threads", "1", "--depth-knobs"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:] + result.stderr[-3000:])
        self.assertEqual(read_depth_knob_curve(prefix + ".xml"), [])
        with open(prefix + ".metrics.json") as fh:
            metrics = json.load(fh)
        self.assertEqual(metrics["depth_knobs"]["curve"], [])
        self.assertIn("depth is a feature", metrics["depth_knobs"]["skipped"])
        self.assertIn("No knob curve: the sample's depth is a feature", result.stdout)

    def test_strains(self):
        _, metrics = self.train("strains")
        rows = {(r["simulated from"], r["fragments"]): r for r in metrics["strains"]["by_fragments"]}
        self.assertEqual(rows[("another genome", "all")]["present"], 40 * 4)
        self.assertEqual(rows[("the representative", "all")]["present"], 40 * 4)
        self.assertIn("strains missed", [r["taxa"] for r in metrics["strains"]["features"]])

    def test_conservation_features_by_class(self):
        _, metrics = self.train("classes")
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
        _, metrics = self.train("no_previous")
        self.assertNotIn("previous_procedure", metrics)
        _, metrics = self.train("no_previous_said", "--no-previous-procedure")
        self.assertNotIn("previous_procedure", metrics)
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


@unittest.skipIf(RandomForestClassifier is None, "needs numpy, pandas and scikit-learn")
class TrainerFalseCallsTest(unittest.TestCase):
    """random_forest_cmdline.py --fdr-calls and the singleton rule: the calibration, the prior adjusted to each sample
    and the calls at a target share of false calls, as protal makes them (context::FalseCallKnob)."""

    @classmethod
    def setUpClass(cls):
        import random_forest_cmdline
        cls.trainer = random_forest_cmdline
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
        result = subprocess.run([sys.executable, os.path.join(HERE, "random_forest_cmdline.py"), "--truth-file", self.table,
                                 "--output-prefix", prefix, "--features", "all", "--ntree", "16", "--evaluation", "basic",
                                 "--threads", "1", "--test-file", self.table, *extra], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:] + result.stderr[-3000:])
        with open(prefix + ".metrics.json") as fh:
            return prefix + ".xml", json.load(fh), result.stdout

    def test_prior_adjustment_as_protal(self):
        # Probabilities whose mean is the prior stay; many unlikely candidates lower the sample's rate and every
        # probability (context::SampleAdjusted).
        q = np.array([0.9, 0.1, 0.5, 0.5])
        adjusted, rate = self.trainer.sample_adjusted(q, 0.5)
        self.assertAlmostEqual(rate, 0.5)
        np.testing.assert_allclose(adjusted, q, atol=1e-9)
        deep = np.concatenate([[0.95, 0.9], np.full(198, 0.05)])
        adjusted, rate = self.trainer.sample_adjusted(deep, 0.3)
        self.assertLess(rate, 0.1)
        self.assertLess(adjusted[0], 0.95)

    def test_calls_at_a_target_as_protal(self):
        self.trainer.SINGLETON_CONGENER = 100  # the rule, off by default, vetoes a row here
        self.addCleanup(setattr, self.trainer, "SINGLETON_CONGENER", 0)
        # context::FalseCallKnob: the highest-scoring taxa while the mean of their 1 - probability is at most the target;
        # every taxon at the last one's score; the singleton rule's taxa never. Checked against that rule sample by sample.
        frame = pd.DataFrame({"meta_sample": ["a"] * 6 + ["b"] * 3,
                              "fragments": [5, 5, 5, 5, 1, 5, 5, 5, 5], "genus_top_fragments": [0, 0, 0, 0, 500, 0, 0, 0, 0],
                              "em_own_share": [1, 1, 1, 1, 0.2, 1, 1, 1, 1], "identity": [0.99] * 9})
        p = np.array([0.99, 0.95, 0.95, 0.4, 0.98, 0.1, 0.9, 0.2, 0.05])
        y = np.array([1, 1, 0, 0, 0, 0, 1, 0, 0])
        curve = [(0.0, 0.0), (1.0, 1.0)]
        prepared = self.trainer.FalseCallSamples(frame, y, p, curve, 0.5)
        for fdr in (0.001, 0.03, 0.1, 0.3, 0.6):
            calls = prepared.calls(fdr)
            self.assertFalse(calls[4], "the singleton rule's taxon")
            for sample in ("a", "b"):
                rows = np.flatnonzero((frame["meta_sample"] == sample).to_numpy() & (frame["fragments"] > 1).to_numpy())
                scores = np.sort(p[rows])[::-1]
                adjusted, _ = self.trainer.sample_adjusted(np.interp(scores, *zip(*curve)), 0.5)
                means = np.cumsum(1 - adjusted) / np.arange(1, len(scores) + 1)
                n = int((means <= fdr).sum())
                knob = scores[n - 1] if n else np.inf
                np.testing.assert_array_equal(calls[rows], p[rows] >= knob, f"fdr {fdr}, sample {sample}")
            tp, fp = prepared.counts(fdr)
            self.assertEqual((tp, fp), (int((calls & (y == 1)).sum()), int((calls & (y == 0)).sum())))
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
        model, metrics, stdout = self.train("fdr", "--fdr-calls", "--depth-knobs", "--singleton-congener", "100")
        calls = read_false_calls(model)
        self.assertIsNotNone(calls)
        self.assertEqual(calls["fdr"], metrics["false_calls"]["fdr"])
        self.assertIn(calls["fdr"], [float(v) for v in self.trainer.FDR_GRID])
        self.assertEqual([tuple(c) for c in metrics["false_calls"]["curve"]], calls["curve"])
        self.assertIn("test_false_calls", metrics)
        self.assertIn("singleton rule (one fragment beside a congener of 100 or more", stdout)
        self.assertIn("FP fdr", metrics["test_by_depth"][0])
        _, plain, _ = self.train("no_fdr")
        self.assertNotIn("false_calls", plain)
        self.assertIsNone(read_false_calls(os.path.join(self.tmp.name, "no_fdr.xml")))



if __name__ == "__main__":
    unittest.main()
