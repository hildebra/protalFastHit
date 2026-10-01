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
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from model_pmml import PmmlForest, float32_split, format_depth_knobs, read_depth_knobs, write_forest  # noqa: E402

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

    def test_depth_knobs(self):
        # Written as an Extension in the header, before the Application (PMML's order), and read back.
        rf = RandomForestClassifier(n_estimators=2, random_state=1).fit(self.X, self.y)
        path = os.path.join(self.tmp.name, "knobs.xml")
        write_forest(rf, self.features, path, depth_knobs={4: 0.37, 2: 0.6})
        with open(path) as fh:
            xml = fh.read()
        self.assertIn('<Extension name="protal_depth_knobs" value="2:0.6,4:0.37"/>\n  <Application ', xml)
        self.assertEqual(read_depth_knobs(path), {2: 0.6, 4: 0.37})
        np.testing.assert_array_equal(PmmlForest(path).predict(self.X), rf.predict_proba(self.X)[:, 1])
        write_forest(rf, self.features, path)
        self.assertEqual(read_depth_knobs(path), {})
        self.assertEqual(format_depth_knobs({6: 0.05, 3: 0.5}), "3:0.5,6:0.05")

    def test_rejects_other_classes(self):
        rf = RandomForestClassifier(n_estimators=2, random_state=1).fit(self.X, np.where(self.y == 1, "yes", "no"))
        with self.assertRaises(ValueError):
            write_forest(rf, self.features, os.path.join(self.tmp.name, "bad.xml"))


@unittest.skipIf(RandomForestClassifier is None, "needs scikit-learn")
class TrainerDepthKnobsTest(unittest.TestCase):
    """random_forest_cmdline.py --depth-knobs on a table of shallow samples (hundreds of fragments, depth bin 2) and
    deep ones (tens of thousands, bin 4), where a present taxon's evidence grows with depth."""

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
                rows.append({"truth": int(present), "taxon": 100 + (taxon + sample) % 60, "taxon_name": "t",
                             "meta_sample": f"s{sample}", "fragments": float(fragments),
                             "x": signal + rng.normal(0, 0.5), "y": rng.normal(0, 1)})
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

    def test_bins_as_protal(self):
        # profiler::DepthKnobBin: the digits of the sample's fragments over all its taxa, less one, 2 to 6.
        frame = pd.DataFrame({"meta_sample": ["a", "a", "b", "c", "d", "e"],
                              "fragments": [60.0, 39.0, 1000.0, 999.0, 3e7, 0.0]})
        self.assertEqual(list(self.trainer.depth_bins(frame)), [2, 2, 3, 2, 6, 2])

    def test_depth_knobs_in_the_model(self):
        model, metrics = self.train("knobs", "--depth-knobs")
        knobs = read_depth_knobs(model)
        self.assertEqual(sorted(knobs), [2, 4])
        self.assertTrue(all(0.05 <= k <= 0.95 for k in knobs.values()), knobs)
        self.assertEqual(metrics["depth_knobs"]["knobs"], {str(b): k for b, k in knobs.items()})
        # The knobs are chosen on the species held out: their F1 there is at least that at the one knob.
        self.assertGreaterEqual(metrics["depth_knobs"]["F1_at_depth_knobs"], metrics["depth_knobs"]["F1_at_knob"])
        plain, metrics = self.train("plain")
        self.assertEqual(read_depth_knobs(plain), {})
        self.assertNotIn("depth_knobs", metrics)


if __name__ == "__main__":
    unittest.main()
