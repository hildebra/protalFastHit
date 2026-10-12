#!/usr/bin/env python3
"""Time one HistGradientBoostingClassifier fit with and without sample weights (0.25 for a fifth of the rows, as the
trainer weighs scenario rows), at a few table sizes, on one thread. Usage: weighted_binning.py ROWS[,ROWS...] FEATURES"""
import os
import sys
import time

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[v] = "1"
import numpy as np  # noqa: E402
import sklearn  # noqa: E402
from sklearn.ensemble import HistGradientBoostingClassifier  # noqa: E402

sizes = [int(s) for s in sys.argv[1].split(",")]
features = int(sys.argv[2])
print("sklearn", sklearn.__version__, "numpy", np.__version__, flush=True)
for n in sizes:
    rng = np.random.default_rng(1)
    X = rng.normal(size=(n, features))
    y = (X[:, 0] + rng.normal(0, 1, n) > 0).astype(int)
    w = np.where(rng.random(n) < 0.2, 0.25, 1.0)
    for weights in (None, w):
        model = HistGradientBoostingClassifier(max_iter=30, max_leaf_nodes=63, min_samples_leaf=20, early_stopping=False,
                                               random_state=1)
        t0 = time.time()
        model.fit(X, y, sample_weight=weights)
        print(f"rows {n} features {features} weights {'no ' if weights is None else 'yes'}: {time.time() - t0:.2f} s",
              flush=True)
