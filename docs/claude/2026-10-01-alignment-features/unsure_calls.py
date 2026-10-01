#!/usr/bin/env python3
"""The unsure calls (KNOB <= p < 0.8) of the base model: how many are absent and present, and how well each
signal separates them (AUC of absent vs present; 0.5 = no separation). Cross-validated training calls and test
calls of forest seed 1."""
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

import exp_lib as L
from run_rules import tables, BASE

pd.set_option("display.width", 220)
SIGNALS = ["congener_share_d0", "congener_share_d1", "congener_share_d2", "other_genus_share_d1", "identity_z",
           "identity", "top_identity", "fragments", "mean_mapq", "low_mapq_share", "mate_concordance", "clipped_share"]
rows = []
for rt in L.READ_TYPES:
    train, test = tables(rt)
    for split, df, p in (("training (CV)", train, L.cv_predict(train, BASE, 1)),
                         ("test", test, L.test_predict(train, test, BASE, 1))):
        u = df[(p >= L.KNOB) & (p < 0.8)]
        y = u["truth"].to_numpy()
        r = {"read_type": rt, "split": split, "calls": int((p >= L.KNOB).sum()), "unsure": len(u),
             "unsure_absent": int((y == 0).sum()), "unsure_present": int((y == 1).sum()),
             "FP_unsure_share": float(((p >= L.KNOB) & (p < 0.8) & (df["truth"] == 0)).sum() / max(1, ((p >= L.KNOB) & (df["truth"] == 0)).sum()))}
        for s in SIGNALS:
            if s in u.columns and len(np.unique(y)) == 2:
                # oriented so that > 0.5 means the signal is higher for absent taxa
                r[s] = roc_auc_score(1 - y, u[s].to_numpy())
        r["median_congener_d1_absent"] = float(u.loc[y == 0, "congener_share_d1"].median()) if (y == 0).any() else np.nan
        r["median_congener_d1_present"] = float(u.loc[y == 1, "congener_share_d1"].median()) if (y == 1).any() else np.nan
        rows.append(r)
out = pd.DataFrame(rows)
print("unsure calls (0.5 <= p < 0.8) and the AUC of each signal for absent vs present among them "
      "(> 0.5: higher for absent taxa; < 0.5: lower)\n")
print(out[["read_type", "split", "calls", "unsure", "unsure_absent", "unsure_present", "FP_unsure_share"]].to_string(index=False))
print()
print(out[["read_type", "split"] + [s for s in SIGNALS if s in out.columns] +
          ["median_congener_d1_absent", "median_congener_d1_present"]].round(3).to_string(index=False))
