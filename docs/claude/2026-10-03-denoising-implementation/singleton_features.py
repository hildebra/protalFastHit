#!/usr/bin/env python3
"""The singleton rule's rows (one fragment beside a congener of 100 or more) of both pipelines' training and test tables:
present against absent, by identity, divergence beyond qualities, the EM share, the close share and the distance."""
import sys
import pandas as pd
OUT = sys.argv[1]
cols = ["identity", "top_identity", "excess_median", "em_own_share", "relative_close_share", "relative_distance", "genus_top_fragments"]
for pipeline in ("base", "new"):
    for kind in ("training", "test"):
        t = pd.read_csv(f"{OUT}/{pipeline}/{kind}/training_data.tsv", sep="\t", low_memory=False)
        truth = t["truth"].astype(str).str.lower().isin(("1", "true"))
        veto = (t["fragments"] <= 1) & (t["genus_top_fragments"] >= 100)
        print(f"\n{pipeline} {kind}: vetoed {int(veto.sum())}, present {int((veto & truth).sum())}")
        for label, sel in (("present", veto & truth), ("absent", veto & ~truth)):
            if sel.any():
                q = t.loc[sel, cols].quantile([0.1, 0.5, 0.9]).round(3)
                print(f"  {label}: " + "; ".join(f"{c} {q[c].tolist()}" for c in cols))
