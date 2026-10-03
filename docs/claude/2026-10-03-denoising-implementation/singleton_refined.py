#!/usr/bin/env python3
"""The singleton rule as is (one fragment beside a congener of 100 or more) against a refined one that also needs the
read to look like the congener's (its EM share below 0.5, or its identity below 0.95): rows vetoed, present among them,
and what a model would have called at 0.5 among them, on both test sets and read types."""
import sys
import pandas as pd
OUT, E = sys.argv[1:3]
for test in ("base", "new"):
    for rt in ("", "_se"):
        t = pd.read_csv(f"{OUT}/{test}/test/training_data{rt}.tsv", sep="\t", low_memory=False)
        p = pd.read_csv(f"{E}/new_normalized-adjacency_{test}{rt}.test_predictions.tsv.gz", sep="\t",
                        usecols=["meta_sample", "taxon", "p"])
        d = t.merge(p, on=["meta_sample", "taxon"])
        truth = d["truth"].astype(str).str.lower().isin(("1", "true"))
        single = (d["fragments"] <= 1) & (d["genus_top_fragments"] >= 100)
        refined = single & ((d["em_own_share"] < 0.5) | (d["identity"] < 0.95))
        called = d["p"] >= 0.5
        for name, veto in (("count only", single), ("refined", refined)):
            print(f"test {test:4s} {rt or 'pe':3s} {name:10s}: vetoed {int(veto.sum()):4d} ({int((veto & truth).sum())} present); "
                  f"of the calls at 0.5 removed {int((veto & called).sum())}: {int((veto & called & ~truth).sum())} false, "
                  f"{int((veto & called & truth).sum())} true")
