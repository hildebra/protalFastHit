#!/usr/bin/env python3
"""What the singleton rule removes on the two test sets: rows it vetoes, how many a model would have called at 0.5
without it, and how many of those are present. Usage: singleton.py OUT EVAL"""
import os
import sys
import numpy as np
import pandas as pd
OUT, E = sys.argv[1:3]
for model in ("new_normalized-adjacency", "new_normalized-adjacency-relatives", "base_normalized-adjacency"):
    for test in ("base", "new"):
        for rt in ("", "_se"):
            t = pd.read_csv(f"{OUT}/{test}/test/training_data{rt}.tsv", sep="\t", low_memory=False,
                            usecols=["meta_sample", "taxon", "truth", "fragments", "genus_top_fragments"])
            p = pd.read_csv(f"{E}/{model}_{test}{rt}.test_predictions.tsv.gz", sep="\t", usecols=["meta_sample", "taxon", "p"])
            d = t.merge(p, on=["meta_sample", "taxon"])
            truth = d["truth"].astype(str).str.lower().isin(("1", "true"))
            veto = (d["fragments"] <= 1) & (d["genus_top_fragments"] >= 100)
            called = d["p"] >= 0.5
            print(f"{model:36s} test {test:4s} {rt or 'pe':3s}: vetoed {int(veto.sum()):5d} rows ({int((veto & truth).sum())} present); "
                  f"called at 0.5 without the rule {int((veto & called).sum())} ({int((veto & called & truth).sum())} present)")
