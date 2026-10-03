#!/usr/bin/env python3
"""Per-read ambiguity of taxa beside a congener with >=10x more fragments: present (minor congeners) against absent
(spillover), test and training tables. Usage: ambiguity.py TABLES_DIR"""
import sys
import pandas as pd
pd.set_option("display.width", 220)
for name in ("training", "test"):
    t = pd.read_csv(f"{sys.argv[1]}/{name}.tsv", sep="\t", low_memory=False)
    t["truth"] = t["truth"].astype(str).str.lower().isin(("1", "true"))
    t["novel"] = t["meta_novel_level"].fillna("").astype(str).replace("nan", "") != ""
    low = t[t["genus_skew"] <= -1]
    print(f"\n{name}: taxa with a congener of >=10x more fragments, medians (quartiles)")
    for label, g in (("present", low[low.truth]), ("absent, near db species", low[~low.truth & ~low.novel]),
                     ("absent, near held-out", low[~low.truth & low.novel])):
        q = lambda c: f"{g[c].median():.3f} ({g[c].quantile(.25):.3f}, {g[c].quantile(.75):.3f})"
        print(f"  {label:26s} n={len(g):6d}  congener_fit_share {q('congener_fit_share')}  identity {q('identity')}  "
              f"excess_median {q('excess_median')}  low_mapq_share {q('low_mapq_share')}  fragments {q('fragments')}")
