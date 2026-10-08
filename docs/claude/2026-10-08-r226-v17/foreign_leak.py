#!/usr/bin/env python3
"""Is foreign_scanned_share the simulation pool in disguise? The foreign scan read the genomes at hand (the simulation
pool, held-out species left out), so a species has scanned copies when its genome was downloaded. Present taxa are in
the pool by construction. Among the absent rows: those of species that are present in some other sample (so in the
pool) against those of species never present in any sample (mostly not in the pool).

    python3 foreign_leak.py --build local/v17 > foreign_leak.txt
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

TABLES = {"pe": "training_data.tsv", "se": "training_data_se.tsv", "pb": "training_data_pb.tsv",
          "ont": "training_data_ont.tsv"}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--build", required=True)
    opts = p.parse_args(argv)
    for rt, name in TABLES.items():
        frames = [pd.read_csv(os.path.join(opts.build, w, name), sep="\t",
                              usecols=["taxon", "truth", "fragments", "foreign_scanned_share", "foreign_copy_share",
                                       "meta_novel_level", "meta_relative_rank"])
                  for w in ("training", "test") if os.path.isfile(os.path.join(opts.build, w, name))]
        t = pd.concat(frames, ignore_index=True)
        t["truth"] = pd.to_numeric(t["truth"], errors="coerce").fillna(0).astype(int)
        ever_present = set(t.loc[t["truth"] == 1, "taxon"])
        a = t[t["truth"] == 0].copy()
        a["pool"] = a["taxon"].isin(ever_present)
        for col in ("foreign_scanned_share", "foreign_copy_share"):
            a[col] = pd.to_numeric(a[col], errors="coerce")
        print(f"{rt}: {len(t)} rows, {len(ever_present)} species ever present; absent rows {len(a)}, of species ever present "
              f"{a['pool'].mean():.3f}")
        g = a.groupby("pool")
        print(f"  foreign_scanned_share of absent rows: species ever present (in the pool) median "
              f"{g['foreign_scanned_share'].median().get(True, float('nan')):.3f}, mean {g['foreign_scanned_share'].mean().get(True, float('nan')):.3f}; "
              f"never present: median {g['foreign_scanned_share'].median().get(False, float('nan')):.3f}, mean "
              f"{g['foreign_scanned_share'].mean().get(False, float('nan')):.3f}")
        novel = a[(a["meta_novel_level"].astype(str) == "species") & (a["meta_relative_rank"].astype(str) == "genus")]
        print(f"  absent novel-congener rows: in the pool {novel['pool'].mean():.3f}; scanned share mean in pool "
              f"{novel.loc[novel['pool'], 'foreign_scanned_share'].mean():.3f}, out of the pool "
              f"{novel.loc[~novel['pool'], 'foreign_scanned_share'].mean():.3f}")
        # Which tells presence better among the absent-vs-present rows: the pool itself, or the feature?
        t["pool"] = t["taxon"].isin(ever_present)
        print(f"  all rows: foreign_scanned_share > 0.5 agrees with 'species in the pool' in "
              f"{np.mean((pd.to_numeric(t['foreign_scanned_share'], errors='coerce') > 0.5) == t['pool']):.4f} of the rows")
    return 0


if __name__ == "__main__":
    sys.exit(main())
