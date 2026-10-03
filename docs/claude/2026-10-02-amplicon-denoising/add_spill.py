#!/usr/bin/env python3
"""spill_divergent = max(0, -genus_skew) * max(0, excess_median): the decades by which a congener in the sample
outnumbers the taxon, counted only as far as the taxon's reads differ from its reference beyond their errors (a
real minor congener's reads fit it; a relative's spill-over does not). Usage: add_spill.py IN OUT"""
import sys
import numpy as np
import pandas as pd
df = pd.read_csv(sys.argv[1], sep="\t", float_precision="round_trip", low_memory=False)
df["spill_divergent"] = np.maximum(0, -df["genus_skew"]) * np.maximum(0, df["excess_median"])
df.to_csv(sys.argv[2], sep="\t", index=False)
