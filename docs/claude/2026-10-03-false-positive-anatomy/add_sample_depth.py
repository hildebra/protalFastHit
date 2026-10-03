#!/usr/bin/env python3
"""sample_log_fragments = log10 of the sample's fragments over all its taxa (rows), at least 1: the depth protal reads
the knob curve at (profiler::DepthKnobAt, random_forest_cmdline.sample_depths), here as a feature of every row.
Usage: add_sample_depth.py IN OUT"""
import sys
import numpy as np
import pandas as pd
df = pd.read_csv(sys.argv[1], sep="\t", float_precision="round_trip", low_memory=False)
totals = df["fragments"].astype(float).groupby(df["meta_sample"].astype(str)).transform("sum")
df["sample_log_fragments"] = np.log10(np.maximum(totals, 1.0))
df.to_csv(sys.argv[2], sep="\t", index=False)
