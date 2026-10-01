#!/usr/bin/env python3
"""Extracts and caches the SAM-derived features (~/fpexp/features/): from the run's own SAMs (primary) and, with
--secondary, from the realignment with -m 3.

    python3 prep_features.py [--secondary] [splits...]
"""
import sys
import time

import exp_lib as L

secondary = "--secondary" in sys.argv
splits = [a for a in sys.argv[1:] if not a.startswith("--")] or ["test", "training"]
for split in splits:
    for rt in L.READ_TYPES:
        t = time.time()
        df = L.secondary_features(split, rt) if secondary else L.primary_features(split, rt)
        cols = [c for c in df.columns if c not in ("meta_sample", "taxon")]
        print(f"{split} {rt}: {len(df)} taxa in {time.time() - t:.0f} s; means: "
              + ", ".join(f"{c} {df[c].mean():.3f}" for c in cols), flush=True)
