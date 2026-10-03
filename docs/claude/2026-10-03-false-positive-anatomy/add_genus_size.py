#!/usr/bin/env python3
"""genus_species = log10 of the number of species of the taxon's genus in the database (internal_taxonomy.dmp): how
speciose the neighbourhood is, a per-species constant the build could store. Usage: add_genus_size.py TAXONOMY IN OUT RFDIR"""
import sys
from collections import Counter
import numpy as np
import pandas as pd
sys.path.insert(0, sys.argv[4])
import lineages as L
lin, _ = L.from_taxonomy(sys.argv[1])
species_per_genus = Counter(l["genus"] for l in lin.values() if "species" in l and "genus" in l)
df = pd.read_csv(sys.argv[2], sep="\t", float_precision="round_trip", low_memory=False)
genus = df["taxon"].astype(str).map(lambda t: lin.get(t, {}).get("genus"))
df["genus_species"] = np.log10(genus.map(lambda g: species_per_genus.get(g, 1) if g else 1).astype(float))
print("genus sizes (species): quantiles", np.quantile(10 ** df["genus_species"], [0.1, 0.5, 0.9]), file=sys.stderr)
df.to_csv(sys.argv[3], sep="\t", index=False)
