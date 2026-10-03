#!/usr/bin/env python3
"""Add abundance-relative-to-relatives features (the UNOISE "skew" / DADA2 "is it more than a more abundant
neighbour's errors explain" idea) to protal training dumps, from the rows of the same sample.

For each taxon k of a sample, with n_k its fragments:
  genus_top       the most fragments of another species of k's genus in the sample (0 if none)
  family_top      the most fragments of a species of another genus of k's family in the sample (0 if none)
  genus_skew      log10((n_k + 1) / (genus_top + 1))
  family_skew     log10((n_k + 1) / (family_top + 1))
  genus_share     n_k over the fragments of k's genus in the sample
  genus_spill     log10((n_k + 0.5) / (lambda_g * genus_sum_other + lambda_f * family_sum_other + 0.5)): the reads seen
                  over the reads expected to spill over from the sample's relatives, at fixed rates per read of the
                  relative (DADA2's n_j * lambda_ji, with lambda by rank only)
All computable by protal at run time: fragments of every taxon of the sample and the database's taxonomy.

Usage: add_relative_features.py TAXONOMY IN.tsv OUT.tsv [LAMBDA_GENUS [LINEAGES_DIR]]
(LAMBDA_GENUS default 1e-3, empty for the default; the family's is a tenth of it; LINEAGES_DIR holds lineages.py.)
"""
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, sys.argv[5] if len(sys.argv) > 5 else ".")
import lineages  # noqa: E402

taxonomy, inp, out = sys.argv[1:4]
lam_g = float(sys.argv[4]) if len(sys.argv) > 4 and sys.argv[4] else 1e-3
lam_f = lam_g / 10

df = pd.read_csv(inp, sep="\t", float_precision="round_trip", low_memory=False)
by_id, _ = lineages.from_taxonomy(taxonomy)
taxa = df["taxon"].astype(str)
genus = taxa.map({t: by_id.get(t, {}).get("genus", "g?" + t) for t in taxa.unique()})
family = taxa.map({t: by_id.get(t, {}).get("family", "f?" + t) for t in taxa.unique()})
n = df["fragments"].astype(float)
sample = df["meta_sample"].astype(str)

key_g = sample + "|" + genus
key_f = sample + "|" + family
# max over the OTHER members: the group's max, unless k is it, then the second largest
def other_max(keys, values):
    frame = pd.DataFrame({"k": keys, "v": values})
    srt = frame.sort_values("v", ascending=False)
    first = srt.groupby("k")["v"].transform(lambda s: s.iloc[0])
    second = srt.groupby("k")["v"].transform(lambda s: s.iloc[1] if len(s) > 1 else 0.0)
    rank = srt.groupby("k").cumcount()
    res = np.where(rank == 0, second, first)
    return pd.Series(res, index=srt.index).reindex(frame.index).to_numpy()

genus_top = other_max(key_g, n)
genus_sum = n.groupby(key_g).transform("sum").to_numpy()
# family, other genera only: per (sample, family) max/sum minus the own genus' contribution
fam_genus = pd.DataFrame({"kf": key_f, "kg": key_g, "v": n})
gsum = fam_genus.groupby(["kf", "kg"])["v"].sum().rename("gsum").reset_index()
gmax = fam_genus.groupby(["kf", "kg"])["v"].max().rename("gmax").reset_index()
g = gsum.merge(gmax, on=["kf", "kg"])
fam_sum = g.groupby("kf")["gsum"].transform("sum")
# the largest genus max in the family other than the own genus
def top_other_genus(frame):
    vals = frame["gmax"].to_numpy()
    order = np.argsort(-vals)
    best, second = vals[order[0]], (vals[order[1]] if len(vals) > 1 else 0.0)
    return pd.Series(np.where(np.arange(len(vals)) == order[0], second, best), index=frame.index)
g["fam_other_top"] = g.groupby("kf", group_keys=False)[["gmax"]].apply(top_other_genus)
g["fam_other_sum"] = fam_sum - g["gsum"]
g = g.set_index(["kf", "kg"])
idx = pd.MultiIndex.from_arrays([key_f, key_g])
family_top = g["fam_other_top"].reindex(idx).to_numpy()
family_sum_other = g["fam_other_sum"].reindex(idx).to_numpy()

df["genus_top"] = genus_top
df["family_top"] = family_top
df["genus_skew"] = np.log10((n + 1) / (genus_top + 1))
df["family_skew"] = np.log10((n + 1) / (family_top + 1))
df["genus_share"] = n / np.maximum(genus_sum, 1)
expected = lam_g * (genus_sum - n) + lam_f * family_sum_other
df["genus_spill"] = np.log10((n + 0.5) / (expected + 0.5))
df.to_csv(out, sep="\t", index=False)
print(f"{out}: {len(df)} rows, {sample.nunique()} samples; lambda genus {lam_g}, family {lam_f}")
