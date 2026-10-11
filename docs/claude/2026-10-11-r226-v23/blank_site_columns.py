"""A copy of v23's training and test tables with the columns that need the ancestry sites set to -1 (unknown) in every
row, for refitting the full feature set without them (ablate_columns.sh): the ancestry group, polymorphic's
ancestry_fixed_gain and ancestry_fixed_agreement, and weights' ancestry_agreement_weighted. A constant column gives the
trees no split, so the refit is the full set less these eight; the other columns of the three groups stay.

usage: blank_site_columns.py <tables dir> <out dir> <read type>
"""
import os
import sys
import pandas as pd

SITE_COLUMNS = ["ancestry_sites_per_record", "ancestry_agreement", "ancestry_congener_share",
                "ancestry_indel_sites_per_record", "ancestry_indel_congener_share",
                "ancestry_fixed_gain", "ancestry_fixed_agreement", "ancestry_agreement_weighted"]

src, out, t = sys.argv[1], sys.argv[2], sys.argv[3]
suf = "" if t == "pe" else "_" + t
for d in ("training", "test"):
    os.makedirs(f"{out}/{d}", exist_ok=True)
    tab = pd.read_csv(f"{src}/{d}/training_data{suf}.tsv", sep="\t", low_memory=False)
    missing = [c for c in SITE_COLUMNS if c not in tab.columns]
    if missing:
        raise SystemExit(f"columns missing: {missing}")
    tab[SITE_COLUMNS] = -1
    tab.to_csv(f"{out}/{d}/training_data{suf}.tsv", sep="\t", index=False)
    print(f"{d} {t}: {len(tab)} rows, {len(SITE_COLUMNS)} columns set to -1")
