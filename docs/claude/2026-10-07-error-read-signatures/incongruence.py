#!/usr/bin/env python3
"""The FP fragments' gene copies against the build's near pairs across genera (model_logs/gene_incongruence.tsv: every
two copies of a gene in different genera within 0.05 k-mer distance, with the suspect verdicts): whether the copies the
beyond-genus false positives land on have a near copy of another genus, how near, and whether the read is any closer
to the copy than that other genus's copy is, against the copies the same-genus FP fragments and the species' own reads
land on. v15's logs archive lacks the table; v14's build has the same reference genes.

    python3 incongruence.py --records ~/v15/err_analysis --incongruence local/v14/model_logs/gene_incongruence.tsv > incongruence.txt
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

import signatures
from fragments import Taxonomy

RANK = {"genus": 0, "family": 1, "order": 2, "class": 3, "phylum": 4, "domain": 5, "none": 6}


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--records", required=True)
    p.add_argument("--incongruence", required=True)
    p.add_argument("--build", required=True, help="the build's local folder (taxonomy, held-out species)")
    p.add_argument("--types", default="pe,se,pb,ont")
    return p.parse_args(argv)


def congeners_in_db(build):
    """{species taxid: the number of other species of its genus in the training database (the taxonomy's species
    less the held-out ones)}."""
    tax = Taxonomy(build)
    held = {tax.name_id[n] for n in tax.held_out if n in tax.name_id}
    genus_of = {t: tax.ancestor(t, "genus") for t, r in tax.rank.items() if r == "species" and t not in held}
    per_genus = pd.Series([g for g in genus_of.values() if g]).value_counts()
    return {t: int(per_genus.get(g, 1)) - 1 for t, g in genus_of.items() if g}


def nearest(path):
    """{(geneid, taxid): (distance to the nearest copy of another genus within 0.05, its shared rank, suspect)}."""
    t = pd.read_csv(path, sep="\t", dtype={"taxid_a": str, "taxid_b": str})
    a = t[["geneid", "taxid_a", "shared_rank", "distance", "suspect_a"]].rename(
        columns={"taxid_a": "taxid", "suspect_a": "suspect"})
    b = t[["geneid", "taxid_b", "shared_rank", "distance", "suspect_b"]].rename(
        columns={"taxid_b": "taxid", "suspect_b": "suspect"})
    both = pd.concat([a, b], ignore_index=True).sort_values("distance")
    first = both.drop_duplicates(["geneid", "taxid"], keep="first")
    sus = both.groupby(["geneid", "taxid"])["suspect"].max()
    out = {}
    for g, x, r, d in zip(first["geneid"], first["taxid"], first["shared_rank"], first["distance"]):
        out[(int(g), x)] = (d, r, int(sus.get((g, x), 0)))
    return out


def main(argv=None):
    opts = parse_args(argv)
    pd.set_option("display.width", 250)
    near = nearest(opts.incongruence)
    print(f"copies with a near copy of another genus within 0.05: {len(near)}; of them suspect: "
          f"{sum(v[2] for v in near.values())}")
    congeners = congeners_in_db(opts.build)
    for rt in opts.types.split(","):
        f = pd.read_pickle(os.path.join(opts.records, f"{rt}.fragments.pkl.gz"))
        f["group"] = signatures.group_of(f)
        c = f[f["counted"] & f["gene"].notna()].copy()
        c["g"] = ["FP, same genus" if (grp == "wrong (FP)" and rel == "genus") else
                  "FP, beyond the genus" if (grp == "wrong (FP)" and rel in ("family", "order", "other")) else
                  "own reads (FN taxa)" if grp == "own (FN)" else "own reads (TP taxa)" if grp == "own (TP)" else ""
                  for grp, rel in zip(c["group"], c["relation"])]
        c = c[c["g"] != ""]
        key = [near.get((int(g), str(t))) for g, t in zip(c["gene"], c["taxon"])]
        c["near_d"] = [k[0] if k else np.nan for k in key]
        c["near_rank"] = [k[1] if k else "none within 0.05" for k in key]
        c["div"] = 1 - c["identity"]
        c["read_not_closer"] = c["div"] >= c["near_d"]  # the read is no nearer the copy than the other genus's copy
        print(f"\n# {rt}\n")
        rows = []
        for g, s in c.groupby("g"):
            per = s.groupby(["sample", "taxon"]).size()
            rows.append({"fragments of": g, "fragments": len(s), "taxa": len(per),
                         "copy has a near copy of another genus (<= 0.05)": round(s["near_d"].notna().mean(), 3),
                         "<= 0.02": round((s["near_d"] <= 0.02).mean(), 3),
                         "<= 0.01": round((s["near_d"] <= 0.01).mean(), 3),
                         "read no closer than that copy": round(s["read_not_closer"].fillna(False).mean(), 3),
                         "median read divergence": round(s["div"].median(), 4),
                         "median near distance (where any)": round(s["near_d"].median(), 4)})
        signatures.show(pd.DataFrame(rows), index=False)
        far = c[c["g"] == "FP, beyond the genus"]
        if len(far):
            print(f"  beyond-genus FP fragments by the near copy's rank: {far['near_rank'].value_counts().to_dict()}")
            n = far["taxon"].map(congeners)
            print(f"  the FP taxon's congeners in the training database (a congruence test needs some): none "
                  f"{(n == 0).mean():.3f}, 1-2 {n.between(1, 2).mean():.3f}, 3 or more {(n >= 3).mean():.3f}; of the "
                  f"fragments on a copy without a near foreign copy: none {(n[far['near_d'].isna()] == 0).mean():.3f}")
            print(f"  by the read's relation to the source: {pd.crosstab(far['relation'], far['near_d'].notna()).to_dict()}")
            # Would a per-copy rule catch them: the read's divergence against the copy's nearest foreign distance
            for m in (1.0, 0.5, 0.25):
                hit = (far["div"] >= m * far["near_d"]).fillna(False)
                own = c[c["g"].str.startswith("own")]
                own_hit = (own["div"] >= m * own["near_d"]).fillna(False)
                print(f"  rule 'read divergence >= {m} x the copy's nearest foreign distance': catches "
                      f"{hit.mean():.3f} of beyond-genus FP fragments, {((c['g'] == 'FP, same genus') & (c['div'] >= m * c['near_d']).fillna(False)).sum() / max(1, (c['g'] == 'FP, same genus').sum()):.3f} of same-genus FP, "
                      f"costs {own_hit.mean():.3f} of own reads")
    return 0


if __name__ == "__main__":
    sys.exit(main())
