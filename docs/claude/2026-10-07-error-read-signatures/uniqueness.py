#!/usr/bin/env python3
"""Would a per-copy uniqueness score (how near the database's closest other copy of the gene is) see the beyond-genus
false positives? A copy that another genus's species shares near-identically shows in a read's ZA tag: that species
aligns too, with few more edits (ZA lists up to 4 taxa within 5 edits of the best). For the counted fragments of the
FP taxa (from the same genus as their source, or from beyond it) and the own reads, the alternatives of their best
records on the taxon: none, in the taxon's genus only, or in another genus (within 0, 1-2 or 3-5 edits more), the
genus of the read's own species among them, and the unique k-mers (ZU) the read has on the taxon.

    python3 uniqueness.py --records ~/v15/err_analysis --build local/v15 > uniqueness.txt
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

import fragments
import signatures


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--records", required=True)
    p.add_argument("--build", required=True)
    p.add_argument("--types", default="pe,se,pb,ont")
    return p.parse_args(argv)


def alternatives(r, f):
    """{(sample, qname): [(taxid, edits more)]} of the records on each fragment's best taxon (primary and
    supplementary: ZA is on those)."""
    key = dict(zip(zip(f["sample"], f["qname"]), f["taxon"]))
    out = {}
    za = r[r["za"].notna() & (r["taxon"].fillna("") != "")]
    for s, q, t, v in zip(za["sample"], za["qname"], za["taxon"].astype(str), za["za"]):
        if key.get((s, q)) != t or v in ("", "*"):
            continue
        lst = out.setdefault((s, q), [])
        for alt in v.split(","):
            a, _, more = alt.partition(":")
            try:
                lst.append((a, int(more)))
            except ValueError:
                pass
    return out


def classify(f, alts, tax):
    rows = []
    for s, q, t, src in zip(f["sample"], f["qname"], f["taxon"], f["source_taxid"]):
        a = alts.get((s, q), [])
        g_t = tax.ancestor(t, "genus")
        g_src = tax.ancestor(src, "genus") if src else None
        other = [m for x, m in a if tax.ancestor(x, "genus") != g_t]
        same = [m for x, m in a if tax.ancestor(x, "genus") == g_t]
        rows.append({
            "any alternative": bool(a),
            "other genus <= 0": any(m <= 0 for m in other),
            "other genus <= 2": any(m <= 2 for m in other),
            "other genus <= 5": bool(other),
            "own genus only": bool(same) and not other,
            "source's genus among them": g_src is not None and g_src != g_t
                                        and any(tax.ancestor(x, "genus") == g_src for x, _ in a)})
    return pd.DataFrame(rows, index=f.index)


def main(argv=None):
    opts = parse_args(argv)
    pd.set_option("display.width", 250)
    tax = fragments.Taxonomy(opts.build)
    for rt in opts.types.split(","):
        f = pd.read_pickle(os.path.join(opts.records, f"{rt}.fragments.pkl.gz"))
        f = f[f["counted"]].copy()
        f["group"] = signatures.group_of(f)
        far = f["relation"].isin(["family", "order", "other"])
        f["class"] = np.select([(f["group"] == "wrong (FP)") & (f["relation"] == "genus"),
                                (f["group"] == "wrong (FP)") & far,
                                f["group"].isin(["own (FN)", "own (TP)"])],
                               ["FP, from the genus", "FP, from beyond the genus", "own reads"], "")
        f = f[f["class"] != ""]
        r = fragments.load(opts.records, rt, "records")[["sample", "qname", "taxon", "za"]]
        c = classify(f, alternatives(r, f), tax)
        out = c.groupby(f["class"]).mean().round(3)
        out.insert(0, "fragments", f["class"].value_counts())
        out["ZU median"] = f.groupby("class")["zu"].median()
        out["identity median"] = f.groupby("class")["identity"].median().round(4)
        print(f"\n# {rt}\n")
        print(out.to_string())
        b = f[f["class"] == "FP, from beyond the genus"]
        if len(b):
            print(f"  beyond the genus, by source: held out {b['source_held_out'].mean():.3f}; another genus fits within "
                  f"2 edits: held-out sources {c.loc[b.index[b['source_held_out']], 'other genus <= 2'].mean():.3f}, "
                  f"sources in the database {c.loc[b.index[~b['source_held_out']], 'other genus <= 2'].mean():.3f}")
            # Recurrence: the same copy (taxon, gene) taking a beyond-genus read in other samples; in the test samples,
            # a copy that took one in a training sample (what a per-copy rate learned from the build's simulations sees)
            b = b.assign(copy=b["taxon"].astype(str) + "_" + b["gene"].astype("Int64").astype(str))
            per = b.groupby("copy")["sample"].nunique()
            print(f"  beyond-genus FP fragments on {len(per)} copies of {b['taxon'].nunique()} taxa; copies in 2 or more "
                  f"samples: {(per >= 2).mean():.3f} (holding {b['copy'].map(per).ge(2).mean():.3f} of the fragments); "
                  f"taxa in 2 or more samples: {b.groupby('taxon')['sample'].nunique().ge(2).mean():.3f}")
            train = set(b.loc[b["set"] == "training", "copy"])
            test = b[b["set"] == "test"]
            if len(test):
                print(f"  test-set beyond-genus FP fragments: {len(test)}; on a copy that took one in a training sample: "
                      f"{test['copy'].isin(train).mean():.3f}; on a taxon that did: "
                      f"{test['taxon'].isin(set(b.loc[b['set'] == 'training', 'taxon'])).mean():.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
