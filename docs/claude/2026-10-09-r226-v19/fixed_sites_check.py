#!/usr/bin/env python3
"""Is the ancestry report's fixed-site signal real, or does it know the strain? ancestry_sites.py (the build's report)
calls a site polymorphic where any of up to 6 other genomes of the species (the full reference) differs from the
representative, and fixed elsewhere. A missed real strain was simulated from one of those genomes, so its differences
from the representative are polymorphic by construction whenever its genome is among the 6; an in-silico strain's
species has one genome, so it has no polymorphic sites at all. This splits the FN taxa's own records by what was
simulated (the training and test tables' meta columns: the representative, a real strain, an in-silico strain) and
asks, per class:
  - how often the species has alleles at all (polymorphic sites on the record's copy);
  - the share of the record's mismatches at polymorphic sites, against the share of its bases there (a strain whose own
    genome is among the alleles has all its true differences there; a novel strain only by chance or by rate);
  - the records whose every mismatch (2 or more) falls on a polymorphic site, the signature of the strain's own copy;
  - the AUC of the report's signals for the class against the FP genus records, per taxon (>= 10 sites).

v19: the tables under work/, and each FN subset by alleles compared with the FP taxa of the same allele status (the
v18 script compared it with every FP taxon).

    python3 fixed_sites_check.py --build local/v19/protal0.7.9_r226_v19 > fixed_sites_check.txt
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

TABLES = {"pe": "training_data.tsv", "se": "training_data_se.tsv", "pb": "training_data_pb.tsv",
          "ont": "training_data_ont.tsv"}
SUMS = ["aligned", "mismatches", "sites1", "agree1", "alt1", "fixed1", "fixed_agree1", "fixed_alt1", "poly_covered",
        "poly_mismatches", "nonpoly_aligned", "nonpoly_mismatches"]


def auc(y, x):
    y, x = np.asarray(y, dtype=bool), np.asarray(x, dtype=float)
    n1, n0 = int(y.sum()), int((~y).sum())
    if n1 == 0 or n0 == 0:
        return float("nan")
    ranks = pd.Series(x).rank().to_numpy()
    return (ranks[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def classes(build, rt):
    """{(set:sample, taxid): class} of the present rows."""
    out = {}
    for which in ("training", "test"):
        t = pd.read_csv(os.path.join(build, "work", which, TABLES[rt]), sep="\t",
                        usecols=["meta_sample", "taxon", "truth", "meta_rep_genome", "meta_insilico_strain"])
        t = t[pd.to_numeric(t["truth"], errors="coerce") == 1]
        rep = pd.to_numeric(t["meta_rep_genome"], errors="coerce").fillna(0) == 1
        ins = pd.to_numeric(t["meta_insilico_strain"], errors="coerce").fillna(0) == 1
        cls = np.select([rep, ins], ["representative", "in-silico strain"], "real strain")
        for s, tx, c in zip(t["meta_sample"], t["taxon"], cls):
            out[(f"{which}:{s}", str(tx))] = c
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--build", required=True)
    p.add_argument("--types", default="pe,se,pb,ont")
    opts = p.parse_args(argv)
    pd.set_option("display.width", 250)
    for rt in opts.types.split(","):
        f = pd.read_csv(os.path.join(opts.build, "model_logs", "ancestry_sites", f"{rt}.fragments.tsv.gz"), sep="\t",
                        dtype={"taxid": str})
        cls = classes(opts.build, rt)
        f["class"] = [cls.get((s, t), "?") if r == "FN own" else f"{r} {rel}"
                      for s, t, r, rel in zip(f["sample"], f["taxid"], f["role"], f["relation"])]
        fn = f[f["role"] == "FN own"]
        print(f"\n# {rt}: {len(f)} records; FN own {len(fn)} by class {fn['class'].value_counts().to_dict()}\n")
        rows = []
        for name, g in f[f["class"].isin(["representative", "real strain", "in-silico strain", "FP genus"])].groupby("class"):
            s = g[SUMS].sum()
            has = g["poly_covered"] > 0
            m2 = g[(g["mismatches"] >= 2) & has]
            rows.append({
                "class": name, "records": len(g),
                "records with polymorphic sites": has.mean(),
                "bases at polymorphic sites": s["poly_covered"] / max(1, s["aligned"]),
                "mismatches at polymorphic sites": s["poly_mismatches"] / max(1, s["mismatches"]),
                "every mismatch polymorphic (>= 2, with sites)": (m2["poly_mismatches"] == m2["mismatches"]).mean(),
                "species' base at fixed sites": s["fixed_agree1"] / max(1, s["fixed1"]),
                "species' base at congener sites": s["agree1"] / max(1, s["sites1"]),
                "identity": 1 - s["mismatches"] / max(1, s["aligned"]),
                "fixed-site identity": 1 - s["nonpoly_mismatches"] / max(1, s["nonpoly_aligned"]),
            })
        print(pd.DataFrame(rows).set_index("class").T.round(4).to_string())
        # Per taxon, as the report does (>= 10 sites against the nearest congener): each FN class against FP genus.
        f["key"] = f["sample"] + "|" + f["taxid"] + "|" + f["class"]
        taxa = f[f["class"].isin(["representative", "real strain", "in-silico strain", "FP genus"])].groupby(
            ["key", "class"])[SUMS].sum().reset_index()
        taxa = taxa[taxa["sites1"] >= 10]
        taxa["identity"] = 1 - taxa["mismatches"] / taxa["aligned"].clip(lower=1)
        taxa["species_base_at_congener_sites"] = taxa["agree1"] / taxa["sites1"].clip(lower=1)
        taxa["species_base_at_fixed_sites"] = taxa["fixed_agree1"] / taxa["fixed1"].clip(lower=1)
        taxa["fixed_site_identity"] = 1 - taxa["nonpoly_mismatches"] / taxa["nonpoly_aligned"].clip(lower=1)
        taxa["has_alleles"] = taxa["poly_covered"] > 0
        fp = taxa[taxa["class"] == "FP genus"]
        out = []
        for name in ("representative", "real strain", "in-silico strain"):
            for alleles in (None, True, False):
                c = taxa[taxa["class"] == name]
                if alleles is not None:
                    c = c[c["has_alleles"] == alleles]
                if len(c) < 20:
                    continue
                f0 = fp if alleles is None else fp[fp["has_alleles"] == alleles]
                both = pd.concat([c, f0])
                y = (both["class"] == name).to_numpy()
                out.append({"FN class": name, "species with alleles": {None: "any", True: "yes", False: "no"}[alleles],
                            "taxa": len(c), **{sig: auc(y, both[sig]) for sig in
                                               ("identity", "species_base_at_congener_sites",
                                                "species_base_at_fixed_sites", "fixed_site_identity")}})
        print(f"\n  per taxon (>= 10 sites), each FN class against {len(fp)} FP genus taxa ({fp['has_alleles'].mean():.3f} "
              f"of them with alleles), AUC (FN high):")
        print(pd.DataFrame(out).round(3).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
