#!/usr/bin/env python3
"""The soil errors in their genus: what else of the taxon's genus the sample holds (soil_errors.py --save rows).

For each false negative: whether its genus has a present taxon of the database called in the same sample (a congener
that takes its reads), how many fragments that congener has against it, and whether a false positive of its genus sits
beside it (the reads went to the wrong species of the right genus). For each false positive: whether its genus has a
present taxon of the database in the sample (TP or FN), else only the novel species whose reads it took
(meta_relative_rank genus). Strains: identity and divergence of real against in-silico strains and representatives.
Genus-level F1: the same calls scored by genus (a genus called if any of its species is called; present if a species
of it is simulated in the sample, which meta_relative_rank says for absent rows too).

    python3 soil_genus_context.py --rows v13_rows.tsv.gz
"""

import argparse

import numpy as np
import pandas as pd


def genus_of(name):
    s = str(name)
    s = s[3:] if s.startswith("s__") else s
    return s.split(" ")[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rows", required=True)
    opts = ap.parse_args()
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    df = pd.read_csv(opts.rows, sep="\t", low_memory=False)
    df["genus"] = df.taxon_name.map(genus_of)
    df["group"] = df.scenario + " " + df.set
    groups = [g for g in ["soil hold-in species held out", "soil hold-out", "soil_shallow hold-in species held out",
                          "soil_shallow hold-out", "(design) design test"] if g in set(df.group)]
    for rt in ["pe", "se", "pb", "ont"]:
        for grp in groups:
            g = df[(df.read_type == rt) & (df.group == grp)].copy()
            if not len(g):
                continue
            key = ["meta_sample", "genus"]
            g["called"] = g.call.astype(bool)
            agg = g.groupby(key).agg(
                tp=("outcome", lambda s: (s == "TP").sum()), fp=("outcome", lambda s: (s == "FP").sum()),
                fn=("outcome", lambda s: (s == "FN").sum()), present=("truth", "sum"),
                rel_genus=("meta_relative_rank", lambda s: (s.astype(str) == "genus").any()),
                max_tp_frag=("fragments", "max"))
            # the most fragments of a called present congener in the genus and sample
            called_present = g[(g.outcome == "TP")].groupby(key).fragments.max().rename("tp_frag_max")
            g = g.join(agg[["tp", "fp", "fn", "rel_genus"]], on=key).join(called_present, on=key)
            fn = g[g.outcome == "FN"]
            fp = g[g.outcome == "FP"]
            print(f"\n=== {rt} {grp}: FN {len(fn)}, FP {len(fp)}")
            if len(fn):
                with_tp = fn.tp > 0
                ratio = (fn.tp_frag_max / fn.fragments.clip(lower=1))
                print(f"FN with a called present congener in the sample: {with_tp.mean():.1%} "
                      f"(its fragments / the FN's: median {ratio[with_tp].median():.1f}, >=10x {(ratio[with_tp] >= 10).mean():.1%}); "
                      f"with a false positive of its genus beside it: {(fn.fp > 0).mean():.1%}; "
                      f"neither (alone in its genus among the called): {((fn.tp == 0) & (fn.fp == 0)).mean():.1%}")
                print("FN by strain kind x a called present congener (count):")
                print(pd.crosstab(fn.strain, with_tp.map({True: "congener called", False: "no congener called"})).to_string())
            if len(fp):
                db_present = (fp.tp + fp.fn) > 0
                print(f"FP with a present database congener in the sample: {db_present.mean():.1%} "
                      f"(called {(fp.tp > 0).mean():.1%}); only a novel species of its genus: "
                      f"{(~db_present & fp.rel_genus).mean():.1%}; neither: {(~db_present & ~fp.rel_genus).mean():.1%}")
            # genus-level scoring: a genus is present if any row says a simulated species shares the genus
            gg = g.groupby(key).agg(called=("called", "any"), present=("truth", "max"),
                                    rel=("meta_relative_rank", lambda s: (s.astype(str).isin(["genus", "species"])).any()))
            gg["present"] = (gg.present > 0) | gg.rel
            tp_, fp_, fn_ = (gg.called & gg.present).sum(), (gg.called & ~gg.present).sum(), (~gg.called & gg.present).sum()
            print(f"genus level: {len(gg)} genera with reads, present {int(gg.present.sum())}, TP {tp_}, FP {fp_}, FN {fn_}, "
                  f"F1 {2 * tp_ / (2 * tp_ + fp_ + fn_):.4f} (genera whose only evidence is a novel species: "
                  f"{int((gg.rel & ~(g.groupby(key).truth.max() > 0)).sum())})")
            pres = g[g.truth == 1]
            cols = ["fragments", "identity", "top_identity", "excess_median", "excess_scaled_median",
                    "variant_sites_per_kb", "low_mapq_share", "congener_fit_share", "cluster_ani_radius",
                    "cluster_min_ani", "cluster_genomes_log10", "p"]
            cols = [c for c in cols if c in pres]
            print("present taxa, medians by strain kind and outcome:")
            print(pres.groupby(["strain", "outcome"])[cols].median().round(4).to_string())
            print("present taxa, FN rate % by strain kind and identity:")
            ib = pd.cut(pres.identity, [0, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995, 1.001])
            print(pres.pivot_table(index=ib, columns="strain", values="outcome", observed=True,
                                   aggfunc=lambda s: f"{100 * (s == 'FN').mean():.1f} ({len(s)})").to_string())
            print("present taxa by species cluster's genomes in GTDB (log10) and strain kind: FN rate % (rows)")
            cb = pd.cut(pres.cluster_genomes_log10, [-0.01, 0.01, 0.5, 1, 2, 5])
            print(pres.pivot_table(index=cb, columns="strain", values="outcome", observed=True,
                                   aggfunc=lambda s: f"{100 * (s == 'FN').mean():.1f} ({len(s)})").to_string())


if __name__ == "__main__":
    main()
