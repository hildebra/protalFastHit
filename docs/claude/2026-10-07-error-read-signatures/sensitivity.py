#!/usr/bin/env python3
"""Whether the FP and FN taxa come from reads the aligner failed on (sensitivity) or from reads it aligned to the
wrong taxon (assignment), from parse_error_records.py's and fragments.py's tables.

FN taxa: the fate of their own fragments. The taxa tables count every fragment with a record (own_best_on_taxon,
own_best_elsewhere, own_unaligned); a read that seeded on nothing has no record. The v15 SAMs hold every one of them
(error_reads.py of 652ed53 had no cap), and their ZF tags tell which of the elsewhere and unaligned fragments had
seeded on the taxon itself and failed to align there: the reads the aligner lost. Paired-end tables also give the read
pairs simulated from the species: own fragments per read pair by its genome (representative, another genome, an
in-silico strain) shows whether divergent genomes lose reads before seeding.

FP taxa: whether the reads' source species is in the database at all, and if so whether they seeded on it and failed
(ZF), aligned to it too (a secondary record or ZA) or never touched it.

    python3 sensitivity.py --records ~/v15/err_analysis --build local/v15 > sensitivity.txt
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

import fragments
from fragments import Taxonomy


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--records", required=True)
    p.add_argument("--build", required=True)
    p.add_argument("--types", default="pe,se,pb,ont")
    return p.parse_args(argv)


def touched(r, keys):
    """{(sample, qname): (taxa aligned to, taxa seeded on and failed)} for these fragments."""
    sub = r[pd.MultiIndex.from_arrays([r["sample"], r["qname"]]).isin(keys)]
    out = {}
    for s, q, t, zf in zip(sub["sample"], sub["qname"], sub["taxon"].fillna("").astype(str), sub["zf"]):
        a, z = out.setdefault((s, q), (set(), set()))
        if t:
            a.add(t)
        if isinstance(zf, str) and zf not in ("", "*"):
            z.update(x.split(":")[0] for x in zf.split(","))
    return out


def own_fates(f, r, fn):
    """Every own fragment of each FN taxon T (fn: the FN rows, taxid as str; the v15 SAMs hold them all, their counts
    are the taxa tables'), with where it went: its fate."""
    own = f[f["xe"].fillna("").str.contains("source:")].copy()
    own["T"] = own["xe"].str.extract(r"source:(\d+)")[0]
    own = own.merge(fn[["sample", "taxid"]].rename(columns={"taxid": "T"}), on=["sample", "T"])
    seen = touched(r, pd.MultiIndex.from_arrays([own["sample"], own["qname"]]))
    own["on_T_any"] = [own_t in seen.get((s, q), (set(), set()))[0]
                       for s, q, own_t in zip(own["sample"], own["qname"], own["T"])]
    own["failed_on_T"] = [own_t in seen.get((s, q), (set(), set()))[1]
                          for s, q, own_t in zip(own["sample"], own["qname"], own["T"])]
    own["fate"] = np.select(
        [own["taxon"] == own["T"],
         (own["taxon"] != "") & own["failed_on_T"], (own["taxon"] != "") & own["on_T_any"], own["taxon"] != "",
         own["failed_on_T"]],
        ["best on the taxon", "elsewhere, failed on the taxon", "elsewhere, also on the taxon",
         "elsewhere, never seeded on the taxon", "nowhere, failed on the taxon"],
        "nowhere, seeded on other taxa only")
    return own


def fn_side(rt, f, r, taxa, tax, build):
    fn = taxa[taxa["error"] == "FN"].copy()
    for c in ("read_pairs", "own_fragments", "own_best_on_taxon", "own_best_on_taxon_mapq4", "own_best_elsewhere",
              "own_unaligned"):
        fn[c] = pd.to_numeric(fn[c], errors="coerce").fillna(0)
    fn["taxid"] = fn["taxid"].astype(str)
    tot = fn[["own_fragments", "own_best_on_taxon", "own_best_on_taxon_mapq4", "own_best_elsewhere",
              "own_unaligned"]].sum()
    print(f"{rt}: {len(fn)} FN taxa; their own fragments with a record {int(tot['own_fragments'])}: best on the taxon "
          f"{int(tot['own_best_on_taxon'])} ({int(tot['own_best_on_taxon_mapq4'])} at MAPQ >= 4), best elsewhere "
          f"{int(tot['own_best_elsewhere'])}, aligned nowhere {int(tot['own_unaligned'])}")
    print(f"  FN taxa without any own fragment with a record: {(fn['own_fragments'] == 0).mean():.3f}; "
          f"without one on the taxon: {(fn['own_best_on_taxon'] == 0).mean():.3f}")

    own = own_fates(f, r, fn)
    print(f"  {len(own)} own fragments of {own.groupby(['sample', 'T']).ngroups} FN taxa; their fate:",
          own["fate"].value_counts().to_dict())
    on = int((own["fate"] == "best on the taxon").sum())
    failed_els = int((own["fate"] == "elsewhere, failed on the taxon").sum())
    failed_now = int((own["fate"] == "nowhere, failed on the taxon").sum())
    print(f"  of the own fragments that reached the aligner on the taxon: {on} aligned to it, {failed_els} failed there "
          f"and aligned elsewhere, {failed_now} failed there and aligned nowhere: alignment sensitivity on the taxon "
          f"{on / max(1, on + failed_els + failed_now):.3f}")
    k = fn.set_index(["sample", "taxid"])
    k.index.names = ["sample", "T"]
    went = own[own["fate"].str.startswith("elsewhere")]
    if len(went):
        print(f"  the own fragments best elsewhere ({len(went)}): relation of that taxon "
              f"{went['relation'].value_counts(normalize=True).round(3).to_dict()}; its role "
              f"{went['role'].value_counts(normalize=True).round(3).to_dict()}; identity there median "
              f"{went['identity'].median():.4f} (>= 0.98: {(went['identity'] >= 0.98).mean():.3f}); counted there "
              f"{went['counted'].mean():.3f}; source genomes {went['source_genome'].value_counts(normalize=True).round(3).to_dict()}")
        on_t = own[own["fate"] == "best on the taxon"]
        print(f"  for comparison, the own fragments best on the taxon ({len(on_t)}): identity median "
              f"{on_t['identity'].median():.4f}; source genomes "
              f"{on_t['source_genome'].value_counts(normalize=True).round(3).to_dict()}")
    fails = own[own["fate"].str.contains("failed")]
    if len(fails):
        print(f"  the failed own fragments: their identity where they aligned elsewhere (median) "
              f"{fails['identity'].median():.3f}; their source genomes {fails['source_genome'].value_counts().to_dict()}")
    reached = own[own["fate"].isin(["best on the taxon", "elsewhere, failed on the taxon",
                                    "nowhere, failed on the taxon"])]
    per = reached.groupby(["sample", "T"]).size().reindex(k.index).fillna(0)
    lost = 1 - k["own_best_on_taxon"] / per.replace(0, np.nan)
    print(f"  per FN taxon with any own fragment that touched it ({int(per.gt(0).sum())}): share lost to failed "
          f"alignments median {lost.median():.3f}, 90% {lost.quantile(.9):.3f}; taxa losing >= half: "
          f"{(lost >= 0.5).mean():.3f}")

    # Where the own fragments that went elsewhere landed (the taxa tables count them all, by taxon name)
    table = fragments.model_tables(build, rt, lambda c: c in ("meta_sample", "taxon_name", "truth"))
    truth = dict(zip(zip(table["meta_sample"], table["taxon_name"]), table["truth"])) if table is not None else {}
    fp = set(zip(taxa.loc[taxa["error"] == "FP", "sample"], taxa.loc[taxa["error"] == "FP", "taxon_name"]))
    rows = []
    for (s, t), text in zip(k.index, k["own_best_elsewhere_on"].fillna("").astype(str)):
        for item in filter(None, text.split(",")):
            name, _, n = item.rpartition(":")
            role = ("FP" if (s, name) in fp else {1: "present (TP)", 0: "absent, not called"}.get(truth.get((s, name)),
                                                                                                  "no row"))
            rows.append({"sample": s, "T": t, "role": role, "n": int(n),
                         "relation": tax.relation(t, tax.name_id.get(name, ""))})
    el = pd.DataFrame(rows)
    if len(el):
        print(f"  all own fragments best elsewhere ({el['n'].sum()}), by that taxon's role: "
              f"{(el.groupby('role')['n'].sum() / el['n'].sum()).round(3).to_dict()}; by its relation: "
              f"{(el.groupby('relation')['n'].sum() / el['n'].sum()).round(3).to_dict()}")
        share = k["own_best_elsewhere"] / (k["own_best_elsewhere"] + k["own_best_on_taxon"]).replace(0, np.nan)
        print(f"  per FN taxon, the share of its aligned own fragments that went elsewhere: median {share.median():.3f}, "
              f"mean {share.mean():.3f}; >= half: {(share >= 0.5).mean():.3f}; all: {(share == 1).mean():.3f}")
        # By the size of the taxon's genus in the database (its species not held out): protal aligns a read against
        # its top 3 anchors only (--align_top), and a crowded genus pushes the own species out of them
        genus_of = {t: tax.ancestor(t, "genus") for t in set(k.index.get_level_values("T"))}
        sizes = pd.Series([tax.ancestor(s, "genus") for s, r_ in tax.rank.items() if r_ == "species"
                           and tax.name.get(s, "") not in tax.held_out]).value_counts()
        k2 = pd.DataFrame({"share": share, "aligned": k["own_best_elsewhere"] + k["own_best_on_taxon"]})
        k2["genus_species"] = [sizes.get(genus_of.get(t), 0) for t in k2.index.get_level_values("T")]
        k2 = k2[k2["aligned"] >= 5]
        k2["bin"] = pd.cut(k2["genus_species"], [0, 1, 3, 10, 30, 100, 1e9], labels=["1", "2-3", "4-10", "11-30",
                                                                                  "31-100", ">100"])
        print(f"  FN taxa with >= 5 aligned own fragments ({len(k2)}), the share that went elsewhere by the species of "
              f"their genus in the database (Spearman {k2['share'].corr(k2['genus_species'], method='spearman'):.3f}):")
        print(k2.groupby("bin", observed=True)["share"].agg(taxa="size", median="median", mean="mean").round(3)
              .to_string())
        top = el.groupby(["sample", "T"])["n"].sum().sort_values(ascending=False)
        print(f"  the 10 FN taxa with most own fragments elsewhere hold {top.head(10).sum() / el['n'].sum():.3f} of "
              f"them: {[(tax.name.get(t, t), int(n)) for (s, t), n in top.head(5).items()]}")

    # Before seeding (paired-end): own fragments with a record per read pair, by the species' genome in the sample
    if (k["read_pairs"] > 0).any():
        kk = k[k["read_pairs"] >= 20].copy()
        kk["genome"] = ["several" if "," in str(g) else "insilico" if str(g).startswith("insilico_") else
                        "rep" if g == tax.rep.get(t) else "strain"
                        for g, t in zip(kk["genomes"], kk.index.get_level_values("T"))]
        idn = f[(f["role"] == "FN") & (f["relation"] == "own")].groupby(["sample", "taxon"])["identity"].mean()
        idn.index.names = ["sample", "T"]
        kk["own_identity"] = idn.reindex(kk.index)
        kk["on_per_pair"] = kk["own_best_on_taxon"] / kk["read_pairs"]
        kk["touched_per_pair"] = per.reindex(kk.index) / kk["read_pairs"]
        out = kk.groupby("genome").agg(taxa=("on_per_pair", "size"), read_pairs=("read_pairs", "median"),
                                       on_taxon_per_1000_pairs=("on_per_pair", lambda x: 1000 * x.median()),
                                       touched_per_1000_pairs=("touched_per_pair", lambda x: 1000 * x.median()),
                                       own_identity=("own_identity", "median"))
        print("  own fragments per 1000 read pairs (FN taxa with >= 20 pairs), by the species' genome in the sample:")
        print(out.round(3).to_string())
        both = kk.dropna(subset=["own_identity"])
        both = both[both["genome"] != "several"]
        if len(both) > 20:
            rho = both[["own_identity", "on_per_pair"]].corr(method="spearman").iloc[0, 1]
            print(f"  Spearman of own identity with fragments on the taxon per read pair ({len(both)} taxa): {rho:.3f}")
            both["bin"] = pd.cut(both["own_identity"], [0, .95, .97, .98, .99, .995, 1.001])
            print(both.groupby("bin", observed=True)["on_per_pair"].agg(
                taxa="size", per_1000_pairs=lambda x: 1000 * x.median()).round(2).to_string())


def fp_side(rt, f, r, taxa):
    w = f[f["counted"] & (f["role"] == "FP")].copy()
    w["kind"] = np.select([w["source"].eq("host"), w["source_held_out"], w["source_taxid"].eq("")],
                          ["host", "held out (not in the database)", "unknown"], "in the database")
    print(f"{rt}: {len(w)} counted FP fragments by their source: {w['kind'].value_counts(normalize=True).round(3).to_dict()}")
    db = w[w["kind"] == "in the database"].copy()
    if len(db):
        seen = touched(r, pd.MultiIndex.from_arrays([db["sample"], db["qname"]]))
        db["also_on_source"] = [sid in seen.get((s, q), (set(), set()))[0]
                                for s, q, sid in zip(db["sample"], db["qname"], db["source_taxid"])]
        db["state"] = np.select([db["zf_source"], db["also_on_source"] | db["za_source"]],
                                ["seeded on its species, failed there", "aligned to its species too (worse or equal)"],
                                "never seeded on its species")
        print(f"  of those from a species in the database ({len(db)}): "
              f"{db['state'].value_counts(normalize=True).round(3).to_dict()}")
        print(f"  their source genomes: {db['source_genome'].value_counts(normalize=True).round(3).to_dict()}; "
              f"relation to the FP taxon: {db['relation'].value_counts(normalize=True).round(3).to_dict()}")
        print(f"  identity on the FP taxon by state: {db.groupby('state')['identity'].median().round(4).to_dict()}")
        fn = set(zip(taxa.loc[taxa["error"] == "FN", "sample"], taxa.loc[taxa["error"] == "FN", "taxid"].astype(str)))
        db["source_is_fn"] = [(s, t) in fn for s, t in zip(db["sample"], db["source_taxid"])]
        print(f"  their species is an FN taxon of the sample: {db['source_is_fn'].mean():.3f}")
    # Per FP taxon: the main kind of its fragments, and the share whose species is in the database and failed there
    per = w.groupby(["sample", "taxon"])["kind"].agg(lambda s: s.value_counts().index[0])
    print(f"  FP taxa by their fragments' main source: {per.value_counts(normalize=True).round(3).to_dict()}")
    if len(db):
        fail = db.groupby(["sample", "taxon"])["state"].agg(lambda s: (s == "seeded on its species, failed there").mean())
        n = w.groupby(["sample", "taxon"]).size()
        frac = (fail.reindex(n.index).fillna(0) * db.groupby(["sample", "taxon"]).size().reindex(n.index).fillna(0) / n)
        print(f"  FP taxa where most fragments seeded on their own (database) species and failed there: "
              f"{(frac > 0.5).mean():.3f} ({int((frac > 0.5).sum())} of {len(n)})")


def main(argv=None):
    opts = parse_args(argv)
    pd.set_option("display.width", 250)
    tax = Taxonomy(opts.build)
    for rt in opts.types.split(","):
        f = pd.read_pickle(os.path.join(opts.records, f"{rt}.fragments.pkl.gz"))
        taxa = fragments.load(opts.records, rt, "taxa")
        r = fragments.load(opts.records, rt, "records")[["sample", "qname", "taxon", "zf"]]
        print(f"\n# {rt}\n\n## FN taxa: their own reads\n", flush=True)
        fn_side(rt, f, r, taxa, tax, opts.build)
        print(f"\n## FP taxa: their reads' species\n", flush=True)
        fp_side(rt, f, r, taxa)
    return 0


if __name__ == "__main__":
    sys.exit(main())
