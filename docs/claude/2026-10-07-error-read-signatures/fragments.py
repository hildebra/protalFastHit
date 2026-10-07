#!/usr/bin/env python3
"""Fragments (paired-end: both mates; long reads: the whole read) of the FP and FN error records, with their source
and features, from parse_error_records.py's tables.

Each fragment: its best record (protal's: the first without flag 0x904) and where it counts (MAPQ >= 4); the read's
source species and how it relates to the taxon it counts for (own, genus, family, order, other; held out of the
training database or not; representative genome, another genome of the species or an in-silico strain); the role of
that taxon in the sample (FP, FN, or another taxon, TP or TN); and features of its records on that taxon:
identity, mismatches, gap opens, clips, gene ends, ZU/ZT, ZA (alternatives, their least edits more, negative ones:
an alternative that fits better than the alignment kept), ZF (failed seeds), mates, ZR, the gene's conservation,
GC and low complexity. Writes OUT/<type>.fragments.pkl.gz.

    python3 fragments.py --records ~/v15/err_analysis --build local/v15
"""
import argparse
import csv
import os
import sys

import numpy as np
import pandas as pd

MIN_MAPQ = 4
TABLES = {"pe": "training_data.tsv", "se": "training_data_se.tsv", "pb": "training_data_pb.tsv",
          "ont": "training_data_ont.tsv"}


def load(records, rt, what):
    """parse_error_records.py's table of one read type (records or taxa), its sample as <set>:<sample>: a scenario's
    hold-in samples (set training) and hold-out samples (set test) share their names."""
    df = pd.read_pickle(os.path.join(records, f"{rt}.{what}.pkl.gz"))
    df["sample"] = df["set"].astype(str) + ":" + df["sample"].astype(str)
    return df


def model_tables(build, rt, usecols):
    """The model's training and test tables of one read type (None without them), meta_sample as <set>:<sample>."""
    frames = []
    for which in ("training", "test"):
        path = os.path.join(build, which, TABLES[rt])
        if os.path.isfile(path):
            t = pd.read_csv(path, sep="\t", usecols=usecols)
            t["meta_sample"] = which + ":" + t["meta_sample"].astype(str)
            frames.append(t)
    return pd.concat(frames, ignore_index=True) if frames else None


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--records", required=True, help="parse_error_records.py's output folder (also the output)")
    p.add_argument("--build", required=True, help="the build's local folder (internal_taxonomy.dmp, "
                                                  "heldout_species.txt, model_logs/gene_congeners.tsv)")
    p.add_argument("--types", default="pe,se,pb,ont")
    return p.parse_args(argv)


class Taxonomy:
    def __init__(self, build):
        self.parent, self.rank, self.rep, self.name_id, self.name = {}, {}, {}, {}, {}
        with open(os.path.join(build, "internal_taxonomy.dmp"), newline="") as fh:
            reader = csv.reader(fh, delimiter="\t")
            next(reader)
            for row in reader:
                self.parent[row[0]], self.rank[row[0]], self.name[row[0]] = row[1], row[4], row[3]
                self.rep[row[0]] = row[6] if len(row) > 6 else ""
                if row[4] == "species":
                    self.name_id.setdefault(row[3], row[0])
        with open(os.path.join(build, "heldout_species.txt")) as fh:
            self.held_out = {line.split("\t")[0].strip() for line in fh if line.strip()}
        self.cache = {}

    def ancestor(self, taxid, rank):
        for _ in range(64):
            if taxid not in self.rank:
                return None
            if self.rank[taxid] == rank:
                return taxid
            if self.parent[taxid] == taxid:
                return None
            taxid = self.parent[taxid]
        return None

    def relation(self, sid, taxid):
        key = (sid, taxid)
        r = self.cache.get(key)
        if r is None:
            if not sid or taxid not in self.rank:
                r = "unknown"
            elif sid == taxid:
                r = "own"
            else:
                r = "other"
                for rank in ("genus", "family", "order"):
                    a = self.ancestor(sid, rank)
                    if a is not None and a == self.ancestor(taxid, rank):
                        r = rank
                        break
            self.cache[key] = r
        return r


def build(r, taxa, tax, genes):
    """The fragments of one read type's records (vectorised over the records; ZA and ZF parsed per record)."""
    key = list(zip(taxa["sample"], taxa["taxid"].astype(str)))
    fp = {k for k, e in zip(key, taxa["error"]) if e == "FP"}
    fn = {k for k, e in zip(key, taxa["error"]) if e == "FN"}
    p_of = dict(zip(key, taxa["p"]))
    r = r.reset_index(drop=True)
    r["taxon"] = r["taxon"].fillna("").astype(str)
    r["frag"] = r.groupby(["sample", "qname"], sort=False).ngroup()
    aligned = r["taxon"] != ""
    prim = aligned & ((r["flag"] & 0x904) == 0)
    best = r[prim].groupby("frag").head(1)
    rest = r[aligned & ~r["frag"].isin(best["frag"])].groupby("frag").head(1)
    best = pd.concat([best, rest]).set_index("frag")

    f = r.groupby("frag").head(1).set_index("frag")[["type", "set", "scenario", "point", "sample", "qname", "xe",
                                                      "xs", "xg"]].copy()
    f["source"] = f["xs"].fillna("").str.replace(" (not in the database)", "", regex=False)
    f["source_taxid"] = [tax.name_id.get(s, "") if s and s != "host" else "" for s in f["source"]]
    f["source_held_out"] = f["source"].isin(tax.held_out)
    f["source_genome"] = ["insilico" if str(g).startswith("insilico_") else "host" if s == "host" else
                          "rep" if sid and g == tax.rep.get(sid) else "strain"
                          for g, s, sid in zip(f["xg"], f["source"], f["source_taxid"])]
    f["records"] = r.groupby("frag").size()
    f["aligned_records"] = aligned.groupby(r["frag"]).sum()
    f["taxa"] = r[aligned].groupby("frag")["taxon"].nunique().reindex(f.index, fill_value=0)

    # ZF: the taxa the fragment seeded on but did not align to (any record)
    zf_sets = {}
    for frag, v in zip(r["frag"], r["zf"]):
        if isinstance(v, str) and v:
            zf_sets.setdefault(frag, set()).update(x.split(":")[0] for x in v.split(","))
    f["zf_n"] = [len(zf_sets.get(i, ())) for i in f.index]
    f["zf_source"] = [bool(sid) and sid in zf_sets.get(i, ()) for i, sid in zip(f.index, f["source_taxid"])]

    f["taxon"] = best["taxon"].reindex(f.index).fillna("")
    f["counted"] = best["mapq"].reindex(f.index).fillna(-1) >= MIN_MAPQ
    f["mapq"] = best["mapq"].reindex(f.index)
    f["gene"] = best["gene"].reindex(f.index)
    f["proper"] = (best["flag"].reindex(f.index).fillna(0).astype(int) & 2) > 0
    bflag = best["flag"].reindex(f.index).fillna(0).astype(int)
    f["mate_unmapped"] = ((bflag & 1) > 0) & ((bflag & 8) > 0)
    f["abs_tlen"] = best["tlen"].reindex(f.index).abs()
    f["read_len"] = (best["clip_left"] + best["clip_right"] + best["matches"] + best["mismatches"]
                     + best["ins"]).reindex(f.index)
    f["role"] = ["unaligned" if not t else "FP" if (s, t) in fp else "FN" if (s, t) in fn else "other"
                 for s, t in zip(f["sample"], f["taxon"])]
    f["p"] = [p_of.get((s, t), np.nan) for s, t in zip(f["sample"], f["taxon"])]
    f["relation"] = [tax.relation(sid, t) if t else "unaligned" for sid, t in zip(f["source_taxid"], f["taxon"])]
    f["zf_self"] = [t in zf_sets.get(i, ()) for i, t in zip(f.index, f["taxon"])]
    f["zf_fn"] = [any((s, x) in fn for x in zf_sets.get(i, ())) for i, s in zip(f.index, f["sample"])]
    f["zf_fp"] = [any((s, x) in fp for x in zf_sets.get(i, ())) for i, s in zip(f.index, f["sample"])]

    # The records on the best record's taxon (not secondary; all of them if only secondary ones)
    on = aligned & (r["taxon"] == r["frag"].map(f["taxon"]))
    use = on & ((r["flag"] & 0x900) == 0)
    missing = ~r["frag"].isin(r.loc[use, "frag"])
    use = use | (on & missing)
    u = r[use].copy()
    u["n"] = u["matches"] + u["mismatches"] + u["ins"] + u["dels"]
    u["ref_end"] = u["pos"] + u["matches"] + u["mismatches"] + u["dels"] - 1
    u["clipped"] = (u["clip_left"] >= 10) | (u["clip_right"] >= 10)
    u["at_edge"] = (u["pos"] == 1) | (u["ref_end"] >= u["gene_len"])
    u["rel_pos"] = (u["pos"] + u["ref_end"]) / 2 / u["gene_len"].clip(lower=1)
    u["aligned"] = u["matches"] + u["mismatches"] + u["ins"]
    g = u.groupby("frag")
    agg = pd.DataFrame({
        "matches": g["matches"].sum(), "n": g["n"].sum(), "mismatches": g["mismatches"].sum(),
        "indels": g["ins"].sum() + g["dels"].sum(), "gap_opens": g["gap_opens"].sum(), "aligned": g["aligned"].sum(),
        "longest_match": g["longest_match"].max(), "clipped": g["clipped"].any(),
        "clip_bases": g["clip_left"].sum() + g["clip_right"].sum(), "at_edge": g["at_edge"].any(),
        "rel_pos": g["rel_pos"].mean(), "zu": g["zu"].sum(), "zt": g["zt"].sum(), "zr": g["zr"].max(),
        "mates": g.size(), "two_genes": g["gene"].nunique() > 1, "gc": g["gc"].mean(),
        "top_trinucleotide": g["top_trinucleotide"].max(), "homopolymer": g["homopolymer"].max()})
    agg["identity"] = agg["matches"] / agg["n"].clip(lower=1)
    f = f.join(agg.drop(columns=["matches", "n"]))
    f["mates_split"] = (f["mates"] == 1) & ((bflag & 1) > 0) & ~f["mate_unmapped"]

    # ZA of those records: alternatives and their edits more (negative: an alternative fits better)
    za = {}
    for frag, v in zip(u["frag"], u["za"]):
        if isinstance(v, str) and v not in ("", "*"):
            for alt in v.split(","):
                t, _, more = alt.partition(":")
                try:
                    za.setdefault(frag, []).append((t, int(more)))
                except ValueError:
                    pass
    alts = [za.get(i, []) for i in f.index]
    f["za_n"] = [len({t for t, _ in a}) for a in alts]
    f["za_min"] = [min((m for _, m in a), default=np.nan) for a in alts]
    f["za_negative"] = [any(m < 0 for _, m in a) for a in alts]
    f["za_tie"] = [any(m == 0 for _, m in a) for a in alts]
    f["za_source"] = [bool(sid) and any(t == sid for t, _ in a) for a, sid in zip(alts, f["source_taxid"])]
    f["za_fn"] = [any((s, t) in fn for t, _ in a) for a, s in zip(alts, f["sample"])]
    f["za_fp"] = [any((s, t) in fp for t, _ in a) for a, s in zip(alts, f["sample"])]

    # The read's other records: secondary, supplementary, on other taxa; its share of aligned bases on the taxon
    sec = ((r["flag"] & 0x100) > 0) & aligned
    sup = ((r["flag"] & 0x800) > 0) & aligned
    f["secondary"] = sec.groupby(r["frag"]).sum()
    f["supplementary"] = sup.groupby(r["frag"]).sum()
    f["other_taxa_records"] = (aligned & ~on).groupby(r["frag"]).sum()
    bases = (r["matches"].fillna(0) + r["mismatches"].fillna(0) + r["ins"].fillna(0)) * aligned
    f["share_on_taxon"] = (bases * on).groupby(r["frag"]).sum() / bases.groupby(r["frag"]).sum().clip(lower=1)

    gf = genes
    for col, name in (("within_factor", "gene_within"), ("between_factor", "gene_between"),
                      ("identical_share", "gene_identical"), ("near_identical_share", "gene_near_identical")):
        f[name] = f["gene"].map(gf[col])
    return f.reset_index(drop=True)


def main(argv=None):
    opts = parse_args(argv)
    tax = Taxonomy(opts.build)
    genes = pd.read_csv(os.path.join(opts.build, "model_logs", "gene_congeners.tsv"), sep="\t").set_index("geneid")
    genes.index = genes.index.astype(float)
    for rt in opts.types.split(","):
        r = load(opts.records, rt, "records")
        taxa = load(opts.records, rt, "taxa")
        f = build(r, taxa, tax, genes)
        f.to_pickle(os.path.join(opts.records, f"{rt}.fragments.pkl.gz"))
        print(f"{rt}: {len(f)} fragments, {int(f['counted'].sum())} counted; roles of the counted: "
              f"{f.loc[f['counted'], 'role'].value_counts().to_dict()}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
