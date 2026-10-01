#!/usr/bin/env python3
"""own_cluster.py [READ_TYPES...] - a prototype of features from a taxon's own read cluster, and their effect on F1.

A present species beside a congener the database lacks gets that congener's reads too: the features averaged over
all its reads (identity, MAPQ, congener fit) then look like those of an absent taxon living on borrowed reads, and
the model misses it (test_errors.py: present taxa with a missing congener in the sample are missed 2-3 times as
often). The reads within 0.04 of the taxon's best (its 98th percentile identity, as protal's MSA margin) are its
own cluster; features of that cluster:
  own_log_fragments       log10(1 + the cluster's records)
  own_share               the cluster's share of the taxon's records
  own_low_mapq_share      the cluster's records with MAPQ < 10 (protal's MAPQ: the score gap to the next candidate)
  own_congener_fit_share  the cluster's records with an alternative on a congener within 1 edit (the ZA tag)
at margins 0.01, 0.02 and 0.04 below the 98th percentile, and the identity distribution's shape: identity_gap (98th
percentile less the median) and identity_99_share (records at 0.99 or more).
From the training and test SAMs of 0.7.1's pipeline ($V071/{training,test}/points/*/protal*/alignments): primary
records, identity M/(M+X+I+D) as protal's AlignmentIdentity. Also, per taxon, the share of its records whose read
was simulated from another species (the read name's genome): how borrowed the reads are, for the report.
Then the forest refitted (exp_lib) with the base features and with these added: CV and test F1 as features_exp.py,
and the test set's misses by whether a missing congener is in the sample. Writes results/own_cluster.md.
"""
import collections
import csv
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
V071 = os.path.expanduser(os.environ.get("V071", "~/bench071/V071"))
os.environ["V2"] = V071
os.environ.setdefault("FPEXP", os.path.expanduser("~/bench071/f1exp"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "2026-10-01-alignment-features"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "..", "scripts"))
import exp_lib as X  # noqa: E402
from model_features import NORMALIZED_FEATURES as BASE  # noqa: E402

READ_TYPES = sys.argv[1:] or ["pe", "se", "pb", "ont"]
SEEDS = [int(s) for s in os.environ.get("SEEDS", "1 2 3").split()]
MARGINS = (0.01, 0.02, 0.04)
OWN = {m: [f"own_log_fragments_{m}", f"own_share_{m}", f"own_low_mapq_share_{m}", f"own_congener_fit_share_{m}"] for m in MARGINS}
SHAPE = ["identity_gap", "identity_99_share"]
CACHE = os.path.join(os.environ["FPEXP"], "own2")


def species_of_genome():
    out = {}
    with open(os.path.join(V071, "genomes.tsv")) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[0]] = f[1].split(";s__")[-1]
    # the world's genomes outside the release (species no database has) too
    world = os.path.expanduser("~/bench071/world/full/simulation/genomes.tsv")
    if os.path.exists(world):
        with open(world) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                out.setdefault(r["accession"], r["gtdb_taxonomy"].split(";s__")[-1])
    return out


def sam_of(split, rt, sample):
    sub = "protal_se" if rt == "se" else "protal"
    point = sample.rsplit("_s_", 1)[0]
    if rt == "se":
        point = point.removesuffix("_se")
    base = os.path.join(V071, split, "points", point, sub, "alignments", sample)
    return next((base + ext for ext in (".sam.zst", ".sam.gz", ".sam") if os.path.exists(base + ext)), base + ".sam.zst")


def features(split, rt):
    path = os.path.join(CACHE, f"{split}_{rt}.tsv")
    if os.path.exists(path):
        return pd.read_csv(path, sep="\t")
    genus = X.db_genera()
    taxid_name = {}
    with open(os.path.join(V071, "training_db", "genome2tiid.tsv")) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4:
                taxid_name[f[1]] = f[3].split(";s__")[-1]
    sp_of = species_of_genome()
    table = X.load(split, rt)
    rows = []
    for sample in sorted(table["meta_sample"].unique()):
        sam = sam_of(split, rt, sample)
        if not os.path.exists(sam):
            print(f"no SAM for {sample}: {sam}", file=sys.stderr)
            continue
        per = collections.defaultdict(list)  # taxid -> [(identity, mapq, congener_fit, borrowed)]
        for f in X.sam_records(sam):
            flag = int(f[1])
            if flag & (4 | 256 | 2048):
                continue
            taxid = f[2].split("_")[0]
            c = X.cigar_counts(f[5])
            aligned = c["M"] + c["X"] + c["I"] + c["D"]
            if aligned == 0:
                continue
            ident = c["M"] / aligned
            za = re.search(r"\tZA:Z:([^\t\n]*)", "\t" + f[11]) if len(f) > 11 else None
            fit = False
            if za and za.group(1) != "*":
                g = genus.get(taxid)
                for alt in za.group(1).split(","):
                    t, _, extra = alt.partition(":")
                    if t != taxid and genus.get(t) == g and extra.isdigit() and int(extra) <= 1:
                        fit = True
                        break
            source = re.sub(r"_contig.*$", "", f[0]).split("/")[0]
            src = sp_of.get(source)
            borrowed = src is not None and src != taxid_name.get(taxid)
            per[taxid].append((ident, int(f[4]), fit, borrowed))
        for taxid, recs in per.items():
            ids = np.array([r[0] for r in recs])
            top = np.quantile(ids, 0.98)
            row = {"meta_sample": sample, "taxon": int(taxid), "borrowed_share": sum(r[3] for r in recs) / len(recs),
                   "identity_gap": top - float(np.median(ids)), "identity_99_share": float((ids >= 0.99).mean())}
            for m in MARGINS:
                own = [r for r in recs if r[0] >= top - m]
                n = max(1, len(own))
                row.update({f"own_log_fragments_{m}": np.log10(1 + len(own)), f"own_share_{m}": len(own) / len(recs),
                            f"own_low_mapq_share_{m}": sum(r[1] < 10 for r in own) / n,
                            f"own_congener_fit_share_{m}": sum(r[2] for r in own) / n,
                            f"own_borrowed_share_{m}": sum(r[3] for r in own) / n})
            rows.append(row)
    df = pd.DataFrame(rows)
    os.makedirs(CACHE, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    return df


def with_own(split, rt):
    t = X.load(split, rt)
    f = features(split, rt)
    t = t.merge(f, on=["meta_sample", "taxon"], how="left")
    for c in [c for cols in OWN.values() for c in cols] + SHAPE:
        t[c] = t[c].fillna(0.0)
    return t


def main():
    lines = ["# Own read cluster features (own_cluster.py)", "",
             f"Seeds {SEEDS}; CV = species held out (5 folds) at 0.5; test = the independent test set; delta = test F1 "
             "against the base with a paired bootstrap (mean, 95% interval).", ""]
    for rt in READ_TYPES:
        train, test = with_own("training", rt), with_own("test", rt)
        lines += [f"## {rt}", "", "| features | CV F1 | test F1 | test FP | test FN | FN with a missing congener in the "
                  "sample | delta (interval) |", "|---|---|---|---|---|---|---|"]
        base_calls = {}
        variants = [("base", BASE)] + [(f"base + own cluster within {m}", BASE + OWN[m]) for m in MARGINS] + \
                   [("base + identity shape", BASE + SHAPE), ("base + own within 0.01 + shape", BASE + OWN[0.01] + SHAPE)]
        for label, feats in variants:
            cv, tf, fp, fn, fnc, diffs = [], [], [], [], [], []
            for seed in SEEDS:
                p_cv = X.cv_predict(train, feats, seed)
                p = X.test_predict(train, test, feats, seed)
                call = p >= X.KNOB
                y = test["truth"].to_numpy()
                cv.append(X.scores(train["truth"].to_numpy(), p_cv >= X.KNOB)["F1"])
                s = X.scores(y, call)
                tf.append(s["F1"]); fp.append(s["FP"]); fn.append(s["FN"])
                fnc.append(int(((y == 1) & ~call & (test["meta_novel_congener"].to_numpy() == 1)).sum()))
                if label == "base":
                    base_calls[seed] = call
                else:
                    diffs.append(X.bootstrap_diff(test, base_calls[seed], call, n=500, seed=seed))
            d = np.array(diffs).mean(axis=0) if diffs else (0, 0, 0)
            lines.append(f"| {label} | {np.mean(cv):.4f} | {np.mean(tf):.4f} | {np.mean(fp):.1f} | {np.mean(fn):.1f} | "
                         f"{np.mean(fnc):.1f} | {d[0]:+.4f} ({d[1]:+.4f}, {d[2]:+.4f}) |")
            print(lines[-1], flush=True)
        pres = test[test["truth"] == 1]
        lines += ["", f"Borrowed reads (simulated from another species) among the present taxa's records, test set: "
                  f"median {pres['borrowed_share'].median():.3f}; with a missing congener in the sample "
                  f"{pres.loc[pres['meta_novel_congener'] == 1, 'borrowed_share'].median():.3f} of all records and "
                  f"{pres.loc[pres['meta_novel_congener'] == 1, 'own_borrowed_share_0.01'].median():.3f} of the own cluster's (0.01); "
                  f"absent taxa: {test.loc[test['truth'] == 0, 'borrowed_share'].median():.3f}.", ""]
    out = os.path.join(HERE, "..", "results", "own_cluster.md" if len(READ_TYPES) == 4 else f"own_cluster_{'_'.join(READ_TYPES)}.md")
    with open(out, "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
