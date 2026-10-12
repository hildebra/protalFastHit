"""Each false positive and false negative of an r226 build (training samples = species held out, and the test set;
the knob as protal calls) given one cause from its reads (model_logs/error_reads/<type>/taxa.tsv.gz) and the calls'
metadata, the first that applies in this order, and the F1 each cause costs (F1 with its errors removed - F1):

False positives, by the main source of the reads on the taxon (fragments_on_taxon_from):
  twin of a held-out species: a species the database lacks, same genus, within 0.02 (meta_novel_distance)
  cloud of a held-out congener: a species the database lacks, same genus, 3 or more fragments
  a held-out congener's stray reads: a species the database lacks, same genus, 1-2 fragments
  a held-out species of another genus
  a present congener's reads (a database species of the genus)
  another genus's reads (a database species of another genus)
False negatives:
  1-2 fragments: at most 2 fragments before and after the MAPQ filter (depth)
  MAPQ filter: at most 2 fragments after it, 3 or more before
  swamped by a held-out congener: most reads on the taxon from a species the database lacks
  swamped by another species: most reads on the taxon from another species in the database
  own reads elsewhere: half or more of its own aligned reads best on another species
  divergent: identity below 0.97
  other

usage: causes.py <build dir> <read type> [<read type> ...]
"""
import re
import sys
import numpy as np
import pandas as pd

D, types = sys.argv[1], sys.argv[2:]
pd.set_option("display.width", 250)


def genus(name):
    m = re.match(r"s__(\S+)", str(name))
    return m.group(1) if m else ""


def parts(s):
    out = []
    for p in str(s).split(","):
        if ":" in p:
            name, n = p.rsplit(":", 1)
            out.append((name.replace(" (not in the database)", ""), int(n), not name.endswith("(not in the database)")))
    return out


def fp_cause(r):
    p = parts(r["fragments_on_taxon_from"])
    if not p:
        return "no source traced"
    name, _, indb = max(p, key=lambda q: q[1])
    same = genus(name) == genus(r["taxon_name"])
    frags = pd.to_numeric(r["fragments_on_taxon"], errors="coerce")
    if not indb and same:
        if r["meta_novel_distance"] < 0.02:
            return "twin of a held-out species (<0.02)"
        return "cloud of a held-out congener (>=3 fragments)" if frags >= 3 else "a held-out congener's stray reads (1-2)"
    if not indb:
        return "a held-out species of another genus"
    return "a present congener's reads" if same else "another genus's reads"


def fn_cause(r):
    after, before = r["fragments"], r["fragments_all"]
    if after <= 2 and before <= 2:
        return "1-2 fragments"
    if after <= 2:
        return "MAPQ filter (<=2 after, >=3 before)"
    p = parts(r["fragments_on_taxon_from"])
    tot = sum(n for _, n, _ in p)
    own = sum(n for name, n, _ in p if name == r["taxon_name"])
    if tot and own / tot < 0.5:
        others = [(name, n, indb) for name, n, indb in p if name != r["taxon_name"]]
        name, _, indb = max(others, key=lambda q: q[1])
        return "swamped by a held-out congener" if not indb else "swamped by another database species"
    on, el = r["own_best_on_taxon"], r["own_best_elsewhere"]
    if on + el >= 3 and el / (on + el) >= 0.5:
        return "own reads best elsewhere"
    if r["identity"] < 0.97:
        return "divergent (identity < 0.97)"
    return "other"


for T in types:
    SUF = "" if T == "pe" else "_" + T
    calls = pd.read_csv(f"{D}/model_logs/trained_model{SUF}.calls.tsv.gz", sep="\t",
                        usecols=["meta_sample", "meta_scenario", "meta_rep_genome", "meta_insilico_strain",
                                 "meta_novel_distance", "taxon_name", "truth", "set", "call"], dtype={"meta_scenario": str})
    tabs = []
    for s in ("training", "test"):
        t = pd.read_csv(f"{D}/work/{s}/training_data{SUF}.tsv", sep="\t",
                        usecols=["meta_sample", "taxon_name", "fragments", "fragments_all", "identity",
                                 "allele_copy_share"]).drop_duplicates(["meta_sample", "taxon_name"])
        tabs.append(calls[calls["set"] == s].merge(t, on=["meta_sample", "taxon_name"], how="left"))
    c = pd.concat(tabs, ignore_index=True)
    c["meta_novel_distance"] = pd.to_numeric(c["meta_novel_distance"], errors="coerce").fillna(1.0)
    rep = pd.to_numeric(c["meta_rep_genome"], errors="coerce") == 1
    ins = pd.to_numeric(c["meta_insilico_strain"], errors="coerce") == 1
    alle = c["allele_copy_share"].fillna(0) > 0.2
    c["genome"] = np.where(rep, "representative", np.where(ins, "in-silico strain",
                                                           np.where(alle, "real strain, alleles", "real strain, none")))
    err = pd.read_csv(f"{D}/model_logs/error_reads/{T}/taxa.tsv.gz", sep="\t", dtype={"scenario": str})
    err = err[err["error"].isin(["FP", "FN"])].merge(
        c[["set", "meta_sample", "taxon_name", "meta_novel_distance", "genome", "fragments", "fragments_all", "identity"]],
        left_on=["set", "sample", "taxon_name"], right_on=["set", "meta_sample", "taxon_name"], how="left")
    for col in ("own_best_on_taxon", "own_best_elsewhere", "fragments", "fragments_all", "identity"):
        err[col] = pd.to_numeric(err[col], errors="coerce").fillna(0 if col != "identity" else 1)
    err["cause"] = [fp_cause(r) if r["error"] == "FP" else fn_cause(r) for _, r in err.iterrows()]
    err["kind"] = err["scenario"].fillna("").replace("", "design")
    print(f"######## {T}")
    for s in ("training", "test"):
        x = c[c["set"] == s]
        tp = int(((x["truth"] == 1) & (x["call"] == 1)).sum())
        fpn = int(((x["truth"] == 0) & (x["call"] == 1)).sum())
        fnn = int(((x["truth"] == 1) & (x["call"] == 0)).sum())
        f1 = 2 * tp / (2 * tp + fpn + fnn)
        e = err[err["set"] == s]
        print(f"== {s} ({'species held out' if s == 'training' else 'final model'}): F1 {f1:.4f}, FP {fpn}, FN {fnn} "
              f"(error_reads: {int((e['error'] == 'FP').sum())} FP, {int((e['error'] == 'FN').sum())} FN)")
        rows = []
        for (er, cause), g in e.groupby(["error", "cause"]):
            n = len(g)
            if er == "FP":
                gain = 2 * tp / (2 * tp + fpn - n + fnn) - f1
            else:
                gain = 2 * (tp + n) / (2 * (tp + n) + fpn + fnn - n) - f1
            rows.append({"error": er, "cause": cause, "n": n, "share": n / (e["error"] == er).sum(),
                         "F1 if removed": f"+{gain:.4f}",
                         "design": int((g["kind"] == "design").sum()), "scenarios": int((g["kind"] != "design").sum())})
        print(pd.DataFrame(rows).to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        if s == "training":
            f = e[e["error"] == "FN"]
            print("   FN causes by genome simulated (training):")
            print(pd.crosstab(f["genome"], f["cause"], margins=True).to_string())
