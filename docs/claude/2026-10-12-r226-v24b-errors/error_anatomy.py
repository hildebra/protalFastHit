"""The false positives and false negatives of an r226 build's models (calls at the knob protal calls at: the calls
table's `call`), with species held out (training samples) and on the test set (final model):

1. where they are: by set and sample kind (design or scenario), with each kind's share of the errors;
2. the false positives' reads (model_logs/error_reads/<type>/taxa.tsv.gz, fragments_on_taxon_from): from a species the
   database lacks or has, of the taxon's genus or another, how many fragments, how far the held-out congener is,
   how sure the model was, and the taxa that recur;
3. the false negatives: by the genome simulated and fragments, and where their own reads went (own_best_on_taxon,
   ..._mapq4, own_best_elsewhere(_on), own_unaligned, own_seeded_not_aligned), against the true positives' features.

usage: error_anatomy.py <build dir> <read type>
"""
import re
import sys
import numpy as np
import pandas as pd

D, T = sys.argv[1], sys.argv[2]
SUF = "" if T == "pe" else "_" + T
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 40)

FEAT = ["fragments", "fragments_all", "em_fragments", "identity", "top_identity", "low_mapq_share", "congener_fit_share",
        "other_genus_fit_share", "allele_copy_share", "allele_explained_share", "hit_gene_fraction",
        "excess_scaled_median", "sample_log_fragments", "db_nearest_congener"]


def genus(name):
    m = re.match(r"s__(\S+)", str(name))
    return m.group(1) if m else ""


calls = pd.read_csv(f"{D}/model_logs/trained_model{SUF}.calls.tsv.gz", sep="\t",
                    dtype={"meta_scenario": str, "meta_relative_rank": str})
calls["kind"] = calls["meta_scenario"].fillna("").replace("", "design")
rows = []
for s in ("training", "test"):
    t = pd.read_csv(f"{D}/work/{s}/training_data{SUF}.tsv", sep="\t", usecols=["meta_sample", "taxon_name"] + FEAT)
    t = t.drop_duplicates(["meta_sample", "taxon_name"])
    c = calls[calls["set"] == s].merge(t, on=["meta_sample", "taxon_name"], how="left")
    rows.append(c)
c = pd.concat(rows, ignore_index=True)
present = c["truth"] == 1
called = c["call"].astype(int) == 1
c["error"] = np.where(present & ~called, "FN", np.where(~present & called, "FP", np.where(present, "TP", "TN")))
rep = pd.to_numeric(c["meta_rep_genome"], errors="coerce") == 1
ins = pd.to_numeric(c["meta_insilico_strain"], errors="coerce") == 1
alle = c["allele_copy_share"].fillna(0) > 0.2
c["genome"] = np.where(rep, "representative", np.where(ins, "in-silico strain",
                                                       np.where(alle, "real strain, alleles", "real strain, none")))
c["frag_bin"] = pd.cut(c["fragments"].fillna(0), [-1, 1, 2, 10, 100, 1e12], labels=["1", "2", "3-10", "11-100", ">100"])
cong = pd.to_numeric(c["meta_novel_congener"], errors="coerce") == 1
dist = pd.to_numeric(c["meta_novel_distance"], errors="coerce")
c["absent_kind"] = np.where(cong, pd.cut(dist.fillna(1), [-1, 0.02, 0.05, 0.15, 9],
                                         labels=["held-out congener <0.02", "held-out congener 0.02-0.05",
                                                 "held-out congener 0.05-0.15", "held-out congener farther"]).astype(str),
                            np.where(c["meta_relative_rank"].fillna("") == "genus", "present congeners only",
                                     "no congener in the sample"))

print(f"######## {T}: {D}")
print("== 1. Where the errors are (knob as protal calls; training = species held out, test = final model)")
g = c.groupby(["set", "kind"])
w = pd.DataFrame({"samples": g["meta_sample"].nunique(), "present": g["truth"].sum(),
                  "absent": g.size() - g["truth"].sum(),
                  "FP": g.apply(lambda x: (x["error"] == "FP").sum(), include_groups=False),
                  "FN": g.apply(lambda x: (x["error"] == "FN").sum(), include_groups=False)})
w["TP"] = w["present"] - w["FN"]
w["F1"] = 2 * w["TP"] / (2 * w["TP"] + w["FP"] + w["FN"])
w["FP/sample"] = w["FP"] / w["samples"]
w["FN rate"] = w["FN"] / w["present"]
w["share of the set's FP"] = w["FP"] / w.groupby(level=0)["FP"].transform("sum")
w["share of the set's FN"] = w["FN"] / w.groupby(level=0)["FN"].transform("sum")
print(w.to_string(float_format=lambda x: f"{x:.4f}"))
for s in ("training", "test"):
    x = w.loc[s]
    tp, fp, fn = x["TP"].sum(), x["FP"].sum(), x["FN"].sum()
    d = x.loc[["design"]] if "design" in x.index else x.iloc[:0]
    print(f"   {s}: F1 {2 * tp / (2 * tp + fp + fn):.4f} over all samples; design samples alone "
          f"{2 * d['TP'].sum() / (2 * d['TP'].sum() + d['FP'].sum() + d['FN'].sum()):.4f}")

# ---- 2. False positives ----
err = pd.read_csv(f"{D}/model_logs/error_reads/{T}/taxa.tsv.gz", sep="\t", dtype={"scenario": str})
fp = err[err["error"] == "FP"].copy()


def sources(s):
    """[(species, fragments, in database)] of fragments_on_taxon_from."""
    out = []
    for part in str(s).split(","):
        if ":" not in part:
            continue
        name, n = part.rsplit(":", 1)
        lacks = name.endswith(" (not in the database)")
        out.append((name.replace(" (not in the database)", ""), int(n), not lacks))
    return out


def fp_source(r):
    src = sources(r["fragments_on_taxon_from"])
    tot = sum(n for _, n, _ in src)
    if not tot:
        return "no source traced"
    g = genus(r["taxon_name"])
    agg = {}
    for name, n, indb in src:
        k = ("a species the database lacks" if not indb else "a database species") + \
            (", same genus" if genus(name) == g else ", other genus")
        agg[k] = agg.get(k, 0) + n
    k, n = max(agg.items(), key=lambda kv: kv[1])
    return k if n / tot >= 0.5 else "mixed"


fp["source"] = fp.apply(fp_source, axis=1)
fp["frag_bin"] = pd.cut(fp["fragments_on_taxon"].fillna(0), [-1, 1, 2, 10, 100, 1e12],
                        labels=["1", "2", "3-10", "11-100", ">100"])
fp["kind"] = fp["scenario"].fillna("design")
fp = fp.merge(c[["set", "meta_sample", "taxon_name", "absent_kind", "identity", "low_mapq_share"]],
              left_on=["set", "sample", "taxon_name"], right_on=["set", "meta_sample", "taxon_name"], how="left")
print(f"\n== 2. False positives: {len(fp)} (training {int((fp['set'] == 'training').sum())}, test "
      f"{int((fp['set'] == 'test').sum())}); traced read sources (error_reads)")
for s in ("training", "test"):
    f = fp[fp["set"] == s]
    print(f"-- {s}: by the main source of the taxon's reads (>= half of them) and fragments on the taxon")
    print(pd.crosstab(f["source"], f["frag_bin"], margins=True).to_string())
    print(f"-- {s}: by what the calls table says is beside it (meta_novel_distance etc.)")
    print(pd.crosstab(f["absent_kind"], f["kind"], margins=True).to_string())
f = fp[fp["set"] == "training"]
print("-- training: identity of the taxon's records (median, 10%, 90%) and score by the main source")
print(f.groupby("source").agg(n=("p", "size"), identity_median=("identity", "median"),
                              identity_10=("identity", lambda x: x.quantile(0.1)),
                              identity_90=("identity", lambda x: x.quantile(0.9)),
                              p_median=("p", "median"), within_0_1_of_knob=("p", lambda x: (x < f["knob"].iloc[0] + 0.1).mean()),
                              mapq4_share=("best_on_taxon_mapq4", "sum")).to_string(float_format=lambda x: f"{x:.3f}"))
print("-- training: the score of the false positives (share above each value)")
print("   " + "  ".join(f">= {q}: {(f['p'] >= q).mean():.2f}" for q in (0.8, 0.9, 0.95, 0.99)))
rec = fp.groupby("taxon_name").agg(FP=("sample", "size"), samples=("sample", "nunique"),
                                   sources=("fragments_on_taxon_from", lambda x: "; ".join(sorted(set(
                                       n for s in x for n, _, _ in sources(s)))[:3])))
absent_rows = c[~present].groupby("taxon_name").size().rename("absent_rows")
rec = rec.join(absent_rows).sort_values("FP", ascending=False)
print(f"-- recurring false positives: {int((rec['FP'] >= 3).sum())} taxa with 3 or more ({int(rec.loc[rec['FP'] >= 3, 'FP'].sum())} "
      f"of {len(fp)} FP); the top 12 (FP, rows where the taxon was absent, the reads' sources)")
print(rec.head(12).to_string())

# ---- 3. False negatives ----
print(f"\n== 3. False negatives")
for s in ("training", "test"):
    x = c[(c["set"] == s) & present]
    tab = x.groupby(["genome", "frag_bin"], observed=True).agg(present=("error", "size"),
                                                               FN=("error", lambda e: (e == "FN").sum()))
    tab["FN rate"] = tab["FN"] / tab["present"]
    tab["share of FN"] = tab["FN"] / tab["FN"].sum()
    print(f"-- {s}: present taxa by genome simulated and fragments (after the MAPQ filter)")
    print(tab.to_string(float_format=lambda v: f"{v:.4f}"))
fn = err[err["error"] == "FN"].copy()
fn = fn.merge(c[["set", "meta_sample", "taxon_name", "genome", "frag_bin", "fragments", "fragments_all", "em_fragments",
                 "identity", "low_mapq_share", "allele_explained_share", "congener_fit_share"]],
              left_on=["set", "sample", "taxon_name"], right_on=["set", "meta_sample", "taxon_name"], how="left")
for col in ("own_fragments", "own_best_on_taxon", "own_best_on_taxon_mapq4", "own_best_elsewhere", "own_unaligned",
            "own_seeded_not_aligned"):
    fn[col] = pd.to_numeric(fn[col], errors="coerce").fillna(0)
f = fn[fn["set"] == "training"]
print(f"-- training: where the false negatives' own reads went ({len(f)} FN, {int(f['own_fragments'].sum())} of their "
      f"fragments in the SAMs), by genome")
fate = f.groupby("genome")[["own_fragments", "own_best_on_taxon", "own_best_on_taxon_mapq4", "own_best_elsewhere",
                            "own_unaligned", "own_seeded_not_aligned"]].sum()
for col in fate.columns[1:]:
    fate[col] = fate[col] / fate["own_fragments"]
fate.insert(0, "FN", f.groupby("genome").size())
print(fate.to_string(float_format=lambda v: f"{v:.3f}"))


def where_else(s, g):
    out = {"same genus": 0, "other genus": 0}
    for part in str(s).split(","):
        if ":" in part:
            name, n = part.rsplit(":", 1)
            out["same genus" if genus(name) == g else "other genus"] += int(n)
    return pd.Series(out)


we = f.apply(lambda r: where_else(r["own_best_elsewhere_on"], genus(r["taxon_name"])), axis=1).sum()
print(f"   their reads that aligned best elsewhere: same genus {int(we['same genus'])}, other genus {int(we['other genus'])}")
print("-- training: FN by what kept them below the knob (fragments after the MAPQ filter, the share of low-MAPQ records)")
f2 = f.copy()
f2["why"] = np.select([f2["fragments"].fillna(0) <= 2,
                       f2["low_mapq_share"].fillna(0) >= 0.5,
                       f2["identity"].fillna(1) < 0.97],
                      ["<= 2 fragments", ">= half the records MAPQ-filtered", "identity < 0.97"], "other")
print(pd.crosstab(f2["genome"], f2["why"], margins=True).to_string())
print("-- training: features of FN against TP with 3 or more fragments, by genome (medians)")
x = c[(c["set"] == "training") & present & (c["fragments"].fillna(0) >= 3)]
med = x.groupby(["genome", "error"])[["fragments", "identity", "low_mapq_share", "allele_explained_share",
                                      "congener_fit_share", "excess_scaled_median", "p"]].median()
med["n"] = x.groupby(["genome", "error"]).size()
print(med.to_string(float_format=lambda v: f"{v:.4f}"))
print("-- training: the score of the false negatives (share below each value)")
print("   " + "  ".join(f"< {q}: {(f['p'] < q).mean():.2f}" for q in (0.05, 0.1, 0.3, 0.5)))
