"""Follow-ups of error_anatomy.py (one read type, training samples = species held out, knob as protal calls):

1. the false negatives' own reads that aligned best to another species: was that species present in the sample (a
   congener that absorbed them), absent and called (they made a false positive), or absent and not called;
2. what the MAPQ filter took: fragments before (fragments_all) and after it, and the EM's share (em_fragments);
3. the false negatives simulated from the representative with 3 or more fragments: whose reads are on the taxon;
4. the false positives' pairs (taxon, the species that lent it its reads): how often each pair recurs, and the
   distance of the nearest held-out congener in the sample (meta_novel_distance: species_clouds.tsv's marker genes);
5. the errors per sample against the sample's present species and depth.

usage: deep_dive.py <build dir> <read type>
"""
import re
import sys
import numpy as np
import pandas as pd

D, T = sys.argv[1], sys.argv[2]
SUF = "" if T == "pe" else "_" + T
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


calls = pd.read_csv(f"{D}/model_logs/trained_model{SUF}.calls.tsv.gz", sep="\t",
                    usecols=["meta_sample", "meta_scenario", "meta_rep_genome", "meta_insilico_strain", "taxon_name",
                             "meta_novel_distance", "truth", "set", "p", "knob", "call"], dtype={"meta_scenario": str})
calls = calls[calls["set"] == "training"]
t = pd.read_csv(f"{D}/work/training/training_data{SUF}.tsv", sep="\t",
                usecols=["meta_sample", "taxon_name", "fragments", "fragments_all", "em_fragments", "identity",
                         "allele_copy_share", "sample_log_fragments"]).drop_duplicates(["meta_sample", "taxon_name"])
c = calls.merge(t, on=["meta_sample", "taxon_name"], how="left")
rep = pd.to_numeric(c["meta_rep_genome"], errors="coerce") == 1
ins = pd.to_numeric(c["meta_insilico_strain"], errors="coerce") == 1
alle = c["allele_copy_share"].fillna(0) > 0.2
c["genome"] = np.where(rep, "representative", np.where(ins, "in-silico strain",
                                                       np.where(alle, "real strain, alleles", "real strain, none")))
state = {(s, n): ("present" if tr == 1 else ("absent, called" if cl == 1 else "absent, not called"))
         for s, n, tr, cl in zip(c["meta_sample"], c["taxon_name"], c["truth"], c["call"])}

err = pd.read_csv(f"{D}/model_logs/error_reads/{T}/taxa.tsv.gz", sep="\t", dtype={"scenario": str})
err = err[err["set"] == "training"]
fn = err[err["error"] == "FN"].merge(c[["meta_sample", "taxon_name", "genome", "fragments", "fragments_all",
                                        "em_fragments", "identity"]],
                                     left_on=["sample", "taxon_name"], right_on=["meta_sample", "taxon_name"], how="left")

print(f"######## {T}: {D} (training samples, species held out)")
print("== 1. The false negatives' own reads that aligned best to another species, by that species' state in the sample")
rows = []
for _, r in fn.iterrows():
    for name, n, _ in parts(r["own_best_elsewhere_on"]):
        rows.append((r["genome"], state.get((r["sample"], name), "no row (not in the sample's table)"),
                     "same genus" if genus(name) == genus(r["taxon_name"]) else "other genus", n))
el = pd.DataFrame(rows, columns=["genome", "state", "genus", "fragments"])
print(el.pivot_table(index=["genome", "genus"], columns="state", values="fragments", aggfunc="sum", fill_value=0).to_string())
own_aligned = fn["own_best_on_taxon"].fillna(0) + fn["own_best_elsewhere"].fillna(0)
share = fn["own_best_elsewhere"].fillna(0) / own_aligned.replace(0, np.nan)
print("   per FN with 3 or more own aligned fragments: median share of them that aligned best elsewhere, by genome")
print(share[own_aligned >= 3].groupby(fn["genome"]).describe()[["count", "25%", "50%", "75%"]].to_string(
    float_format=lambda v: f"{v:.3f}"))

print("\n== 2. What the MAPQ filter took: FN by fragments after it (fragments) and before it (fragments_all)")
x = c[c["truth"] == 1].copy()
x["after"] = pd.cut(x["fragments"].fillna(0), [-1, 2, 10, 1e12], labels=["<=2", "3-10", ">10"])
x["before"] = pd.cut(x["fragments_all"].fillna(0), [-1, 2, 10, 1e12], labels=["<=2", "3-10", ">10"])
x["FN"] = x["call"] == 0
print(pd.crosstab([x["genome"], x["after"]], x["before"], values=x["FN"], aggfunc="sum").fillna(0).astype(int).to_string())
print("   FN with <= 2 fragments after the filter but > 10 before it, and their EM fragments (median):",
      int((x["FN"] & (x["after"] == "<=2") & (x["before"] == ">10")).sum()),
      x.loc[x["FN"] & (x["after"] == "<=2") & (x["before"] == ">10"), "em_fragments"].median())

print("\n== 3. Representatives missed with 3 or more fragments: the reads on the taxon and their sources")
r3 = fn[(fn["genome"] == "representative") & (fn["fragments"].fillna(0) >= 3)].copy()
own = r3["own_best_on_taxon"].fillna(0)
on = pd.to_numeric(r3["fragments_on_taxon"], errors="coerce").fillna(0)
r3["own_share"] = own / on.replace(0, np.nan)
print(f"   {len(r3)} FN; own share of the reads on the taxon: median {r3["own_share"].median():.3f}; "
      f"own reads that aligned elsewhere: {int(r3['own_best_elsewhere'].fillna(0).sum())} of "
      f"{int((r3['own_best_on_taxon'].fillna(0) + r3['own_best_elsewhere'].fillna(0)).sum())} aligned")
print(r3.sort_values("fragments", ascending=False)[["sample", "taxon_name", "p", "fragments", "identity", "own_best_on_taxon",
                                                     "own_best_elsewhere", "fragments_on_taxon_from"]].head(10).to_string(
    index=False, max_colwidth=90))

print("\n== 4. The false positives' (taxon, main source) pairs")
fp = err[err["error"] == "FP"].copy()


def main_source(s):
    p = parts(s)
    if not p:
        return ("", 0, True)
    return max(p, key=lambda q: q[1])


fp["source"], fp["source_fragments"], fp["source_in_db"] = zip(*fp["fragments_on_taxon_from"].map(main_source))
fp = fp.merge(c[["meta_sample", "taxon_name", "meta_novel_distance"]], left_on=["sample", "taxon_name"],
              right_on=["meta_sample", "taxon_name"], how="left")
pairs = fp.groupby(["taxon_name", "source"]).agg(FP=("sample", "size"), in_db=("source_in_db", "first"),
                                                 fragments=("source_fragments", "median"),
                                                 held_out_congener_distance=("meta_novel_distance", "median")).reset_index()
print(f"   {len(pairs)} pairs; {int((pairs['FP'] >= 3).sum())} pairs with 3 or more FP hold "
      f"{int(pairs.loc[pairs['FP'] >= 3, 'FP'].sum())} of {len(fp)} FP")
print(pairs.sort_values("FP", ascending=False).head(15).to_string(index=False))

print("\n== 5. Errors per sample (design and scenario samples)")
c["kind"] = c["meta_scenario"].fillna("").replace("", "design")
per = c.groupby(["kind", "meta_sample"]).agg(present=("truth", "sum"), rows=("truth", "size"),
                                             FP=("call", lambda s: 0), depth=("sample_log_fragments", "first"))
per["FP"] = c[(c["truth"] == 0) & (c["call"] == 1)].groupby(["kind", "meta_sample"]).size()
per["FN"] = c[(c["truth"] == 1) & (c["call"] == 0)].groupby(["kind", "meta_sample"]).size()
per = per.fillna(0)
per["FP per 100 present"] = 100 * per["FP"] / per["present"].replace(0, np.nan)
per["FN rate"] = per["FN"] / per["present"].replace(0, np.nan)
print(per.groupby(level=0).agg(samples=("present", "size"), present_median=("present", "median"),
                               FP_per_100_present=("FP per 100 present", "median"),
                               FN_rate=("FN rate", "median")).to_string(float_format=lambda v: f"{v:.3f}"))
d = per.loc["design"].copy()
d["present_bin"] = pd.cut(d["present"], [0, 30, 60, 100, 1e9], labels=["<=30", "31-60", "61-100", ">100"])
print("   design samples by present species: FP per 100 present, FN rate (medians)")
print(d.groupby("present_bin", observed=True)[["FP per 100 present", "FN rate"]].median().to_string(
    float_format=lambda v: f"{v:.3f}"))
