#!/usr/bin/env python3
"""Traces the false positives of a build_gtdb_database.py run's independent test sets.

For every read type (pe, se, pb, ont) the test set's absent taxa that the model calls present (p >= knob) are
false positives. Each is traced to the reads that made it: the reads aligned to the taxon's genes in the sample's
SAM, grouped by the genome they were simulated from (read names carry it; long-read names carry the genome's
index in the manifest), with the identity of their alignments. The genome's species is in the database or
held out of it (heldout_species.txt).

    python3 trace_false_positives.py BUILD_DIR OUT_DIR [knob]

BUILD_DIR: the build's output folder (trained_model*.test_predictions.tsv.gz, heldout_species.txt, test/).
Needs pandas and numpy, and zstdcat on PATH for .sam.zst. OUT_DIR gets fp_calls.tsv (every false positive with its
closest relatives), fp_reads.tsv (their reads by source genome) and fp_main_source.tsv.
"""
import collections
import glob
import os
import re
import subprocess
import sys

import numpy as np
import pandas as pd

BUILD = os.path.expanduser(sys.argv[1])
OUT = sys.argv[2]
KNOB = float(sys.argv[3]) if len(sys.argv) > 3 else 0.5
RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
TABLES = {"pe": "training_data.tsv", "se": "training_data_se.tsv", "pb": "training_data_pb.tsv",
          "ont": "training_data_ont.tsv"}
PREDS = {"pe": "trained_model.test_predictions.tsv.gz", "se": "trained_model_se.test_predictions.tsv.gz",
         "pb": "trained_model_pb.test_predictions.tsv.gz", "ont": "trained_model_ont.test_predictions.tsv.gz"}
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 40)
pd.set_option("display.max_colwidth", 36)

held = {}
with open(f"{BUILD}/heldout_species.txt") as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if f[0]:
            held[f[0]] = f[1] if len(f) > 1 else "species"

# the paired-end points' manifests: community sample -> species -> lineage, read pairs, abundance
manifests, community_genomes, genome_species = {}, {}, {}
for design in sorted(os.listdir(f"{BUILD}/test/points")):
    path = f"{BUILD}/test/points/{design}/sim/manifest.tsv"
    if not os.path.exists(path):
        continue
    m = pd.read_csv(path, sep="\t")
    for sample, g in m.groupby("sample", sort=False):
        species = {}
        for _, r in g.iterrows():
            s = species.setdefault("s__" + r["species"], {"lineage": r["taxonomy"].split(";"), "pairs": 0})
            s["pairs"] += int(r["read_pairs"])
        manifests[sample] = species
        community_genomes[sample] = list(g["genome"])  # long-read names: g<index>x...
    for _, r in m.iterrows():
        genome_species[r["genome"]] = ("s__" + r["species"], r["taxonomy"])


def community(rt, design, sample):
    if rt == "pe":
        return sample
    if rt == "se":
        return sample[:-3]
    table = pd.read_csv(f"{BUILD}/test/points/{design}/sim/samples.tsv", sep="\t").set_index("sample")
    return table.loc[sample, "community"]


def shared(a, b):
    """Number of leading ranks two lineages (lists of names without prefix) share."""
    n = 0
    for x, y in zip(a, b):
        if x != y:
            break
        n += 1
    return n


# ---- the test tables: every taxon, its call, and the closest simulated relatives of the absent ones ------------
frames, absent = [], []
for rt in TABLES:
    p = pd.read_csv(f"{BUILD}/{PREDS[rt]}", sep="\t", low_memory=False)
    feats = pd.read_csv(f"{BUILD}/test/{TABLES[rt]}", sep="\t", low_memory=False)
    feats = feats[["meta_sample", "taxon_name", "total_hits", "identity", "top_identity", "depth"]]
    feats = feats.drop_duplicates(["meta_sample", "taxon_name"])
    p = p.merge(feats, on=["meta_sample", "taxon_name"], how="left")
    p["read_type"], p["call"] = rt, p["p"] >= KNOB
    frames.append(p)
    for _, r in p[p["truth"] == 0].iterrows():
        lineage = [r["domain"], r["meta_lineage_phylum"][3:], r["meta_lineage_class"][3:], r["meta_lineage_order"][3:],
                   r["meta_lineage_family"][3:], r["meta_lineage_genus"][3:], r["taxon_name"][3:]]
        species = manifests[community(rt, r["meta_design"], r["meta_sample"])]
        best = {"novel": (0, None), "present": (0, None)}  # closest simulated species that the database lacks / has
        for name, s in species.items():
            d = shared(lineage, [x[3:] for x in s["lineage"]])
            kind = "novel" if name in held else "present"
            if d > best[kind][0] or (d == best[kind][0] and d > 0 and s["pairs"] > species[best[kind][1]]["pairs"]):
                best[kind] = (d, name)
        nd, nn = best["novel"]
        pd_, pn = best["present"]
        primary = "novel" if nd >= pd_ and nd > 0 else ("present" if pd_ > 0 else "none")
        absent.append({"read_type": rt, "design": r["meta_design"], "sample": r["meta_sample"], "taxon": r["taxon_name"],
                       "taxid": r["taxon"], "p": r["p"], "call": bool(r["call"]), "primary": primary,
                       "novel_rank": RANKS[nd - 1] if nd else "", "novel_species": nn or "",
                       "novel_heldout": held.get(nn, ""), "novel_pairs": species[nn]["pairs"] if nn else 0,
                       "present_rank": RANKS[pd_ - 1] if pd_ else "", "present_species": pn or "",
                       "present_pairs": species[pn]["pairs"] if pn else 0,
                       "total_hits": r["total_hits"], "identity": r["identity"], "top_identity": r["top_identity"]})
a = pd.concat(frames)
ab = pd.DataFrame(absent)
fp_lineage = {r["taxon_name"]: ["d__" + r["domain"], r["meta_lineage_phylum"], r["meta_lineage_class"],
                                r["meta_lineage_order"], r["meta_lineage_family"], r["meta_lineage_genus"], r["taxon_name"]]
              for _, r in a.drop_duplicates("taxon_name").iterrows()}
fp = ab[ab["call"]].copy()
os.makedirs(OUT, exist_ok=True)
fp.to_csv(f"{OUT}/fp_calls.tsv", sep="\t", index=False)


def depth_of(design):
    return int(re.search(r"_[pb](\d+)", design).group(1))


print(f"## False positives at knob {KNOB}: absent taxa called present\n")
ab["depth"] = ab["design"].map(depth_of)
t = ab.groupby("read_type").agg(absent=("call", "size"), fp=("call", "sum"), samples=("sample", "nunique"))
t["fp_per_sample"] = (t["fp"] / t["samples"]).round(2)
t["fp_rate_%"] = (100 * t["fp"] / t["absent"]).round(2)
print(t.to_string(), "\n")
print("by depth (read pairs, or bases for long reads):")
d = ab.groupby(["read_type", "depth"]).agg(absent=("call", "size"), fp=("call", "sum"))
d["fp_rate_%"] = (100 * d["fp"] / d["absent"]).round(2)
print(d.to_string(), "\n")

print("## Where the false positives sit relative to the sample's simulated species\n")
ab["relation"] = np.where(ab["primary"] == "novel", "closest simulated species is held out of the database, shares " + ab["novel_rank"],
                          np.where(ab["primary"] == "present", "closest simulated species is in the database, shares " + ab["present_rank"],
                                   "no simulated relative"))
r = ab.groupby("relation").agg(absent=("call", "size"), fp=("call", "sum"))
r["fp_rate_%"] = (100 * r["fp"] / r["absent"]).round(2)
print(r.sort_values("absent", ascending=False).to_string(), "\n")

# ---- the reads of every false positive, by the genome they were simulated from ----------------------------------
def sam_path(rt, design, sample):
    base = f"{BUILD}/test/points/{design[:-3] if rt == 'se' else design}"
    hits = [h for h in glob.glob(f"{base}/{'protal_se' if rt == 'se' else 'protal'}/alignments/{sample}.sam*")
            if not h.endswith(".err")]
    return hits[0] if hits else None


def identity(cigar):
    c = collections.Counter()
    for n, op in re.findall(r"(\d+)([MXIDS=])", cigar):
        c[op] += int(n)
    aligned = c["M"] + c["="] + c["X"] + c["I"] + c["D"]
    return (c["M"] + c["="]) / max(1, aligned)


def source_of(qname, rt, comm):
    if rt in ("pe", "se"):
        return re.match(r"(GC[AF]_\d+\.\d+)", qname).group(1)
    return community_genomes[comm][int(re.match(r"g(\d+)x", qname).group(1))]


rows = []
for (rt, design, sample), g in fp.groupby(["read_type", "design", "sample"]):
    path = sam_path(rt, design, sample)
    if path is None:
        print("no SAM for", rt, sample)
        continue
    comm = community(rt, design, sample)
    taxids = {str(r["taxid"]): r for _, r in g.iterrows()}
    cat = "zstdcat" if path.endswith(".zst") else "zcat" if path.endswith(".gz") else "cat"
    awk = ('BEGIN{n=split(T,a," "); for(i=1;i<=n;i++) t[a[i]]=1} !/^@/ {split($3,r,"_"); if (r[1] in t) '
           'print $1"\\t"$3"\\t"$5"\\t"$6}')
    out = subprocess.run(f"{cat} '{path}' | awk -F'\\t' -v T='{' '.join(taxids)}' '{awk}'", shell=True,
                         capture_output=True, text=True).stdout
    per = collections.defaultdict(lambda: collections.defaultdict(list))
    for line in out.splitlines():
        qname, rname, mapq, cigar = line.split("\t")
        per[rname.split("_")[0]][source_of(qname, rt, comm)].append((qname, identity(cigar), int(mapq)))
    for taxid, r in taxids.items():
        sources = per.get(taxid, {})
        fp_reads = sum(len({q for q, _, _ in v}) for v in sources.values())
        for genome, v in sources.items():
            species, lineage = genome_species[genome]
            n = shared(fp_lineage[r["taxon"]], lineage.split(";"))
            rows.append({"read_type": rt, "sample": sample, "fp_taxon": r["taxon"], "p": round(r["p"], 3),
                         "fp_reads": fp_reads, "source_genome": genome, "source_species": species,
                         "source_status": f"held out ({held[species]})" if species in held else "in database",
                         "shared_rank": RANKS[n - 1] if n else "none",
                         "reads_from_source": len({q for q, _, _ in v}),
                         "mean_identity": round(float(np.mean([i for _, i, _ in v])), 3),
                         "same_genus_as_fp": species.split(" ")[0] == r["taxon"].split(" ")[0]})
reads = pd.DataFrame(rows)
reads.to_csv(f"{OUT}/fp_reads.tsv", sep="\t", index=False)
main = (reads.sort_values("reads_from_source", ascending=False).groupby(["read_type", "sample", "fp_taxon"], sort=False)
        .first().reset_index())
agg = reads.groupby(["read_type", "sample", "fp_taxon"]).apply(lambda g: pd.Series({
    "fp_reads": g["fp_reads"].iloc[0], "sources": len(g),
    "reads_from_held_out": int(g.loc[g["source_status"].str.startswith("held"), "reads_from_source"].sum()),
    "reads_from_database_species": int(g.loc[g["source_status"] == "in database", "reads_from_source"].sum()),
    "reads_from_same_genus": int(g.loc[g["same_genus_as_fp"], "reads_from_source"].sum())}), include_groups=False).reset_index()
agg = agg.merge(main[["read_type", "sample", "fp_taxon", "source_species", "source_status", "mean_identity", "p"]],
                on=["read_type", "sample", "fp_taxon"])
agg.to_csv(f"{OUT}/fp_main_source.tsv", sep="\t", index=False)

print("## Reads of the false positives, by the genome they were simulated from\n")
print(reads[["read_type", "sample", "fp_taxon", "p", "fp_reads", "source_species", "source_status", "shared_rank",
             "reads_from_source", "mean_identity"]].to_string(index=False), "\n")
print("## Per false positive: where its reads come from\n")
agg["main_source_is"] = np.where(agg["reads_from_held_out"] >= agg["reads_from_database_species"],
                                 "a species held out of the database", "a species in the database")
agg["main_source_same_genus"] = agg["reads_from_same_genus"] * 2 >= agg["fp_reads"]
print(agg.drop(columns=["main_source_same_genus"]).to_string(index=False), "\n")
print(agg.groupby(["read_type", "main_source_is"]).size().unstack(fill_value=0).to_string(), "\n")
print(f"reads of a false positive: median {agg['fp_reads'].median():.0f}, mean {agg['fp_reads'].mean():.1f}, "
      f"<= 3 reads: {(agg['fp_reads'] <= 3).sum()} of {len(agg)}, >= 20 reads: {(agg['fp_reads'] >= 20).sum()}; "
      f"most reads from the same genus: {int(agg['main_source_same_genus'].sum())} of {len(agg)}")
print("alignment identity of those reads:", reads["mean_identity"].describe().round(3).to_dict(), "\n")

# ---- why the model calls them ------------------------------------------------------------------------------------
print("## Why the model calls them: the evidence behind a call\n")
a["hits_bin"] = pd.cut(a["total_hits"], [0, 1, 3, 8, 48, 10**9], labels=["1", "2-3", "4-8", "9-48", "49+"])
g = a.groupby(["hits_bin", "truth"], observed=True).agg(taxa=("call", "size"), called=("call", "sum")).unstack(fill_value=0)
g.columns = [f"{x}_{'present' if y else 'absent'}" for x, y in g.columns]
g["FP_rate_of_absent_%"] = (100 * g["called_absent"] / g["taxa_absent"]).round(2)
g["sensitivity_%"] = (100 * g["called_present"] / g["taxa_present"]).round(1)
print("read hits on the taxon:\n" + g.to_string(), "\n")
low = a[a["call"] & (a["total_hits"] <= 8)]
print("identity of called taxa with up to 8 hits (0 = absent, 1 = present):")
print(low.groupby("truth")[["identity", "top_identity"]].describe().T.round(3).to_string(), "\n")
print("identity of the true positives by read type and read hits (quantiles):")
tp = a[(a["truth"] == 1) & a["call"]].copy()
tp["hits"] = np.where(tp["total_hits"] >= 20, ">=20", np.where(tp["total_hits"] <= 8, "<=8", "9-19"))
q = tp.groupby(["read_type", "hits"])["identity"].quantile([0.05, 0.25, 0.5]).unstack().round(3)
q["taxa"] = tp.groupby(["read_type", "hits"]).size()
print(q.to_string(), "\n")
print("threshold, false positives and false negatives (all read types):")
rows = []
for th in (0.5, 0.6, 0.7, 0.8, 0.9, 0.95):
    c = a["p"] >= th
    rows.append({"threshold": th, "FP": int(((a["truth"] == 0) & c).sum()), "FN": int(((a["truth"] == 1) & ~c).sum())})
print(pd.DataFrame(rows).to_string(index=False), "\n")
fps = a[(a["truth"] == 0) & a["call"]].copy()
print("p of the false positives:\n" + pd.cut(fps["p"], [0.5, 0.6, 0.7, 0.8, 0.9, 1.0], right=False).value_counts().sort_index().to_string(), "\n")
fps["community"] = np.where(fps["read_type"].isin(["pb", "ont"]), fps["meta_sample"], fps["meta_sample"].str.replace("_se$", "", regex=True))
u = fps.groupby(["community", "taxon_name"])["read_type"].apply(lambda s: ",".join(sorted(s)))
print(f"{len(fps)} calls are {len(u)} distinct (sample, taxon) events; in more than one read type:\n{u[u.str.contains(',')].to_string()}")
