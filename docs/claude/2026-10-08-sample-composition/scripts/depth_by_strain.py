#!/usr/bin/env python3
"""Protal depth against the truth by what the sample holds of a species: its reference genome, another genome, two strains
(db_all run of run_validation.sh; run in the validation folder): results/depth_by_strain.txt."""
import csv, glob, collections, statistics
truth = collections.defaultdict(lambda: collections.defaultdict(lambda: [0.0, 0, []]))
with open("sims/manifest.tsv") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        t = truth[row["sample"]][row["taxonomy"].split(";")[-1]]
        t[0] += float(row["vertical_coverage"]); t[1] += 1; t[2].append(row["genome"])
div = {}
with open("gtdb/simulation/divergence.tsv") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        div[row["accession"]] = float(row["strain_divergence"])
out = collections.defaultdict(list)
for sample in truth:
    comp = next(csv.DictReader(open(glob.glob(f"run_db_all/**/{sample}.profile.composition", recursive=True)[0]), delimiter="\t"))
    share = float(comp["FragmentBaseShare"])
    log = {r["Name"]: (float(r["VCov"]), float(r["LowIdentityShare"])) for r in csv.DictReader(open(glob.glob(f"run_db_all/**/{sample}.profile.log", recursive=True)[0]), delimiter="\t")}
    for name, (cov, n, genomes) in truth[sample].items():
        if name not in log: continue
        reps = [g for g in genomes if g.startswith("GCF")]
        key = ("one strain" if n == 1 else "two strains") + (", the reference" if n == 1 and reps else ", another genome" if n == 1 else "")
        out[key].append((log[name][0] / (cov * share), log[name][1]))
for k, v in sorted(out.items()):
    print(f"{k}: protal depth / truth median {statistics.median(r for r, _ in v):.4f} (n={len(v)}), low-identity share median {statistics.median(s for _, s in v):.4f}")
print("strain divergences:", sorted(set(round(x, 4) for x in div.values()))[:5], "...")
