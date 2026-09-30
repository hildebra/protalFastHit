#!/usr/bin/env python3
"""group_bias.py SAMPLE VARIANT... - where a run's abundance error sits, by kind of species: the true and the
reported relative abundance summed over species simulated from their representative alone, from another
strain alone, and from a mixture of genomes (the representative and others, or several others), and each
kind's share of the Bray-Curtis dissimilarity. Runs as score.py finds them (RUNS, the test points)."""
import collections
import csv
import os
import sys

variants, sample = sys.argv[2:], sys.argv[1]
sys.argv = sys.argv[:1]
import score  # noqa: E402

present, true_ab = score.truth_of(sample)
point = sample.rsplit("_s_", 1)[0]
genomes = collections.defaultdict(set)
with open(os.path.join(score.point_dir(point), "sim", "manifest.tsv")) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        if row["sample"] == sample:
            genomes[row["species"]].add(row["genome"])


def kind(species):
    g = genomes.get(species, set())
    if not g:
        return "absent from the sample"
    if len(g) == 1:
        return "representative alone" if next(iter(g)).startswith("GCF_") else "another strain alone"
    return "mixture with the representative" if any(x.startswith("GCF_") for x in g) else "mixture of other strains"


est = {}
for v in variants:
    got = score.reported(os.path.join(score.RUNS, f"{v}.{sample}"))
    total = sum(got.values()) or 1.0
    est[v] = {s: x / total for s, x in got.items()}
held = score.heldout() if any("ho" in v for v in variants) else set()
if held:
    present -= held
    t = sum(true_ab[s] for s in present)
    true_ab = {s: true_ab[s] / t for s in present}
species = set(true_ab) | set().union(*(set(e) for e in est.values()))
kinds = sorted({kind(s) for s in species})
print(f"{sample}: true / reported abundance by kind of species (Bray-Curtis part)")
print("kind".ljust(34) + "true    " + "".join(f"{v:<22}" for v in variants))
for k in kinds:
    sp = [s for s in species if kind(s) == k]
    t = sum(true_ab.get(s, 0) for s in sp)
    cells = []
    for v in variants:
        r = sum(est[v].get(s, 0) for s in sp)
        bc = 0.5 * sum(abs(true_ab.get(s, 0) - est[v].get(s, 0)) for s in sp)
        cells.append(f"{r:.3f} ({bc:.4f})".ljust(22))
    print(f"{k} ({len(sp)})".ljust(34) + f"{t:.3f}   " + "".join(cells))
