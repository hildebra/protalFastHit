"""The local test's samples: one sample drawn at random from each of ten sets drawn at random among those with a truth
table (seed 20261013). Prints set, sample, read type, truth, group and note.

usage: pick_samples.py <scripts/mock_benchmark/samples.tsv>
"""
import csv
import random
import sys

rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
rng = random.Random(20261013)
sets = sorted({r["set"] for r in rows if r["truth"] != "-"})
for s in sorted(rng.sample(sets, 10)):
    r = rng.choice([r for r in rows if r["set"] == s and r["truth"] != "-"])
    print(r["set"], r["sample"], r["read_type"], r["truth"], r["group"], r["note"][:70], sep="\t")
