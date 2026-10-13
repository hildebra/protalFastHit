#!/usr/bin/env python3
"""A protal truth file (--profile_truth: one internal taxid per line) of a mock benchmark truth table's species,
by the database's genome2tiid.tsv (accession, taxid, representative, lineage). Species the database lacks are
reported on stderr.

usage: truth_taxids.py TRUTH.tsv GENOME2TIID > truth.txt
"""
import sys

truth, genome2tiid = sys.argv[1], sys.argv[2]
taxid_of = {}
with open(genome2tiid) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 4:
            taxid_of[f[3].split(";")[-1]] = f[1]
with open(truth) as fh:
    rows = [line.rstrip("\n").split("\t") for line in fh if not line.startswith("#")]
for species, share in ((r[0], r[1]) for r in rows[1:]):
    if species == "?" or share == "0":
        continue
    if species in taxid_of:
        print(taxid_of[species])
    else:
        print(f"not in the database: {species}", file=sys.stderr)
