#!/usr/bin/env python3
"""Drop given variants from a results*.jsonl and delete their analysis dirs so analyze.py redoes them.
usage: redo_variants.py <run_dir> <tag or ''> <variant,...> [species,...]"""
import json, os, shutil, sys
run, tag, variants = sys.argv[1], sys.argv[2], set(sys.argv[3].split(","))
species = set(sys.argv[4].split(",")) if len(sys.argv) > 4 else None
t = f"_{tag}" if tag else ""
res = os.path.join(run, f"results{t}.jsonl")
keep = []
for line in open(res):
    r = json.loads(line)
    if r["variant"] in variants and (species is None or r["species"] in species):
        d = os.path.join(run, "an" + t, r["species"], r["variant"])
        shutil.rmtree(d, ignore_errors=True)
        continue
    keep.append(line)
open(res, "w").writelines(keep)
print(f"{res}: kept {len(keep)} records")
