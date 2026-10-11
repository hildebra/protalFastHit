"""The report's two tables of section 4 from ablation_<read type>.txt (ablation_table.py's output): each group left out
alone, and the combinations, species held out (and the test set for the combinations), per read type found.

usage: group_tables.py <report dir>
"""
import os
import re
import sys

D = sys.argv[1]
TYPES = [t for t in ("pe", "se", "pb", "ont") if os.path.exists(f"{D}/ablation_{t}.txt")]
ROW = re.compile(r"^\s*(\S+)\s+([\d.]+)\s+([+-][\d.]+ \[[^\]]+\])?\s*([\d.]+)\s+([+-][\d.]+ \[[^\]]+\])?")
GROUPS = ["adjacency", "distance", "depth", "divergence", "unfiltered", "ref", "complexity", "consistency", "shape",
          "neighbourhood", "ancestry", "gaps", "untried", "alleles", "polymorphic", "weights"]
COMBOS = [("no_ancestry", "ancestry"), ("no_anc_dist", "ancestry, distance"),
          ("no_site_columns", "the 8 ancestry-site columns only"),
          ("no_sites", "ancestry, polymorphic, weights"),
          ("no_strain4", "ancestry, alleles, polymorphic, weights"),
          ("no_strain6", "the same + gaps, untried"), ("core", "all but the core set (normalized..unfiltered)")]

res = {}
for t in TYPES:
    for line in open(f"{D}/ablation_{t}.txt"):
        m = ROW.match(line)
        if m and m.group(3):
            res[(t, m.group(1))] = (m.group(3), m.group(5) or "")


def short(v):  # "+0.0002 [-0.0002, +0.0006]" -> "+0.0002 [−0.0002, +0.0006]"
    return v.replace("-", "−") if v else "running"


print("| left out | " + " | ".join(TYPES) + " |")
print("|---|" + "---|" * len(TYPES))
for g in GROUPS:
    print(f"| {g} | " + " | ".join(short(res.get((t, "no_" + g), ("", ""))[0]) for t in TYPES) + " |")
print()
print("| left out | " + " | ".join(f"{t} held out | {t} test" for t in TYPES) + " |")
print("|---|" + "---|" * (2 * len(TYPES)))
for name, label in COMBOS:
    cells = []
    for t in TYPES:
        h, s = res.get((t, name), ("", ""))
        cells += [short(h), short(s)]
    print(f"| {label} | " + " | ".join(cells) + " |")
