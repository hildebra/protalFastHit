"""Where the '-' of each MSA row lie: per gene (partition), the gap count and whether they are at
the gene's ends. msa_gaps.py MSA PARTITIONS ROW..."""
import re
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lrlib import read_fasta_list

rows = dict(read_fasta_list(sys.argv[1]))
parts = []
for line in open(sys.argv[2]):
    m = re.match(r"DNA, gene(\d+) = (\d+)-(\d+)", line.strip())
    if m:
        parts.append((int(m.group(1)), int(m.group(2)) - 1, int(m.group(3))))
for name in sys.argv[3:]:
    row = rows[name]
    out = []
    for gid, s, e in parts:
        seg = row[s:e]
        g = seg.count("-")
        if g:
            lead = len(seg) - len(seg.lstrip("-"))
            trail = len(seg) - len(seg.rstrip("-"))
            out.append(f"gene{gid}:{g}/{len(seg)} (lead {lead}, trail {trail})")
    print(name, "genes with gaps:", len(out))
    for x in out[:15]:
        print("   ", x)
