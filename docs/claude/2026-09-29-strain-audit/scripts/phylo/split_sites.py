#!/usr/bin/env python3
"""For each true split of a species' strain tree, count marker-gene columns whose true tip pattern
supports it exactly (one base on one side, another on the other side), and report the number of
true singleton sites per tip.  usage: split_sites.py <strains/SP dir>"""
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import phylo_utils as pu  # noqa: E402

d = sys.argv[1]
tree = pu.fix_parents(pu.unroot(pu.read_tree(os.path.join(d, "true.nwk"))))
tips = pu.leaf_names(tree)
sp = pu.splits(tree, with_info=True)
seqs = defaultdict(list)
with open(os.path.join(d, "truth_markers.tsv")) as fh:
    next(fh)
    rows = defaultdict(dict)
    for line in fh:
        tip, gid, marker, seq, nins = line.rstrip("\n").split("\t")
        rows[gid][tip] = seq
support = Counter()
single = Counter()
anchor = min(tips)
for gid, bytip in rows.items():
    L = len(next(iter(bytip.values())))
    for p in range(L):
        col = {t: bytip[t][p] for t in tips}
        if "-" in col.values():
            continue
        cnt = Counter(col.values())
        if len(cnt) == 2:
            (b1, n1), (b2, n2) = cnt.most_common()
            side = frozenset(t for t in tips if col[t] == b2)
            side = side if anchor not in side else frozenset(tips) - side
            if len(side) == 1 or len(side) == len(tips) - 1:
                single[next(iter(side if len(side) == 1 else frozenset(tips) - side))] += 1
            else:
                support[side] += 1
print("split\tbranch_length\tsupporting_marker_sites")
for s, (bl, _) in sorted(sp.items(), key=lambda kv: kv[1][0]):
    print(f"{','.join(sorted(s))}\t{bl:.2e}\t{support.get(s, 0)}")
print("tip singleton marker sites:", dict(sorted(single.items())))
