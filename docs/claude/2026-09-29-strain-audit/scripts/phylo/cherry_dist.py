#!/usr/bin/env python3
"""True vs estimated patristic distance of the closest sample pairs (cherries).
usage: cherry_dist.py <run> <tag or ''> <species> <variant>..."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import phylo_utils as pu  # noqa: E402

PHYLO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
run, tag, sp = sys.argv[1:4]
t = pu.fix_parents(pu.read_tree(os.path.join(PHYLO, "strains", sp, "true.nwk")))
pt = pu.patristic(t)
leaves = pu.leaf_names(t)
# cherries of the true tree
pairs = []
for nd in t.nodes():
    if len(nd.children) == 2 and all(c.is_leaf() for c in nd.children):
        pairs.append(tuple(sorted(c.name for c in nd.children)))
for v in sys.argv[4:]:
    p = os.path.join(PHYLO, "runs", run, "an" + (f"_{tag}" if tag else ""), sp, v, "iq.treefile")
    e = pu.fix_parents(pu.read_tree(p))
    pe = pu.patristic(e)
    parts = []
    for a, b in pairs:
        est = pe.get((a, b))
        parts.append(f"{a}-{b} {pt[(a, b)]*1e4:.1f}->{est*1e4:.1f}" if est is not None else f"{a}-{b} dropped")
    print(f"{run} {sp} {v}: " + "; ".join(parts) + "  (1e-4 subst/site, true->est)")
