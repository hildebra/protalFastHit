import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lrlib import World

gtdb, db = sys.argv[1], sys.argv[2]
w = World(gtdb, db)
print("genomes", len(w.genomes), "db genes", len(w.db_genes))
for acc in w.genomes:
    print(acc, w.taxid.get(acc), w.rep[acc], {c: len(s) for c, s in w.genomes[acc].items()}, len(w.markers.get(acc, [])))
same = diff = missing = 0
lens = []
for acc, rep in w.rep.items():
    for m in w.markers.get(acc, []):
        name = f"{w.taxid[acc]}_{w.geneid[m[0]]}"
        g = w.gene_seq(acc, m)
        if rep:
            lens.append(len(g))
            if name not in w.db_genes:
                missing += 1
            elif w.db_genes[name] == g:
                same += 1
            else:
                diff += 1
                if diff < 3:
                    d = w.db_genes[name]
                    print("diff", acc, m, len(g), len(d), g[:30], d[:30])
print("rep genes equal to DB", same, "differ", diff, "missing", missing, "max len", max(lens))
# strains: identity to DB gene
import difflib
for acc, rep in w.rep.items():
    if rep:
        continue
    eq = ln = n = 0
    for m in w.markers.get(acc, []):
        name = f"{w.taxid[acc]}_{w.geneid[m[0]]}"
        g = w.gene_seq(acc, m)
        d = w.db_genes.get(name)
        if d is None:
            continue
        n += 1
        if len(g) == len(d):
            eq += 1
            ln += sum(a != b for a, b in zip(g, d))
    print("strain", acc, "genes", n, "same length", eq, "mismatches in same-length", ln)
