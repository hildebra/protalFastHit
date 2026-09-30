import sys, os
sys.path.insert(0, os.path.expanduser("~/audit5/qc"))
from consistency import read_fasta
d = os.path.expanduser("~/audit5/qc/runs/e5/protal/strains/")
for f in ("s__Mockella_alpha.raw.msa.fna", "s__Mockella_alpha.msa.fna"):
    n, s = read_fasta(d + f)
    idx = [i for i, x in enumerate(n) if x == "sample_2"]
    a, b = s[idx[0]], s[idx[1]]
    print(f, "rows named sample_2:", len(idx), "identical:", a == b, "differing columns:", sum(x != y for x, y in zip(a, b)))
import csv
from collections import Counter
c = Counter((r["sample"], r["gene_id"]) for r in csv.DictReader(open(d + "s__Mockella_alpha.meta.tsv"), delimiter="\t"))
print("meta (sample, gene) keys appearing twice:", sum(1 for v in c.values() if v > 1))
