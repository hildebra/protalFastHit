#!/usr/bin/env python3
"""Simulate a few paired reads from reference.fna genes (100bp, 0.5% errors, fixed seed)."""
import os, random, sys

ref = sys.argv[1]; outdir = sys.argv[2]; prefix = sys.argv[3]
os.makedirs(outdir, exist_ok=True)
COMP = str.maketrans("ACGTacgt", "TGCATGCA")
def rc(s): return s.translate(COMP)[::-1]
genes, name = [], None
with open(ref) as fh:
    for line in fh:
        line = line.strip()
        if line.startswith(">"): name = line[1:].split()[0]
        elif line: genes.append(line.upper())
rng = random.Random(42)
def mutate(seq):
    b = list(seq)
    for i, c in enumerate(b):
        if rng.random() < 0.005: b[i] = rng.choice([x for x in "ACGT" if x != c])
    return "".join(b)
n = 0
with open(os.path.join(outdir, f"{prefix}_R1.fq"), "w") as r1, open(os.path.join(outdir, f"{prefix}_R2.fq"), "w") as r2:
    for gene in genes:
        if len(gene) < 320: continue
        for _ in range(30):
            flen = rng.randint(220, 320); start = rng.randint(0, len(gene) - flen)
            frag = gene[start:start+flen]
            if rng.random() < 0.5: frag = rc(frag)
            s1 = mutate(frag[:100]); s2 = mutate(rc(frag[-100:])); n += 1
            r1.write(f"@{prefix}.{n}/1\n{s1}\n+\n{'I'*100}\n")
            r2.write(f"@{prefix}.{n}/2\n{s2}\n+\n{'I'*100}\n")
print(f"{n} pairs")
