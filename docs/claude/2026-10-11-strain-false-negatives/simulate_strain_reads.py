"""Paired-end reads from the marker genes of the strain world's test strains (halves.tsv: never in a database) and of
congener representatives, in protal map rows, with each sample's truth.

Two sample sets of SAMPLES samples each; sample r holds the test strain of rank r of every focal species (its only
strain of that species) and BACKGROUND congener representatives (other species of the focal genera, which are in
every database as they are), at 20 to 200 pairs each:
  deep     every test strain at DEEP pairs: how its reads align (read level);
  shallow  each test strain at one of SHALLOW pairs (cycled over species and samples): whether it is called.
Fragments: a position drawn over the genome's concatenated marker genes (so by gene length), length N(350, 50) in
[200, 600] cut to the gene, either strand; 150 bp mates (shorter genes give shorter mates, at least 100 bp);
substitution errors from 0.1% at a read's start to 0.6% at its end. Read names: <accession>__<n>.

usage: simulate_strain_reads.py <world dir> <out dir> [seed]
"""
import glob
import gzip
import os
import random
import sys
from collections import defaultdict

W, OUT = sys.argv[1], sys.argv[2]
SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 1
SAMPLES, DEEP, BACKGROUND = 20, 60, 15
SHALLOW = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32]
READ = 150
COMP = str.maketrans("ACGTN", "TGCAN")

rows = [line.rstrip("\n").split("\t") for line in open(f"{W}/taxa.tsv")][1:]
halves = [line.rstrip("\n").split("\t") for line in open(f"{W}/halves.tsv")][1:]
test = defaultdict(dict)  # species -> rank -> accession
for acc, half, sp, st, rank in halves:
    if half == "test":
        test[sp][int(rank)] = acc
species_of = {r[0]: r[2] for r in rows}
congener_reps = sorted(r[0] for r in rows if r[1] == "representative" and r[5] == "congener")
wanted = {a for ranks in test.values() for a in ranks.values()} | set(congener_reps)

genes = defaultdict(list)  # accession -> [sequence]
for path in sorted(glob.glob(f"{W}/genomic_files_all/**/fna/*.fna", recursive=True)):
    acc, seq = None, []
    for line in open(path):
        if line.startswith(">"):
            if acc in wanted and seq:
                genes[acc].append("".join(seq).upper())
            acc, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
    if acc in wanted and seq:
        genes[acc].append("".join(seq).upper())
missing = sorted(a for a in wanted if not genes.get(a))
print(f"{len(genes)} genomes' marker genes read, {len(missing)} without any")

rng = random.Random(SEED)


def mutate(read):
    out = list(read)
    for i in range(len(out)):
        if rng.random() < 0.001 + 0.005 * i / READ:
            out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
    return "".join(out)


def pairs(acc, n, tag):
    seqs = [s for s in genes[acc] if len(s) >= 100]
    total = sum(len(s) for s in seqs)
    for k in range(n):
        x = rng.randrange(total)
        for s in seqs:
            if x < len(s):
                break
            x -= len(s)
        frag = min(max(int(rng.gauss(350, 50)), 200), 600, len(s))
        start = min(x, len(s) - frag)
        f = s[start:start + frag]
        if rng.random() < 0.5:
            f = f.translate(COMP)[::-1]
        mate = min(READ, len(f))
        r1, r2 = f[:mate], f[-mate:].translate(COMP)[::-1]
        name = f"{acc}__{tag}{k}"
        yield name, mutate(r1), mutate(r2)


os.makedirs(OUT, exist_ok=True)
focal = sorted(test)
map_rows = []
with open(f"{OUT}/truth.tsv", "w") as truth:
    truth.write("sample\tset\taccession\tspecies\tkind\tpairs\n")
    for kind in ("deep", "shallow"):
        for r in range(1, SAMPLES + 1):
            sample = f"{kind}_{r:02d}"
            members = []
            for i, sp in enumerate(focal):
                acc = test[sp].get(r) or test[sp][rng.choice(sorted(test[sp]))]
                if not genes.get(acc):
                    continue
                n = DEEP if kind == "deep" else SHALLOW[(i + r) % len(SHALLOW)]
                members.append((acc, sp, "strain", n))
            for acc in rng.sample([a for a in congener_reps if genes.get(a)], BACKGROUND):
                members.append((acc, species_of[acc], "congener", int(round(20 * 10 ** rng.random()))))
            f1, f2 = f"{OUT}/{sample}_1.fq.gz", f"{OUT}/{sample}_2.fq.gz"
            with gzip.open(f1, "wt", compresslevel=1) as o1, gzip.open(f2, "wt", compresslevel=1) as o2:
                for acc, sp, what, n in members:
                    truth.write(f"{sample}\t{kind}\t{acc}\t{sp}\t{what}\t{n}\n")
                    for name, r1, r2 in pairs(acc, n, ""):
                        o1.write(f"@{name}/1\n{r1}\n+\n{'F' * len(r1)}\n")
                        o2.write(f"@{name}/2\n{r2}\n+\n{'F' * len(r2)}\n")
            map_rows.append((sample, f1, f2))
with open(f"{OUT}/samples.tsv", "w") as o:
    o.writelines("\t".join(r) + "\n" for r in map_rows)
print(f"{len(map_rows)} samples written to {OUT}")
