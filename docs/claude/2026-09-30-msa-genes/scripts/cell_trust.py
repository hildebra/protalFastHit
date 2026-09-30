#!/usr/bin/env python3
"""Are the MSA cells of genes with few unique k-mers as trustworthy as the others?

usage: cell_trust.py STRAINS_DIR DB_DIR MANIFEST
Per (sample, gene) cell of every raw MSA: the rate of bases that differ from the reference row, and
of IUPAC codes, over the positions the row writes. In the simulations a strain's substitutions are
spread evenly over its genes, so a cell's rate should be near its row's median; reads of another
species would raise it. Cells are grouped by the gene's unique k-mers (short + long) and by whether
a congener (another species of the genus) is in the sample.
"""
import csv
import os
import re
import sys
from collections import defaultdict

strains, db, manifest = sys.argv[1:4]
BASES = set("ACGT")
IUPAC = set("RYSWKMBDHV")

uniq = {}
with open(os.path.join(db, "unique_kmers.tsv")) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        uniq[(int(f[0]), int(f[1]))] = int(f[2]) + int(f[4])
genera = defaultdict(set)  # sample -> genera with 2+ species present
species_in = defaultdict(set)
with open(manifest) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        species_in[r["sample"]].add(r["species"])
for s, sps in species_in.items():
    count = defaultdict(int)
    for sp in sps:
        count[sp.split()[0]] += 1
    genera[s] = {g for g, c in count.items() if c > 1}


def read_fasta(path):
    names, seqs, cur = [], [], []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if names:
                seqs.append("".join(cur))
            names.append(line[1:])
            cur = []
        else:
            cur.append(line)
    if names:
        seqs.append("".join(cur))
    return names, seqs


cells = []  # (bin, congener, rate, ratio, iupac_rate, written)
with open(os.path.join(strains, "species.tsv")) as fh:
    listed = [(r["species"], int(r["taxid"])) for r in csv.DictReader(fh, delimiter="\t") if r["raw_msa"] != "-"]
for sp, taxid in listed:
    genus = sp[3:].split("_")[0]
    names, seqs = read_fasta(os.path.join(strains, sp + ".raw.msa.fna"))
    ref = seqs[0]
    parts = []
    for line in open(os.path.join(strains, sp + ".raw.partition.txt")):
        m = re.match(r"DNA, gene(\d+) = (\d+)-(\d+)", line)
        if m:
            parts.append((int(m.group(1)), int(m.group(2)) - 1, int(m.group(3))))
    for name, seq in zip(names[1:], seqs[1:]):
        row = []
        for gene, s, e in parts:
            written = diff = iupac = 0
            for c in range(s, e):
                r, b = ref[c], seq[c]
                if r not in BASES:
                    continue
                if b in BASES:
                    written += 1
                    diff += b != r
                elif b in IUPAC:
                    written += 1
                    iupac += 1
            if written >= 200:
                row.append((gene, diff / written, iupac / written, written))
        if len(row) < 10:
            continue
        rates = sorted(x[1] for x in row)
        med = rates[len(rates) // 2]
        congener = genus in genera.get(name, set())
        for gene, rate, irate, written in row:
            u = uniq.get((taxid, gene), 0)
            b = "0" if u == 0 else "1-9" if u < 10 else "10-49" if u < 50 else "50-99" if u < 100 else "100+"
            # a cell is suspect when it holds 3 more differences than the row's median rate predicts,
            # and more than twice that rate
            expected = med * written
            suspect = rate * written > expected + 3 and rate > 2 * med
            cells.append((b, congener, rate, med, irate, suspect))

print("gene's unique k-mers, congener in the sample: cells, mean rate / row median, suspect cells, IUPAC per kb")
for b in ("0", "1-9", "10-49", "50-99", "100+"):
    for congener in (False, True):
        sel = [c for c in cells if c[0] == b and c[1] == congener]
        if not sel:
            continue
        mean_rate = sum(c[2] for c in sel) / len(sel)
        mean_med = sum(c[3] for c in sel) / len(sel)
        susp = sum(c[5] for c in sel)
        iup = 1000 * sum(c[4] for c in sel) / len(sel)
        print(f"  {b:>6} {'congener' if congener else 'alone   '}: {len(sel):6d} cells, rate {mean_rate:.4f} vs median {mean_med:.4f}"
              f" ({mean_rate / mean_med if mean_med else float('nan'):.2f}x), suspect {susp} ({susp / len(sel):.2%}), IUPAC {iup:.3f}/kb")
