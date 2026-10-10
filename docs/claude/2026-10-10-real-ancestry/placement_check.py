"""Which copies does PlaceResidues (ColumnWeightsMsa.h) leave out? Per gene, each representative's protein (its
reference.fna copy translated in frame) against its row of GTDB's alignment (column_msa/): the share of the row's
residues found in order (difflib's matching blocks, the same question PlaceResidues answers), and what differs.

usage: placement_check.py <converted database dir> <unpacked dir (reference.fna)>
"""
import collections
import difflib
import glob
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "scripts"))
import ancestry_oracle as oracle  # noqa: E402

DB, UNPACKED = sys.argv[1], sys.argv[2]


def translate(nt):
    return "".join(oracle.amino_acid(oracle.CODE[nt[i]], oracle.CODE[nt[i + 1]], oracle.CODE[nt[i + 2]])
                   for i in range(0, len(nt) - 2, 3))


copies = {}
with open(f"{UNPACKED}/reference.fna") as f:
    name = None
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            name = tuple(int(x) for x in line[1:].split("_")[:2])
        elif name:
            copies[name] = copies.get(name, "") + line
marker_of = {}
with open(f"{DB}/gene2geneid.tsv") as f:
    for line in f:
        m, g = line.split()[:2]
        marker_of[int(g)] = m
per_gene = collections.defaultdict(list)
for path in sorted(glob.glob(f"{DB}/column_msa/*.faa*")):
    gene = int(os.path.basename(path).split(".")[0])
    text = subprocess.run(["zstd", "-dc", path], capture_output=True, text=True).stdout if path.endswith(".zst") else open(path).read()
    taxid = None
    rows = {}
    for line in text.splitlines():
        if line.startswith(">"):
            taxid = int(line[1:])
            rows[taxid] = ""
        else:
            rows[taxid] += line.strip()
    for taxid, row in rows.items():
        protein = translate(copies.get((taxid, gene), ""))
        residues = row.replace("-", "").replace(".", "").upper()
        if not residues:
            continue
        sm = difflib.SequenceMatcher(None, residues, protein, autojunk=False)
        found = sum(b.size for b in sm.get_matching_blocks())
        per_gene[gene].append((taxid, found / len(residues), len(residues), len(protein), residues[:12], protein[:12]))
for gene in sorted(per_gene):
    rs = per_gene[gene]
    low = [r for r in rs if r[1] < 0.5]
    if low:
        print(f"gene {gene} {marker_of.get(gene)}: {len(low)} of {len(rs)} below 0.5; e.g. taxid {low[0][0]} share "
              f"{low[0][1]:.2f}, row {low[0][2]} residues '{low[0][4]}...', protein {low[0][3]} '{low[0][5]}...'")
print(f"copies below 0.5: {sum(r[1] < 0.5 for rs in per_gene.values() for r in rs)} of {sum(len(rs) for rs in per_gene.values())}")
