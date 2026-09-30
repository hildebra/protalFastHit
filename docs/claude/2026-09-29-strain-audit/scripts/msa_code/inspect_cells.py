#!/usr/bin/env python3
"""List MSA cells of one species row that are N or IUPAC, with a strand-aware SAM pileup and the
genomes the reads come from (read names start with the source genome's accession).

usage: inspect_cells.py <protal_out> <manifest> <species> <sample> <kind: N|IUPAC> [max]
"""
import sys, os, re, gzip, collections, importlib.util
out_dir, manifest, species, sample, kind = sys.argv[1:6]
mx = int(sys.argv[6]) if len(sys.argv) > 6 else 20
sys.argv = [sys.argv[0], out_dir, manifest]
spec = importlib.util.spec_from_file_location("tc", os.path.join(os.path.dirname(os.path.abspath(__file__)), "tc_lib.py"))
tc = importlib.util.module_from_spec(spec); spec.loader.exec_module(tc)

msa = tc.read_fasta(os.path.join(out_dir, "strains", species + ".raw.msa.fna"))
refrow = msa[species + "_reference"]; row = msa[sample]
parts = tc.read_parts(os.path.join(out_dir, "strains", species + ".raw.partition.txt"))
acc = tc.truth[(sample, species)][0]
taxid = tc.tax_of[acc]
cells = []  # (gid, gene pos, cell, rep base, truth base, gene length)
for gid, a, b in parts:
    refseg, rowseg = refrow[a:b], row[a:b]
    repgene = refseg.replace("-", "")
    tg = tc.true_gene(acc, gid)
    p = 0
    for rc, sc in zip(refseg, rowseg):
        if rc == "-": continue
        if (kind == "N" and sc == "N") or (kind == "IUPAC" and sc in tc.IUPAC):
            t = tg[p] if tg is not None and len(tg) == len(repgene) else "?"
            cells.append((gid, p, sc, repgene[p], t, len(repgene)))
        p += 1
print(f"{len(cells)} {kind} cells in {species} {sample} ({acc})")
want = collections.defaultdict(set)
for c in cells[:mx]: want[c[0]].add(c[1])
names = {f"{taxid}_{g}": g for g in want}
pile = collections.defaultdict(lambda: collections.defaultdict(collections.Counter))
with gzip.open(os.path.join(out_dir, "alignments", sample + ".sam.gz"), "rt") as f:
    for line in f:
        if line.startswith("@"): continue
        t = line.split("\t")
        g = names.get(t[2])
        if g is None: continue
        flag = int(t[1]); strand = "-" if flag & 0x10 else "+"
        src = t[0].split("_contig")[0]
        src = "own" if src == acc else src
        pos = int(t[3]) - 1; seq = t[9]; q = 0
        for n, op in re.findall(r"(\d+)([MIDNSHPX=])", t[5]):
            n = int(n)
            if op in "MX=":
                for i in range(n):
                    if pos + i in want[g]: pile[g][pos + i][f"{seq[q + i]}{strand}:{src}"] += 1
                pos += n; q += n
            elif op in "IS":
                if op == "I" and pos in want[g]: pile[g][pos][f"ins{strand}:{src}"] += 1
                q += n
            elif op == "D":
                for i in range(n):
                    if pos + i in want[g]: pile[g][pos + i][f"del{strand}:{src}"] += 1
                pos += n
for gid, p, sc, r, t, L in cells[:mx]:
    print(f"gene {gid} pos {p}/{L}: cell {sc} rep {r} truth {t} pileup {dict(sorted(pile[gid][p].items()))}")
