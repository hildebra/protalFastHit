#!/usr/bin/env python3
"""Per-gene truth breakdown for one species row, and a SAM pileup at positions written as reference
where the truth has a SNP.

usage: inspect_gene.py <protal_out> <manifest> <species s__X_y> <sample> [max_genes]
"""
import sys, os, re, gzip, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
out_dir, manifest, species, sample = sys.argv[1:5]
max_genes = int(sys.argv[5]) if len(sys.argv) > 5 else 3
sys.argv = [sys.argv[0], out_dir, manifest]
import importlib.util
spec = importlib.util.spec_from_file_location("tc", os.path.join(os.path.dirname(os.path.abspath(__file__)), "tc_lib.py"))
tc = importlib.util.module_from_spec(spec); spec.loader.exec_module(tc)

msa = tc.read_fasta(os.path.join(out_dir, "strains", species + ".raw.msa.fna"))
refrow = msa[species + "_reference"]; row = msa[sample]
parts = tc.read_parts(os.path.join(out_dir, "strains", species + ".raw.partition.txt"))
acc = tc.truth[(sample, species)][0]
taxid = tc.tax_of[acc]
bad = {}
for gid, a, b in parts:
    refseg, rowseg = refrow[a:b], row[a:b]
    repgene = refseg.replace("-", "")
    tg = tc.true_gene(acc, gid)
    if tg is None or len(tg) != len(repgene): continue
    p = 0; stats = collections.Counter(); positions = []
    for rc, sc in zip(refseg, rowseg):
        if rc == "-": continue
        t, r = tg[p], repgene[p]
        if t != r:
            stats["snp"] += 1
            if sc == r: stats["ref"] += 1; positions.append(p)
            elif sc == t: stats["ok"] += 1
        p += 1
    if stats["ref"]: bad[gid] = (stats, positions, repgene, tg)
print(f"{species} {sample} {acc} taxid {taxid}: genes with SNPs written as ref: {len(bad)}")
for gid, (s, pos, rg, tg) in sorted(bad.items(), key=lambda kv: -kv[1][0]["ref"])[:10]:
    print(f"  gene {gid}: snps {s['snp']} ok {s['ok']} ref {s['ref']} first positions {pos[:8]}")

sam = os.path.join(out_dir, "alignments", sample + ".sam.gz")
want = {f"{taxid}_{gid}": gid for gid in list(sorted(bad, key=lambda g: -bad[g][0]['ref']))[:max_genes]}
pile = collections.defaultdict(lambda: collections.defaultdict(collections.Counter))
origin = collections.defaultdict(collections.Counter)
with gzip.open(sam, "rt") as f:
    for line in f:
        if line.startswith("@"): continue
        t = line.split("\t")
        if t[2] not in want: continue
        flag = int(t[1])
        if flag & 0x4: continue
        gid = want[t[2]]
        origin[gid][t[0].rsplit("-", 1)[0].rsplit("_contig", 1)[0]] += 1
        pos = int(t[3]) - 1; seq = t[9]; q = 0
        for n, op in re.findall(r"(\d+)([MIDNSHPX=])", t[5]):
            n = int(n)
            if op in "MX=":
                for i in range(n): pile[gid][pos + i][seq[q + i]] += 1
                pos += n; q += n
            elif op in "IS": q += n
            elif op == "D":
                for i in range(n): pile[gid][pos + i]["-"] += 1
                pos += n
for gid in want.values():
    s, pos, rg, tg = bad[gid]
    print(f"gene {gid}: read origins {dict(origin[gid])}")
    for p in pos[:6]:
        print(f"   pos {p}: rep {rg[p]} truth {tg[p]} pileup {dict(pile[gid][p])}")
