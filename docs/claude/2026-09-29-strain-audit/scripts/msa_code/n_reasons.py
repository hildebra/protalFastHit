#!/usr/bin/env python3
"""Why is a cell N? For every N cell of every single-strain row: the pileup's alleles and strands.
Also: of reference-written cells (no variant) how many are covered by one strand only.

usage: n_reasons.py <protal_out> <manifest>
"""
import sys, os, re, gzip, collections, importlib.util
out_dir, manifest = sys.argv[1:3]
spec = importlib.util.spec_from_file_location("tc", os.path.join(os.path.dirname(os.path.abspath(__file__)), "tc_lib.py"))
tc = importlib.util.module_from_spec(spec); spec.loader.exec_module(tc)

strain_dir = os.path.join(out_dir, "strains")
# sample -> taxid -> gene -> list of (pos, kind)
targets = collections.defaultdict(lambda: collections.defaultdict(dict))
cellinfo = {}
for fn in sorted(os.listdir(strain_dir)):
    if not fn.endswith(".raw.msa.fna"): continue
    species = fn[:-len(".raw.msa.fna")]
    msa = tc.read_fasta(os.path.join(strain_dir, fn))
    refrow = msa[species + "_reference"]
    parts = tc.read_parts(os.path.join(strain_dir, species + ".raw.partition.txt"))
    for sample, row in msa.items():
        g = tc.truth.get((sample, species), [])
        if len(g) != 1: continue
        acc = g[0]; taxid = tc.tax_of[acc]
        for gid, a, b in parts:
            refseg, rowseg = refrow[a:b], row[a:b]
            repgene = refseg.replace("-", "")
            tg = tc.true_gene(acc, gid)
            ok = tg is not None and len(tg) == len(repgene)
            p = 0
            for rc, sc in zip(refseg, rowseg):
                if rc == "-": continue
                if sc == "N":
                    targets[sample][taxid].setdefault(gid, set()).add(p)
                    cellinfo[(sample, taxid, gid, p)] = (sc, repgene[p], tg[p] if ok else "?", acc)
                p += 1

reasons = collections.Counter()
examples = collections.defaultdict(list)
for sample in sorted(targets):
    names = {f"{t}_{g}": (t, g) for t in targets[sample] for g in targets[sample][t]}
    pile = collections.defaultdict(collections.Counter)
    with gzip.open(os.path.join(out_dir, "alignments", sample + ".sam.gz"), "rt") as f:
        for line in f:
            if line.startswith("@"): continue
            t = line.split("\t")
            tg = names.get(t[2])
            if tg is None: continue
            want = targets[sample][tg[0]][tg[1]]
            flag = int(t[1]); strand = "-" if flag & 0x10 else "+"
            if flag & 0x100 or int(t[4]) < 4: continue
            pos = int(t[3]) - 1; seq = t[9]; q = 0
            for n, op in re.findall(r"(\d+)([MIDNSHPX=])", t[5]):
                n = int(n)
                if op in "MX=":
                    for i in range(n):
                        if pos + i in want: pile[(tg[0], tg[1], pos + i)][(seq[q + i], strand)] += 1
                    pos += n; q += n
                elif op in "IS": q += n
                elif op == "D":
                    for i in range(n):
                        if pos + i in want: pile[(tg[0], tg[1], pos + i)][("-", strand)] += 1
                    pos += n
    for (taxid, gid, p), c in pile.items():
        sc, rep, truth, acc = cellinfo[(sample, taxid, gid, p)]
        bases = collections.Counter()
        strands = collections.defaultdict(set)
        for (b, s), n in c.items():
            bases[b] += n; strands[b].add(s)
        top = bases.most_common(1)[0][0]
        if len(bases) == 1 and top != rep and len(strands[top]) == 1:
            r = "only alt, one strand"
        elif top != rep and len(strands[top]) == 1:
            r = "alt top, one strand, + others"
        elif sum(bases.values()) < 2:
            r = "depth<2"
        else:
            r = "other"
        r += " (true SNP)" if truth not in ("?", rep) else ""
        reasons[r] += 1
        if len(examples[r]) < 6: examples[r].append((sample, taxid, gid, p, rep, truth, dict(c)))
for r, n in reasons.most_common():
    print(f"{n}\t{r}")
    for e in examples[r]: print("\t", e)
