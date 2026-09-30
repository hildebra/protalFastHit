#!/usr/bin/env python3
"""SAM edits for end-to-end probes of protal's MSA.

  sam_edit.py xm   <in.sam.gz> <out.sam>                 X ops -> M (as a standard aligner writes)
  sam_edit.py n    <in.sam.gz> <out.sam> <gene> <pos0>   base at gene position pos0 -> N in every read
  sam_edit.py snps <protal_out> <manifest> <species> <sample> [n]   true SNPs the row calls correctly
"""
import sys, re, gzip, os, collections, importlib.util

def records(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        for line in f:
            yield line

mode = sys.argv[1]
if mode == "xm":
    with open(sys.argv[3], "w") as out:
        for line in records(sys.argv[2]):
            if line.startswith("@"): out.write(line); continue
            t = line.split("\t")
            ops = [(int(n), "M" if op == "X" else op) for n, op in re.findall(r"(\d+)([MIDNSHPX=])", t[5])]
            merged = []
            for n, op in ops:
                if merged and merged[-1][1] == op: merged[-1] = (merged[-1][0] + n, op)
                else: merged.append((n, op))
            t[5] = "".join(f"{n}{op}" for n, op in merged)
            out.write("\t".join(t))
elif mode == "n":
    gene, target = sys.argv[4], int(sys.argv[5])
    changed = 0
    with open(sys.argv[3], "w") as out:
        for line in records(sys.argv[2]):
            if line.startswith("@"): out.write(line); continue
            t = line.split("\t")
            if t[2] == gene:
                pos, q, seq = int(t[3]) - 1, 0, list(t[9])
                for n, op in re.findall(r"(\d+)([MIDNSHPX=])", t[5]):
                    n = int(n)
                    if op in "MX=":
                        if pos <= target < pos + n:
                            seq[q + target - pos] = "N"; changed += 1
                        pos += n; q += n
                    elif op in "IS": q += n
                    elif op == "D": pos += n
                t[9] = "".join(seq)
            out.write("\t".join(t))
    print(f"set {changed} read bases to N", file=sys.stderr)
elif mode == "snps":
    out_dir, manifest, species, sample = sys.argv[2:6]
    k = int(sys.argv[6]) if len(sys.argv) > 6 else 5
    sys.argv = [sys.argv[0], out_dir, manifest]
    spec = importlib.util.spec_from_file_location("tc", os.path.join(os.path.dirname(os.path.abspath(__file__)), "tc_lib.py"))
    tc = importlib.util.module_from_spec(spec); spec.loader.exec_module(tc)
    msa = tc.read_fasta(os.path.join(out_dir, "strains", species + ".raw.msa.fna"))
    refrow = msa[species + "_reference"]; row = msa[sample]
    acc = tc.truth[(sample, species)][0]; taxid = tc.tax_of[acc]
    found = 0
    for gid, a, b in tc.read_parts(os.path.join(out_dir, "strains", species + ".raw.partition.txt")):
        refseg, rowseg = refrow[a:b], row[a:b]
        repgene = refseg.replace("-", "")
        tg = tc.true_gene(acc, gid)
        if tg is None or len(tg) != len(repgene) or "-" in refseg: continue
        for p in range(20, len(repgene) - 20):
            if tg[p] != repgene[p] and rowseg[p] == tg[p]:
                print(f"{taxid}_{gid}\t{p}\trep {repgene[p]}\ttruth {tg[p]}\tcolumn {a + p + 1}")
                found += 1
                break
        if found >= k: break
