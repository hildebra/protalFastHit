#!/usr/bin/env python3
"""Consistency checks between protal's raw strain outputs for every species in a strains dir."""
import sys, os, glob, re, csv
from collections import defaultdict, Counter

def read_fasta(p):
    names, seqs, cur = [], [], []
    for line in open(p):
        line = line.rstrip("\n")
        if not line:
            continue
        if line[0] == ">":
            if names:
                seqs.append("".join(cur))
            names.append(line[1:]); cur = []
        else:
            cur.append(line)
    if names:
        seqs.append("".join(cur))
    return names, seqs

def parse_part(p):
    out = []
    for line in open(p):
        m = re.search(r"(\S+)\s*=\s*(\d+)-(\d+)", line)
        if m:
            out.append((m.group(1), int(m.group(2)), int(m.group(3))))
    return out

def main(d, min_cov=2):
    for msa in sorted(glob.glob(os.path.join(d, "*.raw.msa.fna"))):
        sp = os.path.basename(msa)[:-len(".raw.msa.fna")]
        names, seqs = read_fasta(msa)
        part = parse_part(os.path.join(d, sp + ".raw.partition.txt"))
        meta = list(csv.DictReader(open(os.path.join(d, sp + ".meta.tsv")), delimiter="\t"))
        stats = list(csv.DictReader(open(os.path.join(d, sp + ".snp_stats.tsv")), delimiter="\t"))
        L = set(len(s) for s in seqs)
        issues = []
        if len(L) != 1:
            issues.append(f"row lengths differ {L}")
        L = max(L)
        # partitions contiguous 1..L
        pos = 1
        for g, s, e in part:
            if s != pos:
                issues.append(f"partition {g} starts {s}, expected {pos}")
            pos = e + 1
        if pos - 1 != L:
            issues.append(f"partitions end {pos-1} != msa len {L}")
        genes = [g[4:] for g, _, _ in part]
        if len(set(genes)) != len(genes):
            issues.append("duplicate genes in partition")
        # names
        if len(set(names)) != len(names):
            issues.append("duplicate row names")
        meta_samples = []
        for r in meta:
            if r["sample"] not in meta_samples:
                meta_samples.append(r["sample"])
        stat_samples = [r["sample"] for r in stats]
        msa_samples = names[1:]
        # meta cell set
        cells = {(r["sample"], r["gene_id"]): r for r in meta}
        meta_genes = set(r["gene_id"] for r in meta)
        # per (sample, gene) non-gap count in MSA vs meta counts_vcov1/2
        mismatch_v2 = 0; mismatch_v1 = 0; cells_checked = 0; nometa_nongap = 0
        hcov_minus_msa = []
        for i, n in enumerate(names):
            if i == 0:
                continue
            for g, s, e in part:
                gid = g[4:]
                sub = seqs[i][s-1:e]
                nongap = sum(1 for c in sub if c not in "-N")
                nongapN = sum(1 for c in sub if c != "-")
                key = (n, gid)
                if key not in cells:
                    if nongapN:
                        nometa_nongap += 1
                    continue
                cells_checked += 1
                r = cells[key]
                glen = int(r["gene_length"])
                if nongapN != int(r["counts_vcov2"]):
                    mismatch_v2 += 1
                if nongapN != int(r["counts_vcov1"]):
                    mismatch_v1 += 1
                hcov_minus_msa.append(float(r["hcov"]) - nongapN / glen)
        # reference row
        ref = seqs[0]
        refgaps = sum(1 for c in ref if c == "-")
        iupac = Counter(c for s in seqs[1:] for c in s if c not in "ACGT-N")
        print(f"== {sp}: {len(names)} rows ({names[0]}), L={L}, {len(part)} partitions, "
              f"meta genes {len(meta_genes)} (partition genes {len(genes)}; in meta not in partition "
              f"{len(meta_genes - set(genes))}), meta samples {meta_samples}, stats samples {stat_samples}")
        print(f"   msa samples {msa_samples}; meta samples not in MSA: {[s for s in meta_samples if s not in msa_samples]}")
        print(f"   cells checked {cells_checked}; MSA non-gap != counts_vcov2: {mismatch_v2}; != counts_vcov1: {mismatch_v1}; "
              f"cells w/o meta but bases: {nometa_nongap}; ref gaps {refgaps}; IUPAC {dict(iupac)}")
        if hcov_minus_msa:
            hm = sorted(hcov_minus_msa)
            print(f"   hcov(meta) - nongap/len(MSA): mean {sum(hm)/len(hm):.4f} max {hm[-1]:.4f}")
        for x in issues:
            print("   ISSUE:", x)

if __name__ == "__main__":
    main(sys.argv[1])
