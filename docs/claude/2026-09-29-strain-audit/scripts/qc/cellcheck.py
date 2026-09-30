#!/usr/bin/env python3
"""Per (sample, gene) cell: meta hcov (depth >= 1) vs the fraction of the gene the raw MSA holds
(non-gap, depth >= snp_min_cov), and meta 'filtered' vs the 'N's the MSA writes."""
import sys, os, csv, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from consistency import read_fasta, parse_part

tot = pass_meta = fail_msa = 0
nf_eq = nf_cells = 0
depth_gate_fail = 0
for msa in sorted(glob.glob(os.path.join(sys.argv[1], "*.raw.msa.fna"))):
    sp = os.path.basename(msa)[:-len(".raw.msa.fna")]
    names, seqs = read_fasta(msa)
    row = dict(zip(names, seqs))
    part = {g[4:]: (s, e) for g, s, e in parse_part(os.path.join(sys.argv[1], sp + ".raw.partition.txt"))}
    for r in csv.DictReader(open(os.path.join(sys.argv[1], sp + ".meta.tsv")), delimiter="\t"):
        if r["sample"] not in row or r["gene_id"] not in part:
            continue
        s, e = part[r["gene_id"]]
        sub = row[r["sample"]][s - 1:e]
        glen = int(r["gene_length"])
        msa_frac = sum(1 for c in sub if c not in "-") / glen
        tot += 1
        if float(r["mean_vcov_nonzero"]) < 1.0:
            depth_gate_fail += 1
        if float(r["hcov"]) >= 0.3:
            pass_meta += 1
            if msa_frac < 0.3:
                fail_msa += 1
        nN = sub.count("N")
        nf_cells += 1
        nf_eq += nN == int(r["filtered"])
print(f"cells {tot}; mean_vcov_nonzero < 1.0: {depth_gate_fail}; hcov>=0.3 per meta: {pass_meta}, of which the MSA holds <30% of the gene: {fail_msa}")
print(f"cells whose MSA 'N' count equals meta 'filtered': {nf_eq} of {nf_cells}")
