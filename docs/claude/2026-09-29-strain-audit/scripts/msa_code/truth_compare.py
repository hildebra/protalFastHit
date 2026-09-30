#!/usr/bin/env python3
"""Compare protal raw MSA cells to the true gene sequences of the simulated genomes.

usage: truth_compare.py <protal_out_dir> <manifest.tsv> [world_dir]
Only genes whose true sequence has the representative's length (substitutions only) are compared.
"""
import gzip, os, sys, collections

out_dir, manifest = sys.argv[1], sys.argv[2]
world = sys.argv[3] if len(sys.argv) > 3 else os.path.expanduser("~/audit5/world")
sim = os.path.join(world, "gtdb_r226/simulation")
IUPAC = {'R': 'AG', 'Y': 'CT', 'W': 'AT', 'S': 'CG', 'M': 'AC', 'K': 'GT', 'B': 'CGT', 'H': 'ACT', 'D': 'AGT', 'V': 'ACG'}
COMP = str.maketrans("ACGTN", "TGCAN")


def read_fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    seqs, name, buf = {}, None, []
    with op(path, "rt") as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None: seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.upper())
    if name is not None: seqs[name] = "".join(buf)
    return seqs


genome_fasta = {}
for i, line in enumerate(open(os.path.join(sim, "genomes.tsv"))):
    if i == 0: continue
    t = line.rstrip("\n").split("\t")
    genome_fasta[t[0]] = t[2]
marker_pos = {}
for i, line in enumerate(open(os.path.join(sim, "marker_positions.tsv"))):
    if i == 0: continue
    acc, marker, contig, start, end, strand = line.rstrip("\n").split("\t")
    marker_pos.setdefault((acc, marker), []).append((contig, int(start), int(end), strand))
geneid_of = {}
for line in open(os.path.join(world, "protal_db/gene2geneid.tsv")):
    m, g = line.split()
    geneid_of[m] = int(g)
marker_of = {v: k for k, v in geneid_of.items()}
tax_of, rep_of = {}, {}
for line in open(os.path.join(world, "protal_db/genome2tiid.tsv")):
    t = line.rstrip("\n").split("\t")
    tax_of[t[0]] = int(t[1]); rep_of[t[0]] = t[2]
    tax_of[t[2]] = int(t[1]); rep_of[t[2]] = t[2]
truth = collections.defaultdict(list)  # (sample, species name) -> genomes
for i, line in enumerate(open(manifest)):
    t = line.rstrip("\n").split("\t")
    if i == 0: hdr = t; continue
    r = dict(zip(hdr, t))
    truth[(r["sample"], "s__" + r["species"].replace(" ", "_"))].append(r["genome"])

genome_cache = {}
def true_gene(acc, geneid):
    locs = marker_pos.get((acc, marker_of[geneid]))
    if not locs or len(locs) != 1: return None
    contig, start, end, strand = locs[0]
    if acc not in genome_cache: genome_cache[acc] = read_fasta(genome_fasta[acc])
    s = genome_cache[acc][contig][start - 1:end]
    if strand == "-": s = s.translate(COMP)[::-1]
    return s


totals = collections.Counter()
per_row = []
strain_dir = os.path.join(out_dir, "strains")
for fn in sorted(os.listdir(strain_dir)):
    if not fn.endswith(".raw.msa.fna"): continue
    species = fn[:-len(".raw.msa.fna")]
    msa = read_fasta(os.path.join(strain_dir, fn))
    ref_name = species + "_reference"
    if ref_name not in msa: continue
    refrow = msa[ref_name]
    parts = []
    for line in open(os.path.join(strain_dir, species + ".raw.partition.txt")):
        name, rng = line.split("=")
        gid = int(name.split("gene")[1])
        a, b = rng.strip().split("-")
        parts.append((gid, int(a) - 1, int(b)))
    assert parts[-1][2] == len(refrow), (species, parts[-1], len(refrow))
    for i in range(1, len(parts)):
        assert parts[i][1] == parts[i - 1][2], (species, parts[i - 1], parts[i])
    for sample, row in msa.items():
        if sample == ref_name: continue
        assert len(row) == len(refrow), (species, sample)
        genomes = truth.get((sample, species), [])
        if len(genomes) != 1: continue
        acc = genomes[0]
        rep = rep_of[acc]
        c = collections.Counter()
        for gid, a, b in parts:
            refseg, rowseg = refrow[a:b], row[a:b]
            repgene = refseg.replace("-", "")
            rg = true_gene(rep, gid)
            if rg is not None and rg != repgene:
                c["rep_gene_mismatch"] += 1
                continue
            tg = true_gene(acc, gid)
            if tg is None or len(tg) != len(repgene):
                c["genes_skipped_indel"] += 1
                continue
            mm = [x != y for x, y in zip(tg, repgene)]
            if any(sum(mm[i:i + 10]) >= 5 for i in range(len(mm))):
                c["genes_skipped_shifted"] += 1
                continue
            c["genes_compared"] += 1
            p = 0
            for rc, sc in zip(refseg, rowseg):
                if rc == "-":
                    c["ins_col_nongap" if sc != "-" else "ins_col_gap"] += 1
                    continue
                t = tg[p]; r = repgene[p]; p += 1
                snp = t != r
                if sc == "-": k = "gap"
                elif sc == "N": k = "N"
                elif sc in IUPAC: k = "iupac_with_truth" if t in IUPAC[sc] else "iupac_without_truth"
                elif sc == t: k = "correct"
                elif sc == r and snp: k = "ref_at_true_snp"
                else: k = "false_snp"
                c[k] += 1
                if snp: c["true_snp_" + k] += 1
        per_row.append((species, sample, acc, c))
        totals.update(c)

for species, sample, acc, c in per_row:
    called = c["correct"] + c["ref_at_true_snp"] + c["false_snp"]
    print(f"{species}\t{sample}\t{acc}\tgenes={c['genes_compared']}\tcorrect={c['correct']}\tref_at_snp={c['ref_at_true_snp']}"
          f"\tfalse_snp={c['false_snp']}\tiupac={c['iupac_with_truth']}+{c['iupac_without_truth']}\tN={c['N']}\tgap={c['gap']}"
          f"\ttrue_snps_called={c['true_snp_correct']}/{sum(v for k, v in c.items() if k.startswith('true_snp_'))}"
          f"\tins_nongap={c['ins_col_nongap']}")
print("TOTAL", dict(sorted(totals.items())))
