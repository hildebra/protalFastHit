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



def read_parts(path):
    parts = []
    for line in open(path):
        name, rng = line.split("=")
        a, b = rng.strip().split("-")
        parts.append((int(name.split("gene")[1]), int(a) - 1, int(b)))
    return parts
