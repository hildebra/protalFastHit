#!/usr/bin/env python3
"""foreign_diagnose.py OUT SAMPLE KIND - why genes of species in a foreign_world.sh sample are flagged foreign.

Profiles OUT/world/reads/<SAMPLE>_<KIND>.fq.gz again with --keep_foreign_genes, keeping the SAM, and for each foreign
gene of a species the sample holds lists the genes its reads put next to it (as MicrobialProfile::FinishLink does: a
long read's best records in read order, at most 3 kb apart on the read; the ends facing each other by orientation),
each partner's taxon, the pairing's verdict in the gene's taxon's clades (gene_neighbours.tsv, as Table::Assess), and
the species' own neighbours of that end in its genomes (gene_positions.tsv).
"""
import collections
import csv
import os
import re
import subprocess
import sys

out, sample, kind = sys.argv[1:4]
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.expanduser("~/protal-hap/src/scripts/mini_db"))
import gene_neighbours as gn  # noqa: E402

db = os.path.join(out, "full_species")
diag = os.path.join(out, "diag", f"{sample}_{kind}")
os.makedirs(diag, exist_ok=True)
sam = os.path.join(diag, f"{sample}.sam.zst")
if not os.path.exists(sam):
    subprocess.run([os.path.expanduser("~/protal-hap/build/protal"), "--db", db, "-1",
                    os.path.join(out, "world", "reads", f"{sample}_{kind}.fq.gz"), "--read_type", kind, "--prefix", sample,
                    "-o", diag, "-t", "6", "--no_strains", "--no_qcmsa", "--keep_foreign_genes"],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

nodes, reps = gn.read_taxonomy(os.path.join(db, "internal_taxonomy.dmp"))
counts, informative = collections.defaultdict(dict), {}
for line in open(os.path.join(db, "gene_neighbours.tsv")):
    if line.startswith("#") or line.startswith("clade"):
        continue
    c, g, e, p, pe, n, inf = (int(x) for x in line.split("\t")[:7])
    counts[(c, g, e)][(p, pe)] = [None] * n
    informative[(c, g, e)] = inf
clades = {c for c, _, _ in informative}


def verdict(taxid, gene, end, partner):
    chain, t, seen = [], taxid, set()
    while t in nodes and t not in seen:
        seen.add(t)
        if t in clades:
            chain.append(t)
        t = nodes[t][0]
    share, populated = gn.smoothed_share(counts, informative, chain, gene, end, partner)
    if share is None:
        return "unknown", None
    if not populated:
        return ("expected" if any(partner in counts.get((k, gene, end), {}) for k in chain) else "unknown"), share
    return ("expected" if share >= 0.2 else "unlikely" if share <= 0.05 else "rare"), share


# The species' own neighbours, from its genomes' positions.
own = collections.defaultdict(lambda: collections.defaultdict(collections.Counter))  # taxid -> (gene, end) -> partners
placements = collections.defaultdict(lambda: collections.defaultdict(list))
for line in open(os.path.join(db, "gene_positions.tsv")):
    if line.startswith("#") or line.startswith("accession"):
        continue
    f = line.rstrip("\n").split("\t")
    placements[int(f[1])][f[0]].append((int(f[5]), f[2], int(f[3]), f[4] == "1", int(f[6]) - 1, int(f[7]), f[8]))
foreign = []
for r in csv.DictReader(open(os.path.join(diag, f"{sample}.profile.genes.log")), delimiter="\t"):
    if r["Foreign"] == "1" and r["Predicted"] == "1":
        foreign.append((int(r["TaxID"]), int(r["GeneID"])))
for taxid in {t for t, _ in foreign}:
    for acc, pl in placements[taxid].items():
        for gene, end, partner, partner_end, gap in gn.neighbour_ends(pl, 3000):
            own[taxid][(gene, end)][(partner, partner_end)] += 1


def ops(cigar):
    return [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHPX=])", cigar)]


# A read's best records (primary and supplementary), with where they lie on the read.
reads = collections.defaultdict(list)
text = subprocess.run(["zstd", "-dcq", sam], check=True, capture_output=True, text=True).stdout
for line in text.splitlines():
    if line.startswith("@"):
        continue
    f = line.split("\t")
    flag = int(f[1])
    if flag & 0x104 or f[2] == "*":
        continue
    taxid, gene = (int(x) for x in f[2].split("_")[:2])
    o = ops(f[5])
    reverse = bool(flag & 0x10)
    lead = o if not reverse else list(reversed(o))
    start = 0
    for n, op in lead:
        if op in "HS":
            start += n
        else:
            break
    length = sum(n for n, op in o if op in "MIX=")
    reads[f[0]].append((start, start + length, taxid, gene, not reverse))
links = collections.defaultdict(collections.Counter)  # (taxid, gene) -> (partner taxid, gene, ends, verdict) -> reads
for name, recs in reads.items():
    recs.sort()
    for a, b in zip(recs, recs[1:]):
        if a[3] == b[3] or b[0] > a[1] + 3000:
            continue
        end_a, end_b = (3 if a[4] else 5), (5 if b[4] else 3)
        for x, y, ex, ey in ((a, b, end_a, end_b), (b, a, end_b, end_a)):
            if (x[2], x[3]) in foreign:
                v, share = verdict(x[2], x[3], ex, (y[3], ey))
                links[(x[2], x[3])][(y[2] == x[2], y[3], ex, ey, v)] += 1
print(f"{len(foreign)} foreign genes of called taxa in {sample} ({kind})")
for key in foreign[:12]:
    taxid, gene = key
    print(f"\ntaxon {taxid} ({nodes[taxid][2]}), gene {gene}; its own neighbours: "
          + "; ".join(f"{e}': {dict(own[taxid][(gene, e)])}" for e in (5, 3)))
    for (same, partner, ex, ey, v), n in links[key].most_common(6):
        print(f"  {n} reads: its {ex}' end next to gene {partner} {ey}' ({'same taxon' if same else 'another taxon'}): {v}")
