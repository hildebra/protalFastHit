#!/usr/bin/env python3
"""What a protal SAM says about some taxa (e.g. a species and its congeners): per taxon, the reads, genes, identity,
MAPQ and the tags of their records (ZR: a long read's gene that its consensus taxon changed (1) or that another taxon
fits clearly better (2); ZA: the other taxa within 5 edits); per gene, each taxon's records and identity; and the reads
that have records of more than one of the taxa, gene by gene.

Identity = matching bases / (matching + mismatching + inserted + deleted), from the CIGAR (protal writes M for
matches and X for mismatches). Records: the primary (flag 0/16), a long read's other genes (2048), other candidates
(256, counted apart).

usage: sam_taxa.py SAM(.zst|.gz) GENOME2TIID TAXID[,TAXID...] [--reads N] > report.txt
"""
import argparse
import collections
import gzip
import re
import statistics
import subprocess
import sys


def open_sam(path):
    if path.endswith(".zst"):
        return subprocess.Popen(["zstd", "-dc", path], stdout=subprocess.PIPE, text=True).stdout
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def identity(cigar):
    counts = collections.Counter()
    for n, op in re.findall(r"(\d+)([MIDNSHPX=])", cigar):
        counts[op] += int(n)
    aligned = counts["M"] + counts["="] + counts["X"] + counts["I"] + counts["D"]
    return (counts["M"] + counts["="]) / aligned if aligned else float("nan"), counts["M"] + counts["="] + counts["X"]


def species_names(path):
    names = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4:
                names[f[1]] = f[3].split(";")[-1]
    return names


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sam")
    ap.add_argument("genome2tiid")
    ap.add_argument("taxids")
    ap.add_argument("--reads", type=int, default=15, help="mixed reads to list (default 15)")
    args = ap.parse_args()
    focus = args.taxids.split(",")
    names = species_names(args.genome2tiid)

    per_taxon = {t: collections.defaultdict(list) for t in focus}  # taxid -> field -> values
    per_gene = collections.defaultdict(lambda: collections.defaultdict(list))  # gene -> taxid -> identities
    reads = collections.defaultdict(list)  # read -> [(taxid, gene, identity, mapq, ZR, ZA, flag)]
    secondary = collections.Counter()
    for line in open_sam(args.sam):
        if line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        flag = int(f[1])
        if flag & 4:
            continue
        taxid, gene = f[2].split("_", 1)
        if taxid not in per_taxon:
            continue
        if flag & 256:
            secondary[taxid] += 1
            continue
        tags = dict(t.split(":", 2)[0::2] for t in f[11:])
        ident, aligned = identity(f[5])
        zr = tags.get("ZR", "0")
        za = tags.get("ZA", "*")
        p = per_taxon[taxid]
        p["identity"].append(ident)
        p["mapq"].append(int(f[4]))
        p["zr"].append(zr)
        p["za_other"].append(za != "*")
        p["genes"].append(gene)
        p["aligned"].append(aligned)
        per_gene[gene][taxid].append(ident)
        reads[f[0]].append((taxid, gene, ident, int(f[4]), zr, za, flag))

    print(f"# {args.sam}")
    print("taxid\tspecies\trecords\treads\tgenes\tmean_identity\tidentity_p10\tmean_mapq\tmapq0\tZR1\tZR2\tZA_other\t"
          "secondary")
    for t in focus:
        p = per_taxon[t]
        if not p["identity"]:
            print(f"{t}\t{names.get(t, '?')}\t0")
            continue
        ids = sorted(p["identity"])
        read_count = sum(1 for r in reads.values() if any(x[0] == t for x in r))
        print(f"{t}\t{names.get(t, '?')}\t{len(ids)}\t{read_count}\t{len(set(p['genes']))}\t"
              f"{statistics.fmean(ids):.4f}\t{ids[len(ids) // 10]:.4f}\t{statistics.fmean(p['mapq']):.1f}\t"
              f"{sum(m == 0 for m in p['mapq'])}\t{p['zr'].count('1')}\t{p['zr'].count('2')}\t"
              f"{sum(p['za_other'])}\t{secondary[t]}")

    print("\n# genes with records of more than one of the taxa: records (mean identity) per taxon")
    print("gene\t" + "\t".join(f"{t} {names.get(t, '?')}" for t in focus))
    for gene in sorted(per_gene, key=lambda g: int(g) if g.isdigit() else g):
        by = per_gene[gene]
        if sum(1 for t in focus if by.get(t)) > 1:
            cells = [f"{len(by[t])} ({statistics.fmean(by[t]):.4f})" if by.get(t) else "-" for t in focus]
            print(f"{gene}\t" + "\t".join(cells))

    print("\n# genes of the first taxon by mean identity (lowest first): records, mean identity")
    first = focus[0]
    genes = {g: per_gene[g][first] for g in per_gene if per_gene[g].get(first)}
    for gene, ids in sorted(genes.items(), key=lambda kv: statistics.fmean(kv[1]))[:15]:
        print(f"{gene}\t{len(ids)}\t{statistics.fmean(ids):.4f}")

    mixed = [(r, recs) for r, recs in reads.items() if len({x[0] for x in recs}) > 1]
    only_other = [(r, recs) for r, recs in reads.items() if first not in {x[0] for x in recs}]
    print(f"\n# reads with records of {first} only: {sum(1 for r in reads.values() if {x[0] for x in r} == {first})}; "
          f"of other taxa only: {len(only_other)}; mixed: {len(mixed)}")
    print("# mixed and other-only reads, gene by gene: taxid:gene identity MAPQ ZR ZA")
    for r, recs in (mixed + only_other)[:args.reads]:
        print(r + "\t" + "  ".join(f"{t}:{g} {i:.3f} q{q} zr{zr} za={za}" for t, g, i, q, zr, za, _ in recs))


if __name__ == "__main__":
    main()
