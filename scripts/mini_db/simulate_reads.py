#!/usr/bin/env python3
"""simulate_reads.py - paired-end reads from a mock community, without ART.

A small, dependency-free read simulator for testing protal on the mini DB:
fragments are drawn uniformly from the chosen genomes, both mates are written
in FR orientation with uniform substitution errors and constant quality.

  --genomes    genome table as written by simulate_gtdb_release.py
               (simulation/genomes.tsv: accession, gtdb_taxonomy, fasta_path, ...)
  --community  TSV "accession<TAB>relative_abundance" ('#' comments allowed).
               Abundances are relative cell abundances (what protal's profile
               estimates): a genome gets read pairs in proportion to
               abundance x genome length. They are normalised to sum to 1.

Writes <out_prefix>_R1.fq, <out_prefix>_R2.fq and <out_prefix>.truth.tsv
(accession, species, gtdb_taxonomy, relative_abundance, genome_length,
read_pairs). Read names are <accession>-<n>, so each read's source genome is
known. The output is fully determined by the inputs and --seed.

Usage:
  simulate_reads.py --genomes genomes.tsv --community community.tsv --out_prefix PREFIX
      [--pairs 30000] [--seed 7] [--read_length 150] [--insert_mean 350]
      [--insert_sd 50] [--error_rate 0.002]
"""

import argparse
import gzip
import math
import os
import random
import sys

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def read_table(path):
    """TSV with a header row -> list of dicts ('#' lines skipped)."""
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip() and not l.startswith("#")]
    header = lines[0].split("\t")
    return [dict(zip(header, l.split("\t"))) for l in lines[1:]]


def read_contigs(path):
    opener = gzip.open if path.endswith(".gz") else open
    contigs, chunks = [], None
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if chunks:
                    contigs.append("".join(chunks))
                chunks = []
            else:
                chunks.append(line.strip().upper())
    if chunks:
        contigs.append("".join(chunks))
    return contigs


def allocate(total, weights):
    """Split total into integers proportional to weights (largest remainder)."""
    s = sum(weights)
    exact = [total * w / s for w in weights]
    counts = [int(x) for x in exact]
    for i in sorted(range(len(exact)), key=lambda i: counts[i] - exact[i])[:total - sum(counts)]:
        counts[i] += 1
    return counts


def add_errors(rng, seq, rate):
    """Uniform substitutions at the given per-base rate (geometric skipping)."""
    if rate <= 0:
        return seq
    s = list(seq)
    log_q = math.log(1 - rate)
    pos = int(math.log(1 - rng.random()) / log_q)
    while pos < len(s):
        s[pos] = rng.choice("ACGT".replace(s[pos], "") if s[pos] in "ACGT" else "ACGT")
        pos += 1 + int(math.log(1 - rng.random()) / log_q)
    return "".join(s)


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--genomes", required=True, help="simulation/genomes.tsv of the mini release")
    ap.add_argument("--community", required=True, help="accession<TAB>relative_abundance")
    ap.add_argument("--out_prefix", required=True)
    ap.add_argument("--pairs", type=int, default=30000, help="total read pairs")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--read_length", type=int, default=150)
    ap.add_argument("--insert_mean", type=float, default=350)
    ap.add_argument("--insert_sd", type=float, default=50)
    ap.add_argument("--error_rate", type=float, default=0.002, help="substitutions per base")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    L = args.read_length
    genomes = {g["accession"]: g for g in read_table(args.genomes)}
    base = os.path.dirname(os.path.abspath(args.genomes))

    community = []
    for row in read_table(args.community):
        acc, abundance = row["accession"], float(row["relative_abundance"])
        if acc not in genomes:
            sys.exit(f"{args.community}: {acc} is not in {args.genomes}")
        if abundance <= 0:
            sys.exit(f"{args.community}: abundance of {acc} must be > 0")
        path = genomes[acc]["fasta_path"]
        if not os.path.exists(path):
            path = os.path.join(base, path)
        contigs = [c for c in read_contigs(path) if len(c) >= args.insert_mean + 3 * args.insert_sd]
        if not contigs:
            sys.exit(f"{path}: no contig is long enough for {args.insert_mean:.0f} bp inserts")
        community.append({"accession": acc, "abundance": abundance, "contigs": contigs,
                          "length": sum(len(c) for c in contigs), "row": genomes[acc]})
    if not community:
        sys.exit(f"{args.community} lists no genomes")

    total_abundance = sum(c["abundance"] for c in community)
    counts = allocate(args.pairs, [c["abundance"] * c["length"] for c in community])

    os.makedirs(os.path.dirname(os.path.abspath(args.out_prefix)), exist_ok=True)
    with open(f"{args.out_prefix}_R1.fq", "w", newline="\n") as r1, \
         open(f"{args.out_prefix}_R2.fq", "w", newline="\n") as r2, \
         open(f"{args.out_prefix}.truth.tsv", "w", newline="\n") as truth:
        truth.write("accession\tspecies\tgtdb_taxonomy\trelative_abundance\tgenome_length\tread_pairs\n")
        qual = "I" * L
        for c, n in zip(community, counts):
            lineage = c["row"]["gtdb_taxonomy"]
            truth.write(f"{c['accession']}\t{lineage.split(';')[-1]}\t{lineage}\t"
                        f"{c['abundance'] / total_abundance:.6f}\t{c['length']}\t{n}\n")
            weights = [len(s) for s in c["contigs"]]
            for i in range(1, n + 1):
                contig = rng.choices(c["contigs"], weights=weights)[0]
                insert = int(round(rng.gauss(args.insert_mean, args.insert_sd)))
                insert = max(L, min(insert, len(contig)))
                start = rng.randint(0, len(contig) - insert)
                fragment = contig[start:start + insert]
                if rng.random() < 0.5:
                    fragment = revcomp(fragment)
                name = f"{c['accession']}-{i}"
                r1.write(f"@{name}\n{add_errors(rng, fragment[:L], args.error_rate)}\n+\n{qual}\n")
                r2.write(f"@{name}\n{add_errors(rng, revcomp(fragment[-L:]), args.error_rate)}\n+\n{qual}\n")

    summary = ", ".join(f"{c['accession']}: {n}" for c, n in zip(community, counts))
    sys.stderr.write(f"Wrote {sum(counts)} read pairs to {args.out_prefix}_R[12].fq ({summary})\n")


if __name__ == "__main__":
    main()
