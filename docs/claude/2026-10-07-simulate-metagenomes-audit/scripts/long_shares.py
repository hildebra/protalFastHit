#!/usr/bin/env python3
"""Each genome's share of a long-read sample's bases against its weight's share, by how fragmented the genome is:
the reads of simulate_metagenomes --long_samples (names g<i>x_<n>), the genomes TSV it read (sample, genome, fasta,
weight, host), and each FASTA's base-weighted mean contig length.

  long_shares.py READS.fq.zst GENOMES.tsv SAMPLE

Prints, by thirds of the genomes' base-weighted contig length, the mean of (bases share / weight share), and the
spread of that ratio over all genomes.
"""

import gzip
import statistics
import subprocess
import sys


def contig_lengths(path):
    lengths, current = [], 0
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                current = 0
            else:
                current += len(line.strip())
    if current:
        lengths.append(current)
    return lengths


def main():
    reads, genomes_tsv, sample = sys.argv[1:4]
    genomes = []
    with open(genomes_tsv) as fh:
        next(fh)
        for line in fh:
            s, name, fasta, weight, host = line.rstrip("\n").split("\t")
            if s == sample:
                genomes.append((name, fasta, float(weight)))
    bases = [0] * len(genomes)
    text = subprocess.run(["zstd", "-dc", reads], check=True, stdout=subprocess.PIPE, text=True).stdout.split("\n")
    for i in range(0, len(text) - 3, 4):
        g = int(text[i][2:text[i].index("x_")])
        bases[g] += len(text[i + 1])
    total_bases, total_weight = sum(bases), sum(w for _, _, w in genomes)
    rows = []
    for (name, fasta, weight), b in zip(genomes, bases):
        lengths = contig_lengths(fasta)
        weighted = sum(x * x for x in lengths) / sum(lengths)  # the contig length a base is in, on average
        rows.append((weighted, (b / total_bases) / (weight / total_weight), weight / total_weight))
    rows = [r for r in rows if r[2] * total_bases > 300000]  # genomes with ~20 reads or more: the ratio is not noise
    rows.sort()
    third = len(rows) // 3
    print(f"genomes with >= 300 kb expected: {len(rows)}")
    for label, part in (("shortest contigs", rows[:third]), ("middle", rows[third:2 * third]), ("longest contigs", rows[2 * third:])):
        print(f"{label}: base-weighted contig length {statistics.median(r[0] for r in part) / 1000:.0f} kb, "
              f"bases share / weight share {statistics.mean(r[1] for r in part):.3f}")
    print(f"ratio over all: median {statistics.median(r[1] for r in rows):.3f}, SD {statistics.pstdev(r[1] for r in rows):.3f}")


if __name__ == "__main__":
    main()
