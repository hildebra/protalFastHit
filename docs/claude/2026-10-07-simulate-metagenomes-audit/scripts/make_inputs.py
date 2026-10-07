#!/usr/bin/env python3
"""Inputs for the simulate_metagenomes benchmarks of this report, from a genome table of real-sized genomes
(name, lineage, FASTA, length; bench_pairs.sh of the Illumina report writes one):

  long WORK TABLE BASES SAMPLES   WORK/long_samples.tsv and WORK/long_genomes.tsv for --long_samples: SAMPLES samples
                                  of every genome of the table, lognormal (sigma 1.5) abundances times length as weights,
                                  BASES template bases each, reads into @OUT@/<sample>.fq.zst (the caller puts a folder
                                  in place of @OUT@)
  big WORK TABLE ROWS             WORK/big_table.tsv: ROWS genomes in ROWS/2 species of ROWS/20 genera (lineages of
                                  their own), each pointing at one of the table's FASTAs, with its length: a design's
                                  cost at GTDB scale (--test reads no genome)
"""

import os
import random
import sys


def read_table(path):
    with open(path) as fh:
        return [line.rstrip("\n").split("\t") for line in fh if line.strip()]


def long_inputs(work, table, bases, samples):
    rows = read_table(table)
    rng = random.Random(7)
    with open(os.path.join(work, "long_samples.tsv"), "w") as s, open(os.path.join(work, "long_genomes.tsv"), "w") as g:
        s.write("sample\tout\tbases\tseed\n")
        g.write("sample\tgenome\tfasta\tweight\thost\n")
        for i in range(samples):
            name = f"long{i + 1}"
            s.write(f"{name}\t@OUT@/{name}.fq.zst\t{bases}\t{101 + i}\n")
            for row in rows:
                g.write(f"{name}\t{row[0]}\t{row[2]}\t{rng.lognormvariate(0, 1.5) * int(row[3])!r}\t0\n")


def big_table(work, table, n):
    rows = read_table(table)
    with open(os.path.join(work, "big_table.tsv"), "w") as fh:
        for i in range(n):
            species, genus = i // 2, i // 20
            source = rows[i % len(rows)]
            fh.write(f"B{i}\td__Bacteria;p__P{genus % 40};c__C{genus % 200};o__O{genus % 600};f__F{genus % 1500};"
                     f"g__G{genus};s__G{genus} sp{species}\t{source[2]}\t{source[3]}\n")


def main():
    kind, work, table = sys.argv[1:4]
    os.makedirs(work, exist_ok=True)
    if kind == "long":
        long_inputs(work, table, int(sys.argv[4]), int(sys.argv[5]))
    elif kind == "big":
        big_table(work, table, int(sys.argv[4]))
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main()
