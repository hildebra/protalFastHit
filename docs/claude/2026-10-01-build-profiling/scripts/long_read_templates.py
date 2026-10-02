#!/usr/bin/env python3
"""Templates for one long-read sample, for pbsim3 --strategy templ (a prototype, for timing).

collect_training_data.py runs pbsim3 --strategy wgs once per genome of a sample. pbsim3 then
computes its read length table (100 bp to 1 Mb) for every contig, and makes reads for a contig
until their bases reach depth x its length, cutting the last one to the rest (at least 100 bp):
at the shallow design points every contig gets one short read. Here the sample's reads are drawn
first, as they arise: read lengths from a gamma distribution with the setup's mean and standard
deviation until the sample's bases are reached; each read from a genome chosen by relative
abundance x genome length, from a contig chosen by length, at a uniform start, on either strand,
cut where the contig ends. pbsim3 --strategy templ then adds its errors and qualities to each
template in one call.

usage: long_read_templates.py --manifest MANIFEST --sample SAMPLE --bases N --length-mean M
       --length-sd S --seed K -o templates.fa
"""

import argparse
import bisect
import csv
import gzip
import random


def read_fasta(path):
    """The contigs of a FASTA (plain or gzipped), as strings."""
    contigs, parts = [], []
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for line in fh:
            if line.startswith(">"):
                if parts:
                    contigs.append("".join(parts))
                parts = []
            else:
                parts.append(line.strip())
    if parts:
        contigs.append("".join(parts))
    return contigs


COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--manifest", required=True, help="the paired-end point's manifest.tsv (the communities)")
    p.add_argument("--sample", required=True, help="the community sample to draw from")
    p.add_argument("--bases", type=float, required=True)
    p.add_argument("--length-mean", type=float, required=True)
    p.add_argument("--length-sd", type=float, required=True)
    p.add_argument("--min-length", type=int, default=100)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("-o", "--out", required=True)
    opts = p.parse_args()
    rng = random.Random(opts.seed)
    with open(opts.manifest) as fh:
        genomes = [row for row in csv.DictReader(fh, delimiter="\t") if row["sample"] == opts.sample]
    weights = [float(g["relative_abundance"]) * float(g["genome_length"]) for g in genomes]
    # Reads until the sample's bases are reached: a genome by weight, a length from a gamma with the setup's
    # mean and sd, a start uniform over the genome's contigs; a read is cut where its contig ends.
    shape = (opts.length_mean / opts.length_sd) ** 2
    scale = opts.length_sd ** 2 / opts.length_mean
    contigs_of, ends_of = {}, {}
    written, total, used = 0, 0, set()
    with open(opts.out, "w") as out:
        while total < opts.bases:
            g = rng.choices(range(len(genomes)), weights=weights)[0]
            if g not in contigs_of:
                contigs_of[g] = [c for c in read_fasta(genomes[g]["fasta_path"]) if len(c) >= opts.min_length]
                ends = []
                for c in contigs_of[g]:
                    ends.append((ends[-1] if ends else 0) + len(c))
                ends_of[g] = ends
            contigs, ends = contigs_of[g], ends_of[g]
            length = max(opts.min_length, int(rng.gammavariate(shape, scale)))
            pos = rng.randrange(ends[-1])
            k = bisect.bisect_right(ends, pos)
            start = pos - (ends[k - 1] if k else 0)
            seq = contigs[k][start:start + length]
            if len(seq) < opts.min_length:  # at a contig's end: the read's end there
                seq = contigs[k][max(0, len(contigs[k]) - length):]
            if rng.random() < 0.5:
                seq = seq.translate(COMPLEMENT)[::-1]
            written += 1
            total += len(seq)
            used.add(g)
            out.write(f">g{g}_r{written}\n{seq}\n")
    print(f"{written} templates, {total} bases, from {len(used)} of {len(genomes)} genomes")


if __name__ == "__main__":
    main()
