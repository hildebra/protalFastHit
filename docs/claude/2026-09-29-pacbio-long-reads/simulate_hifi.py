#!/usr/bin/env python3
"""simulate_hifi.py - PacBio-HiFi-like reads from a mock community of mini-GTDB genomes.

Read lengths are normal (--mean, --sd) within [--min, --max]; a read starts uniformly in a contig
(shorter contigs give shorter reads), on either strand, and gets --error errors per base (40%
substitutions, 30% 1 bp insertions, 30% 1 bp deletions), quality 'I'. A genome gets bases in
proportion to relative abundance x genome length (cell abundance, as protal estimates it).

Writes <out>.fq and <out>.truth.tsv (read, accession, species, contig, start, end, strand).
"""

import argparse
import gzip
import random

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def read_table(path):
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip() and not l.startswith("#")]
    header = lines[0].split("\t")
    return [dict(zip(header, l.split("\t"))) for l in lines[1:]]


def read_contigs(path):
    opener = gzip.open if path.endswith(".gz") else open
    contigs, name, chunks = [], None, []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    contigs.append((name, "".join(chunks).upper()))
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
    if name:
        contigs.append((name, "".join(chunks).upper()))
    return contigs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--community", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--bases", type=float, default=20e6)
    ap.add_argument("--mean", type=int, default=15000)
    ap.add_argument("--sd", type=int, default=3000)
    ap.add_argument("--min", type=int, default=5000)
    ap.add_argument("--max", type=int, default=30000)
    ap.add_argument("--error", type=float, default=0.001)
    ap.add_argument("--seed", type=int, default=5)
    args = ap.parse_args()
    rng = random.Random(args.seed)

    genomes = {g["accession"]: g for g in read_table(args.genomes)}
    community = {row["accession"]: float(row["relative_abundance"]) for row in read_table(args.community)}
    contigs = {acc: read_contigs(genomes[acc]["fasta_path"]) for acc in community}
    lengths = {acc: sum(len(s) for _, s in contigs[acc]) for acc in community}
    weight = {acc: community[acc] * lengths[acc] for acc in community}
    total = sum(weight.values())

    def mutate(seq):
        out = []
        for b in seq:
            r = rng.random()
            if r >= args.error:
                out.append(b)
            elif r < 0.4 * args.error:
                out.append(rng.choice([c for c in "ACGT" if c != b]))
            elif r < 0.7 * args.error:
                out.append(b + rng.choice("ACGT"))
        return "".join(out)

    n = 0
    with open(args.out + ".fq", "w") as fq, open(args.out + ".truth.tsv", "w") as truth:
        truth.write("read\taccession\tspecies\tcontig\tstart\tend\tstrand\n")
        for acc in community:
            species = genomes[acc]["gtdb_taxonomy"].split(";")[-1]
            budget = args.bases * weight[acc] / total
            clen = [len(s) for _, s in contigs[acc]]
            done = 0
            while done < budget:
                contig, seq = rng.choices(contigs[acc], weights=clen)[0]
                length = min(len(seq), max(args.min, min(args.max, int(rng.gauss(args.mean, args.sd)))))
                start = rng.randint(0, len(seq) - length)
                frag = seq[start:start + length]
                strand = "+" if rng.random() < 0.5 else "-"
                if strand == "-":
                    frag = frag[::-1].translate(COMPLEMENT)
                read = mutate(frag)
                n += 1
                name = f"m64001_000000/{n}/ccs"
                fq.write(f"@{name}\n{read}\n+\n{'I' * len(read)}\n")
                truth.write(f"{name}\t{acc}\t{species}\t{contig}\t{start}\t{start + length}\t{strand}\n")
                done += length
    print(f"{n} reads, {args.bases:.0f} bases")


if __name__ == "__main__":
    main()
