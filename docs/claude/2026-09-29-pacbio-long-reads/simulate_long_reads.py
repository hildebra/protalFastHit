#!/usr/bin/env python3
"""simulate_long_reads.py - PacBio- or ONT-like long reads from a mock community of mini-GTDB genomes.

--platform sets the defaults (each can be overridden):
  hifi: lengths normal (mean 15 kb, sd 3 kb) within 5-30 kb; 0.1% errors (40% substitutions, 30%
        1 bp insertions, 30% 1 bp deletions); quality 'I'.
  ont:  lengths log-normal (median 20 kb, sigma 0.6) within 2-150 kb; 1.6% errors (30%
        substitutions, 25% insertions, 45% deletions) plus, in homopolymers of 4 or more, one base
        lost with probability 0.15 (about 2% errors in all, Q17); quality '2' (Q17).
A read starts uniformly in a contig (shorter contigs give shorter reads), on either strand. A genome
gets bases in proportion to relative abundance x genome length (cell abundance, as protal
estimates it).

Writes <out>.fq (or <out>.fa with --fasta) and <out>.truth.tsv (read, accession, species, contig,
start, end, strand).
"""

import argparse
import gzip
import math
import random

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
PLATFORMS = {
    "hifi": dict(lengths="normal", mean=15000, sd=3000, min=5000, max=30000, error=0.001,
                 shares=(0.4, 0.3, 0.3), homopolymer=0.0, quality="I"),
    "ont": dict(lengths="lognormal", mean=20000, sd=0.6, min=2000, max=150000, error=0.016,
                shares=(0.30, 0.25, 0.45), homopolymer=0.15, quality="2"),
}


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
    ap.add_argument("--platform", choices=sorted(PLATFORMS), default="hifi")
    ap.add_argument("--bases", type=float, default=20e6)
    ap.add_argument("--error", type=float, help="per-base error rate (outside homopolymer losses)")
    ap.add_argument("--homopolymer", type=float, help="chance to lose a base of a homopolymer of 4+")
    ap.add_argument("--fasta", action="store_true", help="write FASTA (no qualities)")
    ap.add_argument("--seed", type=int, default=5)
    args = ap.parse_args()
    p = dict(PLATFORMS[args.platform])
    if args.error is not None:
        p["error"] = args.error
    if args.homopolymer is not None:
        p["homopolymer"] = args.homopolymer
    rng = random.Random(args.seed)

    genomes = {g["accession"]: g for g in read_table(args.genomes)}
    community = {row["accession"]: float(row["relative_abundance"]) for row in read_table(args.community)}
    contigs = {acc: read_contigs(genomes[acc]["fasta_path"]) for acc in community}
    lengths = {acc: sum(len(s) for _, s in contigs[acc]) for acc in community}
    weight = {acc: community[acc] * lengths[acc] for acc in community}
    total = sum(weight.values())
    sub, ins, _ = p["shares"]

    def draw_length():
        if p["lengths"] == "normal":
            length = int(rng.gauss(p["mean"], p["sd"]))
        else:
            length = int(math.exp(rng.gauss(math.log(p["mean"]), p["sd"])))
        return max(p["min"], min(p["max"], length))

    def mutate(seq):
        out, i = [], 0
        while i < len(seq):
            b = seq[i]
            run = 1
            while i + run < len(seq) and seq[i + run] == b:
                run += 1
            # a homopolymer of 4 or more may lose a base
            emitted = run - (1 if run >= 4 and rng.random() < p["homopolymer"] else 0)
            for _ in range(emitted):
                r = rng.random()
                if r >= p["error"]:
                    out.append(b)
                elif r < sub * p["error"]:
                    out.append(rng.choice([c for c in "ACGT" if c != b]))
                elif r < (sub + ins) * p["error"]:
                    out.append(b + rng.choice("ACGT"))
            i += run
        return "".join(out)

    n = 0
    suffix = ".fa" if args.fasta else ".fq"
    with open(args.out + suffix, "w") as fq, open(args.out + ".truth.tsv", "w") as truth:
        truth.write("read\taccession\tspecies\tcontig\tstart\tend\tstrand\n")
        for acc in community:
            species = genomes[acc]["gtdb_taxonomy"].split(";")[-1]
            budget = args.bases * weight[acc] / total
            clen = [len(s) for _, s in contigs[acc]]
            done = 0
            while done < budget:
                contig, seq = rng.choices(contigs[acc], weights=clen)[0]
                length = min(len(seq), draw_length())
                start = rng.randint(0, len(seq) - length)
                frag = seq[start:start + length]
                strand = "+" if rng.random() < 0.5 else "-"
                if strand == "-":
                    frag = frag[::-1].translate(COMPLEMENT)
                read = mutate(frag)
                n += 1
                name = f"read{n}"
                if args.fasta:
                    fq.write(f">{name}\n{read}\n")
                else:
                    fq.write(f"@{name}\n{read}\n+\n{p['quality'] * len(read)}\n")
                truth.write(f"{name}\t{acc}\t{species}\t{contig}\t{start}\t{start + length}\t{strand}\n")
                done += length
    print(f"{n} reads, {args.bases:.0f} bases, platform {args.platform}")


if __name__ == "__main__":
    main()
