#!/usr/bin/env python3
"""Reads of known divergence from a protal database's own marker genes, for the gates' loss (F1 audit, B2).

Picks --species species of the reference (reference.fna of an unpacked database: genes named <taxid>_<geneid>) with
at least --min-genes genes, mutates every gene of each at --divergence (substitutions at uniformly drawn positions, each
to another base; no indels), and draws --fragments fragments per species from the mutated genes, by gene length, with
substitution errors at --error. Single-end reads (--mates 1) are a fragment's first --length bases; paired-end reads
(--mates 2) are both ends of a fragment of --insert (+- 10%) bases, the second reverse-complemented. A fragment that
does not fit its gene is cut to the gene. Read names: <taxid>_<geneid>_<n> (and /1, /2). Standard library only.

usage: b2_reads.py REFERENCE_FNA OUT_PREFIX --divergence 0.05 --length 150 --mates 2 [--species 40] [--fragments 500]
       [--seed 1] [--insert 300] [--error 0.001] [--min-genes 80]
"""
import argparse
import random
import sys

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def read_fasta(path):
    genes, name, parts = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    genes[name] = "".join(parts).upper()
                name, parts = line[1:].split()[0], []
            elif line:
                parts.append(line)
    if name:
        genes[name] = "".join(parts).upper()
    return genes


def mutate(seq, rate, rng):
    """`seq` with round(rate * len) positions changed, each to another base."""
    seq = list(seq)
    for pos in rng.sample(range(len(seq)), round(rate * len(seq))):
        seq[pos] = rng.choice([b for b in "ACGT" if b != seq[pos]])
    return "".join(seq)


def with_errors(seq, rate, rng):
    return "".join(rng.choice([b for b in "ACGT" if b != c]) if rng.random() < rate else c for c in seq)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reference")
    ap.add_argument("out")
    ap.add_argument("--divergence", type=float, required=True)
    ap.add_argument("--length", type=int, required=True)
    ap.add_argument("--mates", type=int, choices=(1, 2), required=True)
    ap.add_argument("--species", type=int, default=40)
    ap.add_argument("--fragments", type=int, default=500)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--insert", type=int, default=300)
    ap.add_argument("--error", type=float, default=0.001)
    ap.add_argument("--min-genes", type=int, default=80)
    args = ap.parse_args()

    genes = read_fasta(args.reference)
    by_taxon = {}
    for name, seq in genes.items():
        taxid = name.split("_")[0]
        if len(seq) >= args.length + 20:
            by_taxon.setdefault(taxid, []).append(name)
    eligible = sorted(t for t, names in by_taxon.items() if len(names) >= args.min_genes)
    # The species: the same for every divergence and setup (a seed of their own); the mutations and reads by --seed.
    species = random.Random(1000 + args.seed).sample(eligible, min(args.species, len(eligible)))
    # The mutated genes by --seed and --divergence alone, the same for every read setup; the reads by the setup too.
    mutation_rng = random.Random(f"{args.seed}:{args.divergence}")
    rng = random.Random(f"{args.seed}:{args.divergence}:{args.length}:{args.mates}")
    outs = [open(f"{args.out}_{m}.fq", "w") for m in range(1, args.mates + 1)]
    n = 0
    for taxid in species:
        names = sorted(by_taxon[taxid])
        mutated = {name: mutate(genes[name], args.divergence, mutation_rng) for name in names}
        weights = [len(mutated[name]) for name in names]
        for _ in range(args.fragments):
            name = rng.choices(names, weights)[0]
            seq = mutated[name]
            size = args.length if args.mates == 1 else int(args.insert * rng.uniform(0.9, 1.1))
            size = min(size, len(seq))
            start = rng.randrange(0, len(seq) - size + 1)
            fragment = seq[start:start + size]
            if rng.random() < 0.5:
                fragment = fragment.translate(COMPLEMENT)[::-1]
            reads = [fragment[:args.length]]
            if args.mates == 2:
                reads.append(fragment.translate(COMPLEMENT)[::-1][:args.length])
            n += 1
            for m, (out, read) in enumerate(zip(outs, reads), 1):
                suffix = f"/{m}" if args.mates == 2 else ""
                read = with_errors(read, args.error, rng)
                out.write(f"@{name}_{n}{suffix}\n{read}\n+\n{'I' * len(read)}\n")
    for out in outs:
        out.close()
    print(f"{len(species)} species, {n} fragments", file=sys.stderr)


if __name__ == "__main__":
    sys.exit(main())
