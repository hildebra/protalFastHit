#!/usr/bin/env python3
"""Rank a database's marker genes by how distinctive they are, for a reduced database of the best N.

Usage:
  rank_genes.py --db FOLDER [-o gene_ranking.tsv] [--top N --subset genes.txt]

FOLDER holds the separate files of a built database (protal --build --no_bundle, or
protal --unpack_db on a database.protal): reference.map and unique_kmers.tsv, and if there
gene2geneid.tsv (the markers' GTDB ids, written by the converter), gene_congeners.tsv and
suspect_copies.tsv (reports of the build).

A gene's score is prevalence x unique share:
  prevalence    the share of the database's species whose copy of the gene has k-mers in the
                index (unique_kmers.tsv lists one row per such copy); a gene most species lack
                detects few of them
  unique share  the share of a copy's k-mers unique to its species in the database (short or
                long unique, unique_kmers.tsv), averaged over the species that have the gene;
                a read of a gene with many unique k-mers names its species, one of a gene
                conserved between species does not
The table (-o, default stdout; one line per gene, best first) gives both, the number of species,
the copies' mean length, how the gene differs between congeneric species against within species
(between_factor and identical_share of gene_congeners.tsv, if there) and the gene's suspect copies
(suspect_copies.tsv, if there), so that a subset can be chosen on other grounds too. --top N
--subset FILE writes the N best gene ids as a gene list for scripts/mini_db/gtdb_to_protal_db.py
--genes or protal --build --build_gene_subset (one id per line, the marker and score as comments);
build_gtdb_database.py --n-genes N does all of this in a build-and-train run
(docs/building-a-database.md, reduced marker sets).
"""
import argparse
import collections
import os
import sys

COLUMNS = ("rank", "gene_id", "marker", "score", "prevalence", "unique_share", "species", "mean_length",
           "between_factor", "identical_share", "suspect_copies")


def read_map(path):
    """reference.map -> ({gene id: [lengths]}, species) of the genes' copies."""
    lengths = collections.defaultdict(list)
    species = set()
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 4:
                continue
            species.add(int(f[0]))
            lengths[int(f[1])].append(int(f[3]) - int(f[2]))
    return lengths, species


def read_unique_kmers(path):
    """unique_kmers.tsv -> {gene id: (copies with k-mers in the index, sum of their unique shares)}."""
    rows = collections.defaultdict(lambda: [0, 0.0])
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or not f[0].isdigit():
                continue
            gene, short, long_, total = int(f[1]), int(f[2]), int(f[4]), int(f[8])
            if total == 0:
                continue
            rows[gene][0] += 1
            rows[gene][1] += (short + long_) / total
    return rows


def read_congeners(path):
    """gene_congeners.tsv -> {gene id: (between_factor, identical_share)} (strings as written; NA if not)."""
    table = {}
    if not path or not os.path.isfile(path):
        return table
    with open(path) as fh:
        header = None
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if header is None:
                header = f
                continue
            row = dict(zip(header, f))
            if row.get("geneid", "").isdigit():
                table[int(row["geneid"])] = (row.get("between_factor", "NA"), row.get("identical_share", "NA"))
    return table


def read_suspect_copies(path):
    """suspect_copies.tsv -> {gene id: copies}."""
    counts = collections.Counter()
    if not path or not os.path.isfile(path):
        return counts
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) >= 2 and f[0].isdigit() and f[1].isdigit():
                counts[int(f[1])] += 1
    return counts


def rank(folder):
    """-> [row per gene, best first] of COLUMNS' values, from the folder's files."""
    lengths, species = read_map(os.path.join(folder, "reference.map"))
    unique = read_unique_kmers(os.path.join(folder, "unique_kmers.tsv"))
    if not species or not unique:
        sys.exit(f"{folder}: reference.map or unique_kmers.tsv lists no genes")
    names = {}
    if os.path.isfile(os.path.join(folder, "gene2geneid.tsv")):
        with open(os.path.join(folder, "gene2geneid.tsv")) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2 and f[1].isdigit():
                    names[int(f[1])] = f[0]
    congeners = read_congeners(os.path.join(folder, "gene_congeners.tsv"))
    suspects = read_suspect_copies(os.path.join(folder, "suspect_copies.tsv"))
    rows = []
    for gene in sorted(lengths):
        copies, shares = unique.get(gene, (0, 0.0))
        prevalence = copies / len(species)
        unique_share = shares / copies if copies else 0.0
        between, identical = congeners.get(gene, ("NA", "NA"))
        rows.append([gene, names.get(gene, "NA"), prevalence * unique_share, prevalence, unique_share, copies,
                     sum(lengths[gene]) / len(lengths[gene]), between, identical, suspects.get(gene, 0)])
    rows.sort(key=lambda r: (-r[2], -r[5], r[0]))
    return [[i + 1] + r for i, r in enumerate(rows)]


def read_ranking(path):
    """A ranking table this script wrote -> [[rank, gene id, marker, score], ...] in its order."""
    rows = []
    with open(path) as fh:
        header = None
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if header is None:
                header = f
                for name in ("rank", "gene_id", "marker", "score"):
                    if name not in header:
                        sys.exit(f"{path} is no ranking table of rank_genes.py (no {name} column)")
                continue
            if len(f) < len(header) or not f[header.index("gene_id")].isdigit():
                continue
            row = dict(zip(header, f))
            rows.append([int(row["rank"]), int(row["gene_id"]), row["marker"], float(row["score"])])
    if not rows:
        sys.exit(f"{path} ranks no genes")
    rows.sort(key=lambda r: r[0])
    return rows


def format_row(row):
    out = []
    for value in row:
        out.append(f"{value:.4f}" if isinstance(value, float) and not isinstance(value, bool) else str(value))
    return "\t".join(out)


def write_subset(path, rows, n, source):
    """The gene list of the n best genes of rows ([rank, gene id, marker, score, ...]): one gene id per line, each
    preceded by a comment with its marker and score (protal --build_gene_subset reads only the id lines; the
    converter's --genes takes the first column). `source` says in the header where the ranking is from."""
    chosen = rows[:n]
    with open(path + ".partial", "w", newline="\n") as fh:
        fh.write(f"# the {len(chosen)} most distinctive of {len(rows)} marker genes, ranked from {source} "
                 "(scripts/rank_genes.py: prevalence x unique k-mer share)\n")
        for row in chosen:
            fh.write(f"# rank {row[0]}: {row[2]}, score {row[3]:.4f}\n{row[1]}\n")
    os.replace(path + ".partial", path)
    return chosen


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--db", required=True, help="folder with the built database's separate files (reference.map, "
                                                "unique_kmers.tsv, gene2geneid.tsv, gene_congeners.tsv, suspect_copies.tsv)")
    p.add_argument("-o", "--output", help="the ranking table (default: stdout)")
    p.add_argument("--top", type=int, help="with --subset: how many of the best genes to list")
    p.add_argument("--subset", help="write the --top genes as a gene list (protal --build_gene_subset, "
                                    "gtdb_to_protal_db.py --genes)")
    args = p.parse_args()
    if bool(args.top) != bool(args.subset):
        p.error("--top and --subset go together")
    rows = rank(args.db)
    text = "\t".join(COLUMNS) + "\n" + "".join(format_row(r) + "\n" for r in rows)
    if args.output:
        with open(args.output + ".partial", "w", newline="\n") as fh:
            fh.write(text)
        os.replace(args.output + ".partial", args.output)
    else:
        sys.stdout.write(text)
    if args.subset:
        if args.top < 1 or args.top > len(rows):
            sys.exit(f"--top: {args.top} is not between 1 and the {len(rows)} genes of {args.db}")
        chosen = write_subset(args.subset, rows, args.top, args.db)
        sys.stderr.write(f"{args.subset}: the {len(chosen)} best genes, scores {chosen[0][3]:.3f} down to {chosen[-1][3]:.3f}"
                         f" ({', '.join(str(r[2]) for r in chosen)})\n")


if __name__ == "__main__":
    main()
