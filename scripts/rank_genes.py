#!/usr/bin/env python3
"""Rank a database's marker genes by how distinctive they are, for a reduced database of the best N.

Usage:
  rank_genes.py --db FOLDER [-o gene_ranking.tsv] [--top N [--per-domain M] --subset genes.txt]

FOLDER holds the separate files of a built database (protal --build --no_bundle, or
protal --unpack_db on a database.protal): reference.map, unique_kmers.tsv and
internal_taxonomy.dmp, and if there gene2geneid.tsv (the markers' GTDB ids, written by the
converter; --gene-table names one elsewhere), gene_congeners.tsv and suspect_copies.tsv (reports
of the build).

A gene's score is prevalence x unique share:
  prevalence    the share of the database's species whose copy of the gene has k-mers in the
                index (unique_kmers.tsv lists one row per such copy); a gene most species lack
                detects few of them
  unique share  the share of a copy's k-mers unique to its species in the database (short or
                long unique, unique_kmers.tsv), averaged over the species that have the gene;
                a read of a gene with many unique k-mers names its species, one of a gene
                conserved between species does not
Both are given per domain as well (bacteria_prevalence, bacteria_score, archaea_prevalence,
archaea_score; NA for a domain the database lacks): GTDB's bacterial and archaeal marker sets
share only a few genes, and archaea are a few percent of the species, so a gene of the archaeal
set alone never scores high over all species. --top N therefore reserves --per-domain M of the N
genes for each domain (default a third of N, at least 1): the M best by that domain's score, then
the rest by the overall score; a gene good for both domains counts for both.

The table (-o, default stdout; one line per gene, best first by the overall score) also gives the
number of species, the copies' mean length, how the gene differs between congeneric species
against within species (between_factor and identical_share of gene_congeners.tsv, if there) and
the gene's suspect copies (suspect_copies.tsv, if there), so that a subset can be chosen on other
grounds too. --subset FILE writes the N genes chosen as a gene list for
scripts/mini_db/gtdb_to_protal_db.py --genes or protal --build --build_gene_subset (one id per
line, the marker and scores as comments); build_gtdb_database.py --n-genes N does all of this in a
build-and-train run (docs/building-a-database.md, reduced marker sets).
"""
import argparse
import collections
import os
import sys

DOMAINS = ("bacteria", "archaea")
COLUMNS = ("rank", "gene_id", "marker", "score", "prevalence", "unique_share", "species", "mean_length",
           "bacteria_prevalence", "bacteria_score", "archaea_prevalence", "archaea_score",
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
    """unique_kmers.tsv -> {gene id: {species taxid: unique share}} of the copies with k-mers in the index."""
    shares = collections.defaultdict(dict)
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or not f[0].isdigit():
                continue
            taxid, gene, short, long_, total = int(f[0]), int(f[1]), int(f[2]), int(f[4]), int(f[8])
            if total == 0:
                continue
            shares[gene][taxid] = (short + long_) / total
    return shares


def read_domains(path):
    """internal_taxonomy.dmp -> {species taxid: domain} ("bacteria", "archaea", or the domain's name lowercased
    without its d__)."""
    nodes = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 5 and f[0].isdigit():
                nodes[int(f[0])] = (int(f[1]), f[3], f[4])
    domains = {}
    for taxid, (_, _, rank) in nodes.items():
        if rank != "species":
            continue
        node, seen = taxid, set()
        while node in nodes and node not in seen:
            seen.add(node)
            parent, name, rank = nodes[node]
            if rank == "domain":
                domains[taxid] = name.removeprefix("d__").lower()
                break
            node = parent
    return domains


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


def read_names(path):
    """gene2geneid.tsv -> {gene id: marker}."""
    names = {}
    if path and os.path.isfile(path):
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2 and f[1].isdigit():
                    names[int(f[1])] = f[0]
    return names


def rank(folder, gene_table=None, taxonomy=None):
    """-> [row per gene, best first] of COLUMNS' values (as dicts), from the folder's files; gene2geneid.tsv and
    internal_taxonomy.dmp from `gene_table` and `taxonomy` if given, from the folder otherwise."""
    lengths, species = read_map(os.path.join(folder, "reference.map"))
    shares = read_unique_kmers(os.path.join(folder, "unique_kmers.tsv"))
    if not species or not shares:
        sys.exit(f"{folder}: reference.map or unique_kmers.tsv lists no genes")
    names = read_names(gene_table or os.path.join(folder, "gene2geneid.tsv"))
    domain_of = read_domains(taxonomy or os.path.join(folder, "internal_taxonomy.dmp"))
    by_domain = collections.Counter(domain_of.get(t, "unknown") for t in species)
    congeners = read_congeners(os.path.join(folder, "gene_congeners.tsv"))
    suspects = read_suspect_copies(os.path.join(folder, "suspect_copies.tsv"))
    rows = []
    for gene in sorted(lengths):
        copies = shares.get(gene, {})
        prevalence = len(copies) / len(species)
        unique_share = sum(copies.values()) / len(copies) if copies else 0.0
        row = {"gene_id": gene, "marker": names.get(gene, "NA"), "score": prevalence * unique_share,
               "prevalence": prevalence, "unique_share": unique_share, "species": len(copies),
               "mean_length": sum(lengths[gene]) / len(lengths[gene]), "suspect_copies": suspects.get(gene, 0)}
        row["between_factor"], row["identical_share"] = congeners.get(gene, ("NA", "NA"))
        for domain in DOMAINS:
            mine = [s for t, s in copies.items() if domain_of.get(t) == domain]
            if by_domain.get(domain):
                row[domain + "_prevalence"] = len(mine) / by_domain[domain]
                row[domain + "_score"] = row[domain + "_prevalence"] * (sum(mine) / len(mine) if mine else 0.0)
            else:
                row[domain + "_prevalence"] = row[domain + "_score"] = "NA"
        rows.append(row)
    rows.sort(key=lambda r: (-r["score"], -r["species"], r["gene_id"]))
    for i, row in enumerate(rows, 1):
        row["rank"] = i
    return rows


def per_domain_default(n):
    return max(1, n // 3)


def select(rows, n, per_domain=None):
    """The n genes of a subset: for each domain in the ranking (the one with more species first), the per_domain
    best by its score (default per_domain_default(n)), then the best by the overall score; -> the rows chosen,
    in the order of the ranking, each with "chosen_for" (the domains it was reserved for, or "overall")."""
    if per_domain is None:
        per_domain = per_domain_default(n)
    domains = [d for d in DOMAINS if any(r[d + "_score"] != "NA" for r in rows)]
    if len(domains) < 2:
        per_domain = 0  # one domain: the overall score is its score
    chosen = {}
    for domain in domains:
        best = sorted((r for r in rows if r[domain + "_score"] != "NA" and r[domain + "_score"] > 0),
                      key=lambda r: (-r[domain + "_score"], r["rank"]))
        for row in best[:per_domain]:
            chosen.setdefault(row["gene_id"], []).append(domain)
    for row in rows:
        if len(chosen) >= n:
            break
        chosen.setdefault(row["gene_id"], [])
    result = []
    for row in rows:
        if row["gene_id"] in chosen:
            row = dict(row, chosen_for=",".join(chosen[row["gene_id"]]) or "overall")
            result.append(row)
    return result[:n] if len(result) > n else result


def coverage(chosen, threshold=0.5):
    """How many of the genes chosen each domain has in at least `threshold` of its species: {domain: count}, for
    the domains in the ranking."""
    counts = {}
    for domain in DOMAINS:
        values = [r[domain + "_prevalence"] for r in chosen if r[domain + "_prevalence"] != "NA"]
        if values:
            counts[domain] = sum(v >= threshold for v in values)
    return counts


def read_ranking(path):
    """A ranking table this script wrote -> [row dict] in its order, the numbers as numbers (NA kept)."""
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
            row = {}
            for name, value in zip(header, f):
                if name in ("rank", "gene_id", "species", "suspect_copies"):
                    row[name] = int(value)
                elif name == "marker" or value == "NA":
                    row[name] = value
                else:
                    try:
                        row[name] = float(value)
                    except ValueError:
                        row[name] = value
            for domain in DOMAINS:  # a table of an earlier version: no domain columns
                row.setdefault(domain + "_prevalence", "NA")
                row.setdefault(domain + "_score", "NA")
            rows.append(row)
    if not rows:
        sys.exit(f"{path} ranks no genes")
    rows.sort(key=lambda r: r["rank"])
    return rows


def format_row(row):
    out = []
    for name in COLUMNS:
        value = row[name]
        out.append(f"{value:.4f}" if isinstance(value, float) else str(value))
    return "\t".join(out)


def write_table(path, rows):
    text = "\t".join(COLUMNS) + "\n" + "".join(format_row(r) + "\n" for r in rows)
    if path:
        with open(path + ".partial", "w", newline="\n") as fh:
            fh.write(text)
        os.replace(path + ".partial", path)
    else:
        sys.stdout.write(text)


def write_subset(path, chosen, total, source):
    """The gene list of the genes chosen (select): one gene id per line, each preceded by a comment with its
    marker, scores and what it was reserved for (protal --build_gene_subset reads only the id lines; the
    converter's --genes takes the first column). `source` says in the header where the ranking is from."""
    with open(path + ".partial", "w", newline="\n") as fh:
        fh.write(f"# the {len(chosen)} most distinctive of {total} marker genes, ranked from {source} "
                 "(scripts/rank_genes.py: prevalence x unique k-mer share, per domain too)\n")
        for row in chosen:
            domains = ", ".join(f"{d} {row[d + '_score']:.3f}" for d in DOMAINS if row[d + "_score"] != "NA")
            fh.write(f"# rank {row['rank']}: {row['marker']}, score {row['score']:.4f}"
                     f"{' (' + domains + ')' if domains else ''}, chosen for {row.get('chosen_for', 'overall')}\n"
                     f"{row['gene_id']}\n")
    os.replace(path + ".partial", path)


def describe(chosen):
    """One line on a subset: its genes and how many cover each domain."""
    covered = coverage(chosen)
    return (f"{', '.join(r['marker'] for r in chosen)}; scores {chosen[0]['score']:.3f} down to {chosen[-1]['score']:.3f}" +
            (f"; in half the species or more of: {', '.join(f'{d} {n}' for d, n in covered.items())}" if covered else ""))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--db", required=True, help="folder with the built database's separate files (reference.map, "
                                                "unique_kmers.tsv, internal_taxonomy.dmp, gene2geneid.tsv, gene_congeners.tsv, "
                                                "suspect_copies.tsv)")
    p.add_argument("--gene-table", help="gene2geneid.tsv of the database, if not in the folder")
    p.add_argument("--taxonomy", help="internal_taxonomy.dmp of the database, if not in the folder")
    p.add_argument("-o", "--output", help="the ranking table (default: stdout)")
    p.add_argument("--top", type=int, help="with --subset: how many genes to choose")
    p.add_argument("--per-domain", type=int, help="of them, how many are reserved for the best genes of each domain "
                                                   "(default: a third of --top, at least 1)")
    p.add_argument("--subset", help="write the genes chosen as a gene list (protal --build_gene_subset, "
                                    "gtdb_to_protal_db.py --genes)")
    args = p.parse_args()
    if bool(args.top) != bool(args.subset):
        p.error("--top and --subset go together")
    rows = rank(args.db, args.gene_table, args.taxonomy)
    write_table(args.output, rows)
    if args.subset:
        if args.top < 1 or args.top > len(rows):
            sys.exit(f"--top: {args.top} is not between 1 and the {len(rows)} genes of {args.db}")
        chosen = select(rows, args.top, args.per_domain)
        write_subset(args.subset, chosen, len(rows), args.db)
        sys.stderr.write(f"{args.subset}: {len(chosen)} genes: {describe(chosen)}\n")


if __name__ == "__main__":
    main()
