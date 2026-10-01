#!/usr/bin/env python3
"""gene_neighbours.py - which marker genes lie next to each other in the genomes, per clade, for protal.

Places the marker genes of a converted database (gtdb_to_protal_db.py's reference.fna, reference.map and
internal_taxonomy.dmp) in their species' representative genomes, for the species whose genome is at hand
(a genome table as simulate_metagenomes takes: accession, GTDB taxonomy, FASTA path), and counts per clade
(by default family, order, class, phylum and domain) how often each end of each gene has which neighbour
within --max_gap bases, and how far away. protal --build packs the result, gene_neighbours.tsv, into the
database; protal uses it to follow a fragment or a long read from one gene into the next.

A marker gene of the database is its representative genome's own gene call, so it is found by its exact
sequence, on either strand. One genome per species (the representative's) counts, so that species with
many genomes do not outweigh the others. Each gene has two ends, 5' and 3' in its coding orientation. An end
is informative when the contig goes on for --max_gap bases past it: its partner is then the next placed
marker within --max_gap, or none. An end within --max_gap of its contig's end says nothing (a fragmented
assembly is no evidence that a gene has no neighbour).

Writes <db>/gene_neighbours.tsv (or --output), one line per clade, gene end and partner:

  clade  gene  end  partner  partner_end  species  informative  gap_median  gap_min  gap_max

clade: an internal taxid; gene, partner: gene ids (partner 0: no marker within --max_gap, partner_end 0);
end, partner_end: 5 or 3, the ends that face each other; species: the clade's species with this partner
at this end; informative: its species in which this end was informative; gap: bases between the two genes
(negative: they overlap). Lines starting with '#' are comments.

Usage:
  gene_neighbours.py --db <converted folder> --genome_table <table> [-t 8] [--max_gap 3000]
      [--ranks family,order,class,phylum,domain] [--positions FILE]
"""

import argparse
import bisect
import collections
import gzip
import multiprocessing
import os
import statistics
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
from gtdb_to_protal_db import ACCESSION_RE, normalize_accession  # noqa: E402

FILE_NAME = "gene_neighbours.tsv"
HEADER = "clade\tgene\tend\tpartner\tpartner_end\tspecies\tinformative\tgap_median\tgap_min\tgap_max"
COMPLEMENT = bytes.maketrans(b"ACGTN", b"TGCAN")
SEPARATOR = b"|"  # between contigs: no gene matches across it


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_taxonomy(path):
    """internal_taxonomy.dmp -> ({taxid: (parent, rank, name)}, {representative accession: species taxid})."""
    nodes, reps = {}, {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            nodes[int(f[0])] = (int(f[1]), f[4], f[3])
            if f[4] == "species" and len(f) > 6 and f[6]:
                reps[normalize_accession(f[6])] = int(f[0])
    return nodes, reps


def ancestors(nodes, taxid, ranks):
    """The clades of `ranks` above taxid, {rank: clade taxid}."""
    found, seen = {}, set()
    while taxid in nodes and taxid not in seen:
        seen.add(taxid)
        parent, rank, _ = nodes[taxid]
        if rank in ranks:
            found[rank] = taxid
        taxid = parent
    return found


def read_genome_table(path):
    """accession -> FASTA path, from a genome table (accession, taxonomy, path; a header line is skipped)."""
    paths = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or line.startswith("#") or f[0] == "accession":
                continue
            m = ACCESSION_RE.search(f[0])
            paths.setdefault(m.group(1) if m else f[0], f[2])
    return paths


def read_map(path):
    """reference.map -> {taxid: [(gene id, start byte, end byte)]}."""
    genes = collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) >= 4:
                genes[int(f[0])].append((int(f[1]), int(f[2]), int(f[3])))
    return genes


def read_contigs(path):
    """A genome FASTA (plain or gzip) -> [(name, sequence as upper-case bytes)]."""
    contigs, name, chunks = [], None, []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name is not None:
                    contigs.append((name, b"".join(chunks).upper()))
                name, chunks = line[1:].split()[0].decode() if line[1:].split() else "", []
            else:
                chunks.append(line.strip())
    if name is not None:
        contigs.append((name, b"".join(chunks).upper()))
    return contigs


def place(item):
    """One species: its genes found in its genome. -> (taxid, accession, contigs [(name, length)],
    placements [(gene id, contig, start, end, strand)], genes looked for, genes found more than once)."""
    taxid, accession, fasta, reference, genes = item
    try:
        contigs = read_contigs(fasta)
    except (OSError, EOFError) as error:
        return taxid, accession, None, str(error), len(genes), 0
    genome = SEPARATOR.join(seq for _, seq in contigs)
    offsets, at = [], 0
    for _, seq in contigs:
        offsets.append(at)
        at += len(seq) + 1
    placements, repeated = [], 0
    with open(reference, "rb") as fh:
        for gene, start, _end in genes:
            fh.seek(start)
            seq = fh.readline().strip().upper()
            if not seq:
                continue
            strand, pos = "+", genome.find(seq)
            if pos < 0:
                strand, pos = "-", genome.find(seq.translate(COMPLEMENT)[::-1])
            if pos < 0:
                continue
            other = seq if strand == "+" else seq.translate(COMPLEMENT)[::-1]
            repeated += genome.find(other, pos + 1) >= 0
            contig = bisect.bisect_right(offsets, pos) - 1
            local = pos - offsets[contig]
            placements.append((gene, contig, local, local + len(seq), strand))
    return taxid, accession, [(n, len(s)) for n, s in contigs], placements, len(genes), repeated


def neighbour_ends(contigs, placements, max_gap):
    """The informative ends of the placed genes: [(gene, end, partner, partner end, gap)], partner 0 for
    no marker within max_gap."""
    by_contig = collections.defaultdict(list)
    for gene, contig, start, end, strand in placements:
        by_contig[contig].append((start, end, gene, strand))
    ends = []
    for contig, genes in by_contig.items():
        genes.sort()
        length = contigs[contig][1]
        for i, (start, end, gene, strand) in enumerate(genes):
            right_end = 3 if strand == "+" else 5   # the gene's end at its right in the genome
            left_end = 5 if strand == "+" else 3
            if i + 1 < len(genes) and genes[i + 1][0] - end <= max_gap:
                nstart, _, ngene, nstrand = genes[i + 1]
                ends.append((gene, right_end, ngene, 5 if nstrand == "+" else 3, nstart - end))
            elif i + 1 < len(genes) or length - end >= max_gap:
                ends.append((gene, right_end, 0, 0, 0))
            if i > 0 and start - genes[i - 1][1] <= max_gap:
                _, pend, pgene, pstrand = genes[i - 1]
                ends.append((gene, left_end, pgene, 3 if pstrand == "+" else 5, start - pend))
            elif i > 0 or start >= max_gap:
                ends.append((gene, left_end, 0, 0, 0))
    return ends


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--db", required=True, help="converted database folder (reference.fna, reference.map, "
                                                "internal_taxonomy.dmp; before protal --build packs them)")
    ap.add_argument("--genome_table", required=True, help="accession, GTDB taxonomy, FASTA path (a header line is skipped)")
    ap.add_argument("--output", help=f"default: <db>/{FILE_NAME}")
    ap.add_argument("--max_gap", type=int, default=3000, help="the farthest a neighbour counts, in bases (default 3000)")
    ap.add_argument("--ranks", default="family,order,class,phylum,domain",
                    help="the ranks of the clades to count for (default family,order,class,phylum,domain)")
    ap.add_argument("--min_placed", type=float, default=0.8,
                    help="a genome counts only if this share of its species' genes is found in it (default 0.8)")
    ap.add_argument("--positions", help="also write where each gene was found: accession, taxid, gene id, contig, "
                                        "start (1-based), end, strand")
    ap.add_argument("-t", "--threads", type=int, default=1)
    args = ap.parse_args()

    reference = os.path.join(args.db, "reference.fna")
    for name in ("reference.fna", "reference.map", "internal_taxonomy.dmp"):
        if not os.path.isfile(os.path.join(args.db, name)):
            sys.exit(f"{args.db} has no {name}: run this on the converted folder before protal --build packs it")
    ranks = [r.strip() for r in args.ranks.split(",") if r.strip()]
    nodes, reps = read_taxonomy(os.path.join(args.db, "internal_taxonomy.dmp"))
    genes = read_map(os.path.join(args.db, "reference.map"))
    table = read_genome_table(args.genome_table)
    work = [(taxid, acc, table[acc], reference, genes[taxid]) for acc, taxid in sorted(reps.items(), key=lambda kv: kv[1])
            if acc in table and genes.get(taxid)]
    print(f"{len(work)} of {len(reps)} species have their representative genome in {args.genome_table}", flush=True)

    counts = collections.defaultdict(lambda: collections.defaultdict(list))  # (clade, gene, end) -> partner -> gaps
    informative = collections.Counter()  # (clade, gene, end) -> species
    used = unreadable = sparse = looked = found = repeated = 0
    species_ends = species_adjacent = 0
    gaps_all = []
    positions = open(args.positions, "w", newline="\n") if args.positions else None
    if positions:
        positions.write("accession\ttaxid\tgene\tcontig\tstart\tend\tstrand\n")
    with multiprocessing.Pool(max(1, args.threads)) as pool:
        for taxid, acc, contigs, placements, n_genes, n_repeated in pool.imap(place, work, chunksize=4):
            if contigs is None:
                unreadable += 1
                sys.stderr.write(f"Cannot read the genome of {acc}: {placements}\n")
                continue
            looked += n_genes
            found += len(placements)
            if len(placements) < args.min_placed * n_genes:
                sparse += 1
                continue
            used += 1
            repeated += n_repeated
            if positions:
                for gene, contig, start, end, strand in placements:
                    positions.write(f"{acc}\t{taxid}\t{gene}\t{contigs[contig][0]}\t{start + 1}\t{end}\t{strand}\n")
            clades = ancestors(nodes, taxid, ranks).values()
            for gene, end, partner, partner_end, gap in neighbour_ends(contigs, placements, args.max_gap):
                species_ends += 1
                if partner:
                    species_adjacent += 1
                    gaps_all.append(gap)
                for clade in clades:
                    counts[(clade, gene, end)][(partner, partner_end)].append(gap)
                    informative[(clade, gene, end)] += 1
    if positions:
        positions.close()

    output = args.output or os.path.join(args.db, FILE_NAME)
    rows = 0
    with open(output + ".partial", "w", newline="\n") as fh:
        fh.write(f"# protal gene neighbours: genomes={used} max_gap={args.max_gap} ranks={','.join(ranks)}\n")
        fh.write(HEADER + "\n")
        for key in sorted(counts):
            clade, gene, end = key
            for (partner, partner_end), gaps in sorted(counts[key].items()):
                if partner:
                    median, low, high = int(round(statistics.median(gaps))), min(gaps), max(gaps)
                else:
                    median = low = high = 0
                fh.write(f"{clade}\t{gene}\t{end}\t{partner}\t{partner_end}\t{len(gaps)}\t{informative[key]}\t"
                         f"{median}\t{low}\t{high}\n")
                rows += 1
    os.replace(output + ".partial", output)

    # The snapshot: how much adjacency there is, and how alike it is within clades of each rank.
    print(f"Genes placed: {found} of {looked} ({100 * found / max(1, looked):.1f}%); {used} genomes used, "
          f"{sparse} with fewer than {args.min_placed:.0%} of their genes found, {unreadable} unreadable; "
          f"{repeated} genes found more than once (the first place taken)")
    if species_ends:
        gaps = ""
        if len(gaps_all) > 1:
            q = statistics.quantiles(gaps_all, n=10)
            gaps = f"; gaps: median {statistics.median(gaps_all):.0f}, 10-90% {q[0]:.0f} to {q[-1]:.0f}"
        print(f"Gene ends: {species_ends} informative, {species_adjacent} ({100 * species_adjacent / species_ends:.1f}%) "
              f"with a marker within {args.max_gap} bases{gaps}")
    rank_of = {clade: nodes[clade][1] for clade, _, _ in counts}
    for rank in ranks:
        shares, weights, clades = [], [], set()
        for (clade, gene, end), partners in counts.items():
            if rank_of[clade] != rank:
                continue
            clades.add(clade)
            n = informative[(clade, gene, end)]
            shares.append(max(len(g) for g in partners.values()) / n)
            weights.append(n)
        if weights:
            mean = sum(s * w for s, w in zip(shares, weights)) / sum(weights)
            print(f"  {rank}: {len(clades)} clades; the most common partner of a gene end is that of {100 * mean:.1f}% "
                  f"of the clade's species (weighted by species)")
    print(f"Wrote {output}: {rows} lines")


if __name__ == "__main__":
    main()
