#!/usr/bin/env python3
"""gene_neighbours.py - which marker genes lie next to each other in the genomes, per clade, for protal.

Places the marker genes of a converted database (gtdb_to_protal_db.py's reference.fna, reference.map,
internal_taxonomy.dmp and genome2tiid.tsv) in every genome of their species at hand (a genome table as
simulate_metagenomes takes: accession, GTDB taxonomy, FASTA path; the representatives and the other strains
downloaded to simulate training data), and counts per clade (by default family, order, class, phylum and
domain) how often each end of each gene has which neighbour within --max_gap bases, and how far away. protal
--build packs both results into the database: gene_positions.tsv, where each gene lies in each genome, which
protal keeps but does not read, and gene_neighbours.tsv, the clades' frequencies, which it loads to follow a
fragment or a long read from one gene into the next and to judge whether the genes next to each other on a
taxon's reads are neighbours in its clade.

Placing. A database gene is its species' representative's own gene call: in the representative it is found by
its exact sequence, on either strand. In another strain it differs by the strain's mutations, so it is found by
its k-mer trace: each of the gene's KMER-mers is looked up among the genome's KMER-mers at every STRIDE-th
position, a hit giving a diagonal (genome position less the gene's offset); the band of near diagonals with the
most hits places the gene if it holds MIN_KMER_SHARE or more of the gene's k-mers that can hit there (genes up
to ~9% different). A gene is placed exactly where it can be (its first occurrence, the forward strand first),
else by its trace. A genome counts only if --min_placed of its species' genes are placed in it.

Counting. Each gene has two ends, 5' and 3' in its coding orientation. An end is informative in a genome when
the contig goes on for --max_gap bases past it: its partner is then the next placed marker within --max_gap,
or none. An end within --max_gap of its contig's end says nothing (a fragmented assembly is no evidence that a
gene has no neighbour). A sequence whose header marks it a whole replicon (NCBI's "complete genome" and
"complete sequence", CIRCULAR_MARKS) is circular, and so is a species representative's genome in one sequence
(a closed chromosome, whatever its header says): its last gene faces its first across the origin, and every
end on it is informative. A species counts once: at each gene end its genomes' most common partner (the
representative's breaks a tie), and the end informative if it is in any of its genomes, so that another
strain's complete assembly fills in where the representative's contigs end. A clade's line for a gene end and
partner says in how many of its species that end faces that partner, of how many in which the end is
informative: the observed frequency. A species also gets lines of its own (clade: its taxid, 1 of 1 species) for
the partners any of its genomes shows that would be no expected neighbours by its clades' frequencies alone
(count_clades; --no_species_lines leaves them out), so that protal judges the genes next to each other on the
reads of any of its genomes its neighbours.

Writes into <db> (or --output and --positions):

gene_neighbours.tsv, one line per clade, gene end and partner:
  clade  gene  end  partner  partner_end  species  informative  gap_median  gap_min  gap_max
clade: an internal taxid; gene, partner: gene ids (partner 0: no marker within --max_gap, partner_end 0); end,
partner_end: 5 or 3, the ends that face each other; species: the clade's species with this partner at this
end; informative: its species in which this end was informative; gap: bases between the two genes (negative:
they overlap), per species the median over its genomes.

gene_positions.tsv, one line per gene placed in a genome:
  accession  taxid  contig  contig_length  circular  gene  start  end  strand  placed  kmer_share
circular: 1 for a whole replicon (by its header, or a representative in one sequence); start, end: 1-based,
inclusive; placed: exact or trace; kmer_share: the share of the gene's k-mers that can hit its diagonal band
that do (1 when exact). Lines starting with '#' are comments in both.

--from_positions derives gene_neighbours.tsv from a gene_positions.tsv instead of from the genomes, without
the species of --exclude_taxa (their taxids, one per line): gtdb_to_protal_db.py --from_db does so for a
training database that leaves species out (derive()).

Usage:
  gene_neighbours.py --db <converted folder> --genome_table <table> [-t 8] [--max_gap 3000]
      [--ranks family,order,class,phylum,domain] [--no_species_lines]
  gene_neighbours.py --db <converted folder> --from_positions gene_positions.tsv [--exclude_taxa FILE]
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
POSITIONS_FILE = "gene_positions.tsv"
HEADER = "clade\tgene\tend\tpartner\tpartner_end\tspecies\tinformative\tgap_median\tgap_min\tgap_max"
POSITIONS_HEADER = "accession\ttaxid\tcontig\tcontig_length\tcircular\tgene\tstart\tend\tstrand\tplaced\tkmer_share"
COMPLEMENT = bytes.maketrans(b"ACGTN", b"TGCAN")
SEPARATOR = b"|"  # between contigs: no gene matches across it
KMER, STRIDE = 24, 16  # GenomeIndex: the genome's k-mers at every STRIDE-th position
MIN_KMER_SHARE = 0.1  # a trace places a gene with this share of its k-mers on one diagonal band
MIN_KMER_HITS = 3     # and at least this many
TRACE_BAND = 32       # diagonals this close are one place (a strain's indels shift its gene's diagonal)
PRIOR_SPECIES = 3     # protal's kPriorSpecies: a clade with a gene end informative in fewer species counts less
                      # there than the clades above it, whose share its own is smoothed towards (GeneNeighbours.h)
MIN_INFORMATIVE = 5   # protal's kMinInformative: a top clade with fewer judges no pairing unlikely
EXPECTED_SHARE = 0.2  # protal's kExpectedShare: a pairing of this smoothed share or more is expected
REPEAT_SHARE = 0.5    # a second place with this share of the best's hits: the gene is there twice


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


def accession_of(text):
    m = ACCESSION_RE.search(text)
    return m.group(1) if m else text


def read_genome_table(path):
    """accession -> FASTA path, from a genome table (accession, taxonomy, path[, length]; a header line is
    skipped)."""
    paths = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or line.startswith("#") or f[0] == "accession":
                continue
            paths.setdefault(accession_of(f[0]), f[2])
    return paths


def read_genome_taxa(path):
    """genome2tiid.tsv (accession, species taxid, ...) -> {accession: taxid}: every genome of the release."""
    taxa = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[1].isdigit():
                taxa[accession_of(f[0])] = int(f[1])
    return taxa


def read_map(path):
    """reference.map -> {taxid: [(gene id, start byte, end byte)]}."""
    genes = collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) >= 4:
                genes[int(f[0])].append((int(f[1]), int(f[2]), int(f[3])))
    return genes


# FASTA headers of sequences that are a whole replicon: NCBI's "..., complete genome" (chromosomes) and
# "..., complete sequence" (plasmids, chromosomes), and the topology tag of submissions. Such a sequence is
# taken as circular: its last gene faces its first across the origin, and none of its gene ends is near an end.
CIRCULAR_MARKS = (b"complete genome", b"complete sequence", b"topology=circular", b"[circular]")


def read_contigs(path):
    """A genome FASTA (plain or gzip) -> [(name, sequence as upper-case bytes, circular)], circular if the
    header marks the sequence a whole replicon (CIRCULAR_MARKS)."""
    with open(path, "rb") as fh:
        data = fh.read()
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    contigs = []
    for record in (b"\n" + data).split(b"\n>")[1:]:
        header, _, body = record.partition(b"\n")
        words = header.split()
        circular = any(mark in header.lower() for mark in CIRCULAR_MARKS)
        contigs.append((words[0].decode() if words else "", b"".join(body.split()).upper(), circular))
    return contigs


class GenomeIndex:
    """A genome's KMER-mers at every STRIDE-th position. An occurrence of a sequence at any position p holds the
    indexed k-mer at the next multiple of STRIDE as its own k-mer at offset j = -p mod STRIDE, if the sequence
    has KMER + STRIDE - 1 bases or more: so occurrences() finds every exact occurrence, as bytes.find does,
    without a scan of the genome per gene (1.45 s per GTDB-sized genome with find), and trace() the place of a
    sequence that differs from the genome's by scattered mutations."""

    def __init__(self, genome):
        self.genome = genome
        self.index = {}
        for p in range(0, len(genome) - KMER + 1, STRIDE):
            self.index.setdefault(genome[p:p + KMER], []).append(p)

    def occurrences(self, seq):
        """The positions where seq starts in the genome, ascending (overlapping occurrences too)."""
        genome, n = self.genome, len(seq)
        if n < KMER + STRIDE - 1:
            found, pos = [], genome.find(seq)
            while pos >= 0:
                found.append(pos)
                pos = genome.find(seq, pos + 1)
            return found
        found = set()
        for j in range(STRIDE):
            for p in self.index.get(seq[j:j + KMER], ()):
                if p >= j and genome[p - j:p - j + n] == seq:
                    found.add(p - j)
        return sorted(found)

    def trace(self, seq):
        """Where seq lies by its k-mer trace: (start in the genome, '+' or '-', share of its k-mers that can hit
        there and do, whether a second place has REPEAT_SHARE of the hits), or None. start is the most voted
        diagonal of the best band, the gene's first base if it has no indel before its first hit."""
        if len(seq) < KMER:
            return None
        get = self.index.get
        places = []  # (hits, strand, diagonal)
        for strand, s in (("+", seq), ("-", seq.translate(COMPLEMENT)[::-1])):
            votes = collections.Counter()
            for j in range(len(s) - KMER + 1):
                for p in get(s[j:j + KMER], ()):
                    votes[p - j] += 1
            band, top = None, None
            for d in sorted(votes):
                if band is not None and d - band[1] <= TRACE_BAND:
                    band[1], band[2] = d, band[2] + votes[d]
                    if votes[d] > votes[band[3]]:
                        band[3] = d
                else:
                    if band is not None:
                        places.append((band[2], strand, band[3]))
                    band = [d, d, votes[d], d]
            if band is not None:
                places.append((band[2], strand, band[3]))
        if not places:
            return None
        places.sort(key=lambda p: (-p[0], p[1], p[2]))
        hits, strand, start = places[0]
        possible = -(-(len(seq) - KMER + 1) // STRIDE)
        share = hits / possible
        if hits < MIN_KMER_HITS or share < MIN_KMER_SHARE:
            return None
        repeated = len(places) > 1 and places[1][0] >= REPEAT_SHARE * hits
        return start, strand, min(1.0, share), repeated


def place_genome(contigs, genes):
    """genes [(gene id, sequence)] in a genome's contigs [(name, sequence, circular)] ->
    ([(gene, contig index, start, end, strand, placed, kmer_share)] with start, end 0-based and end exclusive,
    genes found more than once)."""
    genome = SEPARATOR.join(seq for _, seq, _ in contigs)
    index = GenomeIndex(genome)
    offsets, at = [], 0
    for _, seq, _ in contigs:
        offsets.append(at)
        at += len(seq) + 1
    placements, repeated = [], 0
    for gene, seq in genes:
        strand, found = "+", index.occurrences(seq)
        if not found:
            strand, found = "-", index.occurrences(seq.translate(COMPLEMENT)[::-1])
        if found:
            start, placed, share = found[0], "exact", 1.0
            repeated += len(found) > 1
        else:
            traced = index.trace(seq)
            if traced is None:
                continue
            start, strand, share, twice = traced
            placed = "trace"
            repeated += twice
        # The contig of the gene's middle: a strain's copy cut at a contig's start has its diagonal before it.
        contig = max(0, bisect.bisect_right(offsets, start + len(seq) // 2) - 1)
        length = len(contigs[contig][1])
        local = start - offsets[contig]
        placements.append((gene, contig, max(0, local), min(length, local + len(seq)), strand, placed, share))
    return placements, repeated


def place_species(item):
    """One species: its database genes placed in each of its genomes. -> (taxid, [(accession, contigs
    [(name, length, circular)] or None, placements or why the genome could not be read, genes looked for,
    repeated)])."""
    taxid, genomes, reference, genes = item
    sequences = []
    with open(reference, "rb") as fh:
        for gene, start, _end in genes:
            fh.seek(start)
            seq = fh.readline().strip().upper()
            if seq:
                sequences.append((gene, seq))
    out = []
    for accession, path in genomes:
        try:
            contigs = read_contigs(path)
        except (OSError, EOFError, ValueError) as error:
            out.append((accession, None, str(error), len(sequences), 0))
            continue
        placements, repeated = place_genome(contigs, sequences)
        out.append((accession, [(n, len(s), c) for n, s, c in contigs], placements, len(sequences), repeated))
    return taxid, out


def neighbour_ends(placements, max_gap):
    """The informative ends of a genome's placed genes, placements [(gene, contig, contig length, circular,
    start, end, strand)] (0-based, end exclusive): [(gene, end, partner, partner end, gap)], partner 0 for no
    marker within max_gap. On a circular contig the last gene faces the first across the origin, and every end
    is informative."""
    by_contig = collections.defaultdict(list)
    for gene, contig, length, circular, start, end, strand in placements:
        by_contig[contig].append((start, end, gene, strand, length, circular))
    ends = []
    for genes in by_contig.values():
        genes.sort()
        n = len(genes)
        for i, (start, end, gene, strand, length, circular) in enumerate(genes):
            right_end = 3 if strand == "+" else 5   # the gene's end at its right in the genome
            left_end = 5 if strand == "+" else 3
            if i + 1 < n:
                after, gap = genes[i + 1], genes[i + 1][0] - end
            elif circular and n > 1:
                after, gap = genes[0], length - end + genes[0][0]
            else:
                after, gap = None, 0
            if after is not None and gap <= max_gap:
                ends.append((gene, right_end, after[2], 5 if after[3] == "+" else 3, gap))
            elif after is not None or circular or length - end >= max_gap:
                ends.append((gene, right_end, 0, 0, 0))
            if i > 0:
                before, gap = genes[i - 1], start - genes[i - 1][1]
            elif circular and n > 1:
                before, gap = genes[-1], start + length - genes[-1][1]
            else:
                before, gap = None, 0
            if before is not None and gap <= max_gap:
                ends.append((gene, left_end, before[2], 3 if before[3] == "+" else 5, gap))
            elif before is not None or circular or start >= max_gap:
                ends.append((gene, left_end, 0, 0, 0))
    return ends


def species_ends(genomes, representative, max_gap):
    """One species' gene ends from its genomes {accession: placements} (neighbour_ends' input): {(gene, end):
    (partner, partner end, gap)}, the partner most of its genomes informative there show (the representative's
    on a tie, else the smallest), the gap their median."""
    votes = collections.defaultdict(list)  # (gene, end) -> [(partner, partner end, gap, from the representative)]
    for accession, placements in genomes.items():
        for gene, end, partner, partner_end, gap in neighbour_ends(placements, max_gap):
            votes[(gene, end)].append((partner, partner_end, gap, accession == representative))
    chosen = {}
    for key, seen in votes.items():
        counts = collections.Counter((p, pe) for p, pe, _, _ in seen)
        top = max(counts.values())
        tied = sorted(k for k, n in counts.items() if n == top)
        rep = [(p, pe) for p, pe, _, is_rep in seen if is_rep and (p, pe) in tied]
        partner = rep[0] if rep else tied[0]
        gaps = [g for p, pe, g, _ in seen if (p, pe) == partner]
        chosen[key] = (partner[0], partner[1], int(round(statistics.median(gaps))) if partner[0] else 0)
    return chosen


def species_partners(genomes, max_gap):
    """Every partner one species' genomes {accession: placements} show at each gene end: {(gene, end): {(partner,
    partner end): gap}}, the gap the median of the genomes that show it (0 for no partner)."""
    seen = collections.defaultdict(lambda: collections.defaultdict(list))
    for placements in genomes.values():
        for gene, end, partner, partner_end, gap in neighbour_ends(placements, max_gap):
            seen[(gene, end)][(partner, partner_end)].append(gap)
    return {key: {partner: int(round(statistics.median(gaps))) if partner[0] else 0 for partner, gaps in partners.items()}
            for key, partners in seen.items()}


def smoothed_share(counts, informative, chain, gene, end, partner):
    """protal's share of a pairing in a species (gene_neighbours::Table::Assess): over the species' clades from the
    top down (chain: nearest first), the top clade with data on the end has its own share, each one below (species +
    PRIOR_SPECIES x its parent's share) / (informative + PRIOR_SPECIES). -> (share, whether the top clade has the end
    informative in MIN_INFORMATIVE species or more), or (None, False) if no clade has data on the end."""
    share, populated = None, False
    for clade in reversed(chain):
        key = (clade, gene, end)
        n_informative = informative.get(key, 0)
        if n_informative == 0:
            continue
        n = len(counts[key].get(partner, ())) if key in counts else 0
        if share is None:
            share, populated = n / n_informative, n_informative >= MIN_INFORMATIVE
        else:
            share = (n + PRIOR_SPECIES * share) / (n_informative + PRIOR_SPECIES)
    return share, populated


def count_clades(species, nodes, ranks, reps_of, max_gap, species_lines=True):
    """The clades' lines from species {taxid: {accession: placements}}: {(clade, gene, end): {(partner, partner
    end): [gaps of its species]}}, {(clade, gene, end): informative species}, and the species' informative ends'
    gaps (None: no marker within max_gap).

    With species_lines (and no species among the ranks), also lines of the species itself (clade: its taxid, each
    partner 1 of 1 species) at each end where a partner that any of its genomes shows there would be no expected
    neighbour by its clades alone (smoothed_share below EXPECTED_SHARE, in clades with enough species to judge): a
    species whose gene order differs from its family's there, or whose strains differ (a strain's partner, not that
    of most of its genomes, which the clades count). protal reads the species as the nearest clade of its lineage, so
    its own neighbours are expected (a share of at least 1/4 with PRIOR_SPECIES 3) while its clades' stay so (3/4 of
    their share): the genes next to each other on the reads of any of its genomes are judged its neighbours, also
    where its family has them otherwise."""
    counts = collections.defaultdict(lambda: collections.defaultdict(list))
    informative = collections.Counter()
    ends = []
    own = {}
    for taxid in sorted(species):
        clades = ancestors(nodes, taxid, ranks).values()
        own[taxid] = species_ends(species[taxid], reps_of.get(taxid), max_gap)
        for (gene, end), (partner, partner_end, gap) in own[taxid].items():
            ends.append(gap if partner else None)
            for clade in clades:
                counts[(clade, gene, end)][(partner, partner_end)].append(gap)
                informative[(clade, gene, end)] += 1
    if species_lines and "species" not in ranks:
        for taxid in sorted(own):
            chain = list(ancestors(nodes, taxid, ranks).values())  # nearest first
            for (gene, end), partners in sorted(species_partners(species[taxid], max_gap).items()):
                for partner, gap in sorted(partners.items()):
                    share, populated = smoothed_share(counts, informative, chain, gene, end, partner)
                    if share is not None and populated and share < EXPECTED_SHARE:
                        counts[(taxid, gene, end)][partner].append(gap)
                        informative[(taxid, gene, end)] = 1
    return counts, informative, ends


def write_table(path, counts, informative, comment):
    """gene_neighbours.tsv from count_clades; the number of lines."""
    rows = 0
    with open(path + ".partial", "w", newline="\n") as fh:
        fh.write(f"# protal gene neighbours: {comment}\n")
        fh.write(HEADER + "\n")
        for key in sorted(counts):
            clade, gene, end = key
            for (partner, partner_end), gaps in sorted(counts[key].items()):
                median, low, high = (int(round(statistics.median(gaps))), min(gaps), max(gaps)) if partner else (0, 0, 0)
                fh.write(f"{clade}\t{gene}\t{end}\t{partner}\t{partner_end}\t{len(gaps)}\t{informative[key]}\t"
                         f"{median}\t{low}\t{high}\n")
                rows += 1
    os.replace(path + ".partial", path)
    return rows


def read_positions(path, exclude=frozenset()):
    """gene_positions.tsv -> ({species taxid: {accession: placements}} (neighbour_ends' input), {name: value}
    of its comment line), without the species of `exclude`."""
    species = collections.defaultdict(lambda: collections.defaultdict(list))
    settings = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                settings.update(w.split("=", 1) for w in line[1:].split() if "=" in w)
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] == "accession" or len(f) < 9:
                continue
            taxid = int(f[1])
            if taxid in exclude:
                continue
            species[taxid][f[0]].append((int(f[5]), f[2], int(f[3]), f[4] == "1", int(f[6]) - 1, int(f[7]), f[8]))
    return species, settings


def derive(positions, taxonomy, output, exclude=frozenset(), positions_out=None):
    """gene_neighbours.tsv (output) from a gene_positions.tsv without the species of `exclude` (taxids), with the
    max_gap, ranks and species lines its comment names (species lines if it does not say); and, with positions_out,
    the positions file without them. -> (lines, genomes, species)."""
    nodes, reps = read_taxonomy(taxonomy)
    species, settings = read_positions(positions, exclude)
    max_gap = int(settings.get("max_gap", 3000))
    ranks = settings.get("ranks", "family,order,class,phylum,domain").split(",")
    species_lines = settings.get("species_lines", "1") != "0"
    counts, informative, _ = count_clades(species, nodes, ranks, {t: a for a, t in reps.items()}, max_gap, species_lines)
    genomes = sum(len(g) for g in species.values())
    rows = write_table(output, counts, informative, f"genomes={genomes} species={len(species)} max_gap={max_gap} "
                                                     f"ranks={','.join(ranks)} species_lines={int(species_lines)}")
    if positions_out:
        with open(positions) as fin, open(positions_out + ".partial", "w", newline="\n") as fout:
            for line in fin:
                f = line.split("\t")
                if line.startswith("#") or f[0] == "accession" or len(f) < 2 or int(f[1]) not in exclude:
                    fout.write(line)
        os.replace(positions_out + ".partial", positions_out)
    return rows, genomes, len(species)


def summarise(counts, informative, nodes, ranks):
    """How alike each rank's clades are: per rank, its clades, the species-weighted share of a gene end's most
    common partner, and how many of its gene ends lean mostly on the rank above (PRIOR_SPECIES)."""
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
            sparse = sum(1 for w in weights if w < PRIOR_SPECIES)
            print(f"  {rank}: {len(clades)} clades; the most common partner of a gene end is that of {100 * mean:.1f}% "
                  f"of the clade's species (weighted by species); {sparse} of its {len(weights)} gene ends "
                  f"informative in fewer than {PRIOR_SPECIES} species (they lean mostly on the rank above)")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--db", required=True, help="converted database folder (reference.fna, reference.map, "
                                                "internal_taxonomy.dmp, genome2tiid.tsv; before protal --build packs them)")
    ap.add_argument("--genome_table", help="accession, GTDB taxonomy, FASTA path (a header line is skipped)")
    ap.add_argument("--output", help=f"default: <db>/{FILE_NAME}")
    ap.add_argument("--positions", help=f"where each gene was found in each genome; default: <db>/{POSITIONS_FILE}")
    ap.add_argument("--max_gap", type=int, default=3000, help="the farthest a neighbour counts, in bases (default 3000)")
    ap.add_argument("--ranks", default="family,order,class,phylum,domain",
                    help="the ranks of the clades to count for (default family,order,class,phylum,domain)")
    ap.add_argument("--no_species_lines", action="store_true",
                    help="no lines of a species' own where its gene order differs from its clades' (count_clades)")
    ap.add_argument("--min_placed", type=float, default=0.8,
                    help="a genome counts only if this share of its species' genes is placed in it (default 0.8)")
    ap.add_argument("--from_positions", help="derive the table from this gene_positions.tsv instead of the genomes")
    ap.add_argument("--exclude_taxa", help="with --from_positions: species taxids to leave out, one per line")
    ap.add_argument("-t", "--threads", type=int, default=1)
    args = ap.parse_args()

    output = args.output or os.path.join(args.db, FILE_NAME)
    taxonomy = os.path.join(args.db, "internal_taxonomy.dmp")
    if args.from_positions:
        exclude = set()
        if args.exclude_taxa:
            with open(args.exclude_taxa) as fh:
                exclude = {int(w) for w in fh.read().split()}
        rows, genomes, n_species = derive(args.from_positions, taxonomy, output, exclude)
        print(f"Wrote {output}: {rows} lines, from {genomes} genomes of {n_species} species in {args.from_positions}")
        return
    if not args.genome_table:
        ap.error("give --genome_table (or --from_positions)")
    reference = os.path.join(args.db, "reference.fna")
    for name in ("reference.fna", "reference.map", "internal_taxonomy.dmp", "genome2tiid.tsv"):
        if not os.path.isfile(os.path.join(args.db, name)):
            sys.exit(f"{args.db} has no {name}: run this on the converted folder before protal --build packs it")
    ranks = [r.strip() for r in args.ranks.split(",") if r.strip()]
    nodes, reps = read_taxonomy(taxonomy)
    reps_of = {t: a for a, t in reps.items()}
    genes = read_map(os.path.join(args.db, "reference.map"))
    taxa = read_genome_taxa(os.path.join(args.db, "genome2tiid.tsv"))
    by_species = collections.defaultdict(list)
    for accession, path in sorted(read_genome_table(args.genome_table).items()):
        taxid = taxa.get(accession)
        if taxid is not None and genes.get(taxid):
            by_species[taxid].append((accession, path))
    work = [(taxid, sorted(by_species[taxid], key=lambda g: g[0] != reps_of.get(taxid)), reference, genes[taxid])
            for taxid in sorted(by_species)]
    with_rep = sum(1 for taxid, genomes, _, _ in work if any(a == reps_of.get(taxid) for a, _ in genomes))
    print(f"{sum(len(w[1]) for w in work)} genomes of {len(work)} species in {args.genome_table} "
          f"({with_rep} with their representative genome)", flush=True)

    species = collections.defaultdict(dict)  # taxid -> accession -> placements
    used = unreadable = sparse = looked = exact = traced = repeated = by_header = single_reps = 0
    positions_path = args.positions or os.path.join(args.db, POSITIONS_FILE)
    with open(positions_path + ".partial", "w", newline="\n") as positions, \
            multiprocessing.Pool(max(1, args.threads)) as pool:
        positions.write(f"# protal gene positions: max_gap={args.max_gap} ranks={','.join(ranks)} kmer={KMER} "
                        f"stride={STRIDE} min_kmer_share={MIN_KMER_SHARE} min_placed={args.min_placed} "
                        f"species_lines={int(not args.no_species_lines)}\n")
        positions.write(POSITIONS_HEADER + "\n")
        for taxid, results in pool.imap(place_species, work, chunksize=2):
            for accession, contigs, placements, n_genes, n_repeated in results:
                if contigs is None:
                    unreadable += 1
                    sys.stderr.write(f"Cannot read the genome of {accession}: {placements}\n")
                    continue
                looked += n_genes
                if len(placements) < args.min_placed * n_genes:
                    sparse += 1
                    continue
                used += 1
                repeated += n_repeated
                by_header += any(c for _, _, c in contigs)
                if len(contigs) == 1 and accession == reps_of.get(taxid) and not contigs[0][2]:
                    contigs = [(contigs[0][0], contigs[0][1], True)]  # a representative in one sequence: closed
                    single_reps += 1
                rows = []
                for gene, contig, start, end, strand, placed, share in placements:
                    exact += placed == "exact"
                    traced += placed == "trace"
                    name, length, circular = contigs[contig]
                    positions.write(f"{accession}\t{taxid}\t{name}\t{length}\t{int(circular)}\t{gene}\t{start + 1}\t{end}\t"
                                    f"{strand}\t{placed}\t{share:.3f}\n")
                    rows.append((gene, name, length, circular, start, end, strand))
                species[taxid][accession] = rows
    os.replace(positions_path + ".partial", positions_path)

    species_lines = not args.no_species_lines
    counts, informative, ends = count_clades(species, nodes, ranks, reps_of, args.max_gap, species_lines)
    rows = write_table(output, counts, informative, f"genomes={used} species={len(species)} max_gap={args.max_gap} "
                                                     f"ranks={','.join(ranks)} species_lines={int(species_lines)}")

    # The snapshot: how the genes were placed, how much adjacency there is, how alike it is within each rank.
    print(f"Genes placed: {exact} exactly and {traced} by their k-mer trace in the {used} genomes used (of "
          f"{len(species)} species; circular: {by_header} with a whole replicon by its header, {single_reps} "
          f"representatives in one sequence); "
          f"{sparse} genomes with fewer than {args.min_placed:.0%} of their species' genes placed, {unreadable} "
          f"unreadable; {repeated} genes found more than once (the first or best place taken)")
    gaps = [g for g in ends if g is not None]
    if ends:
        spread = ""
        if len(gaps) > 1:
            q = statistics.quantiles(gaps, n=10)
            spread = f"; gaps: median {statistics.median(gaps):.0f}, 10-90% {q[0]:.0f} to {q[-1]:.0f}"
        print(f"Gene ends: {len(ends)} informative in the species, {len(gaps)} ({100 * len(gaps) / len(ends):.1f}%) "
              f"with a marker within {args.max_gap} bases{spread}")
    summarise(counts, informative, nodes, ranks)
    if species_lines:
        own = collections.Counter(clade for clade, _, _ in counts if nodes[clade][1] == "species")
        print(f"  species: {len(own)} of the {len(species)} species have lines of their own, at {sum(own.values())} gene "
              f"ends where their partner would be no expected neighbour by their clades alone")
    print(f"Wrote {output}: {rows} lines; {positions_path}: the genes of {used} genomes")


if __name__ == "__main__":
    main()
