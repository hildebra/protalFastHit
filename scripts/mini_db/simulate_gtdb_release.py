#!/usr/bin/env python3
"""simulate_gtdb_release.py - write a sparse, synthetic GTDB release for local
protal testing.

The output mirrors the layout of a real (extracted) GTDB release, populated with
a handful of made-up species (default: 3). Sequences are random but evolved down
the taxonomy, so related species share diverged copies of the same marker genes:

  <outdir>/
    VERSION.txt
    {bac120,ar53}_taxonomy_r<R>.tsv
    {bac120,ar53}_metadata_r<R>.tsv.gz
    genomic_files_reps/{bac120,ar53}_marker_genes_reps_r<R>/{fna,faa}/<set>_<marker>.{fna,faa}
    genomic_files_reps/gtdb_genomes_reps_r<R>/database/GCF/999/.../<acc>_genomic.fna.gz
    genomic_files_all/{bac120,ar53}_marker_genes_all_r<R>/{fna,faa}/<set>_<marker>.{fna,faa}
    simulation/              (not part of GTDB: truth data for tests)
      genomes.tsv            genome table for simulate_metagenomes (all genomes)
      marker_positions.tsv   where every marker gene sits in its genome
      genomes_nonreps/       FASTA of the non-representative genomes

Every species gets one representative genome (RS_GCF_999SSS001.1) plus
--genomes_per_species - 1 non-representative strains (GB_GCA_999SSSGGG.1).
NCBI has not issued accessions in the 999xxxxxx range, so none of these collide
with a real genome.

The marker universe comes from markers_r226.tsv (GTDB r226 bac120 + ar53 marker
info). Species in d__Archaea carry the ar53 markers, all others bac120.

Usage:
  simulate_gtdb_release.py --outdir DIR [--release 226] [--seed 42]
      [--lineages FILE] [--genomes_per_species 3] [--genome_length 150000]
      [--contigs 3] [--marker_loss 0.02] [--strain_divergence 0.005]

--lineages: one GTDB lineage per line (d__...;s__...); '#' lines are comments.
"""

import argparse
import contextlib
import gzip
import io
import math
import os
import random
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Per-site substitution probability on the branch leading INTO a node of the
# given rank. Two congeneric species differ at ~8% of marker sites, two strains
# of one species at ~1% (2 x --strain_divergence), different phyla are
# essentially unrelated.
BRANCH = {"d": 0.15, "p": 0.10, "c": 0.04, "o": 0.03, "f": 0.03, "g": 0.03, "s": 0.035}
# Fraction of a branch's mutation events that are codon indels (rest: substitutions).
INDEL_FRACTION = 0.02
MAX_INDELS_PER_BRANCH = 3

DEFAULT_LINEAGES = [
    "d__Bacteria;p__Simulatota;c__Simulatia;o__Simulales;f__Simulaceae;g__Mockella;s__Mockella alpha",
    "d__Bacteria;p__Simulatota;c__Simulatia;o__Simulales;f__Simulaceae;g__Mockella;s__Mockella beta",
    "d__Bacteria;p__Fictota;c__Fictia;o__Fictales;f__Fictaceae;g__Fakibacter;s__Fakibacter gamma",
]

BASES = "ACGT"
STOP = {"TAA", "TAG", "TGA"}
COMPLEMENT = str.maketrans("ACGT", "TGCA")
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: _AA[16 * i + 4 * j + k]
         for i, a in enumerate("TCAG") for j, b in enumerate("TCAG") for k, c in enumerate("TCAG")}
SENSE_CODONS = [c for c in sorted(CODON) if c not in STOP]


@contextlib.contextmanager
def gzip_text(path):
    """Text writer for a .gz file whose bytes depend only on the content
    (gzip.open would stamp the current time into the header)."""
    with open(path, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz, \
         io.TextIOWrapper(gz, encoding="ascii", newline="\n") as fh:
        yield fh


def read_lineages(path):
    with open(path) as fh:
        return [l.strip() for l in fh if l.strip() and not l.lstrip().startswith("#")]


def read_markers(path):
    """-> ({marker: length_aa}, {set: [marker, ...]}) from markers_r226.tsv."""
    lengths, by_set = {}, {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("set\t"):
                continue
            mset, marker, _name, length = line.rstrip("\n").split("\t")[:4]
            lengths.setdefault(marker, int(length))
            by_set.setdefault(mset, []).append(marker)
    if "bac120" not in by_set:
        sys.exit(f"Marker table {path} has no bac120 markers")
    by_set.setdefault("ar53", [])
    return lengths, by_set


def random_cds(rng, length_aa):
    """ATG + random sense codons + stop, 1.0-1.15x the HMM length."""
    codons = int(length_aa * (1 + rng.random() * 0.15))
    return "ATG" + "".join(rng.choices(SENSE_CODONS, k=codons)) + "TAA"


def random_dna(rng, length, gc):
    at = (1 - gc) / 2
    return "".join(rng.choices(BASES, weights=(at, gc / 2, gc / 2, at), k=length))


def poisson(rng, lam, cap):
    """Knuth's Poisson sampler, capped (lam is small here)."""
    limit, n, prod = math.exp(-lam), 0, rng.random()
    while prod > limit and n < cap:
        n += 1
        prod *= rng.random()
    return n


def mutate(rng, seq, d, coding=True):
    """Bernoulli(d) substitutions per site (first/last codon kept). For coding
    sequences (in frame 0), substitutions never create a stop codon and a share
    of the events are codon indels instead."""
    if d <= 0:
        return seq
    s = list(seq)
    p_sub = d * (1 - INDEL_FRACTION) if coding else d
    log_q = math.log(1 - p_sub)
    # Geometric skipping: jump straight to the next mutated site.
    pos = 3 + int(math.log(1 - rng.random()) / log_q)
    while pos < len(s) - 3:
        alts = BASES.replace(s[pos], "")
        if coding:
            c = pos - pos % 3
            alts = [b for b in alts if "".join(s[c:pos] + [b] + s[pos + 1:c + 3]) not in STOP]
        if alts:
            s[pos] = rng.choice(alts)
        pos += 1 + int(math.log(1 - rng.random()) / log_q)
    if coding:
        for _ in range(poisson(rng, len(s) * d * INDEL_FRACTION, MAX_INDELS_PER_BRANCH)):
            codons = len(s) // 3
            if codons < 4:
                break
            at = 3 * rng.randint(1, codons - 2)
            if rng.random() < 0.5:
                del s[at:at + 3]
            else:
                s[at:at] = rng.choice(SENSE_CODONS)
    return "".join(s)


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def translate(seq):
    aa = "".join(CODON.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))
    return aa[:-1] if aa.endswith("*") else aa


class Simulator:
    def __init__(self, args, lineages, lengths, markers_by_set):
        self.args = args
        self.rng = random.Random(args.seed)
        self.markers_by_set = markers_by_set
        # Ancestral (root) sequence for every marker of the union of both sets.
        self.root_seq = {m: random_cds(self.rng, lengths[m]) for m in sorted(lengths)}
        self.node_seq = {}  # (lineage prefix, marker) -> sequence, shared by descendants
        self.species = [self._make_species(i + 1, l) for i, l in enumerate(lineages)]
        self.genomes = [g for sp in self.species for g in self._make_genomes(sp)]

    def _make_species(self, index, lineage):
        ranks = lineage.split(";")
        if len(ranks) != 7 or [r[:3] for r in ranks] != ["d__", "p__", "c__", "o__", "f__", "g__", "s__"]:
            sys.exit(f"Lineage needs the 7 ranks d__..s__: {lineage}")
        mset = "ar53" if ranks[0] == "d__Archaea" else "bac120"
        markers = self.markers_by_set[mset]
        if not markers:
            sys.exit(f"No {mset} markers in the marker table for {lineage}")
        return {
            "index": index,
            "lineage": lineage,
            "name": ranks[6][3:],
            "set": mset,
            "markers": markers,
            "seqs": {m: self._evolve(ranks, m) for m in markers},
            "gc": 0.38 + self.rng.random() * 0.26,
        }

    def _evolve(self, ranks, marker):
        """Marker sequence at the species node, evolving (and caching) every ancestor."""
        seq = self.root_seq[marker]
        for depth, rank in enumerate(ranks):
            key = (";".join(ranks[:depth + 1]), marker)
            if key not in self.node_seq:
                self.node_seq[key] = mutate(self.rng, seq, BRANCH[rank[0]])
            seq = self.node_seq[key]
        return seq

    def _make_genomes(self, sp):
        a = self.args
        order = list(sp["markers"])
        self.rng.shuffle(order)
        strand = {m: "+" if self.rng.random() < 0.5 else "-" for m in order}
        spacer = a.genome_length // (len(order) + 1)
        background = [random_dna(self.rng, spacer, sp["gc"]) for _ in range(len(order) + 1)]

        for g in range(1, a.genomes_per_species + 1):
            is_rep = g == 1
            acc = f"GCF_999{sp['index']:03d}001.1" if is_rep else f"GCA_999{sp['index']:03d}{g:03d}.1"
            genes = {}
            for m in order:
                if self.rng.random() >= a.marker_loss:
                    genes[m] = mutate(self.rng, sp["seqs"][m], a.strain_divergence)
            bg = [mutate(self.rng, s, a.strain_divergence, coding=False) for s in background]
            contigs, positions = self._assemble(acc, order, genes, strand, bg)
            yield {
                "accession": acc,
                "gtdb_acc": ("RS_" if is_rep else "GB_") + acc,
                "rep_gtdb_acc": f"RS_GCF_999{sp['index']:03d}001.1",
                "species": sp,
                "is_rep": is_rep,
                "strain": f"SIM-{g}",
                "genes": genes,
                "contigs": contigs,
                "positions": positions,
            }

    def _assemble(self, acc, order, genes, strand, bg):
        """Interleave background spacers and genes; split into contigs between genes."""
        present = [m for m in order if m in genes]
        n_contigs = max(1, min(self.args.contigs, len(present)))
        # n contigs need n - 1 breaks, chosen among the len(present) - 1 gene gaps.
        breaks = set(self.rng.sample(range(len(present) - 1), n_contigs - 1)) if n_contigs > 1 else set()
        contigs, positions, cur = [], [], [bg[0]]
        cur_len = len(bg[0])
        for i, m in enumerate(present):
            seq = genes[m] if strand[m] == "+" else revcomp(genes[m])
            contig = f"{acc}_contig{len(contigs) + 1}"
            positions.append((m, contig, cur_len + 1, cur_len + len(seq), strand[m]))
            cur += [seq, bg[i + 1]]
            cur_len += len(seq) + len(bg[i + 1])
            if i in breaks:
                contigs.append("".join(cur))
                cur, cur_len = [], 0
        contigs.append("".join(cur))
        return contigs, positions


def genome_path(root, acc):
    """GTDB layout: database/GCF/999/001/001/GCF_999001001.1_genomic.fna.gz"""
    db, digits = acc[:3], acc[4:13]
    return os.path.join(root, db, digits[0:3], digits[3:6], digits[6:9], f"{acc}_genomic.fna.gz")


def write_taxonomy(path, genomes):
    with open(path, "w", newline="\n") as fh:
        for g in genomes:
            fh.write(f"{g['gtdb_acc']}\t{g['species']['lineage']}\n")


def write_metadata(path, genomes):
    cols = ["accession", "ambiguous_bases", "checkm2_completeness", "checkm2_contamination",
            "contig_count", "gc_percentage", "genome_size", "gtdb_genome_representative",
            "gtdb_representative", "gtdb_taxonomy", "ncbi_organism_name"]
    with gzip_text(path) as fh:
        fh.write("\t".join(cols) + "\n")
        for g in genomes:
            seq = "".join(g["contigs"])
            gc = 100 * (seq.count("G") + seq.count("C")) / len(seq)
            fh.write("\t".join([
                g["gtdb_acc"], "0", "100.00", "0.00", str(len(g["contigs"])), f"{gc:.2f}",
                str(len(seq)), g["rep_gtdb_acc"], "t" if g["is_rep"] else "f",
                g["species"]["lineage"], f"{g['species']['name']} strain {g['strain']}",
            ]) + "\n")


def write_marker_files(directory, mset, markers, genomes):
    os.makedirs(os.path.join(directory, "fna"), exist_ok=True)
    os.makedirs(os.path.join(directory, "faa"), exist_ok=True)
    for m in markers:
        with open(os.path.join(directory, "fna", f"{mset}_{m}.fna"), "w", newline="\n") as fna, \
             open(os.path.join(directory, "faa", f"{mset}_{m}.faa"), "w", newline="\n") as faa:
            for g in genomes:
                seq = g["genes"].get(m)
                if seq is None:
                    continue
                fna.write(f">{g['gtdb_acc']}\n{seq}\n")
                faa.write(f">{g['gtdb_acc']}\n{translate(seq)}\n")


def write_genome(path, g):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    name, strain = g["species"]["name"], g["strain"]
    with gzip_text(path) as fh:
        for i, seq in enumerate(g["contigs"], 1):
            fh.write(f">{g['accession']}_contig{i} {name} strain {strain}, synthetic contig {i}\n")
            for j in range(0, len(seq), 80):
                fh.write(seq[j:j + 80] + "\n")
    return sum(len(s) for s in g["contigs"])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--outdir", required=True, help="release directory to create")
    ap.add_argument("--release", default="226", help="GTDB release number used in file names")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--markers", default=os.path.join(SCRIPT_DIR, "markers_r226.tsv"))
    ap.add_argument("--lineages", help="file with one GTDB lineage per line (default: 3 built-in species)")
    ap.add_argument("--genomes_per_species", type=int, default=3,
                    help="genomes per species; the first is the GTDB representative")
    ap.add_argument("--genome_length", type=int, default=150_000,
                    help="background (non-marker) DNA per genome")
    ap.add_argument("--contigs", type=int, default=3)
    ap.add_argument("--marker_loss", type=float, default=0.02,
                    help="probability that a genome lacks a given marker")
    ap.add_argument("--strain_divergence", type=float, default=0.005,
                    help="per-site substitution rate of each genome relative to its species")
    args = ap.parse_args()

    if not 1 <= args.genomes_per_species <= 999:
        sys.exit("--genomes_per_species must be in 1..999")
    lineages = read_lineages(args.lineages) if args.lineages else DEFAULT_LINEAGES
    if not 1 <= len(lineages) <= 999:
        sys.exit("Need 1..999 lineages")

    lengths, markers_by_set = read_markers(args.markers)
    sim = Simulator(args, lineages, lengths, markers_by_set)
    out, rel = args.outdir, f"r{args.release}"

    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "VERSION.txt"), "w", newline="\n") as fh:
        fh.write(f"v{args.release}.0\nSynthetic sparse GTDB release written by simulate_gtdb_release.py "
                 f"(seed {args.seed}, {len(sim.species)} species, {len(sim.genomes)} genomes)\n")

    for mset in ("bac120", "ar53"):
        genomes = [g for g in sim.genomes if g["species"]["set"] == mset]
        write_taxonomy(os.path.join(out, f"{mset}_taxonomy_{rel}.tsv"), genomes)
        write_metadata(os.path.join(out, f"{mset}_metadata_{rel}.tsv.gz"), genomes)
        write_marker_files(os.path.join(out, "genomic_files_reps", f"{mset}_marker_genes_reps_{rel}"),
                           mset, markers_by_set[mset], [g for g in genomes if g["is_rep"]])
        write_marker_files(os.path.join(out, "genomic_files_all", f"{mset}_marker_genes_all_{rel}"),
                           mset, markers_by_set[mset], genomes)

    sim_dir = os.path.join(out, "simulation")
    os.makedirs(sim_dir, exist_ok=True)
    reps_db = os.path.join(out, "genomic_files_reps", f"gtdb_genomes_reps_{rel}", "database")
    with open(os.path.join(sim_dir, "genomes.tsv"), "w", newline="\n") as gt, \
         open(os.path.join(sim_dir, "marker_positions.tsv"), "w", newline="\n") as mp:
        gt.write("accession\tgtdb_taxonomy\tfasta_path\tgenome_length\tgtdb_representative\n")
        mp.write("accession\tmarker\tcontig\tstart\tend\tstrand\n")
        for g in sim.genomes:
            fasta = (genome_path(reps_db, g["accession"]) if g["is_rep"] else
                     os.path.join(sim_dir, "genomes_nonreps", f"{g['accession']}_genomic.fna.gz"))
            length = write_genome(fasta, g)
            gt.write(f"{g['accession']}\t{g['species']['lineage']}\t{os.path.abspath(fasta)}\t"
                     f"{length}\t{'t' if g['is_rep'] else 'f'}\n")
            for p in g["positions"]:
                mp.write("\t".join([g["accession"], *map(str, p)]) + "\n")

    n_reps = sum(g["is_rep"] for g in sim.genomes)
    sys.stderr.write(f"Wrote {out}: {len(sim.species)} species, {len(sim.genomes)} genomes "
                     f"({n_reps} representatives)\n")


if __name__ == "__main__":
    main()
