"""Shared helpers for the long-read audit: the mini GTDB's genomes and marker positions, the DB's
genes, read simulation with a per-base truth map, and SAM parsing."""
import gzip
import os
import random
import re
import subprocess

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    seqs, name, buf = {}, None, []
    with op(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def read_fasta_list(path):
    """[(name, seq)] keeping duplicates."""
    out, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(buf)))
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        out.append((name, "".join(buf)))
    return out


class World:
    """genomes.tsv, marker_positions.tsv, gene2geneid.tsv, genome2tiid.tsv and the DB's genes."""

    def __init__(self, gtdb, db):
        sim = os.path.join(gtdb, "simulation")
        self.genomes = {}      # accession -> {contig: seq}
        self.taxonomy = {}
        self.rep = {}
        with open(os.path.join(sim, "genomes.tsv")) as fh:
            next(fh)
            for line in fh:
                acc, tax, path, _, rep = line.rstrip("\n").split("\t")
                self.taxonomy[acc] = tax
                self.rep[acc] = rep == "t"
                self.genomes[acc] = read_fasta(path)
        self.geneid = {}
        with open(os.path.join(db, "gene2geneid.tsv")) as fh:
            for line in fh:
                m, g = line.split()[:2]
                self.geneid[m] = int(g)
        self.taxid = {}
        self.rep_of = {}
        with open(os.path.join(db, "genome2tiid.tsv")) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                self.taxid[f[0]] = int(f[1])
                self.rep_of[f[0]] = f[2]
        self.markers = {}  # accession -> [(marker, contig, start0, end0, strand)]
        with open(os.path.join(sim, "marker_positions.tsv")) as fh:
            next(fh)
            for line in fh:
                acc, marker, contig, start, end, strand = line.rstrip("\n").split("\t")
                self.markers.setdefault(acc, []).append((marker, contig, int(start) - 1, int(end), strand))
        # DB genes: reference.fna of an unpacked DB, named taxid_geneid
        self.db_genes = dict(read_fasta_list(os.path.join(db, "reference.fna")))

    def gene_seq(self, acc, marker_entry):
        marker, contig, s, e, strand = marker_entry
        seq = self.genomes[acc][contig][s:e]
        return seq if strand == "+" else revcomp(seq)


def mutate(seq, rng, error, sub_share, ins_share, homopolymer_loss=0.0):
    """Returns (read, origin): origin[i] is the index in seq of read base i, or -1 (inserted).
    Errors: substitutions, 1 bp insertions and 1 bp deletions in the given shares; with
    homopolymer_loss, one base lost in that share of homopolymers of 4 or more."""
    out, origin = [], []
    sub = error * sub_share
    ins = error * (sub_share + ins_share)
    drop_hp = set()
    if homopolymer_loss > 0:
        for m in re.finditer(r"A{4,}|C{4,}|G{4,}|T{4,}", seq):
            if rng.random() < homopolymer_loss:
                drop_hp.add(m.start())
    for i, b in enumerate(seq):
        if i in drop_hp:
            continue
        r = rng.random()
        if r < sub:
            out.append(rng.choice([c for c in "ACGT" if c != b]))
            origin.append(i)
        elif r < ins:
            out.append(b)
            origin.append(i)
            out.append(rng.choice("ACGT"))
            origin.append(-1)
        elif r < error:
            continue  # deletion
        else:
            out.append(b)
            origin.append(i)
    return "".join(out), origin


def parse_sam(path):
    """Yields records as lists of fields (tags kept as a dict under index 11)."""
    if path.endswith(".zst"):
        p = subprocess.Popen(["zstd", "-dc", path], stdout=subprocess.PIPE, text=True)
        fh = p.stdout
    elif path.endswith(".gz"):
        fh = gzip.open(path, "rt")
    else:
        fh = open(path)
    for line in fh:
        if line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        tags = {}
        for t in f[11:]:
            k, ty, v = t.split(":", 2)
            tags[k] = v
        yield f[:11] + [tags]
    fh.close()


def sam_header(path):
    if path.endswith(".zst"):
        out = subprocess.run(["bash", "-c", f"zstd -dc {path} | grep '^@'"], capture_output=True, text=True).stdout
        return out
    with open(path) as fh:
        return "".join(l for l in fh if l.startswith("@"))


CIG = re.compile(r"(\d+)([MIDNSHPX=])")


def cigar_ops(c):
    return [(int(n), op) for n, op in CIG.findall(c)]
