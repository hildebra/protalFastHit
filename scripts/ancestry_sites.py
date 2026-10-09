#!/usr/bin/env python3
"""Which side the reads behind a model's errors take where the species differs from its congeners, and whether their
mismatches fall where the species' own strains vary: the false positives' reads (mostly a novel congener's) against
the false negatives' own reads (mostly a strain's), from error_reads.py's SAMs. build_gtdb_database.py runs it per
read type once the error reads are taken (model_logs/ancestry_sites/, with --share-logs); the next build's models
also carry the same sites as features (protal's AncestrySites.h: ancestry_sites_per_record, ancestry_agreement,
ancestry_congener_share), which the trainer's reports evaluate.

At the sites where a species T differs from its nearest congener in the database (T's derived states, and the
congener's), a strain of T carries T's base; a species that branched off T's lineage at fraction x of T's branch
carries T's base at about x of them and the congener's base at the rest. With the species' other genomes
(--full-reference, the build's full_reference.fna.zst: every genome's marker genes under the species' taxid), the
sites split further: those where T's own strains vary (polymorphic: a mismatch there is within-species variation,
not evidence against T) and those fixed in T (a difference from the congeners that every known strain shares: a
synapomorphy, strong evidence when the read carries it, strong evidence against when the read carries the
congener's base). Nothing of this is known to the reads' alignments, which put both at the same identity
(docs/claude/2026-10-07-error-read-signatures).

For every FP and FN taxon of the error SAMs: its copies of the genes its counted reads hit, its congeners' copies
(the species of its genus in the training database, held-out species left out, at most --max-congeners), the
differing sites against the nearest congener by identity on that gene (and against the majority of the three
nearest), the polymorphic sites from up to --max-alleles other genomes of the species, and for each counted read on
the taxon the covered sites of each class and the read's base there. Writes OUT.fragments.tsv.gz (one line per
record), OUT.taxa.tsv.gz (per taxon, pooled), OUT.auc.tsv (the AUC of each signal for the FN taxa's own reads
against the FP reads, per taxon, overall and within identity bands) and OUT.summary.txt (the same in words, also on
stdout).

The gene copies are compared along their shared 12-mers (insilico_strains.substitutions), without an alignment, so a
stretch past an indel is not compared.

One run serves several read types: --sams and --out then take a value each per SAM folder, paired in order. The run
reads every folder's records first, makes one pass over --reference and one over --full-reference for the copies and
alleles any folder wants, and then writes each folder's four outputs at its prefix, byte for byte those of a run on
that folder alone (each folder's congeners drawn by a generator seeded afresh, its copies and alleles filtered back to
those it wants); with --threads the folders' analyses run side by side in forked processes. Each folder's lines on
stdout follow a line naming its prefix, in the order given.

    python3 ancestry_sites.py --sams model_logs/error_reads/pe --reference training_db/reference.fna \\
        --full-reference training_db/full_reference.fna.zst --taxonomy work/internal_taxonomy.dmp \\
        --heldout model_logs/heldout_species.txt --out model_logs/ancestry_sites/pe
    (a bundled database: protal --db training_db/database.protal --unpack_db --unpack_dir DIR gives its reference.fna;
    the full reference is one pass of 86 GB at r226, about 10 minutes, and a few GB of memory for pe)
    python3 ancestry_sites.py --sams model_logs/error_reads/{pe,se,pb,ont} --out model_logs/ancestry_sites/{pe,se,pb,ont} \\
        --threads 4 ...  (the four read types over one pass of each reference)
"""
import argparse
import collections
import csv
import glob
import gzip
import io
import multiprocessing
import os
import random
import re
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from insilico_strains import CODE, kmer_codes, substitutions  # noqa: E402
sys.path.insert(0, os.path.join(HERE, "mini_db"))
from gtdb_to_protal_db import allele_genome  # noqa: E402

MIN_MAPQ = 4
CIGAR = re.compile(rb"(\d+)([MIDNSHP=X])")
NEAREST = 3  # congeners whose majority base defines the second set of sites
CONSENSUS = 0.9  # protal's consensus sites (AncestrySites.h kConsensus): of the congeners compared at a position, the
MIN_CONGENERS = 3  # share that carry one other base, where at least this many were compared (kMinCongeners); else the
#                    nearest congener's difference
BANDS = ((0.93, 0.96), (0.96, 0.975), (0.975, 0.99))
MIN_SITES = (3, 10, 30)
COUNTS = ["aligned", "mismatches", "sites1", "agree1", "alt1", "sites3", "agree3", "alt3", "sites4", "agree4", "alt4",
          "fixed1", "fixed_agree1", "fixed_alt1", "poly_covered", "poly_mismatches", "nonpoly_aligned",
          "nonpoly_mismatches"]


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--sams", required=True, nargs="+",
                   help="folder of <set>/<point>/<sample>[.FP|.FN].sam[.zst] (error_reads.py's); several (a read type's "
                        "each) share one pass over each reference, each reported at the --out prefix in its place")
    p.add_argument("--reference", required=True, help="reference.fna (or .zst) of the training database: >taxid_geneid")
    p.add_argument("--full-reference", help="full_reference.fna[.zst]: every genome's marker genes, >taxid_geneid of its "
                                            "species; gives the species' polymorphic sites")
    p.add_argument("--taxonomy", required=True, help="internal_taxonomy.dmp")
    p.add_argument("--heldout", required=True, help="heldout_species.txt (species the training database lacks)")
    p.add_argument("--out", required=True, nargs="+", help="output prefix, one per --sams folder in the same order")
    p.add_argument("--max-congeners", type=int, default=8)
    p.add_argument("--max-alleles", type=int, default=6, help="other genomes' copies kept per species and gene")
    p.add_argument("--allele-genome-share", type=float, default=1.0,
                   help="take alleles only from the genomes that give protal's strain alleles (gtdb_to_protal_db.allele_genome, "
                        "protal --allele_genome_share): build_gtdb_database.py simulates strains from the others, so the "
                        "report's polymorphic sites are those of the alleles protal knows and never a simulated strain's own "
                        "(default 1: every genome)")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("-t", "--threads", type=int, default=1,
                   help="with several --sams folders, analyse up to this many side by side (forked processes, after the "
                        "shared passes over the references; the outputs are the same)")
    opts = p.parse_args(argv)
    if len(opts.sams) != len(opts.out):
        p.error(f"--sams and --out take a value each per SAM folder, paired in order: {len(opts.sams)} SAM folder"
                f"{'s' * (len(opts.sams) != 1)}, {len(opts.out)} output prefix{'es' * (len(opts.out) != 1)}")
    prefixes = [os.path.abspath(out) for out in opts.out]
    if len(set(prefixes)) < len(prefixes):
        p.error("--out names a prefix twice: " + ", ".join(sorted({out for out, a in zip(opts.out, prefixes)
                                                                   if prefixes.count(a) > 1})))
    return opts


class Taxonomy:
    def __init__(self, path, heldout):
        self.parent, self.rank, self.name = {}, {}, {}
        with open(path, newline="") as fh:
            reader = csv.reader(fh, delimiter="\t")
            next(reader)
            for row in reader:
                self.parent[row[0]], self.rank[row[0]], self.name[row[0]] = row[1], row[4], row[3]
        self.name_id = {}
        for t, r in self.rank.items():
            if r == "species":
                self.name_id.setdefault(self.name[t], t)
        held = set()
        if heldout and os.path.isfile(heldout):
            with open(heldout) as fh:
                held = {line.split("\t")[0].strip() for line in fh if line.strip()}
        self.held_out = {self.name_id[n] for n in held if n in self.name_id}
        self.genus_members = collections.defaultdict(list)
        for t, r in self.rank.items():
            if r == "species" and t not in self.held_out:
                g = self.ancestor(t, "genus")
                if g:
                    self.genus_members[g].append(t)

    def ancestor(self, taxid, rank):
        for _ in range(64):
            if taxid not in self.rank:
                return None
            if self.rank[taxid] == rank:
                return taxid
            if self.parent[taxid] == taxid:
                return None
            taxid = self.parent[taxid]
        return None

    def relation(self, sid, taxid):
        if not sid:
            return "unknown"
        if sid == taxid:
            return "own"
        for rank in ("genus", "family", "order"):
            a = self.ancestor(sid, rank)
            if a is not None and a == self.ancestor(taxid, rank):
                return rank
        return "other"

    def congeners(self, taxid, cap, rng):
        g = self.ancestor(taxid, "genus")
        members = [t for t in self.genus_members.get(g, []) if t != taxid] if g else []
        if len(members) > cap:
            members = rng.sample(members, cap)
        return members


def open_text(path):
    """Lines (bytes) of a plain or zstd-compressed file."""
    if path.endswith(".zst"):
        proc = subprocess.Popen(["zstd", "-dcq", path], stdout=subprocess.PIPE)
        yield from proc.stdout
        proc.wait()
        if proc.returncode:
            raise RuntimeError(f"zstd -dc {path} exited {proc.returncode}")
    else:
        with open(path, "rb") as fh:
            yield from fh


def gzip_text(path):
    """A gzip file to write text to, without a time stamp in its header: the same lines give the same bytes, in a run
    on one SAM folder as in a run on several."""
    return io.TextIOWrapper(gzip.GzipFile(path, "wb", mtime=0), encoding="utf-8")


def echo(line):
    print(line, flush=True)


def records_of(path, tax, seen):
    """The counted records of the error taxa in one SAM: (qname, taxid, gene, pos, cigar, seq, role, source taxid,
    source genome) for primary records at MAPQ >= 4 whose taxon is an FP taxon of the read (xe FP:<taxid>) or an FN
    taxon the read belongs to (xe FN:<taxid> and xs the taxon's species); a record in both of a sample's files (FP and
    FN) once (seen). The source genome is the xg tag's accession (accession_of), "" without one."""
    out = []
    for line in open_text(path):
        if line.startswith(b"@"):
            continue
        f = line.rstrip(b"\n").split(b"\t")
        if len(f) < 11:
            continue
        flag, rname, mapq = int(f[1]), f[2], int(f[4])
        if flag & 0x904 or mapq < MIN_MAPQ or rname == b"*":
            continue
        taxid, _, gene = rname.partition(b"_")
        xs = xe = xg = ""
        for tag in f[11:]:
            if tag.startswith(b"xs:Z:"):
                xs = tag[5:].decode()
            elif tag.startswith(b"xe:Z:"):
                xe = tag[5:].decode()
            elif tag.startswith(b"xg:Z:"):
                xg = tag[5:].decode()
        t = taxid.decode()
        reasons = set(xe.split(","))
        source = tax.name_id.get(xs.replace(" (not in the database)", ""), "")
        if f"FP:{t}" in reasons:
            role = "FP"
        elif f"FN:{t}" in reasons and source == t:
            role = "FN own"
        else:
            continue
        key = (f[0], rname, f[3], f[5])
        if key in seen:
            continue
        seen.add(key)
        out.append((f[0].decode(), t, int(gene), int(f[3]), f[5], f[9], role, source,
                    accession_of(xg) if xg and xg != "?" else ""))
    return out


def aligned_pairs(pos, cigar, seq):
    """(gene position 0-based, read base code) for the aligned bases of a record."""
    ref = pos - 1
    read = 0
    pairs = []
    for n, op in CIGAR.findall(cigar):
        n = int(n)
        if op in b"M=X":
            codes = CODE[np.frombuffer(seq[read:read + n], dtype=np.uint8)]
            pairs.append((np.arange(ref, ref + n), codes))
            ref += n
            read += n
        elif op in b"IS":
            read += n
        elif op in b"DN":
            ref += n
    if not pairs:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.uint8)
    return np.concatenate([p[0] for p in pairs]), np.concatenate([p[1] for p in pairs])


ACCESSION_RE = re.compile(r"(GC[AF]_\d{9}\.\d+)")


def accession_of(text):
    """The GTDB accession in a genome's name (GCA_ or GCF_, nine digits, a version), else the name's first token: how
    the full reference's records (the converter, after the name) and the error reads' xg tag name genomes."""
    m = ACCESSION_RE.search(text)
    if m:
        return m.group(1)
    return text.split()[0] if text.split() else text


def fasta_records(path, wanted):
    """(header first token, the header's second token or "", sequence bytes) of the FASTA records whose first token is
    in wanted."""
    name, genome, chunks = None, "", []
    for line in open_text(path):
        if line.startswith(b">"):
            if name is not None:
                yield name, genome, b"".join(chunks)
            tokens = line[1:].split()
            token = tokens[0].decode() if tokens else ""
            name, chunks = (token, []) if token in wanted else (None, [])
            genome = tokens[1].decode() if len(tokens) > 1 else ""
        elif name is not None:
            chunks.append(line.strip())
    if name is not None:
        yield name, genome, b"".join(chunks)


def load_copies(path, wanted):
    """{header first token: codes} of the reference's records whose first token is in wanted."""
    return {name: CODE[np.frombuffer(seq, dtype=np.uint8)] for name, _, seq in fasta_records(path, wanted)}


def load_alleles(path, wanted, reps, cap, share=1.0):
    """{header: [(genome accession or "", codes)]}: up to cap distinct copies per wanted header that differ from the
    representative's, each with its genome (the record's second token, as the converter writes it since 2026-10-08; ""
    from an older full reference, which then leaves no read's own genome out). With share below 1 only the copies of
    genomes that give strain alleles (allele_genome), as protal --build --allele_genome_share takes them; a copy
    without its genome named is kept (an older full reference cannot be split)."""
    alleles = collections.defaultdict(list)
    seen = collections.defaultdict(set)
    for name, genome, seq in fasta_records(path, wanted):
        rep = reps.get(name)
        if rep is None or len(alleles[name]) >= cap or seq in seen[name]:
            continue
        if share < 1 and genome and not allele_genome(accession_of(genome), share):
            continue
        codes = CODE[np.frombuffer(seq, dtype=np.uint8)]
        if len(codes) == len(rep) and (codes == rep).all():
            continue
        seen[name].add(seq)
        alleles[name].append((accession_of(genome) if genome else "", codes))
    return alleles


def diagonal(rep, other, k=12):
    """The main diagonal substitutions() used (positions on rep minus positions on other)."""
    a, b = kmer_codes(rep, k), kmer_codes(other, k)
    ua, ia, ca = np.unique(a, return_index=True, return_counts=True)
    ub, ib, cb = np.unique(b, return_index=True, return_counts=True)
    keep_a, keep_b = (ca == 1) & (ua >= 0), (cb == 1) & (ub >= 0)
    _, pa, pb = np.intersect1d(ua[keep_a], ub[keep_b], assume_unique=True, return_indices=True)
    if not len(pa):
        return 0
    diagonals, counts = np.unique(ia[keep_a][pa] - ib[keep_b][pb], return_counts=True)
    return int(diagonals[np.argmax(counts)])


class Sites:
    """The site classes of one gene copy of T: against the nearest congener (alt1: the congener's base at the sites
    where it differs from T, -1 elsewhere), against the majority of the NEAREST nearest (alt3), protal's consensus over
    every congener compared (alt4: where MIN_CONGENERS or more were compared at the position, the base CONSENSUS of
    them carry if it is not T's; with fewer, the nearest's difference), and the sites where T's own alleles differ
    from the representative (poly; per allele in poly_by, so that a read's own source genome can be left out:
    poly_without). A missed real strain was simulated from a GTDB genome, usually one of the alleles, whose every
    difference from the representative would otherwise be polymorphic by construction (docs/claude/2026-10-08-r226-v18,
    section 3)."""

    def __init__(self, t_copy, congener_copies, alleles):
        n = len(t_copy)
        self.alt1 = np.full(n, -1, dtype=np.int8)
        self.alt3 = np.full(n, -1, dtype=np.int8)
        self.alt4 = np.full(n, -1, dtype=np.int8)
        self.poly = np.zeros(n, dtype=bool)
        self.nearest_identity = None
        self.allele_divergence = 0.0
        ranked = []
        for copy in congener_copies:
            positions, compared, mask = substitutions(t_copy, copy, with_mask=True)
            if compared < n // 2:
                continue
            d = diagonal(t_copy, copy)
            inside = (positions - d >= 0) & (positions - d < len(copy))
            bases = copy[positions[inside] - d]
            inside[inside] = bases < 4  # a site whose congener base is another letter (N) is no site
            ranked.append((1 - len(positions) / compared, positions[inside], copy[positions[inside] - d], mask))
        if ranked:
            ranked.sort(key=lambda r: -r[0])
            self.nearest_identity = ranked[0][0]
            self.alt1[ranked[0][1]] = ranked[0][2]
            k = min(NEAREST, len(ranked))
            votes = np.zeros((n, 4), dtype=np.int16)
            for _, positions, bases, _ in ranked[:k]:
                votes[positions, bases] += 1
            majority = votes.argmax(axis=1)
            count = votes.max(axis=1)
            self.alt3[count * 2 > k] = majority[count * 2 > k]
            votes = np.zeros((n, 4), dtype=np.int16)
            compared_by = np.zeros(n, dtype=np.int16)
            for _, positions, bases, mask in ranked:
                votes[positions, bases] += 1
                compared_by += mask
            majority = votes.argmax(axis=1)
            count = votes.max(axis=1)
            enough = compared_by >= MIN_CONGENERS
            self.alt4[~enough] = self.alt1[~enough]
            consensus = enough & (count > 0) & (count >= CONSENSUS * compared_by)
            self.alt4[consensus] = majority[consensus]
        self.poly_by = []  # (genome accession or "", its allele's polymorphic sites)
        for genome, copy in alleles:
            positions, compared = substitutions(t_copy, copy)
            if compared < n // 2:
                continue
            mask = np.zeros(n, dtype=bool)
            mask[positions] = True
            self.poly_by.append((genome, mask))
            self.poly |= mask
            self.allele_divergence = max(self.allele_divergence, len(positions) / compared)
        self.fixed1 = (self.alt1 >= 0) & ~self.poly  # differs from the nearest congener, fixed in T's alleles

    def has_allele_of(self, genome):
        """Whether one of the alleles is `genome`'s copy."""
        return bool(genome) and any(g == genome for g, _ in self.poly_by)

    def poly_without(self, genome):
        """The polymorphic sites by the alleles of the other genomes than `genome` (a read's own source genome left
        out); all of them when `genome` is not among them."""
        if not self.has_allele_of(genome):
            return self.poly
        poly = np.zeros(len(self.poly), dtype=bool)
        for g, mask in self.poly_by:
            if g != genome:
                poly |= mask
        return poly


def auc(y, x):
    """AUC of x for y == 1 against y == 0 (rank-based, ties averaged), or nan."""
    y, x = np.asarray(y, dtype=bool), np.asarray(x, dtype=float)
    n1, n0 = int(y.sum()), int((~y).sum())
    if n1 == 0 or n0 == 0:
        return float("nan")
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x))
    ranks[order] = np.arange(1, len(x) + 1)
    xs = x[order]
    i = 0
    while i < len(xs):
        j = i
        while j + 1 < len(xs) and xs[j + 1] == xs[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = (i + 1 + j + 1) / 2
        i = j + 1
    return float((ranks[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def count_record(sites, own, ref_pos, bases, genome=""):
    """The site counts of one record on T's copy (COUNTS); the polymorphic and fixed sites by the alleles of other
    genomes than the read's own source `genome`."""
    mism = bases != own[ref_pos]
    out = [len(ref_pos), int(mism.sum())]
    for alt in (sites.alt1, sites.alt3, sites.alt4):
        covered = alt[ref_pos] >= 0
        out += [int(covered.sum()), int((~mism[covered]).sum()), int((bases[covered] == alt[ref_pos][covered]).sum())]
    polymorphic = sites.poly_without(genome)
    fixed1 = (sites.alt1 >= 0) & ~polymorphic
    covered = fixed1[ref_pos]
    out += [int(covered.sum()), int((~mism[covered]).sum()), int((bases[covered] == sites.alt1[ref_pos][covered]).sum())]
    poly = polymorphic[ref_pos]
    out += [int(poly.sum()), int(mism[poly].sum()), int((~poly).sum()), int(mism[~poly].sum())]
    return out


def signals(counts, with_alleles):
    """Per taxon, the signals of its pooled counts: {name: value}."""
    c = dict(zip(COUNTS, counts))
    out = {"identity": 1 - c["mismatches"] / max(1, c["aligned"]),
           "species_base_at_congener_sites": c["agree1"] / max(1, c["sites1"]),
           "species_base_at_majority_sites": c["agree3"] / max(1, c["sites3"]),
           "species_base_at_consensus_sites": c["agree4"] / max(1, c["sites4"])}
    if with_alleles:
        out.update({"species_base_at_fixed_sites": c["fixed_agree1"] / max(1, c["fixed1"]),
                    "fixed_site_identity": 1 - c["nonpoly_mismatches"] / max(1, c["nonpoly_aligned"]),
                    "mismatches_at_polymorphic_sites": c["poly_mismatches"] / max(1, c["mismatches"])})
    return out


def summarize(taxa, with_alleles):
    """-> (the summary's lines, the AUC rows (min_sites, band, signal, taxa, FN, AUC))."""
    lines = ["Per group, pooled over records (identity: median over taxa; shares: pooled):"]
    groups = collections.defaultdict(lambda: np.zeros(len(COUNTS) + 1))
    ids = collections.defaultdict(list)
    for sample, t, group, n, counts in taxa:
        groups[group] += np.concatenate([[n], counts])
        if counts[0]:
            ids[group].append(1 - counts[1] / counts[0])
    for group, v in sorted(groups.items()):
        c = dict(zip(["records"] + COUNTS, v))
        text = (f"  {group:14s} records {int(c['records']):7d}, identity {np.median(ids[group]):.4f}, sites/record "
                f"{c['sites1'] / max(1, c['records']):.2f}, species' base {c['agree1'] / max(1, c['sites1']):.3f}, "
                f"congener's {c['alt1'] / max(1, c['sites1']):.3f} (nearest); {c['agree3'] / max(1, c['sites3']):.3f} / "
                f"{c['alt3'] / max(1, c['sites3']):.3f} (majority of 3); {c['agree4'] / max(1, c['sites4']):.3f} / "
                f"{c['alt4'] / max(1, c['sites4']):.3f} at {c['sites4'] / max(1, c['records']):.2f}/record (consensus, "
                f"protal's)")
        if with_alleles:
            text += (f"; fixed sites/record {c['fixed1'] / max(1, c['records']):.2f}, species' base "
                     f"{c['fixed_agree1'] / max(1, c['fixed1']):.3f}, congener's {c['fixed_alt1'] / max(1, c['fixed1']):.3f}; "
                     f"mismatches at polymorphic sites {c['poly_mismatches'] / max(1, c['mismatches']):.3f} "
                     f"(polymorphic sites cover {c['poly_covered'] / max(1, c['aligned']):.4f} of the aligned bases); "
                     f"fixed-site identity {1 - c['nonpoly_mismatches'] / max(1, c['nonpoly_aligned']):.4f}")
        lines.append(text)
    lines.append("")
    lines.append("Per taxon, FN own (1) against FP genus (0), AUC of each signal (identity: FN low; sites: FN high), taxa "
                 "with at least N sites against the nearest congener:")
    rows = [(g == "FN own", signals(counts, with_alleles), counts) for _, _, g, _, counts in taxa
            if g in ("FN own", "FP genus")]
    aucs = []
    for n_min in MIN_SITES:
        sel = [r for r in rows if r[2][COUNTS.index("sites1")] >= n_min]
        y = np.array([r[0] for r in sel])
        if not len(sel) or not 0 < y.sum() < len(sel):
            continue
        lines.append(f"  N >= {n_min:2d}: {len(sel)} taxa ({int(y.sum())} FN)")
        for name in sel[0][1]:
            a = auc(y, [r[1][name] for r in sel])
            aucs.append((n_min, "all", name, len(sel), int(y.sum()), a))
            line = f"    {name:34s} {a:.3f}"
            for lo, hi in BANDS:
                band = [r for r in sel if lo <= r[1]['identity'] < hi]
                yb = np.array([r[0] for r in band])
                if len(band) and 0 < yb.sum() < len(band):
                    ab = auc(yb, [r[1][name] for r in band])
                    aucs.append((n_min, f"{lo}-{hi}", name, len(band), int(yb.sum()), ab))
                    line += f"   identity {lo}-{hi}: {ab:.3f} ({len(band)})"
            lines.append(line)
    if len(aucs) == 0:
        lines.append("  (no taxa of both kinds with sites: nothing to compare)")
    return lines, aucs


def sams_of(folder):
    """The SAMs of one --sams folder: <set>/<point>/<sample>[.FP|.FN].sam[.zst], or at its top."""
    return sorted(glob.glob(os.path.join(folder, "*", "*", "*.sam*")) + glob.glob(os.path.join(folder, "*.sam*")))


class Folder:
    """One SAM folder's report before its analysis: every sample's counted records (records_of), the taxa they hit and
    the congeners drawn for them, the reference copies it wants (own_keys: its taxa's copies of the genes they hit;
    wanted: those and the congeners' copies of the same genes) and its output prefix. The congeners are drawn by a
    generator seeded afresh (--seed) for each folder, so that they are those of a run on the folder alone. head: the
    summary's first lines, printed as they are known (echo) or, in a run on several folders, kept for the folder's
    block (said)."""

    def __init__(self, sams_dir, out, sams, tax, opts, echoed):
        self.sams_dir, self.out, self.echoed = sams_dir, out, echoed
        self.head, self.said = [], []
        self.records, seen = {}, {}
        for path in sams:
            parts = path.replace("\\", "/").split("/")
            sample = (parts[-3] + ":" if len(parts) >= 3 else "") + os.path.basename(path).split(".")[0]
            self.records.setdefault(sample, []).extend(records_of(path, tax, seen.setdefault(sample, set())))
        n_records = sum(len(r) for r in self.records.values())
        self.tell(f"{len(sams)} SAMs of {len(self.records)} samples, {n_records} counted records on FP taxa or of FN taxa's "
                  f"own reads")
        rng = random.Random(opts.seed)
        self.hits = collections.defaultdict(set)
        for recs in self.records.values():
            for _, t, g, *_ in recs:
                self.hits[t].add(g)
        self.congeners = {t: tax.congeners(t, opts.max_congeners, rng) for t in self.hits}
        self.own_keys = {f"{t}_{g}" for t, genes in self.hits.items() for g in genes}
        self.wanted = set(self.own_keys)
        for t, genes in self.hits.items():
            for g in genes:
                for c in self.congeners[t]:
                    self.wanted.add(f"{c}_{g}")
        self.tell(f"{len(self.hits)} taxa, {len(self.own_keys)} taxon-genes, {len(self.wanted)} copies wanted; taxa without "
                  f"a congener in the training database: {sum(1 for t in self.hits if not self.congeners[t])}")
        self.copies, self.alleles = {}, {}

    def tell(self, line):
        self.head.append(line)
        if self.echoed:
            echo(line)
        else:
            self.said.append(line)

    def take(self, copies, alleles, opts):
        """Keeps of the shared passes' copies and alleles those this folder wants, in the passes' order: a header's copy
        and its alleles (load_alleles caps them per header) depend on that header's records alone, so these are what a
        run on the folder alone reads."""
        self.copies = {k: v for k, v in copies.items() if k in self.wanted}
        self.tell(f"{len(self.copies)} copies read from {opts.reference}")
        if not opts.full_reference:
            return
        self.alleles = {k: v for k, v in alleles.items() if k in self.own_keys}
        named = sum(1 for v in self.alleles.values() for genome, _ in v if genome)
        self.tell(f"alleles: {sum(len(v) for v in self.alleles.values())} copies of {len(self.alleles)} taxon-genes from "
                  f"{opts.full_reference}; taxa with any: {len({k.split('_')[0] for k, v in self.alleles.items() if v})} of "
                  f"{len(self.hits)}; {named} with their genome named (a read's own source genome is left out of its "
                  f"alleles)" +
                  (f"; only from the genomes that give strain alleles (--allele-genome-share {opts.allele_genome_share:g})"
                   if opts.allele_genome_share < 1 else "") +
                  ("" if named else "; NONE named: an older full reference, so no read's own genome is left out and the "
                                    "fixed sites of real strains are circular"))


def report(folder, tax, opts, say):
    """The sites of one folder's taxa and its records' counts at them: writes OUT.fragments.tsv.gz, OUT.taxa.tsv.gz,
    OUT.auc.tsv and OUT.summary.txt at the folder's prefix, and says (say, a line at a time) the line on the records
    whose own genome's allele was left out and the summary."""
    copies, alleles, congeners, records = folder.copies, folder.alleles, folder.congeners, folder.records
    sites = {}
    for t, genes in folder.hits.items():
        for g in genes:
            key = f"{t}_{g}"
            if key in copies:
                sites[(t, g)] = Sites(copies[key], [copies[f"{c}_{g}"] for c in congeners[t] if f"{c}_{g}" in copies],
                                      alleles.get(key, []))

    os.makedirs(os.path.dirname(os.path.abspath(folder.out)) or ".", exist_ok=True)
    per = collections.defaultdict(lambda: [0, np.zeros(len(COUNTS), dtype=np.int64)])
    left_out = collections.Counter()  # records whose own source genome's allele was left out, by role
    with gzip_text(folder.out + ".fragments.tsv.gz") as fh:
        fh.write("sample\tqname\ttaxid\trole\trelation\tgene\tnearest_identity\tallele_divergence\t" + "\t".join(COUNTS) + "\n")
        for sample, recs in records.items():
            for qname, t, g, pos, cigar, seq, role, source, genome in recs:
                key = (t, g)
                own = copies.get(f"{t}_{g}")
                if own is None or key not in sites:
                    continue
                s = sites[key]
                ref_pos, bases = aligned_pairs(pos, cigar, seq)
                inside = ref_pos < len(own)
                left_out[role] += s.has_allele_of(genome)
                counts = count_record(s, own, ref_pos[inside], bases[inside], genome)
                rel = tax.relation(source, t)
                fh.write("\t".join(map(str, [sample, qname, t, role, rel, g,
                                             "" if s.nearest_identity is None else f"{s.nearest_identity:.4f}",
                                             f"{s.allele_divergence:.4f}"] + counts)) + "\n")
                group = role if role == "FN own" else f"FP {rel}"
                v = per[(sample, t, group)]
                v[0] += 1
                v[1] += np.array(counts)

    taxa = [(sample, t, group, n, counts) for (sample, t, group), (n, counts) in per.items()]
    with gzip_text(folder.out + ".taxa.tsv.gz") as fh:
        fh.write("sample\ttaxid\tgroup\trecords\tidentity\t" + "\t".join(COUNTS) + "\n")
        for sample, t, group, n, counts in taxa:
            idn = 1 - counts[1] / counts[0] if counts[0] else float("nan")
            fh.write("\t".join(map(str, [sample, t, group, n, f"{idn:.4f}"] + counts.tolist())) + "\n")
    head = list(folder.head)
    if opts.full_reference:
        head.append("records whose own source genome's allele was left out: " +
                    (", ".join(f"{role} {n}" for role, n in sorted(left_out.items())) or "none"))
        say(head[-1])
    lines, aucs = summarize(taxa, bool(opts.full_reference))
    with open(folder.out + ".auc.tsv", "w") as fh:
        fh.write("min_sites\tidentity_band\tsignal\ttaxa\tfn\tauc\n")
        for n_min, band, name, n, fn, a in aucs:
            fh.write(f"{n_min}\t{band}\t{name}\t{n}\t{fn}\t{a:.4f}\n")
    with open(folder.out + ".summary.txt", "w") as fh:
        fh.write("\n".join(head + [""] + lines) + "\n")
    say("\n" + "\n".join(lines))


_SHARED = None  # (folders, taxonomy, options) of the analyses in forked processes (report_shared)


def report_shared(i):
    """report() on the i-th folder, in a forked process: the lines it says."""
    folders, tax, opts = _SHARED
    said = []
    report(folders[i], tax, opts, said.append)
    return said


def reports(folders, tax, opts):
    """The lines report() says of each folder, in the folders' order: up to --threads folders side by side in forked
    processes, which share the parent's copies and alleles (where the system forks), else one after the other."""
    global _SHARED
    workers = min(opts.threads, len(folders))
    if workers > 1 and "fork" in multiprocessing.get_all_start_methods():
        sys.stdout.flush()  # nothing buffered for the forked processes to print again
        sys.stderr.flush()
        _SHARED = (folders, tax, opts)
        try:
            with multiprocessing.get_context("fork").Pool(workers) as pool:
                yield from pool.imap(report_shared, range(len(folders)))
        finally:
            _SHARED = None
        return
    for folder in folders:
        said = []
        report(folder, tax, opts, said.append)
        yield said


def main(argv=None):
    opts = parse_args(argv)
    for path, what in ((opts.reference, "reference"), (opts.taxonomy, "taxonomy")):
        if not os.path.isfile(path):
            sys.exit(f"ancestry_sites.py: no {what} {path}")
    if opts.full_reference and not os.path.isfile(opts.full_reference):
        print(f"no full reference {opts.full_reference}: the congener sites only", flush=True)
        opts.full_reference = None
    tax = Taxonomy(opts.taxonomy, opts.heldout)
    found = []
    for sams_dir in opts.sams:
        found.append(sams_of(sams_dir))
        if not found[-1]:
            sys.exit(f"ancestry_sites.py: no SAM under {sams_dir}")
    several = len(opts.sams) > 1
    folders = []
    for sams_dir, out, sams in zip(opts.sams, opts.out, found):
        folders.append(Folder(sams_dir, out, sams, tax, opts, not several))
        if several:
            echo(f"{sams_dir}: {folders[-1].head[0]}")

    # One pass over each reference for the copies and alleles any folder wants; each folder keeps its own of them.
    wanted = set().union(*(folder.wanted for folder in folders))
    own_keys = set().union(*(folder.own_keys for folder in folders))
    if several:
        echo(f"{len(folders)} SAM folders: {len(own_keys)} taxon-genes, {len(wanted)} copies wanted in all, read in one "
             f"pass over each reference")
    began = time.time()
    copies = load_copies(opts.reference, wanted)
    if several:
        echo(f"{len(copies)} copies read from {opts.reference} in {time.time() - began:.0f} s")
    alleles = {}
    if opts.full_reference:
        began = time.time()
        alleles = load_alleles(opts.full_reference, own_keys, {k: copies[k] for k in own_keys if k in copies},
                               opts.max_alleles, opts.allele_genome_share)
        if several:
            echo(f"alleles: {sum(len(v) for v in alleles.values())} copies of {len(alleles)} taxon-genes read from "
                 f"{opts.full_reference} in {time.time() - began:.0f} s")
    for folder in folders:
        folder.take(copies, alleles, opts)
    del copies, alleles

    if not several:
        report(folders[0], tax, opts, echo)
        return 0
    for folder, said in zip(folders, reports(folders, tax, opts)):
        echo("")
        echo(f"== {folder.out} (the SAMs of {folder.sams_dir}) ==")
        for line in folder.said + said:
            echo(line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
