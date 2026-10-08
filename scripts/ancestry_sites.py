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

    python3 ancestry_sites.py --sams model_logs/error_reads/pe --reference training_db/reference.fna \\
        --full-reference training_db/full_reference.fna.zst --taxonomy internal_taxonomy.dmp \\
        --heldout heldout_species.txt --out model_logs/ancestry_sites/pe
    (a bundled database: protal --db training_db/database.protal --unpack_db --unpack_dir DIR gives its reference.fna;
    the full reference is one pass of 86 GB at r226, about 10 minutes, and a few GB of memory for pe)
"""
import argparse
import collections
import csv
import glob
import gzip
import os
import random
import re
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from insilico_strains import CODE, kmer_codes, substitutions  # noqa: E402

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
    p.add_argument("--sams", required=True, help="folder of <set>/<point>/<sample>[.FP|.FN].sam[.zst] (error_reads.py's)")
    p.add_argument("--reference", required=True, help="reference.fna (or .zst) of the training database: >taxid_geneid")
    p.add_argument("--full-reference", help="full_reference.fna[.zst]: every genome's marker genes, >taxid_geneid of its "
                                            "species; gives the species' polymorphic sites")
    p.add_argument("--taxonomy", required=True, help="internal_taxonomy.dmp")
    p.add_argument("--heldout", required=True, help="heldout_species.txt (species the training database lacks)")
    p.add_argument("--out", required=True, help="output prefix")
    p.add_argument("--max-congeners", type=int, default=8)
    p.add_argument("--max-alleles", type=int, default=6, help="other genomes' copies kept per species and gene")
    p.add_argument("--seed", type=int, default=1)
    return p.parse_args(argv)


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


def records_of(path, tax, seen):
    """The counted records of the error taxa in one SAM: (qname, taxid, gene, pos, cigar, seq, role, source taxid)
    for primary records at MAPQ >= 4 whose taxon is an FP taxon of the read (xe FP:<taxid>) or an FN taxon the read
    belongs to (xe FN:<taxid> and xs the taxon's species); a record in both of a sample's files (FP and FN) once
    (seen)."""
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
        xs = xe = ""
        for tag in f[11:]:
            if tag.startswith(b"xs:Z:"):
                xs = tag[5:].decode()
            elif tag.startswith(b"xe:Z:"):
                xe = tag[5:].decode()
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
        out.append((f[0].decode(), t, int(gene), int(f[3]), f[5], f[9], role, source))
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


def fasta_records(path, wanted):
    """(header first token, sequence bytes) of the FASTA records whose first token is in wanted."""
    name, chunks = None, []
    for line in open_text(path):
        if line.startswith(b">"):
            if name is not None:
                yield name, b"".join(chunks)
            token = line[1:].split()[0].decode() if line[1:].split() else ""
            name, chunks = (token, []) if token in wanted else (None, [])
        elif name is not None:
            chunks.append(line.strip())
    if name is not None:
        yield name, b"".join(chunks)


def load_copies(path, wanted):
    """{header first token: codes} of the reference's records whose first token is in wanted."""
    return {name: CODE[np.frombuffer(seq, dtype=np.uint8)] for name, seq in fasta_records(path, wanted)}


def load_alleles(path, wanted, reps, cap):
    """{header: [codes]}: up to cap distinct copies per wanted header that differ from the representative's."""
    alleles = collections.defaultdict(list)
    seen = collections.defaultdict(set)
    for name, seq in fasta_records(path, wanted):
        rep = reps.get(name)
        if rep is None or len(alleles[name]) >= cap or seq in seen[name]:
            continue
        codes = CODE[np.frombuffer(seq, dtype=np.uint8)]
        if len(codes) == len(rep) and (codes == rep).all():
            continue
        seen[name].add(seq)
        alleles[name].append(codes)
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
    from the representative (poly)."""

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
        for copy in alleles:
            positions, compared = substitutions(t_copy, copy)
            if compared < n // 2:
                continue
            self.poly[positions] = True
            self.allele_divergence = max(self.allele_divergence, len(positions) / compared)
        self.fixed1 = (self.alt1 >= 0) & ~self.poly  # differs from the nearest congener, fixed in T's alleles


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


def count_record(sites, own, ref_pos, bases):
    """The site counts of one record on T's copy (COUNTS)."""
    mism = bases != own[ref_pos]
    out = [len(ref_pos), int(mism.sum())]
    for alt in (sites.alt1, sites.alt3, sites.alt4):
        covered = alt[ref_pos] >= 0
        out += [int(covered.sum()), int((~mism[covered]).sum()), int((bases[covered] == alt[ref_pos][covered]).sum())]
    covered = sites.fixed1[ref_pos]
    out += [int(covered.sum()), int((~mism[covered]).sum()), int((bases[covered] == sites.alt1[ref_pos][covered]).sum())]
    poly = sites.poly[ref_pos]
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


def main(argv=None):
    opts = parse_args(argv)
    for path, what in ((opts.reference, "reference"), (opts.taxonomy, "taxonomy")):
        if not os.path.isfile(path):
            sys.exit(f"ancestry_sites.py: no {what} {path}")
    if opts.full_reference and not os.path.isfile(opts.full_reference):
        print(f"no full reference {opts.full_reference}: the congener sites only", flush=True)
        opts.full_reference = None
    rng = random.Random(opts.seed)
    tax = Taxonomy(opts.taxonomy, opts.heldout)
    sams = sorted(glob.glob(os.path.join(opts.sams, "*", "*", "*.sam*")) + glob.glob(os.path.join(opts.sams, "*.sam*")))
    if not sams:
        sys.exit(f"ancestry_sites.py: no SAM under {opts.sams}")
    records, seen = {}, {}
    for path in sams:
        parts = path.replace("\\", "/").split("/")
        sample = (parts[-3] + ":" if len(parts) >= 3 else "") + os.path.basename(path).split(".")[0]
        records.setdefault(sample, []).extend(records_of(path, tax, seen.setdefault(sample, set())))
    n_records = sum(len(r) for r in records.values())
    head = [f"{len(sams)} SAMs of {len(records)} samples, {n_records} counted records on FP taxa or of FN taxa's own reads"]
    print(head[0], flush=True)

    hits = collections.defaultdict(set)
    for recs in records.values():
        for _, t, g, *_ in recs:
            hits[t].add(g)
    congeners = {t: tax.congeners(t, opts.max_congeners, rng) for t in hits}
    own_keys = {f"{t}_{g}" for t, genes in hits.items() for g in genes}
    wanted = set(own_keys)
    for t, genes in hits.items():
        for g in genes:
            for c in congeners[t]:
                wanted.add(f"{c}_{g}")
    head.append(f"{len(hits)} taxa, {len(own_keys)} taxon-genes, {len(wanted)} copies wanted; taxa without a congener in "
                f"the training database: {sum(1 for t in hits if not congeners[t])}")
    print(head[-1], flush=True)
    copies = load_copies(opts.reference, wanted)
    head.append(f"{len(copies)} copies read from {opts.reference}")
    print(head[-1], flush=True)
    alleles = {}
    if opts.full_reference:
        alleles = load_alleles(opts.full_reference, own_keys, {k: copies[k] for k in own_keys if k in copies},
                               opts.max_alleles)
        head.append(f"alleles: {sum(len(v) for v in alleles.values())} copies of {len(alleles)} taxon-genes from "
                    f"{opts.full_reference}; taxa with any: {len({k.split('_')[0] for k, v in alleles.items() if v})} of {len(hits)}")
        print(head[-1], flush=True)

    sites = {}
    for t, genes in hits.items():
        for g in genes:
            key = f"{t}_{g}"
            if key in copies:
                sites[(t, g)] = Sites(copies[key], [copies[f"{c}_{g}"] for c in congeners[t] if f"{c}_{g}" in copies],
                                      alleles.get(key, []))

    os.makedirs(os.path.dirname(os.path.abspath(opts.out)) or ".", exist_ok=True)
    per = collections.defaultdict(lambda: [0, np.zeros(len(COUNTS), dtype=np.int64)])
    with gzip.open(opts.out + ".fragments.tsv.gz", "wt") as fh:
        fh.write("sample\tqname\ttaxid\trole\trelation\tgene\tnearest_identity\tallele_divergence\t" + "\t".join(COUNTS) + "\n")
        for sample, recs in records.items():
            for qname, t, g, pos, cigar, seq, role, source in recs:
                key = (t, g)
                own = copies.get(f"{t}_{g}")
                if own is None or key not in sites:
                    continue
                s = sites[key]
                ref_pos, bases = aligned_pairs(pos, cigar, seq)
                inside = ref_pos < len(own)
                counts = count_record(s, own, ref_pos[inside], bases[inside])
                rel = tax.relation(source, t)
                fh.write("\t".join(map(str, [sample, qname, t, role, rel, g,
                                             "" if s.nearest_identity is None else f"{s.nearest_identity:.4f}",
                                             f"{s.allele_divergence:.4f}"] + counts)) + "\n")
                group = role if role == "FN own" else f"FP {rel}"
                v = per[(sample, t, group)]
                v[0] += 1
                v[1] += np.array(counts)

    taxa = [(sample, t, group, n, counts) for (sample, t, group), (n, counts) in per.items()]
    with gzip.open(opts.out + ".taxa.tsv.gz", "wt") as fh:
        fh.write("sample\ttaxid\tgroup\trecords\tidentity\t" + "\t".join(COUNTS) + "\n")
        for sample, t, group, n, counts in taxa:
            idn = 1 - counts[1] / counts[0] if counts[0] else float("nan")
            fh.write("\t".join(map(str, [sample, t, group, n, f"{idn:.4f}"] + counts.tolist())) + "\n")
    lines, aucs = summarize(taxa, bool(opts.full_reference))
    with open(opts.out + ".auc.tsv", "w") as fh:
        fh.write("min_sites\tidentity_band\tsignal\ttaxa\tfn\tauc\n")
        for n_min, band, name, n, fn, a in aucs:
            fh.write(f"{n_min}\t{band}\t{name}\t{n}\t{fn}\t{a:.4f}\n")
    with open(opts.out + ".summary.txt", "w") as fh:
        fh.write("\n".join(head + [""] + lines) + "\n")
    print("\n" + "\n".join(lines), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
