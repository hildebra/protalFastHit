#!/usr/bin/env python3
"""ancestry_oracle.py - protal's ancestry sites, strain alleles and polymorphic sites recomputed in Python from a
database's files and a run's SAM records: the oracle of the true-positive test (ancestry_truth_test.py).

A line-by-line port of the C++ that makes the features, so that the test can tell a feature that does what it says
from one that does not:
  AncestrySites.h    UniqueKmers, Chain, CompareCopies, Consensus, Cache::Get (the congeners: the gene's nearest by
                     congener_gaps.tsv, then species_neighbours.tsv, up to 17), Count, ForEachSite
  StrainAlleles.h    the table, FromSamRecord, Explain, BestAllele, Polymorphism, CountSites
  Profiler.h         NoteRecord (the ancestry counts of every best record), NoteAlleles (the alleles, polymorphic and
                     fixed-site counts of the kept ones), FixedWeight, and the features' formulas
The records are each read's (each mate's) first record without the secondary flag, kept when MAPQ >= 4 and its
alignment, soft clips left out, is longer than 50 (Profiler.h PrepareMAPQ); the suspect copies are not left out (the
test's worlds have none).

    import ancestry_oracle as oracle
    db = oracle.Database(folder)            # reference.fna, species_neighbours.tsv, congener_gaps.tsv,
                                            # strain_alleles.tsv of protal --unpack_db
    taxa = oracle.profile_sam(path, db)     # {taxid: Evidence}
    taxa[taxid].features()                  # {feature: value}, as in <profile>.truth_annotated
"""

import bisect
import collections
import gzip
import os
import re

KMER = 12
CONGENERS = 17
CONSENSUS = 0.9
MIN_CONGENERS = 3
TOLERANT_FROM = 6  # kTolerantFrom: from this many compared, one congener with a base of its own does not block a site
MIN_COMPARED = 0.5
MAX_INDEL = 60
MIN_INDEL_LENGTH = 3
INDEL_TOLERANCE = 4      # AncestrySites.h kIndelTolerance (and StrainAlleles.h kIndelTolerance)
INDEL_MARGIN = 5
MAX_LENGTH = 65535
MIN_MAPQ = 4             # Profiler.h m_min_mapq
MIN_ALIGNMENT = 50       # Profiler.h m_min_alignment_length
UNKNOWN = -1.0

CODE = collections.defaultdict(lambda: 4, {"A": 0, "C": 1, "G": 2, "T": 3, "a": 0, "c": 1, "g": 2, "t": 3})
SUB, INS, DEL = 0, 1, 2
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


# ---- ancestry sites ------------------------------------------------------------------------------------------------

def unique_kmers(s):
    """The 12-mers of s (A/C/G/T only) that occur once, as (value, position), sorted by value."""
    out, value, valid, mask = [], 0, 0, (1 << (2 * KMER)) - 1
    for i, ch in enumerate(s):
        c = CODE[ch]
        if c > 3:
            value = valid = 0
            continue
        value = ((value << 2) | c) & mask
        valid += 1
        if valid >= KMER:
            out.append((value, i + 1 - KMER))
    out.sort()
    unique = []
    i = 0
    while i < len(out):
        j = i
        while j < len(out) and out[j][0] == out[i][0]:
            j += 1
        if j == i + 1:
            unique.append(out[i])
        i = j
    return unique


def chain(anchors, main):
    """AncestrySites.h detail::Chain: anchors (position, diagonal) sorted by position (then diagonal)."""
    first = next((i for i, a in enumerate(anchors) if a[1] == main), None)
    if first is None:
        return []

    def accepts(last, i, forward):
        pos, diag = anchors[i]
        shift = diag - last[1]
        if shift < -MAX_INDEL or shift > MAX_INDEL:
            return False
        own = pos - last[0]
        other = own - shift
        if forward and (own <= 0 or other <= 0):
            return False
        if not forward and (own >= 0 or other >= 0):
            return False
        if shift == 0:
            return True
        nxt = i + 1 if forward else i - 1
        return (i + 1 < len(anchors) if forward else i > 0) and anchors[nxt][1] == diag

    before, last = [], anchors[first]
    for i in range(first - 1, -1, -1):
        if accepts(last, i, False):
            before.append(anchors[i])
            last = anchors[i]
    out = before[::-1] + [anchors[first]]
    last = anchors[first]
    for i in range(first + 1, len(anchors)):
        if accepts(last, i, True):
            out.append(anchors[i])
            last = anchors[i]
    return out


class Sites:
    def __init__(self):
        self.positions, self.bases, self.indels = [], [], []  # indels: (position, length)
        self.congener, self.congeners, self.compared, self.identity = 0, 0, 0, 0.0
        self.covered = None

    def empty(self):
        return not self.positions and not self.indels

    def count(self, begin, end):
        return bisect.bisect_left(self.positions, end) - bisect.bisect_left(self.positions, begin)


def compare_copies(own, own_kmers, other):
    """AncestrySites.h CompareCopies: the sites of `own` against `other` with `covered`, or None."""
    if len(own) > MAX_LENGTH or len(own) < KMER or len(other) < KMER:
        return None
    b = unique_kmers(other)
    anchors, i, j = [], 0, 0
    while i < len(own_kmers) and j < len(b):
        if own_kmers[i][0] < b[j][0]:
            i += 1
        elif b[j][0] < own_kmers[i][0]:
            j += 1
        else:
            anchors.append((own_kmers[i][1], own_kmers[i][1] - b[j][1]))
            i += 1
            j += 1
    if len(anchors) < 2:
        return None
    counts = collections.Counter(d for _, d in anchors)
    main, best = 0, 0
    for d in sorted(counts):
        if counts[d] > best:
            best, main = counts[d], d
    anchors.sort()
    links = chain(anchors, main)
    if len(links) < 2:
        return None
    s = Sites()
    covered = bytearray(len(own))
    n_other = len(other)

    def other_at(i, d):
        k = i - d
        return CODE[other[k]] if 0 <= k < n_other else 4

    def cover(frm, to, d):
        for i in range(frm, min(to, len(own))):
            x, y = CODE[own[i]], other_at(i, d)
            if x > 3 or y > 3:
                continue
            covered[i] = 1
            if x != y:
                s.positions.append(i)
                s.bases.append(y)

    for k in range(len(links) - 1):
        (lp, ld), (rp, rd) = links[k], links[k + 1]
        if ld == rd:
            cover(lp, min(len(own), rp + KMER), ld)
            continue
        cover(lp, lp + KMER, ld)
        cover(rp, min(len(own), rp + KMER), rd)
        i = lp + KMER
        while i < rp and CODE[own[i]] < 4 and CODE[own[i]] == other_at(i, ld):
            covered[i] = 1
            i += 1
        j = rp
        while j > i and CODE[own[j - 1]] < 4 and CODE[own[j - 1]] == other_at(j - 1, rd):
            j -= 1
            covered[j] = 1
        shift = rd - ld
        if abs(shift) >= MIN_INDEL_LENGTH:
            s.indels.append((min(i, MAX_LENGTH), shift))
    n = sum(covered)
    if n < MIN_COMPARED * len(own):
        return None
    s.compared, s.congeners = n, 1
    s.identity = 1.0 - len(s.positions) / n
    s.covered = covered
    return s


def agree(best, other, compared):
    """AncestrySites.h Agree: `best` of `compared` congeners carry one state other than the species', `other` another
    one; CONSENSUS of them, or from TOLERANT_FROM all but one that carries a state of its own."""
    if best == 0:
        return False
    if best >= CONSENSUS * compared:
        return True
    return compared >= TOLERANT_FROM and best + 1 == compared and best + other == compared


def consensus(length, comparisons, outgroup=None):
    """AncestrySites.h Consensus (outgroup: the family's consensus base per position, ColumnWeights.h, or None: with
    fewer than MIN_CONGENERS compared, the nearest congener's difference is a site only where its base is the family's)."""
    out = Sites()
    if not comparisons or length > MAX_LENGTH:
        return out
    nearest = comparisons[0]
    for c in comparisons[1:]:
        if nearest.identity < c.identity:
            nearest = c
    compared = [0] * length
    votes = [[0, 0, 0, 0] for _ in range(length)]
    fallback = [-1] * length
    indel_votes = {}
    order = []
    for c in comparisons:
        for i in range(min(length, len(c.covered))):
            compared[i] += c.covered[i] != 0
        for pos, base in zip(c.positions, c.bases):
            if pos >= length:
                continue
            votes[pos][base & 3] += 1
            if c is nearest:
                fallback[pos] = base & 3
        for indel in c.indels:
            if indel not in indel_votes:
                order.append(indel)
            indel_votes[indel] = indel_votes.get(indel, 0) + 1
    for i in range(length):
        if compared[i] == 0:
            continue
        out.compared += 1
        if compared[i] < MIN_CONGENERS:
            polarised = outgroup is not None and i < len(outgroup) and outgroup[i] < 4
            if fallback[i] >= 0 and (not polarised or outgroup[i] == fallback[i]):
                out.positions.append(i)
                out.bases.append(fallback[i])
            continue
        v = votes[i]
        best = 0
        for b in (1, 2, 3):
            if v[b] > v[best]:
                best = b
        if agree(v[best], sum(v) - v[best], compared[i]):
            out.positions.append(i)
            out.bases.append(best)

    def compared_at(indel):
        pos, length_ = indel
        end = pos + (length_ if length_ > 0 else 0)
        n = 0
        for c in comparisons:
            before = pos > 0 and pos - 1 < len(c.covered) and c.covered[pos - 1]
            after = end < len(c.covered) and c.covered[end]
            n += bool(before or after)
        return n

    for indel in order:
        at = compared_at(indel)
        if at < MIN_CONGENERS:
            site = indel in nearest.indels
        else:
            n = indel_votes[indel]
            other = sum(m for o, m in indel_votes.items() if o[0] == indel[0] and o[1] != indel[1])
            site = agree(n, min(other, at - n if at > n else 0), at)
        if site:
            out.indels.append(indel)
    out.indels.sort()
    out.congener, out.identity = nearest.congener, nearest.identity
    out.congeners = min(len(comparisons), 65535)
    out.compared_at = compared  # the oracle's own: the congeners compared at each position (diagnostics)
    return out


class Counts:
    __slots__ = ("sites", "agree", "congener", "indel_sites", "indel_agree", "indel_congener")

    def __init__(self):
        self.sites = self.agree = self.congener = self.indel_sites = self.indel_agree = self.indel_congener = 0


def parse_cigar(cigar):
    return [(int(n), op) for n, op in CIGAR.findall(cigar)]


def count(sites, cigar, pos, seq):
    """AncestrySites.h Count: a record (1-based pos, SEQ in the reference's orientation) at the sites."""
    c = Counts()
    if sites.empty():
        return c
    bases = bool(sites.positions) and seq not in ("", "*")
    gaps = []
    start = pos - 1 if pos > 0 else 0
    ref = start
    query = 0
    for run, op in parse_cigar(cigar):
        if op in "M=":
            if bases:
                n = sites.count(ref, ref + run)
                c.sites += n
                c.agree += n
            ref += run
            query += run
        elif op == "X":
            if bases:
                k = bisect.bisect_left(sites.positions, ref)
                while k < len(sites.positions) and sites.positions[k] < ref + run:
                    q = query + (sites.positions[k] - ref)
                    if q >= len(seq):
                        break
                    c.sites += 1
                    if CODE[seq[q]] == sites.bases[k]:
                        c.congener += 1
                    k += 1
            ref += run
            query += run
        elif op in "DN":
            if sites.indels:
                gaps.append((ref, run, True))
            ref += run
        elif op == "I":
            if sites.indels:
                gaps.append((ref, run, False))
            query += run
        elif op == "S":
            query += run
    for position, length in sites.indels:
        n = abs(length)
        if n < MIN_INDEL_LENGTH:
            continue
        lo, hi = position, position + (n if length > 0 else 0)
        if lo < start + INDEL_MARGIN or hi + INDEL_MARGIN > ref:
            continue
        matched = other = False
        for g_ref, g_len, deletion in gaps:
            g_hi = g_ref + (g_len if deletion else 0)
            if g_hi + INDEL_TOLERANCE < lo or g_ref > hi + INDEL_TOLERANCE:
                continue
            if deletion == (length > 0) and g_len == n:
                matched = True
            else:
                other = True
        if matched:
            c.indel_sites += 1
            c.indel_congener += 1
        elif not other:
            c.indel_sites += 1
            c.indel_agree += 1
    return c


SPECIES_BASE, CONGENER_BASE, OTHER_BASE = 0, 1, 2


def for_each_site(sites, cigar, pos, seq):
    """AncestrySites.h ForEachSite: (position, outcome) of each base site a record covers."""
    if not sites.positions or seq in ("", "*"):
        return
    ref = pos - 1 if pos > 0 else 0
    query = 0
    for run, op in parse_cigar(cigar):
        if op in "M=X":
            k = bisect.bisect_left(sites.positions, ref)
            while k < len(sites.positions) and sites.positions[k] < ref + run:
                p = sites.positions[k]
                if op != "X":
                    yield p, SPECIES_BASE
                else:
                    q = query + (p - ref)
                    if q >= len(seq):
                        break
                    yield p, CONGENER_BASE if CODE[seq[q]] == sites.bases[k] else OTHER_BASE
                k += 1
            ref += run
            query += run
        elif op in "DN":
            ref += run
        elif op in "IS":
            query += run


# ---- strain alleles ------------------------------------------------------------------------------------------------

class Edit:
    """An edit of the representative (or a read's difference): kind SUB/INS/DEL, base (0-3), length."""
    __slots__ = ("pos", "kind", "base", "length")

    def __init__(self, pos, kind, base=0, length=1):
        self.pos, self.kind, self.base, self.length = pos, kind, base, length

    def __repr__(self):
        return f"Edit({self.pos},{'SID'[self.kind]},{self.base},{self.length})"


def matches(a, r):
    if a.kind != r.kind:
        return False
    if a.kind == SUB:
        return a.pos == r.pos and a.base == r.base
    return abs(a.pos - r.pos) <= INDEL_TOLERANCE and a.length == r.length


class Allele:
    __slots__ = ("begin", "end", "edits", "positions")

    def __init__(self, begin, end, edits):
        self.begin, self.end, self.edits = begin, end, edits
        self.positions = [e.pos for e in edits]


def parse_allele(text):
    rng, _, edits = text.partition(":")
    begin, _, end = rng.partition("-")
    out = []
    for item in filter(None, edits.split(",")):
        m = re.fullmatch(r"(\d+)([ACGTid])(\d*)", item)
        if not m:
            raise ValueError(f"not an edit: {item!r}")
        pos, op, length = int(m.group(1)), m.group(2), m.group(3)
        if op == "i":
            out.append(Edit(pos, INS, 0, int(length)))
        elif op == "d":
            out.append(Edit(pos, DEL, 0, int(length)))
        else:
            out.append(Edit(pos, SUB, CODE[op], 1))
    return Allele(int(begin), int(end), out)


def read_alleles(path):
    """strain_alleles.tsv -> {(taxid, gene): [Allele]}."""
    table = {}
    if not os.path.exists(path):
        return table
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line or line.startswith(("#", "taxid")):
                continue
            taxid, gene, alleles = line.split("\t")
            table[(int(taxid), int(gene))] = [parse_allele(a) for a in alleles.split(";")]
    return table


class ReadDiffs:
    __slots__ = ("begin", "end", "columns", "diffs", "other_base", "positions")

    def __init__(self):
        self.begin = self.end = self.columns = 0
        self.diffs, self.other_base = [], []


def from_sam_record(cigar, pos, seq):
    """StrainAlleles.h FromSamRecord -> ReadDiffs, or None (no SEQ, pos 0 or no column)."""
    if seq in ("", "*") or pos == 0:
        return None
    out = ReadDiffs()
    ref, query = pos - 1, 0
    out.begin = ref
    for run, op in parse_cigar(cigar):
        if op in "=M":
            ref += run
            query += run
            out.columns += run
        elif op == "X":
            for i in range(run):
                if query + i < len(seq):
                    b = CODE[seq[query + i]]
                    out.diffs.append(Edit(min(ref + i, 65535), SUB, 0 if b > 3 else b))
                    out.other_base.append(b > 3)
            ref += run
            query += run
            out.columns += run
        elif op == "I":
            out.diffs.append(Edit(min(ref, 65535), INS, 0, min(run, 4095)))
            out.other_base.append(False)
            query += run
            out.columns += run
        elif op in "DN":
            out.diffs.append(Edit(min(ref, 65535), DEL, 0, min(run, 4095)))
            out.other_base.append(False)
            ref += run
            out.columns += run
        elif op == "S":
            query += run
    out.end = ref
    out.positions = [d.pos for d in out.diffs]
    return out if out.columns > 0 else None


def explain(allele, read, used_out=None):
    """StrainAlleles.h Explain -> (explained, contradicted); used_out gets a flag per read diff."""
    used = [False] * len(read.diffs)
    explained = contradicted = 0
    lo, hi = max(read.begin, allele.begin), min(read.end, allele.end)
    if lo < hi:
        k = bisect.bisect_left(allele.positions, lo)
        while k < len(allele.edits) and allele.edits[k].pos < hi:
            e = allele.edits[k]
            k += 1
            substitution = e.kind == SUB
            if not substitution and (e.pos <= read.begin or e.pos >= read.end):
                continue
            frm = e.pos if substitution else max(e.pos - INDEL_TOLERANCE, 0)
            to = e.pos if substitution else e.pos + INDEL_TOLERANCE
            r = bisect.bisect_left(read.positions, frm)
            found = False
            while r < len(read.diffs) and read.diffs[r].pos <= to:
                if not used[r] and not read.other_base[r] and matches(e, read.diffs[r]):
                    used[r] = found = True
                    break
                r += 1
            if found:
                explained += 1
            else:
                contradicted += 1
    if used_out is not None:
        used_out[:] = used
    return explained, contradicted


def best_allele(alleles, read):
    """StrainAlleles.h BestAllele -> (allele index or -1, shift, explained)."""
    best = (-1, 0, 0)
    for a, allele in enumerate(alleles):
        explained, contradicted = explain(allele, read)
        shift = contradicted - explained
        if shift < best[1]:
            best = (a, shift, explained)
    return best


class Polymorphism:
    """StrainAlleles.h Polymorphism, set for a read's span."""

    def __init__(self, alleles, begin, end):
        sites = {}
        self.ranges = [(a.begin, a.end) for a in alleles]
        for a in alleles:
            for e in a.edits:
                if e.pos >= end:
                    break
                if e.kind == SUB:
                    if e.pos >= begin:
                        bits, indel = sites.get(e.pos, (0, False))
                        sites[e.pos] = (bits | 1 << e.base, indel)
                elif e.kind == INS:
                    if e.pos >= begin:
                        sites[e.pos] = (sites.get(e.pos, (0, False))[0], True)
                else:
                    for p in range(max(e.pos, begin), min(e.pos + e.length, end)):
                        sites[p] = (sites.get(p, (0, False))[0], True)
        self.sites = sorted((p, bits, indel) for p, (bits, indel) in sites.items())
        self.index = {p: (bits, indel) for p, bits, indel in self.sites}

    def cover(self, pos):
        return sum(b <= pos < e for b, e in self.ranges)


def count_sites(poly, read):
    """StrainAlleles.h CountSites -> (sites, known, novel)."""
    sites = known = novel = 0
    diffs, n, lo = read.diffs, len(read.diffs), 0

    def end_of(e):
        return e.pos + (e.length if e.kind == DEL else 1)

    for pos, bits, site_indel in poly.sites:
        if pos < read.begin:
            continue
        if pos >= read.end:
            break
        while lo < n and end_of(diffs[lo]) + INDEL_TOLERANCE <= pos:
            lo += 1
        substitution, indel = -1, False
        j = lo
        while j < n and diffs[j].pos <= pos + INDEL_TOLERANCE:
            r = diffs[j]
            if r.kind == SUB:
                if r.pos == pos:
                    substitution = j
            else:
                over = r.kind == DEL and r.pos <= pos < end_of(r)
                if over or (site_indel and abs(r.pos - pos) <= INDEL_TOLERANCE):
                    indel = True
            j += 1
        if substitution >= 0:
            if read.other_base[substitution]:
                continue
            sites += 1
            if bits >> diffs[substitution].base & 1:
                known += 1
            else:
                novel += 1
        elif indel:
            sites += 1
            if site_indel:
                known += 1
            else:
                novel += 1
        else:
            sites += 1
    return sites, known, novel


def fixed_weight(covering):
    return 0 if covering == 0 else 60 * covering // (covering + 1)


# ---- the database and a run ----------------------------------------------------------------------------------------

NO_CODE, NO_BASE, CONSERVED_CODE = 0, 4, 9  # ColumnWeights.h kNoCode, kNoBase, kConservedCode
AMINO = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"


def amino_acid(b1, b2, b3):
    """ColumnWeights.h AminoAcid: the codon's amino acid ('*' a stop), 'X' with another letter."""
    if b1 > 3 or b2 > 3 or b3 > 3:
        return "X"
    return AMINO[16 * b1 + 4 * b2 + b3]


class Columns:
    """ColumnWeights.h Columns: a copy's codes per position."""

    def __init__(self, length):
        self.within = [NO_CODE] * length
        self.among = [NO_CODE] * length
        self.aa = [NO_CODE] * length
        self.consensus = [NO_BASE] * length


def read_column_weights(path):
    """column_weights.tsv -> ({(family, gene): {within, among, aa, consensus}}, {(taxid, gene): (family, runs)})."""
    families, copies = {}, {}
    if not os.path.exists(path):
        return families, copies
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\r\n").split("\t")
            if f[0] == "F":
                families[(int(f[1]), int(f[2]))] = {
                    "within": [int(c, 16) for c in f[5]], "among": [int(c, 16) for c in f[6]], "aa": [int(c, 16) for c in f[7]],
                    "consensus": [CODE[c] for c in f[8]]}
            elif f[0] == "C":
                runs = [tuple(int(x) for x in item.split(":")) for item in f[4].split(",") if item]
                copies[(int(f[1]), int(f[2]))] = (int(f[3]), runs)
    return families, copies


def count_columns(columns, cigar, pos, seq, copy):
    """ColumnWeights.h Count: (columns, aligned, within_aligned, mismatches, within_mismatches, among_mismatches,
    conserved_mismatches, synonymous, nonsynonymous, nonsynonymous_conserved) of a record."""
    c = [0] * 10
    if not columns or not seq or seq == "*":
        return c
    ref = pos - 1 if pos > 0 else 0
    query = 0
    length = len(columns.within)
    for run, op in parse_cigar(cigar):
        if op in "M=X":
            for i in range(run):
                p = ref + i
                if p >= length:
                    break
                c[0] += 1
                w = columns.within[p]
                if w == NO_CODE:
                    continue
                c[1] += 1
                c[2] += w
                if op != "X":
                    continue
                c[3] += 1
                c[4] += w
                c[5] += columns.among[p]
                c[6] += w >= CONSERVED_CODE
                codon = p - p % 3
                if codon + 3 > len(copy) or query + i >= len(seq):
                    continue
                bases = [CODE[copy[codon]], CODE[copy[codon + 1]], CODE[copy[codon + 2]]]
                before = amino_acid(*bases)
                bases[p - codon] = CODE[seq[query + i]]
                after = amino_acid(*bases)
                if before == "X" or after == "X":
                    continue
                if before == after:
                    c[7] += 1
                else:
                    c[8] += 1
                    c[9] += columns.aa[p] >= CONSERVED_CODE
            ref += run
            query += run
        elif op in "IS":
            query += run
        elif op in "DN":
            ref += run
    return c


class Database:
    """The files of protal --unpack_db that the features read."""

    def __init__(self, folder):
        self.genes = {}
        with open(os.path.join(folder, "reference.fna")) as fh:
            name = None
            for line in fh:
                line = line.rstrip("\r\n")
                if line.startswith(">"):
                    name = line[1:].split()[0]
                elif name:
                    taxid, gene = name.split("_")
                    self.genes[(int(taxid), int(gene))] = line
                    name = None
        self.neighbours = {}
        path = os.path.join(folder, "species_neighbours.tsv")
        if os.path.exists(path):
            with open(path) as fh:
                for line in fh:
                    if line.startswith(("#", "taxid")):
                        continue
                    taxid, _, rest = line.rstrip("\r\n").partition("\t")
                    self.neighbours[int(taxid)] = [int(x.split(":")[0]) for x in rest.split(",") if x]
        self.nearest = {}
        path = os.path.join(folder, "congener_gaps.tsv")
        if os.path.exists(path):
            with open(path) as fh:
                for line in fh:
                    if line.startswith(("#", "taxid")):
                        continue
                    taxid, _, rest = line.rstrip("\r\n").partition("\t")
                    for item in filter(None, rest.split(",")):
                        fields = item.split(":")
                        if len(fields) >= 5:
                            self.nearest[(int(taxid), int(fields[0]))] = int(fields[4])
        self.alleles = read_alleles(os.path.join(folder, "strain_alleles.tsv"))
        self.has_alleles = bool(self.alleles)  # GenomeLoader's table Empty(): the features are then -1
        self.weight_families, self.weight_copies = read_column_weights(os.path.join(folder, "column_weights.tsv"))
        self.has_weights = bool(self.weight_families)
        self._columns = {}
        self._sites = {}

    def columns(self, taxid, gene):
        """ColumnWeights.h Table::Expand through the Cache: the copy's columns (a Columns) or None."""
        key = (taxid, gene)
        if key in self._columns:
            return self._columns[key]
        result = None
        own = self.genes.get(key)
        copy = self.weight_copies.get(key)
        if own and copy and len(own) <= MAX_LENGTH:
            family, runs = copy
            f = self.weight_families.get((family, gene))
            if f is not None:
                result = Columns(len(own))
                for begin, column, length in runs:
                    for i in range(length):
                        p, c = begin + i, column + i
                        if p >= len(own) or c >= len(f["within"]):
                            break
                        result.within[p] = f["within"][c]
                        result.among[p] = f["among"][c]
                        result.aa[p] = f["aa"][c // 3] if c // 3 < len(f["aa"]) else NO_CODE
                        result.consensus[p] = f["consensus"][c]
        self._columns[key] = result
        return result

    def sites(self, taxid, gene):
        """AncestrySites.h Cache::Get."""
        key = (taxid, gene)
        if key in self._sites:
            return self._sites[key]
        result = Sites()
        own = self.genes.get(key)
        if own and KMER <= len(own) <= MAX_LENGTH and (self.neighbours or self.nearest):
            candidates = []
            nearest = self.nearest.get(key, 0)
            if nearest and nearest != taxid:
                candidates.append(nearest)
            for n in self.neighbours.get(taxid, []):
                if n != taxid and n not in candidates:
                    candidates.append(n)
            own_kmers = unique_kmers(own)
            comparisons = []
            for congener in candidates:
                if len(comparisons) >= CONGENERS:
                    break
                other = self.genes.get((congener, gene))
                if other is None:
                    continue
                c = compare_copies(own, own_kmers, other)
                if c is None:
                    continue
                c.congener = congener
                comparisons.append(c)
            if comparisons:
                # The family's consensus polarises a small genus's sites (GenomeLoader::AncestrySitesOf).
                columns = self.columns(taxid, gene) if getattr(self, "has_weights", False) else None
                result = consensus(len(own), comparisons, columns.consensus if columns is not None else None)
        self._sites[key] = result
        return result

    def species_sites(self, taxid, gene):
        """GenomeLoader::SpeciesAncestrySitesOf: the sites the species' alleles share (AncestrySites.h SharedBySpecies:
        without the base sites where an allele has another base or an indel, the indel sites with an allele's indel within
        INDEL_TOLERANCE); the same sites for a copy without alleles."""
        sites = self.sites(taxid, gene)
        alleles = self.alleles.get((taxid, gene)) if self.has_alleles else None
        if sites.empty() or not alleles:
            return sites
        key = ("species", taxid, gene)
        if key in self._sites:
            return self._sites[key]
        poly = Polymorphism(alleles, 0, MAX_LENGTH)
        indels = sorted(p for p, _bits, indel in poly.sites if indel)
        out = Sites()
        out.congener, out.congeners, out.compared, out.identity = sites.congener, sites.congeners, sites.compared, sites.identity
        for p, b in zip(sites.positions, sites.bases):
            if p not in poly.index:
                out.positions.append(p)
                out.bases.append(b)
        for indel in sites.indels:
            lo = bisect.bisect_left(indels, indel[0] - INDEL_TOLERANCE if indel[0] > INDEL_TOLERANCE else 0)
            if not (lo < len(indels) and indels[lo] <= indel[0] + INDEL_TOLERANCE):
                out.indels.append(indel)
        out.covered = sites.covered
        self._sites[key] = out
        return out


class Evidence:
    """A taxon's counts (Profiler.h RecordEvidence's ancestry, allele and polymorphic fields)."""

    def __init__(self):
        self.records = self.kept = 0
        self.ancestry_sites = self.ancestry_agree = self.ancestry_congener = 0
        self.ancestry_indel_sites = self.ancestry_indel_agree = self.ancestry_indel_congener = 0
        self.alleles_known = False
        self.allele_records = self.allele_differences = self.allele_explained = self.allele_gain = 0
        self.allele_aligned = 0
        self.poly_sites = self.poly_known = self.poly_novel = 0
        self.fixed_all = self.fixed_all_agree = self.fixed = self.fixed_agree = 0
        # the oracle's own: the fixed-site agreement unweighted, which no feature gives directly
        self.fixed_plain = self.fixed_plain_agree = 0
        # the column weights (RecordEvidence::cw_*)
        self.weights_known = False
        self.cw = [0] * 10
        self.cw_sites_weight = self.cw_agree_weight = 0

    def features(self):
        def share(a, b, empty=UNKNOWN):
            return empty if b == 0 else a / b
        f = {
            "ancestry_sites_per_record": share(self.ancestry_sites, self.records, 0.0),
            "ancestry_agreement": share(self.ancestry_agree, self.ancestry_sites),
            "ancestry_congener_share": share(self.ancestry_congener, self.ancestry_sites),
            "ancestry_indel_sites_per_record": share(self.ancestry_indel_sites, self.records, 0.0),
            "ancestry_indel_congener_share": share(self.ancestry_indel_congener, self.ancestry_indel_sites),
        }
        if not self.alleles_known:
            for k in ("allele_explained_share", "allele_identity_gain", "allele_copy_share", "polymorphic_known_share",
                      "polymorphic_novel_share", "ancestry_fixed_gain", "ancestry_fixed_agreement", "allele_sites_per_kb",
                      "ancestry_fixed_share"):
                f[k] = UNKNOWN
        else:
            f["allele_explained_share"] = share(self.allele_explained, self.allele_differences, 0.0)
            f["allele_identity_gain"] = share(self.allele_gain, self.allele_aligned, 0.0)
            f["allele_copy_share"] = share(self.allele_records, self.kept, 0.0)
            f["polymorphic_known_share"] = share(self.poly_known, self.poly_sites, 0.0)
            f["polymorphic_novel_share"] = share(self.poly_novel, self.poly_sites, 0.0)
            f["ancestry_fixed_gain"] = (0.0 if self.fixed == 0 or self.fixed_all == 0 else
                                        self.fixed_agree / self.fixed - self.fixed_all_agree / self.fixed_all)
            f["ancestry_fixed_agreement"] = share(self.fixed_agree, self.fixed, 0.0)
            f["allele_sites_per_kb"] = share(1000 * self.poly_sites, self.allele_aligned, 0.0)
            f["ancestry_fixed_share"] = share(self.fixed, self.fixed_all, 0.0)
        f["oracle_fixed_agreement"] = share(self.fixed_plain_agree, self.fixed_plain)
        f["oracle_records"] = self.records
        weights = ["conserved_mismatch_ratio", "conserved_mismatch_rate", "ancestry_agreement_weighted", "nonsynonymous_share",
                   "nonsynonymous_conserved_rate", "column_weight_coverage"]
        if not self.weights_known:
            for k in weights:
                f[k] = UNKNOWN
        else:
            cols, aligned, within_aligned, mism, within_mism, _among, conserved, syn, nonsyn, nonsyn_cons = self.cw
            f["conserved_mismatch_ratio"] = (0.0 if mism == 0 or aligned == 0 or within_aligned == 0 else
                                             (within_mism / mism) / (within_aligned / aligned))
            f["conserved_mismatch_rate"] = share(1000 * conserved, aligned, 0.0)
            f["ancestry_agreement_weighted"] = share(self.cw_agree_weight, self.cw_sites_weight)
            f["nonsynonymous_share"] = share(nonsyn, syn + nonsyn, 0.0)
            f["nonsynonymous_conserved_rate"] = share(1000 * nonsyn_cons, aligned, 0.0)
            f["column_weight_coverage"] = share(aligned, cols, 0.0)
        return f


def note_record(e, db, taxid, gene, cigar, pos, seq, mapq, kept):
    """Profiler.h NoteRecord (the ancestry counts) and NoteAlleles."""
    e.records += 1
    e.kept += kept
    sites = db.species_sites(taxid, gene)  # the ancestry counts: the sites the species' alleles share
    if not sites.empty():
        c = count(sites, cigar, pos, seq)
        e.ancestry_sites += c.sites
        e.ancestry_agree += c.agree
        e.ancestry_congener += c.congener
        e.ancestry_indel_sites += c.indel_sites
        e.ancestry_indel_agree += c.indel_agree
        e.ancestry_indel_congener += c.indel_congener
    if getattr(db, "has_weights", False):
        e.weights_known = True
        columns = db.columns(taxid, gene)
        if columns is not None:
            cw = count_columns(columns, cigar, pos, seq, db.genes.get((taxid, gene), ""))
            for k in range(10):
                e.cw[k] += cw[k]
            if not sites.empty():
                for p, outcome in for_each_site(sites, cigar, pos, seq):
                    w = max(1, columns.among[p]) if p < len(columns.among) else 1
                    e.cw_sites_weight += w
                    if outcome == SPECIES_BASE:
                        e.cw_agree_weight += w
    if not kept or not db.has_alleles:
        return
    e.alleles_known = True
    read = from_sam_record(cigar, pos, seq)
    if read is None:
        return
    e.allele_differences += len(read.diffs)
    e.allele_aligned += read.columns
    alleles = db.alleles.get((taxid, gene))
    if not alleles:
        return
    _, shift, explained = best_allele(alleles, read)
    e.allele_records += 1
    e.allele_explained += explained
    e.allele_gain += -shift
    poly = Polymorphism(alleles, read.begin, read.end)
    n, known, novel = count_sites(poly, read)
    e.poly_sites += n
    e.poly_known += known
    e.poly_novel += novel
    for p, outcome in for_each_site(db.sites(taxid, gene), cigar, pos, seq):  # all the sites (NoteAlleles)
        w = fixed_weight(poly.cover(p))
        if w == 0:
            continue
        agree = outcome == SPECIES_BASE
        e.fixed_all += w
        e.fixed_all_agree += w if agree else 0
        if p in poly.index:
            continue
        e.fixed += w
        e.fixed_agree += w if agree else 0
        e.fixed_plain += 1
        e.fixed_plain_agree += agree


def open_sam(path):
    if path.endswith(".zst"):
        import subprocess
        proc = subprocess.Popen(["zstd", "-dc", path], stdout=subprocess.PIPE, text=True)
        return proc.stdout
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def best_records(path):
    """Each read's (each mate's) first record without the secondary flag: (qname, mate, rname, pos, mapq, cigar, seq,
    flag)."""
    seen = set()
    with open_sam(path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            f = line.rstrip("\n").split("\t", 11)
            if len(f) < 10:
                continue
            flag = int(f[1])
            if flag & 0x4 or flag & 0x100 or f[2] == "*":
                continue
            mate = 1 if flag & 0x40 else 2 if flag & 0x80 else 0
            key = (f[0], mate)
            if key in seen:
                continue
            seen.add(key)
            yield f[0], mate, f[2], int(f[3]), int(f[4]), f[5], f[9], flag


def clipped_length(cigar):
    """AlignmentUtils.h CompressedCigarInfo's clipped_alignment_length: every op but the soft clips."""
    return sum(n for n, op in parse_cigar(cigar) if op != "S")


def profile_sam(path, db, taxa=None):
    """{taxid: Evidence} of a SAM's best records (only of `taxa`, if given)."""
    out = {}
    for _, _, rname, pos, mapq, cigar, seq, _ in best_records(path):
        taxid, _, gene = rname.partition("_")
        taxid, gene = int(taxid), int(gene)
        if taxa is not None and taxid not in taxa:
            continue
        kept = mapq >= MIN_MAPQ and clipped_length(cigar) > MIN_ALIGNMENT
        note_record(out.setdefault(taxid, Evidence()), db, taxid, gene, cigar, pos, seq, mapq, kept)
    return out
