#!/usr/bin/env python3
"""insilico_strains.py - in-silico strains for the species of a simulation genome table that have one genome.

The presence models learn what a strain looks like from species simulated from another genome than the database's
representative. GTDB has other genomes for some species only; in the r226 pool 2,017 of 7,998 species had one, so
they were always simulated from their own reference, and a model given GTDB's cluster sizes (the priors) learned
that a divergent read cloud on a one-genome species is a relative the database lacks: a rule of the simulation,
not of the species (docs/claude/2026-10-04-r226-v10-evaluation). This script gives such species a strain: a copy
of the representative with substitutions as a strain's would be, so that the simulator (which draws a species,
then one of its genomes) simulates them from either.

The mutations:
- Divergence. Each strain gets a marker-gene divergence D (the median marker gene's share of differing bases),
  drawn from the real strains of the table: each other genome's marker copies placed by their k-mer trace
  (gene_neighbours.py's gene_positions.tsv, kmer_share s of K-mers: divergence 1 - s^(1/K)), the median over its
  genes. The genome as a whole differs by D / --marker-scale (marker genes are conserved: at r226 the reads of real
  strains differed from their representative's by 0.43 (median) of 1 - their cluster's mean ANI), at most
  1 - --min-ani. --ani MIN-MAX draws the genome's ANI uniformly instead (marker divergence = (1 - ANI) x scale).
- Genes. Each marker gene differs by D x its conservation factor, estimated like protal's (GeneConservation.h)
  from the same placements: per strain the gene's divergence over the strain's median gene, the median over the
  strains, the median gene scaled to 1, kept within 0.25-4. Other open reading frames (stop to stop, 300 bases or
  more, the longest first where they overlap) diverge at a gamma-distributed rate per frame (shape 2), the rest of
  the genome at 1.3 times the coding rate, scaled so that the genome differs by its divergence.
- Codons. In a coding frame a substitution (transitions --kappa times as likely as each transversion) that keeps
  the amino acid is kept, one that changes it with probability --omega (dN/dS), one that makes a stop codon never;
  the proposals are scaled per gene so that the gene still differs by its target. Most differences fall on third
  codon positions, as a strain's do (protal's third_position_share). --omega auto (the default) measures how many
  of the real strains' substitutions in their marker genes fall on third positions (a sample of the table's real
  strains against their representatives, the copies compared along their shared 12-mers) and takes the omega that
  gives the same share when the representatives' genes are mutated here: at r226 v15 the in-silico strains made
  with 0.15 had 0.06-0.10 less of their differences on third positions than the real strains, a spectrum the
  presence models could learn (docs/claude/2026-10-07-error-read-signatures).
Substitutions only: the strain has the representative's length and its genes' positions.

Writes OUT_DIR/<name>.fna.gz per strain (name: insilico_ and the representative's accession with '_' for '.', so
that no accession pattern takes it for the representative; its contigs <name>_<the representative's contig>, so
that its reads, named by their contig, are told from the representative's), the genome table --output
(the input's rows and one per strain: name, taxonomy, FASTA, length) and OUT_DIR/insilico_strains.tsv (per strain:
representative, species, genome and marker divergence drawn, substitutions made, the coding share of the genome,
the marker divergence reached on the placed genes).

Usage:
  insilico_strains.py --genome-table genomes.tsv --output genomes_simulated.tsv --out-dir insilico_strains
      [--positions gene_positions.tsv --taxonomy internal_taxonomy.dmp] [--share 1] [--ani MIN-MAX] [-t 8]
"""

import argparse
import concurrent.futures
import gzip
import hashlib
import os
import statistics
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "mini_db"))
from gtdb_to_protal_db import normalize_accession  # noqa: E402

PREFIX = "insilico_"  # names of in-silico strains (collect_training_data.py's meta_insilico_strain)
SUMMARY = "insilico_strains.tsv"
KMER = 24  # gene_neighbours.py's KMER: kmer_share is of KMER-mers
MIN_ORF = 300
NONCODING_RATE = 1.3
GAMMA_SHAPE = 2.0
MIN_STRAIN_GENES = 10  # a strain's genes placed, to count its median gene
MIN_TYPICAL = 0.002  # a strain's median gene at least this far from the representative, for the factors
MIN_FACTOR, MAX_FACTOR = 0.25, 4.0
MIN_FACTOR_STRAINS = 3
OMEGA_DEFAULT = 0.15  # without real strains to calibrate on
OMEGA_GRID = (0.01, 0.02, 0.035, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.5, 1.0)
SPECTRUM_STRAINS = 200  # real strains compared with their representatives for the spectrum
SPECTRUM_KMER = 12  # the copies are compared along their shared 12-mers (no alignment)
MIN_SPECTRUM_SUBSTITUTIONS = 2000  # fewer, and the measured share is not trusted

CODE = np.full(256, 4, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    CODE[_b] = _i
    CODE[_b + 32] = _i
BASES = np.frombuffer(b"ACGTN", dtype=np.uint8)
# The standard genetic code, codon index 16 x first + 4 x second + third (A, C, G, T = 0..3); '*' a stop.
AMINO = np.frombuffer(b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF", dtype=np.uint8)
STOP = AMINO == ord("*")


def is_insilico(name):
    """Whether a genome table's name is an in-silico strain of this script."""
    return name.startswith(PREFIX)


def strain_name(accession):
    return PREFIX + accession.replace(".", "_")


def read_table(path):
    with open(path) as fh:
        return [line.rstrip("\n").split("\t") for line in fh if line.strip() and not line.startswith("#")]


def one_genome_species(rows):
    """The rows of the species with one genome (by the last rank of the taxonomy), in the table's order."""
    rows = [r for r in rows if len(r) >= 3 and ";s__" in r[1]]  # not a header
    count = {}
    for r in rows:
        count[r[1].split(";")[-1]] = count.get(r[1].split(";")[-1], 0) + 1
    return [r for r in rows if count[r[1].split(";")[-1]] == 1 and not is_insilico(r[0])]


def representatives(taxonomy):
    """Accessions of the species' representatives, from internal_taxonomy.dmp."""
    with open(taxonomy) as fh:
        header = next(fh).rstrip("\n").split("\t")
        rank, rep = header.index("rank"), header.index("rep_genome")
        return {normalize_accession(f[rep]) for f in (line.rstrip("\n").split("\t") for line in fh)
                if f[rank] == "species" and f[rep]}


def read_positions(path, wanted, reps):
    """gene_positions.tsv -> ({accession: [(contig, gene, start, end, strand)]} of the accessions in `wanted`,
    {strain accession: {gene: divergence}} of the genomes that are no representative)."""
    genes, strains = {}, {}
    with open(path) as fh:
        header = None
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if header is None:
                header = {name: i for i, name in enumerate(f)}
                continue
            acc = normalize_accession(f[header["accession"]])
            if acc in wanted:
                genes.setdefault(acc, []).append((f[header["contig"]], int(f[header["gene"]]), int(f[header["start"]]),
                                                  int(f[header["end"]]), f[header["strand"]]))
            if reps is not None and acc not in reps:
                share = 1.0 if f[header["placed"]] == "exact" else float(f[header["kmer_share"]])
                strains.setdefault(acc, {})[int(f[header["gene"]])] = 1.0 - max(share, 1e-9) ** (1.0 / KMER)
    return genes, strains


def strain_divergences(strains):
    """Each real strain's median marker gene divergence (MIN_STRAIN_GENES genes placed or more)."""
    return [statistics.median(d.values()) for d in strains.values() if len(d) >= MIN_STRAIN_GENES]


def conservation_factors(strains):
    """{gene: factor}: per strain each gene's divergence over the strain's median gene (strains whose median gene is
    MIN_TYPICAL or more from the representative), the median over the strains, scaled so the median gene has 1."""
    ratios = {}
    for d in strains.values():
        if len(d) < MIN_STRAIN_GENES:
            continue
        typical = statistics.median(d.values())
        if typical < MIN_TYPICAL:
            continue
        for gene, div in d.items():
            ratios.setdefault(gene, []).append(div / typical)
    raw = {g: statistics.median(r) for g, r in ratios.items() if len(r) >= MIN_FACTOR_STRAINS}
    if not raw:
        return {}
    middle = statistics.median(raw.values()) or 1.0
    return {g: min(MAX_FACTOR, max(MIN_FACTOR, v / middle)) for g, v in raw.items()}


def gene_copy(contigs, contig, start, end, strand):
    """A marker gene's codes (0-3, 4 other) in reading orientation from a genome's contigs ({name: codes}), or None
    if the placement is outside the contig. start, end: 1-based, inclusive."""
    codes = contigs.get(contig)
    if codes is None or not 1 <= start <= end <= len(codes):
        return None
    copy = codes[start - 1:end]
    if strand == "-":
        copy = np.where(copy < 4, 3 - copy.astype(np.int16), 4).astype(np.uint8)[::-1]
    return copy


def kmer_codes(codes, k):
    """The k-mers of codes as integers (4 ** k values; -1 where a k-mer holds another letter)."""
    if len(codes) < k:
        return np.zeros(0, dtype=np.int64)
    windows = np.lib.stride_tricks.sliding_window_view(codes, k)
    powers = 4 ** np.arange(k - 1, -1, -1, dtype=np.int64)
    values = (np.minimum(windows, 3).astype(np.int64) * powers).sum(axis=1)
    return np.where((windows == 4).any(axis=1), -1, values)


def substitutions(rep, other, k=SPECTRUM_KMER):
    """The substitutions between a representative's gene copy and another genome's (codes in reading orientation),
    without an alignment: the k-mers unique to each copy pair them on their main diagonal, and the bases are compared
    on the stretches between two paired k-mers that lie on that diagonal (an indel moves the diagonal and ends the
    stretch). Returns (positions on the representative's copy that differ, bases compared)."""
    a, b = kmer_codes(rep, k), kmer_codes(other, k)
    if not len(a) or not len(b):
        return np.zeros(0, dtype=np.int64), 0
    ua, ia, ca = np.unique(a, return_index=True, return_counts=True)
    ub, ib, cb = np.unique(b, return_index=True, return_counts=True)
    keep_a = (ca == 1) & (ua >= 0)
    keep_b = (cb == 1) & (ub >= 0)
    shared, pa, pb = np.intersect1d(ua[keep_a], ub[keep_b], assume_unique=True, return_indices=True)
    if len(shared) < 2:
        return np.zeros(0, dtype=np.int64), 0
    pos_a, pos_b = ia[keep_a][pa], ib[keep_b][pb]
    diagonals, counts = np.unique(pos_a - pos_b, return_counts=True)
    d = int(diagonals[np.argmax(counts)])
    on = pos_a - pos_b == d
    pos_a = np.sort(pos_a[on])
    compared = np.zeros(len(rep), dtype=bool)
    for p, q in zip(pos_a[:-1], pos_a[1:]):  # consecutive paired k-mers on the diagonal: the stretch between
        compared[p:q + k] = True
    lo, hi = max(0, d), min(len(rep), len(other) + d)
    if hi <= lo:
        return np.zeros(0, dtype=np.int64), 0
    differ = np.zeros(len(rep), dtype=bool)
    differ[lo:hi] = rep[lo:hi] != other[lo - d:hi - d]
    compared[:lo] = False
    compared[hi:] = False
    compared &= (rep < 4)
    compared[lo:hi] &= other[lo - d:hi - d] < 4  # another letter (N) on either copy: not compared
    differ &= compared
    return np.flatnonzero(differ), int(compared.sum())


def spectrum_pair(job):
    """One real strain against its representative: (substitutions on third codon positions, substitutions, bases
    compared, the representative's gene copies compared (for the calibration))."""
    strain_path, rep_path, strain_genes, rep_genes = job
    strain = {h.split()[0] if h.split() else h: CODE[np.frombuffer(s, dtype=np.uint8)] for h, s in read_fasta(strain_path)}
    rep = {h.split()[0] if h.split() else h: CODE[np.frombuffer(s, dtype=np.uint8)] for h, s in read_fasta(rep_path)}
    by_gene = {g[1]: g for g in strain_genes}
    third = subs = compared = 0
    copies = []
    for contig, gene, s, e, strand in rep_genes:
        other = by_gene.get(gene)
        if other is None:
            continue
        a = gene_copy(rep, contig, s, e, strand)
        b = gene_copy(strain, other[0], other[2], other[3], other[4])
        if a is None or b is None or len(a) < 3 * SPECTRUM_KMER:
            continue
        positions, n = substitutions(a, b)
        if n < len(a) // 2:  # the copies do not pair: another gene, or a poor placement
            continue
        third += int((positions % 3 == 2).sum())
        subs += len(positions)
        compared += n
        if len(copies) < 8:
            copies.append(a[:len(a) - len(a) % 3])
    return third, subs, compared, copies


def real_spectrum(rows, strains, positions_path, reps, rng, threads):
    """The share of the real strains' substitutions in their marker genes that fall on third codon positions, from a
    sample of SPECTRUM_STRAINS strains (gene_positions.tsv's non-representative genomes with MIN_STRAIN_GENES genes
    placed) against their species' representatives. Returns (share or None, substitutions, strains compared, some
    representatives' gene copies)."""
    by_acc = {normalize_accession(r[0]): r for r in rows if len(r) >= 3 and ";s__" in r[1]}
    rep_of = {}
    for acc, r in by_acc.items():
        if acc in reps:
            rep_of.setdefault(r[1].split(";")[-1], acc)
    candidates = sorted(acc for acc, d in strains.items() if len(d) >= MIN_STRAIN_GENES and acc in by_acc
                        and rep_of.get(by_acc[acc][1].split(";")[-1]) not in (None, acc))
    if not candidates:
        return None, 0, 0, []
    sample = [candidates[i] for i in rng.permutation(len(candidates))[:SPECTRUM_STRAINS]]
    pairs = {acc: rep_of[by_acc[acc][1].split(";")[-1]] for acc in sample}
    genes, _ = read_positions(positions_path, set(pairs) | set(pairs.values()), None)
    jobs = [(by_acc[acc][2], by_acc[rep][2], genes.get(acc, []), genes.get(rep, []))
            for acc, rep in pairs.items() if genes.get(acc) and genes.get(rep)]
    third = subs = 0
    copies, used = [], 0
    with concurrent.futures.ProcessPoolExecutor(max(1, threads)) as pool:
        for t, s, c, cp in pool.map(spectrum_pair, jobs, chunksize=4):
            if c:
                used += 1
            third, subs = third + t, subs + s
            if len(copies) < 24:
                copies.extend(cp)
    return (third / subs if subs >= MIN_SPECTRUM_SUBSTITUTIONS else None), subs, used, copies


def simulated_third_share(copies, omega, kappa, rng, rate=0.03):
    """The share of the substitutions mutate() makes on third codon positions, on a contig of these gene copies
    (each a frame on the plus strand at `rate`)."""
    codes = np.concatenate(copies)
    frames, first = [], 0
    for c in copies:
        frames.append((first, first + len(c) - 1, "+", rate))
        first += len(c)
    owner, cpos, minus, fixed, relative = annotate(codes, frames, rng)
    new, *_ = mutate(codes, owner, cpos, minus, fixed, len(frames), rng, kappa, omega)
    changed = new != codes
    return float((changed & (cpos == 2)).sum() / max(1, changed.sum()))


def calibrate_omega(target, copies, kappa, rng, grid=OMEGA_GRID):
    """The omega at which mutate() puts `target` of its substitutions on third codon positions: the share is measured
    on `copies` at each omega of the grid (two draws each) and the target is interpolated in log(omega), clamped to
    the grid."""
    if not copies or target is None:
        return OMEGA_DEFAULT
    shares = [np.mean([simulated_third_share(copies, w, kappa, rng) for _ in range(2)]) for w in grid]
    if target >= shares[0]:
        return grid[0]
    if target <= shares[-1]:
        return grid[-1]
    for (w0, s0), (w1, s1) in zip(zip(grid, shares), zip(grid[1:], shares[1:])):
        if s1 <= target <= s0 and s0 != s1:
            f = (s0 - target) / (s0 - s1)
            return float(np.exp(np.log(w0) + f * (np.log(w1) - np.log(w0))))
    return OMEGA_DEFAULT


def read_fasta(path):
    """[(header line without '>', sequence bytes)]."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rb") as fh:
        data = fh.read()
    records = []
    for chunk in data.split(b"\n>"):
        chunk = chunk.lstrip(b">")
        if not chunk.strip():
            continue
        header, _, body = chunk.partition(b"\n")
        records.append((header.decode().rstrip("\r"), body.replace(b"\n", b"").replace(b"\r", b"")))
    return records


def orfs(codes):
    """Open reading frames of MIN_ORF bases or more, stop to stop, of a contig's codes (0-3, 4 other) on both strands:
    [(first, last, strand)] in contig coordinates (0-based, inclusive), the longest first."""
    n = len(codes)
    found = []
    reverse = np.where(codes < 4, 3 - codes.astype(np.int16), 4).astype(np.uint8)[::-1]
    for strand, seq in (("+", codes), ("-", reverse)):
        for frame in range(3):
            m = (n - frame) // 3
            if m <= 0:
                continue
            c = seq[frame:frame + 3 * m].reshape(m, 3).astype(np.int32)
            valid = (c < 4).all(axis=1)
            index = np.where(valid, 16 * c[:, 0] + 4 * c[:, 1] + c[:, 2], 0)
            stops = np.flatnonzero(~valid | STOP[index])
            edges = np.concatenate(([-1], stops, [m]))
            for a, b in zip(edges[:-1] + 1, edges[1:]):
                if 3 * (b - a) >= MIN_ORF:
                    lo, hi = frame + 3 * a, frame + 3 * b - 1  # in seq coordinates
                    found.append((lo, hi) if strand == "+" else (n - 1 - hi, n - 1 - lo))
                    found[-1] += (strand,)
    found.sort(key=lambda o: o[0] - o[1])
    return found


def annotate(codes, frames, rng):
    """A contig's sites by frame. frames: [(first, last, strand, rate or None)], the marker genes (with their rate)
    before the open reading frames (None: a gamma-distributed rate, mean 1, drawn here, relative to the coding rate);
    the sites of a frame already taken by an earlier one stay with it, and a frame that keeps less than half of its
    sites is left out. Returns (owner: the frame of each site or -1, codon position or -1, on the minus strand,
    fixed rate (the markers'), relative rate (the other frames', NONCODING_RATE elsewhere, 0 for other letters))."""
    n = len(codes)
    fixed = np.zeros(n)
    relative = np.zeros(n)
    cpos = np.full(n, -1, dtype=np.int8)
    minus = np.zeros(n, dtype=bool)
    owner = np.full(n, -1, dtype=np.int64)
    for k, (first, last, strand, r) in enumerate(frames):
        sites = np.arange(first, last + 1)
        free = sites[owner[sites] < 0]
        if len(free) < len(sites) // 2:
            continue
        owner[free] = k
        if r is not None:
            fixed[free] = r
        else:
            relative[free] = rng.gamma(GAMMA_SHAPE, 1.0 / GAMMA_SHAPE)
        cpos[free] = (free - first) % 3 if strand == "+" else (last - free) % 3
        minus[free] = strand == "-"
    relative[(cpos < 0) & (codes < 4)] = NONCODING_RATE
    return owner, cpos, minus, fixed, relative


def mutate(codes, owner, cpos, minus, rate, n_frames, rng, kappa, omega):
    """Substitutes in a contig's codes (annotate's sites, `rate` per site). Returns (new codes, substitutions,
    coding sites, substitutions per frame, coding sites per frame)."""
    n = len(codes)
    coding = cpos >= 0
    noncoding = ~coding & (codes < 4)
    new = codes.copy()
    # Non-coding: every proposal kept.
    hit = np.flatnonzero(noncoding & (rng.random(n) < rate))
    new[hit] = alternative(codes[hit], rng, kappa)
    changed = [hit]
    # Coding: the codon each site is in, in reading orientation.
    sites = np.flatnonzero(coding)
    cp = cpos[sites].astype(np.int64)
    mi = minus[sites]
    start = np.where(mi, sites + cp, sites - cp)  # the codon's first base in reading direction
    step = np.where(mi, -1, 1)
    inside = (start + 2 * step >= 0) & (start + 2 * step < n)
    sites, cp, mi, start, step = sites[inside], cp[inside], mi[inside], start[inside], step[inside]
    trio = np.stack([codes[start], codes[start + step], codes[start + 2 * step]], axis=1).astype(np.int64)
    trio = np.where(mi[:, None] & (trio < 4), 3 - trio, trio)
    valid = (trio < 4).all(axis=1)
    sites, cp, mi, start, trio = sites[valid], cp[valid], mi[valid], start[valid], trio[valid]
    codon = 16 * trio[:, 0] + 4 * trio[:, 1] + trio[:, 2]
    weight = 4 ** (2 - cp)
    old = trio[np.arange(len(cp)), cp]  # the site's base in reading orientation

    p_ts = kappa / (kappa + 2.0)
    acc_ts = kept_share(codon, weight, old, old ^ 2, omega)
    acc_v1 = kept_share(codon, weight, old, old ^ 1, omega)
    acc_v2 = kept_share(codon, weight, old, old ^ 3, omega)
    mean_accept = p_ts * acc_ts + (1 - p_ts) / 2 * (acc_v1 + acc_v2)
    # Per frame, proposals scaled so that the frame's expected substitutions are its rate times its sites.
    frame_of = owner[sites]
    sums = np.bincount(frame_of, weights=mean_accept, minlength=n_frames)
    counts = np.bincount(frame_of, minlength=n_frames)
    frame_accept = np.where(counts > 0, sums / np.maximum(counts, 1), 1.0)
    propose = rate[sites] / np.maximum(frame_accept[frame_of], 1e-6)
    drawn = rng.random(len(sites)) < np.minimum(propose, 1.0)
    idx = np.flatnonzero(drawn)
    alt = alternative(old[idx], rng, kappa)
    keep_p = kept_share(codon[idx], weight[idx], old[idx], alt, omega)
    kept = idx[rng.random(len(idx)) < keep_p]
    # Two substitutions in one codon were each judged against the original codon: the codon both make must not be a
    # stop (else neither is made).
    alt_kept = alternative_lookup(idx, alt, kept)
    # One codon: its first base and strand (a plus- and a minus-strand codon can start at one base where frames on
    # opposite strands meet).
    _, inverse = np.unique(2 * start[kept] + mi[kept], return_inverse=True)
    delta = np.bincount(inverse, weights=weight[kept] * (alt_kept.astype(np.int64) - old[kept])).astype(np.int64)
    fine = ~STOP[codon[kept] + delta[inverse]]
    kept, alt_kept = kept[fine], alt_kept[fine]
    top = np.where(mi[kept], 3 - alt_kept, alt_kept)
    new[sites[kept]] = top
    changed.append(sites[kept])
    per_frame = np.bincount(owner[sites[kept]], minlength=n_frames)
    total = sum(len(c) for c in changed)
    return new, total, int(coding.sum()), per_frame, counts


def kept_share(codon, weight, old, alt, omega):
    """The probability that a substitution of a codon's base `old` (at the place `weight` = 4 ** (2 - codon position))
    by `alt` is kept: 1 if it keeps the amino acid, omega if it changes it, 0 if it makes a stop codon."""
    changed = codon + weight * (alt.astype(np.int64) - old)
    same = AMINO[changed] == AMINO[codon]
    return np.where(STOP[changed], 0.0, np.where(same, 1.0, omega))


def alternative_lookup(idx, alt, kept):
    """The alternative drawn for each of `kept` (a subset of idx, both sorted)."""
    return alt[np.searchsorted(idx, kept)]


def alternative(base, rng, kappa):
    """Another base for each of `base` (0-3): the transition (A<->G, C<->T) with probability kappa / (kappa + 2),
    else one of the two transversions."""
    u = rng.random(len(base))
    p_ts = kappa / (kappa + 2.0)
    return np.where(u < p_ts, base ^ 2, np.where(u < (1 + p_ts) / 2, base ^ 1, base ^ 3)).astype(base.dtype)


def make_strain(job):
    """Writes one strain; returns its summary row."""
    (acc, species, path, out, genes, factors, genome_div, marker_div, seed, kappa, omega) = job
    rng = np.random.default_rng(seed)
    records = read_fasta(path)
    by_contig = {}
    for g in genes:
        by_contig.setdefault(g[0], []).append(g)
    total_sites = sum(len(seq) for _, seq in records)
    contigs = []
    for header, seq in records:
        codes = CODE[np.frombuffer(seq, dtype=np.uint8)]
        name = header.split()[0] if header.split() else header
        markers = [(s - 1, e - 1, strand, marker_div * factors.get(gene, 1.0), gene)
                   for _, gene, s, e, strand in by_contig.get(name, []) if 1 <= s <= e <= len(codes)]
        frames = [(m[0], m[1], m[2], m[3]) for m in markers] + [(a, b, s, None) for a, b, s in orfs(codes)]
        contigs.append((header, codes, markers, len(frames), annotate(codes, frames, rng)))
    # The coding rate, so that the genome differs by genome_div: the marker sites at their rates, the other frames
    # at the coding rate times their gamma draw, the non-coding sites at NONCODING_RATE times it.
    fixed = sum(c[4][3].sum() for c in contigs)
    relative = sum(c[4][4].sum() for c in contigs)
    coding_rate = max(0.0, genome_div * total_sites - fixed) / max(relative, 1e-9)
    tmp = out + ".partial"
    subs = coding = 0
    reached = []
    with gzip.open(tmp, "wb", compresslevel=1) as fh:
        for header, codes, markers, n_frames, (owner, cpos, minus, fixed, relative) in contigs:
            new, n_subs, n_coding, per_frame, counts = mutate(codes, owner, cpos, minus, fixed + coding_rate * relative,
                                                              n_frames, rng, kappa, omega)
            subs += n_subs
            coding += n_coding
            for k, m in enumerate(markers):
                if counts[k]:
                    reached.append(per_frame[k] / counts[k] / max(factors.get(m[4], 1.0), 1e-9))
            fh.write(f">{strain_name(acc)}_{header} in-silico strain of {acc}\n".encode())
            text = BASES[new].tobytes()
            for i in range(0, len(text), 80):
                fh.write(text[i:i + 80] + b"\n")
    os.replace(tmp, out)
    return [strain_name(acc), acc, species, f"{genome_div:.5f}", f"{marker_div:.5f}", str(subs),
            f"{coding / max(total_sites, 1):.3f}", f"{statistics.median(reached):.5f}" if reached else "",
            str(total_sites)]


def parse_ani(text):
    if text in (None, "", "pool"):
        return None
    lo, _, hi = text.partition("-")
    lo, hi = float(lo), float(hi or lo)
    lo, hi = (lo / 100, hi / 100) if hi > 1 else (lo, hi)
    if not 0.5 < lo <= hi < 1:
        raise argparse.ArgumentTypeError(f"--ani {text}: MIN-MAX between 0.5 and 1 (or 50-100)")
    return lo, hi


def seed_of(seed, acc):
    return int.from_bytes(hashlib.sha256(f"{seed}:{acc}".encode()).digest()[:8], "little")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--genome-table", required=True, help="the simulator's genome table (accession, taxonomy, FASTA[, length])")
    ap.add_argument("--output", required=True, help="the genome table with the strains added")
    ap.add_argument("--out-dir", required=True, help="where the strains' FASTAs and insilico_strains.tsv go")
    ap.add_argument("--positions", help="gene_neighbours.py's gene_positions.tsv for the table's genomes: the marker genes' "
                                        "places, and the real strains' divergence and gene factors (without it every "
                                        "gene is an open reading frame and --ani is needed)")
    ap.add_argument("--taxonomy", help="internal_taxonomy.dmp: which placed genomes are representatives (the others are "
                                       "the real strains); with --positions")
    ap.add_argument("--share", type=float, default=1.0, help="of the species with one genome, the share given a strain (default 1)")
    ap.add_argument("--ani", type=parse_ani, help="MIN-MAX: the genome's ANI drawn uniformly (e.g. 95-99), instead of the "
                                                  "real strains' marker divergence (the default)")
    ap.add_argument("--min-ani", type=float, default=0.95, help="no strain further than this ANI (default 0.95)")
    ap.add_argument("--marker-scale", type=float, default=0.45,
                    help="the median marker gene's divergence over the genome's (default 0.45, from r226 v10's strains)")
    ap.add_argument("--omega", default="auto",
                    help="share of amino-acid changes kept (dN/dS): a number, or auto (the default), the omega at which "
                         f"the strains' third-codon-position share of substitutions equals the real strains' ({OMEGA_DEFAULT} "
                         "without real strains)")
    ap.add_argument("--kappa", type=float, default=3.0, help="transition / transversion rate ratio (default 3)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("-t", "--threads", type=int, default=4)
    args = ap.parse_args(argv)

    rows = read_table(args.genome_table)
    rows = [r for r in rows if not is_insilico(r[0])]
    single = one_genome_species(rows)
    rng = np.random.default_rng(args.seed)
    chosen = [r for r in single if rng.random() < args.share]
    wanted = {normalize_accession(r[0]) for r in chosen}
    genes, strains = {}, {}
    if args.positions:
        reps = representatives(args.taxonomy) if args.taxonomy else None
        genes, strains = read_positions(args.positions, wanted, reps)
    divergences = sorted(strain_divergences(strains))
    factors = conservation_factors(strains)
    if args.ani is None and not divergences:
        sys.exit("No real strains' divergence (--positions with --taxonomy, and other genomes than the representatives "
                 "in the table): give --ani MIN-MAX")
    max_genome = 1.0 - args.min_ani
    os.makedirs(args.out_dir, exist_ok=True)
    omega_note = ""
    if args.omega == "auto":
        share, subs, used, copies = real_spectrum(rows, strains, args.positions, reps or set(), rng, args.threads) \
            if args.positions and args.taxonomy else (None, 0, 0, [])
        if share is None:
            omega = OMEGA_DEFAULT
            omega_note = (f"omega {omega} (no real strains to calibrate on: {subs} substitutions of {used} strains)")
        else:
            omega = calibrate_omega(share, copies, args.kappa, rng)
            omega_note = (f"omega {omega:.3f}: {share:.3f} of the real strains' substitutions on third codon positions "
                          f"({subs} substitutions in the marker genes of {used} strains)")
    else:
        omega = float(args.omega)
        omega_note = f"omega {omega} (given)"
    jobs = []
    for r in chosen:
        acc = normalize_accession(r[0])
        if args.ani is None:
            marker = divergences[int(rng.integers(len(divergences)))]
            genome = min(max(marker / args.marker_scale, 0.0005), max_genome)
            marker = genome * args.marker_scale if genome == max_genome else marker
        else:
            genome = 1.0 - rng.uniform(*args.ani)
            marker = genome * args.marker_scale
        out = os.path.join(args.out_dir, strain_name(acc) + ".fna.gz")
        jobs.append((acc, r[1].split(";")[-1], r[2], out, genes.get(acc, []), factors, genome, marker,
                     seed_of(args.seed, acc), args.kappa, omega))
    with concurrent.futures.ProcessPoolExecutor(max(1, args.threads)) as pool:
        summary = list(pool.map(make_strain, jobs, chunksize=4))
    by_name = {s[0]: s for s in summary}
    with open(os.path.join(args.out_dir, SUMMARY), "w") as fh:
        fh.write("strain\trepresentative\tspecies\tgenome_divergence\tmarker_divergence\tsubstitutions\tcoding_share\t"
                 "marker_divergence_reached\tlength\n")
        fh.writelines("\t".join(s) + "\n" for s in summary)
    with open(args.output + ".partial", "w") as fh:
        fh.writelines("\t".join(r) + "\n" for r in rows)
        for (acc, _sp, _p, out, *_), r in zip(jobs, chosen):
            fh.write("\t".join([strain_name(acc), r[1], os.path.abspath(out), by_name[strain_name(acc)][-1]]) + "\n")
    os.replace(args.output + ".partial", args.output)
    placed = sum(1 for j in jobs if j[4])
    genome_divs = [float(s[3]) for s in summary]
    print(f"{len(summary)} in-silico strains of the {len(single)} species with one genome (of {len(rows)} genomes), "
          f"{placed} with their marker genes placed; genome divergence median "
          f"{statistics.median(genome_divs) if genome_divs else 0:.4f} (ANI {100 * (1 - max(genome_divs, default=0)):.1f}-"
          f"{100 * (1 - min(genome_divs, default=0)):.1f}%); "
          + (f"drawn from {len(divergences)} real strains' marker divergence (median {statistics.median(divergences):.4f}), "
             if args.ani is None else f"ANI drawn from {args.ani[0]:.3f}-{args.ani[1]:.3f}, ")
          + f"{len(factors)} gene factors ({min(factors.values(), default=1):.2f}-{max(factors.values(), default=1):.2f}); "
          f"{omega_note}: {args.output}, {os.path.join(args.out_dir, SUMMARY)}", flush=True)


if __name__ == "__main__":
    main()
