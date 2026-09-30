#!/usr/bin/env python3
"""qcmsa.py - QC and filter a strain MSA using protal .meta.tsv quality metrics.

This is the Python port of qcmsa.R (milestone M5). It is dependency-free for the
core filtering (Python >= 3.8 standard library only); matplotlib is imported
lazily and only when --plot is requested.

The filter has two passes:

  1. MRate2 (multi-allelicity) filter  -- identical algorithm to qcmsa.R.
     Iterative Tukey-IQR outlier detection on the *count* of multi-allelic
     (MRate2 > 0) cells per gene and per sample. Because MRate2 is sparse
     (mostly 0) a plain fence on the rates collapses to 0 and fires on any
     non-zero value, so we fence the non-zero counts instead and additionally
     require at least --min-bad bad peers before anything is removed.

  2. Site / sequence cleanup (M5 step 4b) -- second pass over the surviving
     columns: drop columns without any A/C/G/T, optionally drop constant sites
     (--discard-constant) and low-parsimony variable sites
     (--min-parsimony-samples), both off by default because a tree's branch
     lengths need them; optionally mask individual cell outliers, and
     re-apply the per-sequence horizontal-coverage floor after gene removal.

Inputs match protal's output contract:
  <msa>        FASTA (plain or .gz) -- protal's <species>.raw.msa.fna
  <partition>  RAxML-style partition -- "DNA, gene<ID> = <start>-<end>" (1-based inclusive;
               the 0-based files of older protal versions are recognised too)
  <meta.tsv>   protal per-sample x per-gene metrics, WITH a header row

Usage:
  qcmsa.py <msa> <partition> <meta.tsv> [options]
"""

import argparse
import gzip
import math
import os
import re
import sys
from collections import Counter, defaultdict

# Characters treated as "missing" (no informative base) at an MSA position.
MISSING = {"-", "N", "n", "."}
# The unambiguous bases; the site cleanup judges columns on these alone.
BASES = set("ACGTacgt")
# What qcmsa writes next to <prefix>; removed first, so no output of an earlier run survives.
OUTPUT_EXTENSIONS = (".msa.fna", ".partition.txt", ".qcmsa_summary.tsv", ".qc.png")

# --preset -> (iqr_mult, min_bad). Tunes how aggressive the MRate2 fence is.
PRESETS = {
    "strict":    (1.0, 1),
    "default":   (1.5, 2),
    "sensitive": (2.0, 3),
}


# ----------------------------------------------------------------------------
# Small numeric helpers (R quantile type 7, matching qcmsa.R's quantile()).
# ----------------------------------------------------------------------------
def quantile_type7(sorted_vals, p):
    """Linear-interpolation quantile (R type 7 / numpy default) on a sorted list."""
    n = len(sorted_vals)
    if n == 0:
        return float("nan")
    if n == 1:
        return float(sorted_vals[0])
    h = (n - 1) * p
    lo = int(math.floor(h))
    hi = min(lo + 1, n - 1)
    frac = h - lo
    return sorted_vals[lo] + frac * (sorted_vals[hi] - sorted_vals[lo])


def upper_fence(values, iqr_mult):
    """Tukey upper fence: Q3 + iqr_mult * (Q3 - Q1)."""
    s = sorted(values)
    q25 = quantile_type7(s, 0.25)
    q75 = quantile_type7(s, 0.75)
    return q75 + iqr_mult * (q75 - q25)


# ----------------------------------------------------------------------------
# I/O
# ----------------------------------------------------------------------------
def open_maybe_gz(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


# A RAxML-style partition names its blocks with an arbitrary string, so match the
# name rather than assuming protal's `gene<N>`. Any tool that writes a valid
# partition file should be usable here; requiring protal's naming made this
# protal-only in practice while documenting itself as generic.
_PART_RE = re.compile(r"(\S+)\s*=\s*(\d+)\s*-\s*(\d+)")

_GENE_NUM_RE = re.compile(r"gene(\d+)\Z")


def canon_gene(name):
    """Join key shared by the partition and the meta.

    protal names a partition `gene4` and calls it `4` in the meta, so the two are
    reconciled by stripping the prefix. Every other name is its own key, which is
    what lets a partition written against real reference record names join against
    a meta written the same way.
    """
    m = _GENE_NUM_RE.fullmatch(str(name))
    return m.group(1) if m else str(name)


def gene_sort_key(g):
    """Numeric ordering for numeric gene ids, lexicographic for the rest.

    Sorting the keys as plain strings would put gene 10 before gene 9 and silently
    reorder every report protal already produces.
    """
    g = str(g)
    return (0, int(g), "") if g.isdigit() else (1, 0, g)


def parse_partition(path, base="auto"):
    """Return (genes, display, base) where genes is [(key:str, start, end)].

    Coordinates come back 0-based inclusive whatever the input used. `display` maps
    each key back to the name as written, so the output partition is spelled the way
    the input was.

    protal and RAxML -- and so most everything else, rg-msa included -- write 1-based
    inclusive; older protal versions wrote 0-based. They are told apart by the lowest start:
    a 1-based file cannot contain 0, and both tools' partitions begin at the start
    of the alignment.
    """
    raw = []
    display = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            m = _PART_RE.search(line)
            if not m:
                continue
            name = m.group(1)
            key = canon_gene(name)
            display.setdefault(key, name)
            raw.append((key, int(m.group(2)), int(m.group(3))))
    if not raw:
        return [], {}, base

    if base == "auto":
        lo = min(s for _, s, _ in raw)
        if lo == 0:
            detected = 0
        elif lo == 1:
            detected = 1
        else:
            sys.stderr.write(
                f"qcmsa.py: partition starts at {lo}, which is neither 0- nor "
                f"1-based at the alignment start; assuming 1-based (RAxML). Pass "
                f"--partition-base to be explicit.\n"
            )
            detected = 1
    else:
        detected = int(base)

    shift = 1 if detected == 1 else 0
    genes = [(g, s - shift, e - shift) for g, s, e in raw]
    return genes, display, detected


# protal .meta.tsv columns (the header protal writes). We read by header name so
# the extra coverage columns added over time don't shift our indices, and -- the
# bug qcmsa.R had -- the header line is NOT mistaken for data.
META_GENE_COL = "gene_id"
META_SAMPLE_COL = "sample"
META_MRATE2_COL = "multi_rate_vcov2"
META_MULTI_COL = "multi_allelic"      # positions written as an IUPAC code
META_VCOV2_COL = "counts_vcov2"       # positions with >= 2 reads
META_HCOV_COL = "hcov"               # fraction of gene covered (M3 --gene_min_hcov_frac)
META_DEPTH_COL = "mean_vcov_nonzero"  # mean depth over covered positions (M3 --gene_min_mean_depth)


def load_meta(path, gene_whitelist, sample_whitelist=None):
    """Load per-(sample,gene) meta restricted to genes in gene_whitelist and, if given, to
    samples in sample_whitelist (the samples in the MSA: protal lists every sample with reads
    on a gene in the meta, also those whose MSA row it dropped).

    Returns (rows, samples_in_order, genes_sorted, cov) where
      rows = list of (sample:str, gene:int, mrate2:float, multi:int, vcov2:int): the cell's
             multi-allelic rate, its multi-allelic positions and its positions with >= 2 reads
      cov  = {(sample, gene): (hcov:float, mean_depth:float)} (empty if the
             coverage columns are absent).
    """
    rows = []
    cov = {}
    with open_maybe_gz(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            si = header.index(META_SAMPLE_COL)
            gi = header.index(META_GENE_COL)
            mi = header.index(META_MRATE2_COL)
            ni = header.index(META_MULTI_COL)
            vi = header.index(META_VCOV2_COL)
        except ValueError as exc:
            raise SystemExit(
                f"qcmsa.py: meta file '{path}' is missing expected column "
                f"({exc}); header was: {header}"
            )
        hi = header.index(META_HCOV_COL) if META_HCOV_COL in header else None
        di = header.index(META_DEPTH_COL) if META_DEPTH_COL in header else None
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            gene = canon_gene(f[gi])
            if gene not in gene_whitelist:
                continue
            sample = f[si]
            if sample_whitelist is not None and sample not in sample_whitelist:
                continue
            rows.append((sample, gene, float(f[mi]), int(f[ni]), int(f[vi])))
            if hi is not None and di is not None:
                try:
                    cov[(sample, gene)] = (float(f[hi]), float(f[di]))
                except ValueError:
                    pass

    samples_seen = []
    seen = set()
    genes_seen = set()
    for sample, gene, *_ in rows:
        if sample not in seen:
            seen.add(sample)
            samples_seen.append(sample)
        genes_seen.add(gene)
    return rows, samples_seen, sorted(genes_seen, key=gene_sort_key), cov


def coverage_filter(cov, all_genes, hcov_t, depth_t, min_samples):
    """Reproduce protal's M3 coverage gate from the meta coverage columns.

    A (sample,gene) cell passes if hcov >= hcov_t AND mean_depth >= depth_t.
    A gene is dropped if NOT more than min_samples cells pass (protal uses a
    strict '>' on msa_min_samples). Returns (dropped_genes, fail_cells, reason)
    where fail_cells are coverage-failing cells in *surviving* genes (to gap-fill).
    """
    passing = defaultdict(int)
    failing = defaultdict(set)   # gene -> set of failing samples
    for (sample, gene), (h, d) in cov.items():
        if gene not in all_genes:
            continue
        if h >= hcov_t and d >= depth_t:
            passing[gene] += 1
        else:
            failing[gene].add(sample)
    dropped, fail_cells, reason = set(), set(), {}
    for gene in all_genes:
        np = passing.get(gene, 0)
        if np <= min_samples:        # not strictly greater -> dropped (matches protal)
            dropped.add(gene)
            reason[gene] = (np, "M3 coverage: only %d sample(s) pass" % np)
        else:
            for s in failing.get(gene, ()):
                fail_cells.add((s, gene))
    return dropped, fail_cells, reason


def read_fasta(path):
    """Return (names:list, seqs:list) preserving file order."""
    names, seqs = [], []
    cur = []
    with open_maybe_gz(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line[0] == ">":
                if names:
                    seqs.append("".join(cur))
                names.append(line[1:].strip())
                cur = []
            else:
                cur.append(line)
    if names:
        seqs.append("".join(cur))
    return names, seqs


def write_fasta(path, names, seqs, width=80):
    with open(path, "w") as fh:
        for name, seq in zip(names, seqs):
            fh.write(">" + name + "\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


# ----------------------------------------------------------------------------
# Pass 1: multi-allelicity (MRate2) sample/gene filter.
# ----------------------------------------------------------------------------
def mrate2_filter(rows, all_genes, all_samples, min_bad, iqr_mult,
                  gene_abs=0, sample_abs=0, min_rate=0.004):
    """Remove samples, then genes, whose pooled multi-allelic rate is an outlier.

    A sample's pooled rate is its multi-allelic positions (IUPAC codes) over its positions with
    >= 2 reads, summed over genes; a gene's, the same over the samples kept. Unlike a count of
    multi-allelic genes, the rate does not grow with depth, so deep clean samples are not taken
    for mixtures. An item is removed when its rate exceeds both the Tukey upper fence of all
    items' rates (Q3 + iqr_mult * IQR, over >= 4 items) and min_rate, and it has >= min_bad
    multi-allelic genes (samples); or, with sample_abs / gene_abs > 0, when it has at least that
    many. One pass: removals do not move the fence for the rest.

    Returns (kept_genes, kept_samples, filtered_genes, filtered_samples, gene_reason,
    sample_reason), the reasons as {id: (value, text)}.
    """
    def rates(key_index, keep):
        multi, positions, bad = defaultdict(int), defaultdict(int), defaultdict(int)
        for row in rows:
            sample, gene, _, m, v = row
            if not keep(sample, gene):
                continue
            key = row[key_index]
            multi[key] += m
            positions[key] += v
            bad[key] += m > 0
        return multi, positions, bad

    def flag(keys, multi, positions, bad, abs_min, what, peers):
        rate = {k: (multi[k] / positions[k] if positions.get(k) else 0.0) for k in keys}
        fence = upper_fence(list(rate.values()), iqr_mult) if len(rate) >= 4 else None
        threshold = max(fence, min_rate) if fence is not None else None
        flagged, reason = set(), {}
        for k in keys:
            if threshold is not None and rate[k] > threshold and bad[k] >= min_bad:
                flagged.add(k)
                reason[k] = (bad[k], f"multi-allelic rate {rate[k]:.4f} > {threshold:.4f} "
                                     f"(fence {fence:.4f}, floor {min_rate}); {bad[k]} multi-allelic {peers}")
            elif abs_min > 0 and bad[k] >= abs_min:
                flagged.add(k)
                reason[k] = (bad[k], f"multi-allelic in {bad[k]} {peers} >= {abs_min}")
        if fence is None and abs_min <= 0 and any(bad[k] > 0 for k in keys):
            sys.stderr.write(f"qcmsa.py: WARNING: fewer than 4 {what}s, too few for the {what} filter; "
                             f"set --{what}-abs-min-bad to filter them.\n")
        return flagged, reason, threshold

    samples = list(all_samples)
    multi, positions, bad = rates(0, lambda s, g: True)
    bad_samples, sample_reason, s_thr = flag(samples, multi, positions, bad, sample_abs, "sample", "genes")
    kept_samples = set(samples) - bad_samples

    multi, positions, bad = rates(1, lambda s, g: s in kept_samples)
    bad_genes, gene_reason, g_thr = flag(list(all_genes), multi, positions, bad, gene_abs, "gene", "samples")
    kept_genes = set(all_genes) - bad_genes

    fmt = lambda t: "-" if t is None else f"{t:.4f}"
    sys.stderr.write(f"  multi-allelic filter: {len(bad_samples)} sample(s) above rate {fmt(s_thr)}, "
                     f"{len(bad_genes)} gene(s) above rate {fmt(g_thr)}\n")
    return (kept_genes, kept_samples, bad_genes, bad_samples, gene_reason, sample_reason)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def build_argparser():
    p = argparse.ArgumentParser(
        description="QC and filter a strain MSA using protal .meta.tsv metrics."
    )
    p.add_argument("msa", help="MSA FASTA (plain or .gz)")
    p.add_argument("partition", help="RAxML-style partition file")
    p.add_argument("--partition-base", choices=["auto", "0", "1"], default="auto",
                   help="Coordinate base of the partition file. auto (default) "
                        "detects it from the lowest start: protal, RAxML and rg-msa "
                        "write 1-based (older protal versions wrote 0-based).")
    p.add_argument("meta", help="protal .meta.tsv (with header)")
    p.add_argument("--prefix", default=None,
                   help="Output prefix (default: MSA path with .fna/.gz stripped)")

    p.add_argument("--preset", choices=sorted(PRESETS),
                   help="Convenience bundle for --iqr-mult/--min-bad "
                        "(strict=1.0/1, default=1.5/2, sensitive=2.0/3). "
                        "Explicit --iqr-mult/--min-bad override the preset.")
    p.add_argument("--iqr-mult", type=float, default=None,
                   help="Tukey IQR multiplier for the MRate2 fence (default 1.5)")
    p.add_argument("--min-bad", type=int, default=None,
                   help="Min bad peers before a gene/sample is removed (default 2)")
    # Absolute multi-allelicity cutoffs -- catch ANY signal above a clean-zero
    # baseline (e.g. conspecific/mixed strains), which the Tukey fence cannot do
    # because the bad items then ARE the distribution. Default 0 = off.
    p.add_argument("--sample-abs-min-bad", type=int, default=0,
                   help="Remove a sample with multi-allelic (MRate2>0) signal in >= this "
                        "many genes, regardless of the Tukey fence. 0=off. Try 1-2 to catch "
                        "mixed/conspecific strains.")
    p.add_argument("--gene-abs-min-bad", type=int, default=0,
                   help="Remove a gene multi-allelic in >= this many samples, regardless of "
                        "the Tukey fence. 0=off.")
    p.add_argument("--mrate2-min-rate", type=float, default=0.004,
                   help="Floor for the multi-allelic fences: a sample, gene or cell is only "
                        "removed (masked) with a multi-allelic rate above both its Tukey fence "
                        "and this rate (default 0.004: 0.4%% of the positions with >= 2 reads). "
                        "It keeps a fence from collapsing to 0 when most items are clean.")
    p.add_argument("--mrate2-include-zeros", action="store_true",
                   help="No effect; kept for old command lines. The fences are computed on "
                        "every item's rate, floored at --mrate2-min-rate.")

    # Coverage gating -- this is where the gene/sample coverage filtering lives
    # (protal emits a raw MSA). Computed from the meta hcov / mean-depth columns.
    # Defaults are ON; set any to 0 to disable that part.
    p.add_argument("--gene-min-hcov", type=float, default=0.3,
                   help="Min fraction of a gene covered for a (sample,gene) cell to "
                        "pass. Default 0.3; 0 disables.")
    p.add_argument("--gene-min-mean-depth", type=float, default=1.0,
                   help="Min mean depth over covered positions for a cell to pass. "
                        "Default 1.0; 0 disables.")
    p.add_argument("--gene-min-samples", type=int, default=1,
                   help="Drop a gene unless MORE than this many samples pass coverage "
                        "(strict >, like protal's old msa_min_samples). Default 1: a gene "
                        "needs 2 samples, as protal needs 2 samples for an MSA. 0 disables.")
    p.add_argument("--max-mrate2", type=float, default=None,
                   help="Hard per-cell MRate2 cap for the cell-outlier fence "
                        "(default: Tukey fence derived from the data)")

    # M5 step 4b -- site / sequence cleanup
    # Constant (invariant) sites are KEPT by default (they inform branch-length
    # estimation). Pass --discard-constant to drop them. Both site tests judge the
    # bases A/C/G/T only: an IUPAC code is ambiguous, as IQ-TREE reads it.
    p.add_argument("--discard-constant", dest="remove_constant",
                   action="store_true", default=False,
                   help="Discard constant (invariant) sites: columns with at most one of "
                        "A/C/G/T (default: off). A tree from such an MSA needs an "
                        "ascertainment correction (IQ-TREE: +ASC).")
    p.add_argument("--min-parsimony-samples", type=int, default=0,
                   help="Drop variable sites where fewer than N samples (the reference "
                        "row not counted) differ from the majority base (default 0: keep "
                        "every site). A site where only one sample differs is a strain's "
                        "own mutation: dropping those shortens terminal branches to near "
                        "0, so use this for topology-only analyses.")
    p.add_argument("--reapply-hcov", type=int, default=0,
                   help="After gene/site removal, drop sequences with fewer than "
                        "this many valid (non -/N) bases (default 0 = disabled). "
                        "protal passes its --msa_min_hcov here.")
    p.add_argument("--mask-cell-outliers", dest="mask_cells",
                   action="store_true", default=True,
                   help="Mask individual (sample,gene) MRate2 outlier cells with "
                        "'-' (default: on). qcmsa.R only marks these on the plot.")
    p.add_argument("--no-mask-cell-outliers", dest="mask_cells",
                   action="store_false", help="Do not mask cell outliers (qcmsa.R parity)")

    p.add_argument("--plot", action="store_true",
                   help="Also emit a <prefix>.qc.png MRate2 heatmap (needs matplotlib)")
    p.add_argument("--no-summary", dest="summary", action="store_false", default=True,
                   help="Do not write the <prefix>.qcmsa_summary.tsv decision log")
    return p


def resolve_params(args):
    iqr_mult, min_bad = PRESETS["default"]
    if args.preset:
        iqr_mult, min_bad = PRESETS[args.preset]
    if args.iqr_mult is not None:
        iqr_mult = args.iqr_mult
    if args.min_bad is not None:
        min_bad = args.min_bad
    return iqr_mult, min_bad


def main(argv=None):
    args = build_argparser().parse_args(argv)
    iqr_mult, min_bad = resolve_params(args)

    prefix = args.prefix
    if prefix is None:
        # strip ".raw.msa.fna" / ".msa.fna" / ".fna" so the output is <name>.msa.fna
        prefix = re.sub(r"(\.raw)?(\.msa)?\.fna(\.gz)?$", "", args.msa)
    if os.path.abspath(prefix + ".msa.fna") == os.path.abspath(args.msa):
        raise SystemExit("qcmsa.py: output would overwrite the input MSA; pass a "
                         "distinct --prefix (input should be <name>.raw.msa.fna).")
    # Outputs of an earlier run must not stay behind as if this run had written them.
    for ext in OUTPUT_EXTENSIONS:
        if os.path.exists(prefix + ext):
            os.remove(prefix + ext)

    def stop(reason, counts=()):
        """Write no MSA: say why on stderr and in the summary, and return 0."""
        sys.stderr.write(f"qcmsa.py: {reason} - no MSA produced (skipping).\n")
        if args.summary:
            with open(prefix + ".qcmsa_summary.tsv", "w") as fh:
                fh.write("section\tkey\tvalue\treason\n")
                fh.write(f"status\tno_msa\t\t{reason}\n")
                for key, value in counts:
                    fh.write(f"count\t{key}\t{value}\t\n")
        return 0

    # --- inputs ---
    partition, gene_display, part_base = parse_partition(
        args.partition, args.partition_base
    )
    if not partition:
        raise SystemExit(f"qcmsa.py: no genes parsed from partition '{args.partition}'")
    gene_whitelist = {g for g, _, _ in partition}
    sys.stderr.write(f"Partition: {len(partition)} genes\n")

    # --- read MSA: only its samples count (e.g. towards --gene-min-samples) ---
    names, seqs = read_fasta(args.msa)
    if not names:
        raise SystemExit(f"qcmsa.py: empty MSA '{args.msa}'")
    duplicates = sorted(n for n, c in Counter(names).items() if c > 1)
    if duplicates:
        raise SystemExit(f"qcmsa.py: '{args.msa}' names {len(duplicates)} sequence(s) more than once "
                         f"({', '.join(duplicates[:5])}); give every sample its own #SAMPLEID")
    seq_of = dict(zip(names, seqs))
    sys.stderr.write(f"MSA loaded: {len(names)} sequences, {len(seqs[0])} bp\n")

    rows, all_samples, all_genes, cov = load_meta(args.meta, gene_whitelist, set(names))
    all_samples = sorted(all_samples)
    sys.stderr.write(
        f"Loaded meta: {len(all_samples)} samples x {len(all_genes)} genes (in MSA)\n"
    )
    sys.stderr.write(f"Params: iqr_mult={iqr_mult} min_bad={min_bad}"
                     + (f" preset={args.preset}" if args.preset else "") + "\n")

    # --- pass 0: optional coverage gate (reproduces protal M3 from meta) ---
    cov_gate = (args.gene_min_hcov > 0 or args.gene_min_mean_depth > 0
                or args.gene_min_samples > 0)
    cov_dropped_genes, cov_fail_cells, cov_reason = set(), set(), {}
    if cov_gate:
        if not cov:
            sys.stderr.write("qcmsa.py: coverage gating requested but meta has no "
                             "hcov/mean_vcov_nonzero columns; skipping coverage gate.\n")
        else:
            cov_dropped_genes, cov_fail_cells, cov_reason = coverage_filter(
                cov, set(all_genes), args.gene_min_hcov,
                args.gene_min_mean_depth, args.gene_min_samples)
            sys.stderr.write(
                f"Coverage gate (hcov>={args.gene_min_hcov}, depth>="
                f"{args.gene_min_mean_depth}, >{args.gene_min_samples} samples): "
                f"dropped {len(cov_dropped_genes)}/{len(all_genes)} genes, "
                f"gap-filled {len(cov_fail_cells)} cell(s)\n")
    # Genes/cells removed by coverage don't participate in the MRate2 stats (they
    # are gaps in the output), so the adaptive fence is computed on covered data.
    cov_survivor_genes = [g for g in all_genes if g not in cov_dropped_genes]
    rows_for_mrate2 = [row for row in rows
                       if row[1] not in cov_dropped_genes and (row[0], row[1]) not in cov_fail_cells]

    # --- pass 1: MRate2 gene/sample filter ---
    (kept_genes, kept_samples, mr_filtered_genes, filtered_samples,
     gene_reason, sample_reason) = mrate2_filter(
        rows_for_mrate2, cov_survivor_genes, all_samples, min_bad, iqr_mult,
        gene_abs=args.gene_abs_min_bad, sample_abs=args.sample_abs_min_bad,
        min_rate=args.mrate2_min_rate
    )
    filtered_genes = mr_filtered_genes | cov_dropped_genes
    sys.stderr.write(f"Filtered genes: {len(filtered_genes)} / {len(all_genes)}"
                     f" ({len(cov_dropped_genes)} coverage, {len(mr_filtered_genes)} MRate2)\n")
    sys.stderr.write(f"Filtered samples: {len(filtered_samples)} / {len(all_samples)}\n")

    # --- cell-outlier fence ---
    # The fence covers the rates of the kept cells, zeros included, floored at
    # --mrate2-min-rate as the sample and gene fences are (on mostly clean data Q3 is 0). A cell
    # also needs --min-bad multi-allelic positions, so one IUPAC code never masks a gene.
    kept_rows = [row for row in rows_for_mrate2
                 if row[0] not in filtered_samples and row[1] not in filtered_genes]
    if args.max_mrate2 is not None:
        cell_fence = args.max_mrate2
    elif len(kept_rows) >= 4:
        cell_fence = max(upper_fence([mr for _, _, mr, _, _ in kept_rows], iqr_mult),
                         args.mrate2_min_rate)
    else:
        cell_fence = float("inf")
    outlier_cells = {(sample, gene) for sample, gene, mr, multi, _ in kept_rows
                     if mr > cell_fence and multi >= min_bad}
    sys.stderr.write(
        f"Cell fence (MRate2): {cell_fence:.5f} -> {len(outlier_cells)} outlier cell(s)\n"
    )

    sample_set = set(all_samples)
    # Keep references (names not in meta) always; drop filtered samples.
    kept_names = [n for n in names if n not in sample_set or n not in filtered_samples]

    # Degenerate MSA (e.g. only the reference survived protal's row filter): nothing
    # meaningful to filter. Warn and skip gracefully so batch/protal-driven runs continue.
    progress = [("samples_in", len(all_samples)), ("samples_filtered", len(filtered_samples)),
                ("genes_in", len(all_genes)), ("genes_filtered_coverage", len(cov_dropped_genes)),
                ("genes_filtered_mrate2", len(mr_filtered_genes))]
    n_sample_seqs = sum(1 for n in kept_names if n in sample_set)
    if n_sample_seqs < 2:
        return stop(f"only {n_sample_seqs} sample sequence(s) left after the multi-allelic filter", progress)

    # --- surviving genes -> original column ranges (sorted by start) ---
    kept_partition = sorted(
        [(g, s, e) for (g, s, e) in partition if g not in filtered_genes],
        key=lambda t: t[1],
    )
    if not kept_partition:
        return stop("every gene was filtered (coverage gate or multi-allelic filter)", progress)

    # Per surviving column: which gene it belongs to, and its original index.
    col_gene = []
    col_orig = []
    for g, s, e in kept_partition:
        for c in range(s, e + 1):
            col_gene.append(g)
            col_orig.append(c)

    # Subset each kept sequence to the surviving columns; gap-fill coverage-failed
    # cells (M3-equivalent) and, if requested, MRate2 outlier cells.
    msa_rows = []
    for n in kept_names:
        seq = seq_of[n]
        chars = [seq[c] for c in col_orig]
        if n in sample_set:
            for j, g in enumerate(col_gene):
                if (n, g) in cov_fail_cells or (args.mask_cells and (n, g) in outlier_cells):
                    chars[j] = "-"
        msa_rows.append(chars)

    n_cols = len(col_gene)

    # --- reapply-hcov: drop sequences below the valid-base floor AFTER gene removal ---
    # Measured on the gene-filtered alignment (full gene columns), NOT after site
    # cleanup -- otherwise the floor is compared against the tiny informative-only
    # alignment and would drop every sample.
    if args.reapply_hcov > 0:
        keep_idx = []
        for i, (n, row) in enumerate(zip(kept_names, msa_rows)):
            if n in sample_set:  # references are exempt from the floor
                valid = sum(1 for ch in row if ch not in MISSING)
                if valid < args.reapply_hcov:
                    continue
            keep_idx.append(i)
        dropped = len(kept_names) - len(keep_idx)
        if dropped:
            sys.stderr.write(
                f"reapply-hcov={args.reapply_hcov}: dropped {dropped} sequence(s) "
                "below the valid-base floor\n"
            )
        kept_names = [kept_names[i] for i in keep_idx]
        msa_rows = [msa_rows[i] for i in keep_idx]
        n_sample_seqs = sum(1 for n in kept_names if n in sample_set)
        if n_sample_seqs < 2:
            return stop(f"only {n_sample_seqs} sample sequence(s) have {args.reapply_hcov} valid bases "
                        "(--reapply-hcov)", progress)

    # --- pass 2: site cleanup (all-missing / constant / low-parsimony) ---
    # Sites are judged on A/C/G/T: an IUPAC code is an ambiguity, as IQ-TREE reads it, so a
    # column of A and R is constant and a column of only N, '-' and IUPAC codes holds no base.
    # The parsimony count covers the sample rows only: the reference row is not a sample.
    is_sample_row = [n in sample_set for n in kept_names]
    col_keep = [True] * n_cols
    site_removal_reason = Counter()  # reason -> n_sites, for the summary breakdown
    for j in range(n_cols):
        counts = Counter()
        sample_counts = Counter()
        for row, is_sample in zip(msa_rows, is_sample_row):
            ch = row[j]
            if ch in BASES:
                ch = ch.upper()
                counts[ch] += 1
                if is_sample:
                    sample_counts[ch] += 1
        if not counts:
            col_keep[j] = False
            site_removal_reason["all_missing"] += 1
        elif args.remove_constant and len(counts) <= 1:
            col_keep[j] = False
            site_removal_reason["constant"] += 1
        elif args.min_parsimony_samples > 0 and len(sample_counts) > 1:
            minor = sum(sample_counts.values()) - max(sample_counts.values())
            if minor < args.min_parsimony_samples:
                col_keep[j] = False
                site_removal_reason["low_parsimony"] += 1

    surviving = [j for j in range(n_cols) if col_keep[j]]
    n_removed_sites = n_cols - len(surviving)
    if not surviving:
        return stop("every site was removed by the site cleanup", progress + [("sites_in", n_cols)])

    # Final sequences (column subset).
    final_seqs = ["".join(row[j] for j in surviving) for row in msa_rows]

    sys.stderr.write(
        f"Site cleanup: removed {n_removed_sites} / {n_cols} sites "
        f"({len(surviving)} retained) -- "
        f"all_missing={site_removal_reason['all_missing']}, "
        f"constant={site_removal_reason['constant']}, "
        f"low_parsimony={site_removal_reason['low_parsimony']}\n"
    )

    # --- write filtered MSA ---
    msa_out = prefix + ".msa.fna"
    write_fasta(msa_out, kept_names, final_seqs)
    sys.stderr.write(f"Saved: {msa_out}\n")

    # --- write updated partition (recomputed contiguous coordinates) ---
    # Surviving column count per gene, in kept_partition order.
    surv_per_gene = defaultdict(int)
    for j in surviving:
        surv_per_gene[col_gene[j]] += 1
    part_out = prefix + ".partition.txt"
    with open(part_out, "w") as fh:
        new_start = 0
        for g, _, _ in kept_partition:
            length = surv_per_gene.get(g, 0)
            if length == 0:
                continue  # gene lost all its sites in cleanup
            new_end = new_start + length - 1
            # Spelled and based as the input was, so the file round-trips.
            out_shift = 1 if part_base == 1 else 0
            fh.write(
                f"DNA, {gene_display.get(g, g)} = "
                f"{new_start + out_shift}-{new_end + out_shift}\n"
            )
            new_start = new_end + 1
    sys.stderr.write(f"Saved: {part_out}\n")

    # --- machine-readable decision log (drives the strain_report plots) ---
    if args.summary:
        n_seqs_out = len(kept_names)
        n_ref = sum(1 for n in kept_names if n not in sample_set)
        summary_out = prefix + ".qcmsa_summary.tsv"
        with open(summary_out, "w") as fh:
            fh.write("section\tkey\tvalue\treason\n")
            fh.write("status\tmsa\t\t\n")
            fh.write(f"param\tiqr_mult\t{iqr_mult}\t\n")
            fh.write(f"param\tmin_bad\t{min_bad}\t\n")
            fh.write(f"param\tpreset\t{args.preset or ''}\t\n")
            fh.write(f"param\tmrate2_min_rate\t{args.mrate2_min_rate}\t\n")
            fh.write(f"param\tcell_fence\t{cell_fence:.6g}\t\n")
            fh.write(f"count\tsamples_in\t{len(all_samples)}\t\n")
            fh.write(f"count\tsamples_kept\t{len(kept_samples)}\t\n")
            fh.write(f"count\tsamples_filtered\t{len(filtered_samples)}\t\n")
            fh.write(f"count\tgenes_in\t{len(all_genes)}\t\n")
            fh.write(f"count\tgenes_kept\t{len(kept_genes)}\t\n")
            fh.write(f"count\tgenes_filtered\t{len(filtered_genes)}\t\n")
            fh.write(f"count\tgenes_filtered_coverage\t{len(cov_dropped_genes)}\t\n")
            fh.write(f"count\tgenes_filtered_mrate2\t{len(mr_filtered_genes)}\t\n")
            fh.write(f"count\tcoverage_gap_filled_cells\t{len(cov_fail_cells)}\t\n")
            fh.write(f"count\tsites_in\t{n_cols}\t\n")
            fh.write(f"count\tsites_kept\t{len(surviving)}\t\n")
            fh.write(f"count\tsites_removed\t{n_removed_sites}\t\n")
            fh.write(f"count\tsites_removed_all_missing\t{site_removal_reason['all_missing']}"
                     f"\tno A/C/G/T in any kept sequence (only -, N or IUPAC codes)\n")
            fh.write(f"count\tsites_removed_constant\t{site_removal_reason['constant']}"
                     f"\t--discard-constant: at most one of A/C/G/T in the column\n")
            fh.write(f"count\tsites_removed_low_parsimony\t{site_removal_reason['low_parsimony']}"
                     f"\t--min-parsimony-samples {args.min_parsimony_samples}: "
                     f"fewer than that many samples differ from the majority base\n")
            fh.write(f"count\toutlier_cells\t{len(outlier_cells)}\t\n")
            fh.write(f"count\tseqs_out\t{n_seqs_out}\t\n")
            fh.write(f"count\treference_seqs_out\t{n_ref}\t\n")
            for g in sorted(filtered_genes, key=gene_sort_key):
                if g in cov_dropped_genes:
                    nb, txt = cov_reason.get(g, (None, "M3 coverage"))
                    fh.write(f"gene_filtered\t{g}\t{nb}\t{txt}\n")
                    continue
                nb, reason = gene_reason.get(g, (None, "multi-allelic outlier"))
                fh.write(f"gene_filtered\t{g}\t{nb}\t{reason}\n")
            for s in sorted(filtered_samples):
                nb, reason = sample_reason.get(s, (None, "multi-allelic outlier"))
                fh.write(f"sample_filtered\t{s}\t{nb}\t{reason}\n")
            for s, g in sorted(outlier_cells):
                fh.write(f"cell_outlier\t{s}|gene{g}\t\tMRate2 > cell_fence {cell_fence:.3g}\n")
        sys.stderr.write(f"Saved: {summary_out}\n")

    if args.plot:
        try:
            make_plot(prefix, rows, all_genes, all_samples,
                      filtered_genes, filtered_samples, outlier_cells)
        except Exception as exc:  # plotting must never break the filter
            sys.stderr.write(f"qcmsa.py: --plot failed ({exc}); skipping heatmap\n")

    sys.stderr.write("Done.\n")
    return 0


def make_plot(prefix, rows, all_genes, all_samples,
              filtered_genes, filtered_samples, outlier_cells):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    gidx = {g: i for i, g in enumerate(all_genes)}
    sidx = {s: i for i, s in enumerate(all_samples)}
    grid = [[float("nan")] * len(all_genes) for _ in all_samples]
    for sample, gene, mr, *_ in rows:
        if sample in sidx and gene in gidx:
            grid[sidx[sample]][gidx[gene]] = mr

    fig, ax = plt.subplots(figsize=(max(6, len(all_genes) * 0.06 + 1.5),
                                    max(3, len(all_samples) * 0.06 + 1.5)))
    im = ax.imshow(grid, aspect="auto", cmap="Blues", interpolation="nearest")
    ax.set_xlabel("Gene")
    ax.set_ylabel("Sample")
    ax.set_title("MRate2 (red label = filtered, x = outlier cell)")
    ax.set_yticks(range(len(all_samples)))
    ax.set_yticklabels(all_samples, fontsize=4)
    for tick, s in zip(ax.get_yticklabels(), all_samples):
        if s in filtered_samples:
            tick.set_color("red")
    ax.set_xticks([])
    for sample, gene in outlier_cells:
        if sample in sidx and gene in gidx:
            ax.plot(gidx[gene], sidx[sample], marker="x", color="red", markersize=2)
    fig.colorbar(im, ax=ax, shrink=0.5, label="MRate2")
    out = prefix + ".qc.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    sys.stderr.write(f"Saved: {out}\n")


if __name__ == "__main__":
    sys.exit(main())
