#!/usr/bin/env python3
"""hifi_reads.py - PacBio HiFi reads made from templates: the PacBio reads of the training data and benchmarks.

pbsim3 does not simulate HiFi reads: it simulates the subreads of multi-pass sequencing (--pass-num), which PacBio's
ccs turns into HiFi reads; in one pass with its error model it writes every base quality as '!' (Q0). Its reads at
99.9% accuracy were HiFi-like in their errors only, and protal's quality features (excess_*) took every base of
them for an error. Badread has error and quality models trained on HiFi reads (pacbio2021), but runs in Python
with alignments along each read, which its README calls slow. So the collector makes its PacBio reads here, one
read of each template it drew
(collect_training_data.long_read_templates), with numpy:

- each read's quality R (Phred) falls with its length (QUALITY_BY_LENGTH): Q50 up to 5 kb, Q30 at 25 kb, linearly
  between, Q20 at 50 kb and lower still beyond, at that slope; a read's R is that of its length plus a normal
  deviation of the setup's SD (3 by default), kept between Q_MIN and Q_MAX: its bases' mean error probability is
  10^(-R/10);
- a base's error probability is the read's times the base's weight, scaled so that the read's bases average the
  read's: a base of a homopolymer of n >= 3 weighs (n/2)^2 (HiFi's errors are mostly indels in homopolymers), and
  every base has a lognormal factor of its own (CONFIDENCE_SIGMA), so that qualities vary along a read;
- a base is an error with that probability: in a homopolymer, a base more or less of the run; elsewhere a
  substitution, an insertion or a deletion (ERROR_MIX);
- a base's quality is its error probability (Phred, 1 to 93, ccs's limit), an inserted base's that of the base it
  follows: the qualities say how likely each base is wrong, as calibrated HiFi qualities do, so that the
  differences of a read from its genome are, on average, what its qualities expect.

    python3 scripts/hifi_reads.py --templates templates.fa --out reads.fq.gz [--q_sd 3 --seed 1]
"""

import argparse
import gzip

import numpy as np

# Identifies this model in a design point's key (collect_training_data.py): points made by another are simulated again.
MODEL = ("hifi_reads.py v2: read Q by length (Q50 to 5 kb, Q30 at 25 kb, Q20 at 50 kb) + normal SD, base weights "
         "(n/2)^2 in homopolymers of 3+ and lognormal 1, mix 30:35:35")
# A read's mean quality (Phred) by its length: (length, Q) points, linear between them, flat before the first,
# and on at the last two's slope after the last.
QUALITY_BY_LENGTH = ((5_000, 50.0), (25_000, 30.0), (50_000, 20.0))
Q_MIN, Q_MAX = 10.0, 60.0
P_MAX = 0.5                     # no base is more likely wrong than this
CONFIDENCE_SIGMA = 1.0          # the lognormal spread of the bases' error probabilities around their weight
HOMOPOLYMER = 3                 # runs of this many bases or more are homopolymers
ERROR_MIX = (0.30, 0.35, 0.35)  # substitution, insertion, deletion, outside homopolymers
CHUNK_BASES = 4_000_000         # templates simulated together

_BASES = np.frombuffer(b"ACGT", dtype=np.uint8)
_CODE = np.full(256, 4, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _CODE[_b] = _i
    _CODE[ord(chr(_b).lower())] = _i


def length_quality(lengths):
    """The mean read quality (Phred) of reads of these lengths (QUALITY_BY_LENGTH)."""
    x = np.array([p[0] for p in QUALITY_BY_LENGTH], dtype=np.float64)
    y = np.array([p[1] for p in QUALITY_BY_LENGTH], dtype=np.float64)
    lengths = np.asarray(lengths, dtype=np.float64)
    q = np.interp(lengths, x, y)
    beyond = lengths > x[-1]
    q[beyond] = y[-1] + (lengths[beyond] - x[-1]) * (y[-1] - y[-2]) / (x[-1] - x[-2])
    return q


def mutate(seqs, rng, q_sd):
    """HiFi reads of templates seqs [bytes] (each of at least one base). -> ([(read, qualities as Phred+33)],
    {"events": errors per read, "expected": the expected errors of each read by its qualities, "q": each read's
    drawn quality, "homopolymer_substitutions": substitutions in homopolymers (none)})."""
    lengths = np.fromiter((len(s) for s in seqs), dtype=np.int64, count=len(seqs))
    a = np.frombuffer(b"".join(seqs), dtype=np.uint8)
    n = a.size
    starts = np.concatenate(([0], np.cumsum(lengths)[:-1]))
    read_of = np.repeat(np.arange(len(seqs)), lengths)

    # Homopolymer runs, none across two reads.
    first = np.ones(n, dtype=bool)
    first[1:] = a[1:] != a[:-1]
    first[starts] = True
    run_start = np.flatnonzero(first)
    run_length = np.diff(np.append(run_start, n))
    run = np.repeat(run_length, run_length)
    homopolymer = run >= HOMOPOLYMER

    weight = np.where(homopolymer, (run / 2.0) ** 2, 1.0) * rng.lognormal(0.0, CONFIDENCE_SIGMA, n)
    q_read = np.clip(length_quality(lengths) + rng.normal(0.0, q_sd, len(seqs)), Q_MIN, Q_MAX)
    scale = 10.0 ** (-q_read / 10.0) * lengths / np.add.reduceat(weight, starts)
    p = np.minimum(weight * scale[read_of], P_MAX)

    error = rng.random(n) < p
    kind = rng.random(n)
    sub_share, ins_share = ERROR_MIX[0], ERROR_MIX[1]
    acgt = _CODE[a] < 4
    substitution = error & ~homopolymer & (kind < sub_share) & acgt
    insertion = error & ((homopolymer & (kind < 0.5)) | (~homopolymer & (kind >= sub_share) & (kind < sub_share + ins_share)))
    deletion = error & ~substitution & ~insertion

    base = a.copy()
    base[substitution] = _BASES[(_CODE[a[substitution]] + rng.integers(1, 4, int(substitution.sum()))) % 4]
    inserted = np.where(homopolymer, a, _BASES[rng.integers(0, 4, n)])
    quality = (33 + np.clip(np.rint(-10.0 * np.log10(p)), 1, 93)).astype(np.uint8)

    emitted = (~deletion).astype(np.int64) + insertion
    at = np.repeat(np.arange(n), emitted)
    out_base, out_quality = base[at], quality[at]
    out_base[np.cumsum(emitted)[insertion] - 1] = inserted[insertion]  # an inserted base follows its base
    out_lengths = np.add.reduceat(emitted, starts)
    out_starts = np.concatenate(([0], np.cumsum(out_lengths)[:-1]))

    reads = [(out_base[s:s + k].tobytes(), out_quality[s:s + k].tobytes()) for s, k in zip(out_starts, out_lengths)]
    errors = (10.0 ** (-(out_quality.astype(np.float64) - 33) / 10.0))
    events = substitution | insertion | deletion
    stats = {"events": np.add.reduceat(events.astype(np.int64), starts),
             "expected": np.add.reduceat(errors, out_starts) if errors.size else np.zeros(len(seqs)),
             "q": q_read, "homopolymer_substitutions": int((substitution & homopolymer).sum())}
    return reads, stats


def read_fasta(path):
    """(name, sequence) of a FASTA, in file order."""
    name, parts = None, []
    with (gzip.open(path, "rb") if path.endswith(".gz") else open(path, "rb")) as fh:
        for line in fh:
            line = line.rstrip(b"\r\n")
            if line.startswith(b">"):
                if name is not None:
                    yield name, b"".join(parts)
                name, parts = line[1:].split()[0].decode() if len(line) > 1 else "", []
            elif line:
                parts.append(line.upper())
    if name is not None:
        yield name, b"".join(parts)


def simulate(templates, out, q_sd=3.0, seed=1):
    """A HiFi read of each template of the FASTA `templates`, named after it, into the gzipped FASTQ `out`, in
    file order. -> the number of reads."""
    rng = np.random.default_rng(seed)
    reads = 0
    with gzip.open(out, "wb", compresslevel=1) as fh:
        chunk, bases = [], 0

        def flush():
            made, _ = mutate([s for _, s in chunk], rng, q_sd)
            fh.write(b"".join(b"@" + name.encode() + b"\n" + seq + b"\n+\n" + qual + b"\n"
                              for (name, _), (seq, qual) in zip(chunk, made)))
            return len(made)

        for name, seq in read_fasta(templates):
            if not seq:
                continue
            chunk.append((name, seq))
            bases += len(seq)
            if bases >= CHUNK_BASES:
                reads += flush()
                chunk, bases = [], 0
        if chunk:
            reads += flush()
    return reads


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--templates", required=True, help="FASTA: one read of each sequence")
    ap.add_argument("--out", required=True, help="gzipped FASTQ")
    ap.add_argument("--q_sd", type=float, default=3.0, help="the SD of reads' qualities around their length's")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    print(f"{simulate(args.templates, args.out, args.q_sd, args.seed)} reads written to {args.out}")


if __name__ == "__main__":
    main()
