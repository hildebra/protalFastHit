#!/usr/bin/env python3
"""check_results.py - compare a protal run on simulated reads with the truth.

Checks, per species of the mock community:
  detection   every species in the truth is in the profile, nothing else is
  abundance   |profile - truth| <= --max_abundance_error (relative abundance)
  assignment  of the reads protal aligned (primary records), the fraction
              whose reference gene belongs to the read's source species
              >= --min_assignment_accuracy
  alignment   aligned reads / simulated reads >= --min_aligned_fraction (a
              sanity floor: only the marker genes, ~45% of a genome, are indexed)

Prints a table, writes it to --out (TSV) and exits 1 if any check fails.

Usage:
  check_results.py --truth reads.truth.tsv --profile mini.profile --sam mini.sam
                   --taxonomy protal_db/internal_taxonomy.dmp [--out check.tsv]
"""

import argparse
import gzip
import os
import sys
from collections import Counter, defaultdict


def read_truth(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, l.rstrip("\n").split("\t"))) for l in fh if l.strip()]


def read_taxonomy(path):
    """internal taxid -> scientific name."""
    with open(path) as fh:
        next(fh)
        return {int(f[0]): f[3] for f in (l.rstrip("\n").split("\t") for l in fh)}


def read_profile(path):
    """species -> abundance from protal's .profile (rep_genome, lineage, abundance)."""
    abundance = {}
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            lineage = next((f for f in fields if f.startswith("d__")), None)
            if lineage is None:
                continue
            abundance[lineage.split(";")[-1]] = float(fields[-1])
    return abundance


def read_sam(path, taxid_name):
    """-> {source accession: Counter(assigned species)} over primary, mapped records."""
    opener = gzip.open if path.endswith(".gz") else open
    assigned = defaultdict(Counter)
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            qname, flag, rname = line.split("\t", 3)[:3]
            if int(flag) & (0x4 | 0x100 | 0x800):
                continue
            source = qname.rsplit("/", 1)[0].rsplit("-", 1)[0]
            species = taxid_name.get(int(rname.split("_", 1)[0]), f"taxid {rname}")
            assigned[source][species] += 1
    return assigned


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--truth", required=True, help="<reads>.truth.tsv from simulate_reads.py")
    ap.add_argument("--profile", required=True, help="protal .profile")
    ap.add_argument("--sam", required=True, help="protal .sam (or .sam.gz)")
    ap.add_argument("--taxonomy", required=True, help="internal_taxonomy.dmp of the DB")
    ap.add_argument("--out", help="write the result table here")
    ap.add_argument("--max_abundance_error", type=float, default=0.05)
    ap.add_argument("--min_assignment_accuracy", type=float, default=0.95)
    ap.add_argument("--min_aligned_fraction", type=float, default=0.2)
    args = ap.parse_args()

    sam = args.sam if os.path.exists(args.sam) or not os.path.exists(args.sam + ".gz") else args.sam + ".gz"
    truth = read_truth(args.truth)
    profile = read_profile(args.profile)
    assigned = read_sam(sam, read_taxonomy(args.taxonomy))

    # Species-level truth (several genomes of one species are summed).
    expected, reads, aligned, correct = Counter(), Counter(), Counter(), Counter()
    for row in truth:
        sp = row["species"]
        expected[sp] += float(row["relative_abundance"])
        reads[sp] += 2 * int(row["read_pairs"])
        counts = assigned.get(row["accession"], Counter())
        aligned[sp] += sum(counts.values())
        correct[sp] += counts[sp]

    rows, failures = [], []
    for sp in sorted(set(expected) | set(profile)):
        truth_ab, prof_ab = expected.get(sp, 0.0), profile.get(sp)
        problems = []
        if sp not in expected:
            problems.append("false positive")
        elif prof_ab is None:
            problems.append("not detected")
        elif abs(prof_ab - truth_ab) > args.max_abundance_error:
            problems.append(f"abundance off by {abs(prof_ab - truth_ab):.3f}")
        aln_frac = aligned[sp] / reads[sp] if reads[sp] else None
        accuracy = correct[sp] / aligned[sp] if aligned[sp] else None
        if sp in expected:
            if aln_frac is None or aln_frac < args.min_aligned_fraction:
                problems.append("too few reads aligned")
            if accuracy is not None and accuracy < args.min_assignment_accuracy:
                problems.append("reads assigned to the wrong species")
        fmt = lambda v, p=3: "-" if v is None else f"{v:.{p}f}"
        rows.append([sp, fmt(truth_ab), fmt(prof_ab), fmt(aln_frac), fmt(accuracy, 4),
                     "; ".join(problems) or "ok"])
        failures += [f"{sp}: {p}" for p in problems]

    header = ["species", "truth", "profile", "aligned_frac", "assign_acc", "status"]
    widths = [max(len(str(r[i])) for r in rows + [header]) for i in range(len(header))]
    for r in [header] + rows:
        print("  ".join(str(v).ljust(w) for v, w in zip(r, widths)).rstrip())
    if args.out:
        with open(args.out, "w", newline="\n") as fh:
            for r in [header] + rows:
                fh.write("\t".join(r) + "\n")

    if failures:
        print(f"\nFAIL: {len(failures)} check(s) failed")
        for f in failures:
            print(f"  - {f}")
        sys.exit(1)
    print(f"\nPASS: {len(rows)} species detected, abundances within {args.max_abundance_error}, "
          f"read assignment >= {args.min_assignment_accuracy}")


if __name__ == "__main__":
    main()
