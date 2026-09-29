#!/usr/bin/env python3
"""evaluate.py - how well one protal run of simulated long reads did, as one TSV row.

Counted records are those the profiler takes: not secondary (0x100), MAPQ >= 4, more than 50
aligned bases (M, X, I, D). Their reference bases (M, X, D) are split into those on a gene of the
read's own species and those on another species. Expected are the bases of the reads that lie in
marker genes of their genome (marker_positions.tsv). Abundances are each species' VCov in
<profile>.log over the sum, against the community's cell abundances.
"""

import argparse
import re
import sys


def read_table(path):
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip() and not l.startswith("#")]
    header = lines[0].split("\t")
    return [dict(zip(header, l.split("\t"))) for l in lines[1:]]


def seconds(log):
    """Duration of the alignment stage from protal's log ("Aligning reads took 1m 2s 30ms")."""
    total = 0.0
    for line in open(log, errors="replace"):
        if line.startswith("Aligning reads took"):
            for value, unit in re.findall(r"(\d+)(h|ms|m|s)\b", line):
                total += int(value) * {"h": 3600, "m": 60, "s": 1, "ms": 0.001}[unit]
    return total


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for name in ("sam", "truth", "taxonomy", "genomes", "markers", "profile_log", "community", "label", "log"):
        ap.add_argument("--" + name, required=True)
    ap.add_argument("--header", action="store_true")
    args = ap.parse_args()

    names = {}
    for line in open(args.taxonomy):
        fields = line.rstrip("\n").split("\t")
        if fields[0].isdigit():
            names[fields[0]] = fields[3]
    truth = {row["read"]: row for row in read_table(args.truth)}
    markers = {}
    for row in read_table(args.markers):
        markers.setdefault((row["accession"], row["contig"]), []).append((int(row["start"]), int(row["end"])))
    expected = 0
    for row in truth.values():
        start, end = int(row["start"]), int(row["end"])
        for m_start, m_end in markers.get((row["accession"], row["contig"]), []):
            expected += max(0, min(end, m_end) - max(start, m_start))

    records = kept = 0
    right = wrong = dropped = 0
    for line in open(args.sam):
        if line.startswith("@"):
            continue
        r = line.split("\t", 10)
        records += 1
        flag, mapq, cigar = int(r[1]), int(r[4]), r[5]
        if flag & 0x100:
            continue
        ops = re.findall(r"(\d+)([MIDNSHPX=])", cigar)
        aligned = sum(int(n) for n, op in ops if op in "MXID=")
        ref_bases = sum(int(n) for n, op in ops if op in "MXD=")
        if aligned <= 50:
            continue
        if mapq < 4:
            dropped += ref_bases
            continue
        kept += 1
        read = re.sub(r"_p\d+$", "", r[0])
        species = names.get(r[2].split("_")[0], "?")
        if species == truth[read]["species"]:
            right += ref_bases
        else:
            wrong += ref_bases

    community = {row["accession"]: float(row["relative_abundance"]) for row in read_table(args.community)}
    species_of = {row["accession"]: row["gtdb_taxonomy"].split(";")[-1] for row in read_table(args.genomes)}
    true_abundance = {}
    for acc, a in community.items():
        true_abundance[species_of[acc]] = true_abundance.get(species_of[acc], 0) + a
    vcov = {row["Name"]: float(row["VCov"]) for row in read_table(args.profile_log)}
    total_vcov = sum(vcov.values()) or 1
    estimated = {s: v / total_vcov for s, v in vcov.items()}
    error = sum(abs(estimated.get(s, 0) - true_abundance.get(s, 0)) for s in set(estimated) | set(true_abundance)) / 2

    columns = ["approach", "align_s", "records", "counted", "expected_marker_bp", "own_species_bp", "other_species_bp",
               "mapq_below_4_bp", "recovered", "misassigned", "abundance_error", "abundances"]
    values = [args.label, f"{seconds(args.log):.1f}", records, kept, expected, right, wrong, dropped,
              f"{right / expected:.3f}", f"{wrong / max(1, right + wrong):.4f}", f"{error:.4f}",
              ",".join(f"{s.split('__')[-1]}={estimated.get(s, 0):.3f}" for s in sorted(true_abundance))]
    if args.header:
        print("\t".join(columns))
    print("\t".join(str(v) for v in values))


if __name__ == "__main__":
    sys.exit(main())
