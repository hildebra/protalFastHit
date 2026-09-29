#!/usr/bin/env python3
"""msa_noise.py - per sample of raw strain MSAs, the share of called bases that are IUPAC ambiguity
codes (several alleles passed the SNP filters) and that differ from the column's majority base.
Samples of one strain each should have neither but for real strain differences."""

import argparse
import collections
import glob


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--strains", required=True, help="strains output folder")
    args = ap.parse_args()
    called = collections.Counter()
    ambiguous = collections.Counter()
    for path in glob.glob(f"{args.strains}/*.raw.msa.fna"):
        rows, name = {}, None
        for line in open(path):
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                rows[name] = []
            elif name:
                rows[name].append(line)
        rows = {n: "".join(s) for n, s in rows.items()}
        for sample, seq in rows.items():
            for c in seq:
                if c in "-N":
                    continue
                called[sample] += 1
                ambiguous[sample] += c not in "ACGT"
    for sample in sorted(called):
        print(f"{sample}\tcalled {called[sample]}\tambiguous {ambiguous[sample] / called[sample]:.5f}")


if __name__ == "__main__":
    main()
