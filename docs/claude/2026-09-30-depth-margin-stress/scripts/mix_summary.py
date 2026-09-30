#!/usr/bin/env python3
"""mix_summary.py LABEL... - run B's mixtures from evaluate.py's <label>.mix.tsv (in ~/audit5/accuracy/results): at
sites where the two strains differ, the share of called cells that carry both alleles (the IUPAC code of the two)
or the major's base alone; at sites where they agree, wrong calls (another base, or an IUPAC code) per million
called cells. By stage (raw, qc) and the minor strain's share."""
import collections, csv, os, sys
R = os.path.expanduser("~/audit5/accuracy/results")
DIFFER = {"major_alt_only", "minor_alt_only", "both_alt_differ"}
AGREE = {"both_ref", "shared_alt"}
CALLED = {"major_base", "minor_base", "iupac_both", "other", "iupac_other"}
for label in sys.argv[1:]:
    c = collections.defaultdict(collections.Counter)
    for r in csv.DictReader(open(os.path.join(R, label + ".mix.tsv")), delimiter="\t"):
        c[(r["stage"], r["minor_frac"], r["site"] in DIFFER, r["site"] in AGREE)][r["call"]] += int(r["count"])
    print(f"== {label}")
    print("stage minor  differ: both alleles / major only / minor only   agree: wrong per million called")
    for stage in ("raw", "qc"):
        for frac in sorted({k[1] for k in c}, key=float):
            d = c[(stage, frac, True, False)]
            a = c[(stage, frac, False, True)]
            dc = sum(d[x] for x in CALLED) or 1
            ac = sum(a[x] for x in CALLED) or 1
            wrong = sum(a[x] for x in ("other", "iupac_both", "iupac_other", "minor_base"))
            print(f"{stage:5s} {float(frac):5.2f}  {d['iupac_both'] / dc:.3f} / {d['major_base'] / dc:.3f} / {d['minor_base'] / dc:.3f}"
                  f"          {1e6 * wrong / ac:8.0f}")
