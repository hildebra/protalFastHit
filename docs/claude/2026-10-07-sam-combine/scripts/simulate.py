#!/usr/bin/env python3
"""simulate.py - the samples of the SAM-combination experiment (docs/claude/2026-10-07-sam-combine).

Uses the e2e tests' read simulator on the mini database (tests/e2e/test_protal_e2e.py: 100 bp pairs of 220-320 bp
fragments, 0.5% substitutions, 12 pairs on every gene of 320 bp or more, all three species) and lays the samples out
as three separate studies would have them:

  OUT/run1/sa_1.fq, sa_2.fq, sb_1.fq, sb_2.fq        (seeds 1, 2: the e2e module's sa and sb)
  OUT/run2/sc_1.fq, ...      sd_1.fq, ...            (seeds 11, 12)
  OUT/run3/sa_1.fq, sa_2.fq                          (seed 13: another sample that is also named sa)

  PROTAL_TEST_DB=DB PROTAL=protal python3 simulate.py OUT
"""

import argparse
import os
import shutil
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", ".."))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("out", help="folder to write run1/, run2/ and run3/ into")
    parser.add_argument("--root", default=ROOT, help="protal source tree whose tests/e2e to use (default: this one)")
    args = parser.parse_args()

    sys.path.insert(0, os.path.join(args.root, "tests", "e2e"))
    import test_protal_e2e as t  # reads PROTAL_TEST_DB and PROTAL at import

    t.set_up_database()  # unpacks the database, simulates sa (seed 1) and sb (seed 2) into t.READS
    t.simulate_reads("sc", pairs_per_gene=12, seed=11)
    t.simulate_reads("sd", pairs_per_gene=12, seed=12)
    t.simulate_reads("sa_other", pairs_per_gene=12, seed=13)

    layout = {"run1": [("sa", "sa"), ("sb", "sb")], "run2": [("sc", "sc"), ("sd", "sd")], "run3": [("sa_other", "sa")]}
    for run, samples in layout.items():
        os.makedirs(os.path.join(args.out, run), exist_ok=True)
        for source, name in samples:
            # _1/_2: protal_map_utils generate names sa_R1.fq/sa_R2.fq sample "sa_R" (see the report)
            for mate in ("1", "2"):
                shutil.copy(os.path.join(t.READS, f"{source}_R{mate}.fq"), os.path.join(args.out, run, f"{name}_{mate}.fq"))
    for folder in (t.READS, t.UNPACKED):
        if folder:
            shutil.rmtree(folder, ignore_errors=True)


if __name__ == "__main__":
    main()
