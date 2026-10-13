#!/usr/bin/env python3
"""model_fit.py - what a stride pick could gain over the hash sample, predicted from an older probe's rule table.

A read (pair) with n "good" lookups (best flex score >= T) passes rule (c, T, m) of the hash sample with probability
P(Binomial(n, 1/c) >= m): each core is in the sample on its own. From the table's recall of every rule at one T (all c
and m), fit the distribution p(n) of the good lookups among the reads with records (non-negative, summing to 1, with a
light second-difference penalty above n = 5), then predict the share of those reads a gate at m = 1 loses:

  random              the hash sample (the fit's check against the table)
  stride_run          a stride over good lookups that lie in one run of the read: lost only if n < c, with
                      probability 1 - n/c
  stride_scattered    a pick of a fixed number of the read's `--cores` cores (cores/c) with the good lookups scattered
                      among them (hypergeometric)

Real reads lie between the two stride columns. A model, not a measurement: protal's --stage1_probe since 2026-10-13
measures the stride's rules directly.

  python3 model_fit.py local/Perf0.8.2d/stage1_probe/pe1_stage1_probe.tsv --cores 50
  python3 model_fit.py local/Perf0.8.2d/stage1_probe/pb2_stage1_probe.tsv --max-good 160 --cores 2000

docs/claude/2026-10-13-stride-pick.
"""

import argparse
import csv
import math
import sys

import numpy as np


def p_at_least(n, c, m):
    """P(Binomial(n, 1/c) >= m)."""
    if c == 1:
        return 1.0 if n >= m else 0.0
    p = 1.0 / c
    below = sum(math.comb(n, j) * p ** j * (1 - p) ** (n - j) for j in range(min(m, n + 1)))
    return 1.0 - below


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("table", help="a probe's <sample>_stage1_probe.tsv")
    parser.add_argument("--score", type=int, default=14, help="T (14)")
    parser.add_argument("--max-good", type=int, default=80, help="the most good lookups a read has in the fit (80)")
    parser.add_argument("--cores", type=int, default=50, help="cores per read for the scattered case (50: a pair)")
    parser.add_argument("--iterations", type=int, default=400000)
    parser.add_argument("--smooth", type=float, default=1e-3)
    args = parser.parse_args()

    with open(args.table, newline="") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if int(r["min_score"]) == args.score]
    obs = [(int(r["c"]), int(r["min_hits"]), float(r["recall"])) for r in rows]
    ns = np.arange(args.max_good + 1)
    a = np.array([[p_at_least(int(n), c, m) for n in ns] for c, m, _ in obs])
    y = np.array([r for _, _, r in obs])

    # Exponentiated gradient keeps p on the simplex; it starts uniform.
    p = np.full(len(ns), 1.0 / len(ns))
    rate = 0.02
    for _ in range(args.iterations):
        r = a @ p - y
        g = 2 * a.T @ r
        d = p[5:-2] - 2 * p[6:-1] + p[7:]  # second differences at n = 6 .. max-1
        g[5:-2] += 2 * args.smooth * d
        g[6:-1] -= 4 * args.smooth * d
        g[7:] += 2 * args.smooth * d
        p = p * np.exp(-rate * g)
        p /= p.sum()

    print("c\tm\tobserved\tfitted")
    for (c, m, rec), fit in zip(obs, a @ p):
        print(f"{c}\t{m}\t{rec:.4f}\t{fit:.4f}")
    cum = np.cumsum(p)
    print(f"P(n<=4) {cum[4]:.4f}  P(n<=7) {cum[7]:.4f}  P(n<=15) {cum[15]:.4f}  P(n<=31) {cum[min(31, args.max_good)]:.4f}  "
          f"mean {float(ns @ p):.1f}")
    print("c\tobserved_loss\trandom\tstride_run\tstride_scattered")
    observed = {c: 1 - rec for c, m, rec in obs if m == 1}
    for c in (2, 4, 8, 16, 32):
        rand = float(p @ (1 - 1 / c) ** ns)
        run = float(p @ np.maximum(0.0, 1 - ns / c))
        k = args.cores / c
        lo, frac = int(k), k - int(k)

        def none_good(n, draws):
            total = max(args.cores, n)
            v = 1.0
            for i in range(draws):
                v *= max(0, total - n - i) / (total - i)
            return v

        scattered = float(sum(p[n] * ((1 - frac) * none_good(n, lo) + frac * none_good(n, lo + 1)) for n in ns))
        print(f"{c}\t{observed.get(c, float('nan')):.4f}\t{rand:.4f}\t{run:.4f}\t{scattered:.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
