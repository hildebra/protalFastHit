#!/usr/bin/env python3
"""error_breakdown.py [READ_TYPES...] (test_errors.py until 2026-10-06; not a test) - the errors of 0.7.1's model
(refitted as the trainer fits it, seed 1) on its pipeline's
independent test set ($V071/test): false negatives by the taxon's fragments, by whether it was simulated from its
representative (meta_rep_genome), by domain and by whether a congener the database lacks is in the sample;
false positives by how the taxon relates to what was simulated (meta_relative_rank, meta_novel_congener); and
per design point the F1 at 0.5 against the F1-optimal threshold of that point (an oracle: how much a depth-aware threshold could recover at most). Writes results/test_errors.md.
"""
import collections
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
os.environ["V2"] = os.path.expanduser(os.environ.get("V071", "~/bench071/V071"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "2026-10-01-alignment-features"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "..", "scripts"))
import exp_lib as X  # noqa: E402
from model_features import NORMALIZED_FEATURES as BASE  # noqa: E402

READ_TYPES = sys.argv[1:] or ["pe", "se", "pb", "ont"]


def fbin(n):
    return "1" if n <= 1 else "2-3" if n <= 3 else "4-10" if n <= 10 else "11-100" if n <= 100 else ">100"


def main():
    lines = ["# Errors on the independent test set (test_errors.py)", ""]
    for rt in READ_TYPES:
        train, test = X.load("training", rt), X.load("test", rt)
        p = X.test_predict(train, test, BASE, 1)
        y = test["truth"].to_numpy()
        call = p >= X.KNOB
        test = test.assign(p=p, call=call)
        fn, fp = test[(y == 1) & ~call], test[(y == 0) & call]
        present = test[y == 1]
        lines += [f"## {rt}: {len(fn)} false negatives, {len(fp)} false positives, F1 {X.scores(y, call)['F1']:.4f}", ""]
        rows = []
        for label, col in (("fragments", None), ("simulated from the representative", "meta_rep_genome"), ("domain", "meta_domain")):
            if col is None:
                miss = collections.Counter(fbin(n) for n in fn["fragments"])
                allp = collections.Counter(fbin(n) for n in present["fragments"])
            else:
                miss = collections.Counter(str(v) for v in fn[col])
                allp = collections.Counter(str(v) for v in present[col])
            rows.append(f"- missed by {label}: " + ", ".join(f"{k} {miss[k]} of {allp[k]} ({miss[k] / allp[k]:.1%})"
                                                            for k in sorted(allp)))
        for nc, label in ((0, "no congener of theirs missing from the database"), (1, "a missing congener in the sample")):
            g = present[present["meta_novel_congener"] == nc]
            gm = g[~g["call"]]
            rows.append(f"- present taxa with {label}: {len(g)}, missed {len(gm)} ({len(gm) / max(1, len(g)):.1%}), "
                        f"of them with more than 10 fragments {int((gm['fragments'] > 10).sum())}")
        for col in ("meta_relative_rank", "meta_novel_congener"):
            if col in fp:
                c = collections.Counter(str(v) for v in fp[col])
                rows.append(f"- false positives by {col}: " + ", ".join(f"{k} {v}" for k, v in sorted(c.items())))
        rows.append(f"- false positives' fragments: " + ", ".join(f"{k} {v}" for k, v in
                                                                 sorted(collections.Counter(fbin(n) for n in fp["fragments"]).items())))
        lines += rows + ["", "| design point | F1 at 0.5 | best threshold: F1 |", "|---|---|---|"]
        grid = np.arange(0.05, 0.96, 0.01)
        for design, g in test.groupby("meta_design"):
            yy, pp = g["truth"].to_numpy(), g["p"].to_numpy()
            f = [X.scores(yy, pp >= t)["F1"] for t in grid]
            lines.append(f"| {design} | {X.scores(yy, pp >= 0.5)['F1']:.4f} | {grid[int(np.argmax(f))]:.2f}: {max(f):.4f} |")
        lines.append("")
    out = os.path.join(HERE, "..", "results", "test_errors.md")
    with open(out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
