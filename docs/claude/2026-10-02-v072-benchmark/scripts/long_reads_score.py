#!/usr/bin/env python3
"""long_reads_score.py [BENCH_DIR] - scores long_reads.sh's runs (the long-read samples made by 0.7.2's collector,
BENCH_DIR/samples_lr072) as score.py scores the others: BENCH_DIR/results_lr072/runs.tsv, species.tsv, summary.md,
overall.md and paired.md.
"""
import importlib.util
import os
import sys

B = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/bench071")
HERE = os.path.dirname(os.path.abspath(__file__))
# This benchmark's score.py, under another name: it imports the v0.7.1 benchmark's score.py as `score`.
spec = importlib.util.spec_from_file_location("v072_score", os.path.join(HERE, "score.py"))
V = importlib.util.module_from_spec(spec)
spec.loader.exec_module(V)

S = V.S
S.RUNS = os.path.join(B, "runs_lr072")
S.OUT = os.path.join(B, "results_lr072")
S.POINTS = [os.path.join(B, "samples_lr072", "points")]
S.VERSIONS = {"v070": "0.7.0", "v071": "0.7.1", "v072": "0.7.2", "v072k": "0.7.2, --knob 0.5"}
V.PAIRS = [("v071", "v070"), ("v072", "v071"), ("v072k", "v071"), ("v072", "v072k")]

if __name__ == "__main__":
    S.main()
    V.overall()
    V.paired()
