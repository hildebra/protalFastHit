#!/usr/bin/env python3
"""score.py [BENCH_DIR] - this benchmark's runs scored as the v0.7.2 benchmark scores its own
(../../2026-10-02-v072-benchmark/scripts/score.py, which wraps the v0.7.1 benchmark's score.py):
  short  profile.sh's runs, BENCH_DIR/runs_v073 -> BENCH_DIR/results_v073
  lr073  long_reads.sh's runs on the samples of 0.7.3's collector, BENCH_DIR/runs_lr073 -> BENCH_DIR/results_lr073
  lr072  the v0.7.2 benchmark's runs on 0.7.2's collector's samples and 0.7.3's beside them, BENCH_DIR/runs_lr072
         -> BENCH_DIR/results_lr072_v073
Each: runs.tsv, species.tsv, summary.md, overall.md (means per database, read type and version) and paired.md (per
sample, one version less another, with a 95% bootstrap interval).
"""
import importlib.util
import os
import sys

B = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/bench071")
WANTED = sys.argv[2:]  # the v0.7.2 benchmark's score.py sets sys.argv to [script, B]
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "v072_score", os.path.join(HERE, "..", "..", "2026-10-02-v072-benchmark", "scripts", "score.py"))
V = importlib.util.module_from_spec(spec)
spec.loader.exec_module(V)
S = V.S
SHORT_POINTS = list(S.POINTS)

SETS = {
    "short": dict(runs="runs_v073", out="results_v073", points=SHORT_POINTS,
                  versions={"v060": "0.6.0a", "v070": "0.7.0", "v071": "0.7.1", "v072": "0.7.2", "v073": "0.7.3",
                            "v073k": "0.7.3, --knob 0.5"},
                  pairs=[("v070", "v060"), ("v071", "v070"), ("v072", "v071"), ("v073", "v072"), ("v073k", "v072"),
                         ("v073", "v073k")]),
    "lr073": dict(runs="runs_lr073", out="results_lr073", points=[os.path.join(B, "samples_lr073", "points")],
                  versions={"v070": "0.7.0", "v071": "0.7.1", "v072": "0.7.2", "v072k": "0.7.2, --knob 0.5",
                            "v073": "0.7.3", "v073k": "0.7.3, --knob 0.5"},
                  pairs=[("v071", "v070"), ("v072", "v071"), ("v073", "v072"), ("v073k", "v072k"), ("v073", "v073k"),
                         ("v072", "v072k")]),
    "lr072": dict(runs="runs_lr072", out="results_lr072_v073", points=[os.path.join(B, "samples_lr072", "points")],
                  versions={"v070": "0.7.0", "v071": "0.7.1", "v072": "0.7.2", "v072k": "0.7.2, --knob 0.5",
                            "v073": "0.7.3", "v073k": "0.7.3, --knob 0.5"},
                  pairs=[("v073", "v072"), ("v073k", "v072k"), ("v073", "v073k")]),
}

if __name__ == "__main__":
    for name in WANTED or SETS:
        c = SETS[name]
        S.RUNS, S.OUT, S.POINTS, S.VERSIONS = os.path.join(B, c["runs"]), os.path.join(B, c["out"]), c["points"], c["versions"]
        V.PAIRS = c["pairs"]
        print(f"## {name}: {S.RUNS} -> {S.OUT}")
        S.main()
        V.overall()
        V.paired()
