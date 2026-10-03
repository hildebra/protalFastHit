#!/usr/bin/env python3
"""score.py [BENCH_DIR] [SET ...] - this benchmark's runs scored as the v0.7.3 benchmark scores its own (through the
v0.7.2 benchmark's score.py, which wraps the v0.7.1 benchmark's):
  short  profile.sh's runs, BENCH_DIR/runs_v075 -> BENCH_DIR/results_v075 (0.6.0a, 0.7.3, 0.7.5)
  lr     long_reads.sh's runs on the samples of 0.7.3's collector, BENCH_DIR/runs_lr075 -> BENCH_DIR/results_lr075
Each: runs.tsv, species.tsv, summary.md, overall.md (means per database, read type and version) and paired.md (per
sample, one version less another, with a 95% bootstrap interval).
"""
import importlib.util
import os
import sys

B = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/bench071")
WANTED = sys.argv[2:]
HERE = os.path.dirname(os.path.abspath(__file__))

# The v0.7.1 scorer reads the species the 0.7.1 pipeline held out (V071/heldout_species.txt), gone from BENCH_DIR;
# the 0.7.3 pipeline held out the same 290 (the v0.7.3 benchmark checked), so its list stands in. 0.7.5's pipeline
# holds out 349 (its own design), so for its runs against its training database that list is the one that applies:
# the scorer's set is swapped per run.
v071 = os.path.join(B, "V071", "heldout_species.txt")
if not os.path.exists(v071):
    os.makedirs(os.path.dirname(v071), exist_ok=True)
    os.symlink(os.path.join(B, "V073", "heldout_species.txt"), v071)

spec = importlib.util.spec_from_file_location(
    "v072_score", os.path.join(HERE, "..", "..", "2026-10-02-v072-benchmark", "scripts", "score.py"))
V = importlib.util.module_from_spec(spec)
spec.loader.exec_module(V)
S = V.S
SHORT_POINTS = list(S.POINTS)


def heldout_of(pipeline):
    with open(os.path.join(B, pipeline, "heldout_species.txt")) as fh:
        return {S.species_of(line.split("\t")[0]) for line in fh if line.strip()}


HELDOUT = {"v075": heldout_of("V075"), "default": heldout_of("V073")}
_score = S.score


def score_with_the_runs_heldout(run, species_rows):
    S.heldout = HELDOUT.get(run.split(".", 1)[0], HELDOUT["default"])
    return _score(run, species_rows)


S.score = score_with_the_runs_heldout

SETS = {
    "short": dict(runs="runs_v075", out="results_v075", points=SHORT_POINTS,
                  versions={"v060": "0.6.0a", "v073": "0.7.3", "v075": "0.7.5"},
                  pairs=[("v073", "v060"), ("v075", "v073")]),
    "lr": dict(runs="runs_lr075", out="results_lr075", points=[os.path.join(B, "samples_lr073", "points")],
               versions={"v073": "0.7.3", "v075": "0.7.5"},
               pairs=[("v075", "v073")]),
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
