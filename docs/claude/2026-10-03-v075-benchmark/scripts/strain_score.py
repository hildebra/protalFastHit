#!/usr/bin/env python3
"""strain_score.py [BENCH_DIR] - the strain runs of this benchmark (strain_runs.sh, BENCH_DIR/strain_runs_v075: 0.7.3
and 0.7.5, the long reads also without phasing) scored by the v0.7.2 benchmark's strain_score.py:
BENCH_DIR/results_strains_v075/species.tsv, summary.md.
"""
import importlib.util
import os
import sys

B = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/bench071")
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "v072_strain_score", os.path.join(HERE, "..", "..", "2026-10-02-v072-benchmark", "scripts", "strain_score.py"))
M = importlib.util.module_from_spec(spec)
spec.loader.exec_module(M)
M.RUNS = os.path.join(B, "strain_runs_v075")
M.OUT = os.path.join(B, "results_strains_v075")
M.VERSIONS = {"v073": "0.7.3", "v073np": "0.7.3, --no_phasing", "v075": "0.7.5", "v075np": "0.7.5, --no_phasing"}

if __name__ == "__main__":
    # The v0.7.2 scorer names the genes by the 0.7.1 pipeline's gene2geneid.tsv, which is gone from BENCH_DIR; the
    # 0.7.3 pipeline's is the same table (the same release and converter), so it stands in.
    v071 = os.path.join(B, "V071", "protal_db", "gene2geneid.tsv")
    if not os.path.exists(v071):
        os.makedirs(os.path.dirname(v071), exist_ok=True)
        os.symlink(os.path.join(B, "V073", "protal_db", "gene2geneid.tsv"), v071)
    M.main()
