#!/usr/bin/env python3
"""strain_score.py [BENCH_DIR] - the strain benchmark's runs (../../2026-10-02-v072-benchmark/scripts/strain_runs.sh
and this benchmark's strain_runs.sh, all in BENCH_DIR/strain_runs) scored by the v0.7.2 benchmark's strain_score.py,
with 0.7.3 and 0.7.3 without phasing among the versions: BENCH_DIR/results_strains_v073/species.tsv, summary.md.
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
M.OUT = os.path.join(B, "results_strains_v073")
M.VERSIONS = {"v060": "0.6.0a", "v070": "0.7.0", "v071": "0.7.1", "v072": "0.7.2", "v073": "0.7.3",
              "v073np": "0.7.3, --no_phasing"}

if __name__ == "__main__":
    M.main()
