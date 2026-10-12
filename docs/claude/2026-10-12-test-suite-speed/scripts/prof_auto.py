#!/usr/bin/env python3
"""Make TrainerScenariosTest's tables and run the trainer with --features auto under cProfile.
Usage: prof_auto.py SRC OUTDIR"""
import os
import subprocess
import sys

src, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
sys.path.insert(0, os.path.join(src, "scripts"))
import test_model_pmml as t  # noqa: E402

cls = t.TrainerScenariosTest
cls.setUpClass()
print("tables in", cls.tmp.name, flush=True)
cmd = [sys.executable, "-m", "cProfile", "-o", os.path.join(out, "auto.prof"),
       os.path.join(src, "scripts", "machine_learning_cmdline.py"), "--truth-file", cls.training, "--output-prefix",
       os.path.join(out, "auto"), "--features", "auto", "--ntree", "16", "--rounds", "30", "--evaluation", "basic",
       "--threads", "1", "--test-file", cls.test]
import time  # noqa: E402
t0 = time.time()
r = subprocess.run(cmd, capture_output=True, text=True)
print("trainer", r.returncode, round(time.time() - t0, 1), "s", flush=True)
with open(os.path.join(out, "auto.stdout"), "w") as fh:
    fh.write(r.stdout + r.stderr)
import pstats  # noqa: E402
p = pstats.Stats(os.path.join(out, "auto.prof"))
p.sort_stats("cumulative").print_stats(40)
p.sort_stats("tottime").print_stats(25)
