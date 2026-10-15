#!/usr/bin/env python3
"""Audit probe: build_gtdb_database.py end to end on the pipeline test's synthetic release, with the layout of the r226
v28 build that the pipeline test does not cover: paired-end and single-end setups of their own (points pe only, both,
se only), the design's Ultima units, pb and ont, scenarios with Ultima/pb/ont reads, both model kinds, models at other
gates, --test-c-gates, some simulations streamed, error reads and the share archive.

usage: probe_build.py WORKDIR NAME [extra build options...]   (PROTAL, SIMULATE, PROTAL_TRAIN_PYTHON in the environment)
The inputs are made once in WORKDIR; each NAME is a build into WORKDIR/NAME (console in WORKDIR/NAME.log).
"""

import json
import os
import sys
import tempfile
import types

SRC = os.path.expanduser(os.environ.get("PROBE_SRC", "~/bgaudit/src/scripts/mini_db"))
sys.path.insert(0, SRC)
work, name, extra = sys.argv[1], sys.argv[2], sys.argv[3:]
os.makedirs(work, exist_ok=True)


class Kept:
    """A TemporaryDirectory that stays: the probe's outputs are looked at afterwards."""
    def __init__(self, *a, **k):
        self.name = work

    def cleanup(self):
        pass


import test_gtdb_pipeline as pipeline  # noqa: E402

T = pipeline.GtdbBuildTest
if os.path.isdir(os.path.join(work, "inputs")):
    T.python = os.environ.get("PROTAL_TRAIN_PYTHON", sys.executable)
    T.tmp, T.inputs, T.full = Kept(), os.path.join(work, "inputs"), None
else:
    pipeline.tempfile = types.SimpleNamespace(TemporaryDirectory=Kept, mkdtemp=tempfile.mkdtemp)
    T.setUpClass()
definitions, host = T.scenario_inputs(T)
# The scenarios with long reads too (the v28 build's gut and host have pe, Ultima, pb and ont reads).
with open(definitions) as fh:
    scenarios = json.load(fh)
for s in ("gut", "host"):
    kinds = {r["type"] for r in scenarios[s]["reads"]}
    for kind in ("pb", "ont"):
        if kind not in kinds:
            scenarios[s]["reads"].append({"type": kind, "depth": 300000})
with open(definitions, "w") as fh:
    json.dump(scenarios, fh)

scratch = os.path.join(work, name + "_scratch")
options = ["--scratch", scratch, "--read-types", "pe,se,pb,ont",
           "--read-setups", "100:HS20:300:40,150:HSXt:350:50", "--se-read-setups", "150:HSXt:350:50,250:MSv3:550:50",
           "--se-ultima-setup", "ultima:300:50:25:2",
           "--long-read-bases", "100000,400000", "--long-read-samples", "2",
           "--test-long-read-bases", "200000", "--test-long-read-samples", "1",
           "--model", "forest,gbm", "--model-gates", "pe:1/8,se:1/8,pb:8/32,ont:1/16", "--test-c-gates",
           "--scenario-file", definitions, "--scenario-samples", "2", "--scenario-test-samples", "1",
           "--host-genome", host, "--profile-blocks", "0.00001", "--share-logs", *extra]
result = T.build(name, *options, scenarios=True, error_reads="all", foreign=False)
print(f"exit {result.returncode}")
print(result.stdout[-6000:])
sys.exit(result.returncode)
