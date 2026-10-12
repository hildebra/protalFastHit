#!/usr/bin/env python3
"""Run some GtdbBuildTest tests with the class folder kept (tearDownClass replaced); print the folder.
Usage: keep_gtdb.py SRC test_a [test_i ...]"""
import os
import sys
import time
import unittest

src = sys.argv[1]
sys.path.insert(0, os.path.join(src, "scripts", "mini_db"))
import test_gtdb_pipeline as t  # noqa: E402

cls = t.GtdbBuildTest
cls.tearDownClass = classmethod(lambda c: print("kept", c.tmp.name, flush=True))
suite = unittest.TestSuite()
names = [n for n in dir(cls) if n.startswith("test_")]
for prefix in sys.argv[2:]:
    suite.addTests(cls(n) for n in names if n.startswith(prefix))
start = time.time()
result = unittest.TextTestRunner(verbosity=2).run(suite)
print("wall", round(time.time() - start, 1), "s", flush=True)
# keep the TemporaryDirectory from being removed at exit
cls.tmp._finalizer.detach()
