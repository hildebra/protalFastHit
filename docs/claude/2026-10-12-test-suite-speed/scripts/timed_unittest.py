#!/usr/bin/env python3
"""Run unittest modules and write the wall time of every test and of every class (setUpClass included).

Usage: timed_unittest.py OUT.tsv module_or_path [...]
OUT.tsv rows: kind (test|class|module) name seconds outcome.
"""
import sys
import time
import unittest



class TimedResult(unittest.TextTestResult):
    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        self.rows = []
        self.last_event = time.monotonic()
        self.class_start = {}
        self.class_end = {}
        self.current_start = None
        self.outcome = {}

    def startTest(self, test):
        now = time.monotonic()
        cls = type(test).__module__ + "." + type(test).__qualname__
        self.class_start.setdefault(cls, self.last_event)
        self.current_start = now
        super().startTest(test)

    def stopTest(self, test):
        now = time.monotonic()
        cls = type(test).__module__ + "." + type(test).__qualname__
        self.rows.append(("test", test.id(), now - self.current_start, self.outcome.get(test.id(), "ok")))
        self.class_end[cls] = now
        self.last_event = now
        super().stopTest(test)

    def _mark(self, test, what):
        self.outcome[test.id()] = what

    def addFailure(self, test, err):
        self._mark(test, "FAIL")
        super().addFailure(test, err)

    def addError(self, test, err):
        self._mark(test, "ERROR")
        super().addError(test, err)

    def addSkip(self, test, reason):
        self._mark(test, "skip")
        super().addSkip(test, reason)



def main():
    out_path = sys.argv[1]
    names = sys.argv[2:]
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    for name in names:
        if name.endswith(".py"):
            import importlib.util, os
            sys.path.insert(0, os.path.dirname(os.path.abspath(name)))
            mod_name = os.path.splitext(os.path.basename(name))[0]
            spec = importlib.util.spec_from_file_location(mod_name, name)
            mod = importlib.util.module_from_spec(spec)
            sys.modules[mod_name] = mod
            spec.loader.exec_module(mod)
            suite.addTests(loader.loadTestsFromModule(mod))
        else:
            suite.addTests(loader.loadTestsFromName(name))

    runner = unittest.TextTestRunner(verbosity=2, resultclass=TimedResult)
    t0 = time.monotonic()
    result = runner.run(suite)
    total = time.monotonic() - t0
    with open(out_path, "w") as out:
        for row in result.rows:
            out.write("%s\t%s\t%.2f\t%s\n" % row)
        for cls, start in result.class_start.items():
            out.write("class\t%s\t%.2f\t-\n" % (cls, result.class_end[cls] - start))
        out.write("total\t%s\t%.2f\t%s\n" % (" ".join(names), total, "ok" if result.wasSuccessful() else "FAILED"))
    sys.exit(0 if result.wasSuccessful() else 1)


if __name__ == "__main__":
    main()
