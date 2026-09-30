"""Prefix each line of stdin with the seconds since start (monotonic clock: WSL's wall clock jumps
when it resyncs). Usage: protal ... 2>&1 | python3 stamp.py > log"""
import sys
import time

t0 = time.monotonic()
out = sys.stdout
for line in sys.stdin:
    out.write(f"{time.monotonic() - t0:8.2f} {line}")
