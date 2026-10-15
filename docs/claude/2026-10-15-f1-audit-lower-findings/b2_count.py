#!/usr/bin/env python3
"""The fragments of a b2_reads.py sample that protal aligned (any record), and those whose primary record is on their
own species (the taxid in the read's name), from its SAM (.sam, or .sam.zst read through zstdcat).

usage: b2_count.py SAM FRAGMENTS    -> prints: fragments<TAB>aligned<TAB>own_species
"""
import subprocess
import sys


def lines(path):
    if path.endswith(".zst"):
        proc = subprocess.Popen(["zstdcat", path], stdout=subprocess.PIPE, text=True)
        yield from proc.stdout
        proc.wait()
    else:
        with open(path) as fh:
            yield from fh


def main():
    path, fragments = sys.argv[1], int(sys.argv[2])
    aligned, own = set(), set()
    for line in lines(path):
        if line.startswith("@"):
            continue
        f = line.split("\t", 4)
        flag = int(f[1])
        if flag & 0x4 or f[2] == "*":
            continue
        name = f[0].rsplit("/", 1)[0] if f[0].endswith(("/1", "/2")) else f[0]
        aligned.add(name)
        if not flag & 0x900 and f[2].split("_")[0] == name.split("_")[0]:
            own.add(name)
    print(f"{fragments}\t{len(aligned)}\t{len(own)}")


if __name__ == "__main__":
    sys.exit(main())
