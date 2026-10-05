#!/usr/bin/env python3
"""Synthetic genomes of real size for timing the collector's long-read code: random sequence, NCBI's layout (80
bases per line, gzip level 6), lengths and contig counts like GTDB r226 genomes (median ~3.4 Mb, ~80 contigs).

usage: make_genomes.py OUT_DIR --genomes 400 [--threads 6]
Writes OUT_DIR/g<i>.fna.gz and OUT_DIR/genomes.tsv (name, lineage, path, length).
"""

import argparse
import concurrent.futures
import gzip
import os

import numpy as np

LETTERS = np.frombuffer(b"ACGT", dtype=np.uint8)


def one(args):
    out_dir, i = args
    rng = np.random.default_rng(1000 + i)
    length = int(np.clip(rng.lognormal(np.log(3.4e6), 0.35), 1.2e6, 9e6))
    contigs = int(np.clip(rng.lognormal(np.log(78), 0.8), 1, 600))
    cuts = np.sort(rng.choice(np.arange(1, length), contigs - 1, replace=False)) if contigs > 1 else np.array([], int)
    bounds = np.concatenate(([0], cuts, [length]))
    seq = LETTERS[rng.integers(0, 4, length)].tobytes()
    parts = []
    for c in range(contigs):
        s = seq[bounds[c]:bounds[c + 1]]
        parts.append(b">c%d_%d synthetic\n" % (i, c + 1) + b"\n".join(s[j:j + 80] for j in range(0, len(s), 80)) + b"\n")
    path = os.path.join(out_dir, f"g{i}.fna.gz")
    with open(path, "wb") as fh:
        fh.write(gzip.compress(b"".join(parts), compresslevel=6))
    return f"G{i:05d}", path, length


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("out")
    ap.add_argument("--genomes", type=int, default=400)
    ap.add_argument("--threads", type=int, default=6)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    with concurrent.futures.ProcessPoolExecutor(a.threads) as ex:
        rows = list(ex.map(one, [(a.out, i) for i in range(a.genomes)]))
    with open(os.path.join(a.out, "genomes.tsv"), "w") as fh:
        for name, path, length in rows:
            lineage = f"d__Bacteria;p__P;c__C;o__O;f__F;g__G{name};s__G{name} sp"
            fh.write(f"{name}\t{lineage}\t{path}\t{length}\n")
    print(f"{len(rows)} genomes in {a.out}")


if __name__ == "__main__":
    main()
