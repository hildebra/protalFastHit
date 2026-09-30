#!/usr/bin/env python3
"""Helpers for the I/O audit: simulated reads from the database genes (as tests/e2e does), and a
BGZF writer (bgzip is not installed).

  mkreads.py reads <reference.fna> <outprefix> <pairs_per_gene> <seed>
  mkreads.py bgzf <in> <out>            # BGZF with EOF block
  mkreads.py bgzf-noeof <in> <out>      # BGZF without EOF block
"""
import random
import struct
import sys
import zlib


def revcomp(seq):
    return seq[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def genes(path):
    out, name, seq = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    out.append((name, "".join(seq).upper()))
                name, seq = line[1:].split()[0], []
            elif line:
                seq.append(line)
    if name:
        out.append((name, "".join(seq).upper()))
    return out


def reads(ref, prefix, pairs_per_gene, seed):
    rng = random.Random(seed)

    def mutate(seq):
        bases, qual = list(seq), ["I"] * len(seq)
        for i, b in enumerate(bases):
            if rng.random() < 0.005:
                bases[i] = rng.choice([c for c in "ACGT" if c != b])
            if rng.random() < 0.01:
                qual[i] = "#"
        return "".join(bases), "".join(qual)

    n = 0
    with open(f"{prefix}_R1.fq", "w") as r1, open(f"{prefix}_R2.fq", "w") as r2:
        for _, gene in genes(ref):
            if len(gene) < 320:
                continue
            for _ in range(pairs_per_gene):
                flen = rng.randint(220, 320)
                start = rng.randint(0, len(gene) - flen)
                frag = gene[start:start + flen]
                if rng.random() < 0.5:
                    frag = revcomp(frag)
                s1, q1 = mutate(frag[:100])
                s2, q2 = mutate(revcomp(frag[-100:]))
                n += 1
                r1.write(f"@{n}/1\n{s1}\n+\n{q1}\n")
                r2.write(f"@{n}/2\n{s2}\n+\n{q2}\n")
    print(n)


def bgzf_block(data):
    c = zlib.compressobj(6, zlib.DEFLATED, -15)
    comp = c.compress(data) + c.flush()
    bsize = 18 + len(comp) + 8 - 1
    head = bytes([0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, ord('B'), ord('C'), 2, 0]) + struct.pack("<H", bsize)
    return head + comp + struct.pack("<II", zlib.crc32(data) & 0xffffffff, len(data))


def bgzf(src, dst, eof=True):
    data = open(src, "rb").read()
    with open(dst, "wb") as out:
        for i in range(0, len(data), 0xff00):
            out.write(bgzf_block(data[i:i + 0xff00]))
        if eof:
            out.write(bgzf_block(b""))


def seekable(src, dst, chunk):
    """zstd seekable format with frames of `chunk` content bytes (zstd CLI per frame)."""
    import subprocess
    data = open(src, "rb").read()
    entries = []
    with open(dst, "wb") as out:
        for i in range(0, len(data), chunk):
            part = data[i:i + chunk]
            comp = subprocess.run(["zstd", "-q", "-c", "--no-check"], input=part, stdout=subprocess.PIPE, check=True).stdout
            out.write(comp)
            entries.append((len(comp), len(part)))
        body = b"".join(struct.pack("<II", c, d) for c, d in entries)
        table = body + struct.pack("<IBI", len(entries), 0, 0x8F92EAB1)
        out.write(struct.pack("<II", 0x184D2A5E, len(table)) + table)
    print("frames", len(entries))


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "seekable":
        seekable(sys.argv[2], sys.argv[3], int(sys.argv[4]))
    if cmd == "reads":
        reads(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]))
    elif cmd == "bgzf":
        bgzf(sys.argv[2], sys.argv[3])
    elif cmd == "bgzf-noeof":
        bgzf(sys.argv[2], sys.argv[3], eof=False)
