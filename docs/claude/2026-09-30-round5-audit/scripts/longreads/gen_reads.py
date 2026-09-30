"""Simulates long reads from the mini GTDB's genomes with a per-base truth map.

  gen_reads.py community OUT --gtdb G --db D --reads N --platform hifi|ont [--community acc:w,...]
  gen_reads.py sliding OUT --gtdb G --db D --acc ACC --length L --offsets a:b:step [--platform ...]

Writes OUT.fq, OUT.fa (no qualities) and OUT.truth.pkl: per read (name, acc, strand, blocks) where
blocks are (read_start, chrom_start, length) runs of read bases (forward read coordinates) that
come from the genome's virtual chromosome (its contigs concatenated in file order) at chrom_start
onward; inserted bases are in no block. Qualities are random per base (so a QUAL orientation error
shows)."""
import argparse
import math
import os
import pickle
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lrlib import World, mutate, revcomp

PLATFORMS = {
    # error, substitution share, insertion share, homopolymer loss, quality range
    "perfect": (0.0, 0.4, 0.3, 0.0, (30, 40)),
    "hifi": (0.001, 0.4, 0.3, 0.0, (20, 41)),
    "ont": (0.016, 0.30, 0.25, 0.15, (8, 30)),
    "ont5": (0.05, 0.30, 0.25, 0.15, (5, 25)),
    "ont10": (0.10, 0.30, 0.25, 0.15, (3, 20)),
}


def chromosome(world, acc):
    contigs = world.genomes[acc]
    offset, seq, where = 0, [], {}
    for name, s in contigs.items():
        where[name] = offset
        seq.append(s)
        offset += len(s)
    return "".join(seq), where


def make_read(rng, chrom, start, length, strand, platform):
    error, sub, ins, hp, (qlo, qhi) = PLATFORMS[platform]
    piece = chrom[start:start + length]
    read, origin = mutate(piece, rng, error, sub, ins, hp)
    origin = [start + o if o >= 0 else -1 for o in origin]
    if strand == "-":
        read = revcomp(read)
        origin = origin[::-1]
    # blocks: runs where origin increases (forward) or decreases (reverse) by 1
    blocks = []
    i = 0
    n = len(origin)
    while i < n:
        if origin[i] < 0:
            i += 1
            continue
        j = i + 1
        step = 1 if strand == "+" else -1
        while j < n and origin[j] >= 0 and origin[j] == origin[j - 1] + step:
            j += 1
        blocks.append((i, origin[i], j - i))
        i = j
    qual = "".join(chr(33 + rng.randint(qlo, qhi)) for _ in read)
    return read, qual, blocks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode")
    ap.add_argument("out")
    ap.add_argument("--gtdb", required=True)
    ap.add_argument("--db", required=True)
    ap.add_argument("--platform", default="hifi")
    ap.add_argument("--reads", type=int, default=300)
    ap.add_argument("--community", default="GCA_999001002.1:0.5,GCF_999002001.1:0.3,GCA_999003003.1:0.2")
    ap.add_argument("--median", type=float, default=15000)
    ap.add_argument("--sigma", type=float, default=0.2)
    ap.add_argument("--min_len", type=int, default=3000)
    ap.add_argument("--max_len", type=int, default=150000)
    ap.add_argument("--acc")
    ap.add_argument("--length", type=int)
    ap.add_argument("--offsets")
    ap.add_argument("--strands", default="+-")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--name", default="read{}")
    a = ap.parse_args()
    rng = random.Random(a.seed)
    world = World(a.gtdb, a.db)
    chroms = {}
    jobs = []  # (acc, start, length, strand)
    if a.mode == "community":
        accs, weights = [], []
        for item in a.community.split(","):
            acc, w = item.split(":")
            accs.append(acc)
            weights.append(float(w))
        for acc in accs:
            chroms[acc] = chromosome(world, acc)
        # cells: reads proportional to weight x genome length (all similar here) -> just weight
        for _ in range(a.reads):
            acc = rng.choices(accs, weights)[0]
            chrom = chroms[acc][0]
            length = int(min(a.max_len, max(a.min_len, math.exp(rng.gauss(math.log(a.median), a.sigma)))))
            length = min(length, len(chrom))
            start = rng.randint(0, len(chrom) - length)
            jobs.append((acc, start, length, rng.choice("+-")))
    elif a.mode == "sliding":
        chroms[a.acc] = chromosome(world, a.acc)
        lo, hi, step = (int(x) for x in a.offsets.split(":"))
        k = 0
        for off in range(lo, hi, step):
            if off + a.length > len(chroms[a.acc][0]):
                break
            jobs.append((a.acc, off, a.length, a.strands[k % len(a.strands)]))
            k += 1
    else:
        sys.exit("mode?")
    truth = []
    with open(a.out + ".fq", "w") as fq, open(a.out + ".fa", "w") as fa:
        for i, (acc, start, length, strand) in enumerate(jobs, 1):
            read, qual, blocks = make_read(rng, chroms[acc][0], start, length, strand, a.platform)
            name = a.name.format(i)
            fq.write(f"@{name}\n{read}\n+\n{qual}\n")
            fa.write(f">{name}\n{read}\n")
            truth.append({"name": name, "acc": acc, "strand": strand, "start": start, "length": length,
                          "read_length": len(read), "blocks": blocks})
    offsets = {acc: chroms[acc][1] for acc in chroms}
    with open(a.out + ".truth.pkl", "wb") as fh:
        pickle.dump({"reads": truth, "contig_offsets": offsets, "platform": a.platform}, fh)
    print(f"{len(jobs)} reads, {sum(t['read_length'] for t in truth)} bp -> {a.out}.fq")


if __name__ == "__main__":
    main()
