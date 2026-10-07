#!/usr/bin/env python3
"""One table of the FP and FN records that extract_error_records.sh took from a build's model_logs/error_reads.

Reads DIR/<type>/<set>/<point>/<sample>.sam.zst (records tagged xg/xs/xe by error_reads.py) and writes
OUT/<type>.records.pkl.gz: one row per record with its alignment (CIGAR counts, clips, mate), its tags (ZU, ZT, ZA,
ZF, ZR), the read's sequence statistics (GC, the share of its most common trinucleotide, longest homopolymer), the
gene's length and the read's source and reasons; and OUT/<type>.taxa.pkl.gz: the samples' FP and FN rows.

    python3 parse_error_records.py --dir ~/v15/v15_error_records --out ~/v15/err_analysis
"""
import argparse
import collections
import glob
import os
import re
import subprocess
import sys

import pandas as pd

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
HOMOPOLYMER = re.compile(r"A+|C+|G+|T+")


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--dir", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--types", default="pe,se,pb,ont")
    return p.parse_args(argv)


def scenario_of(point):
    m = re.match(r"sc_(soil_shallow|soil|gut|host)_", point)
    return m.group(1) if m else "design"


def cigar(c):
    m = x = i = d = opens = 0
    ops = CIGAR.findall(c)
    left = right = 0
    k = 0
    while k < len(ops) and ops[k][1] in "SH":
        left += int(ops[k][0])
        k += 1
    j = len(ops)
    while j > k and ops[j - 1][1] in "SH":
        right += int(ops[j - 1][0])
        j -= 1
    longest_match = 0
    for n, op in ops[k:j]:
        n = int(n)
        if op in "M=":
            m += n
            longest_match = max(longest_match, n)
        elif op == "X":
            x += n
        elif op == "I":
            i += n
            opens += 1
        elif op == "D":
            d += n
            opens += 1
    return m, x, i, d, opens, left, right, longest_match


def seq_stats(seq):
    """(GC share, share of the most common trinucleotide, longest homopolymer) of a read's bases."""
    if not seq or seq == "*":
        return float("nan"), float("nan"), 0
    gc = (seq.count("G") + seq.count("C")) / len(seq)
    tri = collections.Counter(seq[i:i + 3] for i in range(len(seq) - 2))
    top = max(tri.values()) / max(1, len(seq) - 2) if tri else float("nan")
    hp = max((len(h) for h in HOMOPOLYMER.findall(seq)), default=0)
    return gc, top, hp


def records(path):
    proc = subprocess.Popen(["zstd", "-dcq", path], stdout=subprocess.PIPE, text=True)
    lengths = {}
    for line in proc.stdout:
        if line.startswith("@"):
            if line.startswith("@SQ\t"):
                f = dict(x.split(":", 1) for x in line.rstrip("\n").split("\t")[1:] if ":" in x)
                lengths[f.get("SN")] = int(f.get("LN", 0))
            continue
        yield line.rstrip("\n").split("\t"), lengths
    proc.wait()


def parse_sam(path, rt, which, point, sample):
    rows = []
    scenario = scenario_of(point)
    for f, lengths in records(path):
        tags = {t[:2]: t[5:] for t in f[11:]}
        flag = int(f[1])
        row = {"type": rt, "set": which, "point": point, "scenario": scenario, "sample": sample, "qname": f[0],
               "flag": flag, "zu": int(tags.get("ZU", 0)), "zt": int(tags.get("ZT", 0)), "za": tags.get("ZA", ""),
               "zf": tags.get("ZF", ""), "zr": int(tags.get("ZR", 0)), "xg": tags.get("xg", ""),
               "xs": tags.get("xs", ""), "xe": tags.get("xe", "")}
        if flag & 4 or f[2] == "*":
            row["taxon"] = ""
        else:
            taxon, _, gene = f[2].partition("_")
            m, x, i, d, opens, left, right, longest = cigar(f[5])
            gc, tri, hp = seq_stats(f[9])
            row.update(taxon=taxon, gene=int(gene) if gene.isdigit() else -1, pos=int(f[3]), mapq=int(f[4]),
                       matches=m, mismatches=x, ins=i, dels=d, gap_opens=opens, clip_left=left, clip_right=right,
                       longest_match=longest, gene_len=lengths.get(f[2], 0), rnext=f[6], pnext=int(f[7]),
                       tlen=int(f[8]), seq_len=len(f[9]) if f[9] != "*" else 0, gc=gc, top_trinucleotide=tri,
                       homopolymer=hp)
        rows.append(row)
    return rows


def main(argv=None):
    opts = parse_args(argv)
    os.makedirs(opts.out, exist_ok=True)
    for rt in opts.types.split(","):
        rows, taxa = [], []
        sams = sorted(glob.glob(os.path.join(opts.dir, rt, "*", "*", "*.sam.zst")))
        for path in sams:
            which, point = path.split(os.sep)[-3], path.split(os.sep)[-2]
            sample = os.path.basename(path)[:-len(".sam.zst")]
            rows.extend(parse_sam(path, rt, which, point, sample))
            tpath = path[:-len(".sam.zst")] + ".taxa.tsv"
            if os.path.isfile(tpath):
                t = pd.read_csv(tpath, sep="\t", dtype={"taxid": str})
                t["type"], t["point"], t["scenario"] = rt, point, scenario_of(point)
                taxa.append(t)
        r = pd.DataFrame(rows)
        r.to_pickle(os.path.join(opts.out, f"{rt}.records.pkl.gz"))
        pd.concat(taxa, ignore_index=True).to_pickle(os.path.join(opts.out, f"{rt}.taxa.pkl.gz"))
        print(f"{rt}: {len(sams)} SAMs, {len(r)} records, {r['qname'].nunique()} read names", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
