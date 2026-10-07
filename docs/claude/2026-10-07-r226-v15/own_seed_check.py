#!/usr/bin/env python3
"""Did the unseen species' own reads seed on their own marker genes? Per sample of error_reads.py's output, for each
unseen species of the training database (not held out) whose own fragments aligned nowhere: its own records (xs tag),
and those whose ZF tag (taxa seeded on but not aligned to) names the species itself. A species whose own reads never
seeded on its own genes was missed for want of reads on its markers (depth); one whose own reads did seed there and
failed to align was missed by the aligner.

    python3 own_seed_check.py --heldout local/v15/heldout_species.txt ~/v15/model_logs/error_reads/ont/test/sc_soil_*/*.taxa.tsv
"""
import argparse
import collections
import os
import subprocess

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("taxa", nargs="+", help="<sample>.taxa.tsv files (their .sam.zst beside them)")
    ap.add_argument("--heldout", required=True)
    opts = ap.parse_args()
    with open(opts.heldout) as fh:
        heldout = {line.split("\t")[0].strip() for line in fh if line.strip()}
    out = []
    for path in opts.taxa:
        t = pd.read_csv(path, sep="\t", dtype={"taxid": str})
        u = t[(t["error"] == "unseen") & ~t["taxon_name"].isin(heldout)].copy()
        for c in ("own_fragments", "own_best_on_taxon", "own_best_elsewhere", "own_unaligned"):
            u[c] = pd.to_numeric(u[c], errors="coerce").fillna(0)
        u = u[(u["own_fragments"] > 0) & (u["own_best_on_taxon"] == 0) & (u["own_best_elsewhere"] == 0)]
        want = dict(zip(u["taxon_name"], u["taxid"]))
        records, own_seeded, reads_seeded = collections.Counter(), collections.Counter(), collections.defaultdict(set)
        sam = path[:-len(".taxa.tsv")] + ".sam.zst"
        with subprocess.Popen(["zstd", "-dc", sam], stdout=subprocess.PIPE) as proc:
            for line in proc.stdout:
                if line.startswith(b"@"):
                    continue
                i = line.find(b"\txs:Z:")
                if i < 0:
                    continue
                species = line[i + 6:line.find(b"\t", i + 6) if line.find(b"\t", i + 6) >= 0 else None].rstrip(b"\r\n").decode()
                if species not in want:
                    continue
                records[species] += 1
                z = line.find(b"\tZF:Z:")
                if z >= 0:
                    end = line.find(b"\t", z + 6)
                    failed = line[z + 6:end if end >= 0 else None].rstrip(b"\r\n").split(b",")
                    if want[species].encode() in failed:
                        own_seeded[species] += 1
                        reads_seeded[species].add(line.split(b"\t", 1)[0])
        n = len(want)
        hit = sum(1 for s in want if own_seeded[s] > 0)
        out.append({"sample": os.path.basename(path)[:-len(".taxa.tsv")], "unseen, own reads unaligned": n,
                    "with an own record": sum(1 for s in want if records[s] > 0),
                    "own read seeded on itself": hit, "share": hit / n if n else 0.0,
                    "median own reads seeded on itself": float(pd.Series([len(reads_seeded[s]) for s in want
                                                                          if own_seeded[s] > 0]).median()) if hit else 0.0})
        print(out[-1], flush=True)
    pd.set_option("display.width", 250)
    print(pd.DataFrame(out).to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
