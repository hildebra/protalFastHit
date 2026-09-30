#!/usr/bin/env python3
"""compare.py - compact per-depth comparison of several evaluate.py labels (strain rows only).

Usage: compare.py [--bacteria] [--stage raw|qc] [--nodummya] label [label ...]
Prints one line per (label, depth): called fraction, ACGT accuracy, SNP sensitivity (all / covered),
N rate at covered true-SNP vs invariant sites, false-alt ppm, IUPAC share of calls.
"""
import sys, os
import pandas as pd
sys.path.insert(0, os.path.expanduser("~/audit5/accuracy"))
import summarize as S


def main():
    args = sys.argv[1:]
    bact = "--bacteria" in args
    nod = "--nodummya" in args
    stage = "raw"
    if "--stage" in args:
        stage = args[args.index("--stage") + 1]
        args.remove("--stage"); args.remove(stage)
    labels = [a for a in args if not a.startswith("--")]
    out = []
    for lab in labels:
        df = S.load(lab)
        df = df[(df.stage == stage) & (df.kind == "strain")]
        if bact:
            df = df[df.domain == "bacteria"]
        if nod:
            df = df[df.species != "Dummya solo"]
        t = S.table(df, ["depth"])
        t.insert(0, "label", lab)
        out.append(t)
    t = pd.concat(out)
    cols = ["label", "depth", "called", "acc", "sens", "sensCov", "N_snp", "N_inv", "fSNPppm", "iupac", "refAtSNP"]
    print(f"stage={stage} bacteria_only={bact} no_dummya={nod}")
    print(S.fmt(t[cols]))


if __name__ == "__main__":
    main()
